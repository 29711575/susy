#ifndef PURE_PROJECTOR_FACTOR_SOLVE_H
#define PURE_PROJECTOR_FACTOR_SOLVE_H
#include "types.h"
#include "pure_projector_factor.h"
#include <memory>
#include <vector>

// Apply I plus an exactly supported block without an ambient dense product.
// Kinetic/general matrices retain the dense operation; no entries are cut off.
inline MatType pureFactorApply(const MatType &factor,const MatType &rhs) {
    std::vector<int> active;
    std::vector<bool> used(factor.rows(),false);
    for(int c=0;c<factor.cols();++c)for(int r=0;r<factor.rows();++r)
        if(factor(r,c)!=(r==c?DataType(1,0):DataType(0,0)))used[r]=used[c]=true;
    for(int r=0;r<factor.rows();++r)if(used[r])active.push_back(r);
    if(active.size()>16||active.size()>=std::size_t(factor.rows()))return factor*rhs;
    MatType result=rhs;
    if(active.empty())return result;
    MatType block(active.size(),active.size()),rows(active.size(),rhs.cols());
    for(int r=0;r<int(active.size());++r){
        rows.row(r)=rhs.row(active[r]);
        for(int c=0;c<int(active.size());++c)block(r,c)=factor(active[r],active[c]);
    }
    const MatType product=block*rows;
    for(int r=0;r<int(active.size());++r)result.row(active[r])=product.row(r);
    return result;
}

inline void purePairFactorApplyInto(const PureFactorView &factor,const MatType &rhs,MatType &result){
    result.resize(rhs.rows(),rhs.cols());
    for(const auto &pair:factor.pairBlocks()){
        const MatType block=factor.pairMatrix(pair);
        if(pair.second<0)result.row(pair.first)=block(0,0)*rhs.row(pair.first);
        else{
            result.row(pair.first)=block(0,0)*rhs.row(pair.first)+block(0,1)*rhs.row(pair.second);
            result.row(pair.second)=block(1,0)*rhs.row(pair.first)+block(1,1)*rhs.row(pair.second);
        }
    }
}

inline MatType pureFactorApply(const PureFactorView &factor,const MatType &rhs) {
    if(factor.hasPairBlocks()){
        MatType result;purePairFactorApplyInto(factor,rhs,result);return result;
    }
    if(!factor.hasLocalSupport()||factor.active().size()>=std::size_t(factor.dimension())){
        if(const MatType *dense=factor.denseStorage())return (*dense)*rhs;
        return pureFactorApply(factor.dense(),rhs);
    }
    const auto &active=factor.active();
    MatType result=rhs;
    if(active.empty())return result;
    MatType rows(active.size(),rhs.cols());
    for(int r=0;r<int(active.size());++r)rows.row(r)=rhs.row(active[r]);
    const MatType product=factor.block()*rows;
    for(int r=0;r<int(active.size());++r)result.row(active[r])=product.row(r);
    return result;
}

struct PureFactorApplyWorkspace { MatType rows,product; };
struct PureFactorSolveWorkspace { MatType rows,solved; };

// The destination and workspace must not alias rhs. Preserve the original
// product and row-scatter order, while retaining allocations between calls.
inline void pureFactorApplyInto(const PureFactorView &factor,const MatType &rhs,
                                MatType &result,PureFactorApplyWorkspace &work) {
    if(factor.hasPairBlocks()){purePairFactorApplyInto(factor,rhs,result);return;}
    if(!factor.hasLocalSupport()||factor.active().size()>=std::size_t(factor.dimension())){
        if(const MatType *dense=factor.denseStorage())result=(*dense)*rhs;
        else result=pureFactorApply(factor.dense(),rhs);
        return;
    }
    const auto &active=factor.active();
    result=rhs;
    if(active.empty())return;
    work.rows.resize(active.size(),rhs.cols());
    for(int r=0;r<int(active.size());++r)work.rows.row(r)=rhs.row(active[r]);
    work.product.noalias()=factor.block()*work.rows;
    for(int r=0;r<int(active.size());++r)result.row(active[r])=work.product.row(r);
}

// Reconstruction can update an identity-plus-local-block frame in place:
// all active input rows are packed before scattering any output rows.
inline bool pureFactorApplyInPlace(const PureFactorView &factor,MatType &matrix,
                                  PureFactorApplyWorkspace &work){
    if(!factor.hasLocalSupport()||factor.active().size()>=std::size_t(factor.dimension()))return false;
    const auto &active=factor.active();
    if(active.empty())return true;
    work.rows.resize(active.size(),matrix.cols());
    for(int r=0;r<int(active.size());++r)work.rows.row(r)=matrix.row(active[r]);
    work.product.noalias()=factor.block()*work.rows;
    for(int r=0;r<int(active.size());++r)matrix.row(active[r])=work.product.row(r);
    return true;
}

// Exact local support leaves all other columns unchanged. Pack every input
// column before writing any output, making this operation alias-safe.
inline bool pureFactorApplyTransposedInPlace(const PureFactorView &factor,MatType &matrix,
                                            PureFactorApplyWorkspace &work){
    if(!factor.hasLocalSupport()||factor.active().size()>=std::size_t(factor.dimension()))return false;
    const auto &active=factor.active();
    if(active.empty())return true;
    work.rows.resize(active.size(),matrix.rows());
    for(int r=0;r<int(active.size());++r)work.rows.row(r)=matrix.col(active[r]).transpose();
    work.product.noalias()=factor.block()*work.rows;
    for(int r=0;r<int(active.size());++r)matrix.col(active[r])=work.product.row(r).transpose();
    return true;
}

// Solve a block diagonal matrix A_active + d I_inactive. Support is detected
// by exact comparisons, not numerical truncation. The full-dimensional LU
// pivot threshold is retained, including the exterior diagonal's pivot scale.
// General/dense factors use the original FullPivLU implementation.
class PureFactorSolve {
    struct PairSolve {
        int first,second;
        std::unique_ptr<Eigen::FullPivLU<MatType>> lu;
    };
    std::vector<PairSolve> pairs_;
    bool paired_=false;
    Eigen::Index dimension_=0;
    DataType exterior_;
    std::vector<int> active_;
    std::unique_ptr<Eigen::FullPivLU<MatType>> lu_;
    bool compact_=false,invertible_=false;
public:
    PureFactorSolve(const PureFactorView &factor,DataType exterior)
        :dimension_(factor.dimension()),exterior_(exterior){
        if(factor.hasPairBlocks()&&exterior==DataType(1,0)){
            paired_=true;double maximum=0.;
            for(const auto &pair:factor.pairBlocks()){
                PairSolve entry{pair.first,pair.second,
                    std::unique_ptr<Eigen::FullPivLU<MatType>>(new Eigen::FullPivLU<MatType>(factor.pairMatrix(pair)))};
                maximum=std::max(maximum,entry.lu->maxPivot());pairs_.push_back(std::move(entry));
            }
            if(maximum==0.||!std::isfinite(maximum))return;
            // Match the ambient FullPivLU pivot gate, not a two-dimensional
            // rank threshold for each block considered in isolation.
            for(auto &pair:pairs_){
                const double pivot=pair.lu->maxPivot();
                if(pivot==0.||!std::isfinite(pivot))return;
                pair.lu->setThreshold(Eigen::NumTraits<double>::epsilon()*dimension_*(maximum/pivot));
                if(!pair.lu->isInvertible())return;
            }
            invertible_=true;return;
        }
        if(!factor.hasLocalSupport()||exterior!=DataType(1,0)||
           factor.active().size()>=std::size_t(dimension_/2)){
            if(const MatType *dense=factor.denseStorage())
                *this=PureFactorSolve(*dense,exterior);
            else *this=PureFactorSolve(factor.dense(),exterior);
            return;
        }
        active_=factor.active();compact_=true;
        if(active_.empty()){invertible_=true;return;}
        lu_.reset(new Eigen::FullPivLU<MatType>(factor.block()));
        const double maximum=std::max(std::abs(exterior_),lu_->maxPivot());
        if(lu_->maxPivot()==0)return;
        lu_->setThreshold(Eigen::NumTraits<double>::epsilon()*dimension_*maximum/lu_->maxPivot());
        invertible_=std::abs(exterior_)>Eigen::NumTraits<double>::epsilon()*dimension_*maximum&&lu_->isInvertible();
    }
    PureFactorSolve(const MatType &matrix,DataType exterior,bool allowCompact=true)
        :dimension_(matrix.rows()),exterior_(exterior){
        if(matrix.rows()!=matrix.cols()||matrix.rows()==0)return;
        if(allowCompact){
            std::vector<bool> used(dimension_,false);
            for(int col=0;col<dimension_;++col)for(int row=0;row<dimension_;++row)
                if(matrix(row,col)!=(row==col?exterior_:DataType(0,0)))used[row]=used[col]=true;
            for(int row=0;row<dimension_;++row)if(used[row])active_.push_back(row);
            compact_=active_.size()<std::size_t(dimension_/2)&&active_.size()<=16&&exterior_!=DataType(0,0);
        }
        if(!compact_){lu_.reset(new Eigen::FullPivLU<MatType>(matrix));invertible_=lu_->isInvertible();return;}
        if(active_.empty()){invertible_=std::isfinite(exterior_.real())&&std::isfinite(exterior_.imag());return;}
        MatType block(active_.size(),active_.size());
        for(int c=0;c<block.cols();++c)for(int r=0;r<block.rows();++r)block(r,c)=matrix(active_[r],active_[c]);
        lu_.reset(new Eigen::FullPivLU<MatType>(block));
        const double maximum=std::max(std::abs(exterior_),lu_->maxPivot());
        if(lu_->maxPivot()==0)return;
        lu_->setThreshold(Eigen::NumTraits<double>::epsilon()*dimension_*maximum/lu_->maxPivot());
        invertible_=std::abs(exterior_)>Eigen::NumTraits<double>::epsilon()*dimension_*maximum&&lu_->isInvertible();
    }
    bool isInvertible()const{return invertible_;}
    MatType solve(const MatType &rhs)const{
        if(paired_){MatType result;PureFactorSolveWorkspace work;solveInto(rhs,result,work);return result;}
        if(!compact_)return lu_->solve(rhs);
        MatType result=rhs/exterior_;
        if(active_.empty())return result;
        MatType small(active_.size(),rhs.cols());
        for(int row=0;row<small.rows();++row)small.row(row)=rhs.row(active_[row]);
        const MatType solved=lu_->solve(small);
        for(int row=0;row<small.rows();++row)result.row(active_[row])=solved.row(row);
        return result;
    }
    // Only an identity exterior can leave the untouched columns in place.
    // The local pivoted solve and full-dimensional rank gate are unchanged.
    bool solveTransposedInPlace(MatType &matrix,PureFactorSolveWorkspace &work)const{
        if(!compact_||exterior_!=DataType(1,0))return false;
        if(active_.empty())return true;
        work.rows.resize(active_.size(),matrix.rows());
        for(int r=0;r<int(active_.size());++r)work.rows.row(r)=matrix.col(active_[r]).transpose();
        work.solved=lu_->solve(work.rows);
        for(int r=0;r<int(active_.size());++r)matrix.col(active_[r])=work.solved.row(r).transpose();
        return true;
    }
    void solveInto(const MatType &rhs,MatType &result,PureFactorSolveWorkspace &work)const{
        if(paired_){
            result.resize(rhs.rows(),rhs.cols());
            for(const auto &pair:pairs_){
                work.rows.resize(pair.second<0?1:2,rhs.cols());
                work.rows.row(0)=rhs.row(pair.first);
                if(pair.second>=0)work.rows.row(1)=rhs.row(pair.second);
                work.solved=pair.lu->solve(work.rows);
                result.row(pair.first)=work.solved.row(0);
                if(pair.second>=0)result.row(pair.second)=work.solved.row(1);
            }
            return;
        }
        if(!compact_){result=lu_->solve(rhs);return;}
        result=rhs/exterior_;
        if(active_.empty())return;
        work.rows.resize(active_.size(),rhs.cols());
        for(int row=0;row<work.rows.rows();++row)work.rows.row(row)=rhs.row(active_[row]);
        work.solved=lu_->solve(work.rows);
        for(int row=0;row<work.rows.rows();++row)result.row(active_[row])=work.solved.row(row);
    }
};
#endif
