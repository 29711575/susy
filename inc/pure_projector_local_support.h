#ifndef PURE_PROJECTOR_LOCAL_SUPPORT_H
#define PURE_PROJECTOR_LOCAL_SUPPORT_H

#include "pure_projector_green.h"
#include "pure_projector_factor.h"
#include "skewMatUtils.h"
#include <vector>

// Exact support only: no small matrix elements are dropped. Returning false
// asks the caller to retain the general dense calculation.
inline bool pureLocalSupport(const MatType &matrix,std::vector<int> &active,bool allowFull=false) {
    active.clear();
    if(matrix.rows()!=matrix.cols()||!pureProjectorMatrixFinite(matrix))return false;
    std::vector<bool> used(matrix.rows(),false);
    for(int col=0;col<matrix.cols();++col)for(int row=0;row<matrix.rows();++row)
        if(matrix(row,col)!=DataType(0,0))used[row]=used[col]=true;
    for(int row=0;row<matrix.rows();++row)if(used[row])active.push_back(row);
    return active.size()<=16&&(allowFull||active.size()<std::size_t(matrix.rows()));
}

inline MatType pureLocalBlock(const MatType &matrix,const std::vector<int> &active) {
    MatType block(active.size(),active.size());
    for(int col=0;col<block.cols();++col)for(int row=0;row<block.rows();++row)
        block(row,col)=matrix(active[row],active[col]);
    return block;
}

// An explicitly supported block has the same pivot cutoff as the old
// ambient PureFactorSolve, including the untouched exterior diagonal.
class PureAmbientBlockSolve {
    Eigen::FullPivLU<MatType> lu_;
    bool invertible_=false;
public:
    PureAmbientBlockSolve(const MatType &block,DataType exterior,int dimension):lu_(block) {
        if(block.rows()==0||!pureProjectorMatrixFinite(block)||lu_.maxPivot()==0)return;
        const double maximum=std::max(std::abs(exterior),lu_.maxPivot());
        lu_.setThreshold(Eigen::NumTraits<double>::epsilon()*dimension*maximum/lu_.maxPivot());
        invertible_=std::abs(exterior)>Eigen::NumTraits<double>::epsilon()*dimension*maximum&&lu_.isInvertible();
    }
    bool isInvertible()const{return invertible_;}
    MatType solve(const MatType &rhs)const{return lu_.solve(rhs);}
};

struct PureCompactTransform {
    bool available=false,ok=false;
    std::vector<int> active;
    MatType transform;
};

inline PureCompactTransform pureCompactTransform(const MatType &oldFactor,const MatType &newFactor) {
    PureCompactTransform result;
    if(oldFactor.rows()!=oldFactor.cols()||newFactor.rows()!=oldFactor.rows()||
       newFactor.cols()!=oldFactor.cols())return result;
    const int dimension=oldFactor.rows();
    std::vector<bool> used(dimension,false);
    for(int col=0;col<dimension;++col)for(int row=0;row<dimension;++row){
        const DataType exterior=row==col?DataType(1,0):DataType(0,0);
        if(oldFactor(row,col)!=exterior||newFactor(row,col)!=exterior)used[row]=used[col]=true;
    }
    for(int row=0;row<dimension;++row)if(used[row])result.active.push_back(row);
    // Retain the original ambient fallback when PureFactorSolve would use it.
    if(result.active.size()>16||result.active.size()>=std::size_t(dimension/2))return result;
    result.available=true;
    if(result.active.empty()){result.transform.resize(0,0);result.ok=true;return result;}
    PureAmbientBlockSolve inverse(pureLocalBlock(oldFactor,result.active).transpose(),DataType(1,0),dimension);
    if(!inverse.isInvertible())return result;
    result.transform=inverse.solve(pureLocalBlock(newFactor,result.active).transpose()).transpose();
    result.ok=pureProjectorMatrixFinite(result.transform);
    return result;
}

inline PureCompactTransform pureCompactTransform(const PureFactorMatrix &oldFactor,
                                                 const PureFactorMatrix &newFactor) {
    PureCompactTransform result;
    if(oldFactor.rows()!=newFactor.rows()||!oldFactor.hasLocalSupport()||!newFactor.hasLocalSupport())return result;
    const auto &a=oldFactor.localSupport(),&b=newFactor.localSupport();
    std::set_union(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(result.active));
    const int dimension=oldFactor.rows();
    if(result.active.size()>16||result.active.size()>=std::size_t(dimension/2))return result;
    result.available=true;
    if(result.active.empty()){result.transform.resize(0,0);result.ok=true;return result;}
    // Extract the union, including identity entries and their original signed
    // zeros. The solve and its ambient cutoff are identical to the dense scan.
    MatType oldBlock(result.active.size(),result.active.size()),newBlock(result.active.size(),result.active.size());
    for(int c=0;c<oldBlock.cols();++c)for(int r=0;r<oldBlock.rows();++r){
        oldBlock(r,c)=oldFactor(result.active[r],result.active[c]);
        newBlock(r,c)=newFactor(result.active[r],result.active[c]);
    }
    PureAmbientBlockSolve inverse(oldBlock.transpose(),DataType(1,0),dimension);
    if(!inverse.isInvertible())return result;
    result.transform=inverse.solve(newBlock.transpose()).transpose();
    result.ok=pureProjectorMatrixFinite(result.transform);
    return result;
}

// Wick contractions are transpose-antisymmetric. An oblique projector can
// acquire a symmetric component during long propagation even when its solve
// and involution residuals pass. PFAPACK reads only the upper triangle, so
// passing that component through silently changes the contraction and phase.
// Project both arguments explicitly. A material correction is usable only
// when it preserves det(I+C*G), the squared overlap factor; otherwise the
// caller must rebuild or stop. This gate does not assume a positive weight.
inline PfaffianResult pureSkewProductSign(const MatType &cayley,const MatType &green) {
    PfaffianResult failure;
    const int dimension=cayley.rows();
    if(dimension<=0||(dimension%2)!=0||cayley.cols()!=dimension||
       green.rows()!=dimension||green.cols()!=dimension||
       !pureProjectorMatrixFinite(cayley)||!pureProjectorMatrixFinite(green))return failure;
    const MatType c=(0.5*(cayley-cayley.transpose())).eval();
    const MatType g=(0.5*(green-green.transpose())).eval();
    if(!pureProjectorMatrixFinite(c)||!pureProjectorMatrixFinite(g)){
        failure.status=PfaffianStatus::nonfinite_pivot;return failure;
    }
    const double correction=std::max(
        (c-cayley).norm()/std::max(1.,cayley.norm()),
        (g-green).norm()/std::max(1.,green.norm()));
    if(correction>1e-10){
        const MatType identity=MatType::Identity(dimension,dimension);
        const MatType original=identity+cayley*green,skew=identity+c*g;
        if(!pureProjectorMatrixFinite(original)||!pureProjectorMatrixFinite(skew)){
            failure.status=PfaffianStatus::nonfinite_pivot;return failure;
        }
        const DataType before=original.fullPivLu().determinant();
        const DataType after=skew.fullPivLu().determinant();
        const double scale=std::max(std::abs(before),std::abs(after));
        if(!std::isfinite(scale)||scale==0||
           !std::isfinite(before.real())||!std::isfinite(before.imag())||
           !std::isfinite(after.real())||!std::isfinite(after.imag())||
           std::abs(before-after)>1e-8*scale){
            failure.status=PfaffianStatus::untrusted_condition;return failure;
        }
    }
    return pfaffianForSignOfProductWithStatus(c,g);
}

inline PfaffianResult pureCompactProductSign(const MatType &cayley,const MatType &green,
                                           const std::vector<int> &support) {
    std::vector<int> active;
    if(!pureLocalSupport(cayley,active,true)||(active.size()%2)!=0)return PfaffianResult();
    PfaffianResult result;
    if(active.empty()){
        result.status=PfaffianStatus::success;result.value=DataType(1,0);result.min_pivot=1;
    }else{
        std::vector<int> global;for(int index:active)global.push_back(support[index]);
        result=pureSkewProductSign(pureLocalBlock(cayley,active),pureLocalBlock(green,global));
    }
    if(((green.rows()-int(active.size()))/2)%2)result.value=-result.value;
    return result;
}

// For C supported on A, the exterior rows of [C,-I;I,G] can be eliminated
// exactly. In the original ordering their Pfaffian contributes
// (-1)^((D-|A|)/2). Keep that factor so the existing physical-sign convention
// in pureFastLocalRatioAtCut is unchanged.
inline PfaffianResult pureLocalProductSign(const MatType &cayley,const MatType &green,
                                         bool &usedLocalSupport) {
    usedLocalSupport=false;
    std::vector<int> active;
    if(cayley.rows()!=cayley.cols()||cayley.rows()!=green.rows()||
       green.rows()!=green.cols()||(green.rows()%2)!=0||
       !pureProjectorMatrixFinite(cayley)||!pureProjectorMatrixFinite(green))
        return PfaffianResult();
    if(!pureLocalSupport(cayley,active)||(active.size()%2)!=0)
        return pureSkewProductSign(cayley,green);
    PfaffianResult result;
    if(active.empty()){
        result.status=PfaffianStatus::success;result.value=DataType(1,0);
        result.min_pivot=1;
    }else{
        result=pureSkewProductSign(
            pureLocalBlock(cayley,active),pureLocalBlock(green,active));
    }
    if(((green.rows()-int(active.size()))/2)%2)result.value=-result.value;
    usedLocalSupport=true;
    return result;
}

struct PureLocalGreenResult {
    bool ok=false;
    DataType determinant=DataType(0,0);
    double solve_residual=std::numeric_limits<double>::infinity();
    MatType green;
    MatType update_left,update_right;
};

// P=(I-G)/2, T=I+E*d*E^T, K=I+d*P_AA.
// det(overlap_new)/det(overlap_old)=det(K), and
// G_new=G-2*(I-P)*E*K^{-1}*d*E^T*P.
// The caller supplies a Green already certified at the current cut, and keeps
// the physical-sign, reality, local solve and updated-Green trust gates.
inline PureLocalGreenResult pureLocalGreenUpdateFromBlock(
    const MatType &d,const std::vector<int> &active,const MatType &green,double residualTolerance,
    bool stableDeterminant=false,bool prepareOnly=false) {
    PureLocalGreenResult result;
    if(green.rows()!=green.cols()||!pureProjectorMatrixFinite(green))return result;
    const int dimension=green.rows();
    if(active.empty()){
        result.determinant=DataType(1,0);if(!prepareOnly)result.green=green;
        result.update_left.resize(dimension,0);result.update_right.resize(0,dimension);
        result.solve_residual=0;result.ok=true;return result;
    }
    const int rank=int(active.size());
    const MatType paa=(MatType::Identity(rank,rank)-pureLocalBlock(green,active))*0.5;
    const MatType k=MatType::Identity(rank,rank)+d*paa;
    Eigen::FullPivLU<MatType> solve(k);
    if(!solve.isInvertible())return result;
    MatType pRows(rank,dimension),qColumns(dimension,rank);
    for(int i=0;i<rank;++i){
        pRows.row(i)=-0.5*green.row(active[i]);pRows(i,active[i])+=0.5;
        qColumns.col(i)=0.5*green.col(active[i]);qColumns(active[i],i)+=0.5;
    }
    const MatType rhs=d*pRows;
    const MatType z=solve.solve(rhs);
    result.solve_residual=(k*z-rhs).norm()/std::max(1.0,rhs.norm());
    if(!pureProjectorMatrixFinite(z)||!std::isfinite(result.solve_residual)||
       result.solve_residual>residualTolerance)return result;
    // The expanded 4 x 4 determinant can lose a small proposal through
    // cancellation. Reuse the already computed pivoted LU in production.
    result.determinant=stableDeterminant?solve.determinant():k.determinant();
    result.update_left=-2.0*qColumns;result.update_right=z;
    if(!prepareOnly)result.green=green+result.update_left*result.update_right;
    result.ok=std::isfinite(result.determinant.real())&&
        std::isfinite(result.determinant.imag())&&(prepareOnly||pureProjectorMatrixFinite(result.green));
    return result;
}

inline PureLocalGreenResult pureLocalGreenUpdate(
    const MatType &transform,const MatType &green,double residualTolerance,
    bool stableDeterminant=false,bool prepareOnly=false) {
    if(transform.rows()!=transform.cols()||transform.rows()!=green.rows())return PureLocalGreenResult();
    const MatType delta=transform-MatType::Identity(transform.rows(),transform.cols());
    std::vector<int> active;
    if(!pureLocalSupport(delta,active,true))return PureLocalGreenResult();
    return pureLocalGreenUpdateFromBlock(pureLocalBlock(delta,active),active,green,
                                        residualTolerance,stableDeterminant,prepareOnly);
}

inline void pureMaterializeLocalGreen(const MatType &green,PureLocalGreenResult &local) {
    if(local.update_left.cols()==0)local.green=green;
    else local.green=green+local.update_left*local.update_right;
    local.ok=local.ok&&pureProjectorMatrixFinite(local.green);
}

// The production gate on a candidate G is finiteness (purity is checked at
// stabilization). Certify that gate without forming G+UV for a rejection.
// Ambiguous/extreme inputs still materialize and take the original gate.
inline bool pureLocalGreenFiniteWithoutMaterializing(const MatType &green,const PureLocalGreenResult &local) {
    if(local.update_left.cols()==0)return true;
    auto maximum=[](const MatType &m){return std::max(m.real().cwiseAbs().maxCoeff(),m.imag().cwiseAbs().maxCoeff());};
    const long double g=maximum(green),u=maximum(local.update_left),v=maximum(local.update_right);
    const long double bound=g+2.L*local.update_left.cols()*u*v;
    return std::isfinite(bound)&&bound<std::numeric_limits<double>::max()/16.L;
}

// For G'=G+UV, G'^2-I=(G^2-I)+(GU+UVU)V+U(VG).
// The cached purity norm bounds the first term. The remaining terms cost
// O(D^2 r), with a conservative complex-arithmetic roundoff allowance for the
// cached dense product, the low-rank products and the stored G' addition.
// This is an upper bound, not a relaxed tolerance. If it cannot certify the
// existing gate, the caller evaluates the original dense residual.
inline double pureLocalGreenResidualBound(const PureProjectorGreenResult &trusted,
                                         const PureLocalGreenResult &local) {
    const auto &g=trusted.green;const auto &u=local.update_left;const auto &v=local.update_right;
    if(!local.ok||!trusted.ok()||!std::isfinite(trusted.green_residual)||
       trusted.green_residual<0||u.rows()!=g.rows()||v.cols()!=g.cols()||u.cols()!=v.rows())
        return std::numeric_limits<double>::infinity();
    const double oldNorm=trusted.green_residual*std::max(1.,std::sqrt(double(g.rows())));
    const MatType correction=(g*u+u*(v*u))*v+u*(v*g);
    const double scale=std::max(1.,g.norm()+u.norm()*v.norm());
    const double rounding=128.*std::numeric_limits<double>::epsilon()*
        std::max<Eigen::Index>(1,g.rows())*scale*scale;
    return (oldNorm+correction.norm()+rounding)/std::max(1.,local.green.norm());
}
#endif
