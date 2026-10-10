#ifndef PURE_PROJECTOR_FACTOR_H
#define PURE_PROJECTOR_FACTOR_H
#include "types.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

class PureFactorMatrix;
struct PureFactorPairBlock {
    int first=-1,second=-1;
    MatType matrix;
};
namespace Eigen { namespace internal {
template<> struct traits<PureFactorMatrix> { typedef MatType ReturnType; };
}}

// Immutable shared, lossless storage. Eigen evaluates this expression into an
// ordinary dense matrix before numerical operations: no multiplication order,
// coefficient (including signed zero), or reference precision is changed.
class PureFactorMatrix : public Eigen::ReturnByValue<PureFactorMatrix> {
    struct Entry { Eigen::Index index; DataType value; std::uint64_t skip,sign; };
    struct Storage {
        Eigen::Index rows=0,cols=0;
        MatType dense;
        std::vector<Entry> entries;
        std::uint64_t tail=1,tail_sign=0;
        std::vector<std::uint64_t> zero_signs;
        DataType zero_fill=DataType(0,0);
        bool have_zero=false,uniform_zeros=true;
        // Exact support of A-I, computed once for this immutable value.
        // Dense/general factors retain their existing numerical path.
        bool local=false;
        std::vector<int> active;
        MatType block;
        // Exact connected components of size at most two. No small entries
        // are discarded; any additional coupling retains the dense path.
        bool paired=false;
        std::vector<PureFactorPairBlock> pairs;
    };
    std::shared_ptr<const Storage> storage_;
    static constexpr std::uint64_t prime=1099511628211ULL;
    static constexpr std::uint64_t sign_bit=1ULL<<63;
    static std::uint64_t bits(double value){std::uint64_t result;std::memcpy(&result,&value,sizeof(result));return result;}
    static std::uint64_t power(std::uint64_t n){
        std::uint64_t result=1,p=prime;
        while(n){if(n&1)result*=p;p*=p;n>>=1;}return result;
    }
public:
    PureFactorMatrix()=default;
    PureFactorMatrix(const PureFactorMatrix &other):storage_(other.storage_){}
    PureFactorMatrix(PureFactorMatrix &&other) noexcept:storage_(std::move(other.storage_)){}
    PureFactorMatrix &operator=(const PureFactorMatrix &other){storage_=other.storage_;return *this;}
    PureFactorMatrix &operator=(PureFactorMatrix &&other) noexcept {storage_=std::move(other.storage_);return *this;}
    PureFactorMatrix(const MatType &matrix){assign(matrix);}
    PureFactorMatrix &operator=(const MatType &matrix){assign(matrix);return *this;}
    // Explicit entry point for callers that have already performed their own
    // model-specific finite/Hermitian validation. Generic assignment below
    // keeps all sparse/sign/support scans for legacy density models.
    static PureFactorMatrix fromValidatedDense(const MatType &matrix) {
        if(!matrix.allFinite()) throw std::invalid_argument("validated dense factor is nonfinite");
        PureFactorMatrix result;result.assignValidatedDense(matrix);return result;
    }
    // Transfer a freshly assembled dense matrix into its immutable owner.
    // Existing lvalue callers still copy, and the finite gate is identical.
    static PureFactorMatrix fromValidatedDense(MatType &&matrix) {
        if(!matrix.allFinite()) throw std::invalid_argument("validated dense factor is nonfinite");
        PureFactorMatrix result;result.assignValidatedDense(std::move(matrix));return result;
    }
    template<class Derived> PureFactorMatrix &operator=(const Eigen::MatrixBase<Derived>&matrix){
        MatType evaluated=matrix;assign(evaluated);return *this;
    }
    void assign(const MatType &matrix){
        auto s=std::make_shared<Storage>();s->rows=matrix.rows();s->cols=matrix.cols();
        Eigen::Index previous=0;
        std::uint64_t zeroSign=0;
        const auto limit=std::size_t(matrix.size()/3);
        for(Eigen::Index i=0;i<matrix.size();++i){
            const DataType value=matrix.data()[i];
            const auto re=bits(value.real()),im=bits(value.imag());
            if((re&~sign_bit)||(im&~sign_bit)){
                s->entries.push_back({i,value,power(2*std::uint64_t(i-previous)),zeroSign});
                previous=i+1;zeroSign=0;
                if(s->entries.size()>limit){
                    s->entries.clear();s->entries.shrink_to_fit();
                    s->zero_signs.clear();s->zero_signs.shrink_to_fit();s->dense=matrix;break;
                }
            }else {
                if(!s->have_zero){s->zero_fill=value;s->have_zero=true;}
                else if(re!=bits(s->zero_fill.real())||im!=bits(s->zero_fill.imag()))s->uniform_zeros=false;
                if(re||im){
                if(s->zero_signs.empty())s->zero_signs.resize((2*matrix.size()+63)/64);
                if(re)s->zero_signs[(2*i)/64]|=1ULL<<((2*i)%64);
                if(im)s->zero_signs[(2*i+1)/64]|=1ULL<<((2*i+1)%64);
                zeroSign^=re^im;
                }
            }
        }
        s->tail=power(2*std::uint64_t(matrix.size()-previous));
        s->tail_sign=zeroSign;
        if(matrix.rows()==matrix.cols()){
            std::vector<bool> used(matrix.rows(),false);
            for(Eigen::Index c=0;c<matrix.cols();++c)for(Eigen::Index r=0;r<matrix.rows();++r)
                if(matrix(r,c)!=(r==c?DataType(1,0):DataType(0,0)))used[r]=used[c]=true;
            for(int r=0;r<matrix.rows();++r)if(used[r])s->active.push_back(r);
            s->local=s->active.size()<=16;
            if(s->local){
                s->block.resize(s->active.size(),s->active.size());
                for(int c=0;c<s->block.cols();++c)for(int r=0;r<s->block.rows();++r)
                    s->block(r,c)=matrix(s->active[r],s->active[c]);
            }else {
                s->active.clear();
                std::vector<int> partner(matrix.rows(),-1);
                bool paired=true;
                for(int c=0;c<matrix.cols()&&paired;++c)for(int r=0;r<matrix.rows();++r){
                    if(!std::isfinite(matrix(r,c).real())||!std::isfinite(matrix(r,c).imag())){paired=false;break;}
                    if(r==c||matrix(r,c)==DataType(0,0))continue;
                    if((partner[r]!=-1&&partner[r]!=c)||(partner[c]!=-1&&partner[c]!=r)){
                        paired=false;break;
                    }
                    partner[r]=c;partner[c]=r;
                }
                if(paired){
                    for(int r=0;r<matrix.rows();++r){
                        const int c=partner[r];if(c>=0&&c<r)continue;
                        PureFactorPairBlock pair;pair.first=r;pair.second=c;
                        pair.matrix.resize(c<0?1:2,c<0?1:2);pair.matrix(0,0)=matrix(r,r);
                        if(c>=0){pair.matrix(0,1)=matrix(r,c);pair.matrix(1,0)=matrix(c,r);pair.matrix(1,1)=matrix(c,c);}
                        s->pairs.push_back(std::move(pair));
                    }
                    s->paired=true;
                }
            }
        }
        storage_=std::move(s);
    }
    void assignValidatedDense(const MatType &matrix) {
        auto s=std::make_shared<Storage>();
        s->rows=matrix.rows();s->cols=matrix.cols();s->dense=matrix;
        storage_=std::move(s);
    }
    void assignValidatedDense(MatType &&matrix) {
        auto s=std::make_shared<Storage>();
        s->rows=matrix.rows();s->cols=matrix.cols();s->dense=std::move(matrix);
        storage_=std::move(s);
    }
    // Identity certifies exact immutable coefficients, including signed zeros.
    bool sharesStorageWith(const PureFactorMatrix &other)const {
        return storage_ && storage_==other.storage_;
    }
    Eigen::Index rows()const{return storage_?storage_->rows:0;}
    Eigen::Index cols()const{return storage_?storage_->cols:0;}
    // Check the same coefficients without reconstructing implicit signed zeros.
    bool allFinite()const{
        if(!storage_)return true;
        auto finite=[](DataType x){return std::isfinite(x.real())&&std::isfinite(x.imag());};
        if(storage_->dense.size()){
            for(Eigen::Index i=0;i<storage_->dense.size();++i)
                if(!finite(storage_->dense.data()[i]))return false;
        }else for(const auto &entry:storage_->entries)if(!finite(entry.value))return false;
        return true;
    }
    bool hasLocalSupport()const{return storage_&&storage_->local;}
    const std::vector<int> &localSupport()const{return storage_->active;}
    const MatType &localBlock()const{return storage_->block;}
    // Dense storage is immutable and owned by the factor. Returning a pointer
    // keeps the caller from accidentally extending its lifetime past the
    // owning factor while allowing read-only dense consumers to avoid a copy.
    const MatType *denseStorage()const{
        return storage_&&storage_->dense.size()?&storage_->dense:nullptr;
    }
    bool hasPairBlocks()const{return storage_&&storage_->paired;}
    const std::vector<PureFactorPairBlock> &pairBlocks()const{return storage_->pairs;}
    // Public proposals may be edited after a cache is attached. Verify every
    // coefficient, including signed zero, before trusting cached metadata.
    bool bitwiseMatches(const MatType &matrix)const{
        if(matrix.rows()!=rows()||matrix.cols()!=cols())return false;
        if(!storage_)return matrix.size()==0;
        if(storage_->dense.size())return std::memcmp(matrix.data(),storage_->dense.data(),matrix.size()*sizeof(DataType))==0;
        if(storage_->uniform_zeros){
            const auto re=bits(storage_->zero_fill.real()),im=bits(storage_->zero_fill.imag());
            Eigen::Index previous=0;
            for(const auto &entry:storage_->entries){
                std::uint64_t difference=0;
                for(Eigen::Index i=previous;i<entry.index;++i)
                    difference|=(bits(matrix.data()[i].real())^re)|(bits(matrix.data()[i].imag())^im);
                if(difference||bits(matrix.data()[entry.index].real())!=bits(entry.value.real())||
                   bits(matrix.data()[entry.index].imag())!=bits(entry.value.imag()))return false;
                previous=entry.index+1;
            }
            std::uint64_t difference=0;
            for(Eigen::Index i=previous;i<matrix.size();++i)
                difference|=(bits(matrix.data()[i].real())^re)|(bits(matrix.data()[i].imag())^im);
            return difference==0;
        }
        std::size_t next=0;
        for(Eigen::Index i=0;i<matrix.size();++i){
            if(next<storage_->entries.size()&&storage_->entries[next].index==i){
                const auto value=storage_->entries[next++].value;
                if(bits(matrix.data()[i].real())!=bits(value.real())||bits(matrix.data()[i].imag())!=bits(value.imag()))return false;
            }else{
                const auto &signs=storage_->zero_signs;
                const auto re=signs.empty()?0ULL:((signs[(2*i)/64]>>((2*i)%64))&1)*sign_bit;
                const auto im=signs.empty()?0ULL:((signs[(2*i+1)/64]>>((2*i+1)%64))&1)*sign_bit;
                if(bits(matrix.data()[i].real())!=re||bits(matrix.data()[i].imag())!=im)return false;
            }
        }
        return true;
    }
    DataType operator()(Eigen::Index row,Eigen::Index col)const{
        const Eigen::Index index=row+col*rows();
        if(storage_->dense.size())return storage_->dense(row,col);
        const auto &entries=storage_->entries;
        auto entry=std::lower_bound(entries.begin(),entries.end(),index,
            [](const Entry &e,Eigen::Index i){return e.index<i;});
        if(entry!=entries.end()&&entry->index==index)return entry->value;
        if(storage_->zero_signs.empty())return DataType(0,0);
        const auto &signs=storage_->zero_signs;
        return DataType((signs[(2*index)/64]>>((2*index)%64))&1?-0.0:0.0,
            (signs[(2*index+1)/64]>>((2*index+1)%64))&1?-0.0:0.0);
    }
    template<class Dest> void evalTo(Dest &destination)const{
        destination.resize(rows(),cols());
        if(!storage_)return;
        if(storage_->dense.size()){destination=storage_->dense;return;}
        // Adjoint density factors have a uniform negative imaginary zero.
        // Fill it once instead of restoring O(D^2) individual sign bits.
        if(storage_->uniform_zeros)destination.setConstant(storage_->zero_fill);
        else destination.setZero();
        for(std::size_t w=0;!storage_->uniform_zeros&&w<storage_->zero_signs.size();++w){
            std::uint64_t mask=storage_->zero_signs[w];
            while(mask){
                const unsigned bit=__builtin_ctzll(mask);
                const Eigen::Index word=Eigen::Index(w*64+bit),index=word/2;
                auto &value=destination.coeffRef(index%rows(),index/rows());
                if(word%2)value.imag(-0.0);else value.real(-0.0);
                mask&=mask-1;
            }
        }
        for(const auto &entry:storage_->entries)
            destination.coeffRef(entry.index%rows(),entry.index/rows())=entry.value;
    }
    // Exactly the existing FNV recurrence. A run of k +0 scalar words is
    // multiplication by prime^k modulo 2^64, so no hash schema changes.
    void mixHash(std::uint64_t &hash)const{
        if(!storage_)return;
        auto mix=[&](DataType v){hash^=bits(v.real());hash*=prime;hash^=bits(v.imag());hash*=prime;};
        if(storage_->dense.size()){
            for(Eigen::Index i=0;i<storage_->dense.size();++i)mix(storage_->dense.data()[i]);
        }else{
            for(const auto &entry:storage_->entries){hash*=entry.skip;hash^=entry.sign;mix(entry.value);}
            hash*=storage_->tail;hash^=storage_->tail_sign;
        }
    }
    std::size_t storageBytes()const{
        return storage_?sizeof(Storage)+sizeof(DataType)*(storage_->dense.size()+storage_->block.size())+
            sizeof(int)*storage_->active.capacity()+sizeof(Entry)*storage_->entries.capacity()+sizeof(std::uint64_t)*storage_->zero_signs.capacity()+
            sizeof(PureFactorPairBlock)*storage_->pairs.capacity()+
            [&](){std::size_t bytes=0;for(const auto&p:storage_->pairs)bytes+=p.matrix.size()*sizeof(DataType);return bytes;}():0;
    }
};

// A short-lived view; the owning factor outlives every call using it. Signed
// zeros and all matrix coefficients are preserved by the small transforms.
enum class PureFactorOperation { normal, transpose, adjoint };
class PureFactorView {
    const PureFactorMatrix &factor_;
    PureFactorOperation operation_;
    bool cached_,paired_;
    MatType transformed(const MatType &matrix)const{
        if(operation_==PureFactorOperation::transpose)return matrix.transpose();
        if(operation_==PureFactorOperation::adjoint)return matrix.adjoint();
        return matrix;
    }
public:
    PureFactorView(const PureFactorMatrix &factor,PureFactorOperation operation,bool cached,bool paired=true)
        :factor_(factor),operation_(operation),cached_(cached),paired_(paired){}
    Eigen::Index dimension()const{return factor_.rows();}
    bool hasLocalSupport()const{return cached_&&factor_.hasLocalSupport();}
    const std::vector<int> &active()const{return factor_.localSupport();}
    MatType block()const{return transformed(factor_.localBlock());}
    const MatType *denseStorage()const{
        return operation_==PureFactorOperation::normal?factor_.denseStorage():nullptr;
    }
    MatType dense()const{return transformed(MatType(factor_));}
    bool hasPairBlocks()const{return cached_&&paired_&&factor_.hasPairBlocks();}
    const std::vector<PureFactorPairBlock> &pairBlocks()const{return factor_.pairBlocks();}
    MatType pairMatrix(const PureFactorPairBlock &pair)const{return transformed(pair.matrix);}
};
#endif
