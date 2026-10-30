#ifndef PURE_PROJECTOR_SVD_H
#define PURE_PROJECTOR_SVD_H
#include "types.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

// Values-only overlap SVD for magnitude and signed walks. The workspace is private to
// one chain and retains LAPACK storage across time slices. It never replaces
// the signed Pfaffian/Green evaluation. Small or uncertifiable overlaps keep
// the original Eigen Jacobi SVD and exactly the original trust thresholds.
class PureOverlapSingularValues {
    Eigen::JacobiSVD<MatType> eigen_;
    MatType input_;
    dVecType singular_;
    std::vector<DataType> work_;
    std::vector<double> real_work_;
    lapack_int rows_=0,cols_=0;
    bool accelerated_=false;
    unsigned long long accelerated_calls_=0,reference_calls_=0;

    void reference(const MatType &matrix) {
        eigen_.compute(matrix,0);
        accelerated_=false;
        ++reference_calls_;
    }
public:
    static bool disabled() {
        static const bool value=[] {
            const char *e=std::getenv("PFQMC_NO_OVERLAP_SVD");
            return e&&*e&&*e!='0';
        }();
        return value;
    }
    void compute(const MatType &matrix,bool accelerate,double rankTolerance,
                 double minimumRcond,double zeroTolerance,bool checkDeterminantBoundary=true) {
        const lapack_int m=matrix.rows(),n=matrix.cols(),k=std::min(m,n);
        if(!accelerate||disabled()||k<32||!matrix.allFinite()) {
            reference(matrix);return;
        }
        input_=matrix;singular_.resize(k);
        DataType unusedU,unusedV;
        if(rows_!=m||cols_!=n||work_.empty()) {
            real_work_.resize(std::max<lapack_int>(1,5*k));
            DataType query;
            const lapack_int info=LAPACKE_zgesvd_work(LAPACK_COL_MAJOR,'N','N',m,n,
                input_.data(),m,singular_.data(),&unusedU,1,&unusedV,1,
                &query,-1,real_work_.data());
            if(info!=0||!std::isfinite(query.real())||query.real()<1||
               query.real()>double(std::numeric_limits<lapack_int>::max())) {
                reference(matrix);return;
            }
            work_.resize(static_cast<lapack_int>(query.real()));
            rows_=m;cols_=n;
        }
        const lapack_int info=LAPACKE_zgesvd_work(LAPACK_COL_MAJOR,'N','N',m,n,
            input_.data(),m,singular_.data(),&unusedU,1,&unusedV,1,
            work_.data(),static_cast<lapack_int>(work_.size()),real_work_.data());
        if(info!=0||!singular_.allFinite()||!(singular_(0)>0)) {
            reference(matrix);return;
        }
        // Preserve Eigen as the authority around rank/rcond/zero boundaries.
        // The guard selects a reference decomposition; it does not relax any
        // gate or clip a singular value. Rank-deficient candidates also fall
        // back instead of deciding their status from the accelerated solver.
        const double guard=64.*double(k)*std::numeric_limits<double>::epsilon();
        const double rcond=singular_(k-1)/singular_(0);
        double logDet=0.;
        for(int i=0;i<k;++i) {
            if(!(singular_(i)>0)) {reference(matrix);return;}
            if(checkDeterminantBoundary)logDet+=std::log(singular_(i));
        }
        if(rcond<=std::max(rankTolerance,minimumRcond)+guard||
           (checkDeterminantBoundary&&logDet<=std::log(zeroTolerance)+guard*double(k))) {
            reference(matrix);return;
        }
        accelerated_=true;++accelerated_calls_;
    }
    const dVecType &singularValues()const {
        return accelerated_?singular_:eigen_.singularValues();
    }
    int rank(double threshold) {
        const auto &s=singularValues();
        if(!accelerated_) {
            eigen_.setThreshold(threshold);
            return eigen_.rank();
        }
        if(s.size()==0)return 0;
        const double cutoff=std::max(std::numeric_limits<double>::min(),threshold*s(0));
        return (s.array()>=cutoff).count();
    }
    bool accelerated()const{return accelerated_;}
    unsigned long long acceleratedCalls()const{return accelerated_calls_;}
    unsigned long long referenceCalls()const{return reference_calls_;}
};

// Same singular-value rank/solve contract, with LAPACK for larger live Green
// reconstructions. Reference callers retain Eigen unless explicitly enabled.
class PureGreenSvd {
    std::unique_ptr<Eigen::JacobiSVD<MatType>> eigen_;
    MatType u_,v_,work_,vt_;
    dVecType singular_;
    double threshold_=1e-12;
    bool valid_=false;
    bool compute_vectors_=true;
    bool eigen_active_=false;
    bool threshold_set_=false;
public:
    PureGreenSvd()=default;
    PureGreenSvd(const MatType &matrix,bool accelerated,bool computeVectors=true) {
        compute(matrix,accelerated,computeVectors);
    }
    void compute(const MatType &matrix,bool accelerated,bool computeVectors=true) {
        compute_vectors_=computeVectors;valid_=false;eigen_active_=false;
        if(!accelerated||std::min(matrix.rows(),matrix.cols())<32){
            if(!eigen_)eigen_.reset(new Eigen::JacobiSVD<MatType>);
            eigen_->compute(matrix,compute_vectors_?(Eigen::ComputeThinU|Eigen::ComputeThinV):0);
            if(threshold_set_)eigen_->setThreshold(threshold_);
            eigen_active_=true;return;
        }
        work_=matrix;
        const lapack_int m=matrix.rows(),n=matrix.cols(),k=std::min(m,n);
        singular_.resize(k);
        if(!compute_vectors_){
            DataType unusedU,unusedV;
            const lapack_int info=LAPACKE_zgesdd(LAPACK_COL_MAJOR,'N',m,n,work_.data(),m,
                singular_.data(),&unusedU,1,&unusedV,1);
            valid_=info==0&&singular_.allFinite();
            return;
        }
        u_.resize(m,k);vt_.resize(k,n);
        const lapack_int info=LAPACKE_zgesdd(LAPACK_COL_MAJOR,'S',m,n,work_.data(),m,
            singular_.data(),u_.data(),m,vt_.data(),k);
        valid_=info==0&&singular_.allFinite()&&u_.allFinite()&&vt_.allFinite();
        if(valid_)v_=vt_.adjoint();
    }
    void setThreshold(double threshold){threshold_=threshold;threshold_set_=true;if(eigen_active_)eigen_->setThreshold(threshold);}
    int rank()const{
        if(eigen_active_)return eigen_->rank();
        if(!valid_||singular_.size()==0)return 0;
        return (singular_.array()>std::max(std::numeric_limits<double>::min(),threshold_*singular_(0))).count();
    }
    const dVecType &singularValues()const{return eigen_active_?eigen_->singularValues():singular_;}
    MatType solve(const MatType &rhs)const{
        if(!compute_vectors_)throw std::logic_error("PureGreenSvd values-only decomposition cannot solve");
        if(eigen_active_)return eigen_->solve(rhs);
        const int count=rank();
        return v_.leftCols(count)*singular_.head(count).cwiseInverse().asDiagonal()*
            (u_.leftCols(count).adjoint()*rhs);
    }
};
#endif
