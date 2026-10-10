#ifndef KITAEV_SSH_BIPARTITE_SPECTRUM_H
#define KITAEV_SSH_BIPARTITE_SPECTRUM_H
#include "skewMatUtils.h"
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <vector>

namespace kitaev_ssh {
// K=[0,-iT; iT^T,0], T real, and T=U Sigma V^T. Hence exp(-a K)
// has real diagonal blocks U cosh(a Sigma) U^T / V cosh(a Sigma) V^T
// and imaginary cross block +i U sinh(a Sigma) V^T. No eigenmode,
// bond or matrix coefficient is truncated. Geometry/gauge enter only T.
class BipartiteSpectrum {
    Eigen::MatrixXd input_,u_,vt_,v_,scaled_,aa_,bb_,ab_;
    dVecType singular_,function_;
    std::vector<double> work_;
    std::vector<lapack_int> integer_work_;
    lapack_int workspace_dimension_=0;
    bool valid_=false;
public:
    static bool disabled() {
        static const bool flag=[] {
            const char *e=std::getenv("PFQMC_NO_BIPARTITE_SPECTRUM");
            return e&&*e&&*e!='0';
        }();return flag;
    }
    bool compute(const MatType &k) {
        valid_=false;
        if(k.rows()<=0||k.rows()%2||k.cols()!=k.rows()||!k.allFinite())return false;
        const lapack_int n=k.rows()/2;
        // Exact structure test, not a small-entry cutoff. A generic generator
        // with any extra coupling or real coefficient retains the old solver.
        for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
            if(k(i,j)!=DataType(0)||k(n+i,n+j)!=DataType(0)||
               k(i,n+j).real()!=0.||k(n+j,i)!=std::conj(k(i,n+j)))return false;
        }
        input_=(-k.topRightCorner(n,n).imag()).eval();
        singular_.resize(n);u_.resize(n,n);vt_.resize(n,n);
        if(workspace_dimension_!=n||work_.empty()) {
            integer_work_.resize(8*n);
            double query=0.;
            const lapack_int info=LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,'S',n,n,
                input_.data(),n,singular_.data(),u_.data(),n,vt_.data(),n,
                &query,-1,integer_work_.data());
            if(info!=0||!std::isfinite(query)||query<1.||
               query>double(std::numeric_limits<lapack_int>::max()))return false;
            work_.resize(static_cast<lapack_int>(query));workspace_dimension_=n;
        }
        const lapack_int info=LAPACKE_dgesdd_work(LAPACK_COL_MAJOR,'S',n,n,
            input_.data(),n,singular_.data(),u_.data(),n,vt_.data(),n,
            work_.data(),static_cast<lapack_int>(work_.size()),integer_work_.data());
        if(info!=0||!singular_.allFinite()||!u_.allFinite()||!vt_.allFinite()||
           (singular_.array()<0.).any())return false;
        v_=vt_.transpose();valid_=true;return true;
    }
    bool valid()const{return valid_;}
    void exponential(double a,MatType &out) {
        if(!valid_)throw std::logic_error("invalid bipartite spectrum");
        const int n=singular_.size();
        function_=(a*singular_.array()).cosh();
        scaled_=u_*function_.asDiagonal();aa_.noalias()=scaled_*u_.transpose();
        scaled_=v_*function_.asDiagonal();bb_.noalias()=scaled_*vt_;
        function_=(a*singular_.array()).sinh();
        scaled_=u_*function_.asDiagonal();ab_.noalias()=scaled_*vt_;
        out.resize(2*n,2*n);
        for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
            out(i,j)=DataType(aa_(i,j),0.);out(n+i,n+j)=DataType(bb_(i,j),0.);
            out(i,n+j)=DataType(0.,ab_(i,j));out(n+j,i)=DataType(0.,-ab_(i,j));
        }
    }
    void cayley(double a,MatType &out) {
        if(!valid_)throw std::logic_error("invalid bipartite spectrum");
        const int n=singular_.size();
        function_=(.5*a*singular_.array()).tanh();
        scaled_=u_*function_.asDiagonal();ab_.noalias()=scaled_*vt_;
        out.setZero(2*n,2*n);
        for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
            out(i,n+j)=DataType(0.,-ab_(i,j));out(n+j,i)=DataType(0.,ab_(i,j));
        }
    }
    bool cayleySafe(double a,int dimension)const {
        if(!valid_)return false;
        const double maximum=a*singular_(0);
        const double hi=1.+std::exp(maximum),lo=1.+std::exp(-maximum);
        return std::isfinite(hi)&&lo/hi>
            2.*double(dimension)*std::numeric_limits<double>::epsilon();
    }
};

// Reference callers keep the original full complex spectral construction.
// The accelerated workspace remains owned by the chain, not by this view.
class PhononNodeSpectrum {
    HermitianSpectral reference_;
    BipartiteSpectrum *bipartite_=nullptr;
public:
    bool compute(const MatType &k,bool eligible,BipartiteSpectrum &workspace) {
        bipartite_=nullptr;
        if(eligible&&k.rows()/2>=2&&!BipartiteSpectrum::disabled()&&workspace.compute(k)) {
            bipartite_=&workspace;return true;
        }
        reference_=selfAdjointSpectral(k);return reference_.ok;
    }
    bool reduced()const{return bipartite_!=nullptr;}
    void exponential(double a,MatType &out) {
        if(bipartite_)bipartite_->exponential(a,out);
        else {
            cVecType f=(-a*reference_.D.array()).exp().template cast<DataType>();
            out.noalias()=reference_.V*f.asDiagonal()*reference_.V.adjoint();
        }
    }
    void cayley(double a,MatType &out) {
        if(bipartite_)bipartite_->cayley(a,out);
        else {
            cVecType f=(.5*a*reference_.D.array()).tanh().template cast<DataType>();
            out.noalias()=reference_.V*f.asDiagonal()*reference_.V.adjoint();
        }
    }
    bool cayleySafe(double a,int dimension)const {
        if(bipartite_)return bipartite_->cayleySafe(a,dimension);
        double lo=std::numeric_limits<double>::infinity(),hi=0.;
        for(int e=0;e<reference_.D.size();++e) {
            const double onePlus=1.+std::exp(-a*reference_.D(e));
            if(std::isfinite(onePlus)){lo=std::min(lo,onePlus);hi=std::max(hi,onePlus);}
        }
        const double rcond=hi>0.&&std::isfinite(lo)?lo/hi:0.;
        return rcond>2.*double(dimension)*std::numeric_limits<double>::epsilon();
    }
};
}
#endif
