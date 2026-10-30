#ifndef PURE_PROJECTOR_REAL_H
#define PURE_PROJECTOR_REAL_H
#include "gaussian_trial_state.h"
#include "pure_projector_factor.h"
#include "pure_projector_canonical.h"
#include <cstdlib>
#include <limits>
#include <vector>

// Chain-private, bounded per-slice cache for generic factors. Retain their
// immutable owner as well as the real D^dag B D representation; storage
// addresses cannot be recycled underneath an entry. Model-certified factors
// already own that representation and are borrowed only during evaluation.
class PureRealMagnitudeWorkspace {
    struct Factor {PureFactorMatrix owner;Eigen::MatrixXd value;bool valid=false;};
    std::vector<Factor> factors_;
    std::vector<double> tau_,qr_work_,q_work_,svd_work_;
    Eigen::MatrixXd svd_input_,metric_q_,isotropy_,delta_,gram_;
    lapack_int m_=0,n_=0;
    bool queryWorkspace(int m,int n) {
        if(m_==m&&n_==n&&!qr_work_.empty())return true;
        tau_.resize(n);double query=0.;
        if(LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,m,n,next_q.data(),m,
                              tau_.data(),&query,-1)!=0||!sizeWork(query,qr_work_))return false;
        if(LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,m,n,n,next_q.data(),m,
                              tau_.data(),&query,-1)!=0||!sizeWork(query,q_work_))return false;
        double dummy=0.;singular.resize(n);svd_input_.resize(n,n);
        if(LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'N','N',n,n,svd_input_.data(),n,
                              singular.data(),&dummy,1,&dummy,1,&query,-1)!=0||
           !sizeWork(query,svd_work_))return false;
        m_=m;n_=n;return true;
    }
    static bool sizeWork(double count,std::vector<double> &buffer) {
        if(!std::isfinite(count)||count<1.||count>double(std::numeric_limits<lapack_int>::max()))return false;
        buffer.resize(static_cast<lapack_int>(count));return true;
    }
public:
    Eigen::MatrixXd q,next_q,overlap;
    dVecType singular;
    unsigned long long calls=0,fallbacks=0,cache_hits=0,conversions=0,certified_hits=0;
    static bool disabled() {
        static const bool value=[] {
            const char *e=std::getenv("PFQMC_NO_REAL_MAGNITUDE");return e&&*e&&*e!='0';
        }();return value;
    }
    bool trialMatches(const GaussianTrialState &trial)const {
        const auto &t=trial.real_magnitude_trial;
        return t&&t->source.rows()==trial.Phi.rows()&&t->source.cols()==trial.Phi.cols()&&
            (t->source.array()==trial.Phi.array()).all();
    }
    // A complete candidate rebuild no longer needs the previous generic
    // owners. Keep their real buffers allocated for subsequent conversions.
    void releaseFactorOwners() {
        for(auto &entry:factors_) {entry.owner=PureFactorMatrix{};entry.valid=false;}
    }
    const Eigen::MatrixXd *factor(const PureFactorMatrix &matrix,int index,int count) {
        if(const auto *certified=matrix.certifiedRealStorage()) {
            // The slice owns the immutable factor throughout this weight
            // evaluation. Do not retain proposal owners after a rejection.
            if(index<int(factors_.size())) {
                auto &entry=factors_[index];entry.owner=PureFactorMatrix{};entry.valid=false;
            }
            ++certified_hits;return certified;
        }
        factors_.resize(count);
        auto &entry=factors_[index];
        if(entry.owner.sharesStorageWith(matrix)) {
            ++cache_hits;return entry.valid?&entry.value:nullptr;
        }
        entry.owner=matrix;entry.valid=false;
        const MatType *b=matrix.denseStorage();
        if(!b||b->rows()%2||!b->allFinite())return nullptr;
        const int n=b->rows()/2;
        if(b->cols()!=2*n)return nullptr;
        // Exact structure checks. Even a tiny additional complex coupling
        // remains in the complex reference; no Hamiltonian entry is dropped.
        for(int j=0;j<n;++j)for(int i=0;i<n;++i)
            if((*b)(i,j).imag()!=0.||(*b)(n+i,n+j).imag()!=0.||
               (*b)(i,n+j).real()!=0.||(*b)(n+i,j).real()!=0.)return nullptr;
        entry.value.resize(2*n,2*n);
        for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
            entry.value(i,j)=(*b)(i,j).real();entry.value(n+i,n+j)=(*b)(n+i,n+j).real();
            entry.value(i,n+j)=-(*b)(i,n+j).imag();entry.value(n+i,j)=(*b)(n+i,j).imag();
        }
        entry.valid=true;++conversions;return &entry.value;
    }
    bool orthogonalize(double &detR,double orthTolerance) {
        const int m=next_q.rows(),n=next_q.cols();
        if(!next_q.allFinite()||!queryWorkspace(m,n))return false;
        const auto rawQr=[&](double &scale) {
            if(LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,m,n,next_q.data(),m,tau_.data(),
                                  qr_work_.data(),qr_work_.size())!=0)return false;
            scale=1.;for(int i=0;i<n;++i)scale*=next_q(i,i);
            return LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,m,n,n,next_q.data(),m,tau_.data(),
                                      q_work_.data(),q_work_.size())==0&&next_q.allFinite();
        };
        if(!rawQr(detR))return false;
        auto &stats=pureCanonicalCorrectionStats;++stats.real_calls;
        const auto fail=[&](){++stats.failures;return false;};
        if(m!=2*n)return fail();
        const double normalization=std::sqrt(double(n));
        gram_.noalias()=next_q.transpose()*next_q;
        gram_-=Eigen::MatrixXd::Identity(n,n);
        const double orth=gram_.norm()/normalization;
        if(!std::isfinite(orth)||orth>orthTolerance)return fail();
        // Q_complex=D Q_real, D=diag(I,iI). Hence Q_complex^T Q_complex
        // is Q_real^T J Q_real, J=diag(I,-I), and D^dag conj(D)=J.
        metric_q_=next_q;metric_q_.bottomRows(n)*=-1.;
        isotropy_.noalias()=next_q.transpose()*metric_q_;
        delta_.noalias()=-.5*metric_q_*isotropy_;
        const double relative=delta_.norm()/next_q.norm(),before=isotropy_.norm()/normalization;
        if(std::isfinite(relative))stats.maximum_relative_correction=std::max(stats.maximum_relative_correction,relative);
        if(std::isfinite(before))stats.maximum_input_isotropy=std::max(stats.maximum_input_isotropy,before);
        if(!isotropy_.allFinite()||!delta_.allFinite()||!std::isfinite(relative)||
           relative>pureProjectorCanonicalCorrectionBound(m))return fail();
        next_q+=delta_;
        double correctionDet;
        if(!rawQr(correctionDet))return fail();
        detR*=correctionDet;
        metric_q_=next_q;metric_q_.bottomRows(n)*=-1.;
        gram_.noalias()=next_q.transpose()*metric_q_;
        const double after=gram_.norm()/normalization;
        gram_.noalias()=next_q.transpose()*next_q;
        gram_-=Eigen::MatrixXd::Identity(n,n);
        const double afterOrth=gram_.norm()/normalization;
        if(std::isfinite(after))stats.maximum_output_isotropy=std::max(stats.maximum_output_isotropy,after);
        if(std::isfinite(afterOrth))stats.maximum_output_orthonormality=std::max(stats.maximum_output_orthonormality,afterOrth);
        if(!std::isfinite(after)||!std::isfinite(afterOrth)||!std::isfinite(detR)||
           after>pureProjectorCanonicalCorrectionBound(m)||afterOrth>orthTolerance)return fail();
        return true;
    }
    bool overlapValues() {
        const int n=overlap.rows();svd_input_=overlap;double dummy=0.;
        return LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'N','N',n,n,svd_input_.data(),n,
            singular.data(),&dummy,1,&dummy,1,svd_work_.data(),svd_work_.size())==0&&singular.allFinite();
    }
};
#endif
