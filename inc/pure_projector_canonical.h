#ifndef PURE_PROJECTOR_CANONICAL_H
#define PURE_PROJECTOR_CANONICAL_H

#include "types.h"
#include "pure_projector_qr_workspace.h"
#include <algorithm>
#include <cmath>
#include <limits>

// Roundoff maintenance for an orthonormal Majorana Gaussian frame. This is
// deliberately bounded independently of any user residual tolerance: a
// macroscopically non-Gaussian frame is an error, never a projection request.
inline double pureProjectorCanonicalCorrectionBound(int rows) {
    return 128. * rows * std::numeric_limits<double>::epsilon();
}

struct PureCanonicalCorrectionStats {
    std::uint64_t complex_calls=0,real_calls=0,extended_calls=0,failures=0;
    std::uint64_t extended_green_residual_checks=0;
    double maximum_relative_correction=0.;
    double maximum_input_isotropy=0.,maximum_output_isotropy=0.;
    double maximum_output_orthonormality=0.;
};
inline thread_local PureCanonicalCorrectionStats pureCanonicalCorrectionStats;

using PureCanonicalExtendedComplex=std::complex<long double>;
using PureCanonicalExtendedMatrix=Eigen::Matrix<PureCanonicalExtendedComplex,Eigen::Dynamic,Eigen::Dynamic>;

// This recovery starts from an independently certified double Gaussian frame.
// The step bound is still the original DOUBLE roundoff bound, never an
// overlap-conditioned tolerance. No Green matrix is projected or modified.
inline bool pureProjectorCanonicalMicrocorrectExtended(
    PureCanonicalExtendedMatrix &q,double &relative,double &isotropy,double orthTolerance) {
    auto &stats=pureCanonicalCorrectionStats;++stats.extended_calls;
    const auto fail=[&](){++stats.failures;return false;};
    const int d=q.rows(),n=q.cols();
    if(d<=0||d!=2*n||!q.allFinite())return fail();
    const long double normalization=std::sqrt(static_cast<long double>(n));
    const PureCanonicalExtendedMatrix delta=(-.5L*q.conjugate()*(q.transpose()*q)).eval();
    relative=double(delta.norm()/q.norm());
    if(!delta.allFinite()||!std::isfinite(relative)||relative>pureProjectorCanonicalCorrectionBound(d))return fail();
    const PureCanonicalExtendedMatrix corrected=q+delta;
    Eigen::HouseholderQR<PureCanonicalExtendedMatrix> qr(corrected);
    PureCanonicalExtendedMatrix next=qr.householderQ()*PureCanonicalExtendedMatrix::Identity(d,n);
    isotropy=double((next.transpose()*next).norm()/normalization);
    const double orth=double((next.adjoint()*next-PureCanonicalExtendedMatrix::Identity(n,n)).norm()/normalization);
    if(!next.allFinite()||!std::isfinite(isotropy)||!std::isfinite(orth)||
       isotropy>128.*d*std::numeric_limits<long double>::epsilon()||orth>orthTolerance)return fail();
    q.swap(next);return true;
}

// Evaluate the residual of the exact SAME stored double matrix. Scalar
// long-double accumulation prevents a large oblique G from failing merely
// because its double GEMM has cancellation error. The normalization and
// caller threshold are applied by the original certification policy.
inline double pureProjectorStoredGreenResidualExtended(const MatType &green) {
    ++pureCanonicalCorrectionStats.extended_green_residual_checks;
    const int d=green.rows();long double squared=0.;
    for(int i=0;i<d;++i)for(int j=0;j<d;++j) {
        PureCanonicalExtendedComplex value=0.;
        for(int k=0;k<d;++k)
            value+=PureCanonicalExtendedComplex(green(i,k).real(),green(i,k).imag())*
                   PureCanonicalExtendedComplex(green(k,j).real(),green(k,j).imag());
        if(i==j)value-=1.L;
        squared+=std::norm(value);
    }
    return double(std::sqrt(squared)/std::max(1.L,std::sqrt(static_cast<long double>(d))));
}

inline bool pureProjectorCanonicalFrameTrusted(const MatType &q,double orthTolerance) {
    if(q.rows()<=0 || q.rows()!=2*q.cols() || !q.allFinite())return false;
    const double normalization=std::sqrt(double(q.cols()));
    const double orth=(q.adjoint()*q-MatType::Identity(q.cols(),q.cols())).norm()/normalization;
    const double isotropy=(q.transpose()*q).norm()/normalization;
    return std::isfinite(orth)&&std::isfinite(isotropy)&&
        orth<=orthTolerance&&isotropy<=pureProjectorCanonicalCorrectionBound(q.rows());
}

// If X=Q R and Q+delta=Q_new R2, the maintained propagated frame is
// X_new=(Q+delta)R=Q_new(R2 R). Return det(R2) so every weight consumer
// retains its magnitude AND phase. For a Green-only consumer no scale is
// needed. The span change is first-order in the floating-point defect.
struct PureCanonicalWorkspace {
    PureComplexQRWorkspace qr;
    MatType isotropy,delta,corrected,next,gram;
};

inline bool pureProjectorCanonicalMicrocorrect(
    MatType &q,DataType &correctionDet,double orthTolerance,bool accelerated=false,
    PureCanonicalWorkspace *workspace=nullptr) {
    auto &stats=pureCanonicalCorrectionStats;++stats.complex_calls;
    const auto fail=[&](){++stats.failures;return false;};
    correctionDet=DataType(1.,0.);
    const int m=q.rows(),n=q.cols();
    if(m<=0||m!=2*n||!q.allFinite())return fail();
    const double normalization=std::sqrt(double(n));
    if(workspace&&accelerated&&n>=32) {
        auto &work=*workspace;
        // All products write independent scratch; q is committed only after
        // the same correction, QR determinant and residual gates succeed.
        work.gram.noalias()=q.adjoint()*q;
        work.gram-=MatType::Identity(n,n);
        const double orth=work.gram.norm()/normalization;
        if(!std::isfinite(orth)||orth>orthTolerance)return fail();
        work.isotropy.noalias()=q.transpose()*q;
        work.delta.noalias()=-.5*q.conjugate()*work.isotropy;
        const double relative=work.delta.norm()/q.norm(),before=work.isotropy.norm()/normalization;
        if(std::isfinite(relative))stats.maximum_relative_correction=std::max(stats.maximum_relative_correction,relative);
        if(std::isfinite(before))stats.maximum_input_isotropy=std::max(stats.maximum_input_isotropy,before);
        if(!work.isotropy.allFinite()||!work.delta.allFinite()||!std::isfinite(relative)||
           relative>pureProjectorCanonicalCorrectionBound(m))return fail();
        work.corrected=q+work.delta;work.next=work.corrected;
        if(!work.qr.factorize(work.next))return fail();
        for(int i=0;i<n;++i)correctionDet*=work.next(i,i);
        if(!work.qr.formQ(work.next))return fail();
        work.gram.noalias()=work.next.transpose()*work.next;
        const double after=work.gram.norm()/normalization;
        work.gram.noalias()=work.next.adjoint()*work.next;
        work.gram-=MatType::Identity(n,n);
        const double afterOrth=work.gram.norm()/normalization;
        if(std::isfinite(after))stats.maximum_output_isotropy=std::max(stats.maximum_output_isotropy,after);
        if(std::isfinite(afterOrth))stats.maximum_output_orthonormality=std::max(stats.maximum_output_orthonormality,afterOrth);
        if(!work.next.allFinite()||!std::isfinite(after)||!std::isfinite(afterOrth)||
           !std::isfinite(correctionDet.real())||!std::isfinite(correctionDet.imag())||
           !(std::abs(correctionDet)>0.)||after>pureProjectorCanonicalCorrectionBound(m)||
           afterOrth>orthTolerance)return fail();
        q.swap(work.next);return true;
    }
    const double orth=(q.adjoint()*q-MatType::Identity(n,n)).norm()/normalization;
    if(!std::isfinite(orth)||orth>orthTolerance)return fail();
    const MatType isotropy=q.transpose()*q;
    const MatType delta=(-.5*q.conjugate()*isotropy).eval();
    const double relative=delta.norm()/q.norm(),before=isotropy.norm()/normalization;
    if(std::isfinite(relative))stats.maximum_relative_correction=std::max(stats.maximum_relative_correction,relative);
    if(std::isfinite(before))stats.maximum_input_isotropy=std::max(stats.maximum_input_isotropy,before);
    if(!isotropy.allFinite()||!delta.allFinite()||!std::isfinite(relative)||
       relative>pureProjectorCanonicalCorrectionBound(m))return fail();
    MatType corrected=q+delta,next;
    if(accelerated&&n>=32) {
        cVecType tau(n);next=corrected;
        if(LAPACKE_zgeqrf(LAPACK_COL_MAJOR,m,n,next.data(),m,tau.data())!=0)return fail();
        for(int i=0;i<n;++i)correctionDet*=next(i,i);
        if(LAPACKE_zungqr(LAPACK_COL_MAJOR,m,n,n,next.data(),m,tau.data())!=0)return fail();
    } else {
        Eigen::HouseholderQR<MatType> qr(corrected);
        next=qr.householderQ()*MatType::Identity(m,n);
        for(int i=0;i<n;++i)correctionDet*=qr.matrixQR()(i,i);
    }
    const double after=(next.transpose()*next).norm()/normalization;
    const double afterOrth=(next.adjoint()*next-MatType::Identity(n,n)).norm()/normalization;
    if(std::isfinite(after))stats.maximum_output_isotropy=std::max(stats.maximum_output_isotropy,after);
    if(std::isfinite(afterOrth))stats.maximum_output_orthonormality=std::max(stats.maximum_output_orthonormality,afterOrth);
    if(!next.allFinite()||!std::isfinite(after)||!std::isfinite(afterOrth)||
       !std::isfinite(correctionDet.real())||!std::isfinite(correctionDet.imag())||
       !(std::abs(correctionDet)>0.)||after>pureProjectorCanonicalCorrectionBound(m)||
       afterOrth>orthTolerance)return fail();
    q.swap(next);return true;
}

#endif
