#ifndef PURE_PROJECTOR_ENDPOINT_H
#define PURE_PROJECTOR_ENDPOINT_H

#include "gaussian_trial_state.h"
#include "pure_projector_weight.h"
#include "qmc_core_profile.h"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

enum class PureEndpointRebuildStatus {
    success,
    invalid_configuration,
    nonfinite_scalar_prefactor,
    zero_scalar_prefactor,
    subspace_failure,
    endpoint_overlap_untrusted
};

inline const char *pureEndpointRebuildStatusName(PureEndpointRebuildStatus status) {
    switch (status) {
    case PureEndpointRebuildStatus::success: return "success";
    case PureEndpointRebuildStatus::invalid_configuration: return "invalid_configuration";
    case PureEndpointRebuildStatus::nonfinite_scalar_prefactor:
        return "nonfinite_scalar_prefactor";
    case PureEndpointRebuildStatus::zero_scalar_prefactor: return "zero_scalar_prefactor";
    case PureEndpointRebuildStatus::subspace_failure: return "subspace_failure";
    case PureEndpointRebuildStatus::endpoint_overlap_untrusted:
        return "endpoint_overlap_untrusted";
    }
    return "unknown";
}

struct PureEndpointRebuildResult {
    PureEndpointRebuildStatus status = PureEndpointRebuildStatus::invalid_configuration;
    int cut = -1;
    MatType q_right;
    MatType q_left;
    MatType green;
    int overlap_rank = 0;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    double green_residual = std::numeric_limits<double>::infinity();
    double green_skew_residual = std::numeric_limits<double>::infinity();
    double green_diagonal_residual = std::numeric_limits<double>::infinity();
    double green_residual_double = std::numeric_limits<double>::infinity();
    bool long_double_green_residual = false;
    bool canonical_frame_maintenance = false;
    double canonical_correction_right = 0.0,canonical_correction_left = 0.0;
    double log_abs_weight = -std::numeric_limits<double>::infinity();
    double log_abs_scalar_prefactor = -std::numeric_limits<double>::infinity();
    DataType scalar_prefactor_phase = DataType(0.0,0.0);
    std::string message;
    bool ok() const { return status == PureEndpointRebuildStatus::success; }
};

// Private to the caller/chain. Endpoint results retain independent matrices;
// none of these reusable buffers is exposed through a returned reference.
struct PureEndpointRebuildWorkspace {
    PureComplexQRWorkspace qr;
    PureCanonicalWorkspace canonical;
    PureProjectorGreenWorkspace green;
    MatType propagation;
};

namespace pure_projector_endpoint_detail {

struct ScaledThinQr {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType q;
    double log_abs_determinant = 0.0;
    bool ok() const { return status == PureProjectorStatus::success; }
};

// Diagnostics only; no counter participates in a numerical decision.
inline std::uint64_t endpointAcceleratedQrCalls=0;
inline std::uint64_t endpointReferenceQrCalls=0;
inline thread_local std::uint64_t endpointDirectDenseRightCalls=0;
inline thread_local std::uint64_t endpointDirectDenseLeftCalls=0;

// The input frame and immutable factor never alias this chain-private scratch.
// Preserve the same dense product (including complex adjoint) without expanding
// the factor or allocating an alias-protection temporary for every half-step.
// A null workspace and factors without dense storage keep the reference path.
inline void propagateFactor(const PureFactorMatrix &factor,MatType &frame,bool adjoint,
                            PureEndpointRebuildWorkspace *workspace) {
    if(workspace) {
        if(const MatType *dense=factor.denseStorage()) {
            if(adjoint) {
                workspace->propagation.noalias()=dense->adjoint()*frame;
                ++endpointDirectDenseLeftCalls;
            } else {
                workspace->propagation.noalias()=(*dense)*frame;
                ++endpointDirectDenseRightCalls;
            }
            frame.swap(workspace->propagation);
            return;
        }
    }
    if(adjoint)frame=factor.adjoint()*frame;
    else frame=factor*frame;
}

inline ScaledThinQr scaledThinQr(const MatType &subspace,
                                 const PureProjectorOptions &options, bool canonical = false,
                                 PureEndpointRebuildWorkspace *workspace=nullptr) {
    ScaledThinQr result;
    if (subspace.rows() <= 0 || subspace.cols() <= 0 ||
        subspace.rows() < subspace.cols() || !pureProjectorMatrixFinite(subspace)) {
        result.status = PureProjectorStatus::nonfinite_input;
        return result;
    }
    const auto recordScale=[&](const MatType &packed) {
        for (int i=0;i<subspace.cols();++i) {
            const double magnitude=std::abs(packed(i,i));
            if (!std::isfinite(magnitude) || magnitude <= options.rank_tolerance) {
                result.status = PureProjectorStatus::subspace_rank_deficient;
                return false;
            }
            result.log_abs_determinant += std::log(magnitude);
        }
        return true;
    };
    if(options.accelerated_linear_algebra && subspace.cols()>=32) {
        // Same unpivoted thin QR and diagonal-magnitude contract as reference.
        // Save the scales before zungqr overwrites R. No singular mode is removed.
        ++endpointAcceleratedQrCalls;
        result.q=subspace;
        const lapack_int m=subspace.rows(),n=subspace.cols();
        cVecType tau;
        bool qrOK;
        if(workspace)qrOK=workspace->qr.factorize(result.q);
        else {tau.resize(n);qrOK=LAPACKE_zgeqrf(LAPACK_COL_MAJOR,m,n,result.q.data(),m,tau.data())==0;}
        if(!qrOK || !pureProjectorMatrixFinite(result.q) || !recordScale(result.q)) {
            if(result.status==PureProjectorStatus::invalid_dimension)
                result.status=PureProjectorStatus::solve_residual_exceeded;
            return result;
        }
        const bool qOK=workspace?workspace->qr.formQ(result.q):
            LAPACKE_zungqr(LAPACK_COL_MAJOR,m,n,n,result.q.data(),m,tau.data())==0;
        if(!qOK) {
            result.status=PureProjectorStatus::solve_residual_exceeded;
            return result;
        }
    } else {
        ++endpointReferenceQrCalls;
        Eigen::HouseholderQR<MatType> qr(subspace);
        result.q = qr.householderQ() * MatType::Identity(subspace.rows(),subspace.cols());
        if(!recordScale(qr.matrixQR()))return result;
    }
    if(canonical) {
        DataType correctionDet;
        if(!pureProjectorCanonicalMicrocorrect(result.q,correctionDet,
                options.solve_residual_tolerance,options.accelerated_linear_algebra,
                workspace?&workspace->canonical:nullptr)) {
            result.status=PureProjectorStatus::solve_residual_exceeded;
            return result;
        }
        result.log_abs_determinant+=std::log(std::abs(correctionDet));
    }
    const double orth=(result.q.adjoint()*result.q-
        MatType::Identity(subspace.cols(),subspace.cols())).norm()/
        std::max(1.0,std::sqrt(double(subspace.cols())));
    if (!pureProjectorMatrixFinite(result.q) ||
        !std::isfinite(result.log_abs_determinant) ||
        orth > options.solve_residual_tolerance) {
        result.status = PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    result.status = PureProjectorStatus::success;
    return result;
}

} // namespace pure_projector_endpoint_detail

// Rebuild only the requested cut.  Propagation performs thin QR but never
// constructs a prefix Green function, so an untrustworthy intermediate
// trial-overlap cannot veto a trustworthy physical endpoint.
inline PureEndpointRebuildResult pureProjectorEndpointRebuild(
    const GaussianTrialState &trial,
    const std::vector<PureProjectorSlice> &slices,
    int cut,
    int blockSize,
    const PureProjectorOptions &options = PureProjectorOptions(),
    PureEndpointRebuildWorkspace *workspace=nullptr) {
    PureEndpointRebuildResult result;
    result.cut=cut;
    if (cut<0 || cut>int(slices.size()) || blockSize<=0 ||
        trial.Phi.rows()<=0 || trial.Phi.cols()<=0) {
        result.message="invalid endpoint rebuild configuration";
        return result;
    }
    double logEta=0.0;
    DataType etaPhase(1.0,0.0);
    bool ep_invalid=false;
    profileCore("ep_validate",[&]{
      for (int i=0;i<int(slices.size());++i) {
        const PureProjectorSlice &slice=slices[i];
        if (slice.matrix.rows()!=trial.Phi.rows() ||
            slice.matrix.cols()!=trial.Phi.rows() ||
            !pureProjectorMatrixFinite(slice.matrix)) {
            result.message="invalid factor at canonical index "+std::to_string(i);
            ep_invalid=true;return;
        }
        if (!std::isfinite(slice.eta.real()) || !std::isfinite(slice.eta.imag())) {
            result.status=PureEndpointRebuildStatus::nonfinite_scalar_prefactor;
            result.message="nonfinite scalar prefactor at canonical index "+std::to_string(i);
            ep_invalid=true;return;
        }
        const double magnitude=std::abs(slice.eta);
        if (!(magnitude>0.0)) {
            result.status=PureEndpointRebuildStatus::zero_scalar_prefactor;
            result.message="zero scalar prefactor at canonical index "+std::to_string(i);
            ep_invalid=true;return;
        }
        logEta+=std::log(magnitude);
        etaPhase*=slice.eta/magnitude;
      }
    });
    if(ep_invalid)return result;

    MatType right=trial.Phi,left=trial.Phi;
    double rightScale=0.0,leftScale=0.0;
    bool ep_failed=false;
    profileCore("ep_right",[&]{
      for (int index=0;index<cut;++index) {
        profileCore("ep_gemm",[&]{
            pure_projector_endpoint_detail::propagateFactor(slices[index].matrix,right,false,
                options.accelerated_linear_algebra?workspace:nullptr);
        });
        if (!pureProjectorMatrixFinite(right)) {
            result.status=PureEndpointRebuildStatus::subspace_failure;
            result.message="nonfinite right subspace at canonical index "+std::to_string(index);
            ep_failed=true;return;
        }
        if ((index+1)%blockSize==0 || index+1==cut) {
            pure_projector_endpoint_detail::ScaledThinQr qr;
            profileCore("ep_thin_qr",[&]{
                qr=pure_projector_endpoint_detail::scaledThinQr(right,options,true,workspace);
            });
            if (!qr.ok()) {
                result.status=PureEndpointRebuildStatus::subspace_failure;
                result.message="right thin-QR failed at canonical index "+std::to_string(index);
                ep_failed=true;return;
            }
            right=qr.q;rightScale+=qr.log_abs_determinant;
        }
      }
    });
    if(ep_failed)return result;
    int propagated=0;
    profileCore("ep_left",[&]{
      for (int index=int(slices.size())-1;index>=cut;--index) {
        profileCore("ep_gemm",[&]{
            pure_projector_endpoint_detail::propagateFactor(slices[index].matrix,left,true,
                options.accelerated_linear_algebra?workspace:nullptr);
        });
        ++propagated;
        if (!pureProjectorMatrixFinite(left)) {
            result.status=PureEndpointRebuildStatus::subspace_failure;
            result.message="nonfinite left subspace at canonical index "+std::to_string(index);
            ep_failed=true;return;
        }
        if (propagated%blockSize==0 || index==cut) {
            pure_projector_endpoint_detail::ScaledThinQr qr;
            profileCore("ep_thin_qr",[&]{
                qr=pure_projector_endpoint_detail::scaledThinQr(left,options,true,workspace);
            });
            if (!qr.ok()) {
                result.status=PureEndpointRebuildStatus::subspace_failure;
                result.message="left thin-QR failed at canonical index "+std::to_string(index);
                ep_failed=true;return;
            }
            left=qr.q;leftScale+=qr.log_abs_determinant;
        }
      }
    });
    if(ep_failed)return result;

    PureProjectorGreenResult endpoint;
    PureProjectorOptions gaussianOptions=options;gaussianOptions.fail_on_green_structure=true;
    profileCore("ep_green",[&]{ endpoint=pureProjectorGreen(right,left,gaussianOptions,
        workspace?&workspace->green:nullptr); });    result.q_right=right;result.q_left=left;
    result.overlap_rank=endpoint.overlap_rank;
    result.overlap_rcond=endpoint.overlap_rcond;
    result.solve_residual=endpoint.solve_residual;
    result.green_residual=endpoint.green_residual;
    result.green_skew_residual=endpoint.green_skew_residual;
    result.green_diagonal_residual=endpoint.green_diagonal_residual;
    result.green_residual_double=endpoint.green_residual_double;
    result.long_double_green_residual=endpoint.long_double_green_residual;
    result.canonical_frame_maintenance=endpoint.canonical_frame_maintenance;
    result.canonical_correction_right=endpoint.canonical_correction_right;
    result.canonical_correction_left=endpoint.canonical_correction_left;
    if (!endpoint.ok()) {
        result.status=PureEndpointRebuildStatus::endpoint_overlap_untrusted;
        result.message=std::string("target-cut overlap failed: ")+
            pureProjectorStatusName(endpoint.status);
        return result;
    }
    const MatType overlap=left.adjoint()*right;
    const DataType determinant=overlap.determinant();
    if (!std::isfinite(determinant.real()) || !std::isfinite(determinant.imag()) ||
        !(std::abs(determinant)>0.0)) {
        result.status=PureEndpointRebuildStatus::endpoint_overlap_untrusted;
        result.message="target-cut overlap determinant is zero or nonfinite";
        return result;
    }
    result.green=endpoint.green;
    result.log_abs_scalar_prefactor=logEta;
    result.scalar_prefactor_phase=etaPhase;
    result.log_abs_weight=logEta+0.5*(rightScale+leftScale+
        std::log(std::abs(determinant)));
    result.status=PureEndpointRebuildStatus::success;
    result.message="endpoint-only thin-QR rebuild trusted";
    return result;
}

#endif
