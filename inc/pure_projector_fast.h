#ifndef PURE_PROJECTOR_FAST_H
#define PURE_PROJECTOR_FAST_H

#include "pure_projector_endpoint.h"
#include "pure_projector_real.h"
#include "qmc_core_profile.h"
#include "pure_projector_green.h"
#include "pure_projector_mp.h"
#include "pure_projector_stack.h"
#include "pure_projector_local_support.h"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

enum class PureBranch { Ket, Bra };
enum class PureFastRunMode { AuditLockstep, FastStrict };
enum class PureFastInitializationPolicy {
    SequentialAudit,
    MirroredTheoremZ2Plus,
    EndpointTransportedZ2
};
enum class PureFastRatioStatus {
    success,
    invalid_proposal,
    rank_failure,
    overlap_untrusted,
    pfaffian_untrusted,
    nonfinite_ratio,
    complex_ratio,
    zero_ratio,
    reference_failure,
    terminated,
    positive_ratio_bound_violation
};

inline const char *pureFastRatioStatusName(PureFastRatioStatus status) {
    switch (status) {
    case PureFastRatioStatus::success: return "success";
    case PureFastRatioStatus::invalid_proposal: return "invalid_proposal";
    case PureFastRatioStatus::rank_failure: return "rank_failure";
    case PureFastRatioStatus::overlap_untrusted: return "overlap_untrusted";
    case PureFastRatioStatus::pfaffian_untrusted: return "pfaffian_untrusted";
    case PureFastRatioStatus::nonfinite_ratio: return "nonfinite_ratio";
    case PureFastRatioStatus::complex_ratio: return "complex_ratio";
    case PureFastRatioStatus::zero_ratio: return "zero_ratio";
    case PureFastRatioStatus::positive_ratio_bound_violation: return "positive_ratio_bound_violation";
    case PureFastRatioStatus::reference_failure: return "reference_failure";
    case PureFastRatioStatus::terminated: return "terminated";
    }
    return "unknown";
}

struct PureSliceLocation {
    PureBranch branch = PureBranch::Ket;
    int slice = -1;
    int factor = -1;
    int bond = -1;
    int aux = -1;
};

struct PureFastConfiguration {
    std::vector<PureProjectorSlice> slices;
    std::vector<int> hs_fields;
    std::vector<PureSliceLocation> locations;
};

enum class PureMirroredInitializationStatus {
    success,
    invalid_ket_protocol,
    nonfinite_scalar_prefactor,
    nonpositive_scalar_prefactor
};

struct PureMirroredInitializationResult {
    PureMirroredInitializationStatus status =
        PureMirroredInitializationStatus::invalid_ket_protocol;
    PureFastConfiguration configuration;
    DataType total_scalar_prefactor = DataType(0.0,0.0);
    int initial_z2_sign = 0;
    std::string initialization_policy;
    std::string message;
    bool ok() const { return status == PureMirroredInitializationStatus::success; }
};

// The input is the canonical ket protocol/action order.  Appending the
// adjoints in strict reverse order reuses the same noncommuting factor-order
// convention as PureStaticProjectorContour and gives U_bra=U_ket^dagger.
inline PureMirroredInitializationResult pureProjectorMirroredConfiguration(
    const PureFastConfiguration &ketProtocol,double realityTolerance) {
    PureMirroredInitializationResult result;
    if (ketProtocol.slices.empty() ||
        ketProtocol.slices.size()!=ketProtocol.hs_fields.size() ||
        ketProtocol.slices.size()!=ketProtocol.locations.size()) {
        result.message="invalid ket protocol";
        return result;
    }
    result.configuration=ketProtocol;
    for (int index=int(ketProtocol.slices.size())-1;index>=0;--index) {
        const PureProjectorSlice &source=ketProtocol.slices[index];
        if (!pureProjectorMatrixFinite(source.matrix) ||
            !std::isfinite(source.eta.real()) || !std::isfinite(source.eta.imag())) {
            result.status=PureMirroredInitializationStatus::nonfinite_scalar_prefactor;
            result.message="nonfinite mirrored ket factor or scalar";
            return result;
        }
        PureProjectorSlice mirrored(source.matrix.adjoint().eval(),
            std::conj(source.eta),"mirrored_adjoint:"+source.label);
        PureSliceLocation location=ketProtocol.locations[index];
        location.branch=PureBranch::Bra;
        result.configuration.slices.push_back(std::move(mirrored));
        result.configuration.hs_fields.push_back(ketProtocol.hs_fields[index]);
        result.configuration.locations.push_back(location);
    }
    DataType scalar(1.0,0.0);
    for (const auto &slice:result.configuration.slices) scalar*=slice.eta;
    result.total_scalar_prefactor=scalar;
    if (!std::isfinite(scalar.real()) || !std::isfinite(scalar.imag())) {
        result.status=PureMirroredInitializationStatus::nonfinite_scalar_prefactor;
        result.message="mirrored total scalar prefactor is nonfinite";
        return result;
    }
    if (!(scalar.real()>0.0) ||
        std::abs(scalar.imag())>realityTolerance*std::max(1.0,std::abs(scalar.real()))) {
        result.status=PureMirroredInitializationStatus::nonpositive_scalar_prefactor;
        result.message="mirrored total scalar prefactor is not positive real";
        return result;
    }
    result.status=PureMirroredInitializationStatus::success;
    result.initial_z2_sign=1;
    result.initialization_policy="mirrored_theorem_z2_plus";
    result.message="U_bra=U_ket^dagger and positive scalar prefactor imply W=||U_ket Psi_T||^2>0";
    return result;
}

struct PureFastProposal {
    int index = -1;
    int new_hs = 0;
    MatType new_factor;
    std::shared_ptr<const PureFactorMatrix> cached_factor;
    // An explicit immutable proposal may omit the dense value. A nonempty
    // public dense value always takes precedence, including after mutation.
    bool immutable_factor_only = false;
    PureFastProposal materialized()const{
        auto result=*this;
        if(immutable_factor_only&&new_factor.size()==0&&cached_factor)
            result.new_factor=*cached_factor;
        result.immutable_factor_only=false;
        return result;
    }
    DataType new_eta = DataType(1.0,0.0);
    double uniform = 0.0;
    // Opt-in physical contract, supplied by the model for an HS1 field flip.
    double positive_ratio_log_bound = -1.0;
};

struct PureProposalSnapshot {
    PureBranch branch = PureBranch::Ket;
    int slice = -1;
    int factor = -1;
    int bond = -1;
    int aux = -1;
    int old_hs = 0;
    int new_hs = 0;
    std::uint64_t configuration_hash = 0;
    DataType fast_ratio = DataType(0.0,0.0);
    double uniform = 0.0;
};

struct PureFastRatioResult {
    PureFastRatioStatus status = PureFastRatioStatus::invalid_proposal;
    DataType ratio = DataType(0.0,0.0);
    DataType slow_ratio = DataType(0.0,0.0);
    int low_rank = 0;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    PfaffianStatus pfaffian_status = PfaffianStatus::success;
    double relative_reference_error = std::numeric_limits<double>::infinity();
    double green_reference_error = std::numeric_limits<double>::infinity();
    MatType fast_updated_green;
    PureMpProposalResult mp_reference;
    bool trust_alarm = false;
    bool used_reference = false;
    bool used_trusted_green = false;
    bool used_local_support = false;
    bool used_green_update = false;
    bool ok() const { return status == PureFastRatioStatus::success; }
};

struct PureFastProposalResult {
    PureProposalSnapshot snapshot;
    PureFastRatioResult ratio;
    bool accepted = false;
    bool terminated = false;
    bool reference_accepted = false;
    bool z2_reference_mismatch = false;
    // Stable diagnostic tag for fail-closed exits.  It does not participate
    // in the Markov decision and therefore cannot change the trajectory.
    std::string failure_stage;
    double fast_seconds = 0.0;
    double slow_seconds = 0.0;
};

struct PureFastOptions {
    // Opt-in for callers needing only overlap rank and singular values here.
    // Default preserves legacy-model arithmetic and diagnostics.
    bool reference_singular_values_only = false;
    // Explicit SSH-only magnitude opt-in; signed and generic callers retain
    // their original complex propagation.
    bool real_magnitude = false;
    PureProjectorWeightMode weight_mode = PureProjectorWeightMode::RealZ2;
    double rank_tolerance = 1e-12;
    double minimum_overlap_rcond = 1e-12;
    double residual_tolerance = 1e-10;
    double reality_tolerance = 1e-10;
    double zero_tolerance = 1e-14;
    double decision_margin_tolerance = 1e-12;
    double green_update_tolerance = 1e-8;
    // Local reconstruction policy; strict/reference keeps its existing limit.
    double local_green_rebuild_tolerance = 1e-7;
    std::uint64_t read_only_audit_interval = 0;
    bool mp_same_proposal_fallback = true;
    bool accelerated_linear_algebra = true;
    // Independent/reference callers do not supply a cached live-cut Green.
    bool reuse_trusted_green = true;
    bool compact_local_update = true;
    // Production policy: propagate/update G locally, retain periodic double
    // stabilization, and stop unresolved numerical errors without launching MP.
    // Strict/reference callers keep the historical behavior by default.
    bool stabilized_local_updates = false;
    bool retain_measurement_rebuild_cache = false;
    bool reuse_stabilized_subspaces = true;
    // Exact immutable factor support; independent/reference arithmetic stays dense.
    bool cached_factor_support = true;
    bool pair_block_factors = true;
    bool direct_green_columns = true;
    // Reuse local reconstruction frames and factor-application workspaces
    // across proposals. The exact QR checkpoints and trust gates are kept;
    // only temporary matrix storage is retained between calls.
    bool reuse_rebuild_workspace = true;
    // Local-only lifetime/allocation optimizations; public results still own
    // their matrices and currentWeight().green is synchronized immediately.
    bool borrow_green = true;
    // Buffer reuse reduced allocation volume further, but borrowing alone
    // was faster in both paired L98/L100 runs. Keep reuse for ablation.
    bool reuse_green_buffers = false;
    // Reuse an already certified solve at the same exact QR checkpoint pair.
    bool reuse_rebuild_result = true;
    bool reuse_measurement_green = true;
    bool reuse_long_jump_checkpoints = true;
    // Local updates arriving from the right use the existing index+1 cut;
    // retain the round trip as an ablation/reference option.
    bool avoid_reverse_round_trip = true;
    bool compact_proposal_algebra = true;
    // Four-way L98/L100 timings did not show an additional sweep-time win
    // from deferral. Keep it available for ablation, off in production.
    bool defer_rejected_green = false;
    double local_reality_tolerance = 1e-7;
    // Enable only with a proved positive contour (not for arbitrary HS1 models).
    bool positive_hs1_ratio_contract = false;
};

// Complete state at a signed-overlap stabilization boundary. A propagated Q
// alone is insufficient: continuing the recurrence also needs the accumulated
// determinant and overlap normalization/phase state.
struct PureProjectorSignedPrefixCheckpoint {
    MatType q;
    double log_det_r=0., previous_log_overlap=0.;
    DataType det_r_phase=DataType(1,0), previous_overlap_phase=DataType(1,0);
    double log_abs_weight=0.;
    DataType complex_phase=DataType(1,0);
    int z2_sign=1, overlap_rank=0;
    double overlap_rcond=1., overlap_residual=0.;
    PfaffianStatus pfaffian_status=PfaffianStatus::success;
    std::uint64_t prefix_hash=1469598103934665603ULL;
};

struct PureProjectorSignedPrefixCache {
    std::vector<std::shared_ptr<const PureProjectorSignedPrefixCheckpoint>> checkpoints;
};

inline std::uint64_t pureProjectorSlicePrefixHash(
    std::uint64_t hash,const PureProjectorSlice &slice) {
    slice.matrix.mixHash(hash);
    std::uint64_t bits=0;
    static_assert(sizeof(bits)==sizeof(double),"unexpected double size");
    double real=slice.eta.real(),imag=slice.eta.imag();
    std::memcpy(&bits,&real,sizeof(bits)); hash^=bits; hash*=1099511628211ULL;
    std::memcpy(&bits,&imag,sizeof(bits)); hash^=bits; hash*=1099511628211ULL;
    return hash;
}

struct PureFastDiagnostics {
    std::uint64_t local_double_retries = 0;
    std::uint64_t reverse_cut_shortcuts = 0;
    double maximum_local_ratio_imaginary = 0.0;
    std::uint64_t proposals = 0;
    std::uint64_t rebuild_count = 1;
    std::uint64_t ratio_slow_reference_count = 0;
    std::uint64_t trust_alarm_count = 0;
    std::uint64_t slow_reference_failure_count = 0;
    std::uint64_t mp_fallback_count = 0;
    std::uint64_t mp_fallback_failure_count = 0;
    std::uint64_t mp_endpoint_fallback_count = 0;
    std::uint64_t mp_endpoint_fallback_failure_count = 0;
    std::uint64_t mp_endpoint_rejection_recovery_count = 0;
    std::uint64_t mp_green_double_recovery_count = 0;
    std::uint64_t accepted_factor_center_recovery_count = 0;
    std::uint64_t accepted_factor_center_recovery_failure_count = 0;
    std::uint64_t cut_rebuild_recovery_count = 0;
    std::uint64_t cut_rebuild_reference_count = 0;
    std::uint64_t post_decision_center_recovery_count = 0;
    std::uint64_t post_decision_center_recovery_failure_count = 0;
    PureStackStatus last_accepted_factor_failure_status = PureStackStatus::success;
    std::uint64_t read_only_endpoint_audit_failure_count = 0;
    double minimum_overlap_rcond = 1.0;
    double maximum_green_rebuild_error = 0.0;
    double maximum_ratio_reference_error = 0.0;
    double maximum_endpoint_rebuild_green_residual = 0.0;
    double total_fast_seconds = 0.0;
    double total_reference_seconds = 0.0;
    double total_mp_fallback_seconds = 0.0;
    PureMpPerformanceProfile mp_profile;
    int first_failure_proposal = -1;
};

struct PureEndpointTransportInitializationResult {
    PureProjectorWeightStatus status = PureProjectorWeightStatus::invalid_slice;
    PureProjectorWeightResult weight;
    std::uint64_t transported_factors = 0;
    std::uint64_t mp_fallback_count = 0;
    std::uint64_t mp_failure_count = 0;
    double mirrored_log_abs_weight = -std::numeric_limits<double>::infinity();
    double transported_log_abs_weight = -std::numeric_limits<double>::infinity();
    double endpoint_log_abs_weight = -std::numeric_limits<double>::infinity();
    double log_abs_consistency_error = std::numeric_limits<double>::infinity();
    std::string message;
    bool ok() const { return status == PureProjectorWeightStatus::success; }
};

// Single source of truth for the Green/QR tolerances used by the stable
// signed-overlap reference weight, so a prepared trial QR is built with the
// exact gates of the call that consumes it.
inline PureProjectorOptions pureFastGreenOptions(const PureFastOptions &options) {
    PureProjectorOptions greenOptions;
    greenOptions.rank_tolerance=options.rank_tolerance;
    greenOptions.minimum_overlap_rcond=options.minimum_overlap_rcond;
    greenOptions.solve_residual_tolerance=std::max(options.residual_tolerance,1e-9);
    greenOptions.green_residual_tolerance=options.green_update_tolerance;
    greenOptions.allow_oblique_green=true;
    greenOptions.prefer_lu_green_solve=true;
    greenOptions.accelerated_linear_algebra=options.accelerated_linear_algebra;
    // Scale the relaxed production solve gate by the actual measured rcond,
    // rather than by its configured lower bound. Production uses 1e-9 for
    // rcond <= 0.1; well-conditioned overlaps can use up to 1e-8. The Green
    // residual gate still checks the resulting matrix before it is used.
    greenOptions.well_conditioned_solve_tolerance=1e-8;
    return greenOptions;
}

// One workspace belongs to one projector chain. Eigen's compute() methods
// retain their internal allocations when dimensions and computation options
// stay unchanged; callers still get a fresh local workspace when they do not
// provide one. No numerical thresholds or decomposition options are changed.
struct PureStableWeightWorkspace {
    MatType q,x,next_q,overlap,identity;
    Eigen::HouseholderQR<MatType> qr;
    Eigen::JacobiSVD<MatType> svd;
    PureOverlapSingularValues magnitude_svd;
    PureRealMagnitudeWorkspace real_magnitude;
    std::vector<DataType> qr_tau;
};

// The validated SSH factors retain an immutable dense matrix.  Keep this
// counter separate from the general core profile: it records which dispatch
// path was selected, not an allocator-level claim about Eigen internals.
struct PureProjectorDensePropagationStats {
    unsigned long long stable_direct_dense=0;
    unsigned long long stable_expanded_factor=0;
    unsigned long long magnitude_direct_dense=0;
    unsigned long long magnitude_expanded_factor=0;
};
inline thread_local PureProjectorDensePropagationStats pureProjectorDensePropagationStats;

inline void pureProjectorPropagate(const PureFactorMatrix &factor,const MatType &q,
                                   MatType &x,bool magnitudePath) {
    if(const MatType *dense=factor.denseStorage()) {
        // The pointer is borrowed only for this expression; factor owns the
        // immutable storage and outlives the call.
        x=(*dense)*q;
        if(magnitudePath)++pureProjectorDensePropagationStats.magnitude_direct_dense;
        else ++pureProjectorDensePropagationStats.stable_direct_dense;
    } else {
        // Preserve the generic ReturnByValue path for sparse/local factors.
        x=factor*q;
        if(magnitudePath)++pureProjectorDensePropagationStats.magnitude_expanded_factor;
        else ++pureProjectorDensePropagationStats.stable_expanded_factor;
    }
}

// Diagnostic: how often a supplied spectral Cayley failed its gates and the
// FullPivLU reference solve was used instead (should be ~0 in the tested regime).
inline thread_local unsigned long long spectralCayleyFallbackCount=0;
inline thread_local unsigned long long signedPrefixResumeCount=0;
inline thread_local unsigned long long signedPrefixSkippedSliceCount=0;

// Only this factory can create a record: its provenance is a successful
// reference FullPivLU solve, never a spectral construction. Keeping the factor
// handle alive makes pointer identity safe against storage-address reuse.
class PureProjectorLUCayley {
    PureFactorMatrix factor_;
    MatType value_;
    double zero_,rank_,rcond_,residual_,reality_,green_;
    PureProjectorLUCayley(const PureFactorMatrix &factor,const MatType &value,
                         const PureFastOptions &o):factor_(factor),value_(value),
        zero_(o.zero_tolerance),rank_(o.rank_tolerance),rcond_(o.minimum_overlap_rcond),
        residual_(o.residual_tolerance),reality_(o.reality_tolerance),green_(o.green_update_tolerance) {}
public:
    bool matches(const PureFactorMatrix &factor,const PureFastOptions &o)const {
        return factor_.sharesStorageWith(factor)&&zero_==o.zero_tolerance&&
            rank_==o.rank_tolerance&&rcond_==o.minimum_overlap_rcond&&
            residual_==o.residual_tolerance&&reality_==o.reality_tolerance&&green_==o.green_update_tolerance;
    }
    const MatType &value()const {return value_;}
    static bool solve(const PureFactorMatrix &factor,const PureFastOptions &o,
                      MatType &value,std::shared_ptr<const PureProjectorLUCayley> *record) {
        MatType onePlus=MatType::Identity(factor.rows(),factor.cols())+factor;
        Eigen::FullPivLU<MatType> lu(onePlus);
        if(!lu.isInvertible())return false;
        value=lu.solve(2.0*MatType::Identity(factor.rows(),factor.cols()))-
            MatType::Identity(factor.rows(),factor.cols());
        // Nonfinite results still follow the original downstream failure path;
        // they are never installed as reusable records.
        if(record) {
            record->reset();
            if(pureProjectorMatrixFinite(value))
                *record=std::shared_ptr<const PureProjectorLUCayley>(new PureProjectorLUCayley(factor,value,o));
        }
        return true;
    }
};
using PureProjectorLUCayleyCache=std::vector<std::shared_ptr<const PureProjectorLUCayley>>;
inline std::size_t pureProjectorCayleyBytes(const PureProjectorLUCayleyCache &cache) {
    std::size_t bytes=sizeof(cache)+cache.capacity()*sizeof(cache[0]);
    for(const auto &record:cache)if(record)
        bytes+=sizeof(*record)+std::size_t(record->value().size())*sizeof(DataType);
    return bytes; // excludes shared factor storage and allocator/control blocks
}
inline thread_local unsigned long long cayleyLUSolveCount=0,cayleyCacheHitCount=0,cayleyCacheFallbackCount=0;
inline thread_local unsigned long long signedPropagationSliceCount=0;

inline PureProjectorWeightResult pureProjectorStableReferenceWeight(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice>&slices,
    const PureFastOptions &options,
    const PreparedGreenTrial *preparedTrial=nullptr,
    const std::vector<MatType> *preparedCayley=nullptr,
    PureStableWeightWorkspace *workspace=nullptr,
    const PureProjectorSignedPrefixCache *prefixCache=nullptr,
    PureProjectorSignedPrefixCache *updatedPrefixCache=nullptr,
    int requestedResumeIndex=-1,
    PureProjectorLUCayleyCache *cayleyCache=nullptr) {    const PureProjectorOptions greenOptions=pureFastGreenOptions(options);
    PureProjectorWeightResult result;result.status=PureProjectorWeightStatus::success;
    result.log_abs_weight=0;result.complex_phase=DataType(1,0);result.z2_sign=1;
    result.overlap_rank=trial.Phi.cols();result.overlap_rcond=1;result.overlap_residual=0;
    PureStableWeightWorkspace localWorkspace;
    PureStableWeightWorkspace &work=workspace?*workspace:localWorkspace;
    work.q=trial.Phi;double logDetR=0,previousLogOverlap=0;DataType detRPhase(1,0),previousOverlapPhase(1,0);
    // A1: reuse the cached thin QR of the fixed trial state and the already
    // orthonormal propagated basis; fall back to the full two-sided thin QR
    // whenever the prepared basis is unavailable or fails its residual gate.
    auto greenFor=[&](const MatType &qr)->PureProjectorGreenResult{
        const PureProjectorGreenResult raw=preparedTrial
            ? pureProjectorGreenPreparedRight(qr,trial,*preparedTrial,greenOptions)
            : pureProjectorGreenThinQr(qr,trial.Phi,greenOptions);
        // A large Green norm is expected for an oblique projector.  The
        // dimension-normalized involution residual can therefore reject a
        // numerically valid solve even when the relative residual is tiny.
        // Recheck only this certified solve failure; no matrix or update is
        // changed, and all overlap/solve/structure gates remain intact.
        return raw.ok()?raw:pureProjectorCheckObliqueGreen(raw,greenOptions);
    };
    const int sliceCount=int(slices.size());
    if(cayleyCache)cayleyCache->resize(sliceCount);
    int startIndex=0;
    constexpr std::uint64_t initialPrefixHash=1469598103934665603ULL;
    std::uint64_t prefixHash=initialPrefixHash;
    if(updatedPrefixCache) {
        updatedPrefixCache->checkpoints.clear();
        updatedPrefixCache->checkpoints.reserve(std::size_t(sliceCount)+1);
    }
    auto saveCheckpoint=[&](std::uint64_t hash) {
        if(!updatedPrefixCache)return;
        profileCore("sw_prefix_checkpoint",[&] {
            auto checkpoint=std::make_shared<PureProjectorSignedPrefixCheckpoint>();
            checkpoint->q=work.q;checkpoint->log_det_r=logDetR;
            checkpoint->previous_log_overlap=previousLogOverlap;
            checkpoint->det_r_phase=detRPhase;
            checkpoint->previous_overlap_phase=previousOverlapPhase;
            checkpoint->log_abs_weight=result.log_abs_weight;
            checkpoint->complex_phase=result.complex_phase;
            checkpoint->z2_sign=result.z2_sign;
            checkpoint->overlap_rank=result.overlap_rank;
            checkpoint->overlap_rcond=result.overlap_rcond;
            checkpoint->overlap_residual=result.overlap_residual;
            checkpoint->pfaffian_status=result.pfaffian_status;
            checkpoint->prefix_hash=hash;
            updatedPrefixCache->checkpoints.push_back(std::move(checkpoint));
        });
    };
    bool resumed=false;
    profileCore("sw_prefix_resume",[&] {
        if(prefixCache&&requestedResumeIndex>0&&requestedResumeIndex<=sliceCount&&
           prefixCache->checkpoints.size()>std::size_t(requestedResumeIndex)&&
           prefixCache->checkpoints[requestedResumeIndex]) {
            for(int i=0;i<requestedResumeIndex;++i)
                prefixHash=pureProjectorSlicePrefixHash(prefixHash,slices[i]);
            const auto &checkpoint=*prefixCache->checkpoints[requestedResumeIndex];
            const bool finiteCheckpoint=std::isfinite(checkpoint.log_det_r)&&
                std::isfinite(checkpoint.previous_log_overlap)&&
                std::isfinite(checkpoint.log_abs_weight)&&
                std::isfinite(checkpoint.overlap_rcond)&&
                std::isfinite(checkpoint.overlap_residual)&&
                checkpoint.overlap_rank>=0&&
                checkpoint.pfaffian_status==PfaffianStatus::success&&
                pureProjectorMatrixFinite(checkpoint.q)&&
                std::isfinite(checkpoint.det_r_phase.real())&&
                std::isfinite(checkpoint.det_r_phase.imag())&&
                std::isfinite(checkpoint.previous_overlap_phase.real())&&
                std::isfinite(checkpoint.previous_overlap_phase.imag())&&
                std::isfinite(checkpoint.complex_phase.real())&&
                std::isfinite(checkpoint.complex_phase.imag());
            if(finiteCheckpoint&&checkpoint.prefix_hash==prefixHash&&
               checkpoint.q.rows()==trial.Phi.rows()&&checkpoint.q.cols()==trial.Phi.cols()) {
                work.q=checkpoint.q;logDetR=checkpoint.log_det_r;
                previousLogOverlap=checkpoint.previous_log_overlap;
                detRPhase=checkpoint.det_r_phase;
                previousOverlapPhase=checkpoint.previous_overlap_phase;
                result.log_abs_weight=checkpoint.log_abs_weight;
                result.complex_phase=checkpoint.complex_phase;
                result.z2_sign=checkpoint.z2_sign;
                result.overlap_rank=checkpoint.overlap_rank;
                result.overlap_rcond=checkpoint.overlap_rcond;
                result.overlap_residual=checkpoint.overlap_residual;
                result.pfaffian_status=checkpoint.pfaffian_status;
                startIndex=requestedResumeIndex;resumed=true;
                ++signedPrefixResumeCount;
                signedPrefixSkippedSliceCount+=std::uint64_t(startIndex);
            }
        }
    });
    // A rejected checkpoint must leave both the propagation state and the
    // checkpoint-only hash at the full-reference initial state. Validation is
    // deliberately completed before assigning any cached propagation fields.
    if(!resumed)prefixHash=initialPrefixHash;
    if(updatedPrefixCache) {
        if(resumed) {
            for(int i=0;i<=startIndex;++i)
                updatedPrefixCache->checkpoints.push_back(prefixCache->checkpoints[i]);
        } else {
            auto initial=std::make_shared<PureProjectorSignedPrefixCheckpoint>();
            initial->q=work.q;initial->overlap_rank=result.overlap_rank;
            initial->overlap_rcond=result.overlap_rcond;initial->prefix_hash=prefixHash;
            updatedPrefixCache->checkpoints.push_back(std::move(initial));
        }
    }
    for(int index=startIndex;index<sliceCount;++index){const auto&s=slices[index];
        ++signedPropagationSliceCount;
        if(s.matrix.rows()!=trial.Phi.rows()||s.matrix.cols()!=trial.Phi.rows()||
           (!s.validated_dense&&!pureProjectorMatrixFinite(s.matrix))||
           !std::isfinite(s.eta.real())||!std::isfinite(s.eta.imag())){
            result.status=PureProjectorWeightStatus::invalid_slice;result.first_failing_slice=index;return result;}
        DataType detR(1,0);bool propagation_failure=false;
        profileCore("sw_propagate_qr",[&]{
            pureProjectorPropagate(s.matrix,work.q,work.x,false);
            if(!pureProjectorMatrixFinite(work.x)){propagation_failure=true;return;}
            const int m=work.x.rows(),n=work.x.cols();
            if(options.reference_singular_values_only&&options.accelerated_linear_algebra&&n>=32) {
                // Preserve the complex determinant phase before overwriting
                // R with Q. Basis rephasings cancel against the overlap phase.
                work.next_q=work.x;work.qr_tau.resize(n);
                if(LAPACKE_zgeqrf(LAPACK_COL_MAJOR,m,n,work.next_q.data(),m,
                                  work.qr_tau.data())!=0){propagation_failure=true;return;}
                for(int i=0;i<n;++i)detR*=work.next_q(i,i);
                if(LAPACKE_zungqr(LAPACK_COL_MAJOR,m,n,n,work.next_q.data(),m,
                                  work.qr_tau.data())!=0){propagation_failure=true;return;}
            } else {
                work.qr.compute(work.x);
                work.identity.resize(work.x.rows(),work.x.cols());work.identity.setIdentity();
                work.next_q=work.qr.householderQ()*work.identity;
                const auto &matrixQR=work.qr.matrixQR();
                for(int i=0;i<work.x.cols();++i)detR*=matrixQR(i,i);
            }
        });
        if(propagation_failure){
            result.status=PureProjectorWeightStatus::propagation_failure;result.first_failing_slice=index;return result;}
        if(!std::isfinite(detR.real())||!std::isfinite(detR.imag())||std::abs(detR)<=options.zero_tolerance){
            result.status=PureProjectorWeightStatus::zero_weight;result.first_failing_slice=index;result.z2_sign=0;return result;}
        double nextLogOverlap=0;DataType nextOverlapPhase(1,0);
        profileCore("sw_overlap_svd",[&]{
            logDetR+=std::log(std::abs(detR));detRPhase*=pureProjectorUnitPhase(detR);
            work.overlap=trial.Phi.adjoint()*work.next_q;
            const dVecType *singular=nullptr;
            if(options.reference_singular_values_only) {
                // A small normalized determinant is not a physical zero:
                // signed weights restore the QR scale in logarithmic form.
                work.magnitude_svd.compute(work.overlap,options.accelerated_linear_algebra,
                    options.rank_tolerance,options.minimum_overlap_rcond,options.zero_tolerance,false);
                result.overlap_rank=work.magnitude_svd.rank(options.rank_tolerance);
                singular=&work.magnitude_svd.singularValues();
            } else {
                work.svd.compute(work.overlap,Eigen::ComputeThinU|Eigen::ComputeThinV);
                work.svd.setThreshold(options.rank_tolerance);result.overlap_rank=work.svd.rank();
                singular=&work.svd.singularValues();
            }
            if(result.overlap_rank!=work.overlap.rows()){result.status=PureProjectorWeightStatus::zero_weight;
                result.first_failing_slice=index;result.z2_sign=0;return;}
            const auto &sv=*singular;
            double largest=sv(0),smallest=sv(work.overlap.rows()-1);
            result.overlap_rcond=largest>0?smallest/largest:0;if(result.overlap_rcond<options.minimum_overlap_rcond){
                result.status=PureProjectorWeightStatus::overlap_untrusted;result.first_failing_slice=index;return;}
            // The determinant of Q_trial^dag Q is a product of principal
            // overlaps, not the physical weight: QR scales live in logDetR.
            // Full rank and rcond certify the solve even when that product is
            // tiny. Accumulate its magnitude in log space and its LU phase
            // one pivot at a time, without a dimension-dependent zero gate.
            double logAbsDet=0.;
            for(int i=0;i<sv.size();++i)logAbsDet+=std::log(sv(i));
            Eigen::PartialPivLU<MatType> lu(work.overlap);
            DataType detMPhase(double(lu.permutationP().determinant()),0.);
            for(int i=0;i<work.overlap.rows();++i) {
                const DataType pivot=lu.matrixLU()(i,i);
                if(!std::isfinite(pivot.real())||!std::isfinite(pivot.imag())||std::abs(pivot)==0.) {
                    result.status=PureProjectorWeightStatus::overlap_untrusted;
                    result.first_failing_slice=index;return;
                }
                detMPhase*=pureProjectorUnitPhase(pivot);
            }
            if(!std::isfinite(logAbsDet)) {
                result.status=PureProjectorWeightStatus::overlap_untrusted;
                result.first_failing_slice=index;return;
            }
            nextLogOverlap=logDetR+logAbsDet;nextOverlapPhase=
                pureProjectorUnitPhase(detRPhase*pureProjectorUnitPhase(detMPhase));
        });
        if(result.status!=PureProjectorWeightStatus::success)return result;
        double ratioLog=std::log(std::abs(s.eta))+.5*(nextLogOverlap-previousLogOverlap);
        // ratioLog already includes |eta|; only its unit phase enters the root.
        const DataType etaPhase=pureProjectorUnitPhase(s.eta);
        DataType ratio=std::exp(ratioLog)*std::sqrt(etaPhase*etaPhase*nextOverlapPhase/previousOverlapPhase);
        const double identityResidual=s.identity_residual>=0.0?s.identity_residual:
            (s.matrix-MatType::Identity(s.matrix.rows(),s.matrix.cols())).norm();
        if(identityResidual<=10*options.zero_tolerance)ratio=s.eta;else{
            MatType sliceGreen;bool cayleyPrepared=false;
            // A2: prefer the node spectral Cayley C=U tanh(a eps/2) U^dag, but
            // only after an independent residual certificate against the actual
            // slice matrix; any discrepancy falls back to the FullPivLU solve.
            profileCore("sw_cayley_prepare",[&]{
                if(preparedCayley&&index<int(preparedCayley->size())){
                    const MatType &proposed=(*preparedCayley)[index];
                    const int d=s.matrix.rows();
                    if(proposed.rows()==d&&proposed.cols()==d&&pureProjectorMatrixFinite(proposed)){
                        const MatType eye=MatType::Identity(d,d);
                        const double herm=(proposed-proposed.adjoint()).norm()/
                            std::max(1.0,std::sqrt(double(d)));
                        const MatType residual=(eye+s.matrix)*proposed-(eye-s.matrix);
                        const double res=residual.norm()/std::max(1.0,2.0*std::sqrt(double(d)));
                        if(std::isfinite(herm)&&std::isfinite(res)&&
                           herm<=options.residual_tolerance&&res<=options.residual_tolerance){
                            sliceGreen=proposed;cayleyPrepared=true;
                        }
                    }
                }
            });
            // Spectral data and LU records are deliberately separate. This
            // cache is only consulted on a path which would otherwise solve.
            if(!cayleyPrepared&&cayleyCache) {
                const auto &record=(*cayleyCache)[index];
                if(record&&record->matches(s.matrix,options)) {
                    sliceGreen=record->value();cayleyPrepared=true;++cayleyCacheHitCount;
                } else ++cayleyCacheFallbackCount;
            }
            if(!cayleyPrepared){
                if(preparedCayley)++spectralCayleyFallbackCount;
                profileCore("sw_cayley_solve",[&]{
                    ++cayleyLUSolveCount;
                    if(!PureProjectorLUCayley::solve(s.matrix,options,sliceGreen,
                            cayleyCache?&(*cayleyCache)[index]:nullptr)) {
                        result.status=PureProjectorWeightStatus::pfaffian_untrusted;result.first_failing_slice=index;return;
                    }
                });
            }
            if(result.status!=PureProjectorWeightStatus::success)return result;            PureProjectorGreenResult partial;
            profileCore("sw_partial_green",[&]{
                partial=greenFor(work.q);
            });            if(!partial.ok()){
                result.status=PureProjectorWeightStatus::overlap_untrusted;
                result.first_failing_slice=index;
                result.green_status=partial.status;
                result.green_overlap_rcond=partial.overlap_rcond;
                result.green_solve_residual=partial.solve_residual;
                result.green_residual=partial.green_residual;
                result.green_skew_residual=partial.green_skew_residual;
                result.green_diagonal_residual=partial.green_diagonal_residual;
                return result;
            }            PfaffianResult pf;
            profileCore("sw_pfaffian",[&]{
                pf=pfaffianForSignOfProductWithStatus(sliceGreen,partial.green);result.pfaffian_status=pf.status;
                if(!pf.ok()){result.status=PureProjectorWeightStatus::pfaffian_untrusted;result.first_failing_slice=index;return;}
                // Normalize the Pfaffian ordering at the identity, as in local ratios.
                const double orderingSign=trial.Phi.cols()%2 ? -1.0 : 1.0;
                DataType target=pureProjectorUnitPhase(orderingSign*s.eta*pf.value);if((std::conj(target)*ratio).real()<0)ratio=-ratio;
            });
            if(result.status!=PureProjectorWeightStatus::success)return result;
        }
        if(!std::isfinite(ratio.real())||!std::isfinite(ratio.imag())||std::abs(ratio)<=options.zero_tolerance){
            result.status=PureProjectorWeightStatus::zero_weight;result.first_failing_slice=index;result.z2_sign=0;return result;}
        result.log_abs_weight+=std::log(std::abs(ratio));result.complex_phase=
            pureProjectorUnitPhase(result.complex_phase*pureProjectorUnitPhase(ratio));
        if(std::abs(result.complex_phase.imag())<=options.reality_tolerance)
            result.z2_sign=result.complex_phase.real()>=0?1:-1;else result.z2_sign=0;
        work.q.swap(work.next_q);previousLogOverlap=nextLogOverlap;previousOverlapPhase=nextOverlapPhase;
        // The rolling hash exists only to validate/build checkpoint caches.
        // Avoid scanning immutable matrices when no cache is being generated.
        if(updatedPrefixCache) {
            prefixHash=pureProjectorSlicePrefixHash(prefixHash,s);
            saveCheckpoint(prefixHash);
        }
    }
    PureProjectorGreenResult green;
    profileCore("sw_final_green",[&]{green=greenFor(work.q);});    if(!green.ok()){
        result.status=PureProjectorWeightStatus::overlap_untrusted;
        result.first_failing_slice=slices.empty()?-1:int(slices.size())-1;
        result.green_status=green.status;
        result.green_overlap_rcond=green.overlap_rcond;
        result.green_solve_residual=green.solve_residual;
        result.green_residual=green.green_residual;
        result.green_skew_residual=green.green_skew_residual;
        result.green_diagonal_residual=green.green_diagonal_residual;
        return result;
    }    result.green=green.green;result.overlap_rank=green.overlap_rank;result.overlap_rcond=green.overlap_rcond;
    result.overlap_residual=green.solve_residual;result.green_residual=green.green_residual;
    result.weight=result.log_abs_weight<700?std::exp(result.log_abs_weight)*result.complex_phase:result.complex_phase;
    result.determinant_identity_error=0;return result;
}

// D^dag B D and the certified trial span are real. QR and singular values
// therefore give the same |det overlap| under a change of occupied basis.
// A failure or proximity to ANY original trust boundary replays the entire
// same candidate using complex arithmetic below; no real-path status decides
// a rank/zero rejection. Signed phase and all Green routines remain complex.
inline bool pureProjectorRealMagnitudeWeight(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice> &slices,
    const PureFastOptions &options,PureRealMagnitudeWorkspace &work,
    PureProjectorWeightResult &result) {
    const int n=trial.Phi.cols();
    if(n<16||PureRealMagnitudeWorkspace::disabled()||!work.trialMatches(trial))return false;
    const auto &basis=trial.real_magnitude_trial->basis;
    if(basis.rows()!=trial.Phi.rows()||basis.cols()!=n||!basis.allFinite())return false;
    ++work.calls;work.q=basis;double logDetR=0.,previousLogOverlap=0.;
    const double guard=64.*n*std::numeric_limits<double>::epsilon();
    result=PureProjectorWeightResult{};result.status=PureProjectorWeightStatus::success;
    result.log_abs_weight=0.;result.complex_phase=DataType(1,0);result.z2_sign=0;
    result.overlap_rank=n;result.overlap_rcond=1.;result.overlap_residual=0.;
    for(int index=0;index<int(slices.size());++index) {
        const auto &s=slices[index];
        if(s.matrix.rows()!=basis.rows()||s.matrix.cols()!=basis.rows()||
           !std::isfinite(s.eta.real())||!std::isfinite(s.eta.imag()))return false;
        const auto *factor=work.factor(s.matrix,index,slices.size());
        if(!factor)return false;
        double detR=0.;bool qrOK=false;
        profileCore("sw_propagate_qr",[&] {
            work.next_q.noalias()=(*factor)*work.q;qrOK=work.orthogonalize(detR);
        });
        if(!qrOK||!std::isfinite(detR)||std::abs(detR)<=options.zero_tolerance*(1.+guard))return false;
        double logAbsDet=0.;bool overlapOK=true;
        profileCore("sw_overlap_svd",[&] {
            work.overlap.noalias()=basis.transpose()*work.next_q;
            if(!work.overlap.allFinite()||!work.overlapValues()||!(work.singular(0)>0.)) {
                overlapOK=false;return;
            }
            result.overlap_rcond=work.singular(n-1)/work.singular(0);
            if(result.overlap_rcond<=std::max(options.rank_tolerance,options.minimum_overlap_rcond)+guard) {
                overlapOK=false;return;
            }
            for(int i=0;i<n;++i) {
                if(!(work.singular(i)>0.)){overlapOK=false;return;}
                logAbsDet+=std::log(work.singular(i));
            }
            if(!std::isfinite(logAbsDet))overlapOK=false;
        });
        if(!overlapOK)return false;
        logDetR+=std::log(std::abs(detR));
        const double nextLogOverlap=logDetR+logAbsDet;
        const double identityResidual=s.identity_residual>=0.?s.identity_residual:
            (s.matrix-MatType::Identity(s.matrix.rows(),s.matrix.cols())).norm();
        const double logRatio=std::log(std::abs(s.eta))+
            (identityResidual<=10.*options.zero_tolerance?0.:.5*(nextLogOverlap-previousLogOverlap));
        if(!std::isfinite(logRatio)||logRatio<=std::log(options.zero_tolerance)+guard*n)return false;
        result.log_abs_weight+=logRatio;work.q.swap(work.next_q);previousLogOverlap=nextLogOverlap;
    }
    if(!std::isfinite(result.log_abs_weight))return false;
    result.weight=result.log_abs_weight<700.?DataType(std::exp(result.log_abs_weight),0.):DataType(1,0);
    result.determinant_identity_error=0.;return true;
}

// Validated deferred-phase path: magnitude-only reference weight during
// proposals; the signed phase is reconstructed at measurement.
//
// The production Metropolis test accepts on log_abs_matter alone (see
// PhononSampler::attempt/attemptBond); the Pfaffian complex phase and the
// center Green are consumed only at fixed measurement strides for phase
// reweighting. This routine reproduces the log|w| of
// pureProjectorStableReferenceWeight while omitting the phase-only work:
// no per-slice Cayley solve, no partial Green, no Pfaffian block, no final
// Green and no center-Green rebuild. It is the same QR propagation and the
// same rank/rcond certification gates; log|det overlap| is obtained from
// the SVD singular values (|det|=prod sigma), avoiding the determinant LU and
// the singular-vector computation. The phase value is NOT returned; callers
// must recompute the full signed weight at measurement time.
inline PureProjectorWeightResult pureProjectorMagnitudeWeight(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice>&slices,
    const PureFastOptions &options,PureStableWeightWorkspace *workspace=nullptr) {
    PureProjectorOptions greenOptions=pureFastGreenOptions(options);
    (void)greenOptions;
    PureProjectorWeightResult result;result.status=PureProjectorWeightStatus::success;
    result.log_abs_weight=0;result.complex_phase=DataType(1,0);result.z2_sign=0;
    result.overlap_rank=trial.Phi.cols();result.overlap_rcond=1;result.overlap_residual=0;
    PureStableWeightWorkspace localWorkspace;
    PureStableWeightWorkspace &work=workspace?*workspace:localWorkspace;
    if(options.real_magnitude&&options.accelerated_linear_algebra) {
        if(pureProjectorRealMagnitudeWeight(trial,slices,options,work.real_magnitude,result))return result;
        ++work.real_magnitude.fallbacks;
        result=PureProjectorWeightResult{};result.status=PureProjectorWeightStatus::success;
        result.log_abs_weight=0.;result.complex_phase=DataType(1,0);result.z2_sign=0;
        result.overlap_rank=trial.Phi.cols();result.overlap_rcond=1.;result.overlap_residual=0.;
    }
    work.q=trial.Phi;double logDetR=0,previousLogOverlap=0;
    for(int index=0;index<int(slices.size());++index){const auto&s=slices[index];
        if(s.matrix.rows()!=trial.Phi.rows()||s.matrix.cols()!=trial.Phi.rows()||
           (!s.validated_dense&&!pureProjectorMatrixFinite(s.matrix))||
           !std::isfinite(s.eta.real())||!std::isfinite(s.eta.imag())){
            result.status=PureProjectorWeightStatus::invalid_slice;result.first_failing_slice=index;return result;}
        DataType detR(1,0);bool propagation_failure=false;
        profileCore("sw_propagate_qr",[&]{
            pureProjectorPropagate(s.matrix,work.q,work.x,true);
            if(!pureProjectorMatrixFinite(work.x)){propagation_failure=true;return;}
            const int m=work.x.rows(),n=work.x.cols();
            const bool accelerated=options.accelerated_linear_algebra&&n>=32;
            if(accelerated){
                // Magnitude propagation consumes only |det R| and the Q
                // subspace.  MKL's no-pivot QR preserves that contract while
                // avoiding Eigen's temporary full thin-Q construction.
                work.next_q=work.x;
                work.qr_tau.resize(n);
                const lapack_int info=LAPACKE_zgeqrf(LAPACK_COL_MAJOR,m,n,
                    work.next_q.data(),m,work.qr_tau.data());
                if(info!=0){propagation_failure=true;return;}
                for(int i=0;i<n;++i)detR*=work.next_q(i,i);
                if(LAPACKE_zungqr(LAPACK_COL_MAJOR,m,n,n,work.next_q.data(),m,
                                  work.qr_tau.data())!=0){propagation_failure=true;return;}
            }else{
                work.qr.compute(work.x);
                work.identity.resize(work.x.rows(),work.x.cols());work.identity.setIdentity();
                work.next_q=work.qr.householderQ()*work.identity;
                const auto &matrixQR=work.qr.matrixQR();
                for(int i=0;i<work.x.cols();++i)detR*=matrixQR(i,i);
            }
        });
        if(propagation_failure){
            result.status=PureProjectorWeightStatus::propagation_failure;result.first_failing_slice=index;return result;}
        if(!std::isfinite(detR.real())||!std::isfinite(detR.imag())||std::abs(detR)<=options.zero_tolerance){
            result.status=PureProjectorWeightStatus::zero_weight;result.first_failing_slice=index;result.z2_sign=0;return result;}
        double nextLogOverlap=0;
        profileCore("sw_overlap_svd",[&]{
            logDetR+=std::log(std::abs(detR));
            work.overlap=trial.Phi.adjoint()*work.next_q;
            work.magnitude_svd.compute(work.overlap,options.accelerated_linear_algebra,
                options.rank_tolerance,options.minimum_overlap_rcond,options.zero_tolerance);
            result.overlap_rank=work.magnitude_svd.rank(options.rank_tolerance);
            if(result.overlap_rank!=work.overlap.rows()){result.status=PureProjectorWeightStatus::zero_weight;
                result.first_failing_slice=index;result.z2_sign=0;return;}
            const auto &sv=work.magnitude_svd.singularValues();
            double largest=sv(0),smallest=sv(work.overlap.rows()-1);
            result.overlap_rcond=largest>0?smallest/largest:0;
            if(result.overlap_rcond<options.minimum_overlap_rcond){
                result.status=PureProjectorWeightStatus::overlap_untrusted;result.first_failing_slice=index;return;}
            double logAbsDet=0;for(int i=0;i<sv.size();++i){
                if(!(sv(i)>0)){result.status=PureProjectorWeightStatus::zero_weight;
                    result.first_failing_slice=index;result.z2_sign=0;return;}
                logAbsDet+=std::log(sv(i));
            }
            // A small product of certified nonzero singular values is not a
            // zero overlap. The physical magnitude includes the QR scales.
            if(!std::isfinite(logAbsDet)){
                result.status=PureProjectorWeightStatus::overlap_untrusted;result.first_failing_slice=index;return;}
            nextLogOverlap=logDetR+logAbsDet;
        });
        if(result.status!=PureProjectorWeightStatus::success)return result;
        double ratioLog=std::log(std::abs(s.eta))+.5*(nextLogOverlap-previousLogOverlap);
        const double identityResidual=s.identity_residual>=0.0?s.identity_residual:
            (s.matrix-MatType::Identity(s.matrix.rows(),s.matrix.cols())).norm();
        double logRatio;
        if(identityResidual<=10*options.zero_tolerance) logRatio=std::log(std::abs(s.eta));
        else logRatio=ratioLog;
        if(!std::isfinite(logRatio)||logRatio<=std::log(options.zero_tolerance)){
            result.status=PureProjectorWeightStatus::zero_weight;result.first_failing_slice=index;result.z2_sign=0;return result;}
        result.log_abs_weight+=logRatio;
        work.q.swap(work.next_q);previousLogOverlap=nextLogOverlap;
    }
    result.weight=result.log_abs_weight<700?std::exp(result.log_abs_weight):DataType(1,0);
    result.determinant_identity_error=0;return result;
}

inline std::uint64_t pureFastConfigurationHash(const PureFastConfiguration &configuration) {    std::uint64_t hash=1469598103934665603ULL;
    auto mix=[&](std::uint64_t value){hash^=value;hash*=1099511628211ULL;};
    for(int field:configuration.hs_fields)mix(std::uint64_t(std::int64_t(field)));
    for(const auto &slice:configuration.slices)slice.matrix.mixHash(hash);
    return hash;
}

// Production local update. The Green and its rank/condition/solve certificate
// belong to the manager at this exact cut. The general subspace calculation
// below is used only when this specialization cannot certify the proposal.
inline PureFastRatioResult pureFastLocalRatioFromGreen(
    const GaussianTrialState &trial,const std::vector<PureProjectorSlice> &slices,
    int index,const MatType &newFactor,DataType newEta,const MatType &rightBefore,
    const PureProjectorGreenResult &trusted,const PureFastOptions &options,
    double proposalUniform=-1.,const PureFactorMatrix *cachedNewFactor=nullptr) {
    PureFastRatioResult result;
    const auto &storedOld=slices[index].matrix;
    MatType oldFactor;
    PureCompactTransform compact;
    if(options.stabilized_local_updates&&options.compact_proposal_algebra){
        if(options.cached_factor_support&&cachedNewFactor)
            compact=pureCompactTransform(storedOld,*cachedNewFactor);
        else {oldFactor=storedOld;compact=pureCompactTransform(oldFactor,newFactor);}
    }
    if(!compact.available||!options.stabilized_local_updates)oldFactor=storedOld;
    MatType transform;
    if(compact.available){
        if(!compact.ok)return result;
        transform=std::move(compact.transform);
    }else{
        PureFactorSolve inverseOld(oldFactor.transpose(),DataType(1,0));
        if(!inverseOld.isInvertible())return result;
        transform=inverseOld.solve(newFactor.transpose()).transpose();
    }
    const MatType delta=transform-MatType::Identity(transform.rows(),transform.cols());
    std::vector<int> active;
    if(compact.available){
        std::vector<int> small;
        if(!pureLocalSupport(delta,small,true)||(small.size()%2)!=0)return result;
        for(int i:small)active.push_back(compact.active[i]);
    }else if(!pureLocalSupport(delta,active,options.stabilized_local_updates)||
             (active.size()%2)!=0)return result;
    // Retain the rank threshold on deltaRight without an ambient D x L SVD.
    // Row extraction is exact because delta has no entries outside active.
    result.low_rank=int(active.size());
    if(!options.stabilized_local_updates&&!active.empty()){
        MatType deltaRows(active.size(),oldFactor.cols());
        for(int r=0;r<int(active.size());++r)
            deltaRows.row(r)=newFactor.row(active[r])-oldFactor.row(active[r]);
        const MatType deltaRight=deltaRows*rightBefore;
        Eigen::JacobiSVD<MatType> rankSvd(deltaRight);
        rankSvd.setThreshold(options.rank_tolerance);
        result.low_rank=rankSvd.rank();
        // General/near-rank-deficient updates retain the old truncation path.
        if(result.low_rank!=int(active.size()))return result;
    }
    const bool defer=options.stabilized_local_updates&&options.defer_rejected_green&&
        proposalUniform>=0.&&proposalUniform<1.;
    PureLocalGreenResult local;
    if(compact.available){
        std::vector<int> small;
        for(int a:active)small.push_back(int(std::lower_bound(compact.active.begin(),compact.active.end(),a)-compact.active.begin()));
        local=pureLocalGreenUpdateFromBlock(pureLocalBlock(delta,small),active,trusted.green,
                                            options.residual_tolerance,true,defer);
    }else local=pureLocalGreenUpdate(transform,trusted.green,options.residual_tolerance,
                                     options.stabilized_local_updates,defer);
    if(!local.ok)return result;
    if(defer&&!pureLocalGreenFiniteWithoutMaterializing(trusted.green,local)){
        pureMaterializeLocalGreen(trusted.green,local);
        if(!local.ok)return result;
    }
    const DataType etaRatio=newEta/slices[index].eta;
    const DataType squared=etaRatio*etaRatio*local.determinant;
    if(!std::isfinite(squared.real())||!std::isfinite(squared.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    if(!options.stabilized_local_updates&&std::abs(squared)<=options.zero_tolerance){
        result.status=PureFastRatioStatus::zero_ratio;return result;
    }
    if(squared==DataType(0,0)){
        // A floating-point zero is not a certified physical zero. Route it
        // through the same-proposal fresh rebuild and fail-closed path.
        result.ratio=DataType(0,0);result.status=PureFastRatioStatus::zero_ratio;
        return result;
    }
    MatType cayley;
    if(compact.available){
        if(transform.rows()){
            PureAmbientBlockSolve inverseCayley(MatType::Identity(transform.rows(),transform.cols())+transform,
                                                DataType(2,0),storedOld.rows());
            if(!inverseCayley.isInvertible())return result;
            cayley=inverseCayley.solve(2.0*MatType::Identity(transform.rows(),transform.cols()))-
                MatType::Identity(transform.rows(),transform.cols());
        }else cayley.resize(0,0);
    }else{
        PureFactorSolve inverseCayley(MatType::Identity(transform.rows(),transform.cols())+transform,DataType(2,0));
        if(!inverseCayley.isInvertible())return result;
        cayley=inverseCayley.solve(2.0*MatType::Identity(transform.rows(),transform.cols()))-
            MatType::Identity(transform.rows(),transform.cols());
    }
    const PfaffianResult phase=compact.available?
        pureCompactProductSign(cayley,trusted.green,compact.active):
        pureLocalProductSign(cayley,trusted.green,result.used_local_support);
    if(compact.available)result.used_local_support=true;
    result.pfaffian_status=phase.status;
    if(!phase.ok()){
        result.status=PureFastRatioStatus::pfaffian_untrusted;return result;
    }
    const DataType extraSign=(trial.Phi.cols()%2)?DataType(-1,0):DataType(1,0);
    const DataType target=pureProjectorUnitPhase(etaRatio*phase.value*extraSign);
    const double phaseTolerance=options.stabilized_local_updates?
        options.local_reality_tolerance:options.reality_tolerance;
    if(std::abs(target*target-pureProjectorUnitPhase(squared))>2.*phaseTolerance){
        result.status=PureFastRatioStatus::pfaffian_untrusted;return result;
    }
    DataType ratio=std::sqrt(squared);
    if((std::conj(target)*ratio).real()<0)ratio=-ratio;
    if(!std::isfinite(ratio.real())||!std::isfinite(ratio.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    const double realityTolerance=options.stabilized_local_updates?
        options.local_reality_tolerance:options.reality_tolerance;
    if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
       std::abs(ratio.imag())>realityTolerance*std::max(1.0,std::abs(ratio.real()))){
        result.status=PureFastRatioStatus::complex_ratio;return result;
    }
    if(!options.stabilized_local_updates){
        double residual=pureLocalGreenResidualBound(trusted,local);
        if(!std::isfinite(residual)||residual>options.green_update_tolerance)
            residual=(local.green*local.green-
                MatType::Identity(local.green.rows(),local.green.cols())).norm()/
                std::max(1.0,local.green.norm());
        if(!std::isfinite(residual)||residual>options.green_update_tolerance){
            result.status=PureFastRatioStatus::overlap_untrusted;return result;
        }
    }
    const DataType decisionRatio=options.weight_mode==PureProjectorWeightMode::RealZ2?DataType(ratio.real(),0):ratio;
    if(defer&&proposalUniform<std::min(1.,std::abs(decisionRatio))&&local.green.size()==0){
        pureMaterializeLocalGreen(trusted.green,local);
        if(!local.ok)return result;
    }
    result.ratio=ratio;result.fast_updated_green=std::move(local.green);
    result.overlap_rcond=trusted.overlap_rcond;
    result.solve_residual=std::max(trusted.solve_residual,local.solve_residual);
    result.used_trusted_green=true;result.used_green_update=result.fast_updated_green.size()!=0;
    result.status=PureFastRatioStatus::success;
    return result;
}

inline PureFastRatioResult pureFastLocalRatioAtCut(
    const GaussianTrialState &trial,
    const std::vector<PureProjectorSlice> &slices,
    int index,const MatType &newFactor,DataType newEta,
    const MatType &rightBefore,const MatType &leftAfter,
    const PureFastOptions &options,
    const PureProjectorGreenResult *trustedGreen = nullptr) {
    PureFastRatioResult result;
    PureProjectorOptions localGreenOptions;
    localGreenOptions.accelerated_linear_algebra=options.accelerated_linear_algebra;
    if(index<0||index>=int(slices.size())||newFactor.rows()!=trial.Phi.rows()||
       newFactor.cols()!=trial.Phi.rows()||!pureProjectorMatrixFinite(newFactor))return result;
    if(rightBefore.rows()!=trial.Phi.rows()||leftAfter.rows()!=trial.Phi.rows()||
       rightBefore.cols()!=trial.Phi.cols()||leftAfter.cols()!=trial.Phi.cols())return result;
    // Only the live manager may supply this certificate for exactly index+1.
    // Keep the original local Green thresholds, including its stricter
    // default purity tolerance, even if the stack uses different options.
    const bool canReuse=options.reuse_trusted_green&&trustedGreen&&trustedGreen->ok()&&
        trustedGreen->green.rows()==trial.Phi.rows()&&
        trustedGreen->green.cols()==trial.Phi.rows()&&
        trustedGreen->overlap_rank==trial.Phi.cols()&&
        trustedGreen->overlap_rcond>=localGreenOptions.minimum_overlap_rcond&&
        trustedGreen->overlap_rcond>=options.minimum_overlap_rcond&&
        trustedGreen->overlap_rcond>options.rank_tolerance&&
        std::isfinite(trustedGreen->overlap_rcond)&&
        trustedGreen->solve_residual<=localGreenOptions.solve_residual_tolerance&&
        trustedGreen->solve_residual<=options.residual_tolerance&&
        std::isfinite(trustedGreen->solve_residual)&&
        trustedGreen->green_residual<=localGreenOptions.green_residual_tolerance&&
        std::isfinite(trustedGreen->green_residual)&&
        pureProjectorMatrixFinite(trustedGreen->green);
    if(canReuse&&options.compact_local_update){
        PureFastRatioResult local=pureFastLocalRatioFromGreen(
            trial,slices,index,newFactor,newEta,rightBefore,*trustedGreen,options);
        if(local.ok())return local;
        // Unsupported or untrusted acceleration uses the original double
        // reference; its failure still routes to the unchanged MP protocol.
    }
    const MatType oldFactor=slices[index].matrix;
    const MatType oldRight=oldFactor*rightBefore;
    const MatType deltaRight=(newFactor-oldFactor)*rightBefore;
    Eigen::JacobiSVD<MatType> deltaSvd(deltaRight,Eigen::ComputeThinU|Eigen::ComputeThinV);
    deltaSvd.setThreshold(options.rank_tolerance);
    const int rank=deltaSvd.rank();
    result.low_rank=rank;
    if(rank==0){
        result.ratio=newEta/slices[index].eta;
        if(!std::isfinite(result.ratio.real())||!std::isfinite(result.ratio.imag())){
            result.status=PureFastRatioStatus::nonfinite_ratio;return result;
        }
        if(std::abs(result.ratio)<=options.zero_tolerance){
            result.status=PureFastRatioStatus::zero_ratio;return result;
        }
        if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
           std::abs(result.ratio.imag())>options.reality_tolerance*
               std::max(1.0,std::abs(result.ratio.real()))){
            result.status=PureFastRatioStatus::complex_ratio;return result;
        }
        const PureProjectorGreenResult unchangedGreen=
            canReuse?*trustedGreen:pureProjectorGreenThinQr(oldRight,leftAfter,localGreenOptions);
        if(!unchangedGreen.ok()){
            result.status=PureFastRatioStatus::overlap_untrusted;return result;
        }
        result.status=PureFastRatioStatus::success;
        result.used_trusted_green=canReuse;
        result.fast_updated_green=unchangedGreen.green;
        result.overlap_rcond=unchangedGreen.overlap_rcond;
        result.solve_residual=unchangedGreen.solve_residual;return result;
    }
    MatType u=deltaSvd.matrixU().leftCols(rank);
    MatType c=deltaSvd.singularValues().head(rank).asDiagonal()*
              deltaSvd.matrixV().leftCols(rank).adjoint();
    const MatType overlap=leftAfter.adjoint()*oldRight;
    Eigen::JacobiSVD<MatType> overlapSvd(overlap,Eigen::ComputeThinU|Eigen::ComputeThinV);
    overlapSvd.setThreshold(options.rank_tolerance);
    if(overlapSvd.rank()!=overlap.rows()){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const double largest=overlapSvd.singularValues()(0);
    const double smallest=overlapSvd.singularValues()(overlap.rows()-1);
    result.overlap_rcond=largest>0?smallest/largest:0;
    if(result.overlap_rcond<options.minimum_overlap_rcond){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const MatType a=leftAfter.adjoint()*u;
    const MatType x=overlapSvd.solve(a);
    const MatType y=overlapSvd.solve(leftAfter.adjoint());
    result.solve_residual=std::max(
        (overlap*x-a).norm()/std::max(1.0,a.norm()),
        (overlap*y-leftAfter.adjoint()).norm()/std::max(1.0,leftAfter.norm()));
    if(!std::isfinite(result.solve_residual)||result.solve_residual>options.residual_tolerance){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const MatType k=MatType::Identity(rank,rank)+c*x;
    Eigen::FullPivLU<MatType> klu(k);
    if(!klu.isInvertible()){result.status=PureFastRatioStatus::zero_ratio;return result;}
    const DataType etaRatio=newEta/slices[index].eta;
    const DataType ratioSquared=etaRatio*etaRatio*k.determinant();
    if(!std::isfinite(ratioSquared.real())||!std::isfinite(ratioSquared.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    if(std::abs(ratioSquared)<=options.zero_tolerance){
        result.status=PureFastRatioStatus::zero_ratio;return result;
    }
    DataType ratio=std::sqrt(ratioSquared);

    PureFactorSolve oldTranspose(oldFactor.transpose(),DataType(1,0));
    if(!oldTranspose.isInvertible()){
        result.status=PureFastRatioStatus::rank_failure;return result;
    }
    MatType transform=oldTranspose.solve(newFactor.transpose()).transpose();
    MatType onePlus=MatType::Identity(transform.rows(),transform.cols())+transform;
    PureFactorSolve cayleySolve(onePlus,DataType(2,0));
    if(!cayleySolve.isInvertible()){
        result.status=PureFastRatioStatus::pfaffian_untrusted;return result;
    }
    MatType transformGreen=cayleySolve.solve(2.0*MatType::Identity(
        transform.rows(),transform.cols()))-MatType::Identity(transform.rows(),transform.cols());
    PureProjectorGreenResult oldGreen=canReuse?*trustedGreen:
        pureProjectorGreenThinQr(oldRight,leftAfter,localGreenOptions);
    if(!oldGreen.ok()){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    result.used_trusted_green=canReuse;
    bool usedSupport=false;
    PfaffianResult pf=pureLocalProductSign(transformGreen,oldGreen.green,usedSupport);
    result.pfaffian_status=pf.status;
    if(!pf.ok()){result.status=PureFastRatioStatus::pfaffian_untrusted;return result;}
    const DataType extraSign=((trial.Phi.cols()%2)==0)?DataType(1.0):DataType(-1.0);
    const DataType target=pureProjectorUnitPhase(etaRatio*pf.value*extraSign);
    if(std::abs(target*target-pureProjectorUnitPhase(ratioSquared))>2.*options.reality_tolerance){
        result.status=PureFastRatioStatus::pfaffian_untrusted;return result;
    }
    if((std::conj(target)*ratio).real()<0)ratio=-ratio;
    if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
       std::abs(ratio.imag())>options.reality_tolerance*std::max(1.0,std::abs(ratio.real()))){
        result.status=PureFastRatioStatus::complex_ratio;return result;
    }
    if(!std::isfinite(ratio.real())||!std::isfinite(ratio.imag())){
        result.status=PureFastRatioStatus::nonfinite_ratio;return result;
    }
    result.ratio=ratio;
    const MatType z=klu.solve(c*y);
    result.fast_updated_green=oldGreen.green-2.0*(u-oldRight*x)*z;
    if(!pureProjectorMatrixFinite(result.fast_updated_green)){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    const double greenResidual=(result.fast_updated_green*result.fast_updated_green-
        MatType::Identity(result.fast_updated_green.rows(),result.fast_updated_green.cols())).norm()/
        std::max(1.0,result.fast_updated_green.norm());
    if(!std::isfinite(greenResidual)||greenResidual>options.green_update_tolerance){
        result.status=PureFastRatioStatus::overlap_untrusted;return result;
    }
    result.status=PureFastRatioStatus::success;
    return result;
}

inline PureFastRatioResult pureFastLocalRatio(
    const GaussianTrialState &trial,
    const std::vector<PureProjectorSlice> &slices,
    int index,const MatType &newFactor,DataType newEta,
    const PureFastOptions &options) {
    PureFastRatioResult result;
    PureProjectorStackManager stack(trial,slices,std::max(1,int(slices.size())));
    if(!stack.moveToCut(index)){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    const MatType rightBefore=stack.phiRight();
    if(!stack.moveToCut(index+1)){result.status=PureFastRatioStatus::overlap_untrusted;return result;}
    return pureFastLocalRatioAtCut(
        trial,slices,index,newFactor,newEta,rightBefore,stack.phiLeft(),options);
}

// Initialize an independent bra/ket contour without asking every unphysical
// trial-state prefix to have a well-conditioned overlap.  The physical sign is
// anchored at the mirrored U^dagger U contour and transported to the requested
// bra fields with the same local-ratio and same-proposal MP trust gates used by
// production updates.  No Metropolis decision or random number is consumed.
inline PureEndpointTransportInitializationResult
pureProjectorEndpointTransportInitialization(
    const GaussianTrialState &trial,const PureFastConfiguration &target,
    int blockSize,const PureFastOptions &options) {
    PureEndpointTransportInitializationResult result;
    const int total=int(target.slices.size());
    if(options.weight_mode!=PureProjectorWeightMode::RealZ2||
       total<=0||(total&1)||blockSize<=0||
       target.hs_fields.size()!=target.slices.size()||
       target.locations.size()!=target.slices.size()) {
        result.message="endpoint transport requires a balanced bra/ket contour";
        return result;
    }
    const int center=total/2;
    PureFastConfiguration ket;
    ket.slices.assign(target.slices.begin(),target.slices.begin()+center);
    ket.hs_fields.assign(target.hs_fields.begin(),target.hs_fields.begin()+center);
    ket.locations.assign(target.locations.begin(),target.locations.begin()+center);
    for(int index=0;index<center;++index)if(ket.locations[index].branch!=PureBranch::Ket){
        result.message="endpoint transport found a non-ket factor before the center";
        return result;
    }
    const PureMirroredInitializationResult mirrored=
        pureProjectorMirroredConfiguration(ket,options.reality_tolerance);
    if(!mirrored.ok()){
        result.message=std::string("mirrored anchor failed: ")+mirrored.message;
        return result;
    }
    for(int offset=0;offset<center;++offset){
        const int index=center+offset;
        const int source=center-1-offset;
        const auto &bra=target.locations[index];
        const auto &ketLocation=target.locations[source];
        if(bra.branch!=PureBranch::Bra||ketLocation.branch!=PureBranch::Ket||
           bra.slice!=ketLocation.slice||bra.factor!=ketLocation.factor||
           bra.bond!=ketLocation.bond||bra.aux!=ketLocation.aux){
            result.message="endpoint transport contour is not strict reverse-adjoint ordered";
            return result;
        }
    }

    PureProjectorOptions endpointOptions;
    endpointOptions.rank_tolerance=options.rank_tolerance;
    endpointOptions.minimum_overlap_rcond=options.minimum_overlap_rcond;
    endpointOptions.solve_residual_tolerance=options.residual_tolerance;
    endpointOptions.green_residual_tolerance=options.green_update_tolerance;
    PureFastConfiguration current=mirrored.configuration;
    const PureEndpointRebuildResult anchor=pureProjectorEndpointRebuild(
        trial,current.slices,center,blockSize,endpointOptions);
    if(!anchor.ok()){
        result.status=PureProjectorWeightStatus::overlap_untrusted;
        result.message=std::string("mirrored endpoint anchor failed: ")+anchor.message;
        return result;
    }
    result.mirrored_log_abs_weight=anchor.log_abs_weight;
    double transportedLogAbs=anchor.log_abs_weight;
    DataType transportedPhase(1.0,0.0);
    int transportedZ2=1;
    for(int index=center;index<total;++index){
        const double matrixDifference=(current.slices[index].matrix-
            target.slices[index].matrix).norm();
        const double etaDifference=std::abs(current.slices[index].eta-
            target.slices[index].eta);
        if(matrixDifference<=options.zero_tolerance&&
           etaDifference<=options.zero_tolerance)continue;
        DataType ratio(0.0,0.0);
        PureProjectorStackManager stack(
            trial,current.slices,blockSize,endpointOptions,index);
        PureFastRatioResult fast;
        if(stack.ok()){
            const MatType rightBefore=stack.phiRight();
            if(stack.moveToCut(index+1))fast=pureFastLocalRatioAtCut(
                trial,current.slices,index,target.slices[index].matrix,
                target.slices[index].eta,rightBefore,stack.phiLeft(),options);
        }
        if(fast.ok())ratio=fast.ratio;
        else{
            ++result.mp_fallback_count;
            PureMpProposalOptions mpOptions;
            mpOptions.minimum_endpoint_rcond=options.minimum_overlap_rcond;
            mpOptions.residual_tolerance=options.residual_tolerance;
            mpOptions.reality_tolerance=options.reality_tolerance;
            mpOptions.zero_tolerance=options.zero_tolerance;
            const PureMpProposalResult mp=pureProjectorMpSameProposal(
                trial,current.slices,index,target.slices[index].matrix,
                target.slices[index].eta,mpOptions,nullptr,
                pureFastConfigurationHash(current));
            if(!mp.ok()){
                ++result.mp_failure_count;
                result.status=PureProjectorWeightStatus::overlap_untrusted;
                result.message=std::string("endpoint sign transport MP gate failed: ")+
                    mp.message;
                return result;
            }
            ratio=mp.ratio;
        }
        if(!std::isfinite(ratio.real())||!std::isfinite(ratio.imag())||
           std::abs(ratio)<=options.zero_tolerance){
            result.status=PureProjectorWeightStatus::zero_weight;
            result.message="endpoint sign transport produced a zero/nonfinite ratio";
            return result;
        }
        if(options.weight_mode==PureProjectorWeightMode::RealZ2&&
           std::abs(ratio.imag())>options.reality_tolerance*
                std::max(1.0,std::abs(ratio.real()))){
            result.status=PureProjectorWeightStatus::complex_ratio;
            result.message="endpoint sign transport produced a complex physical ratio";
            return result;
        }
        transportedLogAbs+=std::log(std::abs(ratio));
        transportedPhase=pureProjectorUnitPhase(transportedPhase*ratio);
        transportedZ2*=ratio.real()>=0.0?1:-1;
        current.slices[index]=target.slices[index];
        current.hs_fields[index]=target.hs_fields[index];
        ++result.transported_factors;
    }
    const PureEndpointRebuildResult endpoint=pureProjectorEndpointRebuild(
        trial,target.slices,center,blockSize,endpointOptions);
    if(!endpoint.ok()){
        result.status=PureProjectorWeightStatus::overlap_untrusted;
        result.message=std::string("target endpoint rebuild failed: ")+endpoint.message;
        return result;
    }
    result.transported_log_abs_weight=transportedLogAbs;
    result.endpoint_log_abs_weight=endpoint.log_abs_weight;
    result.log_abs_consistency_error=std::abs(
        transportedLogAbs-endpoint.log_abs_weight);
    if(!std::isfinite(result.log_abs_consistency_error)||
       result.log_abs_consistency_error>
           std::max(options.residual_tolerance,
                    options.green_update_tolerance)){
        result.status=PureProjectorWeightStatus::determinant_identity_failure;
        result.message="transported and endpoint log-weights disagree";
        return result;
    }
    result.weight.status=PureProjectorWeightStatus::success;
    result.weight.log_abs_weight=endpoint.log_abs_weight;
    result.weight.complex_phase=transportedPhase;
    result.weight.z2_sign=transportedZ2;
    result.weight.green=endpoint.green;
    result.weight.overlap_rank=endpoint.overlap_rank;
    result.weight.overlap_rcond=endpoint.overlap_rcond;
    result.weight.overlap_residual=endpoint.solve_residual;
    result.weight.green_residual=endpoint.green_residual;
    result.weight.weight=endpoint.log_abs_weight<700.0?
        std::exp(endpoint.log_abs_weight)*transportedPhase:transportedPhase;
    result.status=PureProjectorWeightStatus::success;
    result.message="mirrored sign anchor, trusted local transport, and endpoint QR passed";
    return result;
}

class PureProjectorFastWalker {
public:
    PureProjectorFastWalker(const GaussianTrialState &trial,
                            PureFastConfiguration configuration,int blockSize,
                            PureFastRunMode mode,
                            PureFastOptions options=PureFastOptions(),
                            PureFastInitializationPolicy initializationPolicy=
                                PureFastInitializationPolicy::SequentialAudit)
        : trial_(trial),configuration_(std::move(configuration)),block_size_(blockSize),
          mode_(mode),options_(options),initialization_policy_(initializationPolicy) {
        if(options_.stabilized_local_updates&&
           (mode_!=PureFastRunMode::FastStrict||
            options_.weight_mode!=PureProjectorWeightMode::RealZ2||
            options_.read_only_audit_interval!=0))
            throw std::invalid_argument("stabilized local policy requires RealZ2 FastStrict without audits");
        if(configuration_.slices.size()!=configuration_.hs_fields.size()||
           configuration_.slices.size()!=configuration_.locations.size())
            throw std::invalid_argument("fast walker configuration size mismatch");
        if(initialization_policy_==PureFastInitializationPolicy::MirroredTheoremZ2Plus){
            const int size=int(configuration_.slices.size());
            if(size==0||size%2)throw std::runtime_error("mirrored initializer requires paired factors");
            for(int i=0;i<size/2;++i){
                const auto &a=configuration_.slices[i],&b=configuration_.slices[size-1-i];
                if(a.eta!=std::conj(b.eta)||a.matrix.rows()!=b.matrix.cols()||a.matrix.cols()!=b.matrix.rows())
                    throw std::runtime_error("mirrored initializer factors do not prove a norm");
                const MatType left=a.matrix,right=b.matrix;
                if(!(left-right.adjoint()).isZero(0.))
                    throw std::runtime_error("mirrored initializer factors do not prove a norm");
            }
            DataType scalar(1.0,0.0);for(const auto&s:configuration_.slices)scalar*=s.eta;
            if(!std::isfinite(scalar.real())||!std::isfinite(scalar.imag())||
               !(scalar.real()>0.0)||std::abs(scalar.imag())>options_.reality_tolerance*
                    std::max(1.0,std::abs(scalar.real())))
                throw std::runtime_error("mirrored initializer scalar-prefactor theorem failed");
            const int center=int(configuration_.slices.size()/2);
            PureProjectorOptions endpointOptions;
            endpointOptions.rank_tolerance=options_.rank_tolerance;
            endpointOptions.minimum_overlap_rcond=options_.minimum_overlap_rcond;
            endpointOptions.solve_residual_tolerance=options_.residual_tolerance;
            endpointOptions.green_residual_tolerance=options_.stabilized_local_updates?
                options_.local_green_rebuild_tolerance:options_.green_update_tolerance;
            const PureEndpointRebuildResult endpoint=pureProjectorEndpointRebuild(
                trial_,configuration_.slices,center,block_size_,endpointOptions);
            if(!endpoint.ok())throw std::runtime_error(
                std::string("mirrored endpoint initializer failed: ")+endpoint.message);
            current_.status=PureProjectorWeightStatus::success;
            current_.log_abs_weight=endpoint.log_abs_weight;
            current_.complex_phase=DataType(1.0,0.0);
            current_.z2_sign=1;current_.green=endpoint.green;
            current_.overlap_rank=endpoint.overlap_rank;
            current_.overlap_rcond=endpoint.overlap_rcond;
            current_.overlap_residual=endpoint.solve_residual;
            current_.green_residual=endpoint.green_residual;
            diagnostics_.maximum_endpoint_rebuild_green_residual=endpoint.green_residual;
            current_.weight=current_.log_abs_weight<700.0?
                DataType(std::exp(current_.log_abs_weight),0.0):DataType(1.0,0.0);
            z2_sign_=1;
        }else if(initialization_policy_==PureFastInitializationPolicy::EndpointTransportedZ2){
            const PureEndpointTransportInitializationResult initialized=
                pureProjectorEndpointTransportInitialization(
                    trial_,configuration_,block_size_,options_);
            if(!initialized.ok())throw std::runtime_error(
                std::string("endpoint-transported initializer failed: ")+
                initialized.message);
            current_=initialized.weight;
            z2_sign_=current_.z2_sign;
            diagnostics_.mp_fallback_count+=initialized.mp_fallback_count;
            diagnostics_.mp_fallback_failure_count+=initialized.mp_failure_count;
            diagnostics_.maximum_endpoint_rebuild_green_residual=
                current_.green_residual;
        }else{
            current_=pureProjectorStableReferenceWeight(trial_,configuration_.slices,options_);
            if(!current_.ok())throw std::runtime_error("fast walker initial weight failed");
            if(options_.weight_mode==PureProjectorWeightMode::RealZ2&&current_.z2_sign==0)
                throw std::runtime_error("RealZ2 initial reference sign unresolved");
            z2_sign_=current_.z2_sign;
        }
        stack_options_.accelerated_linear_algebra=options_.accelerated_linear_algebra;
        stack_options_.pair_block_factors=options_.pair_block_factors;
        stack_options_.direct_green_columns=options_.direct_green_columns;
        stack_options_.reuse_rebuild_workspace=options_.reuse_rebuild_workspace;
        stack_options_.rank_tolerance=options_.rank_tolerance;
        stack_options_.minimum_overlap_rcond=options_.minimum_overlap_rcond;
        stack_options_.solve_residual_tolerance=options_.residual_tolerance;
        stack_options_.green_residual_tolerance=options_.stabilized_local_updates?
            options_.local_green_rebuild_tolerance:options_.green_update_tolerance;
        const int initialCut=initialization_policy_!=
            PureFastInitializationPolicy::SequentialAudit?
            int(configuration_.slices.size()/2):0;
        manager_.reset(new PureProjectorStackManager(
            trial_,configuration_.slices,block_size_,stack_options_,initialCut,
            options_.stabilized_local_updates,true,options_.reuse_stabilized_subspaces,
            options_.cached_factor_support,options_.borrow_green,options_.reuse_green_buffers,
            options_.reuse_rebuild_result));
        if(!manager_->ok())throw std::runtime_error("fast walker initial stack failed");
    }

    std::uint64_t configurationHash() const {
        if(!configuration_hash_valid_){
            configuration_hash_=pureFastConfigurationHash(configuration_);
            configuration_hash_valid_=true;
        }
        return configuration_hash_;
    }
    PureProjectorGreenResult measurementGreen() {
        if(options_.stabilized_local_updates&&options_.reuse_measurement_green){
            auto reused=manager_->readOnlyGreenAtCut(int(configuration_.slices.size()/2));
            if(reused.ok())return reused;
        }
        // Measurement is a read-only center-cut rebuild.  It must not move the
        // live sweep cut, otherwise the next local proposal would require a
        // long inverse propagation and accumulate avoidable round-off.
        PureProjectorStackManager measurement(
            trial_,configuration_.slices,block_size_,stack_options_,
            int(configuration_.slices.size()/2),options_.stabilized_local_updates,
            !options_.stabilized_local_updates||options_.retain_measurement_rebuild_cache,
            options_.reuse_stabilized_subspaces,options_.cached_factor_support,
            options_.borrow_green,options_.reuse_green_buffers,options_.reuse_rebuild_result);
        measurement_scaled_green_rebuilds_+=measurement.scaledGreenRebuilds();
        if(!measurement.ok())
            return PureProjectorGreenResult();
        if(options_.stabilized_local_updates&&options_.borrow_green)
            return std::move(measurement).takeGreen();
        return measurement.green();
    }
    int z2Sign() const {return z2_sign_;}
    std::uint64_t scaledGreenRebuilds() const {
        return measurement_scaled_green_rebuilds_+(manager_?manager_->scaledGreenRebuilds():0);
    }
    std::string stackFailureDetail() const {
        return manager_?manager_->failureDetail():"stack_status=unavailable";
    }
    std::string initializationPolicy() const {
        if(initialization_policy_==PureFastInitializationPolicy::MirroredTheoremZ2Plus)
            return "mirrored_theorem_z2_plus";
        if(initialization_policy_==PureFastInitializationPolicy::EndpointTransportedZ2)
            return "endpoint_transported_z2";
        return "sequential_audit";
    }
    const PureFastConfiguration &configuration() const{return configuration_;}
    const PureProjectorWeightResult &currentWeight() const{return current_;}
    const PureFastDiagnostics &diagnostics() const{return diagnostics_;}

#ifdef PFQMC_PURE_PROJECTOR_TEST_HOOKS
    // Test-only checkpoint restoration. Production builds never expose this
    // hook; the restored state must already have passed an endpoint rebuild.
    void restoreTrustedCheckpointForTest(int z2,double logAbsWeight,DataType phase) {
        if(z2!=1&&z2!=-1)throw std::invalid_argument("invalid checkpoint Z2");
        const int center=int(configuration_.slices.size()/2);
        manager_.reset(new PureProjectorStackManager(
            trial_,configuration_.slices,block_size_,stack_options_,center));
        if(!manager_->ok())throw std::runtime_error("checkpoint endpoint rebuild failed");
        const PureProjectorGreenResult rebuilt=manager_->green();
        if(!rebuilt.ok())throw std::runtime_error("checkpoint Green rebuild failed");
        z2_sign_=z2;current_.status=PureProjectorWeightStatus::success;
        current_.log_abs_weight=logAbsWeight;current_.complex_phase=phase;
        current_.z2_sign=z2;current_.green=rebuilt.green;
        current_.overlap_rank=rebuilt.overlap_rank;
        current_.overlap_rcond=rebuilt.overlap_rcond;
        current_.overlap_residual=rebuilt.solve_residual;
        current_.green_residual=rebuilt.green_residual;
        current_.weight=logAbsWeight<700.0?
            std::exp(logAbsWeight)*phase:phase;
        accepted_since_rebuild_=0;terminated_=false;
    }
    void setDecisionMarginForTest(double tolerance) {
        if(!std::isfinite(tolerance)||tolerance<0)throw std::invalid_argument("invalid test margin");
        options_.decision_margin_tolerance=tolerance;
    }
#endif

    PureFastProposalResult propose(const PureFastProposal &proposal) {
        if(options_.stabilized_local_updates)return proposeLocal(proposal);
        if(proposal.immutable_factor_only&&proposal.new_factor.size()==0&&proposal.cached_factor)
            return propose(proposal.materialized());
        PureFastProposalResult output;
        if(terminated_){output.terminated=true;output.ratio.status=PureFastRatioStatus::terminated;
            output.failure_stage="walker_already_terminated";return output;}
        if(proposal.index<0||proposal.index>=int(configuration_.slices.size())||
           (proposal.new_hs!=-1&&proposal.new_hs!=1)||proposal.uniform<0||proposal.uniform>=1){
            output.ratio.status=PureFastRatioStatus::invalid_proposal;return output;
        }
        const int originalCut=manager_->cut();
        const std::uint64_t proposalNumber=diagnostics_.proposals++;
        const auto &location=configuration_.locations[proposal.index];
        output.snapshot.branch=location.branch;output.snapshot.slice=location.slice;
        output.snapshot.factor=location.factor;output.snapshot.bond=location.bond;
        output.snapshot.aux=location.aux;output.snapshot.old_hs=configuration_.hs_fields[proposal.index];
        output.snapshot.new_hs=proposal.new_hs;output.snapshot.configuration_hash=configurationHash();
        output.snapshot.uniform=proposal.uniform;
        const auto fastStart=std::chrono::steady_clock::now();
        auto moveToTrustedCut=[&](int target){
            if(manager_->moveToCut(target))return true;
            // An intermediate live cut may be ill conditioned even when the
            // requested target cut is sound. Rebuild the target directly so
            // intermediate gauges do not become artificial sweep barriers.
            if(manager_->rebuildToCut(target)){
                ++diagnostics_.cut_rebuild_recovery_count;
                ++diagnostics_.rebuild_count;
                accepted_since_rebuild_=0;
                return true;
            }
            // The caller will route this unresolved cut to the independent
            // same-proposal reference; it is not itself a terminal failure.
            ++diagnostics_.cut_rebuild_reference_count;
            return false;
        };
        MatType rightBeforeForUpdate;
        if(!moveToTrustedCut(proposal.index)){
            output.ratio.status=PureFastRatioStatus::overlap_untrusted;
        }else{
            const MatType rightBefore=manager_->phiRight();
            rightBeforeForUpdate=rightBefore;
            if(!moveToTrustedCut(proposal.index+1))
                output.ratio.status=PureFastRatioStatus::overlap_untrusted;
            else {
                const PureProjectorGreenResult trusted=manager_->green();
                output.ratio=pureFastLocalRatioAtCut(
                    trial_,configuration_.slices,proposal.index,proposal.new_factor,
                    proposal.new_eta,rightBefore,manager_->phiLeft(),options_,&trusted);
            }
        }
        checkPositiveRatio(proposal,output.ratio);
        output.fast_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-fastStart).count();
        diagnostics_.total_fast_seconds+=output.fast_seconds;
        output.snapshot.fast_ratio=output.ratio.ratio;
        if(std::isfinite(output.ratio.overlap_rcond)&&output.ratio.overlap_rcond>0)
            diagnostics_.minimum_overlap_rcond=std::min(
                diagnostics_.minimum_overlap_rcond,output.ratio.overlap_rcond);

        // Materialize a second configuration only for reference/recovery paths.
        // The ordinary trusted proposal changes exactly one factor on acceptance.
        std::unique_ptr<PureFastConfiguration> candidate;
        auto candidateConfiguration=[&]() -> PureFastConfiguration & {
            if(!candidate){
                candidate.reset(new PureFastConfiguration(configuration_));
                candidate->slices[proposal.index].matrix=proposal.new_factor;
                candidate->slices[proposal.index].eta=proposal.new_eta;
                candidate->hs_fields[proposal.index]=proposal.new_hs;
            }
            return *candidate;
        };
        const double fastProbability=output.ratio.ok()?std::min(1.0,std::abs(output.ratio.ratio)):0.0;
        output.ratio.trust_alarm=!output.ratio.ok()||
            std::abs(proposal.uniform-fastProbability)<options_.decision_margin_tolerance;
        if(output.ratio.trust_alarm)++diagnostics_.trust_alarm_count;
        const bool periodicAudit=options_.read_only_audit_interval>0&&
            ((proposalNumber+1)%options_.read_only_audit_interval)==0;
        const bool needReference=mode_==PureFastRunMode::AuditLockstep||
            output.ratio.trust_alarm||periodicAudit;
        bool referenceControlsDecision=mode_==PureFastRunMode::AuditLockstep||
            output.ratio.trust_alarm;
        PureProjectorWeightResult candidateWeight;
        bool mpControlsDecision=false;
        const bool endpointObserverOnly=periodicAudit&&
            mode_==PureFastRunMode::FastStrict&&!output.ratio.trust_alarm;
        if(endpointObserverOnly){
            ++diagnostics_.ratio_slow_reference_count;output.ratio.used_reference=true;
            PureEndpointRebuildResult oldEndpoint=pureProjectorEndpointRebuild(
                trial_,configuration_.slices,int(configuration_.slices.size()/2),
                block_size_,stack_options_);
            PureEndpointRebuildResult newEndpoint=pureProjectorEndpointRebuild(
                trial_,candidateConfiguration().slices,int(candidateConfiguration().slices.size()/2),
                block_size_,stack_options_);
            if(oldEndpoint.ok()&&newEndpoint.ok()&&output.ratio.ok()){
                output.ratio.slow_ratio=std::exp(
                    newEndpoint.log_abs_weight-oldEndpoint.log_abs_weight)*
                    pureProjectorUnitPhase(output.ratio.ratio);
                output.ratio.relative_reference_error=
                    std::abs(output.ratio.ratio-output.ratio.slow_ratio)/
                    std::max(options_.zero_tolerance,std::abs(output.ratio.slow_ratio));
                diagnostics_.maximum_ratio_reference_error=std::max(
                    diagnostics_.maximum_ratio_reference_error,
                    output.ratio.relative_reference_error);
            }else ++diagnostics_.read_only_endpoint_audit_failure_count;
        }
        if(needReference&&!endpointObserverOnly){
            ++diagnostics_.ratio_slow_reference_count;
            output.ratio.used_reference=true;
            const auto slowStart=std::chrono::steady_clock::now();
            if(output.ratio.trust_alarm&&
               options_.weight_mode==PureProjectorWeightMode::RealZ2){
                if(!options_.mp_same_proposal_fallback){
                    output.ratio.status=PureFastRatioStatus::reference_failure;
                    ++diagnostics_.mp_fallback_failure_count;
                    if(diagnostics_.first_failure_proposal<0)
                        diagnostics_.first_failure_proposal=int(proposalNumber);
                    terminated_=true;output.terminated=true;
                    output.failure_stage="mp_fallback_disabled";return output;
                }
                ++diagnostics_.mp_fallback_count;
                PureMpProposalOptions mpOptions;
                mpOptions.minimum_endpoint_rcond=options_.minimum_overlap_rcond;
                mpOptions.residual_tolerance=options_.residual_tolerance;
                mpOptions.reality_tolerance=options_.reality_tolerance;
                mpOptions.zero_tolerance=options_.zero_tolerance;
                mpOptions.endpoint_block_size=block_size_;
                const auto mpStarted=std::chrono::steady_clock::now();
                output.ratio.mp_reference=pureProjectorMpSameProposal(
                    trial_,configuration_.slices,proposal.index,
                    proposal.new_factor,proposal.new_eta,mpOptions,
                    &mp_subspace_cache_,output.snapshot.configuration_hash);
                diagnostics_.total_mp_fallback_seconds+=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-mpStarted).count();
                diagnostics_.mp_profile+=output.ratio.mp_reference.profile;
                if(!output.ratio.mp_reference.ok()||
                   !std::isfinite(output.ratio.mp_reference.ratio.real())||
                   !std::isfinite(output.ratio.mp_reference.ratio.imag())){
                    ++diagnostics_.mp_endpoint_fallback_count;
                    const auto endpointMpStarted=std::chrono::steady_clock::now();
                    const PureMpProposalResult endpointReference=
                        pureProjectorMpEndpointSameProposal(
                            trial_,configuration_.slices,proposal.index,
                            proposal.new_factor,proposal.new_eta,mpOptions);
                    diagnostics_.total_mp_fallback_seconds+=std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-endpointMpStarted).count();
                    diagnostics_.mp_profile+=endpointReference.profile;
                    output.ratio.mp_reference=endpointReference;
                    const bool safeEndpointRejection=
                        !endpointReference.ok()&&endpointReference.ratio_converged&&
                        std::isfinite(endpointReference.ratio.real())&&
                        std::isfinite(endpointReference.ratio.imag())&&
                        proposal.uniform>=std::min(
                            1.0,std::abs(endpointReference.ratio));
                    if((!endpointReference.ok()&&!safeEndpointRejection)||
                       !std::isfinite(endpointReference.ratio.real())||
                       !std::isfinite(endpointReference.ratio.imag())){
                        output.ratio.status=PureFastRatioStatus::reference_failure;
                        ++diagnostics_.mp_fallback_failure_count;
                        ++diagnostics_.mp_endpoint_fallback_failure_count;
                        if(diagnostics_.first_failure_proposal<0)
                            diagnostics_.first_failure_proposal=int(proposalNumber);
                        terminated_=true;output.terminated=true;
                        output.failure_stage="mp_endpoint_same_proposal_untrusted";
                        return output;
                    }
                    if(safeEndpointRejection)
                        ++diagnostics_.mp_endpoint_rejection_recovery_count;
                }
                output.ratio.slow_ratio=output.ratio.mp_reference.ratio;
                output.ratio.ratio=output.ratio.slow_ratio;
                output.ratio.status=PureFastRatioStatus::success;
                output.reference_accepted=proposal.uniform<
                    std::min(1.0,std::abs(output.ratio.slow_ratio));
                mpControlsDecision=true;
            }else{
                candidateWeight=pureProjectorStableReferenceWeight(
                    trial_,candidateConfiguration().slices,options_);
            }
            output.slow_seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-slowStart).count();
            diagnostics_.total_reference_seconds+=output.slow_seconds;
            if(!mpControlsDecision&&(!candidateWeight.ok()||
               (options_.weight_mode==PureProjectorWeightMode::RealZ2&&candidateWeight.z2_sign==0)||
               !std::isfinite(current_.log_abs_weight)||
               std::abs(current_.complex_phase)<=options_.zero_tolerance)){
                output.ratio.status=PureFastRatioStatus::reference_failure;
                ++diagnostics_.slow_reference_failure_count;
                if(diagnostics_.first_failure_proposal<0)
                    diagnostics_.first_failure_proposal=int(proposalNumber);
                terminated_=true;output.terminated=true;
                output.failure_stage="stable_reference_weight_untrusted";return output;}
            if(!mpControlsDecision){
                output.ratio.slow_ratio=std::exp(candidateWeight.log_abs_weight-current_.log_abs_weight)*
                    candidateWeight.complex_phase/current_.complex_phase;
                if(output.ratio.ok())output.ratio.relative_reference_error=
                    std::abs(output.ratio.ratio-output.ratio.slow_ratio)/
                    std::max(1e-14,std::abs(output.ratio.slow_ratio));
                if(std::isfinite(output.ratio.relative_reference_error))
                    diagnostics_.maximum_ratio_reference_error=std::max(
                        diagnostics_.maximum_ratio_reference_error,
                        output.ratio.relative_reference_error);
                const double slowProbability=std::min(1.0,std::abs(output.ratio.slow_ratio));
                output.reference_accepted=proposal.uniform<slowProbability;
                if(output.ratio.fast_updated_green.size()!=0){
                    PureProjectorStackManager candidateStack(
                        trial_,candidateConfiguration().slices,block_size_,stack_options_,proposal.index+1);
                    if(candidateStack.ok()){
                        const PureProjectorGreenResult referenceGreen=candidateStack.green();
                        if(referenceGreen.ok())output.ratio.green_reference_error=
                            (output.ratio.fast_updated_green-referenceGreen.green).norm()/
                            std::max(1.0,referenceGreen.green.norm());
                    }
                }
                if(std::isfinite(output.ratio.green_reference_error))
                    diagnostics_.maximum_green_rebuild_error=std::max(
                        diagnostics_.maximum_green_rebuild_error,
                        output.ratio.green_reference_error);
                if(std::isfinite(output.ratio.green_reference_error)&&
                   output.ratio.green_reference_error>options_.green_update_tolerance&&
                   !output.ratio.trust_alarm){
                    output.ratio.trust_alarm=true;++diagnostics_.trust_alarm_count;
                    referenceControlsDecision=true;
                }
                if(output.ratio.trust_alarm){
                    output.ratio.ratio=output.ratio.slow_ratio;
                    output.ratio.status=PureFastRatioStatus::success;
                }
            }
        }
        if(output.ratio.ok()&&!checkPositiveRatio(proposal,output.ratio)){
            terminated_=true;output.terminated=true;
            output.failure_stage="positive_ratio_reference_unresolved";return output;
        }
        if(!output.ratio.ok()){
            moveToTrustedCut(originalCut);
            synchronizeCurrentGreen();
            return output;
        }
        const double probability=std::min(1.0,std::abs(output.ratio.ratio));
        output.accepted=proposal.uniform<probability;
        if(mode_==PureFastRunMode::AuditLockstep&&output.accepted!=output.reference_accepted){
            output.ratio.status=PureFastRatioStatus::reference_failure;
            terminated_=true;output.terminated=true;
            output.failure_stage="audit_lockstep_decision_mismatch";return output;
        }
        if(output.accepted){
            bool acceptedStackRebuilt=false;
            const int ratioSign=output.ratio.ratio.real()>=0?1:-1;
            if(options_.weight_mode==PureProjectorWeightMode::RealZ2)z2_sign_*=ratioSign;
            if(!mpControlsDecision){PureProjectorSlice replacement(
                    proposal.new_factor,proposal.new_eta,
                    configuration_.slices[proposal.index].label);
                if(!manager_->acceptFactorAtCut(proposal.index,replacement,rightBeforeForUpdate)){
                    const PureStackStatus failedStatus=manager_->status();
                    diagnostics_.last_accepted_factor_failure_status=failedStatus;
                    // The ratio and accept/reject decision are already trusted.
                    // Failure here can be local to the live-cut gauge after an
                    // accepted update, while the same physical configuration
                    // remains well conditioned at the theorem-backed center.
                    // Rebuild that accepted candidate from scratch without
                    // consuming RNG or changing the Markov decision.
                    const int center=int(candidateConfiguration().slices.size()/2);
                    manager_.reset(new PureProjectorStackManager(
                        trial_,candidateConfiguration().slices,block_size_,stack_options_,center));
                    if(!manager_->ok()){
                        ++diagnostics_.accepted_factor_center_recovery_failure_count;
                        terminated_=true;output.terminated=true;
                        output.failure_stage=std::string("accepted_factor_center_recovery_failed_after_")+
                            pureStackStatusName(failedStatus);return output;
                    }
                    ++diagnostics_.accepted_factor_center_recovery_count;
                    ++diagnostics_.rebuild_count;
                    accepted_since_rebuild_=0;
                    acceptedStackRebuilt=true;
                }}
            if(candidate)configuration_=std::move(*candidate);
            else {
                configuration_.slices[proposal.index].matrix=proposal.new_factor;
                configuration_.slices[proposal.index].eta=proposal.new_eta;
                configuration_.hs_fields[proposal.index]=proposal.new_hs;
            }
            configuration_hash_valid_=false;
            diagnostics_.mp_profile+=mp_subspace_cache_.acceptedUpdate(
                output.snapshot.configuration_hash,configurationHash(),
                proposal.index,configuration_.slices);
            if(referenceControlsDecision&&!mpControlsDecision)current_=candidateWeight;
            else{
                current_.log_abs_weight+=std::log(std::abs(output.ratio.ratio));
                current_.complex_phase*=pureProjectorUnitPhase(output.ratio.ratio);
                current_.weight=std::exp(current_.log_abs_weight)*current_.complex_phase;
                current_.z2_sign=z2_sign_;
            }
            if(!mpControlsDecision&&!acceptedStackRebuilt&&
               ++accepted_since_rebuild_>=std::uint64_t(std::max(1,block_size_))){
                if(!manager_->rebuildCurrent()){
                    const int center=int(configuration_.slices.size()/2);
                    manager_.reset(new PureProjectorStackManager(
                        trial_,configuration_.slices,block_size_,stack_options_,center));
                    if(!manager_->ok()){
                        ++diagnostics_.post_decision_center_recovery_failure_count;
                        terminated_=true;output.terminated=true;
                        output.failure_stage="periodic_rebuild_and_center_recovery_failed";
                        return output;
                    }
                    ++diagnostics_.post_decision_center_recovery_count;
                }
                accepted_since_rebuild_=0;++diagnostics_.rebuild_count;
            }
        }else if(!mpControlsDecision&&!moveToTrustedCut(originalCut)){
            const int center=int(configuration_.slices.size()/2);
            manager_.reset(new PureProjectorStackManager(
                trial_,configuration_.slices,block_size_,stack_options_,center));
            if(!manager_->ok()){
                ++diagnostics_.post_decision_center_recovery_failure_count;
                terminated_=true;output.terminated=true;
                output.failure_stage="rejected_cut_restore_and_center_recovery_failed";
            }else{
                ++diagnostics_.post_decision_center_recovery_count;
                ++diagnostics_.rebuild_count;
                accepted_since_rebuild_=0;
            }
        }
        if(mpControlsDecision){
            manager_.reset(new PureProjectorStackManager(
                trial_,configuration_.slices,block_size_,stack_options_,originalCut));
            if(!manager_->ok()){
                const MatType &mpGreen=output.accepted?
                    output.ratio.mp_reference.trusted_post_green:
                    output.ratio.mp_reference.trusted_pre_green;
                const double residual=mpGreen.size()? (mpGreen*mpGreen-
                    MatType::Identity(mpGreen.rows(),mpGreen.cols())).norm()/
                    std::max(1.0,mpGreen.norm()):std::numeric_limits<double>::infinity();
                if(!pureProjectorMatrixFinite(mpGreen)||!std::isfinite(residual)||
                   residual>options_.green_update_tolerance){
                    output.ratio.status=PureFastRatioStatus::reference_failure;
                    ++diagnostics_.mp_fallback_failure_count;
                    terminated_=true;output.terminated=true;
                    output.failure_stage="mp_green_double_recovery_input_untrusted";return output;
                }
                // Re-anchor the same accepted/rejected configuration at the
                // theorem-backed center cut.  This discards only the
                // untrusted live-cut gauge; it does not change HS, Z2, log
                // weight, or RNG.
                const int center=int(configuration_.slices.size()/2);
                manager_.reset(new PureProjectorStackManager(
                    trial_,configuration_.slices,block_size_,stack_options_,center));
                if(!manager_->ok()){
                    output.ratio.status=PureFastRatioStatus::reference_failure;
                    ++diagnostics_.mp_fallback_failure_count;
                    terminated_=true;output.terminated=true;
                    output.failure_stage="mp_center_reanchor_failed";return output;
                }
                ++diagnostics_.mp_green_double_recovery_count;
            }
            const PureProjectorGreenResult restored=manager_->green();
            if(!restored.ok()){
                output.ratio.status=PureFastRatioStatus::reference_failure;
                ++diagnostics_.mp_fallback_failure_count;
                terminated_=true;output.terminated=true;
                output.failure_stage="mp_restored_green_untrusted";return output;
            }
            diagnostics_.maximum_endpoint_rebuild_green_residual=std::max(
                diagnostics_.maximum_endpoint_rebuild_green_residual,
                restored.green_residual);
            ++diagnostics_.rebuild_count;accepted_since_rebuild_=0;
        }
        if(needReference&&candidateWeight.ok()){
            const int expected=output.accepted?
                (candidateWeight.z2_sign==0?z2_sign_:candidateWeight.z2_sign):
                (current_.z2_sign==0?z2_sign_:current_.z2_sign);
            output.z2_reference_mismatch=(expected!=z2_sign_);
        }
        // currentWeight().green is a public diagnostic.  A low-rank accepted
        // update changes the live cut even when the absolute-weight fields are
        // transported incrementally, so never leave the previous Green there.
        // If the live Green cannot be obtained, explicitly invalidate that
        // diagnostic rather than exposing stale state.
        synchronizeCurrentGreen();
        return output;
    }

private:
    bool checkPositiveRatio(const PureFastProposal &proposal,PureFastRatioResult &ratio) const {
        if(!options_.positive_hs1_ratio_contract||!ratio.ok())return ratio.ok();
        const double bound=proposal.positive_ratio_log_bound;
        const double re=ratio.ratio.real(),im=ratio.ratio.imag();
        const double tolerance=options_.reality_tolerance;
        const double phaseTolerance=options_.stabilized_local_updates?
            options_.local_reality_tolerance:options_.reality_tolerance;
        if(!std::isfinite(bound)||bound<0||!std::isfinite(re)||!std::isfinite(im)||re<=0||
           std::abs(im)>phaseTolerance*std::max(1.,std::abs(re))||
           std::abs(std::log(re))>bound+tolerance){
            ratio.status=PureFastRatioStatus::positive_ratio_bound_violation;
            return false;
        }
        return true;
    }

    PureFastProposalResult proposeLocal(const PureFastProposal &proposal) {
        PureFastProposalResult output;
        if(terminated_){output.terminated=true;output.ratio.status=PureFastRatioStatus::terminated;
            output.failure_stage="walker_already_terminated";return output;}
        const bool immutable=proposal.immutable_factor_only&&proposal.new_factor.size()==0&&bool(proposal.cached_factor);
        const bool factorValid=immutable?
            (proposal.cached_factor->rows()==trial_.Phi.rows()&&proposal.cached_factor->cols()==trial_.Phi.rows()&&proposal.cached_factor->allFinite()):
            (proposal.new_factor.rows()==trial_.Phi.rows()&&proposal.new_factor.cols()==trial_.Phi.rows()&&pureProjectorMatrixFinite(proposal.new_factor));
        if(proposal.index<0||proposal.index>=int(configuration_.slices.size())||
           (proposal.new_hs!=-1&&proposal.new_hs!=1)||!std::isfinite(proposal.uniform)||
           proposal.uniform<0||proposal.uniform>=1||
           !factorValid)return output;
        ++diagnostics_.proposals;
        const auto start=std::chrono::steady_clock::now();
        const auto &location=configuration_.locations[proposal.index];
        output.snapshot.branch=location.branch;output.snapshot.slice=location.slice;
        output.snapshot.factor=location.factor;output.snapshot.bond=location.bond;
        output.snapshot.aux=location.aux;
        output.snapshot.old_hs=configuration_.hs_fields[proposal.index];
        output.snapshot.new_hs=proposal.new_hs;output.snapshot.uniform=proposal.uniform;
        auto fail=[&](const char *stage){
            terminated_=true;output.terminated=true;output.failure_stage=stage;
            output.snapshot.configuration_hash=configurationHash();
            output.fast_seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
            return output;
        };
        auto move=[&](int target){
            // Long jumps use the same forward/adjoint QR path, starting from
            // exact valid checkpoints when available. Failed motion or cached
            // reconstruction still takes the original fresh double retry.
            const bool nearby=std::abs(target-manager_->cut())<=block_size_;
            if(nearby&&manager_->moveToCut(target))return true;
            ++diagnostics_.rebuild_count;
            if(!nearby&&options_.reuse_long_jump_checkpoints&&manager_->rebuildToCut(target))return true;
            return manager_->rebuildToCutFresh(target);
        };
        MatType temporaryRightBefore;
        MatType &rightBefore=options_.reuse_green_buffers?right_before_scratch_:temporaryRightBefore;
        bool prepared=false;
        if(options_.avoid_reverse_round_trip&&manager_->cut()>proposal.index){
            // The ratio already uses G at index+1 in either sweep direction.
            // Only R0 is needed from index: solve B_old R0 = R at the live
            // cut instead of moving both frames and G backward then forward.
            prepared=move(proposal.index+1)&&
                manager_->rightBeforeFactorAtCut(proposal.index,rightBefore);
            if(prepared)++diagnostics_.reverse_cut_shortcuts;
        }
        if(!prepared){
            // Preserve the original checked path on forward scans and if
            // the read-only solve fails; move() retains its fresh retry.
            if(move(proposal.index)){
                rightBefore=manager_->phiRight();
                prepared=move(proposal.index+1);
            }
        }
        if(!prepared){
            output.ratio.status=PureFastRatioStatus::overlap_untrusted;
            return fail("local_cut_rebuild_failed");
        }
        PureFactorMatrix cachedNew;
        if(options_.cached_factor_support){
            if(immutable)cachedNew=*proposal.cached_factor;
            else if(proposal.cached_factor&&proposal.cached_factor->bitwiseMatches(proposal.new_factor))
                cachedNew=*proposal.cached_factor;
            else cachedNew=proposal.new_factor;
        }
        // Compact algebra needs only immutable coefficients. General/disabled
        // paths still receive the exact dense factor and original arithmetic.
        MatType denseFallback;
        std::vector<int> unionSupport;
        const auto &oldStored=configuration_.slices[proposal.index].matrix;
        const bool supported=immutable&&cachedNew.hasLocalSupport()&&oldStored.hasLocalSupport();
        if(supported)std::set_union(cachedNew.localSupport().begin(),cachedNew.localSupport().end(),
            oldStored.localSupport().begin(),oldStored.localSupport().end(),std::back_inserter(unionSupport));
        const bool needDense=immutable&&(!options_.cached_factor_support||!options_.compact_proposal_algebra||
            !supported||unionSupport.size()>16||unionSupport.size()>=std::size_t(oldStored.rows()/2));
        if(needDense)denseFallback=*proposal.cached_factor;
        const MatType &newFactor=needDense?denseFallback:proposal.new_factor;
        auto compute=[&](){
            PureProjectorGreenResult snapshot;
            const PureProjectorGreenResult &trusted=options_.borrow_green?
                manager_->greenView():(snapshot=manager_->green());
            auto ratio=pureFastLocalRatioFromGreen(trial_,configuration_.slices,
                proposal.index,newFactor,proposal.new_eta,rightBefore,
                trusted,options_,proposal.uniform,
                options_.cached_factor_support?&cachedNew:nullptr);
            output.snapshot.fast_ratio=ratio.ratio;
            checkPositiveRatio(proposal,ratio);
            return ratio;
        };
        output.ratio=compute();
        if(!output.ratio.ok()){
            ++diagnostics_.local_double_retries;
            // At most one double re-anchor. No changed random draw or MP run.
            if(!manager_->rebuildToCutFresh(proposal.index))return fail("local_ratio_rebuild_failed");
            rightBefore=manager_->phiRight();
            if(!manager_->moveToCut(proposal.index+1))return fail("local_ratio_rebuild_failed");
            output.ratio=compute();
            if(!output.ratio.ok())return fail("local_ratio_unresolved");
        }
        diagnostics_.maximum_local_ratio_imaginary=std::max(
            diagnostics_.maximum_local_ratio_imaginary,
            std::abs(output.ratio.ratio.imag())/std::max(1.,std::abs(output.ratio.ratio.real())));
        output.snapshot.fast_ratio=output.ratio.ratio;
        if(options_.weight_mode==PureProjectorWeightMode::RealZ2)
            output.ratio.ratio=DataType(output.ratio.ratio.real(),0);
        output.accepted=proposal.uniform<std::min(1.,std::abs(output.ratio.ratio));
        if(output.accepted){
            PureProjectorSlice replacement;
            replacement.matrix=options_.cached_factor_support?cachedNew:PureFactorMatrix(newFactor);
            replacement.eta=proposal.new_eta;
            replacement.label=configuration_.slices[proposal.index].label;
            if(!manager_->acceptFactorAtCut(proposal.index,replacement,rightBefore,
                                           &output.ratio.fast_updated_green))
                return fail("local_accepted_update_failed");
            configuration_.slices[proposal.index]=replacement;
            configuration_.hs_fields[proposal.index]=proposal.new_hs;
            configuration_hash_valid_=false;
            if(options_.weight_mode==PureProjectorWeightMode::RealZ2)
                z2_sign_*=output.ratio.ratio.real()>=0?1:-1;
            current_.log_abs_weight+=std::log(std::abs(output.ratio.ratio));
            current_.complex_phase*=pureProjectorUnitPhase(output.ratio.ratio);
            current_.weight=current_.log_abs_weight<700.?
                std::exp(current_.log_abs_weight)*current_.complex_phase:current_.complex_phase;
            current_.z2_sign=z2_sign_;
        }
        synchronizeCurrentGreen();
        output.fast_seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
        diagnostics_.total_fast_seconds+=output.fast_seconds;
        return output;
    }

    void synchronizeCurrentGreen() {
        if(manager_&&manager_->ok()){
            PureProjectorGreenResult snapshot;
            const PureProjectorGreenResult &live=options_.stabilized_local_updates&&options_.borrow_green?
                manager_->greenView():(snapshot=manager_->green());
            if(live.ok()){
                current_.green=live.green;
                current_.overlap_rank=live.overlap_rank;
                current_.overlap_rcond=live.overlap_rcond;
                current_.overlap_residual=live.solve_residual;
                current_.green_residual=live.green_residual;
                return;
            }
        }
        current_.green.resize(0,0);
        current_.overlap_rank=0;
        current_.overlap_rcond=0.0;
        current_.overlap_residual=std::numeric_limits<double>::infinity();
        current_.green_residual=std::numeric_limits<double>::infinity();
    }

    GaussianTrialState trial_;
    PureFastConfiguration configuration_;
    int block_size_;
    PureFastRunMode mode_;
    PureFastOptions options_;
    PureProjectorWeightResult current_;
    std::unique_ptr<PureProjectorStackManager> manager_;
    MatType right_before_scratch_;
    int z2_sign_=1;
    bool terminated_=false;
    PureFastDiagnostics diagnostics_;
    PureProjectorOptions stack_options_;
    PureFastInitializationPolicy initialization_policy_=
        PureFastInitializationPolicy::SequentialAudit;
    std::uint64_t accepted_since_rebuild_=0;
    PureMpSubspaceCache mp_subspace_cache_;
    mutable bool configuration_hash_valid_=false;
    std::uint64_t measurement_scaled_green_rebuilds_=0;
    mutable std::uint64_t configuration_hash_=0;
};

#endif
