#ifndef PURE_PROJECTOR_GREEN_H
#define PURE_PROJECTOR_GREEN_H

#include "operator.h"
#include "types.h"
#include "pure_projector_svd.h"
#include "pure_projector_factor.h"
#include "pure_projector_canonical.h"
#include <vector>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

enum class PureProjectorStatus {
    success,
    invalid_dimension,
    nonfinite_input,
    subspace_rank_deficient,
    overlap_rank_deficient,
    overlap_ill_conditioned,
    solve_residual_exceeded,
    green_residual_exceeded,
    green_structure_exceeded
};

inline const char *pureProjectorStatusName(PureProjectorStatus status) {
    switch (status) {
    case PureProjectorStatus::success: return "success";
    case PureProjectorStatus::invalid_dimension: return "invalid_dimension";
    case PureProjectorStatus::nonfinite_input: return "nonfinite_input";
    case PureProjectorStatus::subspace_rank_deficient: return "subspace_rank_deficient";
    case PureProjectorStatus::overlap_rank_deficient: return "overlap_rank_deficient";
    case PureProjectorStatus::overlap_ill_conditioned: return "overlap_ill_conditioned";
    case PureProjectorStatus::solve_residual_exceeded: return "solve_residual_exceeded";
    case PureProjectorStatus::green_residual_exceeded: return "green_residual_exceeded";
    case PureProjectorStatus::green_structure_exceeded: return "green_structure_exceeded";
    }
    return "unknown";
}

struct PureProjectorOptions {
    double rank_tolerance = 1e-12;
    double minimum_overlap_rcond = 1e-12;
    double solve_residual_tolerance = 1e-10;
    // Optional production tolerance scaled by the measured overlap rcond.
    // Zero preserves the strict solve gate. The original gate remains the
    // floor for ill-conditioned systems, which still receive full recovery.
    double well_conditioned_solve_tolerance = 0.0;
    double green_residual_tolerance = 1e-10;
    // Stable signed weights already accept the relative involution residual
    // of an oblique Green. Apply that same policy before expensive recovery;
    // strict/reference callers retain the dimension-normalized default.
    bool allow_oblique_green = false;
    // Production solves with pivoted LU after the same SVD rank/rcond gate.
    // The SVD then needs singular values only; strict callers keep its solve.
    bool prefer_lu_green_solve = false;
    // Always compute skew/diagonal diagnostics.  Enabling this flag promotes
    // them to an additional fail-closed gate; it is off by default so this
    // review does not silently change the frozen Phase 3C trust policy.
    bool fail_on_green_structure = false;
    // Off by default: endpoint and independent reference arithmetic is frozen.
    bool accelerated_linear_algebra = false;
    // Used only by cached factor views in the local stack.
    bool pair_block_factors = true;
    bool direct_green_columns = true;
    // Opt-in until HPC confirms the parameter-dependent timing benefit.
    bool reuse_rebuild_workspace = false;
};

struct ThinQrResult {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType q;
    int rank = 0;
    double orthonormal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

struct PureProjectorGreenResult {
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    MatType green;
    int overlap_rank = 0;
    double overlap_rcond = 0.0;
    double solve_residual = std::numeric_limits<double>::infinity();
    double initial_solve_residual = std::numeric_limits<double>::infinity();
    int solve_refinement_steps = 0;
    bool extended_precision_solve = false;
    double green_residual = std::numeric_limits<double>::infinity();
    double initial_green_residual = std::numeric_limits<double>::infinity();
    double green_residual_double = std::numeric_limits<double>::infinity();
    bool long_double_green_residual = false;
    bool canonical_frame_maintenance_attempted = false;
    bool canonical_frame_maintenance = false;
    double canonical_correction_right = 0.0,canonical_correction_left = 0.0;
    double canonical_isotropy_right = 0.0,canonical_isotropy_left = 0.0;
    int canonical_overlap_rank = 0;
    double canonical_overlap_rcond = 0.0;
    // Set only by the local oblique-Green recovery; the old diagnostic keeps
    // its dimension-normalized meaning for strict/reference callers.
    double green_relative_residual = std::numeric_limits<double>::infinity();
    double green_skew_residual = std::numeric_limits<double>::infinity();
    double green_diagonal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

struct PureProjectorGreenPrecisionStats {
    std::uint64_t long_double_residual_checks=0,long_double_residual_passes=0;
    std::uint64_t canonical_ld_rebuild_attempts=0,canonical_ld_rebuild_assemblies=0;
    std::uint64_t canonical_ld_rebuild_certified=0;
    double maximum_canonical_relative_correction=0.0;
};
inline thread_local PureProjectorGreenPrecisionStats pureProjectorGreenPrecisionStats;

inline bool pureProjectorMatrixFinite(const PureFactorMatrix &matrix) { return matrix.allFinite(); }

inline bool pureProjectorMatrixFinite(const MatType &matrix) {
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = 0; row < matrix.rows(); ++row) {
            if (!std::isfinite(matrix(row, col).real()) ||
                !std::isfinite(matrix(row, col).imag())) return false;
        }
    }
    return true;
}

inline double pureProjectorSolveTolerance(
    const PureProjectorOptions &options, double overlapRcond) {
    if (!std::isfinite(overlapRcond) || overlapRcond <= 0.0 ||
        !std::isfinite(options.well_conditioned_solve_tolerance))
        return options.solve_residual_tolerance;
    return std::max(options.solve_residual_tolerance,
        options.well_conditioned_solve_tolerance * std::min(1.0, overlapRcond));
}

inline ThinQrResult thinQrSubspace(
    const MatType &subspace,
    const PureProjectorOptions &options = PureProjectorOptions(),
    bool canonical = false) {
    ThinQrResult result;
    if (subspace.rows() <= 0 || subspace.cols() <= 0 ||
        subspace.rows() < subspace.cols()) return result;
    if (!pureProjectorMatrixFinite(subspace)) {
        result.status = PureProjectorStatus::nonfinite_input;
        return result;
    }
    if(options.accelerated_linear_algebra&&subspace.cols()>=32){
        MatType packed=subspace;
        const lapack_int m=subspace.rows(),n=subspace.cols();
        std::vector<lapack_int> pivots(n,0);
        cVecType tau(n);
        const lapack_int info=LAPACKE_zgeqp3(LAPACK_COL_MAJOR,m,n,packed.data(),m,pivots.data(),tau.data());
        if(info!=0||!pureProjectorMatrixFinite(packed)){
            result.status=PureProjectorStatus::subspace_rank_deficient;return result;
        }
        const double maximum=packed.diagonal().cwiseAbs().maxCoeff();
        result.rank=(packed.diagonal().cwiseAbs().array()>maximum*options.rank_tolerance).count();
        if(result.rank!=n){result.status=PureProjectorStatus::subspace_rank_deficient;return result;}
        if(options.reuse_rebuild_workspace)result.q.swap(packed);
        else result.q=packed;
        if(LAPACKE_zungqr(LAPACK_COL_MAJOR,m,n,n,result.q.data(),m,tau.data())!=0){
            result.status=PureProjectorStatus::subspace_rank_deficient;return result;
        }
    }else{
    Eigen::ColPivHouseholderQR<MatType> qr(subspace);
    qr.setThreshold(options.rank_tolerance);
    result.rank = qr.rank();
    if (result.rank != subspace.cols()) {
        result.status = PureProjectorStatus::subspace_rank_deficient;
        return result;
    }
    result.q = qr.householderQ() *
        MatType::Identity(subspace.rows(), subspace.cols());
    }
    if(canonical) {
        DataType correctionDet;
        if(!pureProjectorCanonicalMicrocorrect(result.q,correctionDet,
                options.solve_residual_tolerance,options.accelerated_linear_algebra)) {
            result.status=PureProjectorStatus::solve_residual_exceeded;
            return result;
        }
    }
    result.orthonormal_residual =
        (result.q.adjoint() * result.q -
         MatType::Identity(subspace.cols(), subspace.cols())).norm() /
        std::max(1.0, std::sqrt(double(subspace.cols())));
    if (!pureProjectorMatrixFinite(result.q) ||
        result.orthonormal_residual > options.solve_residual_tolerance) {
        result.status = PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    result.status = PureProjectorStatus::success;
    return result;
}

inline PureProjectorGreenResult pureProjectorGreen(
    const MatType &phiRight, const MatType &phiLeft,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    PureProjectorGreenResult result;
    if (phiRight.rows() <= 0 || phiRight.cols() <= 0 ||
        phiRight.rows() != phiLeft.rows() || phiRight.cols() != phiLeft.cols() ||
        phiRight.rows() < phiRight.cols()) return result;
    if (!pureProjectorMatrixFinite(phiRight) || !pureProjectorMatrixFinite(phiLeft)) {
        result.status = PureProjectorStatus::nonfinite_input;
        return result;
    }

    const MatType overlap = phiLeft.adjoint() * phiRight;
    PureGreenSvd svd(overlap,options.accelerated_linear_algebra,!options.prefer_lu_green_solve);
    svd.setThreshold(options.rank_tolerance);
    result.overlap_rank = svd.rank();
    if (result.overlap_rank != overlap.rows()) {
        result.status = PureProjectorStatus::overlap_rank_deficient;
        return result;
    }
    const double largest = svd.singularValues()(0);
    const double smallest = svd.singularValues()(svd.singularValues().size() - 1);
    result.overlap_rcond = largest > 0.0 ? smallest / largest : 0.0;
    if (!std::isfinite(result.overlap_rcond) ||
        result.overlap_rcond < options.minimum_overlap_rcond) {
        result.status = PureProjectorStatus::overlap_ill_conditioned;
        return result;
    }

    const double solveTolerance = pureProjectorSolveTolerance(options, result.overlap_rcond);

    const MatType rightHandSide = phiLeft.adjoint();
    const double rightHandSideScale = std::max(1.0, rightHandSide.norm());
    Eigen::PartialPivLU<MatType> lu;
    MatType solved;
    if(options.prefer_lu_green_solve) {
        lu.compute(overlap);solved=lu.solve(rightHandSide);
    } else solved=svd.solve(rightHandSide);
    result.solve_residual = (overlap * solved - rightHandSide).norm() /
        rightHandSideScale;
    result.initial_solve_residual = result.solve_residual;
    using ExtendedComplex=std::complex<long double>;
    using ExtendedMatrix=Eigen::Matrix<ExtendedComplex,Eigen::Dynamic,Eigen::Dynamic>;
    const double independentFrameTolerance=std::min(options.solve_residual_tolerance,
        pureProjectorCanonicalCorrectionBound(phiRight.rows()));
    const bool canonicalFrames=
        pureProjectorCanonicalFrameTrusted(phiRight,independentFrameTolerance)&&
        pureProjectorCanonicalFrameTrusted(phiLeft,independentFrameTolerance);
    bool canonicalMaintenanceFailed=false;
    PureProjectorStatus extendedFailureStatus=PureProjectorStatus::solve_residual_exceeded;
    auto extendedSolve=[&]() {
        if(std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits)return false;
        ExtendedMatrix right=phiRight.cast<ExtendedComplex>();
        ExtendedMatrix left=phiLeft.cast<ExtendedComplex>();
        if(canonicalFrames) {
            result.canonical_frame_maintenance_attempted=true;
            ++pureProjectorGreenPrecisionStats.canonical_ld_rebuild_attempts;
            if(!pureProjectorCanonicalMicrocorrectExtended(right,result.canonical_correction_right,
                    result.canonical_isotropy_right,options.solve_residual_tolerance)||
               !pureProjectorCanonicalMicrocorrectExtended(left,result.canonical_correction_left,
                    result.canonical_isotropy_left,options.solve_residual_tolerance)) {
                canonicalMaintenanceFailed=true;return false;
            }
            pureProjectorGreenPrecisionStats.maximum_canonical_relative_correction=std::max(
                pureProjectorGreenPrecisionStats.maximum_canonical_relative_correction,
                std::max(result.canonical_correction_right,result.canonical_correction_left));
        }
        const ExtendedMatrix leftAdj=left.adjoint();
        // Form the overlap at the precision of the solve. Generic subspaces
        // remain unchanged; independently certified canonical frames receive
        // only the bounded roundoff maintenance recorded above.
        const ExtendedMatrix matrix=leftAdj*right;
        if(canonicalFrames) {
            const MatType correctedOverlap=matrix.cast<DataType>();
            PureGreenSvd correctedSvd(correctedOverlap,options.accelerated_linear_algebra,false);
            correctedSvd.setThreshold(options.rank_tolerance);
            result.canonical_overlap_rank=correctedSvd.rank();
            if(result.canonical_overlap_rank!=matrix.rows()) {
                canonicalMaintenanceFailed=true;
                extendedFailureStatus=PureProjectorStatus::overlap_rank_deficient;return false;
            }
            const auto &values=correctedSvd.singularValues();
            result.canonical_overlap_rcond=values(0)>0.?values(values.size()-1)/values(0):0.;
            if(!std::isfinite(result.canonical_overlap_rcond)||
               result.canonical_overlap_rcond<options.minimum_overlap_rcond) {
                canonicalMaintenanceFailed=true;
                extendedFailureStatus=PureProjectorStatus::overlap_ill_conditioned;return false;
            }
            result.overlap_rcond=std::min(result.overlap_rcond,result.canonical_overlap_rcond);
        }
        const double extendedSolveTolerance=std::min(solveTolerance,
            pureProjectorSolveTolerance(options,result.overlap_rcond));
        Eigen::PartialPivLU<ExtendedMatrix> lu(matrix);
        const ExtendedMatrix solution=lu.solve(leftAdj);
        const long double residual=(matrix*solution-leftAdj).norm()/std::max(1.L,leftAdj.norm());
        if(!solution.allFinite()||!std::isfinite(residual)||residual>extendedSolveTolerance)return false;
        const ExtendedMatrix green=ExtendedMatrix::Identity(phiRight.rows(),phiRight.rows())-2.L*right*solution;
        if(!green.allFinite())return false;
        result.green=green.cast<DataType>();result.solve_residual=double(residual);
        result.canonical_frame_maintenance=canonicalFrames;
        if(canonicalFrames)++pureProjectorGreenPrecisionStats.canonical_ld_rebuild_assemblies;
        result.solve_refinement_steps=5;result.extended_precision_solve=true;
        return true;
    };
    // A full-rank overlap can amplify the SVD's decomposition error in the
    // inverse action. Retry the SAME linear system with pivoted LU and mixed
    // precision refinement only after the original solve fails its gate.
    // Rank/rcond still come from the original SVD. Double candidates retain
    // the original residual check and scale; the final extended solve below
    // retains its precision through certification and Green assembly, with
    // independently bounded canonical-frame maintenance when eligible.
    bool extendedGreenReady=false;
    if(pureProjectorMatrixFinite(solved)&&std::isfinite(result.solve_residual)&&
       result.solve_residual>solveTolerance) {
        if(!options.prefer_lu_green_solve)lu.compute(overlap);
        MatType candidate=options.prefer_lu_green_solve?solved:MatType(lu.solve(rightHandSide));
        ExtendedMatrix extendedOverlap,extendedRhs;
        for(int step=1;step<=4;++step) {
            if(!pureProjectorMatrixFinite(candidate))break;
            const double residual=(overlap*candidate-rightHandSide).norm()/
                rightHandSideScale;
            if(std::isfinite(residual)&&residual<result.solve_residual) {
                solved=candidate;result.solve_residual=residual;
                result.solve_refinement_steps=step;
            }
            if(result.solve_residual<=solveTolerance)break;
            if(step==4)break;
            if(extendedOverlap.size()==0) {
                extendedOverlap=overlap.cast<ExtendedComplex>();
                extendedRhs=rightHandSide.cast<ExtendedComplex>();
            }
            const ExtendedMatrix defect=extendedRhs-extendedOverlap*candidate.cast<ExtendedComplex>();
            const MatType correction=lu.solve(defect.cast<DataType>().eval());
            if(!pureProjectorMatrixFinite(correction))break;
            candidate+=correction;
        }
        // Retain precision through overlap formation and Green assembly when
        // double solution rounding cannot pass the original solve gate.
        if(result.solve_residual>solveTolerance)
            extendedGreenReady=extendedSolve();
    }
    if (!pureProjectorMatrixFinite(solved) || !std::isfinite(result.solve_residual) ||
        result.solve_residual > solveTolerance) {
        result.status = canonicalMaintenanceFailed?extendedFailureStatus:
            PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    if(!extendedGreenReady)
        result.green = MatType::Identity(phiRight.rows(), phiRight.rows()) -
                       2.0 * phiRight * solved;
    const double greenScale=std::max(1.0,std::sqrt(double(result.green.rows())));
    bool greenFinite=false;
    auto residualTrusted=[&]() {
        if(result.green_residual<=options.green_residual_tolerance)return true;
        if(!options.allow_oblique_green)return false;
        const double normScale=std::max(greenScale,result.green.norm());
        if(!std::isfinite(normScale))return false;
        result.green_relative_residual=result.green_residual*(greenScale/normScale);
        return std::isfinite(result.green_relative_residual)&&
            result.green_relative_residual<=options.green_residual_tolerance;
    };
    auto certifyGreen=[&]() {
        greenFinite=pureProjectorMatrixFinite(result.green);
        if(!greenFinite)return false;
        result.green_residual=(result.green*result.green-
            MatType::Identity(result.green.rows(),result.green.cols())).norm()/greenScale;
        result.green_residual_double=result.green_residual;
        if(!std::isfinite(result.initial_green_residual))result.initial_green_residual=result.green_residual;
        result.long_double_green_residual=false;
        if(!std::isfinite(result.green_residual))return false;
        if(residualTrusted())return true;
        if(std::numeric_limits<long double>::digits>std::numeric_limits<double>::digits) {
            ++pureProjectorGreenPrecisionStats.long_double_residual_checks;
            result.green_residual=pureProjectorStoredGreenResidualExtended(result.green);
            result.long_double_green_residual=true;
            if(!std::isfinite(result.green_residual))return false;
            const bool trusted=residualTrusted();
            if(trusted)++pureProjectorGreenPrecisionStats.long_double_residual_passes;
            return trusted;
        }
        return false;
    };
    auto certifyStructure=[&]() {
        result.green_skew_residual=(result.green+result.green.transpose()).norm()/greenScale;
        result.green_diagonal_residual=result.green.diagonal().norm()/greenScale;
        return std::isfinite(result.green_skew_residual)&&std::isfinite(result.green_diagonal_residual)&&
            result.green_skew_residual<=options.green_residual_tolerance&&
            result.green_diagonal_residual<=options.green_residual_tolerance;
    };
    bool greenTrusted=certifyGreen();
    bool structureTrusted=greenFinite&&certifyStructure();
    const bool canonicalStructureRecovery=canonicalFrames&&!structureTrusted;
    if((!greenTrusted||canonicalStructureRecovery)&&greenFinite&&
       std::isfinite(result.green_residual)&&!result.extended_precision_solve) {
        if(extendedSolve()) {
            greenTrusted=certifyGreen();structureTrusted=greenFinite&&certifyStructure();
        } else if(canonicalMaintenanceFailed||result.canonical_frame_maintenance_attempted) {
            result.status=canonicalMaintenanceFailed?extendedFailureStatus:
                PureProjectorStatus::solve_residual_exceeded;return result;
        }
    }
    if(canonicalMaintenanceFailed) {
        result.status=extendedFailureStatus;return result;
    }
    if(!greenTrusted) {
        result.status=PureProjectorStatus::green_residual_exceeded;return result;
    }
    if(!std::isfinite(result.green_skew_residual)||!std::isfinite(result.green_diagonal_residual)||
       (options.fail_on_green_structure&&!structureTrusted)) {
        result.status=PureProjectorStatus::green_structure_exceeded;return result;
    }
    if(result.canonical_frame_maintenance)
        ++pureProjectorGreenPrecisionStats.canonical_ld_rebuild_certified;
    result.status = PureProjectorStatus::success;
    return result;
}

inline PureProjectorGreenResult pureProjectorGreenThinQr(
    const MatType &phiRight, const MatType &phiLeft,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    const ThinQrResult rightQr = thinQrSubspace(phiRight, options, true);
    if (!rightQr.ok()) {
        PureProjectorGreenResult result;
        result.status = rightQr.status;
        return result;
    }
    const ThinQrResult leftQr = thinQrSubspace(phiLeft, options, true);
    if (!leftQr.ok()) {
        PureProjectorGreenResult result;
        result.status = leftQr.status;
        return result;
    }
    PureProjectorOptions gaussianOptions=options;gaussianOptions.fail_on_green_structure=true;
    return pureProjectorGreen(rightQr.q, leftQr.q, gaussianOptions);
}

// A1: the LEFT basis is the immutable trial state, identical for every slice
// and every proposal. Its pivoted thin QR is prepared once and reused. The
// Green function I-2 Q_R (Q_L^dag Q_R)^{-1} Q_L^dag is invariant under
// independent unitary rephasings Q->Q O of either Stiefel basis, so caching
// this QR and feeding an already-orthonormal right basis changes no physics;
// every rank/rcond/solve/involution gate of pureProjectorGreen() is retained.
struct PreparedGreenTrial {
    MatType qLeft;
    PureProjectorStatus status = PureProjectorStatus::invalid_dimension;
    int rank = 0;
    double orthonormal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

inline PreparedGreenTrial prepareGreenTrial(
        const GaussianTrialState &trial, const PureProjectorOptions &options) {
    PreparedGreenTrial prepared;
    const ThinQrResult leftQr = thinQrSubspace(trial.Phi, options, true);
    if (!leftQr.ok()) { prepared.status = leftQr.status; return prepared; }
    prepared.qLeft = leftQr.q;
    prepared.status = PureProjectorStatus::success;
    prepared.rank = leftQr.rank;
    prepared.orthonormal_residual = leftQr.orthonormal_residual;
    return prepared;
}

// Right basis is supplied already orthonormal (it is the explicit thin Q of
// the per-slice Householder QR). Verify it with the SAME gate used inside
// thinQrSubspace; any discrepancy, dimension mismatch or failed preparation
// falls back to the full two-sided thin-QR path, preserving old behavior.
inline PureProjectorGreenResult pureProjectorGreenPreparedRight(
        const MatType &qRight, const GaussianTrialState &trial,
        const PreparedGreenTrial &prepared, const PureProjectorOptions &options) {
    if (!prepared.ok() ||
        qRight.rows() != prepared.qLeft.rows() ||
        qRight.cols() != prepared.qLeft.cols() ||
        qRight.rows() <= 0 || qRight.cols() <= 0 ||
        qRight.rows() < qRight.cols() ||
        !pureProjectorMatrixFinite(qRight)) {
        return pureProjectorGreenThinQr(qRight, trial.Phi, options);
    }
    const double orth =
        (qRight.adjoint() * qRight -
         MatType::Identity(qRight.cols(), qRight.cols())).norm() /
        std::max(1.0, std::sqrt(double(qRight.cols())));
    if (!std::isfinite(orth) || orth > options.solve_residual_tolerance ||
        !pureProjectorCanonicalFrameTrusted(qRight,options.solve_residual_tolerance) ||
        !pureProjectorCanonicalFrameTrusted(prepared.qLeft,options.solve_residual_tolerance)) {
        return pureProjectorGreenThinQr(qRight, trial.Phi, options);
    }
    PureProjectorOptions gaussianOptions=options;gaussianOptions.fail_on_green_structure=true;
    const PureProjectorGreenResult preparedResult =
        pureProjectorGreen(qRight, prepared.qLeft, gaussianOptions);
    if (preparedResult.ok()) return preparedResult;

    // The prepared QR is mathematically equivalent to re-orthogonalizing the
    // trial basis, but the two factorizations can make different pivoting
    // decisions near an ill-conditioned overlap.  Retry with the frozen
    // uncached two-sided QR only when the prepared route fails.  A fallback is
    // accepted only if it passes every original rank, solve, and Green gate;
    // otherwise retain the original failure rather than relaxing a criterion.
    const PureProjectorGreenResult uncachedResult =
        pureProjectorGreenThinQr(qRight, trial.Phi, options);
    return uncachedResult.ok() ? uncachedResult : preparedResult;
}
// Only reconsider a purity-residual failure of an otherwise certified solve.
// This changes the check's scale, never the Green matrix or the update formula.
inline PureProjectorGreenResult pureProjectorCheckObliqueGreen(
    const PureProjectorGreenResult& failed,const PureProjectorOptions& options) {
    if(failed.status!=PureProjectorStatus::green_residual_exceeded||
       failed.green.rows()<=0||failed.green.rows()!=failed.green.cols()||
       failed.overlap_rank*2!=failed.green.rows()||
       !pureProjectorMatrixFinite(failed.green)||
       !std::isfinite(failed.overlap_rcond)||failed.overlap_rcond<options.minimum_overlap_rcond||
       !std::isfinite(failed.solve_residual)||
       failed.solve_residual>pureProjectorSolveTolerance(options,failed.overlap_rcond)||
       !std::isfinite(failed.green_residual)||failed.green_residual<0)return failed;
    PureProjectorGreenResult result=failed;
    const double dimensionScale=std::max(1.,std::sqrt(double(result.green.rows())));
    const double greenScale=std::max(dimensionScale,result.green.norm());
    if(!std::isfinite(greenScale))return failed;
    result.green_relative_residual=result.green_residual*(dimensionScale/greenScale);
    if(!std::isfinite(result.green_relative_residual)||
       result.green_relative_residual>options.green_residual_tolerance)return failed;
    result.green_skew_residual=(result.green+result.green.transpose()).norm()/dimensionScale;
    result.green_diagonal_residual=result.green.diagonal().norm()/dimensionScale;
    if(!std::isfinite(result.green_skew_residual)||!std::isfinite(result.green_diagonal_residual)||
       (options.fail_on_green_structure&&(result.green_skew_residual>options.green_residual_tolerance||
                                        result.green_diagonal_residual>options.green_residual_tolerance)))return failed;
    result.status=PureProjectorStatus::success;
    return result;
}

inline bool propagateKet(const MatType &factor, MatType &subspace) {
    if (factor.rows() <= 0 || factor.rows() != factor.cols() ||
        factor.cols() != subspace.rows() || !pureProjectorMatrixFinite(factor) ||
        !pureProjectorMatrixFinite(subspace)) return false;
    subspace = factor * subspace;
    return pureProjectorMatrixFinite(subspace);
}

inline bool propagateBra(const MatType &factor, MatType &subspace) {
    if (factor.rows() <= 0 || factor.rows() != factor.cols() ||
        factor.cols() != subspace.rows() || !pureProjectorMatrixFinite(factor) ||
        !pureProjectorMatrixFinite(subspace)) return false;
    subspace = factor.adjoint() * subspace;
    return pureProjectorMatrixFinite(subspace);
}

inline bool propagateKet(Operator &factor, MatType &subspace) {
    if (!pureProjectorMatrixFinite(subspace)) return false;
    MatType output;
    factor.left_multiply(subspace, output);
    if (output.rows() != subspace.rows() || output.cols() != subspace.cols() ||
        !pureProjectorMatrixFinite(output)) return false;
    subspace.swap(output);
    return true;
}

inline bool propagateBra(Operator &factor, MatType &subspace) {
    if (!pureProjectorMatrixFinite(subspace)) return false;
    MatType output;
    factor.adjoint_left_multiply(subspace, output);
    if (output.rows() != subspace.rows() || output.cols() != subspace.cols() ||
        !pureProjectorMatrixFinite(output)) return false;
    subspace.swap(output);
    return true;
}

#endif
