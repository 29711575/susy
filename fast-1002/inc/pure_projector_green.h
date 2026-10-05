#ifndef PURE_PROJECTOR_GREEN_H
#define PURE_PROJECTOR_GREEN_H

#include "operator.h"
#include "types.h"
#include "pure_projector_svd.h"
#include "pure_projector_factor.h"
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
    double green_residual_tolerance = 1e-10;
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
    double green_residual = std::numeric_limits<double>::infinity();
    // Set only by the local oblique-Green recovery; the old diagnostic keeps
    // its dimension-normalized meaning for strict/reference callers.
    double green_relative_residual = std::numeric_limits<double>::infinity();
    double green_skew_residual = std::numeric_limits<double>::infinity();
    double green_diagonal_residual = std::numeric_limits<double>::infinity();
    bool ok() const { return status == PureProjectorStatus::success; }
};

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

inline ThinQrResult thinQrSubspace(
    const MatType &subspace,
    const PureProjectorOptions &options = PureProjectorOptions()) {
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
    PureGreenSvd svd(overlap,options.accelerated_linear_algebra);
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

    const MatType rightHandSide = phiLeft.adjoint();
    const MatType solved = svd.solve(rightHandSide);
    result.solve_residual = (overlap * solved - rightHandSide).norm() /
        std::max(1.0, rightHandSide.norm());
    if (!pureProjectorMatrixFinite(solved) ||
        result.solve_residual > options.solve_residual_tolerance) {
        result.status = PureProjectorStatus::solve_residual_exceeded;
        return result;
    }
    result.green = MatType::Identity(phiRight.rows(), phiRight.rows()) -
                   2.0 * phiRight * solved;
    result.green_residual =
        (result.green * result.green -
         MatType::Identity(result.green.rows(), result.green.cols())).norm() /
        std::max(1.0, std::sqrt(double(result.green.rows())));
    if (!pureProjectorMatrixFinite(result.green) ||
        result.green_residual > options.green_residual_tolerance) {
        result.status = PureProjectorStatus::green_residual_exceeded;
        return result;
    }
    const double greenScale =
        std::max(1.0, std::sqrt(double(result.green.rows())));
    result.green_skew_residual =
        (result.green + result.green.transpose()).norm() / greenScale;
    result.green_diagonal_residual =
        result.green.diagonal().norm() / greenScale;
    if (!std::isfinite(result.green_skew_residual) ||
        !std::isfinite(result.green_diagonal_residual) ||
        (options.fail_on_green_structure &&
         (result.green_skew_residual > options.green_residual_tolerance ||
          result.green_diagonal_residual > options.green_residual_tolerance))) {
        result.status = PureProjectorStatus::green_structure_exceeded;
        return result;
    }
    result.status = PureProjectorStatus::success;
    return result;
}

inline PureProjectorGreenResult pureProjectorGreenThinQr(
    const MatType &phiRight, const MatType &phiLeft,
    const PureProjectorOptions &options = PureProjectorOptions()) {
    const ThinQrResult rightQr = thinQrSubspace(phiRight, options);
    if (!rightQr.ok()) {
        PureProjectorGreenResult result;
        result.status = rightQr.status;
        return result;
    }
    const ThinQrResult leftQr = thinQrSubspace(phiLeft, options);
    if (!leftQr.ok()) {
        PureProjectorGreenResult result;
        result.status = leftQr.status;
        return result;
    }
    return pureProjectorGreen(rightQr.q, leftQr.q, options);
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
    const ThinQrResult leftQr = thinQrSubspace(trial.Phi, options);
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
    if (!std::isfinite(orth) || orth > options.solve_residual_tolerance) {
        return pureProjectorGreenThinQr(qRight, trial.Phi, options);
    }
    const PureProjectorGreenResult preparedResult =
        pureProjectorGreen(qRight, prepared.qLeft, options);
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
       !std::isfinite(failed.solve_residual)||failed.solve_residual>options.solve_residual_tolerance||
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
