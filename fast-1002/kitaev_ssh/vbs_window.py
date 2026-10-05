#!/usr/bin/env python3
"""Finite-window structure-factor identities for ordinary square PBC.

The production MC estimators remain in vbs_observables.h.  This module is an
offline matrix/statistics helper for the exact relation between an off-grid
finite-window Fourier estimator and the complete PBC Fourier grid.
"""
from __future__ import annotations

import numpy as np


def _validate_bond_length(bond_length: float):
    if not np.isfinite(bond_length) or bond_length <= 0:
        raise ValueError("bond_length must be positive and finite")


def reciprocal_geometry(bond_length: float = 1.0):
    """Return project cell/reciprocal vectors in the stated bond-length unit."""
    _validate_bond_length(bond_length)
    length = float(bond_length)
    a1 = np.array([np.sqrt(3.0) * length / 2, 3.0 * length / 2])
    a2 = np.array([-np.sqrt(3.0) * length / 2, 3.0 * length / 2])
    b1 = np.array([2 * np.pi / (np.sqrt(3.0) * length), 2 * np.pi / (3 * length)])
    b2 = np.array([-2 * np.pi / (np.sqrt(3.0) * length), 2 * np.pi / (3 * length)])
    return a1, a2, b1, b2


def window_factor(L: int, p_fraction, bond_length: float = 1.0) -> complex:
    """W_L(p)=N^-1 sum_(x,y) exp[-2pi*i*(u*x+v*y)].

    `bond_length` is explicit geometry metadata/validation only; this fractional
    phase performs no Cartesian length conversion.
    """
    if int(L) != L or L < 1:
        raise ValueError("L must be a positive integer")
    p = np.asarray(p_fraction, dtype=float)
    if p.shape != (2,) or not np.all(np.isfinite(p)):
        raise ValueError("p_fraction must contain two finite reciprocal fractions")
    _validate_bond_length(bond_length)
    x, y = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
    return complex(np.mean(np.exp(-2j * np.pi * (p[0] * x + p[1] * y))))


def square_window_amplitudes(L: int, Q, bond_length: float = 1.0):
    """Return complex W_L(Q-k) amplitudes on the ordinary square-PBC grid."""
    q = np.asarray(Q, dtype=float)
    if q.shape != (2,) or not np.all(np.isfinite(q)):
        raise ValueError("Q must contain two finite reciprocal fractions")
    if int(L) != L or L < 2:
        raise ValueError("L must be an integer >= 2")
    _validate_bond_length(bond_length)
    amplitudes = np.empty((L, L), dtype=complex)
    for n1 in range(L):
        for n2 in range(L):
            k = np.array([n1 / L, n2 / L], dtype=float)
            amplitudes[n1, n2] = window_factor(L, q - k, bond_length)
    return amplitudes


def square_window_weights(L: int, Q, bond_length: float = 1.0):
    """Weights |W_L(Q-k)|^2 for all allowed square-PBC k, stable grid order."""
    q = np.asarray(Q, dtype=float)
    if q.shape != (2,) or not np.all(np.isfinite(q)):
        raise ValueError("Q must contain two finite reciprocal fractions")
    if int(L) != L or L < 2:
        raise ValueError("L must be an integer >= 2")
    amplitudes = square_window_amplitudes(L, q, bond_length)
    rows = []
    weights = np.abs(amplitudes) ** 2
    for n1 in range(L):
        for n2 in range(L):
            k = np.array([n1 / L, n2 / L], dtype=float)
            w = amplitudes[n1, n2]
            weight = float(weights[n1, n2])
            rows.append({"n1": n1, "n2": n2, "k_fraction": k.tolist(),
                         "W_real": float(w.real), "W_imag": float(w.imag),
                         "weight": weight})
    return weights, rows


def windowed_matrix(S_grid, weights):
    """Apply the weights to 3x3 S(k) matrices before any eigensolve."""
    s = np.asarray(S_grid, dtype=complex)
    w = np.asarray(weights, dtype=float)
    if s.shape != (*w.shape, 3, 3):
        raise ValueError("S_grid must have shape weights.shape + (3,3)")
    if not np.all(np.isfinite(s)) or not np.all(np.isfinite(w)):
        raise ValueError("non-finite structure factor or window weight")
    return np.einsum("xy,xyij->ij", w, s, optimize=True)


def largest_hermitian_eigenvalue(matrix):
    """Largest eigenvalue of the Hermitian part, matching the QMC convention."""
    m = np.asarray(matrix, dtype=complex)
    if m.shape != (3, 3) or not np.all(np.isfinite(m)):
        raise ValueError("matrix must be finite 3x3")
    h = (m + m.conj().T) / 2
    return float(np.linalg.eigvalsh(h)[-1])


def ensemble_connected_matrix(S_raw, mean_F, Ncell: int):
    """Subtract the outer product of the ensemble Fourier mean from raw S."""
    s = np.asarray(S_raw, dtype=complex)
    f = np.asarray(mean_F, dtype=complex)
    if s.shape != (3, 3) or f.shape != (3,) or Ncell <= 0:
        raise ValueError("expected a 3x3 S, 3-vector mean_F and positive Ncell")
    return s - np.outer(f, f.conj()) / float(Ncell)


def translation_averaged_connected_matrix(S_translation_raw, mean_Delta_zero,
                                          W_Q, Ncell: int):
    """Connected matrix after averaging a fixed-origin estimator over shifts.

    ``mean_Delta_zero`` is the phase-reweighted ensemble mean of F(0)/Ncell.
    Spatial translation averaging makes the mean Fourier amplitude at Q equal
    to W_Q*mean_F(0); consequently the connected subtraction is
    ``|W_Q|^2 <F(0)><F(0)>^dagger/Ncell``.  The subtraction is deliberately
    applied to the already-windowed matrix, not separately to each k bin.
    """
    s = np.asarray(S_translation_raw, dtype=complex)
    d0 = np.asarray(mean_Delta_zero, dtype=complex)
    if s.shape != (3, 3) or d0.shape != (3,) or Ncell <= 0:
        raise ValueError("expected a 3x3 matrix, 3-vector zero-mode mean and positive Ncell")
    wq = complex(W_Q)
    if not np.all(np.isfinite(s)) or not np.all(np.isfinite(d0)) or not np.isfinite(wq):
        raise ValueError("non-finite translation-connected input")
    return s - float(Ncell) * abs(wq)**2 * np.outer(d0, d0.conj())


def configuration_spatially_centered_F(F_Q, F_zero, W_Q):
    """Per-configuration uniform-mode projection used by selected `Sprj`."""
    f = np.asarray(F_Q, dtype=complex)
    f0 = np.asarray(F_zero, dtype=complex)
    if f.shape != (3,) or f0.shape != (3,):
        raise ValueError("F_Q and F_zero must be three-component vectors")
    return f - complex(W_Q) * f0


def joint_phase_ratio_window_jackknife(block_phase, block_numerators, weights,
                                       block_first_moments=None,
                                       window_amplitudes=None, Ncell=None):
    """Joint delete-one-block ratio/jackknife for all grid momenta at once.

    `block_numerators[b,k]` is the sum of phase* S_k over block b.  The same
    phase sum and omitted block are used for every momentum, retaining their
    covariance. Each replica first constructs the complete weighted 3x3
    matrix, then takes its leading Hermitian eigenvalue. If first-moment block
    numerators are supplied, they contain phase-weighted Delta(k)=F(k)/Ncell;
    the connected result subtracts Ncell*mean_Delta(Q) mean_Delta(Q)^dagger in
    the same full sample and every leave-one-block replica.
    """
    phase = np.asarray(block_phase, dtype=complex)
    numerator = np.asarray(block_numerators, dtype=complex)
    w = np.asarray(weights, dtype=float)
    if phase.ndim != 1 or len(phase) < 2:
        raise ValueError("at least two phase blocks are required")
    if numerator.shape != (len(phase), *w.shape, 3, 3):
        raise ValueError("block_numerators shape must be (blocks,*weights.shape,3,3)")
    first = None
    amplitudes = None
    if block_first_moments is not None or window_amplitudes is not None or Ncell is not None:
        if block_first_moments is None or window_amplitudes is None or Ncell is None or Ncell <= 0:
            raise ValueError("connected estimator requires first moments, window amplitudes and positive Ncell")
        first = np.asarray(block_first_moments, dtype=complex)
        amplitudes = np.asarray(window_amplitudes, dtype=complex)
        if first.shape != (len(phase), *w.shape, 3) or amplitudes.shape != w.shape:
            raise ValueError("first-moment/window-amplitude dimensions do not match weights")
        if not np.all(np.isfinite(first)) or not np.all(np.isfinite(amplitudes)):
            raise ValueError("non-finite first-moment/window-amplitude data")
    if (not np.all(np.isfinite(phase)) or not np.all(np.isfinite(numerator))
            or not np.all(np.isfinite(w))):
        raise ValueError("non-finite block data")
    total_phase = phase.sum()
    if abs(total_phase) <= 1e-14:
        raise ValueError("full phase denominator is unresolved")
    total_num = numerator.sum(axis=0)
    point_matrices = total_num / total_phase
    matrix = windowed_matrix(point_matrices, w)
    connected_matrix = None
    mean_delta_q = None
    if first is not None:
        mean_delta_grid = first.sum(axis=0) / total_phase
        mean_delta_q = np.einsum("xy,xyi->i", amplitudes, mean_delta_grid, optimize=True)
        connected_matrix = matrix - float(Ncell) * np.outer(mean_delta_q, mean_delta_q.conj())
    replica_matrices = []
    replica_connected_matrices = []
    replica_point_connected_matrices = []
    replica_lambdas = []
    replica_szz = []
    for b in range(len(phase)):
        denominator = total_phase - phase[b]
        if abs(denominator) <= 1e-14:
            raise ValueError("leave-one-block phase denominator is unresolved")
        per_k = (total_num - numerator[b]) / denominator
        wm = windowed_matrix(per_k, w)
        replica_matrices.append(wm)
        replica_lambdas.append(largest_hermitian_eigenvalue(wm))
        replica_szz.append(per_k[..., 2, 2].real.reshape(-1))
        if first is not None:
            delta_grid = (first.sum(axis=0) - first[b]) / denominator
            delta_q = np.einsum("xy,xyi->i", amplitudes, delta_grid, optimize=True)
            wconn = wm - float(Ncell) * np.outer(delta_q, delta_q.conj())
            replica_connected_matrices.append(wconn)
            point_conn = per_k - float(Ncell) * np.einsum(
                "xyi,xyj->xyij", delta_grid, delta_grid.conj(), optimize=True)
            replica_point_connected_matrices.append(point_conn)
    replica_lambdas = np.asarray(replica_lambdas)
    lambda_center = replica_lambdas.mean()
    lambda_error = float(np.sqrt((len(phase) - 1) / len(phase)
                                  * np.sum((replica_lambdas - lambda_center) ** 2)))
    replica_szz = np.asarray(replica_szz)
    centered = replica_szz - replica_szz.mean(axis=0)
    covariance = (len(phase) - 1) / len(phase) * centered.T @ centered
    result = {"matrix": matrix, "point_matrices": point_matrices,
            "leading_eigenvalue": largest_hermitian_eigenvalue(matrix),
            "leading_eigenvalue_jackknife_error": lambda_error,
            "replica_matrices": np.asarray(replica_matrices),
            "replica_leading_eigenvalues": replica_lambdas,
            "szz_real_jackknife_covariance": covariance}
    if first is not None:
        rc = np.asarray(replica_connected_matrices)
        eigs = np.asarray([largest_hermitian_eigenvalue(x) for x in rc])
        center = eigs.mean()
        rpc = np.asarray(replica_point_connected_matrices)
        point_szz = rpc[..., 2, 2].real.reshape(len(phase), -1)
        point_centered = point_szz - point_szz.mean(axis=0)
        point_covariance = (len(phase) - 1) / len(phase) * point_centered.T @ point_centered
        result.update({"connected_matrix": connected_matrix,
                       "point_connected_matrices": point_matrices - float(Ncell) * np.einsum(
                           "xyi,xyj->xyij", mean_delta_grid, mean_delta_grid.conj(), optimize=True),
                       "connected_leading_eigenvalue": largest_hermitian_eigenvalue(connected_matrix),
                       "connected_leading_eigenvalue_jackknife_error": float(np.sqrt(
                           (len(phase)-1)/len(phase)*np.sum((eigs-center)**2))),
                       "mean_Delta_Q": mean_delta_q,
                       "replica_connected_matrices": rc,
                       "replica_point_connected_matrices": rpc,
                       "connected_szz_real_jackknife_covariance": point_covariance,
                       "replica_connected_leading_eigenvalues": eigs})
    return result


def connected_window_association_ratios(block_phase, block_S_numerators,
                                        block_Delta_numerators, point_indices,
                                        L=None, min_blocks=8, sigma_threshold=2.0,
                                        Ncell=None):
    """Joint block estimates of the two directional window structure-factor ratios.

    `point_indices` has `center`, `plus_b1`, `minus_b1`, `plus_b2`, `minus_b2`
    integer slots into the common selected-q table. Every q uses the same
    delete-one block and phase denominator. Negative ratios are retained.
    """
    phase = np.asarray(block_phase, dtype=complex)
    sn = np.asarray(block_S_numerators, dtype=complex)
    dn = np.asarray(block_Delta_numerators, dtype=complex)
    if phase.ndim != 1 or len(phase) < 2 or sn.ndim != 4 or sn.shape[0] != len(phase):
        raise ValueError("invalid joint selected block dimensions")
    if sn.shape[2:] != (3, 3) or dn.shape != sn.shape[:2] + (3,):
        raise ValueError("selected S and Delta blocks must have shapes (B,Nq,3,3)/(B,Nq,3)")
    required = {"center", "plus_b1", "minus_b1", "plus_b2", "minus_b2"}
    if set(point_indices) != required or any(i < 0 or i >= sn.shape[1] for i in point_indices.values()):
        raise ValueError("point_indices must identify center and both +/- neighbors in each direction")
    if Ncell is None:
        if L is None or int(L) != L or L < 2:
            raise ValueError("L must be an integer >= 2")
        nc = int(L) * int(L)
    else:
        if int(Ncell) != Ncell or Ncell < 1:
            raise ValueError("Ncell must be a positive integer")
        nc = int(Ncell)
    total_phase = phase.sum()
    if abs(total_phase) <= 1e-14:
        raise ValueError("full phase denominator is unresolved")
    total_s = sn.sum(axis=0)
    total_d = dn.sum(axis=0)

    def connected_from(snum, dnum, den):
        s = snum / den
        delta = dnum / den
        return s - nc * delta[..., :, None] * delta[..., None, :].conj()

    full = connected_from(total_s, total_d, total_phase)
    values = {}
    replica_values = {"b1": [], "b2": []}
    replica_center = []
    replica_connected_szz = []

    def ratio(conn, axis):
        c = point_indices["center"]
        p = point_indices[f"plus_{axis}"]
        m = point_indices[f"minus_{axis}"]
        # The observable is the real fixed-z diagonal channel.  Retain the
        # complex matrix and its diagnostics separately, but do not let a tiny
        # imaginary diagonal component change R_win.
        denom = float(conn[c, 2, 2].real)
        numerator = float(conn[p, 2, 2].real + conn[m, 2, 2].real)
        return (complex(1.0 - numerator / (2.0 * denom), 0.0)
                if abs(denom) > 1e-14 else complex(np.nan, np.nan))

    for axis in ("b1", "b2"):
        values[axis] = ratio(full, axis)
    replica_center = []
    for b in range(len(phase)):
        den = total_phase - phase[b]
        if abs(den) <= 1e-14:
            raise ValueError("leave-one-block phase denominator is unresolved")
        conn = connected_from(total_s-sn[b], total_d-dn[b], den)
        replica_connected_szz.append(conn[:, 2, 2].real)
        replica_center.append(conn[point_indices["center"], 2, 2])
        for axis in ("b1", "b2"):
            replica_values[axis].append(ratio(conn, axis))

    def complex_jk(samples):
        samples = np.asarray(samples, dtype=complex)
        center = samples.mean()
        d = samples - center
        factor = (len(samples)-1)/len(samples)
        cov = factor * np.array([[np.sum(d.real*d.real),np.sum(d.real*d.imag)],
                                 [np.sum(d.real*d.imag),np.sum(d.imag*d.imag)]])
        return center, np.sqrt(np.maximum(np.diag(cov),0.0)), cov

    center_val, center_se, center_cov = complex_jk(replica_center)
    connected_szz_replicas = np.asarray(replica_connected_szz)
    connected_szz_centered = connected_szz_replicas - connected_szz_replicas.mean(axis=0)
    connected_szz_covariance = ((len(phase)-1)/len(phase)
                                * connected_szz_centered.T @ connected_szz_centered)
    center_est = full[point_indices["center"], 2, 2]
    center_C = float(center_est.real)
    center_C_se = float(center_se[0])
    denom_error = float(np.hypot(*center_se))
    results = {}
    for axis in ("b1", "b2"):
        r = values[axis]
        rmean, rse, rcov = complex_jk(replica_values[axis])
        enough = len(phase) >= min_blocks
        resolved = (np.isfinite(r.real) and np.isfinite(r.imag)
                    and np.isfinite(center_C)
                    and abs(center_C) > sigma_threshold * center_C_se
                    and abs(center_C) > 1e-14)
        results[axis] = {"available": bool(enough and resolved),
                         "unavailable_reason": None if enough and resolved else
                            ("too_few_equal_blocks" if not enough else "connected_center_denominator_unresolved"),
                         "estimate": [float(r.real),float(r.imag)] if np.isfinite(r.real) and np.isfinite(r.imag) else None,
                         "stderr_real": float(rse[0]), "stderr_imag": float(rse[1]),
                         "covariance_real_imag": float(rcov[0,1]),
                         "center_Szz_connected": [float(center_est.real),float(center_est.imag)],
                         "center_C": center_C,
                         "center_C_jackknife_se": center_C_se,
                         "m2_z": center_C / nc,
                         "m2_z_jackknife_se": center_C_se / nc,
                         "center_Szz_complex_jackknife_error": denom_error,
                         "n_blocks": int(len(phase)), "min_blocks": int(min_blocks),
                         "sigma_threshold": float(sigma_threshold),
                         "negative_real_allowed": True}
    return {"ratios": results, "connected_selected_matrices": full,
            "center_Szz_connected_jackknife_center": [float(center_val.real),float(center_val.imag)],
            "center_Szz_connected_jackknife_covariance": center_cov,
            "connected_szz_real_jackknife_covariance_across_selected_q": connected_szz_covariance,
            "replica_connected_szz_by_q": connected_szz_replicas}
