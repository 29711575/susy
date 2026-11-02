#!/usr/bin/env python3
"""Static fixed-background VBS curvature, not a QMC estimator.

E_matter=-sum(singular_values(T))=-Tr|iA|/4 with A_AB=2T.
Real Fourier modes have sum(v_alpha*v_beta)/Ncell=delta_alpha_beta.
Continuous-q references evaluate k+q directly, without snapping to a grid.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from predict_vbs_momenta import (allowed_grid, dirac_nodes, grid_point_for,
                                 momentum_distance, q_predictions,
                                 reciprocal_cartesian, uniform_background_solution)


def _parameters(jeff, g, spring_k=1.):
    jeff = np.asarray(jeff, dtype=float)
    if jeff.shape != (3,) or not np.all(np.isfinite(jeff)):
        raise ValueError("Jeff must contain three finite values")
    if not math.isfinite(g) or not math.isfinite(spring_k) or spring_k <= 0:
        raise ValueError("g must be finite and spring_k finite and positive")
    return jeff


def _geometry(lx, ly, shear):
    if any(not isinstance(v, (int, np.integer)) for v in (lx, ly, shear)):
        raise ValueError("geometry must use integers")
    if lx < 2 or ly < 2 or not 0 <= shear < lx:
        raise ValueError("require Lx,Ly>=2 and 0<=shear<Lx")


def _q_fraction(q):
    q = np.asarray(q, dtype=float)
    if q.shape != (2,) or not np.all(np.isfinite(q)):
        raise ValueError("q must have two finite fractional coordinates")
    return tuple(q % 1.)


def _mode_labels(q):
    phases = ("cos",) if all(abs(2*z-round(2*z)) < 1e-12 for z in q) else ("cos", "sin")
    return [(mu, phase) for mu in range(3) for phase in phases]


def cell_reduce(x, y, lx, ly, shear):
    q, r = divmod(y, ly)
    p, s = divmod(x-q*shear, lx)
    return s, r, p, q


def bond_table(lx, ly, shear):
    _geometry(lx, ly, shear)
    bonds = []
    for x in range(lx):
        for y in range(ly):
            for mu, (dx, dy) in enumerate(((-1, 0), (0, -1), (0, 0))):
                bx, by, _, _ = cell_reduce(x+dx, y+dy, lx, ly, shear)
                # Uniform zero-flux (++ Wilson) gauge, same as the Bloch f(k).
                bonds.append((x, y, mu, bx, by, 1))
    return bonds


def hopping_matrix(lx, ly, shear, jeff, delta=None, g=0.):
    jeff = _parameters(jeff, g)
    n = lx*ly
    delta = np.zeros((n, 3)) if delta is None else np.asarray(delta, dtype=float)
    if delta.shape != (n, 3) or not np.all(np.isfinite(delta)):
        raise ValueError("deltaX must be a finite Ncell by 3 array")
    t = np.zeros((n, n))
    for x, y, mu, bx, by, gauge in bond_table(lx, ly, shear):
        a, b = x*ly+y, bx*ly+by
        t[a, b] += gauge*(jeff[mu]+g*delta[a, mu])
    return t


def finite_matter_energy(lx, ly, shear, jeff, delta=None, g=0.):
    return float(-np.linalg.svd(hopping_matrix(lx, ly, shear, jeff, delta, g),
                              compute_uv=False).sum())


def mode_basis(lx, ly, shear, q):
    """Periodic orthonormal real modes; off-grid q is a reference only."""
    _geometry(lx, ly, shear)
    q = _q_fraction(q)
    point = grid_point_for(q, lx, ly, shear)
    if point is None:
        raise ValueError("finite-torus modulation q must be on the allowed sheared grid")
    q = point["fraction"]
    phases = np.array([2*math.pi*(q[0]*x+q[1]*y)
                       for x in range(lx) for y in range(ly)])
    modes = []
    for mu, name in _mode_labels(q):
        values = np.cos(phases) if name == "cos" else np.sin(phases)
        values = values/np.sqrt(np.mean(values*values))
        arr = np.zeros((lx*ly, 3)); arr[:, mu] = values
        modes.append((mu, name, arr))
    v = np.array([m[2].reshape(-1) for m in modes])
    if np.max(abs(v@v.T/(lx*ly)-np.eye(len(modes)))) > 1e-10:
        raise RuntimeError("Fourier modes failed orthonormality")
    return modes


def finite_difference_hessians(lx, ly, shear, jeff, xbar, g, q,
                              steps, spring_k=1., analytic=None):
    """Independent full mixed-channel energy differences at multiple steps."""
    _parameters(jeff, g, spring_k)
    xbar = np.asarray(xbar, dtype=float)
    if xbar.shape != (3,) or not np.all(np.isfinite(xbar)):
        raise ValueError("Xbar must have three finite values")
    modes = mode_basis(lx, ly, shear, q)
    n, d = lx*ly, len(modes)
    t0 = hopping_matrix(lx, ly, shear, jeff)
    dt = [hopping_matrix(lx, ly, shear, [0., 0., 0.], m[2], g) for m in modes]
    fields = np.array([m[2] for m in modes])
    def energy(coeff):
        t = t0 + np.einsum("i,ijk->jk", coeff, dt)
        x = xbar[None, :] + np.einsum("i,ijk->jk", coeff, fields)
        return float((-np.linalg.svd(t, compute_uv=False).sum()
                      + .5*spring_k*np.sum(x*x))/n)
    base = energy(np.zeros(d))
    mixed = np.arange(1., d+1.)*(-1.)**np.arange(d)
    mixed /= np.linalg.norm(mixed)
    results = []
    for step in steps:
        if not math.isfinite(step) or step <= 0:
            raise ValueError("finite difference steps must be finite and positive")
        h = np.zeros((d, d))
        for i in range(d):
            ei = np.eye(d)[i]*step
            h[i, i] = (energy(ei)-2*base+energy(-ei))/step**2
            for j in range(i):
                ej = np.eye(d)[j]*step
                h[i, j] = h[j, i] = (energy(ei+ej)-energy(ei-ej)
                                     -energy(-ei+ej)+energy(-ei-ej))/(4*step**2)
        mixed_value = (energy(step*mixed)-2*base+energy(-step*mixed))/step**2
        result = {"step": step, "total_hessian": h.tolist(),
                  "mixed_direction": mixed.tolist(), "mixed_curvature": mixed_value}
        if analytic is not None:
            result["max_abs_matrix_error"] = float(np.max(abs(h-analytic)))
            result["mixed_curvature_error"] = float(mixed_value-mixed@analytic@mixed)
        results.append(result)
    return {"method": "central energy differences including all mixed bond/quadrature entries",
            "checks": results,
            "max_matrix_step_change": float(max(
                (np.max(abs(np.array(a["total_hessian"])-np.array(b["total_hessian"])))
                 for a, b in zip(results, results[1:])), default=0.))}


def _response(matter, labels, spring_k, diagnostics):
    total = matter+spring_k*np.eye(len(labels))
    vals, vecs = np.linalg.eigh(total)
    return {"status": "certified", "matter_hessian": matter.tolist(),
            "phonon_hessian": (spring_k*np.eye(len(labels))).tolist(),
            "total_hessian": total.tolist(), "eigenvalues": vals.tolist(),
            "min_eigenvalue": float(vals[0]), "mode_count": len(labels),
            "basis": [{"channel": mu, "phase": phase} for mu, phase in labels],
            "soft_mode": [{"channel": mu, "phase": phase,
                           "coefficient": float(vecs[i, 0])}
                          for i, (mu, phase) in enumerate(labels)],
            "soft_subspace_dimension": int(np.sum(abs(vals-vals[0]) < 1e-10)),
            "soft_subspace_vectors": vecs[:, abs(vals-vals[0]) < 1e-10].T.tolist(),
            "normalization": "sum_(r,mu) v_alpha v_beta / Ncell = delta_alpha_beta",
            "energy_definition": "E_matter[Jeff+g deltaX] + spring_k/2 sum_b (Xbar_b+deltaX_b)^2, per cell",
            **diagnostics}


def static_hessian(lx, ly, shear, jeff, g, q, step=2e-4, spring_k=1.):
    """Occupied-band perturbation, checked independently by energy differences.

    The legacy step argument is used only for fallback energy scans when a
    zero mode prevents a certified Hessian, not for the analytic response.
    """
    jeff = _parameters(jeff, g, spring_k)
    if not math.isfinite(step) or step <= 0:
        raise ValueError("energy scan step must be finite and positive")
    modes = mode_basis(lx, ly, shear, q)
    labels = [(m[0], m[1]) for m in modes]
    n = lx*ly
    t = hopping_matrix(lx, ly, shear, jeff)
    h = np.zeros((2*n, 2*n)); h[:n, n:] = t; h[n:, :n] = t.T
    evals, vecs = np.linalg.eigh(h)
    minimum = float(np.min(abs(evals)))
    threshold = 1e-12*max(1., float(np.max(abs(jeff))))
    diagnostics = {"method": "real-space occupied-band perturbation",
                   "base_matter_energy_per_cell": finite_matter_energy(lx, ly, shear, jeff)/n,
                   "minimum_abs_f": minimum,
                   "minimum_majorana_positive_energy": 2*minimum,
                   "zero_mode_tolerance": threshold,
                   "zero_mode_count": int(np.sum(abs(evals) <= threshold))}
    if minimum <= threshold and g != 0.:
        scans = []
        for mu, phase, field in modes:
            scans.append({"channel": mu, "phase": phase,
                          "matter_energy_per_cell": [
                              {"amplitude": a, "energy": finite_matter_energy(
                                  lx, ly, shear, jeff, a*field, g)/n}
                              for a in (-step, -step/2, 0., step/2, step)]})
        return {"status": "zero_mode_uncertified", "matter_hessian": None,
                "total_hessian": None, "min_eigenvalue": None,
                "phonon_hessian": (spring_k*np.eye(len(modes))).tolist(),
                "mode_count": len(modes), "energy_scans": scans, **diagnostics}
    dt = [hopping_matrix(lx, ly, shear, [0., 0., 0.], m[2], g) for m in modes]
    occupied, empty = np.where(evals < -threshold)[0], np.where(evals > threshold)[0]
    amps = []
    for derivative in dt:
        dh = np.zeros_like(h); dh[:n, n:] = derivative; dh[n:, :n] = derivative.T
        amps.append(vecs[:, occupied].T@dh@vecs[:, empty])
    amps = np.asarray(amps)
    denom = 1./(evals[occupied, None]-evals[None, empty])
    matter = 2*np.einsum("iab,jab,ab->ij", amps, amps, denom)/n
    return _response(matter, labels, spring_k, diagnostics)


def bloch_energy(u, v, jeff):
    return -abs(jeff[2]+jeff[0]*np.exp(-2j*math.pi*u)+jeff[1]*np.exp(-2j*math.pi*v))


def momentum_reference(q, jeff, g, k_points, step=2e-4, spring_k=1.):
    """Both occupied-to-empty vertices for each k,k+q pair, in bounded batches.

    Eigenvectors of [[0,f],[f*,0]] are (-f/|f|,1)/sqrt(2) and
    (f/|f|,1)/sqrt(2). These give the analytic transition amplitudes below.
    Generic normalized cos/sin modes have Fourier coefficients 1/sqrt(2)
    and Hessian factor 2. At q=-q, cosine has coefficient 1 and the repeated
    k,k+q pair needs factor 1. Both vertices and all relative phases remain.
    The legacy step argument is unused because this formula is analytic.
    """
    q = _q_fraction(q)
    jeff = _parameters(jeff, g, spring_k)
    labels = _mode_labels(q)
    self_conjugate = len(labels) == 3
    mu = np.array([label[0] for label in labels])
    zplus = np.array([1. if self_conjugate else
                     (1. if phase == "cos" else 1j)/math.sqrt(2.)
                     for _, phase in labels])
    threshold = 1e-12*max(1., float(np.max(abs(jeff))))
    out = np.zeros((len(labels), len(labels)))
    minimum, npts = math.inf, 0
    points = iter(k_points)
    while True:
        chunk = list(itertools.islice(points, 8192))
        if not chunk:
            break
        uv = np.asarray(chunk, dtype=float)
        if uv.ndim != 2 or uv.shape[1] != 2 or not np.all(np.isfinite(uv)):
            raise ValueError("k points must be finite fractional coordinate pairs")
        phase0 = np.column_stack((np.exp(-2j*math.pi*uv), np.ones(len(uv))))
        phasep = np.column_stack((np.exp(-2j*math.pi*(uv+q)), np.ones(len(uv))))
        f0, fp = phase0@jeff, phasep@jeff
        e0, ep = abs(f0), abs(fp)
        minimum = min(minimum, float(np.min(e0)), float(np.min(ep)))
        if minimum <= threshold and g != 0.:
            bad = int(np.argmin(np.minimum(e0, ep)))
            return {"status": "zero_mode_uncertified", "matter_hessian": None,
                    "total_hessian": None, "min_eigenvalue": None,
                    "zero_mode_at_k": uv[bad].tolist(), "q": list(q),
                    "minimum_abs_f": minimum, "zero_mode_tolerance": threshold,
                    "k_points_evaluated_before_stop": npts+len(chunk)}
        npts += len(chunk)
        if g == 0.:
            continue
        p0, pp = f0/e0, fp/ep
        a = g*zplus[None, :]*phasep[:, mu]
        b = g*np.conj(zplus)[None, :]*phase0[:, mu]
        forward = .5*(-np.conj(p0[:, None])*a + np.conj(b)*pp[:, None])
        reverse = .5*(-np.conj(pp[:, None])*b + np.conj(a)*p0[:, None])
        amps = np.stack((forward, reverse), axis=-1)
        out -= (1. if self_conjugate else 2.)*np.real(np.einsum(
            "kia,kja,k->ij", np.conj(amps), amps, 1./(e0+ep)))
    if not npts:
        raise ValueError("reference k-point list must be nonempty")
    out /= npts
    return _response(out, labels, spring_k,
                     {"method": "Bloch occupied-band analytic response, two Fourier vertices",
                      "k_point_count": npts, "minimum_abs_f": minimum,
                      "zero_mode_tolerance": threshold})


def dense_reference(q, jeff, g, grid=200, step=2e-4, spring_k=1.):
    if not isinstance(grid, (int, np.integer)) or grid < 16 or grid % 2:
        raise ValueError("reference grid must be an even integer >=16")
    points = (np.arange(grid, dtype=float)+.5)/grid
    result = momentum_reference(q, jeff, g, itertools.product(points, repeat=2),
                                step, spring_k)
    result.update(grid=grid, integration_rule="fractional Brillouin-zone midpoint grid")
    return result


def analyze_geometry(lx, ly, shear, qbare, quniform, jeff, xbar, g, ref_grid,
                     step, spring_k=1., reference_cache=None, comparison_qs=()):
    grids = sorted(set([ref_grid] if isinstance(ref_grid, int) else ref_grid))
    if not grids:
        raise ValueError("at least one reference grid required")
    cache = {} if reference_cache is None else reference_cache
    grid = allowed_grid(lx, ly, shear)
    rows = []
    for p in grid:
        if p["index"] == 0:
            continue
        q = p["fraction"]
        h = static_hessian(lx, ly, shear, jeff, g, q, step, spring_k)
        row = {"q": list(q), "n": [p["n1"], p["n2"]], "index": p["index"],
               "distance_to_bare": momentum_distance(qbare, q),
               "distance_to_uniform": momentum_distance(quniform, q),
               "finite": h, "reference": None, "finite_reference_difference": None}
        for name, target in (("bare", qbare), ("uniform", quniform)):
            row["distance_to_"+name+"_pair"] = min(
                (momentum_distance(target, q), momentum_distance((-target[0], -target[1]), q)),
                key=lambda distance: distance["cartesian_distance"])
        finite_k = momentum_reference(q, jeff, g, (p["fraction"] for p in grid),
                                      spring_k=spring_k)
        if h["status"] == finite_k["status"] == "certified":
            row["finite_k_crosscheck_max_abs"] = float(np.max(abs(
                np.array(h["total_hessian"])-np.array(finite_k["total_hessian"]))))
        else:
            row["finite_k_crosscheck_max_abs"] = None
        rows.append(row)
    certified = [r for r in rows if r["finite"]["status"] == "certified"]
    soft_set = []
    if certified:
        best = min(r["finite"]["min_eigenvalue"] for r in certified)
        soft_set = [r for r in certified if r["finite"]["min_eigenvalue"] <= best+1e-10]
    def nearest(key):
        best = min(r[key]["cartesian_distance"] for r in rows)
        return [r for r in rows if r[key]["cartesian_distance"] <= best+1e-10]
    near_b_set, near_u_set = nearest("distance_to_bare"), nearest("distance_to_uniform")
    comparison_indices = {grid_point_for(q, lx, ly, shear)["index"] for q in comparison_qs}
    selected = {r["index"] for r in near_b_set+near_u_set+soft_set} | comparison_indices
    for row in rows:
        if row["index"] not in selected:
            continue
        row["common_q_geometry_comparison"] = row["index"] in comparison_indices
        row["references"] = {}
        for rg in grids:
            key = (tuple(row["q"]), rg, tuple(jeff), g, spring_k)
            if key not in cache:
                cache[key] = dense_reference(row["q"], jeff, g, rg, spring_k=spring_k)
            row["references"][str(rg)] = cache[key]
        row["reference"] = row["references"][str(grids[-1])]
        ref = row["reference"]
        if row["finite"]["status"] == ref["status"] == "certified":
            row["finite_reference_difference"] = row["finite"]["min_eigenvalue"]-ref["min_eigenvalue"]
        row["reference_grid_changes"] = []
        for lo, hi in zip(grids, grids[1:]):
            a, b = row["references"][str(lo)], row["references"][str(hi)]
            change = {"coarse_grid": lo, "fine_grid": hi, "status": "uncertified"}
            if a["status"] == b["status"] == "certified":
                change.update(status="certified",
                              min_eigenvalue_change=b["min_eigenvalue"]-a["min_eigenvalue"],
                              max_abs_matrix_change=float(np.max(abs(
                                  np.array(b["total_hessian"])-np.array(a["total_hessian"])))))
            row["reference_grid_changes"].append(change)
        analytic = row["finite"]["total_hessian"]
        row["finite_difference"] = finite_difference_hessians(
            lx, ly, shear, jeff, xbar, g, row["q"], (step, step/2, step/4),
            spring_k, np.asarray(analytic) if analytic is not None else None)
        row["finite_difference"]["status"] = (
            "analytic_comparison" if analytic is not None else "energy_differences_only_uncertified")
    base = finite_matter_energy(lx, ly, shear, jeff)
    phonon = .5*spring_k*lx*ly*float(np.dot(xbar, xbar))
    bloch = sum(bloch_energy(*p["fraction"], jeff) for p in grid)
    a1, a2 = np.array([math.sqrt(3.)/2, 1.5]), np.array([-math.sqrt(3.)/2, 1.5])
    return {"Lx": lx, "Ly": ly, "shear": shear, "Ncell": lx*ly,
            "period_vectors_cartesian": [(lx*a1).tolist(), (shear*a1+ly*a2).tolist()],
            "wilson_sector": [1, 1], "flux_sector": "zero flux, uniform A-to-B gauge",
            "physical_spin_parity_projection": False,
            "common_Xbar": list(xbar), "common_Jeff": list(jeff),
            "static_energy": {"matter": base, "phonon": phonon, "total": base+phonon,
                              "matter_per_cell": base/(lx*ly), "phonon_per_cell": phonon/(lx*ly),
                              "total_per_cell": (base+phonon)/(lx*ly),
                              "real_space_minus_bloch": base-bloch},
            "all_q": rows, "nearest_bare": near_b_set[0],
            "nearest_bare_candidates": near_b_set, "nearest_uniform": near_u_set[0],
            "nearest_uniform_candidates": near_u_set,
            "softest": soft_set[0] if soft_set else None, "softest_candidates": soft_set,
            "soft_search_complete": len(certified) == len(rows),
            "uncertified_indices": [r["index"] for r in rows if r["finite"]["status"] != "certified"]}


def _plot(results, path, target_q=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(results), figsize=(13, 4), constrained_layout=True)
    colors = [r["finite"]["min_eigenvalue"] for geo in results for r in geo["all_q"]
              if r["finite"]["status"] == "certified"]
    limits = (min(colors), max(colors)) if colors else (0., 1.)
    for ax, geo in zip(np.atleast_1d(axes), results):
        rows = [r for r in geo["all_q"] if r["finite"]["status"] == "certified"]
        if rows:
            xy = np.array([reciprocal_cartesian(*r["q"]) for r in rows])
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=[r["finite"]["min_eigenvalue"] for r in rows],
                            cmap="viridis", s=65, vmin=limits[0], vmax=limits[1])
            fig.colorbar(sc, ax=ax, label="min static curvature")
        for row in geo["softest_candidates"]:
            x, y = reciprocal_cartesian(*row["q"])
            ax.scatter([x], [y], s=130, facecolors="none", edgecolors="red")
        if target_q is not None:
            targets = [reciprocal_cartesian(*_q_fraction(q)) for q in
                       (target_q, (-target_q[0], -target_q[1]))]
            ax.scatter(*np.array(targets).T, marker="x", color="black", s=90, label="uniform Q pair")
            ax.legend(fontsize=8)
        ax.set_aspect("equal", adjustable="box")
        ax.set(title=f"({geo['Lx']},{geo['Ly']},{geo['shear']})", xlabel="q_x", ylabel="q_y")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _write_candidates(results, output_dir):
    """Allowed measurement momenta, with conjugates; never change the model."""
    for geo in results:
        grid = allowed_grid(geo["Lx"], geo["Ly"], geo["shear"])
        selected = {r["index"] for r in geo["all_q"] if r.get("references")}
        for index in list(selected):
            q = grid[index]["fraction"]
            selected.add(grid_point_for((-q[0], -q[1]), geo["Lx"], geo["Ly"], geo["shear"])["index"])
        name = f"q_candidates_Lx{geo['Lx']}_Ly{geo['Ly']}_s{geo['shear']}.txt"
        lines = ["# Allowed static-response candidates and conjugates, label u v.",
                 "# Measurement input for --vbs-q-file with selected/both mode.",
                 "# Static candidates do not establish a QMC peak or a quantum critical point."]
        for index in sorted(selected):
            point = grid[index]
            lines.append(f"grid_{point['n1']}_{point['n2']} {point['fraction'][0]:.17g} {point['fraction'][1]:.17g}")
        (output_dir/name).write_text("\n".join(lines)+"\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--Jx", type=float, default=1.)
    ap.add_argument("--Jy", type=float, default=1.)
    ap.add_argument("--Jz", type=float, default=.5)
    ap.add_argument("--g", type=float, default=.3)
    ap.add_argument("--mass", type=float, default=1.7)
    ap.add_argument("--omega-ph", type=float, default=.8)
    ap.add_argument("--spring-k", type=float, help="optional consistency check: k=mass*omega_ph^2")
    ap.add_argument("--geometries", default="6,6,0;12,3,0;12,3,3")
    ap.add_argument("--uniform-grid", type=int, default=800)
    ap.add_argument("--reference-grids", default="200,400")
    ap.add_argument("--finite-step", type=float, default=1e-3,
                    help="largest energy-difference amplitude; also check step/2 and step/4")
    ap.add_argument("--plot", action="store_true", help="export PNG curvature maps")
    args = ap.parse_args()
    start = time.perf_counter()
    try:
        if any(not math.isfinite(v) or v <= 0 for v in (args.mass, args.omega_ph, args.finite_step)):
            raise ValueError("mass, omega-ph and finite-step must be finite and positive")
        k = args.mass*args.omega_ph**2
        J = _parameters([args.Jx, args.Jy, args.Jz], args.g, k)
        if args.spring_k is not None and not math.isclose(args.spring_k, k, rel_tol=1e-12):
            raise ValueError("spring-k must equal mass*omega-ph^2")
        grids = sorted(set(int(x) for x in args.reference_grids.split(",")))
        if not grids or any(n < 16 or n % 2 for n in grids):
            raise ValueError("reference grids must be even integers >=16")
        if args.uniform_grid < 32 or args.uniform_grid % 2:
            raise ValueError("uniform grid must be an even integer >=32")
        geometries = [tuple(int(x) for x in s.split(",")) for s in args.geometries.split(";")]
        if any(len(geo) != 3 for geo in geometries):
            raise ValueError("geometries require Lx,Ly,shear triples")
        for geo in geometries:
            _geometry(*geo)
        sols = {str(n): uniform_background_solution(
            J, args.g, k, tolerance=1e-13, max_iterations=300, mixing=.5, integration_grid=n)
            for n in sorted(set([args.uniform_grid//2, args.uniform_grid]))}
        if any(not s["converged"] for s in sols.values()):
            raise RuntimeError("uniform background iteration did not converge")
        sol = sols[str(args.uniform_grid)]
        bare = q_predictions(dirac_nodes(*J)["nodes"])
        eff = q_predictions(dirac_nodes(*sol["Jeff"])["nodes"])
        if not bare or not eff:
            raise ValueError("comparison requires isolated bare and uniform Dirac nodes")
        qbare, quniform = (bare[0]["u"], bare[0]["v"]), (eff[0]["u"], eff[0]["v"])
        common = [p["fraction"] for p in allowed_grid(*geometries[0])
                  if p["index"] != 0 and all(grid_point_for(p["fraction"], *geo) is not None
                                             for geo in geometries[1:])]
        comparison_qs = []
        if common:
            distance = min(momentum_distance(quniform, q)["cartesian_distance"] for q in common)
            comparison_qs = [q for q in common
                             if momentum_distance(quniform, q)["cartesian_distance"] <= distance+1e-10]
        cache = {}
        results = [analyze_geometry(*geo, qbare, quniform, sol["Jeff"], sol["Xbar"],
                                    args.g, grids, args.finite_step, k, cache, comparison_qs)
                   for geo in geometries]
    except (ValueError, RuntimeError) as exc:
        ap.error(str(exc))
    output = {"schema": "static-vbs-curvature-v2", "argv": sys.argv,
              "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (Path(__file__), HERE/"predict_vbs_momenta.py")},
              "parameters": {"J": J.tolist(), "g": args.g, "mass": args.mass,
                             "omega_ph": args.omega_ph, "spring_k": k, "lambda": args.g**2/k},
              "uniform_convergence": sols, "common_background_grid": args.uniform_grid,
              "Q_bare": list(qbare), "Q_uniform_effective": list(quniform),
              "reference_grids": grids, "geometries": results,
              "common_q_geometry_comparisons": [
                  {"q": list(q), "geometries": [
                      {"shape": [geo["Lx"], geo["Ly"], geo["shear"]],
                       "finite_min_curvature": row["finite"]["min_eigenvalue"],
                       "reference_min_curvature": row["reference"]["min_eigenvalue"],
                       "finite_minus_reference": row["finite_reference_difference"]}
                      for geo in results for row in geo["all_q"]
                      if grid_point_for(q, geo["Lx"], geo["Ly"], geo["shear"])["index"] == row["index"]]}
                  for q in comparison_qs],
              "elapsed_seconds": time.perf_counter()-start,
              "scope": "fixed uniform static matter+elastic Hessian; no phonon quantum fluctuations, QMC peak or quantum critical point",
              "reference_scope": "same-q BZ quadrature only; no global continuous-q minimum search"}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"static_response.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    _write_candidates(results, args.output_dir)
    with (args.output_dir/"static_response.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Lx", "Ly", "shear", "q", "n", "finite_status", "finite_min_curvature",
                    "reference_min_curvature", "finite_minus_reference", "reference_grid_change",
                    "distance_to_bare", "distance_to_uniform", "finite_k_matrix_error"])
        for geo in results:
            for row in geo["all_q"]:
                ref = row["reference"] or {}
                changes = row.get("reference_grid_changes", [])
                w.writerow([geo["Lx"], geo["Ly"], geo["shear"], row["q"], row["n"],
                            row["finite"]["status"], row["finite"]["min_eigenvalue"],
                            ref.get("min_eigenvalue"), row["finite_reference_difference"],
                            changes[-1].get("min_eigenvalue_change") if changes else None,
                            row["distance_to_bare"]["cartesian_distance"],
                            row["distance_to_uniform"]["cartesian_distance"],
                            row["finite_k_crosscheck_max_abs"]])
    if args.plot:
        _plot(results, args.output_dir/"static_response.png", quniform)
    print(json.dumps({"parameters": output["parameters"], "Q_bare": qbare,
                      "Q_uniform": quniform, "elapsed_seconds": output["elapsed_seconds"],
                      "geometries": [{"shape": [r["Lx"], r["Ly"], r["shear"]],
                                      "softest_q": r["softest"]["q"] if r["softest"] else None,
                                      "min_curvature": r["softest"]["finite"]["min_eigenvalue"] if r["softest"] else None,
                                      "soft_search_complete": r["soft_search_complete"]}
                                     for r in results]}, indent=2))


if __name__ == "__main__":
    main()
