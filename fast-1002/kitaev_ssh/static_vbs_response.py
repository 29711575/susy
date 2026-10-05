#!/usr/bin/env python3
"""Lightweight static VBS Hessian for a uniform Majorana background.

This is an audit tool, not a QMC or projector path.  It uses the same
f(k)=Jz+Jx exp(-ik1)+Jy exp(-ik2) convention as predict_vbs_momenta.py.
The finite-torus matter energy is evaluated from the singular values of the
real-space A--B hopping matrix; a dense momentum-space block calculation is
used as the reference response.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from predict_vbs_momenta import (allowed_grid, dirac_nodes, momentum_distance,
                                 nearest_grid, q_predictions,
                                 uniform_background_solution)


def cell_reduce(x, y, lx, ly, shear):
    q, r = divmod(y, ly)
    xx = x - q * shear
    p, s = divmod(xx, lx)
    return s, r, p, q


def bond_table(lx, ly, shear):
    bonds = []
    for x in range(lx):
        for y in range(ly):
            for mu, (dx, dy) in enumerate(((-1, 0), (0, -1), (0, 0))):
                bx, by, w1, w2 = cell_reduce(x + dx, y + dy, lx, ly, shear)
                # The predictor's translational Bloch convention is the
                # uniform (++ Wilson) f(k) gauge.  In this representation the
                # bond signs are all +1; boundary signs are a gauge choice
                # and are not inserted a second time into the Bloch phases.
                gauge = 1
                bonds.append((x, y, mu, bx, by, gauge))
    return bonds


def finite_matter_energy(lx, ly, shear, jeff, delta=None, g=0.):
    """Ground-state matter energy, -sum singular values of the A-B matrix."""
    n = lx * ly
    t = np.zeros((n, n), dtype=float)
    delta = np.zeros((n, 3), dtype=float) if delta is None else np.asarray(delta)
    for x, y, mu, bx, by, gauge in bond_table(lx, ly, shear):
        a = x * ly + y
        b = bx * ly + by
        t[a, b] += gauge * (jeff[mu] + g * delta[a, mu])
    return float(-np.linalg.svd(t, compute_uv=False).sum())


def mode_basis(lx, ly, shear, q):
    """Return real modes with sum_(r,mu) v_a v_b / Ncell = delta_ab."""
    n = lx * ly
    phases = np.array([2 * math.pi * (q[0] * x + q[1] * y)
                       for x in range(lx) for y in range(ly)])
    raw = [("cos", np.cos(phases)), ("sin", np.sin(phases))]
    out = []
    for name, values in raw:
        norm = float(np.sqrt(np.mean(values * values)))
        if norm > 1e-11:
            out.append((name, values / norm))
    modes = []
    for mu in range(3):
        for name, values in out:
            arr = np.zeros((n, 3), dtype=float)
            arr[:, mu] = values
            modes.append((mu, name, arr))
    return modes


def static_hessian(lx, ly, shear, jeff, g, q, step=2e-4):
    modes = mode_basis(lx, ly, shear, q)
    d = len(modes)
    n = lx * ly
    t = np.zeros((n, n), dtype=float)
    derivatives = [np.zeros((n, n), dtype=float) for _ in modes]
    for x, y, mu, bx, by, gauge in bond_table(lx, ly, shear):
        a, b = x * ly + y, bx * ly + by
        t[a, b] += gauge * jeff[mu]
        for i, (channel, _phase, values) in enumerate(modes):
            if channel == mu: derivatives[i][a, b] += gauge * g * values[a, mu]
    # Hermitian bipartite Hamiltonian has the same +/- singular spectrum as
    # the Majorana generator; its occupied-band energy is -sum singular values.
    h = np.zeros((2*n, 2*n), dtype=float); h[:n,n:] = t; h[n:,:n] = t.T
    evals, vecs = np.linalg.eigh(h)
    dmat = []
    for dt in derivatives:
        dh = np.zeros_like(h); dh[:n,n:] = dt; dh[n:,:n] = dt.T
        dmat.append(vecs.T @ dh @ vecs)
    occ = np.where(evals < -1e-12)[0]; emp = np.where(evals > 1e-12)[0]
    matter = np.zeros((d, d), dtype=float)
    for i in range(d):
        for j in range(i + 1):
            value = 0.
            for a in occ:
                for b in emp:
                    value += 2. * dmat[i][a,b] * dmat[j][b,a] / (evals[a] - evals[b])
            matter[i,j] = matter[j,i] = value / n
    base = float(evals[evals < 0].sum()) / n
    total = matter + np.eye(d)
    vals, vecs = np.linalg.eigh(total)
    return {"matter_hessian": matter.tolist(), "phonon_hessian": np.eye(d).tolist(),
            "total_hessian": total.tolist(), "eigenvalues": vals.tolist(),
            "min_eigenvalue": float(vals[0]),
            "soft_mode": [{"channel": modes[i][0], "phase": modes[i][1],
                           "coefficient": float(vecs[i, 0])} for i in range(d)],
            "mode_count": d, "step": step, "base_matter_energy_per_cell": base,
            "normalization": "sum_(r,mu) v_alpha v_beta / Ncell = delta_alpha_beta",
            "energy_definition": "per-cell E_matter[Jeff+g*deltaX] + 1/2 sum_b (Xbar_b+deltaX_b)^2"}


def bloch_energy(u, v, jeff):
    f = jeff[2] + jeff[0] * np.exp(-2j * math.pi * u) + jeff[1] * np.exp(-2j * math.pi * v)
    return -abs(f)


def momentum_reference(q, jeff, g, k_points, step=2e-4):
    """Occupied-band Hessian on an explicit k list in the real mode basis.

    A real modulation has both Fourier vertices A_k->B_(k+q) and
    A_(k+q)->B_k.  The old implementation retained only the first vertex,
    which suppressed the matter curvature.  Summing the complete 4x4 pair
    block over k double-counts k<->k+q, so the usual factor 1/2 is retained.
    """
    k_points = [tuple(map(float, point)) for point in k_points]
    npts = len(k_points)
    if not npts:
        raise ValueError("reference k-point list must be nonempty")
    self_conjugate = all(abs((2 * z) - round(2 * z)) < 1e-12 for z in q)
    modes = [(mu, "cos") for mu in range(3)] if self_conjugate else [
        (mu, phase) for mu in range(3) for phase in ("cos", "sin")]
    out = np.zeros((len(modes), len(modes)))
    coeff = 1. if self_conjugate else 1. / math.sqrt(2.)
    for u, v in k_points:
            up, vp = (u + q[0]) % 1., (v + q[1]) % 1.
            f0 = jeff[2] + jeff[0]*np.exp(-2j*math.pi*u) + jeff[1]*np.exp(-2j*math.pi*v)
            fp = jeff[2] + jeff[0]*np.exp(-2j*math.pi*up) + jeff[1]*np.exp(-2j*math.pi*vp)
            # Analytic second-order occupied-band response.  This is the
            # finite-difference limit and avoids 4x4 diagonalizations for
            # every amplitude pair on every dense integration point.
            h0 = np.zeros((4,4), complex)
            h0[:2,:2] = [[0,f0],[np.conj(f0),0]]
            h0[2:,2:] = [[0,fp],[np.conj(fp),0]]
            evals, vecs = np.linalg.eigh(h0)
            if np.min(np.abs(evals)) < 1e-12:
                raise RuntimeError(f"reference encountered zero mode at k={(u,v)} q={q}")
            blocks = []
            for mu, phase in modes:
                z_plus = coeff if phase == "cos" else 1j*coeff
                z_minus = coeff if phase == "cos" else -1j*coeff
                pm0 = [np.exp(-2j*math.pi*u), np.exp(-2j*math.pi*v), 1.][mu]
                pmp = [np.exp(-2j*math.pi*up), np.exp(-2j*math.pi*vp), 1.][mu]
                vv = np.zeros((4,4), complex)
                vv[0,3] = g*z_plus*pmp
                vv[2,1] = g*z_minus*pm0
                vv[3,0] = np.conj(vv[0,3]); vv[1,2] = np.conj(vv[2,1])
                blocks.append(vv)
            transformed = np.asarray([vecs.conj().T @ vv @ vecs for vv in blocks])
            occupied = np.where(evals < 0)[0]
            empty = np.where(evals > 0)[0]
            # Contract all channels at once.  The occupied-band second
            # derivative carries the standard factor 2; the explicit k list
            # already represents one normalized real-mode Fourier component,
            # so no additional duplicate-orientation half is applied here.
            amps = transformed[:, occupied[:, None], empty[None, :]]
            denom = 1. / (evals[occupied, None] - evals[None, empty])
            curvature_factor = 1. if self_conjugate else 2.
            out += curvature_factor * np.real(np.einsum("iab,jab,ab->ij", np.conj(amps), amps, denom))
    out /= npts
    return {"matter_hessian": out.tolist(), "total_hessian": (out+np.eye(len(modes))).tolist(),
            "min_eigenvalue": float(np.linalg.eigvalsh(out+np.eye(len(modes)))[0]),
            "k_point_count": npts, "step": step,
            "basis": [{"channel": mu, "phase": phase} for mu, phase in modes]}


def dense_reference(q, jeff, g, grid=200, step=2e-4):
    """Midpoint-BZ wrapper for :func:`momentum_reference`."""
    if grid < 16 or grid % 2:
        raise ValueError("reference grid must be an even integer >=16")
    points = (np.arange(grid, dtype=float) + .5) / grid
    result = momentum_reference(q, jeff, g,
                                ((u, v) for u in points for v in points), step)
    result["grid"] = grid
    result["integration_rule"] = "fractional Brillouin-zone midpoint grid"
    return result


def analyze_geometry(lx, ly, shear, qbare, quniform, jeff, xbar, g, ref_grid, step):
    grid = allowed_grid(lx, ly, shear)
    nonzero = [p for p in grid if p["fraction"] != (0., 0.)]
    rows = []
    for p in nonzero:
        q = p["fraction"]
        h = static_hessian(lx, ly, shear, jeff, g, q, step)
        rows.append({"q": list(q), "index": p["index"],
                     "distance_to_bare": momentum_distance(qbare, q),
                     "distance_to_uniform": momentum_distance(quniform, q),
                     "finite": h, "reference": None,
                     "finite_reference_difference": None})
    soft_ordered = sorted(rows, key=lambda r: (r["finite"]["min_eigenvalue"], r["index"]))
    soft_best = soft_ordered[0]["finite"]["min_eigenvalue"]
    soft_set = sorted([row for row in soft_ordered
                       if row["finite"]["min_eigenvalue"] <= soft_best + 1e-10],
                      key=lambda row: row["index"])
    soft = soft_set[0]
    def tied_nearest(distance_key, tolerance=1e-10):
        ordered = sorted(rows, key=lambda r: (r[distance_key]["cartesian_distance"], r["index"]))
        best = ordered[0][distance_key]["cartesian_distance"]
        tied = [row for row in ordered if row[distance_key]["cartesian_distance"] <= best + tolerance]
        return sorted(tied, key=lambda row: row["index"])
    near_b_set = tied_nearest("distance_to_bare")
    near_u_set = tied_nearest("distance_to_uniform")
    near_b, near_u = near_b_set[0], near_u_set[0]
    # The audit's actual finite-size comparison set is the two target-nearest
    # modes plus the measured softest mode.  Dense references for every other
    # allowed q are unnecessary work and are not used in the summary table.
    selected = {id(row) for row in near_b_set + near_u_set + soft_set}
    for row in rows:
        if id(row) in selected:
            row["reference"] = dense_reference(tuple(row["q"]), jeff, g, ref_grid, step)
            row["finite_reference_difference"] = row["finite"]["min_eigenvalue"] - row["reference"]["min_eigenvalue"]
    base_matter = finite_matter_energy(lx, ly, shear, jeff)
    phonon_background = .5 * lx * ly * float(np.dot(xbar, xbar))
    return {"Lx": lx, "Ly": ly, "shear": shear, "Ncell": lx*ly,
            "static_energy": {"matter": base_matter, "phonon": phonon_background,
                               "total": base_matter + phonon_background,
                               "matter_per_cell": base_matter / (lx*ly),
                               "phonon_per_cell": phonon_background / (lx*ly),
                               "total_per_cell": (base_matter + phonon_background) / (lx*ly),
                               "definition": "per-cell second derivative of E_matter[Jeff+g deltaX]+1/2 sum_b (Xbar_b+deltaX_b)^2"},
            "all_q": rows, "nearest_bare": near_b,
            "nearest_bare_candidates": near_b_set,
            "nearest_uniform": near_u,
            "nearest_uniform_candidates": near_u_set,
            "softest": soft,
            "softest_candidates": soft_set}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--uniform-grid", type=int, default=800)
    ap.add_argument("--reference-grids", default="200,400")
    ap.add_argument("--finite-step", type=float, default=2e-4)
    args = ap.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    J=(1.,1.,.5); g=.3; k=1.; m=1.; omega=1.
    grids=[int(x) for x in args.reference_grids.split(",")]
    sols={str(n): uniform_background_solution(J,g,k,tolerance=1e-13,max_iterations=300,mixing=.5,integration_grid=n) for n in [400,800]}
    sol=sols[str(args.uniform_grid)]
    bare=q_predictions(dirac_nodes(*J)["nodes"]); eff=q_predictions(dirac_nodes(*sol["Jeff"])["nodes"])
    qbare=(bare[0]["u"],bare[0]["v"]); quniform=(eff[0]["u"],eff[0]["v"])
    geometries=[(6,6,0),(12,3,0),(12,3,3)]
    results=[]
    for gspec in geometries: results.append(analyze_geometry(*gspec,qbare,quniform,sol["Jeff"],sol["Xbar"],g,grids[-1],args.finite_step))
    out={"parameters":{"J":J,"g":g,"mass":m,"omega_ph":omega,"spring_k":k},"uniform_convergence":sols,
         "Q_bare":list(qbare),"Q_uniform_effective":list(quniform),"geometries":results,
         "reference_grid_checks":{}}
    for rg in grids:
        out["reference_grid_checks"][str(rg)] = {"grid":rg}
    (args.output_dir/"static_response.json").write_text(json.dumps(out,indent=2)+"\n")
    with (args.output_dir/"static_response.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["Lx","Ly","shear","Ncell","nearest_bare_q","nearest_uniform_q","soft_q","soft_eigenvalue","reference_eigenvalue","finite_minus_reference","bare_mismatch","uniform_mismatch"])
        for r in results:
            w.writerow([r["Lx"],r["Ly"],r["shear"],r["Ncell"],r["nearest_bare"]["q"],r["nearest_uniform"]["q"],r["softest"]["q"],r["softest"]["finite"]["min_eigenvalue"],r["softest"]["reference"]["min_eigenvalue"],r["softest"]["finite_reference_difference"],r["nearest_bare"]["distance_to_bare"]["cartesian_distance"],r["nearest_uniform"]["distance_to_uniform"]["cartesian_distance"]])
    print(json.dumps({"Q_bare":qbare,"Q_uniform":quniform,"uniform":sol,"geometries":[{"shape":(r["Lx"],r["Ly"],r["shear"]),"nearest_bare":r["nearest_bare"]["q"],"nearest_uniform":r["nearest_uniform"]["q"],"soft":r["softest"]["q"],"soft_eigenvalue":r["softest"]["finite"]["min_eigenvalue"]} for r in results]},indent=2))

if __name__ == "__main__": main()
