#!/usr/bin/env python3
"""Export grid and arbitrary-path VBS spectra from one completed coupled_mc run.

Uses the existing full-grid block sums, not new QMC measurements. Off-grid
matrices are spatial-translation-averaged finite-window estimators (handoff
T1/C1/C2), not fixed-origin Ssel. Requires a translation-invariant ensemble
for equality of their physical expectations. No Hamiltonian or RNG changes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from predict_vbs_momenta import dirac_nodes, q_predictions
from vbs_window import connected_window_association_ratios

AXES = "xyz"
ESTIMATORS = ("raw", "connected", "projected")
CHANNELS = ("szz", "trace", "lambda_max")
TOL = 1e-10


def reciprocal_distance(a, b, reciprocal):
    """Physical reciprocal distance modulo reciprocal-lattice periods."""
    d = np.asarray(a, float) - np.asarray(b, float)
    basis = np.asarray(reciprocal, float)
    candidates = [(float(np.linalg.norm((d + np.array([i, j])) @ basis)),
                   [float(d[0] + i), float(d[1] + j)])
                  for i in range(-2, 3) for j in range(-2, 3)]
    return min(candidates, key=lambda x: x[0])


def fixed_vbs_geometry(meta):
    """Return the verified finite-torus reciprocal steps for fixed windows."""
    nx, ny, shear = meta.get("Lx"), meta.get("Ly"), meta.get("shear")
    if any(type(x) is not int for x in (nx, ny, shear)) or nx < 2 or ny < 2:
        raise ValueError("fixed VBS windows require integer Lx/Ly/shear")
    if meta.get("spatial_boundary") != "pbc" or not meta.get("vbs_enabled"):
        raise ValueError("fixed VBS windows require periodic full-grid data")
    pbc_geometry = meta.get("vbs", {}).get("pbc_geometry", "axis_aligned")
    if pbc_geometry != "axis_aligned":
        raise ValueError(f"unsupported finite reciprocal geometry: {pbc_geometry}")
    if shear == 0:
        d1 = np.array([1.0 / nx, 0.0])
        d2 = np.array([0.0, 1.0 / ny])
        source = "ordinary_pbc_b1_over_Lx_b2_over_Ly"
    else:
        # These are the actual dual characters recorded by the sheared torus.
        d1 = np.array([1.0 / nx, -float(shear) / (nx * ny)])
        d2 = np.array([0.0, 1.0 / ny])
        source = "verified_finite_torus_reciprocal_basis"
        momenta = meta.get("vbs", {}).get("momenta", [])
        if momenta:
            recorded = finite_array([p.get("reciprocal_fraction") for p in momenta], (nx * ny, 2))
            if (np.max(abs(wrapped_difference(recorded[ny] - recorded[0], d1))) > TOL or
                    np.max(abs(wrapped_difference(recorded[1] - recorded[0], d2))) > TOL):
                raise ValueError("recorded finite reciprocal basis does not match shear geometry")
    return d1, d2, source


def fixed_vbs_groups(meta, explicit_q=None):
    """Build Q_D and -Q_D groups without snapping their coordinates."""
    d1, d2, step_source = fixed_vbs_geometry(meta)
    if explicit_q is None:
        couplings = finite_array([meta[k] for k in ("Jx", "Jy", "Jz")], (3,))
        nodes = dirac_nodes(*couplings)
        if nodes["status"] != "ok":
            raise ValueError(f"bare Dirac prediction unavailable: {nodes['status']}")
        predicted = q_predictions(nodes["nodes"])
        if not predicted:
            raise ValueError("bare Dirac predictor returned no Q_D")
        q = np.array([predicted[0]["u"], predicted[0]["v"]], float)
        source = "Q_D=K_plus-K_minus from bare Jx,Jy,Jz"
    else:
        q = finite_array(explicit_q, (2,)) % 1.0
        source = "explicit_preselected_Q"
        nodes = {"status": "explicit", "nodes": []}
    q = q % 1.0
    centers = [("Q_D", q), ("-Q_D", (-q) % 1.0)]
    groups = []
    for label, center in centers:
        points = {"center": center,
                  "plus_b1": (center + d1) % 1.0,
                  "minus_b1": (center - d1) % 1.0,
                  "plus_b2": (center + d2) % 1.0,
                  "minus_b2": (center - d2) % 1.0}
        groups.append({"center_label": label, "center": center,
                       "points": points, "step_b1": d1.copy(),
                       "step_b2": d2.copy(), "step_source": step_source,
                       "source": source, "prediction": nodes})
    return groups


def finite_array(values, shape=None):
    out = np.asarray(values, dtype=float)
    if (shape is not None and out.shape != shape) or not np.all(np.isfinite(out)):
        raise ValueError(f"expected finite numeric array of shape {shape}")
    return out


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def wrapped_difference(a, b):
    d = np.asarray(a) - np.asarray(b)
    return d - np.rint(d)


def jackknife_se(replicas):
    """Equal-block delete-one jackknife; never treats momenta as independent."""
    r = np.asarray(replicas)
    if len(r) < 2:
        raise ValueError("at least two equal blocks are required for an error estimate")
    return np.sqrt((len(r) - 1) / len(r) * np.sum((r - r.mean(axis=0)) ** 2, axis=0))


def read_run(results_path, blocks_path=None, reblock=1):
    """Read actual PhaseStatistics CSV field names; reject partial/mismatched runs."""
    results_path = Path(results_path)
    meta = json.loads(results_path.read_text(encoding="utf-8"))
    if meta.get("status") != "complete":
        raise ValueError("results status must be complete; partial/aborted runs are not accepted")
    vbs = meta.get("vbs", {})
    if not meta.get("vbs_enabled") or vbs.get("mode") not in ("all", "both"):
        raise ValueError("full-grid blocks required: run with --vbs 1 --vbs-mode all or both")
    nx, ny, shear = (meta.get(k) for k in ("Lx", "Ly", "shear"))
    if any(type(x) is not int for x in (nx, ny, shear)) or min(nx, ny) < 2 or not 0 <= shear < nx:
        raise ValueError("invalid or missing Lx/Ly/shear")
    nc = nx * ny
    if vbs.get("Nc") != nc or meta.get("spatial_boundary") != "pbc":
        raise ValueError("inconsistent cell count or unsupported boundary")
    reciprocal = finite_array(vbs.get("reciprocal"), (2, 2))
    if abs(np.linalg.det(reciprocal)) < 1e-14:
        raise ValueError("singular reciprocal basis")
    momenta = vbs.get("momenta", [])
    if len(momenta) != nc:
        raise ValueError("metadata does not contain the complete grid")
    grid = finite_array([p.get("reciprocal_fraction") for p in momenta], (nc, 2))
    for i, p in enumerate(momenta):
        n1, n2 = divmod(i, ny)
        expected = [n1 / nx, n2 / ny - shear * n1 / (nx * ny)]
        if p.get("index") != i or p.get("n") != [n1, n2] or np.max(abs(wrapped_difference(grid[i], expected))) > TOL:
            raise ValueError("grid order/fractions do not match the recorded torus")
        if np.max(abs(finite_array(p.get("cartesian"), (2,)) - grid[i] @ reciprocal)) > 1e-8:
            raise ValueError("Cartesian momentum metadata disagrees with reciprocal fractions")
    if blocks_path is None:
        if not results_path.name.endswith("_results.json"):
            raise ValueError("use --blocks for a results filename not ending in _results.json")
        blocks_path = results_path.with_name(results_path.name[:-13] + "_blocks.csv")
    blocks_path = Path(blocks_path)
    snames = [f"S_q{k}_{mu}{nu}" for k in range(nc) for mu in AXES for nu in AXES]
    dnames = [f"Delta_q0_{mu}" for mu in AXES]
    phases, counts, numerators, delta0 = [], [], [], []
    with blocks_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"block", "count", "sum_phase_real", "sum_phase_imag"}
        required.update(f"sum_phase_{name}_{part}" for name in snames + dnames for part in ("real", "imag"))
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not required.issubset(fields):
            raise ValueError("missing/duplicate full-grid CSV columns (selected-only data are insufficient)")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("truncated or malformed CSV row")
            if int(row["block"]) != len(counts) or int(row["count"]) <= 0:
                raise ValueError("blocks must be sequential, starting at zero, with positive counts")
            counts.append(int(row["count"]))
            phases.append(complex(float(row["sum_phase_real"]), float(row["sum_phase_imag"])))
            def summed(name):
                return complex(float(row[f"sum_phase_{name}_real"]), float(row[f"sum_phase_{name}_imag"]))
            numerators.append([summed(name) for name in snames])
            delta0.append([summed(name) for name in dnames])
    if not counts or len(set(counts)) != 1:
        raise ValueError("equal, nonempty blocks required; unequal blocks are not silently pooled")
    if len(counts) != meta.get("blocks") or sum(counts) != meta.get("samples"):
        raise ValueError("CSV counts do not match the completed results JSON")
    if counts[0] != meta.get("block_samples"):
        raise ValueError("CSV block size differs from results metadata")
    phase = np.asarray(phases, complex)
    snum = np.asarray(numerators, complex).reshape(len(counts), nc, 3, 3)
    dnum = np.asarray(delta0, complex)
    if not all(np.all(np.isfinite(a)) for a in (phase, snum, dnum)):
        raise ValueError("non-finite block sums")
    if type(reblock) is not int or reblock < 1 or len(counts) % reblock:
        raise ValueError("--reblock must divide the block count exactly; no remainder is discarded")
    def merge(a):
        return a.reshape(len(counts) // reblock, reblock, *a.shape[1:]).sum(axis=1)
    phase, snum, dnum = (merge(a) for a in (phase, snum, dnum))
    if len(phase) < 2:
        raise ValueError("at least two blocks must remain after reblocking")
    return dict(meta=meta, grid=grid, reciprocal=reciprocal, phase=phase,
                snum=snum, dnum=dnum, block_count=counts[0] * reblock,
                results_path=str(results_path.resolve()), blocks_path=str(blocks_path.resolve()))


def ratio_replicas(run):
    """Index 0 is the full-sample plug-in; remaining indices delete one block."""
    p, s, d = (run[k] for k in ("phase", "snum", "dnum"))
    n = run["block_count"]
    denominators = np.r_[p.sum(), p.sum() - p]
    samples = np.r_[len(p) * n, np.full(len(p), (len(p) - 1) * n)]
    if np.any(abs(denominators) <= 1e-12 * samples):
        raise ValueError("unresolved full/delete-block phase denominator; no division or clipping performed")
    srat = np.concatenate([s.sum(axis=0)[None], s.sum(axis=0)[None] - s], axis=0) / denominators[:, None, None, None]
    drat = np.concatenate([d.sum(axis=0)[None], d.sum(axis=0)[None] - d], axis=0) / denominators[:, None]
    phase_ratios = denominators / samples
    phase_se = np.array([jackknife_se(phase_ratios[1:].real), jackknife_se(phase_ratios[1:].imag)])
    diagnostic = dict(average_phase=[float(phase_ratios[0].real), float(phase_ratios[0].imag)],
                      average_phase_se=phase_se.tolist(),
                      phase_resolved_5sigma=bool(abs(phase_ratios[0]) > 5 * np.linalg.norm(phase_se)),
                      blocks=len(p), block_samples=n, block_independence_verified=False,
                      too_few_blocks_diagnostic=bool(len(p) < 16))
    return srat, drat, diagnostic


def window_weights(q, grid, nx, ny):
    """Exact rectangular-representative window, also valid for sheared grids.

    Uses stable Dirichlet kernels. A grid match is tested against BOTH
    fractional coordinates, never inferred from W(Q)=0.
    """
    q = finite_array(q, (2,))
    delta = wrapped_difference(q, grid)
    matches = np.flatnonzero(np.max(abs(delta), axis=1) < TOL)
    if len(matches):
        w = np.zeros(len(grid)); w[matches[0]] = 1.
        return w, int(matches[0])
    amplitudes = np.ones(len(grid), complex)
    for axis, size in enumerate((nx, ny)):
        t = delta[:, axis]
        amplitudes *= np.exp(-1j * np.pi * (size - 1) * t) * np.sinc(size * t) / np.sinc(t)
    weights = abs(amplitudes) ** 2
    if abs(weights.sum() - 1.) > 5e-10:
        raise ValueError("window completeness check failed; grid/representatives disagree")
    return weights, None


def matrices(srat, drat, weights, nc):
    """Return [replica, estimator(raw, connected, projected), mu,nu]."""
    raw = np.einsum("k,rkij->rij", weights, srat)
    connected = raw - nc * weights[0] * drat[:, :, None] * drat[:, None, :].conj()
    projected = raw - weights[0] * srat[:, 0]
    return np.stack([raw, connected, projected], axis=1)


def scalar_channels(m):
    h = (m + np.swapaxes(m.conj(), -1, -2)) / 2
    return np.stack([m[..., 2, 2].real, np.trace(m, axis1=-2, axis2=-1).real,
                     np.linalg.eigvalsh(h)[..., -1]], axis=-1)


def targets_for_run(meta, qfile=None):
    targets = []
    def add(label, q):
        q = finite_array(q, (2,)) % 1
        for t in targets:
            if np.max(abs(wrapped_difference(q, t["q"]))) < TOL:
                if label not in t["label"].split("|"):
                    t["label"] += "|" + label
                return
        targets.append(dict(label=label, q=q))
    couplings = finite_array([meta[k] for k in ("Jx", "Jy", "Jz")], (3,))
    nodes = dirac_nodes(*couplings)
    if nodes["status"] == "ok":
        for i, q in enumerate(q_predictions(nodes["nodes"])):
            add(f"Qbare_{i}", [q["u"], q["v"]])
    for p in meta["vbs"].get("qpoints", []):
        add(p["label"], p["reciprocal_fraction"])
    if qfile:
        for number, line in enumerate(Path(qfile).read_text(encoding="utf-8").splitlines(), 1):
            line = line.split("#")[0].split("//")[0].replace(",", " ").strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"q-file line {number}: expected label u v")
            add(fields[0], fields[1:])
    return targets, nodes


def segment_images(q, a, b):
    """Representatives of a reciprocal-equivalent point lying on a specified segment."""
    q = np.asarray(q) % 1
    d = b - a
    norm2 = float(d @ d)
    if norm2 <= 1e-24:
        raise ValueError("consecutive path vertices must differ")
    lo = np.ceil(np.minimum(a, b) - q - TOL).astype(int)
    hi = np.floor(np.maximum(a, b) - q + TOL).astype(int)
    sizes = np.maximum(hi - lo + 1, 0)
    if int(sizes[0]) * int(sizes[1]) > 10000:
        raise ValueError("path spans too many reciprocal cells; specify a shorter path")
    for shift in itertools.product(*(range(lo[i], hi[i] + 1) for i in range(2))):
        r = q + shift
        t = float((r - a) @ d / norm2)
        if -TOL <= t <= 1 + TOL and np.linalg.norm(r - a - t*d) < TOL:
            yield float(np.clip(t, 0, 1)), r


def make_path(vertices, reciprocal, samples, grid, targets):
    vertices = finite_array(vertices)
    if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 2:
        raise ValueError("need >=2 finite two-component vertices and >=2 samples per segment")
    if np.isscalar(samples):
        samples = [int(samples)] * (len(vertices) - 1)
    else:
        samples = [int(value) for value in samples]
    if len(samples) != len(vertices) - 1 or any(value < 2 for value in samples):
        raise ValueError("need one sampling density >=2 for each path segment")
    points, ticks, distance = [], [0.], 0.
    for segment_index, (a, b) in enumerate(zip(vertices[:-1], vertices[1:])):
        length = float(np.linalg.norm((b-a) @ reciprocal))
        entries = [(float(t), a+t*(b-a), "") for t in
                   np.linspace(0, 1, samples[segment_index])]
        for q in grid:
            entries.extend((t, r, "") for t, r in segment_images(q, a, b))
        for target in targets:
            entries.extend((t, r, target["label"]) for t, r in segment_images(target["q"], a, b))
        entries.sort(key=lambda x: x[0])
        for t, q, label in entries:
            s = distance + t * length
            if points and abs(points[-1]["distance"] - s) < 1e-11:
                if label:
                    points[-1]["label"] = "|".join(filter(None, [points[-1]["label"], label]))
                    points[-1]["q"] = q  # retain the exact labeled candidate rather than a rounded sample
                continue
            points.append(dict(q=q, distance=s, label=label))
        distance += length
        ticks.append(distance)
    return points, ticks


def evaluate(points, kind, run, srat, drat):
    nx, ny = run["meta"]["Lx"], run["meta"]["Ly"]
    nc = nx * ny
    rows, centers, real_errors, imag_errors = [], [], [], []
    for p in points:
        q = np.asarray(p["q"])
        weights, index = window_weights(q, run["grid"], nx, ny)
        m = matrices(srat, drat, weights, nc)
        scalars = scalar_channels(m)
        errors = jackknife_se(scalars[1:])
        qcart = q @ run["reciprocal"]
        row = dict(kind=kind, label=p.get("label", ""), u=float(q[0]), v=float(q[1]),
                   qx=float(qcart[0]), qy=float(qcart[1]), distance=p.get("distance", ""),
                   grid_index="" if index is None else index, on_grid=int(index is not None),
                   is_q0=int(index == 0), W_abs2=float(weights[0]),
                   largest_window_weight=float(weights.max()),
                   window_weight_sum=float(weights.sum()),
                   window_weight_min=float(weights.min()),
                   window_weight_nonnegative=int(np.min(weights) >= -5e-14))
        for e, estimator in enumerate(ESTIMATORS):
            for c, channel in enumerate(CHANNELS):
                key = estimator + "_" + channel
                row[key] = float(scalars[0, e, c])
                row[key + "_se"] = float(errors[e, c])
                row[key + "_per_cell"] = float(scalars[0, e, c] / nc)
                row[key + "_per_cell_se"] = float(errors[e, c] / nc)
            row[estimator + "_hermitian_residual"] = float(np.linalg.norm(m[0,e] - m[0,e].conj().T))
            row[estimator + "_max_imag_diagonal"] = float(abs(np.diag(m[0,e]).imag).max())
        rows.append(row)
        centers.append(m[0]); real_errors.append(jackknife_se(m[1:].real)); imag_errors.append(jackknife_se(m[1:].imag))
    return rows, centers, real_errors, imag_errors


def evaluate_fixed_associations(run, explicit_q=None, min_blocks=8):
    """Evaluate fixed Q_D/−Q_D and four directional neighbors jointly."""
    groups = fixed_vbs_groups(run["meta"], explicit_q=explicit_q)
    unique = []
    lookup = {}
    for group in groups:
        for role, q in group["points"].items():
            key = tuple(np.round(q % 1.0, 14))
            if key not in lookup:
                lookup[key] = len(unique)
                unique.append(q % 1.0)
            group.setdefault("indices", {})[role] = lookup[key]
    weights, on_grid, w_abs = [], [], []
    for q in unique:
        w, idx = window_weights(q, run["grid"], run["meta"]["Lx"], run["meta"]["Ly"])
        weights.append(w)
        on_grid.append(idx)
        w_abs.append(float(np.sqrt(max(float(w[0]), 0.0))))
    weights = np.asarray(weights, float)
    # Each selected point is reconstructed from the complete allowed momentum
    # block, while its connected mean is W(q)*Delta(0).  The phase of W drops
    # out of the outer product, so the positive magnitude is sufficient here.
    selected_s = np.einsum("qk,bkij->bqij", weights, run["snum"], optimize=True)
    selected_d = np.asarray(w_abs)[None, :, None] * run["dnum"][:, None, :]
    point_indices = {"center": 0, "plus_b1": 1, "minus_b1": 2,
                     "plus_b2": 3, "minus_b2": 4}
    rows = []
    for group in groups:
        # Reorder the unique-point table into the public helper's required slots.
        order = [group["indices"][role] for role in
                 ("center", "plus_b1", "minus_b1", "plus_b2", "minus_b2")]
        local_s = selected_s[:, order]
        local_d = selected_d[:, order]
        result = connected_window_association_ratios(
            run["phase"], local_s, local_d, point_indices,
            Ncell=run["meta"]["vbs"]["Nc"], min_blocks=min_blocks)
        cvals = result["connected_selected_matrices"][:, 2, 2].real
        cse = jackknife_se(result["replica_connected_szz_by_q"])
        row = {"center_label": group["center_label"],
               "source": group["source"],
               "step_source": group["step_source"],
               "center_u": float(group["center"][0]),
               "center_v": float(group["center"][1]),
               "center_qx": float(group["center"] @ run["reciprocal"][:, 0]),
               "center_qy": float(group["center"] @ run["reciprocal"][:, 1]),
               "Ncell": int(run["meta"]["vbs"]["Nc"]),
               "block_samples": int(run["block_count"]),
               "n_blocks": int(len(run["phase"])),
               "phase_resolved_5sigma": bool(abs(run["phase"].sum() /
                   (len(run["phase"]) * run["block_count"])) >
                   5 * np.linalg.norm(ratio_replicas(run)[2]["average_phase_se"])),
               "center_on_grid": bool(on_grid[group["indices"]["center"]] is not None),
               "center_grid_index": "" if on_grid[group["indices"]["center"]] is None else int(on_grid[group["indices"]["center"]]),
               "center_W_abs2": float(weights[group["indices"]["center"], 0]),
               "center_is_q0": bool(on_grid[group["indices"]["center"]] == 0)}
        for role, step_key in (("plus_b1", "b1_plus"), ("minus_b1", "b1_minus"),
                               ("plus_b2", "b2_plus"), ("minus_b2", "b2_minus")):
            idx = group["indices"][role]
            q = group["points"][role]
            row[f"{step_key}_u"] = float(q[0]); row[f"{step_key}_v"] = float(q[1])
            row[f"{step_key}_on_grid"] = bool(on_grid[idx] is not None)
            row[f"{step_key}_grid_index"] = "" if on_grid[idx] is None else int(on_grid[idx])
            row[f"{step_key}_W_abs2"] = float(weights[idx, 0])
            row[f"{step_key}_is_q0"] = bool(on_grid[idx] == 0)
        # Values in the helper's local order, including the center.
        local_c = cvals
        local_se = cse
        labels = {"center": "center", "plus_b1": "b1_plus", "minus_b1": "b1_minus",
                  "plus_b2": "b2_plus", "minus_b2": "b2_minus"}
        for role, name in labels.items():
            j = ("center", "plus_b1", "minus_b1", "plus_b2", "minus_b2").index(role)
            row[f"{name}_C"] = float(local_c[j]); row[f"{name}_C_se"] = float(local_se[j])
            row[f"{name}_m2_z"] = float(local_c[j] / row["Ncell"])
            row[f"{name}_m2_z_se"] = float(local_se[j] / row["Ncell"])
        for axis in ("b1", "b2"):
            rr = result["ratios"][axis]
            applicable = not any(row[f"{axis}_{sign}_is_q0"] for sign in ("plus", "minus")) and not row["center_is_q0"]
            row[f"R_{axis}"] = rr["estimate"][0] if applicable and rr["estimate"] is not None else ""
            row[f"R_{axis}_imag"] = rr["estimate"][1] if applicable and rr["estimate"] is not None else ""
            row[f"R_{axis}_se"] = rr["stderr_real"] if applicable and rr["estimate"] is not None else ""
            row[f"R_{axis}_available"] = bool(applicable and rr["available"])
            row[f"R_{axis}_unavailable_reason"] = (None if row[f"R_{axis}_available"] else
                ("q0_equivalent_in_group" if not applicable else rr["unavailable_reason"]))
        row["center_C_reliable"] = bool(result["ratios"]["b1"]["available"] or result["ratios"]["b2"]["available"])
        row["center_C_jackknife_se"] = float(result["ratios"]["b1"]["center_C_jackknife_se"])
        conjugate = (-group["center"]) % 1.0
        for role, name in (("plus_b1", "b1_plus"), ("minus_b1", "b1_minus"),
                           ("plus_b2", "b2_plus"), ("minus_b2", "b2_minus")):
            row[f"{name}_distance_to_conjugate_center"] = float(
                reciprocal_distance(group["points"][role], conjugate, run["reciprocal"])[0])
        row["delta1_u"] = float(group["step_b1"][0]); row["delta1_v"] = float(group["step_b1"][1])
        row["delta2_u"] = float(group["step_b2"][0]); row["delta2_v"] = float(group["step_b2"][1])
        d1cart = group["step_b1"] @ run["reciprocal"]
        d2cart = group["step_b2"] @ run["reciprocal"]
        row["delta1_qx"] = float(d1cart[0]); row["delta1_qy"] = float(d1cart[1])
        row["delta2_qx"] = float(d2cart[0]); row["delta2_qy"] = float(d2cart[1])
        rows.append(row)
    return rows, groups


def plot_tables(path_rows, grid_rows, ticks, vertices, prefix, estimator, channel, per_cell):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    key = estimator + "_" + channel + ("_per_cell" if per_cell else "")
    label = {"szz": r"$S_{zz}$", "trace": r"$\mathrm{Tr}\,S$", "lambda_max": r"$\lambda_{\max}(S)$"}[channel]
    if per_cell:
        label += r"$/N_c$"
    x = np.array([r["distance"] for r in path_rows])
    y = np.array([r[key] if not r["is_q0"] else np.nan for r in path_rows])
    err = np.array([r[key+"_se"] if not r["is_q0"] else np.nan for r in path_rows])
    fig, ax = plt.subplots(figsize=(8.6, 4.9))
    ax.plot(x, y, label="Translation-averaged window")
    ax.fill_between(x, y-err, y+err, alpha=.2, label="Nominal block-jackknife SE")
    on_grid = np.array([bool(r["on_grid"] and not r["is_q0"]) for r in path_rows])
    ax.plot(x[on_grid], y[on_grid], "o", markersize=4, label="Allowed grid points on path")
    for r in path_rows:
        if "Qbare_" in r["label"]:
            ax.axvline(r["distance"], linestyle="--", linewidth=.8)
            ax.text(r["distance"], .98, r["label"].split("|")[0], transform=ax.get_xaxis_transform(),
                    rotation=90, va="top", ha="right", fontsize=8)
    ax.set_xticks(ticks, [f"({u:g}, {v:g})" for u,v in vertices])
    ax.set_xlabel("Momentum-path distance; ticks are reciprocal fractions")
    ax.set_ylabel(label)
    ax.set_title(f"{estimator}: q=0 omitted; convergence not certified")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(str(prefix)+"_path.png", dpi=170); plt.close(fig)
    rows = [r for r in grid_rows if not r["is_q0"]]
    fig, ax = plt.subplots(figsize=(6.7, 5.4))
    points = ax.scatter([r["qx"] for r in rows], [r["qy"] for r in rows], c=[r[key] for r in rows])
    fig.colorbar(points, ax=ax, label=label)
    ax.set(xlabel=r"$q_x$", ylabel=r"$q_y$", title=f"{estimator}: allowed grid, q=0 omitted")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout(); fig.savefig(str(prefix)+"_grid.png", dpi=170); plt.close(fig)


def export(results, blocks, prefix, vertices, samples=161, reblock=1,
           qfile=None, plot=False, estimator="connected", channel="lambda_max", per_cell=False,
           fixed_q=None, fixed_min_blocks=8):
    run = read_run(results, blocks, reblock)
    srat, drat, diagnostics = ratio_replicas(run)
    targets, nodes = targets_for_run(run["meta"], qfile)
    path, ticks = make_path(vertices, run["reciprocal"], samples, run["grid"], targets)
    segment_samples = ([int(samples)] * (len(vertices) - 1) if np.isscalar(samples)
                       else [int(value) for value in samples])
    association_rows, _ = evaluate_fixed_associations(
        run, explicit_q=fixed_q, min_blocks=fixed_min_blocks)
    sets = [("grid", [dict(q=q) for q in run["grid"]]), ("path", path), ("targets", targets)]
    prefix = Path(prefix)
    suffixes = ["_grid.csv", "_path.csv", "_targets.csv", "_association.csv", "_matrices.npz", "_metadata.json"]
    if plot:
        suffixes += ["_path.png", "_grid.png"]
        import matplotlib  # fail before writing outputs if plotting dependency is missing
    outputs = [Path(str(prefix) + s) for s in suffixes]
    if any(p.exists() for p in outputs):
        raise FileExistsError("refusing to overwrite existing momentum output; choose a new --output-prefix")
    tables, all_rows, matrix_data = {}, [], [[], [], []]
    for kind, points in sets:
        rows, *arrays = evaluate(points, kind, run, srat, drat)
        tables[kind] = rows
        all_rows.extend(rows)
        for out, values in zip(matrix_data, arrays):
            out.extend(values)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fields = list(tables["grid"][0])
    for kind, rows in tables.items():
        with Path(str(prefix)+f"_{kind}.csv").open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
    with Path(str(prefix)+"_association.csv").open("x", newline="", encoding="utf-8") as stream:
        association_fields = list(association_rows[0]) if association_rows else []
        writer = csv.DictWriter(stream, fieldnames=association_fields)
        writer.writeheader(); writer.writerows(association_rows)
    with Path(str(prefix)+"_matrices.npz").open("xb") as stream:
        np.savez_compressed(stream, estimator=np.array(ESTIMATORS),
            kind=np.array([r["kind"] for r in all_rows]), uv=np.array([[r["u"],r["v"]] for r in all_rows]),
            matrix_center=np.array(matrix_data[0]), matrix_se_real=np.array(matrix_data[1]),
            matrix_se_imag=np.array(matrix_data[2]))
    metadata = dict(schema="vbs-momentum-spectrum-v2", source_results=run["results_path"],
        source_blocks=run["blocks_path"],
        input_sha256={"results": sha256_file(run["results_path"]),
                      "blocks": sha256_file(run["blocks_path"])},
        reciprocal_basis=np.asarray(run["reciprocal"]).tolist(),
        parameters={k: run["meta"].get(k) for k in
            ("Lx","Ly","shear","geometry_id","Jx","Jy","Jz","g","mass","omega_ph","k","lambda","intervals","dt","tau","seed")},
        statistics=diagnostics, bare_dirac_prediction=nodes,
        uniform_mean_Delta0=[[float(z.real),float(z.imag)] for z in drat[0]],
        path_vertices=np.asarray(vertices).tolist(), path_vertex_distances=ticks,
        path_points_per_segment=segment_samples,
        rows={k:len(v) for k,v in tables.items()}, estimator_order=list(ESTIMATORS),
        window_association_rows=len(association_rows),
        window_association_file=str(Path(str(prefix)+"_association.csv").name),
        fixed_vbs_association={"groups": association_rows,
                               "definition": "R_win=1-[C(Q+delta)+C(Q-delta)]/[2*C(Q)]",
                               "channel": "connected Szz real; m2_z=C/Nc",
                               "centers_are_Q_D_and_minus_Q_D": True},
        matrix_shape="(points in grid/path/targets CSV order, 3 estimators, 3, 3)",
        normalization="F=sum_R exp(-iq.R) X; S=<F F^dagger>/Nc; *_per_cell=S_scalar/Nc; X is not gX",
        raw_definition="S_T(Q)=sum_k |W(Q-k)|^2 S_raw(k)",
        connected_definition="S_T,c(Q)=S_T(Q)-Nc*|W(Q)|^2 <Delta(0)><Delta(0)>^dagger",
        projected_definition="S_T,prj(Q)=S_T(Q)-|W(Q)|^2 S_raw(0)",
        center="full-sample plug-in; no jackknife bias correction; no negative-value clipping",
        errors="equal-block joint delete-one jackknife, common complex phase denominator, matrix first then eigenvalue",
        q0_handling={"retained_in_csv_and_npz": True, "plot_exclusion_field": "is_q0",
                     "plot_exclusion_rule": "omit rows with is_q0=1",
                     "reason": "q=0 is retained for connected background and excluded from nonzero-momentum plots and peak searches"},
        warnings=["Off-grid values are translation averages reconstructed from the full grid, not fixed-origin Ssel.",
          "Equality to the fixed-origin ensemble expectation assumes translation-invariant Hamiltonian, boundary/gauge treatment and trial density matrix.",
          "q=0 is retained in CSV/NPZ and the background calculation, but excluded from the plots; do not use it in nonzero-momentum peak searches.",
          "No smoothing, interpolation of eigenvalues, covariance independence assumption, or white-noise normalization was applied.",
          "The window and commensurability remain size dependent; denser path sampling is not finer physical momentum resolution.",
          "Error bars do not establish thermalization, independent blocks or projection convergence; a narrow peak alone is not a thermodynamic-order proof.",
          "This exporter processes exactly one chain. Do not concatenate seeds and reblock across seed boundaries."])
    with Path(str(prefix)+"_metadata.json").open("x", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
    if plot:
        plot_tables(tables["path"], tables["grid"], ticks, vertices, prefix, estimator, channel, per_cell)
    return metadata


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results", help="one completed PREFIX_results.json")
    ap.add_argument("--blocks", help="defaults to corresponding PREFIX_blocks.csv")
    ap.add_argument("--output-prefix", required=True, help="new output prefix; never overwrites")
    ap.add_argument("--vertices", default="0,1;1,0", help="semicolon-separated reciprocal fractions, e.g. '0,1;1,0'")
    ap.add_argument("--points-per-segment", default="161",
                    help="one integer for every segment, or comma-separated per-segment integers")
    ap.add_argument("--reblock", type=int, default=1, help="merge this many adjacent blocks within this one chain")
    ap.add_argument("--q-file", help="additional label/u/v candidates; Qbare and existing selected points are included automatically")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--estimator", choices=ESTIMATORS, default="connected", help="plot choice only; all three estimators are exported")
    ap.add_argument("--channel", choices=CHANNELS, default="lambda_max")
    ap.add_argument("--per-cell", action="store_true", help="plot S_scalar/Nc; both normalizations are always exported")
    ap.add_argument("--fixed-q", metavar="U,V", help="explicit preselected Q_D fractional coordinates; -Q_D is added automatically")
    ap.add_argument("--fixed-min-blocks", type=int, default=8,
                    help="minimum equal blocks for an available R_win estimate")
    args = ap.parse_args()
    try:
        vertices = [[float(x) for x in s.split(",")] for s in args.vertices.split(";")]
        fixed_q = None if args.fixed_q is None else [float(x) for x in args.fixed_q.replace(",", " ").split()]
        if fixed_q is not None and len(fixed_q) != 2:
            raise ValueError("--fixed-q expects U,V")
        sample_fields = [int(value.strip()) for value in args.points_per_segment.split(",") if value.strip()]
        if not sample_fields:
            raise ValueError("--points-per-segment must contain an integer")
        samples = sample_fields[0] if len(sample_fields) == 1 else sample_fields
        report = export(args.results, args.blocks, args.output_prefix, vertices,
            samples, args.reblock, args.q_file, args.plot,
            args.estimator, args.channel, args.per_cell, fixed_q, args.fixed_min_blocks)
    except (ValueError, OSError, KeyError, TypeError) as error:
        ap.exit(2, f"error: {error}\n")
    print(f"Exported {report['rows']} to {args.output_prefix}_*; no QMC updates were run.")
    if not report["statistics"]["phase_resolved_5sigma"]:
        print("WARNING: average phase is not resolved at 5 sigma; ratio estimates are only diagnostics.")
    print("Nominal block errors are not a convergence certificate. Inspect metadata and reblocking.")


if __name__ == "__main__":
    main()
