#!/usr/bin/env python3
"""Plot measured P4 output. No synthetic configurations, interpolation or peak alignment.

Requires numpy and matplotlib. All nonlinear errors come from the C++ joint
phase/block jackknife; this reader does not replace it with independent errors.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def normalize_trial_metadata(data):
    # Known old P3-P5 files have no shift field and used the undisplaced vacuum.
    if 'trial_shift' not in data and 'phonon_trial' in data:
        raise ValueError('trial metadata missing displacement')
    d=data.setdefault('trial_shift',[0.,0.,0.])
    if len(d)!=3 or not np.isfinite(np.asarray(d,dtype=float)).all():
        raise ValueError('invalid trial_shift')
    if data.get('phonon_trial','equal_shifted_vacuum_endpoints')!='equal_shifted_vacuum_endpoints':
        raise ValueError('unsupported phonon trial')
    if data.get('matter_trial','free_ground_state_static_X_zero')!='free_ground_state_static_X_zero':
        raise ValueError('unsupported matter trial')
    return data


def normalize_tau_metadata(data):
    """Expose canonical total open-contour ``tau`` while reading old files."""
    candidates = []
    if "tau" in data:
        candidates.append(("tau", float(data["tau"])))
    if "tau_per_branch" in data:
        candidates.append(("tau_per_branch", 2*float(data["tau_per_branch"])))
    # Read-only compatibility for results produced before the tau rename.
    if "theta_per_branch" in data:
        candidates.append(("legacy theta_per_branch", 2*float(data["theta_per_branch"])))
    if "total_time" in data:
        candidates.append(("legacy total_time", float(data["total_time"])))
    if "intervals" in data and "dt" in data:
        candidates.append(("intervals*dt", float(data["intervals"])*float(data["dt"])))
    if not candidates or not all(np.isfinite(value) and value >= 0 for _, value in candidates):
        raise ValueError("missing or invalid imaginary-time metadata")
    tau = candidates[0][1]
    for name, value in candidates[1:]:
        if not np.isclose(value, tau, rtol=0, atol=1e-12):
            raise ValueError(f"inconsistent imaginary-time metadata: {name}")
    data["tau"] = tau
    return data


def normalize_geometry_metadata(data):
    """Accept known historical rectangular files; reject incomplete new geometry."""
    if "shear" not in data:
        if any(key in data for key in ("pbc_geometry", "torus_matrix",
                                       "geometry_id", "momentum_grid_convention")):
            raise ValueError("geometry metadata is incomplete: missing shear")
        data["shear"] = 0
        data["pbc_geometry"] = "axis_aligned"
        data["geometry_id"] = f"honeycomb-pbc-Lx{data['Lx']}-Ly{data['Ly']}-shear0"
        return data
    shear = data["shear"]
    if isinstance(shear, bool) or not isinstance(shear, int) or not 0 <= shear < data["Lx"]:
        raise ValueError("invalid shear metadata")
    expected_kind = "axis_aligned" if shear == 0 else "sheared"
    expected_id = f"honeycomb-pbc-Lx{data['Lx']}-Ly{data['Ly']}-shear{shear}"
    if data.get("pbc_geometry") != expected_kind or data.get("geometry_id") != expected_id:
        raise ValueError("invalid or incomplete geometry identity")
    return data


def value(estimate):
    if any(x is None for x in estimate["value"]):
        raise ValueError("undefined phase ratio; refusing to plot a physical estimate")
    return complex(*estimate["value"])


def error(estimate):
    x = estimate["stderr_real"]
    return np.nan if x is None else x


def load_run(path):
    path = Path(path)
    data = json.loads(path.read_text())
    if data.get("status") != "complete":
        raise ValueError("refusing incomplete or aborted MC results")
    normalize_trial_metadata(data)
    normalize_tau_metadata(data)
    normalize_geometry_metadata(data)
    if data.get("status") != "complete" or not data.get("coupled_mc") or not data.get("vbs_enabled"):
        raise ValueError("need completed actual MC with --vbs 1")
    if not data["ratio_defined"]:
        raise ValueError("phase denominator is numerically zero; raw blocks retained, no physical plot")
    # Check the primary JSON against the saved joint phase/numerator sums.
    blocks_path = path.with_name(path.name.removesuffix("_results.json") + "_blocks.csv")
    with blocks_path.open() as stream:
        blocks = list(csv.DictReader(stream))
    if not blocks or len(blocks) != data["blocks"]:
        raise ValueError("missing/inconsistent block file")
    counts = [int(b["count"]) for b in blocks]
    if len(set(counts)) != 1 or sum(counts) != data["samples"]:
        raise ValueError("unequal or incomplete blocks")
    phase = sum(complex(float(b["sum_phase_real"]), float(b["sum_phase_imag"])) for b in blocks)
    max_difference = 0.
    for name, estimate in data["observables"].items():
        numerator = sum(complex(float(b[f"sum_phase_{name}_real"]), float(b[f"sum_phase_{name}_imag"])) for b in blocks)
        max_difference = max(max_difference, abs(numerator / phase - value(estimate)))
    scale = max(1., *(abs(value(e)) for e in data["observables"].values()))
    if max_difference > 2e-12 * scale:
        raise ValueError("block sums disagree with JSON estimates")
    return data, max_difference


def choose_level(data, samples):
    levels = data["vbs"]["blocking"]
    if samples is not None:
        matches = [x for x in levels if x["block_samples"] == samples]
        if not matches:
            raise ValueError("requested block size unavailable; no samples are discarded to create it")
        return matches[0]
    eligible = [x for x in levels if x["blocks"] >= 16]
    return eligible[-1] if eligible else levels[0]


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def analyze(data, level, candidate_source=None):
    raw, derived = level["raw"], level["derived"]
    nc = data["vbs"]["Nc"]
    s = np.array([[[value(raw[f"S_q{q}_{mu}{nu}"]) for nu in "xyz"] for mu in "xyz"] for q in range(nc)])
    candidates = {"Szz": int(np.argmax(s[1:, 2, 2].real)) + 1,
                  "TrS": int(np.argmax(np.trace(s[1:], axis1=1, axis2=2).real)) + 1}
    selection = "exploratory ensemble maxima, not significant peaks"
    if candidate_source:
        source = json.loads(Path(candidate_source).read_text())
        if (source["parameters"]["Lx"], source["parameters"]["Ly"], source["parameters"].get("shear", 0)) != (data["Lx"], data["Ly"], data["shear"]):
            raise ValueError("fixed candidate requires the same momentum grid")
        for key in ("Jx", "Jy", "Jz", "g", "mass", "omega_ph"):
            if source["parameters"][key] != data[key]:
                raise ValueError("independent candidate check requires the same Hamiltonian parameters")
        if source["parameters"]["seed"] == data["seed"]:
            raise ValueError("candidate checking requires an independent seed")
        candidates = source["candidates"]
        selection = f"fixed from {candidate_source}; independent-seed check, not peak discovery"
    rows = []
    for q, momentum in enumerate(data["vbs"]["momenta"]):
        row = {"q": q, "n1": momentum["n"][0], "n2": momentum["n"][1],
               "b1_fraction": momentum["reciprocal_fraction"][0], "b2_fraction": momentum["reciprocal_fraction"][1],
               "qx": momentum["cartesian"][0], "qy": momentum["cartesian"][1]}
        for mu in "xyz":
            for nu in "xyz":
                estimate = raw[f"S_q{q}_{mu}{nu}"]
                z = value(estimate)
                for key, val in [("real", z.real), ("imag", z.imag), ("stderr_real", estimate["stderr_real"]),
                                 ("stderr_imag", estimate["stderr_imag"]), ("cov_real_imag", estimate["cov_real_imag"])]:
                    row[f"S_{mu}{nu}_{key}"] = val
        for name in ("TrS", "mz2", "leading_eigenvalue", "hermitian_residual"):
            estimate = derived[f"{name}_q{q}"]
            row[name + "_real"] = value(estimate).real
            row[name + "_imag"] = value(estimate).imag
            row[name + "_stderr_real"] = estimate["stderr_real"]
        rows.append(row)
    correlations = []
    for d in range(nc):
        dx, dy = divmod(d, data["Ly"])
        row = {"d": d, "dx": dx, "dy": dy}
        for prefix in ("C", "C_connected"):
            for mu in "xyz":
                for nu in "xyz":
                    e = derived[f"{prefix}_q{d}_{mu}{nu}"]
                    z = value(e)
                    row[f"{prefix}_{mu}{nu}_real"] = z.real
                    row[f"{prefix}_{mu}{nu}_imag"] = z.imag
                    row[f"{prefix}_{mu}{nu}_stderr_real"] = e["stderr_real"]
                    row[f"{prefix}_{mu}{nu}_stderr_imag"] = e["stderr_imag"]
                    row[f"{prefix}_{mu}{nu}_cov_real_imag"] = e["cov_real_imag"]
        correlations.append(row)
    parseval = abs(np.trace(s, axis1=1, axis2=2).sum() - nc * sum(value(raw[f"X2_{mu}"]) for mu in "xyz"))
    parameters = {k: data[k] for k in ("Lx", "Ly", "shear", "pbc_geometry", "geometry_id", "Jx", "Jy", "Jz", "g", "lambda", "mass", "omega_ph", "intervals", "dt", "tau", "seed", "proposal", "trial_shift")}
    summary = {"parameters": parameters, "source": "actual coupled Monte Carlo; center node", "candidates": candidates,
               "candidate_selection": selection, "vbs_claim": "none; no established modulation peak or thermodynamic order",
               "statistics": {k: level[k] for k in ("blocks", "block_samples", "block_sweeps", "jackknife_defined", "phase_resolved_5sigma", "average_phase")},
               "block_independence_verified": False, "size_projection_step_convergence_verified": False,
               "max_hermitian_residual": float(max(np.linalg.norm(m - m.conj().T) for m in s)),
               "max_diagonal_imaginary_part": float(np.max(abs(s[:, range(3), range(3)].imag))),
               "parseval_error": float(parseval), "uniform_X": {mu: raw[f"Delta_q0_{mu}"] for mu in "xyz"},
               "normalization": data["vbs"]["normalization"], "elapsed_seconds": data["elapsed_seconds"], "acceptance": data["acceptance"],
               "raw": raw, "derived": derived}
    if data["g"] == 0:
        variance = 1 / (2 * data["mass"] * data["omega_ph"])
        zscores = []
        for q in range(nc):
            for mu in "xyz":
                for nu in "xyz":
                    e = raw[f"S_q{q}_{mu}{nu}"]
                    for part in ("real", "imag"):
                        err = e[f"stderr_{part}"]
                        exact = variance if part == "real" and mu == nu else 0.
                        if err is not None and err > 1e-14:
                            zscores.append(abs(getattr(value(e), part) - exact) / err)
        summary["g0_max_component_zscore"] = max(zscores) if zscores else None
        summary["g0_variance_exact"] = variance
    return summary, rows, correlations


def plots(data, level, summary, rows, correlations, output):
    raw, derived = level["raw"], level["derived"]
    nc = data["vbs"]["Nc"]
    status = (f"{level['blocks']} blocks x {level['block_sweeps']} sweeps; phase resolved={level['phase_resolved_5sigma']}; "
              "independence NOT established")
    if not level["jackknife_defined"] or level["blocks"] < 16:
        status = "INSUFFICIENT BLOCKS / UNCERTAIN ERRORS; " + status
    proposal = "bond-refresh" if data["proposal"] == "free_whole_bond_refresh" else "local"
    caption = (f"Actual MC {data['Lx']}x{data['Ly']} shear={data['shear']}, J=({data['Jx']},{data['Jy']},{data['Jz']}), g={data['g']}, "
               f"lambda={data['lambda']:.4g}, m={data['mass']}, omega={data['omega_ph']}, "
               f"tau={data['tau']:g}, h={data['dt']:g}, seed={data['seed']}\n"
               f"{proposal}; endpoint d={data['trial_shift']}; " + status)
    q = np.array([m["cartesian"] for m in data["vbs"]["momenta"]])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for ax, name, vals in zip(axes, ("Re Szz(q)", "Re Tr S(q)", "Largest eig[(S+S†)/2]"),
                            ([r["S_zz_real"] for r in rows], [r["TrS_real"] for r in rows], [r["leading_eigenvalue_real"] for r in rows])):
        dots = ax.scatter(q[:, 0], q[:, 1], c=vals, s=240, marker="s", cmap="viridis", edgecolors="black")
        for i, p in enumerate(q):
            ax.annotate(str(tuple(data["vbs"]["momenta"][i]["n"])), p, xytext=(0, 12), textcoords="offset points", ha="center", fontsize=8)
        ax.set(xlabel="qx", ylabel="qy", title=name)
        ax.margins(.3)
        ax.set_aspect("equal")
        fig.colorbar(dots, ax=ax, shrink=.8)
    fig.suptitle(caption + "\nAll allowed discrete points incl. q=0; no interpolation or established VBS peak", fontsize=10)
    fig.savefig(output / "structure_factor.png", dpi=170)
    plt.close(fig)

    a = np.array(data["vbs"]["primitive"])
    positions = []
    for d in range(nc):
        dx, dy = divmod(d, data["Ly"])
        if data["shear"] == 0:
            if dx > data["Lx"] // 2:
                dx -= data["Lx"]
            if dy > data["Ly"] // 2:
                dy -= data["Ly"]
        positions.append(np.array([dx, dy]) @ a)
    positions = np.array(positions)
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    for ax, channel in zip(axes, ("zz", "xz", "yz")):
        vals = np.array([r[f"C_connected_{channel}_real"] for r in correlations])
        # Keep d=0 visible and labeled, but do not let its oscillator variance
        # hide all off-site modulation. No value is clipped or changed.
        lim = max(abs(vals[1:]).max(), 1e-15)
        dots = ax.scatter(positions[1:, 0], positions[1:, 1], c=vals[1:], s=200, cmap="coolwarm", vmin=-lim, vmax=lim, edgecolors="black")
        ax.scatter(positions[0, 0], positions[0, 1], s=200, facecolors="none", edgecolors="black")
        for d, pos in enumerate(positions):
            e = error(derived[f"C_connected_q{d}_{channel}"])
            ax.annotate(f"{vals[d]:.2g}\n±{e:.2g}", pos, xytext=(0, 12), textcoords="offset points", ha="center", fontsize=7)
        ax.set(xlabel="d_x (real space)", ylabel="d_y (real space)", title=f"Re C_connected,{channel}(d)")
        ax.margins(.4)
        ax.set_aspect("equal")
        fig.colorbar(dots, ax=ax, shrink=.8)
    fig.suptitle(caption + "\nInverse FT; subtract <X_mu><X_nu>; joint JK errors; d=0 labeled but excluded from color scale", fontsize=10)
    fig.savefig(output / "connected_correlations.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    levels = [x for x in data["vbs"]["blocking"] if x["blocks"] >= 4 and x["jackknife_defined"]]
    for ax, channel in zip(axes, ("Szz", "TrS")):
        for qi in sorted({0, summary["candidates"][channel]}):
            group, name = ("raw", f"S_q{qi}_zz") if channel == "Szz" else ("derived", f"TrS_q{qi}")
            ax.plot([x["block_sweeps"] for x in levels], [error(x[group][name]) for x in levels], "o-", label=f"q index {qi}")
        ax.set(xscale="log", xlabel="sweeps per block", ylabel="joint JK standard error (real)", title=channel)
        ax.legend()
    fig.suptitle(caption + "\nBlocking diagnostic: few long blocks do not establish an error plateau", fontsize=9)
    fig.savefig(output / "blocking.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-samples", type=int, help="exact saved level; default largest level with >=16 blocks")
    parser.add_argument("--candidate-from", type=Path, help="vbs_analysis.json from a different seed; fixes exploration q")
    args = parser.parse_args()
    data, difference = load_run(args.results)
    level = choose_level(data, args.block_samples)
    summary, rows, correlations = analyze(data, level, args.candidate_from)
    summary["block_csv_ratio_max_error"] = difference
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("plot output directory must be absent or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "vbs_analysis.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv(args.output / "structure_factor.csv", rows)
    write_csv(args.output / "correlations.csv", correlations)
    plots(data, level, summary, rows, correlations, args.output)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("raw", "derived")}, indent=2))


if __name__ == "__main__":
    main()
