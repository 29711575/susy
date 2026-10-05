#!/usr/bin/env python3
"""Fixed-grid seed / projection / step diagnostics, not a convergence test.

The first run fixes all momenta (and the exploratory candidate markers).
Errors are joint phase jackknife errors at a saved block level. Combining
independent-run errors in quadrature is not valid for reused random streams;
this script therefore requires different seeds.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from plot_vbs import analyze, choose_level, error, load_run, value, write_csv
from convergence_analysis import check_trial_comparison
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-samples", type=int)
    args = parser.parse_args()
    runs = [load_run(path)[0] for path in args.results]
    if len(runs) < 2 or len({r["seed"] for r in runs}) != len(runs):
        raise ValueError("need at least two independent seeds")
    for r in runs[1:]:
        check_trial_comparison(runs[0],r)
        for key in ("Lx", "Ly", "Jx", "Jy", "Jz", "mass", "omega_ph", "g"):
            if r[key] != runs[0][key]:
                raise ValueError("seed/tau/h comparison requires same geometry and Hamiltonian parameters")
    levels = [choose_level(r, args.block_samples) for r in runs]
    first, _, _ = analyze(runs[0], levels[0])
    rows = []
    for run, level, path in zip(runs, levels, args.results):
        for q in range(run["vbs"]["Nc"]):
            row = {"run": path.stem, "seed": run["seed"], "trial_shift":run['trial_shift'], "tau": run["tau"], "h": run["dt"],
                   "q": q, "nx": run["vbs"]["momenta"][q]["n"][0], "ny": run["vbs"]["momenta"][q]["n"][1],
                   "block_sweeps": level["block_sweeps"], "blocks": level["blocks"], "elapsed_seconds": run["elapsed_seconds"]}
            for channel, group, name in (("Szz", "raw", f"S_q{q}_zz"), ("TrS", "derived", f"TrS_q{q}")):
                estimate, reference = level[group][name], levels[0][group][name]
                row[channel] = value(estimate).real
                row[channel+"_stderr"] = error(estimate)
                combined = np.hypot(error(estimate), error(reference))
                row[channel+"_difference_reference"] = value(estimate).real-value(reference).real
                row[channel+"_difference_zscore"] = row[channel+"_difference_reference"]/combined if combined > 0 else None
            rows.append(row)
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be absent or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "comparison.csv", rows)
    summary = {"reference": str(args.results[0]), "fixed_exploratory_candidates": first["candidates"],
               "convergence_established": False, "interpretation": "overlapping error bars only bound detectable differences; block independence unproved",
               "runs": [{"path": str(p), "seed": r["seed"], "blocks": l["blocks"], "block_sweeps": l["block_sweeps"],
                         "phase_resolved": l["phase_resolved_5sigma"], "elapsed_seconds": r["elapsed_seconds"], "proposal": r["proposal"]}
                        for p, r, l in zip(args.results, runs, levels)]}
    (args.output / "comparison.json").write_text(json.dumps(summary, indent=2)+"\n")
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    nc = runs[0]["vbs"]["Nc"]
    for ax, channel in zip(axes, ("Szz", "TrS")):
        for i, run in enumerate(runs):
            points = rows[i*nc:(i+1)*nc]
            ax.errorbar(np.arange(nc)+(i-(len(runs)-1)/2)*.07,
                        [p[channel] for p in points], yerr=[p[channel+"_stderr"] for p in points], fmt="o", capsize=3,
                        label=f"seed {run['seed']}, d={run['trial_shift']}, tau={run['tau']:g}, h={run['dt']:g}, "
                              +( "bond-refresh" if run["proposal"]=="free_whole_bond_refresh" else "local"))
        ax.axvline(first["candidates"][channel], color="grey", ls=":", label="candidate fixed from first run")
        ax.set(xticks=range(nc), xticklabels=[str(tuple(m["n"])) for m in runs[0]["vbs"]["momenta"]],
               xlabel="fixed discrete (nx,ny)", ylabel=channel)
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.2))
    p=runs[0]
    status = f"blocks={[l['blocks'] for l in levels]}, phase resolved={[l['phase_resolved_5sigma'] for l in levels]}"
    fig.suptitle(f"Actual MC {p['Lx']}x{p['Ly']}, J=({p['Jx']},{p['Jy']},{p['Jz']}), g={p['g']}, m={p['mass']}, omega={p['omega_ph']}\n"
                 f"Joint JK; {status}\nBlock independence / projection / step convergence NOT established", fontsize=9)
    fig.savefig(args.output / "comparison.png", dpi=170)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
