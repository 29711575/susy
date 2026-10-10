#!/usr/bin/env python3
"""Analyze complete blocks retained from an interrupted fermion scan."""
import argparse
import csv
import datetime
import json
import math
import pathlib
import re
from dataclasses import replace

import numpy as np

from convergence_analysis import Chain, estimates_json, joint_estimate


FERMION = re.compile(r"F_(bare|string)_dx(\d+)_dy(\d+)_([AB])([AB])$")


def read_partial(path, job, plan):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows or [int(row["block"]) for row in rows] != list(range(len(rows))):
        raise ValueError(f"missing or unordered complete blocks: {path}")
    if any(int(row["count"]) != plan["block_samples"] for row in rows):
        raise ValueError(f"incomplete block retained: {path}")
    fields = list(rows[0])
    names = [field.removeprefix("sum_phase_").removesuffix("_real")
             for field in fields
             if field.startswith("sum_phase_") and field.endswith("_real")
             and field != "sum_phase_real"]
    columns = ["sum_phase"] + ["sum_phase_" + name for name in names]
    sums = np.array([[complex(float(row[col + "_real"]), float(row[col + "_imag"]))
                      for col in columns] for row in rows])
    counts = np.array([int(row["count"]) for row in rows])
    metadata = {
        "seed": job["seed"], "Lx": plan["Lx"], "Ly": plan["Ly"],
        "Jx": plan["Jx"], "Jy": plan["Jy"], "Jz": plan["Jz"],
        "g": job["g"], "mass": plan["mass"], "omega_ph": plan["omega_ph"],
        "scope": plan["scope"], "spatial_boundary": "periodic",
        "wilson_loops": plan["wilson_loops"], "time_contour": "symmetric_projector",
        "intervals": job["intervals"], "dt": plan["dt"],
        "trial_shift": plan["trial_shift"], "measure_every": plan["measure_every"],
    }
    diagnostics = {key: sum(int(row[key]) for row in rows)
                   for key in ("positive_count", "negative_count", "complex_count")}
    if sum(diagnostics.values()) != int(counts.sum()):
        raise ValueError(f"phase classification does not cover samples: {path}")
    return Chain(path, metadata, names, sums, counts), diagnostics


def transform(sums, count, names, metadata):
    if abs(sums[0]) <= 1e-12 * count:
        raise ValueError("numerically unresolved phase denominator")
    ratios = sums[1:] / sums[0]
    out = {"phase": sums[0] / count, "phase_magnitude": abs(sums[0] / count)}
    out.update(zip(names, ratios))
    out.update((name + "_magnitude", abs(value))
               for name, value in zip(names, ratios) if name.startswith("F_"))
    return out


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def estimate_rows(chains, block_sweeps, scope, diagnostics):
    estimate = joint_estimate([chain.reblock(block_sweeps) for chain in chains], transform=transform)
    values = estimates_json(estimate)
    meta = chains[0].metadata
    rows = []
    for name, item in values.items():
        if name.endswith("_magnitude") or name in ("phase", "phase_magnitude"):
            continue
        value = complex(*item["value"])
        magnitude = values.get(name + "_magnitude")
        rows.append({
            "scope": scope, "g": meta["g"], "lambda": meta["g"] ** 2,
            "tau_total": meta["intervals"] * meta["dt"], "intervals": meta["intervals"],
            "seeds": ",".join(str(chain.seed) for chain in chains),
            "samples": estimate["samples"],
            "blocks_per_chain": ",".join(str(n) for n in estimate["chain_blocks"]),
            "block_sweeps": block_sweeps, "observable": name,
            "Re": value.real, "Im": value.imag, "magnitude": abs(value),
            "SE_Re": item["stderr_real"], "SE_Im": item["stderr_imag"],
            "SE": math.hypot(item["stderr_real"], item["stderr_imag"]),
            "cov_Re_Im": item["cov_real_imag"],
            "SE_magnitude": magnitude["stderr_real"] if magnitude else None,
            "phase_resolved_5sigma": estimate["phase_resolved_5sigma"],
        })
    phase = values["phase"]
    z = complex(*phase["value"])
    se = math.hypot(phase["stderr_real"], phase["stderr_imag"])
    phase_row = {
        "scope": scope, "g": meta["g"], "lambda": meta["g"] ** 2,
        "tau_total": meta["intervals"] * meta["dt"], "intervals": meta["intervals"],
        "seeds": ",".join(str(chain.seed) for chain in chains),
        "samples": estimate["samples"],
        "blocks_per_chain": ",".join(str(n) for n in estimate["chain_blocks"]),
        "block_sweeps": block_sweeps, "Re": z.real, "Im": z.imag,
        "magnitude": abs(z), "SE_Re": phase["stderr_real"],
        "SE_Im": phase["stderr_imag"], "SE": se,
        "SNR": abs(z) / se if se else None,
        "positive_count": sum(d["positive_count"] for d in diagnostics),
        "negative_count": sum(d["negative_count"] for d in diagnostics),
        "complex_count": sum(d["complex_count"] for d in diagnostics),
        "acceptance": None, "phase_resolved_5sigma": estimate["phase_resolved_5sigma"],
    }
    return rows, phase_row


def analytic_fermion(name):
    match = FERMION.fullmatch(name)
    if not match:
        return None
    kind, dx, dy, a, b = match.groups()
    dx, dy = int(dx), int(dy)
    value = complex(a == b and dx == 0 and dy == 0)
    if a != b:
        value = 0j
        for nx in range(4):
            for ny in range(4):
                qx, qy = 2 * np.pi * nx / 4, 2 * np.pi * ny / 4
                f = .5 + np.exp(-1j * qx) + np.exp(-1j * qy)
                factor = 1j * f if a == "A" else -1j * f.conjugate()
                value += np.exp(-1j * (qx * dx + qy * dy)) * factor / abs(f) / 16
    if kind == "string":
        value *= (-1) ** (dx + dy + int(a == "B" and b == "A"))
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((run / "manifest.json").read_text())
    plan = manifest["plan"]
    entries = []
    for command in manifest["commands"]:
        chain, diagnostics = read_partial(pathlib.Path(command["folder"]) / "run_blocks.csv",
                                          command["job"], plan)
        entries.append((chain, diagnostics, command["job"]))

    all_rows, fermion_rows, phase_rows = [], [], []
    blocking_fermions, blocking_phase = [], []
    for chain, diagnostics, _ in entries:
        rows, phase = estimate_rows([chain], 25, "interrupted_chain", [diagnostics])
        all_rows.extend(rows)
        fermion_rows.extend(row for row in rows if row["observable"].startswith("F_"))
        phase_rows.append(phase)

    pooled_fermions, pooled_phase = [], []
    for group in plan["groups"]:
        selected = [(chain, diag) for chain, diag, job in entries
                    if job["g"] == group["g"] and job["intervals"] == group["intervals"]]
        chains = [item[0] for item in selected]
        diagnostics = [item[1] for item in selected]
        rows, phase = estimate_rows(chains, 25, "interrupted_pooled", diagnostics)
        pooled_fermions.extend(row for row in rows if row["observable"].startswith("F_"))
        pooled_phase.append(phase)
        for block in (50, 100, 200):
            if all(len(chain.counts) % (block // 25) == 0 for chain in chains):
                block_rows, block_phase = estimate_rows(chains, block, "interrupted_pooled", diagnostics)
                blocking_fermions.extend(row for row in block_rows if row["observable"].startswith("F_"))
                blocking_phase.append(block_phase)

    for row in pooled_fermions:
        reference = analytic_fermion(row["observable"]) if row["g"] == 0 else None
        row["analytic_Re"] = reference.real if reference is not None else None
        row["analytic_Im"] = reference.imag if reference is not None else None
        row["analytic_abs_difference"] = abs(complex(row["Re"], row["Im"]) - reference) if reference is not None else None
    save_csv(output / "phase_per_chain_partial.csv", phase_rows)
    save_csv(output / "phase_pooled_partial.csv", pooled_phase)
    save_csv(output / "observables_per_chain_partial.csv", all_rows)
    save_csv(output / "fermions_per_chain_partial.csv", fermion_rows)
    save_csv(output / "fermions_pooled_partial.csv", pooled_fermions)
    if blocking_fermions:
        save_csv(output / "fermion_blocking_sensitivity_partial.csv", blocking_fermions)
    if blocking_phase:
        save_csv(output / "phase_blocking_sensitivity_partial.csv", blocking_phase)

    negatives = sum(row["negative_count"] for row in phase_rows)
    complexes = sum(row["complex_count"] for row in phase_rows)
    samples = sum(row["samples"] for row in phase_rows)
    g0_error = max(row["analytic_abs_difference"] for row in pooled_fermions
                   if row["analytic_abs_difference"] is not None)
    summary = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "interrupted_partial_complete_blocks_only",
        "planned_samples": len(entries) * plan["sweeps"], "retained_samples": samples,
        "negative_count": negatives, "complex_count": complexes,
        "g0_max_analytic_abs_difference": g0_error,
        "phase_per_chain": phase_rows, "phase_pooled": pooled_phase,
        "limitations": [
            "All processes ended without run_results.json or completion.json.",
            "Only fully flushed 25-sweep blocks are analyzed; chain lengths are unequal.",
            "Acceptance and final sampler counters were written only at normal exit and are unavailable.",
            "No exact RNG/state checkpoint exists, so this is not a resumed or completed scan.",
            "Block sensitivity is emitted only when every chain divides exactly; no samples are discarded.",
        ],
    }
    (output / "partial_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Interrupted L=4 fermion scan: retained complete blocks", "",
        f"Generated {summary['created_utc']}. This is a partial analysis of {samples}/{summary['planned_samples']} planned measurements. All ten chains were interrupted before normal exit.", "",
        "| g | lambda | tau_total | seeds | samples | blocks/chain | Re phase | Im phase | abs phase | SE | SNR | negative | complex | acceptance |",
        "|---:|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in pooled_phase:
        lines.append("| {g} | {lambda:.2g} | {tau_total} | {seeds} | {samples} | {blocks_per_chain} | {Re:.10g} | {Im:.10g} | {magnitude:.10g} | {SE:.3g} | {SNR} | {negative_count} | {complex_count} | unavailable |".format(**row))
    lines += [
        "", f"Across retained measurements: negative={negatives}, clearly complex={complexes}. 未观测到负号只描述这些有限且中断的链，不是严格证明无符号。",
        "", f"For g=0, the largest difference from the independent free Bloch result among 128 fermion channels is {g0_error:.3g}.",
        "", "fermions_pooled_partial.csv contains all 128 complex matter-Majorana channels for each parameter group, with Re/Im, component errors, covariance, magnitude, and jackknife magnitude error. These are C_ij=delta_ij-G_ij and specified-path product(u)C_ij, not spin correlations.",
        "", "The estimates reuse the existing chain-preserving joint jackknife on complete 25-sweep blocks. Unequal chain lengths are retained. Larger-block results are provided only where exact reblocking is possible without dropping samples.",
        "", "Acceptance is unavailable because final sampler counters are written only on normal process exit. The abrupt stop also prevents an exact continuation because no chain/RNG checkpoint was written. Treat stationarity, autocorrelation, projection-time convergence, and final precision as unresolved.",
    ]
    (output / "PARTIAL_REPORT.md").write_text("\n".join(lines) + "\n")
    print(output / "PARTIAL_REPORT.md")


if __name__ == "__main__":
    main()
