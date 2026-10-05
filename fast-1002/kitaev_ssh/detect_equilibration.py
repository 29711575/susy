#!/usr/bin/env python3
"""Choose a stable suffix of a QMC block time series.

The production executable may be run with ``--warmup 0``.  This tool then
examines the saved block means after the run and writes a conservative,
chain-local equilibration cut together with a reweighted summary of the
retained suffix.  The cut is a diagnostic heuristic, not a proof of
equilibrium: it requires both a suffix mean compatible with the final window
and a nonsignificant linear drift for every selected observable.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np


def _complex(row: dict[str, str], stem: str) -> complex:
    return complex(float(row[f"{stem}_real"]), float(row[f"{stem}_imag"]))


def read_blocks(path: Path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no blocks in {path}")
    names = [
        field.removeprefix("sum_phase_").removesuffix("_real")
        for field in rows[0]
        if field.startswith("sum_phase_")
        and field.endswith("_real")
        and field != "sum_phase_real"
    ]
    if not names:
        raise ValueError(f"no observables in {path}")
    phase = np.array([_complex(r, "sum_phase") for r in rows], dtype=np.complex128)
    nums = np.array(
        [[_complex(r, f"sum_phase_{name}") for name in names] for r in rows],
        dtype=np.complex128,
    )
    counts = np.array([int(r["count"]) for r in rows], dtype=int)
    if (np.any(counts <= 0) or len(set(counts.tolist())) != 1
            or not np.all(np.isfinite(phase)) or not np.all(np.isfinite(nums))):
        raise ValueError("invalid block count or nonfinite complex sum")
    unresolved = np.flatnonzero(np.abs(phase) <= 1e-12 * counts)
    if unresolved.size:
        raise ValueError("phase denominator unresolved in block(s): "
                         + ", ".join(map(str, unresolved.tolist())))
    values = np.full((len(rows), len(names)), np.nan, dtype=float)
    for i in range(len(rows)):
        if abs(phase[i]) > 0:
            values[i] = (nums[i] / phase[i]).real
    if not np.all(np.isfinite(values)):
        raise ValueError("a block has zero phase sum; cannot diagnose its observable")
    return rows, names, phase, nums, counts, values


def _slope_t(values: np.ndarray) -> float:
    """Absolute t statistic for a linear trend in one block series."""
    n = len(values)
    if n < 3:
        return 0.0
    x = np.arange(n, dtype=float)
    x -= x.mean()
    y = values - values.mean()
    xx = float(x @ x)
    beta = float(x @ y) / xx
    residual = y - beta * x
    dof = n - 2
    s2 = float(residual @ residual) / dof
    if s2 <= 1e-30:
        return 0.0 if abs(beta) <= 1e-15 else float("inf")
    return abs(beta) / math.sqrt(s2 / xx)


def _diagnose(values: np.ndarray, start: int, window: int, z: float) -> dict:
    tail = values[start:]
    final = values[-window:]
    mean_tail = float(tail.mean())
    mean_final = float(final.mean())
    scale_tail = float(tail.std(ddof=1)) if len(tail) > 1 else 0.0
    scale_final = float(final.std(ddof=1)) if len(final) > 1 else 0.0
    se = math.hypot(
        scale_tail / math.sqrt(max(1, len(tail))),
        scale_final / math.sqrt(max(1, len(final))),
    )
    delta = abs(mean_tail - mean_final)
    # A small scale floor avoids declaring a numerically constant channel
    # unstable because of roundoff in the reweighted ratio.
    scale_floor = max(1e-12, 1e-6 * max(1.0, abs(mean_final)))
    compatible = delta <= z * se + scale_floor
    slope_t = _slope_t(tail)
    slope_ok = math.isfinite(slope_t) and slope_t <= z
    return {
        "start_block": int(start),
        "tail_blocks": int(len(tail)),
        "tail_mean": mean_tail,
        "final_window_mean": mean_final,
        "difference": delta,
        "difference_se": se,
        "mean_compatible": bool(compatible),
        "slope_t": slope_t if math.isfinite(slope_t) else None,
        "slope_compatible": bool(slope_ok),
        "stable": bool(compatible and slope_ok),
    }


def choose_cut(values: np.ndarray, names: list[str], selected: Iterable[str] | None,
               window: int | None, min_tail: int | None, z: float) -> tuple[int | None, dict]:
    n = values.shape[0]
    if n < 4:
        raise ValueError("at least four blocks are needed for a stationarity screen")
    if selected:
        unknown = sorted(set(selected) - set(names))
        if unknown:
            raise ValueError("unknown observables: " + ", ".join(unknown))
        indexes = [names.index(name) for name in selected]
    else:
        # Default diagnostics focus on the quantities that track the sampled
        # phonon state and the full-grid z structure factor. In particular,
        # require every Szz(q) channel and the nonlinear all-q maximum when
        # they are available; do not gate on every cross-channel observable.
        preferred = [name for name in names if name in {"X_x", "X_y", "X_z", "majorana_energy"}
                     or re.fullmatch(r"S_q\d+_zz", name)
                     or name == "all_q_max_mz2_block"]
        indexes = [names.index(name) for name in preferred] or list(range(len(names)))
    if window is None:
        window = max(4, min(32, n // 8))
    if min_tail is None:
        min_tail = max(2 * window, 8)
    if window < 2 or min_tail < window or min_tail > n:
        raise ValueError("require 2 <= window <= min_tail <= number of blocks")

    details: dict[str, list[dict]] = {}
    per_observable: dict[str, int | None] = {}
    candidates = range(0, n - min_tail + 1)
    for index in indexes:
        checks = [_diagnose(values[:, index], start, window, z) for start in candidates]
        details[names[index]] = checks
        first = next((d["start_block"] for d in checks if d["stable"]), None)
        per_observable[names[index]] = first
    available = [cut for cut in per_observable.values() if cut is not None]
    cut = max(available) if len(available) == len(per_observable) else None
    return cut, {
        "window_blocks": window,
        "minimum_tail_blocks": min_tail,
        "z_threshold": z,
        "observables": list(per_observable),
        "per_observable_start_block": per_observable,
        "diagnostics": details,
        "recommended_start_block": cut,
    }


def _summary(rows, names, phase, nums, counts, cut):
    p = phase[cut:].sum()
    n = nums[cut:].sum(axis=0)
    total_count = int(counts[cut:].sum())
    if abs(p) <= 1e-12 * total_count:
        raise ValueError("retained suffix has an unresolved total phase denominator")
    value = n / p
    # Delete-one-block jackknife of the same ratio estimator.
    replicas = []
    for i in range(cut, len(rows)):
        pp = p - phase[i]
        nn = n - nums[i]
        replicas.append(nn / pp if abs(pp) > 1e-12 * (total_count - counts[i]) else None)
    replicas = np.asarray(replicas)
    if len(replicas) >= 2 and all(replica is not None for replica in replicas):
        centered = replicas - replicas.mean(axis=0)
        stderr = np.sqrt((len(replicas) - 1) / len(replicas) * np.sum(centered.real ** 2, axis=0))
    else:
        stderr = [None] * len(names)
    return {
        "blocks": len(rows) - cut,
        "samples": int(counts[cut:].sum()),
        "phase_mean": [float((p / counts[cut:].sum()).real), float((p / counts[cut:].sum()).imag)],
        "observables": {
            name: {"value": [float(z.real), float(z.imag)],
                   "stderr_real": float(err) if err is not None and math.isfinite(float(err)) else None}
            for name, z, err in zip(names, value, stderr)
        },
    }


def process(block_path: Path, args: argparse.Namespace) -> dict:
    rows, names, phase, nums, counts, values = read_blocks(block_path)
    results_path = args.results
    if results_path is None:
        candidate = block_path.with_name(block_path.name.removesuffix("_blocks.csv") + "_results.json")
        results_path = candidate if candidate.exists() else None
    metadata = json.loads(results_path.read_text()) if results_path else {}
    if metadata and metadata.get("status") != "complete":
        raise ValueError(f"source run is not complete: {metadata.get('status')}")
    if metadata and int(metadata.get("blocks", len(rows))) != len(rows):
        raise ValueError("block CSV row count disagrees with result metadata")
    q_names = sorted(
        (name for name in names if re.fullmatch(r"S_q\d+_zz", name)),
        key=lambda name: int(re.search(r"q(\d+)", name).group(1)),
    )
    if q_names:
        nc = int(metadata.get("Lx", len(q_names))) * int(metadata.get("Ly", 1))
        if nc != len(q_names):
            raise ValueError(f"all-q Szz channel count {len(q_names)} disagrees with Lx*Ly={nc}")
        q_columns = [names.index(name) for name in q_names]
        q_block_values = (nums[:, q_columns] / phase[:, None]).real / nc
        block_max = q_block_values.max(axis=1)
        names.append("all_q_max_mz2_block")
        values = np.column_stack((values, block_max))
    else:
        nc = None
        q_block_values = None
    cut, diagnostics = choose_cut(values, names, args.observable, args.window_blocks,
                                  args.min_tail_blocks, args.z)
    if cut is not None and abs(phase[cut:].sum()) <= 1e-12 * int(counts[cut:].sum()):
        diagnostics["rejection_reason"] = "retained_suffix_phase_denominator_unresolved"
        diagnostics["recommended_start_block"] = None
        cut = None
    warmup = int(metadata.get("warmup_sweeps", args.warmup))
    stride = int(metadata.get("measure_every", args.measure_every))
    block_count = int(counts[0])
    start_sweep = None if cut is None else warmup + (cut * block_count + 1) * stride
    end_sweep = warmup + len(rows) * block_count * stride
    payload = {
        "method": "posthoc_block_stationarity_screen",
        "source_blocks": str(block_path),
        "source_results": str(results_path) if results_path else None,
        "warmup_sweeps_in_source": warmup,
        "measure_every": stride,
        "block_samples": block_count,
        "measurement_sweeps_available": len(rows) * block_count * stride,
        "measurement_sweep_start": start_sweep,
        "measurement_sweep_end": end_sweep,
        "all_q_max_definition": (
            f"per-block max_q Re(sum_phase_S_q_zz/sum_phase)/{nc}; diagnostic channel only"
            if q_names else None
        ),
        "diagnostics": diagnostics,
        "selected_summary": _summary(rows, names[:len(nums[0])], phase, nums, counts, cut)
        if cut is not None else None,
        "recommended": cut is not None,
        "caveat": "A stable block suffix is evidence of finite-run stationarity, not a proof of equilibrium. Check multiple observables and independent seeds.",
    }
    prefix = args.output or block_path.with_name(block_path.name.removesuffix("_blocks.csv") + "_equilibrated")
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_name(prefix.name + ".json")
    plot_path = prefix.with_name(prefix.name + ".png")
    selected_csv = prefix.with_name(prefix.name + "_blocks.csv")
    outputs = [json_path, plot_path, selected_csv]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not args.force:
        raise FileExistsError("output exists; choose another --output prefix or pass --force: "
                              + ", ".join(existing))
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    _plot_diagnostics(plot_path, values, names, cut, q_block_values, phase, nums, q_names, nc)
    if cut is not None:
        with selected_csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows[cut:])
    elif selected_csv.exists() and args.force:
        selected_csv.unlink()
    print(json.dumps({"recommended_start_block": cut, "measurement_sweep_start": start_sweep,
                      "selected_blocks": None if cut is None else len(rows) - cut,
                      "json": str(json_path), "plot": str(plot_path)}, indent=2))
    return payload


def _plot_diagnostics(path: Path, values: np.ndarray, names: list[str], cut: int | None,
                      q_block_values: np.ndarray | None, phase: np.ndarray,
                      nums: np.ndarray, q_names: list[str], nc: int | None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(values)
    x = np.arange(1, n + 1)
    fig, axes = plt.subplots(2 if q_block_values is not None else 1, 1,
                             figsize=(11, 8 if q_block_values is not None else 5),
                             squeeze=False, layout="constrained")
    ax = axes[0, 0]
    if q_block_values is not None:
        cumulative = []
        for stop in range(1, n + 1):
            denominator = phase[:stop].sum()
            if abs(denominator) <= 1e-12 * countsafe(phase[:stop]):
                cumulative.append(np.nan)
            else:
                numerator = nums[:stop, [names.index(name) for name in q_names]].sum(axis=0)
                cumulative.append(float(np.max((numerator / denominator).real / nc)))
        ax.plot(x, values[:, names.index("all_q_max_mz2_block")], ".", alpha=.45,
                label="block all-q max")
        ax.plot(x, cumulative, "-", linewidth=1.8, label="cumulative all-q max")
        ax.set_ylabel(r"$\max_q\,\mathrm{Re}[S_{zz}(q)]/N_c$")
        lower = axes[1, 0]
        for q in range(q_block_values.shape[1]):
            lower.plot(x, q_block_values[:, q], alpha=.45, linewidth=.8,
                       label=q_names[q] if q < 8 else None)
        lower.set_ylabel("block Szz(q) / Nc")
        lower.set_xlabel("measurement block")
        lower.set_title("Individual momentum channels (first eight shown in legend)")
        if q_block_values.shape[1] <= 8:
            lower.legend(ncol=4, fontsize=8)
        lower.grid(alpha=.2)
    else:
        for i, name in enumerate(names[:min(6, len(names))]):
            ax.plot(x, values[:, i], ".-", label=name)
        ax.set_ylabel("reweighted block mean")
    if cut is not None:
        for target in (ax, axes[1, 0] if q_block_values is not None else ax):
            target.axvline(cut + .5, color="black", linestyle="--", linewidth=1.2,
                           label=f"selected start: block {cut}")
    ax.set_xlabel("measurement block")
    ax.set_title("Post-run stationarity diagnostic; visual and statistical screens")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def countsafe(phase: np.ndarray) -> float:
    """Scale factor for detecting a near-zero average phase in a running sum."""
    return float(len(phase))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blocks", type=Path, nargs="+")
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path, help="output prefix for one input block file")
    parser.add_argument("--observable", action="append", help="diagnostic observable; repeatable")
    parser.add_argument("--window-blocks", type=int)
    parser.add_argument("--min-tail-blocks", type=int)
    parser.add_argument("--z", type=float, default=2.0)
    parser.add_argument("--warmup", type=int, default=0)
    parser.add_argument("--measure-every", type=int, default=1)
    parser.add_argument("--force", action="store_true", help="replace this analysis output prefix")
    args = parser.parse_args()
    if args.z <= 0 or args.warmup < 0 or args.measure_every <= 0:
        parser.error("z must be positive; warmup must be nonnegative; measure-every must be positive")
    if len(args.blocks) > 1 and args.output:
        parser.error("--output is only supported for one input")
    for path in args.blocks:
        process(path, args)


if __name__ == "__main__":
    main()
