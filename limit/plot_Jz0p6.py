#!/usr/bin/env python3
"""Plot the published Jz=0.6 tables offline; never starts a calculation."""
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    default = Path(__file__).resolve().parent / "results/Jz0p6-20261110"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=default)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output = args.output_dir or args.data_dir
    output.mkdir(parents=True, exist_ok=True)
    with (args.data_dir / "summary.csv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    with (args.data_dir / "coverage.csv").open(encoding="utf-8", newline="") as f:
        coverage = list(csv.DictReader(f))
    keys = {(int(r["L"]), round(float(r["lambda_value"]), 8)) for r in rows}
    assert len(keys) == len(rows), "Duplicate parameter points"
    expected = {(int(r["L"]), round(float(r["lambda_value"]), 8))
                for r in coverage if r["status"] == "validated"}
    assert keys == expected, "Coverage and data disagree"
    for r in rows:
        assert r["validated"] == "True" and float(r["Jz"]) == 0.6
        assert all(math.isfinite(float(r[k])) for k in ("m2_z", "chi_zz_max_nonzero"))
        assert float(r["chi_zz_max_nonzero"]) > 0

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    palette = ["#2673ad", "#d37826", "#8c4267", "#777e36"]
    markers = ["o", "s", "^", "D", "v", "P", "X", "h"]
    groups = {L: sorted((r for r in rows if int(r["L"]) == L),
                        key=lambda r: float(r["lambda_value"]))
              for L in range(3, 11)}
    metrics = [("m2_z", r"$m_z^2$"), ("chi_zz_max_nonzero", r"$\chi_{zz}$")]
    missing = [(int(r["L"]), float(r["lambda_value"]))
               for r in coverage if r["status"] != "validated"]
    note = (f"Validated: {len(rows)}/{len(coverage)} points. Missing: "
            + (", ".join(f"L={L}, lambda={lam:g}" for L, lam in missing) or "none"))

    def values(rr, key):
        raw = [float(r[key]) for r in rr]
        return [0.0 if key == "m2_z" and abs(y) < 1e-20 else y for y in raw]

    fig, axes = plt.subplots(4, 4, figsize=(19, 13), sharex=True)
    for index, L in enumerate(range(3, 11)):
        rr = groups[L]
        x = [float(r["lambda_value"]) for r in rr]
        for column, (key, label) in enumerate(metrics):
            ax = axes[index // 2, 2 * (index % 2) + column]
            y = values(rr, key)
            ax.plot(x, y, "o-", color=palette[index % 4], ms=3.5,
                    mfc="white", lw=1.2)
            ax.set_title(f"L={L}: {label} ({len(rr)}/30 points)", fontsize=12)
            ax.set_ylabel(label)
            ax.set_xlim(.04, 3.06)
            ax.set_xticks([.1, 1, 2, 3])
            ax.grid(alpha=.18)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 4))
            if key == "m2_z":
                if max(y) == 0:
                    ax.set_ylim(-1e-4, 1e-4)
                    ax.text(.04, .87, "Numerical zero", transform=ax.transAxes)
                else:
                    ax.set_ylim(bottom=-.04 * max(y))
            for size, lam in missing:
                if size == L:
                    ax.axvspan(lam - .05, lam + .05, color="#dddddd", alpha=.6, zorder=0)
            if index // 2 == 3:
                ax.set_xlabel(r"$\lambda$")
    fig.suptitle(r"Adiabatic limit: $J=(1,1,0.6)$, $k=1$, $L=3$--$10$", fontsize=18)
    fig.text(.5, .94, note + ". Independent y axes; shaded areas are uncomputed.",
             ha="center", fontsize=11)
    fig.subplots_adjust(top=.89, bottom=.11, left=.06, right=.985, hspace=.46, wspace=.30)
    fig.text(.06, .05, "Nonzero-q z displacement channel; stable single-domain susceptibility; fixed Wilson(+,+), physical matter parity.")
    fig.text(.06, .025, "Each point is marked; m2 < 1e-20 is displayed as zero only. Lines guide the eye; no interpolation or global-minimum proof.", fontsize=10)
    fig.savefig(output / "order-and-susceptibility.png", dpi=200, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for index, L in enumerate(range(3, 11)):
        rr = groups[L]
        x = [float(r["lambda_value"]) for r in rr]
        for ax, (key, _) in zip(axes, metrics):
            ax.plot(x, values(rr, key), color=palette[index % 4], marker=markers[index],
                    ms=4, mfc="white", lw=1.15, ls="-" if index < 4 else "--", label=f"L={L}")
    for ax, (_, label) in zip(axes, metrics):
        ax.set_xlabel(r"$\lambda$")
        ax.set_ylabel(label)
        ax.set_xlim(.04, 3.06)
        ax.grid(alpha=.2)
    axes[0].set_ylim(bottom=-.0008)
    axes[0].legend(ncol=2, frameon=False, fontsize=10)
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$\chi_{zz}$ (log scale)")
    fig.suptitle(r"Adiabatic size comparison: $J=(1,1,0.6)$, $k=1$", fontsize=16)
    fig.subplots_adjust(top=.88, bottom=.2, left=.08, right=.98, wspace=.25)
    fig.text(.08, .09, note, fontsize=10)
    fig.text(.08, .04, "Nonzero-q z channel; single-domain response. Logarithmic chi axis. No extrapolation across missing points.", fontsize=10)
    fig.savefig(output / "size-comparison.png", dpi=200, facecolor="white")
    plt.close(fig)
    print(f"Plotted {len(rows)} validated rows in {output}")


if __name__ == "__main__":
    main()
