#!/usr/bin/env python3
"""Plot the published compact static tables; no optimization or QMC."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambda-min", type=float, default=0.1)
    parser.add_argument("--lambda-max", type=float, default=3.5)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    if args.lambda_min > args.lambda_max:
        parser.error("lambda-min must not exceed lambda-max")
    table = Path(__file__).resolve().parent / "results" / "summary.csv"
    with table.open() as handle:
        rows = [r for r in csv.DictReader(handle)
                if args.lambda_min <= float(r["lambda_value"]) <= args.lambda_max]
    if not rows:
        parser.error("No published data in the requested range")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    styles = {6: ("#2878b5", "o"), 7: ("#209e50", "^"),
              8: ("#d53a32", "D"), 9: ("#cc7900", "s")}

    def curves(ax, key):
        for L, (color, marker) in styles.items():
            group = sorted((r for r in rows if int(r["L"]) == L),
                           key=lambda r: float(r["lambda_value"]))
            if not group:
                continue
            values = [float(r[key]) for r in group]
            if key == "m2_z":
                values = [v if v >= 1e-20 else 0.0 for v in values]
            ax.plot([float(r["lambda_value"]) for r in group], values,
                    marker=marker, color=color, linewidth=0.9, markersize=3,
                    markerfacecolor="none" if L in (8, 9) else color,
                    label=f"L={L} ({len(group)} points)")
        ax.set_xlabel(r"$\lambda=g^2/k$")
        ax.set_xlim(args.lambda_min, args.lambda_max)
        ax.grid(alpha=0.2)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), layout="constrained")
    for ax, key, ylabel, title in zip(
            axes, ["m2_z", "order_vector_max_nonzero", "uniform_xy_difference"],
            [r"$\max_{q\ne0}|D_z(q)|^2$", r"$\max_{q\ne0}\|D(q)\|_2$",
             r"$|\overline{X}_x-\overline{X}_y|$"],
            ["z-bond order squared", "Three-bond modulation", "Uniform x/y difference"]):
        curves(ax, key)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(bottom=0)
    axes[0].legend(fontsize=8)
    fig.suptitle("Strict static limit: physical-parity finite tori")
    fig.savefig(args.output_dir / "order-parameters.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.7), layout="constrained")
    curves(ax, "chi_zz_max_nonzero")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\max_{q\ne0}\chi_{X,zz}(q)$")
    ax.set_title("Static single-domain relaxed susceptibility")
    ax.legend(fontsize=9)
    fig.savefig(args.output_dir / "susceptibility-log.png", dpi=180)
    plt.close(fig)
    print(f"Plotted {len(rows)} published points into {args.output_dir}")


if __name__ == "__main__":
    main()
