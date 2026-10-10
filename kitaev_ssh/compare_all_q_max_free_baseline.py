#!/usr/bin/env python3
"""Compare the lambda=0 all-q max with its exact free-phonon baseline."""
from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/sign-size-lambda-L4-tau-3seed-20260918/plots/mz2_all_q_max_by_lambda_tau.csv"
OUTDIR = ROOT / "results/sign-size-lambda-L4-tau-3seed-20260918/plots"
OUTPNG = OUTDIR / "mz2_all_q_max_lambda0_free_ED_baseline.png"
OUTCSV = OUTDIR / "mz2_all_q_max_lambda0_free_ED_baseline.csv"

# These are the parameters recorded in the corresponding MC result files.
MASS = 4.0
OMEGA = 0.5
NC = 16


def main() -> None:
    rows = []
    with INPUT.open(newline="") as f:
        for row in csv.DictReader(f):
            if abs(float(row["lambda"])) < 1e-12:
                rows.append({k: float(row[k]) for k in ("tau", "all_q_max_mz2", "stderr")})
    rows.sort(key=lambda r: r["tau"])
    if not rows:
        raise SystemExit("no lambda=0 rows found")

    # For independent harmonic bond coordinates, S_zz(q)=<|F_q|^2>/Nc=1/(2mω)
    # and m_z^2=S_zz/Nc. This is the exact ED/analytic free-phonon value.
    baseline = 1.0 / (2.0 * MASS * OMEGA * NC)
    with OUTCSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tau", "qmc_all_q_max", "qmc_stderr", "free_exact_mz2"])
        for r in rows:
            w.writerow([r["tau"], r["all_q_max_mz2"], r["stderr"], baseline])

    tau = [r["tau"] for r in rows]
    val = [r["all_q_max_mz2"] for r in rows]
    err = [r["stderr"] for r in rows]
    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=180)
    ax.errorbar(tau, val, yerr=err, fmt="o-", capsize=3, lw=1.6,
                label="QMC: λ=0, all-q max")
    ax.axhline(baseline, color="black", ls="--", lw=1.5,
               label=fr"free exact baseline $1/(2m\omega N_c)={baseline:.6f}$")
    ax.set_xlabel(r"imaginary time $\tau$")
    ax.set_ylabel(r"$\max_q\,m_z^2(q)$")
    ax.set_title("L=4, λ=0: all-q maximum versus exact free-phonon baseline")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPNG)
    plt.close(fig)
    print(OUTPNG)
    print(OUTCSV)
    print(f"baseline={baseline:.12g}")


if __name__ == "__main__":
    main()
