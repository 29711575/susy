#!/usr/bin/env bash
# From the release directory, after all explicitly selected P5 chains complete.
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
P5_ANALYSIS_DIR=${P5_ANALYSIS_DIR:-build-kitaev-p01/p5_runs/final_analysis}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
mkdir -p "$MPLCONFIGDIR"
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
  --group T08_h02=build-kitaev-p01/p5_runs/t08_a_results.json,build-kitaev-p01/p5_runs/t08_b_results.json \
  --group T16_h02=build-kitaev-p01/p5_runs/t16_a_results.json,build-kitaev-p01/p5_runs/t16_b_results.json,build-kitaev-p01/p5_runs/t16_c_results.json \
  --group T32_h02=build-kitaev-p01/p5_runs/t32_a_results.json,build-kitaev-p01/p5_runs/t32_b_results.json,build-kitaev-p01/p5_runs/t32_c_results.json \
  --group T32_h01=build-kitaev-p01/p5_runs/t32_h01_a_results.json,build-kitaev-p01/p5_runs/t32_h01_b_results.json \
  --block-sweeps 100 --looks 2 \
  --compare T08_h02:T16_h02 --compare T16_h02:T32_h02 --compare T32_h02:T32_h01 \
  --step-pair T32_h02:T32_h01 --output "$P5_ANALYSIS_DIR"
