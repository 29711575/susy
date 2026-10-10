#!/bin/bash
# Run in release root, using existing P5 raw chains and the downloaded pilot.
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
P5_FOLLOWUP_DIR=${P5_FOLLOWUP_DIR:-build-kitaev-p01/p5_followup}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
  --group T16_h02=build-kitaev-p01/p5_runs/t16_a_results.json,build-kitaev-p01/p5_runs/t16_b_results.json,build-kitaev-p01/p5_runs/t16_c_results.json \
  --group T32_h02=build-kitaev-p01/p5_runs/t32_a_results.json,build-kitaev-p01/p5_runs/t32_b_results.json,build-kitaev-p01/p5_runs/t32_c_results.json \
  --group "T64_h02=$P5_FOLLOWUP_DIR/results/t64_a_results.json,$P5_FOLLOWUP_DIR/results/t64_b_results.json" \
  --block-sweeps 100 --looks 1 --compare T32_h02:T64_h02 \
  --output "$P5_FOLLOWUP_DIR/analysis"
"$P5_PYTHON" kitaev_ssh/plot_projection_response.py \
  "$P5_FOLLOWUP_DIR/analysis" "$P5_FOLLOWUP_DIR/results/free_diagnostic.json"
