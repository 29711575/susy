#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
SHIFTED_RUN_DIR=${SHIFTED_RUN_DIR:-build-kitaev-p01/shifted_trial}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group U16=build-kitaev-p01/p5_runs/t16_a_results.json,build-kitaev-p01/p5_runs/t16_b_results.json,build-kitaev-p01/p5_runs/t16_c_results.json \
 --group U32=build-kitaev-p01/p5_runs/t32_a_results.json,build-kitaev-p01/p5_runs/t32_b_results.json,build-kitaev-p01/p5_runs/t32_c_results.json \
 --group "D16=$SHIFTED_RUN_DIR/results/shifted_t16_a_results.json,$SHIFTED_RUN_DIR/results/shifted_t16_b_results.json" \
 --group "D32=$SHIFTED_RUN_DIR/results/shifted_t32_a_results.json,$SHIFTED_RUN_DIR/results/shifted_t32_b_results.json" \
 --block-sweeps 100 --looks 1 \
 --compare U16:U32 --compare D16:D32 --compare U16:D16 --compare U32:D32 \
 --output "$SHIFTED_RUN_DIR/analysis"
