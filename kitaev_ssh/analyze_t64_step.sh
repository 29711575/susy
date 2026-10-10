#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
STEP_DIR=${STEP_DIR:-build-kitaev-p01/t64_step}
STEP_ANALYSIS_DIR=${STEP_ANALYSIS_DIR:-$STEP_DIR/analysis}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
coarse=()
for x in a b c d; do coarse+=("build-kitaev-p01/t32_t64_projection/results/t64_${x}_results.json"); done
fine=()
for x in a b c d; do fine+=("$STEP_DIR/results/h01_${x}_results.json"); done
coarse_joined=$(IFS=,; echo "${coarse[*]}")
fine_joined=$(IFS=,; echo "${fine[*]}")
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group "H02=$coarse_joined" --group "H01=$fine_joined" \
 --block-sweeps 200 --looks 1 --compare H02:H01 --step-pair H02:H01 --output "$STEP_ANALYSIS_DIR"
