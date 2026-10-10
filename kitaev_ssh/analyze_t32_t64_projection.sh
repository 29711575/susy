#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
PROJECTION_DIR=${PROJECTION_DIR:-build-kitaev-p01/t32_t64_projection}
PROJECTION_ANALYSIS_DIR=${PROJECTION_ANALYSIS_DIR:-$PROJECTION_DIR/analysis}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
t32=()
for x in a b c d e f g h; do t32+=("build-kitaev-p01/t32_replication8/results/replication8_${x}_results.json"); done
t64=()
for x in a b c d; do t64+=("$PROJECTION_DIR/results/t64_${x}_results.json"); done
t32_joined=$(IFS=,; echo "${t32[*]}")
t64_joined=$(IFS=,; echo "${t64[*]}")
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group "T32=$t32_joined" --group "T64=$t64_joined" \
 --block-sweeps 200 --looks 1 --compare T32:T64 --output "$PROJECTION_ANALYSIS_DIR"
