#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
LONG_RUN_DIR=${LONG_RUN_DIR:-build-kitaev-p01/shifted_long}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
long_groups=(
 --group "L16=$LONG_RUN_DIR/results/long_t16_a_results.json,$LONG_RUN_DIR/results/long_t16_b_results.json"
 --group "L32=$LONG_RUN_DIR/results/long_t32_a_results.json,$LONG_RUN_DIR/results/long_t32_b_results.json"
)
# Primary: only long chains, including 20 blocks/chain at 400 sweeps.
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py "${long_groups[@]}" \
 --block-sweeps 200 --looks 1 --compare L16:L32 --output "$LONG_RUN_DIR/analysis"
# Supplemental: do NOT combine long and short chains into a single group.
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py "${long_groups[@]}" \
 --group S16=build-kitaev-p01/shifted_trial/results/shifted_t16_a_results.json,build-kitaev-p01/shifted_trial/results/shifted_t16_b_results.json \
 --group S32=build-kitaev-p01/shifted_trial/results/shifted_t32_a_results.json,build-kitaev-p01/shifted_trial/results/shifted_t32_b_results.json \
 --block-sweeps 200 --group-block S16=100 --group-block S32=100 --looks 1 \
 --compare S16:L16 --compare S32:L32 --compare S16:S32 \
 --output "$LONG_RUN_DIR/history_comparison"
