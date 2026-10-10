#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
STAT_RUN_DIR=${STAT_RUN_DIR:-build-kitaev-p01/shifted_t32_stationarity}
STAT_ANALYSIS_DIR=${STAT_ANALYSIS_DIR:-$STAT_RUN_DIR/analysis}
STAT_OLD_ANALYSIS_DIR=${STAT_OLD_ANALYSIS_DIR:-$STAT_RUN_DIR/old_seed7002_diagnostic}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group "NEW=$STAT_RUN_DIR/results/t32_stationary_a_results.json,$STAT_RUN_DIR/results/t32_stationary_b_results.json,$STAT_RUN_DIR/results/t32_stationary_c_results.json,$STAT_RUN_DIR/results/t32_stationary_d_results.json" \
 --block-sweeps 200 --looks 1 --output "$STAT_ANALYSIS_DIR"
# Prior seed 7002 is a *separate* diagnostic; do not pool its shorter-warmup run
# into the new-chain stationary group or alter the primary block sensitivity.
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group OLD7002=build-kitaev-p01/shifted_long/results/long_t32_b_results.json \
 --block-sweeps 200 --looks 1 --output "$STAT_OLD_ANALYSIS_DIR"
