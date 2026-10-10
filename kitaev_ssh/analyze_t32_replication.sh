#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
REPLICATION_RUN_DIR=${REPLICATION_RUN_DIR:-build-kitaev-p01/t32_replication}
REPLICATION_AUX_DIR=${REPLICATION_AUX_DIR:-$REPLICATION_RUN_DIR/auxiliary}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py \
 --group "REPLICATION=$REPLICATION_RUN_DIR/results/replication_a_results.json,$REPLICATION_RUN_DIR/results/replication_b_results.json,$REPLICATION_RUN_DIR/results/replication_c_results.json,$REPLICATION_RUN_DIR/results/replication_d_results.json" \
 --block-sweeps 200 --looks 1 --output "$REPLICATION_AUX_DIR"
"$P5_PYTHON" kitaev_ssh/paired_t32_replication.py --run-dir "$REPLICATION_RUN_DIR"
