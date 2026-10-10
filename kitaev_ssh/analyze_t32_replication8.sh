#!/bin/bash
set -euo pipefail
P5_PYTHON=${P5_PYTHON:-python3}
R8_DIR=${R8_DIR:-build-kitaev-p01/t32_replication8}
R8_AUX_DIR=${R8_AUX_DIR:-$R8_DIR/auxiliary}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-p5-mpl}
files=()
for x in a b c d e f g h; do files+=("$R8_DIR/results/replication8_${x}_results.json"); done
joined=$(IFS=,; echo "${files[*]}")
"$P5_PYTHON" kitaev_ssh/convergence_analysis.py --group "R8=$joined" --block-sweeps 200 --looks 1 --output "$R8_AUX_DIR"
"$P5_PYTHON" kitaev_ssh/paired_t32_replication8.py --run-dir "$R8_DIR"
