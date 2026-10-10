#!/bin/bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
python=${P5_PYTHON:-/home/asus/.cache/pfqmc_plot_env/bin/python}
data=${L4_MIXING_DIR:-$root/build-kitaev-p01/l4_mixing}
output=${L4_MIXING_ANALYSIS_DIR:-$data/analysis}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-l4-mixing-matplotlib}
files=
for seed in 7801 7802 7803 7804; do
 file="$data/results/l4_${seed}_results.json"
 test -f "$file"
 if test -n "$files"; then files="$files,$file"; else files="$file"; fi
done
"$python" "$root/kitaev_ssh/convergence_analysis.py" \
 --group "L4_mixing=$files" --block-sweeps 100 --output "$output"
