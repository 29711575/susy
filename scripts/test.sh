#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
build=$(cd "${1:-$root/build}" && pwd)
work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
cd "$work"
"$build/density_sign_regression" "$root/reproduction/pure_state_driven_stage1/density_sign_fixtures"
"$build/cdw_initial_test"
for mode in negative positive eta unknown anomalies; do "$build/hs1_review_test" "$mode"; done
for state in free cdw-even; do
 "$build/pure_state_driven_qmc" --L 18 --boundary pbc --hs 1 --dt 0.01 --parity -1 \
 --burn 1 --measurements 2 --block 16 --seed 98150915 --protocol quench --Vi 0 --Vf 4 --tau 0.1 \
 --initialization-policy mirrored-theorem --update-policy stabilized-local \
 --direct-green-columns true --compact-density-proposals true --reuse-rebuild-workspace false \
 --initial-state "$state" --output "$work/$state.csv"
done
