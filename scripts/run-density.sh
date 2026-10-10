#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
# Pass all physical parameters explicitly; block remains the upstream default 16.
exec "${QMC_BINARY:-$root/build/pure_state_driven_qmc}" --hs 1 --initialization-policy mirrored-theorem \
 --update-policy stabilized-local --direct-green-columns true --compact-density-proposals true \
 --reuse-rebuild-workspace true "$@"
