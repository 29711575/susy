#!/usr/bin/env bash
# Build only the focused correctness targets, then preserve all run evidence.
set -euo pipefail
TASK_ROOT=$(cd -- "$(dirname -- "$0")/.." && pwd)
TASK_BUILD=${1:-"$TASK_ROOT/build-pfqmc-correctness"}
TASK_OUTPUT=${2:-"$TASK_ROOT/audit/pfqmc-check-$(date -u +%Y%m%dT%H%M%SZ)"}
: "${MKLROOT:?Set MKLROOT as for the production build}"
export LD_LIBRARY_PATH="$MKLROOT/lib:$MKLROOT/lib/intel64:${LD_LIBRARY_PATH:-}"
bash "$TASK_ROOT/scripts/build_nodeps.sh" "$TASK_BUILD" \
    kitaev_ssh_total_energy_validate \
    kitaev_ssh_pfqmc_free_phonon_validate \
    kitaev_ssh_pfqmc_fock_validate \
    kitaev_ssh_pfqmc_cache_validate \
    kitaev_ssh_phonon_validate \
    kitaev_ssh_mc_validate \
    kitaev_ssh_vbs_susceptibility_validate \
    kitaev_ssh_numerical_failure_validate \
    kitaev_ssh_green_refinement_validate \
    kitaev_ssh_fermion_validate \
    kitaev_ssh_real_magnitude_validate \
    kitaev_ssh_physical_trial_validate \
    kitaev_ssh_bipartite_validate
python3 "$TASK_ROOT/scripts/run_pfqmc_correctness.py" \
    --build-dir "$TASK_BUILD" --output-dir "$TASK_OUTPUT"
