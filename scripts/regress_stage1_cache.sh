#!/usr/bin/env bash
# Fixed-seed regression for the stage-1 cache/dense/workspace changes.
# Usage: regress_stage1_cache.sh <current_build> <baseline_build> <new_work_dir>
set -euo pipefail

current_build=$(cd -- "${1:?current build directory}" && pwd)
baseline_build=$(cd -- "${2:?baseline build directory}" && pwd)
work_input=${3:?new work directory}
if [ -e "$work_input" ]; then
    echo "refusing to reuse existing regression directory: $work_input" >&2
    exit 2
fi
mkdir -p "$work_input"
work=$(cd "$work_input" && pwd)
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
current_bin="$current_build/kitaev_ssh_mc"
baseline_bin="$baseline_build/kitaev_ssh_mc"
[ -x "$current_bin" ] || { echo "missing $current_bin" >&2; exit 2; }
[ -x "$baseline_bin" ] || { echo "missing $baseline_bin" >&2; exit 2; }

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export LD_LIBRARY_PATH=${MKLROOT:-/opt/intel/oneapi/mkl/2026.1}/lib:${LD_LIBRARY_PATH:-}

sha256sum "$current_bin" "$baseline_bin" > "$work/binary_sha256.txt"
printf '%s\n' \
  "current_build=$current_build" \
  "baseline_build=$baseline_build" \
  "compiler_flags=-O3 -fno-fast-math -ffp-contract=off -DPFQMC_SCALE_SAFE_UDT" \
  "OMP_NUM_THREADS=$OMP_NUM_THREADS" \
  "MKL_NUM_THREADS=$MKL_NUM_THREADS" \
  "MKL_DYNAMIC=$MKL_DYNAMIC" > "$work/environment.txt"

run_case() {
    local label=$1 lx=$2 ly=$3 intervals=$4 dt=$5 coupling=$6 seed=$7 warmup=$8 sweeps=$9
    local stride=${10} blocks=${11} proposal=${12} defer=${13}
    local common=(--Lx "$lx" --Ly "$ly" --Jx 1 --Jy 1 --Jz 0.5
        --intervals "$intervals" --dt "$dt" --mass 1.7 --omega-ph 0.8
        --g "$coupling" --seed "$seed" --warmup "$warmup" --sweeps "$sweeps"
        --measure-every "$stride" --block-samples "$blocks" --qr-block 4
        --weight-engine combined --vbs 0 --proposal "$proposal" --defer-phase "$defer")
    local tag="${label}_${proposal//-/_}_defer${defer}"
    for kind in baseline current; do
        local bin=$baseline_bin
        [ "$kind" = current ] && bin=$current_bin
        /usr/bin/time -v -o "$work/${kind}_${tag}.time" \
            "$bin" "${common[@]}" --output "$work/${kind}_${tag}" \
            > "$work/${kind}_${tag}.stdout" 2> "$work/${kind}_${tag}.stderr"
    done
    python3 "$root/scripts/compare_stage1_outputs.py" \
        --baseline "$work/baseline_${tag}_results.json" \
        --current "$work/current_${tag}_results.json" \
        --baseline-blocks "$work/baseline_${tag}_blocks.csv" \
        --current-blocks "$work/current_${tag}_blocks.csv" \
        --output "$work/compare_${tag}.json"
    python3 - "$work/current_${tag}_results.json" "$proposal" "$intervals" <<'PY'
import json
import sys

data = json.load(open(sys.argv[1]))
proposal = sys.argv[2]
intervals = int(sys.argv[3])
counts = data["cache_diagnostics"]
proposals = data["proposals_warmup"] + data["proposals_measurement"]
if proposal == "local":
    assert counts["local_proposals"] == proposals
    assert counts["bond_refresh_full_rebuilds"] == 0
    assert counts["local_node_rebuilds"] <= proposals + intervals + 1
    assert counts["local_node_reuses"] >= max(0, proposals - 1) * intervals
else:
    assert counts["local_proposals"] == 0
    assert counts["bond_refresh_full_rebuilds"] == proposals
print(f"cache diagnostics PASS: {proposal} proposals={proposals} {counts}")
PY
}

# Base case, zero coupling/small projection, and nonzero larger/longer projection.
for proposal in local bond-refresh; do
    for defer in 0 1; do
        run_case base 3 2 4 0.2 0.3 314159 8 32 2 4 "$proposal" "$defer"
        run_case zero 2 2 2 0.2 0.0 9101 2 8 2 2 "$proposal" "$defer"
        run_case long 4 3 8 0.1 0.6 9102 2 4 2 2 "$proposal" "$defer"
    done
done

printf '%s\n' "REGRESS_STAGE1_CACHE PASS" | tee "$work/status.txt"
