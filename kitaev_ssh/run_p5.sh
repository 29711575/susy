#!/usr/bin/env bash
# Run from the release directory, after sourcing the existing oneAPI environment.
# Distinct output prefixes and seeds preserve chain boundaries; no resume/overwrite.
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_DYNAMIC=FALSE
P5_OUTPUT_DIR=${P5_OUTPUT_DIR:-build-kitaev-p01/p5_runs}
P5_BINARY=${P5_BINARY:-build-kitaev-p01/kitaev_ssh_mc}
mkdir -p "$P5_OUTPUT_DIR"
run_chain() {
    local label=$1 seed=$2 intervals=$3 step=$4
    local sweeps=${5:-8000}
    local trial_args=()
    if [[ $# -ge 6 ]]; then trial_args=(--trial-dx "$6" --trial-dy "$6" --trial-dz "$6"); fi
    local suffix
    for suffix in results.json blocks.csv stdout.log stderr.log; do
        if [[ -e "$P5_OUTPUT_DIR/${label}_$suffix" ]]; then
            echo "Refusing to overwrite $P5_OUTPUT_DIR/${label}_$suffix; choose a new P5_OUTPUT_DIR." >&2
            return 2
        fi
    done
    "$P5_BINARY" \
        --vbs 1 --proposal bond-refresh --Lx 2 --Ly 2 \
        --Jx 1 --Jy 1 --Jz 0.5 --g 0.3 --mass 1.7 --omega-ph 0.8 \
        --intervals "$intervals" --dt "$step" --qr-block 4 --seed "$seed" \
        --warmup 400 --sweeps "$sweeps" --measure-every 1 --block-samples 25 \
        "${trial_args[@]}" \
        --output "$P5_OUTPUT_DIR/$label" \
        > "$P5_OUTPUT_DIR/${label}_stdout.log" \
        2> "$P5_OUTPUT_DIR/${label}_stderr.log"
}
case "${1:-}" in
    short)
        run_chain t08_a 6101 4 0.2
        run_chain t08_b 6102 4 0.2
        ;;
    long)
        run_chain t16_a 6201 8 0.2
        run_chain t16_b 6202 8 0.2
        ;;
    projection_a) run_chain t32_a 6301 16 0.2 ;;
    projection_b) run_chain t32_b 6302 16 0.2 ;;
    supplement)
        run_chain t16_c 6203 8 0.2
        run_chain t32_c 6303 16 0.2
        ;;
    fine_a) run_chain t32_h01_a 6401 32 0.1 ;;
    fine_b) run_chain t32_h01_b 6402 32 0.1 ;;
    followup_a) run_chain t64_a 6601 32 0.2 1600 ;;
    followup_b) run_chain t64_b 6602 32 0.2 1600 ;;
    shifted_t16_a) run_chain shifted_t16_a 6701 8 0.2 1600 0.1378676470588235 ;;
    shifted_t16_b) run_chain shifted_t16_b 6702 8 0.2 1600 0.1378676470588235 ;;
    shifted_t32_a) run_chain shifted_t32_a 6801 16 0.2 1600 0.1378676470588235 ;;
    shifted_t32_b) run_chain shifted_t32_b 6802 16 0.2 1600 0.1378676470588235 ;;
    *) echo 'Usage: bash kitaev_ssh/run_p5.sh short|long|projection_a|projection_b|supplement|fine_a|fine_b|followup_a|followup_b|shifted_t16_a|shifted_t16_b|shifted_t32_a|shifted_t32_b' >&2; exit 2 ;;
esac
