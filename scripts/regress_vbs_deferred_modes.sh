#!/usr/bin/env bash
# Fixed-seed cross-regression for VBS mode selection and deferred phase.
# Usage: bash scripts/regress_vbs_deferred_modes.sh <build_dir> <work_dir>
set -euo pipefail

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export LD_LIBRARY_PATH=${MKLROOT:-/opt/intel/oneapi/mkl/2026.1}/lib:${LD_LIBRARY_PATH:-}

build=$(cd -- "${1:?build directory}" && pwd)
work_input=${2:?work directory}
mkdir -p "$work_input"
work=$(cd "$work_input" && pwd)
bin="$build/kitaev_ssh_mc"
[ -x "$bin" ] || { echo "missing $bin" >&2; exit 2; }

qf="$work/q_points_cross_regress.txt"
printf '%s\n' \
  '# fixed q-file shared by selected and both' \
  'ivbs_Qp 0.5804306 0.4195694' \
  'g_q10 0.5 0.5' \
  'kek_Kp 0.3333333 0.6666667' \
  'qzero 0.0 0.0' > "$qf"

common=(
  --Lx 3 --Ly 2 --Jx 1 --Jy 1 --Jz 0.5
  --intervals 4 --dt 0.2 --mass 1.7 --omega-ph 0.8 --g 0.3
  --seed 314159 --warmup 8 --sweeps 32 --measure-every 2
  --block-samples 4 --qr-block 4 --proposal bond-refresh
  --weight-engine combined --vbs 1
)

for defer in 0 1; do
  for mode in all selected both; do
    tag="$work/${mode}_defer${defer}"
    qarg=()
    if [ "$mode" != all ]; then qarg=(--vbs-q-file "$qf"); fi
    "$bin" "${common[@]}" --defer-phase "$defer" --vbs-mode "$mode" \
      "${qarg[@]}" --output "$tag" >"$tag.stdout" 2>"$tag.stderr"
  done
done

PYTHON=${PYTHON:-/home/asus/.venv/bin/python}
"$PYTHON" - "$work" <<'PY'
import csv
import json
import pathlib
import sys

work = pathlib.Path(sys.argv[1])
modes = ("all", "selected", "both")
defers = (0, 1)
basic = ("X_x", "X_y", "X_z", "X2_x", "X2_y", "X2_z", "majorana_energy")
csv_common = ["block", "count", "sum_phase_real", "sum_phase_imag",
              "positive_count", "negative_count", "complex_count"]
for name in basic:
    csv_common.extend((f"sum_phase_{name}_real", f"sum_phase_{name}_imag"))

def result(mode, defer):
    return json.loads((work / f"{mode}_defer{defer}_results.json").read_text())

def blocks(mode, defer):
    with (work / f"{mode}_defer{defer}_blocks.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return [tuple(row[key] for key in csv_common) for row in rows]

data = {(mode, defer): result(mode, defer) for defer in defers for mode in modes}

def equal(label, left, right):
    if left != right:
        raise AssertionError(label)

def trajectory(d):
    return tuple(d[key] for key in (
        "proposals_warmup", "accepted_warmup", "proposals_measurement",
        "accepted_measurement", "untrusted_rejected_warmup",
        "untrusted_rejected_measurement", "untrusted_rejected_initial"))

def final_state(d):
    equal("center snapshot disagrees with final path", d["final_center_coordinates"],
          d["final_phonon_path"][d["final_center_node"]])
    return (d["final_center_node"], d["final_phonon_path"],
            d["final_center_coordinates"])

def basic_observables(d):
    return {name: d["observables"][name] for name in basic}

def phase_summary(d):
    return (d["average_phase"], d["phase_diagnostics"], d["ratio_defined"],
            d["jackknife_defined"], d["phase_resolved_5sigma"], d["blocks"])

def common_measurement(d, mode, defer):
    return (trajectory(d), final_state(d), basic_observables(d),
            phase_summary(d), blocks(mode, defer))

for defer in defers:
    reference = data[("all", defer)]
    for mode in ("selected", "both"):
        equal(f"{defer}: all/{mode} trajectory/physics/phase mismatch",
              common_measurement(reference, "all", defer),
              common_measurement(data[(mode, defer)], mode, defer))
    equal(f"{defer}: selected/both trajectory/physics/phase mismatch",
          common_measurement(data[("selected", defer)], "selected", defer),
          common_measurement(data[("both", defer)], "both", defer))

for mode in modes:
    eager = data[(mode, 0)]
    deferred = data[(mode, 1)]
    equal(f"{mode}: eager/deferred trajectory/physics/phase mismatch",
          common_measurement(eager, mode, 0),
          common_measurement(deferred, mode, 1))
    equal(f"{mode}: eager/deferred VBS blocking mismatch",
          eager["vbs"]["blocking"], deferred["vbs"]["blocking"])
    if mode != "all":
        equal(f"{mode}: eager/deferred candidate diagnostics mismatch",
              eager["vbs"]["candidates"], deferred["vbs"]["candidates"])
        candidate = eager["vbs"]["candidates"]
        equal(f"{mode}: candidate schema", candidate["schema"], "ivbs-candidates-v2")
        equal(f"{mode}: candidate schema revision", candidate["schema_revision"], 2)
        for state in ("center_defined", "jackknife_defined", "stability_defined"):
            if not candidate[state]:
                raise AssertionError(f"{mode}: candidate {state} is false")

for mode in modes:
    for defer in defers:
        d = data[(mode, defer)]
        equal(f"{mode}/{defer}: PBC boundary", d["spatial_boundary"], "pbc")
        equal(f"{mode}/{defer}: twist flag", d["vbs"]["twist_applied"], False)
        equal(f"{mode}/{defer}: CLI defer metadata", d["defer_phase"], defer)

report = {
    "schema": "vbs_deferred_cross_regression_v1",
    "system": {"Lx": 3, "Ly": 2, "seed": 314159, "q_file": str(work / "q_points_cross_regress.txt")},
    "matrix": [f"{mode}/defer{defer}" for defer in defers for mode in modes],
    "checks": [
        "all_selected_both_fixed_seed_trajectory_physics_phase",
        "eager_deferred_fixed_seed_trajectory_physics_phase",
        "final_phonon_path_and_center_coordinates",
        "blocks_and_candidate_diagnostics",
        "pbc_and_no_twist_metadata",
    ],
    "status": "PASS",
}
(work / "cross_mode_regression.json").write_text(json.dumps(report, indent=2) + "\n")
print("CROSS_VBS_DEFERRED_MODES PASS (all/selected/both x eager/deferred)")
PY
