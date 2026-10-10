#!/usr/bin/env bash
# Regression: with a fixed seed, the all/selected/both VBS measurement modes
# must produce bitwise-identical trajectories. Compares the phase/count
# columns and the 7 basic-observable columns of *_blocks.csv (printed with 17
# significant digits, so identical doubles give identical strings) plus the
# acceptance counters of the JSON summaries. Includes a non-square lattice.
# Usage: bash scripts/regress_vbs_modes.sh <build_dir> <work_dir>
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=${MKLROOT:-/opt/intel/oneapi/mkl/2026.1}/lib:${LD_LIBRARY_PATH:-}
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:?build directory}; shift || true
work=${1:?work directory}
mkdir -p "$work"
work=$(cd "$work" && pwd)
bin=$build/kitaev_ssh_mc
[ -x "$bin" ] || { echo "missing $bin" >&2; exit 2; }
qf=$work/q_points_regress.txt
cat > "$qf" <<'EOF'
# regression q list (subset of q_points_example.txt logic)
ivbs_Qp 0.5804306 0.4195694
g_q10 0.5 0.5
kek_Kp 0.3333333 0.6666667
qzero 0.0 0.0
EOF

run_modes() { # run_modes Lx Ly tag
  local Lx=$1 Ly=$2 tag=$3
  local common="--Lx $Lx --Ly $Ly --Jx 1 --Jy 1 --Jz 0.5 --intervals 4 --dt 0.2 \
    --mass 1.7 --omega-ph 0.8 --g 0.3 --seed 42 --warmup 20 --sweeps 80 \
    --measure-every 2 --block-samples 10 --qr-block 4"
  (cd "$work" \
    && "$bin" $common --vbs 1 --vbs-mode all      --output "${tag}_all"  >"${tag}_all.json" \
    && "$bin" $common --vbs 1 --vbs-mode selected --vbs-q-file "$qf" --output "${tag}_sel" >"${tag}_sel.json" \
    && "$bin" $common --vbs 1 --vbs-mode both     --vbs-q-file "$qf" --output "${tag}_both" >"${tag}_both.json")
}
compare() { # compare tag
  local tag=$1
  "$PYTHON" - "$work" "$tag" "$root" <<'PY'
import json,sys,csv,os
work,tag,root=sys.argv[1],sys.argv[2],sys.argv[3]
def cols(path,n):
    with open(path) as f:
        r=csv.reader(f); head=next(r); rows=[row for row in r]
    # columns: block,count,phase(2),pos,neg,complex, then per-name real/imag
    keep=list(range(0,7))+[7+2*j+k for j in range(n) for k in (0,1)]
    return head,[ [row[i] for i in keep] for row in rows ]
BASIC=7
heads=None;data={}
for mode in ("all","sel","both"):
    h,d=cols(f"{work}/{tag}_{mode}_blocks.csv",BASIC)
    heads=h;data[mode]=d
assert data["all"]==data["sel"]==data["both"], f"{tag}: blocks.csv trajectory columns differ"
j={m:json.load(open(f"{work}/{tag}_{m}.json")) for m in ("all","sel","both")}
acc=j["all"]["acceptance"]
for m in ("sel","both"):
    for k in ("acceptance","proposals_warmup","accepted_warmup","proposals_measurement","accepted_measurement"):
        assert j[m][k]==j["all"][k], (tag,m,k)
    for name in ("X_x","X_y","X_z","X2_x","X2_y","X2_z","majorana_energy"):
        assert j[m]["observables"][name]["value"]==j["all"]["observables"][name]["value"], (tag,m,name)

def block(data):
    assert data["vbs"]["blocking"], "missing vbs blocking"
    return data["vbs"]["blocking"][0]

all_block,sel_block,both_block=(block(j[m]) for m in ("all","sel","both"))
nc=j["all"]["Lx"]*j["all"]["Ly"]
nq=j["sel"]["vbs"]["Nq"]
assert j["all"]["vbs"]["mode"]=="all" and len(j["all"]["vbs"]["momenta"])==nc
for name in ("S_q0_xx","Delta_q0_x","TrS_q0","C_q0_xx","C_connected_q0_xx"):
    assert name in all_block["raw"] or name in all_block["derived"], (tag,"all_old_key",name)
assert len(sel_block["derived"])==8*nq+3, (tag,"selected_derived_length")
assert not any(name.startswith("TrS_q") for name in sel_block["derived"]), (tag,"selected_grid_derived")
assert len(both_block["derived"])==22*nc+8*nq+3, (tag,"both_derived_length")
all_derived_names=list(all_block["derived"])
both_derived_names=list(both_block["derived"])
assert both_derived_names[:len(all_derived_names)]==all_derived_names, (tag,"both_derived_prefix")
all_raw_names=list(all_block["raw"])
both_raw_names=list(both_block["raw"])
assert both_raw_names[:len(all_raw_names)]==all_raw_names, (tag,"both_raw_prefix")
for name in all_raw_names:
    assert both_block["raw"][name]==all_block["raw"][name], (tag,"both_raw_value",name)
for name in all_derived_names:
    assert both_block["derived"][name]==all_block["derived"][name], (tag,"both_derived_value",name)

def check_candidates(data, label):
    v=data["vbs"]; c=v.get("candidates")
    assert c and c["schema"]=="ivbs-candidates-v2" and c["level"]=="finest_blocks", (tag,label,"candidate_schema")
    assert c["center_defined"] and c["jackknife_defined"] and c["stability_defined"], (tag,label,"candidate_states")
    assert c["total_replicas"]==c["valid_replicas"]+c["invalid_replicas"], (tag,label,"candidate_replica_counts")
    names=[g["name"] for g in c["groups"]]
    assert names==["grid_raw","offgrid_raw","offgrid_projected"], (tag,label,"candidate_groups")
    qpoints={q["index"]:q for q in v["qpoints"]}
    grid=sorted(i for i,q in qpoints.items() if q["on_grid"] and q["n"]!=[0,0])
    off=sorted(i for i,q in qpoints.items() if not q["on_grid"])
    groups={g["name"]:g for g in c["groups"]}
    assert sorted(groups["grid_raw"]["participants"])==grid, (tag,label,"grid_partition")
    assert sorted(groups["offgrid_raw"]["participants"])==off, (tag,label,"offgrid_raw_partition")
    assert sorted(groups["offgrid_projected"]["participants"])==off, (tag,label,"offgrid_projected_partition")
    for g in groups.values():
        assert g["available"]==(len(g["participants"])>0), (tag,label,g["name"],"availability")
        assert g["total_replicas"]==c["total_replicas"] and g["valid_replicas"]==c["valid_replicas"] and g["invalid_replicas"]==c["invalid_replicas"], (tag,label,g["name"],"replica_counts")
        if g["available"]:
            assert g["center_defined"] and g["jackknife_defined"] and g["stability_defined"], (tag,label,g["name"],"group_states")
            assert len(g["unique_win_frequency"])==len(g["participants"]), (tag,label,g["name"],"frequency_length")
    for point in c["points"]:
        assert "W_L_abs2" in point, (tag,label,"point_window",point["index"])
        for key in ("defined","jackknife_defined","cov_real_imag"):
            assert key in point["mz2_raw"] and key in point["mz2_projected"], (tag,label,"point_state",point["index"],key)
    derived=block(data)["derived"]
    for name in ("candidate_max_winner_index","candidate_participant_count"):
        assert name in derived and derived[name]["stderr_real"] is None and derived[name]["stderr_imag"] is None, (tag,label,name,"categorical_error")
    summary_csv=f"{work}/{tag}_{label}_candidates_summary.csv"
    summary_json=f"{work}/{tag}_{label}_candidates_summary.json"
    import subprocess
    subprocess.run([sys.executable,os.path.join(root,"kitaev_ssh","summarize_vbs_candidates.py"),
                    f"{work}/{tag}_{label}.json","--csv",summary_csv,"--json",summary_json],
                   check=True)
    summary=json.load(open(summary_json))
    assert summary["status"]=="available" and len(summary["rows"])==v["Nq"], (tag,label,"summary")
    with open(summary_csv) as stream:
        rows=list(csv.DictReader(stream))
    assert len(rows)==v["Nq"] and rows[0]["status"]=="available", (tag,label,"summary_csv")
    return summary

check_candidates(j["sel"],"sel")
check_candidates(j["both"],"both")
print(f"{tag}: trajectory/layout/candidate/postprocess checks PASS (acc={acc:.5f}, Nc={nc}, Nq={nq})")
PY
}
PYTHON=${PYTHON:-/home/asus/.venv/bin/python}
run_modes 4 4 L44
compare L44
run_modes 5 3 L53
compare L53
echo "REGRESS_VBS_MODES PASS"
