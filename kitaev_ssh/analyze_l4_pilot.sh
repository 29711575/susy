#!/bin/bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
python=${P5_PYTHON:-/home/asus/.cache/pfqmc_plot_env/bin/python}
data=${L4_PILOT_DIR:-$root/build-kitaev-p01/l4_pilot}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/kitaev-l4-matplotlib}
"$python" - "$root" "$data" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1])/'kitaev_ssh'))
from plot_vbs import load_run, write_csv
directory=Path(sys.argv[2])
rows=[]
for seed in (7601,7602,7603,7604):
    d, consistency=load_run(directory/'results'/f'l4_{seed}_results.json')
    for level in d['vbs']['blocking']:
        row=dict(seed=seed, tau=d['tau'], elapsed_seconds=d['elapsed_seconds'],
                 seconds_per_sweep=d['elapsed_seconds']/(d['warmup_sweeps']+d['measurement_sweeps']),
                 acceptance=d['acceptance'], numerical_failures=d['numerical_failures'],
                 phase_real=d['average_phase']['value'][0], phase_imag=d['average_phase']['value'][1],
                 blocks=level['blocks'], block_sweeps=level['block_sweeps'], csv_ratio_error=consistency)
        for label, family, name in [('Xz','raw','Delta_q0_z'),('Szz_02','raw','S_q2_zz'),
                                   ('Szz_22','raw','S_q10_zz'),('TrS_02','derived','TrS_q2'),
                                   ('TrS_22','derived','TrS_q10')]:
            e=level[family][name]
            row[label]=e['value'][0]; row[label+'_se']=e['stderr_real']
        rows.append(row)
write_csv(directory/'summary.csv',rows)
print(directory/'summary.csv')
PY
