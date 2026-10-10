#!/usr/bin/env python3
"""Summarize the bounded four-mode benchmark; no images or MC inference."""
import argparse
import csv
import json
from pathlib import Path
from statistics import mean

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args=parser.parse_args();directory=args.directory
    labels=['reference','svd-values','node-merge','combined']
    runs={mode:[json.loads((directory/f'bench_{repeat}_{mode}.json').read_text()) for repeat in (1,2)] for mode in range(4)}
    baseline=mean(r['seconds_per_sweep'] for r in runs[0])
    rows=[]
    max_log=max_phase=0.
    for mode,label in enumerate(labels):
        values=[r['seconds_per_sweep'] for r in runs[mode]]
        rows.append(dict(engine=label,seconds_per_sweep=values,mean_seconds_per_sweep=mean(values),speedup=baseline/mean(values)))
        for i,r in enumerate(runs[mode]):
            old=runs[0][i]
            if r['accepted']!=old['accepted']:
                raise ValueError('benchmark acceptance counts differ; investigate trajectories before comparing')
            max_log=max(max_log,abs(r['final_log_abs']-old['final_log_abs']))
            max_phase=max(max_phase,abs(complex(*r['final_phase'])-complex(*old['final_phase'])))
    blocks={label:list(csv.DictReader((directory/f'smoke_{label}_blocks.csv').open())) for label in labels}
    max_block=0.
    for label in labels[1:]:
        if len(blocks[label])!=len(blocks['reference']):raise ValueError('different block counts')
        for a,b in zip(blocks['reference'],blocks[label]):
            if set(a)!=set(b):raise ValueError('different block schemas')
            max_block=max(max_block,max(abs(float(a[k])-float(b[k])) for k in a))
    if max_log>2e-10 or max_phase>2e-10 or max_block>2e-10:
        raise ValueError('cross-mode end-to-end numerical discrepancy')
    profile=json.loads((directory/'combined_profile.json').read_text())
    out=dict(throughput=rows,max_final_log_error=max_log,max_final_phase_error=max_phase,
             max_smoke_block_error=max_block,accepted_per_30_sweeps=runs[0][0]['accepted'],
             replay=json.loads((directory/'replay.json').read_text()),
             instrumented_combined_seconds_per_sweep=profile['seconds_per_sweep'],
             instrumented_hotspot_percent={k:100*v/profile['elapsed_seconds'] for k,v in profile['timings'].items()},
             timing_scope='30 sweeps incl center observables and statistics, excludes initialization/serialization; no operation timers in primary benchmark')
    target=directory/'benchmark_summary.json'
    if target.exists():raise ValueError('summary exists; preserve earlier result')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
