#!/usr/bin/env python3
"""PNG-only views of recorded results; no refitting or tolerance adjustment."""
import argparse,csv,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/interacting-ed-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('directory');a=p.parse_args()
from pathlib import Path
root=Path(a.directory)
with (root/'QMC-ED-observables.csv').open() as f:rows=list(csv.DictReader(f))
cases=list(dict.fromkeys(r['case'] for r in rows))
labels=[c.replace('bond-guided','guided').replace('production','prod').replace('-tau',' t').replace('-dt',' d').replace('-',' ') for c in cases]
fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
for ax,name in zip(axes,['majorana_energy','chi_cos_q3','B4_q3_z']):
 vals=[next(r for r in rows if r['case']==c and r['observable']==name) for c in cases]
 ax.errorbar([float(r['difference_finite_dt'])/float(r['MC_standard_error']) for r in vals],range(len(vals)),xerr=1,fmt='o',capsize=3)
 ax.axvline(0,color='black',lw=.8);ax.axvspan(-6,6,color='tab:green',alpha=.07)
 ax.set_title(name);ax.set_xlabel('(QMC - matched ED) / MC SE');ax.set_yticks(range(len(cases)),labels if ax is axes[0] else ['']*len(cases));ax.set_xlim(-6.5,6.5);ax.grid(axis='x',alpha=.2)
fig.suptitle('All 28 chains; error bars = 1 MC SE. Local precision incomplete.');fig.savefig(root/'QMC-ED.png',dpi=160);plt.close(fig)
with (root/'blocking.csv').open() as f:blocks=list(csv.DictReader(f))
bad=[r for r in rows if r['status']!='PASS'];fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
for r in bad:
 b=[v for v in blocks if v['case']==r['case'] and v['observable']==r['observable']]
 ax.plot([int(v['block_sweeps']) for v in b],[float(v['stderr']) for v in b],'o-',label=r['observable'])
ax.axhline(.005,color='black',ls='--',label='predeclared absolute SE target');ax.set_xscale('log',base=2);ax.set_xlabel('Sweeps per block (chains kept separate)');ax.set_ylabel('MC standard error');ax.set_title('Local proposal: four precision gates unmet; plateau gates pass');ax.legend();ax.grid(alpha=.2);fig.savefig(root/'local-precision-blocking.png',dpi=160);plt.close(fig)
