from pathlib import Path
import json,gzip,csv,hashlib,collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent;SRC=D/'snapshot.json.gz';Q=[15,20,22,27]
with gzip.open(SRC,'rt') as f:data=json.load(f)
assert data['cut']==4000
groups=collections.defaultdict(list)
for c in data['chains']:
 if c['status']=='captured':groups[float(c['point']['lambda_value']),float(c['point']['tau'])].append(c)
def estimate(a):
 den=complex(a[0],a[1]);assert abs(den)>1e-12
 vals={q:(complex(a[2+8*(q-1)],a[3+8*(q-1)])/den).real/36 for q in Q};q=max(Q,key=lambda q:vals[q]);return vals[q],q,vals
rows=[]
for (lam,tau),g in sorted(groups.items()):
 total=np.sum([c['raw'] for c in g],axis=0);m2,q,values=estimate(total);n=len(g);reps=[estimate(total-np.array(c['raw']))[0] for c in g] if n>1 else [];se=float(np.sqrt((n-1)/n*np.sum((np.array(reps)-np.mean(reps))**2))) if reps else None
 rows.append(dict(lambda_value=lam,tau=tau,seeds=n,m2=m2,m2_error=se,q_peak=q,saved_min=min(c['saved'] for c in g),saved_max=max(c['saved'] for c in g),retained=sum(c['retained'] for c in g),**{'m2_q%d'%q:values[q] for q in Q}))
with (D/'four-q-summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,ax=plt.subplots(figsize=(8.8,6.4),layout='constrained');colors=plt.get_cmap('tab10')
for i,lam in enumerate(sorted(set(r['lambda_value'] for r in rows))):
 rr=sorted([r for r in rows if r['lambda_value']==lam and r['m2']>0],key=lambda r:r['tau']);x=[r['tau'] for r in rr];y=[r['m2'] for r in rr];e=[r['m2_error'] if r['m2_error'] is not None and r['m2_error']<r['m2'] else float('nan') for r in rr]
 ax.errorbar(x,y,yerr=e,color=colors(i),marker='o',markersize=3,lw=1,capsize=2,label=rf'$\lambda={lam:g}$')
anchor=next(r for r in rows if abs(r['lambda_value']-1.8)<1e-8 and abs(r['tau']-8)<1e-8);xx=np.array([4.,40.]);ax.plot(xx,anchor['m2']*(xx/8.)**(2./3.),'k--',lw=1.5,label=r'$2/3$ guide ($\lambda=1.8,\ \tau=8$)')
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'$\tau$');ax.set_ylabel(r'$m_z^2=\max_{q\in\{15,20,22,27\}}\mathrm{Re}\,S^{zz}(q)/36$');ax.set_title(r'$L=6,\ \omega=0.1,\ \Delta\tau=0.4$');ax.grid(True,which='both',alpha=.2);ax.legend(ncol=2,fontsize=9);fig.text(.5,-.015,'Original 3 seeds; cut: 4000 sweeps per chain; seed jackknife errors',ha='center',fontsize=9);fig.savefig(D/'m2-vs-tau-four-q-slope-2over3.png',dpi=180,bbox_inches='tight');plt.close(fig)
meta=dict(snapshot=SRC.name,snapshot_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),captured_utc=data['captured_utc'],cut=4000,momenta=Q,cells=len(rows),seeds=sorted(set(r['seeds'] for r in rows)),anchor={'lambda':1.8,'tau':8,'m2':anchor['m2']},positive_cells=sum(r['m2']>0 for r in rows),definition='phase pool first, max of four momenta, delete-one-seed jackknife with peak re-selection')
(D/'four-q-analysis.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
