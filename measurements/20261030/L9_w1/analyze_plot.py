import json,gzip,csv,collections,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent;allc=[];meta=[]
for port in [20015]:
 with gzip.open(D/('snapshot-%d.json.gz'%port),'rt') as f:d=json.load(f)
 allc.extend(d['chains']);meta.append({'port':port,'captured_utc':d['captured_utc'],'status':dict(collections.Counter(c['status'] for c in d['chains']))})
dedup={}
for c in allc:
 if c['status']!='captured':continue
 p=c['point'];k=(round(float(p['lambda_value']),6),round(float(p['tau']),6),int(p['seed']))
 if k not in dedup or c['retained']>dedup[k]['retained']:dedup[k]=c
groups=collections.defaultdict(list)
for (lam,tau,seed),c in dedup.items():groups[lam,tau].append(c)
def est(a):
 den=complex(a[0],a[1]);assert abs(den)>1e-12
 s=[(complex(a[2+2*i],a[3+2*i])/den).real/81 for i in range(80)];q=int(np.argmax(s));v=s[q];return (math.sqrt(v) if v>=0 else float('nan')),v,q+1
rows=[]
for (lam,tau),g in sorted(groups.items()):
 raw=np.sum([c['raw'] for c in g],axis=0);m,m2,q=est(raw);n=len(g);rep=[est(raw-np.array(c['raw'])) for c in g] if n>1 else [];e=[float(np.sqrt((n-1)/n*np.sum((np.array([v[i] for v in rep])-np.mean([v[i] for v in rep]))**2))) for i in (0,1)] if rep else [None,None]
 rows.append(dict(lambda_value=lam,tau=tau,dt=.1,seeds=n,saved_min=min(c['saved'] for c in g),saved_max=max(c['saved'] for c in g),retained=sum(c['retained'] for c in g),m=m,m_error=e[0],m2=m2,m2_error=e[1],q=q,phase=abs(complex(raw[0],raw[1]))/sum(c['retained'] for c in g),sources=';'.join(sorted(set(c['root'] for c in g)))))
with (D/'summary.csv').open('w',newline='') as f:
 fields=['lambda_value','tau','dt','seeds','saved_min','saved_max','retained','m','m_error','m2','m2_error','q','phase','sources'];w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
for metric,label,out in [('m',r'$m_z$', 'order-vs-tau.png'),('m2',r'$m_z^2$', 'order-squared-vs-tau.png')]:
 fig,ax=plt.subplots(figsize=(8,5.8),layout='constrained')
 for lam in [6,7,8]:
  rr=sorted([r for r in rows if abs(r['lambda_value']-lam)<1e-8 and r[metric]>0],key=lambda r:r['tau'])
  if not rr:continue
  x=[r['tau'] for r in rr];y=[r[metric] for r in rr];er=[r[metric+'_error'] if r[metric+'_error'] is not None and r[metric+'_error']<r[metric] else float('nan') for r in rr];ax.errorbar(x,y,yerr=er,marker='o',markersize=4,lw=1,capsize=2,label=rf'$\lambda={lam:g}$')
 ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'$\tau$');ax.set_ylabel(label);ax.grid(alpha=.2,which='both');ax.legend(fontsize=10);ax.set_title(r'$L=9,\ \omega=1,\ \Delta\tau=0.1$');fig.text(.5,-.02,'Cut: 5000 sweeps per chain; errors: independent-seed jackknife',ha='center',fontsize=9);fig.savefig(D/out,dpi=180,bbox_inches='tight');plt.close(fig)
(D/'analysis.json').write_text(json.dumps({'meta':meta,'selected_chains':len(dedup),'cells':len(rows),'lambda_values':sorted(set(r['lambda_value'] for r in rows))},indent=2)+'\n');print((D/'analysis.json').read_text())
