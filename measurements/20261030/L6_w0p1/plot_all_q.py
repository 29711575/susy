import gzip,json,collections,csv,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent
with gzip.open(D/'snapshot.json.gz','rt') as f:data=json.load(f)
groups=collections.defaultdict(list)
for c in data['chains']:
 if c['status']=='captured':groups[(c['point']['lambda_value'],c['point']['tau'])].append(c)
def estimate(a,count):
 ph=complex(a[0],a[1]);
 def v(i):return complex(a[2+2*i],a[3+2*i])/ph
 s=[v(4*q).real/36 for q in range(35)];chi=[(v(4*q+1)-36*v(4*q+2)*v(4*q+3)).real for q in range(35)];qm=int(np.argmax(s));qc=int(np.argmax(chi));m2=s[qm];return {'m2':m2,'m':math.sqrt(m2) if m2>=0 else float('nan'),'chi':chi[qc],'q_m':qm+1,'q_chi':qc+1,'phase':abs(ph)/count}
def err(g,total,key):
 if len(g)<2:return None
 vals=[estimate(total-np.array(c['raw']),sum(z['retained'] for z in g if z is not c))[key] for c in g];av=np.mean(vals);return float(np.sqrt((len(g)-1)/len(g)*sum((x-av)**2 for x in vals)))
rows=[]
for (lam,tau),g in sorted(groups.items()):
 total=np.sum([c['raw'] for c in g],axis=0);r=estimate(total,sum(c['retained'] for c in g));r.update(lambda_value=lam,tau=tau,seeds=len(g),saved_min=min(c['saved'] for c in g),saved_max=max(c['saved'] for c in g),retained=sum(c['retained'] for c in g),m_error=err(g,total,'m'),m2_error=err(g,total,'m2'),chi_error=err(g,total,'chi'));rows.append(r)
fields=['lambda_value','tau','seeds','saved_min','saved_max','retained','m','m_error','m2','m2_error','chi','chi_error','q_m','q_chi','phase']
with open(D/'summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)
colors=plt.get_cmap('tab10');ls=sorted({r['lambda_value'] for r in rows})
for key,label,out in [('m',r'$m$', 'order-vs-tau.png'),('m2',r'$m^2$', 'order-squared-vs-tau.png'),('chi',r'$\chi$', 'susceptibility-vs-tau.png')]:
 fig,ax=plt.subplots(figsize=(8.4,6),layout='constrained')
 for i,lam in enumerate(ls):
  data1=sorted([r for r in rows if r['lambda_value']==lam and r['seeds']==3 and r[key]>0],key=lambda x:x['tau']); partial=sorted([r for r in rows if r['lambda_value']==lam and r['seeds']<3 and r[key]>0],key=lambda x:x['tau']);col=colors(i)
  if data1:
   x=[r['tau'] for r in data1];y=[r[key] for r in data1];e=[r[key+'_error'] if r[key+'_error'] is not None and r[key+'_error']<r[key] else float('nan') for r in data1];ax.errorbar(x,y,yerr=e,color=col,marker='o',markersize=3,linewidth=1,capsize=2,label=rf'$\lambda={lam:g}$')
  if partial:ax.scatter([r['tau'] for r in partial],[r[key] for r in partial],facecolors='none',edgecolors=[col],s=30)
 if key=='m2':
  anchor=next((r for r in rows if abs(r['lambda_value']-1.8)<1e-9 and abs(r['tau']-8)<1e-9),None)
  if anchor:
   guide=np.array([4.,19.2]);ax.plot(guide,anchor['m2']*(guide/8.)**(2./3.),'k--',linewidth=1.4,label=r'$2/3$ guide (anchored at $\lambda=1.8,\ \tau=8$)')
 ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'$\tau$');ax.set_ylabel(label);ax.set_title(r'$L=6,\ \omega=0.1,\ \Delta\tau=0.4$');ax.grid(True,which='both',alpha=.2);ax.legend(ncol=2,fontsize=9);fig.text(.5,-.015,'Cut: 4000 sweeps; open circles: fewer than 3 seeds',ha='center',fontsize=9);fig.savefig(D/out,dpi=180,bbox_inches='tight');plt.close(fig)
json.dump({'captured_utc':data['captured_utc'],'cut':4000,'statuses':dict(collections.Counter(c['status'] for c in data['chains'])),'cells':len(rows),'complete_seed_cells':sum(r['seeds']==3 for r in rows),'chains_included':sum(r['seeds'] for r in rows)},open(D/'analysis.json','w'),indent=2)
print(json.dumps(json.load(open(D/'analysis.json')),indent=2))
