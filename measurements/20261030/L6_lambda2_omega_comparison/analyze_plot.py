import gzip,json,collections,csv,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parent; data=[]
for port in (20015,20022):
 with gzip.open(D/('snapshot-%d.json.gz'%port),'rt') as f:data.append(json.load(f))
groups=collections.defaultdict(list)
for d in data:
 for c in d['chains']:
  if c['status']=='captured':
   p=c['point'];groups[(float(p['omega']),float(p['dt']),float(p['tau']))].append(c)
def est(a,count):
 ph=complex(a[0],a[1]);
 def v(i):return complex(a[2+2*i],a[3+2*i])/ph
 s=[v(4*q).real/36 for q in range(35)];chis=[(v(4*q+1)-36*v(4*q+2)*v(4*q+3)).real for q in range(35)];q=int(np.argmax(s));z=int(np.argmax(chis));return {'m2':s[q], 'm':math.sqrt(s[q]) if s[q]>0 else float('nan'), 'chi':chis[z], 'q_m':q+1,'q_chi':z+1,'phase':abs(ph)/count}
def jack(g,total,key):
 if len(g)<2:return None
 vals=[est(total-np.array(c['raw']),sum(y['retained'] for y in g if y is not c))[key] for c in g];mean=np.mean(vals);return float(np.sqrt((len(g)-1)/len(g)*sum((v-mean)**2 for v in vals)))
rows=[]
for (om,dt,tau),g in sorted(groups.items()):
 total=np.sum([c['raw'] for c in g],axis=0);r=est(total,sum(c['retained'] for c in g));r.update(omega=om,dt=dt,tau=tau,seeds=len(g),saved_min=min(c['saved'] for c in g),saved_max=max(c['saved'] for c in g),retained=sum(c['retained'] for c in g),m_error=jack(g,total,'m'),m2_error=jack(g,total,'m2'),chi_error=jack(g,total,'chi'));rows.append(r)
fields=['omega','dt','tau','seeds','saved_min','saved_max','retained','m','m_error','m2','m2_error','chi','chi_error','q_m','q_chi','phase']
with open(D/'summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)
for key,label,out in [('m',r'$m$', 'order-comparison.png'),('m2',r'$m^2$', 'order-squared-comparison.png'),('chi',r'$\chi$', 'susceptibility-comparison.png')]:
 fig,ax=plt.subplots(figsize=(8.3,5.7),layout='constrained');colors={.1:'#2878b5',.2:'#ef8636',.3:'#36934a',.5:'#bb486c'}
 for om,dt in [(.1,.2),(.1,.4),(.2,.4),(.3,.4),(.5,.4)]:
  group=sorted([r for r in rows if r['omega']==om and r['dt']==dt and r[key]>0],key=lambda r:r['tau'])
  if not group:continue
  x=[r['tau'] for r in group];y=[r[key] for r in group];e=[r[key+'_error'] if r[key+'_error'] is not None and r[key+'_error']<r[key] else float('nan') for r in group];style='--' if dt==.2 else '-';ax.errorbar(x,y,yerr=e,color=colors[om],linestyle=style,marker='o' if dt==.4 else 's',markersize=4,linewidth=1.4,capsize=2,label=rf'$\omega={om:g},\ dt={dt:g}$')
 ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel(r'$\tau$');ax.set_ylabel(label);ax.set_title(r'$L=6,\ \lambda=2$');ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=9);fig.text(.5,-.02,'Cut: 4000 sweeps; error bars: seed jackknife',ha='center',fontsize=9);fig.savefig(D/out,dpi=180,bbox_inches='tight');plt.close(fig)
print(json.dumps(dict(rows=len(rows),groups=[{'omega':om,'dt':dt,'points':sum(r['omega']==om and r['dt']==dt for r in rows)} for om,dt in [(.1,.2),(.1,.4),(.2,.4),(.3,.4),(.5,.4)]],statuses=[dict(collections.Counter(c['status'] for c in d['chains'])) for d in data]),indent=2))
