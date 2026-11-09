"""Offline reproduction: python analyze.py [--plots]. Requires numpy; plots require matplotlib."""
from pathlib import Path
import csv, gzip, json, math, sys
from collections import defaultdict
import numpy as np
D=Path(__file__).resolve().parent;OUT=D/'derived';OUT.mkdir(exist_ok=True)
def save_csv(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
G=defaultdict(list)
data=json.load(gzip.open(D/'data/seed-moments.json.gz','rt'))
for c in data['chains']:
 p=c['point'];assert p['J']==[1,1,.6] and p['lambda_value']==2
 G[(p['L'],p['omega'],p['dt'],p['tau'])].append(c)
def est(cs,L):
 a=np.sum([c['raw'] for c in cs],axis=0);ph=complex(*a[:2]);z=(a[2::2]+1j*a[3::2])/ph;z=z.reshape(-1,4)
 m=z[:,0].real/L**2;chi=(z[:,1]-L**2*z[:,2]*z[:,3]).real
 return dict(m2=float(max(m)),chi=float(max(chi)),m2_q=int(np.argmax(m))+1,chi_q=int(np.argmax(chi))+1,phase_real=float(ph.real/sum(c['retained'] for c in cs))),m,chi
rows=[];spectra=[]
for (L,w,dt,t),cs in sorted(G.items()):
 assert len({c['point']['seed'] for c in cs})==len(cs)
 v,m,ch=est(cs,L);n=len(cs);jk=[est(cs[:i]+cs[i+1:],L)[0] for i in range(n)] if n>1 else []
 for k in ('m2','chi','phase_real'):
  vals=np.array([r[k] for r in jk]);v[k+'_error']=float(np.sqrt((n-1)/n*np.sum((vals-np.mean(vals))**2))) if n>1 else None
 v.update(L=L,omega=w,dt=dt,tau=t,seeds=n,retained=sum(c['retained'] for c in cs),phase_resolved=abs(v['phase_real'])>5*v['phase_real_error'] if n>1 else None);rows.append(v)
 for q,(mq,cq) in enumerate(zip(m,ch),1):spectra.append(dict(L=L,omega=w,dt=dt,tau=t,q=q,n1=q//L,n2=q%L,m2_q=mq,chi_q=cq))
ref=list(csv.DictReader(open(D/'data/reference-results.csv')))
for r in rows:
 old=next(x for x in ref if all(float(x[k])==r[k] for k in ('L','omega','dt','tau')))
 for k in ('m2','chi','phase_real'):assert abs(r[k]-float(old[k]))<1e-11,(r,k)
save_csv('recomputed-results.csv',rows);save_csv('momentum-spectra.csv',spectra)
def fit(rs,k):
 x=np.log([r['tau'] for r in rs]);y=np.array([r[k] for r in rs]);e=np.array([r[k+'_error'] for r in rs]);X=np.column_stack([np.ones(len(x)),x]);W=np.diag((y/e)**2);V=np.linalg.inv(X.T@W@X);b=V@X.T@W@np.log(y)
 return dict(slope=float(b[1]),slope_SE=float(np.sqrt(V[1,1])),chi2_dof=float(np.sum(((np.log(y)-X@b)*y/e)**2)/(len(x)-2)),unweighted_slope=float(np.polyfit(x,np.log(y),1)[0]))
fits=[]
for w in (.1,.16):
 for lo,hi in [(4,10.4),(10,20),(10,24)]:
  by={L:{r['tau']:r for r in rows if r['L']==L and r['omega']==w and lo<=r['tau']<=hi and r['seeds']==3 and r['phase_resolved']} for L in (6,7)}
  common=sorted(set(by[6])&set(by[7]))
  if len(common)<3:continue
  for L in (6,7):
   for k in ('m2','chi'):fits.append(dict(L=L,omega=w,window=f'{lo:g}-{hi:g}',taus=','.join(f'{t:g}' for t in common),observable=k,**fit([by[L][t] for t in common],k)))
save_csv('matched-window-fits.csv',fits)
bg=[]
for r in rows:
 L,w,t=r['L'],r['omega'],r['tau'];b=w/(2*L**2);cf=1-math.exp(-w*t/2)
 bg.append(dict(L=L,omega=w,tau=t,m2=r['m2'],bare_m2_reference=b,bare_reference_over_m2=b/r['m2'],chi=r['chi'],bare_chi_reference=cf,bare_chi_reference_over_chi=cf/r['chi'],chi_over_tau_L2_m2=r['chi']/(t*L**2*r['m2'])))
save_csv('bare-background-reference.csv',bg)
shell=[]
kappa=np.arccos(-.3)/(2*np.pi);Q=np.array([2*kappa,1-2*kappa])
for L in (6,7):
 fs=[(abs(.6+np.exp(-2j*np.pi*a/L)+np.exp(-2j*np.pi*b/L)),a,b) for a in range(L) for b in range(L)];smin=min(f[0] for f in fs);modes=[(a,b) for s,a,b in fs if abs(s-smin)<1e-10]
 shell.append(dict(L=L,Nc=L**2,free_vacuum_energy=-sum(s for s,a,b in fs),sigma_min=smin,bare_single_quasiparticle_gap=2*smin,min_modes=modes,required_code_parity=1 if L==6 else -1,free_vacuum_code_parity=-1,trial_blocked=L==6,trial_blocking_energy_cost=2*smin if L==6 else 0,valley_transfer_indices=[[(a-c)%L,(b-d)%L] for a,b in modes for c,d in modes if (a,b)!=(c,d)],bare_intervalley_Q=Q.tolist()))
(OUT/'bare-shells.json').write_text(json.dumps(shell,indent=2))
if '--plots' in sys.argv:
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 for k,power in [('m2',2/3),('chi',5/3)]:
  fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
  for ax,w in zip(axes,(.1,.16)):
   for L in (6,7):
    rs=sorted([r for r in rows if r['L']==L and r['omega']==w and r[k]>0],key=lambda r:r['tau'])
    ax.errorbar([r['tau'] for r in rs],[r[k] for r in rs],yerr=[r[k+'_error'] or 0 for r in rs],fmt='o-',color='C0' if L==6 else 'C1',ms=3,capsize=2,label=f'L={L}')
    for r in rs:
     if r['seeds']<3:ax.plot(r['tau'],r[k],'o',mfc='white',mec='black',ms=5)
   anchor=next(r for r in rows if r['L']==7 and r['omega']==w and r['tau']==10.4);x=np.geomspace(4,24,100);ax.plot(x,anchor[k]*.6*(x/10.4)**power,'k--',label=f'slope guide {power:.3f} (not fit)')
   ax.set(xscale='log',yscale='log',xlabel='total tau',ylabel=k,title=f'Jz=.6, lambda=2, omega={w}');ax.grid(which='both',alpha=.2);ax.legend(fontsize=8)
  fig.suptitle('November 9 frozen snapshot; cut4000; errors: seed jackknife; open: <3 seeds')
  fig.savefig(D/'figures'/f'comparison-{k}.png',dpi=180);plt.close(fig)
print(f'Verified {len(rows)} points, {len(data["chains"])} independent selected chains; derived files written.')
