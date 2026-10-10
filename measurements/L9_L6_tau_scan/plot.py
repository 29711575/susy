import os,csv,gzip,json,hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/chi-tau-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
chains=[];excluded=[];sources=[]
for port in [20015,20022]:
 f=P/f'order-snapshot-{port}.json.gz';d=json.load(gzip.open(f,'rt'))
 chains.extend(d['chains']);excluded.extend(d['excluded'])
 sources.append(dict(port=port,start=d['host_utc_started'],finish=d['host_utc_finished'],sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
assert len({c['seed'] for c in chains})==len(chains)
assert all('No retained samples after cutoff' in e['reason'] for e in excluded),excluded
assert len([c for c in chains if c['batch']=='L9-original'])==340
assert len([c for c in chains if c['batch']=='L6-original'])==40

def vector(c):
 a=np.array(c['chi_numerator_sums']).reshape(c['L']**2,3,2)
 return np.r_[complex(*c['phase_sum']),(a[:,:,0]+1j*a[:,:,1]).ravel()]
def spectrum(v,Nc):
 assert abs(v[0])>1e-12
 a=(v[1:]/v[0]).reshape(Nc,3)
 return a[:,2]-Nc*a[:,0]*a[:,1]
def summarize(cs,scope):
 groups=defaultdict(list)
 for c in cs:groups[(c['L'],c['lambda_value'],c['tau'])].append(c)
 rows=[];grid=[];residual=0.
 for (L,lam,tau),group in sorted(groups.items()):
  Nc=L*L;n=len(group);assert n>=2
  v=np.array([vector(c) for c in group]);total=v.sum(axis=0);chi=spectrum(total,Nc);rr=chi.real
  jack=np.array([spectrum(total-x,Nc).real for x in v]);maxima=jack[:,1:].max(axis=1)
  se=np.sqrt((n-1)/n*np.sum((jack-jack.mean(axis=0))**2,axis=0))
  peak_se=np.sqrt((n-1)/n*np.sum((maxima-maxima.mean())**2))
  q=1+rr[1:].argmax();fixed=49 if L==9 else 21
  # Independent expression in terms of each seed's phase-weighted mean.
  counts=np.array([c['count'] for c in group]);phase=v[:,0];weights=counts*(phase/counts)
  means=(v[:,1:]/phase[:,None]).reshape(n,Nc,3)
  pooled=np.sum(weights[:,None,None]*means,axis=0)/weights.sum()
  independent=pooled[:,2]-Nc*pooled[:,0]*pooled[:,1]
  residual=max(residual,float(np.max(np.abs(chi-independent))))
  row=dict(scope=scope,L=L,lambda_value=lam,tau=tau,seeds=n,complete_seeds=sum(c['status']=='complete' for c in group),samples=sum(c['count'] for c in group),samples_min=min(c['count'] for c in group),samples_max=max(c['count'] for c in group),discard_sweeps=group[0]['discard'],chi_max=float(rr[q]),chi_max_stderr=float(peak_se),q_max=int(q),chi_fixed=float(rr[fixed]),chi_fixed_stderr=float(se[fixed]),q_fixed=fixed,batches=';'.join(sorted({c['batch'] for c in group})))
  assert all(c['discard']==row['discard_sweeps'] for c in group)
  rows.append(row)
  for q in range(Nc):grid.append(dict(scope=scope,L=L,lambda_value=lam,tau=tau,q=q,n1=q//L,n2=q%L,chi_real=float(rr[q]),chi_imag=float(chi[q].imag),stderr_real=float(se[q])))
 assert residual<1e-10
 return rows,grid,residual
rows,grid,residual=summarize(chains,'pooled')
original,originalgrid,originalresidual=summarize([c for c in chains if c['batch']=='L9-original'],'L9-original')
assert len([r for r in rows if r['L']==9])==85

def savecsv(name,data):
 with (P/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
savecsv('curves.csv',rows);savecsv('grid-chi.csv',grid);savecsv('L9-original-curves.csv',original)
savecsv('seed-coverage.csv',[{k:c[k] for k in ['L','batch','seed','lambda_value','tau','raw_count','discard','count','sweeps','status','run']} for c in chains])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False})
colors=['#2463a0','#d58b1c','#bd4b35','#728337','#a65087'];markers=['o','s','^','D','v']
outputs=[]
def finish(fig,ax,name,title,xlabel,ylabel,note,log=False):
 ax.set_title(title,loc='left',fontsize=15,pad=15);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
 ax.grid(axis='y',lw=.6,color='#dedede');ax.set_axisbelow(True)
 if log:
  ax.set_xscale('log');ax.set_yscale('log');ax.set_xticks([.2,.4,1,2,4,8,14]);ax.set_xticklabels(['0.2','0.4','1','2','4','8','14'])
 fig.text(.5,.055,note,ha='center',fontsize=9,color='#555555',linespacing=1.5)
 fig.subplots_adjust(left=.12,right=.98,bottom=.22,top=.88)
 fig.savefig(P/name,dpi=220,facecolor='white');plt.close(fig);outputs.append(name)
for fixed in [False,True]:
 for log in [False,True]:
  fig,ax=plt.subplots(figsize=(10.4,6.7));lambdas=sorted({r['lambda_value'] for r in rows if r['L']==9})
  for lam,col,mark in zip(lambdas,colors,markers):
   rr=sorted([r for r in rows if r['L']==9 and r['lambda_value']==lam],key=lambda r:r['tau'])
   y='chi_fixed' if fixed else 'chi_max';e=y+'_stderr'
   ax.errorbar([r['tau'] for r in rr],[r[y] for r in rr],yerr=[r[e] for r in rr],fmt=mark+'-',color=col,ms=5,lw=1.3,capsize=2,mfc='white',label=rf'$\lambda={lam:g}$')
  ax.legend(ncol=5,frameon=False,fontsize=11,loc='upper left')
  if log:ax.set_xticks([.2,.4,1,2,4,8,14]);ax.set_xticklabels(['0.2','0.4','1','2','4','8','14'])
  label=r'$\mathrm{Re}\,\chi_{zz}(q=(5/9,4/9))$' if fixed else r'$\max_{q\ne0}\mathrm{Re}\,\chi_{zz}(q)$'
  finish(fig,ax,'chi-vs-tau-L9'+('-fixed-q49' if fixed else '')+('-loglog' if log else '')+'.png',r'$L=9,\ \omega=0.3$ | connected susceptibility',r'Total imaginary time $\tau$',label,'Phase-weighted raw moments pooled before connected subtraction; jackknife 1 SE.\n4 seeds per point; tau=10,12,14 include up to 6 supplemental seeds.\nDiscard 1000 sweeps per chain for tau>6; no cut for tau<=6. Finite-projection estimates.',log)
for tau in sorted({r['tau'] for r in rows if r['L']==6}):
 fig,ax=plt.subplots(figsize=(9,6.5));rr=sorted([r for r in rows if r['L']==6 and r['tau']==tau],key=lambda r:r['lambda_value'])
 for y,col,mark,label in [('chi_max',colors[0],'o',r'$\max_{q\ne0}\mathrm{Re}\,\chi_{zz}(q)$'),('chi_fixed',colors[1],'s',r'$\mathrm{Re}\,\chi_{zz}(q=(1/2,1/2))$')]:
  ax.errorbar([r['lambda_value'] for r in rr],[r[y] for r in rr],yerr=[r[y+'_stderr'] for r in rr],fmt=mark+'-',color=col,lw=1.4,capsize=3,mfc='white',label=label)
 ax.legend(frameon=False,fontsize=11);ax.set_xlim(1.2,3.1)
 finish(fig,ax,f'chi-vs-lambda-L6-tau{tau:g}.png',rf'$L=6,\ \omega=0.3,\ \tau={tau:g}$ | connected susceptibility',r'$\lambda$',r'$\mathrm{Re}\,\chi_{zz}$','Discard first 3000 sweeps per chain; phase-weighted pooled raw moments.\nOriginal points: 4 seeds each; supplemental tau10 fine scan included where retained data exist.\nError bars: delete-one-seed jackknife 1 SE; finite-projection estimates.')
fig,ax=plt.subplots(figsize=(9,6.5))
for tau,col,mark in zip(sorted({r['tau'] for r in rows if r['L']==6}),colors,markers):
 rr=sorted([r for r in rows if r['L']==6 and r['tau']==tau],key=lambda r:r['lambda_value'])
 ax.errorbar([r['lambda_value'] for r in rr],[r['chi_max'] for r in rr],yerr=[r['chi_max_stderr'] for r in rr],fmt=mark+'-',color=col,lw=1.4,capsize=3,mfc='white',label=rf'$\tau={tau:g}$')
ax.legend(frameon=False)
finish(fig,ax,'chi-vs-lambda-L6-all-tau.png',r'$L=6,\ \omega=0.3$ | connected susceptibility',r'$\lambda$',r'$\max_{q\ne0}\mathrm{Re}\,\chi_{zz}(q)$','All tau curves use the same normalization and discard first 3000 sweeps per chain.\nRaw moments pooled with a common complex-phase denominator; jackknife 1 SE.')
report=dict(sources=sources,chains=len(chains),cells=len(rows),excluded=excluded,retained_sweeps=sum(c['count'] for c in chains),cuts='L9 tau>6:1000, tau<=6:0; L6:3000',definition='chi(q)=<ChiRaw(q)>-Nc<ChiD(q)><ChiIconj(q)>; pool raw sums and common phase before connected subtraction and q maximum; no free-phonon subtraction',error='delete-one-independent-seed jackknife1SE; pool and connected subtraction and peak recomputed',weighted_reconstruction_max_residual=residual,original_scan_weighted_residual=originalresidual,outputs=outputs,scope='OriginalL9 tau scan340 plus add6tau10/12/14 90; originalL6tau10/14 40 plus fineL6tau10 30; controls and later algorithms excluded',limitations='Incomplete/stopped chain data included where past cutoff; no certification of thermalization, autocorrelation or ground-state limit; nonzero-q maxima may have selection bias')
(P/'metadata.json').write_text(json.dumps(report,indent=2)+'\n');(P/'README.txt').write_text(json.dumps(report,indent=2)+'\nReproduce frozen snapshot: python3 plot.py\n')
print(json.dumps(report,indent=2))
