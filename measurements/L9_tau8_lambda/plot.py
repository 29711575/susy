import csv,gzip,json,hashlib,datetime
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tau8-cut500-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
chains=[];excluded=[];sources=[]
for port in [20015,20022]:
 f=P/f'order-snapshot-{port}.json.gz';d=json.load(gzip.open(f,'rt'));chains.extend(d['chains']);excluded.extend(d['excluded'])
 sources.append(dict(port=port,start=d['host_utc_started'],finish=d['host_utc_finished'],sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
assert len(chains)+len(excluded)==279
assert len({c['seed'] for c in chains+excluded})==279
assert all(c['reason']=='No retained samples after 500 sweep cutoff' for c in excluded),excluded
G=defaultdict(list)
for c in chains:
 assert c['tau']==8 and c['discard']==500 and c['count']==c['raw_count']-500>0
 a=np.array(c['chi_numerator_sums']).reshape(81,3,2);v=np.r_[complex(*c['phase_sum']),(a[:,:,0]+1j*a[:,:,1]).ravel()]
 assert np.isfinite(v).all();G[c['lambda_value']].append((c,v))
def spectrum(v):
 assert abs(v[0])>1e-12
 a=(v[1:]/v[0]).reshape(81,3)
 return a[:,2]-81*a[:,0]*a[:,1]
def metrics(v):
 s=spectrum(v)
 return np.array([s[1:].real.max(),s[49].real])
def jk(reps):
 return np.sqrt((len(reps)-1)/len(reps)*np.sum((reps-reps.mean(axis=0))**2,axis=0))
rows=[];grid=[];mergecheck=0
for lam,group in sorted(G.items()):
 n=len(group);assert n>=2
 vectors=np.array([v for c,v in group]);total=vectors.sum(axis=0);val=metrics(total);s=spectrum(total)
 replicas=np.array([metrics(total-v) for v in vectors]);errors=jk(replicas)
 specerr=jk(np.array([spectrum(total-v).real for v in vectors]))
 # Independently reconstruct pooled ratios from per-chain retained means and phase sums.
 phase=np.array([complex(*c['phase_sum']) for c,v in group]);means=np.array([v[1:]/v[0] for c,v in group])
 a=(np.sum(means*phase[:,None],axis=0)/phase.sum()).reshape(81,3)
 check=a[:,2]-81*a[:,0]*a[:,1];mergecheck=max(mergecheck,float(np.max(np.abs(check-s))))
 full=sum(c['status']=='complete' for c,v in group)
 rows.append(dict(L=9,omega=.3,tau=8,lambda_value=lam,chi_max=float(val[0]),chi_max_stderr=float(errors[0]),q_chi_max=int(1+s[1:].real.argmax()),chi_fixed_q49=float(val[1]),chi_fixed_q49_stderr=float(errors[1]),seeds=n,complete_seeds=full,point_complete=full==n,discard_sweeps=500,samples=sum(c['count'] for c,v in group),samples_min=min(c['count'] for c,v in group),samples_max=max(c['count'] for c,v in group)))
 for q in range(81):grid.append(dict(lambda_value=lam,q=q,n1=q//9,n2=q%9,chi_real=float(s[q].real),chi_imag=float(s[q].imag),chi_stderr_real=float(specerr[q])))
assert mergecheck<1e-10
for filename,data in [('chi-vs-lambda-tau8-cut500.csv',rows),('allowed-grid.csv',grid),('seed-coverage.csv',[{k:c[k] for k in ['batch','seed','lambda_value','tau','status','raw_count','discard','count','blocks','root']} for c in chains])]:
 with open(P/filename,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#555555','text.color':'#252525','axes.labelcolor':'#252525','savefig.facecolor':'white'})
for zoom in [False,True]:
 rr=[r for r in rows if not zoom or r['lambda_value']<=3.1]
 fig,ax=plt.subplots(figsize=(11.5,6.8));x=np.array([r['lambda_value'] for r in rr]);y=np.array([r['chi_max'] for r in rr]);err=np.array([r['chi_max_stderr'] for r in rr])
 ax.plot(x,y,color='#2463a0',lw=1.25)
 for done in [False,True]:
  mask=np.array([r['point_complete']==done for r in rr])
  ax.errorbar(x[mask],y[mask],yerr=err[mask],fmt='o',ms=6,mfc='#2463a0' if done else 'white',mec='#2463a0',color='#2463a0',elinewidth=1,capsize=3,ls='none')
 ax.set_xlabel(r'Coupling $\lambda$');ax.set_ylabel(r'$\max_{\mathbf{q}\ne0}\mathrm{Re}\,\chi_{zz}(\mathbf{q})$')
 ax.set_title(r'$L=9,\ \omega=0.3,\ \tau=8$ | VBS static susceptibility',loc='left',fontsize=16,pad=25)
 ax.text(0,1.025,'Discard sweeps 1-500 in every chain; latest complete blocks collected from HPC',transform=ax.transAxes,fontsize=10,color='#666666')
 ax.grid(axis='y',color='#dedede',lw=.6);ax.set_axisbelow(True)
 ax.set_xlim(1.22,3.15 if zoom else 5.08)
 ax.set_xticks([1.3,1.6,1.9,2.2,2.4,2.6,2.8,3,3.1] if zoom else [1.3,1.6,1.9,2.2,2.4,2.6,2.8,3,3.2,3.6,4,4.5,5])
 ax.axhline(0,color='#888888',lw=.6)
 fig.text(.5,.065,'Pool phase-weighted raw moments, then subtract the connected background; maximum over 80 nonzero allowed q.\nError bars: delete-one-independent-seed jackknife (1 SE), reselecting the momentum maximum in every replica.\nOpen: at least one partial chain. Original and supplemental production batches included; lines join sampled points.',ha='center',fontsize=9,color='#555555',linespacing=1.6)
 fig.subplots_adjust(left=.105,right=.985,bottom=.215,top=.855)
 fig.savefig(P/('chi-vs-lambda-tau8-cut500'+('-zoom' if zoom else '')+'.png'),dpi=220);plt.close(fig)
metadata=dict(L=9,omega=.3,tau=8,planned_chains=279,retained_chains=len(chains),excluded_chains=len(excluded),exclusions=excluded,lambda_values=[r['lambda_value'] for r in rows],cells=len(rows),discard=500,total_retained_samples=sum(c['count'] for c in chains),seed_range=[min(r['seeds'] for r in rows),max(r['seeds'] for r in rows)],sources=sources,batch_counts=dict(Counter(c['batch'] for c in chains)),definition='max over 80 allowed nonzero q of Re chi_zz(q), chi=<ChiRaw>-81*<ChiD>*<ChiIconj>; separate conjugate moment preserved, no free-phonon subtraction',error='joint phase-reweighted delete-one-independent-seed jackknife1SE; connected subtraction and momentum maximum reevaluated in every replica',snapshot='Latest read-only snapshots of all original 279 planned chains; retain complete blocks after sweep500, without the old per-chain sample caps.',checks=dict(independent_pooled_ratio_reconstruction_max_residual=mergecheck,unique_chains=True,complete_10_sample_blocks=True),limitations='Includes stopped partial production chains when retained samples exist; no thermalization, autocorrelation or convergence certification.')
(P/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n');(P/'README.txt').write_text(json.dumps(metadata,indent=2)+'\nReproduce frozen snapshot: python3 plot.py\n')
print(json.dumps({k:v for k,v in metadata.items() if k not in ['exclusions','batch_counts','sources']},indent=2))
