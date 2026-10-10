import csv,gzip,json,hashlib,os
from pathlib import Path
from collections import defaultdict
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/chi-full-disc-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
specs=[('tau_scan',P.parent/'L9_L6_tau_scan'),('tau8_lambda_scan',P.parent/'L9_tau8_lambda')]
rows=[];grids=[];seedrows=[];sources=[];identity=0.;reference_residual=0.
def parts(v,N):
 a=(v[1:]/v[0]).reshape(N,3)
 full=a[:,2].real;disc=(N*a[:,0]*a[:,1]).real
 return np.array([full,disc,full-disc])
def error(a):
 n=len(a);return np.sqrt((n-1)/n*np.sum((a-a.mean(axis=0))**2,axis=0))
for campaign,path in specs:
 cs=[]
 for port in [20015,20022]:
  f=path/f'order-snapshot-{port}.json.gz';d=json.load(gzip.open(f,'rt'));cs+=d['chains']
  sources.append(dict(campaign=campaign,path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),start=d['host_utc_started'],finish=d['host_utc_finished']))
 refname='grid-chi.csv' if campaign=='tau_scan' else 'allowed-grid.csv'
 ref={(int(r.get('L',9)),float(r['lambda_value']),float(r.get('tau',8)),int(r['q'])):float(r['chi_real']) for r in csv.DictReader((path/refname).open())}
 groups=defaultdict(list)
 for c in cs:groups[(c.get('L',9),c['lambda_value'],c['tau'])].append(c)
 for (L,lam,tau),group in sorted(groups.items()):
  N=L*L;n=len(group);assert n>=2
  v=[]
  for c in group:
   a=np.array(c['chi_numerator_sums']).reshape(N,3,2)
   v.append(np.r_[complex(*c['phase_sum']),(a[:,:,0]+1j*a[:,:,1]).ravel()])
  v=np.array(v);total=v.sum(axis=0);p=parts(total,N);jack=np.array([parts(total-x,N) for x in v]);se=error(jack)
  # Independent full complex connected formula agrees with separately-real decomposition.
  a=(total[1:]/total[0]).reshape(N,3)
  identity=max(identity,float(np.max(np.abs(p[2]-(a[:,2]-N*a[:,0]*a[:,1]).real))))
  for q in range(N):
   reference_residual=max(reference_residual,abs(p[2,q]-ref[(L,lam,tau,q)]))
   grids.append(dict(campaign=campaign,L=L,lambda_value=lam,tau=tau,q=q,n1=q//L,n2=q%L,C_full=float(p[0,q]),C_full_stderr=float(se[0,q]),C_disc=float(p[1,q]),C_disc_stderr=float(se[1,q]),chi=float(p[2,q]),chi_stderr=float(se[2,q])))
  qfixed=49 if L==9 else 21;qpeak=int(1+p[2,1:].argmax());qraw=int(1+p[0,1:].argmax())
  for selection,q in [('fixed_target',qfixed),('connected_peak',qpeak)]:
   ratio=p[1,q]/p[0,q];ratiojack=jack[:,1,q]/jack[:,0,q]
   rows.append(dict(campaign=campaign,L=L,lambda_value=lam,tau=tau,selection=selection,q=q,n1=q//L,n2=q%L,seeds=n,retained_sweeps=sum(c['count'] for c in group),discard_sweeps=group[0]['discard'],C_full=float(p[0,q]),C_full_stderr=float(se[0,q]),C_disc=float(p[1,q]),C_disc_stderr=float(se[1,q]),chi=float(p[2,q]),chi_stderr=float(se[2,q]),disc_percent_of_full=float(100*ratio),disc_percent_stderr=float(100*error(ratiojack)),raw_peak_q=qraw,connected_peak_q=qpeak))
   for c,x in zip(group,v):
    a=(x[1:]/x[0]).reshape(N,3);seedpart=parts(x,N)[:,q];weight=x[0]/total[0]
    seedrows.append(dict(campaign=campaign,L=L,lambda_value=lam,tau=tau,selection=selection,q=q,seed=c['seed'],retained_sweeps=c['count'],mean_D_real=float(a[q,0].real),mean_D_imag=float(a[q,0].imag),mean_Iconj_real=float(a[q,1].real),mean_Iconj_imag=float(a[q,1].imag),C_full_seed=float(seedpart[0]),C_disc_seed=float(seedpart[1]),chi_seed=float(seedpart[2]),complex_phase_weight_real=float(weight.real),complex_phase_weight_imag=float(weight.imag)))
assert identity<1e-12 and reference_residual<1e-10
for name,data in [('decomposition.csv',rows),('all-q-decomposition.csv',grids),('seed-first-moments.csv',seedrows)]:
 with (P/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
blue='#2463a0';orange='#d58b1c';purple='#a65087'
def pair(ax,rr,x):
 for y,c,style,label in [('C_full',blue,'o--',r'$C_{\rm full}$'),('chi',orange,'s-',r'$\chi=C_{\rm full}-C_{\rm disc}$')]:
  ax.errorbar([r[x] for r in rr],[r[y] for r in rr],yerr=[r[y+'_stderr'] for r in rr],fmt=style,color=c,ms=4,capsize=2,lw=1.2,mfc='white',label=label)
 ax.grid(axis='y',color='#e1e1e1',lw=.6);ax.set_ylabel('Response')
def disc(ax,rr,x):
 ax.errorbar([r[x] for r in rr],[r['C_disc'] for r in rr],yerr=[r['C_disc_stderr'] for r in rr],fmt='^-',color=purple,ms=4,capsize=2,lw=1.2,label=r'$C_{\rm disc}$')
 ax.axhline(0,color='#888888',lw=.6);ax.grid(axis='y',color='#e1e1e1',lw=.6);ax.set_ylabel(r'$C_{\rm disc}$')
fig,axes=plt.subplots(2,5,figsize=(19,8),layout='constrained')
for col,lam in enumerate([1.3,2,2.4,2.6,3]):
 rr=sorted([r for r in rows if r['campaign']=='tau_scan' and r['L']==9 and r['selection']=='fixed_target' and r['lambda_value']==lam],key=lambda r:r['tau'])
 pair(axes[0,col],rr,'tau');disc(axes[1,col],rr,'tau');axes[0,col].set_title(r'$\lambda='+str(lam)+'$');axes[1,col].set_xlabel(r'$\tau$')
 axes[0,col].set_ylim(bottom=0);axes[0,col].set_xticks([0,4,8,12]);axes[1,col].set_xticks([0,4,8,12])
axes[0,0].legend(fontsize=9,frameon=False);axes[1,0].legend(frameon=False)
fig.suptitle(r'$L=9,\ \omega=0.3$ | fixed target $q=(5/9,4/9)$',fontsize=17)
fig.supxlabel('Pooled complex-phase moments, then product of ensemble means; jackknife 1 SE. Cut1000 for tau>6, otherwise no cut.',fontsize=11)
fig.savefig(P/'L9-full-disc-vs-tau-fixed-q49.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for col,selection in enumerate(['fixed_target','connected_peak']):
 rr=sorted([r for r in rows if r['campaign']=='tau8_lambda_scan' and r['selection']==selection],key=lambda r:r['lambda_value'])
 pair(axes[0,col],rr,'lambda_value');disc(axes[1,col],rr,'lambda_value');axes[1,col].set_xlabel(r'$\lambda$')
 axes[0,col].set_title('Fixed q=(5/9,4/9)' if selection=='fixed_target' else 'q of the pooled connected maximum')
axes[0,0].legend(frameon=False);axes[1,0].legend(frameon=False)
fig.suptitle(r'$L=9,\ \omega=0.3,\ \tau=8$ | full and disconnected response',fontsize=17)
fig.supxlabel('Cut500 per chain. Each peak column uses the same q for full, disc, and chi; q may change across lambda. Jackknife 1 SE at selected q.',fontsize=10)
fig.savefig(P/'L9-tau8-full-disc-vs-lambda.png',dpi=200);plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(13,12),layout='constrained')
for col,tau in enumerate([10,14]):
 for idx,selection in enumerate(['fixed_target','connected_peak']):
  rr=sorted([r for r in rows if r['campaign']=='tau_scan' and r['L']==6 and r['tau']==tau and r['selection']==selection],key=lambda r:r['lambda_value'])
  pair(axes[idx,col],rr,'lambda_value');axes[idx,col].set_title(f'tau={tau}: '+('fixed q=(1/2,1/2)' if idx==0 else 'q of connected maximum'))
  axes[2,col].errorbar([r['lambda_value'] for r in rr],[r['C_disc'] for r in rr],yerr=[r['C_disc_stderr'] for r in rr],fmt='o-' if idx==0 else 's--',color=purple if idx==0 else blue,ms=4,capsize=2,label='fixed q' if idx==0 else 'connected-peak q')
 axes[2,col].axhline(0,color='#888888',lw=.6);axes[2,col].set_ylabel(r'$C_{\rm disc}$');axes[2,col].set_xlabel(r'$\lambda$');axes[2,col].legend(frameon=False)
axes[0,0].legend(frameon=False)
fig.suptitle(r'$L=6,\ \omega=0.3$ | full and disconnected response',fontsize=17)
fig.supxlabel('Cut3000 per chain. Pooled raw moments before connected subtraction; 1 SE at selected q. Completed original and fine-scan chains.',fontsize=10)
fig.savefig(P/'L6-full-disc-by-tau.png',dpi=180);plt.close(fig)
summary={}
for campaign,L,taus in [('tau_scan',9,[8,10,12,14]),('tau_scan',6,[10,14]),('tau8_lambda_scan',9,[8])]:
 for selection in ['fixed_target','connected_peak']:
  rr=[r for r in rows if r['campaign']==campaign and r['L']==L and r['tau'] in taus and r['selection']==selection]
  fractions=np.array([abs(r['disc_percent_of_full']) for r in rr]);largest=max(rr,key=lambda r:abs(r['disc_percent_of_full']))
  summary[f'{campaign}_L{L}_{selection}']=dict(cells=len(rr),median_absolute_disc_percent=float(np.median(fractions)),max_absolute_disc_percent=float(max(fractions)),largest=largest)
report=dict(sources=sources,cells=len(rows)//2,identity_max_residual=identity,prior_plotted_chi_max_residual=reference_residual,definition='Cfull=Re(<ChiRaw>)=Nc Re(<D Iconj>); Cdisc=Nc Re(<D><Iconj>); chi=Cfull-Cdisc. Every ensemble mean uses pooled complex-phase sums.',fixed_targets={'L9':'n=(5,4), q=(5/9,4/9)','L6':'n=(3,3), q=(1/2,1/2)'},error='Delete-one-seed jackknife, recompute pooled means and their product. For decomposition at selected connected-peak q, q is held fixed across replicas to keep the additive identity; prior peak curves recomputed their maximum per replicate.',summary=summary,limitations='Frozen current plotted snapshots, not a fresh HPC collection; differing cutoff and chain scopes kept as separate campaigns. Small pooled first moments may reflect sector cancellation; does not establish equilibration or a physical low-energy mechanism. No autocorrelation estimate from aggregates. Cdisc need not be clipped positive for separately complex-phase-reweighted Iconj.')
(P/'metadata.json').write_text(json.dumps(report,indent=2)+'\n');(P/'README.txt').write_text(json.dumps(report,indent=2)+'\nReproduce: python3 analyze.py\n')
print(json.dumps({'identity_residual':identity,'plot_residual':reference_residual,'summary':summary},indent=2))
