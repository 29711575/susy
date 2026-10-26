import os,csv,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/pcn-grid-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
P=Path.cwd()/'analysis'
def read(file):
 with open(file)as f:return list(csv.DictReader(f))
rows=read(P/'observables.csv');fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for ax,n in zip(axes,['X_z','BondXK_z','B4_q3_z']):
 for alg in ['pcn','guided','local']:
  r=sorted([r for r in rows if r['stage']=='grid'and r['algorithm']==alg and r['observable']==n],key=lambda x:float(x['dt']));x=[float(v['dt'])for v in r];y=[float(v['estimated_CPU_seconds_to_common_error'])if v['estimated_CPU_seconds_to_common_error']else np.nan for v in r];ax.plot(x,y,'o-',label=alg)
 ax.set_title(n);ax.set_xlabel('dt; fixed physical tau=2');ax.set_ylabel('Estimated CPU seconds to common error');ax.set_yscale('log');ax.grid(alpha=.2)
axes[0].legend();fig.suptitle('Time estimates require resolved blocking/phase diagnostics; blank points are inconclusive');fig.savefig(P/'grid-common-error-time.png',dpi=160);plt.close(fig)
r=read(P/'projection-convergence.csv');fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
for n in ['majorana_energy','X_z','X2_z','BondK_z','BondXK_z','S_q3_zz','B4_q3_z']:
 v=[v for v in r if v['observable']==n];ax.plot([float(v['tau'])for v in v],[max(1e-8,abs(float(v['finite_projection_bias']))/float(v['ground_convergence_tolerance']))for v in v],'o-',label=n)
ax.axhline(1,color='black',ls='--');ax.set_yscale('log');ax.set_xlabel('Total projection tau');ax.set_ylabel('|continuous ED - same-parity ground| / declared tolerance');ax.set_title('Physical projection bias; separate from dt and MC error');ax.legend(fontsize=8);ax.grid(alpha=.2);fig.savefig(P/'projection-bias.png',dpi=160);plt.close(fig)
r=read(P/'pcn/observables.csv');fig,axes=plt.subplots(3,2,figsize=(10,10),layout='constrained')
for i,beta in enumerate([.3,.6,1.]):
 for col,source in enumerate(['center_F','center_F2_minus_variance']):
  v=[v for v in r if float(v['beta'])==beta and v['source']==source and v['observable'].startswith('rho_lag')];lags=[int(v['observable'][7:])for v in v];ax=axes[i,col];ax.errorbar(lags,[float(v['value'])for v in v],yerr=[float(v['MC_SE'])for v in v],fmt='o',capsize=2,label='production 1SE');k=np.arange(1,129);a=np.sqrt(1-beta*beta);ax.plot(k,a**(k*(1 if col==0 else 2)),label='analytic');ax.set_xscale('log',base=2);ax.set_xlabel('Lag in selected-path updates');ax.set_title(f'beta={beta}, {source}');ax.grid(alpha=.2)
axes[0,0].legend();fig.savefig(P/'pcn-analytic-acf.png',dpi=160);plt.close(fig)
