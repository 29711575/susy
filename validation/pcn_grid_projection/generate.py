"""Freeze pilot/pure-pCN definitions; formal budgets determined once by planner."""
import argparse,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();out=a.directory;out.mkdir(parents=True,exist_ok=True)
BASE=['--Lx','2','--Ly','2','--shear','0','--Jx','1','--Jy','1','--Jz','.5','--trial-dx','0','--trial-dy','0','--trial-dz','0','--measure-every','1','--qr-block','4','--q-mode-every','0','--q-mode-n1','1','--q-mode-n2','1','--q-mode-beta','1','--q-mode-bonds','z','--numerics','production','--weight-engine','combined','--defer-phase','1','--prefix-cache','1','--cayley-cache','1','--center-green-fast','1','--acceptance-every','1024']
def arguments(case,seed,warmup,sweeps,pcn=False):
 v=BASE+['--g',str(0 if pcn else .4),'--mass',str(100/9 if pcn else 1),'--omega-ph',str(.3 if pcn else 1),'--intervals',str(round(case['tau']/case['dt'])),'--dt',str(case['dt']),'--seed',str(seed),'--warmup',str(warmup),'--sweeps',str(sweeps),'--block-samples','32','--fermions','0','--phonon-diagnostics','0','--measure-energy',str(int(not pcn)),'--bond-diagnostics',str(int(not pcn)),'--vbs',str(int(not pcn)),'--vbs-amplitude',str(int(not pcn)),'--vbs-susceptibility','0','--vbs-amplitude-response','0','--validation-series',str(int(not pcn))]
 if not pcn:v+=['--vbs-mode','all']
 else:v+=['--path-series-bond','2']
 alg=case['algorithm'];v+=['--proposal','local'if alg=='local'else'bond-guided']
 if alg!='local':v+=['--bond-beta',str(case.get('beta',.6)),'--bond-guide','bare'if alg=='guided'else'zero']
 return v
cases=[dict(name=f'grid-{alg}-dt{dt}',stage='grid',tau=2,dt=dt,algorithm=alg,beta=.6)for alg in ['pcn','guided','local']for dt in [.2,.1,.05]]
cases += [dict(name=f'projection-t{tau}-dt{dt}',stage='projection',tau=tau,dt=dt,algorithm='guided',beta=.6)for tau in [4,8,16]for dt in [.1,.05]]
cases.sort(key=lambda c:-(c['tau']/c['dt']*(3 if c['algorithm']=='local'else 1)))
pilot=[]
for i,c in enumerate(cases):
 for j in range(2):
  seed=202610191001+10*i+j;pilot.append(dict(index=len(pilot),case=c,seed=seed,sweeps=8192,warmup=1024,args=arguments(c,seed,1024,8192)))
pcn=[]
for beta in [.3,.6,1.]:
 c=dict(name=f'pure-pcn-beta{beta}',stage='pcn',tau=2,dt=.1,algorithm='pcn',beta=beta)
 for j in range(4):
  seed=202610193001+10*len(pcn)+j;pcn.append(dict(index=len(pcn),case=c,seed=seed,sweeps=65536,warmup=512,args=arguments(c,seed,512,65536,True)))
for name,jobs in [('pilot',pilot),('pcn',pcn)]:
 f=out/(name+'-plan.json');assert not f.exists();f.write_text(json.dumps(dict(tier=name,jobs=jobs),indent=2)+'\n')
config=dict(physics=dict(L=2,J=[1,1,.5],g=.4,mass=1,omega=1,k=1,trial_shift=[0,0,0],gauge='uAB+1,Wilson+1+1',tau_definition='intervals*dt, two-sided total contour'),cases=cases,pilot=dict(seeds='explicit pilot-plan.json; disjoint formal seeds',warmup=1024,sweeps=8192,replicas=2),formal_budget_rule=dict(replicas=4,minimum_sweeps_per_chain=32768,maximum_sweeps_per_chain=262144,maximum_CPU_seconds_per_chain=14400,variance_safety_factor=2,unresolved_pilot_safety_factor=4,warmup_minimum=2048,warmup_maximum=32768,warmup_integrated_time_multiplier=20,round_sweeps='powers of two, capped by measured runtime; no formal optional stopping',formal_seeds='202610192001+10*case_index+replica'),gates=dict(sigma=6,absolute_floor=2e-9,ED_cutoff_absolute_relative=2e-5,ED_B4_cutoff_absolute=.0002,ED_eigen_residual=1e-8,ED_propagation_residual=1e-9,minimum_blocks_per_seed=20,block_sizes=[32,64,128,256,512,1024,2048,4096,8192],plateau_ratio=[.7,1.3],minimum_block_in_tau_int=10,split_Rhat_maximum=1.05,phase_resolution_sigma=5),pcn_gates=dict(absolute_reference_tolerance=2e-9,sigma=6,acf_absolute_SE_maximum=.02,IAT_relative_SE_maximum=.10,covariance_diagonal_relative_SE=.02,endpoint_kurtosis_SE=.05,jackknife_block_sweeps=2048,acf_lags=[1,2,4,8,16,32,64,128],IAT_definition='tau_int=1/2+sum_{k>=1}rho(k), variance inflation=2*tau_int; measurement every sweep, one selected-path pCN per sweep',IAT_window='fixed from analytic a: max(16,ceil(log(1e-5)/log(a))); beta1 window16'),projection_ground_gates=dict(matter_energy_absolute=.01,mean_or_bond_absolute=.002,second_relative=.01,fourth_relative=.02,B4_absolute=.02,near_zero_absolute=.002),scope='No default sampling policy/Hamiltonian/trial/measurement definition changes. Pure controls disable q explicitly. All original assessments preserved.')
(out/'config.json').write_text(json.dumps(config,indent=2)+'\n');print('pilot',len(pilot),'pcn',len(pcn))
