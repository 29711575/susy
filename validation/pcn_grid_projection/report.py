"""Finite-dt ED compatibility, efficiency and separate projection convergence."""
import csv,json,math,time
from pathlib import Path
import numpy as np
from moments import RAW,NAMES,derive,read_series,grouped,jackknife,influences,acf,rhat,target
from analyze_pcn import analyze as analyze_pcn
P=Path.cwd();OUT=P/'analysis';OUT.mkdir(exist_ok=True);cfg=json.loads((P/'config.json').read_text());g=cfg['gates'];ED=Path('/home/sunxr/quick1010-pcn-projection-ed-20261019-20022/results')
def write(file,rows):
 if rows:
  with (OUT/file).open('w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def discrete(c,cut):
 name=f"g0.4-tau{c['tau']}-dt{c['dt']}-cut{cut}.json"
 for p in [ED/name,P/'references'/name]:
  if p.exists():return json.loads(p.read_text())
 raise FileNotFoundError(name)
def continuous(cut,tau):return next(r for r in json.loads((ED/f'continuous-cut{cut}.json').read_text())['results']if r['tau']==tau)
def cutoff(values,name):
 coarse=abs(values[1][name]-values[0][name]);fine=abs(values[2][name]-values[1][name]);tol=g['ED_B4_cutoff_absolute']if name.startswith('B4')else g['ED_cutoff_absolute_relative']*max(1,abs(values[2][name]));return fine<=tol and (fine<=coarse or fine<1e-10),2*fine
pcn=analyze_pcn(P);plan=json.loads((P/'formal-plan.json').read_text());rows=[];curves=[];chains_report=[];autocorr=[];cases=[]
for name in dict.fromkeys(r['case']['name']for r in plan['jobs']):
 jobs=[r for r in plan['jobs']if r['case']['name']==name];c=jobs[0]['case'];aa=[];js=[];problems=[]
 for r in jobs:
  pre=P/'formal/runs'/f"chain{r['index']:03d}";j=None
  try:
   j=json.loads(Path(str(pre)+'_results.json').read_text());assert j['status']=='complete'and j['numerical_failures']==0,(j['status'],j['numerical_failures'],j.get('error'));assert j['seed']==r['seed']and j['samples']==r['sweeps'];assert j['g']==.4 and j['mass']==1 and j['omega_ph']==1 and j['trial_shift']==[0,0,0]and j['wilson_loops']==[1,1];assert abs(j['dt']-c['dt'])<1e-14 and j['intervals']*j['dt']==c['tau'];assert j['q_mode']['every']==0
   assert [j[k]for k in ['weight_engine','defer_phase','prefix_cache','cayley_cache','center_green_fast']]==['combined',1,1,1,1]
   a=read_series(str(pre)+'_validation_series.csv');assert len(a)==r['sweeps'];energy_identity=a[:,2+RAW.index('majorana_energy')]-4*(a[:,2+RAW.index('BondK_x')]+a[:,2+RAW.index('BondK_y')]+.5*a[:,2+RAW.index('BondK_z')]+.4*(a[:,2+RAW.index('BondXK_x')]+a[:,2+RAW.index('BondXK_y')]+a[:,2+RAW.index('BondXK_z')]));identity_error=float(np.max(abs(energy_identity)));assert identity_error<2e-9; t=a.sum(axis=0);err=max(abs(t[k+2]/t[0]-complex(*j['observables'][n]['value']))for k,n in enumerate(RAW));assert err<2e-10
   aa.append(a);js.append(j);chains_report.append(dict(case=name,index=r['index'],seed=r['seed'],samples=len(a),status='COMPLETE',numerical_failures=j['numerical_failures'],raw_reconstruction_error=float(err),matter_energy_BondK_XK_identity_max_error=identity_error,proposal=j['proposal'],guided_bond=j['guided_bond'],numerics={k:j[k]for k in ['numerics_profile','weight_engine','defer_phase','prefix_cache','cayley_cache','center_green_fast']},execution=json.loads(Path(str(pre)+'.execution.json').read_text())))
  except Exception as e:problems.append(dict(seed=r['seed'],kind='NUMERICAL_FAILURE'if j and j.get('numerical_failures',0)else'EXECUTION_OR_MISSING_DATA',error=repr(e)))
 if problems:cases.append(dict(case=name,status='BLOCKED_OR_NUMERICAL_ERROR',problems=problems));continue
 refs=[discrete(c,k)for k in [4,6,8]];cr=continuous(8,c['tau']);ref=refs[-1]['observables'];ct=cr['observables'];ground=cr['ground_observables'];budgetref=continuous(8,2)['observables']if c['stage']=='grid'else ref
 valid=[]
 for b in g['block_sizes']:
  th,se,si=jackknife(grouped(aa,b));eligible=min(math.ceil(len(a)/b)for a in aa)>=g['minimum_blocks_per_seed']
  if eligible:valid.append((b,th,se,si))
  for i,n in enumerate(NAMES):curves.append(dict(case=name,observable=n,block_sweeps=b,value=float(th[i].real),MC_SE=float(se[i]),MC_imag_SE=float(si[i]),minimum_blocks_per_seed=min(math.ceil(len(a)/b)for a in aa),eligible_for_main_error=eligible))
 theta=derive(sum(a.sum(axis=0)for a in aa));se=np.max([v[2]for v in valid],axis=0);si=np.max([v[3]for v in valid],axis=0);h=influences(aa);taus=[];rhats=[]
 for i,n in enumerate(NAMES):
  ts=[]
  for r,v in zip(jobs,h):
   rho,tau,window=acf(v[:,i]);ts.append(tau);autocorr.append(dict(case=name,seed=r['seed'],observable=n,tau_int=tau,variance_inflation=2*tau,IPS_window=window,rho=rho.tolist(),units='one recorded production sweep',estimator='common-phase nonlinear-observable linear influence, centered within seed'))
  taus.append(max(ts));rhats.append(rhat([v[:,i]for v in h]))
 phase=sum(a[:,0].sum()for a in aa)/sum(len(a)for a in aa);pse=0.
 for b,_,_,_ in valid:
  zz=grouped(aa,b)[:,:2];pw=zz[:,1].real/zz[:,1].real.sum();pv=zz[:,0]/zz[:,1];pse=max(pse,float(np.sqrt((pw*abs(pv-phase)**2).sum()/(len(pw)-1))))
 phaseok=abs(phase)>g['phase_resolution_sigma']*pse;cpu=sum(j['timing']['measurement_seconds']for j in js);samples=sum(len(a)for a in aa);tps=cpu/samples;statuses=[]
 for i,n in enumerate(NAMES):
  edok,cuterr=cutoff([r['observables']for r in refs],n);numeric=cr['ground_residual']<=g['ED_eigen_residual']and cr['propagation_segmentation_residual']<=g['ED_propagation_residual'];ratio=valid[-1][2][i]/valid[-2][2][i]if len(valid)>1 and valid[-2][2][i]>1e-12 else 1.;plateau=len(valid)>1 and g['plateau_ratio'][0]<=ratio<=g['plateau_ratio'][1]and valid[-1][0]>=g['minimum_block_in_tau_int']*taus[i];stationary=rhats[i]<=g['split_Rhat_maximum'];limit=target(n,budgetref[n],c['stage']=='projection');precise=se[i]<=limit;delta=float(theta[i].real-ref[n]);compatible=abs(delta)<=g['sigma']*se[i]+cuterr+g['absolute_floor']and abs(theta[i].imag)<=g['sigma']*si[i]+g['absolute_floor']
  if not edok or not numeric:status='INCONCLUSIVE_ED'
  elif not phaseok:status='INCONCLUSIVE_PHASE'
  elif not plateau or not stationary:status='INCONCLUSIVE_AUTOCORRELATION'
  elif not precise:status='INCONCLUSIVE_PRECISION'
  elif not compatible:status='INCOMPATIBLE'
  else:status='PASS'
  statuses.append(status);reliable=edok and numeric and phaseok and plateau and stationary;cost=(float(se[i])/limit)**2*cpu if reliable else None
  rows.append(dict(case=name,stage=c['stage'],algorithm=c['algorithm'],tau=c['tau'],dt=c['dt'],observable=n,ED_finite_dt=ref[n],ED_continuous=ct[n],ED_ground=ground[n],QMC=float(theta[i].real),QMC_imag=float(theta[i].imag),MC_SE=float(se[i]),MC_imag_SE=float(si[i]),difference_QMC_minus_discrete_ED=delta,ED_cutoff_error_estimate=cuterr,ED_solver_residual=cr['ground_residual'],ED_propagation_residual=cr['propagation_segmentation_residual'],ED_time_discretization_bias=ref[n]-ct[n],ED_finite_projection_bias=ct[n]-ground[n],fixed_precision_target=limit,mean_compatible_six_SE=bool(compatible),precision_met=bool(precise),blocking_plateau=bool(plateau),last_block_SE_ratio=float(ratio),largest_eligible_block=valid[-1][0],max_tau_int_across_seeds=taus[i],split_Rhat=rhats[i],phase_resolved=bool(phaseok),acceptance=sum(j['accepted_measurement']for j in js)/sum(j['proposals_measurement']for j in js),measurement_seconds_per_sweep=tps,measurement_CPU_seconds=cpu,warmup_CPU_seconds=sum(j['timing']['warmup_seconds']for j in js),estimated_CPU_seconds_to_common_error=cost,estimated_four_core_wall_seconds_to_common_error=cost/4 if cost else None,time_estimate_scope='stationary variance scaling estimate, not an executed stopping time'if reliable else'not certified: unresolved autocorrelation/phase/reference',status=status))
 case_status='INCOMPATIBLE'if'INCOMPATIBLE'in statuses else'PASS'if all(s=='PASS'for s in statuses)else'INCONCLUSIVE';cases.append(dict(case=name,status=case_status,samples=samples,observables_by_status={s:statuses.count(s)for s in sorted(set(statuses))},average_phase=[float(phase.real),float(phase.imag)]))
 # Save incrementally in case report processing is interrupted; no partial PASS.
 write('observables.csv',rows);write('blocking.csv',curves)
# Relative time cost within each dt and common absolute precision goal.
ratios=[]
for r in rows:
 if r['stage']!='grid':continue
 baseline=next((v for v in rows if v['stage']=='grid'and v['algorithm']=='guided'and v['dt']==r['dt']and v['observable']==r['observable']),None)
 if baseline:ratios.append(dict(case=r['case'],dt=r['dt'],observable=r['observable'],CPU_time_ratio_to_guided_same_dt=r['estimated_CPU_seconds_to_common_error']/baseline['estimated_CPU_seconds_to_common_error']if r['estimated_CPU_seconds_to_common_error']and baseline['estimated_CPU_seconds_to_common_error']else None,status='ESTIMATED'if r['estimated_CPU_seconds_to_common_error']and baseline['estimated_CPU_seconds_to_common_error']else'INCONCLUSIVE'))
time_differences=[]
for tau in [2,4,8,16]:
 for alg in ['pcn','guided','local']:
  dts=[.2,.1,.05]if tau==2 else[.1,.05]
  for large,small in zip(dts,dts[1:]):
   for n in NAMES:
    left=next((r for r in rows if r['tau']==tau and r['algorithm']==alg and r['dt']==large and r['observable']==n),None);right=next((r for r in rows if r['tau']==tau and r['algorithm']==alg and r['dt']==small and r['observable']==n),None)
    if left is None or right is None:continue
    delta=left['QMC']-right['QMC'];ed=left['ED_finite_dt']-right['ED_finite_dt'];se=math.hypot(left['MC_SE'],right['MC_SE']);reliable=left['status']=='PASS'and right['status']=='PASS';compatible=abs(delta-ed)<=g['sigma']*se+left['ED_cutoff_error_estimate']+right['ED_cutoff_error_estimate']+g['absolute_floor'];time_differences.append(dict(tau=tau,algorithm=alg,dt_large=large,dt_small=small,observable=n,QMC_difference=delta,independent_seed_MC_SE=se,ED_discretization_difference=ed,difference_of_differences=delta-ed,status='INCONCLUSIVE'if not reliable else'PASS'if compatible else'INCOMPATIBLE'))
write('dt-comparisons.csv',time_differences)
projection=[]
for tau in [2,4,8,16]:
 rr=[continuous(k,tau)for k in [4,6,8]]
 for n in NAMES:
  values=[r['observables']for r in rr];gs=[r['ground_observables']for r in rr];edok,te=cutoff(values,n);gok,ge=cutoff(gs,n);edok=edok and all(r['ground_residual']<=g['ED_eigen_residual']and r['propagation_segmentation_residual']<=g['ED_propagation_residual']for r in rr);a=values[-1][n];b=gs[-1][n];rules=cfg['projection_ground_gates']
  if n=='majorana_energy':tol=rules['matter_energy_absolute']
  elif n.startswith('B4'):tol=rules['B4_absolute']
  elif n.startswith('m4'):tol=rules['fourth_relative']*abs(b)
  elif n.startswith(('X2_','m2','S_')):tol=rules['second_relative']*abs(b)
  else:tol=rules['mean_or_bond_absolute']
  tol=max(tol,rules['near_zero_absolute']if abs(b)<1e-8 else 0);bound=abs(a-b)+te+ge;projection.append(dict(tau=tau,observable=n,ED_continuous=a,ED_ground=b,finite_projection_bias=a-b,finite_ED_cutoff_error_estimate=te,ground_ED_cutoff_error_estimate=ge,ground_convergence_tolerance=tol,status='INCONCLUSIVE_ED'if not(edok and gok)else'CONVERGED_AT_DECLARED_PRECISION'if bound<=tol else'PHYSICALLY_UNCONVERGED'))
for r in rows:
 pr=next(v for v in projection if v['tau']==r['tau']and v['observable']==r['observable']);r['ground_convergence_tolerance']=pr['ground_convergence_tolerance'];r['QMC_ground_distance_bound']=abs(r['QMC']-r['ED_ground'])+2*r['MC_SE']+abs(r['ED_time_discretization_bias'])+r['ED_cutoff_error_estimate'];r['QMC_ground_resolution']='RESOLVED_AT_DECLARED_PRECISION'if r['status']=='PASS'and r['QMC_ground_distance_bound']<=pr['ground_convergence_tolerance']else'NOT_RESOLVED';r['ED_projection_status']=pr['status']
write('observables.csv',rows);write('efficiency-ratios.csv',ratios);write('projection-convergence.csv',projection)
(OUT/'autocorrelation.json').write_text(json.dumps(autocorr)+'\n');(OUT/'chain-records.json').write_text(json.dumps(chains_report,indent=2)+'\n')
guided_forces=[r['guided_bond']['force_by_type']for r in chains_report if r['guided_bond']['enabled']and r['guided_bond']['reference']=='bare']if chains_report and 'reference'in chains_report[0]['guided_bond']else[]
if guided_forces:assert np.max(abs(np.array(guided_forces)-guided_forces[0]))<1e-12
summary=dict(pcn=pcn['status'],cases=cases,projection={str(tau):'CONVERGED_AT_DECLARED_PRECISION'if all(v['status']=='CONVERGED_AT_DECLARED_PRECISION'for v in projection if v['tau']==tau)else'NOT_CERTIFIED'for tau in [2,4,8,16]},definitions=dict(tau_int='1/2+sum rho; IPS per seed for interacting chains, variance inflation2*tau_int',matter_energy='excludes oscillator kinetic and zero-point energy',efficiency='same absolute error per observable across dt; statistics and full joint Binder phase jackknife',physical_convergence='continuous-time ED to same-parity ground, every listed observable must meet its predeclared tolerance; QMC finite-tau compatibility is separate'),original_history='PASS/INCONCLUSIVE reports unchanged; pilot data not pooled with formal',status='INCOMPATIBLE'if pcn['status']=='INCOMPATIBLE'or any(c['status']=='INCOMPATIBLE'for c in cases)else'COMPLETE_WITH_INCONCLUSIVE'if pcn['status']!='PASS'or any(c['status']!='PASS'for c in cases)else'PASS')
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');lines=['PCN / fixed-physics grid / projection validation',time.strftime('%Y-%m-%d %H:%M:%S UTC',time.gmtime()),'Pure pCN: '+pcn['status']]+[c['case']+': '+c['status']for c in cases]+['Projection ED at declared per-observable tolerances: '+str(summary['projection']),'matter_energy is matter only. All systematic biases and MC errors are separate columns.','QMC ground-distance resolution (two SE plus dt/cut bounds) is a separate per-observable column, not implied by ED state convergence.',
'No inference of thermodynamic criticality, full spin projection, universal update correctness or SUSY.'];(OUT/'REPORT.txt').write_text('\n'.join(lines)+'\n');print(json.dumps(summary,indent=2))

import sys
sys.exit(1 if summary['status']=='INCOMPATIBLE' else 0 if summary['status']=='PASS' else 2)
