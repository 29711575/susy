"""One-time variance/runtime-based allocation; no formal results are consulted."""
import json,math,hashlib,subprocess,sys,time
from pathlib import Path
import numpy as np
from moments import read_series,grouped,jackknife,influences,acf,NAMES,target
P=Path.cwd();cfg=json.loads((P/'config.json').read_text());pil=json.loads((P/'pilot-plan.json').read_text());budget=cfg['formal_budget_rule'];cases=cfg['cases'];groups=[];decisions=[]
def loadref(c,cut):
 name=f"g0.4-tau{c['tau']}-dt{c['dt']}-cut{cut}.json"
 for f in [Path('/home/sunxr/quick1010-pcn-projection-ed-20261019-20022/results')/name,P/'references'/name]:
  if f.exists():return json.loads(f.read_text())
 raise FileNotFoundError(name)
for ci,c in enumerate(cases):
 rows=[r for r in pil['jobs']if r['case']['name']==c['name']];chains=[];timing=[];hashes=[]
 for r in rows:
  pre=P/'pilot/runs'/f"chain{r['index']:03d}";j=json.loads(Path(str(pre)+'_results.json').read_text());assert j['status']=='complete'and j['numerical_failures']==0 and j['seed']==r['seed'];a=read_series(str(pre)+'_validation_series.csv');assert len(a)==r['sweeps'];chains.append(a);timing.append(j['elapsed_seconds']/(r['sweeps']+r['warmup']));hashes.append(hashlib.sha256(Path(str(pre)+'_validation_series.csv').read_bytes()).hexdigest())
 refs=[loadref(c,k)['observables']for k in [4,6,8]]
 for n0 in NAMES:
  coarse=abs(refs[1][n0]-refs[0][n0]);fine=abs(refs[2][n0]-refs[1][n0]);tol=cfg['gates']['ED_B4_cutoff_absolute']if n0.startswith('B4')else cfg['gates']['ED_cutoff_absolute_relative']*max(1,abs(refs[2][n0]));assert fine<=tol and (fine<=coarse or fine<1e-10),('ED cutoff unresolved',c['name'],n0,coarse,fine)
 ref=refs[-1];ref=next(r for r in json.loads(Path('/home/sunxr/quick1010-pcn-projection-ed-20261019-20022/results/continuous-cut8.json').read_text())['results']if r['tau']==2)['observables']if c['stage']=='grid'else ref;se=np.max([jackknife(grouped(chains,f))[1]for f in [64,128,256,512,1024]],axis=0);h=influences(chains);taus=[acf(v[:,i])[1]for v in h for i in range(len(NAMES))];tau=max(taus);unresolved=max(taus)>8192/50
 ratio=max(float(se[i]/max(target(n,ref[n],c['stage']=='projection'),1e-12))for i,n in enumerate(NAMES));mult=budget['unresolved_pilot_safety_factor']if unresolved else budget['variance_safety_factor'];wanted=max(budget['minimum_sweeps_per_chain'],mult*8192*2/4*ratio**2)
 tps=max(timing);warm=min(budget['warmup_maximum'],max(budget['warmup_minimum'],int(math.ceil(budget['warmup_integrated_time_multiplier']*tau))));cap=min(budget['maximum_sweeps_per_chain'],int(budget['maximum_CPU_seconds_per_chain']/tps)-warm);cap=max(1024,2**int(math.floor(math.log2(max(1024,cap)))));n=min(cap,2**int(math.ceil(math.log2(wanted))))
 decisions.append(dict(case=c,pilot_series_sha256=hashes,pilot_SE=se.tolist(),observable_order=NAMES,pilot_max_tau_int=tau,pilot_tail_unresolved=unresolved,required_sweeps_before_cap=wanted,frozen_sweeps_per_seed=n,frozen_warmup=warm,seconds_per_sweep_bound=tps,resource_capped=n<wanted))
 jobs=[]
 for si in range(4):
  seed=202610192001+10*ci+si;args=rows[0]['args'].copy()
  for flag,v in [('--seed',seed),('--warmup',warm),('--sweeps',n)]:args[args.index(flag)+1]=str(v)
  jobs.append(dict(case=c,seed=seed,warmup=warm,sweeps=n,args=args))
 groups.append(dict(jobs=jobs,cost=(warm+n)*tps))
groups.sort(key=lambda g:-g['cost']);jobs=[];ranges=[]
for g in groups:
 start=len(jobs)
 for row in g['jobs']:row['index']=len(jobs);jobs.append(row)
 ranges.append((start,len(jobs)-1))
plan=dict(tier='formal',jobs=jobs,allocation_decisions=decisions,policy='Frozen before first formal sample; no pilot samples pooled; no conditional extension or tolerance relaxation',input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in [P/'config.json',P/'pilot-plan.json',P/'moments.py',P/'plan_formal.py',P/'report.py',P/'analyze_pcn.py',P/'run.py',P/'run_formal.py',P/'maybe_submit_report.py']},created_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
f=P/'formal-plan.json';assert not f.exists();f.write_text(json.dumps(plan,indent=2)+'\n');(P/'formal-plan.sha256').write_text(hashlib.sha256(f.read_bytes()).hexdigest()+'  formal-plan.json\n')
submissions=[]
for k,(start,end)in enumerate(ranges):
 cmd=['qsub','-N',f'pcnFormal{k:02d}','-v','PLAN=formal-plan.json','-J',f'{start}-{end}','run-formal.pbs'];r=subprocess.run(cmd,text=True,capture_output=True);submissions.append(dict(command=cmd,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr));(P/'formal-submissions.json').write_text(json.dumps(submissions,indent=2)+'\n')
 if r.returncode:sys.exit(r.returncode)
# Last terminal formal worker submits a single report without stale-array dependencies.
import runpy
runpy.run_path(str(P/'maybe_submit_report.py'),run_name='__main__')
