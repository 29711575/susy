"""Report ED cutoff, finite dt, finite projection and solver errors separately."""
import csv,json
from pathlib import Path
from moments import NAMES
P=Path.cwd();E=Path('/home/sunxr/quick1010-pcn-projection-ed-20261019-20022/results');O=P/'analysis/ed';O.mkdir(parents=True,exist_ok=True);cfg=json.loads((P/'config.json').read_text());g=cfg['gates'];rows=[];discrete=[];cases=[]
def cutoff(a,n):
 c=abs(a[1][n]-a[0][n]);f=abs(a[2][n]-a[1][n]);tol=g['ED_B4_cutoff_absolute']if n.startswith('B4')else g['ED_cutoff_absolute_relative']*max(1,abs(a[2][n]));return f<=tol and(f<=c or f<1e-10),2*f,c,f
cts=[json.loads((E/f'continuous-cut{k}.json').read_text())for k in [4,6,8]]
for tau in [2,4,8,16]:
 r=[next(v for v in d['results']if v['tau']==tau)for d in cts];num=all(v['ground_residual']<=g['ED_eigen_residual']and v['propagation_segmentation_residual']<=g['ED_propagation_residual']for v in r);statuses=[]
 for n in NAMES:
  vals=[v['observables']for v in r];gs=[v['ground_observables']for v in r];ok,te,c,f=cutoff(vals,n);gok,ge,_,_=cutoff(gs,n);a=vals[-1][n];b=gs[-1][n];rules=cfg['projection_ground_gates'];tol=rules['matter_energy_absolute']if n=='majorana_energy'else rules['B4_absolute']if n.startswith('B4')else rules['fourth_relative']*abs(b)if n.startswith('m4')else rules['second_relative']*abs(b)if n.startswith(('X2_','S_','m2'))else rules['mean_or_bond_absolute'];tol=max(tol,rules['near_zero_absolute']if abs(b)<1e-8 else 0);status='INCONCLUSIVE_ED'if not(ok and gok and num)else'CONVERGED_AT_DECLARED_PRECISION'if abs(a-b)+te+ge<=tol else'PHYSICALLY_UNCONVERGED';statuses.append(status);rows.append(dict(tau=tau,observable=n,ED_continuous=a,ED_ground=b,finite_projection_bias=a-b,ED_truncation_error_estimate=te,ground_truncation_error_estimate=ge,cut4_to_cut6=c,cut6_to_cut8=f,eigensolver_residual=r[-1]['ground_residual'],propagation_residual=r[-1]['propagation_segmentation_residual'],boundary_probability=r[-1]['truncation']['boundary_probability'],near_boundary_probability=r[-1]['truncation']['near_boundary_probability'],ground_convergence_tolerance=tol,status=status))
 cases.append(dict(tau=tau,status='CONVERGED_AT_DECLARED_PRECISION'if all(s=='CONVERGED_AT_DECLARED_PRECISION'for s in statuses)else'INCONCLUSIVE_ED'if'INCONCLUSIVE_ED'in statuses else'PHYSICALLY_UNCONVERGED',counts={s:statuses.count(s)for s in set(statuses)}))
for c in cfg['cases']:
 r=[]
 for k in [4,6,8]:
  name=f"g0.4-tau{c['tau']}-dt{c['dt']}-cut{k}.json";p=E/name
  if not p.exists():p=P/'references'/name
  r.append(json.loads(p.read_text()))
 ct=next(v for v in cts[-1]['results']if v['tau']==c['tau'])
 for n in NAMES:
  ok,te,coarse,fine=cutoff([v['observables']for v in r],n);discrete.append(dict(case=c['name'],tau=c['tau'],dt=c['dt'],observable=n,ED_finite_dt=r[-1]['observables'][n],ED_continuous=ct['observables'][n],time_discretization_bias=r[-1]['observables'][n]-ct['observables'][n],truncation_error_estimate=te,cut4_to_cut6=coarse,cut6_to_cut8=fine,status='PASS'if ok else'INCONCLUSIVE_ED'))
for f,data in [('projection.csv',rows),('finite-dt.csv',discrete)]:
 with (O/f).open('w')as h:w=csv.DictWriter(h,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
summary=dict(projection=cases,discrete_reference_status='PASS'if all(r['status']=='PASS'for r in discrete)else'INCONCLUSIVE',scope='all12 phonons, odd trial parity, independent continuous and identical-split ED; physical nonconvergence is not QMC error',ground_total_energy_ED_only=cts[-1]['results'][0]['ground_energy'],trial_ground_overlap_squared=cts[-1]['results'][0]['trial_ground_overlap_squared']);(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
