#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--continuous',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--config',type=Path,default=Path(__file__).with_name('full-config.json'));a=p.parse_args();gates=json.loads(a.config.read_text())['gates'];rows=[];cases=[]
for g,tau,dt in [(.4,2,.2),(.4,2,.1),(.4,2,.05),(.4,4,.1),(.2,2,.1)]:
 ds=[json.loads((a.directory/f'g{g}-tau{tau}-dt{dt}-cut{n}.json').read_text()) for n in [4,6,8]];ct=next(r for r in json.loads((a.continuous/f"ed-g{str(g).replace('.','p')}-cut8.json").read_text())['results'] if r['tau']==tau);bad=[]
 for name,v in ds[-1]['observables'].items():
  coarse=abs(ds[1]['observables'][name]-ds[0]['observables'][name]);fine=abs(v-ds[1]['observables'][name]);tol=gates['ED_B4_truncation_absolute'] if name.startswith('B4') else gates['ED_truncation_abs_relative']*max(1,abs(v));ok=fine<tol and (fine<=coarse or fine<1e-10)
  if not ok:bad.append(name)
  rows.append(dict(g=g,tau=tau,dt=dt,observable=name,ED_discrete=v,ED_continuous=ct['observables'][name],discretization_shift=v-ct['observables'][name],N4_to_N6=coarse,N6_to_N8=fine,estimated_cutoff_error=2*fine,status='PASS' if ok else 'INCONCLUSIVE'))
 cases.append(dict(g=g,tau=tau,dt=dt,status='PASS' if not bad else 'INCONCLUSIVE',unconverged_observables=bad,boundary_probability=ds[-1]['truncation']['boundary_probability'],near_boundary_probability=ds[-1]['truncation']['near_boundary_probability'],source_tangent_finite_difference_error=ds[0].get('tangent_difference')))
a.output.mkdir(parents=True,exist_ok=True)
with (a.output/'discrete-ED-observables.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(a.output/'discrete-ED-summary.json').write_text(json.dumps(dict(status='PASS' if all(c['status']=='PASS' for c in cases) else 'INCONCLUSIVE',cases=cases,scope='All twelve quantum phonons, Ncut4/6/8, same symmetric split and center cut as production. Source derivative via number-basis exponential divided differences.'),indent=2)+'\n')
print([(c['g'],c['tau'],c['dt'],c['status']) for c in cases])
