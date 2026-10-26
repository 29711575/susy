#!/usr/bin/env python3
"""Deterministic ED checks, cutoff convergence and separate finite-tau bias."""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
from ed_reference import Bosons,Reference

def self_test():
 b=Bosons(2);ext=Bosons(4);v=np.zeros(ext.size);row=np.flatnonzero((b.occ[:,0]==2)&(b.total==2))[0];v[ext.embedding(b)[row]]=1
 x=ext.X[0]@v;xx=ext.X[0]@x
 assert abs(np.vdot(x,x)-2.5)<1e-13 and abs(np.vdot(xx,xx)-9.75)<1e-13
 small=np.zeros(b.size);small[row]=1
 assert abs(np.linalg.norm(b.X[0]@small)**2-2.5)>.1
 r=Reference(0,0);psi,_=r.project(2);o,_=r.measure(psi)
 assert abs(o['X2_z']-.5)<1e-13 and abs(o['m4_q3_z']-3/64)<1e-13 and abs(o['B4_q3_z']-3)<1e-13
 return dict(status='PASS',boundary_number_state_n2_X2=2.5,boundary_number_state_n2_X4=9.75,wrong_truncated_X_squared=float(np.linalg.norm(b.X[0]@small)**2),vacuum_cut0_P_X2_P=.5,vacuum_cut0_P_D4_P=3/64)

def analyze(root,out,config):
 gates=config['gates'];rows=[];cases=[]
 for g,cuts in [(0,[2,4,6]),(.2,[4,6,8]),(.4,[4,6,8])]:
  data=[json.loads((root/f"ed-g{str(g).replace('.','p')}-cut{n}.json").read_text()) for n in cuts]
  for ti in range(2):
   rs=[d['results'][ti] for d in data];bad=[]
   for name,value in rs[-1]['observables'].items():
    coarse=abs(rs[1]['observables'][name]-rs[0]['observables'][name]);fine=abs(value-rs[1]['observables'][name]);tol=gates['ED_B4_truncation_absolute'] if name.startswith('B4') else gates['ED_truncation_abs_relative']*max(1,abs(value))
    contraction=fine<=coarse or fine<1e-10;passed=fine<=tol and contraction
    rows.append(dict(g=g,tau=rs[-1]['tau'],observable=name,ED_value=value,coarse_cutoff=cuts[1],fine_cutoff=cuts[2],coarse_change=coarse,fine_change=fine,truncation_error_estimate=2*fine,cutoff_tolerance=tol,cutoff_status='PASS' if passed else 'INCONCLUSIVE',projection_difference_from_same_parity_ground=rs[-1]['finite_tau_minus_ground'].get(name),eigensolver_residual=rs[-1]['ground_residual'],propagation_residual=rs[-1]['propagation_segmentation_residual']))
    if not passed:bad.append(name)
   h=rs[-1]['response_finite_differences'];halfdiff=abs(h['0.005']['derivative']-h['0.01']['derivative']);numok=rs[-1]['ground_residual']<gates['ED_eigen_residual'] and rs[-1]['propagation_segmentation_residual']<gates['ED_propagation_segmentation_residual'];fieldok=halfdiff<gates['ED_field_halving_abs_relative']
   # Ground quantities also need cutoff convergence, independently from the
   # finite-time vectors. Their difference is NOT added to finite-tau errors.
   gd={k:abs(rs[-1]['ground_observables'][k]-rs[-2]['ground_observables'][k]) for k in rs[-1]['ground_observables']}
   groundbad=[k for k,v in gd.items() if v>(gates['ED_B4_truncation_absolute'] if k.startswith('B4') else gates['ED_truncation_abs_relative']*max(1,abs(rs[-1]['ground_observables'][k])))]
   cases.append(dict(g=g,tau=rs[-1]['tau'],status='PASS' if not bad and numok and fieldok and not groundbad else 'INCONCLUSIVE',cutoffs=cuts,imprecise_cutoff_observables=bad,ground_imprecise_cutoff_observables=groundbad,max_ground_cutoff_change=max(gd.values()),ground_energy=rs[-1]['ground_energy'],trial_parity=data[-1]['metadata']['trial_parity'],trial_ground_overlap_squared=rs[-1]['trial_ground_overlap_squared'],boundary_probability=rs[-1]['truncation']['boundary_probability'],near_boundary_probability=rs[-1]['truncation']['near_boundary_probability'],field_halving_change=halfdiff,response_field_error_estimate=halfdiff/3,eigensolver_residual=rs[-1]['ground_residual'],propagation_residual=rs[-1]['propagation_segmentation_residual'],projection_matter_energy_difference=rs[-1]['finite_tau_minus_ground']['majorana_energy']))
 out.mkdir(parents=True,exist_ok=True)
 with (out/'ED-observables.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary=dict(status='PASS' if all(c['status']=='PASS' for c in cases) else 'INCONCLUSIVE',self_test=self_test(),cases=cases,ground_projection_status='finite tau 2/4 not converged for interacting observables; differences are reported, not used as finite-tau comparison error',scope='all 12 quantum oscillators; physical 2x2 honeycomb; fixed odd trial sector')
 (out/'ED-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

def main():
 p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--config',type=Path,default=Path(__file__).with_name('full-config.json'));a=p.parse_args();analyze(a.directory,a.output,json.loads(a.config.read_text()))
if __name__=='__main__':main()
