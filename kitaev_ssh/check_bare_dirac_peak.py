import json,csv
from pathlib import Path
import numpy as np
from convergence_analysis import read_chain,joint_estimate,full_grid_observables,estimates_json
root=Path(__file__).resolve().parents[1]
out=root/'results/phase-tau-20260917/dirac_peak_check';out.mkdir(exist_ok=True)
m=json.loads((root/'runs/phase-tau-20260917/manifest.json').read_text())
cs=[read_chain(Path(x['folder'])/'run_results.json') for x in m['commands']]
theta=np.arccos(-.25)
K=np.array([theta,-theta])/(2*np.pi);Q=(2*K)%1
B=np.array([[2*np.pi/np.sqrt(3),2*np.pi/3],[-2*np.pi/np.sqrt(3),2*np.pi/3]])
def dist(a,b):
 return min(np.linalg.norm((a-b+np.array([i,j]))@B) for i in range(-2,3) for j in range(-2,3))
grid=np.array([[i/4,j/4] for i in range(4) for j in range(4)])
ds=np.array([min(dist(x,Q),dist(x,(-Q)%1)) for x in grid])
near=np.flatnonzero(np.isclose(ds,ds.min())).tolist()
print('K fractions',K,'theta',theta,'Q fractions',Q,'nearest indices',near,'distance',ds.min())
rows=[];summary=[]
for intervals in m['plan']['intervals']:
 cc=[c for c in cs if c.metadata['intervals']==intervals]
 est=estimates_json(joint_estimate(cc,full_grid_observables))
 vals=np.array([est[f'S_q{i}_zz']['value'][0] for i in range(16)])
 peak=int(np.argmax(vals));pn=int(1+np.argmax(vals[1:]))
 tau=round(intervals*.1,12)
 for q in range(16):
  e=est[f'S_q{q}_zz']
  rows.append(dict(tau=tau,index=q,nx=q//4,ny=q%4,Szz=vals[q],SE=e['stderr_real'],distance_to_bare_valley_transfer=ds[q]))
 summary.append(dict(tau=tau,global_peak=peak,global_Szz=vals[peak],global_SE=est[f'S_q{peak}_zz']['stderr_real'],nonzero_peak=pn,nonzero_Szz=vals[pn],nonzero_SE=est[f'S_q{pn}_zz']['stderr_real'],Szz_nearest=vals[10],SE_nearest=est['S_q10_zz']['stderr_real']))
 print(summary[-1])
for name,data in [('full_grid.csv',rows),('peaks.csv',summary)]:
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
assert abs(.5+np.exp(-1j*theta)+np.exp(1j*theta))<1e-12
info=dict(reference='bare g=0, X=0; J=(1,1,0.5); k_i=k dot a_i',theta=float(theta),K_plus_fraction=K.tolist(),K_minus_fraction=(-K).tolist(),Q_mod_reciprocal=Q.tolist(),nearest_grid_indices=near,nearest_distance=float(ds.min()),conclusion='All measured global Szz peaks at q=0; not equal to bare valley transfer. Bare valley transfer is off the L4 grid. Nonzero maxima are exploratory noisy estimates, not established order.')
(out/'comparison.json').write_text(json.dumps(info,indent=2))

