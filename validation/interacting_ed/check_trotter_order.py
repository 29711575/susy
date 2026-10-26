#!/usr/bin/env python3
"""Independent symmetric splitting diagnostic; not the converged C reference."""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from ed_reference import Reference
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
r=Reference(4,.4);ct,_=r.project(2);base=r.measure(ct)[0];rows=[]
for dt in [.2,.1,.05]:
 M=round(2/dt);v=r.psi0.copy();logs=[]
 for j in range(M//2):
  v=expm_multiply((-dt/2)*r.HM,v,traceA=float((-dt/2)*r.HM.diagonal().sum()))
  v*=np.repeat(np.exp(-dt*r.w*(r.bos.total+6)),r.nf)
  v=expm_multiply((-dt/2)*r.HM,v,traceA=float((-dt/2)*r.HM.diagonal().sum()))
  norm=np.linalg.norm(v);v/=norm;logs.append(float(math.log(norm)))
 o,_=r.measure(v);rows.append(dict(dt=dt,M=M,observables=o,normalization_logs=logs))
order=[]
for name in ['majorana_energy','X_z','X2_z','S_q3_zz','m4_q3_z','B4_q3_z']:
 errs=[abs(x['observables'][name]-base[name]) for x in rows]
 orders=[math.log(errs[i]/errs[i+1],2) for i in [0,1]] if min(errs)>1e-12 else None
 order.append(dict(observable=name,errors=errs,observed_orders=orders))
a.output.write_text(json.dumps(dict(status='PASS' if all(x['observed_orders'] is None or all(1.8<z<2.2 for z in x['observed_orders']) for x in order) else 'INCONCLUSIVE',Ncut=4,g=.4,tau=2,definition='S(dt)=exp(-dt HM/2) exp(-dt Hph) exp(-dt HM/2), center after M/2 complete S; finite-cut diagnostic only',continuous_observables=base,discrete=rows,orders=order),indent=2)+'\n')
