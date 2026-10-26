#!/usr/bin/env python3
"""Independent split exponential numerical consistency at maximum cutoff."""
import argparse,json
import numpy as np
from scipy.sparse.linalg import expm_multiply
from ed_reference import Reference
p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
r=Reference(8,.4);A=-.1*r.HM;trace=float(A.diagonal().sum())
v=np.random.default_rng(2026101701).normal(size=r.dim);v/=np.linalg.norm(v)
y=expm_multiply(A,v,traceA=trace)
z=expm_multiply(A/2,v,traceA=trace/2);z=expm_multiply(A/2,z,traceA=trace/2)
error=float(np.linalg.norm(y-z)/np.linalg.norm(y));assert error<1e-10
with open(a.output,'w') as f:json.dump(dict(status='PASS',Ncut=8,dimension=r.dim,g=.4,matter_half_step=.1,random_seed=2026101701,relative_semigroup_residual=error,tolerance=1e-10,scope='numerical exponential consistency, not cutoff or MC error'),f,indent=2)
print('PASS',error)
