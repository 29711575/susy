#!/usr/bin/env python3
"""Synthetic analysis tests; these are not physical Pfaffian sign witnesses."""
import json
import numpy as np
from compare_qmc import RAW,NAMES,derive,jk

def main():
    phase=np.array([1, -1, 1j, .6+.8j, -.8+.6j, 1],complex)
    d=np.array([1+2j,2-1j,-1+.5j,3+.2j,.2-2j,-.1+.4j])
    integ=np.array([.3-1j,1+2j,2-.4j,-.2+.9j,1-.7j,.4+1j])
    values=np.ones((6,len(RAW)),complex)
    def put(n,x):values[:,RAW.index(n)]=x
    for q in [0,3]:
        put(f'AmpA2_q{q}_z',abs(d)**2);put(f'AmpA4_q{q}_z',abs(d)**4)
        put(f'ChiDRe_q{q}_z',d.real)
    put('ChiIRe_q3_z',integ.real);put('ChiCosRaw_q3_zz',4*d.real*integ.real)
    put('ChiD_q3_z',d);put('ChiIconj_q3_z',integ.conj());put('ChiRaw_q3_zz',4*d*integ.conj())
    blocks=np.column_stack([phase,np.ones(6),phase[:,None]*values])
    weighted=lambda x:sum(phase*x)/sum(phase)
    out=derive(blocks.sum(axis=0))
    expected=4*(weighted(d*integ.conj())-weighted(d)*weighted(integ.conj()))
    assert abs(out[NAMES.index('chi_complex_q3')]-expected)<1e-12
    binder=weighted(abs(d)**4)/weighted(abs(d)**2)**2
    assert abs(out[NAMES.index('B4_q3_z')]-binder)<1e-12
    assert abs(expected-4*(weighted(d*integ.conj())-weighted(d)*weighted(integ).conjugate()))>.01
    # A global phase cancels from every numerator and the common denominator.
    rotated=blocks.copy();rotated[:,0]*=np.exp(.37j);rotated[:,2:]*=np.exp(.37j)
    assert np.max(abs(derive(rotated.sum(axis=0))-out))<1e-11
    # Distinct chains must not be combined into one reblocking group.
    a=blocks[:3].copy();b=blocks[3:].copy()
    result=jk([a,b],2)
    assert len(result[3])==4 and np.isfinite(result[1]).all()
    # With unit phases, the linear component jackknife equals the usual SE.
    unit=blocks.copy();unit[:,0]=1;unit[:,2:]=values
    theta,se,_,_=jk([unit[:3],unit[3:]],1)
    i=NAMES.index('Dmean_q3_z')
    assert abs(theta[i]-d.real.mean())<1e-12
    assert abs(se[i]-d.real.std(ddof=1)/np.sqrt(6))<1e-12
    print(json.dumps({'status':'PASS','scope':'synthetic postprocessing, not a physical phase witness','checks':['complex common denominator','configuration-level conjugation','nonlinear fourth moment ratio','global phase invariance','chain boundary remainder groups','linear jackknife standard error']}))
if __name__=='__main__':main()
