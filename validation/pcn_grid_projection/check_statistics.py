"""Deterministic algebra check; no reference action/sampler is called."""
import numpy as np
from moments import RAW,NAMES,derive,jackknife,influences,acf
from analyze_pcn import precision
x=np.array([.2,.3,.7,.4,1.1,.5]);phase=np.array([1,-.2+.1j,1j,.6+.8j,-.8+.6j,1]);r=np.ones((6,len(RAW)),complex)
for q in [0,3]:r[:,RAW.index(f'AmpA2_q{q}_z')]=x*x;r[:,RAW.index(f'AmpA4_q{q}_z')]=x**4
z=np.column_stack([phase,np.ones(6),phase[:,None]*r]);value=derive(z.sum(axis=0));expected=(phase@x**4/phase.sum())/(phase@x**2/phase.sum())**2;assert abs(value[NAMES.index('B4_q3_z')]-expected)<1e-12
z[:,0]=1;z[:,2:]=r;z[:,2+RAW.index('X_y')]=x;theta,se,_=jackknife(z);assert abs(se[NAMES.index('X_y')]-x.std(ddof=1)/np.sqrt(6))<1e-12
Q=precision(20,.1,100/9,.3);C=np.linalg.inv(Q);ref=.15*np.exp(-.03*abs(np.arange(21)[:,None]-np.arange(21)));assert np.max(abs(C-ref))<2e-9
for b in [.3,.6,1.]:a=np.sqrt(1-b*b);assert np.max(abs(a*a*C+b*b*C-C))<1e-14
print('PASS complex common-phase Binder, linear jackknife, Mehler Q and pCN stationary covariance')
