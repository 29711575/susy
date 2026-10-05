"""Exact displaced harmonic-oscillator path, no QMC weight code reused.
Single commuting Majorana bond: H=(J+gX)Q+Hph, trial Q=-1.
Trapezoidal matter action gives p=N(Cgw,C), q=N(0,C).
"""
from pathlib import Path
import numpy as np
from scipy.special import ndtr
import json
ROOT=Path(__file__).resolve().parent
out=[]
for omega in [.4,.2,.1,.05]:
 M=int(round(80/omega)); dt=.1; mass=1/omega**2; g=np.sqrt(2.4); var=omega/2; decay=np.exp(-omega*dt)
 weights=np.full(M+1,dt); weights[[0,-1]]*=.5
 # Apply C to weights with two linear recurrences, independently of core code.
 left=np.zeros(M+1);right=left.copy()
 for j in range(M+1):left[j]=weights[j]+(decay*left[j-1] if j else 0)
 for j in range(M,-1,-1):right[j]=weights[j]+(decay*right[j+1] if j<M else 0)
 mean=g*var*(left+right-weights)
 D=float(g*np.dot(weights,mean))
 rng=np.random.default_rng(38100+int(omega*1000))
 # AR draws may be batched across independent proposals but MH remains serial.
 def draws(count):
  noise=rng.normal(size=(count,M+1))*np.sqrt(var*(1-decay**2));noise[:,0]=rng.normal(size=count)*np.sqrt(var)
  for j in range(1,M+1):noise[:,j]+=decay*noise[:,j-1]
  return noise
 proposals=draws(5000)
 records={}
 for init in ['free','exact_equilibrium']:
  x=draws(1)[0]+(mean if init=='exact_equilibrium' else 0)
  logw=g*weights@x
  accepted=0; trace=[]
  for y in proposals:
   ly=g*weights@y
   if np.log(rng.random())<min(0,ly-logw):x=y.copy();logw=ly;accepted+=1
   trace.append(float(x[M//2]))
  trace=np.array(trace);records[init]={'accepted':accepted,'sweeps':5000,'mean_center':float(trace.mean()),'first_quarter_mean':float(trace[:1250].mean()),'last_quarter_mean':float(trace[-1250:].mean()),'final_center':float(trace[-1])}
 rec={'omega':omega,'mass':mass,'M':M,'tau':M*dt,'g':float(g),'exact_mean_center':float(mean[M//2]),'exact_variance_center':var,'continuous_mean_center':float(g*(1-np.exp(-omega*M*dt/2))),'whitened_mean_distance_squared':D,'KL_p_q':D/2,'exact_stationary_acceptance_independence_MH':float(2*ndtr(-np.sqrt(D/2))),'MC':records}
 out.append(rec)
 print(json.dumps(rec),flush=True)
(ROOT/'analytic_mixing.json').write_text(json.dumps(out,indent=2)+'\n')
