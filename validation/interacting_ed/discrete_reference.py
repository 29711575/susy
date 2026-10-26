#!/usr/bin/env python3
"""Independent number-basis Strang reference and exact source tangent.
S=exp(-dt HM/2) exp[-dt(Hph-h B)] exp(-dt HM/2), B=Nc Dcos.
The derivative of the middle exponential is evaluated by its spectral divided
 difference in the occupation basis, not by a production Mehler-score helper.
"""
import argparse,hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import expm_multiply
from ed_reference import Reference,Bosons

def projected(r,tau,dt):
 M=round(tau/dt);assert M%2==0 and abs(M*dt-tau)<1e-12
 v=np.zeros((r.dim,2));v[:,0]=r.psi0
 energies=r.w*(r.bos.total+6);d=np.exp(-dt*energies)
 B=(4*r.D[3]).tocoo();en=energies[B.row];em=energies[B.col]
 kernel=B.copy();kernel.data=B.data*(np.exp(-dt*em)-np.exp(-dt*en))/(en-em);kernel=kernel.tocsr()
 trace=float((-dt/2)*r.HM.diagonal().sum());logs=[]
 for step in range(M//2):
  v=expm_multiply((-dt/2)*r.HM,v,traceA=trace)
  f=v.reshape(r.bos.size,r.nf,2).copy();old=f[:,:,0].copy()
  f*=d[:,None,None];f[:,:,1]+=kernel@old
  v=expm_multiply((-dt/2)*r.HM,f.reshape(r.dim,2),traceA=trace)
  norm=np.linalg.norm(v[:,0]);v/=norm;logs.append(float(math.log(norm)))
 psi=v[:,0];tangent=v[:,1];Dpsi=(r.D[3]@psi.reshape(-1,r.nf)).reshape(-1)
 mean=np.vdot(psi,Dpsi).real;chi=2*(np.vdot(tangent,Dpsi).real-np.vdot(psi,tangent).real*mean)
 return psi,float(chi),dict(normalization_step_logs=logs,half_log_norm=sum(logs),full_log_norm=2*sum(logs),tangent_norm=float(np.linalg.norm(tangent)),source_kernel='B_nm*(exp(-dt E_m)-exp(-dt E_n))/(E_n-E_m)')

def field_check(r,tau,dt,h):
 from scipy import sparse
 middle=sparse.diags(r.w*(r.bos.total+6))-h*4*r.D[3]
 v=r.psi0.copy();M=round(tau/dt)
 for step in range(M//2):
  v=expm_multiply((-dt/2)*r.HM,v,traceA=float((-dt/2)*r.HM.diagonal().sum()))
  v=expm_multiply(-dt*middle,v.reshape(-1,r.nf),traceA=float(-dt*middle.diagonal().sum())).reshape(-1)
  v=expm_multiply((-dt/2)*r.HM,v,traceA=float((-dt/2)*r.HM.diagonal().sum()));v/=np.linalg.norm(v)
 return float(np.vdot(v.reshape(-1,r.nf),r.D[3]@v.reshape(-1,r.nf)).real)

def main():
 p=argparse.ArgumentParser();p.add_argument('--cut',type=int,required=True);p.add_argument('--g',type=float,required=True);p.add_argument('--tau',type=float,required=True);p.add_argument('--dt',type=float,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--field-check',action='store_true');a=p.parse_args();assert not a.output.exists();t=time.monotonic()
 r=Reference(a.cut,a.g);psi,chi,norm=projected(r,a.tau,a.dt);obs,bound=r.measure(psi);obs['chi_cos_q3']=chi
 record=dict(schema='independent-full-quantum-Strang-ed-v1',Ncut=a.cut,g=a.g,tau=a.tau,dt=a.dt,M=round(a.tau/a.dt),nb=12,metadata=r.meta,dimension=r.dim,observables=obs,truncation=bound,normalization=norm,seconds=time.monotonic()-t,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 if a.field_check:
  values={str(h):(field_check(r,a.tau,a.dt,h)-field_check(r,a.tau,a.dt,-h))/(2*h) for h in [.01,.005]};record['field_finite_differences']=values;record['tangent_difference']=abs((4*values['0.005']-values['0.01'])/3-chi);assert record['tangent_difference']<2e-8
 a.output.write_text(json.dumps(record,indent=2)+'\n');print('COMPLETE',a.output,record['seconds'],flush=True)
if __name__=='__main__':main()
