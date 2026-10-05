"""Independent full matter+phonon reference for the 2x1 torus.
JW Fock operators and sparse oscillator matrices; no QMC code is imported.
In the trial even-parity sector the six bond operators span two Pauli
components. Orthogonal rotation of six identical oscillators leaves two
coupled modes and four exactly decoupled vacua. All six bonds retained in X.
"""
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,kron,eye,diags
from scipy.sparse.linalg import expm_multiply,eigsh
import json
ROOT=Path(__file__).resolve().parent
ann=[]
for bit in range(2):
 a=np.zeros((4,4),complex)
 for state in range(4):
  if state&(1<<bit):a[state^(1<<bit),state]=(-1)**((state&((1<<bit)-1)).bit_count())
 ann.append(a)
c=[a+a.conj().T for a in ann]+[1j*(a-a.conj().T) for a in ann]
Q=[]
for cell in range(2):
 for bcell in [1-cell,cell,cell]:Q.append(1j*c[cell]@c[2+bcell])
Q=np.array(Q); J=np.array([1.,1.,.5,1.,1.,.5]);H0=np.einsum('b,bij->ij',J,Q)
e,v=np.linalg.eigh(H0);ground=v[:,0]
parity=np.diag([1.,-1.,-1.,1.]);P0=float((ground.conj()@parity@ground).real)
indices=np.where(np.diag(parity)==round(P0))[0];Q=Q[:,indices,:][:,:,indices];H0=H0[indices,:][:,indices];ground=ground[indices]
pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
coeff=np.array([[np.trace(q@p).real/2 for p in pauli]for q in Q])
O,s,vh=np.linalg.svd(coeff,full_matrices=False);rank=int(np.sum(s>1e-10));assert rank==2
O=O[:,:rank];B=np.einsum('bn,bij->nij',O,Q)
assert np.max(abs(Q-np.einsum('bn,nij->bij',O,B)))<1e-12
results=[]
for omega in [.4,.1]:
 for dt in [.1,.05]:
  for cutoff in [8,12,16]:
   g=.3; tau=4.; steps=int(round(tau/(2*dt)));n=cutoff
   a=diags(np.sqrt(np.arange(1,n)),1,shape=(n,n),format='csr');x=np.sqrt(omega/2)*(a+a.T)
   Id=eye(n,format='csr');activeX=[kron(x,Id,format='csr'),kron(Id,x,format='csr')]
   num=np.add.outer(np.arange(n),np.arange(n)).ravel();hpdiag=omega*(num+1)
   HM=kron(csr_matrix(H0),eye(n*n,format='csr'),format='csr')
   for k in range(2):HM+=g*kron(csr_matrix(B[k]),activeX[k],format='csr')
   HP=diags(np.tile(hpdiag,2),format='csr');H=HM+HP
   psi=np.zeros(2*n*n,complex);psi[0]=ground[0];psi[n*n]=ground[1]
   for step in range(steps):
    psi=expm_multiply(-dt/2*HM,psi,traceA=-dt/2*HM.diagonal().sum())
    psi*=np.tile(np.exp(-dt*hpdiag),2)
    psi=expm_multiply(-dt/2*HM,psi,traceA=-dt/2*HM.diagonal().sum());psi/=np.linalg.norm(psi)
   Em=float(np.vdot(psi,HM@psi).real)
   AX=[kron(eye(2,format='csr'),xx,format='csr')for xx in activeX]
   mean=np.array([np.vdot(psi,xx@psi).real for xx in AX])
   second=np.array([[np.vdot(xx@psi,yy@psi).real for yy in AX]for xx in AX])
   mean_b=O@mean
   x2_b=np.einsum('bi,ij,bj->b',O,second,O)+(omega/2)*(1-np.sum(O**2,axis=1))
   resh=psi.reshape(2,n,n);tail=float(np.sum(abs(resh[:,-2:,:])**2)+np.sum(abs(resh[:,:,-2:])**2))
   rec={'omega':omega,'mass':1/omega**2,'dt':dt,'tau':tau,'g':g,'lambda':g*g,'cutoff':cutoff,'active_modes':2,'decoupled_modes':4,'trial_parity':P0,'Em':Em,'mean_X':mean_b.reshape(2,3).mean(axis=0).tolist(),'mean_X2':x2_b.reshape(2,3).mean(axis=0).tolist(),'cutoff_last_2_levels_population_upper_bound':tail}
   results.append(rec);print(json.dumps(rec),flush=True)
(ROOT/'references.json').write_text(json.dumps({'Lx':2,'Ly':1,'J':J[:3].tolist(),'basis_free_matter_eigenvalues':e.tolist(),'operator_span_singular_values':s.tolist(),'oscillator_rotation':O.tolist(),'protocol':'same finite total contour and symmetric split as production, not ground-state comparison','results':results},indent=2)+'\n')
