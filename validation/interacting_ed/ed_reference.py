#!/usr/bin/env python3
"""Independent full 2x2 quantum phonon-Majorana ED, all twelve oscillators.
No imports from production model, generator, Pfaffian or Green code.
Majoranas are built in the full 16-dimensional Fock space BEFORE even
bilinears are restricted to the independently identified trial parity.
"""
import argparse,hashlib,itertools,json,math,os,platform,time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh,expm_multiply

NB=12;NC=4;J=(1.,1.,.5)
def fermions():
    c=np.zeros((8,16,16),complex)
    for j in range(4):
        for n in range(16):
            sign=(-1)**((n&((1<<j)-1)).bit_count());target=n^(1<<j)
            c[j,target,n]=sign;c[4+j,target,n]=1j*(sign if n&(1<<j) else -sign)
    cliff=max(np.max(abs(c[a]@c[b]+c[b]@c[a]-2*(a==b)*np.eye(16))) for a in range(8) for b in range(8))
    assert cliff==0
    edges=[]
    for x in range(2):
        for y in range(2):
            r=2*x+y
            for mu,(bx,by) in enumerate([((x-1)%2,y),(x,(y-1)%2),(x,y)]):
                edges.append(dict(index=3*r+mu,a=r,b=4+2*bx+by,mu=mu,u=1,cell=[x,y],J=J[mu]))
    K=np.array([1j*c[e['a']]@c[e['b']] for e in edges]);assert np.max(abs(K.imag))==0;K=K.real
    assert np.max(abs(K-K.transpose(0,2,1)))==0
    for t in [-1.3,.7]:assert np.max(abs(np.linalg.eigvalsh(t*K[0])-np.array([-abs(t)]*8+[abs(t)]*8)))<2e-14
    hf=np.einsum('b,bij->ij',[e['J'] for e in edges],K);vals,vec=np.linalg.eigh(hf);trial=vec[:,0]
    parity=np.array([(-1)**n.bit_count() for n in range(16)])
    trial_parity=int(round(np.dot(abs(trial)**2,parity)));assert abs(np.dot(abs(trial)**2,parity)-trial_parity)<1e-13
    selected=np.flatnonzero(parity==trial_parity)
    assert np.linalg.norm(trial[parity!=trial_parity])<1e-13
    return K[:,selected][:,:,selected],hf[np.ix_(selected,selected)],trial[selected],dict(edges=edges,A_order='A(0,0),A(0,1),A(1,0),A(1,1); then B in same order',boundary='periodic shear=0',gauge='all u_AB=+1; zero flux Wilson(+1,+1)',full_Clifford_max_error=float(cliff),trial_parity=trial_parity,parity_basis=selected.tolist(),free_energy=float(vals[0]),free_gap=float(vals[1]-vals[0]),single_bond_check='PASS',even_operators_projected_after_full_space_construction=True)

def compositions(total,width):
    if width==1:yield(total,);return
    for n in range(total+1):
        for tail in compositions(total-n,width-1):yield (n,)+tail

class Bosons:
    def __init__(self,cut,m=1.,w=1.):
        self.cut=cut;self.m=m;self.w=w
        occ=np.array([row for n in range(cut+1) for row in compositions(n,NB)],dtype=np.uint8)
        base=cut+1 if cut else 2;power=np.array([base**b for b in range(NB)],dtype=np.int64)
        codes=occ.astype(np.int64)@power;order=np.argsort(codes);self.occ=occ[order];self.codes=codes[order];self.power=power
        self.total=self.occ.sum(axis=1);self.size=len(self.occ);assert self.size==math.comb(NB+cut,cut)
        self.X=[]
        for b in range(NB):
            src=np.flatnonzero(self.total<cut);dst=np.searchsorted(self.codes,self.codes[src]+power[b])
            val=np.sqrt(self.occ[src,b].astype(float)+1)/math.sqrt(2*m*w)
            self.X.append(sparse.csr_matrix((np.r_[val,val],(np.r_[dst,src],np.r_[src,dst])),shape=(self.size,self.size)))
    def embedding(self,smaller):
        return np.searchsorted(self.codes,smaller.occ.astype(np.int64)@self.power)

class Reference:
    def __init__(self,cut,g,m=1.,w=1.):
        t=time.monotonic();self.K,self.hf,self.trial,self.meta=fermions();self.bos=Bosons(cut,m,w);self.g=g;self.m=m;self.w=w
        self.nf=len(self.trial);self.dim=self.nf*self.bos.size
        ident=sparse.eye(self.nf,format='csr');self.HM=sparse.kron(sparse.eye(self.bos.size,format='csr'),sparse.csr_matrix(self.hf),format='csr')
        for b in range(NB):
            if g:self.HM=self.HM+g*sparse.kron(self.bos.X[b],sparse.csr_matrix(self.K[b]),format='csr')
        self.H=self.HM+sparse.kron(sparse.diags(w*(self.bos.total+NB/2)),ident,format='csr');self.H.sum_duplicates();self.H.sort_indices()
        herm=self.H-self.H.T;assert herm.nnz==0 or np.max(abs(herm.data))<1e-14
        self.psi0=np.zeros(self.dim);self.psi0[:self.nf]=self.trial
        self.shift=self.meta['free_energy']+NB*w/2
        self.build_seconds=time.monotonic()-t
        self.D={q:sum(((-1)**(n1*e['cell'][0]+n2*e['cell'][1]))*self.bos.X[e['index']]/NC for e in self.meta['edges'] if e['mu']==2) for q,(n1,n2) in [(0,(0,0)),(3,(1,1))]}
    def shifted(self,h=0.,q=3):
        hh=self.H-self.shift*sparse.eye(self.dim,format='csr')
        if h:hh=hh-h*NC*sparse.kron(self.D[q],sparse.eye(self.nf,format='csr'),format='csr')
        return hh
    def project(self,tau,h=0.,q=3,step=.5):
        A=self.shifted(h,q);psi=self.psi0.copy();lognorm=0.;records=[]
        n=max(1,math.ceil((tau/2)/step));dt=(tau/2)/n
        # Trace supplied explicitly. Each step is independently normalized;
        # log norm is accumulated rather than hiding a normalization change.
        for k in range(n):
            v=expm_multiply((-dt)*A,psi,traceA=float((-dt)*A.diagonal().sum()));norm=np.linalg.norm(v)
            assert np.isfinite(norm) and norm>0;lognorm+=math.log(norm);psi=v/norm;records.append(float(math.log(norm)))
        return psi,dict(energy_shift=self.shift,step=dt,steps=n,half_state_log_norm=lognorm,full_physical_log_norm=2*lognorm-tau*self.shift,normalization_step_logs=records)
    def measure(self,psi,extended=None):
        psi=psi/np.linalg.norm(psi);v=psi.reshape(-1,self.nf);o={}
        # X2 and D4 are P O P, evaluated by exact ladder insertions through
        # Ncut+1 / Ncut+2, NOT by powers of the truncated X matrix.
        ext=extended or Bosons(self.bos.cut+2,self.m,self.w);ve=np.zeros((ext.size,self.nf));ve[ext.embedding(self.bos)]=v
        k=[];xk=[];xm=[];x2=[]
        for b in range(NB):
            xv=ext.X[b]@ve
            xm.append(float(np.vdot(ve,xv).real));x2.append(float(np.vdot(xv,xv).real))
            kv=v@self.K[b].T;k.append(float(np.vdot(v,kv).real))
            xk.append(float(np.vdot(v,self.bos.X[b]@kv).real))
            o['BondK_b'+str(b)]=k[-1];o['BondXK_b'+str(b)]=xk[-1]
        for mu,axis in enumerate('xyz'):
            for stem,vals in [('X_',xm),('X2_',x2),('BondK_',k),('BondXK_',xk)]:o[stem+axis]=float(np.mean(vals[mu::3]))
        o['majorana_energy']=float(np.vdot(psi,self.HM@psi).real);o['E_total']=float(np.vdot(psi,self.H@psi).real)
        for q,(n1,n2) in [(0,(0,0)),(3,(1,1))]:
            De=sum(((-1)**(n1*e['cell'][0]+n2*e['cell'][1]))*ext.X[e['index']]/NC for e in self.meta['edges'] if e['mu']==2)
            dv=De@ve;ddv=De@dv;m2=float(np.vdot(dv,dv).real);m4=float(np.vdot(ddv,ddv).real)
            o.update({f'Dmean_q{q}_z':float(np.vdot(ve,dv).real),f'S_q{q}_zz':NC*m2,f'm2_q{q}_z':m2,f'm4_q{q}_z':m4,f'B4_q{q}_z':m4/(m2*m2)})
        # One explicit connected spatial phonon correlation, in addition to S.
        xv=ext.X[2]@ve;yv=ext.X[5]@ve
        o['Xz0_Xz1']=float(np.vdot(xv,yv).real)
        prob=np.sum(abs(v)**2,axis=1)
        return o,dict(boundary_probability=float(prob[self.bos.total==self.bos.cut].sum()),near_boundary_probability=float(prob[self.bos.total>=max(0,self.bos.cut-1)].sum()),total_phonon_number=float(prob@self.bos.total),operator_intermediate_cut=self.bos.cut+2,number_shell_probability=[float(prob[self.bos.total==n].sum()) for n in range(self.bos.cut+1)])
    def solve(self,taus,fields=True,ground=True):
        ext=Bosons(self.bos.cut+2,self.m,self.w);out=[]
        e0=gs=None
        if ground:
            start=np.random.default_rng(20261017+self.bos.cut).normal(size=self.dim);start/=np.linalg.norm(start)
            vals,vec=eigsh(self.H,k=1,which='SA',v0=start,tol=1e-11,maxiter=1200,ncv=24);e0=float(vals[0]);gs=vec[:,0];res=float(np.linalg.norm(self.H@gs-e0*gs));gobs,gbound=self.measure(gs,ext)
        for tau in taus:
            t=time.monotonic();psi,norm=self.project(tau);obs,bound=self.measure(psi,ext)
            # Independent time segmentation is a numerical propagation check,
            # not an estimate of statistical or physical projection error.
            check,norm2=self.project(tau,step=.25);overlap=abs(np.vdot(check,psi));prop_res=float(np.linalg.norm(check-np.sign(np.vdot(check,psi))*psi))
            resp={}
            if fields:
                for h in [.01,.005]:
                    values=[]
                    for sign in [1,-1]:
                        ph,_=self.project(tau,h=sign*h);d=float(np.vdot(ph.reshape(-1,self.nf),self.D[3]@ph.reshape(-1,self.nf)).real);values.append(d)
                    resp[str(h)]=dict(plus=values[0],minus=values[1],derivative=(values[0]-values[1])/(2*h))
                obs['chi_cos_q3']=(4*resp['0.005']['derivative']-resp['0.01']['derivative'])/3
            out.append(dict(tau=tau,observables=obs,truncation=bound,normalization=norm,propagation_segmentation_residual=prop_res,response_finite_differences=resp,seconds=time.monotonic()-t,ground_energy=e0,ground_residual=res if ground else None,ground_start_seed=20261017+self.bos.cut,trial_ground_overlap_squared=float(abs(np.vdot(self.psi0,gs))**2) if ground else None,ground_observables=gobs if ground else None,ground_truncation=gbound if ground else None,finite_tau_minus_ground={k:obs[k]-gobs[k] for k in obs if k in gobs} if ground else None))
        return out

def estimate(cuts):
    out=[]
    for n in cuts:
        p=math.comb(NB+n,n);dim=8*p
        # Conservative approximate CSR capacity (diagonal+12 fermion terms+
        # 24 ladder links), Krylov24 and at most 6 working vectors.
        out.append(dict(Ncut=n,boson_states=p,full_dimension=16*p,trial_parity_dimension=dim,dense_H_GiB=dim*dim*8/2**30,conservative_sparse_and_work_GiB=(dim*37*12+dim*(24+6)*8)/2**30))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--cut',type=int,default=4);p.add_argument('--g',type=float,default=.4);p.add_argument('--taus',type=float,nargs='+',default=[2.,4.]);p.add_argument('--output',type=Path,required=True);p.add_argument('--no-fields',action='store_true');p.add_argument('--estimate',action='store_true');a=p.parse_args()
    if a.estimate:
        a.output.write_text(json.dumps({'estimates':estimate([2,4,6,8,10]),'per_mode_nmax1_dimension':16*2**12,'per_mode_nmax2_dimension':16*3**12,'all_12_oscillators_active':True},indent=2)+'\n');return
    if a.output.exists():raise FileExistsError(a.output)
    t=time.monotonic();ref=Reference(a.cut,a.g);data=dict(schema='independent-full-quantum-phonon-ed-v1',g=a.g,mass=1,omega=1,k=1,Ncut=a.cut,dimension=ref.dim,nb=12,metadata=ref.meta,build_seconds=ref.build_seconds,nnz=ref.H.nnz,csr_bytes=sum(x.nbytes for x in [ref.H.data,ref.H.indices,ref.H.indptr]),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=platform.python_version(),numpy=np.__version__,threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']})
    data['results']=ref.solve(a.taus,fields=not a.no_fields);data['seconds']=time.monotonic()-t
    a.output.write_text(json.dumps(data,indent=2)+'\n');print('COMPLETE',a.output,'dim',ref.dim,'seconds',data['seconds'],flush=True)
if __name__=='__main__':main()
