"""Deterministic independent Fock/derivative/geometry tests, no MC."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,sys
from pathlib import Path
import numpy as np
from physical_model import PhysicalStaticModel,parity_rule,orthogonal_descent


def kron_all(items):
    out=np.array([[1]],complex)
    for a in items:out=np.kron(out,a)
    return out


def fock_checks():
    rng=np.random.default_rng(106290400);z=np.diag([1.,-1.]);x=np.array([[0.,1.],[1.,0.]]);y=np.array([[0.,-1j],[1j,0.]])
    rows=[]
    for n in [2,3,4]:
        aa=[kron_all([z]*i+[x]+[np.eye(2)]*(n-i-1)) for i in range(n)]
        bb=[kron_all([z]*i+[-y]+[np.eye(2)]*(n-i-1)) for i in range(n)]
        parity=kron_all([z]*n).diagonal().real
        for repeat in range(12):
            t=rng.normal(size=(n,n));u,s,vh=np.linalg.svd(t);det=int(np.linalg.slogdet(t)[0]);ham=sum(1j*t[i,j]*aa[i]@bb[j] for i in range(n) for j in range(n))
            for required in [-1,1]:
                idx=np.flatnonzero(parity==required);exact=float(np.linalg.eigvalsh(ham[np.ix_(idx,idx)])[0]);formula=float(-s.sum()+(2*s[-1] if ((-1)**n)*det!=required else 0.));err=abs(exact-formula);assert err<2e-12
                rows.append(dict(n=n,required_parity=required,detT=det,error=err))
    return rows


def main():
    rows=fock_checks();rng=np.random.default_rng(106290411);other=[]
    for L in [2,3,6,9]:
        rule=parity_rule(L,L)
        if L==6:assert rule['required_code_parity']==1 and rule['target_detT_sign']==1
        if L==9:assert rule['required_code_parity']==-1 and rule['target_detT_sign']==1
        m=PhysicalStaticModel(L,2.4);x=rng.normal(size=m.nb)*.3;d=rng.normal(size=m.nb);d/=np.linalg.norm(d);step=1e-5
        e,g=m.objective(x);ep,gp=m.objective(x+step*d);em,gm=m.objective(x-step*d);ge=abs((ep-em)/(2*step)-g@d);assert ge<2e-7
        h,method,sym=m.hessian(x);he=float(abs((gp-gm)/(2*step)-h@d).max());assert he<2e-6
        o=m.physical_polar(x);reduced,_=m.reduced(o);completed=m.objective(m.g*o[m.a,m.b])[0];assert completed<=reduced+1e-10
        assert reduced>=m.lower_bound()-1e-9
        other.append(dict(L=L,parity=rule,gradient_error=float(ge),hessian_error=he,hessian_method=method))
    m=PhysicalStaticModel(3,2.4);seed,_=orthogonal_descent(m,rng,maxiter=80);assert np.isfinite(seed).all()
    report=dict(status='PASS',Fock_cases=len(rows),maximum_Fock_energy_error=max(r['error'] for r in rows),Fock=rows,derivatives=other,orthogonal_search='orthogonality, parity, monotonic energy and analytic lower bound checked',scope='Independent finite Fock matrices and deterministic derivative checks; not MC or global-minimum proof')
    out=Path(sys.argv[1]);out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['Fock','derivatives']},indent=2))
if __name__=='__main__':main()
