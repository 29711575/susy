"""Explicit physical-parity static model, with an independent SO(N) search.

No full gauge/flux-sector search is implied. The numerical base is the verified
frozen StaticModel; this adapter removes the ambiguous public 'trial' policy.
"""
import numpy as np
from static_limit_response import StaticModel


def parity_rule(lx,ly,shear=0,gauge_product=1):
    nc=lx*ly
    geom=-1 if (lx+ly+shear*(lx-shear))%2 else 1
    required=geom*(-1 if nc%2 else 1)*int(gauge_product)
    # In the fixed A/B ordering: P_ground=(-1)^Nc det(T).
    target_det=required*(-1 if nc%2 else 1)
    return dict(required_code_parity=required,target_detT_sign=target_det,
                geometry_sign=geom,gauge_product=int(gauge_product),Nc=nc)


class PhysicalStaticModel(StaticModel):
    def __init__(self,L,lam):
        super().__init__(L,L,0,[1.,1.,.5],np.sqrt(lam),1.,'trial')
        self.bare_detT_sign=self.trial_sign
        self.physical_rule=parity_rule(L,L)
        self.trial_sign=self.physical_rule['target_detT_sign']
        self.lambda_value=float(lam)

    def physical_polar(self,x):
        u,s,vh,sign,correction=self.spectral(x)
        d=np.ones(self.nc)
        if correction:d[-1]=-1
        o=(u*d)@vh
        assert abs(np.linalg.det(o)-self.trial_sign)<1e-8
        return o

    def reduced(self,o):
        values=o[self.a,self.b]
        energy=-np.dot(self.j[self.mu],values)-self.lambda_value/2*np.dot(values,values)
        grad=np.zeros_like(o)
        np.add.at(grad,(self.a,self.b),-self.j[self.mu]-self.lambda_value*values)
        return float(energy),grad

    def lower_bound(self):
        # Completing the square in X and using sum_(bond j in row i) O_ij^2<=1:
        # E >= -Nc*(sqrt(Jx^2+Jy^2+Jz^2)+lambda/2).
        return -self.nc*(float(np.linalg.norm(self.j))+self.lambda_value/2)


def orthogonal_descent(model,rng,maxiter=600,initial=None):
    n=model.nc
    if initial is None:
        o,r=np.linalg.qr(rng.normal(size=(n,n)))
        o*=np.where(np.diag(r)<0,-1.,1.)
        if int(round(np.linalg.det(o)))!=model.trial_sign:o[:,-1]*=-1
    else:o=initial.copy()
    first,_=model.reduced(o);energy=first;accepted=0
    for iteration in range(maxiter):
        energy,grad=model.reduced(o);a=grad@o.T-o@grad.T
        norm2=float(np.square(a).sum())/2
        if norm2<1e-16:break
        step=1./(2.+model.lambda_value)
        for backtrack in range(24):
            trial=np.linalg.solve(np.eye(n)+step*a/2,(np.eye(n)-step*a/2)@o)
            new_energy,_=model.reduced(trial)
            if new_energy<=energy-1e-4*step*norm2:
                o=trial;energy=new_energy;accepted+=1;break
            step*=.5
        else:break
        if (iteration+1)%100==0:
            u,s,vh=np.linalg.svd(o,full_matrices=False);o=u@vh
            if int(round(np.linalg.det(o)))!=model.trial_sign:o[:,-1]*=-1
    err=float(abs(o.T@o-np.eye(n)).max())
    assert err<1e-9 and int(round(np.linalg.det(o)))==model.trial_sign
    energy,_=model.reduced(o)
    assert energy<=first+1e-9 and energy>=model.lower_bound()-1e-8
    x=model.g*o[model.a,model.b]
    return x,dict(initial_energy=first,final_energy=energy,iterations=iteration+1,
                  accepted=accepted,orthogonality_error=err,detT_sector=model.trial_sign)
