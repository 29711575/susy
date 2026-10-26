"""Common-phase sufficient statistics; no production Hamiltonian imports."""
import csv,math
import numpy as np
RAW=['X_x','X_y','X_z','X2_x','X2_y','X2_z','majorana_energy','BondK_b2','BondXK_b2','BondK_x','BondK_y','BondK_z','BondXK_x','BondXK_y','BondXK_z','S_q0_zz','S_q3_zz','AmpA2_q0_z','AmpA4_q0_z','AmpA2_q3_z','AmpA4_q3_z']
NAMES=RAW[:15]+[s for q in [0,3]for s in [f'S_q{q}_zz',f'm2_q{q}_z',f'm4_q{q}_z',f'B4_q{q}_z']]
def derive(t):
 r={name:t[k+2]/t[0]for k,name in enumerate(RAW)};v=[r[n]for n in RAW[:15]]
 for q in [0,3]:
  m2=r[f'AmpA2_q{q}_z'];m4=r[f'AmpA4_q{q}_z'];v.extend([r[f'S_q{q}_zz'],m2,m4,m4/m2**2])
 return np.asarray(v)
def read_series(file):
 with open(file)as f:
  r=list(csv.DictReader(f))
 z=np.array([complex(float(v['phase_real']),float(v['phase_imag']))for v in r]);x=np.array([[complex(float(v[n+'_real']),float(v[n+'_imag']))for n in RAW]for v in r]);assert all(int(v['sweep'])==i+1 for i,v in enumerate(r));assert np.isfinite(x).all()
 return np.column_stack([z,np.ones(len(z)),z[:,None]*x])
def grouped(chains,f):return np.array([a[k:k+f].sum(axis=0)for a in chains for k in range(0,len(a),f)])
def jackknife(groups,fun=derive,count_column=1):
 t=groups.sum(axis=0);theta=fun(t);w=groups[:,count_column].real/t[count_column].real;deleted=np.array([fun(t-g)for g in groups]);pseudo=(theta-(1-w[:,None])*deleted)/w[:,None];mu=(w[:,None]*pseudo).sum(axis=0);se=np.sqrt((w[:,None]*(pseudo.real-mu.real)**2).sum(axis=0)/(len(w)-1));si=np.sqrt((w[:,None]*(pseudo.imag-mu.imag)**2).sum(axis=0)/(len(w)-1));return theta,se,si

def influences(chains):
 t=sum(a.sum(axis=0)for a in chains);r=t[2:]/t[0];p=t[0]/t[1];out=[]
 for a in chains:
  h=(a[:,2:]-a[:,0,None]*r)/p;v=[h[:,k]for k in range(15)]
  for q in [0,3]:
   ix=[RAW.index(n)for n in [f'S_q{q}_zz',f'AmpA2_q{q}_z',f'AmpA4_q{q}_z']];s,m2,m4=ix;v.extend([h[:,s],h[:,m2],h[:,m4],h[:,m4]/r[m2]**2-2*r[m4]*h[:,m2]/r[m2]**3])
  out.append(np.array(v).T.real)
 return out

def acf(y,maxlag=None):
 z=y-y.mean();n=len(z);var=float(z@z/n)
 if var<1e-25:return np.r_[1.,np.zeros(min(n//4,4096))],.5,0
 size=1<<(2*n-1).bit_length();f=np.fft.rfft(z,n=size);c=np.fft.irfft(f*f.conj(),n=size)[:n]/n;rho=c[:min(n//4,4096)+1]/c[0];tau=.5;window=0
 for k in range(1,len(rho)-1,2):
  pair=rho[k]+rho[k+1]
  if pair<=0:break
  tau+=pair;window=k+1
 return rho,float(tau),window

def rhat(values):
 n=min(len(v)//2 for v in values);v=np.array([p for a in values for p in [a[:n],a[-n:]]]);w=np.var(v,axis=1,ddof=1).mean();b=n*np.var(v.mean(axis=1),ddof=1)
 return float(np.sqrt(((n-1)*w/n+b/n)/w))if w>1e-25 else 1.
def target(name,ref,projection=False):
 if name.startswith('B4'):return .03 if projection else .05
 if name.startswith('m4'):return .03*abs(ref) if projection else .04*abs(ref)
 if name.startswith(('m2','S_','X2_')):return (.01 if projection else .02)*abs(ref)
 if name=='majorana_energy':return .005 if projection else .01
 return .002 if projection else .005
