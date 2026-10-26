"""Production selected-path AR(1) test; analytic Mehler Q independent of C++."""
import csv,json,math
from pathlib import Path
import numpy as np
from moments import jackknife

def precision(M,dt,m,w):
 Q=np.zeros((M+1,M+1));diag=m*w/math.tanh(w*dt);off=-m*w/math.sinh(w*dt)
 for j in range(M):Q[j,j]+=diag;Q[j+1,j+1]+=diag;Q[j,j+1]+=off;Q[j+1,j]+=off
 Q[0,0]+=m*w;Q[-1,-1]+=m*w
 return Q

def analyze(P):
 cfg=json.loads((P/'config.json').read_text())['pcn_gates'];plan=json.loads((P/'pcn-plan.json').read_text());out=P/'analysis/pcn';out.mkdir(parents=True,exist_ok=True);rows=[];covrows=[];curves=[];cases=[]
 M=20;m=100/9;w=.3;dt=.1;C=np.linalg.inv(precision(M,dt,m,w));analytic=np.exp(-w*dt*abs(np.arange(M+1)[:,None]-np.arange(M+1)))/(2*m*w);res=float(np.max(abs(C-analytic)));assert res<=cfg['absolute_reference_tolerance'];B=64;lags=cfg['acf_lags'];d=M+1
 for beta in [.3,.6,1.]:
  cc=[r for r in plan['jobs']if r['case']['beta']==beta];a=math.sqrt(1-beta*beta);window=max(16,int(math.ceil(math.log(1e-5)/math.log(a))))if a else 16;L=max(window,max(lags));sources={'path':[],'noise':[]};corrgroups=[];problems=[]
  coeff=[np.eye(d)[M//2],np.eye(d)[0],np.eye(d)[0]-np.eye(d)[-1]];labels=['center','endpoint','endpoint_difference'];vars=[float(c@C@c)for c in coeff]
  for r in cc:
   pre=P/'pcn/runs'/f"chain{r['index']:03d}"
   try:
    j=json.loads(Path(str(pre)+'_results.json').read_text());assert j['status']=='complete'and j['numerical_failures']==0;assert j['samples']==r['sweeps']and j['seed']==r['seed'];assert j['g']==0 and abs(j['mass']-m)<1e-12 and abs(j['omega_ph']-w)<1e-12 and j['intervals']==M and j['dt']==dt and j['trial_shift']==[0,0,0]and j['wilson_loops']==[1,1];assert j['path_series_bond']==2 and j['proposal']=='fixed_guided_bond_pcn';assert j['accepted_measurement']==j['proposals_measurement']==12*r['sweeps'];assert j['q_mode']['every']==0 and j['guided_bond']['force_by_type']==[0,0,0];assert abs(j['guided_bond']['beta']-beta)<1e-14
    with open(str(pre)+'_validation_series.csv')as f:v=list(csv.DictReader(f))
    z=np.array([complex(float(t['phase_real']),float(t['phase_imag']))for t in v]);assert np.max(abs(z-1))<2e-9;assert all(int(t['sweep'])==i+1 for i,t in enumerate(v));x=np.array([[float(t[f'PathX_j{k}'])for k in range(d)]for t in v]);assert len(x)==r['sweeps'];noise=(x[1:]-a*x[:-1])/beta
    for key,data,zz in [('path',x,z),('noise',noise,z[1:])]:
     seedgroups=[]
     for start in range(0,len(data),B):
      u=data[start:start+B];ph=zz[start:start+B];sums=(ph[:,None]*u).sum(axis=0);products=np.einsum('n,ni,nj->ij',ph,u,u).ravel();four=np.array([ph@u[:,k]**4 for k in [0,d-1]]);seedgroups.append(np.r_[ph.sum(),len(u),sums,products,four])
     sources[key].append(seedgroups)
    fns=[]
    for c,var in zip(coeff,vars):
     f=x@c;fns.extend([f,f*f-var])
    features=[]
    for f in fns:
     v=f[:len(f)-L];features.append(v*v)
     features.extend(v*f[k:k+len(v)]for k in lags);total=np.zeros(len(v))
     for k in range(1,window+1):total+=v*f[k:k+len(v)]
     features.append(total)
    features=np.array(features).T;ph=z[:len(features)];cg=[]
    for start in range(0,len(features),B):
     u=features[start:start+B];zz=ph[start:start+B];cg.append(np.r_[zz.sum(),len(u),(zz[:,None]*u).sum(axis=0)])
    corrgroups.append(cg)
   except Exception as e:problems.append(dict(seed=r['seed'],error=repr(e)))
  if problems:cases.append(dict(beta=beta,status='BLOCKED',problems=problems));continue
  failed=imprecise=0
  def distribution(t):
   q=t[2:]/t[0];mean=q[:d];cov=q[d:d+d*d].reshape(d,d)-mean[:,None]*mean[None,:];v=np.r_[mean,cov.ravel(),q[-2]/q[d+0]/q[d+0],q[-1]/q[d+(d-1)*d+(d-1)]**2];return v
  def reblock(chains,factor):return np.asarray([np.sum(c[k:k+factor],axis=0)for c in chains for k in range(0,len(c),factor)])
  def evaluate(chains,fun,source):
   levels=[]
   for block in [64,128,256,512,1024,2048,4096]:
    th,se,si=jackknife(reblock(chains,block//64),fun);eligible=all(len(c)*64/block>=20 for c in chains)
    for k in range(len(th)):
     if source in ['path','noise']:
      label=f'mean_j{k}'if k<d else f'cov_j{(k-d)//d}_j{(k-d)%d}'if k<d+d*d else'left_endpoint_kurtosis'if k==d+d*d else'right_endpoint_kurtosis'
     else:
      fn=k//(len(lags)+1);local=k%(len(lags)+1);label=labels[fn//2]+('_F'if fn%2==0 else'_F2_minus_variance')+('_tau_int'if local==len(lags)else f'_rho{lags[local]}')
     curves.append(dict(beta=beta,source=source,index=k,observable=label,block_sweeps=block,value=float(th[k].real),MC_SE=float(se[k]),eligible_for_main_error=eligible))
    if eligible:levels.append((th,se,si))
   return levels[0][0],np.max([r[1]for r in levels],axis=0),np.max([r[2]for r in levels],axis=0)
  for source,groups in sources.items():
   th,se,si=evaluate(groups,distribution,source)
   for k in range(d):
    delta=float(th[k].real);ok=abs(delta)<=cfg['sigma']*se[k]+2e-9;failed+=not ok;rows.append(dict(beta=beta,source=source,observable=f'mean_j{k}',reference=0,value=delta,MC_SE=float(se[k]),status='PASS'if ok else'INCOMPATIBLE'))
   for j in range(d):
    for k in range(d):
     ix=d+j*d+k;value=float(th[ix].real);err=float(se[ix]);ok=abs(value-C[j,k])<=cfg['sigma']*err+2e-9;precise=j!=k or err/C[j,k]<=cfg['covariance_diagonal_relative_SE'];failed+=not ok;imprecise+=not precise;covrows.append(dict(beta=beta,source=source,j=j,k=k,reference=C[j,k],value=value,MC_SE=err,status='INCONCLUSIVE_PRECISION'if not precise else'PASS'if ok else'INCOMPATIBLE'))
   for q,label in [(-2,'left_endpoint_kurtosis'),(-1,'right_endpoint_kurtosis')]:
    value=float(th[q].real);err=float(se[q]);ok=abs(value-3)<=cfg['sigma']*err+2e-9;precise=err<=cfg['endpoint_kurtosis_SE'];failed+=not ok;imprecise+=not precise;rows.append(dict(beta=beta,source=source,observable=label,reference=3,value=value,MC_SE=err,status='INCONCLUSIVE_PRECISION'if not precise else'PASS'if ok else'INCOMPATIBLE'))
  stride=len(lags)+2
  def correlation(t):
   q=t[2:]/t[0];out=[]
   for i in range(6):
    f=q[i*stride:(i+1)*stride];out.extend(f[1:-1]/f[0]);out.append(.5+f[-1]/f[0])
   return np.array(out)
  th,se,_=evaluate(corrgroups,correlation,'ACF_IAT');strideout=len(lags)+1
  for i in range(6):
   power=1 if i%2==0 else 2;b=a**power
   for k,lag in enumerate(lags):
    ix=i*strideout+k;ref=b**lag;err=float(se[ix]);ok=abs(th[ix].real-ref)<=cfg['sigma']*err+2e-9;precise=err<=cfg['acf_absolute_SE_maximum'];failed+=not ok;imprecise+=not precise;rows.append(dict(beta=beta,source=labels[i//2]+('_F'if power==1 else'_F2_minus_variance'),observable=f'rho_lag{lag}',reference=ref,value=float(th[ix].real),MC_SE=err,status='INCONCLUSIVE_PRECISION'if not precise else'PASS'if ok else'INCOMPATIBLE'))
   ix=i*strideout+len(lags);ref=.5+sum(b**k for k in range(1,window+1));err=float(se[ix]);ok=abs(th[ix].real-ref)<=cfg['sigma']*err+2e-9;precise=err/ref<=cfg['IAT_relative_SE_maximum'];failed+=not ok;imprecise+=not precise;rows.append(dict(beta=beta,source=labels[i//2]+('_F'if power==1 else'_F2_minus_variance'),observable=f'tau_int_window{window}',reference=ref,value=float(th[ix].real),MC_SE=err,status='INCONCLUSIVE_PRECISION'if not precise else'PASS'if ok else'INCOMPATIBLE'))
  cases.append(dict(beta=beta,CLI_parameter='--bond-beta',a=a,status='INCOMPATIBLE'if failed else'INCONCLUSIVE_PRECISION'if imprecise else'PASS',failed_metrics=failed,imprecise_metrics=imprecise,noise_covariance_reference='independent Mehler Q inverse',tau_int_F_infinite=(1+a)/(2*(1-a)),tau_int_F2_infinite=(1+a*a)/(2*(1-a*a)),finite_IAT_window=window,MC_measurement_interval_sweeps=1,selected_path_updates_per_sweep=1,bond_paths_updated_independently_per_sweep=12,proposals_and_trials='g0 zero-force pCN only, q/local/refresh disabled'))
 for file,data in [('observables.csv',rows),('full-time-covariance.csv',covrows),('blocking.csv',curves)]:
  if data:
   with (out/file).open('w')as f:c=csv.DictWriter(f,fieldnames=list(data[0]));c.writeheader();c.writerows(data)
 result=dict(status='PASS'if all(r['status']=='PASS'for r in cases)else'INCOMPATIBLE'if any(r['status']=='INCOMPATIBLE'for r in cases)else'INCONCLUSIVE',Mehler_inverse_maximum_error=res,cases=cases,IAT_definition='1/2+sum rho, variance inflation2*tau_int. No thinning; one selected path pCN per recorded sweep.',scope='one selected z bond path; all12 paths receive same independent pCN once per sweep atg0; other update types disabled')
 (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result
