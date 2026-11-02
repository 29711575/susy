#!/usr/bin/env python3
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import sys,json,csv,time,hashlib,traceback,argparse
from pathlib import Path
import numpy as np
from physical_model import PhysicalStaticModel,orthogonal_descent
from static_limit_response import response_grid,validation,write_csv
ROOT=Path(os.environ.get('STATIC_TASK_ROOT',str(Path(__file__).resolve().parents[1])))
PLAN=json.loads((ROOT/'plan.json').read_text());POINTS=PLAN['points']

def save(path,obj):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_name(path.name+'.tmp.%d'%os.getpid());temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(str(temp),str(path))

def model(p):
 m=PhysicalStaticModel(p['L'],p['lambda_value']);assert m.physical_rule['required_code_parity']==p['required_code_parity'] and m.trial_sign==p['target_detT_sign'];return m

def signature(m,x):
 f=x.reshape(m.nc,3);d=np.fft.fft2(f.reshape(m.lx,m.ly,3),axes=(0,1))/m.nc;d[0,0]=0
 return tuple(np.round(np.r_[f.mean(0),f.std(0),np.max(abs(d),axis=(0,1))],5))

def select_candidates(m,records,xs,limit=4):
 selected=[];seen=set()
 for i in sorted(range(len(records)),key=lambda i:records[i].get('energy',1e100)):
  if records[i].get('energy') is None or not np.isfinite(xs[i]).all():continue
  s=signature(m,xs[i])
  if s in seen:continue
  seen.add(s);selected.append(i)
  if len(selected)==limit:break
 assert selected;return np.array([xs[i] for i in selected]),[records[i] for i in selected]

def seed_starts(m,p,batch):
 rng=np.random.default_rng(106290000+p['index']*10000+batch*1000);n=PLAN['random_starts_per_batch'];base=m.relax(np.zeros(m.nb)).x;starts=[]
 for i in range(n):
  amp=(.12,.3,.6,1.,1.8,3.,.9,1.4)[i%8]*m.g;center=base if i%3 else np.zeros(m.nb)
  starts.append(('random_%d'%i,center+rng.normal(0.,amp,m.nb)))
 structured=[('zero',np.zeros(m.nb)),('uniform',np.ones(m.nb)*m.g/np.sqrt(3)),('uniform_stationary',base)]
 for mu in range(3):
  x=np.zeros((m.nc,3));x[:,mu]=m.g;structured.append(('uniform_bond_%d'%mu,x.ravel()))
 # Cover every inversion-inequivalent allowed momentum, each bond component and real quadrature.
 for a in range(m.lx):
  for b in range(m.ly):
   if (a,b)==(0,0) or (a,b)>((-a)%m.lx,(-b)%m.ly):continue
   phase=2*np.pi*(m.coords@np.array([a/m.lx,b/m.ly]))
   for quad,profile in [('cos',np.cos(phase)),('sin',np.sin(phase))]:
    if np.linalg.norm(profile)<1e-9:continue
    for mu in range(3):
     x=base.reshape(m.nc,3).copy();x[:,mu]+=.9*m.g*profile;structured.append(('q%d_%d_%s_mu%d'%(a,b,quad,mu),x.ravel()))
 for i,item in enumerate(structured):
  if i%PLAN['random_batches']==batch:starts.append(item)
 library=np.load(ROOT/'seed-library.npz');meta=json.loads((ROOT/'seed-library.json').read_text())
 for entry in meta:
  if entry['L']==p['L']:starts.append(('archived_'+entry['key'],library[entry['key']]*np.sqrt(p['lambda_value']/entry['lambda_value'])))
 return starts,rng

def relax_one(m,label,x0,long=False):
 try:
  r=m.relax(x0,maxiter=12000 if long else 3500);e,g=m.objective(r.x)
  if abs(g).max()>1e-6:
   r=m.relax(r.x,maxiter=18000);e,g=m.objective(r.x)
  return dict(label=label,energy=e,gradient_max=float(abs(g).max()),optimizer_success=bool(r.success),iterations=int(r.nit),detT=m.spectral(r.x)[3],error=''),r.x
 except Exception as e:return dict(label=label,energy=1e100,gradient_max=1e100,optimizer_success=False,iterations=0,detT=0,error=str(e)),np.full(m.nb,np.nan)

def search(index):
 batch=index%PLAN['random_batches'];p=POINTS[index//PLAN['random_batches']];m=model(p);out=ROOT/'search'/p['name']/('batch%d'%batch);out.mkdir(parents=True,exist_ok=False);start=time.time();starts,rng=seed_starts(m,p,batch);records=[];xs=[]
 for label,x0 in starts:
  r,x=relax_one(m,label,x0);records.append(r);xs.append(x)
  if len(records)%32==0:print(json.dumps(dict(point=p['name'],batch=batch,done=len(records),total=len(starts))),flush=True)
 # A second formulation eliminates X and searches SO(Nc), preserving the physical parity.
 so=[]
 for i in range(PLAN['reduced_SO_starts_per_batch']):
  x,diag=orthogonal_descent(m,rng,maxiter=600);so.append(diag);r,x=relax_one(m,'SO_random_%d'%i,x);records.append(r);xs.append(x)
 write_csv(out/'all-starts.csv',records);np.savez_compressed(out/'all-configurations.npz',x=np.array(xs));chosen,cr=select_candidates(m,records,xs,limit=4);np.savez_compressed(out/'candidates.npz',x=chosen)
 best=min(records,key=lambda r:r['energy']);assert best['energy']<1e90
 result=dict(status='COMPLETE',point=p,batch=batch,starts=len(starts),SO_starts=len(so),candidate_records=cr,best=best,SO_diagnostics=so,seconds=time.time()-start,software=dict(numpy=np.__version__))
 save(out/'result.json',result);print(json.dumps(dict(stage='search',point=p['name'],batch=batch,best=best,seconds=result['seconds'])),flush=True)

def collect_start_candidates(L,round_index):
 seeds=[]
 for p in POINTS:
  if p['L']!=L:continue
  candidates=[]
  if round_index==0:
   for batch in range(PLAN['random_batches']):
    root=ROOT/'search'/p['name']/('batch%d'%batch);r=json.loads((root/'result.json').read_text());assert r['status']=='COMPLETE';x=np.load(root/'candidates.npz')['x']
    candidates.extend((rr['energy'],xx,rr['label']) for rr,xx in zip(r['candidate_records'],x))
  else:
   root=ROOT/('cross%d'%(round_index-1))/p['name'];r=json.loads((root/'result.json').read_text());assert r['status']=='COMPLETE';x=np.load(root/'candidates.npz')['x'];candidates.extend((rr['energy'],xx,rr['label']) for rr,xx in zip(r['candidate_records'],x))
  candidates.sort(key=lambda a:a[0]);seen=set();n=0;mod=model(p)
  for e,x,label in candidates:
   sig=signature(mod,x)
   if sig in seen:continue
   seen.add(sig);seeds.append((p['name']+':'+label,p['lambda_value'],x));n+=1
   if n==2:break
 return seeds

def cross(index,round_index):
 p=POINTS[index];m=model(p);out=ROOT/('cross%d'%round_index)/p['name'];out.mkdir(parents=True,exist_ok=False);start=time.time();seeds=collect_start_candidates(p['L'],round_index);records=[];xs=[]
 for label,lam,x in seeds:
  r,x=relax_one(m,label,x*np.sqrt(p['lambda_value']/lam),long=True)
  if not r['error']:
   try:x=m.polish(x);r['energy'],g=m.objective(x);r['gradient_max']=float(abs(g).max())
   except Exception as err:r['polish_warning']=str(err)
  records.append(r);xs.append(x)
 write_csv(out/'all-cross-starts.csv',[{k:r.get(k,'') for k in ['label','energy','gradient_max','optimizer_success','iterations','detT','error','polish_warning']} for r in records]);np.savez_compressed(out/'all-cross-configurations.npz',x=np.array(xs))
 chosen,cr=select_candidates(m,records,xs,limit=4);np.savez_compressed(out/'candidates.npz',x=chosen);best=min(records,key=lambda r:r['energy']);bi=records.index(best);np.save(out/'equilibrium.npy',xs[bi])
 # If a lower unresolved candidate exists, keep it visible rather than replacing it by a higher convenient minimum.
 stable=[r for r in records if r['gradient_max']<2e-7 and not r['error']];assert stable
 best_stable=min(stable,key=lambda r:r['energy']);lower_unresolved=[r for r in records if r['energy']<best_stable['energy']-1e-8 and r['gradient_max']>=2e-7]
 result=dict(status='COMPLETE',point=p,round=round_index,candidate_records=cr,best=best,best_stable=best_stable,lower_unresolved=lower_unresolved,starts=len(seeds),seconds=time.time()-start)
 save(out/'result.json',result);print(json.dumps(dict(stage='cross',point=p['name'],round=round_index,best=best,lower_unresolved=len(lower_unresolved))),flush=True)

def fixed_q_check(m,x,hinv,q,hval):
 phase=2*np.pi*(m.coords@np.array([q[0]/m.lx,q[1]/m.ly]));checks=[]
 for label,values in [('cos',np.cos(phase)),('minus_sin',-np.sin(phase))]:
  source=np.zeros(m.nb);source[2::3]=values;expected=float(source@hinv@source/m.nc);dx=hinv@source
  plus=m.polish(m.relax(x+hval*dx,hval*source).x,hval*source);minus=m.polish(m.relax(x-hval*dx,-hval*source).x,-hval*source)
  measured=float(source@(plus-minus)/(2*hval*m.nc));checks.append(dict(quadrature=label,h=hval,expected=expected,measured=measured,error=abs(measured-expected),gradient_plus=float(abs(m.objective(plus,hval*source)[1]).max()),gradient_minus=float(abs(m.objective(minus,-hval*source)[1]).max())))
 return checks

def certify(index):
 p=POINTS[index];m=model(p);out=ROOT/'certified'/p['name'];out.mkdir(parents=True,exist_ok=False);last=PLAN['cross_rounds']-1;crossdir=ROOT/('cross%d'%last)/p['name'];r=json.loads((crossdir/'result.json').read_text());assert not r['lower_unresolved'],r['lower_unresolved'];x=np.load(crossdir/'equilibrium.npy');e,g=m.objective(x)
 previous=json.loads((ROOT/('cross%d'%(last-1))/p['name']/'result.json').read_text());improvement=(previous['best']['energy']-e)/m.nc
 assert abs(improvement)<=PLAN['cross_energy_tolerance_per_cell'],('SEARCH_NOT_STABLE',improvement)
 assert abs(g).max()<PLAN['gradient_tolerance']
 h,method,sym=m.hessian(x);eig=np.linalg.eigvalsh(h);assert eig[0]>PLAN['minimum_hessian_eigenvalue'];hinv=np.linalg.inv(h);rows,matrices=response_grid(m,x,hinv);peak=max(rows[1:],key=lambda z:z['chi_zz']);fixedidx=p['fixed_q'][0]*m.ly+p['fixed_q'][1]
 attempts=[];passed=None
 for fields in [(1e-4,5e-5,2.5e-5,5e-6,2.5e-6),(2.5e-6,1e-6,5e-7,2.5e-7),(2.5e-7,1e-7,5e-8,2.5e-8)]:
  c=validation(m,x,h,hinv,peak,106296000+index,field_strengths=fields);c['fixed_production_q_checks']=fixed_q_check(m,x,hinv,p['fixed_q'],min(fields));attempts.append(c)
  fixederror=max(z['error'] for z in c['fixed_production_q_checks']);fixedgrad=max(max(z['gradient_plus'],z['gradient_minus']) for z in c['fixed_production_q_checks'])
  if c['status']=='PASS' and fixederror<2e-5 and fixedgrad<1e-9:passed=c;break
 save(out/'validation-attempts.json',attempts);assert passed is not None,('RESPONSE_NOT_CERTIFIED',attempts[-1])
 field=x.reshape(m.nc,3);fft=np.fft.fft2(field.reshape(m.lx,m.ly,3),axes=(0,1))/m.nc;allrows=[];ferr=0;cerr=0
 for z,mat in zip(rows,matrices):
  d=np.array(mat['D_real'])+1j*np.array(mat['D_imag']);qchi=np.array(mat['chi_real'])+1j*np.array(mat['chi_imag']);ferr=max(ferr,float(abs(d-fft[z['n1'],z['n2']]).max()));other=matrices[((-z['n1'])%m.lx)*m.ly+(-z['n2'])%m.ly];oc=np.array(other['chi_real'])+1j*np.array(other['chi_imag']);cerr=max(cerr,float(abs(qchi.conj()-oc).max()));assert np.linalg.eigvalsh(qchi)[0]>0
  z.update(D_x_abs=float(abs(d[0])),D_y_abs=float(abs(d[1])),D_z_abs=float(abs(d[2])),D_vector_norm=float(np.linalg.norm(d)),Szz=float(m.nc*abs(d[2])**2),g_D_z_abs=float(m.g*abs(d[2])),g2_D_z_abs2=float(m.g**2*abs(d[2])**2),g2_chi_zz=float(m.g**2*z['chi_zz']));allrows.append(z)
 assert ferr<1e-12 and cerr<1e-8
 polar=m.physical_polar(x);reduced,_=m.reduced(polar);stationarity=float(abs(x-m.g*polar[m.a,m.b]).max());assert abs(reduced-e)<1e-8 and stationarity<2e-7
 lower=m.lower_bound();assert e>=lower-1e-9
 ground_code_parity=((-1)**m.nc)*m.spectral(x)[3];selected_code_parity=ground_code_parity*(-1 if m.spectral(x)[4] else 1);assert selected_code_parity==m.physical_rule['required_code_parity']
 translation=[]
 for a in range(m.lx):
  for b in range(m.ly):translation.append(m.objective(np.roll(np.roll(field.reshape(m.lx,m.ly,3),a,axis=0),b,axis=1).ravel())[0])
 tspread=max(translation)-min(translation);assert tspread<1e-8
 order=max(rows[1:],key=lambda z:z['D_z_abs']);vector=max(rows[1:],key=lambda z:z['D_vector_norm']);fixed=rows[fixedidx];mean=field.mean(0);std=field.std(0);branch='z_modulated' if order['D_z_abs']>1e-5 else 'xy_modulated' if vector['D_vector_norm']>1e-5 else 'uniform_xy_broken' if abs(mean[0]-mean[1])>1e-5 else 'uniform_symmetric'
 phase=np.exp(-2j*np.pi*(m.coords@np.array([peak['u'],peak['v']])));dz=complex(matrices[peak['index']]['D_real'][2],matrices[peak['index']]['D_imag'][2]);domains=dz*phase
 domain=dict(translation_energy_spread=tspread,cos_D_range=[float(domains.real.min()),float(domains.real.max())],minus_sin_D_range=[float(domains.imag.min()),float(domains.imag.max())],global_zero_field_response='domain_selection_not_a_unique_finite_scalar' if abs(dz)>1e-8 else 'no_translation_domain_jump_for_this_z_mode',scope='Single stable branch response; symmetry-restored/domain-selected susceptibility is not identified with this Hessian value')
 summary=dict(status='NUMERICALLY_VALIDATED_LOCAL_MINIMUM',L=p['L'],lambda_value=p['lambda_value'],energy_total=e,energy_per_cell=e/m.nc,physical_rule=m.physical_rule,ground_code_parity=ground_code_parity,selected_code_parity=selected_code_parity,equilibrium_detT_sign=m.spectral(x)[3],parity_excitation_applied=bool(m.spectral(x)[4]),mean_X=mean.tolist(),std_X=std.tolist(),branch=branch,order_z_max_nonzero=order['D_z_abs'],order_vector_max_nonzero=vector['D_vector_norm'],order_z_fixed=fixed['D_z_abs'],chi_zz_max_nonzero=peak['chi_zz'],chi_zz_fixed=fixed['chi_zz'],fixed_q=p['fixed_q'],order_q=[order['n1'],order['n2']],chi_peak_q=[peak['n1'],peak['n2']],mixed_chi_max=max(z['chi_eigenvalue_max'] for z in rows[1:]),gradient_max=float(abs(g).max()),hessian_minimum=float(eig[0]),hessian_method=method,inverse_residual=float(abs(h@hinv-np.eye(m.nb)).max()),stationarity_reduced_error=stationarity,independent_reduced_energy_error=float(abs(reduced-e)),last_cross_improvement_per_cell=improvement,rigorous_global_lower_bound=lower,rigorous_lower_bound_gap=e-lower,global_minimum_proven=False,response_validation='PASS',finest_field=passed['finest_field'],source_field_error=passed['finest_field_max_absolute_error'],fixed_q_field_error=max(z['error'] for z in passed['fixed_production_q_checks']),FFT_error=ferr,chi_conjugacy_error=cerr,domain_scope=domain)
 np.save(out/'equilibrium.npy',x);np.save(out/'hessian.npy',h);write_csv(out/'all-q.csv',allrows);write_csv(out/'equilibrium-bonds.csv',[dict(x=int(c[0]),y=int(c[1]),Xx=float(v[0]),Xy=float(v[1]),Xz=float(v[2])) for c,v in zip(m.coords,field)]);save(out/'result.json',dict(summary=summary,all_q_matrices=matrices,domain=domain));save(out/'validation.json',passed);print(json.dumps(summary),flush=True)

def aggregate():
 rows=[];missing=[];failures=[]
 for p in POINTS:
  path=ROOT/'certified'/p['name']/'result.json'
  if not path.exists():missing.append(p['name']);continue
  r=json.loads(path.read_text())['summary'];assert r['status']=='NUMERICALLY_VALIDATED_LOCAL_MINIMUM' and r['response_validation']=='PASS';rows.append(r)
 save(ROOT/'collection-status.json',dict(expected=82,certified=len(rows),missing=missing,all_local_gates_pass=not missing,global_optimality_claimed=False))
 if missing:raise RuntimeError('Incomplete certification; no final curve: '+str(missing))
 flat=[]
 for r in rows:
  r=dict(r)
  for key in ['physical_rule','mean_X','std_X','fixed_q','order_q','chi_peak_q','domain_scope']:r[key]=json.dumps(r[key],separators=(',',':'))
  flat.append(r)
 write_csv(ROOT/'scan-summary.csv',flat);print('ALL82_LOCAL_CERTIFICATIONS_COMPLETE',flush=True)


def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['search','cross','certify','aggregate']);ap.add_argument('--index',type=int);ap.add_argument('--round',type=int,default=0);args=ap.parse_args()
 try:
  if args.stage=='search':search(args.index)
  elif args.stage=='cross':cross(args.index,args.round)
  elif args.stage=='certify':certify(args.index)
  else:aggregate()
 except Exception as err:
  failure=dict(stage=args.stage,index=args.index,round=args.round,error=str(err),traceback=traceback.format_exc(),pbs_jobid=os.environ.get('PBS_JOBID'));save(ROOT/'failures'/('%s-r%d-i%s.json'%(args.stage,args.round,args.index)),failure);raise
if __name__=='__main__':main()
