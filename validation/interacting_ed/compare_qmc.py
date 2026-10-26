#!/usr/bin/env python3
"""Joint complex phase-ratio block jackknife: production CLI versus full ED."""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np

RAW=['X_'+a for a in 'xyz']+['X2_'+a for a in 'xyz']+['majorana_energy']
RAW += [stem+str(b) for b in range(12) for stem in ['BondK_b','BondXK_b']]
RAW += [stem+a for stem in ['BondK_','BondXK_'] for a in 'xyz']
for q in [0,3]:RAW += [f'S_q{q}_zz',f'AmpA2_q{q}_z',f'AmpA4_q{q}_z',f'ChiDRe_q{q}_z']
RAW += ['ChiIRe_q3_z','ChiCosRaw_q3_zz','ChiD_q3_z','ChiIconj_q3_z','ChiRaw_q3_zz']
NAMES=RAW[:37]+[s for q in [0,3] for s in [f'S_q{q}_zz',f'm2_q{q}_z',f'm4_q{q}_z',f'B4_q{q}_z',f'Dmean_q{q}_z']]+['chi_cos_q3','chi_complex_q3']
assert len(RAW[:37])==37 and RAW[36]=='BondXK_z'

def derive(total):
    r={name:total[i+2]/total[0] for i,name in enumerate(RAW)};out=[r[n] for n in RAW[:37]]
    for q in [0,3]:
        m2=r[f'AmpA2_q{q}_z'];m4=r[f'AmpA4_q{q}_z'];out += [r[f'S_q{q}_zz'],m2,m4,m4/(m2*m2),r[f'ChiDRe_q{q}_z']]
    out += [r['ChiCosRaw_q3_zz']-4*r['ChiDRe_q3_z']*r['ChiIRe_q3_z'],r['ChiRaw_q3_zz']-4*r['ChiD_q3_z']*r['ChiIconj_q3_z']]
    return np.array(out)

def jk(chains,factor):
    groups=np.array([a[i:i+factor].sum(axis=0) for a in chains for i in range(0,len(a),factor)]);t=groups.sum(axis=0);theta=derive(t)
    w=groups[:,1].real/t[1].real;deleted=np.array([derive(t-b) for b in groups]);pseudo=(theta[None,:]-(1-w[:,None])*deleted)/w[:,None];mu=np.sum(w[:,None]*pseudo,axis=0)
    se=np.sqrt(np.sum(w[:,None]*(pseudo.real-mu.real)**2,axis=0)/(len(groups)-1));si=np.sqrt(np.sum(w[:,None]*(pseudo.imag-mu.imag)**2,axis=0)/(len(groups)-1))
    return theta,se,si,pseudo

def precision_limit(name,value,gates):
    if name.startswith('X2_') or name.startswith(('S_','m2_','chi_')):return gates['MC_second_and_response_relative_se']*abs(value)
    if name.startswith('B4'):return gates['MC_B4_absolute_se']
    if name.startswith('m4'):return gates['MC_fourth_relative_se']*abs(value)
    if name.startswith('BondXK'):return gates['MC_bond_XK_absolute_se']
    if name.startswith('BondK'):return gates['MC_bond_K_absolute_se']
    if name=='majorana_energy':return gates['MC_matter_energy_absolute_se']
    return gates['MC_mean_X_absolute_se']

def acf(y):
    z=y-y.mean();n=len(z);variance=np.dot(z,z)/n
    if variance<1e-24:return [1.]+[0.]*min(n//4,100),1.
    rho=[float(np.dot(z[:n-k],z[k:])/(n-k)/variance) for k in range(min(n//4,100)+1)]
    integ=1.
    for k in range(1,len(rho)-1,2):
        pair=rho[k]+rho[k+1]
        if pair<=0:break
        integ+=2*pair
    return rho,integ

def analyze(package,edroot,continuous,out):
    plan=json.loads((package/'plan.json').read_text());gates=plan['config']['gates'];out.mkdir(parents=True,exist_ok=True);rows=[];blocking=[];auto=[];cases=[];records=[];estimates={};comparegroups={}
    for case in dict.fromkeys(c['case'] for c in plan['chains']):
        cc=[c for c in plan['chains'] if c['case']==case];pt=cc[0]['case_parameters'];g,tau,dt=pt['g'],pt['tau'],pt['dt'];chains=[];problems=[]
        for c in cc:
            pre=package/'runs'/f"chain{c['index']:03d}";rp=Path(str(pre)+'_results.json');bp=Path(str(pre)+'_blocks.csv');rec=dict(index=c['index'],case=case,seed=c['seed'])
            try:
                j=json.loads(rp.read_text())
                if j['status']=='running':raise FileNotFoundError('chain still running')
                assert j['status']=='complete' and j['numerical_failures']==0
                for key in ['Lx','Ly','shear','g','mass','Jx','Jy','Jz','intervals','dt','seed']:
                    assert abs(j[key]-c['parameters'][key])<1e-12
                assert j['omega_ph']==1 and j['trial_shift']==[0,0,0] and j['wilson_loops']==[1,1]
                assert j['samples']==c['parameters']['sweeps'] and j['block_samples']==20 and j['bond_diagnostics']==1
                cfg={k:j[k] for k in ['weight_engine','defer_phase','prefix_cache','cayley_cache','center_green_fast','numerics_profile']}
                expected=['reference',0,0,0,0,'reference-v1'] if pt['numerics']=='reference' else ['combined',1,1,1,1,'production-v1'];assert list(cfg.values())==expected
                a=[]
                with bp.open() as f:
                    for v in csv.DictReader(f):a.append([complex(float(v['sum_phase_real']),float(v['sum_phase_imag'])),float(v['count'])]+[complex(float(v['sum_phase_'+name+'_real']),float(v['sum_phase_'+name+'_imag'])) for name in RAW])
                a=np.array(a);assert np.isfinite(a).all() and len(a)==j['samples']//20
                t=a.sum(axis=0);err=max(abs(t[i+2]/t[0]-complex(*j['observables'][name]['value'])) for i,name in enumerate(RAW));assert err<2e-10
                rec.update(status='PASS',samples=j['samples'],numerics=cfg,guided_reference=j['guided_bond'],q_mode=j['q_mode'],average_phase=[float((t[0]/t[1]).real),float((t[0]/t[1]).imag)],raw_reconstruction_error=float(err),numerical_failures=j['numerical_failures'],acceptance=j['acceptance'],elapsed=j['elapsed_seconds'],result_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),blocks_sha256=hashlib.sha256(bp.read_bytes()).hexdigest());chains.append(a)
            except Exception as exc:
                rec.update(status='BLOCKED' if isinstance(exc,FileNotFoundError) else 'FAIL',error=repr(exc));problems.append(rec)
            records.append(rec)
        if problems:cases.append(dict(case=case,status='FAIL' if any(p['status']=='FAIL' for p in problems) else 'BLOCKED',problems=problems));continue
        # Full quantum finite-dt ED, with a cutoff check independent of QMC.
        er=[]
        try:
            for cut in [4,6,8]:er.append(json.loads((edroot/f'g{g}-tau{tau}-dt{dt}-cut{cut}.json').read_text()))
        except FileNotFoundError as exc:cases.append(dict(case=case,status='BLOCKED',error=str(exc)));continue
        ct=json.loads((continuous/f"ed-g{str(g).replace('.','p')}-cut8.json").read_text());cr=next(r for r in ct['results'] if r['tau']==tau)
        result=[]
        for f in gates['reblock_factors']:
            z,se,si,pseudo=jk(chains,f);result.append((z,se,si,pseudo))
            for i,name in enumerate(NAMES):blocking.append(dict(case=case,observable=name,block_sweeps=20*f,value=float(z[i].real),stderr=float(se[i]),stderr_imag=float(si[i])))
        theta=result[0][0];se=np.max([s[1] for s in result],axis=0);si=np.max([s[2] for s in result],axis=0)
        # Complex phase resolution, using the same independent-chain groups.
        phases=np.array([b[0]/b[1] for a in chains for b in a]);phase=phases.mean();phase_errors=[]
        for f in gates['reblock_factors']:
            pg=np.array([a[k:k+f].sum(axis=0)[:2] for a in chains for k in range(0,len(a),f)]);pw=pg[:,1].real/pg[:,1].real.sum();pv=pg[:,0]/pg[:,1]
            phase_errors.append(float(np.sqrt(np.sum(pw*abs(pv-phase)**2)/(len(pg)-1))))
        phase_se=max(phase_errors);phaseok=abs(phase)>gates['phase_resolution_sigma']*phase_se
        for name in ['majorana_energy','m2_q3_z','chi_cos_q3']:
            i=NAMES.index(name);start=0
            for c,a in zip(cc,chains):
                influence=result[0][3][start:start+len(a),i].real;rho,tint=acf(influence);start+=len(a)
                auto.append(dict(case=case,seed=c['seed'],observable=name,base_block_sweeps=20,integrated_correlation_in_base_blocks=tint,acf=rho))
        failed=imprecise=unresolved=edbad=0
        for i,name in enumerate(NAMES):
            ename='chi_cos_q3' if name=='chi_complex_q3' else name
            ev=er[-1]['observables'][ename];ec=cr['observables'][ename];fine=abs(ev-er[-2]['observables'][ename]);coarse=abs(er[-2]['observables'][ename]-er[-3]['observables'][ename]);etol=gates['ED_B4_truncation_absolute'] if name.startswith('B4') else gates['ED_truncation_abs_relative']*max(1,abs(ev));edok=fine<=etol and (fine<=coarse or fine<1e-10)
            te=2*fine;num=cr['propagation_segmentation_residual'];bound=gates['MC_sigma']*se[i]+te+gates['MC_absolute_floor'];delta=float(theta[i].real-ev);consistent=abs(delta)<=bound and abs(theta[i].imag)<=gates['MC_sigma']*si[i]+gates['MC_absolute_floor']
            precise=se[i]<=precision_limit(name,ev,gates);ratio=result[-1][1][i]/result[-2][1][i] if result[-2][1][i]>1e-12 else 1
            plateau=gates['plateau_ratio'][0]<=ratio<=gates['plateau_ratio'][1] and min(math.ceil(len(a)/gates['reblock_factors'][-1]) for a in chains)>=gates['minimum_final_groups_per_chain']
            status='INCONCLUSIVE' if not(edok and phaseok and precise and plateau) else 'PASS' if consistent else 'FAIL'
            failed+=status=='FAIL';imprecise+=not precise;unresolved+=not plateau;edbad+=not edok
            rows.append(dict(case=case,g=g,tau=tau,dt=dt,observable=name,ED_finite_dt=ev,ED_continuous=ec,QMC=float(theta[i].real),QMC_imag=float(theta[i].imag),MC_standard_error=float(se[i]),MC_imag_standard_error=float(si[i]),difference_finite_dt=delta,difference_continuous=float(theta[i].real-ec),ED_truncation_error_estimate=te,ED_continuous_propagation_residual=num,ED_ground_eigensolver_residual=cr['ground_residual'],ED_finite_dt_minus_continuous=ev-ec,finite_tau_minus_ground=cr['finite_tau_minus_ground'].get(ename),statistical_acceptance_tolerance=bound,precision_limit=precision_limit(name,ev,gates),last_SE_ratio=float(ratio),ED_cutoff_converged=edok,MC_precision_met=precise,MC_blocking_plateau=plateau,phase_resolved=bool(phaseok),six_sigma_consistent=bool(consistent),status=status))
        status='FAIL' if failed else 'INCONCLUSIVE' if imprecise or unresolved or edbad or not phaseok else 'PASS'
        cases.append(dict(case=case,status=status,failed_metrics=failed,imprecise_metrics=imprecise,plateau_unresolved_metrics=unresolved,ED_cutoff_unresolved_metrics=edbad,samples=sum(c['parameters']['sweeps'] for c in cc),average_phase=[float(phase.real),float(phase.imag)]))
        estimates[case]=(pt,theta,se)
    # Compare algorithms/policies at identical physics and dt using disjoint seeds.
    cross=[]
    available=list(estimates)
    for ia,left in enumerate(available):
        lp,lv,ls=estimates[left]
        for right in available[ia+1:]:
            rp,rv,rs=estimates[right]
            if any(lp[k]!=rp[k] for k in ['g','tau','dt']):continue
            reliable=all(next(c for c in cases if c['case']==case)['status']=='PASS' for case in [left,right])
            for i,name in enumerate(NAMES):
                delta=float(lv[i].real-rv[i].real);error=float(math.hypot(ls[i],rs[i]));ok=abs(delta)<=gates['MC_sigma']*error+gates['MC_absolute_floor']
                cross.append(dict(left=left,right=right,observable=name,difference=delta,independent_chain_standard_error=error,status='INCONCLUSIVE' if not reliable else 'PASS' if ok else 'FAIL'))
    # Independent-seed weighted extrapolation. Error of nonlinear quantities
    # has already been propagated by joint delete-block phase jackknife.
    extrap=[]
    chosen=[(case,*v) for case,v in estimates.items() if v[0]['g']==.4 and v[0]['tau']==2 and v[0]['proposal']=='bond-guided' and v[0]['numerics']=='production']
    if len(chosen)==3:
        chosen=sorted(chosen,key=lambda v:v[1]['dt']);x=np.array([v[1]['dt']**2 for v in chosen]);A=np.column_stack([np.ones(3),x])
        ct=json.loads((continuous/'ed-g0p4-cut8.json').read_text())['results'][0]
        for i,name in enumerate(NAMES):
            y=np.array([v[2][i].real for v in chosen]);errors=np.maximum([v[3][i] for v in chosen],1e-12);W=np.diag(1/errors**2);cov=np.linalg.inv(A.T@W@A);beta=cov@A.T@W@y
            quart=np.linalg.solve(np.column_stack([np.ones(3),x,x*x]),y);small=np.linalg.solve(A[:2],y[:2]);ref=ct['observables']['chi_cos_q3' if name=='chi_complex_q3' else name]
            extrap.append(dict(observable=name,ED_continuous=ref,QMC_continuous_intercept=float(beta[0]),MC_extrapolation_standard_error=float(math.sqrt(cov[0,0])),difference=float(beta[0]-ref),dt_fit_model_spread=float(max(abs(beta[0]-quart[0]),abs(beta[0]-small[0]))),quadratic_dt_fit_chi2=float(np.sum(((A@beta-y)/errors)**2)),fit_dof=1,note='model spread includes sampling fluctuations; finite-dt ED correction is tabulated separately'))
    for filename,values in [('QMC-ED-observables.csv',rows),('blocking.csv',blocking),('continuum-extrapolation.csv',extrap),('algorithm-comparisons.csv',cross)]:
        if values:
            with (out/filename).open('w') as f:w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
    (out/'autocorrelation.json').write_text(json.dumps(auto,indent=2)+'\n')
    statuses=[c['status'] for c in cases]+[r['status'] for r in cross];summary=dict(status='FAIL' if 'FAIL' in statuses else 'BLOCKED' if 'BLOCKED' in statuses else 'INCONCLUSIVE' if 'INCONCLUSIVE' in statuses else 'PASS',cases=cases,chains=records,expected_chains=len(plan['chains']),completed_chains=sum(r['status']=='PASS' for r in records),scope='real production sampler and CLI; all twelve quantum phonons in independent ED reference; common complex phase denominator before connected subtraction and B4; no total-energy QMC estimator claimed')
    (out/'C-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='chains'},indent=2));return summary

def main():
 p=argparse.ArgumentParser();p.add_argument('package',type=Path);p.add_argument('--discrete-ed',type=Path,required=True);p.add_argument('--continuous-ed',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();s=analyze(a.package,a.discrete_ed,a.continuous_ed,a.output);return 1 if s['status']=='FAIL' else 0 if s['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
