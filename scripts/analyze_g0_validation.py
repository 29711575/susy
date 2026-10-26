#!/usr/bin/env python3
"""Independent CLI block recovery, joint phase-ratio jackknife and fixed g=0 gates.
Reblocking is always inside each chain. Every declared seed is required.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np

def energy(L):
    return -sum(abs(.5+np.exp(-2j*np.pi*x/L)+np.exp(-2j*np.pi*y/L)) for x in range(L) for y in range(L))

def schema(c):
    L,M=c['L'],c['M'];nc=L*L;w=c['omega'];var=1/(2*(100/9)*w);chi=-math.expm1(-w)/((100/9)*w*w)
    qs=[0,L*(L//2)+L//2]+([L] if L==4 else [])
    raw=[];spec=[]
    def add(name,ref,precision=None):spec.append((name,ref,precision))
    for q in qs:
        selfq=all((2*n)%L==0 for n in divmod(q,L))
        for mu in 'xyz':
            raw += [f'AmpA2_q{q}_{mu}',f'AmpA4_q{q}_{mu}',f'ChiD_q{q}_{mu}',f'ChiIconj_q{q}_{mu}']
            v=var/nc
            add(f'm2_q{q}_{mu}',v,'relative');add(f'm4_q{q}_{mu}',(3 if selfq else 2)*v*v)
            add(f'B4_q{q}_{mu}',3 if selfq else 2,'binder');add(f'Dmean_q{q}_{mu}',0)
        for mu in 'xyz':
            for nu in 'xyz':
                raw += [f'S_q{q}_{mu}{nu}',f'ChiRaw_q{q}_{mu}{nu}']
                add(f'S_q{q}_{mu}{nu}',var if mu==nu else 0,'relative' if mu==nu else None)
                add(f'chi_q{q}_{mu}{nu}',chi if mu==nu else 0,'relative' if mu==nu else None)
        raw += [f'ChiDRe_q{q}_z',f'ChiDIm_q{q}_z',f'ChiIRe_q{q}_z',f'ChiIIm_q{q}_z',f'ChiCosRaw_q{q}_zz',f'ChiSinRaw_q{q}_zz']
        add(f'cos_q{q}',chi if selfq else chi/2)
        add(f'sin_q{q}',0 if selfq else chi/2)
    for j in range(M+1):
        raw += [f'PhDiag{s}_j{j}' for s in ['X','X2','CenterC','OtherBondC']]
        add(f'node_X_j{j}',0);add(f'node_X2_j{j}',var)
        add(f'node_center_C_j{j}',var*math.exp(-w*c['dt']*abs(j-M//2)))
        add(f'node_other_bond_C_j{j}',0)
    raw += ['majorana_energy','PhDiagLogMatter']
    add('matter_energy',energy(L));add('log_matter',-2*energy(L))
    return list(dict.fromkeys(raw)),spec,qs

def derive(t,names,c,qs):
    # Common COMPLEX denominator first; all nonlinear operations happen after.
    z=t[0];r={k:t[i+2]/z for i,k in enumerate(names)};nc=c['L']**2;out=[]
    for q in qs:
        for mu in 'xyz':
            m2=r[f'AmpA2_q{q}_{mu}'];m4=r[f'AmpA4_q{q}_{mu}']
            out += [m2,m4,m4/(m2*m2),r[f'ChiD_q{q}_{mu}']]
        for mu in 'xyz':
            for nu in 'xyz':
                out += [r[f'S_q{q}_{mu}{nu}'],r[f'ChiRaw_q{q}_{mu}{nu}']-nc*r[f'ChiD_q{q}_{mu}']*r[f'ChiIconj_q{q}_{nu}']]
        out += [r[f'ChiCosRaw_q{q}_zz']-nc*r[f'ChiDRe_q{q}_z']*r[f'ChiIRe_q{q}_z'],
                r[f'ChiSinRaw_q{q}_zz']-nc*r[f'ChiDIm_q{q}_z']*r[f'ChiIIm_q{q}_z']]
    for j in range(c['M']+1):
        out += [r[f'PhDiag{s}_j{j}'] for s in ['X','X2','CenterC','OtherBondC']]
    out += [r['majorana_energy'],r['PhDiagLogMatter']]
    return np.asarray(out)

def jackknife(chains,names,c,qs,factor):
    groups=[]
    for a in chains:
        groups += [a[i:i+factor].sum(axis=0) for i in range(0,len(a),factor)]
    groups=np.array(groups);total=groups.sum(axis=0);theta=derive(total,names,c,qs)
    # Weighted delete-group pseudovalues retain short remainder groups instead
    # of throwing away observations at coarser blocking levels.
    weights=groups[:,1].real/total[1].real
    deleted=np.array([derive(total-b,names,c,qs) for b in groups])
    pseudo=(theta[None,:]-(1-weights[:,None])*deleted)/weights[:,None]
    avg=(weights[:,None]*pseudo).sum(axis=0)
    se=np.sqrt(((weights[:,None]*(pseudo.real-avg.real)**2).sum(axis=0))/(len(groups)-1))
    sei=np.sqrt(((weights[:,None]*(pseudo.imag-avg.imag)**2).sum(axis=0))/(len(groups)-1))
    return theta,se,sei,len(groups)

def self_test():
    c=dict(L=2,M=2,dt=1.,omega=.3);names,spec,qs=schema(c)
    # Complex signed finite ensemble tests common denominator and nonlinear
    # connection. Distinct configurations distinguish true fourth moments.
    weights=np.array([1+.2j,-.13+.07j,.8-.11j,.3+.04j]);d=np.array([.2,.7,1.3,.4]);integ=np.array([.9,.4,1.7,.6])
    a=np.zeros((4,len(names)+2),complex);a[:,0]=weights;a[:,1]=1
    for q in qs:
        for mu in 'xyz':
            for key,val in [(f'AmpA2_q{q}_{mu}',d*d),(f'AmpA4_q{q}_{mu}',d**4),(f'ChiD_q{q}_{mu}',d),(f'ChiIconj_q{q}_{mu}',integ)]:a[:,2+names.index(key)]=weights*val
            for nu in 'xyz':a[:,2+names.index(f'ChiRaw_q{q}_{mu}{nu}')]=weights*4*d*integ
    val=derive(a.sum(axis=0),names,c,qs);idx={s[0]:i for i,s in enumerate(spec)};z=weights.sum()
    assert abs(val[idx['B4_q0_z']]-(weights@d**4/z)/(weights@d**2/z)**2)<1e-13
    assert abs(val[idx['chi_q0_zz']]-4*((weights@(d*integ))/z-(weights@d/z)*(weights@integ/z)))<1e-13
    assert abs((weights@d**4/z)-(weights@d**2/z)**2)>.01
    th,se,_,_=jackknife([a[:2],a[2:]],names,c,qs,1)
    reps=np.array([derive(a.sum(axis=0)-b,names,c,qs) for b in a]);ref=np.sqrt(3/4*((reps.real-reps.real.mean(axis=0))**2).sum(axis=0))
    assert np.max(abs(se-ref))<1e-11
    # A three-row chain followed by a one-row chain must give three groups,
    # retaining both short remainders without merging across the boundary.
    a=np.ones((4,len(names)+2),complex);a[:,2+names.index('majorana_energy')]=[2,4,8,16]
    th,se,_,ng=jackknife([a[:3],a[3:]],names,c,qs,2)
    ei=idx['matter_energy']
    assert ng==3 and abs(th[ei]-7.5)<1e-13 and abs(se[ei]-math.sqrt(14.125))<1e-12
    return 'PASS'

def analyze(root,out):
    plan=json.loads((root/'plan.json').read_text());gates=plan['gates'];out.mkdir(parents=True,exist_ok=True)
    report={'statistics_self_test':self_test(),'gates':gates,'cases':[],'chain_records':[]}
    detail=[];levels=[];data_by_case={}
    for c in plan['cases']:
        names,spec,qs=schema(c);chains=[];problems=[]
        for ch in [ch for ch in plan['chains'] if ch['case']==c['name']]:
            pre=root/'runs'/('chain%03d'%ch['index']);rp=Path(str(pre)+'_results.json');bp=Path(str(pre)+'_blocks.csv')
            record={'index':ch['index'],'case':c['name'],'seed':ch['seed']}
            try:
                j=json.loads(rp.read_text())
                if j['status']=='running':raise FileNotFoundError('chain still running')
                assert j['status']=='complete';assert j['numerical_failures']==0
                assert j['seed']==ch['seed'] and j['measurement_sweeps']==c['sweeps'] and j['samples']==c['sweeps']
                for k in ['Lx','Ly','g','mass','Jx','Jy','Jz','shear','intervals','dt']:
                    assert abs(j[k]-ch['parameters'][k])<1e-12
                assert j['trial_shift']==[0,0,0] and j['wilson_loops']==[1,1] and abs(j['omega_ph']-c['omega'])<1e-14
                eff={k:j[k] for k in ['weight_engine','defer_phase','prefix_cache','cayley_cache','center_green_fast','numerics_profile']}
                expected=['reference',0,0,0,0,'reference-v1'] if c['numerics']=='reference' else ['combined',1,1,1,1,'production-v1']
                assert list(eff.values())==expected
                assert j['accepted_measurement']==j['proposals_measurement'] and j['accepted_warmup']==j['proposals_warmup']
                assert j['q_mode']['every']==c['qon'] and abs(j['q_mode']['beta']-c['beta'])<1e-14
                for stage in j['move_counts'].values():assert stage['accepted']==stage['proposed'] and stage['first_failures']==0
                if c['proposal']=='bond-guided':assert abs(j['guided_bond']['beta']-c['beta'])<1e-14
                if c['proposal']=='bond-guided':assert j['guided_bond']['force_by_type']==[0,0,0]
                a=[]
                with bp.open() as f:
                    for row in csv.DictReader(f):
                        z=complex(float(row['sum_phase_real']),float(row['sum_phase_imag']));n=int(row['count'])
                        assert abs(z-n)<2e-9*n
                        a.append([z,n]+[complex(float(row['sum_phase_'+k+'_real']),float(row['sum_phase_'+k+'_imag'])) for k in names])
                a=np.asarray(a);assert len(a)==c['sweeps']//20;assert np.isfinite(a).all()
                total=a.sum(axis=0);max_recovery=max(abs(total[i+2]/total[0]-complex(*j['observables'][name]['value'])) for i,name in enumerate(names))
                assert max_recovery<2e-10
                record.update(status='PASS',numerics=eff,raw_recovery_max_error=float(max_recovery),max_block_phase_error=float(np.max(abs(a[:,0]/a[:,1]-1))),accepted_proposals=j['accepted_measurement'],proposed_proposals=j['proposals_measurement'],q_proposals=j['move_counts']['q_measurement']['proposed'],elapsed=j['elapsed_seconds'],
                  result_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),blocks_sha256=hashlib.sha256(bp.read_bytes()).hexdigest())
                chains.append(a)
            except Exception as exc:
                record.update(status='BLOCKED' if isinstance(exc,FileNotFoundError) else 'FAIL',error=repr(exc));problems.append(record)
            report['chain_records'].append(record)
        if problems:
            report['cases'].append({'name':c['name'],'tier':c['tier'],'status':'FAIL' if any(p['status']=='FAIL' for p in problems) else 'BLOCKED','problems':problems});continue
        all_stats=[]
        for factor in gates['reblock_factors']:
            th,se,sei,nb=jackknife(chains,names,c,qs,factor);all_stats.append((th,se,sei))
            for i,(name,ref,precision) in enumerate(spec):levels.append(dict(case=c['name'],factor=factor,block_sweeps=20*factor,blocks=nb,observable=name,value=float(th[i].real),stderr=float(se[i]),stderr_imag=float(sei[i])))
        th=all_stats[0][0];err=np.max([s[1] for s in all_stats],axis=0);erri=np.max([s[2] for s in all_stats],axis=0)
        nfail=nprecision=nplateau=0;max_rel=0;max_b4=0
        for i,(name,ref,precision) in enumerate(spec):
            realpass=abs(th[i].real-ref)<=gates['sigma']*err[i]+gates['deterministic_absolute_floor']
            imagpass=abs(th[i].imag)<=gates['sigma']*erri[i]+gates['deterministic_absolute_floor']
            precise=True;plateau=True
            ratio=float(all_stats[-1][1][i]/all_stats[-2][1][i]) if all_stats[-2][1][i]>1e-12 else 1.
            if precision:
                plateau=gates['plateau_min_ratio']<=ratio<=gates['plateau_max_ratio'] and min(math.ceil(len(a)/gates['reblock_factors'][-1]) for a in chains)>=gates['minimum_final_blocks_per_chain']
                precise=err[i]<=(gates['B4_absolute_se'] if precision=='binder' else gates['relative_se']*abs(ref))
                if precision=='binder':max_b4=max(max_b4,float(err[i]))
                else:max_rel=max(max_rel,float(err[i]/abs(ref)))
            status='FAIL' if not(realpass and imagpass) else 'INCONCLUSIVE' if not(precise and plateau) else 'PASS'
            nfail+=status=='FAIL';nprecision+=not precise;nplateau+=not plateau
            detail.append(dict(case=c['name'],observable=name,reference=ref,value=float(th[i].real),imag=float(th[i].imag),stderr=float(err[i]),stderr_imag=float(erri[i]),z=float(abs(th[i].real-ref)/max(err[i],1e-30)),status=status,precision_met=precise,plateau_met=plateau,last_se_ratio=ratio))
        status='FAIL' if nfail else 'INCONCLUSIVE' if nprecision or nplateau or c['tier']=='smoke' else 'PASS'
        report['cases'].append(dict(name=c['name'],status=status,tier=c['tier'],samples=c['sweeps']*4,failed_metrics=nfail,imprecise_metrics=nprecision,plateau_unresolved_metrics=nplateau,max_main_relative_se=max_rel,max_B4_se=max_b4))
        data_by_case[c['name']]=(names,chains)
    # Matched reference/production chain streams: compare raw phase numerators
    # as well as final coordinates without relying only on aggregate error bars.
    comparisons=[]
    for left,right in [('complex-reference','complex-production'),('fine-reference','fine-production')]:
        if left in data_by_case and right in data_by_case:
            nl,cl=data_by_case[left];nr,cr=data_by_case[right];assert nl==nr
            err=max(float(np.max(abs(a-b))) for a,b in ((a[:min(len(a),len(b))],b[:min(len(a),len(b))]) for a,b in zip(cl,cr)))
            comparisons.append(dict(left=left,right=right,max_raw_block_sum_difference=err,status='PASS' if err<2e-8 else 'FAIL'))
    report['matched_numerics_comparisons']=comparisons
    report['cli_recovery_status']='PASS' if len(report['chain_records'])==len(plan['chains']) and all(r['status']=='PASS' for r in report['chain_records']) else 'BLOCKED' if any(r['status']=='BLOCKED' for r in report['chain_records']) else 'FAIL'
    statuses=[c['status'] for c in report['cases']]+[c['status'] for c in comparisons]
    report['main_production_status']='FAIL' if any(c['status']=='FAIL' for c in report['cases'] if c.get('tier')=='full') else 'BLOCKED' if any(c['status']=='BLOCKED' for c in report['cases'] if c.get('tier')=='full') else 'INCONCLUSIVE' if any(c['status']=='INCONCLUSIVE' for c in report['cases'] if c.get('tier')=='full') else 'PASS'
    report['production_status']='FAIL' if 'FAIL' in statuses else 'BLOCKED' if 'BLOCKED' in statuses else 'INCONCLUSIVE' if 'INCONCLUSIVE' in statuses else 'PASS'
    for filename,rows in [('metrics.csv',detail),('blocking.csv',levels)]:
        if rows:
            with (out/filename).open('w') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['chain_records','gates']},indent=2))
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument('package',type=Path,nargs='?');p.add_argument('--output',type=Path);p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:print(self_test());return 0
    if a.package is None:p.error('package is required unless --self-test is used')
    r=analyze(a.package,a.output or a.package/'analysis')
    return 1 if r['production_status']=='FAIL' or r['cli_recovery_status']=='FAIL' else 2 if r['production_status']!='PASS' or r['cli_recovery_status']!='PASS' else 0
if __name__=='__main__':raise SystemExit(main())
