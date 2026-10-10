#!/usr/bin/env python3
"""Read-only complex fermion analysis using the existing joint block jackknife."""
import argparse,datetime,json,math,pathlib
from dataclasses import replace
import numpy as np
from convergence_analysis import joint_estimate,estimates_json,fixed_observables
from analyze_phase_tau_scan import read_phase_chain,phase_transform,chain_row,save_csv,jsonable

def observables(sums,count,names,metadata):
    out=fixed_observables(sums,count,names,metadata)
    out['phase_magnitude']=abs(out['phase'])
    for name,value in zip(names,sums[1:]/sums[0]):
        if name.startswith('F_'):out[name]=value;out[name+'_magnitude']=abs(value)
    return out

def analytic(c):
    dx,dy=c['cell_displacement'];a,b=c['sublattices'];z=complex(int(a==b and dx==dy==0))
    if a!=b:
        z=0j
        for nx in range(4):
            for ny in range(4):
                qx,qy=2*np.pi*nx/4,2*np.pi*ny/4
                f=.5+np.exp(-1j*qx)+np.exp(-1j*qy)
                z+=np.exp(-1j*(qx*dx+qy*dy))*(1j*f if a=='A' else -1j*f.conjugate())/abs(f)/16
    if c['kind']=='string':z*=(-1)**(dx+dy+int(a=='B' and b=='A'))
    return z

def rows_for(chains,block,scope,output,tag):
    blocked=[c.reblock(block) for c in chains]
    unresolved=''
    try:est=joint_estimate(blocked,transform=observables)
    except ValueError as exc:
        if 'unresolved phase denominator' not in str(exc):raise
        unresolved=str(exc);est=joint_estimate(blocked,transform=phase_transform)
    ee=estimates_json(est);d=chains[0].metadata
    channels={x['name']:x for x in d['fermions']['channels']}
    np.savez_compressed(output/f'{tag}_block{block}_covariance.npz',names=np.array(est['names']),values=est['values'],covariance=est['covariance'],phase_covariance=est['phase_covariance'])
    rows=[]
    for name,e in ee.items():
        if name.endswith('_magnitude'):continue
        z=complex(*e['value']);mag=ee.get(name+'_magnitude')
        se=math.hypot(e['stderr_real'],e['stderr_imag'])
        row=dict(scope=scope,g=d['g'],tau_total=d['intervals']*d['dt'],intervals=d['intervals'],seeds=','.join(str(c.seed) for c in chains),block_sweeps=block,observable=name,Re=z.real,Im=z.imag,magnitude=abs(z),SE_Re=e['stderr_real'],SE_Im=e['stderr_imag'],SE=se,cov_Re_Im=e['cov_real_imag'],SE_magnitude=mag['stderr_real'] if mag else None,phase_resolved_5sigma=est['phase_resolved_5sigma'],samples=est['samples'],blocks_per_chain=','.join(map(str,est['chain_blocks'])),kind='',dx=None,dy=None,sublattices='',analytic_Re=None,analytic_Im=None,analytic_abs_difference=None,physical_status=unresolved or 'defined')
        if name in channels:
            c=channels[name];row.update(kind=c['kind'],dx=c['cell_displacement'][0],dy=c['cell_displacement'][1],sublattices=c['sublattices'])
            if d['g']==0:
                ref=analytic(c);row.update(analytic_Re=ref.real,analytic_Im=ref.imag,analytic_abs_difference=abs(z-ref))
        rows.append(row)
    phase=ee['phase'];z=complex(*phase['value']);se=math.hypot(phase['stderr_real'],phase['stderr_imag'])
    diag=dict(g=d['g'],tau_total=d['intervals']*d['dt'],intervals=d['intervals'],seeds=','.join(str(c.seed) for c in chains),block_sweeps=block,Re=z.real,Im=z.imag,magnitude=abs(z),SE=se,SNR=abs(z)/se if se else None,positive_count=sum(c.metadata['phase_diagnostics']['positive_count'] for c in chains),negative_count=sum(c.metadata['phase_diagnostics']['negative_count'] for c in chains),complex_count=sum(c.metadata['phase_diagnostics']['complex_count'] for c in chains),acceptance=sum(c.metadata['accepted_measurement'] for c in chains)/sum(c.metadata['proposals_measurement'] for c in chains),physical_status=unresolved or 'defined')
    return rows,diag

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=pathlib.Path,required=True);ap.add_argument('--output-dir',type=pathlib.Path);args=ap.parse_args()
    root=pathlib.Path(__file__).resolve().parents[1];run=args.run_dir.resolve();manifest=json.loads((run/'manifest.json').read_text());plan=manifest['plan']
    chains=[read_phase_chain(pathlib.Path(c['folder'])/'run_results.json') for c in manifest['commands']]
    assert len(chains)==len(plan['jobs'])
    for c,job in zip(chains,plan['jobs']):
        d=c.metadata
        for k,v in [('Lx',4),('Ly',4),('Jx',1),('Jy',1),('Jz',.5),('g',job['g']),('mass',4),('omega_ph',.5),('dt',.1),('intervals',job['intervals']),('seed',job['seed']),('warmup_sweeps',plan['warmup']),('measurement_sweeps',plan['sweeps']),('weight_engine','combined')]:assert d[k]==v,(k,d[k],v)
        assert d['fermions_enabled'] and d['numerical_failures']==0 and len(d['fermions']['channels'])==128
        assert d['trial_shift']==[0,0,0]
    output=args.output_dir or root/'results'/plan['name'];output.mkdir(parents=True,exist_ok=True)
    per=[];pool=[];blocking=[];phase=[];stationarity=[]
    base=plan['block_samples']
    for c in chains:
        rr,dd=rows_for([c],base,'chain',output,f'seed{c.seed}');per+=rr
    for group in plan['groups']:
        cc=[c for c in chains if c.metadata['g']==group['g'] and c.metadata['intervals']==group['intervals']];assert len(cc)==2
        assert cc[0].metadata['fermions']==cc[1].metadata['fermions']
        tag=f"g{group['g']:.1f}_m{group['intervals']}"
        for block in plan['block_sensitivity']:
            rr,dd=rows_for(cc,block,'pooled',output,tag);blocking+=rr;phase.append(dd)
            if block==base:pool+=rr
        for c in cc:
            cut=len(c.counts)//2
            halves=[replace(c,sums=c.sums[:cut],counts=c.counts[:cut]),replace(c,sums=c.sums[cut:],counts=c.counts[cut:])]
            try:ea,eb=[estimates_json(joint_estimate([h],transform=observables)) for h in halves]
            except ValueError as exc:
                if 'unresolved phase denominator' not in str(exc):raise
                continue
            for name in ea:
                if name.endswith('_magnitude'):continue
                for component,index in [('real',0),('imag',1)]:
                    delta=eb[name]['value'][index]-ea[name]['value'][index];se=math.hypot(ea[name]['stderr_'+component],eb[name]['stderr_'+component])
                    stationarity.append(dict(g=group['g'],tau_total=group['tau_total'],seed=c.seed,observable=name,component=component,second_minus_first=delta,SE_screen=se,z_screen=delta/se if se>1e-14 else None,caveat='screen only; cross-half covariance unestimated; multiple channels; roundoff-scale SE excluded'))
    phaseper=[dict(g=c.metadata['g'],**chain_row(c)) for c in chains]
    save_csv(output/'phase_per_chain.csv',phaseper);save_csv(output/'phase_pooled.csv',phase)
    save_csv(output/'observables_per_chain.csv',per);save_csv(output/'observables_pooled.csv',pool)
    save_csv(output/'fermions_per_chain.csv',[r for r in per if r['observable'].startswith('F_')]);save_csv(output/'fermions_pooled.csv',[r for r in pool if r['observable'].startswith('F_')])
    save_csv(output/'blocking_sensitivity.csv',blocking)
    if stationarity:save_csv(output/'stationarity.csv',stationarity)
    (output/'fermion_paths.json').write_text(json.dumps(chains[0].metadata['fermions'],indent=2)+'\n')
    flags=[r for r in stationarity if r['z_screen'] is not None and abs(r['z_screen'])>=3]
    summary={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest':manifest,'phase_per_chain':phaseper,'phase_pooled':phase,'stationarity_flags':flags,'g0_max_absolute_error':max((r['analytic_abs_difference'] for r in pool if r['analytic_abs_difference'] is not None),default=None)}
    (output/'summary.json').write_text(json.dumps(jsonable(summary),indent=2,allow_nan=False)+'\n')
    lines=['# L=4 coupling scan with complex matter-Majorana correlations','',f"Completed {len(chains)} independent chains. Warmup {plan['warmup']}, measurements {plan['sweeps']} per chain; stride 1.",'','|g|tau_total|seeds|Re phase|Im phase|abs phase|SE|SNR|negative|complex|acceptance|','|---|---|---|---|---|---|---|---|---|---|---|']
    for p in phase:
        if p['block_sweeps']==base:lines.append('|'+ '|'.join(str(p[k]) for k in ['g','tau_total','seeds','Re','Im','magnitude','SE','SNR','negative_count','complex_count','acceptance'])+'|')
    lines+=['','Fermionic estimates: fermions_per_chain.csv and fermions_pooled.csv. Full complex Re/Im, component SE, covariance, modulus and modulus SE. Raw joint phase/numerator blocks remain in each run directory. NPZ files retain full covariance among reported channels. VBS/phonon observables are in observables_*.csv.','', 'C_ij=delta_ij-G_ij, translational mean over 16 source cells at the center equal-time slice. Gauge-string channels multiply the specified directed-u path; paths and PBC conventions are in fermion_paths.json. These are matter-Majorana, not spin or electron correlations. Different paths may differ; no radial averaging or long-distance exponent fit is implied.','',f"Independent analytic g=0 Bloch reference: largest absolute difference = {summary['g0_max_absolute_error']}.",'', 'Blocking 25/50/100/200 is a sensitivity check, not an independence proof. First/second-half screens are not corrected for multiple channels and do not estimate cross-half covariance. A stable phase alone does not establish equilibrium, projection-time, time-step or size convergence.', '', '未观测到负号（若 negative=0）不等于严格证明无符号。浮点虚部与巨大 SNR 不能排除罕见负号。An unresolved denominator yields only phase diagnostics, not physical ratios. No numerical failure is converted to rejection.','',f"Stationarity screens with abs(z)>=3: {len(flags)}; see stationarity.csv. Exact parameters, binary/source hashes, elapsed times and return codes are in run manifest/status/completion files."]
    (output/'FERMION_GSCAN_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('REPORT',output/'FERMION_GSCAN_REPORT.md',flush=True)
if __name__=='__main__':main()
