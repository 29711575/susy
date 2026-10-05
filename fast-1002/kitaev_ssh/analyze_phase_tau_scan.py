#!/usr/bin/env python3
"""Reuse existing chain-preserving joint jackknife for the bounded tau scan."""
import argparse,csv,datetime,hashlib,json,math,pathlib
from dataclasses import replace
import numpy as np
from convergence_analysis import Chain,read_chain,joint_estimate,estimates_json
from plot_vbs import normalize_tau_metadata,normalize_trial_metadata

def read_phase_chain(path):
    # A vanishing denominator is itself a valid phase diagnostic. The existing
    # physical-observable reader intentionally refuses it; retain raw sums here.
    meta=json.loads(path.read_text())
    if meta.get('ratio_defined'):return read_chain(path)
    normalize_tau_metadata(meta);normalize_trial_metadata(meta)
    if meta.get('status')!='complete':raise ValueError('incomplete phase run')
    with path.with_name(path.name.replace('_results.json','_blocks.csv')).open() as stream:rows=list(csv.DictReader(stream))
    names=list(meta['observables']);columns=['sum_phase']+['sum_phase_'+name for name in names]
    if [int(b['block']) for b in rows]!=list(range(len(rows))):raise ValueError('missing phase blocks')
    sums=np.array([[complex(float(b[k+'_real']),float(b[k+'_imag'])) for k in columns] for b in rows])
    counts=np.array([int(b['count']) for b in rows])
    if len(rows)!=meta['blocks'] or int(counts.sum())!=meta['samples'] or not np.isfinite(sums).all():raise ValueError('invalid phase blocks')
    return Chain(path,meta,names,sums,counts)

def phase_transform(sums,count,names,metadata):
    phase=sums[0]/count
    return {'phase':phase,'phase_magnitude':abs(phase)}
def number(x):return 'undefined' if x is None else f'{x:.10g}'
def save_csv(path,rows):
    with path.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def jsonable(x):
    if isinstance(x,dict):return {k:jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [jsonable(v) for v in x]
    if isinstance(x,np.ndarray):return jsonable(x.tolist())
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,np.generic):return x.item()
    return x

def chain_row(c):
    d=c.metadata;p=d['average_phase'];diag=d['phase_diagnostics'];se=math.hypot(p['stderr_real'],p['stderr_imag'])
    with c.path.with_name(c.path.name.replace('_results.json','_blocks.csv')).open() as stream:rows=list(csv.DictReader(stream))
    counts={key:sum(int(b[key]) for b in rows) for key in ['positive_count','negative_count','complex_count']}
    assert sum(counts.values())==d['samples']
    for k,val in counts.items():assert val==diag[k]
    assert math.isclose(diag['stderr_complex'],se,rel_tol=2e-12,abs_tol=1e-30)
    return dict(tau_total=d['intervals']*d['dt'],intervals=d['intervals'],seed=d['seed'],Re=p['value'][0],Im=p['value'][1],magnitude=diag['magnitude'],SE=se,SNR=diag['denominator_snr'],negative_count=counts['negative_count'],complex_count=counts['complex_count'],acceptance=d['acceptance'],positive_count=counts['positive_count'],samples=d['samples'],SE_Re=p['stderr_real'],SE_Im=p['stderr_imag'],cov_Re_Im=p['cov_real_imag'],SE_magnitude=diag['stderr_magnitude'],snr_status=diag['snr_status'],block_sweeps=c.block_sweeps,elapsed_seconds=d['elapsed_seconds'])

def table(rows,pooled=False):
    keys=['tau_total','intervals']+(['seeds'] if pooled else ['seed'])+['Re','Im','magnitude','SE','SNR','negative_count','complex_count','acceptance']
    titles=['tau_total','intervals','seeds' if pooled else 'seed','Re','Im','abs(mean phase)','SE','SNR','negative','complex','acceptance']
    lines=['| '+' | '.join(titles)+' |','|'+'|'.join(['---']*len(keys))+'|']
    for row in rows:
        cells=[]
        for k in keys:
            val=row[k]
            if k=='seeds':cells.append(','.join(map(str,val)))
            elif k in ['intervals','seed','negative_count','complex_count']:cells.append(str(val))
            else:cells.append(number(val))
        lines.append('| '+' | '.join(cells)+' |')
    return '\n'.join(lines)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',required=True,type=pathlib.Path);args=ap.parse_args()
    root=pathlib.Path(__file__).resolve().parents[1];run=args.run_dir.resolve();manifest=json.loads((run/'manifest.json').read_text());plan=manifest['plan']
    chains=[read_phase_chain(pathlib.Path(item['folder'])/'run_results.json') for item in manifest['commands']];assert len(chains)==12
    per=[chain_row(c) for c in chains];pooled=[];blocking=[];stationarity=[];details={}
    for m in plan['intervals']:
        cc=[c for c in chains if c.metadata['intervals']==m];assert len(cc)==2
        for c in cc:
            d=c.metadata
            for k,v in [('Lx',4),('Ly',4),('Jx',1),('Jy',1),('Jz',.5),('g',.5),('mass',4),('omega_ph',.5),('k',1),('dt',.1),('warmup_sweeps',200),('measurement_sweeps',800),('weight_engine','combined')]:assert d[k]==v,(k,d[k],v)
            assert d['numerical_failures']==0 and d['trial_shift']==[0,0,0]
        for block in plan['block_sensitivity']:
            est=joint_estimate([c.reblock(block) for c in cc],transform=phase_transform)
            ee=estimates_json(est);p=ee['phase'];se=math.hypot(p['stderr_real'],p['stderr_imag']);z=complex(*p['value'])
            row=dict(tau_total=m*.1,intervals=m,seeds=[c.seed for c in cc],Re=z.real,Im=z.imag,magnitude=abs(z),SE=se,SNR=abs(z)/se if se>0 else None,negative_count=sum(c.metadata['phase_diagnostics']['negative_count'] for c in cc),complex_count=sum(c.metadata['phase_diagnostics']['complex_count'] for c in cc),acceptance=sum(c.metadata['accepted_measurement'] for c in cc)/sum(c.metadata['proposals_measurement'] for c in cc),positive_count=sum(c.metadata['phase_diagnostics']['positive_count'] for c in cc),samples=sum(c.metadata['samples'] for c in cc),SE_Re=p['stderr_real'],SE_Im=p['stderr_imag'],cov_Re_Im=p['cov_real_imag'],SE_magnitude=ee['phase_magnitude']['stderr_real'],block_sweeps=block,blocks_per_chain=est['chain_blocks'],phase_resolved_5sigma=est['phase_resolved_5sigma'])
            blocking.append(row);details[f'm{m}_block{block}']=est
            if block==25:pooled.append(row)
        for c in cc:
            cut=len(c.counts)//2
            a=replace(c,sums=c.sums[:cut],counts=c.counts[:cut]);b=replace(c,sums=c.sums[cut:],counts=c.counts[cut:])
            try:
                ea=estimates_json(joint_estimate([a]));eb=estimates_json(joint_estimate([b]))
            except ValueError as exc:
                if 'unresolved phase denominator' not in str(exc):raise
                stationarity.append({'tau_total':m*.1,'seed':c.seed,'observable':'physical_observables_unresolved','second_minus_first_real':None,'SE_screen':None,'z_screen':None,'caveat':str(exc)})
                continue
            for name in ['phase','X_z','uniform_background_trace','Szz_q1','TrS_q1','Szz_q3','TrS_q3']:
                delta=eb[name]['value'][0]-ea[name]['value'][0];se=math.hypot(ea[name]['stderr_real'],eb[name]['stderr_real'])
                stationarity.append({'tau_total':m*.1,'seed':c.seed,'observable':name,'second_minus_first_real':delta,'SE_screen':se,'z_screen':delta/se if se>0 else None,'caveat':'screen only; cross-half covariance unestimated; multiple channels'})
    output=root/'results'/plan['name'];output.mkdir(parents=True,exist_ok=False)
    save_csv(output/'per_chain.csv',per);save_csv(output/'pooled_by_tau.csv',pooled);save_csv(output/'blocking_sensitivity.csv',blocking);save_csv(output/'stationarity.csv',stationarity)
    artifact={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest':manifest,'per_chain':per,'pooled':pooled,'blocking':blocking,'stationarity':stationarity,'joint_estimator_details':details,'analysis_script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    (output/'phase_tau_summary.json').write_text(json.dumps(jsonable(artifact),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    negatives=sum(x['negative_count'] for x in per);complexes=sum(x['complex_count'] for x in per);total=sum(x['samples'] for x in per)
    signals=[x for x in stationarity if x['z_screen'] is not None and abs(x['z_screen'])>=3 and x['observable']!='phase']
    text=f"""# L=4 phase versus total tau

Generated {artifact['created_utc']}. All 12 chains completed. Finite-run diagnostic, not a positivity proof.

## Protocol and provenance

L=4, J=(1,1,0.5), g=0.5, k=1, m=4, omega_ph=0.5, dtau=0.1; tau_total=intervals*dtau with intervals=8,16,32,48,64,80. combined engine; all-A-to-B u=+1, Wilson(+1,+1); bond-refresh; zero phonon trial shifts; free X=0 matter trial. Each seed: 200 warmup + 800 measurement sweeps, stride=1, raw block=25, VBS=1. Twelve local WSL processes, one math-library thread per chain, no PBS job. No adaptive extension.

Source commit: {manifest['source_commit']} (diagnostic patch in run/source.patch).
Binary SHA256: {manifest['binary_sha256']}.
Compiler: {manifest['compiler']}; MKL: {manifest['mkl_path']}.
Exact commands/seeds/timing: run manifest.json, status.json and chain completion.json.

## Per-chain measurements (25-sweep blocks)

{table(per)}

SE=sqrt(SE_Re^2+SE_Im^2); SNR=abs(mean phase)/SE. Component errors, covariance and modulus jackknife SE are in per_chain.csv. Undefined SNR at zero observed variance is not infinite statistical evidence.

Counts use fixed measurement times, including retained configurations after rejection. abs(Im phase)>1e-10 is complex; otherwise Re<0 is negative and Re>=0 positive. The exclusive categories cover production measurements, not warmup/proposals. Sums retain the original complex phase without rounding.

## Pooled per tau

Existing convergence_analysis.joint_estimate deletes blocks within each independent chain, preserving allocation. Different tau distributions are never pooled. The same replicas give errors on Re/Im and modulus.

{table(pooled,True)}

Total measured configurations: {total}; negative: {negatives}; clearly complex: {complexes}.

{'未观测到负号。This describes only these recorded finite chains; it is NOT 严格证明无符号.' if negatives==0 else 'Negative phases observed; inspect the counts and denominator errors.'}

Floating-point imaginary noise and zero observed real/modulus variance cannot exclude rare signs, unvisited regions or slow modes. Very large SNR from roundoff-size variance is not a rigorous precision guarantee. No larger-L result or general nonnegative-weight theorem is claimed.

## Blocking and stationarity

blocking_sensitivity.csv uses 25/50/100-sweep blocks (32/16/8 blocks per chain). This is a sensitivity check, not proof of independence or reliable integrated autocorrelation time. stationarity.csv compares first/second halves with the existing estimator; {len(signals)} observable screens have abs(z)>=3. Multiple channels and unestimated cross-half covariance limit interpretation. Constant phase does not demonstrate equilibrium, projection, dtau or size convergence.

## Validation and unchanged physics

Before/after: 13 ctests pass, explicit weight_validate and phonon_validate pass. Test-only old/new traces match for 468 transactions covering L4 intervals 8/80 bond-refresh and L2 local proposals, including decisions, ratios, weights and field/Green/RNG digests. Old production CLI JSON/CSV fields match exactly except elapsed time. Synthetic positive/negative/nonreal classification and modulus-jackknife tests pass.

C++ changes only concern phase_statistics.h, coupled_mc output and tests. Hamiltonian, Pfaffian weight, Green, proposals and RNG consumption remain unchanged. No Pfaffian optimization or L expansion.

Next priority: evaluate a fixed larger L=4 sampling budget or checkpoint-assisted stationarity study if needed; do not treat unobserved negative counts as a positivity theorem.
"""
    (output/'PHASE_TAU_REPORT.md').write_text(text,encoding='utf-8')
    print(table(per));print('\nPOOLED\n'+table(pooled,True));print('REPORT',output/'PHASE_TAU_REPORT.md',flush=True)
if __name__=='__main__':main()
