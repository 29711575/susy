"""Autocorrelation of SAVED block means; not per-sweep sample ACF."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from convergence_analysis import read_chain,fixed_momentum_mapping
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/phase-tau-20260917/autocorrelation'
OUT.mkdir(exist_ok=True)
manifest=json.loads((ROOT/'runs/phase-tau-20260917/manifest.json').read_text())
names=['m2_q1','m2_q3','X_z','majorana_energy']
labels=[r'$S_{zz}(q_1)/N_c$',r'$S_{zz}(q_3)/N_c$',r'$X_z$', 'Matter energy']
rows=[]; traces=[]; curves={}
for command in manifest['commands']:
    c=read_chain(Path(command['folder'])/'run_results.json')
    mp=fixed_momentum_mapping(c.metadata)
    nc=c.metadata['Lx']*c.metadata['Ly']
    sources=[f"S_q{mp['q1']['index']}_zz",f"S_q{mp['q3']['index']}_zz",'X_z','majorana_energy']
    assert len(c.counts)==32 and np.all(c.counts==25)
    tau=round(c.metadata['intervals']*c.metadata['dt'],12)
    for obs,source in zip(names,sources):
        x=(c.sums[:,c.names.index(source)+1]/c.sums[:,0]).real
        if obs.startswith('m2'): x=x/nc
        z=x-x.mean(); ss=z@z
        if ss<=0: raise ValueError('Constant series: ACF undefined')
        # Standard biased sample ACF: common normalization at every lag.
        ac=np.array([z[:len(z)-lag]@z[lag:]/ss for lag in range(9)])
        curves[tau,c.seed,obs]=ac
        for i,value in enumerate(x):
            traces.append(dict(tau_total=tau,seed=c.seed,observable=obs,block=i,
                               measurement_sweep_start=25*i+1,measurement_sweep_end=25*(i+1),value=value))
        for lag,rho in enumerate(ac):
            rows.append(dict(tau_total=tau,seed=c.seed,observable=obs,lag_blocks=lag,
                             lag_sweeps=lag*25,pairs=32-lag,acf=rho))
def save(name,data):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
save('acf_per_chain.csv',rows);save('block_time_series.csv',traces)
taus=sorted({k[0] for k in curves}); pooled=[]
for tau in taus:
    for obs in names:
        aa=np.array([v for (t,s,o),v in curves.items() if t==tau and o==obs])
        assert aa.shape==(2,9)
        for lag in range(9):
            pooled.append(dict(tau_total=tau,observable=obs,lag_blocks=lag,lag_sweeps=25*lag,
                               mean_acf=aa[:,lag].mean(),seed_min=aa[:,lag].min(),seed_max=aa[:,lag].max()))
save('acf_mean_of_seeds.csv',pooled)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,obs,label in zip(axes.flat,names,labels):
    for tau in taus:
        rr=[r for r in pooled if r['tau_total']==tau and r['observable']==obs]
        ax.plot([r['lag_sweeps'] for r in rr],[r['mean_acf'] for r in rr],'o-',ms=3,label=f'tau={tau:g}')
    ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.2)
    ax.set(xlabel='Block separation (sweeps; 25 sweeps/block)',ylabel='Block-mean ACF',title=label,ylim=(-.65,1.05))
    ax.legend(ncol=3,fontsize=8)
fig.suptitle('Saved block means: ACF averaged over two independently centered seeds\n32 blocks/seed; lines guide the eye; not per-sweep ACF')
fig.savefig(OUT/'acf_vs_interval.png',dpi=200);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for ax,tau in zip(axes.flat,taus):
    seeds=sorted({s for t,s,o in curves if t==tau})
    for obs,color,lab in zip(names[:2],['#0072B2','#D55E00'],['q1','q3']):
        for seed,ls in zip(seeds,['-','--']):
            ax.plot(np.arange(9)*25,curves[tau,seed,obs],ls,color=color,marker='o',ms=3,label=f'{lab}, seed {seed}')
    ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.2)
    ax.set(title=f'tau={tau:g}',xlabel='Block separation (sweeps)',ylabel='Block-mean ACF',ylim=(-.65,1.05))
    ax.legend(fontsize=8)
fig.suptitle('VBS candidate channels: each seed separately (25-sweep block means)')
fig.savefig(OUT/'order_acf_per_seed.png',dpi=200);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,obs,label in zip(axes.flat,names,labels):
    for ti,tau in enumerate(taus):
        seeds=sorted({s for t,s,o in curves if t==tau})
        for si,seed in enumerate(seeds):
            rr=[r for r in traces if r['tau_total']==tau and r['seed']==seed and r['observable']==obs]
            ax.plot([r['measurement_sweep_end']-12 for r in rr],[r['value'] for r in rr],
                    color=f'C{ti}',ls='-' if si==0 else '--',alpha=.8,lw=1,label=f'tau={tau:g}' if si==0 else None)
    ax.set(title=label,xlabel='Measurement sweep (block center)',ylabel='Block mean');ax.grid(alpha=.2)
    ax.legend(fontsize=8,ncol=3)
fig.suptitle('Saved block time series; solid/dashed = independent seeds')
fig.savefig(OUT/'block_time_series.png',dpi=200);plt.close(fig)
lines=['# Saved-block autocorrelation diagnostic','',
'Each chain: 800 measurements saved as 32 consecutive 25-sweep blocks. No per-sweep series is available.',
'We compute Re[sum(phase*O)/sum(phase)] for each block; complex sums are retained before division.',
'q1 and q3 are the fixed momentum mapping from convergence_analysis; Nc=16.',
'ACF(k)=sum_{t=0}^{N-k-1}[(x_t-mean)(x_{t+k}-mean)] / sum_t[(x_t-mean)^2].',
'Each seed is centered and normalized separately; reported mean ACF is the arithmetic mean of two seed ACFs, never a concatenated chain.',
'Block separation is NOT thinning of the original measurements. Within-block correlations at 1--24 sweeps are not recoverable.',
'Only lags 0--8 (0--200 sweeps) are shown, with 32--24 pairs per seed. Finite-chain centering biases ACF and negative estimates need not be physical anticorrelation.',
'No confidence bands, integrated autocorrelation time or effective sample size are claimed: 32 block means per seed cannot resolve slow tails reliably. Seed ranges in CSV are NOT confidence intervals.',
'Stationarity is not established; drift can also produce autocorrelation. Inspect block_time_series.png and existing stationarity.csv.',
'Phase real part is constant, so its normalized ACF is undefined; imaginary roundoff is not a meaningful mixing diagnostic.','',
'## Mean ACF at one block (25 sweeps)','',
'| tau | m2 q1 | m2 q3 | X_z | matter energy |','|---|---|---|---|---|']
for tau in taus:
    vals=[next(r['mean_acf'] for r in pooled if r['tau_total']==tau and r['observable']==obs and r['lag_blocks']==1) for obs in names]
    lines.append('| '+f'{tau:g} | '+' | '.join(f'{v:.4f}' for v in vals)+' |')
(OUT/'README.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[-8:]))
print(OUT)

