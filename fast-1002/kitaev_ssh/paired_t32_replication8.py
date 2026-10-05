#!/usr/bin/env python3
"""Fixed 8-chain paired endpoint analysis; chain is the independent unit."""
import argparse, json, math
from pathlib import Path
import numpy as np
from convergence_analysis import read_chain, checked_group, fixed_observables

# 0.975 quantile, Student t with 7 degrees of freedom.
T_CRITICAL_975 = {3: 3.1824463052837046, 7: 2.3646242510102993}

def contrasts(chains):
    checked_group(chains)
    rows=[]
    for c in chains:
        rb=c.reblock(200)
        assert len(rb.counts)==40 and int(rb.counts.sum())==8000
        vals=[]
        for sl in (slice(0,10),slice(30,40)):
            sums=rb.sums[sl].sum(axis=0); count=int(rb.counts[sl].sum())
            vals.append(fixed_observables(sums,count,c.names,c.metadata)['Szz_q1'])
        rows.append(dict(seed=c.seed,first=[vals[0].real,vals[0].imag],
                         fourth=[vals[1].real,vals[1].imag],difference=float((vals[1]-vals[0]).real)))
    x=np.array([r['difference'] for r in rows]); mean=float(x.mean()); sd=float(x.std(ddof=1));se=sd/math.sqrt(len(x))
    df=len(x)-1; critical=T_CRITICAL_975[df]
    return dict(chains=rows,n=len(x),degrees_of_freedom=df,mean_difference=mean,
                chain_difference_sd=sd,stderr=se,critical_t95=critical,
                ci95=[mean-critical*se,mean+critical*se],
                tolerance=.02,complete_interval_within_tolerance=bool(abs(mean)+critical*se<.02),
                limitation=f'Student-t assumes approximately normal independent chain differences; n={len(x)} remains small; blocks are not independent repeats')

def main():
    p=argparse.ArgumentParser();p.add_argument('--run-dir',type=Path,default=Path('build-kitaev-p01/t32_replication8'));a=p.parse_args()
    root=Path(__file__).resolve().parents[1];plan=json.loads(Path(__file__).with_name('t32_replication8_plan.json').read_text())
    new=[read_chain(a.run_dir/'results'/f'replication8_{s}_results.json') for s in 'abcdefgh']
    assert [c.seed for c in new]==plan['seeds']
    for c in new:
        assert c.metadata['status']=='complete' and c.metadata['numerical_failures']==0
        assert c.metadata['warmup_sweeps']==3200 and c.metadata['measurement_sweeps']==8000
        assert c.metadata['block_samples']==25 and c.metadata['trial_shift']==plan['trial_shift']
    old1=[read_chain(root/'build-kitaev-p01/shifted_t32_stationarity/results'/f't32_stationary_{s}_results.json') for s in 'abcd']
    old2=[read_chain(root/'build-kitaev-p01/t32_replication/results'/f'replication_{s}_results.json') for s in 'abcd']
    out=dict(plan=plan,new_primary=contrasts(new),historical_7101_7104=contrasts(old1),historical_7201_7204=contrasts(old2),historical_not_pooled=True)
    target=a.run_dir/'paired_replication8.json';target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('Saved',target)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    r=out['new_primary'];ax[0].plot([x['seed'] for x in r['chains']],[x['difference'] for x in r['chains']],'o');ax[0].axhline(0,color='grey',lw=.8)
    ax[0].set(xlabel='new seed',ylabel='Szz(q1): Q4-Q1',title='8 independent paired contrasts')
    for i,(label,z) in enumerate((('historical 7101-7104',out['historical_7101_7104']),('historical 7201-7204',out['historical_7201_7204']),('new 7301-7308',r))):
        critical=z['critical_t95'];ax[1].errorbar(z['mean_difference'],i,xerr=critical*z['stderr'],fmt='o',capsize=4)
    ax[1].axvspan(-.02,.02,color='grey',alpha=.15);ax[1].axvline(0,color='grey',lw=.8);ax[1].set(yticks=range(3),yticklabels=['historical 7101-7104','historical 7201-7204','new 7301-7308'],xlabel='mean paired difference; batch-specific 95% t interval')
    fig.suptitle('tau=3.2 shifted trial: new batch primary; historical batches not pooled');fig.savefig(a.run_dir/'paired_replication8.png',dpi=160,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':main()
