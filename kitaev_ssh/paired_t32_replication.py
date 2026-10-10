#!/usr/bin/env python3
"""Fixed four-chain replication; independent unit = paired chain contrast."""
import argparse
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
from convergence_analysis import read_chain, checked_group, fixed_observables, joint_estimate


def student3_cdf(x):
    # Integral of the df=3 Student density, avoiding an additional dependency.
    z = x/math.sqrt(3)
    return .5+(math.atan(z)+z/(1+z*z))/math.pi


def student3_critical():
    lo, hi = 0., 100.
    for _ in range(80):
        mid = (lo+hi)/2
        if student3_cdf(mid) < .975:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def summarize(chains):
    checked_group(chains)
    assert len(chains) == 4 and len({c.seed for c in chains}) == 4
    rows, endpoints, phase_resolved = [], [[], []], True
    for c in chains:
        rb = c.reblock(200)
        assert len(rb.counts) == 40 and sum(rb.counts) == 8000
        means = []
        for endpoint, sl in enumerate((slice(0,10), slice(30,40))):
            part = replace(rb, sums=rb.sums[sl], counts=rb.counts[sl])
            endpoints[endpoint].append(part)
            e = joint_estimate([part])
            phase_resolved &= e['phase_resolved_5sigma']
            value = fixed_observables(part.sums.sum(axis=0), int(part.counts.sum()),
                                      c.names, c.metadata)['Szz_q1']
            means.append(value)
        rows.append(dict(seed=c.seed, first_quarter=[means[0].real, means[0].imag],
                         last_quarter=[means[1].real, means[1].imag],
                         paired_difference=float((means[1]-means[0]).real)))
    first = np.array([r['first_quarter'][0] for r in rows])
    last = np.array([r['last_quarter'][0] for r in rows])
    differences = last-first
    estimate = float(np.mean(differences))
    sd = float(np.std(differences, ddof=1))
    se = sd/2
    paired_covariance = float(np.cov(first,last,ddof=1)[0,1])
    variance_identity = float((np.var(first,ddof=1)+np.var(last,ddof=1)-2*paired_covariance)/4)
    assert abs(variance_identity-se*se) < 1e-14
    critical = student3_critical()
    ci = [estimate-critical*se, estimate+critical*se]
    pooled = []
    for parts in endpoints:
        e = joint_estimate(parts)
        pooled.append(e['values'][e['names'].index('Szz_q1')])
    return dict(independent_replicates=4, degrees_of_freedom=3, chains=rows,
                mean_difference=estimate, chain_difference_sd=sd, stderr=se,
                critical_t95=critical, ci95=ci,
                two_sided_t_p=2*(1-student3_cdf(abs(estimate/se))) if se else None,
                covariance_of_first_last_chain_means=paired_covariance,
                variance_of_mean_from_paired_identity=variance_identity,
                pooled_complex_ratio_difference=[(pooled[1]-pooled[0]).real,(pooled[1]-pooled[0]).imag],
                average_chain_vs_pooled_difference=float(estimate-(pooled[1]-pooled[0]).real),
                all_quarter_phases_resolved=bool(phase_resolved),
                numerical_failures=sum(c.metadata['numerical_failures'] for c in chains),
                sample_normality_caveat='four independent chain contrasts; t3 assumes approximate normality; '
                    'no blocks treated as independent replicates and no assumption of zero within-chain endpoint covariance')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=Path('build-kitaev-p01/t32_replication'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    plan = json.loads(Path(__file__).with_name('t32_replication_plan.json').read_text())
    new_chains = [read_chain(args.run_dir/'results'/f'replication_{s}_results.json') for s in 'abcd']
    assert [c.seed for c in new_chains] == plan['seeds']
    for c in new_chains:
        assert c.metadata['status'] == 'complete'
        assert c.metadata['warmup_sweeps'] == 3200 and c.metadata['measurement_sweeps'] == 8000
        assert c.metadata['intervals'] == 16 and c.metadata['dt'] == .2
        assert c.metadata['trial_shift'] == plan['trial_shift']
        for field, expected in dict(Lx=2,Ly=2,Jx=1.,Jy=1.,Jz=.5,g=.3,mass=1.7,omega_ph=.8,
                                    measure_every=1,block_samples=25).items():
            assert c.metadata[field] == expected, field
    new = summarize(new_chains)
    old_chains = [read_chain(root/'build-kitaev-p01/shifted_t32_stationarity/results'/
                           f't32_stationary_{s}_results.json') for s in 'abcd']
    old = summarize(old_chains)
    tolerance = plan['primary']['practical_tolerance']
    new['positive_change_independently_supported'] = new['ci95'][0] > 0
    new['complete_interval_within_tolerance'] = new['ci95'][0] > -tolerance and new['ci95'][1] < tolerance
    new['ready_for_projection_comparison_by_primary_rule'] = bool(
        new['complete_interval_within_tolerance'] and new['all_quarter_phases_resolved'] and not new['numerical_failures'])
    output = dict(plan=plan, new_independent_replication=new, historical_only=old,
                  historical_and_new_not_pooled=True,
                  auxiliary_diagnostics='auxiliary/; no additional sigma-based vetoes')
    target = args.run_dir/'paired_replication.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    axes[0].plot([r['seed'] for r in new['chains']],
                 [r['paired_difference'] for r in new['chains']], 'o')
    axes[0].axhline(0,color='grey',linewidth=.8)
    axes[0].set(xticks=[r['seed'] for r in new['chains']],xlabel='new independent seed',
                ylabel='Szz(q1): Q4 minus Q1',title='Four paired chain contrasts')
    axes[1].axvspan(-tolerance,tolerance,color='grey',alpha=.15,label='prespecified +/-0.02')
    for i, (label, result) in enumerate((('old 7101-7104',old),('new 7201-7204',new))):
        axes[1].errorbar(result['mean_difference'],i,xerr=result['critical_t95']*result['stderr'],
                         fmt='o',capsize=4)
    axes[1].axvline(0,color='grey',linewidth=.8)
    axes[1].set(yticks=[0,1],yticklabels=['old: historical only','new: independent replication'],
                xlabel='mean paired difference; 95% t(3) interval')
    axes[1].legend(fontsize=8)
    fig.suptitle('tau=3.2, h=0.2, shifted trial fixed; 4 chains per batch, no pooling across batches')
    fig.savefig(args.run_dir/'paired_replication.png',dpi=160,bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(output,indent=2))
    print('Saved',target)


if __name__ == '__main__':
    main()
