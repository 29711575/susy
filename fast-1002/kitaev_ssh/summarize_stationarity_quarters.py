#!/usr/bin/env python3
"""Existing tau=3.2 data only: chain-preserving quarter summary and multiplicity."""
import json
import math
from dataclasses import replace
from pathlib import Path

from convergence_analysis import read_chain, joint_estimate, STATIONARITY


def contrast(chains):
    # Match the original screen: separate first/last-quarter joint ratio JKs.
    # Their same-chain cross-quarter covariance is not estimated here.
    ends = []
    for sl in (slice(0, 10), slice(30, 40)):
        parts = []
        for chain in chains:
            rb = chain.reblock(200)
            assert len(rb.counts) == 40
            parts.append(replace(rb, sums=rb.sums[sl], counts=rb.counts[sl]))
        ends.append(joint_estimate(parts))
    i = ends[0]['names'].index('Szz_q1')
    d = float((ends[1]['values'][i]-ends[0]['values'][i]).real)
    se = math.sqrt(ends[0]['covariance'][i, i]+ends[1]['covariance'][i, i])
    return dict(seeds=[c.seed for c in chains], difference=d, stderr=se,
                ci95=[d-1.96*se, d+1.96*se], z=d/se,
                normal_two_sided_p=math.erfc(abs(d/se)/math.sqrt(2)))


def adjusted(rows):
    n = len(rows)
    order = sorted(range(n), key=lambda i: rows[i]['p'])
    output = [dict(r, bonferroni=min(1., n*r['p'])) for r in rows]
    running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., (n-rank)*rows[i]['p']))
        output[i]['holm'] = running
    return output


def main():
    root = Path(__file__).resolve().parents[1]
    run = root/'build-kitaev-p01/shifted_t32_stationarity'
    group = json.loads((run/'analysis/groups.json').read_text())['NEW']
    chains = [read_chain(run/'results'/f't32_stationary_{letter}_results.json')
              for letter in 'abcd']
    single = [contrast([c]) for c in chains]
    pooled = contrast(chains)
    leave_one_out = [dict(omitted=c.seed, **contrast([x for x in chains if x.seed != c.seed]))
                     for c in chains]
    families = {}
    for key in ('split_chain_diagnostics', 'quarter_first_last_diagnostics', 'between_chain_checks'):
        rows = []
        for r in group[key]:
            if r['observable'] not in STATIONARITY:
                continue
            z = r.get('z_diagnostic', r.get('z'))
            rows.append(dict(family=key, id=r.get('seed', r.get('comparison')),
                             observable=r['observable'], z=z,
                             p=math.erfc(abs(z)/math.sqrt(2)) if z is not None else 1.))
        families[key] = adjusted(rows)
    all_rows = adjusted([r for rows in families.values() for r in rows])
    output = dict(
        method='200-sweep blocks; pooled complex phase/numerator sums and stratified joint JK; '
               'same-chain first/last-quarter cross covariance omitted as in original screen; '
               'normal approximate p/intervals; leave-one-out is sensitivity, not data removal',
        individual=single, pooled=pooled, leave_one_out=leave_one_out,
        elementary_family_counts={key: len(rows) for key, rows in families.items()},
        family_adjustments=families, all_84_adjustments=all_rows,
        caveat='84-test union is a supplemental family choice, not a replacement for the original '
               'decision rule. The newly requested pooled/leave-one-out contrasts were not '
               'prespecified tests and are not added retrospectively to the original 84.')
    target = run/'analysis/quarter_summary_multiplicity.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k: output[k] for k in ('individual','pooled','leave_one_out','elementary_family_counts')}, indent=2))
    for key, rows in {**families, 'all84': all_rows}.items():
        print(key, 'min:', min(rows, key=lambda r: r['p']),
              'Holm<.05:', sum(r['holm']<.05 for r in rows))
    print('Saved', target)


if __name__ == '__main__':
    main()
