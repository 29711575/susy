#!/usr/bin/env python3
"""Fixed-mode zero-field amplitude scan. Raw common-phase sums -> joint ratios.

Point errors reuse the chain-preserving joint jackknife. Scan replicas resample
whole reblocked rows WITHIN each independent chain, with fixed chain allocation.
Local degree <=2 Binder fits are repeated on every joint bootstrap replica.
Neither Gamma selection nor data-dependent momentum maxima are used.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from convergence_analysis import read_chain, checked_group, joint_estimate, estimates_json
from plot_vbs import write_csv


def mode(metadata, n1, n2):
    lx, ly = metadata['Lx'], metadata['Ly']
    if not (0 <= n1 < lx and 0 <= n2 < ly) or (n1 == 0 and n2 == 0):
        raise ValueError('select one explicit nonzero allowed momentum')
    grid = metadata.get('vbs', {}).get('momenta', [])
    if len(grid) != lx*ly or metadata['vbs'].get('Nc') != lx*ly:
        raise ValueError('missing full allowed momentum grid')
    matches = [m for m in grid if m['n'] == [n1, n2]]
    if len(matches) != 1:
        raise ValueError('momentum absent or duplicated')
    return matches[0]


def amplitude_chain(path, n1, n2, bond):
    chain = read_chain(path)
    meta = chain.metadata
    amp = meta.get('vbs_amplitude', {})
    if amp.get('schema') != 'vbs-amplitude-v1':
        raise ValueError('missing true amplitude moments; old S/chi blocks cannot reconstruct A or A4')
    if meta.get('vbs_field', {}).get('h', 0) != 0 or amp.get('source_epsilon') != 0:
        raise ValueError('this analysis requires zero physical source')
    if meta['g'] < 0:
        raise ValueError('lambda path requires g=+sqrt(k*lambda)')
    q = mode(meta, n1, n2)['index']
    stems = ['AmpA', 'AmpA2', 'AmpA4']
    if amp['response_enabled']:
        stems += ['AmpIA', 'AmpAIA']
    names = [f'{s}_q{q}_{bond}' for s in stems]
    if any(name not in chain.names for name in names):
        raise ValueError('incomplete amplitude moment schema')
    indices = [0]+[1+chain.names.index(name) for name in names]
    return replace(chain, names=names, sums=chain.sums[:, indices])


def amplitude_observables(sums, count, names, metadata):
    if abs(sums[0]) <= 1e-12*count:
        raise ValueError('numerically unresolved phase denominator')
    r = sums[1:]/sums[0]
    a, a2, a4 = r[:3]
    if abs(a2) <= 1e-14*np.sqrt(abs(a4)) or a2 == 0:
        raise ValueError('undefined Binder denominator')
    nc = metadata['Lx']*metadata['Ly']
    out = dict(phase=sums[0]/count, A=a, A2=a2, A4=a4,
               VA=nc*(a2-a*a), B4=a4/(a2*a2))
    if len(r) == 5:
        out.update(IA=r[3], chiA=nc*(r[4]-a*r[3]))
    return out


def scan_groups(chains, n1, n2):
    if not chains or len({str(c.path.resolve()) for c in chains}) != len(chains):
        raise ValueError('empty or duplicated input chain files')
    # Different lambda points must be statistically independent for this
    # bootstrap; common-random-number scans need a different replica scheme.
    if len({c.seed for c in chains}) != len(chains):
        raise ValueError('all input chains need distinct independent seeds, including across lambda')
    fixed = ('Lx','Ly','shear','geometry_id','Jx','Jy','Jz','mass','omega_ph',
             'intervals','dt','trial_shift','scope','spatial_boundary','wilson_loops','time_contour')
    reference = chains[0].metadata
    for c in chains:
        for key in fixed:
            if c.metadata[key] != reference[key]:
                raise ValueError(f'lambda scan changed fixed parameter: {key}')
        if c.names != chains[0].names:
            raise ValueError('amplitude response schema changed within scan')
        if mode(c.metadata, n1, n2) != mode(reference, n1, n2):
            raise ValueError('fixed momentum changed within scan')
    grouped = {}
    for c in chains:
        k = c.metadata['mass']*c.metadata['omega_ph']**2
        lam = c.metadata['g']**2/k
        grouped.setdefault(lam, []).append(c)
    for group in grouped.values():
        checked_group(group)
    return dict(sorted(grouped.items()))


def fit_slopes(x, binder, window, degree=2):
    """Signed local slopes on interior centers only; no extrapolated endpoints."""
    x = np.asarray(x, dtype=float)
    binder = np.asarray(binder)
    if window < 3 or window % 2 != 1 or degree not in (1, 2) or degree >= window:
        raise ValueError('odd window >=3 and fit degree 1 or 2 required')
    if len(x) < window or np.any(np.diff(x) <= 0):
        raise ValueError('need enough strictly increasing lambda points')
    half = window//2
    centers = np.arange(half, len(x)-half)
    slopes = []
    for i in centers:
        delta = x[i-half:i+half+1]-x[i]
        scale = max(abs(delta))
        design = np.vander(delta/scale, degree+1, increasing=True)
        if np.linalg.cond(design) > 1e8:
            raise ValueError('ill-conditioned local lambda fit')
        weights = np.linalg.pinv(design)[1]/scale
        slopes.append(binder[..., i-half:i+half+1] @ weights)
    return centers, np.stack(slopes, axis=-1)


def bootstrap_curves(groups, replicas, rng):
    first = next(iter(groups.values()))
    labels = list(amplitude_observables(sum(c.sums.sum(axis=0) for c in first),
                                       sum(int(c.counts.sum()) for c in first),
                                       first[0].names,first[0].metadata))
    curves = []
    failures = 0
    for group in groups.values():
        chain = group[0]
        count = sum(int(c.counts.sum()) for c in group)
        values = []
        for _ in range(replicas):
            sums = sum(c.sums[rng.integers(0, len(c.counts), len(c.counts))].sum(axis=0) for c in group)
            try:
                result = amplitude_observables(sums, count, chain.names, chain.metadata)
                labels = list(result)
                values.append(list(result.values()))
            except ValueError:
                failures += 1
                # Keep failed replicas in the array; do not condition the
                # reported quantiles on successful phase denominators.
                width = 8 if len(chain.names) == 5 else 6
                values.append([complex(np.nan, np.nan)]*width)
        curves.append(values)
    return labels, np.transpose(np.asarray(curves), (1, 0, 2)), failures


def quality(estimate):
    reasons = []
    if not estimate['phase_resolved_5sigma']:
        reasons.append('phase unresolved at 5 sigma')
    if min(estimate['chain_blocks']) < 16:
        reasons.append('fewer than 16 blocks in at least one chain')
    n = len(estimate['names'])
    j = estimate['names'].index('A2')
    err = np.sqrt(estimate['covariance'][j,j]+estimate['covariance'][n+j,n+j])
    if abs(estimate['values'][j]) <= 5*err:
        reasons.append('Binder second moment unresolved at 5 sigma')
    for j, name in enumerate(estimate['names']):
        if name == 'phase':
            continue
        err = np.sqrt(estimate['covariance'][n+j,n+j])
        if abs(estimate['values'][j].imag) > 5*err+1e-12*max(1, abs(estimate['values'][j])):
            reasons.append(f'{name} imaginary residual exceeds 5 sigma')
    return reasons


def peak_scan(x, binder, replicas, windows, degree, qualified):
    rows, peaks = [], []
    if not np.isfinite(replicas).all():
        return rows, [dict(status='bootstrap phase/Binder failures; slope intervals withheld')]
    for window in windows:
        if len(x) < window:
            continue
        centers, signed = fit_slopes(x, binder, window, degree)
        _, boot_signed = fit_slopes(x, replicas, window, degree)
        magnitude = np.abs(signed)
        boot = np.abs(boot_signed)
        intervals = np.quantile(boot, [.025,.975], axis=0)
        for j, i in enumerate(centers):
            rows.append(dict(lambda_value=x[i], window=window, degree=degree,
                             signed_slope=signed[j], abs_slope=magnitude[j],
                             abs_slope_low=intervals[0,j], abs_slope_high=intervals[1,j]))
        if len(centers) < 3:
            peaks.append(dict(window=window,degree=degree,
                              status='fewer than three interior fit centers; peak location undefined'))
            continue
        k = int(np.argmax(magnitude))
        indices = centers[np.argmax(boot, axis=1)]
        boundary = k in (0, len(centers)-1)
        peaks.append(dict(window=window, degree=degree, lambda_grid_peak=x[centers[k]],
                          lambda_peak_95_interval=np.quantile(x[indices],[.025,.975]).tolist(),
                          lambda_grid_cell=[float((x[centers[k]-1]+x[centers[k]])/2),
                                            float((x[centers[k]]+x[centers[k]+1])/2)],
                          bootstrap_fit_boundary_fraction=float(np.mean((indices==centers[0])|(indices==centers[-1]))),
                          status='truncated fit domain; extend lambda scan' if boundary else
                                 ('candidate only; block independence and projection/size convergence unverified' if qualified else
                                  'insufficient point statistics; candidate unqualified')))
    return rows, peaks


def run_analysis(chains, n1, n2, windows, degree, replicas, seed, block_sweeps=None):
    groups = scan_groups(chains, n1, n2)
    if block_sweeps is None:
        block_sweeps = int(np.lcm.reduce([c.block_sweeps for c in chains]))
    groups = {lam:[c.reblock(block_sweeps) for c in group] for lam,group in groups.items()}
    estimates = [joint_estimate(group, amplitude_observables) for group in groups.values()]
    points = []
    for (lam, group), e in zip(groups.items(), estimates):
        points.append(dict(lambda_value=lam, samples=e['samples'], seeds=e['seeds'],
                           chain_blocks=e['chain_blocks'], estimates=estimates_json(e),
                           joint_covariance_real_then_imag=e['covariance'].tolist(),
                           average_phase=[e['average_phase'].real,e['average_phase'].imag],
                           quality_flags=quality(e), files=[str(c.path.resolve()) for c in group]))
    labels, curves, failures = bootstrap_curves(groups, replicas, np.random.default_rng(seed))
    x = np.array(list(groups))
    binder = np.array([p['estimates']['B4']['value'][0] for p in points])
    bi = labels.index('B4')
    slopes, peaks = peak_scan(x, binder, curves[:,:,bi].real, windows, degree,
                             all(not p['quality_flags'] for p in points))
    amplitude_peaks = []
    if len(x) >= 3 and np.isfinite(curves).all():
        for name in ['VA']+(['chiA'] if 'chiA' in labels else []):
            j = labels.index(name)
            values = np.array([p['estimates'][name]['value'][0] for p in points])
            k = int(np.argmax(values))
            indices = np.argmax(curves[:,:,j].real,axis=1)
            amplitude_peaks.append(dict(observable=name,lambda_grid_peak=x[k],
                lambda_peak_95_interval=np.quantile(x[indices],[.025,.975]).tolist(),
                lambda_grid_cell=[float(x[0] if k==0 else (x[k-1]+x[k])/2),
                                  float(x[-1] if k==len(x)-1 else (x[k]+x[k+1])/2)],
                bootstrap_boundary_fraction=float(np.mean((indices==0)|(indices==len(x)-1))),
                status='scan-boundary maximum; extend lambda range' if k in (0,len(x)-1) else
                       ('candidate only; convergence unverified' if all(not p['quality_flags'] for p in points) else
                        'insufficient point statistics; candidate unqualified')))
    sensitivity = []
    for (lam, group) in groups.items():
        if all(len(c.counts)%2==0 and len(c.counts)>=4 for c in group):
            e = joint_estimate([c.reblock(2*block_sweeps) for c in group], amplitude_observables)
            sensitivity.append(dict(lambda_value=lam, block_sweeps=2*block_sweeps,
                                    estimates=estimates_json(e), quality_flags=quality(e)))
        else:
            sensitivity.append(dict(lambda_value=lam, status='doubling blocks unavailable without discarding samples'))
    # Check the error-bar dependence of slopes using the doubled within-chain
    # jackknife point variances. No independent moment error propagation.
    doubled_slope_errors = []
    if all('estimates' in p for p in sensitivity):
        bse = np.array([p['estimates']['B4']['stderr_real'] for p in sensitivity])
        for w in windows:
            if len(x) < w:
                continue
            centers, coefficients = fit_slopes(x, np.eye(len(x)), w, degree)
            se = np.sqrt(np.sum((coefficients*bse[:,None])**2, axis=0))
            doubled_slope_errors += [dict(lambda_value=x[i],window=w,signed_slope_stderr=se[j]) for j,i in enumerate(centers)]
    fixed_keys=('Lx','Ly','shear','geometry_id','Jx','Jy','Jz','mass','omega_ph','intervals','dt','tau',
                'trial_shift','scope','spatial_boundary','wilson_loops','time_contour')
    return dict(schema='vbs-amplitude-scan-v1', mode=mode(chains[0].metadata,n1,n2),
                fixed_parameters={k:chains[0].metadata[k] for k in fixed_keys if k in chains[0].metadata},
                bond=chains[0].names[0].rsplit('_',1)[1], block_sweeps=block_sweeps,
                bootstrap=dict(replicas=replicas,seed=seed,failed_point_replicas=failures,
                               scheme='joint rows resampled within independent chains; fixed chain sample allocation'),
                points=points,amplitude_peak_candidates=amplitude_peaks,binder_slopes=slopes,binder_peak_candidates=peaks,
                doubled_block_sensitivity=sensitivity,doubled_block_signed_slope_errors=doubled_slope_errors,
                interpretation='zero-field finite-dt finite-projection data; grid peaks are pseudocritical candidates; no lambda_c or SUSY certification')


def plots(result, output):
    import matplotlib.pyplot as plt
    points = result['points']
    names = ['VA','B4']+(['chiA'] if 'chiA' in points[0]['estimates'] else [])
    fig, axes = plt.subplots(1,len(names)+1,figsize=(4*(len(names)+1),3.5))
    x = [p['lambda_value'] for p in points]
    for ax,name in zip(axes,names):
        ax.errorbar(x,[p['estimates'][name]['value'][0] for p in points],
                    yerr=[p['estimates'][name]['stderr_real'] for p in points],fmt='o-')
        ax.set(xlabel=r'$\lambda$',ylabel=name)
    for w in sorted(set(p['window'] for p in result['binder_slopes'])):
        rows = [p for p in result['binder_slopes'] if p['window']==w]
        xx = [p['lambda_value'] for p in rows]
        axes[-1].plot(xx,[p['abs_slope'] for p in rows],'o-',label=f'window {w}')
        axes[-1].fill_between(xx,[p['abs_slope_low'] for p in rows],[p['abs_slope_high'] for p in rows],alpha=.15)
    axes[-1].set(xlabel=r'$\lambda$',ylabel=r'$|\partial_\lambda B_4|$')
    if result['binder_slopes']:
        axes[-1].legend()
    flagged = any(p['quality_flags'] for p in points)
    params=result['fixed_parameters']
    geometry=f"{params['Lx']}x{params['Ly']}, shear={params['shear']}, tau={params['intervals']*params['dt']:.4g}, dt={params['dt']:.4g}"
    fig.suptitle(('INSUFFICIENT STATISTICS; ' if flagged else '')+geometry+f"; fixed q={result['mode']['n']}, bond={result['bond']}; convergence unverified")
    fig.tight_layout();fig.savefig(output/'amplitude_scan.png',dpi=180);plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',nargs='+',type=Path)
    p.add_argument('--n1',required=True,type=int);p.add_argument('--n2',required=True,type=int)
    p.add_argument('--bond',choices=['x','y','z'],default='z')
    p.add_argument('--block-sweeps',type=int)
    p.add_argument('--windows',type=int,nargs='+',default=[3,5,7])
    p.add_argument('--degree',type=int,choices=[1,2],default=2)
    p.add_argument('--replicas',type=int,default=1000);p.add_argument('--seed',type=int,default=102019)
    p.add_argument('--output',required=True,type=Path)
    args = p.parse_args()
    if args.replicas < 100 or any(w<3 or w%2==0 for w in args.windows):
        p.error('at least 100 replicas and odd windows >=3 required')
    try:
        chains = [amplitude_chain(path,args.n1,args.n2,args.bond) for path in args.results]
        result = run_analysis(chains,args.n1,args.n2,args.windows,args.degree,args.replicas,args.seed,args.block_sweeps)
        args.output.mkdir(parents=True,exist_ok=True)
        if (args.output/'amplitude_scan.json').exists():
            raise ValueError('refusing to overwrite existing analysis')
        (args.output/'amplitude_scan.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        if result['binder_slopes']:
            write_csv(args.output/'binder_slopes.csv',result['binder_slopes'])
        plots(result,args.output)
        print(f"wrote {args.output}; {len(result['points'])} lambda points, fixed q=({args.n1},{args.n2})")
    except (ValueError,KeyError) as e:
        p.error(str(e))


if __name__=='__main__':
    main()
