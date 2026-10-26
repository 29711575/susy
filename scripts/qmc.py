#!/usr/bin/env python3
"""One entry point for offline final analysis, diagnostics, and PBS inventory."""
import argparse
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import subprocess
from zoneinfo import ZoneInfo

# Each chain/analysis process uses one math thread, including NumPy BLAS.
for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[name] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qmc-matplotlib-'+str(os.getuid()))

from qmc_tools import VERSION

BASE = Path(__file__).resolve().parents[1]


def clean(value):
    if isinstance(value, dict):
        return {k:clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name+'.tmp')
    with temp.open('x') as stream:
        json.dump(clean(value), stream, indent=2, allow_nan=False, ensure_ascii=False)
        stream.write('\n')
    temp.replace(path)


def output_dir(path, name):
    now = datetime.now(ZoneInfo('Asia/Shanghai'))
    out = Path(path) if path else BASE/'audit'/now.strftime('%m%d')/'qmc-tools'/(name+'-'+now.strftime('%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    return out


def code_manifest():
    files = [Path(__file__)]+sorted(Path(__file__).with_name('qmc_tools').glob('*.py'))
    return {str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def load_inputs(args):
    from qmc_tools.blocks import checked_group, read_chain
    paths = list(args.inputs)
    if args.selection:
        paths += json.loads(Path(args.selection).read_text())['local_blocks']
    if not paths:
        raise ValueError('provide input paths or --selection from registry fetch')
    chains = [read_chain(p, args.metadata, args.allow_prefix) for p in paths]
    checked_group(chains)
    return [c.window(args.cut_sweeps, args.end_sweeps) for c in chains]


def write_estimates_csv(path, levels):
    import csv
    with path.open('w', newline='') as stream:
        fields = ['factor', 'samples', 'observable', 'real', 'imag', 'error_real', 'error_imag', 'covariance_real_imag', 'defined']
        writer = csv.DictWriter(stream, fields); writer.writeheader()
        for level in levels:
            for name, estimate in level.get('estimates', {}).items():
                writer.writerow(dict(factor=level['factor'], samples=level['samples'], observable=name, **clean(estimate)))


def analyze(args):
    from qmc_tools.observables import transform_for
    from qmc_tools.statistics import joint
    chains = load_inputs(args)
    transform = (lambda r: {}) if args.raw_only else transform_for(chains[0].names, chains[0].metadata)
    levels = []
    factors = args.reblock or [1, 2, 4, 8, 16]
    if len(set(factors)) != len(factors) or min(factors) < 1:
        raise ValueError('reblock factors must be distinct positive integers')
    for factor in factors:
        try:
            result = joint([c.reblock(factor) for c in chains], transform)
            result.update(factor=factor, status='estimated')
        except ValueError as error:
            result = dict(factor=factor, status='unavailable', reason=str(error))
        levels.append(result)
    out = output_dir(args.out, 'analysis')
    report = dict(schema=VERSION+'-analysis', status='analysis_finished',
                  inputs=[c.provenance for c in chains],
                  metadata=[c.metadata for c in chains],
                  window=dict(cut_sweeps=args.cut_sweeps, end_sweeps=args.end_sweeps),
                  windows=[dict(path=str(c.path), start_sweep=c.start_sweep, end_sweep=int(c.ends[-1]), samples=int(c.counts.sum())) for c in chains],
                  phase_counts=[dict(path=str(c.path), positive=int(c.classifications[:, 0].sum()),
                                     negative=int(c.classifications[:, 1].sum()), nonreal=int(c.classifications[:, 2].sum())) for c in chains],
                  analysis_does_not_establish_sampling_completion=True,
                  block_independence_verified=False, equilibration_certified=False,
                  derived_omissions=getattr(transform, 'omissions', []), levels=levels,
                  statistics='complex common-phase ratios; stratified joint weighted delete-group pseudovalues; never reblock across chains',
                  code_sha256=code_manifest(), command=sys.argv)
    write_json(out/'analysis.json', report)
    write_estimates_csv(out/'estimates.csv', levels)
    (out/'README.txt').write_text('从冻结的完整原始块重新分析；未运行蒙特卡洛。\n'
        '每个分块层均保留全部输入样本；短尾留在原链，按实际样本比例作加权删组估计。\n'
        '中途/失败任务的可用前缀不等同于计划已完成；检查 analysis.json 的原始状态、窗口和限制。\n'
        '误差为有限样本联合 jackknife，未认证平衡、块独立性、有限投影或尺寸收敛。\n')
    valid = sum(level['status'] == 'estimated' for level in levels)
    print(json.dumps(dict(output=str(out.resolve()), estimated_levels=valid, unavailable_levels=len(levels)-valid)))
    return 0 if valid else 2


def diagnostic_plot(out, chains, selected):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(selected), 2, squeeze=False, figsize=(12, 3*len(selected)))
    for row, name in enumerate(selected):
        for c in chains:
            lookup = {n:i+1 for i, n in enumerate(c.names)}
            if name in {'Xplus', 'Xminus'}:
                num = (c.sums[:, lookup['X_x']]+(1 if name == 'Xplus' else -1)*c.sums[:, lookup['X_y']])/(2 if name == 'Xplus' else 1)
            else:
                num = c.sums[:, lookup[name]]
            phase = c.sums[:, 0]
            raw = np.full(len(num), np.nan, dtype=complex)
            np.divide(num, phase, out=raw, where=np.abs(phase)>1e-12*c.counts)
            cumulative = np.cumsum(num); denom = np.cumsum(phase)
            curve = np.full(len(num), np.nan, dtype=complex)
            np.divide(cumulative, denom, out=curve, where=np.abs(denom)>1e-12*np.cumsum(c.counts))
            label = str(c.metadata.get('seed', c.path.name))
            axes[row, 0].plot(c.ends, raw.real, lw=.6, alpha=.65, label=label)
            axes[row, 1].plot(c.ends, curve.real, lw=1, label=label)
        axes[row, 0].set_ylabel(name+' (real)')
        axes[row, 0].set_title('Stored-block phase ratios')
        axes[row, 1].set_title('Cumulative common-phase ratio')
        for ax in axes[row]:
            ax.set_xlabel('measurement sweep'); ax.grid(alpha=.2)
    axes[0, 1].legend(title='seed', fontsize=7)
    fig.suptitle('Descriptive trajectories; equilibration is not certified')
    fig.tight_layout(rect=(0, 0, 1, .975)); fig.savefig(out/'trajectories.png', dpi=160); plt.close(fig)


def diagnose_command(args):
    import csv
    from qmc_tools.diagnostics import diagnose
    chains = load_inputs(args)
    result, traces = diagnose(chains, args.observable)
    out = output_dir(args.out, 'diagnostics')
    result.update(schema=VERSION+'-diagnostics', inputs=[c.provenance for c in chains],
                  code_sha256=code_manifest(), command=sys.argv)
    write_json(out/'diagnostics.json', result)
    for name, rows in [('rhat-ess.csv', result['diagnostics']), ('influence-traces.csv', traces)]:
        with (out/name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, list(rows[0])); writer.writeheader(); writer.writerows(clean(rows))
    selected = list(dict.fromkeys(r['observable'] for r in result['diagnostics']))
    if not args.no_plot:
        diagnostic_plot(out, chains, selected)
    # Reuse exactly the same analysis engine and requested window for blocking
    # and per-seed statistics, with no second read of growing source CSV files.
    from qmc_tools.statistics import joint
    levels = []
    for factor in [1, 2, 4, 8, 16]:
        try:
            levels.append(dict(factor=factor, **joint([c.reblock(factor) for c in chains])))
        except ValueError as error:
            levels.append(dict(factor=factor, status='unavailable', reason=str(error)))
    write_json(out/'blocking.json', levels)
    write_estimates_csv(out/'blocking.csv', levels)
    write_json(out/'per-seed.json', [dict(seed=c.metadata.get('seed'), path=str(c.path), **joint([c])) for c in chains])
    cuts = []
    for cut in args.compare_cuts or []:
        try:
            subset = [c.window(cut) for c in chains]
            diagnostic, _ = diagnose(subset, selected)
            cuts.append(dict(cut_sweeps=cut, status='descriptive_only',
                             statistics=joint(subset), diagnostics=diagnostic))
        except ValueError as error:
            cuts.append(dict(cut_sweeps=cut, status='unavailable', reason=str(error)))
    write_json(out/'cut-sensitivity.json', cuts)
    print(json.dumps(dict(output=str(out.resolve()), independent_chains=len(chains),
                          common_block_records=result['common_block_records'], equilibration_certified=False)))
    return 0


def registry_command(args):
    from qmc_tools import registry as reg
    if args.action == 'refresh':
        from concurrent.futures import ThreadPoolExecutor
        ports = args.port or [20022, 20015]
        previous = json.loads(Path(args.previous).read_text()) if args.previous else None
        with ThreadPoolExecutor(max_workers=min(2, len(ports))) as pool:
            snapshots = list(pool.map(lambda port: reg.collect_remote(port, args.host, args.ssh_key), ports))
        result = reg.merge_snapshots(snapshots, previous)
        out = output_dir(args.out, 'registry')
        for snap in snapshots:
            write_json(out/f"snapshot-{snap['port']}.json", snap)
        write_json(out/'registry.json', result); reg.export_csv(out/'tasks.csv', result['records'])
        reg.export_cells(out/'cells.csv', result['records'])
        print(json.dumps(dict(output=str(out.resolve()), records=len(result['records']), issues=len(result['issues']))))
    elif args.action == 'import':
        result = reg.import_submission(args.submission, args.points)
        if args.previous:
            previous = json.loads(Path(args.previous).read_text())
            by_key = {(r['port'], r['job_id']):r for r in result['records']}
            for record in previous['records']:
                key = (record['port'], record['job_id'])
                by_key[key] = {**by_key.get(key, {}), **record}
            result['records'] = list(by_key.values())
        out = output_dir(args.out, 'registry-import')
        write_json(out/'registry.json', result); reg.export_csv(out/'tasks.csv', result['records'])
        reg.export_cells(out/'cells.csv', result['records'])
        print(json.dumps(dict(output=str(out.resolve()), records=len(result['records']))))
    else:
        result = json.loads(Path(args.registry).read_text())
        records = reg.select(result, args)
        if args.action == 'fetch':
            out = output_dir(args.out, 'snapshot')
            files = reg.fetch(records, out, args.host, args.ssh_key)
            write_json(out/'selection.json', dict(records=records, local_blocks=files))
            print(json.dumps(dict(output=str(out.resolve()), chains=len(files))))
            return 0
        if args.out:
            out = output_dir(args.out, 'selection')
            write_json(out/'registry.json', dict(schema=result['schema'], records=records))
            reg.export_csv(out/'tasks.csv', records)
        print('port\tjob\tstate\tL\tomega\ttau\tlambda\tseed\tsweeps\tobserved_utc')
        for r in records:
            state = r['scheduler_state'] if r.get('currently_listed') else f"{r['scheduler_state']}(historical)"
            print('\t'.join(str(v if v is not None else '?') for v in [r['port'], r['job_id'], state, r.get('Lx'), r.get('omega_ph'), r.get('tau'), r.get('lambda'), r.get('seed'), r.get('stored_sweeps'), r['observed_utc']]))
        print(f'{len(records)} records', file=sys.stderr)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='command', required=True)
    for name in ['analyze', 'diagnose']:
        p = subs.add_parser(name)
        p.add_argument('inputs', nargs='*', help='raw blocks CSV, result JSON, prefix, or run directory')
        p.add_argument('--selection', help='selection.json produced by registry fetch')
        p.add_argument('--metadata', help='explicit metadata JSON for incomplete legacy records')
        p.add_argument('--allow-prefix', action='store_true', help='explicitly analyze complete records from running/failed/unknown jobs')
        p.add_argument('--cut-sweeps', type=int, default=0)
        p.add_argument('--end-sweeps', type=int)
        p.add_argument('--out')
        if name == 'analyze':
            p.add_argument('--reblock', nargs='+', type=int)
            p.add_argument('--raw-only', action='store_true')
        else:
            p.add_argument('--observable', action='append', help='raw observable name, Xplus, or Xminus; repeatable')
            p.add_argument('--no-plot', action='store_true')
            p.add_argument('--compare-cuts', nargs='+', type=int, help='additional block-aligned discard windows; no automatic selection')
    p = subs.add_parser('registry')
    commands = p.add_subparsers(dest='action', required=True)
    refresh = commands.add_parser('refresh')
    refresh.add_argument('--port', type=int, action='append', help='repeatable; default 20022 and 20015')
    refresh.add_argument('--host', default='sunxr@remote3.ginpie.com')
    refresh.add_argument('--ssh-key', default='/home/asus/.ssh/id_ed25519_hpc')
    refresh.add_argument('--previous')
    refresh.add_argument('--out')
    imp = commands.add_parser('import')
    imp.add_argument('submission'); imp.add_argument('--points'); imp.add_argument('--previous'); imp.add_argument('--out')
    for action in ['query', 'fetch']:
        query = commands.add_parser(action)
        query.add_argument('--registry', required=True)
        query.add_argument('--port', type=int); query.add_argument('--job'); query.add_argument('--state')
        query.add_argument('--L', type=int); query.add_argument('--omega', type=float); query.add_argument('--tau', type=float)
        query.add_argument('--lambda', dest='lambda_value', type=float); query.add_argument('--seed', type=int); query.add_argument('--out')
        if action == 'fetch':
            query.add_argument('--host', default='sunxr@remote3.ginpie.com')
            query.add_argument('--ssh-key', default='/home/asus/.ssh/id_ed25519_hpc')
    args = parser.parse_args()
    if getattr(args, 'cut_sweeps', 0) < 0:
        parser.error('cut-sweeps must be nonnegative')
    try:
        return analyze(args) if args.command == 'analyze' else diagnose_command(args) if args.command == 'diagnose' else registry_command(args)
    except (ValueError, OSError, KeyError, subprocess.TimeoutExpired) as error:
        print('qmc: '+str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
