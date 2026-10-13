"""Check implicit production defaults against the previously exported binary."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import hashlib
import json
import os
import subprocess

ROOT = Path('/home/asus/newsusy/quick1010')
AUDIT = ROOT / 'audit/1013/core-numerics-defaults'
NEW = AUDIT / 'build/kitaev_ssh_mc'
OLD = ROOT / 'audit/1013/github-latest-source-publication/build/kitaev_ssh_mc'
RUNS = AUDIT / 'runs'
RUNS.mkdir(exist_ok=True)
ENV = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_DYNAMIC='FALSE', LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2026.1/lib')
for key in ['PFQMC_NO_A1', 'PFQMC_NO_A2', 'PFQMC_NO_BIPARTITE_SPECTRUM', 'PFQMC_NO_OVERLAP_SVD', 'PFQMC_NO_REAL_MAGNITUDE']:
    ENV.pop(key, None)
PRODUCTION = ['--weight-engine', 'combined', '--defer-phase', '1', '--prefix-cache', '1', '--cayley-cache', '1', '--center-green-fast', '1']
REFERENCE = ['--weight-engine', 'reference', '--defer-phase', '0', '--prefix-cache', '0', '--cayley-cache', '0', '--center-green-fast', '0']
common = ['--intervals', '4', '--dt', '.1', '--mass', str(100/9), '--omega-ph', '.3', '--g', str(3**.5), '--seed', '106130901', '--warmup', '1', '--sweeps', '2', '--measure-every', '1', '--block-samples', '1', '--vbs', '1', '--vbs-mode', 'all', '--vbs-susceptibility', '1', '--measure-energy', '1', '--fermions', '0']
cases = {
    'legacy_bond': ['--Lx', '6', '--Ly', '6', '--proposal', 'bond-refresh', '--q-mode-every', '0'],
    'guided_q': ['--Lx', '6', '--Ly', '6', '--proposal', 'bond-guided', '--bond-guide', 'bare', '--bond-beta', '1', '--q-mode-every', '1', '--q-mode-n1', '3', '--q-mode-n2', '3'],
    'local_reference': ['--Lx', '3', '--Ly', '2', '--shear', '1', '--proposal', 'local', '--q-mode-every', '0'],
}
jobs = []
for case, settings in cases.items():
    variants = [('old', OLD, REFERENCE if case == 'local_reference' else PRODUCTION),
                ('new', NEW, ['--numerics', 'reference'] if case == 'local_reference' else []),
                ('explicit', NEW, REFERENCE if case == 'local_reference' else PRODUCTION)]
    if case == 'legacy_bond':
        variants += [('named', NEW, ['--numerics', 'production']),
                     ('old_center0', OLD, PRODUCTION[:-1] + ['0']),
                     ('new_center0', NEW, ['--center-green-fast', '0'])]
    for variant, binary, options in variants:
        prefix = RUNS / f'{case}-{variant}'
        jobs.append((case, variant, [str(binary), *common, *settings, *options, '--output', str(prefix)], prefix))

def run(job):
    case, variant, argv, prefix = job
    result = subprocess.run(argv, env=ENV, capture_output=True, text=True)
    prefix.with_suffix('.stdout').write_text(result.stdout)
    prefix.with_suffix('.stderr').write_text(result.stderr)
    assert result.returncode == 0, (case, variant, result.stderr[-2000:])
    data = json.loads(Path(str(prefix) + '_results.json').read_text())
    assert data['status'] == 'complete'
    return (case, variant), data

with ThreadPoolExecutor(max_workers=2) as pool:
    results = dict(pool.map(run, jobs))

comparisons = []
fields = ['weight_engine', 'defer_phase', 'prefix_cache', 'cayley_cache', 'center_green_fast', 'proposal', 'samples', 'move_counts', 'numerical_failures', 'final_phonon_path']
for case, settings in cases.items():
    baseline = results[case, 'old']
    for variant in ['new', 'explicit'] + (['named'] if case == 'legacy_bond' else []):
        data = results[case, variant]
        for field in fields:
            assert data[field] == baseline[field], (case, variant, field)
        a = Path(str(RUNS / f'{case}-old') + '_blocks.csv').read_bytes()
        b = Path(str(RUNS / f'{case}-{variant}') + '_blocks.csv').read_bytes()
        assert a == b, (case, variant, 'raw blocks differ')
        comparisons.append({'case': case, 'variant': variant, 'raw_blocks_byte_identical': True, 'profile': data['numerics_profile']})
    if case != 'local_reference':
        assert results[case, 'new']['numerics_profile'] == 'production-v1'
    else:
        assert results[case, 'new']['numerics_profile'] == 'reference-v1'
    if case == 'legacy_bond':
        for field in fields:
            assert results[case, 'old_center0'][field] == results[case, 'new_center0'][field], field
        assert Path(str(RUNS / f'{case}-old_center0') + '_blocks.csv').read_bytes() == Path(str(RUNS / f'{case}-new_center0') + '_blocks.csv').read_bytes()
        comparisons.append({'case': case, 'variant': 'center0 override', 'raw_blocks_byte_identical': True})

help_result = subprocess.run([str(NEW), '--help'], env=ENV, capture_output=True, text=True, check=True)
(AUDIT / 'help.txt').write_text(help_result.stdout)
assert '--numerics reference' in help_result.stdout
assert '--center-green-fast 1:' not in help_result.stdout
invalid = subprocess.run([str(NEW), *common, '--numerics', 'typo', '--output', str(RUNS / 'invalid-profile')], env=ENV, capture_output=True, text=True)
assert invalid.returncode != 0 and '--numerics must be' in invalid.stderr
assert not Path(str(RUNS / 'invalid-profile') + '_blocks.csv').exists()
summary = {'comparisons': comparisons, 'successful_cli_runs': len(jobs), 'invalid_profile_rejected': True,
           'old_binary_sha256': hashlib.sha256(OLD.read_bytes()).hexdigest(), 'new_binary_sha256': hashlib.sha256(NEW.read_bytes()).hexdigest(),
           'scope': 'Short deterministic fixed-seed compatibility checks; no production submission or sampling-efficiency certification.'}
(AUDIT / 'cli-checks.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(summary, ensure_ascii=False, indent=2))
