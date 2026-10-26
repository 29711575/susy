"""Port-qualified immutable observations, queries, and explicit manifest imports."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

from .blocks import identity


def now():
    return datetime.now(timezone.utc).isoformat()


def normalize(job, port, observed):
    meta = identity(job.get('complete_metadata', {}), job.get('argv'), job.get('point'))
    record = dict(port=port, job_id=job['job_id'], index=job.get('index'), root=job.get('root'),
                  scheduler_state=job.get('job_state', 'unknown'), observed_utc=observed,
                  state_source='live_qstat', output_prefix=job.get('output_prefix'),
                  run_status=job.get('run_status', 'unknown'), stored_samples=job.get('stored_samples'),
                  identity=meta, source=job.get('source', {}), inventory_error=job.get('inventory_error'),
                  acceptance_latest=job.get('acceptance_latest'), checkpoint=job.get('checkpoint'),
                  result_is_startup_placeholder=job.get('result_is_startup_placeholder', False),
                  resume=job.get('resume', {}))
    record['recorded_run_prefixes'] = job.get('recorded_run_prefixes', [])
    record['mapping_source'] = job.get('mapping_source')
    record['global_index'] = job.get('global_index', record['index'])
    for key in ['Lx', 'Ly', 'lambda', 'omega_ph', 'tau', 'dt', 'seed', 'replica']:
        record[key] = meta.get(key)
    record['identity_missing'] = [k for k in ['Lx', 'Ly', 'lambda', 'omega_ph', 'tau', 'dt', 'seed'] if record.get(k) is None]
    if record['stored_samples'] is not None and 'measure_every' in meta:
        record['stored_sweeps'] = record['stored_samples']*meta['measure_every']
    return record


def collect_remote(port, host, key, timeout=180):
    script = Path(__file__).with_name('remote_inventory.py').read_bytes()
    command = ['ssh', '-F', '/dev/null', '-o', 'BatchMode=yes', '-o', 'LogLevel=ERROR',
               '-o', 'ConnectTimeout=15', '-p', str(port)]
    if key:
        command += ['-i', str(Path(key).expanduser())]
    result = subprocess.run(command+[host, 'python -'], input=script, capture_output=True, timeout=timeout)
    if result.returncode:
        raise ValueError(f'port {port}: inventory SSH failed: '+result.stderr.decode(errors='replace')[-2000:])
    data = json.loads(result.stdout)
    data['port'] = port
    return data


def merge_snapshots(snapshots, previous=None):
    records = {(r['port'], r['job_id']): dict(r, state_source='previous_observation', currently_listed=False)
               for r in (previous or {}).get('records', [])}
    for snap in snapshots:
        for job in snap['jobs']:
            r = normalize(job, snap['port'], snap['snapshot_finished_utc'])
            r['currently_listed'] = True
            records[(r['port'], r['job_id'])] = r
    issues = []
    groups = {}
    for record in records.values():
        if record.get('currently_listed') and record.get('scheduler_state') in {'R', 'Q', 'H', 'E', 'W', 'S', 'T'} and record.get('seed') is not None:
            physical = {k:record['identity'].get(k) for k in ['Lx', 'Ly', 'shear', 'g', 'mass', 'omega_ph', 'tau', 'dt', 'trial_shift', 'J']}
            key = (record['seed'], json.dumps(physical, sort_keys=True), json.dumps(record['source'], sort_keys=True))
            groups.setdefault(key, []).append(record)
    for group in groups.values():
        if len(group) > 1:
            issues.append(dict(kind='repeated_seed_requires_review', jobs=[f"{r['port']}:{r['job_id']}" for r in group]))
    return dict(schema='qmc-task-registry-v1', generated_utc=now(),
                records=sorted(records.values(), key=lambda r: (r['port'], r['job_id'])), issues=issues,
                note='A missing job is not assumed finished/cancelled. Previous observations retain their timestamps.')


def indices(text):
    out = []
    for item in text.split(','):
        match = re.fullmatch(r'(\d+)(?:-(\d+)(?::(\d+))?)?', item.strip())
        if not match:
            raise ValueError('invalid PBS index range: '+item)
        start = int(match[1]); stop = int(match[2] or match[1]); step = int(match[3] or 1)
        if step < 1 or stop < start:
            raise ValueError('invalid PBS range')
        out.extend(range(start, stop+1, step))
    if len(out) != len(set(out)):
        raise ValueError('overlapping PBS index ranges')
    return out


def import_submission(path, points_path=None):
    path = Path(path).resolve()
    raw = path.read_bytes(); data = json.loads(raw)
    points_path = Path(points_path) if points_path else path.parent/'package/points.json'
    points = json.loads(points_path.read_text())
    by_index = {p['index']: p for p in points}
    if len(by_index) != len(points):
        raise ValueError('duplicated point indices')
    if 'port' not in data:
        raise ValueError('submission must record its port')
    jobs = []
    for group in data.get('groups', []):
        if group.get('exit', 0) != 0:
            continue
        parent = group['job']
        if '[]' not in parent:
            raise ValueError('expected an explicitly recorded array parent')
        for index in indices(group['indices']):
            if index not in by_index:
                raise ValueError('submitted index not found in points.json')
            jobs.append(dict(job_id=parent.replace('[]', f'[{index}]'), index=index, point=by_index[index],
                             root=data.get('root'), job_state='unknown'))
    snap = dict(port=data['port'], jobs=jobs, snapshot_finished_utc=data.get('utc', now()))
    registry = merge_snapshots([snap])
    for r in registry['records']:
        r.update(state_source='submission_record', currently_listed=False,
                 evidence=dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest()))
    return registry


def select(registry, args):
    out = []
    for r in registry['records']:
        match = True
        for attr, key in [('port', 'port'), ('L', 'Lx'), ('omega', 'omega_ph'), ('tau', 'tau'), ('lambda_value', 'lambda'), ('seed', 'seed')]:
            value = getattr(args, attr, None)
            if value is not None and (r.get(key) is None or abs(float(r[key])-float(value)) > 1e-10):
                match = False
        job = getattr(args, 'job', None)
        if job and not (r['job_id'] == job or r['job_id'].split('[')[0].split('.')[0] == job):
            match = False
        state = getattr(args, 'state', None)
        if state and (not r.get('currently_listed') or r['scheduler_state'] != state):
            match = False
        if match:
            out.append(r)
    return out


def export_csv(path, records):
    fields = ['port', 'job_id', 'scheduler_state', 'currently_listed', 'observed_utc', 'Lx', 'Ly', 'omega_ph', 'tau', 'dt', 'lambda', 'seed', 'stored_sweeps', 'run_status', 'output_prefix', 'state_source', 'inventory_error', 'identity_missing']
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(records)


def export_cells(path, records):
    cells = {}
    physical_keys = ['Lx', 'Ly', 'shear', 'Jx', 'Jy', 'Jz', 'g', 'mass', 'omega_ph', 'tau', 'dt', 'trial_shift', 'matter_trial_policy', 'arm', 'kind']
    for record in records:
        signature = {k:record['identity'].get(k) for k in physical_keys}
        signature['source'] = record.get('source', {})
        key = json.dumps(signature, sort_keys=True)
        cells.setdefault(key, []).append(record)
    rows = []
    for group in cells.values():
        first = group[0]
        row = {k:first.get(k) for k in ['Lx', 'Ly', 'omega_ph', 'tau', 'dt', 'lambda']}
        row.update(arm=first['identity'].get('arm', ''), source=json.dumps(first.get('source', {}), sort_keys=True),
                   records=len(group), active_seeds=len({r['seed'] for r in group if r.get('currently_listed') and r['scheduler_state'] in {'R', 'Q', 'H', 'W', 'S', 'E', 'T'}}),
                   R=0, Q=0, H=0, X=0, F=0, other=0, historical=0)
        for record in group:
            state = record['scheduler_state']
            if not record.get('currently_listed'):
                row['historical'] += 1
            else:
                row[state if state in {'R', 'Q', 'H', 'X', 'F'} else 'other'] += 1
        rows.append(row)
    with Path(path).open('w', newline='') as stream:
        fields = ['Lx', 'Ly', 'omega_ph', 'tau', 'dt', 'lambda', 'arm', 'records', 'active_seeds', 'R', 'Q', 'H', 'X', 'F', 'other', 'historical', 'source']
        writer = csv.DictWriter(stream, fields); writer.writeheader(); writer.writerows(rows)


def fetch(records, out, host, key):
    if not records or any(not r.get('output_prefix') for r in records):
        raise ValueError('fetch requires a nonempty selection with an unambiguous recorded output prefix')
    groups = {}
    for i, record in enumerate(records):
        copy = dict(record, snapshot_label=f'chain-{i:04d}-p{record["port"]}-'+record['job_id'].split('.')[0])
        groups.setdefault(record['port'], []).append(copy)
    for port, group in groups.items():
        config = 'import json\nCONFIG = json.loads('+repr(json.dumps(group))+')\n'
        script = config.encode()+Path(__file__).with_name('remote_fetch.py').read_bytes().replace(b'from __future__ import print_function\n', b'')
        command = ['ssh', '-F', '/dev/null', '-o', 'BatchMode=yes', '-o', 'LogLevel=ERROR',
                   '-o', 'ConnectTimeout=15', '-p', str(port)]
        if key:
            command += ['-i', str(Path(key).expanduser())]
        archive_path = out/f'snapshot-{port}.tar.gz'
        with archive_path.open('xb') as stream:
            result = subprocess.run(command+[host, 'python -'], input=script, stdout=stream,
                                    stderr=subprocess.PIPE, timeout=300)
        if result.returncode:
            raise ValueError('snapshot SSH failed: '+result.stderr.decode(errors='replace')[-2000:])
        with tarfile.open(archive_path) as archive:
            for member in archive:
                target = (out/member.name).resolve()
                if not member.isfile() or not target.is_relative_to(out.resolve()):
                    raise ValueError('unsafe snapshot archive entry')
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as stream:
                    stream.write(archive.extractfile(member).read())
        for record in group:
            folder = out/record['snapshot_label']
            evidence = json.loads((folder/'snapshot.json').read_text())
            for filename, receipt in evidence['files'].items():
                actual = hashlib.sha256((folder/filename).read_bytes()).hexdigest()
                if actual != receipt['sha256']:
                    raise ValueError('snapshot checksum mismatch: '+filename)
    return sorted(str(p) for p in out.glob('chain-*/mc_blocks.csv'))
