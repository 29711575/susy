"""Python 2/3 remote reader. Executed over SSH stdin; never changes jobs/files."""
from __future__ import print_function
import datetime
import glob
import json
import os
import re
import subprocess
import sys


def command(args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return p.returncode, out.decode('utf-8', 'replace'), err.decode('utf-8', 'replace')


def fields(text):
    jobs = []
    current, key = None, None
    for line in text.splitlines():
        if line.startswith('Job Id:'):
            current = {'job_id': line.split(':', 1)[1].strip()}
            jobs.append(current); key = None
        elif current is not None:
            match = re.match(r'^    ([A-Za-z_][\w.]*) = (.*)$', line)
            if match:
                key = match.group(1); current[key] = match.group(2)
            elif line.startswith('\t') and key:
                current[key] += line.strip()
    return jobs


def load(path):
    try:
        with open(path) as stream:
            return json.load(stream)
    except (IOError, ValueError):
        return None


def last_record(path):
    try:
        with open(path, 'rb') as stream:
            stream.seek(0, 2); cap = stream.tell()
            size = min(cap, 65536)
            while True:
                stream.seek(cap-size); data = stream.read(size)
                lines = data.split(b'\n')
                if len(lines) >= 3 or size == cap:
                    break
                size = min(cap, size*2)
            if len(lines) < 2:
                return None
            row = lines[-2].decode('utf-8').split(',')
            return row
    except (IOError, ValueError):
        return None


def root_data(root):
    starts, points, folders, offsets = {}, {}, {}, {}
    for f in glob.glob(root+'/runs/*/*.start.json')+glob.glob(root+'/chains/*/*.start.json'):
        item = load(f)
        if isinstance(item, dict) and item.get('pbs_jobid'):
            starts.setdefault(item['pbs_jobid'], []).append(f[:-11])
    for f in [root+'/points.json', root+'/plan.json']:
        item = load(f)
        if isinstance(item, dict):
            item = item.get('chains', item.get('points', []))
        if isinstance(item, list):
            for point in item:
                if isinstance(point, dict) and 'index' in point:
                    points[point['index']] = point
    for f in glob.glob(root+'/runs/*/point.json')+glob.glob(root+'/chains/*/point.json'):
        point = load(f)
        if isinstance(point, dict) and 'index' in point:
            folders.setdefault(point['index'], []).append(os.path.dirname(f))
    for f in glob.glob(root+'/*.pbs'):
        try:
            with open(f) as stream:
                text = stream.read()
            name = re.search(r'^#PBS\s+-N\s+(\S+)\s*$', text, re.M)
            offset = re.search(r'^export CHAIN_OFFSET=(\d+)\s*$', text, re.M)
            if name and offset:
                offsets.setdefault(name.group(1), []).append((int(offset.group(1)), f))
        except IOError:
            pass
    source = load(root+'/manifest.json') or load(root+'/runtime-origin.json') or {}
    return starts, points, source, folders, offsets


def collect():
    started = datetime.datetime.utcnow().isoformat()+'Z'
    user = os.environ.get('USER', 'sunxr')
    rc, table, err = command(['qstat', '-t', '-u', user])
    if rc:
        raise RuntimeError('qstat list failed: '+err)
    ids = [line.split()[0] for line in table.splitlines()
           if line.split() and re.match(r'^\d+(?:\[\d+\])?\.[\w.-]+$', line.split()[0])]
    details, errors = [], []
    for start in range(0, len(ids), 50):
        rc, text, err = command(['qstat', '-f']+ids[start:start+50])
        details.extend(fields(text))
        if rc:
            errors.append(err.strip())  # A job may end between the two reads.
    cache = {}
    for job in details:
        var = job.get('Variable_List', '')
        match = re.search(r'(?:^|,)PBS_O_WORKDIR=([^,]+)', var)
        root = match.group(1) if match else None
        if not root:
            job['inventory_error'] = 'PBS_O_WORKDIR unavailable'; continue
        job['root'] = root
        if root not in cache:
            cache[root] = root_data(root)
        starts, points, source, folders, offsets = cache[root]
        job['source'] = dict((k, source[k]) for k in ['binary_sha256', 'source_archive_sha256', 'origin_runtime'] if k in source)
        prefixes = starts.get(job['job_id'], [])
        job['recorded_run_prefixes'] = prefixes
        primary = [p for p in prefixes if os.path.basename(p) == 'mc']
        # Historical wrappers can record a control/preflight and production in
        # the same PBS child. Preserve all candidates, select only named mc.
        if len(primary) == 1:
            prefixes = primary
        index = re.search(r'\[(\d+)\]', job['job_id'])
        job['index'] = int(index.group(1)) if index else None
        global_index = job['index']
        matches = offsets.get(job.get('Job_Name'), [])
        job['mapping_source'] = 'PBS job id in run start record' if prefixes else 'points index'
        if not prefixes and len(matches) == 1 and global_index is not None:
            global_index += matches[0][0]
            job['mapping_source'] = 'PBS script name plus literal CHAIN_OFFSET: '+matches[0][1]
        elif not prefixes and any(p.get('array_index', p['index']) != p['index'] for p in points.values()):
            job['inventory_error'] = 'array/global index mapping unresolved'
            global_index = None
        job['global_index'] = global_index
        job['point'] = points.get(global_index, {})
        if not prefixes and len(folders.get(global_index, [])) == 1:
            folder = folders[global_index][0]
            candidate = folder+'/mc'
            if os.path.isfile(candidate+'_results.json'):
                prefixes = [candidate]
                job['mapping_source'] += '; unique recorded point directory'
        if len(prefixes) > 1:
            job['inventory_error'] = 'multiple recorded run prefixes for this PBS job'
        if len(prefixes) == 1:
            prefix = prefixes[0]
            job['output_prefix'] = prefix
            job['argv'] = load(prefix+'.argv.json') or load(os.path.dirname(prefix)+'/argv.json') or load(os.path.dirname(prefix)+'/planned.argv.json') or []
            job['point'] = load(os.path.dirname(prefix)+'/point.json') or job['point']
            result = load(prefix+'_results.json') or {}
            if not job['argv']:
                job['argv'] = result.get('argv', [])
            job['run_status'] = result.get('status', 'unknown')
            job['result_numerical_failures'] = result.get('numerical_failures')
            job['result_is_startup_placeholder'] = result.get('status') == 'running'
            job['resume'] = result.get('resume', {})
            if result.get('status') == 'complete':
                job['complete_metadata'] = dict((k, result[k]) for k in [
                    'Lx', 'Ly', 'shear', 'g', 'mass', 'omega_ph', 'lambda', 'tau', 'dt', 'intervals',
                    'seed', 'measure_every', 'matter_trial', 'matter_trial_sector'] if k in result)
            record = last_record(prefix+'_blocks.csv')
            if record and len(record) >= 7 and record[0].isdigit():
                job['complete_block_records'] = int(record[0])+1
                job['stored_samples'] = (int(record[0])+1)*int(record[1])
            acceptance = last_record(prefix+'_acceptance.csv')
            if acceptance and len(acceptance) >= 11 and acceptance[1].isdigit():
                job['acceptance_latest'] = dict(stage=acceptance[0], sweep=int(acceptance[1]),
                    proposed=int(acceptance[2]), accepted=int(acceptance[3]),
                    first_failures=int(acceptance[6]), reference_failed=int(acceptance[8]),
                    elapsed_seconds=float(acceptance[-1]))
            job['checkpoint'] = load(prefix+'_checkpoint.receipt.json')
    observed = set(j['job_id'] for j in details)
    # Environment values are needed only to locate PBS_O_WORKDIR. Do not retain
    # the user's complete login environment in inventory artifacts.
    for job in details:
        job.pop('Variable_List', None)
    return dict(schema='qmc-remote-inventory-v1', snapshot_started_utc=started,
                snapshot_finished_utc=datetime.datetime.utcnow().isoformat()+'Z',
                jobs=details, qstat_errors=errors, disappeared_during_snapshot=sorted(set(ids)-observed),
                access='read_only', scheduler_listing=table)


if __name__ == '__main__':
    print(json.dumps(collect(), separators=(',', ':')))
