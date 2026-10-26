"""Read-only Python 2/3 snapshot transport. CONFIG is supplied by the caller."""
from __future__ import print_function
import datetime
import hashlib
import io
import json
import os
import sys
import tarfile


def add_bytes(archive, name, data):
    info = tarfile.TarInfo(name); info.size = len(data); info.mode = 0o600
    archive.addfile(info, io.BytesIO(data))


def complete_size(stream):
    stream.seek(0, 2); cap = stream.tell()
    size = min(cap, 65536)
    while size:
        stream.seek(cap-size); tail = stream.read(size)
        last = tail.rfind(b'\n')
        if last >= 0:
            return cap, cap-size+last+1
        if size == cap:
            break
        size = min(cap, size*2)
    raise ValueError('raw CSV has no complete record')


class HashReader(object):
    def __init__(self, stream):
        self.stream = stream; self.hash = hashlib.sha256()

    def read(self, size):
        data = self.stream.read(size); self.hash.update(data); return data


def main(config):
    stdout = getattr(sys.stdout, 'buffer', sys.stdout)
    archive = tarfile.open(fileobj=stdout, mode='w|gz')
    for record in config:
        prefix = record['output_prefix']; label = record['snapshot_label']
        root = os.path.realpath(record['root'])
        if not os.path.realpath(prefix).startswith(root+'/'):
            raise ValueError('recorded output is outside its PBS work directory')
        evidence = dict(port=record['port'], job_id=record['job_id'], remote_prefix=prefix,
                        snapshot_started_utc=datetime.datetime.utcnow().isoformat()+'Z', files={})
        # Read status before freezing blocks: a running producer may complete
        # during collection, but this snapshot never upgrades its old status.
        for suffix, target in [('_results.json', 'mc_results.json'), ('.argv.json', 'mc.argv.json'),
                               ('_matter_trial.json', 'mc_matter_trial.json')]:
            path = prefix+suffix
            if os.path.isfile(path):
                with open(path, 'rb') as stream:
                    data = stream.read()
                json.loads(data.decode('utf-8'))
                add_bytes(archive, label+'/'+target, data)
                evidence['files'][target] = dict(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
        point_path = os.path.dirname(prefix)+'/point.json'
        point = record.get('identity', {})
        if os.path.isfile(point_path):
            with open(point_path) as stream:
                point = json.load(stream)
        add_bytes(archive, label+'/point.json', json.dumps(point).encode('utf-8'))
        add_bytes(archive, label+'/source.json', json.dumps(record.get('source', {})).encode('utf-8'))
        with open(prefix+'_blocks.csv', 'rb') as stream:
            cap, used = complete_size(stream); stream.seek(0)
            reader = HashReader(stream)
            info = tarfile.TarInfo(label+'/mc_blocks.csv'); info.size = used; info.mode = 0o600
            archive.addfile(info, reader)
        evidence['files']['mc_blocks.csv'] = dict(sha256=reader.hash.hexdigest(), bytes_cap=cap, bytes_used=used,
                                                incomplete_tail_bytes=cap-used)
        evidence['snapshot_finished_utc'] = datetime.datetime.utcnow().isoformat()+'Z'
        add_bytes(archive, label+'/snapshot.json', json.dumps(evidence, indent=2).encode('utf-8'))
    archive.close()


if __name__ == '__main__':
    main(CONFIG)
