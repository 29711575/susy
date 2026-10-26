#!/usr/bin/env python3
"""Request a checkpoint from one new CLI/PBS child and wait for its receipt.

Does not qdel/qsub, signal a process group, or modify an existing executable.
Requires compute-node SSH from the login host and a matching control record.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


# Executed with the cluster's Python (including Python2.7 login nodes).
REMOTE = r'''
from __future__ import print_function
import json,os,re,subprocess,time,signal,socket,sys
job,control,stop,timeout=PARAMETERS
assert re.match(r'^\d+(?:\[\d+\])?(?:\.mgt)?$',job), 'Specify one job/array child, never an entire array'
assert os.path.isabs(control) and control.endswith('_checkpoint.control.json')
expected=job if job.endswith('.mgt') else job+'.mgt'
s=subprocess.check_output(['qstat','-f',expected]).decode()
assert re.search(r'job_state = R\b',s), 'Job is not running'
m=re.search(r'exec_host = ([A-Za-z0-9._-]+)/',s);assert m,'Missing execution host'
host=m.group(1)
c=json.load(open(control));assert c['pbs_jobid']==expected,'Control record belongs to another PBS job'
assert c['hostname'].split('.')[0]==host.split('.')[0], 'Control record belongs to another node'
assert c['pid']>1 and str(c['process_start_ticks']).isdigit()
receipt=control[:-len('_checkpoint.control.json')]+'_checkpoint.receipt.json'
old=json.load(open(receipt)).get('generation',0) if os.path.exists(receipt) else 0
code="""import os,json,signal,socket
c=json.loads(%r)
assert socket.gethostname().split('.')[0]==c['hostname'].split('.')[0]
pid=c['pid'];text=open('/proc/%%d/stat'%%pid).read();fields=text[text.rfind(')')+2:].split()
assert fields[19]==str(c['process_start_ticks']),'PID reused; refusing signal'
env=open('/proc/%%d/environ'%%pid,'rb').read().split(b'\\0')
assert (b'PBS_JOBID='+c['pbs_jobid'].encode()) in env,'PID is not part of the requested job'
os.kill(pid,signal.SIGUSR2 if %r else signal.SIGUSR1)
"""%(json.dumps(c),stop)
# The destination is from qstat, not a caller-supplied shell fragment.
try:
 from shlex import quote
except ImportError:
 from pipes import quote
subprocess.check_call(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15',host,'python -c '+quote(code)])
print(json.dumps(dict(signal_sent=True,job=expected,pid=c['pid'],host=host,stop=stop)));sys.stdout.flush()
deadline=time.time()+timeout
while time.time()<deadline:
 if os.path.exists(receipt):
  r=json.load(open(receipt))
  if r['pid']==c['pid'] and r['generation']>old and (not stop or r['stopped']):
   print(json.dumps(dict(checkpoint_saved=True,receipt=r)));sys.stdout.flush();sys.exit(0)
 time.sleep(1)
raise RuntimeError('No new successful checkpoint receipt before timeout; do not qdel or assume a checkpoint was saved')
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job', required=True, help='One PBS job or array child, e.g. 12345[7].mgt')
    p.add_argument('--control', required=True, help='Absolute shared-filesystem PREFIX_checkpoint.control.json')
    p.add_argument('--port', type=int, default=20022)
    p.add_argument('--host', default='sysu-hpc', help='SSH host/alias (default sysu-hpc)')
    p.add_argument('--ssh-config', default=str(Path.home()/'.ssh/config'))
    p.add_argument('--stop', action='store_true', help='Save then stop with exit75; otherwise save and continue')
    p.add_argument('--timeout', type=int, default=600, help='Receipt wait; one sweep may take minutes')
    a = p.parse_args()
    if a.timeout <= 0 or not os.path.isabs(a.control):
        p.error('--timeout must be positive and --control must be absolute')
    parameters = (a.job, a.control, a.stop, a.timeout)
    code = 'PARAMETERS='+repr(parameters)+'\n'+REMOTE
    command = ['ssh', '-F', a.ssh_config, '-p', str(a.port), '-o', 'BatchMode=yes',
               '-o', 'ConnectTimeout=15', a.host, 'python -u -c '+shlex.quote(code)]
    return subprocess.call(command)


if __name__ == '__main__':
    sys.exit(main())
