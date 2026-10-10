#!/usr/bin/env python3
"""Bounded local tau scan; runs the unchanged MC proposal/weight algorithm."""
import argparse,concurrent.futures,datetime,hashlib,json,os,pathlib,shutil,subprocess,sys,time

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=pathlib.Path,required=True);args=ap.parse_args()
    root=pathlib.Path(__file__).resolve().parents[1];run=args.run_dir.resolve();plan=json.loads((root/'kitaev_ssh/phase_tau_scan_plan.json').read_text())
    run.mkdir(parents=True,exist_ok=True)
    if (run/'manifest.json').exists(): raise RuntimeError('refusing to overwrite existing scan')
    # Submission gate: completed tests and bit-identical old/new transaction trace.
    validation=run/'validation'
    for tag in ('before','after'):
        if '100% tests passed' not in (validation/f'ctest_{tag}.log').read_text():raise RuntimeError('ctest gate')
        json.loads((validation/f'weight_{tag}.json').read_text())
        if 'PASS' not in (validation/f'phonon_{tag}.log').read_text():raise RuntimeError('phonon gate')
    if (validation/'trace_before/transactions.txt').read_bytes()!=(validation/'trace_after/transactions.txt').read_bytes():raise RuntimeError('transaction trace mismatch')
    if not json.loads((validation/'cli_comparison.json').read_text())['passed']:raise RuntimeError('CLI regression gate')
    binary=run/'kitaev_ssh_mc';shutil.copy2(root/'build/kitaev_ssh_mc',binary)
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    env['MKL_DYNAMIC']='FALSE';env['PYTHONDONTWRITEBYTECODE']='1'
    commands=[]
    for job in plan['jobs']:
        folder=run/f"m{job['intervals']:03d}_seed{job['seed']}";folder.mkdir(exist_ok=False)
        options={'Lx':4,'Ly':4,'Jx':1,'Jy':1,'Jz':.5,'g':.5,'mass':4,'omega-ph':.5,'dt':.1,'intervals':job['intervals'],'seed':job['seed'],'warmup':plan['warmup'],'sweeps':plan['sweeps'],'measure-every':1,'block-samples':25,'qr-block':4,'vbs':1,'proposal':'bond-refresh','trial-dx':0,'trial-dy':0,'trial-dz':0,'weight-engine':'combined','output':str(folder/'run')}
        command=[str(binary)]+[str(x) for k,v in options.items() for x in ('--'+k,v)]
        commands.append({'job':job,'folder':str(folder),'argv':command})
    source_files=['phase_statistics.h','coupled_mc.cpp','phonon_mc.h','phonon_projector.h','phonon_path.h','vbs_observables.h']
    manifest={'started_utc':now(),'host':os.uname().nodename,'pid':os.getpid(),'execution':'local WSL, up to 12 processes; no PBS job','plan':plan,'binary_sha256':digest(binary),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'source_hashes':{n:digest(root/'kitaev_ssh'/n) for n in source_files},'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'mkl_path':str(pathlib.Path(env['MKLROOT']).resolve()),'threads':{k:env[k] for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_DYNAMIC')},'commands':commands}
    (run/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(run/'source.patch').write_bytes(subprocess.check_output(['git','diff'],cwd=root))
    statuses={str(c['job']['seed']):{'state':'queued',**c['job']} for c in commands}
    def save_status():
        temp=run/'status.tmp';temp.write_text(json.dumps({'updated_utc':now(),'jobs':statuses},indent=2));temp.replace(run/'status.json')
    save_status()
    def execute(c):
        folder=pathlib.Path(c['folder']);start=time.monotonic()
        with (folder/'stdout.json').open('w') as out,(folder/'stderr.log').open('w') as err:
            proc=subprocess.Popen(c['argv'],stdout=out,stderr=err,env=env,cwd=root)
            (folder/'process.json').write_text(json.dumps({'pid':proc.pid,'started_utc':now()}))
            code=proc.wait()
        result={'exit_code':code,'elapsed_wall_seconds':time.monotonic()-start,'ended_utc':now()}
        if code==0:
            data=json.loads((folder/'run_results.json').read_text());counts=data['phase_diagnostics']
            if sum(counts[k] for k in ['positive_count','negative_count','complex_count'])!=data['samples']:raise RuntimeError('phase count conservation failed')
            result['phase_diagnostics']=counts
        (folder/'completion.json').write_text(json.dumps(result,indent=2));return result
    print('START',now(),'12 chains, workers=',plan['max_workers'],flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=plan['max_workers']) as pool:
        futures={pool.submit(execute,c):c for c in commands}
        for c in commands:statuses[str(c['job']['seed'])]['state']='running'
        save_status()
        for future in concurrent.futures.as_completed(futures):
            c=futures[future];seed=str(c['job']['seed'])
            try:
                result=future.result();statuses[seed].update(result,state='complete' if result['exit_code']==0 else 'failed')
            except Exception as exc:statuses[seed].update(state='failed',error=str(exc))
            save_status();print('FINISHED',seed,statuses[seed]['state'],now(),flush=True)
    if any(s['state']!='complete' for s in statuses.values()):return 1
    python=pathlib.Path('/home/asus/.cache/pfqmc_plot_env/bin/python')
    return subprocess.call([str(python),str(root/'kitaev_ssh/analyze_phase_tau_scan.py'),'--run-dir',str(run)],env=env,cwd=root)
if __name__=='__main__':sys.exit(main())
