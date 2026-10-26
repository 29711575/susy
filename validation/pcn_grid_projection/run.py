"""Execute one frozen plan row through the actual production CLI."""
import os,json,subprocess,time,hashlib,socket,sys
from pathlib import Path
P=Path.cwd();planfile=P/os.environ['PLAN'];plan=json.loads(planfile.read_text());row=plan['jobs'][int(os.environ['PBS_ARRAY_INDEX'])];tier=plan['tier'];folder=P/tier/'runs';folder.mkdir(parents=True,exist_ok=True);prefix=folder/f"chain{row['index']:03d}"
build=Path(os.environ['VALIDATION_BUILD']);assert (build/'build.exit').read_text().strip()=='0';assert json.loads((build/'check-cli.json').read_text())['status']=='PASS'
with open(str(prefix)+'.integrity','w')as f:subprocess.run(['sha256sum','-c','binary.sha256'],cwd=build/'runtime',stdout=f,check=True)
env=os.environ.copy();env.update({k:'1'for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']});env['MKL_DYNAMIC']='FALSE';env['LD_LIBRARY_PATH']=str(build/'runtime/lib')+':/opt/ohpc/pub/apps/intel/oneapi/mkl/2021.2.0/lib/intel64:'+env.get('LD_LIBRARY_PATH','')
cmd=[str(build/'runtime/kitaev_ssh_mc')]+row['args']+['--output',str(prefix)];t=time.monotonic()
with open(str(prefix)+'.stdout','w')as f,open(str(prefix)+'.stderr','w')as e:r=subprocess.run(cmd,env=env,stdout=f,stderr=e)
record=dict(command=cmd,exit_code=r.returncode,seconds=time.monotonic()-t,plan_sha256=hashlib.sha256(planfile.read_bytes()).hexdigest(),host=socket.gethostname(),PBS_JOBID=os.environ.get('PBS_JOBID'),seed=row['seed'],threads={k:env[k]for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']},binary_sha256=hashlib.sha256((build/'runtime/kitaev_ssh_mc').read_bytes()).hexdigest());Path(str(prefix)+'.execution.json').write_text(json.dumps(record,indent=2)+'\n');Path(str(prefix)+'.exit').write_text(str(r.returncode)+'\n');sys.exit(r.returncode)
