"""Use unchanged CLI runner, then mark completion and trigger one final report."""
from pathlib import Path
import os,sys,json,subprocess,runpy
P=Path.cwd();r=subprocess.run([sys.executable,str(P/'run.py')]);plan=json.loads((P/os.environ['PLAN']).read_text());row=plan['jobs'][int(os.environ['PBS_ARRAY_INDEX'])];prefix=P/plan['tier']/'runs'/f"chain{row['index']:03d}";prefix.parent.mkdir(parents=True,exist_ok=True)
if not Path(str(prefix)+'.exit').exists():
 Path(str(prefix)+'.exit').write_text(str(r.returncode or 1)+'\n');Path(str(prefix)+'.wrapper-failure.json').write_text(json.dumps(dict(worker_return_code=r.returncode,reason='worker did not produce its normal terminal record; preserve PBS log and all seed data'))+'\n')
runpy.run_path(str(P/'maybe_submit_report.py'),run_name='__main__');sys.exit(r.returncode)
