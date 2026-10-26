"""Gate on complete successful files, avoiding already-finished-array PBS holds."""
from pathlib import Path
import json,time,runpy,subprocess
P=Path.cwd();plan=json.loads((P/'pilot-plan.json').read_text());ed=Path('/home/sunxr/quick1010-pcn-projection-ed-20261019-20022');start=time.monotonic()
while True:
 files=[P/'pilot/runs'/f"chain{r['index']:03d}.exit"for r in plan['jobs']]+[ed/f'task{i:03d}.exit'for i in range(18)];missing=[]
 for file in files:
  if not file.exists()or not file.read_text().strip():missing.append(str(file));continue
  if file.read_text().strip()!='0':raise RuntimeError('prerequisite failed, no formal samples submitted: '+str(file))
 (P/'gate-state.json').write_text(json.dumps(dict(waiting_for=missing,wait_seconds=time.monotonic()-start),indent=2)+'\n')
 if not missing:break
 if time.monotonic()-start>21600:raise RuntimeError('prerequisite completion not established within6h; no formal samples submitted')
 time.sleep(30)
runpy.run_path(str(P/'plan_formal.py'),run_name='__main__')
