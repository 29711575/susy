"""All registered seeds, including failures, must terminate before reporting."""
from pathlib import Path
import json,subprocess,time
P=Path.cwd()
def maybe():
 files=[]
 for tier in ['formal','pcn']:
  for r in json.loads((P/(tier+'-plan.json')).read_text())['jobs']:files.append(P/tier/'runs'/f"chain{r['index']:03d}.exit")
 missing=[f for f in files if not f.exists()or not f.read_text().strip()]
 if missing:return len(missing)
 lock=P/'report-submit.lock'
 try:lock.mkdir()
 except FileExistsError:return 0
 cmd=['qsub','report.pbs'];r=subprocess.run(cmd,text=True,capture_output=True);(P/'report-submission.json').write_text(json.dumps(dict(command=cmd,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr,policy='all terminal markers including failed seeds; no after-completed-array dependency'),indent=2)+'\n');return 0
for attempt in range(10):
 missing=maybe()
 if missing==0 or missing>4:break
 time.sleep(3)
