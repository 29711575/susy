from pathlib import Path
import json,hashlib
out=Path(__file__).resolve().parent
results=[]
for status in json.loads((out/'run_status.json').read_text()):
 case={'name':status['name'],'exit':status['exit'],'observables':[],'phase':[]}
 for line in (out/(status['name']+'.log')).read_text().splitlines():
  tag,*items=line.split()
  if tag in ['METADATA','COUNTS','PHASE','OBS','TIMING']:
   fields={}
   for item in items:
    k,v=item.split('=',1)
    try: v=float(v) if any(c in v for c in '.e') else int(v)
    except ValueError: pass
    fields[k]=v
   if tag=='OBS':case['observables'].append(fields)
   elif tag=='PHASE':case['phase'].append(fields)
   else:case[tag.lower()]=fields
  elif tag in ['real','user','sys']:case[tag+'_seconds']=float(items[0])
 results.append(case)
(out/'generic_summary.json').write_text(json.dumps({'scope':'validation-only generic four-Majorana six-phonon core; production HoneycombLattice 2x1 CLI remains unsupported','cases':results},indent=2)+'\n')
files=['generic_driver.cpp','generic_driver','build_driver.sh','run_cases.py','collect_summary.py']
(out/'driver-sha256.txt').write_text(''.join(f'{hashlib.sha256((out/f).read_bytes()).hexdigest()}  {f}\n' for f in files))
for c in results:
 print(c['name'],c['real_seconds'], 'phase',c['phase'][0]['real'],'numerical_failures',c['counts']['first_failures'])
 for j in range(7):
  a=c['observables'][j];b=c['observables'][j+7]
  print(a['name'],a['real'],a['error_real'],b['error_real'])
