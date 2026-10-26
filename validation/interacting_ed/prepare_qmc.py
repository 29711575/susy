#!/usr/bin/env python3
"""Generate fixed smoke/full production CLI commands, never submit implicitly."""
import argparse,json,shlex
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--config',type=Path,default=Path(__file__).with_name('full-config.json'));p.add_argument('--output',type=Path,required=True);p.add_argument('--smoke',action='store_true');a=p.parse_args();cfg=json.loads(a.config.read_text());r=a.output;r.mkdir(parents=True,exist_ok=False);(r/'commands').mkdir();chains=[]
 # Longest estimates first. Exact execution placement is fixed after timing
 # pilots using the same frozen build; data/seeds/budgets stay unchanged.
 cases=cfg['qmc']['cases'];cases=sorted(cases,key=lambda c:((2 if c['proposal']=='local' else 1)*c['tau']/c['dt']),reverse=True)
 for ci,c in enumerate(cases):
  sweeps=40 if a.smoke else cfg['qmc']['local_sweeps_per_seed'] if c['proposal']=='local' else cfg['qmc']['reference_sweeps_per_seed'] if c['numerics']=='reference' else cfg['qmc']['guide_sweeps_per_seed']
  for si,seed0 in enumerate(cfg['qmc']['seeds'][:2] if a.smoke else cfg['qmc']['seeds']):
   seed=seed0+ci*100+(90000 if a.smoke else 0);idx=len(chains);M=round(c['tau']/c['dt']);assert M%2==0
   params={'Lx':2,'Ly':2,'shear':0,'Jx':1,'Jy':1,'Jz':.5,'mass':1,'omega-ph':1,'g':c['g'],'intervals':M,'dt':c['dt'],
    'trial-dx':0,'trial-dy':0,'trial-dz':0,'seed':seed,'warmup':20 if a.smoke else cfg['qmc']['warmup'],'sweeps':sweeps,'measure-every':1,'block-samples':20,
    'qr-block':4,'proposal':c['proposal'],'q-mode-every':1,'q-mode-n1':1,'q-mode-n2':1,'q-mode-beta':1,'q-mode-bonds':'z',
    'measure-energy':1,'fermions':0,'bond-diagnostics':1,'phonon-diagnostics':0,'vbs':1,'vbs-mode':'all','vbs-susceptibility':1,'vbs-amplitude':1,'vbs-amplitude-response':0,
    'vbs-field-h':0,'vbs-field-n1':1,'vbs-field-n2':1,'vbs-field-phase':0,'acceptance-every':200,'numerics':c['numerics']}
   if c['proposal']=='bond-guided':params.update({'bond-beta':1,'bond-guide':'bare'})
   label=f"g{c['g']}-tau{c['tau']}-dt{c['dt']}-{c['proposal']}-{c['numerics']}";args=[s for k,v in params.items() for s in ['--'+k,str(v)]]+['--output',f'runs/chain{idx:03d}']
   chains.append(dict(index=idx,case=label,case_parameters=c,seed=seed,parameters=params,args=args,tier='smoke' if a.smoke else 'full'))
   script='''#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export LD_LIBRARY_PATH="$PWD/runtime/lib:/opt/ohpc/pub/apps/intel/oneapi/mkl/2021.2.0/lib/intel64:${LD_LIBRARY_PATH:-}"
(cd runtime && sha256sum -c binary.sha256) > runs/chain%03d.integrity
set +e
./runtime/kitaev_ssh_mc %s > runs/chain%03d.stdout 2> runs/chain%03d.stderr
rc=$?
printf '%%s\\n' "$rc" > runs/chain%03d.exit
exit "$rc"
'''%(idx,shlex.join(args),idx,idx,idx)
   (r/'commands'/f'{idx:03d}.sh').write_text(script)
 (r/'plan.json').write_text(json.dumps(dict(config=cfg,chains=chains,tier='smoke' if a.smoke else 'full'),indent=2)+'\n')
 (r/'points.tsv').write_text('index\tcase\tseed\tg\ttau\tdt\tproposal\tnumerics\twarmup\tsweeps\n'+''.join('\t'.join(map(str,[c['index'],c['case'],c['seed'],c['case_parameters']['g'],c['case_parameters']['tau'],c['case_parameters']['dt'],c['case_parameters']['proposal'],c['case_parameters']['numerics'],c['parameters']['warmup'],c['parameters']['sweeps']]))+'\n' for c in chains))
 print(len(chains),'chains')
if __name__=='__main__':main()
