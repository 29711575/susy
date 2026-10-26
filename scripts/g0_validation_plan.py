#!/usr/bin/env python3
"""Freeze a finite, non-adaptive g=0 CLI validation campaign (no submission)."""
import argparse, hashlib, json, shlex
from pathlib import Path

def main():
    a=argparse.ArgumentParser();a.add_argument('directory',type=Path);x=a.parse_args();root=x.directory
    if (root/'plan.json').exists():raise FileExistsError('refusing to overwrite an existing campaign plan')
    root.mkdir(parents=True,exist_ok=True);(root/'commands').mkdir(exist_ok=True)
    # Runtime descending within q-off; q-on is submitted after q-off succeeds.
    cases=[
      ('local-coarse',2,.8,2,'local','production',0,1.,60000,'full'),
      # Supplementary same-seed prefix control; L2 reference retains full precision.
      ('complex-reference',4,.3,10,'bond-guided','reference',0,1.,200,'smoke'),
      ('complex-production',4,.3,10,'bond-guided','production',0,1.,10000,'full'),
      ('fine-reference',2,.3,20,'bond-guided','reference',0,1.,10000,'full'),
      ('fine-production',2,.3,20,'bond-guided','production',0,1.,10000,'full'),
      ('refresh-regression',2,.3,20,'bond-refresh','production',0,1.,10000,'full'),
      ('fine-local-smoke',2,.3,20,'local','production',0,1.,2000,'smoke'),
      ('coarse-production',2,.3,10,'bond-guided','production',0,1.,10000,'full'),
      ('omega08-fine',2,.8,20,'bond-guided','production',0,1.,10000,'full'),
      ('omega08-coarse',2,.8,2,'bond-guided','production',0,1.,10000,'full'),
      ('complex-q-on',4,.3,10,'bond-guided','production',1,1.,10000,'full'),
      ('fine-pcn-q-on',2,.3,20,'bond-guided','production',1,.6,10000,'full'),
    ]
    manifest={'scope':'g=0, fixed zero-flux Wilson(+1,+1), normalized identical free trials; not finite-g validation',
      'gates':{'sigma':6,'deterministic_absolute_floor':2e-10,'relative_se':.01,'B4_absolute_se':.05,
               'block_samples':20,'reblock_factors':[1,2,4,8,16], 'plateau_max_ratio':1.3,
               'plateau_min_ratio':.7,'minimum_final_blocks_per_chain':20,
               'nonzero_imaginary_gate_scale':'corresponding real reference',
               'decision_order':'FAIL for out-of-gate or numerical error; otherwise INCONCLUSIVE for precision/plateau/smoke; otherwise PASS'},
      'budget':'exact declared sweeps; no seed deletion, replacement, optional stopping, or tolerance adaptation',
      'cases':[],'chains':[]}
    for ci,(name,L,w,M,proposal,numerics,qon,beta,sweeps,tier) in enumerate(cases):
        case=dict(name=name,L=L,omega=w,M=M,dt=2./M,proposal=proposal,numerics=numerics,qon=qon,beta=beta,sweeps=sweeps,tier=tier)
        manifest['cases'].append(case)
        for si in range(4):
            idx=len(manifest['chains']);seed=202610170000+ci*100+si
            # Use matched seeds for reference/production comparisons.
            if name=='complex-production':seed=202610170100+si
            if name=='fine-production':seed=202610170300+si
            params={'Lx':L,'Ly':L,'shear':0,'Jx':1,'Jy':1,'Jz':.5,'g':0,'mass':100./9,'omega-ph':w,
              'intervals':M,'dt':2./M,'trial-dx':0,'trial-dy':0,'trial-dz':0,'warmup':200 if proposal!='local' else 2000,
              'sweeps':sweeps,'measure-every':1,'block-samples':20,'qr-block':4,'seed':seed,
              'proposal':proposal,'q-mode-every':qon,'q-mode-n1':1,'q-mode-n2':1 if L==2 else 0,'q-mode-beta':beta,
              'q-mode-bonds':'z','measure-energy':1,'fermions':0,'vbs':1,'vbs-mode':'all','vbs-susceptibility':1,
              'vbs-amplitude':1,'vbs-amplitude-response':0,'phonon-diagnostics':1,
              'acceptance-every':200,'vbs-field-h':0,'vbs-field-n1':0,'vbs-field-n2':0,'vbs-field-phase':0}
            if numerics=='reference':params['numerics']='reference'
            # For production deliberately exercise omitted numerical options.
            if proposal=='bond-guided':params.update({'bond-guide':'zero','bond-guide-x':0,'bond-guide-y':0,'bond-guide-z':0,'bond-beta':beta})
            args=[v for k,val in params.items() for v in ['--'+k,str(val)]]+['--output','runs/chain%03d'%idx]
            manifest['chains'].append(dict(index=idx,case=name,seed=seed,parameters=params,args=args))
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
            (root/'commands'/('%03d.sh'%idx)).write_text(script)
    (root/'plan.json').write_text(json.dumps(manifest,indent=2)+'\n')
    rows=['index\tcase\tseed\tL\tomega\tM\tdt\tproposal\tnumerics\tqon\tbeta\tsweeps\ttier']
    for ch in manifest['chains']:
        c=next(c for c in manifest['cases'] if c['name']==ch['case'])
        rows.append('\t'.join(map(str,[ch['index'],c['name'],ch['seed']]+[c[k] for k in ['L','omega','M','dt','proposal','numerics','qon','beta','sweeps','tier']])))
    (root/'points.tsv').write_text('\n'.join(rows)+'\n')
    print(len(manifest['chains']),'chains; qoff 0-39, qon 40-47')
if __name__=='__main__':main()
