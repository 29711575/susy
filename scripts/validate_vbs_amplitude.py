#!/usr/bin/env python3
"""Native fixed-seed compatibility and amplitude integration checks (<5 min)."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'kitaev_ssh'))
from analyze_vbs_amplitude import amplitude_chain, run_analysis


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True,type=Path)
    p.add_argument('--baseline',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    env.update(LD_LIBRARY_PATH='/opt/intel/oneapi/mkl/2026.1/lib:'+env.get('LD_LIBRARY_PATH',''),
               OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',BLIS_NUM_THREADS='1',MKL_DYNAMIC='FALSE')
    base=['--Lx','3','--Ly','2','--Jx','1','--Jy','1','--Jz','.5',
          '--mass','4','--omega-ph','.5','--g','.3','--intervals','4','--dt','.17',
          '--warmup','4','--sweeps','32','--measure-every','1','--block-samples','4',
          '--vbs','1','--vbs-mode','all','--weight-engine','combined','--defer-phase','1',
          '--prefix-cache','1','--cayley-cache','1','--proposal','bond-guided','--bond-guide','bare',
          '--bond-beta','1','--q-mode-every','1','--q-mode-beta','1','--q-mode-n1','1','--q-mode-n2','0']
    def override(options,updates):
        selected=dict(zip(options[::2],options[1::2]));selected.update(updates)
        return [value for pair in selected.items() for value in pair]
    records=[]
    def run(label,binary,options,extra=()):
        prefix=args.output/label
        cmd=[str(binary.resolve())]+options+list(extra)+['--output',str(prefix.resolve())]
        r=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
        prefix.with_suffix('.log').write_text(r.stdout)
        records.append(dict(label=label,argv=cmd,exit=r.returncode))
        if r.returncode:
            raise RuntimeError(f'{label} failed: {r.stdout[-2000:]}')
        meta=json.loads(Path(str(prefix)+'_results.json').read_text())
        with Path(str(prefix)+'_blocks.csv').open() as f:
            rows=list(csv.DictReader(f))
        return meta,rows,prefix
    for case,(shear,chi,response) in enumerate([(0,0,0),(1,0,1),(0,1,1)]):
        options=base+['--seed',str(52010+case),'--shear',str(shear),'--vbs-susceptibility',str(chi)]
        if chi:
            options+=['--vbs-field-h','0','--vbs-field-n1','1','--vbs-field-n2','0','--vbs-field-phase','.31']
        old,oldrows,_=run(f'case{case}-baseline',args.baseline,options)
        disabled,disabledrows,_=run(f'case{case}-disabled',args.binary,options)
        enabled,enabledrows,prefix=run(f'case{case}-enabled',args.binary,options,
                                     ['--vbs-amplitude','1','--vbs-amplitude-response',str(response)])
        assert oldrows==disabledrows,'disabled interface altered legacy raw blocks'
        assert len(enabledrows)==len(oldrows)
        for a,b in zip(oldrows,enabledrows):
            assert all(b[k]==v for k,v in a.items()),'measurement altered existing phase/raw observable columns'
        assert old['observables']==disabled['observables']
        assert all(enabled['observables'][k]==v for k,v in old['observables'].items())
        assert all(a['derived']==b['derived'] for a,b in zip(old['vbs']['blocking'],enabled['vbs']['blocking']))
        if chi:
            for row in enabledrows:
                for name,field in [('AmpA2','Field_absD2'),('AmpA4','Field_absD4')]:
                    for part in ['real','imag']:
                        a=float(row[f'sum_phase_{name}_q2_z_{part}']);b=float(row[f'sum_phase_{field}_{part}'])
                        assert abs(a-b)<2e-13*max(1,abs(a),abs(b)),'fixed-q field moment mismatch'
        c=amplitude_chain(Path(str(prefix)+'_results.json'),1,0,'z')
        analysis=run_analysis([c],1,0,[3],2,100,843)
        Path(str(prefix)+'_amplitude-analysis.json').write_text(json.dumps(analysis,indent=2,allow_nan=False)+'\n')
        assert analysis['points'][0]['quality_flags'],'short smoke must remain statistically unqualified'
    # L6 exercises production dimensions with new moments and old static
    # susceptibility simultaneously. This is an interface smoke, not a scan.
    opts=override(base,{'--Lx':'6','--Ly':'6','--seed':'621003','--warmup':'1','--sweeps':'2','--block-samples':'1',
                       '--q-mode-n1':'3','--q-mode-n2':'3','--vbs-susceptibility':'1'})
    run('L6-smoke',args.binary,opts,['--vbs-amplitude','1','--vbs-amplitude-response','1'])
    for label,extra in [('no-vbs',['--vbs','0','--vbs-amplitude','1']),
                        ('response-only',['--vbs-amplitude-response','1'])]:
        options=override(base,dict(zip(extra[::2],extra[1::2])))
        r=subprocess.run([str(args.binary.resolve())]+options,env=env,capture_output=True,text=True,timeout=30)
        assert r.returncode!=0,'invalid switch combination accepted'
        assert 'requires --vbs' in r.stderr,'rejection did not enforce the intended switch dependency'
        records.append(dict(label=label,exit=r.returncode,expected_rejection=True,error=r.stderr))
    out=dict(status='PASS',binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
             baseline_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
             cases=records,legacy_common_columns_bit_identical=True,physics_convergence_certified=False)
    (args.output/'native-integration.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: three fixed-seed native legacy/disabled/enabled comparisons, field-moment agreement, L6 smoke and invalid CLI checks')


if __name__=='__main__':
    main()
