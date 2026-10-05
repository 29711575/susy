#!/usr/bin/env python3
"""Read-only audit of the completed tau scan and its bare Dirac reference."""
import csv, hashlib, json, os, subprocess
from pathlib import Path
from statistics import NormalDist
import numpy as np
from convergence_analysis import read_chain, joint_estimate, estimates_json, fixed_momentum_mapping
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/phase-tau-20260917'
OUT=ROOT/'results/phase-tau-20260917/peak_audit'
OUT.mkdir(exist_ok=True)
manifest=json.loads((RUN/'manifest.json').read_text())
chains=[read_chain(Path(c['folder'])/'run_results.json') for c in manifest['commands']]
checks={}
def check(label,error,tolerance=2e-11):
    checks[label]=max(checks.get(label,0.),float(error))
    if not np.isfinite(error) or error>tolerance:
        raise AssertionError((label,error,tolerance))
hashes={name:hashlib.sha256((ROOT/'kitaev_ssh'/name).read_bytes()).hexdigest()==expected
        for name,expected in manifest['source_hashes'].items()}
for name in ['vbs_observables.h','phonon_mc.h','phonon_projector.h','phonon_path.h']:
    assert hashes[name], name
assert hashlib.sha256((RUN/'kitaev_ssh_mc').read_bytes()).hexdigest()==manifest['binary_sha256']
# Identities on EVERY saved block, before any pooling.
for c in chains:
    assert len(c.counts)==32 and np.all(c.counts==25)
    check('JSON_block_ratio_max_error',c.csv_ratio_error)
    r={name:c.sums[:,i+1]/c.sums[:,0] for i,name in enumerate(c.names)}
    meta=c.metadata
    nc=meta['Lx']*meta['Ly']
    assert nc==16 and fixed_momentum_mapping(meta)['q3']['index']==10
    for mom in meta['vbs']['momenta']:
        q=mom['index'];nx,ny=mom['n']
        assert q==4*nx+ny
        expect=np.array([nx/4,ny/4])@np.array(meta['vbs']['reciprocal'])
        check('momentum_cartesian_mapping',np.max(abs(expect-mom['cartesian'])))
        for mu in 'xyz':
            check('diagonal_imaginary',np.max(abs(r[f'S_q{q}_{mu}{mu}'].imag)))
            partner=((4-nx)%4)*4+(4-ny)%4
            check('S_diagonal_plus_minus_q',np.max(abs(r[f'S_q{q}_{mu}{mu}']-r[f'S_q{partner}_{mu}{mu}'])))
            for nu in 'xyz':
                check('S_hermitian',np.max(abs(r[f'S_q{q}_{mu}{nu}']-r[f'S_q{q}_{nu}{mu}'].conj())))
    for mu in 'xyz':
        check('Parseval_per_channel',np.max(abs(sum(r[f'S_q{q}_{mu}{mu}'] for q in range(nc))/nc-r[f'X2_{mu}'])))
        check('uniform_X_equals_Delta0',np.max(abs(r[f'X_{mu}']-r[f'Delta_q0_{mu}'])))
        check('block_variance_nonnegative',max(0.,float(np.max(nc*r[f'Delta_q0_{mu}'].real**2-r[f'S_q0_{mu}{mu}'].real))))
# Bare Bloch dispersion versus an independent real-space matrix construction.
theta=np.arccos(-.25)
K=np.array([theta,-theta])/(2*np.pi)
Q=(2*K)%1
check('Dirac_zero',abs(.5+np.exp(-1j*theta)+np.exp(1j*theta)))
A=np.zeros((32,32))
for x in range(4):
    for y in range(4):
        a=x*4+y
        for dx,dy,j in [(-1,0,1),(0,-1,1),(0,0,.5)]:
            b=16+((x+dx)%4)*4+(y+dy)%4
            A[a,b]+=2*j;A[b,a]-=2*j
grid=np.array([(x/4,y/4) for x in range(4) for y in range(4)])
absf=np.abs(.5+np.exp(-2j*np.pi*grid[:,0])+np.exp(-2j*np.pi*grid[:,1]))
check('bare_realspace_Bloch_spectrum',np.max(abs(np.linalg.eigvalsh(1j*A)-np.sort(np.r_[2*absf,-2*absf]))))
B=np.array(chains[0].metadata['vbs']['reciprocal'])
def distance(a,b):
    return min(np.linalg.norm((a-b+np.array([i,j]))@B) for i in range(-2,3) for j in range(-2,3))
distance_Q=np.array([min(distance(q,Q),distance(q,(-Q)%1)) for q in grid])
nearest=np.flatnonzero(np.isclose(distance_Q,distance_Q.min())).tolist()
assert nearest==[10]
# Full complex phase ratio; derived contrasts share every jackknife replica.
def transform(sums,count,names,meta):
    r=dict(zip(names,sums[1:]/sums[0]));nc=16
    s=[r[f'S_q{q}_zz'] for q in range(nc)]
    d0=r['Delta_q0_z']
    bg=nc*d0*d0  # X is real; preserve the original complex ratio algebra.
    other=sum(s[q] for q in range(1,16) if q!=10)/14
    result={'S0':s[0],'uniform_background':bg,'S0_connected':s[0]-bg,
            'X_z':d0,'X2_z':r['X2_z'],'S_Qnearest':s[10],
            'nonzero_mean':sum(s[1:])/15,'Qnearest_minus_other_mean':s[10]-other}
    for q in range(16):
        result[f'S_q{q}']=s[q]
        result[f'S_q{q}_minus_Qnearest']=s[q]-s[10]
        result[f'S_q{q}_minus_X2']=s[q]-r['X2_z']
    return result
rows=[];contrasts=[];peaks=[]
for m in manifest['plan']['intervals']:
    cc=[c for c in chains if c.metadata['intervals']==m]
    tau=round(m*.1,12)
    for block in (25,50,100):
        est=estimates_json(joint_estimate([c.reblock(block) for c in cc],transform))
        for name,e in est.items():
            rows.append(dict(tau=tau,block_sweeps=block,observable=name,real=e['value'][0],imag=e['value'][1],SE=e['stderr_real']))
        for q in range(1,16):
            if q==10:continue
            e=est[f'S_q{q}_minus_Qnearest'];v=e['value'][0];se=e['stderr_real']
            contrasts.append(dict(tau=tau,block_sweeps=block,q=q,nx=q//4,ny=q%4,difference=v,SE=se,z=v/se))
        if block==25:
            s=np.array([est[f'S_q{q}']['value'][0] for q in range(16)])
            nz=1+int(np.argmax(s[1:]))
            equal=[q for q in range(1,16) if abs(s[q]-s[nz])<1e-12]
            peak=dict(tau=tau,global_peak=int(np.argmax(s)),nonzero_peak_indices=';'.join(map(str,equal)),
                      S0=est['S0']['value'][0],uniform_background=est['uniform_background']['value'][0],
                      S0_connected=est['S0_connected']['value'][0],S0_connected_SE=est['S0_connected']['stderr_real'],
                      S_Qnearest=est['S_Qnearest']['value'][0],S_Qnearest_SE=est['S_Qnearest']['stderr_real'],
                      X2_z=est['X2_z']['value'][0],nonzero_mean=est['nonzero_mean']['value'][0],
                      Qnearest_minus_other_mean=est['Qnearest_minus_other_mean']['value'][0],
                      Qnearest_minus_other_mean_SE=est['Qnearest_minus_other_mean']['stderr_real'])
            peaks.append(peak)
def save(name,data):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
save('derived_blocking.csv',rows);save('fixed_q_contrasts.csv',contrasts);save('peak_decomposition.csv',peaks)
# Existing validation executables only: no new MC production chains.
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    env[k]='1'
env['MKL_DYNAMIC']='FALSE'
testresults={}
for name in ['kitaev_ssh_validate','kitaev_ssh_phonon_validate','kitaev_ssh_mc_validate','kitaev_ssh_vbs_validate','kitaev_ssh_weight_validate']:
    p=subprocess.run([str(ROOT/'build'/name)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=300)
    (OUT/(name+'.log')).write_text(p.stdout+p.stderr)
    testresults[name]=p.returncode
    print('VALIDATION',name,p.returncode,flush=True)
    assert p.returncode==0, (name,p.stderr)
critical=NormalDist().inv_cdf(1-.05/(2*6*14))
maximum=max(abs(x['z']) for x in contrasts if x['block_sweeps']==25)
audit=dict(check_max_errors=checks,scan_source_hash_matches=hashes,scan_binary_hash_match=True,
           current_validation_exit_codes=testresults,bare_K_plus_fractions=K.tolist(),bare_Q_fractions=Q.tolist(),
           nearest_Q_grid_indices=nearest,bare_finite_grid_min_positive_energy=float(2*absf.min()),
           fixed_contrasts_84_tests_bonferroni_z=float(critical),largest_abs_z_at_block25=maximum,
           scope='Existing saved block audit plus existing deterministic/test validation; no production rerun.',
           caveats=['Source coupled_mc.cpp has later fermion changes; audit uses saved scan outputs and preserved binary hash.',
                    'Error bars conditional on adequate blocking/stationarity; max peak selection is exploratory.',
                    'No stored individual configurations: no independent Fourier recomputation on original MC configurations.',
                    'A q=0 uniform displacement is allowed and is not finite-q VBS order.',
                    'No theorem requires a weak-coupling equal-time structure-factor maximum at the bare Dirac difference.'])
(OUT/'audit_summary.json').write_text(json.dumps(audit,indent=2))
lines=['# Structure-factor peak audit','',
'No measurement/coordinate/normalization failure found in the tested scope. All five existing validation executables passed.',
'Checked 12 chains / 384 saved blocks. No MC production rerun. Existing Hamiltonian, Green, proposal and RNG paths unchanged by this audit.',
'Fourier convention: F=sum exp(-iqR)X; S=<F F*>/Nc. Parseval checked independently per saved block and channel.',
'Bare f=.5+exp(-ik1)+exp(-ik2); k_i=k dot a_i. Q fractions=(.58043062,.41956938); nearest L4 grid=(.5,.5), flattened q=10.',
'Uniform contribution: Szz(0)=Nc*<Delta_z(0)>^2+connected Szz(0); only q=0 is altered by subtracting the uniform ensemble mean.',
'Bare phonon baseline Szz(q)=1/(2*m*omega)=.25; this is a reference, not a predicted exact value at g=.5.','',
'| tau | Szz(0) | uniform contribution | connected Szz(0) | nearest-Q Szz | nonzero mean |',
'|---|---|---|---|---|---|']
for p in peaks:
    lines.append(f"| {p['tau']:g} | {p['S0']:.6f} | {p['uniform_background']:.6f} | {p['S0_connected']:.6f} +- {p['S0_connected_SE']:.6f} | {p['S_Qnearest']:.6f} +- {p['S_Qnearest_SE']:.6f} | {p['nonzero_mean']:.6f} |")
lines+=['',f'Across 84 predefined contrasts (14 nonzero q against q=10, six tau), max |z|={maximum:.4f} at block25; nominal Bonferroni 95% threshold={critical:.4f}.',
'These normal-approximation screens are conditional on adequate blocking and stationarity, not proof of no finite-q response. 50/100 sweep sensitivity is retained.',
'No reconstructed or reweighted interacting Dirac point is substituted for the requested bare reference.',
'No original per-configuration fields were saved, so this audit cannot independently recompute Fourier transforms for the original 9600 configurations.',
'The previous all-q maximum plot included the uniform q=0 background; it is not a verified VBS order parameter.',
'The previous claim that the maximum must match the bare valley transfer was too strong. A valley-coupling instability is a hypothesis, not a constraint on a raw structure factor at every coupling.',
'For measurement definitions, code still uses all three independently fluctuating bond types with J_mu+g X_mu. No subtraction of the matter bond mean or constraint sum_R X_R=0 is present. Such a change would change the model and is not part of this audit.',
'Current test logs are in this directory; original scan validation remains in runs/phase-tau-20260917/validation.']
(OUT/'PEAK_AUDIT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(audit,indent=2),flush=True)
print('\n'.join(lines[8:17]),flush=True)

