#!/usr/bin/env python3
"""Deterministic tests only. Fixtures use the real C++ output schema, not QMC data."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from export_vbs_momentum import (AXES, export, jackknife_se, make_path, matrices,
    ratio_replicas, read_run, scalar_channels, targets_for_run, window_weights,
    evaluate_fixed_associations, fixed_vbs_geometry, fixed_vbs_groups)

BASIS = np.array([[2*np.pi/np.sqrt(3), 2*np.pi/3],[-2*np.pi/np.sqrt(3), 2*np.pi/3]])
Q = np.array([.5804306232551663,.4195693767448337])


def geometry(nx, ny, shear):
    grid = np.array([[n1/nx, n2/ny-shear*n1/(nx*ny)] for n1 in range(nx) for n2 in range(ny)])
    coords = np.array([[x,y] for x in range(nx) for y in range(ny)])
    return grid, coords


def fixture(root, nblocks=4, block_samples=2, phases=None):
    """Create a labeled synthetic fixture in the current coupled_mc schema."""
    nx=ny=4; nc=nx*ny
    grid, coords = geometry(nx,ny,0)
    rng=np.random.default_rng(271828)
    fields=rng.normal(size=(nblocks,block_samples,nc,3)) + np.array([.2,.3,.4])
    f=np.einsum('kr,bsri->bski',np.exp(-2j*np.pi*(grid@coords.T)),fields)
    s=f[..., :,None]*f[...,None,:].conj()/nc
    d0=f[:,:,0]/nc
    if phases is None:
        phases=np.exp(1j*rng.uniform(-.15,.15,size=(nblocks,block_samples)))
    phase=phases.sum(axis=1)
    snum=np.einsum('bs,bskij->bkij',phases,s)
    dnum=np.einsum('bs,bsi->bi',phases,d0)
    meta=dict(status="complete",vbs_enabled=True,spatial_boundary="pbc",Lx=nx,Ly=ny,shear=0,
        Jx=1,Jy=1,Jz=.5,g=.3,mass=1,omega_ph=1,k=1,dt=.1,intervals=2,tau=.2,seed=0,
        samples=nblocks*block_samples,blocks=nblocks,block_samples=block_samples,
        vbs=dict(mode="all",Nc=nc,reciprocal=BASIS.tolist(),momenta=[dict(index=i,n=[i//ny,i%ny],
            reciprocal_fraction=q.tolist(),cartesian=(q@BASIS).tolist()) for i,q in enumerate(grid)]))
    results=Path(root)/'synthetic_ONLY_results.json'
    blocks=Path(root)/'synthetic_ONLY_blocks.csv'
    results.write_text(json.dumps(meta))
    names=[f'S_q{k}_{mu}{nu}' for k in range(nc) for mu in AXES for nu in AXES]
    names += [f'Delta_q0_{mu}' for mu in AXES]
    columns=['block','count','sum_phase_real','sum_phase_imag']
    columns += [f'sum_phase_{name}_{part}' for name in names for part in ('real','imag')]
    with blocks.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=columns); writer.writeheader()
        for b in range(nblocks):
            row=dict(block=b,count=block_samples,sum_phase_real=phase[b].real,sum_phase_imag=phase[b].imag)
            for name,z in zip(names,np.r_[snum[b].reshape(-1),dnum[b]]):
                row[f'sum_phase_{name}_real']=z.real; row[f'sum_phase_{name}_imag']=z.imag
            writer.writerow(row)
    return results,blocks,meta,phase,snum,dnum


class ExportMomentumTests(unittest.TestCase):
    def test_fixed_vbs_steps_and_inversion_keep_exact_fractional_coordinates(self):
        meta=dict(status="complete", vbs_enabled=True, spatial_boundary="pbc",
                  Lx=6, Ly=3, shear=0, vbs=dict(pbc_geometry="axis_aligned"),
                  Jx=1, Jy=1, Jz=.5)
        d1,d2,source=fixed_vbs_geometry(meta)
        np.testing.assert_allclose(d1,[1/6,0]); np.testing.assert_allclose(d2,[0,1/3])
        self.assertIn("b1_over_Lx",source)
        groups=fixed_vbs_groups(meta, explicit_q=[.37,.21])
        np.testing.assert_allclose(groups[1]["center"],[-.37 % 1,-.21 % 1])
        np.testing.assert_allclose(groups[0]["points"]["plus_b1"],[(.37+1/6)%1,.21])

    def test_fixed_association_flat_connected_spectrum_has_zero_ratios(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,p,s,d=fixture(root,nblocks=8)
            run=read_run(results)
            run["snum"][:] = 0
            run["snum"][:] = np.eye(3)[None,None,:,:] * 2.0
            run["dnum"][:] = 0
            rows,_=evaluate_fixed_associations(run, explicit_q=[.37,.21], min_blocks=2)
            for row in rows:
                self.assertTrue(row["R_b1_available"])
                self.assertTrue(row["R_b2_available"])
                self.assertAlmostEqual(float(row["R_b1"]),0.0,places=13)
                self.assertAlmostEqual(float(row["R_b2"]),0.0,places=13)

    def test_fixed_association_marks_zero_mode_group_unavailable(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,p,s,d=fixture(root,nblocks=8)
            run=read_run(results)
            rows,_=evaluate_fixed_associations(run, explicit_q=[0,0], min_blocks=2)
            self.assertEqual(len(rows),2)
            self.assertTrue(all(not row["R_b1_available"] and
                                row["R_b1_unavailable_reason"] == "q0_equivalent_in_group"
                                for row in rows))
    def test_weights_match_explicit_sum_square_rectangle_and_shear(self):
        for nx,ny,s in [(6,6,0),(12,3,0),(12,3,3),(5,4,2)]:
            grid,coords=geometry(nx,ny,s)
            for q in (Q, [.137,.271],grid[3]):
                w,index=window_weights(q,grid,nx,ny)
                direct=abs(np.exp(-2j*np.pi*((np.asarray(q)-grid)@coords.T)).mean(axis=1))**2
                np.testing.assert_allclose(w,direct,atol=2e-14,rtol=0)
                self.assertAlmostEqual(w.sum(),1.,places=13)
            self.assertEqual(index,3)

    def test_zero_window_does_not_imply_allowed_momentum(self):
        grid,_=geometry(12,3,0)
        w,index=window_weights([5/12,7/12],grid,12,3)
        self.assertIsNone(index)
        self.assertLess(w[0],1e-28)

    def test_uniform_constant_is_removed_in_both_subtractions(self):
        grid,_=geometry(6,6,0); nc=36
        mean=np.array([.25,-.7,1.2],complex)
        s=np.zeros((3,nc,3,3),complex)
        s[:,0]=nc*np.outer(mean,mean.conj())
        d=np.tile(mean,(3,1))
        w,_=window_weights(Q,grid,6,6)
        out=matrices(s,d,w,nc)
        self.assertGreater(np.linalg.norm(out[0,0]),.01)
        np.testing.assert_allclose(out[:,1:],0,atol=2e-14,rtol=0)

    def test_connected_keeps_uniform_fluctuation_projected_removes_it(self):
        grid,_=geometry(6,6,0); nc=36
        d=np.array([[.2,.3,.4]],complex)
        covariance=np.diag([.1,.2,.3])
        s=np.zeros((1,nc,3,3),complex)
        s[0,0]=nc*(np.outer(d[0],d[0].conj())+covariance)
        w,_=window_weights(Q,grid,6,6)
        out=matrices(s,d,w,nc)[0]
        np.testing.assert_allclose(out[1]-out[2],w[0]*nc*covariance,atol=2e-14,rtol=0)
        np.testing.assert_allclose(out[2],0,atol=2e-14,rtol=0)

    def test_translation_orbit_identity_including_sheared_quotient(self):
        for nx,ny,shear in [(4,4,0),(4,3,1)]:
            grid,coords=geometry(nx,ny,shear); nc=len(grid)
            field=np.random.default_rng(190).normal(size=(nc,3))
            f=np.exp(-2j*np.pi*(grid@coords.T))@field
            sg=f[:,:,None]*f[:,None,:].conj()/nc
            w,_=window_weights(Q,grid,nx,ny)
            expected=matrices(sg[None],(f[0]/nc)[None],w,nc)[0]
            explicit=np.zeros((3,3),complex); explicit_prj=explicit.copy()
            phase=np.exp(-2j*np.pi*(coords@Q)); W=phase.mean()
            for rx,ry in coords:
                yy=coords[:,1]-ry; wy=np.floor_divide(yy,ny); yy-=wy*ny
                xx=(coords[:,0]-rx-shear*wy)%nx
                shifted=field[xx*ny+yy]
                fq=phase@shifted
                fp=fq-W*f[0]
                explicit+=np.outer(fq,fq.conj())/(nc*nc)
                explicit_prj+=np.outer(fp,fp.conj())/(nc*nc)
            np.testing.assert_allclose(expected[0],explicit,atol=2e-13,rtol=0)
            np.testing.assert_allclose(expected[2],explicit_prj,atol=2e-13,rtol=0)

    def test_nonzero_grid_is_unchanged_and_eigenvalues_follow_matrix_average(self):
        grid,_=geometry(2,2,0); nc=4
        s=np.zeros((1,nc,3,3),complex)
        s[0,1]=np.diag([4,0,0]); s[0,2]=np.diag([0,4,0])
        d=np.array([[100,200,300]],complex)
        w,index=window_weights(grid[1],grid,2,2)
        out=matrices(s,d,w,nc)
        for e in range(3): np.testing.assert_allclose(out[0,e],s[0,1])
        w=np.array([0,.5,.5,0])
        mixed=scalar_channels(matrices(s,d,w,nc))[0,0,2]
        self.assertEqual(mixed,2.)
        self.assertNotEqual(mixed,4.)  # average of the two leading eigenvalues is wrong

    def test_joint_complex_reweighting_and_jackknife(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,_,p,s,d=fixture(root)
            run=read_run(results)
            sr,dr,diag=ratio_replicas(run)
            np.testing.assert_allclose(sr[0],s.sum(axis=0)/p.sum())
            w,_=window_weights(Q,run['grid'],4,4)
            actual=scalar_channels(matrices(sr,dr,w,16))
            manual=[]
            for i in range(len(p)):
                pi=p.sum()-p[i]
                sg=(s.sum(axis=0)-s[i])/pi
                bg=(d.sum(axis=0)-d[i])/pi
                raw=np.einsum('k,kij->ij',w,sg)
                con=raw-16*w[0]*np.outer(bg,bg.conj())
                prj=raw-w[0]*sg[0]
                manual.append(scalar_channels(np.stack([raw,con,prj])))
            np.testing.assert_allclose(actual[1:],manual,atol=2e-13,rtol=0)
            np.testing.assert_allclose(jackknife_se(actual[1:]),jackknife_se(manual),atol=2e-13,rtol=0)
            self.assertEqual(diag['blocks'],4)
            self.assertFalse(diag['block_independence_verified'])

    def test_bare_Q_is_exactly_present_and_path_uses_physical_distance(self):
        meta=dict(Jx=1,Jy=1,Jz=.5,vbs={})
        targets,nodes=targets_for_run(meta)
        self.assertEqual(nodes['status'],'ok')
        np.testing.assert_allclose(targets[0]['q'],Q,atol=2e-15,rtol=0)
        grid,_=geometry(6,6,0)
        points,ticks=make_path([[0,1],[1,0]],BASIS,17,grid,targets)
        self.assertTrue(any(np.linalg.norm(p['q']-Q)<2e-15 for p in points))
        self.assertAlmostEqual(ticks[-1],np.linalg.norm(np.array([1,-1])@BASIS))
        self.assertTrue(all(points[i]['distance']<points[i+1]['distance'] for i in range(len(points)-1)))

    def test_schema_end_to_end_and_matrix_order(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,*_=fixture(root)
            prefix=Path(root)/'export'
            report=export(results,None,prefix,[[0,1],[1,0]],samples=17)
            self.assertEqual(report['rows']['grid'],16)
            self.assertEqual(report['rows']['targets'],2)
            with Path(str(prefix)+'_grid.csv').open() as f: rows=list(csv.DictReader(f))
            self.assertEqual(sum(int(r['is_q0']) for r in rows),1)
            for r in rows[1:]:
                self.assertAlmostEqual(float(r['raw_szz']),float(r['connected_szz']))
                self.assertAlmostEqual(float(r['connected_szz_per_cell'])*16,float(r['connected_szz']))
            data=np.load(str(prefix)+'_matrices.npz',allow_pickle=False)
            self.assertEqual(data['matrix_center'].shape,(sum(report['rows'].values()),3,3,3))
            self.assertEqual(list(data['estimator']),['raw','connected','projected'])
            self.assertTrue(np.all(np.isfinite(data['matrix_center'])))
            with self.assertRaises(FileExistsError):
                export(results,None,prefix,[[0,1],[1,0]],samples=17)

    def test_reblock_matches_manual_sums_and_rejects_remainder(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,p,s,d=fixture(root,nblocks=6)
            run=read_run(results,reblock=2)
            np.testing.assert_allclose(run['phase'],p.reshape(3,2).sum(axis=1))
            np.testing.assert_allclose(run['snum'],s.reshape(3,2,16,3,3).sum(axis=1))
            self.assertEqual(run['block_count'],4)
            with self.assertRaises(ValueError): read_run(results,reblock=4)

    def test_rejects_incomplete_and_selected_only_results(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,*_=fixture(root)
            for status in ('aborted','running'):
                meta['status']=status; results.write_text(json.dumps(meta))
                with self.assertRaises(ValueError): read_run(results)
            meta['status']='complete'; meta['vbs']['mode']='selected'; results.write_text(json.dumps(meta))
            with self.assertRaises(ValueError): read_run(results)

    def test_rejects_mismatched_counts_and_nonfinite_data(self):
        with tempfile.TemporaryDirectory() as root:
            results,blocks,meta,*_=fixture(root)
            meta['samples']+=1; results.write_text(json.dumps(meta))
            with self.assertRaises(ValueError): read_run(results)
            meta['samples']-=1; results.write_text(json.dumps(meta))
            lines=blocks.read_text().splitlines(); fields=lines[1].split(','); fields[4]='nan'
            lines[1]=','.join(fields); blocks.write_text('\n'.join(lines)+'\n')
            with self.assertRaises(ValueError): read_run(results)

    def test_zero_phase_denominator_is_not_silently_regularized(self):
        with tempfile.TemporaryDirectory() as root:
            phase=np.array([[1,1],[1,1],[-1,-1],[-1,-1]],complex)
            results,*_=fixture(root,phases=phase)
            with self.assertRaises(ValueError): ratio_replicas(read_run(results))

    def test_plots_smoke_only(self):
        with tempfile.TemporaryDirectory() as root:
            results,*_=fixture(root)
            prefix=Path(root)/'synthetic_not_QMC'
            export(results,None,prefix,[[0,1],[1,0]],samples=17,plot=True,per_cell=True)
            for ending in ('_path.png','_grid.png'):
                self.assertGreater(Path(str(prefix)+ending).stat().st_size,5000)


if __name__=='__main__':
    unittest.main()
