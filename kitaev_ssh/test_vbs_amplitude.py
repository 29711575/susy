import unittest
from pathlib import Path
from dataclasses import replace
import numpy as np
from convergence_analysis import Chain, joint_estimate
from analyze_vbs_amplitude import (amplitude_observables, scan_groups, bootstrap_curves,
                                  fit_slopes, peak_scan, run_analysis, mode)


def chain(seed=1, lam=.3, blocks=20, factor=1.):
    meta=dict(Lx=3,Ly=2,shear=0,geometry_id='honeycomb-pbc-Lx3-Ly2-shear0',
              Jx=1.,Jy=1.,Jz=.5,g=np.sqrt(lam),mass=4.,omega_ph=.5,
              scope='fixed_gauge_matter_sector',spatial_boundary='pbc',wilson_loops=[1,1],
              time_contour='open',intervals=4,dt=.2,seed=seed,measure_every=1,
              trial_shift=[0.,0.,0.],vbs=dict(Nc=6,momenta=[dict(index=2*x+y,n=[x,y],
                  reciprocal_fraction=[x/3,y/2]) for x in range(3) for y in range(2)]))
    a = factor*(.2+.7*np.arange(blocks)/(blocks-1))
    ia = .4+.3*a
    rows = np.column_stack([np.ones(blocks),a,a*a,a**4,ia,a*ia])*2
    return Chain(Path(f'fake_{seed}_{lam}'),meta,[f'{s}_q2_z' for s in ['AmpA','AmpA2','AmpA4','AmpIA','AmpAIA']],
                 rows.astype(complex),np.full(blocks,2))


class AmplitudeTests(unittest.TestCase):
    def test_complex_phase_ratios_use_no_conjugated_means(self):
        a=np.array([.2,.7,1.1]);ia=np.array([.4,.9,.6]);w=np.array([1+.2j,-.13+.1j,.7-.3j])
        r=np.array([w.sum(),w@a,w@(a*a),w@(a**4),w@ia,w@(a*ia)])
        out=amplitude_observables(r,3,[],dict(Lx=3,Ly=2))
        self.assertAlmostEqual(abs(out['VA']-6*(w@(a*a)/w.sum()-(w@a/w.sum())**2)),0,places=13)
        self.assertAlmostEqual(abs(out['chiA']-6*(w@(a*ia)/w.sum()-(w@a/w.sum())*(w@ia/w.sum()))),0,places=13)
        self.assertGreater(abs(out['VA']-6*(w@(a*a)/w.sum()-abs(w@a/w.sum())**2)),.001)

    def test_true_fourth_and_connected_fluctuation(self):
        c=chain();r=amplitude_observables(c.sums.sum(axis=0),40,c.names,c.metadata)
        self.assertGreater(abs(r['A4']-r['A2']**2),.01)
        self.assertAlmostEqual(r['chiA'].real,.3*r['VA'].real,places=13)

    def test_joint_jackknife_unequal_chain_allocation(self):
        chains=[chain(1,blocks=20),chain(2,blocks=12,factor=1.3)]
        e=joint_estimate(chains,amplitude_observables)
        cov=np.zeros_like(e['covariance']);total=sum(c.sums.sum(axis=0) for c in chains)
        for c in chains:
            b=len(c.counts);values=[]
            for row in c.sums:
                replica=total+(c.sums.sum(axis=0)-b*row)/(b-1)
                values.append(list(amplitude_observables(replica,64,c.names,c.metadata).values()))
            v=np.asarray(values);v-=v.mean(axis=0)
            v=np.concatenate([v.real,v.imag],axis=1)
            cov+=(b-1)/b*(v.T@v)
        np.testing.assert_allclose(cov,e['covariance'],atol=1e-14)
        labels=e['names'];i=labels.index('VA');j=labels.index('chiA')
        self.assertAlmostEqual(e['covariance'][j,j],.09*e['covariance'][i,i],places=13)

    def test_reblock_preserves_chain_and_no_remainder_loss(self):
        c=chain(blocks=20);b=c.reblock(4)
        np.testing.assert_allclose(b.sums.sum(axis=0),c.sums.sum(axis=0))
        self.assertEqual(len(b.counts),10)
        with self.assertRaisesRegex(ValueError,'incomplete'):
            c.reblock(6)

    def test_zero_phase_and_binder_rejected(self):
        with self.assertRaisesRegex(ValueError,'phase'):
            amplitude_observables(np.array([0,1,2,3]),4,[],dict(Lx=2,Ly=2))
        with self.assertRaisesRegex(ValueError,'Binder'):
            amplitude_observables(np.array([1,0,0,0]),1,[],dict(Lx=2,Ly=2))

    def test_fixed_geometry_and_seeds(self):
        c,d=chain(1,.3),chain(2,.7)
        self.assertEqual(len(scan_groups([c,d],1,0)),2)
        d.metadata['shear']=1
        with self.assertRaisesRegex(ValueError,'shear'):
            scan_groups([c,d],1,0)
        with self.assertRaisesRegex(ValueError,'seeds'):
            scan_groups([chain(1,.3),chain(1,.7)],1,0)
        with self.assertRaisesRegex(ValueError,'nonzero'):
            mode(c.metadata,0,0)

    def test_joint_bootstrap_preserves_moment_relationship(self):
        groups=scan_groups([chain(1,.3),chain(2,.3),chain(3,.7)],1,0)
        names,v,failed=bootstrap_curves(groups,120,np.random.default_rng(73))
        self.assertEqual(failed,0)
        np.testing.assert_allclose(v[:,:,names.index('chiA')],.3*v[:,:,names.index('VA')],atol=2e-15)

    def test_nonuniform_quadratic_fit_exact_derivative(self):
        x=np.array([.2,.29,.4,.58,.72,.93,1.3])
        for w in [3,5,7]:
            centers,s=fit_slopes(x,3*x*x-2*x+4,w)
            np.testing.assert_allclose(s,6*x[centers]-2,atol=1e-13)

    def test_local_peak_and_window_dependence(self):
        x=np.linspace(0,2,41);b=2-1/(1+np.exp(-(x-1)*7))
        boot=np.tile(b,(200,1))
        rows,peaks=peak_scan(x,b,boot,[3,5,7],2,True)
        self.assertEqual(len(rows),39+37+35)
        for peak in peaks:
            self.assertAlmostEqual(peak['lambda_grid_peak'],1.,places=13)
            self.assertIn('candidate',peak['status'])

    def test_boundary_and_failed_bootstrap_not_certified(self):
        x=np.linspace(0,1,7);b=x*x
        _,peaks=peak_scan(x,b,np.tile(b,(100,1)),[3],2,True)
        self.assertIn('truncated',peaks[0]['status'])
        boot=np.tile(b,(100,1));boot[0,0]=np.nan
        rows,peaks=peak_scan(x,b,boot,[3],2,True)
        self.assertFalse(rows);self.assertIn('withheld',peaks[0]['status'])
        _,peaks=peak_scan(x,b,np.tile(b,(100,1)),[7],2,True)
        self.assertIn('undefined',peaks[0]['status'])
        self.assertNotIn('lambda_peak_95_interval',peaks[0])

    def test_full_analysis_and_doubled_blocks(self):
        chains=[chain(i+1,lam=.2+.1*i) for i in range(7)]
        out=run_analysis(chains,1,0,[3,5,7],2,120,531)
        self.assertEqual(len(out['points']),7)
        self.assertEqual(len(out['binder_slopes']),5+3+1)
        self.assertEqual(len(out['doubled_block_signed_slope_errors']),9)
        for row in out['doubled_block_sensitivity']:
            self.assertEqual(row['block_sweeps'],4)
            self.assertTrue(row['quality_flags']) # 10 blocks; deliberately short fixture.

    def test_free_real_modulus_response_dt_convergence(self):
        # Independent continuous stationary Gaussian covariance of |D|:
        # v*(2/pi)*(sqrt(1-rho^2)+rho*asin(rho)-1), rho=exp(-omega*t).
        nc=6;variance=1/(2*4*.5*nc);omega=.5;tau=4.
        def cov(t):
            rho=np.exp(-omega*np.asarray(t))
            return variance*2/np.pi*(np.sqrt(1-rho*rho)+rho*np.arcsin(rho)-1)
        nodes,weights=np.polynomial.legendre.leggauss(160)
        reference=nc*tau/2*np.dot(weights,cov((nodes+1)*tau/4))
        errors=[]
        for intervals in [20,40,80,160]:
            dt=tau/intervals;t=np.arange(intervals+1)*dt
            w=np.full(intervals+1,dt);w[[0,-1]]*=.5
            response=nc*np.dot(w,cov(abs(t-tau/2)))
            errors.append(abs(response-reference))
        self.assertTrue(all(errors[j+1]<.28*errors[j] for j in range(3)))
        self.assertLess(errors[-1],3e-5)


if __name__=='__main__':
    unittest.main()
