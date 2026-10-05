#!/usr/bin/env python3
"""Only the new offline chain handling; reuses real P4 outputs for regression."""
import unittest
from pathlib import Path
import numpy as np
from convergence_analysis import Chain, checked_group, comparison, fixed_observables, joint_estimate, read_chain, PRIMARY, EPSILON, second_order_residual, add_block_sensitivity


def metadata(seed=1):
    return dict(Lx=2, Ly=2, Jx=1., Jy=1., Jz=.5, g=.3, mass=1.7, omega_ph=.8,
                scope='fixed_gauge_matter_sector', spatial_boundary='pbc', wilson_loops=[1, 1],
                time_contour='open', intervals=4, dt=.2, seed=seed, measure_every=1)


def make_chain(block_values, seed=1, count=2):
    means=np.array(block_values, dtype=complex)
    sums=np.column_stack([np.ones(len(means))*count, means*count])
    return Chain(Path(f'chain_{seed}'), metadata(seed), ['x'], sums, np.full(len(means), count))


def transform(sums, count, names, meta):
    z=sums[1]/sums[0]
    return {'phase': sums[0]/count, 'x': z, 'nonlinear': z*z}


class ConvergenceTests(unittest.TestCase):
    def test_l4_momentum_mapping_normalization_and_phase(self):
        from convergence_analysis import fixed_momentum_mapping, full_grid_observables
        meta=metadata(); meta.update(Lx=4,Ly=4)
        mapping=fixed_momentum_mapping(meta)
        self.assertEqual(mapping['q1']['index'],2)
        self.assertEqual(mapping['q3']['index'],10)
        values={f'S_q{q}_{mu}{nu}':complex(q+1 if mu==nu else 0) for q in range(16) for mu in 'xyz' for nu in 'xyz'}
        values.update(Delta_q0_x=.2,Delta_q0_y=.3,Delta_q0_z=.4,X2_z=.7)
        phase=7*np.exp(.43j)
        sums=np.array([phase]+[phase*v for v in values.values()])
        r=fixed_observables(sums,7,list(values),meta)
        self.assertAlmostEqual(abs(r['Szz_q1']-3),0,places=13)
        self.assertAlmostEqual(abs(r['TrS_q3']-33),0,places=13)
        self.assertAlmostEqual(abs(r['uniform_background_trace']-16*(.04+.09+.16)),0,places=13)
        self.assertAlmostEqual(abs(r['Szz_q3_minus_other_nonzero_mean']-(11-np.mean([q+1 for q in range(1,16) if q!=10]))),0,places=13)
        grid=full_grid_observables(sums,7,list(values),meta)
        self.assertEqual(len(grid),145)
        self.assertAlmostEqual(abs(grid['S_q15_zz']-16),0,places=13)
        meta['vbs']={'Nc':4}
        with self.assertRaisesRegex(ValueError,'Nc'):
            fixed_momentum_mapping(meta)

    def test_tau_metadata_and_legacy_compatibility(self):
        from plot_vbs import normalize_tau_metadata
        new = {'tau': 3.2, 'intervals': 16, 'dt': .2}
        self.assertEqual(normalize_tau_metadata(new)['tau'], 3.2)
        old = {'theta_per_branch': 1.6, 'intervals': 16, 'dt': .2}
        self.assertEqual(normalize_tau_metadata(old)['tau'], 3.2)
        with self.assertRaisesRegex(ValueError, 'inconsistent imaginary-time metadata'):
            normalize_tau_metadata({'tau': 3.1, 'intervals': 16, 'dt': .2})

    def test_four_quarters_preserve_raw_ratio_and_chain_boundaries(self):
        from convergence_analysis import group_analysis, fixed_observables, STATIONARITY
        root=Path(__file__).resolve().parents[1]/'build-kitaev-p01/shifted_long/results'
        a=read_chain(root/'long_t32_a_results.json')
        b=read_chain(root/'long_t32_b_results.json')
        summary,_=group_analysis('test',[a,b],200)
        self.assertEqual(next(z for z in summary['between_chain_checks'] if z['observable']=='X_z')['family_size'],6)
        self.assertEqual(len(summary['quarter_chain_diagnostics']),8)
        self.assertEqual(len(summary['quarter_first_last_diagnostics']),12)
        quarter=next(z for z in summary['quarter_chain_diagnostics'] if z['seed']==a.seed and z['quarter']==2)
        rb=a.reblock(200)
        sums=rb.sums[10:20].sum(axis=0)
        expected=fixed_observables(sums,2000,a.names,a.metadata)
        self.assertEqual(quarter['blocks'],10)
        self.assertEqual(quarter['samples'],2000)
        self.assertEqual((quarter['measurement_sweep_start'],quarter['measurement_sweep_end']),(2000,4000))
        self.assertEqual(set(quarter['estimates']),set(STATIONARITY))
        for name in STATIONARITY:
            self.assertLess(abs(complex(*quarter['estimates'][name]['value'])-expected[name]),2e-14)
            self.assertGreater(quarter['estimates'][name]['stderr_real'],0)
        first=next(z for z in summary['quarter_chain_diagnostics'] if z['seed']==a.seed and z['quarter']==1)
        last=next(z for z in summary['quarter_chain_diagnostics'] if z['seed']==a.seed and z['quarter']==4)
        for name in STATIONARITY:
            contrast=next(z for z in summary['quarter_first_last_diagnostics'] if z['seed']==a.seed and z['observable']==name)
            self.assertAlmostEqual(contrast['fourth_minus_first'],last['estimates'][name]['value'][0]-first['estimates'][name]['value'][0],places=14)

    def test_split_uniform_background(self):
        from convergence_analysis import group_analysis
        c=read_chain(Path(__file__).resolve().parents[1]/'build-kitaev-p01/p4_runs/refresh_a_results.json')
        summary,_=group_analysis('test',[c],100)
        row=next(r for r in summary['split_chain_diagnostics'] if r['observable']=='uniform_background_trace')
        rb=c.reblock(100); mid=len(rb.counts)//2
        def background(sums):
            ratios=sums[1:]/sums[0]
            return 4*sum(ratios[c.names.index('Delta_q0_'+mu)]**2 for mu in 'xyz')
        expected=background(rb.sums[mid:].sum(axis=0))-background(rb.sums[:mid].sum(axis=0))
        self.assertAlmostEqual(row['second_half_minus_first'],expected.real,places=14)
        self.assertGreater(row['difference_se'],0)

    def test_shifted_trial_guards(self):
        from convergence_analysis import check_trial_comparison
        a,b=make_chain([1,2,3,4],1),make_chain([1,2,3,4],2)
        checked_group([a,b])
        self.assertEqual(a.metadata['trial_shift'],[0,0,0])
        b.metadata['trial_shift']=[.1,.2,.3]
        with self.assertRaisesRegex(ValueError,'different distributions'):
            checked_group([a,b])
        check_trial_comparison(a.metadata,b.metadata)
        with self.assertRaisesRegex(ValueError,'trial shift'):
            check_trial_comparison(a.metadata,b.metadata,step=True)
        b.metadata['intervals']=8
        with self.assertRaisesRegex(ValueError,'trial shift'):
            check_trial_comparison(a.metadata,b.metadata)
        b.metadata['trial_shift']=[0,0,0];b.metadata['dt']=.15
        with self.assertRaisesRegex(ValueError,'length and step'):
            check_trial_comparison(a.metadata,b.metadata)
    def test_prespecified_targets(self):
        import json
        plan=json.loads(Path(__file__).with_name('p5_plan.json').read_text())
        self.assertEqual(plan['primary_observables'],PRIMARY)
        self.assertEqual(plan['absolute_tolerances'],EPSILON)

    def test_one_chain_complex_joint_reference(self):
        chain=make_chain([1+.3j, -.2+.8j, 2-.4j, 1-.9j])
        phases=np.exp(1j*np.array([.2, -.6, .4, 1.1]))
        chain.sums *= phases[:, None]
        result=joint_estimate([chain], transform)
        total=chain.sums.sum(axis=0)
        replicas=[]
        for block in chain.sums:
            z=(total[1]-block[1])/(total[0]-block[0])
            replicas.append([(total[0]-block[0])/6, z, z*z])
        r=np.array(replicas);r-=r.mean(axis=0)
        ri=np.concatenate([r.real,r.imag],axis=1)
        np.testing.assert_allclose(result['covariance'], .75*ri.T@ri, atol=3e-14, rtol=0)
        np.testing.assert_allclose(result['values'], list(transform(total,8,None,None).values()), atol=1e-14, rtol=0)

    def test_unequal_chain_lengths_and_weights(self):
        c1=make_chain([1,2,4,8], seed=1, count=2)
        c2=make_chain([0,1,2,3,4,5], seed=2, count=4)
        result=joint_estimate([c1,c2],transform)
        weight=np.array([8,24])/32
        exact_variance=weight[0]**2*np.var([1,2,4,8],ddof=1)/4+weight[1]**2*np.var(range(6),ddof=1)/6
        self.assertAlmostEqual(result['values'][1].real, weight[0]*3.75+weight[1]*2.5, places=14)
        self.assertAlmostEqual(result['covariance'][1,1], exact_variance, places=14)
        self.assertEqual(result['chain_samples'],[8,24])

    def test_reblocking_does_not_cross_chains(self):
        a,b=make_chain([1,2,3],1,1),make_chain([4,5,6],2,1)
        for c in (a,b):
            with self.assertRaisesRegex(ValueError,'remainder'):
                c.reblock(2)
        c=make_chain([1,2,3,4],1,1).reblock(2)
        np.testing.assert_array_equal(c.sums,[[2,3],[2,7]])
        np.testing.assert_array_equal(c.counts,[2,2])

    def test_target_and_seed_guards(self):
        a,b=make_chain([1,2,3,4],1),make_chain([1,2,3,4],2)
        for key,new in [('dt',.1),('g',.4),('intervals',8),('wilson_loops',[-1,1])]:
            old=b.metadata[key];b.metadata[key]=new
            with self.assertRaisesRegex(ValueError,'different distributions'):
                checked_group([a,b])
            b.metadata[key]=old
        with self.assertRaisesRegex(ValueError,'distinct seeds'):
            checked_group([a,a])
        e=joint_estimate([a],transform)
        with self.assertRaisesRegex(ValueError,'independent seed'):
            comparison(e,e,'same_chain')

    def test_ratio_not_chain_ratio_average(self):
        a,b=make_chain([1,2,3,4],1),make_chain([2,4,6,8],2)
        a.sums *= np.exp(.7j)
        b.sums *= np.exp(-.4j)
        result=joint_estimate([a,b],transform)
        expected=(a.sums[:,1].sum()+b.sums[:,1].sum())/(a.sums[:,0].sum()+b.sums[:,0].sum())
        self.assertAlmostEqual(abs(result['values'][1]-expected),0,places=14)
        self.assertGreater(abs(expected-3.75),.1)

    def test_unresolved_phase_is_reported(self):
        c=make_chain(range(8),1,2)
        phases=.02+1j*np.r_[np.ones(4),-np.ones(4)]*np.sqrt(1-.02**2)
        c.sums*=phases[:,None]
        self.assertFalse(joint_estimate([c],transform)['phase_resolved_5sigma'])

    def test_conditional_step_sign_and_interval(self):
        row={'observable':'x','difference_b_minus_a':.03,'difference_se':.006,
             'ci95_low':.018,'ci95_high':.042,'look_safeguard95_low':.01,'look_safeguard95_high':.05}
        r=second_order_residual([row])[0]
        self.assertAlmostEqual(r['fine_minus_continuum_conditional'],-.01)
        self.assertAlmostEqual(r['stderr_conditional'],.002)
        self.assertAlmostEqual(r['ci95_low_conditional'],-.014)
        self.assertAlmostEqual(r['ci95_high_conditional'],-.006)

    def test_look_and_block_safeguard_never_narrows_interval(self):
        a=joint_estimate([make_chain([1,2,3,4],1)],transform)
        b=joint_estimate([make_chain([2,4,6,8],2)],transform)
        rows=comparison(a,b,'b-a',primary=['x'],epsilon={'x':2},looks=2)
        selected=next(r for r in rows if r['observable']=='x')
        self.assertLess(selected['look_safeguard95_low'],selected['simultaneous95_low'])
        group={'estimates':{r['observable']:{'stderr_real':.5} for r in rows},
               'blocking':[{'chain_blocks':[32], 'estimates':{r['observable']:{'stderr_real':2.} for r in rows}}]}
        add_block_sensitivity(rows,group,group)
        self.assertGreaterEqual(selected['max_block_sensitivity_se'],selected['difference_se'])

    def test_real_p4_regression(self):
        root=Path(__file__).resolve().parents[1]/'build-kitaev-p01/p4_runs'
        for name in ('seed_a','refresh_a','refresh_long'):
            c=read_chain(root/f'{name}_results.json')
            for level in c.metadata['vbs']['blocking']:
                if level['blocks']<2:
                    continue
                e=joint_estimate([c.reblock(level['block_sweeps'])])
                for label,group,p4name in [('Szz_q1','raw','S_q1_zz'),('TrS_q1','derived','TrS_q1'),
                                          ('X_z','raw','X_z')]:
                    i=e['names'].index(label);n=len(e['names']);old=level[group][p4name]
                    self.assertLess(abs(e['values'][i]-complex(*old['value'])),2e-13)
                    self.assertLess(abs(np.sqrt(e['covariance'][i,i])-old['stderr_real']),2e-13)
                    self.assertLess(abs(np.sqrt(e['covariance'][n+i,n+i])-old['stderr_imag']),2e-13)


if __name__=='__main__':
    unittest.main(verbosity=2)
