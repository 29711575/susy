"""Deterministic recovery/statistics checks; no MC or scheduler submissions."""
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qmc_tools.blocks import Chain, checked_group, read_chain
from qmc_tools.diagnostics import rank_diagnostic
from qmc_tools.observables import transform_for
from qmc_tools.registry import import_submission, merge_snapshots, indices
from qmc_tools.statistics import joint

BASE = Path(__file__).resolve().parents[2]


class StatisticsTests(unittest.TestCase):
    def chain(self, values, seed=1, phase=None):
        phase = np.ones(len(values), complex) if phase is None else np.array(phase, complex)
        return Chain(Path(f'/tmp/seed{seed}_blocks.csv'), ['O'], np.ones(len(values), dtype=np.int64),
                     np.array([phase, phase*np.array(values)]).T,
                     np.zeros((len(values), 3), dtype=np.int64), dict(seed=seed, measure_every=1), {})

    def test_signed_complex_common_phase(self):
        phase = np.exp(1j*np.array([0., .2, .4, 2.7, -.1]))
        values = np.array([2., 4., 1., 7., 3.])
        c = self.chain(values, phase=phase)
        result = joint([c])['estimates']['O']
        ratio = (phase @ values)/phase.sum()
        delete = np.array([(phase @ values-p*v)/(phase.sum()-p) for p, v in zip(phase, values)])
        reference = np.sqrt(4/5*np.sum((delete.real-delete.real.mean())**2))
        self.assertAlmostEqual(result['real'], ratio.real)
        self.assertAlmostEqual(result['imag'], ratio.imag)
        self.assertAlmostEqual(result['error_real'], reference)

    def test_unequal_independent_chains_keep_original_allocation(self):
        a, b = self.chain([1, 4, 8], 1), self.chain([10, 12, 17, 19, 26], 2)
        result = joint([a, b])['estimates']['O']
        self.assertAlmostEqual(result['real'], 97/8)
        reference = sum((len(c.counts)/8)**2*np.var(c.sums[:, 1].real, ddof=1)/len(c.counts) for c in [a, b])
        self.assertAlmostEqual(result['error_real']**2, reference)

    def test_no_cross_chain_or_remainder_discard(self):
        a, b = self.chain([1, 2, 3], 1), self.chain([4], 2)
        grouped = a.reblock(2)
        self.assertEqual(grouped.counts.tolist(), [2, 1])
        self.assertEqual(grouped.sums[:, 1].real.tolist(), [3, 3])
        self.assertEqual(joint([a, b])['samples'], 4)
        self.assertFalse(joint([a, b])['jackknife_defined'])

    def test_weighted_tail_groups_match_independent_pseudovalue_reference(self):
        c = self.chain([1, 1, 1])
        c.counts = np.array([2, 1, 1])
        c.sums = np.array([[2, 6], [1, 8], [1, 16]], dtype=complex)
        result = joint([c])['estimates']['O']
        self.assertAlmostEqual(result['real'], 7.5)
        self.assertAlmostEqual(result['error_real']**2, 14.125)

    def test_duplicate_seeds_and_ancestry_rejected(self):
        a, b = self.chain([1, 2]), self.chain([3, 4])
        with self.assertRaises(ValueError):
            checked_group([a, b])
        b.metadata['seed'] = 2
        b.metadata['resume'] = dict(ancestor_output_prefix=str(a.path)[:-11])
        with self.assertRaises(ValueError):
            checked_group([a, b])

    def test_phase_zero_is_undefined(self):
        with self.assertRaisesRegex(ValueError, 'phase'):
            joint([self.chain([1, 4], phase=[1, -1])])

    def test_rank_diagnostics_detect_chain_disagreement(self):
        rng = np.random.default_rng(591)
        x = rng.normal(size=(4, 2000))
        self.assertLess(rank_diagnostic(x)['rhat'], 1.01)
        x[0] += 3
        bad = rank_diagnostic(x)
        self.assertGreater(bad['rhat'], 1.1)
        self.assertIsNone(bad['ess_bulk_blocks'])
        self.assertIsNone(rank_diagnostic(np.ones((4, 20)))['rhat'])


class RecoveryTests(unittest.TestCase):
    def test_existing_cpp_complete_result_parity(self):
        path = BASE/'audit/1023/signal-checkpoint/tests-final/new-local_results.json'
        if not path.exists():
            self.skipTest('historical reference fixture not present')
        old = json.loads(path.read_text())
        c = read_chain(path)
        result = joint([c], transform_for(c.names, c.metadata))['estimates']
        for name, ref in {**old['observables'], **old['vbs']['blocking'][0]['derived']}.items():
            got = result[name]
            self.assertAlmostEqual(got['real'], ref['value'][0], delta=1e-10, msg=name)
            self.assertAlmostEqual(got['imag'], ref['value'][1], delta=1e-10, msg=name)
            for key, target in [('error_real', 'stderr_real'), ('error_imag', 'stderr_imag'), ('covariance_real_imag', 'cov_real_imag')]:
                if ref[target] is not None:
                    self.assertAlmostEqual(got[key], ref[target], delta=1e-10, msg=name+key)

    def test_interacting_cpp_result_parity(self):
        path = BASE/'audit/1023/signal-checkpoint/tests-final/interacting-full_results.json'
        if not path.exists():
            self.skipTest('historical interacting reference fixture absent')
        old = json.loads(path.read_text()); c = read_chain(path)
        transform = transform_for(c.names, c.metadata)
        for level, reference in enumerate(old.get('vbs', {}).get('blocking', [])):
            result = joint([c.reblock(2**level)], transform)['estimates']
            for name, ref in {**reference['raw'], **reference['derived']}.items():
                self.assertAlmostEqual(result[name]['real'], ref['value'][0], delta=1e-10, msg=name)
                if ref['stderr_real'] is not None:
                    self.assertAlmostEqual(result[name]['error_real'], ref['stderr_real'], delta=1e-10, msg=name)

    def test_aborted_final_summary_recovered_without_mc(self):
        path = BASE/'audit/1023/signal-checkpoint/edge-tests/baseline-vbs-combination_results.json'
        if not path.exists():
            self.skipTest('historical failed-finalization fixture not present')
        with self.assertRaises(ValueError):
            read_chain(path)
        c = read_chain(path, metadata=dict(shear=0), allow_prefix=True)
        self.assertEqual(c.metadata['status'], 'aborted')
        result = joint([c], transform_for(c.names, c.metadata))
        self.assertEqual(result['samples'], 20)
        self.assertIn('C_connected_q0_zz', result['estimates'])
        self.assertTrue(result['jackknife_defined'])

    def test_corrupt_full_record_and_partial_tail_are_distinct(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/'mc_blocks.csv'
            header = 'block,count,sum_phase_real,sum_phase_imag,positive_count,negative_count,complex_count,sum_phase_O_real,sum_phase_O_imag\n'
            p.write_text(header+'0,1,1,0,1,0,0,2,0\n1,1,1')
            c = read_chain(p, dict(measure_every=1), allow_prefix=True)
            self.assertEqual(len(c.counts), 1)
            self.assertGreater(c.provenance['incomplete_tail_bytes'], 0)
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                read_chain(p)
            p.write_text(header+'0,1,1,0,1,0,0,2,0\n1,1,1\n')
            with self.assertRaisesRegex(ValueError, 'malformed'):
                read_chain(p, dict(measure_every=1), allow_prefix=True)

    def test_block_window_requires_exact_boundaries(self):
        c = StatisticsTests().chain([1, 2, 3, 4])
        c.counts[:] = 10
        self.assertEqual(c.window(10, 30).counts.sum(), 20)
        with self.assertRaises(ValueError):
            c.window(5)

    def test_nonlinear_response_uses_ratios_before_connection(self):
        meta = dict(Lx=2, Ly=2)
        names = [f'{s}_q0_{a}' for s in ['ChiD', 'ChiIconj'] for a in 'xyz']+[f'ChiRaw_q0_{a}{b}' for a in 'xyz' for b in 'xyz']
        derive = transform_for(names, meta)
        r = dict.fromkeys(names, 1+0j)
        r.update(ChiD_q0_z=2+.5j, ChiIconj_q0_z=3-.1j, ChiRaw_q0_zz=30+1j)
        self.assertAlmostEqual(derive(r)['chi_q0_zz'], (30+1j)-4*(2+.5j)*(3-.1j))

    def test_energy_components_and_center_location_are_distinct(self):
        names = ['energy_bare_matter', 'energy_matter_phonon', 'energy_phonon_kinetic',
                 'energy_phonon_potential', 'energy_external_source', 'energy_total', 'X2_x', 'X2_y', 'X2_z']
        r = dict(zip(names, [-5., -2., 1., 3., -.2, -3.2, 2., 3., 4.]))
        derived = transform_for(names, dict(Lx=2, Ly=2, k=1))(r)
        self.assertAlmostEqual(derived['energy_sum_residual'], 0.)
        self.assertAlmostEqual(derived['energy_total_per_cell'], -.8)
        self.assertEqual(derived['energy_phonon_potential_center_reconstructed'], 18.)
        self.assertEqual(r['energy_phonon_potential'], 3.)
        with self.assertRaises(ValueError):
            transform_for(names[:1], dict(Lx=2, Ly=2))

    def test_center_matter_energy_uses_recorded_gauge(self):
        names = [f'{s}_b{b}' for s in ['BondK', 'BondXK'] for b in range(3)]
        r = dict(BondK_b0=1., BondK_b1=2., BondK_b2=3., BondXK_b0=4., BondXK_b1=5., BondXK_b2=6.)
        meta = dict(Lx=1, Ly=1, Jx=1., Jy=2., Jz=.5, g=.4, gauge_u=[-1, 1, -1])
        derived = transform_for(names, meta)(r)
        self.assertEqual(derived['energy_bare_matter_center_reconstructed'], -1+4-1.5)
        self.assertAlmostEqual(derived['energy_matter_phonon_center_reconstructed'], .4*(-4+5-6))
        self.assertNotIn('energy_total', derived)
        with self.assertRaises(ValueError):
            transform_for(names, {**meta, 'gauge_u':[0, 1, 1]})


class RegistryTests(unittest.TestCase):
    def test_submission_indices_and_migration_subset(self):
        p = BASE/'audit/1023/queued-migrate-20022-tau1p2-supplement/submission.json'
        if not p.exists():
            self.skipTest('submission fixture absent')
        registry = import_submission(p)
        self.assertEqual(len(registry['records']), 59)
        self.assertEqual(len({r['seed'] for r in registry['records']}), 59)
        self.assertTrue(all(not r['currently_listed'] for r in registry['records']))
        self.assertEqual(indices('1-3,7,10-14:2'), [1, 2, 3, 7, 10, 12, 14])
        with self.assertRaises(ValueError):
            indices('1-3,3')

    def test_missing_live_job_not_assumed_finished(self):
        snap = dict(port=20022, snapshot_finished_utc='old', jobs=[dict(job_id='1[0].mgt', job_state='R')])
        previous = merge_snapshots([snap])
        new = merge_snapshots([dict(port=20022, snapshot_finished_utc='new', jobs=[])], previous)
        record = new['records'][0]
        self.assertFalse(record['currently_listed'])
        self.assertEqual(record['observed_utc'], 'old')
        self.assertEqual(record['scheduler_state'], 'R')

    def test_remote_qstat_parser_and_wide_tail(self):
        spec = importlib.util.spec_from_file_location('remote_inventory', BASE/'scripts/qmc_tools/remote_inventory.py')
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        rows = mod.fields('Job Id: 1[0].mgt\n    job_state = R\n    Variable_List = PBS_O_WORKDIR=/abc,\n\tFOO=bar\n')
        self.assertEqual(rows[0]['Variable_List'], 'PBS_O_WORKDIR=/abc,FOO=bar')
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/'wide.csv'
            p.write_bytes(b'header\n0,10,'+b'x'*90000+b'\n1,10,incomplete')
            self.assertEqual(mod.last_record(str(p))[:2], ['0', '10'])

    def test_gamma_offset_and_nonstandard_record_layout(self):
        spec = importlib.util.spec_from_file_location('remote_inventory', BASE/'scripts/qmc_tools/remote_inventory.py')
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'points.json').write_text(json.dumps([dict(index=24, array_index=0, seed=999)]))
            (root/'tau10.pbs').write_text('#PBS -N gamma48T10\nexport CHAIN_OFFSET=24\n')
            folder = root/'runs/index-24'; folder.mkdir(parents=True)
            (folder/'point.json').write_text(json.dumps(dict(index=24, seed=999)))
            _, points, _, folders, offsets = mod.root_data(str(root))
            self.assertEqual(offsets['gamma48T10'][0][0], 24)
            self.assertEqual(points[24]['seed'], 999)
            self.assertEqual(folders[24], [str(folder)])


if __name__ == '__main__':
    unittest.main()
