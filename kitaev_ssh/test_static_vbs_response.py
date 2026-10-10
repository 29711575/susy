#!/usr/bin/env python3
import unittest
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np

from static_vbs_response import (allowed_grid, bloch_energy, dense_reference,
                                 finite_matter_energy, mode_basis, momentum_reference,
                                 static_hessian, finite_difference_hessians,
                                 hopping_matrix)


class StaticVBSResponseTests(unittest.TestCase):
    JEFF = [1.0541670812557202, 1.0541670812557202, 0.5285238690991194]
    G = .3
    XBAR = [0.18055693751906726, 0.1805569375190672, 0.09507956366373145]

    def test_finite_energy_matches_bloch_sum(self):
        J = [1., 1., .5]
        for geom in ((6, 6, 0), (12, 3, 0), (12, 3, 3)):
            real = finite_matter_energy(*geom, J)
            bloch = sum(bloch_energy(*p["fraction"], J) for p in allowed_grid(*geom))
            self.assertAlmostEqual(real, bloch, places=11)

    def test_g_zero_hessian_is_identity(self):
        h = static_hessian(6, 6, 0, [1., 1., .5], 0., (.5, 1./3.), 2e-4)
        self.assertTrue(np.allclose(h["matter_hessian"], 0., atol=1e-11))
        self.assertTrue(np.allclose(h["total_hessian"], np.eye(h["mode_count"]), atol=1e-11))

    def test_q_minus_q_spectra_and_self_conjugate_sine(self):
        a = static_hessian(6, 6, 0, [1., 1., .5], .3, (.5, 1./3.), 2e-4)
        b = static_hessian(6, 6, 0, [1., 1., .5], .3, (0.5, 2./3.), 2e-4)
        self.assertTrue(np.allclose(a["eigenvalues"], b["eigenvalues"], atol=1e-10))
        self.assertEqual(len(mode_basis(6, 6, 0, (.5, .5))), 3)

    def test_dense_reference_is_finite_and_has_unit_phonon_shift(self):
        ref = dense_reference((.5, 1./3.), [1.054, 1.054, .529], .3, grid=40)
        self.assertTrue(np.isfinite(ref["min_eigenvalue"]))
        self.assertEqual(len(ref["total_hessian"]), 6)

    def test_small_amplitude_second_difference_is_stable(self):
        q = (.5, 1./3.)
        modes = mode_basis(6, 6, 0, q)
        mode = modes[0][2]
        base = finite_matter_energy(6, 6, 0, [1.054, 1.054, .529])
        values = []
        for h in (1e-3, 5e-4, 2.5e-4):
            ep = finite_matter_energy(6, 6, 0, [1.054, 1.054, .529], mode*h, .3)
            em = finite_matter_energy(6, 6, 0, [1.054, 1.054, .529], -mode*h, .3)
            values.append((ep - 2*base + em) / (36*h**2) + 1.)
        self.assertLess(max(values) - min(values), 2e-6)

    def test_allowed_k_reference_matches_real_space_elementwise(self):
        cases = [((6, 6, 0), (.5, 1./3.)),
                 ((12, 3, 0), (.25, 1./3.)),
                 ((12, 3, 3), (5./12., 7./12.)),
                 ((6, 6, 0), (.5, .5))]
        for geom, q in cases:
            finite = np.asarray(static_hessian(*geom, self.JEFF, self.G, q)["matter_hessian"])
            klist = [p["fraction"] for p in allowed_grid(*geom)]
            reference = np.asarray(momentum_reference(q, self.JEFF, self.G, klist)["matter_hessian"])
            self.assertTrue(np.allclose(finite, reference, atol=2e-11), (geom, q, np.max(abs(finite-reference))))

    def test_mixed_bond_mode_energy_second_difference(self):
        geom = (6, 6, 0); q = (.5, 1./3.)
        modes = mode_basis(*geom, q)
        coeff = np.array([.31, -.22, .17, .29, -.13, .41])
        coeff /= np.linalg.norm(coeff)
        delta = sum((coeff[i] * modes[i][2] for i in range(len(modes))), np.zeros((36, 3)))
        base = finite_matter_energy(*geom, self.JEFF)
        x2 = np.dot(self.XBAR, self.XBAR) * 36
        cross = float(np.sum(np.asarray(self.XBAR)[None, :] * delta))
        values = []
        for h in (1e-3, 5e-4):
            norm2 = float(np.sum(delta*delta))
            ep = finite_matter_energy(*geom, self.JEFF, h*delta, self.G) + .5*(x2 + 2*h*cross + h*h*norm2)
            em = finite_matter_energy(*geom, self.JEFF, -h*delta, self.G) + .5*(x2 - 2*h*cross + h*h*norm2)
            values.append((ep - 2*(base + .5*x2) + em) / (36*h**2))
        expected = coeff @ np.asarray(static_hessian(*geom, self.JEFF, self.G, q)["total_hessian"]) @ coeff
        self.assertAlmostEqual(values[0], expected, places=6)
        self.assertAlmostEqual(values[1], expected, places=6)

    def test_phonon_k_and_majorana_energy_normalization(self):
        for geom in ((6, 6, 0), (12, 3, 3)):
            q = allowed_grid(*geom)[5]["fraction"]
            for k in (.37, 1.088, 2.3):
                h = static_hessian(*geom, self.JEFF, 0., q, spring_k=k)
                self.assertTrue(np.allclose(h["total_hessian"], k*np.eye(h["mode_count"])))
                ref = momentum_reference(q, self.JEFF, 0.,
                                         [p["fraction"] for p in allowed_grid(*geom)], spring_k=k)
                self.assertTrue(np.allclose(ref["total_hessian"], h["total_hessian"]))
            t = hopping_matrix(*geom, self.JEFF)
            a = np.block([[np.zeros_like(t), 2*t], [-2*t.T, np.zeros_like(t)]])
            self.assertAlmostEqual(-.25*np.sum(abs(np.linalg.eigvalsh(1j*a))),
                                   finite_matter_energy(*geom, self.JEFF), places=11)

    def test_full_energy_difference_matrix_three_steps(self):
        for geom, q in (((6, 6, 0), (.5, 1./3.)),
                        ((12, 3, 3), (5./12., 7./12.)),
                        ((6, 6, 0), (.5, .5))):
            analytic = np.array(static_hessian(*geom, self.JEFF, self.G, q,
                                               spring_k=1.088)["total_hessian"])
            checks = finite_difference_hessians(*geom, self.JEFF, self.XBAR,
                                                self.G, q, (1e-3, 5e-4, 2.5e-4),
                                                spring_k=1.088, analytic=analytic)
            self.assertEqual(len(checks["checks"]), 3)
            for check in checks["checks"]:
                self.assertLess(check["max_abs_matrix_error"], 5e-7)
                self.assertLess(abs(check["mixed_curvature_error"]), 5e-7)

    def test_periodicity_and_mode_gram_on_all_three_grids(self):
        for geom in ((6, 6, 0), (12, 3, 0), (12, 3, 3)):
            for point in allowed_grid(*geom):
                modes = mode_basis(*geom, point["fraction"])
                v = np.array([m[2].reshape(-1) for m in modes])
                self.assertTrue(np.allclose(v@v.T/(geom[0]*geom[1]), np.eye(len(modes)), atol=1e-12))
        with self.assertRaisesRegex(ValueError, "allowed sheared grid"):
            mode_basis(12, 3, 3, (.4, .6))

    def test_zero_modes_flagged_with_energy_scan_and_no_hessian(self):
        result = static_hessian(3, 3, 0, [1., 1., 1.], .3, (1./3., 2./3.))
        self.assertEqual(result["status"], "zero_mode_uncertified")
        self.assertIsNone(result["min_eigenvalue"])
        self.assertIsNone(result["matter_hessian"])
        self.assertGreater(result["zero_mode_count"], 0)
        self.assertEqual(len(result["energy_scans"]), 6)
        ref = momentum_reference((1./3., 2./3.), [1., 1., 1.], .3,
                                 [p["fraction"] for p in allowed_grid(3, 3, 0)])
        self.assertEqual(ref["status"], "zero_mode_uncertified")
        free = static_hessian(3, 3, 0, [1., 1., 1.], 0., (1./3., 2./3.), spring_k=1.088)
        self.assertTrue(np.allclose(free["total_hessian"], 1.088*np.eye(6)))

    def test_continuous_q_periodicity_and_conjugate_spectrum(self):
        q = (.58043, .41957)
        for other in ((q[0]+1, q[1]-1), (1-q[0], 1-q[1])):
            a = dense_reference(q, self.JEFF, self.G, grid=32, spring_k=1.088)
            b = dense_reference(other, self.JEFF, self.G, grid=32, spring_k=1.088)
            self.assertTrue(np.allclose(a["eigenvalues"], b["eigenvalues"], atol=1e-12))

    def test_cli_records_real_reference_changes_and_arbitrary_uniform_grid(self):
        script = Path(__file__).with_name("static_vbs_response.py")
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
            subprocess.run([sys.executable, str(script), "--output-dir", tmp,
                            "--uniform-grid", "64", "--reference-grids", "16,32",
                            "--geometries", "3,3,0"], check=True, capture_output=True, env=env)
            data = json.loads((Path(tmp)/"static_response.json").read_text())
            self.assertAlmostEqual(data["parameters"]["spring_k"], 1.088)
            self.assertEqual(data["common_background_grid"], 64)
            self.assertEqual(set(data["uniform_convergence"]), {"32", "64"})
            self.assertTrue(data["common_q_geometry_comparisons"])
            geo = data["geometries"][0]
            self.assertEqual(len(geo["all_q"]), 8)
            self.assertTrue(geo["soft_search_complete"])
            for row in geo["all_q"]:
                self.assertLess(row["finite_k_crosscheck_max_abs"], 1e-11)
                if row["reference"]:
                    self.assertEqual(set(row["references"]), {"16", "32"})
                    change = row["reference_grid_changes"][0]
                    self.assertAlmostEqual(change["min_eigenvalue_change"],
                                           row["references"]["32"]["min_eigenvalue"]-
                                           row["references"]["16"]["min_eigenvalue"])
                    self.assertEqual(len(row["finite_difference"]["checks"]), 3)
            qfile = Path(tmp)/"q_candidates_Lx3_Ly3_s0.txt"
            points = [line.split() for line in qfile.read_text().splitlines() if not line.startswith("#")]
            coords = {(float(p[1]), float(p[2])) for p in points}
            self.assertNotIn((0., 0.), coords)
            for u, v in coords:
                self.assertTrue(any(np.allclose([(1-u)%1, (1-v)%1], q) for q in coords))
            rejected = subprocess.run([sys.executable, str(script), "--output-dir", tmp,
                                       "--spring-k", "1"], capture_output=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn(b"must equal mass*omega-ph^2", rejected.stderr)


if __name__ == "__main__":
    unittest.main()
