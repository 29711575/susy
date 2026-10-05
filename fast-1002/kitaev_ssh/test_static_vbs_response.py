#!/usr/bin/env python3
import unittest
import numpy as np

from static_vbs_response import (allowed_grid, bloch_energy, dense_reference,
                                 finite_matter_energy, mode_basis, momentum_reference,
                                 static_hessian)


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


if __name__ == "__main__":
    unittest.main()
