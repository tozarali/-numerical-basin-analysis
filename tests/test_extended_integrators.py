"""Tests for the tilted well, the generic discrete-gradient and the rho-Bathe integrators."""
import unittest

import numpy as np

from basinlab import bathe, dg_mdof
from basinlab.beam2mode import Beam2Mode
from basinlab.integrators_extra import discrete_gradient_step
from basinlab.mdof_genalpha import integrate_to_masked
from basinlab.model import classify as classify_sym
from basinlab.reference_mdof import rk4
from basinlab.well import Well


class Linear:
    """V = w^2 q^2 / 2 (d=1)."""
    def __init__(self, w=1.0):
        self.w = w

    def potential(self, q):
        return 0.5*self.w**2*q[..., 0]**2

    def force_tangent(self, q):
        return self.w**2*q, np.full(q.shape + (1,), self.w**2)


class WellTests(unittest.TestCase):
    def test_symmetric_limit(self):
        w = Well(0.0)
        self.assertAlmostEqual(w.barrier, 0.25)
        self.assertAlmostEqual(w.saddle, 0.0, places=12)
        x = np.array([[-1.0], [0.3], [1.0], [1.9]]); v = np.array([[0.1], [0.0], [-0.1], [0.5]])
        np.testing.assert_array_equal(w.classify(x, v), classify_sym(x[:, 0], v[:, 0]))

    def test_tilt_shifts_saddle_and_barrier(self):
        w = Well(0.1)
        e = w.equilibria
        self.assertEqual(len(e), 3)
        np.testing.assert_allclose(w.force(e[:, None]).ravel(), 0.0, atol=1e-12)
        self.assertGreater(w.potential(np.array([[e[0]]]))[0], w.potential(np.array([[e[2]]]))[0])
        self.assertLess(w.saddle, 0.0)

    def test_merged_wells_rejected(self):
        with self.assertRaises(ValueError):
            Well(0.5)


class DGTests(unittest.TestCase):
    def test_matches_scalar_dg_for_symmetric_well(self):
        rng = np.random.default_rng(1)
        x = rng.uniform(-2, 2, 30); v = rng.uniform(-2, 2, 30)
        x1, v1 = discrete_gradient_step(x, v, 0.4, 0.15)
        m = np.eye(1); c = np.array([[0.15]])
        q1, w1, failed = dg_mdof.step_masked(x[:, None], v[:, None], 0.4, m, c, Well(0.0))
        self.assertFalse(failed.any())
        np.testing.assert_allclose(q1[:, 0], x1, atol=1e-9)
        np.testing.assert_allclose(w1[:, 0], v1, atol=1e-8)

    def test_energy_conserved_undamped_in_2dof(self):
        mdl = Beam2Mode(P=2.0, alpha=1.0)
        rng = np.random.default_rng(2)
        q = rng.uniform(-1, 1, (12, 2)); v = rng.uniform(-1, 1, (12, 2))
        e0 = mdl.energy(q, v)
        qT, vT, nfail = dg_mdof.integrate_to_masked(q, v, 20.0, 0.5, np.eye(2), np.zeros((2, 2)), mdl)
        self.assertEqual(nfail, 0)
        np.testing.assert_allclose(mdl.energy(qT, vT), e0, rtol=0, atol=1e-7)

    def test_second_order(self):
        mdl = Well(0.1)
        q = np.array([[1.3], [-0.4]]); v = np.array([[0.2], [0.9]])
        m, c = np.eye(1), np.array([[0.15]])
        ref, _ = rk4(q, v, 4.0, 1e-3, np.ones(1), c, mdl)[0], None
        errs = []
        for h in (0.2, 0.1):
            qT, _, _ = dg_mdof.integrate_to_masked(q, v, 4.0, h, m, c, mdl)
            errs.append(np.max(np.abs(qT - ref)))
        self.assertGreater(errs[0]/errs[1], 3.3)


class BatheTests(unittest.TestCase):
    def test_coefficients(self):
        g, c = bathe.coefficients(1.0)
        self.assertAlmostEqual(g, 0.5)
        np.testing.assert_allclose(c, (0.25, 0.5, 0.25))
        g0, c0 = bathe.coefficients(0.0)
        self.assertAlmostEqual(g0, 2 - np.sqrt(2))
        for rho in (0.0, 0.3, 0.7, 1.0):
            self.assertAlmostEqual(sum(bathe.coefficients(rho)[1]), 1.0, places=12)
        # standard Bathe (gamma = 1/2, rho = 0) has all three coefficients 1/3
        _, c = bathe.coefficients(0.0, gamma=0.5)
        np.testing.assert_allclose(c, (1/3, 1/3, 1/3))

    def amplification(self, rho, wh):
        mdl, m, c = Linear(1.0), np.eye(1), np.zeros((1, 1))
        cols = []
        for x0, v0 in ((1.0, 0.0), (0.0, 1.0)):
            q = np.array([[x0]]); v = np.array([[v0]]); a = -mdl.force_tangent(q)[0]
            q1, v1, a1, failed = bathe.step_masked(q, v, a, wh, m, c, mdl, rho)
            self.assertFalse(failed.any())
            cols.append([q1[0, 0], v1[0, 0]])
        return np.array(cols).T

    def test_spectral_radius_at_large_step(self):
        for rho in (0.0, 0.3, 0.6, 1.0):
            A = self.amplification(rho, 1e4)
            self.assertAlmostEqual(max(abs(np.linalg.eigvals(A))), rho, delta=2e-3)

    def test_unconditionally_stable(self):
        for rho in (0.0, 0.5):
            for wh in (0.1, 1.0, 10.0, 100.0):
                self.assertLessEqual(max(abs(np.linalg.eigvals(self.amplification(rho, wh)))), 1 + 1e-9)

    def test_second_order_accuracy(self):
        mdl, m, c = Linear(1.0), np.eye(1), np.zeros((1, 1))
        for rho in (0.0, 0.5, 1.0):
            errs = []
            for h in (0.2, 0.1):
                q, v, nf = bathe.integrate_to_masked(np.array([[1.0]]), np.array([[0.0]]), 4.0, h, m, c, mdl, rho)
                errs.append(abs(q[0, 0] - np.cos(4.0)))
            self.assertGreater(errs[0]/errs[1], 3.5, (rho, errs))

    def test_nonlinear_agrees_with_reference(self):
        mdl = Well(0.1)
        q = np.array([[1.3], [-0.4], [0.2]]); v = np.array([[0.2], [0.9], [-0.5]])
        m, c = np.eye(1), np.array([[0.15]])
        ref = rk4(q, v, 6.0, 1e-3, np.ones(1), c, mdl)[0]
        qT, _, nf = bathe.integrate_to_masked(q, v, 6.0, 0.05, m, c, mdl, 0.5)
        self.assertEqual(nf, 0)
        np.testing.assert_allclose(qT, ref, atol=2e-3)


if __name__ == "__main__":
    unittest.main()


class RegistryTests(unittest.TestCase):
    def test_new_experiments_registered(self):
        from basinlab.run_all import BY_NAME
        for n in ("asym_well", "beam2mode_controls", "fe_dg_control", "fe_crosscheck_openseespy"):
            self.assertIn(n, BY_NAME)
            self.assertTrue(BY_NAME[n].quick_ok)
