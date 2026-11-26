import unittest

import numpy as np
from scipy.integrate import solve_ivp

from basinlab import (beam2mode, diagnostic, dtheory, fe_arch, integrators,
                      integrators_extra, linear_analysis, mdof_genalpha, stats)
from basinlab.model import energy, force


class LinearAnalysis(unittest.TestCase):
    def test_rho_infinity_matches_family_map(self):
        for fam, target in [("generalized-alpha", .4), ("HHT", .7), ("WBZ", .6), ("Newmark", .9)]:
            p, rho = integrators.parameters_for_spectral_radius(fam, target)
            self.assertAlmostEqual(linear_analysis.spectral_radius(p, 1e8), rho, places=3)

    def test_trapezoid_closed_form(self):
        p = integrators.parameters("generalized-alpha", 1.0)
        for w in (0.1, 0.5, 1.0):
            r = linear_analysis.dissipation_and_period(p, w)
            self.assertAlmostEqual(r["xi_bar"], 0.0, places=10)
            self.assertAlmostEqual(r["period_error"], w/(2*np.arctan(w/2)) - 1, places=10)

    def test_physical_damping_recovered_small_step(self):
        p = integrators.parameters("generalized-alpha", .5)
        r = linear_analysis.dissipation_and_period(p, 1e-3, xi=0.05)
        self.assertAlmostEqual(r["xi_bar"], 0.05, places=5)


class DTheory(unittest.TestCase):
    def test_small_energy_period(self):
        self.assertAlmostEqual(dtheory.conservative_period(1e-6), 2*np.pi/np.sqrt(2), places=4)

    def test_dsep_tends_to_one(self):
        p = integrators.parameters("generalized-alpha", .5)
        self.assertAlmostEqual(dtheory.d_sep(p, 0.01, 0.15)["D_sep"], 1.0, places=3)

    def test_dtheory_tends_to_one(self):
        p = integrators.parameters("generalized-alpha", .5)
        self.assertAlmostEqual(dtheory.d_theory(p, 1e-3, 0.15, 0.5), 1.0, places=3)

    def test_crossing(self):
        self.assertAlmostEqual(dtheory.crossing([0, 1], [2, 0]), 0.5)


class Stats(unittest.TestCase):
    def test_parabola_and_bootstrap(self):
        r = np.linspace(0, 1, 21)
        rng = np.random.default_rng(1)
        prob = 0.1 + 0.4*(r - 0.43)**2
        errors = rng.random((21, 20000)) < prob[:, None]
        out = stats.bootstrap_optimum(r, errors, n_boot=200)
        self.assertLess(abs(out["rho_star"] - 0.43), 0.05)
        self.assertLess(out["ci95"][0], out["rho_star"])


class ExtraIntegrators(unittest.TestCase):
    def test_discrete_gradient_conserves_energy(self):
        x, v = np.array([0.3, 1.5]), np.array([1.0, 0.2])
        e0 = energy(x, v)
        for _ in range(500):
            x, v = integrators_extra.discrete_gradient_step(x, v, 0.2, 0.0)
        np.testing.assert_allclose(energy(x, v), e0, rtol=1e-9)

    def test_bdf2_second_order(self):
        x0, v0, c, T = 0.5, 0.3, 0.15, 2.0
        ref = solve_ivp(lambda t, y: (y[1], -c*y[1] - force(y[0])), (0, T), (x0, v0),
                        method="DOP853", rtol=1e-12, atol=1e-14).y[:, -1]
        errs = [abs(integrators_extra.bdf2_integrate(x0, v0, T, h, c)[0] - ref[0])
                for h in (0.02, 0.01)]
        self.assertAlmostEqual(np.log2(errs[0]/errs[1]), 2.0, delta=0.3)


class Beam(unittest.TestCase):
    def test_force_tangent_fd_and_equilibria(self):
        m = beam2mode.Beam2Mode(P=2.0, alpha=1.0)
        q = np.array([0.3, -0.2])
        f, K = m.force_tangent(q)
        e = 1e-6
        for j in range(2):
            dq = np.zeros(2); dq[j] = e
            np.testing.assert_allclose((m.potential(q+dq) - m.potential(q-dq))/(2*e), f[j], rtol=1e-6)
            np.testing.assert_allclose((m.force_tangent(q+dq)[0] - m.force_tangent(q-dq)[0])/(2*e),
                                       K[:, j], rtol=1e-6, atol=1e-8)
        for qe in m.equilibria():
            np.testing.assert_allclose(m.force_tangent(qe)[0], 0, atol=1e-12)
            self.assertTrue(np.all(np.linalg.eigvalsh(m.force_tangent(qe)[1]) > 0))

    def test_mdof_genalpha_converges(self):
        m = beam2mode.Beam2Mode()
        M, C = np.eye(2), np.diag(m.damping)
        q0, v0, T = np.array([[0.4, 0.1]]), np.array([[0.0, 0.3]]), 3.0
        def rhs(t, y):
            f, _ = m.force_tangent(y[:2]); return np.r_[y[2:], -C @ y[2:] - f]
        ref = solve_ivp(rhs, (0, T), np.r_[q0[0], v0[0]], method="DOP853",
                        rtol=1e-12, atol=1e-14).y[:2, -1]
        p = integrators.parameters("generalized-alpha", .6)
        errs = [np.linalg.norm(mdof_genalpha.integrate_to(q0, v0, T, h, M, C, m, p)[0][0] - ref)
                for h in (0.02, 0.01)]
        self.assertAlmostEqual(np.log2(errs[0]/errs[1]), 2.0, delta=0.3)


# Arch tests moved to tests/test_new_modules.py (batched FE API, 2026-11-24)


class Diagnostic(unittest.TestCase):
    def test_evaluate(self):
        cases = [{"true_error": t, "disagreement": t/2} for t in (0.01, 0.03, 0.08, 0.2)]
        out = diagnostic.evaluate(cases)
        self.assertAlmostEqual(out["spearman"], 1.0)
        self.assertEqual(out["hit_rate"], 1.0)


if __name__ == "__main__":
    unittest.main()
