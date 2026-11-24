import os
import tempfile
import unittest

import numpy as np

from basinlab import beam2mode, mdof_genalpha, metrics_ext, reference_mdof
from basinlab.experiments import _common, confirmatory_family_comparison as cfc
from basinlab.fe_arch import Arch
from basinlab.integrators_extra import discrete_gradient_step
from basinlab.model import classify, energy


class MetricsExt(unittest.TestCase):
    def test_signed_shift_and_net_transfer(self):
        n = 40
        ref = np.where(np.arange(n)[None, :] < 20, -1, 1).repeat(n, axis=0).astype(np.int8)
        cand = np.where(np.arange(n)[None, :] < 17, -1, 1).repeat(n, axis=0).astype(np.int8)
        s = metrics_ext.signed_boundary_shift(ref, cand, (1.0, 1.0))
        self.assertAlmostEqual(s["mean_abs"], 3.0)
        self.assertGreater(s["mean_signed"], 0)          # + basin grew
        nt = metrics_ext.net_basin_transfer(ref, cand)
        self.assertAlmostEqual(nt["net"], 3/40)

    def test_uncertainty_exponent_smooth_boundary(self):
        rng = np.random.default_rng(0)
        ics = rng.uniform(-1, 1, size=(20000, 2))
        cls = lambda p: np.where(p[:, 0] > 0, 1, -1)
        out = metrics_ext.uncertainty_exponent(cls, ics, [0.01, 0.02, 0.04, 0.08])
        self.assertAlmostEqual(out["alpha"], 1.0, delta=0.15)


class DiscreteGradient(unittest.TestCase):
    def test_damped_energy_decreases_and_grid_converges(self):
        x, v = _common.grid(30)
        e0 = energy(x, v)
        x1, v1 = discrete_gradient_step(x, v, 0.8, 0.15)
        self.assertTrue(np.all(energy(x1, v1) <= e0 + 1e-12))


class FEArch(unittest.TestCase):
    def setUp(self):
        self.a = Arch(n_el=8)

    def test_rigid_body_motion_is_force_free(self):
        a = Arch(n_el=8, imperfection=0.05)
        th = 0.3
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        moved = a.X @ R.T + np.array([0.2, -0.1])
        d = np.zeros(a.ndof)
        d[0::3], d[1::3], d[2::3] = (moved - a.X)[:, 0], (moved - a.X)[:, 1], th
        np.testing.assert_allclose(a.internal_force(d), 0, atol=1e-8)

    def test_tangent_and_energy_consistency(self):
        a = Arch(n_el=6, imperfection=0.05, end_shortening=0.002)
        rng = np.random.default_rng(1)
        q = 0.01*rng.normal(size=len(a.free))
        f, K = a.force_tangent(q)
        eps = 1e-6
        for j in range(0, len(q), 3):
            dq = np.zeros_like(q); dq[j] = eps
            np.testing.assert_allclose((a.force(q + dq) - a.force(q - dq))/(2*eps), K[:, j],
                                       rtol=1e-5, atol=1e-4)
            self.assertAlmostEqual((a.potential(q + dq) - a.potential(q - dq))/(2*eps), f[j],
                                   delta=1e-5*max(1, abs(f[j])))

    def test_euler_buckling_and_post_buckling(self):
        a = Arch(n_el=16)
        dcr = a.euler_shortening
        a.end_shortening = 0.97*dcr
        self.assertGreater(a.flat_state()[2][0], 0)
        a.end_shortening = 1.05*dcr
        self.assertLess(a.flat_state()[2][0], 0)
        a.end_shortening = 1.6*dcr
        st = a.buckled_states()
        wm = st[1][0][a.mid_index]
        pred = np.sqrt(4*0.6*dcr/np.pi**2)
        self.assertAlmostEqual(wm/pred, 1.0, delta=0.03)
        self.assertAlmostEqual(st[-1][0][a.mid_index], -wm, places=8)
        self.assertGreater(st[1][2][0], 0)

    def test_constrained_midspan_and_classify(self):
        a = Arch(n_el=8)
        a.end_shortening = 1.6*a.euler_shortening
        up = a.buckled_states()[1][0]
        wm = up[a.mid_index]
        Q, _ = a.constrained_midspan(np.array([-wm, 0.5*wm, wm]))
        np.testing.assert_allclose(Q[:, a.mid_index], [-wm, 0.5*wm, wm])
        np.testing.assert_allclose(Q[2], up, atol=1e-7)
        barrier = a.potential(a.flat_state()[0])
        lab = a.classify(Q, np.zeros_like(Q), barrier)
        self.assertEqual(lab.tolist(), [-1, 1, 1])


class Beam2ModeReduction(unittest.TestCase):
    def test_q2_zero_reduces_to_sdof(self):
        from basinlab.integrators import parameters
        m = beam2mode.Beam2Mode(P=2.0, alpha=1.0, damping=(0.15, 0.15))
        X, V = _common.grid(10)
        q = np.column_stack((X.ravel(), np.zeros(X.size)))
        v = np.column_stack((V.ravel(), np.zeros(X.size)))
        p = parameters("generalized-alpha", 0.5)
        qT, vT, nf = mdof_genalpha.integrate_to_masked(q, v, 20.0, 0.4, np.eye(2),
                                                        np.diag(m.damping), m, p)
        sdof = classify(*_common.alpha_state(X.ravel(), V.ravel(), 20.0, 0.4, 0.15, p))
        self.assertEqual(nf, 0)
        np.testing.assert_array_equal(m.classify(qT, vT), sdof)

    def test_mdof_rk4_matches_sdof_rk4(self):
        from basinlab.reference import rk4 as rk4_sdof
        m = beam2mode.Beam2Mode()
        q = np.array([[0.3, 0.0], [-1.2, 0.0]]); v = np.array([[0.5, 0.0], [0.1, 0.0]])
        qT, vT = reference_mdof.rk4(q, v, 2.0, 0.01, np.ones(2), np.diag(m.damping), m)
        xs, vs = rk4_sdof(q[:, 0], v[:, 0], 2.0, 0.01, 0.15)
        np.testing.assert_allclose(qT[:, 0], xs, atol=1e-12)


class Confirmatory(unittest.TestCase):
    def test_unique_crossing(self):
        r = np.array([0, .1, .2, .3])
        self.assertAlmostEqual(cfc.unique_crossing(r, [1.2, 1.1, .9, .8]), 0.15)
        self.assertIsNone(cfc.unique_crossing(r, [1.2, .9, 1.1, .8]))
        self.assertIsNone(cfc.unique_crossing(r, [1.2, 1.1, 1.05, 1.01]))

    def test_grids_match_preregistered_counts(self):
        counts = {f: len(g) for f, g in cfc.rho_grids().items()}
        self.assertEqual(counts, {"generalized-alpha": 41, "HHT": 21, "WBZ": 39, "Newmark": 17})
        for f, g in cfc.rho_grids().items():
            for r in g:
                self.assertIsNotNone(_common.family_parameters(f, r), (f, r))

    def test_draft_protocol_is_locked(self):
        self.assertFalse(cfc.lock_status()["ok"])


class Parts(unittest.TestCase):
    def test_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = _common.Parts(d)
            self.assertFalse(p.has("k"))
            p.save("k", {"a": 1}, x=np.arange(3))
            meta, arr = p.load("k")
            self.assertEqual(meta, {"a": 1})
            np.testing.assert_array_equal(arr["x"], np.arange(3))


if __name__ == "__main__":
    unittest.main()


class Runner(unittest.TestCase):
    def test_registry_dependencies_and_smoke_flag(self):
        from basinlab import run_all
        names = [e.name for e in run_all.REGISTRY]
        for i, e in enumerate(run_all.REGISTRY):
            for d in e.deps:
                self.assertIn(d, names[:i], f"{e.name} depends on later/unknown {d}")
            self.assertTrue(run_all.source_file(e).exists(), e.module)
        self.assertTrue(run_all.BY_NAME["confirmatory_smoke"].always_quick)
        self.assertFalse(run_all.BY_NAME["confirmatory_family_comparison"].always_quick)
