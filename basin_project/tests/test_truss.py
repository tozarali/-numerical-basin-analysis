import unittest

import numpy as np

from basinlab.truss import Truss, dop853, generalized_alpha, rk4


class TrussTests(unittest.TestCase):
    def setUp(self):
        self.truss = Truss()

    def test_force_and_tangent_are_potential_derivatives(self):
        q = np.array([.22, .48])
        force, tangent = self.truss.force_tangent(q)
        eps = 1e-6
        for i in range(2):
            perturb = np.eye(2)[i] * eps
            derivative = (self.truss.potential(q+perturb)
                          - self.truss.potential(q-perturb))/(2*eps)
            self.assertAlmostEqual(force[i], derivative, places=8)
            fplus, _ = self.truss.force_tangent(q+perturb)
            fminus, _ = self.truss.force_tangent(q-perturb)
            np.testing.assert_allclose((fplus-fminus)/(2*eps), tangent[:, i],
                                       rtol=1e-7, atol=1e-9)

    def test_equilibria_and_barrier_classification(self):
        q = np.array([[0, .7], [0, -.7], [0, 0], [.2, .7]])
        v = np.zeros_like(q)
        f, _ = self.truss.force_tangent(q)
        np.testing.assert_allclose(f[:3], 0, atol=1e-15)
        np.testing.assert_array_equal(self.truss.classify(q[:3], v[:3]), [1, -1, 0])
        self.assertAlmostEqual(self.truss.potential(q[2]), self.truss.barrier)
        yzero = np.column_stack((np.linspace(-3, 3, 601), np.zeros(601)))
        self.assertGreaterEqual(float(np.min(self.truss.potential(yzero))),
                                self.truss.barrier-1e-14)

    def test_rk4_and_generalized_alpha_converge_to_dop853(self):
        q, v = np.array([.2, .4]), np.array([.1, .6])
        reference = np.asarray(dop853(self.truss, q, v, 4))
        reference_rk = np.asarray(rk4(self.truss, q, v, 4, .01))
        np.testing.assert_allclose(reference_rk, reference, rtol=1e-7, atol=1e-8)
        coarse = np.asarray(generalized_alpha(self.truss, q, v, 4, .2, 1))
        fine = np.asarray(generalized_alpha(self.truss, q, v, 4, .1, 1))
        self.assertLess(np.linalg.norm(fine-reference),
                        np.linalg.norm(coarse-reference)/3)

    def test_constrained_model_matches_symmetric_full_model(self):
        q, v = np.array([[0., .2], [0., -.6]]), np.array([[0., .8], [0., -.4]])
        for solve in (rk4, generalized_alpha):
            args = (self.truss, q, v, 4, .1)
            kwargs = {"rho": .4} if solve is generalized_alpha else {}
            np.testing.assert_allclose(solve(*args, constrained=True, **kwargs),
                                       solve(*args, constrained=False, **kwargs),
                                       rtol=0, atol=1e-12)

    def test_invalid_states_raise(self):
        with self.assertRaises(ArithmeticError):
            self.truss.force_tangent(np.array([1., 0.]))
        with self.assertRaises(ValueError):
            rk4(self.truss, [.1, .2], [0, 0], 1, .1, constrained=True)
        with self.assertRaises(ValueError):
            generalized_alpha(self.truss, [0, .5], [0, 0], 1, 0, .5)

    def test_reference_energy_decreases_with_damping(self):
        q, v = np.array([.6, .4]), np.array([.0, .5])
        initial = self.truss.energy(q, v)
        q, v = rk4(self.truss, q, v, 4, .005)
        self.assertLess(float(self.truss.energy(q, v)), float(initial))


if __name__ == "__main__":
    unittest.main()
