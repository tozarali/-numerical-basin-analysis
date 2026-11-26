import unittest

import numpy as np

from basinlab import integrators, linear_analysis


class LinearExploratoryTests(unittest.TestCase):
    def test_matrix_matches_direct_linear_update_for_all_families(self):
        x, v, a = 0.3, -0.2, 0.7
        h, xi = 0.4, 0.1
        for family in ("generalized-alpha", "HHT", "WBZ", "Newmark"):
            p, _ = integrators.parameters_for_spectral_radius(family, 0.8)
            am, af, gamma, beta = p
            xp = x + h*v + h*h*(0.5-beta)*a
            vp = v + h*(1-gamma)*a
            a1 = (-am*a - 2*xi*((1-af)*vp + af*v)
                  - ((1-af)*xp + af*x)) / (
                      1-am + (1-af)*(2*xi*gamma*h + beta*h*h))
            expected = [xp + beta*h*h*a1, vp + gamma*h*a1, a1]
            np.testing.assert_allclose(
                linear_analysis.amplification_matrix(p, h, xi) @ [x, v, a],
                expected, rtol=1e-13, atol=1e-13)

    def test_damped_trapezoid_matches_bilinear_complex_root(self):
        p = integrators.parameters("generalized-alpha", 1.0)
        h, xi = 0.4, 0.1
        continuous_root = -xi + 1j*np.sqrt(1-xi*xi)
        expected = (1+h*continuous_root/2)/(1-h*continuous_root/2)
        actual = linear_analysis.principal_pair(p, h, xi)
        self.assertAlmostEqual(actual.real, expected.real, places=12)
        self.assertAlmostEqual(actual.imag, expected.imag, places=12)

    def test_zero_and_nonoscillatory_cases_have_no_period(self):
        p = integrators.parameters("generalized-alpha", 0.5)
        for h, xi in ((0, 0), (0.4, 1), (0.4, 2)):
            self.assertIsNone(linear_analysis.principal_pair(p, h, xi))
            self.assertTrue(np.isnan(linear_analysis.dissipation_and_period(p, h, xi)["period_error"]))
            self.assertTrue(np.isnan(linear_analysis.per_cycle_energy_ratio(p, h, xi)))
        for h, xi in ((-0.1, 0), (0.1, -0.1), (float("nan"), 0)):
            with self.assertRaises(ValueError):
                linear_analysis.amplification_matrix(p, h, xi)


if __name__ == "__main__":
    unittest.main()
