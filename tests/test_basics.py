import unittest

import numpy as np

from basinlab.integrators import (
    integrate, parameters, parameters_for_spectral_radius,
    spectral_radius_infinity, step,
)
from basinlab.experiments.energy_loss_retest import exact_step_grid
from basinlab.metrics import boundary_distance, entropy, mismatch
from basinlab.model import classify, energy
from basinlab.reference import adaptive, rk4


class BasinTests(unittest.TestCase):
    def test_capture_requires_energy_below_barrier(self):
        np.testing.assert_array_equal(classify([-1, 1, 0, 1], [0, 0, 0, 1]),
                                      [-1, 1, 0, 0])
        self.assertAlmostEqual(float(energy(0, 0)), 0.25)

    def test_entropy_and_missing_boxes(self):
        labels = np.array([[-1, 1, -1, -1], [-1, 1, -1, -1]])
        result = entropy(labels, 2)
        self.assertAlmostEqual(result["Sb"], np.log(2) / 2)
        self.assertAlmostEqual(result["Sbb"], np.log(2))
        labels[0, 2] = 0
        self.assertEqual(entropy(labels, 2)["excluded_boxes"], 1)
        self.assertEqual(mismatch(labels, labels)["excluded"], 1)
        self.assertAlmostEqual(mismatch(np.array([-1, 1]), np.array([1, 1]))["fraction"], .5)

    def test_boundary_distance_known_shift(self):
        a = np.tile([-1, -1, 1, 1], (3, 1))
        b = np.tile([-1, 1, 1, 1], (3, 1))
        result = boundary_distance(a, b, (1, 1))
        self.assertAlmostEqual(result["mean_symmetric"], 1)
        self.assertAlmostEqual(result["hausdorff"], 1)

    def test_reference_and_integrator_converge(self):
        x, v = np.array([1.2, -0.3]), np.array([0.3, 0.7])
        exact = adaptive(x, v, 4, .15, rtol=1e-12, atol=1e-14)
        coarse = integrate(x, v, 4, .2, .15, parameters("generalized-alpha", 1))
        fine = integrate(x, v, 4, .1, .15, parameters("generalized-alpha", 1))
        err = lambda state: np.linalg.norm(np.array(state) - np.array(exact))
        self.assertLess(err(fine), err(coarse) / 3)
        self.assertLess(err(rk4(x, v, 4, .01, .15)), 1e-7)

    def test_exact_step_grid_uses_short_final_step(self):
        x, v, damping, h, end = np.array([.8]), np.array([.4]), .15, .06, .1
        p = parameters("generalized-alpha", .4)
        actual = exact_step_grid(x, v, end, h, damping, p)
        a = -damping * v - (-x + x**3)
        x, v, a = step(x, v, a, .06, damping, p)
        expected = step(x, v, a, .04, damping, p)[:2]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)

    def test_bad_inputs_are_not_silently_accepted(self):
        with self.assertRaises(ValueError):
            integrate([1], [0], 1, .3, .15, parameters("generalized-alpha", 1))
        with self.assertRaises(ValueError):
            boundary_distance(np.array([[0, 1]]), np.array([[-1, 1]]), (1, 1))
        with self.assertRaises(ValueError):
            parameters("unknown", 0)

    def test_spectral_radius_family_matching(self):
        for family in ("HHT", "WBZ", "Newmark"):
            p, measured = parameters_for_spectral_radius(family, .8)
            self.assertAlmostEqual(measured, .8, places=10)
            self.assertAlmostEqual(spectral_radius_infinity(p), .8, places=10)
        with self.assertRaises(ValueError):
            parameters_for_spectral_radius("HHT", .4)


if __name__ == "__main__":
    unittest.main()
