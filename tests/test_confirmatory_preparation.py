import unittest

import numpy as np

from basinlab.experiments.confirmatory_preparation import _differences, smoke_population
from basinlab.model import energy


class PreparationTests(unittest.TestCase):
    def test_population_is_filtered_24_by_24_grid(self):
        x, v = smoke_population()
        self.assertGreater(len(x), 0)
        self.assertLess(len(x), 24*24)
        e = energy(x, v)
        self.assertTrue(np.all((e > 0.6) & (e < 2)))

    def test_unresolved_labels_are_not_silently_discarded(self):
        result = _differences(np.array([1, 0, -1]), np.array([-1, 1, -1]))
        self.assertEqual(result, {"different_labels": 2, "unresolved_a": 1,
                                  "unresolved_b": 0})


if __name__ == "__main__":
    unittest.main()
