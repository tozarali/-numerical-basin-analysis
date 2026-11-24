import unittest

import numpy as np

from basinlab.experiments.dsep_exploratory import crossings, loss_ratio, reference_events
from basinlab.integrators import parameters


class DsepExploratoryTests(unittest.TestCase):
    def test_reference_capture_prevents_cycle_at_requested_damping(self):
        for delta in (0.02, 0.05, 0.1):
            _, ref = reference_events(delta, 0.15)
            self.assertEqual(ref["cycle_status"], "no_cycle_captured_below_barrier")
            self.assertIsNone(ref["cycle_time"])
            self.assertLess(ref["capture_time"], ref["turn_time"])
            self.assertLess(abs(ref["final_energy_balance_residual"]), 1e-9)

    def test_low_damping_can_complete_real_cycle(self):
        _, ref = reference_events(0.1, 0.001)
        self.assertEqual(ref["cycle_status"], "completed")
        self.assertGreater(ref["cycle_time"], ref["turn_time"])

    def test_matched_time_ratio_converges(self):
        sol, ref = reference_events(0.05, 0.15)
        p = parameters("generalized-alpha", 0.5)
        coarse = loss_ratio(sol, ref["e0"], ref["v0"], ref["turn_time"], 0.4, 0.15, p)
        fine = loss_ratio(sol, ref["e0"], ref["v0"], ref["turn_time"], 0.005, 0.15, p)
        self.assertLess(abs(fine - 1), abs(coarse - 1))
        self.assertLess(abs(fine - 1), 1e-4)

    def test_all_crossings_reported(self):
        np.testing.assert_allclose(crossings([0, 0.5, 1], [0.9, 1.1, 0.9]), [0.25, 0.75])


if __name__ == "__main__":
    unittest.main()
