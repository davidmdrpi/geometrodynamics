"""Small numerical regressions for the clocks, physical ports and negative controls."""
import json
import unittest
from pathlib import Path

import numpy as np

from experiments.closure_ledger.antipodal_network import clock, gravity, network, receiver, refinement

ARCHIVE = (Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger" /
           "runs" / "20261117_antipodal_network")


class AntipodalClockNetworkTests(unittest.TestCase):
    def test_clock_selection_and_action_scaling(self):
        result = clock.verify()
        self.assertEqual(result["modes_n_m"][:5], [[2, 0], [3, 1], [7, 3], [18, 8], [47, 21]])
        self.assertAlmostEqual(result["energy_ratio_at_0p37"], 0.37**2)

    def test_pure_tt_seed_fails_quadratic_condition(self):
        result = gravity.verify()
        self.assertAlmostEqual(result["flat_cover_quadratic_TT_mean"], 9.0)

    def test_free_tube_exchange_and_port_work(self):
        result = network.run(layers=4, stop=0.8)
        self.assertAlmostEqual(result["initial_tube_energy"], 1.0)
        self.assertLess(result["max_tube_port_balance_error"], 1e-10)
        self.assertGreater(result["max_cumulative_out_at_ports"][1], 0.1)

    def test_exact_static_tube_and_passive_receiver(self):
        (tube, exterior), _, _ = refinement.domain(0, full_ports=True)
        np.testing.assert_allclose(
            tube(0.0),
            receiver.rim_area() / 0.95 * np.array([[1, -1], [-1, 1]]),
            rtol=1e-10, atol=1e-12,
        )
        measurement = receiver.solve(tube(1), exterior(1), 1, 0.2375)
        self.assertGreater(measurement["tube_power_to_B"], 0)
        self.assertAlmostEqual(measurement["receiver_power"], measurement["source_power"])

    def test_archives_preserve_failed_predictions(self):
        roots = json.loads((ARCHIVE / "refinement.json").read_text())
        self.assertEqual([r["cover_tetrahedra"] for r in roots], [3072, 24576, 196608])
        self.assertLess(roots[0]["coarse_clock_cross_reaction"], 1e-9)
        self.assertGreater(roots[-1]["coarse_clock_cross_reaction"],
                           roots[-1]["coarse_zero_clock_cross_reaction"])
        packet = json.loads((ARCHIVE / "transport.json").read_text())
        self.assertEqual(packet[0]["columns"][2], "net_out_B")
        self.assertGreater(packet[0]["history"][-1][2], 0.1)
        self.assertLess(packet[-1]["history"][-1][2], -0.1)


if __name__ == "__main__":
    unittest.main()
