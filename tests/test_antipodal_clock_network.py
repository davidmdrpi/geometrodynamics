"""Small numerical regressions for the clocks, physical ports and negative controls."""
import json
import unittest
from pathlib import Path

import numpy as np

from experiments.closure_ledger.antipodal_network import clock, gravity, network, receiver, refinement, round_exterior

ARCHIVE = (Path(__file__).resolve().parents[1] / "experiments" / "closure_ledger" /
           "runs" / "20261117_antipodal_network")


class AntipodalClockNetworkTests(unittest.TestCase):
    def test_clock_selection_and_action_scaling(self):
        result = clock.verify()
        self.assertEqual(result["modes_n_m"][:5], [[2, 0], [3, 1], [7, 3], [18, 8], [47, 21]])
        self.assertAlmostEqual(result["energy_ratio_at_0p37"], 0.37**2)

    def test_clock_generator_composition_and_fixed_orbits(self):
        controls = [clock.clock_control(q) for q in range(4)]
        self.assertEqual([(c["descends"], c["fixed_orbit"]) for c in controls],
                         [(True, True), (False, None), (True, False), (False, None)])

    def test_round_chamber_reproduces_a_spherical_harmonic(self):
        # Y=x0²-x1² solves -ΔY=8Y on S³ and vanishes on the exchange plane.
        from scipy.sparse.linalg import spsolve
        radius = 0.31
        errors = []
        for n in (8, 16, 32):
            k, m = round_exterior.chamber(n, radius)
            alpha = np.linspace(0, np.pi/2, n+1)
            limit = np.arctan2(1, np.cos(alpha))
            r = radius + np.linspace(0, 1, n+1)[:, None]*(limit-radius)
            exact = (np.cos(r)**2 - (np.sin(r)*np.cos(alpha))**2).ravel()
            interior = np.arange(n+1, n*(n+1))
            z = k-8*m
            residual = z @ exact
            error = spsolve(z[interior][:, interior], residual[interior])
            errors.append(np.sqrt(error @ (m[interior][:, interior] @ error)))
        self.assertTrue(all(3.7 < a/b < 4.3 for a, b in zip(errors, errors[1:])))

    def test_round_response_and_positive_receiver_work(self):
        area = receiver.rim_area()
        ext = round_exterior.RoundExterior(12, area)
        r = ext.radius
        self.assertAlmostEqual(ext.volume, np.pi**2 - 4*np.pi*(r-np.sin(r)*np.cos(r)), places=6)
        np.testing.assert_allclose(ext.ports(0) @ np.ones(2), 0, atol=1e-10)
        row = receiver.solve(round_exterior.tube_ports(area, 1), ext.ports(1), 1, 0.2375)
        self.assertGreater(row["tube_power_to_B"], 0)
        self.assertAlmostEqual(row["source_power"], row["receiver_power"])

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

    def test_finer_grids_use_the_coarse_clock_calibration(self):
        rows = json.loads((ARCHIVE / "round_exterior.json").read_text())
        by_size = {row["subdivisions"]: row["switch"] for row in rows}
        predicted = (4*by_size[64]["fitted_root"] - by_size[32]["fitted_root"]) / 3
        for n in (128, 256):
            setting = by_size[n]
            self.assertAlmostEqual(setting["fixed_frequency"], predicted, places=10)
            self.assertAlmostEqual(setting["fixed_shift"], np.pi/predicted, places=10)
            self.assertGreater(setting["zero_shift_power"], 1)
            self.assertLess(setting["fixed_power_ratio"], 1e-5)
            self.assertAlmostEqual(setting["fixed_power_ratio"],
                                   setting["fixed_dark_power"]/setting["zero_shift_power"])

    def test_archives_preserve_mesh_dependent_results(self):
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
