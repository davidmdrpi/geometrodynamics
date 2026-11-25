import json
import unittest
from pathlib import Path

import numpy as np

from experiments.closure_ledger.antipodal_network import two_detector_records as records


class TwoDetectorRecordsTests(unittest.TestCase):
    def test_pulse_writes_one_or_two_actual_records(self):
        q, p = records.preparation(1, 363, 0.)
        rows, codes, residual = records.evolve([[0., 0.], [1., 0.], [1., 1.]], q, p, 0., .025)
        np.testing.assert_array_equal(codes[:, 0], [0, 1, 3])
        np.testing.assert_allclose(rows[1]["mean_position"], [1., 0.], atol=1e-7)
        self.assertLess(residual, 1e-5)

    def test_first_exit_threshold_is_strict(self):
        q = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
        barrier = records.potential(records.SADDLE)
        thresholds = np.array([barrier, barrier, barrier, barrier+1])
        np.testing.assert_array_equal(records.certify(q, thresholds-.01), [0, 1, 2, 3])
        np.testing.assert_array_equal(records.certify(q, thresholds), [4, 4, 4, 4])

    def test_symmetric_motion_is_independent_of_link_inertia(self):
        q, p = records.preparation(1, 363, 0.)
        rows, codes, _ = records.evolve([[.5, .5], [1., 1.]], q, p, .25, .025)
        np.testing.assert_array_equal(codes[:, 0], [0, 3])
        self.assertLess(rows[0]["final_energy_range"][1], 1e-8)

    def test_energy_cap_and_independent_balanced_ports(self):
        e = .06
        largest_q = np.sqrt(2*e/18.)
        bound = records.potential(-largest_q)+e
        self.assertLess(2*bound + records.PULSE_ENERGY, 2.)
        self.assertLess(bound + records.PULSE_ENERGY/2, records.potential(records.SADDLE))
        q, p = records.preparation(32, 363, e)
        _, labels, _ = records.evolve([[.5, .5]], q, p, 0., .025)
        self.assertTrue(np.all(labels == 0))

    def test_archived_gates_use_every_prepared_trial(self):
        folder = Path(__file__).resolve().parents[1] / "experiments/closure_ledger/runs"
        archive = json.loads((folder / "20261125_two_detector_records/probe.json").read_text())
        for run in archive["runs"].values():
            errors = []
            all_exclusive = True
            self.assertEqual(len(run["rows"]), len(records.THETAS))
            for theta, row in zip(records.THETAS, run["rows"]):
                counts = row["counts"]
                self.assertEqual(set(counts), set(records.LABELS))
                self.assertTrue(all(isinstance(n, int) and n >= 0 for n in counts.values()))
                self.assertEqual(sum(counts.values()), run["count_per_angle"])
                p_right = counts["right_only"]/run["count_per_angle"]
                self.assertAlmostEqual(row["theta"], theta)
                self.assertAlmostEqual(row["right_only_probability"], p_right)
                errors.append(abs(p_right - np.cos(theta)**2))
                all_exclusive &= counts["right_only"] + counts["left_only"] == run["count_per_angle"]
            self.assertAlmostEqual(run["max_born_error"], max(errors))
            self.assertEqual(run["born_window"], max(errors) <= .01)
            self.assertEqual(run["all_exclusive"], all_exclusive)
        coarse = [r["counts"] for r in archive["runs"]["primary"]["rows"]]
        fine = [r["counts"] for r in archive["runs"]["refined"]["rows"]]
        self.assertEqual(archive["checks"]["refined_counts_equal"], coarse == fine)
