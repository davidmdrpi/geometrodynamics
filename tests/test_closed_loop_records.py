import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from experiments.closure_ledger.antipodal_network import closed_loop_records as loop


class ClosedLoopTests(unittest.TestCase):
    def test_asymmetric_echo_against_method_of_steps(self):
        # Reconstruct the incident first echo from an adaptive, dense solution
        # of the open endpoint equation, independent of the delay buffer.
        q, p = loop.base.preparation(2, 363, .02)
        q, p = q[:1], p[:1]
        weights = [.8, .2]
        v = loop.base.inverse_mass(p, .25)
        y0 = np.r_[q[0], v[0]]
        initial = lambda t: np.sqrt(weights)*loop.base.incoming(t)

        def solve(begin, end, y, incident):
            def rhs(t, state):
                x, v = state[:2], state[2:]
                force = loop.base.GAMMA*(2*incident(t)-v)-loop.base.slope(x)
                return np.r_[v, loop.base.inverse_mass(force, .25)]
            out = solve_ivp(rhs, (begin, end), y, method="DOP853", dense_output=True,
                            rtol=1e-11, atol=1e-13, max_step=.1)
            self.assertTrue(out.success)
            return out

        before = solve(0., 4., y0, initial)
        def echo(t):
            return (before.sol(t-4)[2:]-initial(t-4))[::-1]
        after = solve(4., 8., before.y[:, -1], echo)
        errors = []
        for dt in (.025, .0125):
            row, track, _ = loop.evolve(weights, q, p, dt=dt)
            t = np.arange(round(4/dt), round(8/dt)+1)*dt
            errors.append(np.max(abs(track[round(4/dt):round(8/dt)+1]-after.sol(t)[:2].T)))
            self.assertLess(row["max_energy_residual"], 1e-8)
        self.assertLess(errors[1], errors[0]/3)
        self.assertLess(errors[1], .002)

    def test_crossing_invalidates_persistence(self):
        low = np.array([[.9, .9, -.1, .8, -.1], [-.1, .9, .9, -.1, -.1]])
        high = np.array([[1.1, 1.1, .1, 1.1, .1], [.1, 1.1, 1.1, .1, .1]])
        np.testing.assert_array_equal(loop.persistent(low, high), [1, 3, 2, 4, 0])

    def test_dark_and_symmetric_ports(self):
        cold = np.zeros((1, 2))
        for w in (0., .5):
            row, trace, _ = loop.evolve([w, w], cold, cold)
            np.testing.assert_array_equal(trace[:, 0], trace[:, 1])
            self.assertEqual(row["exclusive"], 0.)
            if w == 0:
                np.testing.assert_equal(trace, 0.)

    def test_archive_from_counts(self):
        folder = Path(__file__).resolve().parents[1]/"experiments/closure_ledger/runs"
        report = json.loads((folder/"20261125_closed_loop_records/probe.json").read_text())
        for name in ("primary", "fine", "repeat"):
            rows = report[name]
            self.assertEqual([r["index"] for r in rows], list(range(9)))
            expected_errors, exclusives = [], []
            for i, row in enumerate(rows):
                counts = row["persistent_counts"]
                self.assertEqual(set(counts), set(loop.LABELS))
                self.assertTrue(all(isinstance(v, int) and v >= 0 for v in counts.values()))
                size = 1024 if name == "repeat" else 512
                self.assertEqual(sum(counts.values()), size)
                prob = counts["right_only"]/size
                expected_errors.append(abs(prob-np.cos(loop.base.THETAS[i])**2))
                exclusives.append((counts["right_only"]+counts["left_only"])/size)
                self.assertEqual(prob, row["right_only"])
                self.assertEqual(exclusives[-1], row["exclusive"])
            open_row = report["repeat_open_balanced"] if name == "repeat" else report["open"][4]
            baseline = sum(open_row["persistent_counts"][k] for k in ("right_only", "left_only"))
            baseline /= open_row["launches"]
            self.assertEqual(report["gates"][name], dict(
                above_half_and_open=exclusives[4] > max(.5, baseline),
                near_unit_all_angles=min(exclusives) >= .99, born_window=max(expected_errors) <= .01,
                max_born_error=max(expected_errors)))
        for ctl in report["length_controls"]:
            delay = ctl["closed"]["delay"]
            self.assertLess(ctl["pre_echo_max"], 1e-10)
            self.assertGreater(ctl["first_deviation"], delay)
            self.assertLessEqual(ctl["first_deviation"], delay+ctl["closed"]["dt"]+1e-12)
        all_rows = [row for key in ("primary", "fine", "repeat", "open", "extra_refinement")
                    for row in report[key]]
        all_rows.extend(report[key] for key in ("repeat_open_balanced", "cold_balanced",
                                                "cold_dark", "balanced_no_mutual_inertia"))
        all_rows.extend(ctl[key] for ctl in report["length_controls"] for key in ("closed", "open"))
        budget = max(r["max_energy_residual"] for r in all_rows)
        self.assertEqual(report["checks"]["max_energy_residual"], budget)
        self.assertLess(budget, 1e-8)
        for coarse, fine in zip(report["primary"], report["fine"]):
            self.assertEqual(coarse["persistent_counts"], fine["persistent_counts"])

    def test_fluence_archive_preserves_snapshot_and_persistent_counts(self):
        path = (Path(__file__).resolve().parents[1] / "experiments/closure_ledger/runs"
                / "20261125_closed_loop_records" / "fluence.json")
        rows = json.loads(path.read_text())
        cases = {(row["energy_ratio"], row["closed"], row["delay"]) for row in rows}
        expected = {(1+i/8, closed, delay) for i in range(11)
                    for closed, delay in ((False, 4), (True, 3), (True, 4), (True, 6))}
        self.assertEqual(cases, expected)
        self.assertEqual(len(cases), len(rows))
        for row in rows:
            singles = row["persistent"]["right_only"] + row["persistent"]["left_only"]
            self.assertEqual(sum(row["persistent"].values()), 512)
            self.assertEqual(sum(row["endpoint"].values()), 512)
            self.assertEqual(row["exclusive"], singles/512)
            if row["closed"]:
                self.assertEqual(singles, 0)
            if row["energy_ratio"] == 1.75 and row["delay"] == 4:
                end = row["endpoint"]
                self.assertEqual(end["right_only"]+end["left_only"], 18 if row["closed"] else 236)
