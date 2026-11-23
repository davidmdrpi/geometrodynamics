import unittest

import numpy as np

from experiments.closure_ledger.antipodal_network import timelike_current as current


class TimelikeCurrentTests(unittest.TestCase):
    def test_full_tensor_matches_modal_flux_including_shift_and_trace(self):
        rng = np.random.default_rng(359)
        u, pi, grad = rng.normal(size=(3,2)) + 1j*rng.normal(size=(3,2))
        for a in (0.,-.3,.3):
            _,_,energy,flux = current.currents(u[:,None],pi[:,None],grad[:,None],8.,a)
            expected = current.direct_fiber_contraction(u,pi,grad,8.,a)
            np.testing.assert_allclose([energy[0],flux[0]],expected,rtol=1e-13,atol=1e-13)

    def test_exact_midpoint_and_split_potential(self):
        self.assertEqual(current.exact_control(offset=0.)["energy_velocity"],0.)
        for gap in (.2,-.2,.002):
            case = current.exact_control(offset=gap/2)
            self.assertGreater(case["energy_density"],0.)
            self.assertEqual(case["norm_velocity"],0.)
            self.assertAlmostEqual(case["energy_velocity"],case["analytic"],places=14)
            self.assertGreater(abs(case["energy_velocity"]),1e-6)

    def test_reference_field_and_bulk_evolution(self):
        weight = .85
        lift = current.archived_field_check(1.,8.,weight)
        result = current.run(8.,.02,512,weight)
        drive = next(row for row in result["snapshots"] if row["t"]==1.)
        self.assertGreater(drive["rms_energy_minus_norm"],.005)
        self.assertLess(abs(lift["rms_energy_minus_358"]-drive["rms_energy_minus_norm"]),1e-4)
        self.assertLess(result["free_energy_relative_drift"],1e-9)
        self.assertLess(result["drive_budget_relative"],5e-6)
        self.assertTrue(all(s["midpoint_energy"]>0 for s in result["snapshots"]))

    def test_number_current_control_and_transport(self):
        from experiments.closure_ledger.antipodal_network import branch_current_shift

        control = current.exact_control()
        self.assertEqual(control["number_velocity"], 0.)
        self.assertGreater(control["number_density"], 1.)
        run = branch_current_shift.run(8., .02, 1024)
        result = run["comparison"]
        self.assertLess(result["max_number_error_from_final_mass"], 5e-6)
        self.assertLess(result["max_relative_right_work_residual"], 5e-6)
        self.assertGreater(result["max_shift_from_same_pair"], 0.0009)
        self.assertAlmostEqual(run["rows"][4]["transported"]["energy"], .5, places=8)
