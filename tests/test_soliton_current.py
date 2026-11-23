import unittest

import numpy as np

from experiments.closure_ledger import soliton_current as soliton


class SolitonDiscriminatorTests(unittest.TestCase):
    def test_rest_profile_first_integral(self):
        f, slope = soliton.profile(np.linspace(-4, 4, 129))
        np.testing.assert_allclose(slope**2, f**2-2*f**4+.5*f**6, atol=1e-14)

    def test_free_boost(self):
        rows = soliton.run(n=512, dt=.002, stride=10, force=0., boost=.5)
        for row in rows:
            self.assertLess(abs(row["spectral_peak"]-.5*row["t"]), 2e-7)
            for key in ("number_path", "energy_path"):
                self.assertLess(abs(row[key]-.5*row["t"]), 5e-6)

    def test_driven_core_and_independent_charge_median(self):
        rows = soliton.run(n=512, dt=.002, stride=5)
        result = soliton.summarize(rows)
        self.assertLess(result["max_work_residual"], 5e-6)
        self.assertLess(result["median_number_path_error"], 2e-5)
        self.assertGreater(result["predictor_separation"], .0005)
        self.assertGreater(min(result["core_energy_error"], result["core_number_error"]), .1)
        self.assertLess(result["max_secondary"], .01)
