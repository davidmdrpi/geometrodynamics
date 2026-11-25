import unittest
from math import tan

from experiments.closure_ledger import alice_ground_split as splitting


class GroundSplittingTests(unittest.TestCase):
    def test_m_one_prefactor(self):
        self.assertAlmostEqual(splitting.prediction(.24, 1) / tan(.12)**4, 6.)

    def test_sub_roundoff_splitting(self):
        measured = splitting.measure(.08, 4)["split"]
        self.assertGreater(measured, 0.)
        self.assertLess(abs(measured / splitting.independent(.08, 4) - 1), 1e-7)
