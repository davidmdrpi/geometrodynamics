import unittest

import numpy as np

from experiments.closure_ledger.antipodal_network import elliptic_pointer as pointer
from experiments.closure_ledger.antipodal_network import elliptic_pointer_ring as ring


class EllipticPointerTests(unittest.TestCase):
    def test_orientation_and_relative_phase_set_channel_populations(self):
        winding = np.array([1, 1j])/np.sqrt(2)
        for angle, expected in [(0, 0), (np.pi/2, -1), (-np.pi/2, 1)]:
            plus, minus = pointer.channels(winding, angle)
            self.assertAlmostEqual(abs(plus)**2-abs(minus)**2, expected)

    def test_individual_balanced_records_and_separatrix(self):
        initial = np.array([-0.06, 0., 0.06])
        amplitude = np.sqrt(0.5)
        result = pointer.evolve(initial, amplitude, amplitude)
        np.testing.assert_array_equal(result["records"], [1, 0, -1])
        np.testing.assert_allclose(result["q_switch"], -pointer.DECAY*initial, atol=1e-10)
        np.testing.assert_allclose(result["v_switch"], 0, atol=1e-10)
        self.assertLess(result["ledger_error"].max(), 1e-8)
        # At z=0, integrating the spring equation gives the exact field phase.
        integral = pointer.DAMPING*(1+pointer.DECAY)*initial
        phase = np.exp(1j*pointer.GAIN*integral)
        np.testing.assert_allclose(result["plus_end"], amplitude*phase, atol=1e-10, rtol=0)
        np.testing.assert_allclose(result["minus_end"], amplitude*phase.conj(), atol=1e-10, rtol=0)

    def test_changing_only_the_apparatus_changes_the_law(self):
        for density in ("uniform", "circle_projection"):
            r = pointer.ensemble(32, pointer.MATCHED_WIDTH, density)
            self.assertEqual(r["basin_disagreements"], 0)
            self.assertEqual(sum(r["no_record"]), 0)
            if density == "uniform":
                self.assertEqual(r["max_born_error"], 0)
            else:
                self.assertGreater(r["max_born_error"], 0.08)

    def test_ring_force_and_projected_coefficient(self):
        basis = np.stack([ring.WINDING, ring.WINDING.conj()], axis=1)
        step = 1e-5
        high = ring.ring_terms(step, basis)[0]
        low = ring.ring_terms(-step, basis)[0]
        projected = basis.conj().T @ (high-low)/(2*step)
        expected = -pointer.GAIN*np.array([[0, 1], [1, 0]])
        np.testing.assert_allclose(projected, expected, atol=1e-10)
        field = ((basis[:, 0]+0.5j*basis[:, 1])/np.sqrt(1.25))[:, None]
        q = 0.11
        force = ring.ring_terms(q, field)[2]
        derivative = (ring.ring_terms(q+step, field)[1]-ring.ring_terms(q-step, field)[1])/(2*step)
        np.testing.assert_allclose(force, -derivative, atol=1e-10)

    def test_ring_with_leakage_still_records_both_basins(self):
        r = ring.run(16)
        self.assertEqual(r["positive"][4], 8)
        self.assertEqual(sum(r["no_record"]), 0)
        self.assertGreater(r["max_leakage_at_switch"], 1e-7)
        self.assertLess(r["max_energy_residual"], 1e-8)
