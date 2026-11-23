import unittest

import numpy as np

from experiments.closure_ledger.antipodal_network import dilatonic_waist as waist
from experiments.closure_ledger.antipodal_network import elliptic_transport as transport
from experiments.closure_ledger.antipodal_network import natural_priors as priors


class PointerAndWaistTests(unittest.TestCase):
    def test_release_momentum_changes_uniform_phase_records(self):
        rows = priors.run(128)
        self.assertEqual(rows["sphere_projection_at_rest"]["max_error"], 0.)
        self.assertGreater(rows["phase_at_rest"]["max_error"], 0.09)
        self.assertGreater(rows["oscillator_phase"]["max_error"], 0.2)
        self.assertLess(rows["thermal_rms_factor_1"]["positive"][-1], 128)

    def test_mixer_axis_and_duration(self):
        for theta in transport.THETAS:
            p, leakage = transport.split_weights(theta,.02)
            self.assertAlmostEqual(transport.split_weights(theta,.02,enabled=False)[0], 0.5)
            self.assertLess(abs(p-np.cos(theta)**2), 1e-4)
            self.assertLess(leakage, 2e-4)
        self.assertAlmostEqual(transport.split_weights(0,.02,axis=0)[0], .5, places=10)

    def test_trajectory_equivariance_and_nonequilibrium(self):
        ranks = (np.arange(256)+0.5)/256
        eq = transport.transport(.75,ranks,1.)
        cold = transport.transport(.75,ranks,.6)
        self.assertLess(eq["max_rank_error"], 1e-8)
        self.assertLess(abs(eq["right_fraction"]-eq["right_mass"]), 1/256)
        self.assertGreater(cold["right_fraction"]-eq["right_fraction"], .09)

    def test_einstein_characteristic_and_scalar_decoupling(self):
        riem = waist.curvature(0.,0.,0.)
        for omega, wave in ((0.,1.),(1.,0.),(1.,1.)):
            tensor, scalar = waist.shear_symbol(riem,np.zeros((5,5)),0.,0.,omega,wave)
            self.assertAlmostEqual(tensor[3,4], (wave**2-omega**2)/2)
            self.assertAlmostEqual(scalar, 0.)

    def test_transverse_angular_shears(self):
        jet = waist.audit()["jet"]
        rr = waist.curvature(jet["X"], jet["Y"], jet["Z"])
        hh = np.diag([-13/48, 5/48, 0, 0, 0])
        for freq, radial, angular in [(1., 0., 0.), (0., 0., 1.), (.2, -.4, .7)]:
            out, scalar = waist.shear_symbol(rr,hh,0.,1.,freq,radial,angular)
            self.assertAlmostEqual(2*out[3,4], -7/12*freq**2-radial**2/12-angular**2/2)
            self.assertAlmostEqual(scalar,0.)
        for alpha in (-.1, -.25, -.5):
            out, _ = waist.shear_symbol(rr,np.zeros((5,5)),alpha,0.,.7,0.,.3,(1,3))
            expected = np.zeros((5,5))
            expected[1,3] = expected[3,1] = (-(1+4*alpha)*.7**2+.3**2)/2
            np.testing.assert_allclose(out,expected,atol=1e-12)

    def test_on_shell_jet_and_independent_tensor_contractions(self):
        result = waist.audit()
        self.assertLess(result["random_contraction_error"],1e-10)
        jet = result["jet"]
        self.assertLess(max(jet["metric_residual"], abs(jet["scalar_residual"])),1e-12)
        self.assertAlmostEqual(jet["radial_speed_squared"],-1/7)
