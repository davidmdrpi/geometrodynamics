import json
import math
import unittest

import numpy as np

from experiments.closure_ledger.gravity import cold_loop_mass as cold
from experiments.closure_ledger.gravity.charged_nariai import metric, radial_force
from experiments.closure_ledger.gravity.nariai_loop import gaussian_rhs


class ColdMassTests(unittest.TestCase):
    def test_reference_horizons_and_surface_gravities(self):
        kappas = []
        for epsilon in cold.GAPS:
            family = cold.cold_family(epsilon)
            q, mu = family['charge'], family['mu']
            roots = np.sqrt(family['squared_horizons'])
            np.testing.assert_allclose(metric(roots, mu, q), 0., atol=1e-14)
            np.testing.assert_allclose(np.abs(radial_force(roots, mu, q)), family['kappa'])
            self.assertTrue(0 < roots[0] < roots[1] < roots[2])
            kappas.append(family['kappa'][0])
        self.assertTrue(np.all(np.diff(kappas) < 0))

    def test_radial_integral_against_synchronous_einstein_equation(self):
        q = .27
        r, u, p = .61, -.7, 1.1
        state = np.tile([r, .9, u, .2, 0., p, 0.], (16, 1)).T
        du = gaussian_rhs(state, 2*math.pi/16, q)[2,0]
        mass_rate = 2*r*u*(1+u*u-r*r) + 2*r*r*u*du - 2*q*q*u/r**3
        self.assertAlmostEqual(mass_rate, -r**3*p*p*u/3)

    def test_primary_replay_at_fixed_areal_radius(self):
        rows = json.loads(cold.OUTPUT.read_text())['runs']
        primary = [row for row in rows if row['epsilon'] is None]
        for row in primary:
            point = next(p for p in row['crossings'] if p['target'] == .197)
            self.assertLess(point['u'], 0)
            self.assertLess(point['y_acc'], 0)
        point = next(p for p in primary[-1]['crossings'] if p['target'] == .197)
        self.assertAlmostEqual(point['y_acc'], -990.7, delta=1.)
        self.assertAlmostEqual(point['mass'], point['integrated_mass'], delta=.001)

    def test_live_spatial_return_and_refinement(self):
        saved = json.loads(cold.OUTPUT.read_text())['spatial']
        fine = saved[-1]['profiles'][-1]
        whole_loop = 2*math.pi
        self.assertAlmostEqual(fine['null_grid_time'], whole_loop)
        self.assertTrue(all(p['valid'] and p['u'] < 0 for p in fine['samples']))
        self.assertGreater(fine['samples'][0]['mass'], 3e5)
        masses = np.array([row['profiles'][-1]['samples'][0]['mass'] for row in saved])
        self.assertLess(abs(masses[2]-masses[1]), .3*abs(masses[1]-masses[0]))
        q = cold.cold_family(.01)['charge']
        direct = cold.propagate(128, charge=q, radius=math.sqrt(.4), amplitude=.0002,
                                targets=(), profile_times=(0.,), time_limit=.05, audit_mass=True)
        self.assertAlmostEqual(direct['mass_profiles'][0]['samples'][0]['mass'],
                               saved[0]['profiles'][0]['samples'][0]['mass'])
        self.assertAlmostEqual(direct['frames'][-1]['mass'],
                               direct['frames'][-1]['integrated_mass'], places=5)

    def test_cold_packet_replay_and_mass_gap(self):
        rows = json.loads(cold.OUTPUT.read_text())['runs']
        old = next(row for row in rows if row['epsilon'] == .01 and row['n'] == 512)
        replay = cold.simulate((.01, 512))
        self.assertEqual(replay['center_minima'], [])
        self.assertEqual(len(replay['crossings']), len(cold.SAMPLES))
        for actual, archived in zip(replay['crossings'], old['crossings']):
            self.assertAlmostEqual(actual['mass'], archived['mass'], places=9)
            self.assertAlmostEqual(actual['integrated_mass'], archived['integrated_mass'], places=9)
        q = replay['charge']
        cold_root = min(s.real for s in np.roots([1., -1., 0., q*q]) if 0 < s.real < 2/3)
        cold_mass = 2*cold_root-1.5*cold_root**2
        initial_mass = replay['frames'][0]['mass']
        inward = replay['crossings'][1]
        self.assertAlmostEqual(inward['r']**2, replay['reference']['squared_horizons'][0])
        self.assertGreater((inward['mass']-initial_mass)/(initial_mass-cold_mass), 18.)
