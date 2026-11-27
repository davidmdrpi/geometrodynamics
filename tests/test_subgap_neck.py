import json
import unittest

import numpy as np
from scipy.integrate import quad

from experiments.closure_ledger.gravity import subgap_neck as neck


class SubgapNeckTests(unittest.TestCase):
    def test_fixed_neck_and_constraint_solved_load(self):
        ref = neck.reference(.01)
        for eps in (.005,.05,1.):
            pos,plus,minus,diag = neck.prepare(4096,.01,eps,'standing')
            self.assertAlmostEqual(pos[0,0]**2,.4)
            self.assertAlmostEqual(plus[0,0]+minus[0,0],0.)
            self.assertAlmostEqual(diag['recovered_load']/ref['mass_gap'],eps,places=8)
            self.assertLess(diag['mass_error'],1e-14)
            fine = neck.prepare(8192,.01,eps,'standing')[3]
            self.assertLess(fine['ham_max'],.12*diag['ham_max'])
            self.assertLess(fine['ham_max'],1e-7)
            self.assertEqual(diag['mom_max'],0.)
            self.assertEqual(diag['scalar_shift_charge'],0.)

    def test_chiral_controls_and_quotient(self):
        diag = [neck.prepare(4096,.01,.05,m)[3] for m in ('standing','one_way','paired')]
        energies = np.array([d['energy_each_direction'] for d in diag])
        np.testing.assert_allclose(energies.sum(axis=1),energies[0].sum(),rtol=1e-12)
        self.assertEqual(energies[1,1],0.)
        self.assertGreater(diag[1]['null_generation_rate'],.05)
        self.assertGreater(diag[1]['signed_half_shift_h_error'],1e-3)
        for i in (0,2):
            self.assertLess(diag[i]['signed_half_shift_h_error'],1e-12)
        for d in diag:
            self.assertLess(abs(d['scalar_shift_charge']),1e-13)
            self.assertLess(d['mom_max'],1e-16)

    def test_vacuum_turn_and_optical_feet(self):
        ref=neck.reference(.01)
        times=ref['turn_tau']*np.array([.9,.99,1.,1.005])
        orbit=neck.exact_clock(ref,times)
        self.assertAlmostEqual(orbit[2][0]**2,.396,places=12)
        self.assertTrue(orbit[1][1]<0<orbit[3][1])
        self.assertAlmostEqual(neck.vacuum_foot(ref,orbit[0][0]),-58.9471,places=3)
        self.assertLess(neck.vacuum_foot(ref,orbit[1][0]),-130)

    def test_causal_collar_with_independent_finite_radius_quadrature(self):
        ref=neck.reference(.01)
        for eps in (0.,.005,.05,1.):
            d=neck.cosmological_collar(ref,eps)
            self.assertGreater(d['collar_margin'],.21)
            mu,q=d['outer_mass'],ref['charge']
            def inverse_metric(r):
                return 1/(1-mu/r**2+q*q/r**4-r*r/2)
            root=d['limiting_initial_foot']
            for R,foot in zip((2.,10.,100.),d['finite_radius_feet']):
                inside=quad(inverse_metric,root,foot)[0]
                exterior=-quad(inverse_metric,R,np.inf)[0]
                self.assertAlmostEqual(inside,exterior,places=8)

    def test_focusing_step_preserves_ingoing_sign(self):
        n=4096
        pos,plus,minus,_=neck.prepare(n,.01,.05,'one_way',True)
        active,z,p,m=neck.focused_step(pos,plus,minus,np.ones(n,dtype=bool),16/n,
                                        neck.reference(.01)['charge'])
        j=int(np.flatnonzero(active==n//2)[0])
        self.assertTrue(p[0,j]<0 and m[0,j]<0)
        self.assertTrue(np.isfinite(z[:,j]).all())

    def test_archive_and_live_mixed_gauge_replay(self):
        rows=json.loads(neck.OUT.read_text())['runs']
        self.assertEqual(len(rows),36)
        live=neck.run((2048,.05,'standing',.01))
        saved=next(r for r in rows if not r['optical'] and r['n']==2048
                   and r['epsilon']==.05 and r['mode']=='standing')
        self.assertAlmostEqual(live['frames'][-1]['tau'],saved['frames'][-1]['tau'])
        opt=[r for r in rows if r['optical'] and r['n']==16384]
        self.assertTrue(all(not r['turns'] for r in opt))
        for row in opt:
            half=next(s for s in row['samples']
                      if abs(s['tau']/row['reference']['turn_tau']-.5)<1e-10)
            self.assertLess(half['u'],0.)
            self.assertLess(abs(half['mass']-row['reference']['mu']-half['stress_gain']),4e-10)
