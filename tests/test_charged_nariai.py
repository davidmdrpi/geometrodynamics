import json
from pathlib import Path
import unittest

import numpy as np
from scipy.optimize import brentq

from experiments.closure_ledger.gravity import charged_nariai as charged
from experiments.closure_ledger.gravity import charged_characteristics as rays
from experiments.closure_ledger.gravity import nariai_loop as loop


ARCHIVE = (Path(__file__).resolve().parents[1]/"experiments/closure_ledger/runs"
           /"20261125_charged_nariai"/"probe.json")


class ChargedNariaiTests(unittest.TestCase):
    def test_product_solves_both_gauges_and_constraint(self):
        q = .14
        s = brentq(lambda z: z*z*(1-z)-q*q,2/3,1,xtol=1e-15)
        r = np.sqrt(s)
        k = 2-4*q*q/r**6
        self.assertGreater(k,0.)
        self.assertAlmostEqual(3/r**2,3+3*q*q/r**6)
        n = 64
        dx = 2*np.pi/n
        time = .31
        # dS2 proper-time scale a=cosh(sqrt(k)*t).
        a,h = np.cosh(np.sqrt(k)*time),np.sqrt(k)*np.tanh(np.sqrt(k)*time)
        state = np.tile([r,a,0,h,0,0,0],(n,1)).T
        out = loop.gaussian_rhs(state,dx,q)
        np.testing.assert_allclose(out[2],0.,atol=1e-12)
        np.testing.assert_allclose(out[3],k-h*h,atol=1e-12)
        # Conformal derivative sigma_t = a*h, sigma_tt = k*a².
        state = np.tile([r,np.log(a),0,0,a*h,0],(n,1)).T
        np.testing.assert_allclose(loop.rhs(state,dx,q)[3:],
                                   np.tile([0,k*a*a,0],(n,1)).T,atol=2e-12)

    def test_packet_gauss_parity_and_constraints(self):
        q=.14
        with self.assertRaises(ValueError):
            loop.initial(128,.05,charge=.25)
        y=loop.initial(256,.05,charge=q)
        np.testing.assert_array_equal(y,loop.initial(256,.05,charge=-q))
        np.testing.assert_allclose(np.roll(y,128,axis=1),
                                   np.array([1,1,-1,-1,-1,1])[:,None]*y,atol=1e-12)
        neutral=loop.measure(y,2*np.pi/256,0.)["null_constraint"]
        charged_err=loop.measure(y,2*np.pi/256,0.,q)["null_constraint"]
        self.assertLess(charged_err,2e-6)
        self.assertGreater(neutral,1e-2)
        with self.subTest("signed flux leaves geometry unchanged"):
            np.testing.assert_array_equal(loop.rhs(y,2*np.pi/256,q),loop.rhs(y,2*np.pi/256,-q))

    def test_timelike_turn_at_and_below_inner_horizon(self):
        r=charged.exact()
        saved=json.loads(ARCHIVE.read_text())["electrovac"]
        self.assertLess(r["max_mass_defect"],1e-12)
        self.assertLess(r["J_velocity_error"],1e-12)
        self.assertTrue(r["paired_velocities"][0]<0<r["paired_velocities"][1])
        for t,h in zip(r["trajectories"],saved["trajectories"]):
            self.assertAlmostEqual(t["radius"],h["radius"],places=10)
            self.assertAlmostEqual(charged.metric(t["radius"],r["mu"],r["charge"]),
                                   t["killing_energy"]**2,places=9)
            self.assertGreater(t["outward_velocity_after"],0.)
            self.assertGreater(t["force"],0.)
            self.assertLess(t["energy_residual"],1e-10)
        zero,moving=r["trajectories"]
        self.assertAlmostEqual(zero["radius"],r["horizons"][0],places=11)
        self.assertLess(moving["radius"],r["horizons"][0])
        np.testing.assert_allclose(r["neck_clocks"],zero["time_from_outer"]+np.array([-.25,.25]),
                                   rtol=1e-10)
        # Independent leading RN curvature divergence.
        mu,q=r["mu"],r["charge"]
        self.assertAlmostEqual(charged.curvature(1e-4,mu,q)*1e-48/q**4,508.,places=3)

    def test_packet_archive_and_independent_crossings(self):
        a=json.loads(ARCHIVE.read_text())
        self.assertEqual([p['charge'] for p in a['packets']],[.14,.14,.14,0.,.08,.17])
        primary=charged.packet(256,.14)
        row=a["packets"][1]
        self.assertEqual(primary['gaussian_stop'],row['gaussian_stop'])
        self.assertAlmostEqual(primary['gaussian_end']['t'],row['gaussian_end']['t'],places=7)
        self.assertLess(max(a['packets'][2]['cross_gauge_error']),1e-9)
        for r in a['homogeneous']:
            actual=charged.homogeneous(.14,r['scale_floor'])
            self.assertAlmostEqual(actual['r'],r['r'],places=8)
            self.assertLess(actual['u'],0.)
            self.assertGreater(actual['peak_y_acc'],0.)
        raw=rays.propagate(2048)
        saved=a['characteristics'][0]
        for new,ref in zip(raw['crossings'],saved['crossings'],strict=True):
            self.assertAlmostEqual(new['mass'],ref['mass'],places=8)
        by_n={r['n']:r for r in a['characteristics'] if r['remapped']}
        for n in (8192,16384,32768):
            run=by_n[n]
            self.assertEqual(run['center_minima'],[])
            self.assertEqual(run['stop'],'valid_stencil_limit')
            pts={c['target']:c for c in run['crossings']}
            for radius in (.23,.205,.2,.198,.197):
                c=pts[radius]
                self.assertAlmostEqual(c['r'],radius)
                self.assertLess(max(c['theta_plus'],c['theta_minus']),0.)
                self.assertAlmostEqual(c['u'],radius*(c['theta_plus']+c['theta_minus'])/6,delta=.003)
                self.assertEqual(c['y_acc']>0,radius>=.198)
        for radius in (.23,.205,.2,.198):
            cells=[next(c for c in by_n[n]['crossings'] if c['target']==radius)
                   for n in (8192,16384)]
            self.assertLess(abs(cells[0]['u']-cells[1]['u']),.002)
            self.assertLess(cells[1]['null_constraint'],cells[0]['null_constraint'])

    def test_null_relabelling_and_maxwell_null_stress(self):
        x=np.linspace(-np.pi,np.pi,101)
        b,jac,deriv=rays.null_map(x,2.292,.015,.18)
        h=1e-6
        bp,jp,_=rays.null_map(x+h,2.292,.015,.18)
        bm,jm,_=rays.null_map(x-h,2.292,.015,.18)
        np.testing.assert_allclose((bp-bm)/(2*h),jac,rtol=1e-7)
        np.testing.assert_allclose((jp-jm)/(2*h),deriv,atol=1e-6,rtol=1e-7)
        np.testing.assert_allclose(b,-b[::-1],atol=1e-12)
        shifted=rays.null_map(x+np.pi,2.292,.015,.18)[0]
        np.testing.assert_allclose(shifted-b,np.pi*1.36)
        # At reflection, conformal constraints contain no direct Maxwell null source:
        # the shift in the null residual is exactly the orbit acceleration shift.
        y=loop.initial(128,.05,charge=.14)
        diff=loop.rhs(y,2*np.pi/128,.14)-loop.rhs(y,2*np.pi/128)
        np.testing.assert_allclose(diff[3],2*.14**2/y[0]**5)
        np.testing.assert_allclose(diff[4],-7*.14**2/y[0]**6)
