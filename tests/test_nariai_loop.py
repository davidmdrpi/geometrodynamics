import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from experiments.closure_ledger.gravity import nariai_loop as loop


class NariaiLoopTests(unittest.TestCase):
    def test_exact_product_and_vacuum_collapse(self):
        for t in (-.7,0.,.9):
            r,a,v,w = loop.vacuum(t,np.linspace(0,2*np.pi,20),0.)
            np.testing.assert_allclose(r,1)
            np.testing.assert_allclose(a,np.sqrt(2)*np.cosh(np.sqrt(2)*t))
            np.testing.assert_allclose(v,0)
            np.testing.assert_allclose(w,0)
        for eps in (.05,.1,.2):
            mu=(1-eps**2)/2
            run=loop.vacuum_probe(eps)
            for c in run['centers']:
                def rhs(t,y):
                    return [y[1], y[0]/2-mu/y[0]**3]
                def event(t,y):
                    return y[0]-.1
                event.terminal=True
                sol=solve_ivp(rhs,(0,5),[c['initial_radius'],c['normal_velocity']],events=event,
                              rtol=2e-11,atol=1e-13)
                self.assertTrue(sol.success)
                self.assertAlmostEqual(sol.t_events[0][0]+run['near_pinch'][-1]['proper_time_left'],
                                       c['proper_time_to_future_pinch'],places=8)
            fine=loop.vacuum_probe(eps,size=512)
            for key in ('max_hamiltonian','max_momentum'):
                self.assertLess(fine[key],run[key]/12)
            self.assertLess(run['max_mass_defect'],1e-12)
            self.assertLess(run['max_longitudinal_pair_error'],1e-12)

    def test_full_J_data_and_evolution_constraints(self):
        n=256
        dx=2*np.pi/n
        y=loop.initial(n,.05)
        parity=np.array([1,1,-1,-1,-1,1])[:,None]
        np.testing.assert_allclose(np.roll(y,n//2,axis=1),parity*y,atol=1e-12)
        plus=loop.rk4(y,dx,.002,loop.rhs)
        minus=loop.rk4(y,dx,-.002,loop.rhs)
        np.testing.assert_allclose(np.roll(minus,n//2,axis=1),parity*plus,atol=1e-12)
        coarse=loop.measure(y,dx,0.)['null_constraint']
        fine=loop.measure(loop.initial(2*n,.05),dx/2,0.)['null_constraint']
        self.assertLess(fine,coarse/12)
        # Keeping the pulse but deleting its supporting momentum flux fails.
        y[2]=0
        self.assertGreater(loop.measure(y,dx,0.)['null_constraint'],.05)

    def test_homogeneous_scalar_exact_flux_and_independent_ode(self):
        # Conserved homogeneous scalar flux eliminates the wave variables.
        p0=np.sqrt(6*(1/.98**2-1))
        def f(t,y):
            r,a,u,h=y
            p=p0*.98**3/(a*r**3)
            acc=r-(1+u*u)/r-r*p*p/6
            return [u,a*h,acc,2-p*p-3*acc/r-h*h]
        def event(t,y):
            return y[0]-.2
        event.terminal=True
        s=solve_ivp(f,(0,4),[.98,1,0,0],events=event,rtol=1e-11,atol=1e-13)
        calc=loop.gaussian(256,0.,cut=.2)
        self.assertTrue(s.success)
        self.assertEqual(calc['stop'],'radius_cutoff')
        self.assertAlmostEqual(calc['frames'][-1]['t'],s.t_events[0][0],places=5)

    def test_archive_and_cutoff_replay(self):
        path=Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20261125_nariai_loop/probe.json'
        saved=json.loads(path.read_text())
        fresh=loop.probe()
        for data in (saved,fresh):
            conf=data['conformal']
            self.assertEqual([r['n'] for r in conf[:3]],[128,256,512])
            self.assertTrue(all(r['stop']=='lapse_cutoff' for r in conf[:3]))
            self.assertTrue(all(f['ray']['affine_velocity']<0 for f in conf[2]['frames']))
            self.assertTrue(np.all(np.diff([f['ray']['affine_velocity'] for f in conf[2]['frames']])<=1e-8))
            self.assertLess(conf[2]['antipode_maximum']['r'],1.)
            gs=data['gaussian']
            self.assertEqual([r['cut'] for r in gs[:5]],[.1,.1,.1,.3,.2])
            for r in gs[:5]:
                last=r['frames'][-1]
                self.assertEqual(r['stop'],'radius_cutoff')
                self.assertAlmostEqual(last['minimum']['radius'],r['cut'],places=8)
                self.assertLess(max(last['minimum'][k] for k in ('theta_plus','theta_minus')),0)
                self.assertLess(last['light_positions'][0],last['light_positions'][1])
            last=[r['frames'][-1] for r in gs[:3]]
            self.assertLess(abs(last[2]['t']-last[1]['t']),1e-5)
            self.assertLess(last[2]['max_momentum'],last[1]['max_momentum']/10)
            self.assertLess(last[2]['max_hamiltonian'],1e-6)
            self.assertLess(max(data['cross_gauge'][-1]['radius_discrepancies']),2e-8)
            self.assertEqual(gs[-1]['stop'],'caustic_guard')
            for run in (gs[2],*data['short_loop_followup']):
                bound = run['center_comparison_time']
                velocity = -run['amplitude']
                def upper(t):
                    return (1+(.98**2-1)*np.cosh(np.sqrt(2)*t)
                            +np.sqrt(2)*.98*velocity*np.sinh(np.sqrt(2)*t))
                self.assertAlmostEqual(upper(bound),0.)
                for f in run['frames']:
                    self.assertLessEqual(f['centers'][0]['radius']**2,upper(f['t'])+1e-9)
            for run in data['short_loop_followup']:
                self.assertEqual(run['stop'],'radius_cutoff')
                self.assertLess(max(run['center_signal_arrivals']),.61)
                self.assertGreater(run['frames'][-1]['t'],2.4)
                self.assertAlmostEqual(run['frames'][-1]['centers'][0]['radius'],.1,places=8)
        for section in ('gaussian','conformal','short_loop_followup'):
            self.assertEqual(len(saved[section]),len(fresh[section]))
            for a,b in zip(saved[section],fresh[section]):
                self.assertEqual(a['stop'],b['stop'])
                self.assertAlmostEqual(a['frames'][-1]['t'],b['frames'][-1]['t'],places=8)
