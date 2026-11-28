import json
import unittest

import numpy as np

from experiments.closure_ledger.gravity import inner_horizon_flux as flux
from experiments.closure_ledger.gravity import inner_horizon_coupled as coupled


class InnerHorizonTests(unittest.TestCase):
    def test_background_and_nonzero_bifurcation_lapse(self):
        ref=flux.reference(.01)
        a,b,c=ref['squared_horizons']
        h=np.array([-4.,-.5,0.,.25,.6,1.])
        r,ratio=flux.background(ref,np.arcsinh(h))
        s=r*r
        np.testing.assert_allclose((s-a)*(s-b)*(c-s)/(2*s*s),h*np.exp(ratio),
                                   atol=6e-16,rtol=2e-12)
        self.assertAlmostEqual(s[2],b)
        self.assertTrue(np.all(np.isfinite(ratio)))

    def test_weak_transport_path_balance_and_fixed_clock(self):
        runs=json.loads(flux.OUT.read_text())['runs']
        fine=runs[3]
        self.assertEqual(len(fine['samples']),len(flux.FRACTIONS))
        for s in fine['samples']:
            self.assertLess(abs(s['tau']-s['integrated_tau']),4e-7)
            for v in s['values']:
                for k in ('from_right_going','from_left_going','geometric_variation'):
                    self.assertLess(abs(v['mass_gain']-v[k]),5e-12)
        self.assertEqual(fine['samples'][0]['values'][0]['mass_gain'],0.)
        for index in (0,1):
            gains=[run['samples'][-1]['values'][index]['mass_gain'] for run in runs]
            self.assertLess(abs(gains[2]-gains[3]),.3*abs(gains[1]-gains[2]))
            self.assertLess(abs(gains[4]/gains[2]-1),1e-5)
        last=fine['samples'][-1]
        self.assertGreater(last['values'][0]['delta_u'],abs(last['vacuum_u']))

    def test_amplitude_and_live_characteristic_transport(self):
        first=flux.transport(.024,epsilon=.002,table_step=.001)
        second=flux.transport(.024,epsilon=.004,table_step=.001)
        for a,b in zip(first['samples'],second['samples']):
            for va,vb in zip(a['values'],b['values']):
                self.assertAlmostEqual(vb['mass_gain'],2*va['mass_gain'],places=16)

    def test_coupled_subgap_sign_and_comparison_to_parent(self):
        rows=json.loads(coupled.OUT.read_text())
        parent=json.loads((coupled.OUT.parents[1]/'20261126_subgap_neck/probe.json').read_text())['runs']
        for e in (.005,.05):
            run=next(r for r in rows if r['epsilon']==e and r['step']==.003)
            self.assertEqual(run['turns'],[])
            self.assertEqual(run['stop'],'valid_stencil_limit')
            other=next(r for r in parent if r['mode']=='standing' and r['epsilon']==e
                       and r['n']==16384 and r['optical'])
            row=run['samples'][1]
            archived=min(other['samples'],key=lambda s:abs(row['tau']-s['tau']))
            self.assertLess(abs(row['u']-archived['u']),2e-7)
            self.assertLess(abs(row['mass_gain']/archived['stress_gain']-1),1e-4)
            self.assertLess(abs(row['u_integral']-row['u']),1e-7)
            self.assertTrue(row['theta_plus']<0 and row['theta_minus']<0)
            finer=next(r for r in rows if r['epsilon']==e and r['step']==.0015)
            event=finer['samples'][2]
            self.assertLess(abs(event['geometric_gain']/event['stress_gain']-1),.008)
            self.assertLess(abs(event['u']-event['u_integral']),1e-6)
            self.assertGreater(finer['end']['tau'],1.)
            self.assertGreater(abs(finer['end']['geometric_gain']-finer['end']['stress_gain']),
                               abs(finer['end']['stress_gain']))

    def test_live_vacuum_and_matter_control(self):
        for e in (0.,.005):
            run=coupled.coupled(.02,e,end_fraction=.5)
            point=run['samples'][-1]
            self.assertAlmostEqual(point['tau'],.5*flux.reference(.01)['turn_tau'])
            self.assertLess(abs(point['geometric_gain']-point['stress_gain']),2e-10)
            self.assertLess(point['u'],0)
