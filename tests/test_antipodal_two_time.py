import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from experiments.closure_ledger.antipodal_network import antipodal_two_time as bvp


class AntipodalBoundaryTests(unittest.TestCase):
    def test_stationary_vertical_control_is_not_minimal(self):
        x = np.array([1., 0., 0., 0.])
        t = np.linspace(0,np.pi,9)[:,None]
        nodes = np.cos(t)*x + np.sin(t)*bvp.hopf(x)
        rng = np.random.default_rng(365)
        perturb = np.zeros_like(nodes)
        perturb[1:-1] = rng.normal(size=(7,4))
        def varied(eps):
            return bvp.arc_action(bvp.unit(nodes+eps*perturb), 2.)
        self.assertAlmostEqual(varied(0),np.pi**2)
        self.assertAlmostEqual((varied(1e-6)-varied(-1e-6))/2e-6,0.,places=7)
        horizontal = np.cos(t)*x + np.sin(t)*np.array([0.,0.,1.,0.])
        self.assertAlmostEqual(bvp.arc_action(horizontal,2.),np.pi**2/2)
        self.assertEqual(bvp.rank(nodes[4]),1.)

    def test_point_load_trial_breaks_complementarity(self):
        start = np.array([1.,0.,1.,0.])/np.sqrt(2)
        tangent = np.array([-1.,0.,1.,0.])/np.sqrt(2)
        for row in bvp.midpoint_loading():
            delta,k = row["trial_shift"], row["load"]
            phase = np.array([0,np.pi/2+delta,np.pi])[:,None]
            nodes = np.cos(phase)*start + np.sin(phase)*tangent
            actual = bvp.arc_action(nodes,2.)-k*bvp.rank(nodes[1])
            improvement = np.pi**2/2-k*.5-actual
            self.assertAlmostEqual(improvement,row["action_improvement"])
            self.assertGreater(improvement,0)
            self.assertAlmostEqual(improvement/abs(k),row["guaranteed_absolute_defect"])

    def test_crossing_times_against_newton_equation(self):
        for energy,barrier in ((.7,.3),(.501,.5),(.5,.7)):
            def rhs(t,state):
                q,v = state
                return [v, 4*barrier*q*(1-q*q) if abs(q)<1 else 0.]
            def crossed(t,state):
                return state[0]-1
            crossed.terminal = True
            result = solve_ivp(rhs,(0,20),[-1,np.sqrt(2*energy)],events=crossed,
                               rtol=1e-11,atol=1e-13,max_step=.02)
            self.assertTrue(result.success)
            hits = result.t_events[0]
            if energy>barrier:
                self.assertAlmostEqual(hits[0],bvp.flight_time(energy,barrier),places=7)
            else:
                self.assertEqual(len(hits),0)
                self.assertLess(result.y[0,-1],-1)
        self.assertTrue(np.isinf(bvp.flight_time(.5,.5)))

    def test_paths_and_unconditional_fractions(self):
        path = (Path(__file__).resolve().parents[1]/"experiments/closure_ledger/runs"
                /"20261125_antipodal_two_time/probe.json")
        saved = json.loads(path.read_text())
        replay = bvp.probe()
        np.testing.assert_allclose(
            [[r[k] for k in ("load","trial_shift","action_improvement","guaranteed_absolute_defect")]
             for r in saved["midpoint_load_followup"]],
            [[r[k] for k in ("load","trial_shift","action_improvement","guaranteed_absolute_defect")]
             for r in replay["midpoint_load_followup"]], atol=1e-12)
        for report in (saved,replay):
            self.assertEqual(len(report["paths"]),9)
            for row in report["paths"]:
                self.assertAlmostEqual(row["action"],min(row["lam"],1)*np.pi**2/2,places=8)
                u,v = row["initial_rank"],row["midpoint_rank"]
                if row["lam"]>1:
                    self.assertAlmostEqual(u+v,1.,places=5)
                if row["lam"]<1:
                    self.assertAlmostEqual(u,v,places=5)
            e = report["ensemble"]
            self.assertLess(e["horizontal_max_defect"],1e-12)
            self.assertAlmostEqual(e["round_mean_defect"],.25,delta=.004)
            for k in ("round_balanced_misses","round_balanced_doubles"):
                self.assertAlmostEqual(e[k],.5-1/np.pi,delta=.004)
            self.assertAlmostEqual(e["vertical_mean_defect"],.5,delta=.004)
            self.assertLess(e["uniform_basin_max_error"],.005)
            self.assertLess(e["biased_basin_max_error_from_p_squared"],.005)
            self.assertEqual([r["index"] for r in report["rows"]],list(range(9)))
            for row in report["rows"]:
                p = np.cos(row["index"]*np.pi/16)**2
                self.assertAlmostEqual(row["eventual"]["right_only"],p)
                self.assertAlmostEqual(row["biased_right"],p*p)
                self.assertEqual([r["time"] for r in row["finite"]],list(bvp.TIMES))
                for finite in row["finite"]:
                    for key,delay in (("equal_elapsed",0),("common_slice",.5)):
                        f = finite[key]
                        self.assertEqual(f["both"],0.)
                        self.assertGreater(f["neither"],0.)
                        self.assertAlmostEqual(sum(f.values()),1.)
                        for prob,energy,t in ((f["right_only"],p,finite["time"]),
                                              (f["left_only"],1-p,finite["time"]-delay)):
                            if prob:
                                self.assertAlmostEqual(bvp.flight_time(energy,prob),t,places=6)
                            else:
                                self.assertGreaterEqual(bvp.flight_time(energy,0.)+1e-12,t)
        for recorded,new in zip(saved["rows"],replay["rows"]):
            for f,g in zip(recorded["finite"],new["finite"]):
                for key in ("equal_elapsed","common_slice"):
                    np.testing.assert_allclose(list(f[key].values()),list(g[key].values()),atol=1e-11)
