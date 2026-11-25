import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad, solve_ivp

from experiments.closure_ledger.antipodal_network import expanding_detector_loop as model


class ExpandingLoopTests(unittest.TestCase):
    def test_arrival_clock_and_packet_energy(self):
        # Compute the null distance independently in detector proper time.
        D, tau = model.D, model.base.DURATION
        def a(t):
            return (1+3*max(t-tau, 0.)/D)**(1/3)
        for s in (0., .7, 1.6, 2.4):
            arrival = model.proper_time(D+s, 1.)
            length = quad(lambda t: 1/a(t), s, arrival, points=[tau])[0]
            self.assertAlmostEqual(length, D, places=10)
        def arrived_voltage(t):
            emission = tau+D*(a(t)**2-1)/2-D
            return model.base.incoming(emission)/a(t)
        lo, hi = model.proper_time(D,1.), model.proper_time(D+tau,1.)
        measured = quad(lambda t: arrived_voltage(t)**2, lo, hi)[0]
        expected = quad(lambda s: model.base.incoming(s)**2/a(model.proper_time(s+D,1.)), 0, tau)[0]
        self.assertAlmostEqual(measured, expected, places=10)
        self.assertLess(measured, model.base.PULSE_ENERGY/model.base.GAMMA)

    def test_entire_first_echo_against_proper_time_ode(self):
        D, tau = model.D, model.base.DURATION
        pair = [.8,.2]
        q,p = model.base.preparation(2,363,.02)
        y0 = np.r_[q[0],model.base.inverse_mass(p,.25)[0]]
        def solve(incident, begin, end, state):
            def rhs(t,y):
                force = model.base.GAMMA*(2*incident(t)-y[2:])-model.base.slope(y[:2])
                return np.r_[y[2:],model.base.inverse_mass(force,.25)]
            out = solve_ivp(rhs,(begin,end),state,method="DOP853",dense_output=True,
                            rtol=1e-10,atol=1e-12,max_step=.08)
            self.assertTrue(out.success)
            return out
        def incident(t):
            return np.sqrt(pair)*model.base.incoming(t)
        def a(t):
            return (1+3*max(t-tau,0.)/D)**(1/3)
        first, second = model.proper_time(D,1.), model.proper_time(2*D,1.)
        opened = solve(incident,0.,first,y0)
        def echo(t):
            emission_eta = tau+D*(a(t)**2-1)/2-D
            emitted = model.proper_time(emission_eta,1.)
            return (a(emitted)/a(t))*(opened.sol(emitted)[2:]-incident(emitted))[::-1]
        reference = solve(echo,first,second,opened.y[:,-1])
        errors = []
        for limit in (.05,.025):
            row, trace = model.evolve(pair,1.,limit=limit,count=2)
            k = row["slots"]
            times = model.proper_time(np.arange(k,2*k+1)*D/k,1.)
            errors.append(np.max(abs(trace[k:] - reference.sol(times)[:2].T)))
        self.assertLess(errors[1],.002)
        self.assertLess(errors[1],errors[0]/3)

    def test_archive_counts_and_gates(self):
        root = Path(__file__).resolve().parents[1]
        data = json.loads((root/"experiments/closure_ledger/runs/20261125_expanding_loop/probe.json").read_text())
        for row in data["rows"]:
            counts = row["persistent"]
            self.assertEqual(sum(counts.values()),512)
            self.assertEqual(sum(row["final"].values()),512)
            self.assertTrue(0 <= row["unbroken_and_durable"] <= row["unbroken_singles"] <= row["first_singles"] <= 512)
            self.assertLessEqual(row["unbroken_and_durable"],counts["left_only"]+counts["right_only"])
            self.assertLess(row["max_budget_residual"],1e-8)
            self.assertGreaterEqual(row["lost_fraction_mean"],-1e-10)
        for label,limit in (("coarse_gates",.025),("fine_gates",.0125)):
            primary = sorted((r for r in data["rows"] if r["rate"]==1 and r["mode"]=="redshift"
                              and r["limit"]==limit and r["input"]!="bright"),key=lambda r:int(r["input"]))
            self.assertEqual([r["input"] for r in primary],list(map(str,range(9))))
            c = [r["persistent"] for r in primary]
            right = np.array([v["right_only"]/512 for v in c])
            single = right + np.array([v["left_only"]/512 for v in c])
            opened = next(r["persistent"] for r in data["rows"] if r["mode"]=="open" and r["input"]=="4")
            error = max(abs(right-np.cos(np.arange(9)*np.pi/16)**2))
            self.assertEqual(data[label],dict(max_born_error=error,born_window=error<=.01,
                              near_unit_every_angle=min(single)>=.99,
                              balanced_above_half_and_open=single[4]>max(.5,(opened["right_only"]+opened["left_only"])/512)))
