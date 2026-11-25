import json
from pathlib import Path
import unittest

import numpy as np

from experiments.closure_ledger.antipodal_network import symmetric_two_visit as model


class SymmetricVisitTests(unittest.TestCase):
    def test_arc_gradient_tangential_variations(self):
        rng = np.random.default_rng(368)
        nodes = model.unit(rng.normal(size=(9,4)))
        variation = rng.normal(size=nodes.shape)
        variation -= (nodes*variation).sum(axis=1,keepdims=True)*nodes
        eps = 1e-6
        diff = (model.arc_action(model.unit(nodes+eps*variation),2.)
                -model.arc_action(model.unit(nodes-eps*variation),2.))/(2*eps)
        self.assertAlmostEqual(diff,float((model.arc_gradient(nodes)*variation).sum()),places=6)

    def test_certificates_use_actual_two_arc_paths(self):
        for k in model.LOADS:
            d = k/4
            for floating in (False,True):
                # At zero offset both readings are one half. In the floating
                # trial the two directions straddle the second coordinate axis.
                a = np.pi/4-d/2 if floating else np.pi/4
                angles = np.array([a,a+np.pi/2+d,a+np.pi])
                nodes = np.zeros((3,4))
                nodes[:,0],nodes[:,2] = np.cos(angles),np.sin(angles)
                value = model.arc_action(nodes,2.)-k*sum(model.rank(nodes[:2]))
                gain = model.BARE-k-value
                self.assertAlmostEqual(gain,model.direct_certificate(k,floating)["trial_gain"])

    def test_invariant_orbit_discriminator(self):
        rng = np.random.default_rng(236)
        x = model.unit(rng.normal(size=(150000,4)))
        w = np.stack((-x[:,2],x[:,3],x[:,0],-x[:,1]),axis=1)
        phi = rng.uniform(0,2*np.pi,size=(len(x),1))
        w = w*np.cos(phi)+model.hopf(w)*np.sin(phi)
        for chi in (0.,4.,4.1,6.):
            pred = model.pointer_prediction(chi)
            c = pred["overlap"]
            y = c*x+np.sqrt(1-c*c)*w
            u,v = model.rank(x),model.rank(y)
            self.assertAlmostEqual(abs(u+v-1).mean(),pred["invariant_mean_defect"],delta=.004)
            self.assertAlmostEqual(((u<.5)&(v<.5)).mean(),pred["balanced_doubles"],delta=.004)
            self.assertAlmostEqual(((u>.5)&(v>.5)).mean(),pred["balanced_misses"],delta=.004)

    def test_archive_and_free_node_replay(self):
        root = Path(__file__).resolve().parents[1]
        path = root/"experiments/closure_ledger/runs/20261125_symmetric_two_visit/probe.json"
        for report in (json.loads(path.read_text()),model.probe()):
            direct = report["direct"]
            self.assertEqual([(r["kappa"],r["floating"]) for r in direct],
                             [(k,f) for k in model.LOADS for f in (False,True)])
            for r in direct:
                k,f = r["kappa"],r["floating"]
                cert = model.direct_certificate(k,f)
                self.assertEqual(len(r["fits"]),3)
                self.assertAlmostEqual(cert["trial_gain"],r["certificate"]["trial_gain"])
                for fit in r["fits"]:
                    defect = fit["u"]+fit["v"]-1
                    if k:
                        self.assertGreaterEqual(k*defect,cert["trial_gain"]-1e-9)
                    else:
                        self.assertLess(abs(defect),1e-5)
            self.assertEqual([r["chi"] for r in report["pointer"]],list(model.SUSCEPTIBILITIES))
            for r in report["pointer"]:
                chi = r["chi"]
                pred = model.pointer_prediction(chi)
                for k in pred:
                    self.assertAlmostEqual(r["prediction"][k],pred[k],places=10)
                self.assertEqual(len(r["fits"]),3)
                for fit in r["fits"]:
                    self.assertAlmostEqual(fit["action"],pred["action"],places=8)
                    self.assertAlmostEqual(fit["pointer"],np.sqrt(chi)*fit["real_overlap"],places=5)
                    self.assertLess(abs(abs(fit["real_overlap"])-pred["overlap"]),
                                    .003 if chi==4 else 2e-5)
            pair = report["ready_pair_followup"]
            gain = model.BARE+1/8-(model.BARE+2/16**2+(1-np.sin(1/16))**2/8)
            self.assertAlmostEqual(pair["action_improvement"],gain)
            self.assertAlmostEqual(pair["minimum_defect_bound"],1-np.sqrt(1-8*gain))
            for fit in pair["fits"]:
                self.assertLessEqual(fit["action"],model.BARE+1/8-gain+1e-9)
                self.assertLessEqual(fit["u"]+fit["v"],1-pair["minimum_defect_bound"]+1e-9)
                np.testing.assert_allclose(fit["paired_readout"],np.array([fit["u"],fit["v"]])/2,
                                           atol=1e-6)
