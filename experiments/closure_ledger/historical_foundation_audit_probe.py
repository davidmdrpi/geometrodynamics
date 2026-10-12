"""Targeted dependency audit; old producers and archives are not rewritten."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.closure_ledger import pair_production_threshold_probe as pair
from experiments.closure_ledger import antipodal_focusing_threshold_probe as focus
from geometrodynamics.bell.bulk_identity import make_bulk_pair, bulk_amplitude

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261012_historical_foundation_audit/result.json'


def threshold_audit():
    minimum = pair.test_T1_lowest_stable_configuration()
    threshold = pair.test_T3_threshold_2mc2()
    ratio = minimum['E_at_min_J']/pair.M_E_C2_J
    return dict(legacy_minimum_test_pass=minimum['pass'],
                legacy_threshold_test_pass=threshold['pass'],
                minimum_to_electron_rest_energy=ratio,
                analytic_ratio=3*pair.ALPHA/4,
                computed_pair_minimum_MeV=2*minimum['E_at_min_J']*pair.MEV_PER_J,
                asserted_pair_threshold_MeV=threshold['pair_threshold_MeV'],
                conclusion='The two passing legacy checks do not test their claimed energy equality')


def finite_band_focus_audit():
    n = np.arange(1, 17)
    coefficients = np.exp(-n*n/40.)
    omega = focus.conformal_frequencies(len(n))
    grid = np.linspace(.03, np.pi-.03, 301)
    reflected = focus.physical_field(coefficients, np.pi-grid, 0., omega)
    evolved = focus.physical_field(coefficients, grid, np.pi, omega)
    pole = float(np.sum(n*coefficients))
    values = [float(focus.physical_field(coefficients, np.array([h]), 0., omega)[0])
              for h in [.001, .0005, .00025]]
    return dict(mode_count=len(n), exact_north_pole_limit=pole,
                near_pole_values=values,
                near_pole_errors=[abs(v-pole) for v in values],
                antipodal_identity_error=float(np.max(abs(evolved+reflected))),
                conclusion='Finite harmonic sum has a finite pole limit; exact linear refocus survives')


def readout_audit():
    pair_state = make_bulk_pair()
    def probabilities(a,b,power):
        weights = np.array([[abs(bulk_amplitude(pair_state,a,b,sa,sb))**power
                             for sb in [1,-1]] for sa in [1,-1]])
        return weights/weights.sum()
    records=[]
    for power in [2,4]:
        def E(a,b):
            return float(np.sum(probabilities(a,b,power)*np.array([[1,-1],[-1,1]])))
        settings = [0., np.pi/2, np.pi/4, -np.pi/4]
        a0,a1,b0,b1=settings
        S=abs(E(a0,b0)+E(a0,b1)+E(a1,b0)-E(a1,b1))
        error=max(float(np.max(abs(probabilities(a,b,power).sum(axis=axis)-.5)))
                  for a in np.linspace(0,np.pi,7) for b in np.linspace(0,np.pi,9) for axis in [0,1])
        records.append(dict(amplitude_power=power,chsh=S,max_marginal_error=error))
    return dict(records=records,
                scope='Alternative normalized readout on the same amplitudes; not a quantum or GR realization',
                conclusion='Binary outcomes, amplitudes and no-signaling do not by themselves select the squared-amplitude readout')


def measure():
    threshold, focus_result, bell = threshold_audit(), finite_band_focus_audit(), readout_audit()
    sources=['experiments/closure_ledger/pair_production_threshold_probe.py',
             'experiments/closure_ledger/antipodal_focusing_threshold_probe.py',
             'experiments/closure_ledger/self_consistent_throat_radius_probe.py',
             'experiments/closure_ledger/throat_dirac_spinor_probe.py',
             'experiments/closure_ledger/charge_conjugation_swap_probe.py',
             'experiments/closure_ledger/pair_creation_cobordism_capstone_probe.py',
             'experiments/closure_ledger/born_rule_equivariance_probe.py',
             'experiments/closure_ledger/guidance_law_from_5d_probe.py',
             'geometrodynamics/bell/bulk_identity.py','geometrodynamics/history/closure.py',
             'geometrodynamics/blackhole/interior.py',
             'experiments/closure_ledger/historical_foundation_audit_probe.py']
    return dict(label='HISTORICAL_FOUNDATION_DEPENDENCIES_REOPENED',
                baseline_commit='d42cb5e695f689dee254b51dbeabed3c21b6ad3c',
                scope='Foundation/event chain; not an exhaustive audit of all mass, flavor, QCD or QFT results',
                pair_threshold=threshold, finite_band_focus=focus_result, bell_readout=bell,
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=RUN)
    args=parser.parse_args();result=measure();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
