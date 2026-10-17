"""Checks of the independent response solver and first-crossing termination."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'experiments/closure_ledger'))
import mouth_exchange_probe as probe
from geometrodynamics.waves.antipodal_kick import Parameters
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs


def test_bilinear_reference_noether_symmetry_and_uncoupled_null():
    for coupling in [1., 0.]:
        model = MovingMouthPairs(Parameters(cutoff=8, mass=.1, coupling=coupling))
        ref = probe.fixed_reference(model, .7)
        pair = np.asarray(ref['pair_mixed_impulse'])
        assert ref['max_mixed_ledger'] < 1e-9
        np.testing.assert_allclose(pair[:, 0, 1], -pair[:, 1, 1], atol=1e-10)
        if coupling == 0:
            np.testing.assert_array_equal(pair, 0.)
        else:
            assert abs(pair[-1, 0, 1]) > 1e-5


def test_actual_collision_stop_retains_last_admissible_state():
    model = MovingMouthPairs(Parameters(cutoff=4, mass=.1, coupling=0))
    s = model.initial(empty=True)
    s.angular[0, 1] = .1
    run = probe.moving(model, s, .02)
    lo, hi = run['stats']['stop']
    assert lo < np.pi/2-.65 <= hi
    last = run['records'][-1]
    assert last['t'] == lo
    x = np.asarray(last['x'])
    assert np.arccos(x[0]@x[1]) > .65
    assert hi-lo <= .02000001


def test_reference_tolerance_refinement():
    model = MovingMouthPairs(Parameters(cutoff=8, mass=.1))
    coarse = probe.fixed_reference(model, 0., rtol=2e-8, atol=2e-10)
    fine = probe.fixed_reference(model, 0.)
    np.testing.assert_allclose(coarse['pair_mixed_impulse'], fine['pair_mixed_impulse'], atol=1e-9, rtol=1e-5)
