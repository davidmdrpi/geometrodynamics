from dataclasses import replace

import numpy as np
import pytest

from experiments.closure_ledger.coupled_mouth_throat_review import (
    FreeRotationCharge, radius_diagnostics, window_mean,
)
from geometrodynamics.waves.coupled_mouth_throat import CoupledMouthThroat, Parameters


def test_window_mean_integrates_nonuniform_linear_samples():
    t = np.array([0., .07, .4, .9, 1.8])
    assert window_mean(t, 2*t+3, .13, 1.12) == pytest.approx(4.25)
    with pytest.raises(ValueError, match='complete window'):
        window_mean(t, t, .5, 2.)


def test_full_period_mean_separates_oscillation_from_contraction():
    t = np.linspace(0, 4, 401)
    d = radius_diagnostics(t, 1+.1*np.cos(2*np.pi*t), 1.)
    assert d['instantaneous_minimum'] == pytest.approx(.9)
    assert d['minimum_full_radial_period_mean'] == pytest.approx(1., abs=1e-12)
    with pytest.raises(ValueError, match='complete positive-period'):
        radius_diagnostics(t, t, 5.)


def test_free_rotation_charge_has_canonical_sign_when_untilted():
    m = CoupledMouthThroat(replace(Parameters(), lmax=4, offset=0., coupling=0.))
    state = m.initial()
    assert FreeRotationCharge(m)(state) == pytest.approx(m.angular(state), abs=1e-12)


def test_tilted_free_field_correction_is_not_a_mechanical_support():
    m = CoupledMouthThroat(replace(Parameters(), lmax=6, offset=-.3, coupling=0.))
    _, states, _ = m.evolve(horizon=.5)
    assert np.max(abs(states[:, 2*m.n+3:2*m.n+6])) == 0
    j = np.array([m.angular(s) for s in states])
    assert np.ptp(j) > 1e-5
    assert np.ptp(states[:, -1]) > 1e-5
    assert np.ptp(j+states[:, -1]) < 1e-8
