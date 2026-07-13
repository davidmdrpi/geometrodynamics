"""
Unit tests for the wormhole-network traversal module (PR #216).

These verify the clock algebra, the derived transfer factor
U_BA = e^{i w Delta_BA}, local future-directedness alongside global
pastward reach, the closure offset, the advanced projection of the
traversal kernel, and frequency elasticity.  A toy greybody amplitude
stands in for the PR #215 Tangherlini solver (which the companion
probe uses in full).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.transaction import (
    NetworkMouth,
    NetworkThroat,
    closure_offset,
    emergent_frequency,
    network_confirmation,
    projected_kernel,
    traverse_throat,
)


def make_throat(delta=-7.08, tau_th=0.8, t_amp=0.9 + 0.1j,
                rate_A=1.0, rate_B=1.0, orient=(1, 1), phases=(0.0, 0.0)):
    a = NetworkMouth("A", psi=math.pi, link_id="L1", clock_rate=rate_A,
                     clock_offset=0.0, orientation=orient[0],
                     transfer_phase=phases[0])
    b = NetworkMouth("B", psi=math.pi, link_id="L1", clock_rate=rate_B,
                     clock_offset=delta, orientation=orient[1],
                     transfer_phase=phases[1])
    return NetworkThroat(a, b, tau_th=tau_th, t_AB=lambda w: t_amp)


def test_clock_algebra_roundtrip():
    m = NetworkMouth("A", psi=0.3, link_id="L", clock_rate=0.7,
                     clock_offset=-2.5)
    for t in (-3.0, 0.0, 4.2):
        assert m.global_time(m.local_time(t)) == pytest.approx(t)


def test_traversal_exit_time_and_local_duration():
    th = make_throat(delta=-7.08, tau_th=0.8)
    leg = traverse_throat(th, w=1.3, t_entry=math.pi)
    assert leg.t_end == pytest.approx(math.pi + 0.8 - 7.08)
    assert leg.local_duration == pytest.approx(0.8)
    assert leg.local_duration > 0            # future-directed in its clock
    assert leg.t_end < leg.t_start           # yet globally pastward


def test_derived_transfer_factor_is_clock_offset_phase():
    # the global-frame factor of the traversal must equal
    # t_AB * e^{-i w tau_th}; re-expressed against the global wave the
    # extra factor is U_BA = e^{i w Delta_BA}
    th = make_throat(delta=-5.0, tau_th=0.6, t_amp=0.8 + 0.0j)
    for w in (0.5, 1.7, 3.2):
        leg = traverse_throat(th, w, t_entry=0.0)
        # global-frame transfer: factor * e^{i w (t_exit - t_entry)}
        global_frame = leg.factor * np.exp(1j * w * (leg.t_end - leg.t_start))
        expected = th.t_AB(w) * th.U_BA(w)
        assert abs(global_frame - expected) < 1e-12
        assert abs(abs(th.U_BA(w)) - 1.0) < 1e-12   # unitary transfer


def test_full_confirmation_is_locally_forward_and_globally_advanced():
    d = math.pi
    tau = 0.8
    th = make_throat(delta=closure_offset(d, d, tau), tau_th=tau)
    tr = network_confirmation(th, w=2.0, t_emit=0.0, d_A=d, d_B=d)
    assert tr.locally_forward
    assert tr.globally_advanced
    assert tr.t_emerge < tr.t_emit < tr.t_absorb
    # time closure: the return lands exactly on the emission event
    assert tr.t_return == pytest.approx(0.0, abs=1e-12)


def test_closure_offset_formula():
    assert closure_offset(math.pi, math.pi, 0.8) == pytest.approx(
        -(2 * math.pi + 0.8))


def test_projected_kernel_is_advanced_with_greybody_weight():
    d_B = math.pi
    tau = 0.5
    th = make_throat(delta=-8.0, tau_th=tau, t_amp=0.7 + 0.2j)
    for w in (0.6, 1.9):
        pk = projected_kernel(th, w, d_B)
        assert pk['interval'] < 0                       # pastward segment
        # kernel = confirmation_weight * advanced_phase, |weight| = |t_AB|
        recon = pk['confirmation_weight'] * pk['advanced_phase']
        assert abs(recon - pk['kernel']) < 1e-12
        assert abs(abs(pk['confirmation_weight'])
                   - pk['greybody_magnitude']) < 1e-12


def test_phase_closure_makes_weight_real_positive():
    # with real t_AB and Delta chosen so w*Delta = 0 mod 2pi, the
    # confirmation weight is real positive: the coherent point
    w = 1.0
    delta = -2 * math.pi * 3          # w*Delta = -6 pi = 0 mod 2 pi
    th = make_throat(delta=delta, tau_th=0.4, t_amp=0.9 + 0.0j)
    pk = projected_kernel(th, w, d_B=1.0)
    weight = pk['confirmation_weight']
    assert abs(weight.imag) < 1e-12
    assert weight.real > 0


def test_frequency_elastic_iff_rates_match():
    th = make_throat(rate_A=1.0, rate_B=1.0)
    assert emergent_frequency(th, 2.4) == pytest.approx(2.4)
    th2 = make_throat(rate_A=1.0, rate_B=0.5)   # B deep in a well
    assert emergent_frequency(th2, 2.4) == pytest.approx(4.8)


def test_orientation_composition():
    # a (+1, -1) mouth pair flips the sign of the traversal factor
    th_pp = make_throat(orient=(1, 1), t_amp=1.0 + 0j)
    th_pm = make_throat(orient=(1, -1), t_amp=1.0 + 0j)
    w = 1.1
    f_pp = traverse_throat(th_pp, w, 0.0).factor
    f_pm = traverse_throat(th_pm, w, 0.0).factor
    assert abs(f_pp + f_pm) < 1e-12
