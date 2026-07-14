"""
Unit tests for the wormhole-network traversal module (PR #216).

These verify the clock algebra, the derived transfer factor
U_BA = e^{i w Delta_BA}, local future-directedness alongside global
pastward reach, the closure offset, the advanced projection of the
traversal kernel, frequency elasticity, and the TWO-PORT throat: the
Fabry-Perot loop composition, its transparent-port reduction, composite
unitarity for unitary ports, resonant transmission, and the echo train.
Toy port amplitudes stand in for the PR #215 Tangherlini solver (which
the companion probe uses in full).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.transaction import (
    MouthPort,
    NetworkMouth,
    NetworkThroat,
    closure_offset,
    effective_green,
    emergence_train,
    emergent_frequency,
    loop_eigenvalue,
    network_confirmation,
    projected_kernel,
    transparent_port,
    traverse_throat,
)


def unitary_port(T=0.8, delta=0.3):
    """A symmetric-barrier unitary port: t = sqrt(T) e^{i d},
    r_out = r_in = i sqrt(1-T) e^{i d} (satisfies |t|^2 + |r|^2 = 1 and
    r t* + t r* = 0)."""
    tt = math.sqrt(T) * np.exp(1j * delta)
    rr = 1j * math.sqrt(1 - T) * np.exp(1j * delta)
    return MouthPort(t=lambda w: tt, r_out=lambda w: rr,
                     r_in=lambda w: rr)


def flat_port(t_amp):
    """A port with prescribed transmission and no interior reflection
    (no loops) - the v1 single-interface behavior."""
    return MouthPort(t=lambda w: t_amp, r_out=lambda w: 0.0 + 0j,
                     r_in=lambda w: 0.0 + 0j)


def make_throat(delta=-7.08, tau_th=0.8, port_A=None, port_B=None,
                rate_A=1.0, rate_B=1.0, orient=(1, 1), phases=(0.0, 0.0)):
    a = NetworkMouth("A", psi=math.pi, link_id="L1", clock_rate=rate_A,
                     clock_offset=0.0, orientation=orient[0],
                     transfer_phase=phases[0])
    b = NetworkMouth("B", psi=math.pi, link_id="L1", clock_rate=rate_B,
                     clock_offset=delta, orientation=orient[1],
                     transfer_phase=phases[1])
    return NetworkThroat(a, b, tau_th=tau_th,
                         port_A=port_A or flat_port(0.9 + 0.1j),
                         port_B=port_B or transparent_port())


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
    th = make_throat(delta=-5.0, tau_th=0.6, port_A=flat_port(0.8 + 0j))
    for w in (0.5, 1.7, 3.2):
        leg = traverse_throat(th, w, t_entry=0.0)
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
    th = make_throat(delta=-8.0, tau_th=tau, port_A=flat_port(0.7 + 0.2j))
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
    th = make_throat(delta=delta, tau_th=0.4, port_A=flat_port(0.9 + 0j))
    pk = projected_kernel(th, w, d_B=1.0)
    weight = pk['confirmation_weight']
    assert abs(weight.imag) < 1e-12
    assert weight.real > 0


def test_frequency_elastic_iff_rates_match():
    th = make_throat(rate_A=1.0, rate_B=1.0)
    assert emergent_frequency(th, 2.4) == pytest.approx(2.4)
    # B deep in a well (slow clock): the wave climbs out and REDSHIFTS
    th2 = make_throat(rate_A=1.0, rate_B=0.5)
    assert emergent_frequency(th2, 2.4) == pytest.approx(1.2)


def test_clock_rate_correct_exit_time():
    # equal rates rho: global traversal duration is tau_th / rho
    th = make_throat(delta=-9.0, tau_th=0.8, rate_A=0.5, rate_B=0.5,
                     port_A=unitary_port(), port_B=unitary_port())
    leg = traverse_throat(th, w=1.3, t_entry=2.0)
    assert leg.t_end == pytest.approx(2.0 + 0.8 / 0.5 - 9.0)
    assert leg.local_duration == pytest.approx(0.8)   # proper time


def test_loop_eigenvalue_and_effective_green():
    d, tau = math.pi, 0.8
    th = make_throat(delta=closure_offset(d, d, tau), tau_th=tau,
                     port_A=unitary_port(0.6, 0.2),
                     port_B=unitary_port(0.6, 0.2))
    for w in (0.7, 1.3, 2.9):
        lam = loop_eigenvalue(th, w, d, d)
        # value transport: at time closure the eigenvalue is the
        # composite throat amplitude alone - the carrier closes on the
        # throat's own scattering phase (delays live in arrival times)
        assert abs(lam - th.t_AB(w)) < 1e-12
        # passivity: the two-port throat cannot self-amplify
        assert abs(lam) <= 1.0 + 1e-12
        # the self-consistent resolvent
        g = effective_green(th, w, d, d)
        assert abs(g - 1.0 / (1.0 - lam)) < 1e-12
    # off time closure the displacement phase enters the eigenvalue
    th_off = make_throat(delta=closure_offset(d, d, tau) + 0.3,
                         tau_th=tau, port_A=unitary_port(0.6, 0.2),
                         port_B=unitary_port(0.6, 0.2))
    w = 1.3
    lam_off = loop_eigenvalue(th_off, w, d, d)
    assert abs(lam_off - th_off.t_AB(w) * np.exp(1j * w * 0.3)) < 1e-12


def test_orientation_composition():
    # a (+1, -1) mouth pair flips the sign of the traversal factor
    th_pp = make_throat(orient=(1, 1), port_A=flat_port(1.0 + 0j))
    th_pm = make_throat(orient=(1, -1), port_A=flat_port(1.0 + 0j))
    w = 1.1
    f_pp = traverse_throat(th_pp, w, 0.0).factor
    f_pm = traverse_throat(th_pm, w, 0.0).factor
    assert abs(f_pp + f_pm) < 1e-12


# ── the two-port throat: repeated-loop physics ──────────────────────────────

def test_transparent_ports_reduce_to_free_throat():
    th = make_throat(port_A=transparent_port(), port_B=transparent_port())
    for w in (0.5, 2.0):
        assert abs(th.t_AB(w) - 1.0) < 1e-14
        assert abs(th.r_AA(w)) < 1e-14


def test_loop_expansion_converges_to_composite():
    th = make_throat(port_A=unitary_port(0.5, 0.2),
                     port_B=unitary_port(0.7, -0.4), tau_th=0.9)
    w = 1.3
    partial = sum(th.loop_expansion(w, 60))
    assert abs(partial - th.t_AB(w)) < 1e-10


def test_composite_unitarity_for_unitary_ports():
    # flux conservation of the two-port composition off resonance and on
    th = make_throat(port_A=unitary_port(0.4, 0.15),
                     port_B=unitary_port(0.4, 0.15), tau_th=1.1)
    for w in np.linspace(0.3, 3.3, 25):
        flux = abs(th.t_AB(w)) ** 2 + abs(th.r_AA(w)) ** 2
        assert abs(flux - 1.0) < 1e-12


def test_resonant_transmission_identical_ports():
    # identical unitary ports: on interior resonance |t_net|^2 = 1
    # exactly, even for low single-port transmission
    T_single = 0.15
    th = make_throat(port_A=unitary_port(T_single, 0.3),
                     port_B=unitary_port(T_single, 0.3), tau_th=1.0)
    # resonance: arg(r_inA r_inB) + 2 w tau = 0 mod 2pi
    # r_in phases: pi/2 + 0.3 each -> resonance at 2w = 2pi - (pi + 0.6)
    w_res = (2 * math.pi - (math.pi + 0.6)) / 2.0
    assert abs(abs(th.t_AB(w_res)) ** 2 - 1.0) < 1e-12
    # far off resonance the transmission collapses well below T_single
    w_off = w_res + math.pi / 2
    assert abs(th.t_AB(w_off)) ** 2 < T_single


def test_echo_train_timing_and_damping():
    th = make_throat(delta=-9.0, tau_th=0.7,
                     port_A=unitary_port(0.6, 0.0),
                     port_B=unitary_port(0.6, 0.0))
    train = emergence_train(th, w=1.2, t_entry=2.0, kmax=3)
    damp = math.sqrt(0.4) * math.sqrt(0.4)     # |r_inA r_inB| per loop
    for k, leg in enumerate(train):
        assert leg.local_duration == pytest.approx((2 * k + 1) * 0.7)
        assert leg.t_end == pytest.approx(2.0 + (2 * k + 1) * 0.7 - 9.0)
        assert leg.local_duration > 0          # every echo future-directed
        if k:
            ratio = abs(leg.factor) / abs(train[k - 1].factor)
            assert ratio == pytest.approx(damp, rel=1e-9)
