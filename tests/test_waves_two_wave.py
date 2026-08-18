"""Tests for the dynamical two-wave invariant."""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.waves.field_solve import image_field
from geometrodynamics.waves.throat_operator import MouthPair, gamma_at
from geometrodynamics.waves.two_source import geodesic, mouth_positions
from geometrodynamics.waves.two_wave import (
    GaussianPulse, RetardedGrid, SOURCE_GAP, OBSERVER_REACH, TWO_PI,
    TwoWaveSetup, arrival_directions, circle_point, contract_stress,
    gamma_omega, green_omega, green_omega_derivatives,
    measure_multipath_destroys_the_collinear_null,
    measure_the_arrivals_are_the_branch_ledger_with_maslov_signs,
    measure_the_caustic_is_where_wkb_stops,
    measure_the_improved_stress_tensor_is_traceless,
    measure_the_low_frequency_limit_recovers_the_tomography,
    measure_the_only_tail_is_the_throats,
    measure_the_solved_field_satisfies_the_conformal_wave_equation,
    measure_the_solver_reproduces_the_closed_form_free_field,
    measure_the_wkb_collinear_head_on_result_is_recovered,
    cross_stress_tensor,
    measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth,
    measure_the_interference_tensor_is_largest_where_the_invariant_is_null,
    normalized_invariant, orthonormal_frame, radial_frame_data, solve_field,
    stress_tensor, superpose, trace_of, two_leg_channels, wkb_invariant,
    working_pair, zero_like)


def _grid() -> RetardedGrid:
    return RetardedGrid(n=1 << 15, span=300.0, eps=0.05)


# ── the kernel ──────────────────────────────────────────────────────────────
def test_the_complex_green_function_matches_the_real_one():
    for chi in (0.4, 1.3, 2.7):
        for lam in (0.3, 0.8, 2.5, 6.0):
            w = math.sqrt(lam)
            assert green_omega(chi, complex(w, 0.0)) == pytest.approx(
                float(gamma_at(lam, chi)[0, 1].real), abs=1e-14)


def test_gamma_omega_matches_gamma_at_on_the_real_axis():
    for lam in (0.2, 0.9, 3.1, 7.7):
        w = math.sqrt(lam)
        assert np.abs(gamma_omega(complex(w, 0.0), 1.3)
                      - gamma_at(lam, 1.3)).max() < 1e-14


def test_the_green_derivatives_are_the_analytic_ones():
    """Checked against a high-order finite difference in χ, once."""
    w = complex(3.7, 0.05)
    for chi in (0.5, 1.6, 2.9):
        g0, g1, g2 = green_omega_derivatives(chi, w)
        h = 1e-5
        d1 = (green_omega(chi + h, w) - green_omega(chi - h, w)) / (2 * h)
        d2 = (green_omega(chi + h, w) - 2 * green_omega(chi, w)
              + green_omega(chi - h, w)) / h ** 2
        assert abs(g0 - green_omega(chi, w)) < 1e-15
        assert abs(g1 - d1) < 1e-6 * max(1.0, abs(d1))
        assert abs(g2 - d2) < 1e-4 * max(1.0, abs(d2))


def test_the_antipode_is_a_removable_limit():
    w = complex(2.3, 0.05)
    near = green_omega(math.pi - 1e-7, w)
    at = green_omega(math.pi, w)
    assert abs(near - at) < 1e-8
    assert abs(at - w / (4.0 * math.pi * np.sin(math.pi * w))) < 1e-12


def test_the_pulse_spectrum_is_array_aware():
    p = GaussianPulse(carrier=5.0, width=0.2)
    om = np.array([0.0, 1.0 + 0.05j, -3.0 + 0.05j])
    got = p.spectrum(om)
    assert got.shape == om.shape
    for k, w in enumerate(om):
        assert got[k] == pytest.approx(complex(p.spectrum(complex(w))),
                                       abs=1e-15)


# ── geometry ────────────────────────────────────────────────────────────────
def test_the_radial_gradient_is_a_unit_vector_pointing_away():
    x = circle_point(1.1)
    y = circle_point(0.0)
    frame = orthonormal_frame(x)
    chi, n = radial_frame_data(x, y, frame)
    assert chi == pytest.approx(geodesic(x, y), abs=1e-14)
    assert np.linalg.norm(n) == pytest.approx(1.0, abs=1e-12)
    # moving along +n increases the distance from y
    step = 1e-5
    moved = math.cos(step) * x + math.sin(step) * (frame.T @ n)
    assert geodesic(moved, y) > chi


def test_the_collinear_and_head_on_directions_are_exact():
    obs_c = circle_point(SOURCE_GAP + OBSERVER_REACH)
    obs_h = circle_point(SOURCE_GAP - OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    assert arrival_directions(obs_c, a, b) == pytest.approx(1.0, abs=1e-12)
    assert arrival_directions(obs_h, a, b) == pytest.approx(-1.0, abs=1e-12)


def test_the_source_circle_avoids_the_mouths():
    c1, c2 = mouth_positions(1.3)
    for theta in (0.0, SOURCE_GAP, SOURCE_GAP + OBSERVER_REACH,
                  SOURCE_GAP - OBSERVER_REACH):
        p = circle_point(theta)
        assert geodesic(p, c1) > 0.5
        assert geodesic(p, c2) > 0.5


# ── the solve ───────────────────────────────────────────────────────────────
def test_the_free_solve_is_the_image_sum():
    width, chi = 0.06, 1.7
    amp = 1.0 / (width * math.sqrt(TWO_PI))
    pulse = GaussianPulse(amplitude=amp, carrier=0.0, width=width)
    obs = circle_point(chi)
    setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                         source_b=circle_point(0.0), observer=obs,
                         pulse_a=pulse, pulse_b=pulse, grid=_grid(),
                         with_throat=False)
    sol = solve_field(setup, setup.source_a, pulse)
    for t_arr in (chi, TWO_PI - chi, TWO_PI + chi):
        i = int(round(t_arr / setup.grid.dt))
        assert sol["phi"][i] == pytest.approx(
            image_field(chi, float(setup.grid.times[i]), width=width,
                        n_images=10), abs=1e-12)


def test_the_inversion_uses_the_forward_kernel():
    """``ifft`` would evaluate the transform at ``−t`` and give ~0 for a
    retarded field — the sign that was actually caught this way."""
    grid = _grid()
    pulse = GaussianPulse(carrier=0.0, width=0.06)
    obs = circle_point(1.7)
    setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                         source_b=circle_point(0.0), observer=obs,
                         pulse_a=pulse, pulse_b=pulse, grid=grid,
                         with_throat=False)
    sol = solve_field(setup, setup.source_a, pulse)
    ts = grid.times
    # well before the light cone at χ = 1.7 — "well" because the Gaussian
    # source is not compactly supported, so its own wing legitimately precedes
    # the geometric arrival by a few widths and that is the source's doing,
    # not the propagator's
    before = np.abs(sol["phi"][(ts > 0.05) & (ts < 1.0)]).max()
    after = np.abs(sol["phi"][(ts > 1.6) & (ts < 1.8)]).max()
    assert before < 1e-6 * after


def test_the_wave_equation_holds_with_and_without_the_throat():
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    for throat in (False, True):
        pulse = GaussianPulse(carrier=16.0, width=0.10)
        setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                             source_b=circle_point(SOURCE_GAP), observer=obs,
                             pulse_a=pulse, pulse_b=pulse, grid=grid,
                             with_throat=throat)
        sol = solve_field(setup, setup.source_a, pulse)
        ts = grid.times
        sl = (ts > 0.5) & (ts < 9.0)
        resid = sol["dtt"][sl] - (sol["laplacian"][sl] - sol["phi"][sl])
        assert np.abs(resid).max() < 1e-12 * np.abs(sol["dtt"][sl]).max()


def test_the_throat_changes_the_field():
    grid = _grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    pulse = GaussianPulse(carrier=8.0, width=0.12)
    sols = []
    for throat in (False, True):
        setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                             source_b=circle_point(SOURCE_GAP), observer=obs,
                             pulse_a=pulse, pulse_b=pulse, grid=grid,
                             with_throat=throat)
        sols.append(solve_field(setup, setup.source_a, pulse))
    ts = grid.times
    sl = ts < 8.0
    assert np.abs(sols[1]["phi"][sl] - sols[0]["phi"][sl]).max() > 1e-4


# ── the stress tensor ───────────────────────────────────────────────────────
def test_the_stress_tensor_is_traceless_on_the_solved_field():
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    pulse = GaussianPulse(carrier=16.0, width=0.10)
    setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                         source_b=circle_point(SOURCE_GAP), observer=obs,
                         pulse_a=pulse, pulse_b=pulse, grid=grid,
                         with_throat=True)
    sol = solve_field(setup, setup.source_a, pulse)
    i = int(round(geodesic(obs, setup.source_a) / grid.dt))
    t = stress_tensor(sol, i)
    assert abs(trace_of(t)) < 1e-12 * np.abs(t).max()
    assert t[0, 0] > 0.0                              # positive energy density


def test_the_trace_would_be_vacuous_on_shell():
    """``T^μ_μ = φ(□φ − φ)``: substituting the on-shell ``□φ`` makes it vanish
    for *any* input, which is why the solved value is used instead."""
    grid = _grid()
    obs = circle_point(1.4)
    pulse = GaussianPulse(carrier=8.0, width=0.12)
    setup = TwoWaveSetup(pair=working_pair(), source_a=circle_point(0.0),
                         source_b=circle_point(0.0), observer=obs,
                         pulse_a=pulse, pulse_b=pulse, grid=grid,
                         with_throat=False)
    sol = solve_field(setup, setup.source_a, pulse)
    i = int(round(1.4 / grid.dt))
    phi = float(sol["phi"][i])
    box = float(sol["laplacian"][i] - sol["dtt"][i])
    assert trace_of(stress_tensor(sol, i)) == pytest.approx(
        phi * (box - phi), abs=1e-13)


def test_contract_stress_is_symmetric():
    rng = np.random.default_rng(3)
    a = rng.normal(size=(4, 4))
    b = rng.normal(size=(4, 4))
    a = a + a.T
    b = b + b.T
    assert contract_stress(a, b) == pytest.approx(contract_stress(b, a),
                                                  rel=1e-12)


def test_the_wkb_invariant_is_zero_collinear_and_four_head_on():
    e = 3.0
    k1 = np.array([e, e, 0.0, 0.0])
    k_par = np.array([e, e, 0.0, 0.0])
    k_anti = np.array([e, -e, 0.0, 0.0])
    assert wkb_invariant(1.0, k1, 1.0, k_par) == pytest.approx(0.0, abs=1e-20)
    assert wkb_invariant(1.0, k1, 1.0, k_anti) == pytest.approx(
        (2.0 * e * e) ** 2, rel=1e-12)


# ── the physics ─────────────────────────────────────────────────────────────
def test_the_collinear_configuration_is_null_and_head_on_is_four():
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    for theta, target, tol in ((SOURCE_GAP + OBSERVER_REACH, 0.0, 1e-4),
                               (SOURCE_GAP - OBSERVER_REACH, 4.0, 5e-2)):
        obs = circle_point(theta)
        a, b = circle_point(0.0), circle_point(SOURCE_GAP)
        chi_a, chi_b = geodesic(obs, a), geodesic(obs, b)
        t_star = 3.0
        setup = TwoWaveSetup(
            pair=working_pair(), source_a=a, source_b=b, observer=obs,
            pulse_a=GaussianPulse(carrier=24.0, width=0.10,
                                  t0=t_star - chi_a),
            pulse_b=GaussianPulse(carrier=24.0, width=0.10,
                                  t0=t_star - chi_b),
            grid=grid, with_throat=False)
        got = normalized_invariant(setup, t_star, 0.2)["invariant"]
        assert abs(got - target) < tol


def test_the_winding_image_turns_collinear_into_head_on():
    """The same two sources at the same event, with ``A`` on its long-way
    branch: the arrival direction reverses and the invariant goes 0 → 4."""
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    chi_a, chi_b = geodesic(obs, a), geodesic(obs, b)
    t_star = 3.0
    setup = TwoWaveSetup(
        pair=working_pair(), source_a=a, source_b=b, observer=obs,
        pulse_a=GaussianPulse(carrier=24.0, width=0.10,
                              t0=t_star - (TWO_PI - chi_a)),
        pulse_b=GaussianPulse(carrier=24.0, width=0.10, t0=t_star - chi_b),
        grid=grid, with_throat=False)
    assert normalized_invariant(setup, t_star, 0.2)["invariant"] == (
        pytest.approx(4.0, abs=1e-2))


def test_the_throat_arrival_matches_the_mouth_geometry():
    grid = RetardedGrid(n=1 << 17, span=600.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    c1, c2 = mouth_positions(1.3)
    legs = [(geodesic(b, cj) + geodesic(ci, obs), i)
            for i, ci in enumerate((c1, c2)) for cj in (c1, c2)]
    delay, which = min(legs)
    t_star = 3.0
    setup = TwoWaveSetup(
        pair=working_pair(), source_a=a, source_b=b, observer=obs,
        pulse_a=GaussianPulse(carrier=24.0, width=0.10,
                              t0=t_star - geodesic(obs, a)),
        pulse_b=GaussianPulse(carrier=24.0, width=0.10, t0=t_star - delay),
        grid=grid, with_throat=True)
    got = normalized_invariant(setup, t_star, 0.2)["invariant"]
    predicted = (1.0 - arrival_directions(obs, a, (c1, c2)[which])) ** 2
    assert got == pytest.approx(predicted, rel=5e-3)
    assert 0.1 < got < 3.0                     # neither null nor head-on


def test_the_free_field_has_no_tail():
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src = circle_point(0.0)
    pulse = GaussianPulse(carrier=16.0, width=0.10)
    setup = TwoWaveSetup(pair=working_pair(), source_a=src, source_b=src,
                         observer=obs, pulse_a=pulse, pulse_b=pulse,
                         grid=grid, with_throat=False)
    sol = solve_field(setup, src, pulse)
    ts = grid.times
    chi = geodesic(obs, src)
    mask = (ts > 0.2) & (ts < 5.0)
    for t_arr in (chi, TWO_PI - chi, TWO_PI + chi):
        mask &= np.abs(ts - t_arr) > 0.6
    peak = np.abs(sol["phi"][ts < 6.0]).max()
    assert np.abs(sol["phi"][mask]).max() < 1e-6 * peak


def test_the_caustic_saturates_instead_of_diverging():
    for c in (8.5, 16.5, 32.5):
        w = complex(c, 0.0)
        assert abs(green_omega(math.pi, w)) == pytest.approx(
            c / (4.0 * math.pi), rel=1e-12)
        # and the ratio to WKB depends only on ω·e
        for scaled in (0.5, 1.0, 2.0):
            e = scaled / c
            ratio = (abs(green_omega(math.pi - e, w))
                     * 4.0 * math.pi * math.sin(e))
            assert ratio == pytest.approx(abs(math.sin(scaled)), rel=1e-10)


# ── the measurements ────────────────────────────────────────────────────────
def test_measure_the_solver_reproduces_the_closed_form_free_field():
    r = measure_the_solver_reproduces_the_closed_form_free_field()
    assert r["the_two_constructions_agree"]
    assert r["the_signs_alternate"]


def test_measure_the_solved_field_satisfies_the_conformal_wave_equation():
    r = measure_the_solved_field_satisfies_the_conformal_wave_equation()
    assert r["the_equation_holds"]
    assert r["loewner_margin"] > 0.3


def test_measure_the_improved_stress_tensor_is_traceless():
    r = measure_the_improved_stress_tensor_is_traceless()
    assert r["the_tensor_is_traceless"]


def test_measure_the_wkb_collinear_head_on_result_is_recovered():
    r = measure_the_wkb_collinear_head_on_result_is_recovered()
    assert r["the_directions_are_exactly_parallel"]
    assert r["the_directions_are_exactly_antiparallel"]
    assert r["head_on_converges_to_four"]
    assert r["collinear_converges_to_zero"]
    assert r["the_contour_needs_eps_above_the_spacing"]


def test_measure_multipath_destroys_the_collinear_null():
    r = measure_multipath_destroys_the_collinear_null()
    assert r["the_direct_pair_is_null"]
    assert r["the_winding_image_reads_head_on"]
    assert r["the_throat_matches_its_geometry"]
    assert r["the_control_has_no_second_arrival"]
    assert r["loewner_margin"] > 0.3


def test_measure_the_arrivals_are_the_branch_ledger_with_maslov_signs():
    r = measure_the_arrivals_are_the_branch_ledger_with_maslov_signs()
    assert r["the_free_signs_alternate"]
    assert r["the_free_arrivals_are_sharp"]
    assert r["the_throat_onset_is_causal"]
    assert r["the_throat_arrivals_are_new"]


def test_measure_the_only_tail_is_the_throats():
    r = measure_the_only_tail_is_the_throats()
    assert r["the_free_field_has_no_tail"]
    assert r["the_throat_has_one"]


def test_measure_the_caustic_is_where_wkb_stops():
    r = measure_the_caustic_is_where_wkb_stops()
    assert r["the_saturation_is_linear_in_omega"]
    assert r["the_ratio_collapses_in_omega_times_e"]


def test_measure_the_low_frequency_limit_recovers_the_tomography():
    r = measure_the_low_frequency_limit_recovers_the_tomography()
    assert r["the_bridge_closes"]
    assert r["W_error"] < 1e-3


# ── the (i,j) audit and the β = 0 control ───────────────────────────────────
def test_the_two_leg_channels_are_labelled_by_the_exit_mouth():
    """All four ``(i,j)`` paths, and the predicted invariant depending only on
    the exit mouth ``i``."""
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    src = circle_point(SOURCE_GAP)
    chans = two_leg_channels(obs, src, 1.3)
    assert len(chans) == 4
    assert sorted(c["delay"] for c in chans) == [c["delay"] for c in chans]
    by_exit = {}
    for c in chans:
        by_exit.setdefault(c["exit_mouth"], set()).add(
            round(c["predicted_invariant"], 12))
    assert set(by_exit) == {1, 2}
    for vals in by_exit.values():
        assert len(vals) == 1           # depends only on the exit mouth
    assert len({next(iter(v)) for v in by_exit.values()}) == 2


def test_the_field_picks_the_predicted_exit_mouth():
    """Two channels with *different* predicted values, both matched."""
    grid = RetardedGrid(n=1 << 17, span=600.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    chans = two_leg_channels(obs, b, 1.3)
    width, carrier, t_star = 0.035, 60.0, 3.0
    for ch in (chans[0], chans[-1]):
        setup = TwoWaveSetup(
            pair=working_pair(), source_a=a, source_b=b, observer=obs,
            pulse_a=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - geodesic(obs, a)),
            pulse_b=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - ch["delay"]),
            grid=grid, with_throat=True)
        got = normalized_invariant(setup, t_star, 2.0 * width)["invariant"]
        assert got == pytest.approx(ch["predicted_invariant"], rel=3e-3)


def test_disconnected_mouths_give_the_same_invariant():
    """``β = 0`` is the control PR #258's review taught this arc to run: the
    invariant is set by the exit mouth, not by whether the mouths are joined."""
    grid = RetardedGrid(n=1 << 17, span=600.0, eps=0.05)
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    ch = two_leg_channels(obs, b, 1.3)[0]
    width, carrier, t_star = 0.035, 60.0, 3.0
    vals = []
    for beta in (0.0, 0.06, 0.26):
        setup = TwoWaveSetup(
            pair=MouthPair(1.3, 0.30, 0.35, beta), source_a=a, source_b=b,
            observer=obs,
            pulse_a=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - geodesic(obs, a)),
            pulse_b=GaussianPulse(carrier=carrier, width=width,
                                  t0=t_star - ch["delay"]),
            grid=grid, with_throat=True)
        vals.append(normalized_invariant(setup, t_star,
                                         2.0 * width)["invariant"])
    spread = max(vals) - min(vals)
    assert spread < 1e-5
    # …and five orders below what separates the two exit mouths
    chans = two_leg_channels(obs, b, 1.3)
    sep = (max(c["predicted_invariant"] for c in chans)
           - min(c["predicted_invariant"] for c in chans))
    assert spread < 1e-4 * sep


# ── the interference stress tensor ──────────────────────────────────────────
def test_superposition_is_exact():
    grid = _grid()
    obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    pulse = GaussianPulse(carrier=8.0, width=0.12)
    setup = TwoWaveSetup(pair=working_pair(), source_a=a, source_b=b,
                         observer=obs, pulse_a=pulse, pulse_b=pulse, grid=grid,
                         with_throat=True)
    sa = solve_field(setup, a, pulse)
    sb = solve_field(setup, b, pulse)
    tot = superpose(sa, sb)
    for k in ("phi", "dt", "dtt", "grad", "hess"):
        assert np.abs(tot[k] - (sa[k] + sb[k])).max() == 0.0


def test_the_cross_stress_tensor_is_bilinear_and_traceless():
    grid = RetardedGrid(n=1 << 16, span=400.0, eps=0.05)
    obs = circle_point(SOURCE_GAP - OBSERVER_REACH)
    a, b = circle_point(0.0), circle_point(SOURCE_GAP)
    t_star = 3.0
    setup = TwoWaveSetup(
        pair=working_pair(), source_a=a, source_b=b, observer=obs,
        pulse_a=GaussianPulse(carrier=24.0, width=0.10,
                              t0=t_star - geodesic(obs, a)),
        pulse_b=GaussianPulse(carrier=24.0, width=0.10,
                              t0=t_star - geodesic(obs, b)),
        grid=grid, with_throat=False)
    sa = solve_field(setup, a, setup.pulse_a)
    sb = solve_field(setup, b, setup.pulse_b)
    i = int(round(t_star / grid.dt))
    dt = cross_stress_tensor(sa, sb, i)
    assert np.abs(dt).max() > 1e-3
    assert abs(trace_of(dt)) < 1e-12 * np.abs(dt).max()
    # switching a source off kills it exactly
    assert np.abs(cross_stress_tensor(sa, zero_like(sb), i)).max() == 0.0
    assert np.abs(cross_stress_tensor(zero_like(sa), sb, i)).max() == 0.0
    # and it is the cross term of the same functional
    total = stress_tensor(superpose(sa, sb), i)
    assert np.abs(total - stress_tensor(sa, i) - stress_tensor(sb, i)
                  - dt).max() < 1e-12


def test_the_interference_is_maximal_where_the_invariant_is_null():
    """The collinear configuration nulls ``T_A:T_B`` and *maximises* ``ΔT``."""
    r = measure_the_interference_tensor_is_largest_where_the_invariant_is_null()
    assert abs(r["collinear_invariant"]) < 1e-5
    assert r["collinear_interference"] == pytest.approx(2.0, abs=1e-2)
    assert r["head_on_invariant"] == pytest.approx(4.0, abs=1e-2)
    assert r["head_on_interference"] < 0.7 * r["collinear_interference"]


def test_measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth():
    r = measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth()
    assert r["the_prediction_depends_only_on_the_exit_mouth"]
    assert r["the_field_picks_the_right_one"]
    assert r["the_invariant_is_beta_independent"]
    assert r["every_sweep_point_is_inside_the_cone"]


def test_measure_the_interference_tensor_is_largest_where_the_invariant_is_null():
    r = measure_the_interference_tensor_is_largest_where_the_invariant_is_null()
    assert r["delta_T_is_traceless"]
    assert r["delta_T_vanishes_when_a_source_is_removed"]
    assert r["the_interference_is_maximal_where_the_invariant_is_null"]
