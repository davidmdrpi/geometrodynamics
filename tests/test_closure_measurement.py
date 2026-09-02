"""The closed-history measurement-dependence round, against
`docs/closure_measurement_dependence_prereg.md` (P1-P7 and the controls)."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import closure_measurement as cm


def test_P3_quadrature_matches_the_hand_derived_closed_form():
    for g in (0.3, 0.25 * math.pi, 1.0, 0.5 * math.pi, 2.0):
        assert cm.correlation(g) == pytest.approx(float(cm.closed_form_correlation(g)), abs=1e-4)
        W, Wc = cm.closure_weights(g), cm.closed_form_weights(g)
        for k in W:
            assert W[k] == pytest.approx(Wc[k], abs=1e-4)


def test_P3_pre_registered_values_and_not_cos():
    targets = {0.3: 0.82095, 0.25 * math.pi: 0.53557, 1.0: 0.39850, 0.5 * math.pi: 0.0, 2.0: -0.30350}
    for g, t in targets.items():
        E = float(cm.closed_form_correlation(g))
        assert E == pytest.approx(t, abs=1e-4)
        if g not in (0.5 * math.pi,):
            assert abs(E - math.cos(g)) > 0.1
    for g in (0.3, 1.0, 1.4):
        assert float(cm.closed_form_correlation(math.pi - g)) == pytest.approx(-float(cm.closed_form_correlation(g)), abs=1e-12)


def test_P2_no_signalling_marginals_are_exactly_one_half():
    for g in (0.2, 1.0, 2.5):
        for variant in ("abs", "signed", "pos"):
            m = cm.marginals(g, variant)
            assert m["P(A=+)"] == pytest.approx(0.5, abs=1e-12)
            assert m["P(B=+)"] == pytest.approx(0.5, abs=1e-12)


def test_P4_bell_violation_at_the_standard_angles_and_below_tsirelson():
    assert cm.chsh() == pytest.approx(2.1423, abs=2e-3)
    r = cm.chsh_max(n_grid=25)
    assert r["S_max"] > 2.0 and r["below_tsirelson"]


def test_P5_signed_density_is_the_bargmann_invariant_and_gives_cos():
    r = cm.bargmann_identity()
    assert r["max_residual"] < 1e-12
    assert r["signed_variant_is_cos"]
    assert cm.chsh(variant="signed") == pytest.approx(2.0 * math.sqrt(2.0), abs=1e-3)
    for g, frac in r["negative_weight_fraction_of_circle_like_outcomes"].items():
        assert frac == pytest.approx(g / (2 * math.pi), abs=2e-3)


def test_P6_strict_zero_variant_is_still_not_cos():
    assert cm.correlation(1.0, "pos") == pytest.approx(0.46495, abs=1e-4)
    assert cm.chsh(variant="pos") == pytest.approx(2.4649, abs=2e-3)


def test_P7_window_monte_carlo_converges_to_the_coarea_value():
    target = float(cm.closed_form_correlation(1.0))
    errs = [abs(cm.window_monte_carlo(1.0, eps, n=1000000) - target) for eps in (0.4, 0.2, 0.1)]
    assert errs[0] > errs[1] > 5e-4 or errs[2] < 5e-3
    assert errs[-1] < 5e-3


def test_P1_measurement_dependence_non_coplanar_supports_differ_and_coplanar_densities_differ():
    non = cm.measurement_dependence((0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (1.0, 0.0))
    assert not non["coplanar"] and non["total_variation"] == 1.0
    cop = cm.measurement_dependence((0.0, 0.0), (1.0, 0.0), (0.5, 0.0), (1.0, 0.0))
    assert cop["coplanar"] and cop["total_variation"] > 1e-2
    same = cm.measurement_dependence((0.0, 0.0), (1.0, 0.0), (0.0, 0.0), (1.0, 0.0))
    assert same["total_variation"] < 1e-12


def test_controls():
    assert cm.two_leg_loop_control()["closure_automatic"]
    loc = cm.local_detector_control(n=200000)
    assert loc["local_bound_respected"]
    w = cm.gaussian_width_control(n=200000)
    assert w["depends_on_sigma"]


def test_verdict_rule():
    assert cm.verdict(0.15, 2.14, 0.0, 1.0) == \
        "CLOSURE_INDUCES_SETTING_DEPENDENT_SOURCE_MEASURE_NO_SIGNALLING_NOT_BORN"
    assert cm.verdict(1e-6, 2.828, 0.0, 1.0) == "CLOSURE_REPRODUCES_QUANTUM_CORRELATIONS"
    assert cm.verdict(0.15, 2.14, 1e-3, 1.0) == "CLOSURE_SIGNALS"
    assert cm.verdict(0.15, 2.14, 0.0, 0.0) == "CLOSURE_INDUCES_NO_MEASUREMENT_DEPENDENCE"
