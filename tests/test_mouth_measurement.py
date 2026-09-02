"""The classical Born-rule round, against `docs/classical_born_prereg.md`."""

import math

import numpy as np
import pytest

from geometrodynamics.bulk import mouth_measurement as mm

TH = np.linspace(0.0, math.pi, 361)


def test_H1_symmetry_and_fibre_haar_do_not_select_born():
    c = mm.classification_theorem()
    assert c["all_satisfy_reversal"]
    assert c["every_f_is_realisable"]


def test_relative_configuration_shifts_psi_by_minus_two_phi_and_fixes_theta():
    q = np.array([0.3, 0.5, -0.2, 0.7]); q /= np.linalg.norm(q)
    phi = 0.4
    g = np.array([math.cos(phi), math.sin(phi), 0.0, 0.0])
    r0 = mm.relative_configuration([0, 0, 1.0], q)
    r1 = mm.relative_configuration([0, 0, 1.0], mm.qmul(g, q))
    assert r1["theta"] == pytest.approx(r0["theta"], abs=1e-12)
    shift = (r1["psi"] - r0["psi"] + math.pi) % (2 * math.pi) - math.pi
    assert abs(abs(shift) - 2 * phi) < 1e-10        # magnitude 2 phi; sign is convention


def test_C1_linear_threshold_family_never_reaches_born():
    r = mm.linear_family_best_fit()
    assert r["best_max_miss"] > 0.10
    assert r["closed_form_vs_arc_measure"] < 2e-4
    assert not r["reaches_born"]


def test_C2_intensity_detector_is_the_step():
    f = mm.induced_probability(mm.intensity_detector(), TH, n=500)
    assert np.all(f[TH < math.pi / 2] == 1.0) and np.all(f[TH > math.pi / 2] == 0.0)
    assert np.max(np.abs(f - mm.born(TH))) == pytest.approx(0.5, abs=1e-2)


def test_C3_two_harmonic_natural_weightings_miss_born_by_half():
    for row in mm.two_harmonic_natural_weightings():
        assert row["max_miss"] > 0.45


def test_C4_repository_winding_detector_is_phase_independent_step():
    r = mm.repository_winding_detector()
    assert r["phase_independent"] and r["induced_f_is_step"]


def test_C5_archimedes_is_born_exactly_at_kappa_one_and_only_there():
    assert mm.archimedes_uniformity()["uniform"]
    b = mm.born(TH)
    assert np.max(np.abs(mm.archimedes_probability(1.0, TH) - b)) < 1e-12
    for k, miss in ((0.9, 0.050), (1.1, 0.046), (0.5, 0.25), (2.0, 0.25)):
        assert np.max(np.abs(mm.archimedes_probability(k, TH) - b)) == pytest.approx(miss, abs=2e-3)
    for t in (0.3, 1.0, 2.0):
        assert mm.archimedes_monte_carlo(1.0, t) == pytest.approx(math.cos(t / 2) ** 2, abs=3e-3)


def test_basin_control_reproduces_born_but_is_tuned():
    f = mm.induced_probability(mm.tuned_born_basin(), TH[1:-1])
    assert np.max(np.abs(f - mm.born(TH[1:-1]))) < 2e-4


def test_measure_control_non_invariant_weight_changes_the_answer():
    assert mm.measure_control()["measure_matters"]


def test_reversal_control_on_every_reported_curve():
    curves = {"born": mm.born(TH), "arch": mm.archimedes_probability(1.0, TH),
              "lin": mm.linear_threshold_closed_form(0.8, TH)}
    for v in mm.reversal_control(curves).values():
        assert v < 1e-12


def test_verdict_rule():
    v = mm.verdict(0.109, [0.5, 0.5], 0.0, [0.05, 0.046, 0.25, 0.25])
    assert v == "BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW"
    assert mm.verdict(1e-4, [0.5], 0.0, [0.05]) == "BORN_RULE_DERIVED_FROM_CLASSICAL_BAM_MEASURE"
