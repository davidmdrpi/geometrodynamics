"""Klein time-attachment tests. Freeze ``bb6ff2c``; corrected after review."""

import numpy as np
import pytest

from geometrodynamics.bulk import klein_time as kt
from experiments.closure_ledger import klein_time_probe as probe


@pytest.fixture(scope="module")
def lat():
    return kt.Lattice(kt.N_FROZEN)


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def test_unit_courant_lattice_has_integer_frequencies(lat):
    """Leapfrog at ``dt = dchi`` has eigenphases exactly ``k pi / N``."""
    phases = np.sort(np.abs(np.angle(np.linalg.eigvals(lat.step))))
    expected = np.sort(np.repeat(np.arange(1, lat.n) * np.pi / lat.n, 2))
    assert np.allclose(phases, expected, atol=1e-9)


def test_mode_parity_under_the_antipode(lat):
    for k in range(1, lat.n):
        v = lat.mode(k)
        assert np.allclose(lat.flip @ v, (-1) ** (k + 1) * v, atol=1e-12)


def test_wronskian_is_conserved(lat):
    w = lat.wronskian
    assert np.abs(lat.step.T @ w @ lat.step - w).max() < 1e-12


@pytest.mark.parametrize("m", [0, 1, 2, 7, 24, 31])
@pytest.mark.parametrize("eta", [1, -1])
def test_admissible_fields_are_lagrangian(lat, m, eta):
    r = kt.halving_report(lat, m, eta)
    assert r["involution_error"] < 1e-9
    assert r["antisymplectic_error"] < 1e-9
    assert r["dim_plus"] == r["dim_minus"] == lat.dim
    assert r["omega_on_plus"] < 1e-9 and r["omega_on_minus"] < 1e-9
    assert r["pairing_min_singular"] > 1e-6


def test_orientable_quotient_stays_symplectic(lat):
    c = kt.orientable_control(lat)
    assert c["dim_fixed"] == 24 and c["restricted_rank"] == 24


def test_halving_is_not_an_artefact_of_the_frozen_size():
    for n in (9, 16, 25):
        small = kt.Lattice(n)
        r = kt.halving_report(small, 5, 1)
        assert r["dim_plus"] == r["dim_minus"] == n - 1
        assert r["omega_on_plus"] < 1e-9


@pytest.mark.parametrize("m", [3, 7, 24])
@pytest.mark.parametrize("eta", [1, -1])
def test_marching_reproduces_the_identification(lat, m, eta):
    """Independent of the state-space deck matrix: loop, then compare levels."""
    rng = np.random.default_rng(kt.SEED + m)
    r = kt.direct_admissibility(lat, m, eta, rng)
    assert r["admissible_defect"] < 1e-9
    assert r["generic_defect"] > 0.1


def test_march_agrees_with_matrix_power(lat):
    rng = np.random.default_rng(kt.SEED)
    state = rng.standard_normal(2 * lat.dim)
    levels = kt.march(lat, state, 11)
    final = np.linalg.matrix_power(lat.step, 10) @ state
    assert np.allclose(np.concatenate([levels[10], levels[11]]), final,
                       atol=1e-9)


@pytest.mark.parametrize("m", [0, 7, 24])
def test_deck_is_a_conjugation(lat, m):
    r = kt.conjugation_report(lat, m, 1)
    assert r["j_squared_error"] < 1e-9
    assert r["antilinear_error"] < 1e-9
    assert r["j_plus_in_minus"] < 1e-9
    assert r["metric_definite_min"] > 1e-6


def test_orientable_deck_is_complex_linear(lat):
    """Control: the time-preserving antipode commutes with ``j``."""
    j = kt.complex_structure(lat)
    p = kt.orientable_deck(lat)
    assert np.abs(p @ j - j @ p).max() < 1e-9


@pytest.mark.parametrize("m", sorted(kt.TABLE_P3))
def test_two_cap_modes_match_the_frozen_table(lat, m):
    r = kt.two_cap_report(lat, m)
    untwisted, twisted = kt.TABLE_P3[m]
    assert r["pp"]["modes"] == untwisted
    assert r["pm"]["modes"] == twisted
    assert r["mm"] == r["pp"] and r["mp"] == r["pm"]


@pytest.mark.parametrize("m", sorted(kt.TABLE_P3))
def test_loop_phase_is_zero_or_pi(lat, m):
    """The surviving modes are exactly those with loop phase 0 or pi."""
    r = kt.two_cap_report(lat, m)
    for k in range(1, lat.n):
        phase = kt.loop_phase(k, m, lat.n)
        closes = min(phase, 2 * np.pi - phase) < 1e-9
        antic = abs(phase - np.pi) < 1e-9
        assert (k in r["pp"]["modes"]) == closes
        assert (k in r["pm"]["modes"]) == antic


def test_half_transit_history_splits_even_and_odd(lat):
    r = kt.two_cap_report(lat, lat.n)
    assert all(k % 2 == 0 for k in r["pp"]["modes"])
    assert all(k % 2 == 1 for k in r["pm"]["modes"])
    assert len(r["pp"]["modes"]) + len(r["pm"]["modes"]) == lat.dim


@pytest.mark.parametrize("m", [5, 7, 11, 13])
def test_incommensurate_history_has_no_free_field(lat, m):
    for etas in ((1, 1), (1, -1)):
        assert kt.two_cap_space(lat, m, *etas).shape[1] == 0


@pytest.mark.parametrize("m", sorted(kt.TABLE_P3))
def test_periodic_operator_is_an_independent_count(lat, m):
    k = kt.kernel_report(lat, m)
    assert k["nullity"] == 2 * len(kt.TABLE_P3[m][0])
    assert k["invariant_nullity"] == len(kt.TABLE_P3[m][0])


def test_generic_history_is_slaved_to_its_source(lat):
    s = kt.slaved_field_report(lat, 7)
    assert s["smallest_singular_value"] > 1e-6
    assert s["green_symmetry_error"] < 1e-9
    assert s["green_deck_commutator"] < 1e-9
    assert s["field_before_source_norm"] > 0.1
    assert s["before_after_ratio_off_cap"] == pytest.approx(1.0, abs=1e-9)


def test_resonant_history_fredholm(lat):
    f = kt.fredholm_report(lat, 24, np.random.default_rng(kt.SEED))
    assert f["residual_orthogonal_source"] < 1e-8
    assert f["residual_generic_source"] == pytest.approx(
        f["kernel_component_norm"], abs=1e-8)


# ---- Review of 5b79975 ----------------------------------------------------

def test_fixed_geometry_refinement_finds_survivors_the_cutoff_hid():
    """R1. Holding tau = 7 pi / 48 fixed, refinement admits homogeneous modes."""
    rows = {r["n"]: r for r in kt.refinement_report()}
    assert (rows[24]["untwisted"], rows[24]["twisted"]) == ((), ())
    assert rows[48]["twisted"] == (24,)
    assert rows[72]["untwisted"] == (48,)
    assert rows[96]["twisted"] == (24, 72)
    for r in rows.values():
        assert r["untwisted"] == r["continuum_untwisted"]
        assert r["twisted"] == r["continuum_twisted"]


def test_reviewer_counterexample_satisfies_both_caps():
    """R1. ``sin(48 chi) sin(48 t)`` is admissible at tau = 7 pi / 48."""
    c = kt.continuum_counterexample(7, 48, twisted=False)
    assert (c["k"], c["quadrature"]) == (48, "sin")
    assert c["amplitude"] > 0.5
    assert c["cap_0_defect"] < 1e-9 and c["cap_tau_defect"] < 1e-9


def test_counterexample_check_can_fail():
    """A mode that should not survive is rejected by the same marching test."""
    lat = kt.Lattice(96)
    chi = np.arange(1, 96) * np.pi / 96
    state = np.concatenate([np.zeros(95), np.sin(47 * chi)
                            * np.sin(47 * np.pi / 96)])
    deck = kt.klein_deck(lat, 28, 1)
    assert np.abs(deck @ state - state).max() > 0.1


@pytest.mark.parametrize("p,q", [(1, 3), (2, 5), (3, 7), (5, 11), (7, 48)])
def test_every_rational_spacing_has_an_untwisted_mode(p, q):
    c = kt.continuum_counterexample(p, q, twisted=False)
    assert c["exists"] and c["k"] == q
    assert c["cap_0_defect"] < 1e-9 and c["cap_tau_defect"] < 1e-9


@pytest.mark.parametrize("p,q,exists", [(1, 3, False), (2, 5, False),
                                        (3, 8, True), (7, 48, True)])
def test_twisted_modes_need_an_even_denominator(p, q, exists):
    c = kt.continuum_counterexample(p, q, twisted=True)
    assert c["exists"] is exists
    if exists:
        assert c["k"] == q // 2
        assert c["cap_0_defect"] < 1e-9 and c["cap_tau_defect"] < 1e-9


def test_continuum_modes_agree_with_the_lattice_formula():
    for n, m in ((24, 7), (24, 12), (24, 24), (48, 14), (96, 28)):
        for twisted in (False, True):
            assert kt.continuum_modes(m, 2 * n, n, twisted) == \
                kt.predicted_modes(n, m, twisted)


def test_divisor_diagnostic_separates_spacings():
    """Diagnostic only: exact zero at 7/48, a floor at the golden ratio."""
    assert kt.continuum_gap(7 / 48, 1000)["min_gap"] == 0.0
    golden = kt.continuum_gap((5 ** 0.5 - 1) / 2, 100_000)
    # Corrected after review of c433f91: k = 1, q = 0 gives exactly 1.
    assert golden["min_gap"] == pytest.approx(1.0)
    assert (golden["at_k"], golden["at_q"]) == (1, 0)
    assert kt.continuum_gap(float(np.pi - 3), 100_000)["min_gap"] < 0.1


@pytest.mark.parametrize("m", sorted(kt.TABLE_P3))
def test_orientable_closed_time_selects_the_same_harmonics(lat, m):
    """R2. The resonance condition is temporal closure; the caps halve it."""
    c = kt.closure_control_report(lat, m)
    untwisted, twisted = kt.TABLE_P3[m]
    assert c["orientable_periodic_modes"] == untwisted
    assert c["orientable_antiperiodic_modes"] == twisted
    assert c["orientable_periodic_dim"] == 2 * c["klein_untwisted_dim"]
    assert c["orientable_antiperiodic_dim"] == 2 * c["klein_twisted_dim"]


def test_time_symmetric_response_does_not_need_the_attachment(lat):
    """R2. Unique, reciprocal, time-symmetric with no deck constraint."""
    r = kt.response_control_report(lat)
    assert not r["source_is_deck_invariant"]
    assert r["unrestricted_green_symmetry"] < 1e-9
    assert r["time_mirror_defect_without_deck"] < 1e-9
    assert r["conditions_unrestricted"] == 2 * r["conditions_on_invariant_sources"]
    assert r["deck_invariant_source_dim"] == 81


@pytest.mark.parametrize("m", [0, 1, 7, 24])
@pytest.mark.parametrize("eta", [1, -1])
def test_complex_structure_leaves_every_twist_sector(lat, m, eta):
    """R3. ``j`` maps no nonzero admissible field to an admissible field."""
    assert kt.sector_action_report(lat, m, eta)["dim_plus_meet_j_plus"] == 0


def test_verdict_separates_missing_failed_and_counterexample():
    """R4. A missed tolerance is not a refutation."""
    groups = {"phase_space": ["a"], "complex_structure": ["b"],
              "phase_closure": ["c"], "field_freedom": ["d"],
              "continuum_free_field": ["e"]}
    v = kt.verdict({"a": True, "b": True, "c": False, "e": True}, groups,
                   {"continuum_free_field": "explicit mode"})
    assert v["phase_space"] == kt.CONFIRMED["phase_space"]
    assert v["phase_closure"] == kt.UNRESOLVED_NUMERICALLY
    assert v["field_freedom"] == kt.INCOMPLETE
    assert v["continuum_free_field"] == kt.COUNTEREXAMPLE
    assert v["failed_checks"] == ["c"] and v["missing_checks"] == ["d"]
    assert "REFUTED" not in v.values()
    assert v["born_rule"] == "NOT_DERIVED"


def test_probe_passes_and_claims_nothing_quantum(report):
    assert report["checks_passed"]
    assert len(report["checks"]) == 28
    v = report["verdict"]
    for field, value in kt.CONFIRMED.items():
        assert v[field] == value
    for field, value in kt.NOT_CLAIMED.items():
        assert v[field] == value
    assert v["continuum_free_field"] == kt.COUNTEREXAMPLE
    assert v["field_freedom"] != kt.FROZEN_LABELS["field_freedom"]
    assert set(kt.VERDICT_FIELDS) <= set(v)


def test_continuum_gap_is_the_true_minimum_over_both_candidates():
    """The nearest integer to k alpha is not always the minimiser."""
    for alpha in ((5 ** 0.5 - 1) / 2, 2 ** 0.5 - 1, 0.3, 7 / 48):
        got = kt.continuum_gap(alpha, 400)["min_gap"]
        brute = min(abs(k * k - (q / alpha) ** 2)
                    for k in range(1, 401) for q in range(0, 400))
        assert got == pytest.approx(brute, abs=1e-9)
