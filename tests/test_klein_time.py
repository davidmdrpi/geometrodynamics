"""Klein time-attachment tests. Freeze ``bb6ff2c``."""

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


def test_verdict_fails_on_a_missing_check():
    groups = {"phase_space": ["a"], "complex_structure": ["b"],
              "phase_closure": ["c"], "field_freedom": ["d"]}
    v = kt.verdict({"a": True, "b": True, "c": False}, groups)
    assert v["phase_space"] == kt.CONFIRMED["phase_space"]
    assert v["phase_closure"] == "REFUTED"
    assert v["field_freedom"] == "REFUTED"
    assert v["born_rule"] == "NOT_DERIVED"


def test_probe_passes_and_claims_nothing_quantum(report):
    assert report["checks_passed"]
    assert len(report["checks"]) == 20
    v = report["verdict"]
    for field, value in kt.CONFIRMED.items():
        assert v[field] == value
    for field, value in kt.NOT_CLAIMED.items():
        assert v[field] == value
