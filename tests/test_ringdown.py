"""Settling PR #270's ringdown cross-validation.

Four kinds of check. The *exact* ones pin the `D = 5` identities that make the
answer checkable at all — the `−1/r` tortoise tail and the Bessel potential
tail — at the precision they actually hold. The *asymptotic* ones pin the
eikonal invariants, which are exact and are what every solver is judged
against. The *convergence* ones pin that the independent solver earns its
number. And the *verdict* ones guard the conclusion: Kerr–Schild confirmed,
tortoise damping excluded, and the two honest negatives kept as negatives.
"""

import math

import numpy as np
import pytest

from geometrodynamics.tangherlini import ringdown as rd


# ── exact: the D=5 identities ───────────────────────────────────────────────

def test_the_tortoise_correction_is_minus_one_over_r_not_a_log():
    """`r* − r → −1/r`. A log would grow; this decays, so no Coulomb phase."""
    for r in (50.0, 200.0, 1000.0, 5000.0):
        assert (rd.tortoise(r) - r) * r == pytest.approx(-1.0, abs=2e-4)
    # and the deviation must fall like 1/r^2
    devs = [abs((rd.tortoise(r) - r) * r + 1.0) for r in (50.0, 200.0, 1000.0)]
    assert devs == sorted(devs, reverse=True)
    assert devs[0] / devs[-1] > 100.0


def test_the_potential_tail_is_exactly_bessel():
    """`V → [(ℓ+1)² − ¼]/r²` — the identity PR #271 used to fix the operator."""
    for ell in (0, 1, 2, 3, 5):
        limit = (ell + 1) ** 2 - 0.25
        assert float(rd.potential(1000.0, ell)) * 1000.0 ** 2 == \
            pytest.approx(limit, rel=1e-5)


def test_the_tortoise_map_inverts():
    for r in (1.0001, 1.01, 1.5, 3.0, 50.0, 500.0):
        assert rd.radius_of_tortoise(rd.tortoise(r)) == pytest.approx(r, rel=1e-9)


def test_the_deep_horizon_inverse_matches_the_closed_form_where_it_is_resolvable():
    """`r − 1 → 2e^{2(r*−1)}`, checked where a double can still hold the offset."""
    for rs in (-6.0, -10.0, -14.0):
        r = rd.radius_of_tortoise(rs)
        assert r > 1.0
        assert r - 1.0 == pytest.approx(2.0 * math.exp(2.0 * (rs - 1.0)), rel=1e-5)


def test_the_deep_horizon_inverse_saturates_at_the_horizon():
    """Below `r* ≈ −18` the offset is under machine epsilon, so `1 + offset` *is*
    `1.0`. The asymptotic branch exists to avoid a failed bracket, not to buy
    precision that a double cannot hold — and the limit it lands on is the right
    one, since the potential vanishes at the horizon anyway."""
    for rs in (-25.1, -30.0, -60.0, -200.0):
        r = rd.radius_of_tortoise(rs)
        assert r == 1.0
        assert float(rd.potential(r, 1)) == 0.0
    # the two branches meet without a jump across the r* = −25 switch
    assert rd.radius_of_tortoise(-24.9) == rd.radius_of_tortoise(-25.1)


def test_the_asymptotics_measurement_agrees_with_the_direct_checks():
    result = rd.measure_the_background_asymptotics_are_exact()
    assert result["no_logarithmic_tail"]
    assert result["the_tail_is_exactly_bessel"]
    assert result["the_tail_is_exactly_minus_one_over_r"]


# ── the exact eikonal asymptote ─────────────────────────────────────────────

def test_the_photon_sphere_is_exactly_root_two():
    e = rd.eikonal_limit(0)
    assert e["r_photon"] ** 2 == pytest.approx(2.0, abs=1e-12)


def test_the_orbital_frequency_is_exactly_one_half():
    """`Ω_c = √f(r_ph)/r_ph = √(1/2)/√2 = 1/2` — exact at `n = 3`."""
    assert rd.eikonal_limit(0)["omega_c"] == pytest.approx(0.5, abs=1e-12)


def test_the_lyapunov_exponent_is_one_over_root_two():
    assert rd.eikonal_limit(0)["lyapunov"] == pytest.approx(
        1.0 / math.sqrt(2.0), abs=1e-6)


def test_the_eikonal_frequency_scales_with_ell_plus_one():
    for ell in (0, 1, 5, 20):
        w = rd.eikonal_limit(ell)["omega"]
        assert w.real == pytest.approx(0.5 * (ell + 1), abs=1e-12)
        assert w.imag == pytest.approx(-0.353553, abs=1e-5)


# ── the independent solver ──────────────────────────────────────────────────

def test_the_characteristic_solver_converges_in_step_size():
    result = rd.measure_the_characteristic_scheme_converges()
    assert result["all_converging"]
    for row in result["rows"]:
        deltas = row["successive_differences"]
        assert len(deltas) == 2
        assert deltas[1] < deltas[0], "successive differences must shrink"


def test_the_solver_needs_no_spatial_boundary_conditions():
    """The null diamond's domain of dependence is why this is immune to the
    excision question that broke the other two approaches."""
    note = rd.measure_the_characteristic_scheme_converges()[
        "no_boundary_conditions_are_applied"]
    assert "domain of dependence" in note
    assert "null diamond" in note


def test_the_fundamental_modes_track_the_eikonal_asymptote():
    """The pattern a correct solver must show, and does."""
    result = rd.measure_the_fundamental_modes()
    assert result["every_real_part_sits_above_the_eikonal"]
    assert result["every_damping_within_10_percent_of_the_asymptote"]
    rows = {r["ell"]: r for r in result["rows"] if r["omega"] is not None}
    # the real part must approach 0.5(l+1) as l grows
    gaps = [abs(rows[l]["omega"][0] - 0.5 * (l + 1)) for l in (1, 2, 3)
            if l in rows]
    assert gaps == sorted(gaps, reverse=True), "must converge to the eikonal"


def test_the_ell_one_fundamental_is_where_the_disagreement_was():
    w = rd.fundamental_mode(1, step=0.05)
    assert w.real == pytest.approx(1.0162, abs=1e-3)
    assert w.imag == pytest.approx(-0.3624, abs=1e-3)


# ── the verdict ─────────────────────────────────────────────────────────────

def test_kerr_schild_is_confirmed_and_the_tortoise_damping_excluded():
    """The deliverable. #270's own prime suspect was the wrong code."""
    v = rd.measure_the_cross_validation_verdict()
    assert v["kerr_schild_is_confirmed"]
    assert v["tortoise_damping_is_excluded"]
    assert v["gap_to_kerr_schild"]["imag_percent"] < 0.1
    assert v["gap_to_tortoise"]["imag_percent"] > 30.0


def test_the_reported_gap_reproduces_pr_270s_own_37_percent():
    """#270 measured the two codes as 37% apart. So does this round."""
    v = rd.measure_the_cross_validation_verdict()
    assert v["gap_to_tortoise"]["imag_percent"] == pytest.approx(37.0, abs=1.5)


def test_four_independent_lines_all_land_near_minus_point_three_six():
    v = rd.measure_the_cross_validation_verdict()
    lines = v["damping_lines_of_evidence"]
    agreeing = [val for name, val in lines.items() if "tortoise" not in name]
    assert len(agreeing) == 4
    assert max(agreeing) - min(agreeing) < 0.02
    assert all(-0.39 < val < -0.34 for val in agreeing)
    assert lines["tortoise (PR #270)"] > -0.30, "the outlier, on the other side"


def test_the_verdict_does_not_claim_an_autopsy_it_cannot_do():
    """Neither #270 code was landed, so 'which line' is not answerable."""
    v = rd.measure_the_cross_validation_verdict()
    assert "cannot" in v["what_this_round_cannot_do"] or \
        "no autopsy" in v["what_this_round_cannot_do"]
    assert "landed" in v["what_this_round_cannot_do"]


# ── the negatives stay negative ─────────────────────────────────────────────

def test_the_frequency_domain_failure_is_reported_not_hidden():
    n = rd.measure_what_did_not_work()["frequency_domain_shooting"]
    assert "NON-CONVERGENCE" in n["status"]
    assert n["so_pr_270s_diagnosis_stands"]
    assert len(n["roots_across_epsilon"]) >= 3
    assert "ill-conditioned" in n["why"]


def test_the_sixth_order_wkb_failure_is_reported_not_hidden():
    n = rd.measure_what_did_not_work()["sixth_order_wkb"]
    assert "UNUSABLE" in n["status"]
    assert "DIVERGE" in n["why"] or "diverge" in n["why"]


def test_first_order_wkb_accuracy_is_stated_not_assumed():
    """It is good on damping and poor on the real part at low ℓ. Both are said."""
    n = rd.measure_what_did_not_work()["first_order_wkb_accuracy"]
    assert "0.4%" in n["damping_at_ell_1"]
    assert "13%" in n["real_part_at_ell_1"]
    assert "not a discrepancy" in n["reading"].lower()
    # and the claim must be checkable
    w1 = rd.wkb_fundamental(1)
    exact = rd.fundamental_mode(1, step=0.05)
    assert abs(w1.imag - exact.imag) / abs(exact.imag) < 0.01
    assert abs(w1.real - exact.real) / exact.real > 0.10


def test_the_ledger_marks_ell_zero_as_less_certain():
    entries = {e["claim"]: e for e in rd.measure_the_ringdown_ledger()["entries"]}
    zero = next(e for k, e in entries.items() if "l = 0" in k)
    assert zero["verdict"] == "NO"
    quotable = next(e for k, e in entries.items() if "can now be quoted" in k)
    assert "l = 1, 2, 3" in quotable["verdict"]


def test_the_ledger_records_the_wrong_suspect():
    entries = {e["claim"]: e for e in rd.measure_the_ringdown_ledger()["entries"]}
    suspect = next(e for k, e in entries.items() if "prime suspect" in k)
    assert suspect["verdict"] == "WRONG SUSPECT"


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import ringdown_cross_validation_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
