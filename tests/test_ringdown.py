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

def test_the_tortoise_correction_is_exactly_minus_artanh_one_over_r():
    """The *exact* closed form. `−1/r` is only its leading term."""
    for r in (1.5, 3.0, 50.0, 200.0, 1000.0):
        assert rd.tortoise(r) - r == pytest.approx(-math.atanh(1.0 / r), rel=1e-12)


def test_the_leading_tail_is_minus_one_over_r_not_a_log():
    """`r* − r → −1/r`. A log would grow; this decays, so no Coulomb phase."""
    for r in (50.0, 200.0, 1000.0, 5000.0):
        assert (rd.tortoise(r) - r) * r == pytest.approx(-1.0, abs=2e-4)


def test_the_next_series_coefficient_is_one_third_as_predicted():
    """`r(r*−r) + 1 → −1/(3r²)`. Checking the coefficient, not just that it
    shrinks — a wrong tail could also shrink."""
    for r in (50.0, 200.0, 1000.0):
        deviation = (rd.tortoise(r) - r) * r + 1.0
        assert deviation == pytest.approx(-1.0 / (3.0 * r ** 2), rel=1e-3)


def test_the_potential_tail_is_asymptotically_bessel():
    """`V → [(ℓ+1)² − ¼]/r²` — the identity PR #271 used to fix the operator.
    Asymptotically: the exact `V` also carries `1/r⁴` and `1/r⁶` terms."""
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
    assert result["the_tail_is_asymptotically_bessel"]
    assert result["the_leading_tail_is_minus_one_over_r"]
    assert result["the_next_series_coefficient_is_confirmed"]
    # the prose must not overclaim an equality the series contradicts
    assert "NOT an exact equality" in result["the_exact_closed_form"]


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
    assert v["gap_to_kerr_schild"]["imag_percent_of_reference"] < 0.1
    assert v["gap_to_tortoise"]["imag_percent_of_reference"] > 20.0


def test_both_denominators_are_reported_and_neither_is_bare():
    """`27.1%` against truth, `37.3%` against the tortoise value. The module
    used to quote one in one place and one in another."""
    d = rd.measure_the_cross_validation_verdict()["the_denominator_is_named"]
    assert d["tortoise_relative_error_against_published"] == \
        pytest.approx(27.1, abs=0.5)
    assert d["correct_damping_is_larger_than_tortoise_by"] == \
        pytest.approx(37.3, abs=0.5)
    # the two conventions must not be confused for each other
    assert d["tortoise_relative_error_against_published"] < \
        d["correct_damping_is_larger_than_tortoise_by"]
    assert "denominator" in d["which_to_quote"]


def test_five_independent_lines_all_land_near_minus_point_three_six():
    v = rd.measure_the_cross_validation_verdict()
    lines = v["damping_lines_of_evidence"]
    agreeing = [val for name, val in lines.items() if "tortoise" not in name]
    assert len(agreeing) == 5, "published, characteristic, KS, WKB, eikonal"
    assert max(agreeing) - min(agreeing) < 0.02
    assert all(-0.39 < val < -0.34 for val in agreeing)
    assert lines["tortoise (PR #270)"] > -0.30, "the outlier, on the other side"


# ── the external oracle ─────────────────────────────────────────────────────

def test_the_published_reference_is_the_scaled_frequency_converted_correctly():
    """`ω = ω̃ · T_H`, `T_H = 1/(2π)`. The paper tabulates `ω̃`."""
    scaled = complex(6.38382253011, -2.27657411582)   # the paper's l = 1 entry
    converted = scaled / (2.0 * math.pi)
    published = rd.PUBLISHED_FUNDAMENTAL[1]
    assert converted.real == pytest.approx(published.real, rel=1e-11)
    assert converted.imag == pytest.approx(published.imag, rel=1e-11)


def test_the_characteristic_solver_reproduces_the_published_spectrum():
    """The strongest check in the round: an independent implementation landing
    on a spectrum computed by continued fractions and Hill determinants."""
    result = rd.measure_against_the_published_spectrum()
    assert result["every_mode_within_0p3_percent"]
    assert result["ell_1_and_2_within_0p05_percent"]
    assert "arXiv:2107.04815" in result["source"]


def test_the_published_values_confirm_the_verdict_independently():
    """Kerr–Schild lands on the reference; the tortoise damping does not."""
    published = rd.PUBLISHED_FUNDAMENTAL[1]
    ks_error = abs(rd.KERR_SCHILD_ELL_1.imag - published.imag) / abs(published.imag)
    tort_error = abs(rd.TORTOISE_ELL_1.imag - published.imag) / abs(published.imag)
    assert ks_error < 1e-3, "Kerr-Schild agrees with the published damping"
    assert tort_error > 0.2, "the tortoise damping does not"


def test_the_published_reference_marks_ell_zero_as_the_loosest():
    """Which validates having given `ℓ = 0` a visibly wider uncertainty."""
    assert rd.measure_against_the_published_spectrum()["ell_0_is_the_loosest"]


def test_step_refinement_understated_this_solvers_error():
    """The lesson only an external reference could supply: self-convergence
    measures the error it is refining away, not the total error."""
    r = rd.measure_against_the_published_spectrum()["refinement_versus_truth"]
    assert r["understatement_factor"] > 2.0
    assert r["the_finest_step_is_not_the_closest"], \
        "h = 0.05 lands closer to the published value than h = 0.025"


def test_the_round_does_not_claim_to_be_the_most_accurate_code():
    """#270's Kerr–Schild is ~6x better. Saying otherwise would misread it."""
    note = rd.measure_against_the_published_spectrum()["who_is_most_accurate"]
    assert "Kerr-Schild" in note
    assert "not the most" in note


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


# ── the discretization the arbitrator rests on ──────────────────────────────

def test_the_potential_is_sampled_at_the_diamond_centre():
    """For `N = (i,j)`, `S = (i−1,j−1)` the centre is `u_c = (i−½)h`,
    `v_c = (j−½)h`, so `r*_c = (j−i)h/2` — the half-steps cancel. An earlier
    version sampled `r*_c − h/4`, contradicting its own docstring."""
    import inspect
    source = inspect.getsource(rd.characteristic_evolution)
    assert "0.5 * step * offsets" in source
    assert "offsets - 0.5" not in source, "the h/4 shift must not come back"


def test_the_first_diamond_samples_the_potential_at_the_origin():
    """`i = j = 1` sits on `r* = 0`, not `−h/4`."""
    step = 0.05
    count = int(400.0 / step)
    offsets = np.arange(-count, count + 1)
    centres = 0.5 * step * offsets
    assert centres[count] == 0.0                       # the j − i = 0 entry
    assert centres[count + 1] == pytest.approx(0.5 * step)


def test_the_exact_potential_carries_more_than_the_inverse_square_term():
    """`V = L/r² + (9/4−L)/r⁴ − (9/4)/r⁶`, which is why the tail is only
    asymptotically Bessel."""
    for ell in (0, 1, 2, 3):
        limit = (ell + 1) ** 2 - 0.25
        for r in (5.0, 20.0, 100.0):
            expected = (limit / r ** 2 + (2.25 - limit) / r ** 4
                        - 2.25 / r ** 6)
            assert float(rd.potential(r, ell)) == pytest.approx(expected, rel=1e-12)


# ── where the error floor actually lives ────────────────────────────────────

def test_the_extraction_window_dominates_the_error_floor():
    """The claim that the residual is extraction systematics, measured."""
    result = rd.measure_the_extraction_systematics()
    assert result["the_window_dominates"]
    # varying the extraction moves the answer more than refining the step does
    assert result["extraction_band_exceeds_step_refinement"]
    assert result["how_many_times_larger"] > 2.0


def test_t_max_does_not_move_the_answer_at_all():
    """The extraction window sits well inside `t_max`, so it cannot."""
    assert rd.measure_the_extraction_systematics()["t_max_is_irrelevant"]


def test_the_round_no_longer_claims_only_an_external_check_could_find_it():
    """That claim was too strong — this internal scan finds it with no
    external value. The reference supplies the anchor, not the discovery."""
    note = rd.measure_the_extraction_systematics()["what_this_corrects"]
    assert "too strong" in note
    assert "anchor" in note


def test_no_extrapolated_limit_is_claimed():
    """The quoted number is a raw value at a stated step, not an `h → 0` limit."""
    v = rd.measure_the_cross_validation_verdict()
    assert v["this_round_step_size"] == 0.025
    note = v["this_round_is_a_raw_value_not_an_extrapolated_limit"]
    assert "not an h -> 0 limit" in note.lower().replace("NOT", "not")
    assert "Richardson" in note


def test_the_failed_shoot_does_not_claim_a_sole_cause():
    """Boundary truncation is a second confounder this round did not separate."""
    n = rd.measure_what_did_not_work()["frequency_domain_shooting"]
    assert "not a demonstrated sole cause" in n["a_second_confounder_not_separated"].lower()
    assert "series order" in n["what_would_separate_them"]
