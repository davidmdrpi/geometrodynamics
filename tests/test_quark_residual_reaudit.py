"""The quark residual sector under the corrected radial scalar operator.

Same three kinds of check as `test_scalar_operator_audit.py`. The *structural*
ones (which knobs the ladder consumes, what the derivation formulas are) are
asserted against closed forms and the module's own definitions. The *impact*
ones are asserted as directions and orderings, not as frozen digits — their job
is to characterise a change. And the *scoping* ones pin the statements that
must not silently drift back, including the README number this round corrects.
"""

import math

import numpy as np
import pytest

from geometrodynamics.constants import R_OUTER
from geometrodynamics.qcd import residual_audit as ra
from geometrodynamics.qcd.quark_spectrum import LOCKED_QUARK_PARAMS
from geometrodynamics.tangherlini.operator_audit import OPERATORS
from geometrodynamics.tangherlini.radial import (V_scalar_tangherlini,
                                                 V_tangherlini_legacy)


# ── structural: the derivations are what they claim to be ───────────────────

def test_the_locked_constants_match_the_frozen_v3_quark_lock():
    """The audit must read the lock, not a copy that can drift from it."""
    assert ra.LOCKED_PINHOLE == LOCKED_QUARK_PARAMS.pinhole
    assert ra.LOCKED_TRANSPORT == LOCKED_QUARK_PARAMS.transport
    assert ra.LOCKED_RESISTANCE == LOCKED_QUARK_PARAMS.resistance


def test_the_transport_object_is_the_off_diagonal_not_a_diagonal_expectation():
    """Both eigenvectors must appear, or it is not a transport element.

    Guarded because exactly this bug was found in PR #271's review: a
    `dot(v**2, dV)` had been standing in for `<u_l1|dV|u_l2>`.
    """
    _, per_pair = ra.derive_transport(V_scalar_tangherlini)
    assert set(per_pair) == set(ra.TRANSPORT_PAIRS)
    for (l1, l2) in ra.TRANSPORT_PAIRS:
        assert l1 != l2, "a same-ell pair is a diagonal expectation, not transport"

    # A genuine off-diagonal cannot be reproduced by either diagonal.
    profiles, grid = ra.mode_profiles(V_scalar_tangherlini)
    from geometrodynamics.tangherlini.radial import r_to_rstar
    rstar = np.array([r_to_rstar(float(v), 1.0) for v in grid])
    order = np.argsort(rstar)
    xs = rstar[order]
    delta = (np.asarray(V_scalar_tangherlini(grid, 3, 1.0), dtype=float)
             - np.asarray(V_scalar_tangherlini(grid, 1, 1.0), dtype=float))[order]
    for ell in (1, 3):
        u = profiles[ell]["funcs"][0]["u_half"][order]
        u = u / math.sqrt(float(np.trapezoid(u * u, xs)))
        diagonal = abs(float(np.trapezoid(u * u * delta, xs)))
        assert not math.isclose(diagonal, per_pair[(1, 3)], rel_tol=1e-6)


def test_the_cross_ell_operator_itself_is_exactly_invariant():
    """`V_l2 - V_l1` cannot move: the correction `3A^2/4r^2` carries no ell."""
    grid, _, _ = ra.chebyshev_grid(80)
    for (l1, l2) in ra.TRANSPORT_PAIRS:
        legacy = (np.asarray(V_tangherlini_legacy(grid, l2, 1.0), dtype=float)
                  - np.asarray(V_tangherlini_legacy(grid, l1, 1.0), dtype=float))
        correct = (np.asarray(V_scalar_tangherlini(grid, l2, 1.0), dtype=float)
                   - np.asarray(V_scalar_tangherlini(grid, l1, 1.0), dtype=float))
        assert np.max(np.abs(legacy - correct)) < 1e-12


def test_alpha_q_is_normalised_to_one_at_ell_equals_one():
    for V in OPERATORS.values():
        table = ra.derive_alpha_q_table(V)
        assert table[1] == pytest.approx(1.0, abs=1e-12)
        assert sorted(table) == [1, 2, 3, 4, 5]
        # the throat flux ratio must rise with ell, or the log is meaningless
        assert all(table[k] < table[k + 1] for k in range(1, 5))


def test_the_winning_pinhole_construction_still_wins_under_the_correction():
    """The anti-numerology guardrail: the ranking must survive, not just the winner.

    `docs/quark_axioms.md` rules out eight alternative constructions. If the
    correction had promoted one of them, the "clean structural reading" would be
    a coincidence of the buggy operator. It does not — and the gap to the
    runner-up widens.
    """
    from geometrodynamics.tangherlini.operator_audit import vmax_sum

    def candidates(V):
        grid, _, _ = ra.chebyshev_grid(80)
        vmax = {l: float(np.max(np.asarray(V(grid, l, 1.0), dtype=float)))
                for l in range(1, 7)}
        return {
            "sum_1_5": sum(vmax[l] for l in range(1, 6)),
            "sum_1_6": sum(vmax[l] for l in range(1, 7)),
            "odd_only": vmax[1] + vmax[3] + vmax[5],
            "twice_odd": 2 * (vmax[1] + vmax[3] + vmax[5]),
            "degeneracy_weighted": sum((2 * l + 1) * vmax[l] for l in range(1, 6)),
            "ell_weighted": sum(l * vmax[l] for l in range(1, 6)),
            "vmax_5_only": vmax[5],
            "vmax_5_times_pi": vmax[5] * math.pi,
            "vmax_5_times_e": vmax[5] * math.e,
        }

    gaps = {}
    for name, V in OPERATORS.items():
        cand = candidates(V)
        off = {k: abs(v - ra.LOCKED_PINHOLE) / ra.LOCKED_PINHOLE
               for k, v in cand.items()}
        assert min(off, key=off.get) == "sum_1_5", \
            f"{name}: a different construction won"
        assert cand["sum_1_5"] == pytest.approx(vmax_sum(V, 1, 5), rel=1e-12)
        runner_up = min((k for k in off if k != "sum_1_5"), key=off.get)
        assert runner_up == "vmax_5_times_e"
        gaps[name] = off[runner_up] - off["sum_1_5"]

    assert gaps["scalar_correct"] > gaps["legacy"], "the ranking gap must widen"


def test_the_tortoise_and_raw_grid_pinholes_are_different_quantities():
    """The README quotes the tortoise value; a raw-r linspace gives another.

    Kept explicit so the two are never quietly swapped: the difference is
    quadrature domain, not physics, and it is larger than the correction.
    """
    for V in OPERATORS.values():
        tortoise = ra.derive_pinhole(V)
        raw = ra.derive_pinhole_raw_grid(V)
        assert abs(tortoise - raw) > 0.15


# ── impact: the headline, asserted as a direction ───────────────────────────

def test_all_three_quark_residuals_move_toward_their_locked_values():
    """The headline. Opposite to the lepton sector, and the reason for the round."""
    result = ra.measure_the_three_residuals_under_both_operators()
    assert result["all_three_move_toward_the_lock"]
    for row in result["rows"]:
        assert abs(row["corrected_rel_pct"]) < abs(row["legacy_rel_pct"]), \
            f"{row['residual']} did not improve"
        assert abs(row["corrected_rel_pct"]) < 1.0


def test_the_pinhole_residual_changes_sign_under_the_correction():
    """Legacy undershoots the fitted 22.25, corrected overshoots it."""
    rows = ra.measure_the_three_residuals_under_both_operators()["rows"]
    pinhole = next(r for r in rows if r["residual"] == "pinhole")
    assert pinhole["legacy_rel_pct"] < 0.0 < pinhole["corrected_rel_pct"]


def test_the_composite_reverses_the_per_knob_ordering():
    """Each corrected knob is closer to its lock; the corrected triple is worse.

    The legacy triple's errors partially cancel. This is the round's main
    caution and must not be smoothed over.
    """
    lad = ra.measure_the_ladder_without_retuning()
    assert lad["corrected_pinhole_alone_is_better"]
    assert lad["corrected_composite_is_worse"]
    assert lad["legacy_composite_beats_its_own_pinhole"]


def test_neither_derived_composite_reaches_the_fitted_lock():
    """'Each knob is derived to 1%' was never 'the derived knobs fit the ladder'."""
    lad = ra.measure_the_ladder_without_retuning()
    assert lad["neither_composite_reaches_the_lock"]
    assert lad["locked_max_rel_err"] < 0.02
    for value in lad["composite_max_rel_err"].values():
        assert value > 0.03


def test_the_quark_ladder_is_far_less_stiff_than_the_lepton_chain():
    """+4.8 against -17.5: the whole reason the two sectors disagree."""
    two = ra.measure_the_two_sectors_side_by_side()
    assert two["elasticity_ratio"] > 3.0
    assert two["quark"]["elasticity"] > 0.0
    assert two["lepton"]["elasticity"] < 0.0
    assert abs(two["quark"]["elasticity"]) < abs(two["lepton"]["elasticity"])


def test_the_local_derivative_and_the_secant_are_reported_separately():
    """PR #271's review insisted on this distinction; it holds here too."""
    ela = ra.measure_the_pinhole_elasticity()
    local = ela["local_d_ln_m_d_ln_pinhole_at_the_lock"]
    secant = ela["secant_across_the_two_derived_pinholes"]
    assert set(local) == set(secant) == set(ra.FITTED_SPECIES)
    for s in ra.FITTED_SPECIES:
        assert local[s] != secant[s], "reported as one quantity, not two"
        assert math.copysign(1, local[s]) == math.copysign(1, secant[s])


def test_the_ladder_optimum_is_not_the_fitted_lock():
    """The knob was fitted jointly, so scanning it alone lands elsewhere."""
    ela = ra.measure_the_pinhole_elasticity()
    assert ela["ladder_error_at_what_it_wants"] < ela["ladder_error_at_the_v3_lock"]
    assert ela["pinhole_the_ladder_wants"] != ra.LOCKED_PINHOLE
    assert ela["correction_halves_the_miss"]


# ── the cross-sector bracket, and its caveat ────────────────────────────────

def test_only_the_corrected_operator_brackets_the_canonical_r_outer():
    brk = ra.measure_the_cross_sector_r_outer_bracket()
    assert brk["only_the_corrected_operator_brackets_1_26"]
    assert brk["correction_narrows_the_split"]
    corrected = brk["per_operator"]["scalar_correct"]
    assert corrected["quark_r_outer"] < R_OUTER < corrected["lepton_r_outer"]
    legacy = brk["per_operator"]["legacy"]
    assert R_OUTER < min(legacy["quark_r_outer"], legacy["lepton_r_outer"])


def test_the_bracket_position_is_measured_against_its_own_width():
    """Each operator's position must use *that* operator's bracket width.

    Guarded because the prose first quoted `1.07` bracket-widths for the legacy
    miss, which is the legacy distance divided by the *corrected* width. The
    correct figure is `0.81`.
    """
    per_op = ra.measure_the_cross_sector_r_outer_bracket()["per_operator"]
    for name, b in per_op.items():
        lo = min(b["lepton_r_outer"], b["quark_r_outer"])
        hi = max(b["lepton_r_outer"], b["quark_r_outer"])
        assert b["canonical_position_in_bracket"] == pytest.approx(
            (R_OUTER - lo) / (hi - lo), rel=1e-12)
    assert per_op["legacy"]["canonical_position_in_bracket"] == \
        pytest.approx(-0.81, abs=0.02)
    docs = _read("docs/quark_residual_reaudit.md")
    assert "1.07 bracket-widths" not in docs


def test_the_bracket_is_recorded_as_suggestive_and_not_as_a_derivation():
    """A 0.9% window round 1.26 is not a selection of 1.26."""
    brk = ra.measure_the_cross_sector_r_outer_bracket()
    assert "does not derive it" in brk["caveat"]
    entry = next(e for e in ra.measure_the_quark_ledger()["entries"]
                 if "bracketed by two independent sectors" in e["claim"])
    assert entry["verdict"] == "NEW, AND WEAK"


def test_the_single_sector_fixed_point_is_not_quietly_restored():
    """PR #271 reopened it; this round adds different evidence, not a repair."""
    entry = next(e for e in ra.measure_the_quark_ledger()["entries"]
                 if "single-sector fixed point" in e["claim"])
    assert entry["verdict"] == "INTERPRETATION CHANGED"
    assert "not restored" in entry["evidence"]


# ── scoping: what must not drift back ───────────────────────────────────────

def test_the_readme_nine_percent_muon_error_is_corrected_to_the_measured_value():
    """The number this round fixes.

    #271's README row applied a generic half-percent illustration to the actual
    -0.75% residual. The measured value is ~15%, and it must be the one quoted.
    """
    two = ra.measure_the_two_sectors_side_by_side()
    assert two["the_readme_9_percent_is_wrong"]
    assert two["lepton"]["measured_error_pct"] > 14.0
    assert two["corrected_to_pct"] == two["lepton"]["measured_error_pct"]
    # the linearisation is close to the measured value, and neither is 9%
    assert abs(two["lepton"]["linearised_error_pct"]
               - two["lepton"]["measured_error_pct"]) < 2.0
    assert two["lepton"]["linearised_error_pct"] > 12.0

    readme = _read("README.md")
    assert "makes that a **9 %** muon error" not in readme


def test_the_lepton_resistance_conclusion_survives_but_its_reason_does_not():
    sel = ra.measure_the_lepton_resistance_selector()
    # under legacy the REJECTED candidate fitted better
    assert sel["competitor_beat_the_closed_form_under_legacy"]
    # under the correction the closed form wins on proximity too
    assert sel["closed_form_wins_under_the_correction"]
    entry = next(e for e in ra.measure_the_quark_ledger()["entries"]
                 if "7pi/100" in e["claim"])
    assert entry["verdict"] == "INTERPRETATION CHANGED"


def test_n_still_drifts_so_beta_remains_the_last_fit_knob():
    """The v3 lock's own conclusion, unchanged by the correction."""
    nst = ra.measure_the_n_stability()
    assert nst["N_drifts_under_every_mode"]
    for mode in ("baseline", "legacy", "corrected"):
        assert nst["per_mode"][mode]["width"] > 20


def test_the_topological_quark_claims_are_marked_invariant_not_re_run():
    """Proximity is not dependence — the same exclusion PR #271 made."""
    entries = ra.measure_the_quark_ledger()["entries"]
    shell = next(e for e in entries if "shell-index axioms" in e["claim"])
    assert shell["verdict"] == "EXACTLY INVARIANT"
    assert "k_5" in shell["evidence"]


def _read(rel: str) -> str:
    from pathlib import Path
    return (Path(__file__).resolve().parents[1] / rel).read_text()


def test_the_status_docs_carry_the_quark_reaudit_outcome():
    """The README must not still call the quark residuals plainly 'Verified'.

    They are verified *and improved*, but the composite claim beside them is
    not, and the row has to say which is which.
    """
    readme = _read("README.md")
    assert "quark_residual_reaudit" in readme
    # the three per-knob rows must carry the corrected residuals
    for fragment in ("+0.36%", "+0.70%", "−0.02%"):
        assert fragment in readme, f"missing corrected residual {fragment}"


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import quark_residual_reaudit_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
