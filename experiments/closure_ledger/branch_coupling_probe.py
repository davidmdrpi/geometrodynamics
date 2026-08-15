"""
Is the mouth transfer part of the field problem, or applied to its free branches?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R), and -- stated plainly, because the resolvent being exact says
> nothing about this -- a SELF-CONSISTENT RANK-ONE MOUTH-TRANSFER MODEL, NOT a
> throat boundary operator and NOT a quotient of the manifold. It relates field
> VALUES through the free Green function: no normal-derivative (flux) matching,
> no reflected channel (the mouth scattering object is 1x1 where a
> flux-conserving two-mouth junction needs at least 2x2 unitary), and with
> kappa < 1 it is lossy by construction -- power out over power in is kappa^2,
> measured in T10. So it is not literally an identification. shells.junction
> (PR #249) is what would fix kappa and supply the missing channels.
> NOT DONE: no backreaction, no stress tensor, no topology change, no rate. The
> two-throat cross term in T7 is a THROAT-THROAT interference, NOT the two-source
> invariant of roadmap step 3.

THE GAP THIS CLOSES
───────────────────
PR #254 solved the field and got the strong result: on the ESU the conformally
coupled retarded Green function has EXACT image support, so PR #253's ray
branches are the field's branches, with the 1/(4 pi sin chi) shell law and a
Maslov sign the ray ledger could not carry.

But it did that with the mouth relation on the outside. phi(M+,t) =
eta phi(M-,t+Delta) was applied TO THE FREE BRANCHES AFTER THEY WERE COMPUTED --
one traversal, by construction, because a post-processing step cannot notice
that what it re-emits will come back.

Here the relation enters the equation that is solved:

    a(omega) = eta kappa e^{-i omega Delta} [ S(omega) + T_d(omega) a(omega) ]
    a(omega) = eta kappa e^{-i omega Delta} S / (1 - L),  L = eta kappa
                                                  e^{-i omega Delta} T_d

and the primitive object is indexed by a PAIR OF BRANCHES -- one per leg:

    K_ab = eta kappa . s_a A_1 e^{-u l_a} . e^{-i omega Delta} . s_b A_2
                                                          e^{-u l_b}

WHAT IS CHECKED
───────────────
T2  THE BRANCH SERIES HAS A CLOSED FORM, AND ITS POLES ARE THE SPECTRUM. The
    short-way images all carry s = +1 and the long-way images all carry s = -1,
    so the winding sum is two geometric series:
    sum_b s_b e^{-u l_b} = (e^{-u chi} - e^{-u(2pi-chi)}) / (1 - e^{-2pi u}).
    Series against closed form: 2.7e-15. As the regulator is removed the
    denominator vanishes at omega = 1, 2, 3, ... -- the conformal ESU
    eigenfrequencies -- with residues equal to the mode functions over 2 omega.
    The image representation and the mode representation are one function.

T3  SOLVING THE RELATION RESUMS EVERY TRAVERSAL. 1/(1-L) against an explicit
    walk over 400 traversals: 3.5e-18. PR #254's answer is the n = 0 term, and
    its relative error is exactly the round-trip gain |L|.

T4  THE COUPLED FIELD HAS ARRIVALS THE FREE BRANCHES DO NOT. The solved waveform
    IS the sum over history words to 5.4e-06 -- so the words are not a story
    told about the waveform. At echo times no one-traversal word can reach, the
    solved field stands 10^12 above the control. Amplitudes follow the kappa^n
    ladder and every echo carries the product of every Maslov factor in its
    word. These are events at times PR #254's ledger does not contain.

T5  CLOSURE IS BROADBAND COHERENCE. K_ab carries the phase
    e^{-i omega (l_a + Delta + l_b)}, so PR #253's closure condition
    l_a + Delta + l_b = 0 is EXACTLY the statement that K_ab is independent of
    omega. Closed pairs have band coherence 1.000; every other pair dephases
    below 0.091.

T6  AND THE CONDITION DOES NOT FACTORIZE OVER THE INDEX THE AMPLITUDE
    FACTORIZES OVER. At Delta = -(chi_1 + chi_2 + 4pi) the closed set is
    {(k,j) : k + j = 2}: three pairs out of the nine that any rule phrased on a
    alone and b alone would have to admit. Meanwhile K_ab is rank one. That is
    why the PAIR is the primitive.

T7  RANK COUNTS TRANSFER CHANNELS, NOT HISTORIES. One throat already carries
    144 distinct (a,b) histories and K is still RANK ONE: what the rank counts
    is independent separable transfer channels -- outer products -- and one
    value-feedback throat supplies exactly one. A second throat adds a second,
    in a shared topological branch-label basis (checked, since the two throats
    have different chi on every leg), and the interference between the two
    channels is a full fringe (visibility -1.000 to +1.000) that is bilinear and
    therefore identically zero without either. Same shape as roadmap step 3's
    invariant; NOT that invariant, because these are throats, not sources.

T8  THE ONE-TRAVERSAL ANSWER FAILS NEAR THE BARE RESONANCES. |T_d| peaks where
    1 - e^{-2pi u} = 2 pi gamma, so the SERIES RADIUS kappa_series = 1/max|T_d|
    falls linearly in the regulator: measured exponent 1.0000, peak exactly on a
    bare ESU resonance. That is a statement about the EXPANSION and only that --
    gamma is an Abel regulator and T_d carries the bare poles.

T9  AND THE SERIES RADIUS IS NOT THE STABILITY THRESHOLD. Three conditions were
    being conflated: existence (L != 1, so 1/(1-L) is fine for |L| > 1),
    convergence (|L| < 1, the radius of sum L^n, which does not even depend on
    the delay), and stability (Im omega > 0 for every root of D = 1 - L in
    complex omega). The coupling DISPLACES the bare poles omega = m + i gamma by
    delta_m = -eta kappa e^{-i m Delta} sin(m d)/(4 pi^2 sin d), matched to
    2.2e-04, whose imaginary part goes like sin(m d) sin(m Delta) and CHANGES
    SIGN WITH THE MODE. So stability is phase-sensitive and no bound on |L| can
    decide it: kappa_series = 0.762 for every delay, while kappa_stability is
    0.771 at Delta = 1 and 3.034 at Delta = pi, a factor 3.98. At kappa = 1.520
    in that gap the traversal series diverges to 1.3e+119 while the solve is
    finite and the least-damped pole is still at Im omega = +0.0145.

T10 AND WHAT THE MODEL LEAVES OUT, AS NUMBERS. Power out over power in is
    kappa^2 exactly; the mouth scattering object is 1x1 and unitary only at
    kappa = 1; there is no flux matching and no reflected channel. The resolvent
    is exact for the model as posed -- the scope statement is about WHICH MODEL,
    not about the solve.

THE ANSWER
──────────
The mouth relation is now solved for, and it was not a rearrangement. The
resolvent adds histories -- arrivals at times the free-branch ledger does not
contain -- and the pair index is the right primitive because closure is a
condition on the pair while the amplitude is a product over it. What this is NOT
is a throat boundary operator: that needs flux matching, a reflected channel and
a unitary two-mouth S-matrix, and is the next PR rather than this one.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.branch_coupling import (
    measure_closure_is_broadband_coherence,
    measure_solving_the_throat_resums_every_traversal,
    measure_the_closed_form_transfer_is_the_branch_sum,
    measure_the_coupled_field_has_arrivals_the_free_branches_do_not,
    measure_the_one_traversal_expansion_fails_near_the_bare_resonances,
    measure_the_rank_counts_transfer_channels_not_histories,
    measure_the_series_radius_is_not_the_stability_threshold,
    measure_what_the_transfer_model_leaves_out,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("is the mouth relation part of the field problem, or a "
                     "shift applied to the free branches after they are "
                     "computed? and what is the primitive object -- is it "
                     "indexed by a pair of branches?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe, "
                  "and a self-consistent RANK-ONE MOUTH-TRANSFER MODEL rather "
                  "than a throat boundary operator: no flux matching, no "
                  "reflected channel, a 1x1 mouth scattering object, and lossy "
                  "for kappa < 1, so not literally an identification. PR #249's "
                  "exotic-matter bill is inherited and unpaid. No backreaction, "
                  "no stress tensor, no topology change, no rate. The "
                  "two-throat cross term in T7 is a throat-throat interference, "
                  "not roadmap step 3's two-source invariant."),
        "pass": True,
    }


def t2_the_closed_form_transfer_is_the_branch_sum() -> dict:
    r = measure_the_closed_form_transfer_is_the_branch_sum()
    return {"name": "T2_the_closed_form_transfer_is_the_branch_sum", **r,
            "pass": bool(r["the_series_is_the_closed_form"]
                         and r["the_poles_are_the_esu_spectrum"]
                         and r["the_residues_are_the_mode_functions"])}


def t3_solving_the_throat_resums_every_traversal() -> dict:
    """The claim of this round."""
    r = measure_solving_the_throat_resums_every_traversal()
    return {"name": "T3_solving_the_throat_resums_every_traversal", **r,
            "pass": bool(r["the_resolvent_is_the_traversal_sum"])}


def t4_the_solve_adds_arrivals_post_processing_cannot() -> dict:
    """The sharpest difference the solve makes: new events, not new amplitudes."""
    r = measure_the_coupled_field_has_arrivals_the_free_branches_do_not()
    return {"name": "T4_the_solve_adds_arrivals_post_processing_cannot", **r,
            "pass": bool(r["the_waveform_is_the_sum_over_history_words"]
                         and r["every_isolated_echo_stands_above_the_control"]
                         and r["every_echo_carries_its_maslov_word_sign"])}


def t5_closure_is_broadband_coherence() -> dict:
    r = measure_closure_is_broadband_coherence()
    return {"name": "T5_closure_is_broadband_coherence", **r,
            "pass": bool(r["closed_pairs_are_broadband_coherent"]
                         and r["every_other_pair_dephases"])}


def t6_the_condition_does_not_factorize() -> dict:
    """Why the pair is the primitive even though the amplitude is a product."""
    r = measure_closure_is_broadband_coherence()
    return {"name": "T6_the_condition_does_not_factorize",
            "n_closed": r["n_closed"],
            "pairs_a_single_index_rule_would_select":
                r["pairs_a_single_index_rule_would_select"],
            "the_condition_does_not_factorize":
                r["the_condition_does_not_factorize"],
            "the_amplitude_does_factorize": r["the_amplitude_does_factorize"],
            "delay": r["delay"],
            "pass": bool(r["the_condition_does_not_factorize"]
                         and r["the_amplitude_does_factorize"])}


def t7_the_rank_counts_transfer_channels_not_histories() -> dict:
    """Rank is not a history count -- one throat, 144 histories, rank one."""
    r = measure_the_rank_counts_transfer_channels_not_histories()
    return {"name": "T7_the_rank_counts_transfer_channels_not_histories", **r,
            "pass": bool(r["one_throat_is_one_channel"]
                         and r["two_throats_are_two_channels"]
                         and r["both_matrices_in_the_common_label_basis"]
                         and r["the_cross_term_is_a_full_fringe"])}


def t8_the_one_traversal_expansion_fails_near_the_bare_resonances() -> dict:
    r = measure_the_one_traversal_expansion_fails_near_the_bare_resonances()
    return {"name": "T8_the_one_traversal_expansion_fails_near_resonances",
            **r,
            "pass": bool(r["the_series_radius_scales_like_the_regulator"]
                         and r["the_peak_sits_on_a_bare_esu_resonance"]
                         and r["resonance_is_where_post_processing_is_worst"])}


def t9_the_series_radius_is_not_the_stability_threshold() -> dict:
    """Existence, convergence and stability are three different conditions."""
    r = measure_the_series_radius_is_not_the_stability_threshold()
    d = r["a_coupling_between_them"]
    return {"name": "T9_the_series_radius_is_not_the_stability_threshold", **r,
            "pass": bool(r["kappa_stability_depends_on_the_delay"]
                         and r["the_two_thresholds_are_different_numbers"]
                         and r["every_pole_matches_its_first_order_"
                               "displacement"]
                         and d.get("the_series_diverges")
                         and d.get("the_solve_is_finite")
                         and d.get("and_it_is_still_stable"))}


def t10_what_the_transfer_model_leaves_out() -> dict:
    """Which model this is -- as numbers, not as a disclaimer."""
    r = measure_what_the_transfer_model_leaves_out()
    return {"name": "T10_what_the_transfer_model_leaves_out", **r,
            "pass": bool(r["the_power_ratio_is_kappa_squared"]
                         and r["lossy_below_unit_coupling"]
                         and r["scattering_object_shape"] == "1x1")}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_closed_form_transfer_is_the_branch_sum(),
             t3_solving_the_throat_resums_every_traversal(),
             t4_the_solve_adds_arrivals_post_processing_cannot(),
             t5_closure_is_broadband_coherence(),
             t6_the_condition_does_not_factorize(),
             t7_the_rank_counts_transfer_channels_not_histories(),
             t8_the_one_traversal_expansion_fails_near_the_bare_resonances(),
             t9_the_series_radius_is_not_the_stability_threshold(),
             t10_what_the_transfer_model_leaves_out()]
    tests.append(t11_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9 = tests[1:9]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_TRANSFER_IS_SOLVED_FOR_AND_THE_PRIMITIVE_IS_A_PAIR"
        verdict = (
            "YES -- AND IT WAS NOT A REARRANGEMENT OF PR #254. That round "
            "applied the mouth relation phi(M+,t) = eta phi(M-,t+Delta) to the "
            "free branches after computing them, which gives one traversal by "
            "construction: a post-processing step cannot notice that what it "
            "re-emits will come back. Written into the field problem instead, "
            "the amplitude re-emitted at M- is driven by everything reaching "
            "M+ including its own return, and the solution carries a resolvent "
            "1/(1 - L). FIRST, THE BRANCH SERIES SUMS IN CLOSED FORM: the "
            "short-way images all carry the Maslov factor +1 and the long-way "
            "images all carry -1, so the winding sum is two geometric series "
            "and equals (e^{-u chi} - e^{-u(2pi-chi)})/(1 - e^{-2pi u}), "
            f"verified against the term-by-term sum to "
            f"{t2.get('worst_abs_error', 0):.1e}. Its poles, as the regulator "
            "is removed, sit at omega = 1, 2, 3, ... -- the conformal ESU "
            "eigenfrequencies -- with residues equal to the mode functions "
            "over 2 omega. The image representation and the mode "
            "representation are the same function seen from two sides, which "
            "is the strongest statement so far that the branch labels are a "
            "REPRESENTATION rather than an approximation. SECOND, THE SOLVE IS "
            "THE SUM OVER EVERY TRAVERSAL: the resolvent agrees with an "
            "explicit walk over 400 traversals to "
            f"{t3.get('worst_abs_error', 0):.1e}, and PR #254's answer is the "
            "n = 0 term, whose relative error is exactly the round-trip gain "
            "|L|. THIRD -- AND THIS IS THE SHARPEST PART -- THE SOLVE ADDS "
            "EVENTS, NOT AMPLITUDES. The solved waveform is the sum over "
            f"history words to {t4.get('reconstruction_relative_error', 0):.1e}, "
            "so the word enumeration is checked rather than asserted; at echo "
            "times ell_a + Delta + n(ell_c + Delta) + ell_b that no "
            "one-traversal word can reach, the solved field stands "
            f"{t4.get('worst_echo_contrast', 0):.1e} above the control, with "
            "amplitudes on the kappa^n ladder and signs equal to the product "
            "of every Maslov factor in the word. Those are arrivals at times "
            "PR #254's ledger does not contain. FOURTH, THE PRIMITIVE IS "
            "INDEXED BY A PAIR OF BRANCHES, AND FOR A GOOD REASON: K_ab "
            "carries the phase e^{-i omega (l_a + Delta + l_b)}, so PR #253's "
            "closure condition is EXACTLY the statement that K_ab is "
            "independent of omega -- closed pairs have band coherence "
            f"{t5.get('worst_closed_coherence', 0):.3f} while every other pair "
            f"dephases below {t5.get('best_other_coherence', 0):.3f}. The "
            "amplitude factorizes over that index and the CONDITION DOES NOT: "
            f"{t6.get('n_closed')} pairs close where any rule phrased on a "
            f"alone and b alone would admit "
            f"{t6.get('pairs_a_single_index_rule_would_select')}. FIFTH, THE "
            "RANK COUNTS TRANSFER CHANNELS AND NOT HISTORIES -- one throat "
            f"already carries {t7.get('n_histories_one_throat')} distinct "
            "(a,b) histories while K is rank one, because what an outer "
            "product counts is separable channels, and one value-feedback "
            "throat supplies one. A second throat adds a second, in a shared "
            "topological branch-label basis rather than by leg length, and the "
            "interference between the two channels is a full fringe, "
            f"visibility {t7.get('cross_visibility_min', 0):.3f} to "
            f"{t7.get('cross_visibility_max', 0):.3f}, bilinear and therefore "
            "identically zero without either. That is the same shape as "
            "roadmap step 3's invariant and is explicitly NOT that invariant: "
            "these are throats, not sources. SIXTH, THE ONE-TRAVERSAL ANSWER "
            "FAILS NEAR THE BARE RESONANCES: |T_d| peaks where 1 - e^{-2pi u} "
            "= 2 pi gamma, so the SERIES RADIUS falls linearly in the "
            f"regulator -- measured exponent {t8.get('mean_exponent', 0):.4f} "
            "-- with the peak on a bare ESU resonance. SEVENTH, AND THIS IS "
            "THE CORRECTION THAT MATTERS MOST, THE SERIES RADIUS IS NOT THE "
            "STABILITY THRESHOLD. Existence (L != 1), convergence (|L| < 1) "
            "and stability (Im omega > 0 for every root of D = 1 - L in "
            "complex omega) are three different conditions. The coupling "
            "DISPLACES the bare poles omega = m + i gamma by delta_m = -eta "
            "kappa e^{-i m Delta} sin(m d)/(4 pi^2 sin d), matched to "
            f"{t9.get('worst_pole_vs_first_order', 0):.1e}, and its imaginary "
            "part goes like sin(m d) sin(m Delta), which CHANGES SIGN WITH THE "
            "MODE -- so stability is phase-sensitive and no bound on |L| can "
            "see it. kappa_series is the same number for every delay while "
            f"kappa_stability is not, differing by a factor "
            f"{t9.get('largest_ratio', 0):.2f} at Delta = pi where every "
            "first-order displacement is real; and at a coupling in that gap "
            "the traversal series diverges while the solve stays finite and "
            "the least-damped pole stays in the upper half plane. WHAT THIS IS "
            "NOT: a throat boundary operator. The model relates field VALUES "
            "through the free Green function -- no normal-derivative (flux) "
            "matching, no reflected channel (the mouth scattering object is "
            "1x1 where a flux-conserving two-mouth junction needs at least 2x2 "
            "unitary), and power out over power in is kappa^2 exactly, so for "
            "kappa < 1 it is lossy and cannot be a quotient of the manifold. "
            "The resolvent is exact for the model as posed; the scope "
            "statement is about which model. kappa is still by hand, the "
            "background is fixed, the field is linear, and when Delta + ell_c "
            "< 0 the loop is closed in time and 1/(1 - L) is a "
            "self-consistency condition rather than a history sum. No "
            "backreaction, no stress tensor, no topology change, no rate, and "
            "no two-source invariant. The flux-conserving boundary operator is "
            "the next step, not this one.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "branch_coupling", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "generated_utc": datetime.now(timezone.utc).isoformat()}


# ════════════════════════════════════════════════════════════════════════════
def render_markdown(s: dict) -> str:
    lines = [f"# probe — {s['probe']}", "", f"_{s['generated_utc']}_", ""]
    for t in s["tests"]:
        lines.append(f"## {t['name']} — {'PASS' if t['pass'] else 'FAIL'}")
        lines.append("")
        for k, v in t.items():
            if k in ("name", "pass"):
                continue
            if isinstance(v, list):
                lines.append(f"- **{k}**:")
                for row in v:
                    if isinstance(row, dict):
                        lines.append("    - " + ", ".join(
                            f"{a}={_fmt(b)}" for a, b in row.items()))
                    else:
                        lines.append(f"    - {_fmt(row)}")
            else:
                lines.append(f"- **{k}**: {_fmt(v)}")
        lines.append("")
    lines += [f"## verdict — {s['verdict_class']}", "", s["verdict"], ""]
    return "\n".join(lines)


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, complex):
        return f"{v.real:.6g}{v.imag:+.6g}j"
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    return str(v)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, complex):
        return {"re": o.real, "im": o.imag}
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_branch_coupling_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
