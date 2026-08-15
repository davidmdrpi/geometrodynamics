"""
Is the throat part of the field problem, or a shift applied to its free branches?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is still an IDENTIFICATION MAP with a coupling
> kappa put in by hand -- shells.junction (PR #249) priced that throat and the
> bill is inherited, unpaid. What is new is that the identification is now
> SOLVED FOR rather than applied after the fact. NOT DONE: no backreaction, no
> stress tensor, no topology change, no rate. The two-throat cross term measured
> in T7 is a THROAT-THROAT interference, NOT the two-source invariant of
> roadmap step 3.

THE GAP THIS CLOSES
───────────────────
PR #254 solved the field and got the strong result: on the ESU the conformally
coupled retarded Green function has EXACT image support, so PR #253's ray
branches are the field's branches, with the 1/(4 pi sin chi) shell law and a
Maslov sign the ray ledger could not carry.

But it did that with the throat on the outside. phi(M+,t) = eta phi(M-,t+Delta)
was applied TO THE FREE BRANCHES AFTER THEY WERE COMPUTED -- one traversal, by
construction, because a post-processing step cannot notice that what it re-emits
will come back.

Here the identification enters the equation that is solved:

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

T3  SOLVING THE THROAT RESUMS EVERY TRAVERSAL. 1/(1-L) against an explicit walk
    over 400 traversals: 3.5e-18. PR #254's answer is the n = 0 term, and its
    relative error is exactly the round-trip gain |L|.

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

T7  ONE THROAT IS RANK ONE; TWO THROATS ARE RANK TWO. The second throat adds a
    second outer product, and the interference between them is a full fringe
    (visibility -1.000 to +1.000) that is bilinear -- one factor from each
    throat -- so it is identically zero without either. Same shape as roadmap
    step 3's invariant; NOT that invariant, because these are throats, not
    sources.

T8  AND THE POST-PROCESSED ANSWER FAILS OUTRIGHT AT THE EIGENFREQUENCIES.
    |T_d| peaks where 1 - e^{-2pi u} = 2 pi gamma, so the critical coupling
    kappa_c = 1/max|T_d| falls LINEARLY in the regulator: measured exponent
    1.0000. As gamma is removed, every coupling is critical at some frequency,
    and the peak sits exactly on an ESU eigenfrequency. A one-traversal answer
    is not the leading term of a convergent expansion there.

THE ANSWER
──────────
The throat is now solved for, and it was not a rearrangement. The resolvent adds
histories -- arrivals at times the free-branch ledger does not contain -- the
pair index is the right primitive because closure is a condition on the pair
while the amplitude is a product over it, and the expansion PR #254 implicitly
truncated has an unbounded expansion parameter at the eigenfrequencies.
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
    measure_the_expansion_fails_at_the_eigenfrequencies,
    measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("is the throat identification part of the field problem, "
                     "or a shift applied to the free branches after they are "
                     "computed? and what is the primitive object -- is it "
                     "indexed by a pair of branches?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is still an identification map with a coupling "
                  "kappa put in by hand, and PR #249's exotic-matter bill is "
                  "inherited and unpaid. No backreaction, no stress tensor, no "
                  "topology change, no rate. The two-throat cross term in T7 is "
                  "a throat-throat interference, not roadmap step 3's "
                  "two-source invariant."),
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


def t7_one_throat_is_rank_one_and_two_are_rank_two() -> dict:
    r = measure_the_primitive_is_rank_one_for_one_throat_and_not_for_two()
    return {"name": "T7_one_throat_is_rank_one_and_two_are_rank_two", **r,
            "pass": bool(r["one_throat_is_rank_one"]
                         and r["two_throats_are_rank_two"]
                         and r["the_cross_term_is_a_full_fringe"])}


def t8_the_expansion_fails_at_the_eigenfrequencies() -> dict:
    r = measure_the_expansion_fails_at_the_eigenfrequencies()
    return {"name": "T8_the_expansion_fails_at_the_eigenfrequencies", **r,
            "pass": bool(r["kappa_critical_scales_like_damping"]
                         and r["the_peak_sits_on_an_esu_eigenfrequency"]
                         and r["resonance_is_where_post_processing_is_worst"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_closed_form_transfer_is_the_branch_sum(),
             t3_solving_the_throat_resums_every_traversal(),
             t4_the_solve_adds_arrivals_post_processing_cannot(),
             t5_closure_is_broadband_coherence(),
             t6_the_condition_does_not_factorize(),
             t7_one_throat_is_rank_one_and_two_are_rank_two(),
             t8_the_expansion_fails_at_the_eigenfrequencies()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_THROAT_IS_SOLVED_FOR_AND_THE_PRIMITIVE_IS_A_PAIR"
        verdict = (
            "YES -- AND IT WAS NOT A REARRANGEMENT OF PR #254. That round "
            "applied the identification phi(M+,t) = eta phi(M-,t+Delta) to the "
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
            "amplitude factorizes over that index (K is rank one) and the "
            f"CONDITION DOES NOT: {t6.get('n_closed')} pairs close where any "
            f"rule phrased on a alone and b alone would admit "
            f"{t6.get('pairs_a_single_index_rule_would_select')}. FIFTH, the "
            "rank counts histories -- one throat rank one, two throats rank "
            "two -- and the throat-throat cross term is a full fringe, "
            f"visibility {t7.get('cross_visibility_min', 0):.3f} to "
            f"{t7.get('cross_visibility_max', 0):.3f}, bilinear and therefore "
            "identically zero without either throat. That is the same shape as "
            "roadmap step 3's invariant and is explicitly NOT that invariant: "
            "these are throats, not sources. SIXTH, THERE IS A REGIME WHERE "
            "POST-PROCESSING IS NOT AN APPROXIMATION AT ALL: |T_d| peaks where "
            "1 - e^{-2pi u} = 2 pi gamma, so the critical coupling falls "
            f"linearly in the regulator -- measured exponent "
            f"{t8.get('mean_exponent', 0):.4f} -- and the peak sits exactly on "
            "an ESU eigenfrequency. As the regulator is removed every coupling "
            "is critical somewhere, and there the one-traversal answer is the "
            "first term of a divergent series. WHAT IS STILL PUT IN: the "
            "throat remains an identification map with kappa by hand, the "
            "background is fixed, the field is linear, and when Delta + ell_c "
            "< 0 the loop is closed in time and 1/(1 - L) is a "
            "self-consistency condition rather than a history sum -- it has a "
            "unique solution exactly when the branch-resolved loop gain is "
            "subcritical, which is a bound on kappa and not a derivation of "
            "it. No backreaction, no stress tensor, no topology change, no "
            "rate, and no two-source invariant.")
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
