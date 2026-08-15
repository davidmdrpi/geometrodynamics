"""
Does a solved wave field with a throat identification reproduce the branch ledger?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R), with the throat still an IDENTIFICATION MAP rather than a
> solution -- shells.junction (PR #249) priced that throat and the bill is
> inherited, unpaid. What is new is that propagation is SOLVED rather than
> enumerated. NOT DONE: no backreaction, no stress tensor, no topology change,
> no rate, and NO TWO-SOURCE INVARIANT yet. The self-consistent sum over
> repeated throat traversals is PR #251's fixed point and is not redone here;
> only single-traversal arrivals are computed.

THE QUESTION
────────────
PR #253 built a RAY ledger -- short way, long way, winding -- and ended by
conceding that rank counting had reached its limit. Before building the
two-source invariant that is supposed to replace it:

    does a solved wave field, with the throat imposed as the same
    identification, reproduce the exact branch ledger -- its arrival times, its
    amplitudes, and phases the ray picture never had?

WHY THE ANSWER IS EXACT AND NOT ASYMPTOTIC
──────────────────────────────────────────
On S3 the scalar Laplacian has eigenvalues -n(n+2), and R = 6, so the
CONFORMALLY coupled massless field has omega^2 = n(n+2) + 1 = (n+1)^2, i.e.

    omega_n = n + 1,  INTEGER.

The retarded Green function is therefore exactly periodic and its closed form is
a sum of images:

    G(chi,t) = 1/(4 pi sin chi) [ sum_k delta(t - chi - 2pi k)
                                - sum_k delta(t + chi - 2pi k) ]

so the geometric-optics branches are the EXACT SUPPORT of the field, not the
leading term of a stationary-phase expansion. Nothing has to be argued
asymptotically.

WHAT IS CHECKED
───────────────
T2  THE SOLVE MATCHES THE IMAGE SUM. A truncated eigenmode sum (which never sees
    an image) against the closed-form image sum (which never sees a mode):
    8.3e-13. The solve IS the image sum.

T3  THE FIELD'S SUPPORT IS THE RAY LEDGER. Peaks read off the SOLVED field land
    on PR #253's branch times to 3.0e-04 -- half a grid cell, so grid-limited
    rather than a disagreement -- with every branch matched by a peak and no
    peak spurious.

T4  AND THE FIELD SUPPLIES PHASES THE RAY LEDGER COULD NOT. Every arrival
    carries a sign, and it is (-1)^m with m the number of focal crossings: the
    antipode at t = pi, the source point again at t = 2pi, and so on. That is
    the MASLOV INDEX. 12 of 12 signs agree. Path-length bookkeeping gives
    arrival times and has no way to produce a sign -- this is the first thing in
    the arc that the ray picture structurally could not have said.

T5  THE AMPLITUDE IS PR #251's SHELL LAW, DERIVED RATHER THAN IMPOSED. That
    round set A ~ 1/sin chi by conserving energy across a shell of area
    4 pi sin^2 chi. Here peak * sin(chi) is the same constant at every chi to
    7.0e-16, and equals 1/(4 pi sin chi) exactly.

T6  AND THE LEDGER BELONGS TO THE CONFORMAL FIELD SPECIFICALLY. The MINIMALLY
    coupled field has omega = sqrt(n(n+2)), irrational, so no images and no
    sharp branches: 63% of the peak amplitude sits BETWEEN the arrivals, against
    4.0e-08 for the conformal field. PR #253 never said which field its ledger
    belonged to, because rays cannot tell the two apart.

T7  THE THROAT REPRODUCES THE CLOSURE CONDITION. A through-throat contribution
    lands at l_1 + Delta + l_2; setting Delta to minus a branch-pair sum -- which
    is exactly PR #253's closure condition -- puts an arrival back on the
    emission event. 9 of 9. And the field adds the sign: eta times the two
    Maslov factors.

THE ANSWER
──────────
Yes, and more sharply than asked. The branches are exact support rather than
stationary-phase contributions; the arrival times and amplitudes reproduce PR
#251 and PR #253 without being told to; and the phases are new information the
ray ledger structurally could not carry.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.field_solve import (
    measure_minimal_coupling_has_no_branch_structure,
    measure_the_amplitude_is_the_shell_law,
    measure_the_field_support_is_the_branch_ledger,
    measure_the_phases_are_the_maslov_index,
    measure_the_spectral_solve_matches_the_image_sum,
    measure_the_throat_reproduces_the_closure_condition,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("does a solved wave field, with the throat imposed as the "
                     "same identification, reproduce the exact branch ledger "
                     "of PR #253 -- arrival times, amplitudes, and phases the "
                     "ray picture never had?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is still an identification map, not a solution, "
                  "with PR #249's exotic-matter bill inherited and unpaid. No "
                  "backreaction, no stress tensor, no topology change, no rate, "
                  "and no two-source invariant yet. Repeated throat traversals "
                  "are PR #251's fixed point and are not redone here."),
        "pass": True,
    }


def t2_the_solve_matches_the_image_sum() -> dict:
    r = measure_the_spectral_solve_matches_the_image_sum()
    return {"name": "T2_the_solve_matches_the_image_sum", **r,
            "pass": bool(r["the_two_constructions_agree"])}


def t3_the_field_support_is_the_branch_ledger() -> dict:
    r = measure_the_field_support_is_the_branch_ledger()
    return {"name": "T3_the_field_support_is_the_branch_ledger", **r,
            "pass": bool(r["so_the_field_reproduces_the_ray_ledger"])}


def t4_the_phases_are_the_maslov_index() -> dict:
    """The part the ray ledger structurally could not supply."""
    r = measure_the_phases_are_the_maslov_index()
    return {"name": "T4_the_phases_are_the_maslov_index", **r,
            "pass": bool(r["every_sign_is_the_maslov_factor"])}


def t5_the_amplitude_is_the_shell_law() -> dict:
    r = measure_the_amplitude_is_the_shell_law()
    return {"name": "T5_the_amplitude_is_the_shell_law", **r,
            "pass": bool(r["the_product_is_constant"]
                         and r["matches_one_over_four_pi_sin_chi"])}


def t6_the_ledger_belongs_to_the_conformal_field() -> dict:
    r = measure_minimal_coupling_has_no_branch_structure()
    return {"name": "T6_the_ledger_belongs_to_the_conformal_field", **r,
            "pass": bool(r["the_ledger_belongs_to_the_conformal_field"])}


def t7_the_throat_reproduces_the_closure_condition() -> dict:
    r = measure_the_throat_reproduces_the_closure_condition()
    return {"name": "T7_the_throat_reproduces_the_closure_condition", **r,
            "pass": bool(r["closure_puts_an_arrival_on_the_emission_event"]
                         and r["every_closing_sign_is_eta_times_the_maslov"
                               "_factors"])}


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_solve_matches_the_image_sum(),
             t3_the_field_support_is_the_branch_ledger(),
             t4_the_phases_are_the_maslov_index(),
             t5_the_amplitude_is_the_shell_law(),
             t6_the_ledger_belongs_to_the_conformal_field(),
             t7_the_throat_reproduces_the_closure_condition()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_FIELD_REPRODUCES_THE_LEDGER_AND_ADDS_ITS_PHASES"
        verdict = (
            "YES -- AND THE BRANCHES ARE EXACT SUPPORT, NOT STATIONARY-PHASE "
            "CONTRIBUTIONS. On the Einstein static universe the scalar "
            "Laplacian has eigenvalues -n(n+2) and R = 6, so the CONFORMALLY "
            "coupled massless field has omega^2 = n(n+2) + 1 = (n+1)^2: the "
            "frequencies are exactly the integers. The retarded Green function "
            "is therefore exactly periodic, and in closed form it is a sum of "
            "images, G = 1/(4 pi sin chi) [sum_k delta(t - chi - 2pi k) - sum_k "
            "delta(t + chi - 2pi k)]. The geometric-optics branches are the "
            "EXACT SUPPORT of the solved field; nothing has to be argued "
            "asymptotically. A truncated eigenmode sum -- which never sees an "
            "image -- agrees with the closed-form image sum -- which never sees "
            f"a mode -- to {t2.get('worst_abs_error', 0):.1e}. THE FIELD'S "
            "SUPPORT IS PR #253's RAY LEDGER: peaks read off the SOLVED field "
            f"land on the branch times to "
            f"{t3.get('worst_time_error', 0):.1e}, which is half a grid cell "
            "and therefore grid-limited rather than a disagreement, with every "
            "branch matched by a peak and no peak spurious. AND THE FIELD "
            "SUPPLIES WHAT THE RAY LEDGER STRUCTURALLY COULD NOT: every arrival "
            "carries a sign, and it is (-1)^m with m the number of focal "
            "crossings -- the antipode at t = pi, the source point again at "
            f"t = 2pi, and so on. That is the MASLOV INDEX, and "
            f"{len(t4.get('rows', []))} of {len(t4.get('rows', []))} signs "
            "agree. Path-length bookkeeping gives arrival times and has no way "
            "to produce a phase; this is the first quantity in the arc that the "
            "ray picture could not in principle have carried. THE AMPLITUDE IS "
            "PR #251's SHELL LAW, NOW DERIVED RATHER THAN IMPOSED: that round "
            "set A ~ 1/sin chi by conserving energy across a shell of area "
            "4 pi sin^2 chi, and here peak * sin(chi) is the same constant at "
            f"every chi to {t5.get('relative_spread', 0):.1e}, equal to "
            "1/(4 pi sin chi). AND THE LEDGER BELONGS TO THE CONFORMAL FIELD "
            "SPECIFICALLY, which PR #253 had no way to notice: the MINIMALLY "
            "coupled field has omega = sqrt(n(n+2)), irrational, so no images "
            f"and no sharp branches -- {100 * t6['minimal']['ratio']:.0f}% of "
            "the peak amplitude sits BETWEEN the arrivals, against "
            f"{t6['conformal']['ratio']:.1e} for the conformal field. Rays "
            "cannot tell the two apart. FINALLY THE THROAT REPRODUCES THE "
            "CLOSURE CONDITION: a through-throat contribution lands at "
            "l_1 + Delta + l_2, and setting Delta to minus a branch-pair sum -- "
            "exactly PR #253's closure condition -- puts an arrival back on the "
            "emission event, with the field adding the sign as eta times the "
            "two Maslov factors. WHAT IS STILL PUT IN: the throat remains an "
            "identification map rather than a solution, the background is fixed, "
            "the field is linear, repeated traversals are PR #251's fixed point "
            "and are not redone, and the two-source invariant that vanishes "
            "when a source is removed is NOT built here -- it is the next step, "
            "and the branch structure established above is why it will need a "
            "branch-resolved definition.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "field_solve", "tests": tests,
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
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    return str(v)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_field_solve_probe"
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
