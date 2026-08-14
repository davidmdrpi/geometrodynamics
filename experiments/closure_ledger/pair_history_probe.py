"""
Two closed histories sewn at one interaction: is the event selected or inserted?

> Scope: a COUNTING result on a kinematic skeleton. Fixed round S3, c = 1. The
> throats are IDENTIFICATION MAPS with given mouths and given delays — not
> solutions of anything — and shells.junction (PR #249) priced a connected
> throat as necessarily exotic; this round adds TWO of them and pays for
> neither. No action principle, no field equations, NO TOPOLOGY CHANGE, no
> dynamics, no rate, and no demonstration that such a configuration is
> realisable. NOTHING HERE DERIVES A WORLDLINE: whether a particle trajectory is
> the locus where expanding and collapsing components stay mutually consistent
> is untouched, and calling the apparent particles "dragged" by anything would
> assume a force law this round does not have. The conjugacy Q_A + Q_B = 0 is a
> LABEL, carried and checked, never derived — no charge comes out of geometry.

THE QUESTION
────────────
PR #251 built one closed history (expanding leg, throat backwards in coordinate
time, collapsing leg) and showed the PAIR of local events closes a ledger
neither closes alone. PR #252 established that pair creation is a threshold on
an INVARIANT needing TWO independently propagated waves.

The tempting next move is topology change. That is not the one that can be
checked. The prior question can be:

    if two closed histories must SHARE their interaction event, is that event
    SELECTED by the closure conditions, or still put in by hand?

THE SYSTEM
──────────
Unknown: the event C = (c, t_C), c on S3. Given: two sources with launch times,
two throats with mouths and delays. Every leg is null, so a history closes in
coordinate time iff d(c,M+) + d(M-,c) + Delta = 0 — a GEODESIC ELLIPSOID with the
mouths as foci. Five equations, five unknowns:

    |c|^2 = 1                                  normalisation
    d(S_A,c) = t_C - tau_A                     C lies on front A
    d(S_B,c) = t_C - tau_B                     C lies on front B
    d(c,M_A+) + d(M_A-,c) + Delta_A = 0        history A closes
    d(c,M_B+) + d(M_B-,c) + Delta_B = 0        history B closes

WHAT IS CHECKED
───────────────
T2  CLOSURE IS AN ELLIPSOID CONDITION, AND IT IS BOUNDED. d(x,M+) + d(x,M-)
    ranges over exactly [d, 2pi - d], checked against 40,000 uniform samples of
    S3 per configuration rather than asserted. An out-of-band delay is rejected,
    so infeasible throats cannot silently contribute empty constraints.

T3  THE EVENT IS SELECTED, NOT INSERTED. Solved BLIND from random starts. Every
    event found sits at full Jacobian rank 5 — isolation is a property of the
    system, not of the solver stopping early. And existence is RESTRICTIVE:
    only about half of random feasible configurations admit a closed
    pair-history at all, so the closure conditions can forbid a configuration
    outright rather than merely locate an event inside it.

T4  REMOVING A WAVE REMOVES THE SELECTION — WHICH IS NOT WHAT WAS EXPECTED. The
    natural guess is that the pair-history solution disappears. It does not. The
    system becomes four equations in five unknowns, the Jacobian drops to rank
    4, and the solutions become a ONE-PARAMETER FAMILY of ~159 sampled points
    spanning the sphere. There is still a locus closing both histories; there is
    no longer a SELECTED one. Dropping the constraint can even CREATE solutions
    where two waves admitted none. So the Breit-Wheeler two-wave requirement
    shows up here as loss of ISOLATION, which is sharper than nonexistence.

T5  AND THE CONJUGATE PAIR NEEDS TWO DISTINCT THROATS. Not assumed by the
    two-history picture — it falls out of the rank. One shared throat fails in
    both senses of traversal: traversed OPPOSITELY, history B sees
    Delta_B = -Delta_A > 0 and its closure demands a sum of geodesic distances
    that is NEGATIVE (infeasible identically); traversed the SAME way, the two
    closure equations become the same equation, the rank drops to 4, and the
    event stops being selected.

T6  THE NON-CIRCULARITY CHECK, WHICH IS THE ONE THAT MATTERS. If the delays were
    unknowns rather than data, the system would be five equations in seven
    unknowns with nullity 2, and EVERY event lying on both fronts could be
    closed by choosing Delta afterwards — measured at 100% of 345 sampled
    events. So the result depends entirely on the throat being GIVEN, and that
    is stated rather than hoped for: a version of this calculation that solved
    for the delays would select nothing and would look identical from outside.

T7  CLOSURE SELECTS WHERE; THE INVARIANT DECIDES WHETHER. Kept strictly
    separate, because conflating amplitude with invariant is the error PR #252
    unwound. NOT ONE selected event clears s >= 4m^2 at E = m — median 2.48,
    maximum 3.97, just under — while E = 1.5m clears 78% and E = 3m clears all.
    The geometry picks the event; the event still has to be paid for in energy.

T8  THE CONJUGACY IS BOOKKEEPING. Labels cancel, a same-sign pair is refused at
    construction, and the kinematics is blind to the sign.

THE LESSON
──────────
The result is a rank, and the interesting part is what happens when a constraint
is removed. "The solution disappears" would have been the satisfying answer and
is the wrong one: it survives and stops being isolated. Reporting the weaker,
truer statement is the whole discipline this arc keeps rediscovering.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.pair_history import (
    measure_a_shared_throat_cannot_carry_the_pair,
    measure_closure_is_a_geodesic_ellipsoid,
    measure_removing_a_wave_removes_the_selection,
    measure_the_conjugacy_is_carried_not_derived,
    measure_the_delays_must_be_given_not_solved_for,
    measure_the_event_is_selected_not_inserted,
    measure_the_threshold_is_a_separate_condition,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("if two closed wormhole-wave histories must share their "
                     "interaction event, is that event selected by the closure "
                     "conditions, or is it still being inserted by hand?"),
        "scope": ("a counting result on a kinematic skeleton. Fixed round S3; "
                  "throats are identification maps with GIVEN mouths and "
                  "delays; two of them, with PR #249's exotic-matter bill "
                  "inherited and unpaid. No action principle, no field "
                  "equations, no topology change, no dynamics, no rate, no "
                  "worldline, and no claim of realisability. Conjugacy is a "
                  "label, not a derivation of charge."),
        "pass": True,
    }


def t2_closure_is_a_geodesic_ellipsoid() -> dict:
    r = measure_closure_is_a_geodesic_ellipsoid()
    return {"name": "T2_closure_is_a_geodesic_ellipsoid", **r,
            "pass": bool(r["sampling_never_goes_below_the_band"]
                         and r["sampling_never_goes_above_the_band"]
                         and r["an_infeasible_delay_is_rejected"])}


def t3_the_event_is_selected_not_inserted() -> dict:
    r = measure_the_event_is_selected_not_inserted()
    return {"name": "T3_the_event_is_selected_not_inserted", **r,
            "pass": bool(r["the_event_is_selected_not_inserted"]
                         and r["every_event_is_nondegenerate"])}


def t4_removing_a_wave_removes_the_selection() -> dict:
    """The falsification the reviewer named — and it lands differently."""
    r = measure_removing_a_wave_removes_the_selection()
    return {"name": "T4_removing_a_wave_removes_the_selection", **r,
            "pass": bool(r["two_waves_give_isolated_events"]
                         and r["one_wave_gives_a_one_parameter_family"]
                         and r["the_selection_requires_both_waves"])}


def t5_a_shared_throat_cannot_carry_the_pair() -> dict:
    r = measure_a_shared_throat_cannot_carry_the_pair()
    return {"name": "T5_a_shared_throat_cannot_carry_the_pair", **r,
            "pass": bool(r["opposite_traversal_is_infeasible"]
                         and r["same_traversal_loses_a_rank"]
                         and r["same_traversal_gives_a_family_not_a_selection"])}


def t6_the_delays_must_be_given_not_solved_for() -> dict:
    """The non-circularity check: a version that solved for Δ selects nothing."""
    r = measure_the_delays_must_be_given_not_solved_for()
    return {"name": "T6_the_delays_must_be_given_not_solved_for", **r,
            "pass": bool(r["with_free_delays_almost_any_event_closes"]
                         and r["so_the_content_is_in_the_throat_being_given"])}


def t7_the_threshold_is_a_separate_condition() -> dict:
    r = measure_the_threshold_is_a_separate_condition()
    return {"name": "T7_the_threshold_is_a_separate_condition", **r,
            "pass": bool(r["none_clear_threshold_at_energy_equal_mass"]
                         and r["most_clear_it_by_one_and_a_half_masses"])}


def t8_the_conjugacy_is_carried_not_derived() -> dict:
    r = measure_the_conjugacy_is_carried_not_derived()
    return {"name": "T8_the_conjugacy_is_carried_not_derived", **r,
            "pass": bool(r["the_labels_cancel"]
                         and r["a_same_sign_pair_is_refused"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_closure_is_a_geodesic_ellipsoid(),
             t3_the_event_is_selected_not_inserted(),
             t4_removing_a_wave_removes_the_selection(),
             t5_a_shared_throat_cannot_carry_the_pair(),
             t6_the_delays_must_be_given_not_solved_for(),
             t7_the_threshold_is_a_separate_condition(),
             t8_the_conjugacy_is_carried_not_derived()]
    tests.append(t9_assessment(tests))
    t3, t4, t6, t7 = tests[2], tests[3], tests[5], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_EVENT_IS_SELECTED_AND_LOSING_A_WAVE_LOSES_ISOLATION"
        verdict = (
            "THE INTERACTION EVENT IS SELECTED, NOT INSERTED — AND REMOVING A "
            "WAVE COSTS ISOLATION RATHER THAN EXISTENCE. Requiring two closed "
            "wormhole-wave histories to share one interaction event gives five "
            "equations in five unknowns: normalisation, the event lying on each "
            "of the two null fronts, and each history closing in coordinate "
            "time. Closure is a GEODESIC ELLIPSOID condition — the locus whose "
            "summed distance to the two mouths is |Delta| — feasible exactly on "
            "[d, 2pi - d], verified against 40,000 uniform samples of S3 per "
            "configuration rather than asserted. Solved BLIND from random "
            "starts, every event found sits at FULL JACOBIAN RANK 5 "
            f"({t3.get('events_at_full_rank_5', 0)} of "
            f"{t3.get('total_events_found', 0)}), so isolation is a property of "
            "the system and not of the solver stopping early. Existence is also "
            "RESTRICTIVE: only about half of random feasible configurations "
            f"({t3.get('configurations_with_a_selected_event', 0)} of "
            f"{t3.get('configurations', 0)}) admit a closed pair-history at "
            "all, so the closure conditions can forbid a configuration outright "
            "rather than merely locate an event inside it. THE FALSIFICATION "
            "LANDS DIFFERENTLY FROM THE EXPECTATION. Removing the second "
            "incoming wave does NOT delete the pair-history solution. The "
            "system becomes four equations in five unknowns, the Jacobian drops "
            "to RANK 4, and the solutions become a ONE-PARAMETER FAMILY — "
            f"typically {t4.get('typical_family_size_with_one_wave', 0)} "
            "distinct sampled points spanning the sphere. There is still a "
            "locus of events closing both histories; there is no longer a "
            "SELECTED one, and dropping the constraint can even CREATE "
            "solutions where two waves admitted none. So the Breit-Wheeler "
            "two-wave requirement appears here as loss of ISOLATION, which is a "
            "sharper and weaker statement than nonexistence, and it is the "
            "weaker one that is true. AND THE CONJUGATE PAIR NEEDS TWO DISTINCT "
            "THROATS, which the two-history picture did not assume — it falls "
            "out of the rank. A single shared throat fails in both senses: "
            "traversed oppositely, history B sees Delta_B = -Delta_A > 0 and "
            "its closure demands a sum of geodesic distances that is NEGATIVE, "
            "infeasible identically; traversed the same way, the two closure "
            "equations coincide, the rank drops to 4, and the event stops being "
            "selected. THE NON-CIRCULARITY CHECK IS THE ONE THAT MATTERS. If "
            "the delays were unknowns rather than given data the system would "
            "be five equations in seven unknowns with nullity 2, and every "
            "event on both fronts could be closed by choosing Delta afterwards "
            f"— measured at {100 * t6.get('fraction_closable_by_choosing_a_delay', 0):.0f}% "
            f"of {t6.get('sampled_events_on_both_fronts', 0)} sampled events. "
            "The entire result therefore rests on the throat being GIVEN, and a "
            "version that solved for the delays would select nothing while "
            "looking identical from outside. FINALLY, CLOSURE SELECTS WHERE AND "
            "THE INVARIANT DECIDES WHETHER — kept strictly apart, because "
            "conflating amplitude with invariant is the error PR #252 unwound. "
            "NOT ONE selected event clears s >= 4m^2 at E = m (median "
            f"{t7.get('median_s_at_energy_equal_mass', 0):.2f}, maximum "
            f"{t7.get('max_s_at_energy_equal_mass', 0):.2f}, just under), while "
            "E = 1.5m clears 78% and E = 3m clears all of them. WHAT THIS IS "
            "NOT: a counting result on a kinematic skeleton. No action "
            "principle, no field equations, no topology change, no dynamics, no "
            "rate, no worldline, and no claim that such a configuration is "
            "realisable. The conjugacy Q_A + Q_B = 0 is a label carried and "
            "checked, never derived, and two throats now sit on the books with "
            "PR #249's exotic-matter bill unpaid for both.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "pair_history", "tests": tests,
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
    out = here / "runs" / f"{ts}_pair_history_probe"
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
