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

THE SYSTEM, AND ITS BRANCH SCOPE
────────────────────────────────
Unknown: the event C = (c, t_C), c on S3. Given: two sources with launch times,
two throats with mouths and delays. Every leg is null, so a history closes in
coordinate time iff l_1 + l_2 + Delta = 0. Five equations, five unknowns:

    |c|^2 = 1                                  normalisation
    d(S_A,c) = t_C - tau_A                     C lies on front A
    d(S_B,c) = t_C - tau_B                     C lies on front B
    l_1 + l_2 + Delta_A = 0                    history A closes
    l_1 + l_2 + Delta_B = 0                    history B closes

d is the PRINCIPAL geodesic distance. On the principal branch (short-way, no
winding) the closure locus is a geodesic ellipsoid, but a null leg on a closed
S3 may also take the long way (2pi - d) or wind (+2pi k), and a MIXED branch
fixes the DIFFERENCE of the two distances — a hyperboloid, not an ellipsoid. So
the branch is part of every statement below, and is stated first rather than
discovered later.

WHAT IS CHECKED
───────────────
T2  CLOSURE IS AN ELLIPSOID CONDITION ON THE PRINCIPAL BRANCH, AND IT IS
    BOUNDED. d(x,M+) + d(x,M-) ranges over exactly [d, 2pi - d], checked against
    40,000 uniform samples of S3 per configuration. An out-of-band delay is
    rejected, so infeasible throats cannot silently contribute empty
    constraints.

T3  AND THE BRANCH SCOPE IS LOAD-BEARING. Inside the principal delay band the
    principal branch is the ONLY feasible one — which means every other
    measurement here is principal-branch BY CONSTRUCTION OF ITS PRIOR, not by
    argument. Sampling the whole delay axis opens up to four branches, of both
    sum and DIFFERENCE type. What survives branching is that discreteness is a
    PER-BRANCH property: on any fixed branch the system is still 5x5 and its
    roots are still isolated at full rank. Branching multiplies the number of
    candidate events and changes the existence rate, not the local structure.

T4  THE ALLOWED EVENTS ARE DISCRETE AND LOCALLY ISOLATED. Solved BLIND from
    random starts; every root found sits at full Jacobian rank 5. Stated
    carefully: multi-start root-finding plus full rank shows each root FOUND is
    locally isolated — NOT that all roots were found, and NOT that the event is
    unique. "The event is selected" was an overstatement and has been withdrawn.

T5  REMOVING AN EQUATION IS A DIMENSIONALITY CONTROL, NOT PHYSICS. Dropping
    wave B removes one scalar equation from a square nondegenerate system, so
    rank 5 -> 4 and a one-parameter family is exactly the implicit-function
    result — for ANY deleted equation. The probe therefore also deletes a
    CLOSURE equation instead and gets the identical drop. This establishes
    nondegeneracy and is NOT evidence that pair creation needs two photons;
    that content lives in the invariant s, which needs two independent momenta.
    What survives as interesting is only the direction: the solutions do not
    vanish, they stop being isolated.

T6  IN THIS MODEL THE CONJUGATE PAIR NEEDS TWO DISTINCT THROATS. Traversed
    OPPOSITELY, history B's closure demands a NEGATIVE sum of leg lengths, which
    is infeasible on EVERY branch — the one conclusion here independent of the
    branch scope. Traversed the SAME way on the SAME branch, the two closure
    equations coincide and the rank drops to 4. That second half is scoped to
    the minimal single-pass model and is SCANNED rather than argued: at winding
    <= 1 every distinct branch pair either reduces to the identical constraint
    or is jointly inconsistent, and no counterexample was found. Not found is
    not impossible, and a different gluing is not excluded.

T7  THE NON-CIRCULARITY CHECK, WHICH IS THE ONE THAT MATTERS. With the delays as
    unknowns the system is five equations in seven, and the nullity is MEASURED
    on the actual 5x7 Jacobian (rank 5, nullity 2) rather than counted. Every
    event lying on both fronts can then be closed by choosing Delta afterwards —
    100% of 345 sampled events, with feasibility checked for BOTH throats. So
    the result depends entirely on the throat being GIVEN; a version that solved
    for the delays would select nothing and look identical from outside.

T8  CLOSURE CONSTRAINS WHERE; THE INVARIANT DECIDES WHETHER. Kept strictly
    separate. But the numbers carry two warnings: that no event clears
    s >= 4m^2 at E = m is FORCED, not measured — s = 2E^2(1-cos theta) <= 4E^2
    with equality only at exactly head-on, a measure-zero set — and every
    fraction reported is conditioned on this module's arbitrary prior over
    mouths, delays and launch times. They are regression diagnostics, NOT
    predictions.

T9  THE CONJUGACY IS BOOKKEEPING. Labels cancel, a same-sign pair is refused at
    construction, and the kinematics is blind to the sign.

THE SURVIVING CLAIM
───────────────────
With fixed mouth data and delays, and on a fixed propagation branch,
intersecting two null fronts with two independent closure hypersurfaces
generically produces locally isolated candidate events; removing one front
constraint restores a continuous degree of freedom. That is useful without
promoting the control into pair-production dynamics.

THE LESSON
──────────
Every correction this round needed was a SCOPE correction, not an arithmetic
one: which branch, which prior, which equation, and what a rank actually proves.
The numbers never changed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.pair_history import (
    measure_a_shared_throat_cannot_carry_the_pair,
    measure_the_results_are_scoped_to_the_principal_branch,
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


def t2b_the_results_are_scoped_to_the_principal_branch() -> dict:
    """Which branch everything else is about — stated first, not last."""
    r = measure_the_results_are_scoped_to_the_principal_branch()
    return {"name": "T3_the_results_are_scoped_to_the_principal_branch", **r,
            "pass": bool(
                r["inside_the_band_only_the_principal_branch_is_feasible"]
                and r["outside_it_more_branches_open"]
                and r["off_branch_loci_are_difference_type"]
                and r["discreteness_survives_on_a_fixed_off_branch"])}


def t3_the_event_is_selected_not_inserted() -> dict:
    r = measure_the_event_is_selected_not_inserted()
    return {"name": "T4_the_events_are_discrete_and_locally_isolated", **r,
            "pass": bool(r["events_are_discrete_and_locally_isolated"]
                         and r["every_event_is_nondegenerate"])}


def t4_removing_a_wave_removes_the_selection() -> dict:
    """The falsification the reviewer named — and it lands differently."""
    r = measure_removing_a_wave_removes_the_selection()
    return {"name": "T5_removing_an_equation_is_a_dimensionality_control", **r,
            "pass": bool(r["two_waves_give_isolated_events"]
                         and r["one_wave_gives_a_one_parameter_family"]
                         and r["deleting_a_closure_instead_drops_the_rank"
                               "_the_same_way"]
                         and r["the_square_system_behaves_nondegenerately"])}


def t5_a_shared_throat_cannot_carry_the_pair() -> dict:
    r = measure_a_shared_throat_cannot_carry_the_pair()
    return {"name": "T6_a_shared_throat_cannot_carry_the_pair", **r,
            "pass": bool(r["opposite_traversal_is_infeasible"]
                         and r["same_traversal_loses_a_rank"]
                         and r["no_branch_pair_rescues_a_shared_throat"]
                         and r["same_traversal_gives_a_family_not_a_selection"])}


def t6_the_delays_must_be_given_not_solved_for() -> dict:
    """The non-circularity check: a version that solved for Δ selects nothing."""
    r = measure_the_delays_must_be_given_not_solved_for()
    return {"name": "T7_the_delays_must_be_given_not_solved_for", **r,
            "pass": bool(r["with_free_delays_almost_any_event_closes"]
                         and r["the_nullity_is_measured_not_counted"]
                         and r["so_the_content_is_in_the_throat_being_given"])}


def t7_the_threshold_is_a_separate_condition() -> dict:
    r = measure_the_threshold_is_a_separate_condition()
    return {"name": "T8_the_threshold_is_a_separate_condition", **r,
            "pass": bool(r["none_clear_threshold_at_energy_equal_mass"]
                         and r["zero_percent_at_E_equals_m_is_forced_not"
                               "_measured"])}


def t8_the_conjugacy_is_carried_not_derived() -> dict:
    r = measure_the_conjugacy_is_carried_not_derived()
    return {"name": "T9_the_conjugacy_is_carried_not_derived", **r,
            "pass": bool(r["the_labels_cancel"]
                         and r["a_same_sign_pair_is_refused"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T10_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_closure_is_a_geodesic_ellipsoid(),
             t2b_the_results_are_scoped_to_the_principal_branch(),
             t3_the_event_is_selected_not_inserted(),
             t4_removing_a_wave_removes_the_selection(),
             t5_a_shared_throat_cannot_carry_the_pair(),
             t6_the_delays_must_be_given_not_solved_for(),
             t7_the_threshold_is_a_separate_condition(),
             t8_the_conjugacy_is_carried_not_derived()]
    tests.append(t9_assessment(tests))
    t2b, t3, t4 = tests[2], tests[3], tests[4]
    t6, t7 = tests[6], tests[7]

    if all(t["pass"] for t in tests):
        verdict_class = "DISCRETE_ON_A_FIXED_BRANCH_WITH_GIVEN_THROAT_DATA"
        verdict = (
            "WITH FIXED THROAT DATA AND ON A FIXED PROPAGATION BRANCH, TWO "
            "CLOSED HISTORIES SHARING AN EVENT LEAVE THE ALLOWED EVENTS "
            "DISCRETE. Requiring two closed wormhole-wave histories to share "
            "one interaction event gives five equations in five unknowns: "
            "normalisation, the event lying on each of the two null fronts, and "
            "each history closing in coordinate time. Solved BLIND from random "
            "starts, every root found sits at FULL JACOBIAN RANK 5 "
            f"({t3.get('events_at_full_rank_5', 0)} of "
            f"{t3.get('total_events_found', 0)}). Said precisely: multi-start "
            "root-finding plus full rank shows each root FOUND is locally "
            "isolated — not that all roots were found, and not that the event "
            "is unique. An earlier draft called this 'the event is selected'; "
            "that was an overstatement and is withdrawn. THE BRANCH SCOPE IS "
            "LOAD-BEARING AND IS STATED FIRST. d is the PRINCIPAL geodesic "
            "distance, and inside the principal delay band it is the only "
            "feasible branch — so every other measurement here is "
            "principal-branch BY CONSTRUCTION OF ITS PRIOR rather than by "
            "argument. Sampling the whole delay axis opens up to four branches, "
            "of both sum and DIFFERENCE type, and a mixed branch fixes the "
            "difference of the two distances: a hyperboloid, not an ellipsoid. "
            "What survives is that discreteness is a PER-BRANCH property — on "
            "any fixed branch the system is still 5x5 with isolated full-rank "
            "roots — so branching multiplies the candidate count and shifts the "
            "existence rate without touching the local structure. REMOVING A "
            "WAVE IS A DIMENSIONALITY CONTROL, NOT PHYSICS. Dropping wave B "
            "deletes one scalar equation from a square nondegenerate system, so "
            "rank 5 -> 4 and a one-parameter family is exactly what the "
            "implicit function theorem gives for ANY deleted equation. The "
            "probe therefore deletes a CLOSURE equation instead and gets the "
            "identical drop. This is evidence of nondegeneracy and NOT evidence "
            "that pair creation needs two photons; that content lives in the "
            "invariant s, which needs two independent momenta. What survives as "
            "interesting is only the direction of the surprise: the solutions "
            "do not vanish, they stop being isolated. IN THIS MODEL THE "
            "CONJUGATE PAIR NEEDS TWO DISTINCT THROATS. Traversed oppositely, "
            "history B's closure demands a NEGATIVE sum of leg lengths, "
            "infeasible on EVERY branch — the one conclusion here independent "
            "of branch scope. Traversed the same way on the same branch, the "
            "two closure equations coincide and the rank drops to 4; that half "
            "is scoped to the minimal single-pass model and SCANNED rather than "
            "argued, with every distinct branch pair at winding <= 1 either "
            "reducing to the identical constraint or jointly inconsistent, and "
            "no counterexample found. Not found is not impossible, and a "
            "different gluing is not excluded. THE NON-CIRCULARITY CHECK IS THE "
            "ONE THAT MATTERS: with the delays as unknowns the nullity is "
            f"MEASURED on the actual 5x7 Jacobian as "
            f"{t6.get('nullity_with_delays_free', 0)}, and "
            f"{100 * t6.get('fraction_closable_by_choosing_a_delay', 0):.0f}% "
            f"of {t6.get('sampled_events_on_both_fronts', 0)} sampled events on "
            "both fronts can then be closed by choosing Delta afterwards, with "
            "feasibility checked for BOTH throats. The entire result therefore "
            "rests on the throat being GIVEN. FINALLY, CLOSURE CONSTRAINS WHERE "
            "AND THE INVARIANT DECIDES WHETHER, and the numbers carry two "
            "warnings: that no event clears s >= 4m^2 at E = m is FORCED rather "
            "than measured, since s = 2E^2(1 - cos theta) <= 4E^2 with equality "
            "only at exactly head-on, a measure-zero set; and every fraction "
            "reported is conditioned on an arbitrary prior over mouths, delays "
            "and launch times, so they are regression diagnostics and NOT "
            "predictions. WHAT THIS IS NOT: a counting result on a kinematic "
            "skeleton. No action principle, no field equations, no topology "
            "change, no dynamics, no rate, no worldline, and no claim of "
            "realisability. Conjugacy is a label carried and checked, never "
            "derived, and two throats sit on the books with PR #249's "
            "exotic-matter bill unpaid for both.")
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
