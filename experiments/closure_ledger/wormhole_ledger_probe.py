"""
One conserved wave, seen in pieces: an S³ wormhole loop and its ledger

> Scope: kinematics and bookkeeping on a fixed round S³.  No Einstein equation,
> no backreaction, no throat construction, and NO CLAIM that such an
> identification is dynamically realisable.  The wormhole is an identification
> map — a pair of mouths and a time offset — and kappa and the delay are
> PARAMETERS.  shells.junction (PR #249) priced what such a throat costs: a
> minimal surface has sigma < 0 identically, so a connected throat needs exotic
> matter, and this round INHERITS that bill without paying it.  Amplitudes are
> scalar and the "momentum exchange" is a bookkeeping label on the flux, not a
> computed cross-section.  c = 1, R_S3 = 1.

THE SCENE
─────────
An emitter fires a shell.  While it expands, a collapsing shell sweeps past it
and lands on a receiver, which recoils.  Locally: two unrelated events.
Globally: one object.

    E --expand--> M_future --throat, dt < 0--> M_past --collapse--> R

THE RESULT
──────────
The staging is geometry and the self-consistency is arithmetic, and neither is
the interesting part.  The interesting part is that a closed timelike loop with
a LINEAR wave has exactly one amplitude, A = A_src/(1 - kappa), for EVERY
kappa != 1 — so nothing has to be tuned for consistency and no paradox is
available.  That is a result about linear equations, not about time travel, and
the probe demonstrates the difference rather than asserting it.

WHAT IS CHECKED
───────────────
T2  THE DESTINATION IS GEOMETRY, NOT STAGING.  A geodesic sphere at distance chi
    has area 4 pi sin^2(chi), so a shell refocuses exactly at the antipode.
    Checked against uniform sampling of S³ through the enclosed VOLUME — the
    area's integral — and scored against the estimator's own BINOMIAL STANDARD
    ERROR rather than a fixed percentage.  Two earlier drafts failed here by
    asserting a percentage, which at small chi measures the sample size.

T3  AND THE ARRIVING SHELL IS ACTUALLY COLLAPSING.  The same fact used twice:
    the receiver sits at the past mouth's antipode, which is the only place
    dA/dchi = 4 pi sin(2 chi) is negative all the way in.  Run against a
    receiver displaced to chi = 1.2, where the same wave is still EXPANDING
    when it lands.  A first draft placed the receiver generically and kept the
    word "collapse" anyway.

T4  THE DRAWN OBJECT IS THE SHELL, AND ITS DRAWN SIZE IS THE AREA.  Every point
    of the rendered family sits at geodesic distance chi to 6.7e-15; the
    orthographic shadow never leaves the unit ball; and the shadow's screen
    extent is proportional to sin(chi) with ONE constant to 3.6e-16 — that is
    sqrt(A/4pi), so the apparent size in the picture IS the area law.

    The projection had to be changed to get this.  A stereographic chart is
    unbounded at its own pole, and a shell launched from a point sweeps all of
    S³, so whatever pole is chosen the shell crosses it once: the radius grows
    as 2/epsilon, verified constant to 4e-6 across four decades, and never
    converges.  The first renderer projected from q3 = +1 — the emitter's own
    position — so the emitter was a division by zero and never got drawn.

T5  ONE SELF-CONSISTENT AMPLITUDE.  The closed form against brute iteration of
    2000-4000 round trips, which solves nothing and just keeps sending the wave
    around: 7e-13 even at kappa = 0.99, where the amplification is x100.

T6  AND ONLY THE RESONANCE OBSTRUCTS.  A sweep of the complex plane including
    |kappa| > 1, where the ITERATION diverges but the fixed point is still there
    and still unique.  Divergence of a summation method is not non-existence of
    a solution.

T7  WHAT WOULD ACTUALLY BREAK IT IS NONLINEARITY.  Replace the linear return
    with A = S + kappa A^2 and the SAME closed loop has two solutions below a
    threshold and NONE above it.  Reported because the uniqueness would
    otherwise read as a statement about wormholes when it is a statement about
    linear equations.

T8  THE LEDGER CLOSES — AND THAT IS THE ASSUMPTION, NOT THE RESULT.  Flux
    conservation through the throat is put in when the mouths are identified;
    the residual checks the arithmetic and nothing else.  Said plainly, because
    a number closing to 1e-16 invites the opposite reading.

T9  THE DELAY DECIDES THE STORY AND NO CONSERVED QUANTITY.  Whether the receiver
    is struck before the emitter fires flips between delay = -2 and -12; the
    amplitude spread across four decades of delay is 0.0.

T10 TWO LOCAL EVENTS, ONE CONSERVED WAVE.  Neither half conserves anything
    alone — the emitter loses flux, the receiver gains it — and the pair closes
    to 1.1e-16.  That balance is the only thing making them one object.

THE LESSON
──────────
Both of this round's corrections are representational: a receiver placed where
the caption said "collapse" but the geometry said "arrive", and a projection
that sent the very point the scene is about to infinity.  Neither was a
numerical error and neither would have been caught by refining anything.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.wormhole_ledger import (
    measure_nonlinearity_is_what_would_break_it,
    measure_only_the_resonance_obstructs,
    measure_the_delay_does_not_enter_the_ledger,
    measure_the_drawn_shell_is_the_geodesic_level_set,
    measure_the_ledger_closes,
    measure_the_loop_has_one_self_consistent_solution,
    measure_the_receiver_is_struck_by_a_collapsing_shell,
    measure_the_shell_focuses_at_the_antipode,
    measure_two_local_events_one_conserved_wave,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("previous rounds visualised wormhole networks whose "
                     "mouths supply antipodal wavefronts; can the emitted "
                     "shell, the passing collapsing shell, the receiver's "
                     "recoil and a time-delayed past mouth be shown to be ONE "
                     "conserved wave rather than several that happen to "
                     "balance?"),
        "scope": ("kinematics and bookkeeping on a fixed round S3. The "
                  "wormhole is an IDENTIFICATION MAP, not a solution; kappa "
                  "and the delay are parameters; shells.junction (PR #249) "
                  "priced the throat and this round inherits that bill "
                  "without paying it. No Einstein equation, no backreaction, "
                  "no cross-section."),
        "pass": True,
    }


def t2_the_destination_is_geometry() -> dict:
    r = measure_the_shell_focuses_at_the_antipode()
    return {"name": "T2_the_destination_is_geometry_not_staging", **r,
            "pass": bool(r["the_area_law_holds"]
                         and r["derivative_of_volume_is_the_area"]
                         and r["total_volume_is_two_pi_squared"]
                         and r["area_vanishes_at_both_poles"])}


def t3_the_arriving_shell_is_collapsing() -> dict:
    """The word "collapse" made into the sign of dA/dchi, against a control."""
    r = measure_the_receiver_is_struck_by_a_collapsing_shell()
    return {"name": "T3_the_arriving_shell_is_actually_collapsing", **r,
            "pass": bool(r["default_arrival_is_the_antipode"]
                         and r["the_shell_collapses_only_at_the_refocus"]
                         and r["the_focused_arrival_has_vanishing_area"]
                         and r["the_displaced_arrival_does_not"])}


def t4_the_drawn_object_is_the_shell() -> dict:
    """And its drawn size is sqrt(A/4pi) — the picture measures rather than
    illustrates.  Includes why the stereographic chart had to go."""
    r = measure_the_drawn_shell_is_the_geodesic_level_set()
    return {"name": "T4_the_drawn_object_is_the_shell_and_its_size_is_the_area",
            **r,
            "pass": bool(r["the_drawn_points_are_the_level_set"]
                         and r["the_shadow_never_leaves_the_unit_ball"]
                         and r["the_drawn_size_is_sqrt_of_the_area"]
                         and r["the_stereographic_radius_diverges_as_two_over"
                               "_epsilon"]
                         and r["it_does_not_converge_under_refinement"]
                         and r["while_the_shadow_stays_in_the_unit_ball"])}


def t5_one_self_consistent_amplitude() -> dict:
    r = measure_the_loop_has_one_self_consistent_solution()
    return {"name": "T5_one_self_consistent_amplitude", **r,
            "pass": bool(r["the_fixed_point_is_what_the_loop_does"])}


def t6_only_the_resonance_obstructs() -> dict:
    r = measure_only_the_resonance_obstructs()
    return {"name": "T6_only_the_resonance_obstructs", **r,
            "pass": bool(r["the_fixed_point_exists_everywhere_but_one_point"]
                         and r["the_resonance_is_refused"])}


def t7_nonlinearity_is_what_would_break_it() -> dict:
    """The control that stops T5 reading as a result about time travel."""
    r = measure_nonlinearity_is_what_would_break_it()
    return {"name": "T7_nonlinearity_is_what_would_break_it", **r,
            "pass": bool(r["two_solutions_below_threshold"]
                         and r["none_above_it"]
                         and r["so_uniqueness_is_a_linearity_result"]
                         and r["every_reported_root_solves_the_loop"])}


def t8_the_ledger_closes_and_that_is_the_assumption() -> dict:
    r = measure_the_ledger_closes()
    return {"name": "T8_the_ledger_closes_and_that_is_the_assumption", **r,
            "pass": bool(r["every_ledger_closes"])}


def t9_the_delay_decides_the_story_only() -> dict:
    r = measure_the_delay_does_not_enter_the_ledger()
    return {"name": "T9_the_delay_decides_the_story_and_nothing_conserved", **r,
            "pass": bool(r["the_delay_changes_no_conserved_quantity"]
                         and r["but_it_decides_the_ordering"])}


def t10_two_local_events_one_wave() -> dict:
    r = measure_two_local_events_one_conserved_wave()
    return {"name": "T10_two_local_events_one_conserved_wave", **r,
            "pass": bool(r["it_is_exactly_the_antipode"]
                         and r["the_pair_conserves"]
                         and not r["emitter_alone_conserves"]
                         and not r["receiver_alone_conserves"])}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_destination_is_geometry(),
             t3_the_arriving_shell_is_collapsing(),
             t4_the_drawn_object_is_the_shell(),
             t5_one_self_consistent_amplitude(),
             t6_only_the_resonance_obstructs(),
             t7_nonlinearity_is_what_would_break_it(),
             t8_the_ledger_closes_and_that_is_the_assumption(),
             t9_the_delay_decides_the_story_only(),
             t10_two_local_events_one_wave()]
    tests.append(t11_assessment(tests))
    t2, t3, t4, t5 = tests[1], tests[2], tests[3], tests[4]
    t7, t9, t10 = tests[6], tests[8], tests[9]

    if all(t["pass"] for t in tests):
        verdict_class = "ONE_WAVE_AND_LINEARITY_IS_WHY"
        verdict = (
            "ONE CONSERVED WAVE, AND LINEARITY IS WHY IT COSTS NOTHING. The "
            "scene's staging is geometry rather than choice: a geodesic sphere "
            "at distance chi on S3 has area 4 pi sin^2 chi, so an "
            "energy-conserving shell refocuses EXACTLY at the antipode. "
            "Checked against uniform sampling of S3 through the enclosed "
            "volume — the area's integral — and scored against the "
            "estimator's own binomial standard error rather than a fixed "
            f"percentage: worst z = {t2.get('worst_z_score', 0):.2f}. That "
            "fact is used TWICE. The future mouth sits at the emitter's "
            "antipode, and the receiver sits at the past mouth's, which is "
            "the only place the arriving shell is genuinely COLLAPSING "
            "(dA/dchi = 4 pi sin 2chi < 0 all the way in) rather than merely "
            "arriving; against a receiver displaced to chi = "
            f"{t3.get('displaced_arrival_chi', 0):.1f} the same wave is still "
            "expanding when it lands. THE CONTENT IS THE SELF-CONSISTENCY. A "
            "closed timelike loop carrying a LINEAR wave has exactly one "
            "amplitude, A = A_src/(1 - kappa), matching brute iteration of "
            "the loop — which solves nothing and just keeps sending the wave "
            f"around — to {t5.get('worst_abs_error', 0):.1e}, including at "
            "kappa = 0.99 where the amplification is x100. The fixed point "
            "exists and is unique for EVERY kappa != 1, the single resonance "
            "being the only obstruction in the whole complex plane; it exists "
            "even where |kappa| > 1 and the iteration diverges, because "
            "divergence of a summation method is not absence of a solution. "
            "Nothing is tuned for consistency and no paradox is available. "
            "THAT IS A STATEMENT ABOUT LINEAR EQUATIONS, NOT ABOUT TIME "
            "TRAVEL, and the probe shows the difference rather than asserting "
            "it: the same loop with a quadratic return has two solutions "
            f"below a source threshold of {t7.get('threshold_source', 0):.4f} "
            "and NONE above it. What each observer sees is two unrelated "
            "events, and neither conserves anything alone — the emitter loses "
            f"flux ({t10.get('flux_out_of_emitter', 0):.4f}), the receiver "
            f"gains it ({t10.get('flux_into_receiver', 0):.4f}) — while the "
            f"pair closes to {t10.get('the_pair_closes_to', 0):.1e}. That "
            "balance is the only thing making them one object rather than "
            "two. THE PICTURE MEASURES RATHER THAN ILLUSTRATES: every drawn "
            "point sits at geodesic distance chi from its centre to "
            f"{t4.get('worst_geodesic_distance_error', 0):.1e}, and the "
            "shadow's screen extent is proportional to sin chi with one "
            f"constant to {t4.get('extent_ratio_spread', 0):.1e} — that is "
            "sqrt(A/4pi), so the apparent size in the figure IS the area law. "
            "Getting there required changing the projection, which is this "
            "round's real correction. A stereographic chart is unbounded at "
            "its own pole and a shell launched from a point sweeps ALL of S3, "
            "so whatever pole is chosen the shell crosses it once: the radius "
            "grows as 2/epsilon across four decades and never converges. The "
            "first renderer projected from q3 = +1, which is the emitter's own "
            "position — the emitter was a division by zero and never got "
            "drawn. WHAT IS PUT IN, PLAINLY. The wormhole is an identification "
            "map, not a solution. Flux conservation through the throat is an "
            "ASSUMPTION made when the mouths are identified, so the ledger "
            "residual checks the arithmetic and is not evidence for the model. "
            "PR #249 priced this throat — a minimal surface has sigma < 0 "
            "identically, so a connected throat needs exotic matter — and this "
            "round inherits that bill without paying it. And the delay, which "
            "decides the entire story, changes NO conserved quantity: the "
            "amplitude spread across delays from -2 to -1e4 is "
            f"{t9.get('amplitude_spread', 0):.1e}, while whether the receiver "
            "is struck before the emitter fires flips inside that range.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = ("INCONCLUSIVE. Failed checks: " + ", ".join(failed))

    return {"probe": "wormhole_ledger", "tests": tests,
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
    out = here / "runs" / f"{ts}_wormhole_ledger_probe"
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
