"""Metric backreaction: does A + B do something rescaling A or B cannot?

> Scope: a LINEAR conformally coupled scalar on a FIXED Einstein static
> universe with PR #257's POINT throat, and the metric response computed to
> FIRST ORDER in the lowest (n = 3, homogeneous) TRANSVERSE-TRACELESS harmonic.
> This is backreaction as a linear response, not a solved coupled system, and
> the mouths are not the resolved ones of PR #261/#262.

THE QUESTION
────────────
PR #260 gated stationary action and backreaction on a growing mode. PR #261
removed it by resolving the mouth; PR #262 removed the balls and made the
answer a theorem. The gate is closed, so this round asks the first GR question
the roadmap actually named -- and deliberately not "does spacetime pinch off?":

    DOES A + B PRODUCE A METRIC RESPONSE THAT RESCALING A OR B ALONE CANNOT?

The structure makes this a linear-algebra question with a number attached. The
field equation is linear so phi_{A+B} = phi_A + phi_B; T is quadratic so
T[A+B] = T[A] + T[B] + dT with dT bilinear; linearized Einstein is linear in T
so beta[A+B] = beta[A] + beta[B] + beta[dT]. Rescaling A -> cA sends
beta[A] -> c^2 beta[A], so EVERYTHING REACHABLE BY RESCALING IS THE
TWO-PARAMETER FAMILY {c^2 beta_A + d^2 beta_B}, and the question is exactly
whether beta[dT] lies in it. That is a projection residual.

THE ANSWER IS YES -- most of the interference response is unreachable.

THE CHANNEL, AND WHY IT IS THE ONLY HONEST ONE
──────────────────────────────────────────────
The transverse-traceless sector, lowest harmonic: the shear of the universe.
The reason is not convenience. THE ESU IS HELD STATIC BY A FLUID THIS ARC NEVER
SPECIFIES. A perfect fluid carries no anisotropic stress -- in an orthonormal
frame its T_ab = diag(rho,p,p,p) whatever the anisotropy -- so it contributes
nothing to the traceless spatial part, and neither does Lambda. This is the
ONLY channel whose answer does not depend on what was never put in. The scalar
sector depends on the sound speed, and it is also where the Eddington mode
lives, so an A/B/A+B comparison there would measure the background.

TT perturbations are also gauge-invariant, so no gauge choice enters; the mode
is the softest tensor channel, so if anything responds this does; and the
response reduces to a five-component driven oscillator, which makes the
measurement a statement about the SOURCE rather than about a PDE solver.

WHAT IS CHECKED
───────────────
T2  THE RESPONSE, DERIVED RATHER THAN RECALLED. Cartan about the ESU in the
    homogeneous anisotropy gives dG^TT_ij = beta''_ij + (8/a^2) beta_ij, so
    w3 = 2 sqrt(2)/a and w3^2 > 0: THE TENSOR SECTOR IS STABLE, which is the
    gate PR #260 taught the arc to check before computing on a background. The
    connection is obtained by SOLVING the torsion-free condition, not from a
    remembered formula -- a first draft quoted one and produced a G_00
    containing a'', which is impossible. Three validations with known answers:
    round S^3 (Ric = 2/a^2, R = 6/a^2), ESU (G = diag(3,-1,-1,-1)/a^2, which
    independently reproduces two_wave's _EINSTEIN), closed FRW
    (G_00 = 3(adot^2+1)/a^2). Two more facts fall out and are checked: the
    first-order piece is AUTOMATICALLY TRACELESS, and dG_0i = 0 identically --
    the momentum constraint -- so the mode is genuinely TT.

T3  THE SPLIT THAT MAKES IT LINEAR ALGEBRA. T[A+B] = T[A] + T[B] + dT measured
    rather than assumed: the self term scales as c^2 and the cross term as c to
    machine precision, and the cross term is EXACTLY zero with one source off.

T4  *** THE CONTROL, AND NOTHING HERE MEANS ANYTHING WITHOUT IT. *** The source
    is an integral over all of S^3 of an integrand with 1/chi^4 singularities at
    EIGHT points -- two sources, two mouths, four antipodes. A FIRST ATTEMPT AT
    THIS ROUND REPORTED 0.982 UNREACHABLE AND IT WAS PURE QUADRATURE NOISE:
    independent rules for the same quantity correlated at -0.04. Nothing about
    the run looked wrong. Two fixes were needed -- the mouths in the singular
    set (two_source puts the first at (1,0,0,0), exactly the natural quadrature
    axis, so the grid's own pole sat on a singularity and refining made it
    WORSE), and a smooth partition of unity instead of excision (excision leaves
    the bulk integrand discontinuous and the rule stalls at 1e-03). The rule
    now is that two refinement levels must agree in correlation AND magnitude
    before any residual is quoted.

T5  *** THE ANSWER. *** Most of the interference response lies outside the span
    of the two single-wave responses, and it is comparable to them in size.
    Reported off the full linear SPAN, which strictly contains the physical cone
    {c^2 beta_A + d^2 beta_B}, so the figure is CONSERVATIVE -- the true
    unreachable fraction is at least this large. It MOVES with the time window
    (0.88 to 1.00 over the four reported), which is why T6 scans the carrier as
    well and the result is quoted as a range rather than a constant.

T6  *** THE RESONANCE TEST, DONE ON THE COUPLED SYSTEM -- WHICH REVERSES IT. ***
    An earlier version of this probe argued the tensor channel was off resonance
    BY CONSTRUCTION: the conformal scalar on the ESU has spectrum w_n = n+1,
    INTEGERS; T is quadratic and integers are closed under sums and differences,
    so a source built from the FREE field rings on integers; and w3 = 2 sqrt(2)
    is irrational, 0.172 from the nearest. All of that is true -- OF THE
    UNCOUPLED AMBIENT -- and the conclusion does not survive coupling. The
    throat rings where det(A - Gamma(w)) vanishes, and those zeros are NOT
    integers: at the working point 0.875, 1.854, 1.872, 2.706, 2.878, 3.698.
    The nearest sits 0.050 from w3, and near w3 they are spaced only ~0.17
    apart, so the mode cannot be far from one. Across SIXTEEN throat
    configurations the nearest is always within 0.086, and at d = 0.9 it is
    0.001 -- resonant to the width of the scan. So the corrected statement is a
    WORKING-POINT one and it points the other way: off resonance with the free
    ambient, generically NEAR-RESONANT with the coupled system. That is also the
    mechanism the earlier version lacked for the headline's carrier sensitivity.
    The species of error is worth naming: an argument established for one system
    carried over to a system that differs precisely in the thing being argued
    about. The same scan is this round's honesty check on its own headline: the
    unreachable fraction MOVES with the carrier and with the time window, so it
    is reported as a range and not as a constant.

T7  WHICH BRANCHES WERE THERE. PR #257 measured the same configuration giving
    an invariant of anything from 0 to 4 depending on the arrival branches
    present, and said any quantity built by integrating over the field has to
    name them. Backreaction is such a quantity, so the throat is switched off as
    a control and both answers are reported rather than one being called THE
    answer.

WHAT IS STILL PUT IN
────────────────────
The n = 3 harmonic only -- the homogeneous shear, not the full TT tower. A FIXED
ESU background, POINT sources, and PR #257's POINT throat rather than the
resolved mouths of PR #261/#262. The response is LINEAR: first-order response,
not a solved coupled system. And the sources are points, so their own near
fields are in the integrand -- finite here, because the traceless angular
average kills the 1/chi^4 exactly, but a resolved source is the obvious next
refinement and it is the same finite-support move PR #261/#262 made.

    python -m experiments.closure_ledger.backreaction_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.backreaction import (
    TENSOR_MODE_FREQUENCY, ShearQuadrature, WORKING_BACKREACTION,
    measure_the_answer_needs_the_branches,
    measure_the_coupled_spectrum_is_near_resonant,
    measure_the_interference_response_is_unreachable,
    measure_the_quadrature_converges,
    measure_the_stress_tensor_splits_bilinearly,
    measure_the_tensor_sector_is_stable_and_fluid_free,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("PR #260's gate is closed by #261 and #262, so ask the "
                     "roadmap's first GR question: does A + B produce a metric "
                     "response that rescaling A or B alone cannot reproduce? "
                     "Everything reachable by rescaling is the two-parameter "
                     "family {c^2 beta_A + d^2 beta_B}, so the question is a "
                     "projection residual."),
        "scope": ("a linear conformally coupled scalar on a FIXED ESU with PR "
                  "#257's POINT throat; the metric response to FIRST ORDER in "
                  "the lowest (n=3, homogeneous) transverse-traceless harmonic. "
                  "Linear response, not a solved coupled system."),
        "pass": True,
    }


def t2_the_response_is_derived() -> dict:
    r = measure_the_tensor_sector_is_stable_and_fluid_free()
    return {"name": "T2_the_response_is_derived", **r,
            "pass": bool(r["the_validations_pass"]
                         and r["the_first_order_piece_is_traceless"]
                         and r["the_momentum_constraint_holds"]
                         and r["the_frequency_is_eight_over_a_squared"]
                         and r["the_sector_is_stable"])}


def t3_the_split_is_bilinear() -> dict:
    r = measure_the_stress_tensor_splits_bilinearly()
    return {"name": "T3_the_split_is_bilinear", **r,
            "pass": bool(r["the_self_term_is_quadratic"]
                         and r["the_cross_term_is_linear"]
                         and r["it_vanishes_with_one_source_off"])}


def t4_the_quadrature_control() -> dict:
    r = measure_the_quadrature_converges()
    return {"name": "T4_the_quadrature_control", **r,
            "pass": bool(r["every_component_is_converged"]
                         and r["the_residual_is_stable"])}


def t5_the_interference_is_unreachable() -> dict:
    r = measure_the_interference_response_is_unreachable()
    return {"name": "T5_the_interference_is_unreachable", **r,
            "pass": bool(r["most_of_it_is_unreachable"]
                         and r["the_two_singles_are_independent"]
                         and r["spread_over_windows"] < 0.2)}


def t6_the_coupled_spectrum() -> dict:
    r = measure_the_coupled_spectrum_is_near_resonant()
    return {"name": "T6_the_coupled_spectrum", **r,
            "pass": bool(r["the_free_ambient_is_off_resonance"]
                         and r["the_coupled_spectrum_is_not_the_free_one"]
                         and r["the_mode_is_near_resonant_everywhere"]
                         and r["it_is_a_working_point_quantity"])}


def t7_which_branches() -> dict:
    r = measure_the_answer_needs_the_branches()
    return {"name": "T7_which_branches", **r,
            "pass": bool(r["unreachable_both_ways"])}


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_response_is_derived(), t3_the_split_is_bilinear(),
             t4_the_quadrature_control(), t5_the_interference_is_unreachable(),
             t6_the_coupled_spectrum(), t7_which_branches()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t5, t6, t7 = tests[1:7]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_INTERFERENCE_RESPONSE_IS_NOT_REACHABLE_BY_RESCALING"
        verdict = (
            "*** A + B DOES SOMETHING RESCALING A OR B CANNOT. *** PR #260 "
            "gated backreaction on a growing mode; #261 and #262 closed that "
            "gate, so this round asks the roadmap's first GR question -- and "
            "not 'does spacetime pinch off?', but whether the two-wave "
            "configuration has a metric response outside the family a single "
            "wave can produce. The structure makes it a linear-algebra "
            "question: the field equation is linear so the fields add, T is "
            "quadratic so the cross term is bilinear, linearized Einstein is "
            "linear so the responses add, and rescaling A -> cA sends "
            "beta_A -> c^2 beta_A. Everything reachable is therefore "
            "{c^2 beta_A + d^2 beta_B}, and the measurement is the residual of "
            "beta[dT] off it. MEASURED: "
            f"{t5.get('unreachable_fraction', 0):.3f} of the interference "
            "response is unreachable, with the cross term "
            f"{t5.get('cross_over_single', 0):.2f} times the size of a "
            "single-wave response and the two singles independent "
            f"(|cos| = {t5.get('cos_between_the_singles', 0):.2f}). The "
            "residual is quoted off the full linear SPAN, which strictly "
            "contains the physical cone, so it is conservative; and it moves by "
            f"only {t5.get('spread_over_windows', 0):.3f} across time windows. "
            "*** THE RESPONSE IS DERIVED, NOT RECALLED. *** Cartan about the "
            "ESU in the homogeneous anisotropy gives "
            f"dG^TT = beta'' + (8/a^2) beta, so w3 = {TENSOR_MODE_FREQUENCY:.4f} "
            "and w3^2 > 0 -- THE TENSOR SECTOR IS STABLE, which is the gate PR "
            "#260 taught the arc to check first. The connection comes from "
            "SOLVING the torsion-free condition rather than a remembered "
            "formula, after a first draft's formula produced a G_00 containing "
            "a''. Three validations with known answers pass, one of them "
            "independently reproducing two_wave's _EINSTEIN; the first-order "
            "piece is automatically traceless; and the momentum constraint "
            "dG_0i = 0 holds identically, so the mode is genuinely TT. THAT IS "
            "ALSO WHY THIS CHANNEL: the ESU is held static by a fluid this arc "
            "never specifies, a perfect fluid carries no anisotropic stress, "
            "and neither it nor Lambda touches the traceless spatial part -- so "
            "this is the only channel whose answer does not depend on what was "
            "never put in. *** AND THE CONTROL IS THE POINT. *** A first "
            "attempt reported 0.982 unreachable and it was PURE QUADRATURE "
            "NOISE -- independent rules for the same quantity correlated at "
            "-0.04, with nothing about the run looking wrong. The integrand has "
            "1/chi^4 singularities at EIGHT points, and the two that were "
            "missing were the MOUTHS; two_source puts the first at (1,0,0,0), "
            "exactly the natural quadrature axis, so the grid's own pole sat on "
            "a singularity and refining made the answer worse. With all eight "
            "handled by a smooth partition of unity (excision leaves the bulk "
            "integrand discontinuous and stalls at 1e-03), the refinement "
            f"control now gives worst correlation "
            f"{t4.get('worst_correlation', 0):.4f} and worst magnitude drift "
            f"{t4.get('worst_magnitude_drift', 0):.4f}, with the residual "
            f"moving {t4.get('residual_drift', 0):.4f} between levels. THE "
            "CARRIER IS A MEASUREMENT TOO: the channel's transfer function is "
            "1/(w3^2 - w^2). *** AND THE RESONANCE TEST, DONE PROPERLY, "
            "REVERSES ITSELF. *** An earlier version of this probe proved the "
            "tensor mode incommensurate with the FREE ambient's integer "
            "spectrum -- w_n = n+1, and w3 = 2 sqrt(2) irrational, "
            f"{t6.get('free_ambient_gap_to_an_integer', 0):.3f} from the "
            "nearest -- and concluded the channel was off resonance BY "
            "CONSTRUCTION. That is true of the UNCOUPLED ambient and false of "
            "the coupled system: the throat rings where det(A - Gamma(w)) "
            "vanishes, and those zeros are not integers. Across sixteen throat "
            "configurations the nearest coupled resonance to w3 is always "
            f"within {t6.get('farthest_coupled_resonance', 0):.3f}, and at best "
            f"{t6.get('closest_coupled_resonance', 0):.4f} -- resonant to the "
            "width of the scan. So the channel is off resonance with the free "
            "ambient and generically NEAR-RESONANT with the coupled one, which "
            "is also the mechanism the earlier version lacked for the "
            "headline's carrier sensitivity. *** THE NUMBER IS NOT A UNIVERSAL "
            "CONSTANT, AND IS NOT REPORTED AS ONE. *** Across carriers the "
            f"unreachable fraction runs "
            f"{t6.get('unreachable_range',[0,0])[0]:.3f} to "
            f"{t6.get('unreachable_range',[0,0])[1]:.3f}, and across time "
            f"windows {t5.get('spread_over_windows', 0):.3f} wide -- large "
            "everywhere, constant nowhere. AND THE BRANCHES ARE NAMED: PR "
            "#257's rule that an integrated quantity must state which arrivals "
            "were present is why the throat is switched off as a control -- the "
            f"conclusion survives it "
            f"({t7.get('no_throat', {}).get('unreachable', 0):.3f} unreachable "
            "without the throat), so this is a statement about two waves rather "
            "than about the throat. WHAT IS STILL PUT IN: the n = 3 harmonic "
            "only, a FIXED ESU, POINT sources and PR #257's POINT throat rather "
            "than the resolved mouths of #261/#262, and a LINEAR response "
            "rather than a solved coupled system.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "backreaction", "tests": tests,
            "verdict_class": verdict_class, "verdict": verdict,
            "working_setup": {
                "source_gap": WORKING_BACKREACTION.source_gap,
                "carrier": WORKING_BACKREACTION.carrier,
                "width": WORKING_BACKREACTION.width,
                "tensor_mode_frequency": TENSOR_MODE_FREQUENCY},
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
                for row in v[:30]:
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
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, np.ndarray):
        return np.array2string(v, precision=5)
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in list(v)[:12]) + "]"
    return str(v)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_backreaction_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
