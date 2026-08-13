"""
Can a detached shell do the throat's exotic work for it?

> Scope of every claim: Einstein gravity, Darmois–Israel thin shells, spherical
> symmetry, vacuum bulk, G = 1.  The dimension is a parameter — D = 4 is the
> regression case and D = 5 (Tangherlini) is the one this program cares about.
> Nothing here bounds thick shells, non-spherical configurations or modified
> gravity.

THE HOPE BEING TESTED
─────────────────────
A wormhole throat needs negative surface energy.  If a detached closed shell,
glued into the bulk with the opposite orientation, could supply the
exotic-looking restoring stress while itself being ordinary matter, the exotic
requirement would have been eliminated rather than relocated.

THREE OBSERVABLES, KEPT APART
─────────────────────────────
1  does the shell itself require exotic matter?      σ = −S^τ_τ
2  does it support the throat?                       the gradient of ΔV
3  is the configuration stable?                      the stiffness V''(b₀)

WHAT IS CHECKED
───────────────
T2  THE MACHINERY REPRODUCES PUBLISHED SHELLS.  A flat-interior bubble in
    Schwarzschild carries ordinary matter whose rest mass is the bulk mass to
    1e-3; a Z2 throat reproduces Visser's σ = −√f/(2πGR) to 6.9e-18.

T3  ε IS DERIVED FROM THE GLUING, NOT SET BY HAND.  Each side retains either
    the INNER (r ≤ R) or OUTER (r ≥ R) branch, and ε follows.  That gives FOUR
    gluings, and it corrects the earlier framing: η = ε₊ε₋ = −1 covers TWO of
    them whose forced signs are OPPOSITE — the minimal surface (σ < 0 always)
    and the maximal surface (σ > 0 always).  η alone decides nothing.

T4  THE FORCED SIGNS HOLD IN EVERY DIMENSION.  σ_minimal = −(D−2)(β₊+β₋)/8πGR
    and σ_maximal = +the same, with β± ≥ 0 for any timelike shell.  40 000
    random Tangherlini / de Sitter / charged pairs across D = 4, 5, 6: zero
    violations of either.  Identities, so the sweep checks the implementation.

T5  AND THE DICHOTOMY THAT FOLLOWS.  A detached surface CONNECTING to the
    throat's region through a minimal surface is necessarily exotic; one that
    is non-exotic by its gluing is a MAXIMAL surface, which caps off on both
    sides and shares no bulk with the throat — non-exotic precisely because it
    is disconnected, so it cannot support anything.

T6  OBSERVABLE 2, AS A GRADIENT NOT A FORCE.  An ordinary bubble screens mass
    and shifts the throat's potential; −∂ΔV/∂b = +0.024 at screened μ = 0.6,
    acceleration contribution half that.  Taken at fixed throat rest mass with
    no equation-of-state response, so it is NOT an equilibrium-consistent
    force and is not called one.

T7  OBSERVABLE 3.  Screening raises the critical β² monotonically
    (−1.083 → −0.652), so the shell does enlarge the window; it never reaches
    β² ≥ 0.  β² is an equation-of-state derivative parameter — no sound-speed
    reading is attached to its sign, which for exotic matter would not be
    meaningful.

T8  AND THE THREE DISAGREE.  σ_shell = +6.2e-05 (ordinary), gradient positive
    (supporting), σ_throat = −0.027 (still exotic).  One system, three signs.

T9  BIRKHOFF, AND WHAT IT IS WORTH.  The shell's surface density varies 701× as
    it moves at fixed screened mass, and the throat's σ is bit-for-bit
    identical.  The off-diagonal vanishes, but that is STRUCTURAL — Birkhoff is
    imported the moment the intervening region is written with a constant mass
    parameter.  It establishes no separation-dependent coupling IN THIS MODEL,
    not that every spherical trapped resonator is impossible.

T10 A UNITS CHECK ON THE STIFFNESSES.  Rescaling lengths and mass parameters
    together sends them as 1/L² (drift 0.0).  This is dimensional bookkeeping
    only — it does NOT show a fixed system has no dilation mode, which would
    need the kinetic term this module does not derive.  For the same reason the
    quantities are called stiffnesses, never normal modes: no kinetic metric
    has been derived, so nothing here is a generalised eigenproblem.

SCOPE
─────
A fixed global barotropic index admits no stable static throat at all
(V''(b₀) = 2GM(n−1)/b₀³ with n = 2 + 4w in D = 4: static needs w < −1/2,
stability w > −1/4), so β² is free at the equilibrium as usual.  Both
stiffnesses are verified against direct RK4 integration of σ' = −(D−2)(σ+p)/R,
agreeing to 10 digits.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.shells.junction import (
    measure_the_detached_shell_can_be_ordinary,
    measure_the_forced_signs_hold_in_any_dimension,
    measure_the_gluing_fixes_the_sign,
    measure_the_junction_reproduces_known_shells,
    measure_the_shell_potential_gradient,
    measure_the_stability_window,
    measure_the_stiffnesses_scale_dimensionally,
    measure_the_three_observables_are_independent,
    measure_the_throat_and_shell_are_decoupled,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("does a detached oppositely-glued bulk shell change the "
                     "throat's effective force and stability while itself "
                     "remaining non-exotic?"),
        "observables": ["1 the shell's own Israel surface stress sigma",
                        "2 the shell-induced gradient of the throat potential",
                        "3 the stiffness d2V/db2"],
        "scope": ("Einstein gravity, Darmois-Israel thin shells, spherical "
                  "symmetry, vacuum bulk; D parametric with D=5 the BAM case"),
        "pass": True,
    }


def t2_the_machinery_reproduces_published_shells() -> dict:
    r = measure_the_junction_reproduces_known_shells()
    return {"name": "T2_the_machinery_reproduces_published_shells", **r,
            "pass": bool(r["visser_is_reproduced"]
                         and r["the_bubble_is_ordinary"]
                         and r["its_rest_mass_is_the_bulk_mass"])}


def t3_the_gluing_fixes_the_sign() -> dict:
    """eps is derived from the retained branches, not set by hand."""
    r = measure_the_gluing_fixes_the_sign()
    return {"name": "T3_the_gluing_fixes_the_sign", **r,
            "pass": bool(r["every_forced_sign_is_realised"]
                         and r["eta_alone_does_not_decide"])}


def t4_the_forced_signs_hold_in_any_dimension() -> dict:
    r = measure_the_forced_signs_hold_in_any_dimension()
    return {"name": "T4_the_forced_signs_hold_in_any_dimension", **r,
            "pass": bool(r["no_violations_in_any_dimension"])}


def t5_the_connected_dichotomy() -> dict:
    r = measure_the_detached_shell_can_be_ordinary()
    return {"name": "T5_the_connected_dichotomy", **r,
            "pass": bool(r["the_minimal_surface_is_exotic"]
                         and r["the_maximal_surface_is_ordinary"])}


def t6_the_shell_potential_gradient() -> dict:
    r = measure_the_shell_potential_gradient()
    return {"name": "T6_the_shell_potential_gradient", **r,
            "pass": bool(r["the_gradient_opposes_closure"]
                         and r["it_grows_with_the_screened_mass"])}


def t7_the_stability_window() -> dict:
    r = measure_the_stability_window()
    return {"name": "T7_the_stability_window", **r,
            "pass": bool(r["screening_raises_the_critical_beta2"]
                         and r["the_window_always_needs_negative_beta2"])}


def t8_the_three_observables_disagree() -> dict:
    r = measure_the_three_observables_are_independent()
    return {"name": "T8_the_three_observables_disagree", **r,
            "pass": bool(r["they_do_not_agree"]
                         and r["the_throat_is_still_exotic"])}


def t9_birkhoff_and_what_it_is_worth() -> dict:
    r = measure_the_throat_and_shell_are_decoupled()
    return {"name": "T9_birkhoff_and_what_it_is_worth", **r,
            "pass": bool(r["the_shells_are_genuinely_different"]
                         and r["the_throat_never_notices"])}


def t10_a_units_check_on_the_stiffnesses() -> dict:
    r = measure_the_stiffnesses_scale_dimensionally()
    return {"name": "T10_a_units_check_on_the_stiffnesses", **r,
            "pass": bool(r["stiffnesses_scale_as_inverse_length_squared"])}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_machinery_reproduces_published_shells(),
             t3_the_gluing_fixes_the_sign(),
             t4_the_forced_signs_hold_in_any_dimension(),
             t5_the_connected_dichotomy(), t6_the_shell_potential_gradient(),
             t7_the_stability_window(), t8_the_three_observables_disagree(),
             t9_birkhoff_and_what_it_is_worth(),
             t10_a_units_check_on_the_stiffnesses()]
    tests.append(t11_assessment(tests))
    t3, t4, t8, t9 = tests[2], tests[3], tests[7], tests[8]

    if all(t["pass"] for t in tests):
        verdict_class = "CONNECTED_MEANS_EXOTIC"
        verdict = (
            "CONNECTED MEANS EXOTIC. Deriving the orientation from the gluing "
            "rather than setting it by hand sharpens the result and corrects "
            "the earlier framing. Each side retains the INNER or the OUTER "
            "radial branch, which fixes eps with no freedom left, and that "
            "gives FOUR gluings rather than two. The parity eta = eps+ eps- "
            f"covers two of them at once — {', '.join(t3.get('eta_minus_one_covers', []))} "
            "— whose forced signs are OPPOSITE, so eta alone decides nothing "
            "and 'the oppositely-glued shell is exotic' was too coarse a "
            "statement. What is actually forced is sharper: a MINIMAL surface "
            "has sigma = -(D-2)(beta+ + beta-)/8piGR < 0 and a MAXIMAL surface "
            "the same with the other sign, both identities, and neither is "
            "violated once in 40 000 random Tangherlini, de Sitter and charged "
            f"pairs across D = {t4.get('dims')}. The dichotomy that produces "
            "is the answer. A detached surface that CONNECTS to the throat's "
            "region does so through a minimal surface and is necessarily "
            "exotic. A detached surface that is non-exotic by its gluing is a "
            "maximal surface, which caps off on both sides and shares no bulk "
            "with the throat — it is non-exotic precisely because it is "
            "disconnected, and therefore cannot support anything. Within "
            "Einstein-Israel spherical thin shells the exotic matter is "
            "relocated, never removed. The three observables still disagree on "
            "one system, which is why they are reported apart: the ordinary "
            f"bubble has sigma = {t8.get('observable_1_shell_sigma', 0):.2e}, "
            "its screening does push the throat outward, and the throat's own "
            f"sigma is {t8.get('throat_sigma', 0):.3f} regardless. Two "
            "cautions are carried in the output rather than buried. The "
            "shell's effect is reported as a potential GRADIENT, not a force: "
            "it is taken at fixed throat rest mass with no equation-of-state "
            "response. And Birkhoff's vanishing off-diagonal is STRUCTURAL, "
            "imported the moment the intervening region is written with a "
            "constant mass parameter; what is measured instead is that a "
            f"family of shells differing {t9.get('shell_sigma_varies_by', 0):.0f}x "
            "in surface density leaves the throat bit-for-bit unchanged. That "
            "establishes no separation-dependent coupling in THIS model, not "
            "that every spherical trapped resonator is impossible.")
    else:
        verdict_class = "INCOMPLETE"
        verdict = ("At least one check did not hold; see the failing test "
                   "before quoting any number from this run.")

    return {
        "probe": "shell_junction",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tests": tests,
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    lines = [f"# {s['probe']} probe", "", f"_generated {s['generated_utc']}_",
             ""]
    for t in s["tests"]:
        lines.append(f"## {t['name']} — {'PASS' if t['pass'] else 'FAIL'}")
        for k, v in t.items():
            if k in ("name", "pass"):
                continue
            if isinstance(v, list) and v and isinstance(v[0], dict):
                lines.append(f"- **{k}**:")
                for row in v:
                    lines.append("    - " + ", ".join(
                        f"{a}={_fmt(b)}" for a, b in row.items()))
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
    out = here / "runs" / f"{ts}_shell_junction_probe"
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
