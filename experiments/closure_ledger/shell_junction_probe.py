"""
Can a detached shell do the throat's exotic work for it?

> Framing: static thin shells in spherical symmetry, Israel junction
> conditions, G = 1.  No radiation, no ℓ ≥ 2 structure, no backreaction on the
> shells' equations of state.

THE HOPE BEING TESTED
─────────────────────
A wormhole throat needs negative surface energy.  If a *detached* closed shell,
glued into the bulk with the opposite orientation, could supply the
exotic-looking restoring stress while itself being ordinary matter, the exotic
requirement would have been eliminated rather than relocated.

THREE OBSERVABLES, KEPT APART
─────────────────────────────
A single number cannot answer this, because three different questions are being
asked at once and they do not agree:

1  does the shell itself require exotic matter?      σ = −S^τ_τ
2  does it support the throat?                       F = −∂V_shell/∂b
3  is the supported configuration stable?            V''(b₀), and the modes

WHAT IS CHECKED
───────────────
T2  THE MACHINERY REPRODUCES PUBLISHED SHELLS.  A flat-interior bubble in
    Schwarzschild carries ordinary matter whose rest mass is the bulk mass to
    1e-3; a Z2 throat reproduces Visser's σ = −√f/(2πGR) to 6.9e-18.  Run
    before anything new is asked of it.

T3  OBSERVABLE 1 IS A THEOREM.  For an anti-aligned shell the two extrinsic
    curvature terms ADD: σ = −(β₊ + β₋)/4πGR with both roots non-negative for
    any timelike shell.  So EVERY oppositely-glued shell is exotic, whatever
    the bulk.  200 000 random Schwarzschild / de Sitter / Reissner–Nordström
    pairs: zero counterexamples, worst σ = −6.7e-04.  The aligned control is
    positive 50.1% of the time, so the sweep can find a positive σ when one
    exists.  The same identity applies to the throat, because a throat IS a
    minimal surface — no bulk arrangement can relieve it.

T4  AN ALIGNED SHELL CAN BE ORDINARY.  Same regions, same radius, only the
    gluing differs: σ = +6.2e-05 aligned against −9.5e-02 anti-aligned.

T5  OBSERVABLE 2 IS POSITIVE.  The shell screens mass, shifting the throat's
    potential by 2GΔM/b and pushing outward with F = 2GΔM/b², matched to 1e-6.
    Real support, from non-exotic matter.

T6  OBSERVABLE 3 IS POSITIVE TOO, BUT NOT FOR FREE.  Both normal modes can be
    stable at once — diag(0.151, 0.022) — but only because the throat is given
    β² = dp/dσ < 0.  Screening RAISES the critical β² monotonically
    (−1.083 → −0.652 as the interior mass falls 1.0 → 0.5), so the shell does
    enlarge the window; it never reaches β² ≥ 0.

T7  AND THE THROAT IS STILL EXOTIC.  σ_throat = −0.027 in exactly the
    configuration where the shell is ordinary and supporting.  Three
    observables, three different signs, on one system.

T8  BIRKHOFF DECOUPLES THEM.  The shell's surface density varies 701× as it is
    moved from a = 8 to a = 200 at fixed screened mass — genuinely different
    shells — and the throat's σ is bit-for-bit identical, spread exactly 0.
    The off-diagonal Hessian entry vanishes, but that is STRUCTURAL: Birkhoff
    is assumed the moment the region between is written as Schwarzschild with a
    constant mass.  Reporting it as a measurement would be circular, and it is
    labelled as such.

T9  NO FLAT DIRECTION.  Under a dilation the spectrum scales as 1/L² exactly
    (drift 0.0) and the smallest scaled eigenvalue stays at 5.4e-04, so "both
    eigenvalues positive" is not an artefact of whatever fixed the scale.

SCOPE
─────
β² is a free parameter at the equilibrium, in the standard treatment — a fixed
global barotropic index is too rigid to admit any stable static throat at all
(V''(b₀) = 2GM(n−1)/b₀³ with n = 2 + 4w gives static only for w < −1/2 and
unstable there).  Both stiffnesses are verified against direct RK4 integration
of σ' = −(2/b)(σ + p).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.shells.junction import (
    measure_any_minimal_surface_is_exotic,
    measure_the_detached_shell_can_be_ordinary,
    measure_the_hessian_has_no_flat_direction,
    measure_the_junction_reproduces_known_shells,
    measure_the_shell_force_on_the_throat,
    measure_the_stability_window,
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
        "observables": [
            "1 shell's own Israel surface stress sigma",
            "2 shell-induced force on the throat, -dV_shell/db",
            "3 stiffness d2V/db2 and the coupled normal modes",
        ],
        "why_three": ("a positive stiffness alone would mean restoring and "
                      "would not establish that the shell supplied it"),
        "pass": True,
    }


def t2_the_machinery_reproduces_published_shells() -> dict:
    r = measure_the_junction_reproduces_known_shells()
    return {"name": "T2_the_machinery_reproduces_published_shells", **r,
            "pass": bool(r["visser_is_reproduced"]
                         and r["the_bubble_is_ordinary"]
                         and r["its_rest_mass_is_the_bulk_mass"])}


def t3_every_minimal_surface_is_exotic(samples: int = 200_000) -> dict:
    r = measure_any_minimal_surface_is_exotic(samples=samples)
    return {"name": "T3_every_minimal_surface_is_exotic", **r,
            "pass": bool(r["every_minimal_surface_is_exotic"]
                         and r["the_sweep_can_find_positive_sigma"])}


def t4_an_aligned_shell_can_be_ordinary() -> dict:
    r = measure_the_detached_shell_can_be_ordinary()
    return {"name": "T4_an_aligned_shell_can_be_ordinary", **r,
            "pass": bool(r["the_aligned_shell_is_ordinary"]
                         and r["the_anti_aligned_shell_is_exotic"])}


def t5_the_shell_supports_the_throat() -> dict:
    r = measure_the_shell_force_on_the_throat()
    return {"name": "T5_the_shell_supports_the_throat", **r,
            "pass": bool(r["the_force_opposes_closure"]
                         and r["it_grows_with_the_screened_mass"]
                         and r["it_matches_2_G_dM_over_b_squared"])}


def t6_the_stability_window_and_what_screening_does() -> dict:
    r = measure_the_stability_window()
    return {"name": "T6_the_stability_window_and_what_screening_does", **r,
            "pass": bool(r["screening_raises_the_critical_beta2"]
                         and r["the_window_always_needs_negative_beta2"])}


def t7_the_three_observables_disagree() -> dict:
    r = measure_the_three_observables_are_independent()
    return {"name": "T7_the_three_observables_disagree", **r,
            "pass": bool(r["they_do_not_agree"]
                         and r["the_throat_is_still_exotic"])}


def t8_birkhoff_decouples_them() -> dict:
    r = measure_the_throat_and_shell_are_decoupled()
    return {"name": "T8_birkhoff_decouples_them", **r,
            "pass": bool(r["the_shells_are_genuinely_different"]
                         and r["the_throat_never_notices"])}


def t9_no_flat_direction() -> dict:
    r = measure_the_hessian_has_no_flat_direction()
    return {"name": "T9_no_flat_direction", **r,
            "pass": bool(r["no_flat_direction"]
                         and r["eigenvalues_scale_as_inverse_length_squared"])}


def t10_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T10_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe(samples: int = 200_000) -> dict:
    tests = [t1_goal(), t2_the_machinery_reproduces_published_shells(),
             t3_every_minimal_surface_is_exotic(samples),
             t4_an_aligned_shell_can_be_ordinary(),
             t5_the_shell_supports_the_throat(),
             t6_the_stability_window_and_what_screening_does(),
             t7_the_three_observables_disagree(),
             t8_birkhoff_decouples_them(), t9_no_flat_direction()]
    tests.append(t10_assessment(tests))
    t3, t6, t7, t8 = tests[2], tests[5], tests[6], tests[7]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_EXOTIC_MATTER_MOVES_BUT_DOES_NOT_LEAVE"
        verdict = (
            "THE EXOTIC MATTER MOVES BUT DOES NOT LEAVE. Of the three "
            "observables the hope needed, two come back positive and the "
            "decisive one comes back negative. An ANTI-ALIGNED detached shell "
            "— the oppositely-glued one the proposal is about — is "
            "necessarily exotic: its two extrinsic-curvature terms add rather "
            "than cancel, so sigma = -(beta+ + beta-)/4piGR with both roots "
            f"non-negative. Over {t3.get('samples', 0)} random Schwarzschild, "
            "de Sitter and Reissner-Nordstrom pairs there is not one "
            "counterexample, and there cannot be, because it is an identity "
            "and not a statistic. The same identity applies to the throat, "
            "which IS a minimal surface, so no arrangement of bulk content "
            "can relieve it. What an ORDINARY aligned shell can do is real "
            "but different: it screens mass, pushing the throat outward with "
            f"F = 2G dM/b^2, and it raises the critical beta2 monotonically "
            "as it screens more, enlarging the throat's stability window. "
            "Both normal modes can then be positive at once. But the window "
            "never reaches beta2 >= 0, and in exactly the configuration where "
            "the shell is ordinary and supporting, the throat's own sigma is "
            f"{t7.get('throat_sigma', 0):.3f} — still exotic. Three "
            "observables, three different signs, on one system: which is why "
            "they had to be reported separately rather than collapsed into a "
            "single stability verdict. The last finding is the one that "
            "shapes what comes next. Birkhoff decouples the two surfaces "
            "exactly: the shell's surface density varies "
            f"{t8.get('shell_sigma_varies_by', 0):.0f}x as it is moved at "
            "fixed screened mass, and the throat's sigma does not change in "
            "its last bit. Spherical symmetry has no radiative channel — the "
            "same l = 0 fact wave_constraints found for the scalar — so a "
            "genuine two-mode trapped resonator cannot exist here at all, and "
            "the l >= 2 internal modes are not a later refinement but the "
            "only place such a coupling could live.")
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
