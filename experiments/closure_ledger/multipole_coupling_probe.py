"""
ℓ ≥ 2: the coupling Birkhoff forbids at ℓ = 0, and what it costs to have it

> Scope: the STATIC NEWTONIAN (Laplace) multipole problem — the weak-field limit
> of the junction problem, with interior/exterior solutions r^ℓ and
> r^{−(ℓ+D−3)}.  Not a Regge–Wheeler/Zerilli treatment, and no claim about
> ℓ ≥ 2 quasinormal frequencies or radiative dynamics.  G = 1.

WHAT THIS CLOSES
────────────────
shell_junction found that two concentric surfaces in a vacuum spherical model
cannot talk, and imported Birkhoff to say why.  The natural next guess was that
ℓ ≥ 2 modes are "the only place a coupling could live".  This probe makes that
precise, and the answer is sharper than the guess in both directions.

THE RESULT
──────────
For two concentric shells at b < a, each deformed by δR = α R P_ℓ,

    ∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)

The prefactor is ℓ(ℓ+1) — the Laplacian eigenvalue on the sphere.  So the
decoupling of the previous round is NOT a separate theorem about spheres: it is
the ℓ = 0 case of a one-line multipole fact, and it vanishes because the
constant mode has zero eigenvalue.

WHAT IS CHECKED
───────────────
T2  THE CLOSED FORM AGAINST BRUTE FORCE.  Direct double integration over both
    deformed surfaces, never expanding in multipoles, agrees to 9e-06 for
    ℓ = 1…4 and gives exactly 0.000000000 at ℓ = 0.

T3  BIRKHOFF IS THE ℓ = 0 CASE.  The ℓ = 0 mutual stiffness is exactly zero at
    every separation, while ℓ = 1, 2, 3 all couple.  The previous round measured
    this zero and imported a theorem for it; here it is a vanishing eigenvalue.

T4  A TRAP CAUGHT FIRST — THE PURE DIPOLE IS NOT A TRANSLATION.  The second
    variation of the area of r = R(1+αP_ℓ) is [2+ℓ(ℓ+1)]/(2ℓ+1) exactly —
    2, 4/3, 8/5, 2, 22/9, 32/11 — which does NOT vanish at ℓ = 1.  A naive
    translation-invariance test built on a pure P₁ deformation would have
    reported a stiffness that is not there.

T5  AND TRANSLATION INVARIANCE IS EXACT ANYWAY.  A rigid displacement is
    r = R + dP₁ − d²/3R + (d²/3R)P₂ + O(d³): the induced ℓ = 0 and ℓ = 2 pieces
    are what preserve the area.  The exact displaced sphere is area-preserving
    to 4e-16 at every d; the truncated family beats the pure dipole by O(d²).

T6  THE AREA COST IS THE EXACT RATIONAL at every ℓ checked.

T7  BUT THE COUPLING IS SCREENED.  The same formula suppresses it as (b/a)^ℓ,
    so at b/a = 0.4 the mutual stiffness falls 544× from ℓ = 1 to ℓ = 8.  The
    multipoles that carry a spin-2 signal are the ones separation screens
    hardest.  Both halves are the answer.

T8  AND THE SHEAR RESPONSE IS AN INPUT.  A perfect-fluid shell is
    S_ij = diag(−σ, p, p): it resists area change and nothing else.  Resisting
    SHAPE change at fixed area needs an elastic modulus no equation of state
    supplies, so μ_s is carried explicitly and is zero unless someone sets it.

SCOPE
─────
What is established is the multipole structure of the coupling and where
Birkhoff sits inside it.  What is NOT established: ℓ ≥ 2 dynamics on the
Schwarzschild background, quasinormal frequencies, or that any particular shell
would resonate.  The constitutive gap is named rather than closed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.shells.multipole import (
    measure_a_fluid_shell_has_no_shear_modulus,
    measure_birkhoff_is_the_ell_zero_case,
    measure_the_area_cost_of_a_deformation,
    measure_the_coupling_is_screened_by_separation,
    measure_the_mutual_coupling_is_the_laplacian_eigenvalue,
    measure_the_pure_dipole_is_not_a_translation,
    measure_translation_invariance_is_exact,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("shell_junction found no throat-shell coupling at l = 0 "
                     "and imported Birkhoff to explain it; does l >= 2 supply "
                     "the coupling, and on what terms?"),
        "scope": ("static Newtonian (Laplace) multipole problem, the weak-field "
                  "limit of the junction problem; not Regge-Wheeler/Zerilli"),
        "pass": True,
    }


def t2_the_closed_form_against_brute_force() -> dict:
    r = measure_the_mutual_coupling_is_the_laplacian_eigenvalue()
    return {"name": "T2_the_closed_form_against_brute_force", **r,
            "pass": bool(r["the_closed_form_is_confirmed"])}


def t3_birkhoff_is_the_ell_zero_case() -> dict:
    r = measure_birkhoff_is_the_ell_zero_case()
    return {"name": "T3_birkhoff_is_the_ell_zero_case", **r,
            "pass": bool(r["ell_zero_is_exactly_zero_at_every_separation"]
                         and r["every_other_ell_couples"])}


def t4_the_pure_dipole_is_not_a_translation() -> dict:
    r = measure_the_pure_dipole_is_not_a_translation()
    return {"name": "T4_the_pure_dipole_is_not_a_translation", **r,
            "pass": bool(r["the_closed_form_is_confirmed"]
                         and r["the_dipole_area_cost_is_not_zero"])}


def t5_translation_invariance_is_exact() -> dict:
    r = measure_translation_invariance_is_exact()
    return {"name": "T5_translation_invariance_is_exact", **r,
            "pass": bool(r["the_exact_sphere_is_area_preserving"]
                         and r["the_truncated_family_is_preserving_to_order_d4"]
                         and r["the_pure_dipole_is_not"])}


def t6_the_area_cost_is_the_exact_rational() -> dict:
    r = measure_the_area_cost_of_a_deformation()
    return {"name": "T6_the_area_cost_is_the_exact_rational", **r,
            "pass": bool(r["every_value_is_the_exact_rational"])}


def t7_the_coupling_is_screened_by_separation() -> dict:
    r = measure_the_coupling_is_screened_by_separation()
    return {"name": "T7_the_coupling_is_screened_by_separation", **r,
            "pass": bool(r["the_coupling_exists_for_every_ell_at_least_one"]
                         and r["but_it_falls_geometrically"])}


def t8_the_shear_response_is_an_input() -> dict:
    r = measure_a_fluid_shell_has_no_shear_modulus()
    return {"name": "T8_the_shear_response_is_an_input", **r,
            "pass": bool(r["a_fluid_shell_has_no_shear_response"]
                         and r["an_elastic_one_needs_an_extra_input"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_closed_form_against_brute_force(),
             t3_birkhoff_is_the_ell_zero_case(),
             t4_the_pure_dipole_is_not_a_translation(),
             t5_translation_invariance_is_exact(),
             t6_the_area_cost_is_the_exact_rational(),
             t7_the_coupling_is_screened_by_separation(),
             t8_the_shear_response_is_an_input()]
    tests.append(t9_assessment(tests))
    t2, t4, t7 = tests[1], tests[3], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "BIRKHOFF_IS_A_VANISHING_EIGENVALUE"
        verdict = (
            "BIRKHOFF IS A VANISHING EIGENVALUE. The previous round measured "
            "that two concentric surfaces in a vacuum spherical model cannot "
            "talk, and imported Birkhoff's theorem to say why. The multipole "
            "form shows that import was not needed: the mutual stiffness of "
            "two concentric deformed shells is G m_b m_a l(l+1) (b/a)^l / "
            "(a (2l+1)^2), with the LAPLACIAN EIGENVALUE l(l+1) as its "
            "prefactor, so it vanishes identically at l = 0 and nowhere else. "
            "The decoupling is a property of the constant mode, not of "
            f"spheres. The closed form agrees with brute-force double "
            f"integration over both deformed surfaces — which never expands in "
            f"multipoles — to {t2.get('worst_relative_error', 0):.1e}, and "
            "gives exactly zero at l = 0. So l >= 2 does supply the coupling "
            "the earlier round found forbidden, and the answer has a second "
            "half that matters as much as the first: the same formula screens "
            f"it as (b/a)^l, a factor of {t7.get('suppression_from_ell_1_to_ell_8', 0):.0f} "
            "from l = 1 to l = 8 at b/a = 0.4. The multipoles that carry a "
            "spin-2 signal are precisely the ones separation suppresses "
            "hardest, so having the channel and being able to use it are "
            "different statements. Two things had to be got right before any "
            "of this could be trusted. A pure P1 deformation is NOT a "
            "translation past linear order: the area second variation is "
            f"[2+l(l+1)]/(2l+1) exactly, which at l = 1 is "
            f"{t4.get('dipole_second_variation', 0):.4f} and not zero, so a "
            "naive translation-invariance check would have reported a "
            "stiffness that is not there. Translation invariance does hold "
            "exactly — the rigid displacement carries induced l = 0 and l = 2 "
            "pieces, and along that family the exact displaced sphere is "
            "area-preserving to 4e-16 at every displacement. And the shear "
            "response is an INPUT: a perfect-fluid shell resists area change "
            "and nothing else, so resisting shape change at fixed area needs "
            "an elastic modulus no equation of state supplies. It is carried "
            "explicitly and never fitted. The coupling l >= 2 restores is "
            "real; what a shell would do with it is a constitutive choice "
            "spherical symmetry never had to make.")
    else:
        verdict_class = "INCOMPLETE"
        verdict = ("At least one check did not hold; see the failing test "
                   "before quoting any number from this run.")

    return {
        "probe": "multipole_coupling",
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
    out = here / "runs" / f"{ts}_multipole_coupling_probe"
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
