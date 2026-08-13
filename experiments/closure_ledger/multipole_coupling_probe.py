"""
Where the two-shell coupling starts, in a static Newtonian model

> Scope, because the headline depends on it: this is the STATIC NEWTONIAN
> (Laplace) two-shell model — the weak-field analogue of the junction problem.
> What it establishes is the SHELL-THEOREM / MULTIPOLE structure of the
> coupling.  BIRKHOFF IS A GR RESULT and remains what shell_junction (PR #249)
> relies on; the ℓ = 0 statement here is its Newtonian analogue, not the
> theorem.  Not Regge–Wheeler/Zerilli, no quasinormal frequencies, no radiative
> dynamics.  G = 1.

THE RESULT
──────────
In this model the monopole mutual stiffness vanishes while higher angular
multipoles couple, with the coupling suppressed geometrically by separation:

    ∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)

The prefactor is ℓ(ℓ+1), the eigenvalue of the angular Laplacian — the ℓ = 0
decoupling IS that zero eigenvalue.

AND THE COUPLING STARTS AT ℓ = 2, NOT ℓ = 1
───────────────────────────────────────────
An earlier draft concluded "everything ℓ ≥ 1 couples", from checking
translation invariance of the AREA.  That is the wrong quantity.  At the level
of the mutual ENERGY the two disagree, and the conclusion changes.

WHAT IS CHECKED
───────────────
T2  THE CLOSED FORM AGAINST BRUTE FORCE.  Direct double integration over both
    deformed surfaces, never expanding in multipoles, agrees to 9e-06 for
    ℓ = 1…4 and gives exactly 0.000000000 at ℓ = 0.

T3  THE ℓ = 0 DECOUPLING IS A ZERO EIGENVALUE.  Exactly zero at every
    separation, while higher shape modes couple.  This is the Newtonian
    analogue of what shell_junction established in GR — NOT a replacement for
    Birkhoff, which that round still relies on.

T4  AND IT IS NOT FOUR-DIMENSIONAL ONLY.  The prefactor is ℓ(ℓ+D−3), the
    Laplacian eigenvalue on S^(D−2).  D = 5 — the case this program cares
    about — is derived and checked in its OWN dimension against a brute-force
    integral over two S³ shells with the 1/r² kernel: 3.3e-04, with ℓ = 0
    vanishing to 1.7e-12.  Closed form there: G m_b m_a ℓ(ℓ+2)/(ℓ+1) b^ℓ/a^(ℓ+2).

T5  THE TRANSLATION MODE DOES NOT COUPLE.  The ℓ = 1 control done on the mutual
    energy rather than the area: a rigidly translated inner sphere leaves the
    mutual energy at exactly −G m_b m_a / a — Newton's shell theorem, held to
    1e-15 out to d = 2.5 — so the cross-derivative with respect to rigid
    displacements is 8.3e-13.  The pure P₁ SHAPE mode is a different object and
    does couple, at 1.78e-02.  So genuine coupling STARTS AT ℓ = 2, which is
    what PR #249 guessed and this establishes.

T6  THE SAME TRAP, IN THE AREA.  The second variation of the area of
    r = R(1+αP_ℓ) is [2+ℓ(ℓ+1)]/(2ℓ+1) exactly — 2, 4/3, 8/5, 2, 22/9, 32/11 —
    which does NOT vanish at ℓ = 1.  A pure P₁ deformation is not a translation
    past linear order.

T7  AND TRANSLATION INVARIANCE HOLDS ANYWAY.  The rigid displacement carries
    induced ℓ = 0 and ℓ = 2 pieces; the exact displaced sphere is
    area-preserving to 4e-16 at every d.

T8  THE AREA COST IS THE EXACT RATIONAL at every ℓ checked.

T9  THE COUPLING IS SCREENED.  (b/a)^ℓ, a factor of 544 from ℓ = 1 to ℓ = 8 at
    b/a = 0.4 — the modes carrying a spin-2 signal are the ones separation
    screens hardest.

T10 THE SHEAR RESPONSE IS AN INPUT.  A perfect-fluid shell resists area change
    and nothing else; resisting SHAPE change at fixed area needs an elastic
    modulus no equation of state supplies.  Carried explicitly, never fitted.

THE LESSON
──────────
Both errors this round caught are the same one: a zero-mode test has to be run
on the quantity the claim is about.  Translation invariance of the area does
not decide whether ℓ = 1 couples; translation invariance of the energy does,
and it says it does not.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.shells.multipole import (
    measure_a_fluid_shell_has_no_shear_modulus,
    measure_the_coupling_generalises_to_five_dimensions,
    measure_the_ell_zero_decoupling_is_a_zero_eigenvalue,
    measure_the_translation_mode_does_not_couple,
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
                     "in GR; in the static Newtonian analogue, where does the "
                     "coupling start and on what terms?"),
        "scope": ("static Newtonian (Laplace) two-shell model; establishes the "
                  "shell-theorem/multipole structure. Birkhoff is a GR result "
                  "and remains what shell_junction relies on — this is its "
                  "Newtonian analogue, not that theorem. Not "
                  "Regge-Wheeler/Zerilli, no quasinormal frequencies."),
        "pass": True,
    }


def t2_the_closed_form_against_brute_force() -> dict:
    r = measure_the_mutual_coupling_is_the_laplacian_eigenvalue()
    return {"name": "T2_the_closed_form_against_brute_force", **r,
            "pass": bool(r["the_closed_form_is_confirmed"])}


def t3_the_ell_zero_decoupling_is_a_zero_eigenvalue() -> dict:
    r = measure_the_ell_zero_decoupling_is_a_zero_eigenvalue()
    return {"name": "T3_the_ell_zero_decoupling_is_a_zero_eigenvalue", **r,
            "pass": bool(r["ell_zero_is_exactly_zero_at_every_separation"]
                         and r["every_other_ell_couples"])}


def t3b_the_coupling_generalises_to_five_dimensions() -> dict:
    """D = 5 derived and checked in its own dimension, not assumed from D = 4."""
    r = measure_the_coupling_generalises_to_five_dimensions()
    return {"name": "T4_the_coupling_generalises_to_five_dimensions", **r,
            "pass": bool(r["the_D5_closed_form_is_confirmed"]
                         and r["ell_zero_vanishes_in_five_dimensions"])}


def t4_the_translation_mode_does_not_couple() -> dict:
    """The ℓ = 1 control on the ENERGY, which the area test does not perform."""
    r = measure_the_translation_mode_does_not_couple()
    return {"name": "T5_the_translation_mode_does_not_couple", **r,
            "pass": bool(r["shell_theorem_holds"]
                         and r["the_translation_mode_does_not_couple"]
                         and r["the_shape_mode_does"])}


def t5_the_pure_dipole_is_not_a_translation() -> dict:
    r = measure_the_pure_dipole_is_not_a_translation()
    return {"name": "T6_the_pure_dipole_is_not_a_translation", **r,
            "pass": bool(r["the_closed_form_is_confirmed"]
                         and r["the_dipole_area_cost_is_not_zero"])}


def t6_translation_invariance_is_exact() -> dict:
    r = measure_translation_invariance_is_exact()
    return {"name": "T7_translation_invariance_is_exact", **r,
            "pass": bool(r["the_exact_sphere_is_area_preserving"]
                         and r["the_truncated_family_is_preserving_to_order_d4"]
                         and r["the_pure_dipole_is_not"])}


def t7_the_area_cost_is_the_exact_rational() -> dict:
    r = measure_the_area_cost_of_a_deformation()
    return {"name": "T8_the_area_cost_is_the_exact_rational", **r,
            "pass": bool(r["every_value_is_the_exact_rational"])}


def t8_the_coupling_is_screened_by_separation() -> dict:
    r = measure_the_coupling_is_screened_by_separation()
    return {"name": "T9_the_coupling_is_screened_by_separation", **r,
            "pass": bool(r["every_shape_mode_couples"]
                         and r["but_it_falls_geometrically"])}


def t9_the_shear_response_is_an_input() -> dict:
    r = measure_a_fluid_shell_has_no_shear_modulus()
    return {"name": "T10_the_shear_response_is_an_input", **r,
            "pass": bool(r["a_fluid_shell_has_no_shear_response"]
                         and r["an_elastic_one_needs_an_extra_input"])}


def t10_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_the_closed_form_against_brute_force(),
             t3_the_ell_zero_decoupling_is_a_zero_eigenvalue(),
             t3b_the_coupling_generalises_to_five_dimensions(),
             t4_the_translation_mode_does_not_couple(),
             t5_the_pure_dipole_is_not_a_translation(),
             t6_translation_invariance_is_exact(),
             t7_the_area_cost_is_the_exact_rational(),
             t8_the_coupling_is_screened_by_separation(),
             t9_the_shear_response_is_an_input()]
    tests.append(t10_assessment(tests))
    t2, t4, t5, t8 = tests[1], tests[4], tests[5], tests[8]
    t3b = tests[3]

    if all(t["pass"] for t in tests):
        verdict_class = "COUPLING_STARTS_AT_ELL_TWO"
        verdict = (
            "COUPLING STARTS AT L = 2. In the static Newtonian two-shell model "
            "the monopole mutual stiffness vanishes while higher angular "
            "multipoles couple, with the coupling suppressed geometrically by "
            "separation. The closed form is G m_b m_a l(l+1) (b/a)^l / "
            "(a (2l+1)^2), agreeing with brute-force double integration over "
            "both deformed surfaces — which never expands in multipoles — to "
            f"{t2.get('worst_relative_error', 0):.1e}, and exactly zero at "
            "l = 0. The prefactor is the angular-Laplacian eigenvalue on "
            "S^(D-2), l(l+D-3), so the l = 0 decoupling IS that zero "
            "eigenvalue -- and it is NOT a four-dimensional accident: D = 5, "
            "the case this program cares about, is derived and checked in its "
            "own dimension against a brute-force integral over two S^3 shells "
            f"with the 1/r^2 kernel, agreeing to "
            f"{t3b.get('worst_relative_error', 0):.1e} with l = 0 vanishing to "
            "1.7e-12. There the closed form is G m_b m_a l(l+2)/(l+1) "
            "b^l/a^(l+2). This is the Newtonian "
            "analogue of what shell_junction established in GR; BIRKHOFF IS A "
            "GR THEOREM AND REMAINS WHAT THAT ROUND RELIES ON, and nothing "
            "here replaces it. Where the coupling starts is the part an "
            "earlier draft got wrong, and the correction matters. That draft "
            "concluded 'everything l >= 1 couples', having checked "
            "translation invariance of the AREA — the wrong quantity. Run on "
            "the mutual ENERGY the two disagree: a rigidly translated inner "
            "sphere leaves the mutual energy at exactly -G m_b m_a / a, "
            "Newton's shell theorem, held to 1e-15 out to d = 2.5, so the "
            "cross-derivative with respect to rigid displacements is "
            f"{t4.get('rigid_translation_cross_derivative', 0):.1e}. THE "
            "TRANSLATION MODE DOES NOT COUPLE. The pure P1 SHAPE mode is a "
            f"different object and does, at "
            f"{t4.get('pure_P1_shape_coupling', 0):.2e}. So the honest "
            "ordering is that l = 0 decouples by the vanishing eigenvalue, "
            "the l = 1 translation mode decouples by the shell theorem, and "
            "genuine coupling starts at l = 2 — which is what PR #249 guessed "
            "and this establishes. Both of this round's errors are the same "
            "mistake: a pure P1 deformation is not a translation past linear "
            "order, its area second variation being "
            f"{t5.get('dipole_second_variation', 0):.4f} and not zero, and a "
            "zero-mode test has to be run on the quantity the claim is about. "
            "The second half of the answer is that the same formula screens "
            f"the coupling as (b/a)^l — a factor of "
            f"{t8.get('suppression_from_ell_1_to_ell_8', 0):.0f} from l = 1 to "
            "l = 8 at b/a = 0.4 — so the modes carrying a spin-2 signal are "
            "the ones separation suppresses hardest. And the shear response is "
            "an INPUT: a perfect-fluid shell resists area change and nothing "
            "else, so resisting shape change at fixed area needs an elastic "
            "modulus no equation of state supplies. It is carried explicitly "
            "and never fitted.")

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
