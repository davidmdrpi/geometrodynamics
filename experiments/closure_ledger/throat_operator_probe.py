"""
Does a flux-conserving two-mouth throat behave differently from PR #255's model?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is a POINT-SUPPORTED self-adjoint extension: it
> has no interior, no proper length, and therefore no delay -- the Delta of PRs
> #251-#255 is not a parameter of a point extension and does not survive into
> one. The boundary matrix A is a CHOICE (four real parameters), not a
> derivation; shells.junction (PR #249) is what would fix it from a matter
> model, and nothing here computes the exotic-matter bill. NOT DONE: no
> backreaction, no stress tensor, no topology change, no rate, no two-source
> invariant.

THE GAP THIS CLOSES
───────────────────
PR #255 solved a mouth relation self-consistently and said plainly what it was
not: a relation between field VALUES carried by the free Green function, with no
normal-derivative matching, no reflected channel, a 1x1 mouth object where a
conserving junction needs 2x2 unitary, and kappa^2 power throughput. It also
found its poles OFF the real axis and had to separate three thresholds --
existence, convergence, stability -- to say what that meant.

This round replaces the relation with the real object. A point-supported throat
is a SELF-ADJOINT EXTENSION of the Laplacian on S3 minus two points, and von
Neumann's theorem says those are parametrized by a unitary matrix between the
deficiency spaces -- here U(2) -- equivalently, by Krein's formula, a Hermitian
2x2 matrix A:

    M(omega) q = phi_in|_mouths,     M = A - Gamma(omega)
    Gamma = [[g(omega), G_d(omega)], [G_d(omega), g(omega)]]

WHAT IS CHECKED
───────────────
T2  THE GREEN FUNCTION HAS A CLOSED FORM AND A FINITE PART. Fourier-transforming
    PR #254's image sum gives G(chi,omega) = sin(omega(pi-chi))/(4 pi sin chi
    sin(pi omega)) -- REAL on the real axis, poles exactly at omega = n+1 --
    checked against PR #255's branch series to 6.3e-12, an independent
    construction. Its short-distance split is 1/(4 pi chi) + g(omega) + O(chi)
    with g = -(omega/4pi) cot(pi omega), and the remainder is first order in chi
    (ratio 10.0 per decade). That is what makes a point interaction definable:
    the divergence is the universal Coulomb one.

T3  THE BOUNDARY OPERATOR IS UNITARY, WITH BOTH CHANNELS. The Cayley transform
    S = (A - ic)(A + ic)^-1 of any Hermitian A is unitary to 4.4e-16 and inverts
    back to 2.0e-16, so the self-adjoint two-mouth conditions ARE U(2): four real
    parameters. Diagonal is REFLECTION, off-diagonal TRANSMISSION, and
    |r|^2 + |t|^2 = 1 at each mouth to 4.4e-16. A real beta is reciprocal; a
    complex one is not. PR #255's model in the same language has r = 0 and
    |t| = kappa, column norm kappa^2 -- outside U(2) unless kappa = 1.

T4  FLUX CONSERVATION IS EXACTLY HERMITICITY. The radial current through a small
    sphere at mouth j is Im(q_j* phi_j^reg), independent of the sphere, so the
    total absorbed is Im(q^dag A q). For Hermitian A that is zero for EVERY q --
    1.8e-16 over 200 random draws, not on average and not to leading order -- and
    for a purely off-diagonal A what one mouth absorbs the other emits to
    1.7e-16. The directional control's median net flux is 0.54.

T5  AND THEREFORE THE SPECTRUM IS REAL, FOR EVERY COUPLING. Gamma is real
    symmetric on the axis, so M is Hermitian there and det M is a real function
    (imaginary part 1.5e-15 relative). Newton from a grid of COMPLEX seeds --
    the same method PR #255 used to find its poles -- converges only onto the
    real axis: 0 off-axis roots, worst |Im omega| = 4.5e-18. The directional
    control gives 9 roots, ALL 9 off-axis, 2 of them growing, worst
    |Im omega| = 0.68 -- and it is unstable even at kappa = 1, so it is the
    DIRECTIONALITY and not the loss that does it. PR #255's instability was the
    non-conservation, not a throat.

T6  THE COUPLED SPECTRUM INTERLACES THE FREE ONE. For an exchange-symmetric pair
    the secular equation factorizes into g + G_d = alpha + beta and
    g - G_d = alpha - beta; both left-hand sides run monotonically from -infinity
    to +infinity across every unit gap, so each contributes exactly one root:
    EXACTLY TWO coupled frequencies strictly between consecutive free ones, over
    8 gaps.

T7  AND THE FREE SPECTRUM RETURNS WHEN THE THROAT IS SWITCHED OFF -- which is
    ||A|| -> infinity, not A -> 0, because the diagonal of A is an INVERSE
    scattering length. Shift to the nearest free frequency falls like 1/||A||,
    measured exponent 0.999, down to 6.2e-04.

T8  WHERE PR #255 SITS. Its relation is exactly the strictly lower-triangular
    boundary matrix A(omega) = [[0,0],[1/(eta kappa e^{-i omega Delta}),0]]: no
    self-energy hence no reflection; one direction only hence anti-Hermitian with
    the anti-Hermitian part equal to the whole matrix; and frequency-dependent,
    which a boundary condition is not allowed to be. Three defects, each with a
    number, and none of them touches that round's resolvent -- which was exact
    for the model it posed.

THE ANSWER
──────────
Yes, and the difference is the one that mattered. A flux-conserving throat has a
real spectrum for every coupling, interlacing the ESU spectrum two per gap and
returning it on decoupling. The instability PR #255 measured was its own
non-conservation. What is still put in is the boundary matrix itself.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.throat_operator import (
    measure_self_adjointness_makes_the_spectrum_real,
    measure_the_boundary_operator_is_unitary_with_both_channels,
    measure_the_coupled_spectrum_interlaces_the_free_one,
    measure_the_directional_model_is_what_pr255_solved,
    measure_the_flux_balance_is_exactly_hermiticity,
    measure_the_green_function_has_a_finite_part,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("does a genuinely flux-conserving two-mouth throat -- a "
                     "self-adjoint extension, with reflection, transmission "
                     "and a unitary boundary operator -- behave differently "
                     "from PR #255's directional mouth relation? and what is "
                     "its spectrum?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is POINT-SUPPORTED: no interior, no proper "
                  "length, and no delay -- Delta is not a parameter of a point "
                  "extension. The boundary matrix A is a choice of four real "
                  "parameters, not a derivation; PR #249 is what would fix it. "
                  "No backreaction, no stress tensor, no topology change, no "
                  "rate, no two-source invariant."),
        "pass": True,
    }


def t2_the_green_function_has_a_finite_part() -> dict:
    r = measure_the_green_function_has_a_finite_part()
    return {"name": "T2_the_green_function_has_a_finite_part", **r,
            "pass": bool(r["the_closed_form_is_the_branch_series"]
                         and r["the_remainder_is_first_order_in_chi"])}


def t3_the_boundary_operator_is_unitary_with_both_channels() -> dict:
    r = measure_the_boundary_operator_is_unitary_with_both_channels()
    return {"name": "T3_the_boundary_operator_is_unitary_with_both_channels",
            **r,
            "pass": bool(r["the_cayley_transform_is_unitary"]
                         and r["every_mouth_conserves"]
                         and r["both_channels_are_present"]
                         and r["pr255_is_outside_U2_unless_kappa_is_one"])}


def t4_flux_conservation_is_exactly_hermiticity() -> dict:
    r = measure_the_flux_balance_is_exactly_hermiticity()
    return {"name": "T4_flux_conservation_is_exactly_hermiticity", **r,
            "pass": bool(r["flux_is_conserved_identically"]
                         and r["what_one_mouth_absorbs_the_other_emits"]
                         and r["the_control_does_not_conserve"])}


def t5_self_adjointness_makes_the_spectrum_real() -> dict:
    """The result this round exists for."""
    r = measure_self_adjointness_makes_the_spectrum_real()
    return {"name": "T5_self_adjointness_makes_the_spectrum_real", **r,
            "pass": bool(r["the_secular_function_is_real_on_the_axis"]
                         and r["nothing_off_the_axis"]
                         and r["the_control_fails_both_tests"]
                         and r["and_the_control_is_unstable_even_at_unit"
                               "_transmission"])}


def t6_the_coupled_spectrum_interlaces_the_free_one() -> dict:
    r = measure_the_coupled_spectrum_interlaces_the_free_one()
    return {"name": "T6_the_coupled_spectrum_interlaces_the_free_one", **r,
            "pass": bool(r["exactly_two_per_gap"]
                         and r["every_root_strictly_between_free_ones"]
                         and r["both_channel_functions_are_monotone"])}


def t7_the_free_spectrum_returns_on_decoupling() -> dict:
    r = measure_the_coupled_spectrum_interlaces_the_free_one()
    return {"name": "T7_the_free_spectrum_returns_on_decoupling",
            "decoupling": r["decoupling"],
            "asymptotic_exponents": r["asymptotic_exponents"],
            "the_shift_goes_like_one_over_the_boundary_norm":
                r["the_shift_goes_like_one_over_the_boundary_norm"],
            "free_spectrum_recovered": r["free_spectrum_recovered"],
            "off_is_large_A_not_small_A": r["off_is_large_A_not_small_A"],
            "pass": bool(r["the_shift_goes_like_one_over_the_boundary_norm"]
                         and r["free_spectrum_recovered"])}


def t8_where_pr255_sits() -> dict:
    r = measure_the_directional_model_is_what_pr255_solved()
    return {"name": "T8_where_pr255_sits", **r,
            "pass": bool(r["no_reflection_channel"]
                         and r["not_hermitian_at_any_frequency"]
                         and r["anti_hermitian_part_is_the_whole_matrix"]
                         and r["the_boundary_matrix_depends_on_frequency"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_green_function_has_a_finite_part(),
             t3_the_boundary_operator_is_unitary_with_both_channels(),
             t4_flux_conservation_is_exactly_hermiticity(),
             t5_self_adjointness_makes_the_spectrum_real(),
             t6_the_coupled_spectrum_interlaces_the_free_one(),
             t7_the_free_spectrum_returns_on_decoupling(),
             t8_where_pr255_sits()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "A_CONSERVING_THROAT_HAS_A_REAL_SPECTRUM"
        ctl = t5.get("control_directional", {})
        verdict = (
            "YES -- AND THE DIFFERENCE IS THE ONE THAT MATTERED. PR #255 owed a "
            "boundary operator and this is it: a point-supported throat is a "
            "SELF-ADJOINT EXTENSION of the Laplacian on S3 minus two points, "
            "and von Neumann's theorem parametrizes those by a unitary between "
            "the deficiency spaces -- U(2) -- equivalently, by Krein's formula, "
            "a Hermitian 2x2 boundary matrix A. Everything below follows from "
            "that one substitution. FIRST, IT IS DEFINABLE AT ALL: the free "
            "Green function has the closed form sin(omega(pi-chi))/(4 pi sin "
            "chi sin(pi omega)), real on the real axis with poles exactly at "
            "the free spectrum omega = n+1, agreeing with PR #255's branch "
            f"series to {t2.get('worst_abs_error', 0):.1e}; and its "
            "short-distance split is 1/(4 pi chi) + g(omega) + O(chi) with "
            "g = -(omega/4 pi) cot(pi omega), the remainder first order in chi. "
            "The divergence is the universal Coulomb one, so the subtraction is "
            "not a choice. SECOND, THE BOUNDARY OPERATOR IS UNITARY AND HAS "
            "BOTH CHANNELS: the Cayley transform of any Hermitian A is unitary "
            f"to {t3.get('worst_unitarity_defect', 0):.1e} and inverts back, "
            "with reflection on the diagonal, transmission off it, and "
            f"|r|^2 + |t|^2 = 1 at each mouth to "
            f"{t3.get('worst_sum_of_squares_defect', 0):.1e}. PR #255's model "
            "in the same language has r = 0 and |t| = kappa: outside U(2) "
            "unless kappa = 1. THIRD, FLUX CONSERVATION IS EXACTLY HERMITICITY. "
            "The current through a small sphere at mouth j is Im(q_j* "
            "phi_j^reg), independent of the sphere, so the total absorbed is "
            "Im(q^dag A q) -- zero for every q when A = A^dag, measured at "
            f"{t4.get('worst_relative_net_flux', 0):.1e} over "
            f"{t4.get('n_draws')} random draws, with a purely off-diagonal "
            "throat moving flux from one mouth to the other exactly. FOURTH -- "
            "AND THIS IS WHAT THE ROUND EXISTS FOR -- THE SPECTRUM IS REAL FOR "
            "EVERY COUPLING. Gamma is real symmetric on the axis, so M = A - "
            "Gamma is Hermitian there and det M is a real function; Newton from "
            "a grid of complex seeds, the same method PR #255 used to find its "
            f"poles, converges only onto the axis: 0 off-axis roots, worst "
            f"|Im omega| = "
            f"{t5.get('worst_abs_imaginary_over_all_seeds', 0):.1e}. The "
            f"directional control gives {ctl.get('n_roots')} roots of which "
            f"{ctl.get('n_off_axis')} are off-axis and {ctl.get('n_growing')} "
            f"growing, worst |Im omega| = "
            f"{ctl.get('worst_abs_imaginary', 0):.3f} -- and it is unstable "
            "even at unit transmission, so it is the DIRECTIONALITY and not the "
            "loss. PR #255's instability was its own non-conservation, and the "
            "three thresholds that round had to separate collapse here to one "
            "statement: a conserving throat cannot ring up. FIFTH, THE COUPLED "
            "SPECTRUM INTERLACES THE FREE ONE. For an exchange-symmetric pair "
            "the secular equation splits into g + G_d = alpha + beta and "
            "g - G_d = alpha - beta, both monotone from -infinity to +infinity "
            "across every unit gap, so there are EXACTLY TWO coupled "
            "frequencies strictly between consecutive free ones -- verified "
            f"over {len(t6.get('roots_per_gap', {}))} gaps. SIXTH, THE FREE "
            "SPECTRUM RETURNS WHEN THE THROAT IS SWITCHED OFF, and off is "
            "||A|| -> infinity rather than A -> 0, because the diagonal of A is "
            "an INVERSE scattering length and alpha = 0 is a resonant throat "
            "rather than no throat. The shift to the nearest free frequency "
            f"falls like 1/||A||, exponent "
            f"{(t7.get('asymptotic_exponents') or [0])[-1]:.4f}. WHAT IS STILL "
            "PUT IN: the boundary matrix. A is four real numbers chosen, not "
            "derived, and PR #249 is what would fix them from a matter model; "
            "nothing here computes the exotic-matter bill. The throat is "
            "POINT-SUPPORTED, so it has no interior and no proper length, and "
            "the delay Delta of PRs #251-#255 is not a parameter of a point "
            "extension and does not survive into one -- a real loss of "
            "structure relative to those rounds, stated rather than hidden. No "
            "backreaction, no stress tensor, no topology change, no rate, and "
            "no two-source invariant, which is the next step.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "throat_operator", "tests": tests,
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
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, complex):
        return f"{v.real:.6g}{v.imag:+.6g}j"
    if isinstance(v, np.ndarray):
        return np.array2string(v, precision=5)
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    return str(v)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return [_json_default(x) for x in o.tolist()]
    if isinstance(o, complex):
        return {"re": o.real, "im": o.imag}
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_throat_operator_probe"
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
