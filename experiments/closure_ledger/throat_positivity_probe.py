"""
For which Hermitian A is the self-adjoint point-throat operator non-negative?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is a POINT-SUPPORTED self-adjoint extension: no
> interior, no proper length, no delay. The boundary data is four real numbers
> CHOSEN, not derived; shells.junction (PR #249) is what would fix them from a
> matter model, and nothing here computes the exotic-matter bill. NOT DONE: no
> backreaction, no stress tensor, no topology change, no rate, no two-source
> invariant.

THE QUESTION PR #256 LEFT
─────────────────────────
That round established that a point-supported throat is a self-adjoint extension
parametrized by a Hermitian 2x2 boundary matrix A, that Hermiticity is exactly
flux conservation -- and that it does NOT imply stability: lambda = omega^2 is
real but need not be positive, and a negative lambda means omega = +/- i
sqrt|lambda| with a growing mode. It mapped the stable region only for the
two-parameter exchange-symmetric slice, by scanning.

The full four-parameter answer is not a scan. It is one inequality.

    IN THE FINITE-A CHART:
    THE OPERATOR IS NON-NEGATIVE  <=>  A >= Gamma(0)   (Loewner order)

    Gamma(0) = [[g0, G0],[G0, g0]],  g0 = -1/(4 pi^2),
                                     G0 = (pi - d)/(4 pi^2 sin d)

and in general, on the whole U(2) family, by restriction to the allowed-charge
subspace:  A_eff >= P^dag Gamma(0) P.

WHAT IS CHECKED
───────────────
T2  THE CRITERION IS EXACT, AND THE REGION IS A LIGHT CONE. For 200 random
    Hermitian A -- all with complex beta, unequal mouths -- the criterion is
    compared against an actual negative-lambda root scan: 0 mismatches, with 19
    stable and 181 not, so both verdicts occur. Writing
    A - Gamma(0) = x0 I + x.sigma, positive semidefiniteness is x0 >= |x|: the
    stable set is a FORWARD LIGHT CONE in the four dimensions of Hermitian
    boundary data, apex at A = Gamma(0). The reason is one line:
    dGamma/dlambda is positive definite below threshold, so every eigenvalue of
    A - Gamma(lambda) is decreasing in lambda while both run to +infinity as
    lambda -> -infinity; an eigenvalue crosses zero below lambda* iff it is
    already negative at lambda*.

T3  AND THE SAME ARGUMENT COUNTS THEM. Nothing in it is special to lambda* = 0:

        #{mouth-active eigenvalues < lambda*}
              = #{negative eigenvalues of A - Gamma(lambda*)}

    for any lambda* below the free ground state lambda = 1. A Krein-type
    INERTIA THEOREM, 0 mismatches in 160 random tests at lambda* = -2, 0, 0.5,
    0.9. Stability is the lambda* = 0 case, and the COUNT comes out of it, not
    just the yes/no.

T4  THE BOUNDARY IS A ZERO MODE, NOT A CONVENTION. On the null surface
    x0 = |x| the matrix A - Gamma(0) is rank one, so lambda = 0 is in the
    spectrum: a static solution supported by the throat, sitting below the free
    ground state. Located independently by root-finding at 1.4e-14, with the
    secular function vanishing to 1.8e-17. At the APEX A = Gamma(0) there are
    TWO zero modes, which is what makes the apex a distinguished point rather
    than an artefact of coordinates.

T5  AND THE INSTABILITY TURNS ON LIKE A SQUARE ROOT. Step a distance eps past
    the boundary and lambda ~ -eps/mu' -- LINEAR, measured
    lambda/eps -> -7.374476 -- so sigma = sqrt|lambda| rises with exponent
    0.50001. The coefficient is not fitted: mu' is the slope of the relevant
    eigenvalue of Gamma at threshold, computed independently, giving -7.374433.

T6  PR #256's WEDGE IS THE x2 = x3 = 0 SLICE. Setting alpha1 = alpha2 and beta
    real is exactly two of the four dimensions, and on that slice the wedge
    alpha +/- beta >= g0 +/- G0 reproduces the cone at all 143 sampled points.
    Applied to GENERAL boundary data by averaging the mouths and dropping
    Im beta, the same rule gets 65 of 400 draws wrong -- which is the practical
    reason the general form was needed.

T7  WHERE THE APEX SITS -- AND WHAT HAPPENS AT THE ANTIPODE. tr Gamma(0) =
    2 g0 = -1/(2 pi^2), INDEPENDENT of the mouth separation; its eigenvalues are
    exactly PR #256's two channel thresholds. For 0 < d < pi,
    det Gamma(0) = g0^2 - G0^2 < 0, so Gamma(0) is indefinite and A = 0 is
    unstable. THE EXACT ANTIPODAL ENDPOINT IS A DIFFERENT STATEMENT: G_d has a
    REMOVABLE singularity at d = pi, with G_pi(w) = w/(4 pi sin pi w) and
    G_pi(0) = +1/(4 pi^2) = -g0, so Gamma(0) = g0[[1,-1],[-1,1]] has eigenvalues
    (2 g0, 0) -- negative SEMIdefinite, not indefinite. At antipodes A = 0 is
    therefore MARGINALLY non-negative, sitting on the cone's boundary with a
    zero mode in the symmetric channel. Measured at d = pi exactly.

T9  THE MONOTONICITY IS A GRAM MATRIX, NOT A SAMPLE. Gamma_ij(lambda) =
    <delta_i, (H0-lambda)^-1 delta_j> up to a lambda-independent coincidence
    subtraction, so dGamma_ij/dlambda = <delta_i, (H0-lambda)^-2 delta_j> is the
    Gram matrix of (H0-lambda)^-1 delta_j: PSD for free, positive definite for
    distinct mouths. That is the proof. Computed independently from the S3
    addition theorem with its analytic tail and compared with the closed form's
    own derivative: 8.1e-12, positive definite at every sampled lambda and
    separation INCLUDING the antipode. The random root scans elsewhere are
    regression checks, not the argument.

T10 AND THE CRITERION EXTENDS BEYOND THE CHART. phi_reg = A q is the pair
    B = I, C = A, so a finite Hermitian A needs B invertible; in
    B = i(U-I), C = U+I that is 1 not in spec U, and on an eigenvector with
    U v = e^{i theta} v the A-eigenvalue is -cot(theta/2), which runs to
    -/+ infinity as theta -> 0. The missing strata are DIRICHLET DIRECTIONS
    (q = 0 in a subspace), not reachable by any finite A. The criterion extends
    by restriction: A_eff >= P^dag Gamma(0) P on the allowed-charge subspace.
    k = 2 is the chart, k = 1 drops one direction and leaves a scalar
    inequality, k = 0 is the free operator with no mouth-active spectrum. Agreed
    with the cone on 60/60 chart draws and with a root scan on 60/60 stratum
    draws; the reduction's own assumptions -- that the row space of B is the
    allowed subspace and that A_eff is Hermitian -- checked at 9.6e-16 and
    1.8e-12.

THE ANSWER
──────────
Non-negative iff A >= Gamma(0). The positive sector is a forward light cone with
apex Gamma(0), its boundary is where a zero mode enters, the instability outside
turns on like sqrt(distance), and the number of growing modes is the number of
negative eigenvalues of A - Gamma(0). The next round has a stated region to work
inside and a count of what goes wrong outside it.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.throat_positivity import (
    cone_fraction,
    measure_the_criterion_extends_to_the_boundary_strata,
    measure_the_monotonicity_is_a_gram_matrix,
    measure_the_boundary_of_the_cone_is_a_zero_mode,
    measure_the_exchange_symmetric_wedge_is_a_slice,
    measure_the_growth_rate_turns_on_with_a_square_root,
    measure_the_inertia_theorem_counts_modes_below_any_threshold,
    measure_the_positive_sector_is_a_shifted_psd_cone,
    measure_where_the_apex_sits_as_the_mouths_separate,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("for which Hermitian A is the self-adjoint point-throat "
                     "operator non-negative? PR #256 showed flux conservation "
                     "does not give stability and mapped a two-parameter slice "
                     "by scanning; what is the answer on all four parameters?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is point-supported: no interior, no proper "
                  "length, no delay. The boundary data is four real numbers "
                  "chosen, not derived; PR #249 is what would fix them. No "
                  "backreaction, no stress tensor, no topology change, no "
                  "rate, no two-source invariant."),
        "pass": True,
    }


def t2_the_positive_sector_is_a_shifted_psd_cone() -> dict:
    """The answer: one inequality, checked against an actual scan."""
    r = measure_the_positive_sector_is_a_shifted_psd_cone()
    return {"name": "T2_the_positive_sector_is_a_shifted_psd_cone", **r,
            "pass": bool(r["the_criterion_is_exact"]
                         and r["both_verdicts_occur"]
                         and r["the_light_cone_form_agrees"])}


def t3_the_inertia_theorem() -> dict:
    """And it counts, at any threshold below the free ground state."""
    r = measure_the_inertia_theorem_counts_modes_below_any_threshold()
    return {"name": "T3_the_inertia_theorem", **r,
            "pass": bool(r["the_inertia_theorem_holds"]
                         and r["d_gamma_d_lambda_is_positive_definite"])}


def t4_the_boundary_is_a_zero_mode() -> dict:
    r = measure_the_boundary_of_the_cone_is_a_zero_mode()
    return {"name": "T4_the_boundary_is_a_zero_mode", **r,
            "pass": bool(r["every_boundary_point_has_a_zero_mode"]
                         and r["the_secular_function_vanishes_there"]
                         and r["the_marginal_mode_sits_at_lambda_zero"]
                         and r["no_boundary_point_is_strictly_unstable"]
                         and r["the_apex_carries_two_zero_modes"])}


def t5_the_growth_rate_turns_on_with_a_square_root() -> dict:
    r = measure_the_growth_rate_turns_on_with_a_square_root()
    return {"name": "T5_the_growth_rate_turns_on_with_a_square_root", **r,
            "pass": bool(r["exponent_is_one_half"]
                         and r["lambda_is_linear_in_epsilon"]
                         and r["the_coefficient_matches_the_eigenvalue_slope"])}


def t6_the_wedge_is_a_slice() -> dict:
    r = measure_the_exchange_symmetric_wedge_is_a_slice()
    return {"name": "T6_the_wedge_is_a_slice", **r,
            "pass": bool(r["the_wedge_is_exactly_the_slice"]
                         and r["the_slice_really_is_x2_equals_x3_equals_zero"]
                         and r["the_slice_rule_is_not_enough_in_general"])}


def t7_where_the_apex_sits() -> dict:
    r = measure_where_the_apex_sits_as_the_mouths_separate()
    return {"name": "T7_where_the_apex_sits", **r,
            "pass": bool(r["trace_is_separation_independent"]
                         and r["trace_matches_minus_one_over_two_pi_squared"]
                         and r["the_apex_is_indefinite_away_from_the_antipode"]
                         and r["the_zero_matrix_is_unstable_away_from"
                               "_the_antipode"]
                         and r["the_apex_is_negative_semidefinite_at"
                               "_the_antipode"]
                         and r["the_antipodal_endpoint_is_marginal"]
                         and r["at_the_antipode_A_zero_sits_on_the_boundary"]
                         and r["eigenvalues_are_the_channel_thresholds"])}


def t8_how_big_is_it() -> dict:
    """A cone is unbounded, so the number only means something with its box."""
    r = cone_fraction(1.3, half_width=0.2, n_draws=4000)
    return {"name": "T8_how_big_is_it", **r,
            "pass": bool(0.0 < r["fraction"] < 1.0)}


def t9_the_monotonicity_is_a_gram_matrix() -> dict:
    """The proof, not a sample of it."""
    r = measure_the_monotonicity_is_a_gram_matrix()
    return {"name": "T9_the_monotonicity_is_a_gram_matrix", **r,
            "pass": bool(r["the_gram_sum_is_the_closed_form_derivative"]
                         and r["positive_definite_everywhere"]
                         and r["including_at_the_antipode"])}


def t10_the_criterion_extends_beyond_the_chart() -> dict:
    """phi_reg = A q is a chart of U(2); here is the rest of it."""
    r = measure_the_criterion_extends_to_the_boundary_strata()
    return {"name": "T10_the_criterion_extends_beyond_the_chart", **r,
            "pass": bool(r["the_general_form_agrees_with_the_cone_on_the_chart"]
                         and r["the_general_form_agrees_with_the_scan_on_the"
                               "_strata"]
                         and r["the_reduction_is_legitimate"]
                         and r["the_free_stratum_is_non_negative"])}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_positive_sector_is_a_shifted_psd_cone(),
             t3_the_inertia_theorem(),
             t4_the_boundary_is_a_zero_mode(),
             t5_the_growth_rate_turns_on_with_a_square_root(),
             t6_the_wedge_is_a_slice(),
             t7_where_the_apex_sits(),
             t8_how_big_is_it(),
             t9_the_monotonicity_is_a_gram_matrix(),
             t10_the_criterion_extends_beyond_the_chart()]
    tests.append(t11_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9, t10 = tests[1:10]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_POSITIVE_SECTOR_IS_A_LIGHT_CONE_AT_GAMMA_ZERO"
        verdict = (
            "NON-NEGATIVE IF AND ONLY IF A >= Gamma(0), IN THE LOEWNER ORDER -- "
            "in the finite-A chart, with the general form on the whole U(2) "
            "family got by restricting to the allowed-charge subspace. "
            "PR #256 left the question open on three of its four parameters and "
            "answered the fourth by scanning; the general answer is a single "
            "inequality, and the reason is one line. Gamma(lambda) has "
            "dGamma/dlambda POSITIVE DEFINITE below threshold, so every "
            "eigenvalue of M(lambda) = A - Gamma(lambda) is strictly decreasing "
            "in lambda, while as lambda -> -infinity Gamma -> -(sigma/4pi) I "
            "and both eigenvalues run to +infinity. An eigenvalue therefore "
            "crosses zero somewhere below threshold if and only if it is "
            "already negative AT threshold. Checked against an actual "
            f"negative-lambda root scan on {t2.get('n_draws')} random Hermitian "
            "A -- every one with complex beta and unequal mouths, so all four "
            f"parameters are exercised -- with {t2.get('n_mismatches')} "
            f"mismatches and {t2.get('n_stable')} of them stable, so both "
            "verdicts occur. THE GEOMETRY IS A LIGHT CONE: Hermitian 2x2 "
            "matrices are R^4 under A - Gamma(0) = x0 I + x.sigma, and positive "
            "semidefiniteness is x0 >= |x|, so the stable set is the forward "
            "light cone with apex at A = Gamma(0) -- convex, closed under "
            "positive scaling from the apex, and four-dimensional. AND THE SAME "
            "ARGUMENT COUNTS: #{mouth-active eigenvalues < lambda*} = "
            "#{negative eigenvalues of A - Gamma(lambda*)} for any lambda* "
            f"below the free ground state, {t3.get('n_mismatches')} mismatches "
            f"in {t3.get('n_tested')} tests at lambda* = -2, 0, 0.5 and 0.9. "
            "That is a Krein-type inertia theorem, and stability is its "
            "lambda* = 0 case. THE BOUNDARY IS DETECTABLE, NOT CONVENTIONAL: on "
            "the null surface A - Gamma(0) is rank one, so lambda = 0 enters "
            "the spectrum as a genuine ZERO MODE -- a static solution supported "
            "by the throat, below the free ground state -- located independently "
            f"by root-finding at {t4.get('worst_marginal_lambda', 0):.1e} with "
            "the secular function vanishing to "
            f"{t4.get('worst_secular_at_zero', 0):.1e}. At the apex there are "
            "TWO. AND THE INSTABILITY OUTSIDE TURNS ON LIKE A SQUARE ROOT: "
            f"lambda is linear in the distance eps past the boundary "
            f"({t5.get('lambda_over_epsilon_limit', 0):.6f}, against "
            f"{t5.get('predicted_from_the_eigenvalue_slope', 0):.6f} predicted "
            "from the eigenvalue slope rather than fitted), so the growth rate "
            f"rises with exponent {t5.get('asymptotic_exponent', 0):.5f}. "
            "PR #256's WEDGE IS THE x2 = x3 = 0 SLICE, reproduced exactly at "
            f"all {t6.get('n_slice_points')} sampled points -- and applied to "
            "general boundary data by averaging the mouths and dropping Im "
            f"beta it gets {t6.get('n_the_wedge_rule_gets_wrong')} of "
            f"{t6.get('n_general_draws')} draws wrong, which is why the general "
            "form was needed rather than a wider scan. FINALLY, WHERE THE APEX "
            "SITS: tr Gamma(0) = 2 g0 = -1/(2 pi^2) at EVERY mouth separation, "
            "its eigenvalues are exactly PR #256's two channel thresholds, and "
            "and its determinant g0^2 - G0^2 is negative for 0 < d < pi -- so "
            "Gamma(0) is indefinite there and A = 0 is unstable. THE EXACT "
            "ANTIPODE IS A DIFFERENT STATEMENT, and for this geometry it is the "
            "load-bearing one: G_d has a REMOVABLE singularity at d = pi, with "
            "G_pi(0) = +1/(4 pi^2) = -g0, so Gamma(0) = g0[[1,-1],[-1,1]] has "
            "eigenvalues (2 g0, 0) -- negative SEMIdefinite, not indefinite -- "
            "and A = 0 is MARGINALLY non-negative there, sitting on the cone's "
            "boundary with a zero mode in the symmetric channel. TWO FURTHER "
            "TIGHTENINGS. First, the monotonicity everything rests on is not a "
            "sampled fact: Gamma_ij(lambda) = <delta_i,(H0-lambda)^-1 delta_j> "
            "up to a lambda-independent subtraction, so dGamma/dlambda is the "
            "GRAM MATRIX of (H0-lambda)^-1 delta_j -- PSD for free, positive "
            "definite for distinct mouths -- computed independently from the S3 "
            f"addition theorem and agreeing with the closed form to "
            f"{t9.get('worst_abs_error', 0):.1e}, antipode included. Second, "
            "A >= Gamma(0) is the criterion IN A CHART: phi_reg = A q needs B "
            "invertible, and the strata it misses are Dirichlet directions, "
            "reached as ||A|| -> infinity and not represented by any finite "
            "Hermitian A. The general criterion is A_eff >= P^dag Gamma(0) P on "
            "the allowed-charge subspace, agreeing with the cone on "
            f"{t10.get('chart_draws')} chart draws and with a root scan on "
            f"{t10.get('stratum_draws')} stratum draws, with the reduction's own "
            "assumptions checked rather than assumed. HOW BIG THE REGION IS "
            "depends on the box, and "
            f"the box is stated: {t8.get('fraction', 0):.3f} of a uniform draw "
            f"over |alpha_j|, |Re beta|, |Im beta| <= {t8.get('half_width')}. "
            "WHAT IS STILL PUT IN: the boundary data itself, four real numbers "
            "chosen and not derived, with PR #249 still the thing that would "
            "fix them from matter. The throat is point-supported -- no "
            "interior, no proper length, no delay. No backreaction, no stress "
            "tensor, no topology change, no rate, and no two-source invariant; "
            "what this round buys the next one is a stated region to work "
            "inside and a count of what goes wrong outside it.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "throat_positivity", "tests": tests,
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
    if isinstance(v, complex):
        return f"{v.real:.6g}{v.imag:+.6g}j"
    if isinstance(v, np.ndarray):
        return np.array2string(v, precision=5)
    if isinstance(v, dict):
        return ", ".join(f"{a}={_fmt(b)}" for a, b in v.items())
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in list(v)[:12]) + "]"
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
    out = here / "runs" / f"{ts}_throat_positivity_probe"
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
