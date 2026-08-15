"""Does a two-source invariant distinguish a throat from two scatterers?

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is a POINT-SUPPORTED self-adjoint extension: no
> interior, no proper length, no delay. The boundary data is four real numbers
> CHOSEN, not derived; shells.junction (PR #249) is what would fix them from a
> matter model. Everything below is evaluated at a boundary matrix STRICTLY
> INSIDE PR #257's positive cone, with its Loewner margin quoted, and the exact
> antipodal endpoint is tested SEPARATELY rather than approached. NOT DONE: no
> backreaction, no stress tensor, no topology change, no rate.

THE QUESTION PR #253 LEFT
─────────────────────────
Rank counting ended with an explicit statement of what it could not supply: a
quantity that VANISHES when a source is removed rather than merely becoming
underdetermined. Deleting a scalar equation costs a dimension whatever was
deleted -- a theorem about square systems, not about photons. The replacement
has to be a FIELD quantity, and PRs #254-#257 built the field it has to be
written in: an exactly solvable point throat with a known positive sector.

THE OBJECT
──────────
Superposition is exact for a linear field, so no LINEAR functional carries
two-source information. A QUADRATIC one does, in its cross term:

    C  =  Q[phi_A + phi_B] - Q[phi_A] - Q[phi_B]

identically zero if either source is switched off. For static sources Q is the
interaction energy and C is the throat's Green function between the two source
points:

    C(y_A, y_B) = G(y_A,y_B) + Re sum_ij G(y_A,c_i) R_ij G(c_j,y_B),
    R = (C - B Gamma)^-1 B          (= (A - Gamma)^-1 in the finite-A chart)

which is exactly PR #255's requested index: a MATRIX IN A PAIR OF BRANCHES, the
branch being which mouth the field entered and which it left, plus the extra
"neither mouth" channel.

WHAT IS CHECKED
───────────────
T2  IT VANISHES WHEN A SOURCE IS REMOVED. Bilinear in the two source strengths,
    so switching one off gives exactly zero -- not an underdetermined system.
    Non-vacuous: with both present the smallest sampled value is 3.7e-2.

T3  THE THROAT CHANNEL IS RANK TWO AT ANY SOURCE COUNT. The N x N table of
    throat-mediated cross terms is V^T S V with V of shape 2 x N, so its rank is
    at most two however many sources are placed, while the direct table has full
    rank N. Off the chart the statement needs care: rank R = rank B, so PR
    #257's Dirichlet strata are rank one -- but static sources see only Re R,
    and the real part of a rank-one COMPLEX Hermitian matrix is generically rank
    two. Complex and real Dirichlet directions are drawn separately; they give 2
    and 1.

T4  ANISOTROPY IS NOT THE SIGNATURE. Holding the geodesic separation fixed and
    moving one source over the sphere of that radius, the free interaction is
    constant to 8e-17 and the throat's varies by 66% of its mean. Real effect,
    real measurement -- and TWO DISCONNECTED SCATTERERS produce a 69% variation,
    so it detects structure at the mouths, not a connection between them.

T5  THE DISCONNECTED FAMILY IS A SURFACE. The static invariant determines the
    three entries of S = Re R; two independent scatterers have only two knobs.
    Their image satisfies S12 = G0 det S EXACTLY, so

        W  =  S12 / det S  -  G0

    is zero on it and nonzero off it. 0 defect on 200 disconnected draws.

T6  AND ON REAL BETA THE DEFECT IS THE COUPLING ITSELF: W = -beta, to 5.0e-16
    over 120 random draws, independent of the self-energies, the mouth
    separation and the LOEWNER MARGIN. That last one answers PR #255's caution
    that a test built from a resummed field measures the pole rather than the
    source: driven from margin 0.4 to 0.004 the invariant grows 3.8x and W does
    not move at all (drift 2.1e-17).

T7  AND IT IS A PROTOCOL, NOT A FORMULA. An observer who measures interaction
    energies, knows the background and knows where the mouths are -- but is not
    told the boundary data -- recovers S by least squares from source placements
    and then forms W. Three placements suffice; 24 give W to 1.1e-16.

T8  AGAINST THE ROUND: A ONE-FREQUENCY TEST HAS A BLIND SPOT. For complex beta,
    W = -Re beta - (G_d - Re beta)(Im beta)^2 / P, which vanishes away from
    beta = 0 on two branches. PR #257's gate excludes one of them -- Re beta >
    G_d has det(A - Gamma) < 0 and is unstable -- and does NOT exclude the
    other: Re beta < 0 gives INVISIBLE CONNECTED THROATS with |beta| up to 0.25
    sitting strictly inside the cone at margin 0.034. A single-frequency
    two-source test cannot falsify those.

T9  THE REPAIR, MEASURED. Gamma depends on lambda, so the blind surface moves
    with it: two frequencies give six equations for four parameters and the
    boundary matrix comes back exactly (1.2e-15), the blind family included. The
    only thing left unobservable is the SIGN of Im beta, which PR #256 already
    established is a time reversal rather than a gap.

T10 THE ANTIPODAL ENDPOINT, ON ITS OWN. At d = pi, Gamma(0) is negative
    semidefinite with a zero in the symmetric channel, so (A - Gamma(0))^-1 is
    singular as A -> 0 and the invariant DIVERGES like 1/eps. The discriminator
    does not move: W stays 0 to 3.6e-15 across four decades of divergence, and
    W = -beta still holds exactly for a connected antipodal throat. The loudest
    possible two-source signal, carrying no information about whether the mouths
    are connected.

WHAT IS PUT IN
──────────────
The background, the boundary data, and the mouth positions -- the observer is
assumed to know where the mouths are, and G0 enters W explicitly. Nothing here
derives A from matter.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.two_source import (
    measure_anisotropy_is_not_the_signature,
    measure_the_antipodal_endpoint_on_its_own,
    measure_the_blind_spot_of_a_single_frequency_test,
    measure_the_defect_is_the_mouth_mixing_amplitude,
    measure_the_invariant_is_recoverable_from_observations,
    measure_the_invariant_vanishes_when_a_source_is_removed,
    measure_the_throat_channel_has_the_rank_of_the_boundary_condition,
    measure_two_disconnected_scatterers_lie_on_a_surface,
    measure_two_frequencies_reconstruct_the_boundary_matrix,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("PR #253 ended rank counting by naming what it could not "
                     "supply: a quantity that vanishes when a source is "
                     "removed rather than becoming underdetermined. Build it "
                     "as a field quantity, and ask whether it distinguishes a "
                     "throat from two ordinary scatterers."),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is point-supported: no interior, no proper "
                  "length, no delay. Everything is evaluated at a boundary "
                  "matrix strictly inside PR #257's cone with its Loewner "
                  "margin quoted, and the exact antipodal endpoint is tested "
                  "separately. No backreaction, no stress tensor, no topology "
                  "change, no rate."),
        "pass": True,
    }


def t2_it_vanishes_when_a_source_is_removed() -> dict:
    """The property PR #253 said rank counting structurally could not have."""
    r = measure_the_invariant_vanishes_when_a_source_is_removed()
    return {"name": "T2_it_vanishes_when_a_source_is_removed", **r,
            "pass": bool(r["it_vanishes_exactly"]
                         and r["it_is_not_vacuous"]
                         and r["is_bilinear"]
                         and r["inside_the_cone"])}


def t3_the_throat_channel_is_rank_two() -> dict:
    r = measure_the_throat_channel_has_the_rank_of_the_boundary_condition()
    return {"name": "T3_the_throat_channel_is_rank_two", **r,
            "pass": bool(r["the_throat_table_is_rank_two_in_the_chart"]
                         and r["the_direct_table_has_full_rank"]
                         and r["the_complex_response_has_the_rank_of_B"]
                         and r["a_complex_dirichlet_direction_still_fills"
                               "_both_channels"]
                         and r["a_real_dirichlet_direction_drops_the_table"
                               "_to_one"]
                         and r["the_response_is_hermitian_off_the_chart"])}


def t4_anisotropy_is_not_the_signature() -> dict:
    """The obvious observable, and why it decides nothing."""
    r = measure_anisotropy_is_not_the_signature()
    return {"name": "T4_anisotropy_is_not_the_signature", **r,
            "pass": bool(r["the_free_field_is_isotropic"]
                         and r["both_break_it"]
                         and r["anisotropy_does_not_discriminate"])}


def t5_the_disconnected_family_is_a_surface() -> dict:
    r = measure_two_disconnected_scatterers_lie_on_a_surface()
    return {"name": "T5_the_disconnected_family_is_a_surface", **r,
            "pass": bool(r["the_disconnected_family_is_a_surface"]
                         and r["connected_throats_are_off_it"])}


def t6_the_defect_is_the_coupling() -> dict:
    """W = -beta, and it does not move with the resonance."""
    r = measure_the_defect_is_the_mouth_mixing_amplitude()
    return {"name": "T6_the_defect_is_the_coupling", **r,
            "pass": bool(r["W_is_minus_beta"]
                         and r["every_row_is_inside_the_cone"]
                         and r["the_discriminator_does_not_see_the_resonance"])}


def t7_it_is_a_protocol() -> dict:
    r = measure_the_invariant_is_recoverable_from_observations()
    return {"name": "T7_it_is_a_protocol", **r,
            "pass": bool(r["the_protocol_closes"] and r["W_error"] < 1e-9)}


def t8_the_blind_spot() -> dict:
    """Against the round, and the half of it the stability gate does remove."""
    r = measure_the_blind_spot_of_a_single_frequency_test()
    return {"name": "T8_the_blind_spot", **r,
            "pass": bool(r["the_blind_family_is_not_empty"]
                         and r["the_upper_branch_is_excluded_by_the"
                               "_stability_gate"]
                         and r["the_lower_branch_survives_it"]
                         and r["they_are_invisible_at_lambda_zero"]
                         and r["they_are_visible_at_another_frequency"])}


def t9_two_frequencies_reconstruct_the_throat() -> dict:
    r = measure_two_frequencies_reconstruct_the_boundary_matrix()
    return {"name": "T9_two_frequencies_reconstruct_the_throat", **r,
            "pass": bool(r["the_boundary_matrix_is_reconstructed"]
                         and r["even_the_blind_family_is_reconstructed"])}


def t10_the_antipodal_endpoint() -> dict:
    """d = pi tested as itself, not as a limit."""
    r = measure_the_antipodal_endpoint_on_its_own()
    return {"name": "T10_the_antipodal_endpoint", **r,
            "pass": bool(r["the_antipodal_value_is_minus_g0"]
                         and r["the_invariant_diverges_like_one_over_epsilon"]
                         and r["the_defect_stays_zero"]
                         and r["the_identity_survives_the_endpoint"])}


def t11_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T11_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_it_vanishes_when_a_source_is_removed(),
             t3_the_throat_channel_is_rank_two(),
             t4_anisotropy_is_not_the_signature(),
             t5_the_disconnected_family_is_a_surface(),
             t6_the_defect_is_the_coupling(),
             t7_it_is_a_protocol(),
             t8_the_blind_spot(),
             t9_two_frequencies_reconstruct_the_throat(),
             t10_the_antipodal_endpoint()]
    tests.append(t11_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9, t10 = tests[1:10]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_TWO_SOURCE_INVARIANT_MEASURES_THE_MOUTH_MIXING"
        verdict = (
            "THE INVARIANT IS THE CROSS TERM OF A QUADRATIC FUNCTIONAL, AND ITS "
            "DISCONNECTION DEFECT IS MINUS THE MOUTH-MIXING AMPLITUDE. PR #253 "
            "ended rank counting by naming what it could not give: something "
            "that VANISHES when a source is removed rather than becoming "
            "underdetermined. Superposition makes every linear functional "
            "additive, so the object has to be quadratic, and its cross term is "
            "the throat's Green function between the two source points -- "
            "bilinear in the sources, exactly zero when either is switched off, "
            f"and not vacuous ({t2.get('smallest_value_with_both_present', 0):.3g} "
            "at its smallest with both present). Written out it is PR #255's "
            "requested index, A MATRIX IN A PAIR OF BRANCHES: which mouth the "
            "field entered, which it left, plus the channel that used neither. "
            "THE THROAT CHANNEL IS RANK TWO AT ANY SOURCE COUNT -- the N x N "
            "table of throat-mediated cross terms is V^T S V with V of shape "
            f"2 x N, rank {t3.get('chart_throat_rank')} against a direct table "
            f"of rank {t3.get('chart_direct_rank')} for "
            f"{t3.get('n_sources')} sources -- and off the chart rank R = "
            "rank B, though static sources see only Re R, whose rank is two "
            "even for a COMPLEX rank-one boundary condition and one for a real "
            "one. TWO THINGS THAT LOOK LIKE THE SIGNATURE AND ARE NOT. The "
            "cross term being nonzero is interference. And the interaction "
            "being ANISOTROPIC -- depending on more than the geodesic "
            "separation, which the free field on this background cannot do at "
            f"all ({t4.get('free_spread', 0):.1e} spread) -- is real, "
            f"{t4.get('throat_relative_spread', 0):.2f} of the mean, and two "
            "DISCONNECTED scatterers produce "
            f"{t4.get('disconnected_relative_spread', 0):.2f}. It detects "
            "structure at the mouths, not a connection between them. WHAT DOES "
            "DISCRIMINATE IS A PARAMETER COUNT. The static invariant determines "
            "three numbers, the entries of S = Re R; two independent scatterers "
            "have two knobs, so their image is a SURFACE with the exact "
            "equation S12 = G0 det S, satisfied to "
            f"{t5.get('worst_defect_on_the_disconnected_family', 0):.1e} on "
            f"{t5.get('n_draws')} draws. The defect W = S12/det S - G0 is "
            "therefore the discriminator -- and on real beta it is not merely "
            "nonzero but EQUAL TO THE COUPLING: W = -beta to "
            f"{t6.get('worst_error_in_W_plus_beta', 0):.1e}, independent of the "
            "self-energies, the separation, and the Loewner margin. That last "
            "independence answers PR #255's caution that a resummed field "
            "measures the pole instead of the source: driven from margin 0.4 to "
            f"0.004 the invariant grows "
            f"{t6.get('invariant_growth_toward_the_boundary', 0):.1f}x while W "
            f"drifts {t6.get('worst_defect_drift_across_margins', 0):.1e}. AND "
            "IT IS A PROTOCOL: an observer who measures interaction energies "
            "and knows the background and the mouth positions, but is not told "
            "the boundary data, recovers S by least squares from source "
            f"placements and gets W to {t7.get('W_error', 0):.1e} at "
            f"{t7.get('n_observations')} observations, condition number "
            f"{t7.get('condition_number', 0):.1f}. AGAINST THE ROUND: A ONE-"
            "FREQUENCY TEST IS BLIND ON A ONE-PARAMETER FAMILY. For complex "
            "beta, W = -Re beta - (G_d - Re beta)(Im beta)^2/P vanishes away "
            "from beta = 0 on two branches. PR #257's gate excludes one of them "
            "-- Re beta > G_d has det(A - Gamma) < 0 -- and leaves the other: "
            f"connected throats with |beta| up to "
            f"{t8.get('largest_stable_invisible_coupling', 0):.3f}, strictly "
            f"inside the cone at margin "
            f"{t8.get('smallest_stable_invisible_margin', 0):.3f}, that the "
            "static invariant cannot see at all. Not fine-tuned and not "
            "unstable. THE REPAIR IS MEASURED RATHER THAN ASSERTED: Gamma "
            "depends on lambda so the blind surface moves, two frequencies give "
            "six equations for four parameters, and the boundary matrix comes "
            f"back to {t9.get('worst_parameter_error', 0):.1e} -- the blind "
            "family included. What is still not observable is the SIGN of "
            "Im beta, which PR #256 established is a time reversal. FINALLY, "
            "THE ANTIPODAL ENDPOINT ON ITS OWN, because PR #257 showed it is a "
            "different statement rather than a limit: at d = pi, Gamma(0) is "
            "negative semidefinite with a zero in the symmetric channel, so the "
            "static response is singular as A -> 0 and the invariant DIVERGES "
            f"like 1/eps (residue {t10.get('residue_of_the_divergence', 0):.3g}, "
            "flat across four decades). W stays exactly zero through all of it, "
            "and W = -beta still holds for a connected antipodal throat. The "
            "loudest available two-source signal carries no information about "
            "whether the mouths are connected -- size is not evidence. WHAT IS "
            "STILL PUT IN: the background, the mouth positions, and the "
            "boundary data itself, four real numbers chosen and not derived, "
            "with PR #249 still the thing that would fix them from matter. The "
            "throat is point-supported -- no interior, no proper length, no "
            "delay. No backreaction, no stress tensor, no topology change, no "
            "rate. What this round hands the next one is a quantity that is "
            "zero without a second source, equal to the non-local part of the "
            "operator when the throat is time-reversal invariant, and a stated "
            "blind spot with a stated repair.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "two_source", "tests": tests,
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
    out = here / "runs" / f"{ts}_two_source_probe"
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
