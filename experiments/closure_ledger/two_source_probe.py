"""Static two-source throat tomography -- NOT yet the two-wave invariant.

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is a POINT-SUPPORTED self-adjoint extension: no
> interior, no proper length, no delay. The boundary data is four real numbers
> CHOSEN, not derived. Everything is evaluated at a boundary matrix STRICTLY
> INSIDE PR #257's positive cone with its Loewner margin quoted, and the exact
> antipodal endpoint is tested SEPARATELY. NOT DONE: no backreaction, no stress
> tensor, no topology change, no rate, and NOT the two-wave collision invariant.

WHAT THIS IS NOT
────────────────
This is NOT the roadmap's two-wave invariant. The object built here is a STATIC
source-interaction kernel at a fixed spectral parameter. It carries NO local
null momenta, so it cannot distinguish equal-energy collinear from
counterpropagating waves -- which was the entire load-bearing control behind
C = I_A I_B (k_A . k_B)^2. The dynamical object, built from T_A^{mu nu}
T^B_{mu nu} and resolved on the geodesic/winding branches of PRs #253-#255, is
still owed and is the next PR.

The index pair (i,j) here labels MOUTH CHANNELS -- which mouth the field entered
and which it left -- and NOT the branches of #253-#255, which are geodesic
histories with winding numbers and Maslov signs.

WHAT IS BUILT
─────────────
PR #253 ended rank counting by naming what it could not supply: a quantity that
VANISHES when a source is removed rather than merely becoming underdetermined.
Superposition makes every linear functional additive, so the object has to be
quadratic, and its cross term is the throat's Green function between the two
source points:

    C(y_A,y_B) = G(y_A,y_B) + Re sum_ij G(y_A,c_i) R_ij G(c_j,y_B),
    R = (C - B Gamma)^-1 B      (= (A - Gamma)^-1 in the finite-A chart)

WHAT IS CHECKED
───────────────
T2  IT VANISHES WHEN A SOURCE IS REMOVED -- computed the honest way. The
    quadratic functional is built with explicit source strengths and its own
    self-energy terms, and the cross term is Q[a,b] - Q[a,0] - Q[0,b] from three
    separate evaluations, matching a*b*C to 2.8e-17. Removing a source is
    evaluating the same functional at b = 0, not multiplying the answer by zero.

T3  THE THROAT CHANNEL IS RANK TWO AT ANY SOURCE COUNT. The N x N table is
    V^T S V with V of shape 2 x N -- rank 2 against a direct table of rank 12 at
    12 sources. Off the chart rank R = rank B, but static sources see only Re R,
    and the real part of a rank-one COMPLEX Hermitian matrix is generically rank
    two: complex and real Dirichlet directions give 2 and 1, drawn separately.

T4  ANISOTROPY IS NOT THE SIGNATURE. Free interaction constant to 8e-17 on a
    fixed-separation sphere, throat varies 66% -- and TWO DISCONNECTED
    SCATTERERS vary 69%. It detects structure at the mouths, not a connection.
    The off-diagonal response block is the same trap: Gamma fills it for
    diagonal boundary data too, so it is a CROSS-MOUTH channel, not "through the
    throat".

T5  THE DISCONNECTED FAMILY IS A SURFACE. Three observables (the entries of
    S = Re R), two knobs, so the image satisfies S12 = G0 det S exactly -- 0 to
    1.4e-16 on 200 draws. W = S12/det S - G0 is its defining function.

T6  AND ON REAL BETA, W = -beta to 5.0e-16, independent of the self-energies,
    the separation and the LOEWNER MARGIN: driven from margin 0.4 to 0.004 the
    invariant grows 3.8x and W drifts 2.1e-17, answering PR #255's caution.
    SCOPE, exactly: W detects OFF-DIAGONAL MOUTH-BOUNDARY MIXING relative to the
    diagonal two-scatterer null model, inside this point-interaction model.

T7  AND IT IS A PROTOCOL. An observer who measures interaction energies and
    knows the background and the mouth positions, but not the boundary data,
    recovers S by least squares and gets W to 1.1e-16 from 24 placements.

T8  THE BLIND FAMILY, AND THE TWO THINGS THAT REMOVE IT. For complex beta,
    W = -Re beta - (G_d - Re beta)(Im beta)^2/P vanishes away from beta = 0 on
    two branches. PR #257's gate removes one (Re beta > G_d has
    det(A - Gamma) < 0). Reality of the field removes the rest: every blind
    point needs Im beta != 0. The stable half has couplings COMPARABLE TO AND
    SMALLER THAN the self-energies (0.215-0.254 against 0.25-0.40), which an
    earlier draft had backwards.

T9  WHICH FIELD IS BEING SOLVED. PR #254's is a REAL scalar, and a real solution
    needs the self-adjoint domain to be conjugation-invariant: A = A*, hence
    beta real. Measured, not argued -- with complex beta a real unit static
    source produces a field with Im = -2.4e-03. So the blind family belongs to a
    deliberately time-reversal-breaking COMPLEX-scalar extension, and for the
    arc's actual field content W = -beta settles the question at one spectral
    parameter.

T10 AND EVEN THERE, THE LIMITATION IS THE PROTOCOL. Real static sources see only
    Re R -- three numbers for four parameters. PHASE-SENSITIVE COMPLEX SOURCES
    see both quadratures, hence the full complex R, and then A = Gamma + R^-1 at
    ONE spectral parameter: reconstructed to 3.9e-15. The multi-parameter
    requirement is a restriction of the real-static-source protocol, not of the
    operator.

T11 THE REAL-STATIC-SOURCE RECONSTRUCTION, at two POSITIVE spectral parameters
    below the free ground state (lambda = omega^2, so a negative lambda is an
    imaginary frequency, not a second driving frequency). Six equations, four
    unknowns, boundary matrix back to 1.1e-16 -- with several starts, because a
    single start does land in a local minimum and the reported residual is what
    catches it.

T12 THE ANTIPODAL ENDPOINT, ON ITS OWN. At d = pi the static response is
    singular as A -> 0 and the invariant DIVERGES like 1/eps (residue 0.003375,
    flat across four decades) while W stays 0 to 3.6e-15. The loudest available
    signal carries no information about whether the mouths are connected.

WHAT IS PUT IN
──────────────
The background, the mouth positions, and the boundary data. Nothing here derives
A from matter.
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
    measure_a_real_field_forces_beta_real,
    measure_phase_sensitive_sources_need_only_one_spectral_parameter,
    measure_two_spectral_parameters_reconstruct_the_boundary_matrix,
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
        "what_this_is_not": ("NOT the roadmap's two-wave collision invariant. "
                             "This is a static source-interaction kernel with "
                             "no local null momenta, so it cannot separate "
                             "collinear from counterpropagating waves. The "
                             "index (i,j) labels mouth channels, not the "
                             "geodesic/winding branches of PRs #253-#255."),
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
    """The blind family, and the two conditions that shrink it to nothing."""
    r = measure_the_blind_spot_of_a_single_frequency_test()
    return {"name": "T8_the_blind_spot", **r,
            "pass": bool(r["the_blind_family_is_not_empty"]
                         and r["the_upper_branch_is_excluded_by_the"
                               "_stability_gate"]
                         and r["the_lower_branch_survives_the_stability_gate"]
                         and r["but_no_blind_point_is_real_field_compatible"]
                         and r["they_are_invisible_at_lambda_zero"]
                         and r["they_are_visible_at_a_second_spectral"
                               "_parameter"]
                         and r["every_stable_coupling_is_smaller_than_its"
                               "_self_energies"])}


def t9_a_real_field_forces_beta_real() -> dict:
    """Which field is being solved, and what it costs the blind family."""
    r = measure_a_real_field_forces_beta_real()
    return {"name": "T9_a_real_field_forces_beta_real", **r,
            "pass": bool(r["a_real_beta_gives_a_real_field"]
                         and r["a_complex_beta_does_not"]
                         and r["so_for_PR254s_field_there_is_no_blind_family"])}


def t10_one_spectral_parameter_suffices_with_phase() -> dict:
    """The multi-lambda requirement is the protocol's, not the operator's."""
    r = measure_phase_sensitive_sources_need_only_one_spectral_parameter()
    return {"name": "T10_one_spectral_parameter_suffices_with_phase", **r,
            "pass": bool(r["the_quadratures_give_the_kernel"]
                         and r["one_spectral_parameter_suffices"])}


def t11_two_spectral_parameters_reconstruct_the_throat() -> dict:
    r = measure_two_spectral_parameters_reconstruct_the_boundary_matrix()
    return {"name": "T11_two_spectral_parameters_reconstruct_the_throat", **r,
            "pass": bool(r["the_boundary_matrix_is_reconstructed"]
                         and r["even_the_blind_family_is_reconstructed"])}


def t12_the_antipodal_endpoint() -> dict:
    """d = pi tested as itself, not as a limit."""
    r = measure_the_antipodal_endpoint_on_its_own()
    return {"name": "T12_the_antipodal_endpoint", **r,
            "pass": bool(r["the_antipodal_value_is_minus_g0"]
                         and r["the_invariant_diverges_like_one_over_epsilon"]
                         and r["the_defect_stays_zero"]
                         and r["the_identity_survives_the_endpoint"])}


def t13_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T13_assessment", "n_passed": n, "n_total": len(tests),
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
             t9_a_real_field_forces_beta_real(),
             t10_one_spectral_parameter_suffices_with_phase(),
             t11_two_spectral_parameters_reconstruct_the_throat(),
             t12_the_antipodal_endpoint()]
    tests.append(t13_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 = tests[1:12]

    if all(t["pass"] for t in tests):
        verdict_class = "STATIC_THROAT_TOMOGRAPHY_MEASURES_THE_MOUTH_MIXING"
        verdict = (
            "THIS IS STATIC TWO-SOURCE THROAT TOMOGRAPHY, AND ITS DISCONNECTION "
            "DEFECT IS MINUS THE MOUTH-MIXING AMPLITUDE. FIRST, WHAT IT IS NOT: "
            "not the roadmap's two-wave collision invariant. The object is a "
            "STATIC source-interaction kernel at a fixed spectral parameter, it "
            "carries NO local null momenta, and it therefore cannot distinguish "
            "equal-energy collinear from counterpropagating waves -- the "
            "load-bearing control behind I_A I_B (k_A.k_B)^2. The index (i,j) "
            "labels MOUTH CHANNELS, not the geodesic/winding branches of PRs "
            "#253-#255. The dynamical object built from T_A^{mu nu} T^B_{mu nu} "
            "is still owed and is the next PR; the roadmap entry stays open. "
            "WHAT IS BUILT: PR #253 ended rank counting by naming what it could "
            "not give, something that VANISHES when a source is removed rather "
            "than becoming underdetermined. Superposition makes every linear "
            "functional additive, so the object has to be quadratic, and its "
            "cross term is the throat's Green function between the two source "
            "points. Computed the honest way -- the quadratic functional carries "
            "its own self-energy terms and the cross term is "
            "Q[a,b] - Q[a,0] - Q[0,b] from three separate evaluations, matching "
            f"a*b*C to {t2.get('worst_error_in_Q_minus_Q_minus_Q', 0):.1e}; "
            "removing a source is evaluating the same functional at b = 0, not "
            "multiplying the answer by zero. THE THROAT CHANNEL IS RANK TWO AT "
            f"ANY SOURCE COUNT: rank {t3.get('chart_throat_rank')} against a "
            f"direct table of rank {t3.get('chart_direct_rank')} at "
            f"{t3.get('n_sources')} sources, because the table is V^T S V with V "
            "of shape 2 x N. Off the chart rank R = rank B, but static sources "
            "see only Re R, whose rank is two even for a COMPLEX rank-one "
            "boundary condition and one for a real one. TWO THINGS THAT LOOK "
            "LIKE THE SIGNATURE AND ARE NOT. The cross term being nonzero is "
            "interference. And ANISOTROPY -- the interaction depending on more "
            "than the geodesic separation, which the free field on this "
            f"background cannot do at all ({t4.get('free_spread', 0):.1e}) -- is "
            f"real at {t4.get('throat_relative_spread', 0):.2f} of the mean, and "
            "two DISCONNECTED scatterers give "
            f"{t4.get('disconnected_relative_spread', 0):.2f}. It detects "
            "structure at the mouths, not a connection between them, and the "
            "off-diagonal response block is the same trap one level down: Gamma "
            "fills it for diagonal boundary data too, so it is a CROSS-MOUTH "
            "channel and not 'through the throat'. WHAT DOES DISCRIMINATE IS A "
            "PARAMETER COUNT. The static invariant determines three numbers, the "
            "entries of S = Re R; two independent scatterers have two knobs, so "
            "their image is a SURFACE with the exact equation S12 = G0 det S, "
            f"satisfied to "
            f"{t5.get('worst_defect_on_the_disconnected_family', 0):.1e} on "
            f"{t5.get('n_draws')} draws. The defect W = S12/det S - G0 is its "
            "defining function, and on real beta it EQUALS THE COUPLING: "
            f"W = -beta to {t6.get('worst_error_in_W_plus_beta', 0):.1e}, "
            "independent of the self-energies, the separation and the Loewner "
            f"margin -- driven from margin 0.4 to 0.004 the invariant grows "
            f"{t6.get('invariant_growth_toward_the_boundary', 0):.1f}x while W "
            f"drifts {t6.get('worst_defect_drift_across_margins', 0):.1e}, which "
            "answers PR #255's caution that a resummed field measures the pole "
            "instead of the source. STATED EXACTLY, and this is the claim the "
            "round proves: W DETECTS OFF-DIAGONAL MOUTH-BOUNDARY MIXING RELATIVE "
            "TO THE DIAGONAL TWO-SCATTERER NULL MODEL, inside this "
            "point-interaction model -- not topology, not a traversable "
            "interior. AND IT IS A PROTOCOL: an observer who measures "
            "interaction energies and knows the background and the mouth "
            "positions, but not the boundary data, recovers S by least squares "
            f"and gets W to {t7.get('W_error', 0):.1e} at "
            f"{t7.get('n_observations')} placements. THE BLIND FAMILY, AND THE "
            "TWO THINGS THAT REMOVE IT. For complex beta, "
            "W = -Re beta - (G_d - Re beta)(Im beta)^2/P vanishes away from "
            "beta = 0 on two branches. PR #257's gate removes one, Re beta > "
            "G_d, where det(A - Gamma) < 0. REALITY OF THE FIELD REMOVES THE "
            "REST: PR #254 solves a REAL scalar, a real solution needs the "
            "self-adjoint domain conjugation-invariant, and that is exactly "
            "A = A*, hence beta real -- measured rather than argued, since with "
            "complex beta a real unit static source produces a field with "
            f"Im = {t9.get('rows', [{}, {}])[1].get('imaginary_part_of_the_field', 0):.2e}. "
            "Every blind point needs Im beta != 0, so the family belongs to a "
            "deliberately time-reversal-breaking COMPLEX-scalar extension and "
            "not to the arc's field. Its stable couplings are COMPARABLE TO AND "
            f"SMALLER THAN the self-energies (largest "
            f"{t8.get('largest_stable_invisible_coupling', 0):.3f} at margin "
            f"{t8.get('smallest_stable_invisible_margin', 0):.3f}); an earlier "
            "draft said larger, which is false. AND EVEN INSIDE THAT EXTENSION "
            "THE LIMITATION IS THE PROTOCOL, NOT THE OPERATOR: real static "
            "sources see only Re R, three numbers for four parameters, while "
            "PHASE-SENSITIVE COMPLEX SOURCES see both quadratures and give the "
            "full complex R, whence A = Gamma + R^-1 at ONE spectral parameter, "
            f"reconstructed to {t10.get('worst_boundary_error', 0):.1e}. The "
            "real-static-source reconstruction needs two spectral parameters -- "
            f"both POSITIVE and below the free ground state, since lambda = "
            "omega^2 makes a negative lambda an imaginary frequency rather than "
            "a second driving frequency -- and returns the boundary matrix to "
            f"{t11.get('worst_parameter_error', 0):.1e}, blind family included, "
            "using several starts because a single start does land in a local "
            "minimum and the reported residual is what catches it. FINALLY, THE "
            "ANTIPODAL ENDPOINT ON ITS OWN, because PR #257 showed it is a "
            "different statement rather than a limit: at d = pi the static "
            "response is singular as A -> 0 and the invariant DIVERGES like "
            f"1/eps (residue {t12.get('residue_of_the_divergence', 0):.3g}, flat "
            "across four decades) while W stays exactly zero throughout. The "
            "loudest available two-source signal carries no information about "
            "whether the mouths are connected -- size is not evidence. WHAT IS "
            "STILL PUT IN: the background, the mouth positions, and the boundary "
            "data itself, four real numbers chosen and not derived, with PR #249 "
            "still the thing that would fix them from matter. Everything is "
            "STATIC-SOURCE, so this is an interaction-energy statement and not a "
            "scattering one. The throat is point-supported -- no interior, no "
            "proper length, no delay. No backreaction, no stress tensor, no "
            "topology change, no rate, AND NOT THE TWO-WAVE INVARIANT: what this "
            "round hands the next one is a measured non-locality, W = -beta, and "
            "a stated reason the dynamical object still has to be built.")
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
