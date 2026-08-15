"""
A flux-conserving two-mouth throat: what self-adjointness buys, and what it does not.

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R). The throat is a POINT-SUPPORTED self-adjoint extension: no
> interior, no proper length, and therefore no delay -- the Delta of PRs
> #251-#255 is not a parameter of a point extension and does not survive into
> one. The boundary data is a CHOICE (four real parameters), not a derivation;
> shells.junction (PR #249) is what would fix it from a matter model, and
> nothing here computes the exotic-matter bill. NOT DONE: no backreaction, no
> stress tensor, no topology change, no rate, no two-source invariant.

THE GAP THIS CLOSES
───────────────────
PR #255 solved a mouth relation self-consistently and said plainly what it was
not: a relation between field VALUES carried by the free Green function, with no
normal-derivative matching, no reflected channel, a 1x1 mouth object, and
kappa^2 power throughput.

A point-supported throat is a SELF-ADJOINT EXTENSION of the Laplacian on S3
minus two points. Von Neumann parametrizes the extensions by a unitary between
the deficiency spaces -- U(2) -- equivalently, by Krein's formula, a Hermitian
2x2 matrix. Writing phi = phi_in + sum_j G(chi_j) q_j, the regular part at mouth
j is phi_j^reg = (phi_in)_j + (Gamma q)_j, and a linear boundary condition is a
PAIR of matrices:

    B phi^reg = C q     =>     (C - B Gamma) q = B phi_in|_mouths

so the mouth-active spectrum is det(C - B Gamma) = 0. The extension is
self-adjoint iff rank[B|C] = 2 and B C^dag is Hermitian. The general (B,C) form
is needed because PR #255's relation is NOT of the form phi^reg = A q.

WHAT IS CHECKED
───────────────
T2  THE GREEN FUNCTION HAS A CLOSED FORM AND A FINITE PART.
    G(chi,omega) = sin(omega(pi-chi))/(4 pi sin chi sin(pi omega)) -- real on the
    real axis, poles exactly at omega = n+1, and FINITE at the antipode because
    the numerator's zero cancels sin chi. Checked against PR #255's branch series
    to 6.3e-12, an independent construction. Short-distance split is
    1/(4 pi chi) + g(omega) + O(chi) with g = -(omega/4pi) cot(pi omega),
    remainder first order in chi (ratio 10.0 per decade). The divergence is the
    universal Coulomb one, so the subtraction is forced rather than chosen.

T3  HERMITICITY IS FLUX CONSERVATION -- AND THE CAYLEY ENTRIES ARE NOT
    AMPLITUDES. The current through a small sphere at mouth j is
    Im(q_j* phi_j^reg), independent of the sphere, so the total absorbed is
    Im(q^dag A q): zero for EVERY q when A = A^dag (1.8e-16 over 200 draws), with
    a purely off-diagonal throat moving flux from one mouth to the other exactly.
    A non-Hermitian control has median net flux 0.54. The Cayley transform
    S = (A - ic)(A + ic)^-1 is unitary to 2.2e-16 for every reference scale c --
    a real fact about the PARAMETRIZATION -- but its entry magnitudes swing from
    0.955 to 0.713 as c goes 0.05 -> 0.2 for the SAME A, so they are
    boundary-mixing coefficients, not reflection and transmission. An earlier
    version of this module read them as amplitudes; that is corrected here.

T4  SELF-ADJOINTNESS MAKES lambda REAL, NOT omega. Gamma is real symmetric for
    real lambda = omega^2 of EITHER sign, so the secular function is real
    (2.9e-15) and the eigenvalues of the spatial operator are real. They are not
    thereby POSITIVE. Two of the three boundary matrices this module previously
    advertised have lambda < 0: (0.2,-0.13,0.15+0.07i) gives sigma = 2.470532
    and (-0.4,0.07,-0.09+0.31i) gives sigma = 7.090982, i.e. omega = +/- i sigma
    with one member of each pair growing like e^{sigma t}. They were missed
    because the earlier search seeded only Re omega in [1.1,6.9] and discarded
    roots leaving that window -- a search that structurally could not see a root
    on the imaginary axis. THE CLAIM "real spectrum for every coupling, so a
    conserving throat cannot ring up" IS WITHDRAWN.

T5  AND THIS IS WHAT REPLACES IT: THE STABILITY REGION, IN CLOSED FORM. Along
    the imaginary axis both channel functions fall monotonically from their
    lambda = 0 values to -infinity, so a growing mode exists iff alpha +/- beta
    drops below the threshold. The stable set is the wedge

        alpha + beta >= g0 + G0 = +0.02308202
        alpha - beta >= g0 - G0 = -0.07374262
        g0 = -1/(4 pi^2),  G0 = (pi - d)/(4 pi^2 sin d)

    verified against a negative-lambda scan at every one of 221 grid points, 0
    mismatches, with 56 stable and 165 not. Positivity is a condition on the
    boundary data, separate from self-adjointness.

T6  det(C - B Gamma) = 0 IS THE RANK-TWO MOUTH-ACTIVE SECTOR, NOT THE SPECTRUM.
    Level n has degeneracy (n+1)^2 and only two combinations can move, so
    (n+1)^2 - 2 modes per level stay exactly at the free eigenvalue -- 23 of 25
    at level 4. Within the sector there is also a mode BELOW the free ground
    state (0 < lambda < 1) that an omega-scan starting above 1 cannot see, and
    then two per interlacing gap. And the convenient claim that both channel
    functions run -infinity -> +infinity across EVERY gap is false: the m = 1
    pole cancels in the antisymmetric channel, because the constant mode is
    equal at both mouths, so its first-gap endpoint is finite (-0.0383) and a
    root there is conditional on alpha - beta. Scanned.

T7  PR #255's RELATION EMBEDS EXACTLY -- AS A NON-SELF-ADJOINT CONDITION. That
    round set q_1 = 0 and q_2 = gain . phi_1^reg, which is B = [[0,0],[gain,0]],
    C = I, giving det(C - B Gamma) = 1 - gain . G_d: EXACTLY its own 1 - L,
    matched to 3.5e-18. It is maximal (rank[B|C] = 2) but B C^dag = B is not
    Hermitian, and no finite Hermitian A reproduces it -- it needs the singular
    B. An earlier version of this module used A = [[0,0],[1/gain,0]] as the
    control; that gives g^2 - G_d^2 + G_d/gain, a DIFFERENT function (off by
    1.44), so the old control was not PR #255's model and nothing was concluded
    from comparing against it.

T8  AND THE PHASE OF beta IS PHYSICAL, BUT THERE IS NO PREFERRED DIRECTION. The
    secular function is (a1-g)(a2-g) - |beta - G_d|^2 with G_d real, so it
    depends on Re beta: the mouths are joined through the bulk as well as
    through the throat, and that fixes the relative phase. It is invariant under
    beta -> conj(beta), which is time reversal. An earlier version called a
    complex beta "non-reciprocal", reading the Cayley entries as amplitudes;
    withdrawn.

THE ANSWER
──────────
What survives: two point removals on S3 admit a U(2) self-adjoint-extension
family; the regularized Green/Krein matrix is well defined; Hermitian boundary
data conserve the point-boundary flux exactly; and the rank-two mouth-active
sector is described by det(C - B Gamma) = 0. What does NOT survive: that
self-adjointness implies stability. Positivity is separate, it has a closed form
here, and most of the sampled boundary family is outside it.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.throat_operator import (
    measure_hermiticity_is_flux_conservation,
    measure_self_adjointness_makes_lambda_real_not_omega,
    measure_the_green_function_has_a_finite_part,
    measure_the_mouth_active_sector_is_rank_two,
    measure_the_pr255_boundary_condition_is_not_self_adjoint,
    measure_the_spectrum_is_conjugation_symmetric_in_beta,
    measure_the_stability_region_in_the_boundary_family,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("what does making the throat a genuine self-adjoint "
                     "extension buy -- and what does it not? in particular, "
                     "does flux conservation imply stability?"),
        "scope": ("a linear scalar field on a fixed Einstein static universe. "
                  "The throat is POINT-SUPPORTED: no interior, no proper "
                  "length, and no delay -- Delta is not a parameter of a point "
                  "extension. The boundary data is a choice of four real "
                  "parameters, not a derivation; PR #249 is what would fix it. "
                  "No backreaction, no stress tensor, no topology change, no "
                  "rate, no two-source invariant."),
        "pass": True,
    }


def t2_the_green_function_has_a_finite_part() -> dict:
    r = measure_the_green_function_has_a_finite_part()
    return {"name": "T2_the_green_function_has_a_finite_part", **r,
            "pass": bool(r["the_closed_form_is_the_branch_series"]
                         and r["the_remainder_is_first_order_in_chi"]
                         and r["the_antipodal_focus_is_finite"])}


def t3_hermiticity_is_flux_conservation() -> dict:
    """And the Cayley entries are not amplitudes -- a correction."""
    r = measure_hermiticity_is_flux_conservation()
    return {"name": "T3_hermiticity_is_flux_conservation", **r,
            "pass": bool(r["flux_is_conserved_identically"]
                         and r["what_one_mouth_absorbs_the_other_emits"]
                         and r["the_control_does_not_conserve"]
                         and r["the_cayley_transform_is_unitary"]
                         and r["the_cayley_entries_are_not_physical"
                               "_amplitudes"])}


def t4_self_adjointness_makes_lambda_real_not_omega() -> dict:
    """The withdrawn claim, and why the earlier search could not see it."""
    r = measure_self_adjointness_makes_lambda_real_not_omega()
    return {"name": "T4_self_adjointness_makes_lambda_real_not_omega", **r,
            "pass": bool(r["the_secular_function_is_real_in_lambda"]
                         and r["hermiticity_gives_real_lambda"]
                         and r["hermiticity_does_not_give_positivity"])}


def t5_the_stability_region_in_the_boundary_family() -> dict:
    """What replaces it: positivity, mapped."""
    r = measure_the_stability_region_in_the_boundary_family()
    return {"name": "T5_the_stability_region_in_the_boundary_family", **r,
            "pass": bool(
                r["the_channel_functions_are_monotone_on_the_imaginary_axis"]
                and r["the_closed_form_agrees_with_every_probe"]
                and r["the_closed_form_matches_everywhere"]
                and r["both_signs_are_represented"])}


def t6_the_mouth_active_sector_is_rank_two() -> dict:
    r = measure_the_mouth_active_sector_is_rank_two()
    return {"name": "T6_the_mouth_active_sector_is_rank_two", **r,
            "pass": bool(r["at_most_two_modes_move_per_level"]
                         and r["there_is_a_sector_below_the_ground_state"]
                         and r["two_per_interlacing_gap"]
                         and r["the_first_gap_antisymmetric_endpoints"
                               "_are_finite"]
                         and r["existence_in_the_first_gap_is_conditional"])}


def t7_the_pr255_embedding() -> dict:
    r = measure_the_pr255_boundary_condition_is_not_self_adjoint()
    return {"name": "T7_the_pr255_embedding", **r,
            "pass": bool(r["the_embedding_is_exact"]
                         and r["the_old_control_was_a_different_model"]
                         and r["every_embedding_is_maximal"]
                         and r["none_is_self_adjoint"])}


def t8_the_phase_of_beta() -> dict:
    r = measure_the_spectrum_is_conjugation_symmetric_in_beta()
    return {"name": "T8_the_phase_of_beta", **r,
            "pass": bool(r["the_phase_of_beta_is_physical"]
                         and r["the_spectrum_is_conjugation_symmetric"])}


def t9_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T9_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_green_function_has_a_finite_part(),
             t3_hermiticity_is_flux_conservation(),
             t4_self_adjointness_makes_lambda_real_not_omega(),
             t5_the_stability_region_in_the_boundary_family(),
             t6_the_mouth_active_sector_is_rank_two(),
             t7_the_pr255_embedding(),
             t8_the_phase_of_beta()]
    tests.append(t9_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8 = tests[1:8]

    if all(t["pass"] for t in tests):
        verdict_class = "SELF_ADJOINTNESS_IS_CONSERVATION_NOT_STABILITY"
        th = t5.get("thresholds", {})
        verdict = (
            "SELF-ADJOINTNESS BUYS CONSERVATION AND A REAL lambda; POSITIVITY "
            "IS A SEPARATE CONDITION, AND IT IS THE ONE THAT DECIDES "
            "STABILITY. A point-supported throat is a self-adjoint extension of "
            "the Laplacian on S3 minus two points, parametrized by U(2), and "
            "writing the boundary condition as the PAIR B phi^reg = C q makes "
            "the mouth-active spectrum det(C - B Gamma) = 0. FIRST, IT IS "
            "DEFINABLE AT ALL: G(chi,omega) = sin(omega(pi-chi))/(4 pi sin chi "
            "sin(pi omega)), real on the axis, poles exactly at omega = n+1, "
            "finite at the antipode because the numerator's zero cancels sin "
            f"chi, and agreeing with PR #255's branch series to "
            f"{t2.get('worst_abs_error', 0):.1e}; its short-distance split is "
            "1/(4 pi chi) + g(omega) + O(chi), remainder first order in chi, so "
            "the subtraction a point interaction needs is forced rather than "
            "chosen. SECOND, HERMITICITY IS EXACTLY FLUX CONSERVATION: the "
            "current through a small sphere at mouth j is Im(q_j* phi_j^reg), "
            "so the total absorbed is Im(q^dag A q), zero for EVERY q when "
            f"A = A^dag -- {t3.get('worst_relative_net_flux', 0):.1e} over "
            f"{t3.get('n_draws')} random draws, against a median "
            f"{t3.get('control_median_net_flux', 0):.2f} for a non-Hermitian "
            "control. THIRD -- AND THIS IS A CORRECTION -- THE CAYLEY ENTRIES "
            "ARE NOT AMPLITUDES. The transform is unitary for every reference "
            "scale c, which is a real fact about the parametrization, but the "
            "same A gives entry magnitudes spread by "
            f"{t3.get('cayley_diagonal_spread_over_the_reference_scale', 0):.3f} "
            "across c, so they are boundary-mixing coefficients; a closed "
            "universe has no asymptotic region to normalize a scattering matrix "
            "against. FOURTH -- AND THIS IS THE MAIN CORRECTION -- "
            "SELF-ADJOINTNESS MAKES lambda = omega^2 REAL AND NOTHING MORE. "
            "Gamma is real symmetric for real lambda of either sign, so the "
            "secular function is real; but lambda can be NEGATIVE, and then "
            "omega = +/- i sqrt|lambda| with one member of the pair growing. "
            f"{t4.get('n_unstable_examples')} of the three boundary matrices "
            "this module previously advertised do exactly that -- sigma = "
            "2.470532 and 7.090982 -- and they were missed because the earlier "
            "search seeded only Re omega in [1.1, 6.9] and discarded roots "
            "leaving that window, so a root on the imaginary axis was outside "
            "its reach by construction. The claim 'real spectrum for every "
            "coupling, so a conserving throat cannot ring up' is WITHDRAWN. "
            "FIFTH, WHAT REPLACES IT IS A STABILITY REGION WITH A CLOSED FORM: "
            "both channel functions fall monotonically along the imaginary axis "
            "from their lambda = 0 values, so a growing mode exists iff alpha + "
            f"beta < {th.get('symmetric_threshold', 0):+.8f} or alpha - beta < "
            f"{th.get('antisymmetric_threshold', 0):+.8f}, verified against a "
            f"negative-lambda scan at all {t5.get('grid_points')} grid points "
            f"with {t5.get('grid_mismatches')} mismatches and only "
            f"{t5.get('grid_stable')} of them stable. SIXTH, SCOPE: det(C - B "
            "Gamma) = 0 is the RANK-TWO MOUTH-ACTIVE SECTOR, not the spectrum "
            "-- level n has degeneracy (n+1)^2 and only two combinations can "
            f"move, so {t6.get('untouched_modes_at_level_4')} of 25 modes at "
            "level 4 never leave the free eigenvalue. Inside the sector there "
            "is a mode below the free ground state that an omega-scan starting "
            "above 1 cannot see, then two per interlacing gap; and the "
            "convenient claim that both channels run -infinity to +infinity "
            "across every gap is false, since the m = 1 pole cancels in the "
            "antisymmetric channel and a first-gap root is conditional on alpha "
            "- beta. SEVENTH, PR #255's RELATION EMBEDS EXACTLY, as q_1 = 0 "
            "with q_2 = gain . phi_1^reg, i.e. B = [[0,0],[gain,0]], C = I, "
            f"giving det(C - B Gamma) = 1 - gain . G_d to "
            f"{t7.get('worst_embedding_error', 0):.1e} against that round's own "
            "expression. It is maximal but not self-adjoint, and needs the "
            "singular B -- no finite Hermitian A reproduces it. The previous "
            "version of this control was a different function entirely, off by "
            f"{t7.get('worst_old_control_error', 0):.2f}, so nothing is "
            "concluded from it: this is a classification of PR #255's boundary "
            "condition, NOT a diagnosis of its off-axis poles, because a "
            "self-adjoint throat can have growing modes too. EIGHTH, the phase "
            "of beta is physical because Gamma is not diagonal -- the mouths "
            "are joined through the bulk as well as the throat -- while the "
            "spectrum is invariant under beta -> conj(beta), which is time "
            "reversal; the earlier 'non-reciprocal' reading is withdrawn. WHAT "
            "IS STILL PUT IN: the boundary data, four real numbers chosen and "
            "not derived, with PR #249 still the thing that would fix them. The "
            "throat is POINT-supported, so it has no interior and no proper "
            "length, and the delay Delta of PRs #251-#255 does not survive into "
            "a point extension -- a real loss of structure, stated rather than "
            "hidden. No backreaction, no stress tensor, no topology change, no "
            "rate, and no two-source invariant.")
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
                for row in v[:40]:
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
        return "[" + ", ".join(_fmt(x) for x in v[:12]) + "]"
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
