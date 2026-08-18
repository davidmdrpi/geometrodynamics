"""The dynamical two-wave invariant, measured against its known WKB limit.

> Scope: a LINEAR scalar field on a FIXED round background (the Einstein static
> universe S3 x R), conformally coupled. The throat is a POINT-SUPPORTED
> self-adjoint extension, taken STRICTLY INSIDE PR #257's stability cone with
> its Loewner margin quoted. NO BACKREACTION: the stress tensor is computed FROM
> the field and never fed back. No topology change, no rate, and the boundary
> data is still four real numbers chosen rather than derived.

WHAT THIS IS FOR
────────────────
PR #258 built a STATIC source-interaction kernel and said plainly it was not the
two-wave invariant: no local null momenta, so it could not tell equal-energy
collinear from counterpropagating waves. This round builds the time-dependent
object and applies exactly that control.

The point is NOT to re-derive the WKB identity. C = A_A^2 A_B^2 (k_A.k_B)^2 with
k_A.k_B = -E_A E_B (1 - cos theta) is known: zero collinear, -2 E_A E_B head-on.
The research content is the DIFFERENCE between the exact time-dependent,
multipath, throat-coupled field and that limit.

WHAT IS SOLVED
──────────────
    phi_s(x,w) = s(w) [ G(chi_xs,w) + sum_ij G(chi_xi,w) R_ij(w) G(chi_js,w) ],
    R(w) = (A - Gamma(w^2))^-1

inverted along the RETARDED CONTOUR Im w = eps, which is exact:
phi(t) = e^{eps t} (1/2pi) int du e^{-iut} phihat(u + i eps). Derivatives are
analytic, not differenced: the four-gradient and Hessian come from radial
derivatives of G plus the sphere's own geometry, grad chi = -(y-(x.y)x)/sin chi
and Hess chi = cot chi (delta_ab - n_a n_b).

WHAT IS CHECKED
───────────────
T2  THE SOLVER IS THE IMAGE SUM. The contour integral reproduces PR #254's
    closed-form winding-image sum to 3.3e-16 at four arrival times, INCLUDING
    the alternating Maslov signs. Two constructions sharing no code.

T3  THE SOLVED FIELD OBEYS THE CONFORMAL WAVE EQUATION. d_t^2 phi = lap phi -
    phi, residual 4e-16 relative, with and without the throat. Nothing is
    finite-differenced, so this tests the whole derivative construction at once.

T4  AND THE IMPROVED STRESS TENSOR IS TRACELESS, 1.9e-15 relative. Not a
    tautology: box phi is taken FROM THE SOLVE rather than substituted on shell,
    so the trace equals phi(box phi - phi) -- substituting on shell would make
    it vanish identically for any input, measuring nothing.

T5  THE WKB RESULT IS RECOVERED, AS A LIMIT. Head-on -> 4 (error 5.0e-05 at
    omega = 48); collinear -> 0 (1.8e-10). The arriving directions are exactly
    parallel/antiparallel by construction, to 1e-12. The collinear null turns
    out to be far STRONGER than a leading-order statement -- the two wavefronts
    share their normal exactly, so the residue falls faster than any fixed power
    and the exponent steepens. Convergence is part of the measurement: the same
    run with eps at the frequency spacing is wrong by 4-5 orders of magnitude.

T6  *** MULTIPATH DESTROYS THE COLLINEAR NULL. *** Same two sources, same
    observation point, three different answers, selected by WHICH BRANCH has
    arrived:

        A direct     + B direct              N = 1.9e-07   (WKB 0)
        A long-way   + B direct              N = 3.998     (WKB 4)
        A direct     + B via a mouth         N = 0.5650    (geometry 0.5637)

    The winding image propagates the other way round the sphere, so its arrival
    direction is REVERSED and a collinear pair reads head-on. The cross-mouth
    leg emerges from a mouth, at an angle set by that mouth's position, and the
    exact field agrees with that prediction to 0.24% -- a prediction from
    positions, not a fit. The free-propagation control at the same instant has
    NO second arrival at all (energy product 4e-29 against 1.2e-02), which is
    why the comparison is stated as amplitudes and not as a ratio there. That
    control says the MOUTHS create the branch; T7's beta = 0 control says their
    CONNECTION does not.

    So the collinear null is not destroyed by curvature corrections, which are
    1e-7 here. It is destroyed by MULTIPATH, at O(1). That is the
    branch-resolved invariant PR #255 said was needed.

T7  THE CROSS-MOUTH CHANNELS, AUDITED EXPLICITLY -- AND THE beta = 0 CONTROL.
    All four two-leg paths (i,j) are enumerated rather than minimised over: j is
    the mouth the source drives, i the mouth the signal leaves from, and the
    predicted invariant depends ONLY ON i. That gives two distinct predictions,
    0.651935 and 0.563669, so the test discriminates rather than matching one
    number; the field picks the right one at both resolvable delays, to 8.1e-04
    relative.

    And then the control PR #258's review taught this arc to run first: the same
    measurement with beta = 0, two DISCONNECTED mouths. The invariant does not
    move -- swept over beta in [0, 0.26], all inside PR #257's cone, N shifts by
    6.2e-07, a part in 10^6 and FIVE ORDERS BELOW the 0.088 that separates the
    two exit mouths, while the channel's weight moves by 0.6%.

    So the honest scope is the dynamical version of #258's: THIS OBSERVABLE SEES
    STRUCTURE AT THE MOUTHS, NOT THE CONNECTION BETWEEN THEM. What sees the
    connection is still W = -beta, from the low-frequency limit of the same
    solve (T12). The multipath result stands -- a second arrival direction
    destroys the collinear null -- but the throat's non-locality is not what
    supplies it.

T8  T[phi_A + phi_B] AND THE INTERFERENCE TENSOR. T is quadratic, so the
    two-wave content of the TOTAL stress tensor is the bilinear cross term
    dT = T[phi_A+phi_B] - T[phi_A] - T[phi_B], obtained from three evaluations
    of the same functional. It is traceless (1.8e-15) and EXACTLY zero when
    either source is switched off -- PR #253's missing property, now at tensor
    level.

    And it disagrees with T_A:T_B completely. Normalised,
    dT^00/sqrt(T_A^00 T_B^00) reaches its MAXIMUM 2.000 in the COLLINEAR
    configuration -- two parallel waves add coherently -- which is precisely
    where T_A:T_B vanishes. Head-on, where the invariant is maximal, the
    interference energy is 1.044, roughly half. A backreaction estimate driven
    by C = T_A:T_B would look at the collinear case, see nothing, and be wrong
    about its own source by the size of the whole effect.

T9  THE ARRIVALS ARE THE BRANCH LEDGER. Free arrivals land at chi, 2pi-chi,
    2pi+chi to 1.3e-03 with signs + - +, out of a solve that never saw the
    ledger. The throat adds arrivals the free ledger does not contain, at the
    two-leg times chi(y,c_j) + chi(c_i,x); checked at the CAUSAL ONSET rather
    than at a peak, because R(w) has poles and a throat arrival rings up instead
    of pulsing.

T10 AND THAT RINGING IS THE ONLY TAIL. S3 x R is conformally flat, so the
    conformal scalar obeys Huygens EXACTLY: between geometric arrivals the free
    field is 1.4e-08 of its peak, which is the Gaussian source's own wing. With
    the throat it is 8.1e-02 -- a factor of 5.7e+06. Every bit of tail in this
    model belongs to the throat.

T11 THE CAUSTIC IS WHERE WKB STOPS, WITH A SCALE. Geometric optics gives
    1/(4 pi sin chi), divergent at the antipode; the exact kernel is finite,
    G(pi,w) = w/(4 pi sin pi w), and LINEAR IN w. In between everything depends
    on the single combination w*e with e = pi - chi: the exact/WKB ratio is
    |sin(w e)|, identical across carriers to 6.6e-15. The caustic is cut off at
    e* ~ 1/w.

T12 AND THE LOW-FREQUENCY LIMIT RECOVERS PR #258. int dt phi(t) = phihat(0), so
    the DC content of the solved time series IS the static kernel #258 did its
    tomography on. Running that protocol on numbers produced by the dynamic
    solver -- least squares for S = Re R, then W = S12/det S - G0 -- returns
    -0.060010 against -beta = -0.06. The route goes through the whole contour
    integral, and the O(eps^2) contour bias is Richardson-extrapolated with both
    the raw and corrected numbers reported.

WHAT IS PUT IN
──────────────
The background, the mouth positions, and the boundary data. Nothing derives A
from matter; PR #249 is still what would.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.two_wave import (
    measure_multipath_destroys_the_collinear_null,
    measure_the_arrivals_are_the_branch_ledger_with_maslov_signs,
    measure_the_caustic_is_where_wkb_stops,
    measure_the_improved_stress_tensor_is_traceless,
    measure_the_low_frequency_limit_recovers_the_tomography,
    measure_the_only_tail_is_the_throats,
    measure_the_solved_field_satisfies_the_conformal_wave_equation,
    measure_the_solver_reproduces_the_closed_form_free_field,
    measure_the_wkb_collinear_head_on_result_is_recovered,
    measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth,
    measure_the_interference_tensor_is_largest_where_the_invariant_is_null,
)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "question": ("build the time-dependent two-wave invariant on the "
                     "stable point-throat background and measure it against "
                     "the known WKB collinear/head-on result -- not to "
                     "re-derive that identity, but to quantify how the exact "
                     "multipath throat-coupled field differs from it."),
        "scope": ("a linear conformally coupled scalar on a fixed Einstein "
                  "static universe. The throat is point-supported and taken "
                  "strictly inside PR #257's cone with its Loewner margin "
                  "quoted. No backreaction: the stress tensor is computed from "
                  "the field and never fed back. No topology change, no rate, "
                  "and the boundary data is chosen rather than derived."),
        "pass": True,
    }


def t2_the_solver_is_the_image_sum() -> dict:
    r = measure_the_solver_reproduces_the_closed_form_free_field()
    return {"name": "T2_the_solver_is_the_image_sum", **r,
            "pass": bool(r["the_two_constructions_agree"]
                         and r["the_signs_alternate"])}


def t3_the_wave_equation_holds() -> dict:
    r = measure_the_solved_field_satisfies_the_conformal_wave_equation()
    return {"name": "T3_the_wave_equation_holds", **r,
            "pass": bool(r["the_equation_holds"])}


def t4_the_stress_tensor_is_traceless() -> dict:
    r = measure_the_improved_stress_tensor_is_traceless()
    return {"name": "T4_the_stress_tensor_is_traceless", **r,
            "pass": bool(r["the_tensor_is_traceless"])}


def t5_the_wkb_limit_is_recovered() -> dict:
    r = measure_the_wkb_collinear_head_on_result_is_recovered()
    return {"name": "T5_the_wkb_limit_is_recovered", **r,
            "pass": bool(r["the_directions_are_exactly_parallel"]
                         and r["the_directions_are_exactly_antiparallel"]
                         and r["head_on_converges_to_four"]
                         and r["collinear_converges_to_zero"]
                         and r["the_contour_needs_eps_above_the_spacing"])}


def t6_multipath_destroys_the_collinear_null() -> dict:
    """The result: the invariant is branch-resolved."""
    r = measure_multipath_destroys_the_collinear_null()
    return {"name": "T6_multipath_destroys_the_collinear_null", **r,
            "pass": bool(r["the_direct_pair_is_null"]
                         and r["the_winding_image_reads_head_on"]
                         and r["the_throat_matches_its_geometry"]
                         and r["the_control_has_no_second_arrival"])}


def t7_the_cross_mouth_audit() -> dict:
    """All four (i,j), and the beta = 0 control that scopes the claim."""
    r = measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth()
    return {"name": "T7_the_cross_mouth_audit", **r,
            "pass": bool(r["the_prediction_depends_only_on_the_exit_mouth"]
                         and r["the_field_picks_the_right_one"]
                         and r["the_invariant_is_beta_independent"]
                         and r["every_sweep_point_is_inside_the_cone"])}


def t8_the_interference_tensor() -> dict:
    """T[phi_A + phi_B], and where it disagrees with T_A:T_B."""
    r = measure_the_interference_tensor_is_largest_where_the_invariant_is_null()
    return {"name": "T8_the_interference_tensor", **r,
            "pass": bool(r["delta_T_is_traceless"]
                         and r["delta_T_vanishes_when_a_source_is_removed"]
                         and r["the_interference_is_maximal_where_the"
                               "_invariant_is_null"])}


def t9_the_arrivals_are_the_branch_ledger() -> dict:
    r = measure_the_arrivals_are_the_branch_ledger_with_maslov_signs()
    return {"name": "T9_the_arrivals_are_the_branch_ledger", **r,
            "pass": bool(r["the_free_signs_alternate"]
                         and r["the_free_arrivals_are_sharp"]
                         and r["the_throat_onset_is_causal"]
                         and r["the_throat_arrivals_are_new"])}


def t10_the_only_tail_is_the_throats() -> dict:
    r = measure_the_only_tail_is_the_throats()
    return {"name": "T10_the_only_tail_is_the_throats", **r,
            "pass": bool(r["the_free_field_has_no_tail"]
                         and r["the_throat_has_one"])}


def t11_the_caustic_is_where_wkb_stops() -> dict:
    r = measure_the_caustic_is_where_wkb_stops()
    return {"name": "T11_the_caustic_is_where_wkb_stops", **r,
            "pass": bool(r["the_saturation_is_linear_in_omega"]
                         and r["the_ratio_collapses_in_omega_times_e"])}


def t12_the_low_frequency_limit_recovers_the_tomography() -> dict:
    r = measure_the_low_frequency_limit_recovers_the_tomography()
    return {"name": "T12_the_low_frequency_limit_recovers_the_tomography",
            **r,
            "pass": bool(r["the_bridge_closes"])}


def t13_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T13_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(),
             t2_the_solver_is_the_image_sum(),
             t3_the_wave_equation_holds(),
             t4_the_stress_tensor_is_traceless(),
             t5_the_wkb_limit_is_recovered(),
             t6_multipath_destroys_the_collinear_null(),
             t7_the_cross_mouth_audit(),
             t8_the_interference_tensor(),
             t9_the_arrivals_are_the_branch_ledger(),
             t10_the_only_tail_is_the_throats(),
             t11_the_caustic_is_where_wkb_stops(),
             t12_the_low_frequency_limit_recovers_the_tomography()]
    tests.append(t13_assessment(tests))
    t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 = tests[1:12]

    if all(t["pass"] for t in tests):
        verdict_class = "THE_TWO_WAVE_INVARIANT_IS_BRANCH_RESOLVED"
        verdict = (
            "THE EXACT TWO-WAVE INVARIANT IS BRANCH-RESOLVED, AND THAT IS WHERE "
            "IT DEPARTS FROM WKB. The known limit is not re-derived here; it is "
            "used as the control. Two sources and an observer on one great "
            "circle give arriving directions that are exactly parallel or "
            "exactly antiparallel -- to 1e-12, by construction rather than by "
            "fitting -- so WKB's N = (1 - cos theta)^2 predicts 0 and 4, and "
            f"the exact field returns {t5.get('head_on_at_the_largest_carrier', 0):.6f} "
            f"head-on (error {t5.get('head_on_error', 0):.1e}) and "
            f"{t5.get('collinear_at_the_largest_carrier', 0):.1e} collinear. "
            "THE COLLINEAR NULL IS STRONGER THAN LEADING ORDER: on this geometry "
            "the two wavefronts share their normal exactly, so amplitude "
            "gradients cannot tilt either k, the residue falls faster than any "
            "fixed power, and the measured exponent steepens with omega. "
            "*** AND THAT IS WHAT MAKES THE MULTIPATH RESULT LARGE. *** Holding "
            "the sources and the observation point fixed and changing only WHICH "
            "BRANCH has arrived, the same pair gives "
            f"{t6.get('collinear_floor', 0):.1e} on the direct branches, "
            f"{t6.get('long_way_value', 0):.4f} when A arrives on its LONG-WAY "
            "WINDING IMAGE -- that branch runs the other way round the sphere, "
            "so its arrival direction is reversed and a collinear pair reads "
            f"head-on -- and {t6.get('through_the_throat_value', 0):.4f} when B "
            "arrives VIA A MOUTH, against "
            f"{t6.get('through_the_throat_prediction', 0):.4f} predicted from "
            "the mouth's position alone, a "
            f"{100 * t6.get('throat_relative_error', 0):.2f}% agreement with a "
            "number that was never fitted. The free-propagation control at the "
            "same instant has no second arrival at all -- energy product "
            f"{t6.get('free_control_energy_product', 0):.1e} against "
            f"{t6.get('throat_energy_product', 0):.1e} -- so the mouths are "
            "CREATING the second arrival, not merely bending one. "
            "AND THE (i,j) AUDIT SCOPES THAT. All four two-leg paths are "
            "enumerated rather than minimised over -- j the mouth the source "
            "drives, i the mouth the signal leaves from -- and the predicted "
            "invariant depends ONLY ON i, giving two distinct predictions "
            f"{t7.get('distinct_predictions', [0, 0])[0]:.6f} and "
            f"{t7.get('distinct_predictions', [0, 0])[-1]:.6f}, both matched by "
            f"the field to {t7.get('worst_relative_error', 0):.1e} relative. "
            "Then the control PR #258's review taught this arc to run first: "
            "beta = 0, two DISCONNECTED mouths. The invariant does not move -- "
            f"swept over beta in [0, 0.26], all inside the cone, N shifts by "
            f"{t7.get('beta_sweep_spread', 0):.1e}, which is "
            f"{t7.get('beta_spread_over_the_signal', 0):.1e} of the "
            f"{t7.get('exit_mouth_separation', 0):.4f} that separates the two "
            "exit mouths, while the channel's WEIGHT moves by "
            f"{100 * t7.get('the_weight_moves_instead', 0):.1f}%. So the honest "
            "scope is the dynamical version of #258's: THIS OBSERVABLE SEES "
            "STRUCTURE AT THE MOUTHS, NOT THE CONNECTION BETWEEN THEM. What "
            "sees the connection is W = -beta, from the low-frequency limit of "
            "the same solve. AND T[phi_A + phi_B] IS A DIFFERENT DIAGNOSTIC "
            "ENTIRELY. T is quadratic, so the two-wave content of the TOTAL "
            "stress tensor is the bilinear cross term dT = T[phi_A+phi_B] - "
            "T[phi_A] - T[phi_B], built from three evaluations of the same "
            f"functional: traceless to {t8.get('worst_trace', 0):.1e} and "
            "EXACTLY zero when either source is switched off, which is PR "
            "#253's missing property at tensor level. Normalised, "
            "dT^00/sqrt(T_A^00 T_B^00) reaches its MAXIMUM "
            f"{t8.get('collinear_interference', 0):.3f} in the COLLINEAR "
            "configuration -- two parallel waves add coherently -- precisely "
            f"where T_A:T_B vanishes; head-on it is only "
            f"{t8.get('head_on_interference', 0):.3f}. A backreaction estimate "
            "driven by C = T_A:T_B would look at the collinear case, see "
            "nothing, and be wrong about its own source by the size of the "
            "whole effect. THE CONCLUSION: the collinear null is not spoiled "
            "by curvature corrections, which are 1e-7 here; it is spoiled by "
            "MULTIPATH, at "
            "O(1). The invariant has to carry the branch index PR #255 named, "
            "and a single-branch WKB formula is not merely approximate on this "
            "background -- it is answering a different question. THE SOLVER "
            "EARNS THAT. Its free part reproduces PR #254's closed-form "
            f"winding-image sum to {t2.get('worst_difference', 0):.1e} including "
            "the alternating Maslov signs, two constructions sharing no code; "
            "the solved field satisfies the conformal wave equation to 4e-16 "
            "relative with and without the throat, with nothing "
            "finite-differenced; and the improved stress tensor is traceless to "
            f"{t4.get('worst_relative_trace', 0):.1e} -- which is a real test "
            "because box phi is taken from the solve rather than substituted on "
            "shell, so the trace equals phi(box phi - phi) instead of vanishing "
            "algebraically. Convergence is measured too: with the contour offset "
            "at the frequency spacing the collinear value comes out "
            f"{t5.get('under_resolved_contour_value', 0):.1e} instead of "
            f"{t5.get('converged_value_there', 0):.1e}, four orders wrong, so "
            "eps must sit well above 2pi/span. THE OTHER CORRECTIONS, "
            "QUANTIFIED. Arrivals: the free ones land at chi, 2pi-chi, 2pi+chi "
            f"to {t9.get('worst_free_offset', 0):.1e} with signs + - +, and the "
            "throat adds two-leg arrivals the free ledger does not contain, "
            "checked at the causal onset because R(w) has poles and a throat "
            "arrival rings up rather than pulsing. TAIL: S3 x R is conformally "
            "flat, so the conformal scalar obeys Huygens exactly -- between "
            f"geometric arrivals the free field is {t10.get('free_ratio', 0):.1e} "
            f"of its peak against {t10.get('throat_ratio', 0):.1e} with the "
            f"throat, a factor of {t10.get('amplification', 0):.1e}. Every bit of "
            "tail in this model is the throat's. CAUSTIC: geometric optics gives "
            "1/(4 pi sin chi), divergent at the antipode, where the exact kernel "
            "is finite and LINEAR in omega; in between the exact/WKB ratio is "
            "|sin(omega e)|, a function of omega*e alone, identical across "
            f"carriers to {t11.get('worst_collapse_spread', 0):.1e}, so the "
            "caustic is cut off at e* ~ 1/omega. AND THE ROUND CLOSES BACK ON "
            "PR #258: the DC content of the solved time series is exactly the "
            "static kernel that round did its tomography on, and running that "
            "protocol on numbers from the dynamic solver returns W = "
            f"{t12.get('W_from_the_time_dependent_solve', 0):.6f} against "
            f"-beta = {t12.get('minus_beta', 0):.2f}, error "
            f"{t12.get('W_error', 0):.1e}, through the whole contour integral "
            "with the O(eps^2) contour bias Richardson-extrapolated and both "
            "numbers reported. WHAT IS STILL PUT IN: the background, the mouth "
            "positions, and the boundary data -- four real numbers chosen and "
            "not derived, with PR #249 still the thing that would fix them from "
            "matter. NO BACKREACTION: the stress tensor is computed from the "
            "field and never fed back, which is the next step, and it now has a "
            "concrete object to feed.")
    else:
        verdict_class = "INCONCLUSIVE"
        failed = [t["name"] for t in tests if not t["pass"]]
        verdict = "INCONCLUSIVE. Failed checks: " + ", ".join(failed)

    return {"probe": "two_wave", "tests": tests,
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
    out = here / "runs" / f"{ts}_two_wave_probe"
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
