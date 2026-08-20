"""The signed areal response: does the interference metric deform toward a neck?

> Scope: an INITIAL-DATA constraint solve on a maximal slice, linearized in the
> conformal factor, on the RESOLVED neck of PR #261/#262 with the balls removed.
> The source is the interference part of a LINEAR conformally coupled scalar's
> stress tensor on a FIXED Einstein static universe, from PR #263's quadrature.
> This is a constraint, not an evolution: there is no sound speed and no
> Eddington mode in it, and no fluid response is modelled.

THE QUESTION
────────────
PR #263 asked whether A + B produces a metric response that rescaling A or B
alone cannot, and answered yes -- but in the transverse-traceless sector, and
then PROVED that sector cannot give a geometric verdict: dA/A = -<h_nn>/2 and
that average vanishes identically for any TT field. Tracelessness kills the
isotropic part of <n_i n_j f(k.n)> and transversality kills the rest. Building
more harmonics adds more exact zeros. So the question that was actually asked --

    DOES THE INTERACTION METRIC DEFORM TOWARD A NECK, AWAY FROM ONE, OR
    MERELY OSCILLATE?

-- had to move to the SCALAR sector, and as an INITIAL-DATA problem rather than
an evolution, which is what removes the two things PR #263 avoided that sector
for. On a maximal slice the K terms in the Hamiltonian constraint are quadratic,
so at the order the field's stress enters, dR3 = 16 pi G drho with no time
derivatives at all.

THE ANSWER IS: TOWARD A NECK. BOTH MOUTHS CLOSE.

    dA/A = (-2.06e-03, -1.88e-03)   in units of 2 pi G

Negative at both mouths, in all eight control combinations. And the mechanism
is not the obvious one: the interference energy ALONE would open the mouths --
U(c_j) > 0 at both -- and the throat's own monopole layers overshoot that and
invert it. The neck closes because the throat cannot support the conformal
factor the energy piles up around its mouths.

WHAT IS CHECKED
───────────────
T1  THE CONSTRAINT, AND WHY THE CONFORMAL ANSATZ IS NOT A RESTRICTION. With
    g = psi^4 ghat and psi = 1 + u, the constraint is nabla^2 u + 3u = -2 pi G
    drho and dA/A = 4u. A TT piece contributes nothing to dR3 and a longitudinal
    piece is a diffeomorphism of a constant-curvature background, so the
    conformal factor is the WHOLE of what this equation sees.

T2  THE DEGENERACY, TWICE. nabla^2 u = -3u is -nabla^2 u + u = 4u and 4 =
    (n+1)^2 at n = 1: the operator sits EXACTLY on the ESU's dipole level. Its
    kernel is the four k=1 harmonics x^A. About a mouth those split into ONE
    l=0 member cos chi and THREE l=1 members sin chi, and the probe checks both
    -- because every closed-form check in PR #262 was l=0, which is exactly how
    an l(l+2) error survived that round's entire suite.

T3  THE PROJECTOR IS THE GREEN FUNCTION'S OWN TAIL. L G_perp = -delta +
    (2/pi^2) cos chi, and (2/pi^2) cos chi IS the integral kernel of the
    projector onto span{x^A}, because sum_A x^A y^A = cos chi(x,y) and
    ||cos chi||^2 = pi^2/2. That identity is what makes the solvability
    condition sum_j A_j c_j + sum_j D_j = S_sigma an IDENTITY rather than a
    modelling choice, so it is checked rather than asserted.

T4  *** THE ASSEMBLY AGAINST EXACT SOLVES, IN BOTH SECTORS. *** The 12x12
    matching is checked against one-dimensional boundary-value solves on the
    punctured sphere that share no code with it. Agreement is 4e-10 or better
    at every radius in both sectors. GETTING THERE TOOK ONE CORRECTION WORTH
    RECORDING: the first assembly agreed at 1e-06 and no better, at every
    radius, and that number DID NOT MOVE when the reference's quadrature,
    stencil and BVP tolerances were each tightened by four orders. So it was
    not the reference's floor but a systematic error of the assembly -- the
    two-point stencil for the mouth-sphere radial derivative, whose RELATIVE
    truncation error on a 1/chi field is exactly step^2. A discrepancy that
    refuses to move when you refine the OTHER side is the other side telling
    you it is not the problem.

T5  THE DIPOLE LAYERS ARE REQUIRED -- AND THEN INVISIBLE. A_1 c_1 + A_2 c_2
    sweeps only the plane of the two mouths, so two monopoles can meet at most
    two of the four solvability equations; measured, the monopole-only condition
    fails by 62.5% of the obstruction, so WITHOUT the l=1 layers there is no
    solution at all. THIS WAS EXPECTED TO BE THE HEART OF THE ROUND AND IT IS
    NOT: the response to the off-plane obstruction is 6e-17, which is zero. The
    layers deposit it in the kernel elements x^2 and x^3, which vanish at both
    mouths because both mouths lie in the (x^0, x^1) plane. A first draft of the
    module docstring said the l=1 sector IS the calculation. It is required for
    the calculation to EXIST and contributes nothing to its VALUE.

T6  *** WHAT THE ANSWER IS ACTUALLY MADE OF, AND WHY THAT IS LUCKY. *** dA/A is
    to 0.09% a linear functional of the obstruction S_sigma alone -- of its
    in-plane part alone, exactly. Deleting the obstruction and keeping the local
    source data leaves a response a THOUSAND times smaller. Scaling the l=1
    source moments by three, by zero, or replacing them with noise moves the
    answer by 5e-04. That matters because those moments are the worst-converged
    input the source quadrature produces: 41% drift between levels, against 1.5%
    for the obstruction. The signed answer rests on the single best-converged
    number available, which is not how these rounds usually go.

T7  THE CONTROLS. Two quadrature levels, two mouth radii, two gluings, plus a
    full 2 pi twist of the transverse frames. The sign is the same in all of
    them. The magnitude is stable to 2.2% across quadrature levels at fixed
    radius. THE MOUTH RADIUS IS NOT A REGULATOR: the source goes as 1/chi^4 at
    the mouths, a is a parameter of the throat, and there is no a -> 0 limit to
    converge to -- that is the singular point PR #262 removed. Doubling a moves
    the answer by a factor 1.75, almost exactly the factor by which the source's
    own mouth-weighted moment moves.

T8  THE SIGN IS A STATEMENT ABOUT A THROAT. The tube's l=0 constraint channel is
    d_s^2 + 4 pi / area -- a CAVITY. At kL = n pi the response has a pole and the
    sign flips; the scan finds flips at 3.133 and 6.260 against closed-form poles
    at pi and 2 pi. The working throat is at kL = 0.9, inside the first cell.
    Beyond the first pole the mouths can even move in OPPOSITE directions.

WHAT IS STILL PUT IN
────────────────────
The gluing of the two mouths' transverse frames through the tube is a modelling
choice -- but a measurably harmless one: a full twist moves dA/A by less than
1e-12, because the areal response sees the dipoles only through sum_j D_j, which
the solvability rows pin. The source's channel data come from the local expansion
of U about each mouth centre, exact through O(a^2) in l=0 and leading order in
l=1 -- and T6 bounds that cost from above by showing the answer barely depends on
those data at all. The fluid holding the ESU static is held RIGID; that is
consistent because the scalar's stress tensor is separately conserved, but a
responsive fluid is the obvious next refinement. And the source is PR #263's,
computed on a FIXED background with POINT sources.

    python -m experiments.closure_ledger.areal_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.areal import (
    INTERFERENCE_MOMENTS, MOUTHS, WORKING_TUBE,
    measure_the_dipole_layers_are_required_not_optional,
    measure_the_kernel_projector_is_the_green_functions_own_tail,
    measure_the_matching_reproduces_an_exact_radial_solve,
    measure_the_obstruction_carries_the_answer,
    measure_the_signed_areal_response,
    measure_the_throat_is_a_resonant_cavity,
)
from geometrodynamics.waves.initial_data import (
    CONSTRAINT_EIGENVALUE, areal_response_from_conformal_factor,
    constraint_operator_eigenvalue,
    measure_the_constraint_operator_is_degenerate_on_the_sphere,
    measure_removing_the_balls_lifts_the_degeneracy,
)


def run_probe() -> dict:
    checks: List[dict] = []

    eig = constraint_operator_eigenvalue()
    checks.append({
        "id": "T1", "name": "the constraint and the areal response",
        "detail": {
            "operator": "nabla^2 u + 3u = -2 pi G drho",
            "eigenvalue": eig,
            "areal_response_is_four_u": float(
                areal_response_from_conformal_factor(0.25)[()]
                if np.ndim(areal_response_from_conformal_factor(0.25)) == 0
                else areal_response_from_conformal_factor(0.25)),
        },
        "pass": bool(abs(eig - CONSTRAINT_EIGENVALUE) < 1e-12),
    })

    deg = measure_the_constraint_operator_is_degenerate_on_the_sphere()
    lift = measure_removing_the_balls_lifts_the_degeneracy()
    checks.append({
        "id": "T2", "name": "the degeneracy, and both sectors of it",
        "detail": {"on_the_sphere": deg, "lifted_by_the_balls": {
            k: v for k, v in lift.items() if not isinstance(v, list)}},
        "pass": bool(deg.get("it_is_degenerate", True)
                     and lift.get("the_dipole_partners_are_also_degenerate",
                                  True)),
    })

    proj = measure_the_kernel_projector_is_the_green_functions_own_tail()
    checks.append({
        "id": "T3", "name": "the projector is the Green function's own tail",
        "detail": proj, "pass": bool(proj["it_is_the_projector"]),
    })

    ref = measure_the_matching_reproduces_an_exact_radial_solve()
    checks.append({
        "id": "T4", "name": "the assembly against exact solves, both sectors",
        "detail": ref, "pass": bool(ref["both_sectors_agree"]),
    })

    dip = measure_the_dipole_layers_are_required_not_optional()
    checks.append({
        "id": "T5", "name": "the dipole layers: required, then invisible",
        "detail": dip,
        "pass": bool(dip["monopoles_alone_cannot_close_it"]
                     and dip["the_dipoles_close_it"]
                     and dip["and_then_they_do_not_move_the_answer"]),
    })

    obs = measure_the_obstruction_carries_the_answer()
    checks.append({
        "id": "T6", "name": "what the answer is made of",
        "detail": obs, "pass": bool(obs["the_obstruction_is_the_answer"]),
    })

    head = measure_the_signed_areal_response()
    checks.append({
        "id": "T7", "name": "*** the signed areal response ***",
        "detail": head, "pass": bool(head["every_variant_agrees_in_sign"]),
    })

    cav = measure_the_throat_is_a_resonant_cavity()
    checks.append({
        "id": "T8", "name": "the sign is a statement about a throat",
        "detail": cav,
        "pass": bool(cav["flips_land_on_the_poles"]
                     and cav["the_working_throat_is_off_resonance"]),
    })

    return {
        "probe": "areal",
        "question": "does the interaction metric deform toward a neck, away "
                    "from one, or merely oscillate?",
        "answer": "toward a neck: dA/A < 0 at both mouths",
        "areal_response": head["areal_response"],
        "units": "2 pi G, with the source being PR #263's interference dT00",
        "mouths": [list(c) for c in MOUTHS],
        "tube": {"area": WORKING_TUBE.area, "length": WORKING_TUBE.length,
                 "wavenumber": WORKING_TUBE.wavenumber()},
        "moment_levels": [{"radius": m.radius, "points": m.points}
                          for m in INTERFERENCE_MOMENTS],
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    lines = [
        "# Areal probe — the signed ΔA/A on the resolved neck",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']} — "
        f"ΔA/A = ({summary['areal_response'][0]:+.4e}, "
        f"{summary['areal_response'][1]:+.4e}) in units of {summary['units']}.",
        "",
        f"Throat: area {summary['tube']['area']:.4f}, length "
        f"{summary['tube']['length']}, wavenumber "
        f"{summary['tube']['wavenumber']:.4f} — phase "
        f"{summary['tube']['wavenumber'] * summary['tube']['length']:.2f}, "
        "inside the first cavity cell.",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")
    head = next(c for c in summary["checks"] if c["id"] == "T7")["detail"]
    lines += [
        "",
        "## The controls",
        "",
        "| radius | points | gluing | ΔA/A mouth 1 | ΔA/A mouth 2 |",
        "|--------|--------|--------|--------------|--------------|",
    ]
    for r in head["rows"]:
        lines.append(f"| {r['radius']:.2f} | {r['points']} | {r['gluing']} | "
                     f"{r['areal_response'][0]:+.4e} | "
                     f"{r['areal_response'][1]:+.4e} |")
    lines += [
        "",
        f"Sign agrees in every variant: {head['every_variant_agrees_in_sign']}. "
        f"Quadrature spread at fixed radius: "
        f"{head['quadrature_spread_at_fixed_radius']:.2%}. "
        f"Worst condition number: {head['worst_condition_number']:.2e}, "
        f"worst residual: {head['worst_residual']:.1e}.",
        "",
        f"**What it means.** {head['the_answer']}",
        "",
    ]
    return "\n".join(lines)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_areal_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
