"""Which throat is physical -- and the sign reverses on the one that is forced.

> Scope: the same INITIAL-DATA constraint solve as PR #264, on the same maximal
> slice with the same PR #263 source.  What changes is the throat: instead of a
> product tube whose area and length were free parameters, the throat that the
> constraint and a smooth gluing FORCE.  Still LINEAR in the conformal factor,
> still a FIXED background for the field, still a RIGID fluid.

THE QUESTION
────────────
PR #264 measured dA/A on a product tube of area 4 pi and length 0.9 and found
the mouths CLOSE. Then, asked to match the tube's area to its own mouths', it
found the sign REVERSES. So the throat's geometry stopped being decoration:

    A and L were free parameters. WHICH VALUES ARE PHYSICAL?

THEY WERE NEVER FREE
────────────────────
On a maximal slice the BACKGROUND constraint is R_hat = 16 pi G rho, so a
profile does not choose its matter -- the matter is whatever the profile
implies. A product tube of area A is S^2(r) x R with R_hat = 8 pi / A, so

    A = 4 pi         (PR #261-#264)   ->  rho_tube = rho_ambient / 3
    A = 4 pi sin^2 a (the matched one) ->  rho_tube = 133 rho_ambient

Neither was chosen for a reason, and neither is the ambient's own fluid.

THE ANSWER: THE VACUUM THROAT, AND IT HAS NO FREE PARAMETER AT ALL.

Ask for a throat that needs NO matter (R_hat = 0) and glues on with NO surface
layer. R_hat = 0 integrates once to f'^2 = 1 + C/f; a neck fixes C = -f0; and
smooth gluing to the round S^3 at mouth radius a needs f' = cos a where
f = sin a, which forces f0 = sin^3 a. Both conditions used, nothing left over.
Length and resistance then follow in closed form, checked against quadrature:

    L = 2[sin^3 a . arccosh(1/sin a) + sin a cos a]   ~  2a
    I = int ds/f^2 = 4 cos a / sin^3 a

and the conductance is EXACTLY a quarter of the exterior's own monopole
stiffness, 4 pi / I = N0(a,4)/4, at every mouth radius.

AND THE THROAT HAS A NAME. R_hat = 0, K = 0 and a spherical neck do not merely
permit an Einstein-Rosen bridge -- they ARE one: f'^2 = 1 - f0/f is exactly the
time-symmetric Schwarzschild slice with f0 = 2M. So the forced neck radius is
twice a mass, and

    M = sin^3(a) / 2

is the throat's MASS DERIVED FROM THE SIZE OF THE EXCISED MOUTH, with nothing
left to choose. Three quasi-local masses agree on it exactly, and the gluing
condition turns out to BE the continuity of the Hawking mass across the seam.

AND THE SIGN REVERSES.  dA/A = (+6.64, +8.58) in units of 2 pi G: on the throat
that is forced rather than chosen, the interference OPENS both mouths.

WHAT IS CHECKED
───────────────
T1  THE CURVATURE, DERIVED FROM THE METRIC. R_hat = 2/f^2 - 4f''/f - 2f'^2/f^2,
    obtained by computing Christoffels and contracting rather than by quoting,
    and checked against two cases whose answers are known independently: the
    round S^3 (f = sin s) must give exactly 6, and the vacuum profile must give
    exactly 0. Both hold to machine precision.

T2  THE GLUING FORCES EVERYTHING. f0 = sin^3 a, with the closed forms for L and
    I verified against quadrature to 1e-12 and the conductance ratio to N0
    exactly 1/4 at every radius. THE THROAT IS SHORT: L ~ 2a, not the 0.9 that
    PR #261-#264 carried as a free parameter.

T3  WHAT THE OTHER THROATS WOULD NEED. rho/3 and 133 rho. And the one product
    area that DOES carry the ambient fluid, A = 4 pi/3, cannot be glued without
    a surface layer -- two regions of equal constant R_hat joined smoothly are
    one region of that R_hat, and the only rotationally symmetric one is the
    round sphere, which has no neck.

T4  *** THERE IS NO CAVITY. *** The constraint operator is nabla^2 + R_hat/2, so
    R_hat = 0 leaves the PLAIN LAPLACIAN: (f^2 u')' = 0 in l=0, solutions
    A + B int ds/f^2, MONOTONE. No standing wave, no poles, nothing for a sign
    to flip across. PR #264's cavity, its resonances at kL = n pi and its sign
    flips were properties of MATTER IN THE TUBE, not of a throat.

T5  *** THE MECHANISM: ONE NUMBER DECIDES THE SIGN. *** Decompose a symmetric
    two-port as Y = G [[-1,1],[1,-1]] + shunt I -- a conductance THROUGH the
    tube and a shunt INTO it. The shunt is the flux a UNIFORM potential drives
    in, and for a vacuum tube (f^2 u')' = 0 makes it vanish IDENTICALLY: there
    is nowhere to put that flux. Scanned separately, the conductance moves the
    answer over EIGHT ORDERS and never changes its sign; the shunt passes
    through a pole near 2e-03 and flips it. PR #264's tube sat at 6.07.

T6  *** THE THROAT IS AN EINSTEIN-ROSEN BRIDGE, AND ITS MASS IS DERIVED. ***
    f0 = 2M, so M = sin^3(a)/2. The Schwarzschild parameter, the irreducible
    mass sqrt(A/16pi) and the Hawking mass agree to 1e-13, the neck area is
    16 pi M^2 exactly, and the gluing condition IS Hawking-mass continuity --
    (f/2)(1 - f'^2) is f0/2 on the throat and sin^3(chi)/2 on the ambient.
    FOUR THINGS IT DOES NOT SAY, asserted in the tests because the claim is
    strong enough to be worth not overstating: there is no asymptotic region,
    so no ADM mass -- the derived mass is quasi-local, unambiguous only because
    the Hawking mass is constant on the vacuum piece; it is DIMENSIONLESS, M/R
    against the ESU curvature radius, which is the only kind the scale-modulus
    theorem of #52 allows; both ends sew into the SAME S^3, so it is a handle
    of Misner's kind; and the neck is a minimal surface, which on a K = 0 slice
    is an APPARENT HORIZON. That last is the vacuum face of #7's result that a
    traversable connection needs exotic matter: the throat with none in it is
    the one that is not traversable.

T7  *** THE ANSWER, WITH ITS CONTROLS. *** Both mouths open, in all eight
    combinations of two quadrature levels, two mouth radii and two gluings, and
    across the whole vacuum family -- four orders in the neck radius, with the
    smooth-gluing point in the middle of it, so the answer does not depend on
    hitting the gluing condition exactly.

T8  *** WHERE THE RICCATI SOLVE STOPS RESOLVING. *** (PR #267, release
    hardening.) The static problem on this profile has a CLOSED FORM:
    f = f0 cosh^2 x with ds = 2f dx turns (f^2 u')' = l(l+1)u into
    y'' = (2l+1)^2 y, a constant-coefficient equation, and the half-length
    X = arcosh(1/sin a) has e^-X = tan(a/2) exactly. So with k = 2l+1 and
    q = tan^2k(a/2),
        D_l = -2pi sin(a) [ k(1+q^2)/(1-q^2) - cos a ]
        C_l = +4pi sin(a) k q / (1-q^2)  .
    That is now the PRODUCTION admittance. The Riccati solve is RETAINED, as
    an independent validator, and demoted -- because its last step forms
    Y12 = (s - t)/2 from two eigenchannel values that agree to more digits than
    the solver carries. At a = 0.05 it is right to 1e-13 at l = 0, to four
    figures at l = 1, and at l = 2 it returns -1.17e-14 for a true +3.00e-16:
    wrong sign, and larger than the answer. The DIAGONAL was never affected --
    the two routes agree to 1e-14 in every channel. The floor is around 1e-12
    in |Y12|. This probe does NOT pin the 39x factor: that is one solver's step
    sequence in one build of SciPy. It pins the boundary.
    The l = 0 channel is special-cased to (pi sin^3 a / cos a) [[-1,1],[1,-1]],
    which avoids a DIFFERENT near-cancellation -- k coth(2kX) - tanh X is a
    difference of two quantities both tending to 1 as a -> 0 -- and makes the
    zero-shunt identity Y(1,1)^T = 0 exact to the last bit rather than to a
    tolerance. #265's answer is UNCHANGED: solve_matching uses only l = 0 and
    l = 1, both above the floor, and moves in the 13th digit.

T9  *** THE MOUTH-TO-MOUTH HIERARCHY: C_l ~ a^(4l+3). *** From the closed
    form, C_l -> 4pi(2l+1) sin(a) tan^(4l+2)(a/2) as a -> 0 -- verified to
    2e-16 relative at a = 0.05, and the fitted exponents are 3.000, 7.000,
    11.000, 15.000. Each unit of angular momentum costs FOUR powers of the
    mouth radius. The ESU kernel is the four n = 1 harmonics x^A, degenerate
    on the round S^3; a throat cut at one point splits them by LOCAL angular
    momentum about that point, 1 (X^0 = cos chi) + 3 (X^i = sin chi n^i), and
    the two pieces cross with C0/C1 = 8.5e+05 at a = 0.05.
    THE STATEMENT THIS SUPPORTS, and the only one: the static scalar Laplacian
    on this scalar-flat spatial throat suppresses the local l = 1 mouth-to-
    mouth channel by ~1e-09 at a = 0.05, while preserving a much stronger
    monopole channel. It is NOT a statement that orientation cannot cross the
    throat. One operator, one slice, no lapse chosen, and the l = 1 channel is
    SMALL rather than zero.

WHAT IT COSTS, AND WHY IT IS SAID OUT LOUD
──────────────────────────────────────────
The response is 3000x PR #264's and grows as a^-3, and the matching system's
condition number grows with it. That is not noise -- it is the physics of a
throat with zero shunt BY IDENTITY, which therefore does almost nothing to lift
the constraint's exact degeneracy, so the linear response sits close to a mode
the operator nearly annihilates. THE SIGN IS ROBUST. Whether the linearisation
was legitimate at that amplitude is the question this round leaves open, and it
is now the binding one.

WHAT IS STILL PUT IN
────────────────────
The source is PR #263's, on a FIXED background with POINT sources. The ESU's
fluid is held RIGID. The exterior is the round S^3 with two balls removed, so
the throat is glued to an ambient that does not itself deform -- the C^1 gluing
is exact for the profile but the ambient's own response to the handle is not
modelled. And a vacuum tube is a choice of the SIMPLEST acceptable matter, not
the only one: any rho(s) >= 0 that glues smoothly is admissible, and the vacuum
one is distinguished by needing no matter at all and no free parameter.

    python -m experiments.closure_ledger.physical_throat_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.waves.areal import (WORKING_TUBE,
                                          measure_the_signed_areal_response)
from geometrodynamics.waves.physical_throat import (
    WORKING_VACUUM_THROAT,
    measure_the_curvature_formula_is_derived_not_quoted,
    measure_the_gluing_forces_the_neck_radius,
    measure_the_product_tubes_need_anomalous_matter,
    measure_the_shunt_decides_the_sign,
    measure_the_signed_response_on_the_physical_throat,
    measure_the_mouth_to_mouth_hierarchy,
    measure_the_throat_is_an_einstein_rosen_bridge,
    measure_the_vacuum_throat_has_no_cavity,
    measure_where_the_riccati_solve_stops_resolving,
)


def run_probe() -> dict:
    checks: List[dict] = []

    curv = measure_the_curvature_formula_is_derived_not_quoted()
    checks.append({"id": "T1", "name": "the curvature, derived from the metric",
                   "detail": curv, "pass": bool(curv["both_are_exact"])})

    glue = measure_the_gluing_forces_the_neck_radius()
    checks.append({
        "id": "T2", "name": "the gluing forces the neck radius, L and I",
        "detail": glue,
        "pass": bool(glue["there_is_no_free_parameter"]
                     and glue["the_conductance_is_a_quarter_of_the_exteriors"]
                     and glue["worst_length_drift"] < 1e-9)})

    matter = measure_the_product_tubes_need_anomalous_matter()
    checks.append({
        "id": "T3", "name": "what the other throats would need",
        "detail": matter,
        "pass": bool(matter["neither_used_area_is_the_ambient_fluid"])})

    cav = measure_the_vacuum_throat_has_no_cavity()
    checks.append({
        "id": "T4", "name": "*** there is no cavity ***",
        "detail": cav,
        "pass": bool(cav["it_matches_the_closed_form"]
                     and cav["the_symmetric_channel_is_exactly_dead"]
                     and not cav["vacuum_throat_has_resonances"])})

    mech = measure_the_shunt_decides_the_sign()
    checks.append({
        "id": "T5", "name": "*** the shunt decides the sign ***",
        "detail": mech,
        "pass": bool(mech["conductance_never_changes_the_sign"]
                     and mech["the_shunt_does"]
                     and mech["the_shunt_is_the_tubes_matter"])})

    bridge = measure_the_throat_is_an_einstein_rosen_bridge()
    checks.append({
        "id": "T6",
        "name": "*** it is an Einstein-Rosen bridge, and M is derived ***",
        "detail": bridge,
        "pass": bool(bridge["it_is_an_einstein_rosen_bridge"]
                     and bridge["three_masses_agree"] < 1e-12
                     and bridge["the_gluing_is_hawking_mass_continuity"] < 1e-12)})

    head = measure_the_signed_response_on_the_physical_throat()
    checks.append({
        "id": "T7", "name": "*** the answer, and it reverses ***",
        "detail": head,
        "pass": bool(head["every_variant_agrees_in_sign"]
                     and head["sign"] == ["opens", "opens"])})

    floor = measure_where_the_riccati_solve_stops_resolving()
    checks.append({
        "id": "T8",
        "name": "*** where the Riccati solve stops resolving (PR #267) ***",
        "detail": floor,
        "pass": bool(floor["the_cross_term_fails_at_ell_two"]
                     and floor["the_diagonal_was_never_affected"]
                     and floor["riccati_is_trustworthy_through_ell"] == 1)})

    hier = measure_the_mouth_to_mouth_hierarchy()
    checks.append({
        "id": "T9",
        "name": "*** the mouth-to-mouth hierarchy C_l ~ a^(4l+3) ***",
        "detail": hier,
        "pass": bool(hier["the_exponent_law_holds"]
                     and hier["the_asymptotic_is_the_leading_term"])})

    wide = measure_the_signed_areal_response()
    return {
        "probe": "physical_throat",
        "question": "the throat's area and length were free parameters. "
                    "which values are physical?",
        "answer": "the vacuum throat, glued with no surface layer -- an "
                  "Einstein-Rosen bridge with no free parameter at all: "
                  "f0 = sin^3 a is forced, and f0 = 2M, so M = sin^3(a)/2",
        "mass_law": bridge["mass_law"],
        "mass": WORKING_VACUUM_THROAT.mass(),
        "areal_response": head["areal_response"],
        "previous_round": wide["areal_response"],
        "sign": head["sign"],
        "previous_sign": wide["sign"],
        "units": "2 pi G, with the source being PR #263's interference dT00",
        "throat": {"mouth_radius": WORKING_VACUUM_THROAT.mouth_radius,
                   "neck_radius": WORKING_VACUUM_THROAT.neck_radius(),
                   "length": WORKING_VACUUM_THROAT.length(),
                   "resistance": WORKING_VACUUM_THROAT.resistance(),
                   "shunt": WORKING_VACUUM_THROAT.shunt(),
                   "mass": WORKING_VACUUM_THROAT.mass()},
        "previous_throat": {"area": WORKING_TUBE.area,
                            "length": WORKING_TUBE.length,
                            "shunt": WORKING_TUBE.shunt()},
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    t, p = summary["throat"], summary["previous_throat"]
    lines = [
        "# Physical-throat probe — which throat, and the sign it gives",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']}",
        "",
        "| | this round (vacuum throat) | PR #264 (product tube) |",
        "|--|--|--|",
        f"| neck / area | `f₀ = {t['neck_radius']:.5e}` | `𝒜 = {p['area']:.4f}` |",
        f"| length | `{t['length']:.5f}` (forced) | `{p['length']}` (chosen) |",
        f"| shunt | `{t['shunt']:.1f}` — zero by identity | `{p['shunt']:.4f}` |",
        f"| mass | `M = {t['mass']:.5e}` — **derived** | not defined |",
        f"| `ΔA/A` | `({summary['areal_response'][0]:+.4e}, "
        f"{summary['areal_response'][1]:+.4e})` | "
        f"`({summary['previous_round'][0]:+.4e}, "
        f"{summary['previous_round'][1]:+.4e})` |",
        f"| verdict | **{' / '.join(summary['sign'])}** | "
        f"{' / '.join(summary['previous_sign'])} |",
        "",
        f"In units of {summary['units']}.  Mass law: `{summary['mass_law']}`.",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")
    mech = next(c for c in summary["checks"] if c["id"] == "T5")["detail"]
    lines += ["", "## The mechanism", "",
              "| conductance | shunt | `ΔA/A` mouth 1 | verdict |",
              "|--|--|--|--|"]
    for c in mech["corners"]:
        lines.append(f"| `{c['conductance']:.4e}` | `{c['shunt']:.4f}` | "
                     f"`{c['areal_response'][0]:+.5e}` | {'/'.join(c['sign'])} |")
    lines += [
        "",
        "The conductance is scanned over eight orders and never changes the "
        "sign. The shunt passes through a pole near `2e-03` and flips it. "
        "A vacuum tube has zero shunt **by identity** — `(f²u')' = 0` — so "
        "there is nowhere for monopole flux to go.",
        "",
        f"**What it costs.** {next(c for c in summary['checks'] if c['id'] == 'T7')['detail']['what_it_costs']}.",
        "",
    ]

    floor = next(c for c in summary["checks"] if c["id"] == "T8")["detail"]
    lines += [
        "## The two-port in closed form, and where the solve stopped resolving",
        "",
        "`f = f₀cosh²x` with `ds = 2f dx` turns `(f²u')' = ℓ(ℓ+1)u` into "
        "`y'' = (2ℓ+1)²y`, and `e^{−X} = tan(a/2)` exactly.  With `k = 2ℓ+1` "
        "and `q = tan^{2k}(a/2)`:",
        "",
        "> `D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]`  ,  "
        "`C_ℓ = +4π sin a · kq/(1−q²)`",
        "",
        f"That is now the production admittance at `a = {floor['mouth_radius']}`.  "
        "The Riccati solve is kept as an independent validator, not deleted.",
        "",
        "| `ℓ` | `C_ℓ` closed form | `C_ℓ` Riccati | rel. error | signs | diagonal rel. error |",
        "|--|--|--|--|--|--|",
    ]
    for r in floor["rows"]:
        lines.append(
            f"| {r['ell']} | `{r['cross_closed_form']:+.5e}` | "
            f"`{r['cross_riccati']:+.5e}` | `{r['cross_relative_error']:.2e}` | "
            f"{'agree' if r['signs_agree'] else '**opposite**'} | "
            f"`{r['diagonal_relative_error']:.1e}` |")
    lines += [
        "",
        "The cross term is a *difference* of two eigenchannel values in the "
        "solve and a *product* of small factors in the closed form.  The "
        f"floor is around `{floor['the_floor_is_about']:.0e}`; the diagonal "
        "was never affected.",
        "",
        "The `39×` factor is deliberately **not** pinned: it is one solver's "
        "step sequence in one build of SciPy and would move under any of "
        "them.  What is pinned is the boundary — that the sign is wrong and "
        "the magnitude is more than an order of magnitude out, in a channel "
        "whose honest size is below the solver's floor.",
        "",
    ]

    hier = next(c for c in summary["checks"] if c["id"] == "T9")["detail"]
    lines += [
        "## The hierarchy",
        "",
        f"`{hier['law']}` — fitted exponents "
        + ", ".join(f"`{f['fitted_exponent']:.6f}`" for f in hier["fits"])
        + " against `3, 7, 11, 15`.",
        "",
        "| `a` | `C₀` | `C₁` | `C₀/C₁` | `C₁/D₁` |",
        "|--|--|--|--|--|",
    ]
    for r in hier["rows"]:
        lines.append(f"| `{r['mouth_radius']:.2f}` | `{r['C0']:.6e}` | "
                     f"`{r['C1']:.6e}` | `{r['C0_over_C1']:.4e}` | "
                     f"`{r['ell_one_transmission']:.3e}` |")
    lines += [
        "",
        f"{hier['kernel_split']}, so the ESU kernel splits `1 ⊕ 3` and the "
        "two pieces cross four powers of `a` apart.",
        "",
        f"**The statement this supports:** {hier['the_narrow_statement']}.",
        "",
        f"**What is not claimed:** {hier['what_is_not_claimed']}.",
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
    out = here / "runs" / f"{ts}_physical_throat_probe"
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
