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
    measure_the_throat_is_an_einstein_rosen_bridge,
    measure_the_vacuum_throat_has_no_cavity,
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
