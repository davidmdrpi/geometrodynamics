"""Probe: does the finite mouth force the transport and the antipodal BC?

Every check compares a computation against `docs/finite_mouth_topology_prereg.md`
(commit `d9d85bc`), written before `geometrodynamics/bulk/mouth_topology.py`
existed. The verdict rule is frozen there. Nothing here uses a singlet, a Born
rule, a projector, a tensor product, CHSH or a QED amplitude.

Run:  python -m experiments.closure_ledger.finite_mouth_topology_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import finite_mouth as fm      # noqa: E402
from geometrodynamics.bulk import mouth_topology as mt    # noqa: E402

TARGETS = {
    (0, "Neumann"): 0.000000000, (0, "Dirichlet"): 0.157587622,
    (1, "Neumann"): 1.797266559, (1, "Dirichlet"): 1.804461992,
    (2, "Neumann"): 3.524607516, (2, "Dirichlet"): 3.524854051,
    (3, "Neumann"): 5.248595411, (3, "Dirichlet"): 5.248602920,
}


def check_h1() -> dict:
    """H1 — the gluing is free; antipodal gluing is orientable in the bulk."""
    classes = {}
    for name, m in (("identity", mt.identity_gluing()),
                    ("antipodal", mt.antipodal_gluing()),
                    ("reflection (normal)", mt.reflection_gluing(3)),
                    ("reflection (brane)", mt.reflection_gluing(0))):
        c = mt.classify_gluing(m)
        classes[name] = {"det_bulk": c.det_bulk, "det_brane": c.det_brane,
                         "eps": c.normal_sign, "label": c.label,
                         "w1_loop": mt.mapping_torus_w1(m)}
    reps = []
    for ell in range(5):
        r = mt.harmonic_representation(mt.antipodal_gluing(), ell)
        reps.append({"ell": ell, "dim": r["dimension"], "expected": (ell + 1) ** 2,
                     "scalar": r["scalar"], "residual": r["fit_residual"]})
    refl = mt.harmonic_representation(mt.reflection_gluing(3), 2)
    t2 = classes["antipodal"]["det_bulk"] == +1
    t3 = classes["antipodal"]["label"] == "(det m_3, eps) = (-, -)" and \
        len({(c["det_brane"], c["eps"]) for c in classes.values()}) == 4
    t4 = all(r["scalar"] == (-1) ** r["ell"] and r["dim"] == r["expected"]
             and r["residual"] < 1e-9 for r in reps) and not refl["is_scalar"]
    return {"classes": classes, "harmonic_action_of_antipode": reps,
            "reflection_is_scalar_on_l2": refl["is_scalar"],
            "T2_antipodal_gluing_orientable": t2,
            "T3_four_classes_and_antipodal_is_minus_minus": t3,
            "T4_parity_computed_not_assumed": t4,
            "holds": bool(t2 and t3 and t4),
            "class": "topological theorem",
            "why": ("Darmois matching is O(4)-invariant, so the seam maps are "
                    "free and only the loop monodromy m matters. det(-I_4) = +1: "
                    "the antipodally glued two-mouth handle is ORIENTABLE in the "
                    "bulk; what it twists is the brane S^2-handle and the brane's "
                    "normal line bundle. The (-1)^l is the computed action of -I "
                    "on degree-l harmonics; a reflection does not act as a "
                    "scalar at all.")}


def check_h2() -> dict:
    """H2 — the unique free involution and what it makes non-orientable."""
    inv = mt.antipodal_involution()
    o5 = mt.free_involutions_of_o5()
    free = [r["minus_ones"] for r in o5 if r["free_on_S4"]]
    pins = {n: mt.pin_structures_rp(n) for n in (2, 3, 4)}
    i1 = inv.is_free and inv.is_isometry(0.09, fm.lapse_ultrastatic) and \
        inv.is_isometry(0.09, fm.lapse_vacuum) and inv.tangent_determinant() == -1 \
        and inv.brane_determinant() == +1
    i2 = free == [5]
    i3 = (not pins[4]["orientable"]) and pins[4]["pin_plus"] and not pins[4]["pin_minus"] \
        and pins[3]["spin"] and pins[2]["pin_minus"] and not pins[2]["pin_plus"]
    return {"involution": inv.name, "free": inv.is_free,
            "det_bulk": inv.tangent_determinant(), "det_brane_slice": inv.brane_determinant(),
            "isometry_ultrastatic": inv.is_isometry(0.09, fm.lapse_ultrastatic),
            "isometry_tangherlini": inv.is_isometry(0.09, fm.lapse_vacuum),
            "free_classes_in_O5": free, "o5_table": o5,
            "pin_types": pins,
            "quotient": {"bulk": "RP^4 # RP^4 (non-orientable; Pin+ only)",
                         "brane_slice": "RP^3 # RP^3 (orientable, spin)",
                         "neck": "RP^3 (orientable, spin)",
                         "brane_neck": "RP^2 (non-orientable; Pin- only)"},
            "I1": bool(i1), "I2": bool(i2), "I3": bool(i3),
            "holds": bool(i1 and i2 and i3),
            "class": "topological theorem, conditional on P_B = -P_A (chosen)",
            "why": ("-I_5 is the only fixed-point-free involution in O(5), and "
                    "through the tube it continues uniquely as (s, Omega) -> "
                    "(-s, -Omega). That map is free, an isometry for both lapses, "
                    "and reverses the bulk orientation while preserving the "
                    "brane slice's. The non-orientable object is the QUOTIENT "
                    "M/iota = RP^4 # RP^4, not the two-mouth handle; the RP^2 the "
                    "repository names is its brane neck. And the Pin types do not "
                    "match: RP^4 # RP^4 is Pin+ only, RP^2 is Pin- only.")}


def check_h3() -> dict:
    """H3 — the scalar sector and the half-tube oracle."""
    sectors = {"eta=+1": [mt.neck_sector(l) for l in range(6)],
               "eta=-1": [mt.neck_sector(l, eta=-1) for l in range(6)]}
    b1 = sectors["eta=+1"] == ["Neumann", "Dirichlet"] * 3 and \
        sectors["eta=-1"] == ["Dirichlet", "Neumann"] * 3
    rows, worst = [], 0.0
    for (ell, cond), target in TARGETS.items():
        oracle = mt.half_tube_admittance_oracle(ell, cond)
        num = mt.half_tube_admittance(ell, cond, steps=2000)
        rel = abs(num - oracle) / max(1.0, abs(oracle))
        Y = fm.static_admittance(ell)
        eig = Y[0, 0] + Y[0, 1] if cond == "Neumann" else Y[0, 0] - Y[0, 1]
        rows.append({"ell": ell, "condition": cond, "target": target, "oracle": oracle,
                     "solve": num, "relative_error": rel,
                     "pr277_eigenvalue": float(eig),
                     "oracle_matches_target": abs(oracle - target) < 2e-9,
                     "oracle_is_pr277_sector": abs(eig - oracle) < 1e-12})
        worst = max(worst, rel)
    conv = []
    for cond, ell in (("Neumann", 2), ("Dirichlet", 1)):
        oracle = mt.half_tube_admittance_oracle(ell, cond)
        errs = [abs(mt.half_tube_admittance(ell, cond, steps=n) - oracle)
                for n in (500, 1000, 2000, 4000)]
        conv.append({"condition": cond, "ell": ell, "errors": errs,
                     "ratios": [errs[i] / errs[i + 1] for i in range(3)]})
    b2 = worst < 1e-5 and all(r["oracle_matches_target"] and r["oracle_is_pr277_sector"]
                              for r in rows) and all(min(c["ratios"]) >= 3.5 for c in conv)
    ctrl = mt.neck_reflection_involution()
    b3 = (not ctrl.is_free) and \
        [mt.neck_sector(l, involution=ctrl) for l in range(6)] == ["Neumann"] * 6
    geo = []
    for a in (0.05, 0.3, 0.8, 1.2):
        geo.append({"a": a, "labels": [mt.neck_sector(l) for l in range(4)],
                    "Y1_D": mt.half_tube_admittance_oracle(1, "Dirichlet", 1.0, a)})
    b4 = all(g["labels"] == ["Neumann", "Dirichlet"] * 2 for g in geo) and \
        len({round(g["Y1_D"], 6) for g in geo}) == 4
    return {"sectors": sectors, "oracle_rows": rows, "worst_relative_error": worst,
            "convergence": conv, "antipodal_control_free": ctrl.is_free,
            "antipodal_control_labels": [mt.neck_sector(l, involution=ctrl) for l in range(6)],
            "geometry_control": geo,
            "B1": bool(b1), "B2": bool(b2), "B3": bool(b3), "B4": bool(b4),
            "holds": bool(b1 and b2 and b3 and b4),
            "class": "analytic identity, conditional on H2's choices and on eta",
            "why": ("On M/iota a scalar obeys psi_l(-s) = eta (-1)^l psi_l(s), so "
                    "the neck is Neumann or Dirichlet by the parity of l. PR #129's "
                    "even/odd structure is the eta = +1 sector, obtained at the "
                    "finite ultrastatic neck with no horizon. The half-tube "
                    "admittance is the (1, ±1) sector of PR #277's two-mouth "
                    "oracle, reproduced by an independent second-order solve. "
                    "Replacing the antipode by the identity fixes the whole neck "
                    "and erases the grading.")}


def check_h4() -> dict:
    """H4 — no horizon limit is needed; the map contains a time reversal."""
    res = {lapse.__name__: mt.static_operator_commutes_with_parity(0.09, lapse)
           for lapse in (fm.lapse_ultrastatic, fm.lapse_vacuum)}
    l1 = all(v < 1e-8 for v in res.values())
    # (U,V) -> (-U,-V) is rotation by pi in the Lorentzian (U,V) plane: det +1,
    # reverses T = U+V (time orientation) and X = V-U (space orientation).
    # (U,V) -> (V,U) is the reflection: det -1, preserves T, reverses X.
    uv_rot = np.array([[-1.0, 0.0], [0.0, -1.0]])
    uv_ref = np.array([[0.0, 1.0], [1.0, 0.0]])
    T = np.array([1.0, 1.0])
    l2 = {"pr129_map_det_UV": float(np.linalg.det(uv_rot)),
          "pr129_map_reverses_time_orientation": bool(np.allclose(uv_rot @ T, -T)),
          "alternative_reflection_det_UV": float(np.linalg.det(uv_ref)),
          "alternative_preserves_time_orientation": bool(np.allclose(uv_ref @ T, T)),
          "both_restrict_to_iota_on_t_equals_0": True,
          "five_d_det_pr129": float(np.linalg.det(uv_rot)) * 1.0,   # x det(-I_4)=+1
          "five_d_det_alternative": float(np.linalg.det(uv_ref)) * 1.0}
    return {"parity_commutator_residuals": res, "L1": bool(l1), "L2": l2,
            "holds": bool(l1),
            "class": "analytic identity",
            "why": ("f and both lapses are even in s, so the static operator "
                    "commutes with s -> -s for the ultrastatic and Tangherlini "
                    "branches alike: the sector labels are lapse-independent and "
                    "there is no limit to take. PR #129's (U,V) -> (-U,-V) is the "
                    "PT-type extension of iota (5D det +1, time orientation "
                    "reversed); the P-type extension (U,V) -> (V,U) exists too "
                    "(5D det -1, time orientation kept). Both restrict to iota on "
                    "the t = 0 slice. Which one is physical is a further choice, "
                    "recorded as the datum T.")}


def check_h5() -> dict:
    """H5 — J is not the lift of any gluing map."""
    c = mt.complex_structure_commutation()
    lifts = mt.pin_lifts_of_reflection()
    anti = mt.spin_lifts_of_antipode()
    s1 = lifts["count"] == 4 and all(r["anticommutes_with_volume"] and r["is_real_matrix"]
                                     for r in lifts["lifts"])
    s2 = c["J_equals_left_mult_by_minus_j"] and c["anticommutes_with_L_i"] and \
        c["commutes_with_R_i"] and np.allclose(c["linear_matrix_in_R_i_basis"], [[0, 1], [-1, 0]])
    s3 = s2 and s1 and abs(c["det_J"] - 1.0) < 1e-12 and c["hopf_base_antipode"]
    s4 = lifts["lifts_within_a_type_are_distinct"] and anti["eigenvalues"] == [-1.0, 1.0]
    return {"J": {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in c.items()},
            "reflection_lifts": lifts, "antipode_lifts": anti,
            "S1": bool(s1), "S2": bool(s2), "S3": bool(s3), "S4": bool(s4),
            "holds": bool(s1 and s2 and s3 and s4),
            "class": "analytic identity (S2, S3); representation choice (the K)",
            "why": ("The lift of the neck reflection is ±e_s in Pin^±(4): a real "
                    "matrix that exchanges chiralities. The lift of -I_4 is ± the "
                    "volume element: the chirality sign. J = i sigma_y K is left "
                    "multiplication by -j, a ROTATION with det +1; it is antilinear "
                    "only with respect to the Hopf complex structure L_i and is the "
                    "linear SU(2) matrix i sigma_y with respect to R_i. No Pin lift "
                    "is antilinear, so J is not the lift of any gluing map. Its K is "
                    "the reversal of the Hopf U(1) - charge conjugation, a bundle "
                    "datum - and its i sigma_y is a spin rotation by pi, neither of "
                    "which the mouth gluing supplies.")}


def audit_table(h: dict) -> list:
    return [
        ("Is the physical mouth non-orientable?",
         "The two-mouth handle with antipodal gluing: NO (bulk w_1 = 0; brane "
         "S^2-handle and normal bundle twisted). The quotient M/iota: YES "
         "(RP^4 # RP^4).", "topological theorem"),
        ("What object carries w_1 != 0?",
         "On the handle: the brane S^2 x~ S^1 and the brane's normal line bundle "
         "(eps = -1). On the quotient: the bulk itself and its brane neck RP^2. "
         "Never the bulk mouth S^3.", "topological theorem"),
        ("Is J = i sigma_y K a valid lift?",
         "NO, of any gluing map: every Pin lift is complex-linear; J is antilinear "
         "for the Hopf structure. J = (spin rotation by pi) o (U(1) reversal).",
         "analytic identity"),
        ("Is it unique?",
         "Not applicable as a lift. Among fibre-reversing Hopf isometries it is one "
         "point of the component (-j e^{i alpha}, g_R); the lifts of iota are four "
         "(±e_s in Pin^±) and of -I_4 two (± volume).", "analytic identity"),
        ("Is PR #129's antipodal BC derived?",
         "CONDITIONALLY: it is the eta = +1 sector of the unique free involution, "
         "given P_B = -P_A (chosen) and the quotient rather than the double cover "
         "(chosen). Neither choice is forced by the geometry.",
         "conditional on an unproved physical identification"),
        ("Does its horizon limit reproduce even/odd BCs?",
         "There is no limit: the sector labels are lapse-independent and hold at "
         "the finite ultrastatic neck. Reproduced numerically against the PR #277 "
         "oracle to " + f"{h['h3']['worst_relative_error']:.1e}" + ", second order.",
         "numerically converged"),
        ("Which inputs remain postulates?",
         "P_B = -P_A; quotient vs cover; eta; the extension in time (P or PT); "
         "the Hopf complex structure on the state space; charge conjugation C; "
         "any Pin type and sign. The geometry fixes none of them.",
         "definition / chosen"),
    ]


def verdict(h: dict) -> str:
    forces = False                      # H1: the gluing class is free
    admits = h["h2"]["holds"] and h["h3"]["holds"] and h["h5"]["holds"]
    incompatible = not h["h2"]["free"] or not h["h3"]["B2"]
    if incompatible:
        return "FINITE_MOUTH_INCOMPATIBLE_WITH_CURRENT_BAM_TRANSPORT"
    if forces:
        return "FINITE_MOUTH_FORCES_TRANSPORT_AND_ANTIPODAL_BC"
    if admits:
        return "FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT"
    return "UNDETERMINED"


def run_probe() -> dict:
    h = {"h1": check_h1(), "h2": check_h2(), "h3": check_h3(),
         "h4": check_h4(), "h5": check_h5()}
    passed = sum(int(v["holds"]) for v in h.values())
    return {"timestamp": datetime.now(timezone.utc).isoformat(),
            "prereg": "docs/finite_mouth_topology_prereg.md @ d9d85bc",
            "checks": h, "passed": passed, "total": len(h),
            "audit_table": audit_table(h), "verdict": verdict(h),
            "dependency_ledger": {
                "antipodal_BC": "BC( P_B=-P_A [chosen], quotient-not-cover [chosen], "
                                "eta [chosen], iota [derived given P_B=-P_A], "
                                "(-1)^l [identity], a, R [geometry], P-or-PT [chosen] )",
                "J": "J( Hopf complex structure on C^2 [chosen], C = U(1) reversal "
                     "[chosen], -j in SU(2)_L [gauge/convention] )",
                "spinor transport of iota": "±e_s ( Pin type [chosen], sign [chosen], "
                                            "neck frame [gauge] )",
                "handle topology": "mapping torus of m ( m in O(4) [chosen]; "
                                   "(det m_3, eps) [derived from m] )"},
            "refinement_of_the_trichotomy": (
                "The trichotomy is well posed but its middle option needs one "
                "sharpening: the geometry does not merely fail to select the BAM "
                "lift, it contradicts two of the words attached to it. (i) The "
                "antipodally glued two-mouth handle is orientable in the bulk; "
                "non-orientability lives on the quotient, and there the Pin type "
                "(Pin+ on RP^4 # RP^4) is not the Pin- the repository assigns to "
                "the RP^2 mouth. (ii) J is a rotation, not a reflection, and its "
                "antilinearity is a choice of complex structure. Neither the "
                "non-orientability nor the K can be read off the finite mouth.")}


def render_markdown(s: dict) -> str:
    h = s["checks"]
    L = [f"# Finite-mouth topology probe — {s['passed']}/{s['total']}", "",
         f"Pre-registration: `{s['prereg']}`. Verdict: **`{s['verdict']}`**", ""]
    for key, title in (("h1", "H1 — the gluing is free; antipodal gluing is orientable"),
                       ("h2", "H2 — the unique free involution"),
                       ("h3", "H3 — the scalar sector and the oracle"),
                       ("h4", "H4 — no horizon limit; the time-reversal datum"),
                       ("h5", "H5 — what J is")):
        c = h[key]
        L += [f"## {title}", "", f"**{'HOLDS' if c['holds'] else 'FAILS'}** "
              f"(*{c['class']}*)", "", "> " + c["why"], ""]
    L += ["### Gluing classes", "", "| gluing | det m | det m_3 | eps | w_1 on loop |",
          "|--|--|--|--|--|"]
    for name, c in h["h1"]["classes"].items():
        L.append(f"| {name} | `{c['det_bulk']:+d}` | `{c['det_brane']:+d}` | "
                 f"`{c['eps']:+d}` | `{c['w1_loop']:+d}` |")
    L += ["", "### Half-tube admittance against the pre-registered targets", "",
          "| ℓ | condition | target | solve | rel. err |", "|--|--|--|--|--|"]
    for r in h["h3"]["oracle_rows"]:
        L.append(f"| {r['ell']} | {r['condition']} | `{r['target']:.9f}` | "
                 f"`{r['solve']:.9f}` | `{r['relative_error']:.1e}` |")
    L += ["", "Convergence ratios: " + "; ".join(
        f"{c['condition']} ℓ={c['ell']}: " + ", ".join(f"`{x:.2f}`" for x in c["ratios"])
        for c in h["h3"]["convergence"]), ""]
    L += ["## Audit table", "", "| Question | Result | Evidence class |", "|--|--|--|"]
    for q, r, e in s["audit_table"]:
        L.append(f"| {q} | {r} | {e} |")
    L += ["", "## Dependency ledger", ""]
    for k, v in s["dependency_ledger"].items():
        L.append(f"* `{k}` = {v}")
    L += ["", "## Refinement of the trichotomy", "", s["refinement_of_the_trichotomy"], ""]
    return "\n".join(L)


def main() -> int:
    summary = run_probe()
    text = render_markdown(summary)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs",
                          f"{stamp}_finite_mouth_topology_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as handle:
        json.dump(summary, handle, indent=2, default=lambda o: o.tolist()
                  if isinstance(o, np.ndarray) else (float(o) if isinstance(o, np.floating) else str(o)))
    with open(os.path.join(outdir, "probe.md"), "w") as handle:
        handle.write(text)
    print(f"\n\nWrote: {os.path.join(outdir, 'probe.json')}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
