"""Can any classical field already in BAM keep the finite mouth open?

THE QUESTION, AND WHY IT IS BINARY
──────────────────────────────────
``finite_mouth`` showed that a smooth traversable neck costs

    8 pi G_5 T_ab k^a k^b |_0  =  -3/b^2  =  -3/(R^2 sin^4 a)

for **every** lapse with ``N(0) > 0``. This module asks whether any classical
field the repository already contains can pay it, by building each candidate's
null-contracted stress tensor from its actual action rather than from prose.

THE REQUIREMENT IS A FLARE-OUT THEOREM, NOT A FEATURE OF N = 1
──────────────────────────────────────────────────────────────
The radial null congruence has a three-dimensional screen in ``D = 5``, so
``theta = 3 f'/f`` and Raychaudhuri reads

    dtheta/dlam = -theta^2/3 - sigma^2 - R_ab k^a k^b

At the neck ``theta = 0`` and ``sigma = 0`` by spherical symmetry, while
``dtheta/dlam = 3 f''/f = 3/b^2 > 0`` — the flare-out. Hence
``R_kk = -3/b^2``, and the ``Lambda`` term drops under null contraction, so
``8 pi G_5 T_kk = -3/b^2``. That reproduces the component result with no
reference to ``p_s``, ``N = 1``, or reflection symmetry.

**This is not a new theorem.** It is the 5D specialisation of the standard
traversable-wormhole flare-out requirement (Morris-Thorne 1988) in the local
form of Hochberg-Visser, who proved NEC violation at a throat defined as a
marginally anti-trapped surface without symmetry or asymptotic-flatness
assumptions. Recovering it here is a **validation of the finite-mouth
construction**, and is reported that way rather than as a discovery.

WHAT IS AUDITED
───────────────
Ten candidates, each reduced to its null contraction. Nine are manifestly
non-negative for structural reasons that survive any field configuration; the
tenth, a nonminimally coupled scalar, is the only one whose sign is not fixed
a priori — and it closes at a node for a reason that holds in every dimension.

WHAT IS NOT CLAIMED
───────────────────
The theorem's premises are classical Einstein gravity in 5D with the listed
matter. Higher-curvature gravity (Gauss-Bonnet is the natural ``D = 5`` term),
ghost fields, quantum stress, and a classical Dirac field are all genuine
escapes, and none is refuted here. See ``docs/source_audit_prereg.md``, frozen
before this module existed.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable, Dict, List, Optional

import numpy as np

from geometrodynamics.bulk.finite_mouth import (
    BULK_RADIUS,
    MOUTH_ANGLE,
    geometry,
    throat_radius,
)

__all__ = [
    "CONFORMAL_COUPLING_5D",
    "required_null_stress",
    "conformal_coupling",
    "one_minus_two_xi_conformal",
    "null_stress_minimal_scalar",
    "null_stress_complex_order_field",
    "null_stress_maxwell",
    "null_stress_perfect_fluid",
    "null_stress_nonminimal_at_node",
    "sign_flip_threshold",
    "raychaudhuri_expansion",
    "measure_the_flare_out_requirement",
    "measure_every_candidate_stress",
    "measure_the_nonminimal_loophole",
    "measure_the_dynamic_escape_fails",
    "measure_the_source_audit_ledger",
]

#: Conformal coupling in five spacetime dimensions, ``(D-2)/(4(D-1))``.
CONFORMAL_COUPLING_5D = 3.0 / 16.0


def required_null_stress(bulk_radius: float = BULK_RADIUS,
                         mouth_angle: float = MOUTH_ANGLE) -> float:
    """``8 pi G_5 T_kk`` the neck demands: ``-3/b^2 = -3/(R^2 sin^4 a)``."""
    return -3.0 / (bulk_radius * math.sin(mouth_angle) ** 2) ** 2


def conformal_coupling(dimension: int) -> float:
    """``xi_c = (D-2)/(4(D-1))``."""
    return (dimension - 2.0) / (4.0 * (dimension - 1.0))


def one_minus_two_xi_conformal(dimension: int) -> float:
    """``1 - 2 xi_c = D/(2(D-1))``, positive in every dimension ``D >= 2``.

    This closed form is why conformal coupling cannot rescue a node: the sign
    of the null stress there is ``sign(1 - 2 xi)``, and the conformal value
    never reaches the ``xi > 1/2`` a flip requires.
    """
    return dimension / (2.0 * (dimension - 1.0))


def sign_flip_threshold() -> float:
    """The ``xi`` at which a nonminimal scalar's node contribution changes sign."""
    return 0.5


# ── the candidate null stresses, from their actual stress tensors ───────────

def null_stress_minimal_scalar(gradient_along_k, potential=0.0):
    """``T_kk = (k^a d_a psi)^2``.

    Every ``g_ab`` term — the kinetic trace and the potential alike — drops
    under a null contraction, so ``V`` is irrelevant to the NEC *even where it
    is negative*. That is why the GL quartic's symmetry-breaking region cannot
    help.
    """
    return np.asarray(gradient_along_k, dtype=float) ** 2


def null_stress_complex_order_field(gradient_along_k, stiffness: float = 1.0):
    """``T_kk = kappa |k^a D_a q|^2`` for the repository's GL order field.

    The gradient term in the repository's throat-order energy is the ordinary
    positive one, ``(kappa/2)|grad q|^2``. Promoted covariantly with that sign,
    the null contraction is a modulus squared.
    """
    magnitude = np.abs(np.asarray(gradient_along_k, dtype=complex)) ** 2
    return stiffness * magnitude


def null_stress_maxwell(field_strength: np.ndarray, metric: np.ndarray,
                        null_vector: np.ndarray) -> float:
    """``T_kk = |F_{ab}k^b|^2``, computed from ``F`` rather than asserted.

    ``V_a = F_{ab}k^b`` satisfies ``V_a k^a = 0`` identically by antisymmetry,
    so ``V`` is orthogonal to a null vector and therefore spacelike or null:
    its square is non-negative. The ``-(1/4) g_ab F^2`` piece drops.
    """
    inverse = np.linalg.inv(metric)
    v_lower = field_strength @ null_vector
    return float(v_lower @ inverse @ v_lower)


def null_stress_perfect_fluid(density, pressure, flow_dot_k) -> np.ndarray:
    """``T_kk = (rho + p)(u^a k_a)^2``: non-negative exactly when the NEC holds."""
    return (np.asarray(density, dtype=float)
            + np.asarray(pressure, dtype=float)) * np.asarray(
                flow_dot_k, dtype=float) ** 2


def null_stress_nonminimal_at_node(gradient_along_k, coupling: float):
    """``R_kk ∝ (1 - 2 xi)(dq/dlam)^2`` at a zero of the field.

    The full null-contracted equation for
    ``S = -(1/2) int [ (grad q)^2 + xi R q^2 + 2V ]`` is

        (1 - 8 pi G xi q^2) R_kk = 8 pi G [ (k.grad q)^2 - xi d^2(q^2)/dlam^2 ]

    At a node the prefactor is exactly ``1``, and ``d^2(q^2)/dlam^2`` reduces to
    ``2 (dq/dlam)^2`` because the ``q q''`` term carries a factor ``q = 0``.
    BAM places the defect core at ``q = 0``, so this is the relevant point.
    """
    return (1.0 - 2.0 * coupling) * np.asarray(gradient_along_k, dtype=float) ** 2


# ── the geometry side ───────────────────────────────────────────────────────

def raychaudhuri_expansion(s, neck_radius: float):
    """``theta = 3 f'/f`` for the radial null congruence (screen dimension 3)."""
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    return 3.0 * (s / f) / f


@lru_cache(maxsize=8)
def measure_the_flare_out_requirement(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """S0 — the requirement, re-derived from Raychaudhuri rather than ``p_s``."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    s = np.linspace(-g.half_length, g.half_length, 601)
    theta = raychaudhuri_expansion(s, b)
    # dtheta/dlam from the closed form, and R_kk = -3 b^2/f^4
    f = throat_radius(s, b)
    dtheta = 3.0 * (b * b / f ** 3) / f - 3.0 * (s / f) ** 2 / (f * f)
    ricci_kk = -3.0 * b * b / f ** 4
    residual = float(np.max(np.abs(dtheta - (-theta ** 2 / 3.0 - ricci_kk))))
    mid = len(s) // 2
    return {
        "neck_radius": b,
        "theta_at_neck": float(theta[mid]),
        "dtheta_at_neck": float(dtheta[mid]),
        "expected_dtheta": 3.0 / (b * b),
        "ricci_kk_at_neck": float(ricci_kk[mid]),
        "required_null_stress": required_null_stress(bulk_radius, mouth_angle),
        "raychaudhuri_residual": residual,
        "screen_dimension": 3,
        "identity_holds": bool(residual < 1e-9),
        "flare_out_is_positive": bool(dtheta[mid] > 0.0),
        "matches_the_component_result": bool(
            abs(float(ricci_kk[mid]) - required_null_stress(bulk_radius, mouth_angle))
            < 1e-9),
        "why_this_is_stronger": (
            "Derived from theta and Raychaudhuri alone, with no reference to "
            "p_s, to N = 1, or to reflection symmetry. Smooth radial flare-out "
            "plus Einstein gravity forces T_kk < 0, whatever the lapse and "
            "whatever the matter."),
        "normalisation": (
            "THE SIGN IS INVARIANT; THE MAGNITUDE IS NOT. R_ab k^a k^b scales "
            "as lambda^2 under k -> lambda k, so -3/b^2 is quoted in the local "
            "normalisation k^{hat t} = 1 at the neck. theta = 3 f'/f is "
            "likewise the ultrastatic/locally normalised form; for an affinely "
            "parametrised radial null ray in a general lapse it is "
            "theta = 3(E/N) f'/f with E the conserved energy, and the neck "
            "value -3/b^2 is recovered by setting k^{hat t} = 1 there. The "
            "flare-out theorem depends only on the sign, which no rescaling "
            "can move."),
        "attribution": (
            "NOT a new theorem: this is the 5D specialisation of the standard "
            "traversable-wormhole flare-out requirement (Morris-Thorne 1988) in "
            "the local form of Hochberg-Visser, who proved NEC violation at a "
            "throat defined as a marginally anti-trapped surface. Recovering it "
            "here VALIDATES the finite-mouth construction; it is not a "
            "discovery, and the repository has too few external anchors to "
            "waste one by miscrediting it."),
    }


@lru_cache(maxsize=4)
def measure_every_candidate_stress() -> Dict[str, object]:
    """S1 — each candidate's null stress, computed rather than asserted."""
    rng = np.random.default_rng(20260831)
    gradients = rng.normal(size=64) * 3.0
    rows: List[dict] = []

    rows.append({
        "id": "C1", "candidate": "minimally coupled scalar psi",
        "null_stress_min": float(np.min(null_stress_minimal_scalar(gradients))),
        "sign": "NON-NEGATIVE",
        "reason": "T_kk = (k.grad psi)^2; all g_ab terms drop under the null "
                  "contraction"})

    complex_grad = rng.normal(size=64) + 1j * rng.normal(size=64)
    rows.append({
        "id": "C2", "candidate": "complex GL throat-order field q",
        "null_stress_min": float(np.min(
            null_stress_complex_order_field(complex_grad, stiffness=1.7))),
        "sign": "NON-NEGATIVE",
        "reason": "T_kk = kappa |k.Dq|^2 with the repository's positive "
                  "gradient term; this is the field introduced to REPRESENT "
                  "the throat, and it cannot pay for it"})

    # C3: contract the FULL scalar stress tensor, potential included, with a
    # genuine null vector -- rather than calling a routine that ignores V.
    g = geometry()
    metric = np.diag([-1.0, 1.0, g.neck_radius ** 2,
                      g.neck_radius ** 2, g.neck_radius ** 2])
    inverse = np.linalg.inv(metric)
    null_vector = np.array([1.0, 1.0, 0.0, 0.0, 0.0])
    assert abs(null_vector @ metric @ null_vector) < 1e-15, "k must be null"
    potential_spread = 0.0
    for potential in (-9.0, -1.0, 0.0, 4.0):
        covector = rng.normal(size=5)
        kinetic = float(covector @ inverse @ covector)
        stress = (np.outer(covector, covector)
                  - 0.5 * metric * (kinetic + 2.0 * potential))
        contracted = float(null_vector @ stress @ null_vector)
        potential_spread = max(potential_spread,
                               abs(contracted - (covector @ null_vector) ** 2))
    rows.append({
        "id": "C3", "candidate": "GL potential V(q), incl. symmetry-breaking",
        "null_stress_min": potential_spread,
        "sign": "EXACTLY ZERO",
        "reason": f"contracting the full T_ab for V = -9, -1, 0, +4 leaves "
                  f"(k.grad q)^2 to {potential_spread:.1e}: the potential drops "
                  "identically, so V < 0 regions are irrelevant to the NEC"})

    # Maxwell: build a random antisymmetric F in the actual metric
    worst_maxwell, worst_ortho = math.inf, 0.0
    for _ in range(200):
        raw = rng.normal(size=(5, 5))
        strength = raw - raw.T
        worst_maxwell = min(worst_maxwell,
                            null_stress_maxwell(strength, metric, null_vector))
        worst_ortho = max(worst_ortho,
                          abs(float((strength @ null_vector) @ null_vector)))
    rows.append({
        "id": "C4", "candidate": "Maxwell / KK gauge field",
        "null_stress_min": worst_maxwell, "sign": "NON-NEGATIVE",
        "reason": f"T_kk = |F.k|^2 with V.k = 0 to {worst_ortho:.1e} over 200 "
                  "random field strengths, so V is spacelike or null"})

    rows.append({
        "id": "C5", "candidate": "cosmological constant",
        "null_stress_min": 0.0, "sign": "EXACTLY ZERO",
        "reason": "T_ab ∝ g_ab, and g_ab k^a k^b = 0"})

    fluid = null_stress_perfect_fluid(rng.uniform(0.1, 5.0, 64),
                                      rng.uniform(-0.03, 2.0, 64),
                                      rng.normal(size=64))
    rows.append({
        "id": "C6", "candidate": "perfect fluid obeying the NEC",
        "null_stress_min": float(np.min(fluid)), "sign": "NON-NEGATIVE",
        "reason": "T_kk = (rho+p)(u.k)^2, zero only for rho + p = 0"})

    rows.append({
        "id": "C7", "candidate": "classical 5D gravitational waves",
        "null_stress_min": 0.0, "sign": "EXACTLY ZERO",
        "reason": "R_ab = 0 in vacuum, so R_kk = 0 identically; the Isaacson "
                  "effective energy is positive"})

    rows.append({
        "id": "C8", "candidate": "non-orientable identification / wrap sign",
        "null_stress_min": 0.0, "sign": "NO CONTRIBUTION",
        "reason": "changes global boundary data, not the local R_kk. "
                  "hopf/spinor.py is SU(2) holonomy transport with no stress "
                  "tensor, so U_spin is a transport object, not a source"})

    rows.append({
        "id": "C9", "candidate": "projected bulk Weyl stress (#167/#168)",
        "null_stress_min": 0.0, "sign": "NO CONTRIBUTION TO THE 5D EQUATION",
        "reason": "there T^(5)_ab = 0 exactly; the projection is an effective "
                  "4D BRANE source and cannot appear on the right-hand side of "
                  "the 5D Einstein equation, which is what the neck needs"})

    negative = [r for r in rows if r["null_stress_min"] < -1e-12]
    return {
        "rows": rows,
        "candidates_with_negative_null_stress": [r["id"] for r in negative],
        "all_non_negative": bool(not negative),
        "the_order_field_cannot_pay_its_own_bill": (
            "C2 is the field the repository introduced AS the throat-order "
            "degree of freedom. With the ordinary positive gradient term its "
            "null stress is a modulus squared, so it cannot support the object "
            "it was introduced to represent. That is stronger than saying its "
            "constants were never derived."),
    }


@lru_cache(maxsize=4)
def measure_the_nonminimal_loophole() -> Dict[str, object]:
    """S2 — the one candidate whose sign is not fixed a priori, and it closes."""
    gradients = np.array([0.3, 1.0, 2.5])
    couplings = {
        "minimal xi = 0": 0.0,
        "conformal 4D xi = 1/6": conformal_coupling(4),
        "conformal 5D xi = 3/16": conformal_coupling(5),
        "conformal 6D xi = 1/5": conformal_coupling(6),
        "threshold xi = 1/2": 0.5,
        "beyond threshold xi = 0.8": 0.8,
    }
    rows = []
    for name, xi in couplings.items():
        value = null_stress_nonminimal_at_node(gradients, xi)
        rows.append({"coupling": name, "xi": xi,
                     "one_minus_two_xi": 1.0 - 2.0 * xi,
                     "null_stress_sign": ("negative" if value[0] < -1e-15
                                          else "zero" if abs(value[0]) < 1e-15
                                          else "positive")})
    dimensions = [{"D": d, "xi_conformal": conformal_coupling(d),
                   "one_minus_two_xi": one_minus_two_xi_conformal(d),
                   "closed_form_matches": bool(abs(
                       (1.0 - 2.0 * conformal_coupling(d))
                       - one_minus_two_xi_conformal(d)) < 1e-15)}
                  for d in (3, 4, 5, 6, 10, 26)]
    return {
        "rows": rows, "dimensions": dimensions,
        "threshold": sign_flip_threshold(),
        "conformal_5d": CONFORMAL_COUPLING_5D,
        "one_minus_two_xi_5d": 1.0 - 2.0 * CONFORMAL_COUPLING_5D,
        "conformal_never_flips": bool(
            all(d["one_minus_two_xi"] > 0.0 for d in dimensions)),
        "closed_form": "1 - 2 xi_c = D/(2(D-1)) > 0 for every D >= 2",
        "why_the_node_matters": (
            "At q = 0 the prefactor 1 - 8 pi G xi q^2 is exactly 1 and "
            "d^2(q^2)/dlam^2 = 2 (dq/dlam)^2, because the q q'' term carries a "
            "factor q. BAM places the defect core AT the node, so the sign "
            "there is sign(1 - 2 xi)."),
        "the_loophole_closes": (
            "A flip needs xi > 1/2. Conformal coupling gives D/(2(D-1)), which "
            "is 5/8 in five dimensions and never reaches zero in ANY dimension. "
            "So even a conformally improved throat-order field remains "
            "NEC-positive at a smooth q = 0 core."),
    }


@lru_cache(maxsize=4)
def measure_the_dynamic_escape_fails(samples: int = 4000) -> Dict[str, object]:
    """S3 — nonzero ``K_ij`` cannot rescue it either.

    Integrates ``dtheta/dlam = -theta^2/3 - sigma^2 - R_kk`` forward for a
    battery of **non-negative** ``R_kk`` profiles and shear, starting from a
    converging ray, and looks for any turning point. This is a falsification
    attempt: a single ``theta: - -> +`` would mean a sign error upstream.
    """
    profiles = {
        "vacuum R_kk = 0": lambda x: np.zeros_like(x),
        "constant R_kk = 0.5": lambda x: 0.5 * np.ones_like(x),
        "gaussian bump": lambda x: 2.0 * np.exp(-((x - 1.0) ** 2) / 0.2),
        "oscillating, non-negative": lambda x: 1.0 + np.cos(6.0 * x) ** 2,
        "decaying tail": lambda x: 3.0 / (1.0 + x ** 2),
    }
    rows, turning, late = [], [], []
    for name, profile in profiles.items():
        for theta0 in (-0.05, -0.5, -2.0):
            # The focusing theorem bounds the caustic by (D-2)/|theta0| = 3/|theta0|
            # in vacuum, and non-negative R_kk or shear can only shorten it. An
            # arbitrary fixed window would report a false negative here: an
            # earlier draft used lam <= 6 and missed the theta0 = -0.05 case,
            # whose analytic caustic sits at lam = 60.
            bound = 3.0 / abs(theta0)
            lam = np.linspace(0.0, 1.2 * bound, samples)
            h = lam[1] - lam[0]
            theta, previous = theta0, theta0
            turned, caustic = False, None
            for index, value in enumerate(profile(lam)):
                shear = 0.3 * abs(math.sin(value))          # sigma^2 >= 0
                theta = theta + h * (-theta ** 2 / 3.0 - shear - value)
                if previous < 0.0 and theta > 0.0:
                    turned = True
                if not np.isfinite(theta) or theta < -1e6:
                    caustic = float(lam[index])
                    break
                previous = theta
            rows.append({"profile": name, "theta_initial": theta0,
                         "focusing_bound": bound,
                         "measured_caustic": caustic,
                         "turned_positive": turned,
                         "focused_within_bound": bool(
                             caustic is not None and caustic <= bound * 1.05)})
            if turned:
                turning.append(f"{name} @ theta0={theta0}")
            if not rows[-1]["focused_within_bound"]:
                late.append(f"{name} @ theta0={theta0}")
    return {
        "rows": rows,
        "any_turning_point": bool(turning),
        "turning_cases": turning,
        "cases_not_focused_within_the_bound": late,
        "all_focused_within_the_analytic_bound": bool(not late),
        "focusing_bound_formula": "lam_caustic <= (D-2)/|theta_0| = 3/|theta_0|",
        "why_the_window_matters": (
            "An earlier draft integrated to a fixed lam <= 6 and reported that "
            "not every ray focused -- which was an artefact of the window, not "
            "a failure of the theorem: at theta0 = -0.05 the vacuum caustic "
            "sits at lam = 3/|theta0| = 60. Each ray is now integrated past its "
            "own analytic bound. The no-turning-point result never depended on "
            "the window; the focusing sub-claim did."),
        "why": (
            "With T_kk >= 0 every term on the right of Raychaudhuri is "
            "non-positive, so theta is non-increasing: a ray that enters "
            "converging cannot flare back out. Changing K_ij moves where the "
            "throat sits in a slice and whether the neck is trapped, but it "
            "cannot make a smooth radial null bundle turn around in ordinary "
            "Einstein vacuum."),
        "the_stronger_statement": (
            "Not just static BAM: ANY smooth two-way traversable BAM throat in "
            "classical Einstein gravity requires null-convergence violation "
            "somewhere. A numerical-relativity round looking for a vacuum "
            "dynamic rescue is looking for something the null equations "
            "forbid."),
    }


@lru_cache(maxsize=4)
def measure_the_source_audit_ledger() -> Dict[str, object]:
    """S4 — the verdict, and the branches it leaves open."""
    requirement = measure_the_flare_out_requirement()
    candidates = measure_every_candidate_stress()
    loophole = measure_the_nonminimal_loophole()
    dynamic = measure_the_dynamic_escape_fails()
    verdict_is_no = bool(candidates["all_non_negative"]
                         and loophole["conformal_never_flips"]
                         and not dynamic["any_turning_point"])
    return {
        "required_null_stress": requirement["required_null_stress"],
        "verdict": ("NO -- the current classical BAM field content cannot "
                    "support a traversable particle throat" if verdict_is_no
                    else "SOMETHING SUPPLIES IT -- see the rows"),
        "verdict_is_no": verdict_is_no,
        "entries": [
            {"claim": "the neck NEC price is an artefact of N = 1 or of symmetry",
             "verdict": "NO -- IT IS A FLARE-OUT THEOREM",
             "evidence": "derived from theta and Raychaudhuri alone, residual "
                         f"{requirement['raychaudhuri_residual']:.1e}, with no "
                         "reference to p_s"},
            {"claim": "that theorem is new",
             "verdict": "NO -- IT IS MORRIS-THORNE / HOCHBERG-VISSER IN 5D",
             "evidence": "recovering a known theorem validates the finite-mouth "
                         "construction; it is a second external anchor, not a "
                         "discovery"},
            {"claim": "the GL throat-order field q can support the throat",
             "verdict": "NO",
             "evidence": "T_kk = kappa |k.Dq|^2 >= 0 with the repository's own "
                         "positive gradient term; the field introduced to "
                         "represent the throat cannot pay for it"},
            {"claim": "the GL potential's symmetry-breaking region helps",
             "verdict": "NO -- IT DROPS IDENTICALLY",
             "evidence": "every g_ab term vanishes under null contraction, so "
                         "V < 0 is irrelevant to the NEC"},
            {"claim": "a conformally improved scalar rescues it at the core",
             "verdict": "NO -- 1 - 2 xi_c = D/(2(D-1)) > 0 IN EVERY DIMENSION",
             "evidence": f"5/8 in five dimensions; a flip needs xi > 1/2, and "
                         "BAM puts the defect core exactly at the q = 0 node"},
            {"claim": "the #167/#168 Weyl mechanism can keep this throat open",
             "verdict": "NO -- WRONG EQUATION",
             "evidence": "there T^(5)_ab = 0; a projected Weyl tensor is an "
                         "effective 4D brane source and cannot appear in the 5D "
                         "Einstein equation the neck lives in"},
            {"claim": "nonzero K_ij (gravitational waves) can rescue it",
             "verdict": "NO",
             "evidence": f"{len(dynamic['rows'])} Raychaudhuri integrations "
                         "with non-negative R_kk and shear produce zero turning "
                         "points; theta is non-increasing and focuses"},
        ],
        "what_this_forces": (
            "If every existing BAM action gives T_kk >= 0, the honest headline "
            "is a negative result about the CURRENT field content, not about "
            "wormholes in general."),
        "the_remaining_branches": {
            "1 accept the horizon": "take the Tangherlini branch N(0) = 0 as the "
                                    "particle, and abandon MTY traversability",
            "2 add a ghost": "a wrong-sign field buys T_kk < 0 and brings an "
                             "obvious stability problem",
            "3 leave Einstein gravity": "higher curvature -- in D = 5 the "
                                        "natural term is Gauss-Bonnet, which is "
                                        "exactly where the theorem's premise "
                                        "fails; Lovelock wormholes can satisfy "
                                        "the matter NEC",
            "4 quantum stress": "Casimir-type support, meaning the particle "
                                "geometry is no longer derived from classical GR",
            "5 reinterpret": "particle exchange does not require a traversable "
                             "throat at all",
        },
        "what_the_audit_does_not_do": (
            "It does not choose among those five, and it refutes none of them. "
            "Each is a premise of the theorem rather than a consequence."),
    }
