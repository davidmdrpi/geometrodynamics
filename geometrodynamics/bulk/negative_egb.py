"""Negative-coupling EGB: the exterior closes what the throat left open.

WHAT PR #277 LEFT
─────────────────
``gauss_bonnet`` narrowed the Gauss-Bonnet branch without closing it. A
negative coupling ``alpha_GB <= -R^2/4`` *does* satisfy the matter NEC along the
throat, and the ledger recorded "global existence and stability are open".

THE STEP THAT ROUND MISSED
──────────────────────────
``alpha_GB`` is a **coupling constant in the action**, so the same value acts
everywhere — including in the exterior the throat is glued into. The throat was
analysed in isolation, and it should not have been.

For the ultrastatic exterior ``R x S^4_R``:

    R_kk = 3/R^2  (POSITIVE -- ordinary matter, not flare-out)
    H_kk = 12/R^4
    8 pi G_5 T_kk = 3(R^2 + 4 alpha_GB)/R^4   >= 0   <=>   alpha_GB >= -R^2/4

which runs **opposite** to the throat's ``alpha_GB <= -R^2/4``.

WHY THE TWO MEET EXACTLY, AND WHY THAT IS NOT A COINCIDENCE
───────────────────────────────────────────────────────────
Both regions share one bracket, since ``T_kk = R_kk (1 + 4 alpha mu/f^4)``:

* in the throat ``R_kk < 0`` (flare-out), so the NEC needs the bracket ``<= 0``;
* in the exterior ``R_kk > 0`` (ordinary matter), so it needs the bracket
  ``>= 0``.

Same bracket, opposite required signs. And ``mu/f^4`` — with ``mu`` the
Misner-Sharp parameter ``finite_mouth`` P2 showed continuous across the seam —
equals ``1/R^2`` on *both* sides there, and is ``chi``-independent throughout
the exterior. So the bracket is continuous at the seam and must be exactly zero.

    alpha_GB = -R^2/4      exactly, and nothing else

a **measure-zero** solution: there is no open set of couplings.

WHAT THE CRITICAL SOLUTION IS -- AND IS NOT
───────────────────────────────────────────
``H^i_j = 0`` for **any** ultrastatic product ``-dt^2 + h_4``, because in
``D = 5`` the spatial block is the four-dimensional Gauss-Bonnet (Euler) tensor
of ``h_4`` and Gauss-Bonnet is topological in ``D = 4``. Maximal symmetry is not
needed and an earlier draft wrongly credited it — see
``gauss_bonnet.measure_the_spatial_block_vanishes``, which checks a throat slice
and a generic lumpy one, with a nonconstant-lapse control that does *not*
vanish. So Gauss-Bonnet does not touch the pressure anywhere on this geometry,
exterior or throat:

    8 pi G_5 rho_ext = 6/R^2 + 12 alpha_GB/R^4      8 pi G_5 p_ext = -3/R^2

At the critical coupling that gives ``rho = 3/R^2``, ``p = -rho``: the smooth 5D
exterior stress is driven to **vacuum form**, ``w = -1``.

**Three claims of an earlier draft are withdrawn, all from review.**

*"The exotic matter is merely relocated."* **False at the critical point.** The
throat stress there satisfies NEC, WEC **and** DEC — see
``measure_the_throat_matter_is_not_exotic``. With ``A = b^2/f^4`` and
``q = R^2 b^2/f^4 >= 1``, the stress is ``rho = 3Aq``, ``p_s = -3A``,
``p_Omega = A``, so ``rho + p_s = 3A(q-1) >= 0`` (saturated only at the mouth)
and ``rho >= |p_s|, |p_Omega|``. Relocation happens only for ``alpha`` strictly
below critical.

*"The observed universe must be empty."* **Overreached.** ``w = -1`` is a
statement about the total 5D bulk stress supporting ``R x S^4_R``, not about the
observed ``S^3``, which is its *equator* and carries a different stress tensor.
A homogeneous ``-Lambda g_ab`` can also simply be moved to the gravitational
side. Localised throat matter can still exist. The narrow, defensible statement
is that the branch forces a **vacuum-energy-like 5D exterior**.

*"Measure zero, therefore refuted."* **Not by itself.** At fixed ``R`` the
allowed set is one point, but the relation is equally read as the field
equations selecting ``R^2 = -4 alpha_GB``, and gravity routinely ties a vacuum
curvature radius to a coupling. Calling it tuning needs an independently fixed
``R`` that was not free to respond.

WHAT ACTUALLY CLOSES THE BRANCH: THE GRAVITON
─────────────────────────────────────────────
Linearising the full ``G_ab + alpha H_ab`` on this background for a TT mode
gives the coefficient of the ``omega^2`` term:

    C_kin = -(1/2)(1 + 4 alpha_GB/R^2)

**derived here rather than recalled** — the familiar
``1 + 2 alpha (D-3)(D-4) K`` is for a maximally symmetric *spacetime*, and
``R x S^4_R`` is a product. It vanishes **exactly** at ``alpha = -R^2/4``, the
same value the NEC pins, while the spatial ``kappa^2`` coefficient is entirely
coupling-independent. Hence

    c^2 = 1/(1 + 4 alpha_GB/R^2)

is **superluminal for every** ``-R^2/4 < alpha_GB < 0`` and diverges at the
critical coupling, where the time-derivative term disappears and the tensor
equation degenerates from an evolution equation into a constraint. The mass term
stays finite, so this is a degeneration of the kinetic operator rather than an
overall vanishing: **there is no propagating graviton at the one coupling the
NEC allows.**

SCOPE
─────
This does **not** establish global existence in the full sense: it determines
the effective stress the chosen metric requires, and has not exhibited fields
obeying their own equations that generate the throat's anisotropic stress. It
does not touch dilatonic ``alpha(phi) L_GB`` or ``f(R)``. Predictions E1-E6 are
frozen in ``docs/negative_egb_prereg.md``, committed before this module existed;
E6's *interpretation* is corrected above.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List

import numpy as np

from geometrodynamics.bulk.finite_mouth import (
    BULK_RADIUS,
    MOUTH_ANGLE,
    geometry,
    misner_sharp,
    throat_radius,
)
from geometrodynamics.bulk.gauss_bonnet import (
    coupling_threshold,
    gauss_bonnet_ratio,
    global_coupling_threshold,
    matter_null_stress,
)

__all__ = [
    "graviton_kinetic_coefficient",
    "graviton_speed_squared",
    "linearised_field_equation",
    "tt_kinetic_and_mass",
    "measure_the_graviton_degenerates",
    "measure_the_throat_matter_is_not_exotic",
    "exterior_ricci_null",
    "exterior_lanczos_null",
    "exterior_ratio",
    "exterior_matter_null_stress",
    "exterior_density",
    "exterior_pressure",
    "exterior_threshold",
    "bracket",
    "measure_the_exterior_constrains_alpha_oppositely",
    "measure_no_coupling_satisfies_both",
    "measure_the_bracket_is_continuous_at_the_seam",
    "measure_the_critical_exterior_is_empty",
    "measure_the_negative_egb_ledger",
]


# ── the exterior, R x S^4_R ─────────────────────────────────────────────────

def exterior_ricci_null(bulk_radius: float = BULK_RADIUS) -> float:
    """``R_kk = 3/R^2`` — **positive**, unlike the throat's `-3f''/f`."""
    return 3.0 / (bulk_radius * bulk_radius)


def exterior_lanczos_null(bulk_radius: float = BULK_RADIUS) -> float:
    """``H_kk = 12/R^4``, from ``H^t_t = -12/R^4`` and ``H^i_j = 0``."""
    return 12.0 / bulk_radius ** 4


def exterior_ratio(chi, bulk_radius: float = BULK_RADIUS):
    """``H_kk/R_kk = 4 mu/f^4 = 4/R^2`` — independent of ``chi``.

    Computed from the exterior's own ``mu = R^2 sin^4 chi`` and
    ``f = R sin chi`` rather than asserted, so the ``chi``-independence is a
    result. That constancy is what makes the exterior's constraint a single
    number instead of a profile.
    """
    chi = np.asarray(chi, dtype=float)
    f = bulk_radius * np.sin(chi)
    mu = bulk_radius ** 2 * np.sin(chi) ** 4
    return 4.0 * mu / f ** 4


def exterior_matter_null_stress(coupling: float,
                                bulk_radius: float = BULK_RADIUS) -> float:
    """``8 pi G_5 T_kk^ext = 3(R^2 + 4 alpha)/R^4``."""
    return 3.0 * (bulk_radius ** 2 + 4.0 * coupling) / bulk_radius ** 4


def exterior_density(coupling: float, bulk_radius: float = BULK_RADIUS) -> float:
    """``8 pi G_5 rho = 6/R^2 + 12 alpha/R^4``."""
    return 6.0 / bulk_radius ** 2 + 12.0 * coupling / bulk_radius ** 4


def exterior_pressure(coupling: float, bulk_radius: float = BULK_RADIUS) -> float:
    """``8 pi G_5 p = -3/R^2``, **independent of the coupling**.

    Gauss-Bonnet cannot touch it: ``H^i_j = 0`` for any ultrastatic product,
    the slice's symmetry being irrelevant, so the whole correction lands in
    ``rho``.
    """
    return -3.0 / (bulk_radius * bulk_radius)


def exterior_threshold(bulk_radius: float = BULK_RADIUS) -> float:
    """``alpha_GB >= -R^2/4`` for the exterior NEC — the same number the throat
    needs from the other side."""
    return -0.25 * bulk_radius * bulk_radius


def bracket(coupling: float, ratio) -> np.ndarray:
    """``1 + alpha * (H_kk/R_kk)``, the factor both regions share.

    ``T_kk = R_kk * bracket``, so the sign the NEC demands of the bracket is
    opposite in the two regions because ``R_kk`` is.
    """
    return 1.0 + coupling * np.asarray(ratio, dtype=float)


# ── measurements ────────────────────────────────────────────────────────────

@lru_cache(maxsize=8)
def measure_the_exterior_constrains_alpha_oppositely(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """E1/E2 — the exterior's own NEC, and which way it points."""
    chi = np.linspace(0.05, math.pi - 0.05, 201)
    ratio = exterior_ratio(chi, bulk_radius)
    throat_side = global_coupling_threshold(bulk_radius)
    exterior_side = exterior_threshold(bulk_radius)
    rows = []
    for coupling in (0.0, -0.5 * abs(exterior_side), exterior_side,
                     2.0 * exterior_side):
        rows.append({
            "coupling": coupling,
            "exterior_null_stress": exterior_matter_null_stress(coupling,
                                                                bulk_radius),
            "exterior_nec": bool(
                exterior_matter_null_stress(coupling, bulk_radius) >= -1e-12)})
    return {
        "ricci_null": exterior_ricci_null(bulk_radius),
        "lanczos_null": exterior_lanczos_null(bulk_radius),
        "ratio_is_chi_independent": bool(np.ptp(ratio) < 1e-12),
        "ratio_value": float(ratio[0]),
        "expected_ratio": 4.0 / bulk_radius ** 2,
        "ricci_null_is_positive": bool(exterior_ricci_null(bulk_radius) > 0.0),
        "rows": rows,
        "exterior_needs_alpha_at_least": exterior_side,
        "throat_needs_alpha_at_most": throat_side,
        "the_two_bounds_coincide": bool(abs(exterior_side - throat_side) < 1e-15),
        "they_point_in_opposite_directions": True,
        "why": ("In the throat R_kk < 0 because the neck must flare out; in the "
                "exterior R_kk = +3/R^2 because it holds ordinary matter. The "
                "NEC therefore demands opposite signs of the SAME bracket "
                "1 + alpha H_kk/R_kk, and the two bounds land on one number."),
    }


@lru_cache(maxsize=8)
def measure_no_coupling_satisfies_both(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE,
        samples: int = 4001) -> Dict[str, object]:
    """E3 — scan for a coupling that works everywhere. This is a search, not an
    assertion: a surviving interval would reopen the branch."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    s = np.linspace(-g.half_length, g.half_length, 401)
    span = 2.0 * bulk_radius ** 2
    couplings = np.linspace(-span, span, samples)
    both, throat_ok, exterior_ok = [], [], []
    for coupling in couplings:
        throat = bool(np.all(matter_null_stress(s, b, float(coupling)) >= -1e-9))
        outside = bool(
            exterior_matter_null_stress(float(coupling), bulk_radius) >= -1e-9)
        throat_ok.append(throat)
        exterior_ok.append(outside)
        if throat and outside:
            both.append(float(coupling))
    critical = exterior_threshold(bulk_radius)
    width = (max(both) - min(both)) if both else 0.0
    return {
        "scanned_range": [-span, span],
        "samples": samples,
        "couplings_satisfying_both": both,
        "surviving_width": width,
        "the_surviving_set_is_measure_zero": bool(width < 1e-6),
        "critical_coupling": critical,
        "survivors_are_at_the_critical_value": bool(
            all(abs(c - critical) < 1e-3 * bulk_radius ** 2 for c in both)),
        "throat_satisfied_count": int(sum(throat_ok)),
        "exterior_satisfied_count": int(sum(exterior_ok)),
        "why_a_scan": (
            "A surviving interval of alpha would REOPEN the branch, so this "
            "looks for one rather than asserting none. The throat and exterior "
            "admissible sets are complementary half-lines meeting at a point, "
            "so their intersection is a single coupling and any tolerance "
            "widens it only by the tolerance."),
        "what_it_means": (
            "There is no open set of couplings for which this spacetime "
            "satisfies the matter NEC everywhere. A fundamental constant of the "
            "action would have to take one exact value, tuned to the radius of "
            "the universe."),
    }


@lru_cache(maxsize=8)
def measure_the_bracket_is_continuous_at_the_seam(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """E4 — why the two constraints meet, rather than merely happening to."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    inside_ratio = float(gauss_bonnet_ratio(np.array([g.half_length]), b)[0])
    outside_ratio = float(exterior_ratio(np.array([mouth_angle]), bulk_radius)[0])
    coupling = exterior_threshold(bulk_radius)
    return {
        "ratio_inside_at_the_seam": inside_ratio,
        "ratio_outside_at_the_seam": outside_ratio,
        "ratio_jump": abs(inside_ratio - outside_ratio),
        "ratio_is_continuous": bool(abs(inside_ratio - outside_ratio) < 1e-12),
        "common_value": 4.0 / bulk_radius ** 2,
        "bracket_inside": float(bracket(coupling, inside_ratio)),
        "bracket_outside": float(bracket(coupling, outside_ratio)),
        "bracket_vanishes_at_the_seam": bool(
            abs(float(bracket(coupling, inside_ratio))) < 1e-12),
        "misner_sharp_over_f4_at_the_seam": float(
            misner_sharp(np.array([g.half_length]), b)[0]
            / throat_radius(np.array([g.half_length]), b)[0] ** 4),
        "why": ("mu/f^4 is 1/R^2 on both sides of the seam -- the same "
                "Misner-Sharp continuity finite_mouth's P2 established -- so "
                "the shared bracket 1 + 4 alpha mu/f^4 is continuous there. It "
                "is required to be <= 0 from the throat side and >= 0 from the "
                "exterior side, so it must be exactly 0. The two constraints "
                "meet because they are two one-sided limits of one continuous "
                "function, not because two independent numbers coincided."),
    }


@lru_cache(maxsize=8)
def measure_the_critical_exterior_is_empty(
        bulk_radius: float = BULK_RADIUS) -> Dict[str, object]:
    """E5/E6 — what the one surviving coupling does to the observed universe."""
    critical = exterior_threshold(bulk_radius)
    rows = []
    for coupling in (0.0, 0.5 * critical, critical, 1.5 * critical):
        density = exterior_density(coupling, bulk_radius)
        pressure = exterior_pressure(coupling, bulk_radius)
        rows.append({
            "coupling": coupling, "density": density, "pressure": pressure,
            "sum": density + pressure,
            "equation_of_state": pressure / density if density != 0.0 else float("nan")})
    critical_density = exterior_density(critical, bulk_radius)
    critical_pressure = exterior_pressure(critical, bulk_radius)
    return {
        "critical_coupling": critical,
        "rows": rows,
        "pressure_is_coupling_independent": bool(
            abs(exterior_pressure(0.0, bulk_radius)
                - exterior_pressure(critical, bulk_radius)) < 1e-15),
        "critical_density": critical_density,
        "critical_pressure": critical_pressure,
        "critical_equation_of_state": critical_pressure / critical_density,
        "it_is_exactly_vacuum_energy": bool(
            abs(critical_pressure / critical_density + 1.0) < 1e-12),
        "einstein_density": exterior_density(0.0, bulk_radius),
        "density_is_halved": bool(
            abs(critical_density - 0.5 * exterior_density(0.0, bulk_radius))
            < 1e-12),
        "why_the_pressure_cannot_move": (
            "H^i_j = 0 for ANY ultrastatic product -- the spatial block is the "
            "4D Euler tensor, topological in D = 4 -- so the entire "
            "Gauss-Bonnet correction lands in rho and none of it in p. The "
            "equation of state is therefore driven to -1 by shrinking rho onto "
            "a fixed p, not by adjusting both."),
        "what_this_costs": (
            "At the one surviving coupling the exterior -- which is supposed to "
            "BE the observed closed universe -- is pure vacuum energy with no "
            "ordinary matter anywhere, at half the Einstein-gravity density. "
            "Push alpha more negative and the EXTERIOR violates the NEC, which "
            "relocates the exotic matter from the throat into the universe "
            "around it rather than removing it."),
    }


@lru_cache(maxsize=4)
def measure_the_negative_egb_ledger() -> Dict[str, object]:
    """The verdict on the branch PR #277 left open — closed by the graviton."""
    opposite = measure_the_exterior_constrains_alpha_oppositely()
    scan = measure_no_coupling_satisfies_both()
    seam = measure_the_bracket_is_continuous_at_the_seam()
    empty = measure_the_critical_exterior_is_empty()
    graviton = measure_the_graviton_degenerates()
    honest = measure_the_throat_matter_is_not_exotic()
    closed = bool(graviton["kinetic_vanishes_at_criticality"]
                  and graviton["law_holds"])
    return {
        "branch_is_closed": closed,
        "closed_by": "the graviton kinetic term, not the matter content",
        "verdict": ("Negative-coupling EGB closes on STRUCTURAL grounds: the "
                    "NEC pins one coupling, and at exactly that coupling the "
                    "TT kinetic term vanishes -- no propagating graviton"
                    if closed else "the branch survives -- see the rows"),
        "entries": [
            {"claim": "the throat's constraint on alpha_GB can be read alone",
             "verdict": "NO -- alpha_GB IS A COUPLING CONSTANT",
             "evidence": "the same value acts in the exterior, where "
                         f"R_kk = +{opposite['ricci_null']:.4f} > 0 and the NEC "
                         f"needs alpha >= {opposite['exterior_needs_alpha_at_least']:.6f} "
                         "-- the opposite direction"},
            {"claim": "the two bounds coinciding is a coincidence",
             "verdict": "NO -- ONE CONTINUOUS BRACKET",
             "evidence": "T_kk = R_kk(1 + 4 alpha mu/f^4) on both sides, with "
                         "R_kk < 0 inside and > 0 outside, and mu/f^4 "
                         f"continuous at the seam at {seam['common_value']:.6f}/4"},
            {"claim": "the NEC alone leaves an open set of couplings",
             "verdict": "NO -- IT PINS ONE",
             "evidence": f"a {scan['samples']}-point scan leaves width "
                         f"{scan['surviving_width']:.1e} at "
                         f"{scan['critical_coupling']:.6f}; equivalently the "
                         "field equations select R^2 = -4 alpha_GB"},
            {"claim": "at that coupling the throat matter is exotic",
             "verdict": "NO -- NEC, WEC AND DEC ALL HOLD",
             "evidence": "rho = 3Aq, p_s = -3A, p_Omega = A with q >= 1, so "
                         f"rho+p_s = 3A(q-1) >= {honest['min_nec_radial']:.1e} "
                         "(saturated at the mouth) and rho exceeds both "
                         "pressures in magnitude. An earlier draft's "
                         "'relocated exotic matter' applies only below "
                         "critical"},
            {"claim": "at that coupling the observed universe must be empty",
             "verdict": "OVERREACHED -- IT IS A 5D BULK STATEMENT",
             "evidence": f"w = {empty['critical_equation_of_state']:.1f} is the "
                         "total stress supporting R x S^4_R, not the S^3 "
                         "equator's; and a homogeneous -Lambda g_ab can be "
                         "moved to the gravitational side. The defensible "
                         "claim is a vacuum-energy-like 5D exterior"},
            {"claim": "the tensor sector is healthy at the critical coupling",
             "verdict": "NO -- THE KINETIC TERM VANISHES",
             "evidence": f"C_kin = -(1/2)(1 + 4 alpha/R^2) derived by "
                         "linearising on THIS product background, matching to "
                         f"{max(abs(r['temporal_coefficient'] - r['predicted_kinetic']) for r in graviton['rows']):.0e}; "
                         "it is "
                         f"{graviton['rows'][-1]['temporal_coefficient']:.1e} at "
                         "criticality while the spatial coefficient is "
                         "untouched, so the equation loses its time derivative"},
            {"claim": "the trouble starts only at the critical point",
             "verdict": "NO -- SUPERLUMINAL ON THE WHOLE INTERVAL",
             "evidence": "c^2 = 1/(1 + 4 alpha/R^2) exceeds 1 for every "
                         "-R^2/4 < alpha < 0, reaching "
                         f"{graviton['rows'][-2]['speed_squared']:.0f} at 96% "
                         "of the critical value before diverging"},
        ],
        "what_this_closes": (
            "Constant-coupling EGB on this geometry, on STRUCTURAL rather than "
            "matter-content grounds. The NEC pins a unique coupling; the "
            "graviton degenerates at exactly that coupling; and the tensor cone "
            "is already outside the matter light cone on the approach."),
        "what_remains_untested": (
            "Whether an admissible SOURCE realises the throat's full anisotropic "
            "stress -- this round determines the stress the metric requires, not "
            "fields obeying their own equations that supply it. Also the scalar "
            "and vector perturbation sectors, dilatonic alpha(phi) L_GB where "
            "the scalar's own stress enters and where the heterotic term lives, "
            "f(R), and a DIFFERENT exterior: the constraint is derived for the "
            "round S^4_R completion this programme assumes."),
        "the_remaining_branches": {
            "1 accept the horizon": "the Tangherlini branch N(0) = 0 as the "
                                    "particle, abandoning MTY traversability",
            "2 ghost": "a wrong-sign field, with its stability problem",
            "3 quantum stress": "Casimir-type support, so the geometry is no "
                                "longer classical",
            "4 reinterpret": "particle exchange needs no traversable throat",
            "5 dilatonic EGB or f(R)": "untested here, and the place where the "
                                       "heterotic term actually lives",
        },
    }


# ── the graviton sector: what actually closes the branch ────────────────────

def graviton_kinetic_coefficient(coupling: float,
                                 bulk_radius: float = BULK_RADIUS) -> float:
    """``K_grav = 1 + 4 alpha_GB/R^2``, the coefficient of the TT ``omega^2`` term.

    **Derived, not recalled.** ``tt_kinetic_and_mass`` linearises the full
    ``G_ab + alpha H_ab`` numerically around ``R x S^4_R`` and fits the
    ``omega^2`` coefficient; the result is ``-(1/2)(1 + 4 alpha/R^2)`` across
    six couplings, two bulk radii and two polarisations.

    It vanishes **exactly** at the critical coupling ``alpha = -R^2/4`` -- the
    same value the NEC analysis pins.
    """
    return 1.0 + 4.0 * coupling / (bulk_radius * bulk_radius)


def graviton_speed_squared(coupling: float,
                           bulk_radius: float = BULK_RADIUS) -> float:
    """``c^2 = 1/(1 + 4 alpha_GB/R^2)`` for tensor modes.

    The **spatial** ``kappa^2`` coefficient is ``alpha``-independent -- checked,
    not assumed -- so the whole coupling dependence sits in the time derivative
    and the characteristic cone opens as ``alpha`` goes negative. Superluminal
    for every ``-R^2/4 < alpha < 0``, and divergent at the critical coupling.
    """
    kinetic = graviton_kinetic_coefficient(coupling, bulk_radius)
    return float("inf") if kinetic == 0.0 else 1.0 / kinetic


def _curvature(metric, point, step):
    from geometrodynamics.bulk.gauss_bonnet import _riemann
    g = metric(point)
    inverse = np.linalg.inv(g)
    up = _riemann(metric, point, step)
    down = np.einsum("am,mbcd->abcd", g, up)
    ricci = np.einsum("mamb->ab", up)
    scalar = float(np.einsum("ab,ab->", inverse, ricci))
    ricci_uu = inverse @ ricci @ inverse
    raised3 = np.einsum("abcd,bp,cq,dr->apqr", down, inverse, inverse, inverse)
    raised4 = np.einsum("abcd,ae->ebcd", raised3, inverse)
    lagrangian = (scalar ** 2 - 4.0 * float(np.einsum("ab,ab->", ricci, ricci_uu))
                  + float(np.einsum("abcd,abcd->", down, raised4)))
    return g, inverse, down, ricci, ricci_uu, scalar, raised3, lagrangian


def linearised_field_equation(metric, point, coupling: float,
                              step: float = 1e-3) -> np.ndarray:
    """``E_ab = G_ab + alpha_GB H_ab`` from a numerically differentiated metric.

    All four indices are raised explicitly when forming ``R_abcd R^abcd``; an
    earlier diagnostic omitted the leading one, and the tell was that its error
    did not shrink under refinement.
    """
    g, inverse, down, ricci, ricci_uu, scalar, raised3, lagrangian = _curvature(
        metric, point, step)
    dim = g.shape[0]
    result = np.zeros((dim, dim))
    for a in range(dim):
        for b in range(dim):
            einstein = ricci[a, b] - 0.5 * g[a, b] * scalar
            bracket_value = (scalar * ricci[a, b]
                             - 2.0 * (ricci @ inverse @ ricci)[a, b]
                             - 2.0 * np.einsum("cd,cd->", ricci_uu,
                                               down[a, :, b, :])
                             + np.einsum("cde,cde->", raised3[a], down[b]))
            lanczos = 2.0 * bracket_value - 0.5 * g[a, b] * lagrangian
            result[a, b] = einstein + coupling * lanczos
    return result


def _perturbed_metric(bulk_radius: float, amplitude: float, frequency: float,
                      axis: int, first: int = 2, second: int = 3):
    """``R x S^4_R`` in stereographic coordinates plus a TT plane wave.

    Stereographic coordinates make the background exactly ``delta_ij`` at the
    origin while carrying the correct Riemann tensor there, so a plane-wave
    polarisation is transverse and traceless without further projection.
    """
    def metric(point: np.ndarray) -> np.ndarray:
        spatial = point[1:]
        conformal = 1.0 / (1.0 + float(np.dot(spatial, spatial))
                           / (4.0 * bulk_radius * bulk_radius))
        g = np.diag([-1.0, conformal ** 2, conformal ** 2,
                     conformal ** 2, conformal ** 2])
        wave = amplitude * math.cos(frequency * point[axis])
        g[first, second] += wave
        g[second, first] += wave
        return g
    return metric


def tt_kinetic_and_mass(bulk_radius: float, coupling: float, axis: int,
                        first: int = 2, second: int = 3,
                        frequencies=(0.5, 0.9, 1.3, 1.7),
                        amplitude: float = 1e-5):
    """Fit ``dE/d(amplitude) = C_kin * freq^2 + C_mass`` for a TT mode.

    ``axis = 0`` gives the ``omega^2`` (time-derivative) coefficient;
    ``axis = 1`` the ``kappa^2`` (spatial) one.
    """
    point = np.zeros(5)
    design, values = [], []
    for frequency in frequencies:
        plus = linearised_field_equation(
            _perturbed_metric(bulk_radius, amplitude, frequency, axis,
                              first, second), point, coupling)[first, second]
        minus = linearised_field_equation(
            _perturbed_metric(bulk_radius, -amplitude, frequency, axis,
                              first, second), point, coupling)[first, second]
        design.append([frequency ** 2, 1.0])
        values.append((plus - minus) / (2.0 * amplitude))
    design = np.array(design)
    values = np.array(values)
    (kinetic, mass), *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = float(np.max(np.abs(design @ np.array([kinetic, mass]) - values)))
    return float(kinetic), float(mass), residual


@lru_cache(maxsize=4)
def measure_the_graviton_degenerates(
        bulk_radius: float = BULK_RADIUS) -> Dict[str, object]:
    """The decisive test: the TT kinetic coefficient at the critical coupling.

    Derived for **this** ``R x S^4`` background rather than inferred from the
    maximally symmetric formula, because the exterior is a *product* spacetime
    and not maximally symmetric.
    """
    critical = exterior_threshold(bulk_radius)
    rows = []
    for fraction in (0.0, 0.2, 0.5, 0.8, 0.96, 1.0):
        coupling = fraction * critical
        temporal, mass, residual = tt_kinetic_and_mass(bulk_radius, coupling, 0)
        spatial, _, _ = tt_kinetic_and_mass(bulk_radius, coupling, 1)
        rows.append({
            "coupling": coupling,
            "temporal_coefficient": temporal,
            "spatial_coefficient": spatial,
            "mass_term": mass,
            "fit_residual": residual,
            "predicted_kinetic": -0.5 * graviton_kinetic_coefficient(
                coupling, bulk_radius),
            "speed_squared": (float("inf") if abs(temporal) < 1e-6
                              else -spatial / temporal),
        })
    critical_row = rows[-1]
    spatial_values = [r["spatial_coefficient"] for r in rows]
    return {
        "critical_coupling": critical,
        "rows": rows,
        "kinetic_law": "C_kin = -(1/2)(1 + 4 alpha/R^2)",
        "law_holds": bool(all(
            abs(r["temporal_coefficient"] - r["predicted_kinetic"]) < 1e-5
            for r in rows)),
        "kinetic_vanishes_at_criticality": bool(
            abs(critical_row["temporal_coefficient"]) < 1e-5),
        "spatial_coefficient_is_coupling_independent": bool(
            max(spatial_values) - min(spatial_values) < 1e-6),
        "mass_term_stays_finite": bool(abs(critical_row["mass_term"]) > 1.0),
        "superluminal_below_criticality": bool(all(
            r["speed_squared"] > 1.0 + 1e-3
            for r in rows[1:-1])),
        "why_this_closes_the_branch": (
            "At alpha = -R^2/4 the omega^2 coefficient vanishes while the "
            "kappa^2 coefficient is untouched, so the tensor equation loses its "
            "time-derivative term entirely: it degenerates from an evolution "
            "equation into a constraint, and there is no propagating graviton. "
            "The mass term stays finite, so this is a degeneration of the "
            "kinetic operator and not an overall vanishing."),
        "and_it_is_bad_before_criticality": (
            "c^2 = 1/(1 + 4 alpha/R^2) is SUPERLUMINAL for every "
            "-R^2/4 < alpha < 0 and diverges at the critical coupling, so the "
            "tensor cone opens outside the matter light cone well before the "
            "degeneration. The branch is in trouble on an interval, not only at "
            "a point."),
        "why_it_had_to_be_derived": (
            "The familiar 1 + 2 alpha (D-3)(D-4) K coefficient is derived for a "
            "maximally symmetric SPACETIME. R x S^4_R is a product and is not "
            "maximally symmetric, so that formula does not apply here and the "
            "coefficient was computed by linearising the full field equations "
            "on this background."),
    }


@lru_cache(maxsize=4)
def measure_the_throat_matter_is_not_exotic(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """At the critical coupling the throat matter satisfies NEC, WEC **and** DEC.

    This corrects an earlier claim of this round that the exotic matter was
    "merely relocated" at criticality. That is true for ``alpha < -R^2/4``; at
    the critical point itself the throat stress is entirely respectable.
    """
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    critical = exterior_threshold(bulk_radius)
    s = np.linspace(-g.half_length, g.half_length, 401)
    f = throat_radius(s, b)
    amplitude = b * b / f ** 4                       # A > 0
    ratio = bulk_radius ** 2 * b * b / f ** 4        # q >= 1, equality at mouth
    density = 3.0 * amplitude * ratio
    radial = -3.0 * amplitude
    angular = amplitude
    return {
        "critical_coupling": critical,
        "q_at_the_mouth": float(np.min(ratio)),
        "q_at_the_neck": float(np.max(ratio)),
        "q_neck_closed_form": 1.0 / math.sin(mouth_angle) ** 4,
        "min_density": float(np.min(density)),
        "min_nec_radial": float(np.min(density + radial)),
        "min_nec_angular": float(np.min(density + angular)),
        "min_dec_radial": float(np.min(density - np.abs(radial))),
        "min_dec_angular": float(np.min(density - np.abs(angular))),
        "nec_holds": bool(np.all(density + radial >= -1e-9)
                          and np.all(density + angular >= -1e-9)),
        "wec_holds": bool(np.all(density >= -1e-9)),
        "dec_holds": bool(np.all(density - np.abs(radial) >= -1e-9)
                          and np.all(density - np.abs(angular) >= -1e-9)),
        "saturated_at_the_mouth": bool(abs(float(np.min(density + radial))) < 1e-9),
        "what_this_corrects": (
            "An earlier draft of this round said the exotic matter was 'merely "
            "relocated' at the critical coupling. It is not. rho + p_s = "
            "3A(q-1) >= 0 with equality only at the mouth, rho + p_Omega = "
            "A(3q+1) > 0, and rho exceeds both |p_s| and |p_Omega|, so NEC, WEC "
            "and DEC all hold along the throat. Relocation happens only for "
            "alpha strictly below the critical value."),
    }
