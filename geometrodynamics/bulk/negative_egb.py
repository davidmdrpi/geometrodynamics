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

WHAT SURVIVES THERE IS NOT USABLE
─────────────────────────────────
``H^i_j = 0`` on a maximally symmetric slice, so Gauss-Bonnet does not touch the
exterior pressure at all:

    8 pi G_5 rho_ext = 6/R^2 + 12 alpha_GB/R^4      8 pi G_5 p_ext = -3/R^2

At the critical coupling that gives ``rho = 3/R^2``, ``p = -rho``, so
``w = -1`` **exactly**: the exterior is forced to be pure vacuum energy with no
ordinary matter anywhere. And pushing ``alpha_GB`` more negative makes the
*exterior* violate the NEC — relocating the exotic matter from the throat into
the universe around it rather than removing it.

SCOPE
─────
This settles **global existence** for constant-coupling EGB on this geometry. It
does not address stability or the graviton kinetic term at that coupling, and it
does not touch dilatonic ``alpha(phi) L_GB`` or ``f(R)``. Predictions E1-E6 are
frozen in ``docs/negative_egb_prereg.md``, committed before this module existed.
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

    Gauss-Bonnet cannot touch it: ``H^i_j = 0`` on a maximally symmetric
    spatial slice, so the whole correction lands in ``rho``.
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
            "H^i_j = 0 on a maximally symmetric spatial slice, so the entire "
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
    """The verdict on the branch PR #277 left open."""
    opposite = measure_the_exterior_constrains_alpha_oppositely()
    scan = measure_no_coupling_satisfies_both()
    seam = measure_the_bracket_is_continuous_at_the_seam()
    empty = measure_the_critical_exterior_is_empty()
    closed = bool(scan["the_surviving_set_is_measure_zero"]
                  and empty["it_is_exactly_vacuum_energy"])
    return {
        "branch_is_closed": closed,
        "verdict": ("Negative-coupling EGB closes on PHYSICAL grounds: it "
                    "survives at one exact coupling, and there the observed "
                    "universe must be empty of ordinary matter" if closed
                    else "the branch survives -- see the rows"),
        "entries": [
            {"claim": "the throat's constraint on alpha_GB can be read alone",
             "verdict": "NO -- alpha_GB IS A COUPLING CONSTANT",
             "evidence": "the same value acts in the exterior, where "
                         f"R_kk = +{opposite['ricci_null']:.4f} > 0 and the NEC "
                         f"needs alpha >= {opposite['exterior_needs_alpha_at_least']:.6f} "
                         "-- the opposite direction"},
            {"claim": "some open interval of alpha_GB works everywhere",
             "verdict": "NO -- THE SURVIVING SET IS MEASURE ZERO",
             "evidence": f"a {scan['samples']}-point scan over "
                         f"{scan['scanned_range']} leaves a set of width "
                         f"{scan['surviving_width']:.1e}, at the single value "
                         f"{scan['critical_coupling']:.6f}"},
            {"claim": "the two bounds coinciding is a coincidence",
             "verdict": "NO -- ONE CONTINUOUS BRACKET",
             "evidence": "T_kk = R_kk(1 + 4 alpha mu/f^4) on both sides, with "
                         "R_kk < 0 inside and > 0 outside, and mu/f^4 "
                         f"continuous at the seam at {seam['common_value']:.6f}/4"},
            {"claim": "Gauss-Bonnet can adjust the exterior pressure",
             "verdict": "NO -- H^i_j = 0",
             "evidence": "on a maximally symmetric slice the whole correction "
                         "lands in rho; p = -3/R^2 whatever alpha is"},
            {"claim": "the surviving coupling is physically usable",
             "verdict": "NO -- IT EMPTIES THE UNIVERSE",
             "evidence": f"w = {empty['critical_equation_of_state']:.1f} exactly, "
                         f"rho = {empty['critical_density']:.6f} = half the "
                         "Einstein value: pure vacuum energy, no ordinary "
                         "matter anywhere"},
        ],
        "what_this_closes": (
            "Global existence for constant-coupling EGB on this geometry. The "
            "previous round narrowed the branch by showing a negative coupling "
            "was needed; this one closes it by noting that the coupling is not "
            "the throat's to choose."),
        "what_remains_untested": (
            "Stability and the graviton kinetic term at the critical coupling "
            "-- the second half of PR #277's open item, and not addressed here. "
            "Also dilatonic alpha(phi) L_GB, where the scalar's own stress "
            "enters the null contraction and known 5D solutions exist, and "
            "f(R). And a DIFFERENT exterior: the constraint is derived for the "
            "round S^4_R completion this programme assumes, so a different bulk "
            "completion would need its own version of this calculation."),
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
