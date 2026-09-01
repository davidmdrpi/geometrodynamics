"""Does Gauss–Bonnet reopen the finite-mouth throat? No — it reinforces it.

THE BRANCH, AND WHY IT WAS WORTH COMPUTING
──────────────────────────────────────────
``source_audit`` closed with five open branches, each a *premise* of the
flare-out theorem rather than a consequence. Gauss-Bonnet is the one to compute
first: in ``D = 5`` the invariant is **dynamical** rather than topological, it
is the only escape that keeps both a classical geometry and a traversable
throat, and it is a calculation instead of a philosophical choice.

    G_ab + alpha_GB H_ab = 8 pi G_5 T_ab

    H_ab = 2[ R R_ab - 2 R_ac R^c_b - 2 R^{cd} R_acbd + R_a^{cde} R_bcde ]
           - (1/2) g_ab L_GB

Under null contraction ``g_ab L_GB`` drops and ``G_kk = R_kk``, so

    R_kk + alpha_GB H_kk = 8 pi G_5 T_kk

**Raychaudhuri does not care about the gravitational action.** Flare-out is
pure geometry, so a neck still forces ``R_kk = -3 f''/f``. Only the question of
*which tensor supplies it* changes. The branch's hope was that
``alpha_GB H_kk`` supplies the negative part geometrically, leaving the matter
NEC intact.

WHAT IS FOUND
─────────────
The opposite. At **any** neck, for **any** lapse,

    R_kk = -3 f''/f_0 ,    H_kk = -12 f''/f_0^3 ,    H_kk/R_kk = 4/f_0^2 > 0

Gauss-Bonnet has the *same sign* as Einstein and reinforces the violation. More
generally ``H_kk/R_kk = 4(1-f'^2)/f^2 = 4 mu/f^4`` with ``mu = f^2(1-f'^2)``
the Misner-Sharp parameter — the same quantity ``finite_mouth`` P2 showed is
continuous across the seam — so the reinforcement holds along the whole throat
wherever ``mu > 0``.

Satisfying the matter NEC therefore needs

    alpha_GB <= -f_0^2/4

which fails three ways at once: the sign is opposite to the heterotic value
``alpha' /8 > 0``; the magnitude ties the Gauss-Bonnet length to the throat
radius; and at that magnitude ``alpha_GB H_kk/R_kk = -1`` exactly, so the
"correction" equals the leading term and truncating Lovelock at Gauss-Bonnet is
no longer justified.

HOW IT IS CHECKED
─────────────────
``lanczos_null_contraction`` builds the full Riemann tensor by numerical
differentiation and contracts it, sharing no algebra with the closed forms. It
is validated where the answer is known: in ``D = 4`` the Gauss-Bonnet term is
topological, so ``H_ab`` must vanish identically — and it does, for a general
non-vacuum metric as well as for Schwarzschild, so the zero is not an artefact
of vacuum structure.

SCOPE
─────
Constant-coupling Einstein-Gauss-Bonnet only. A *dilatonic* ``alpha(phi)L_GB``,
the full Lovelock tower, and ``f(R)`` are separate premises and are **not**
refuted here. Predictions G1-G6 are frozen in ``docs/gauss_bonnet_prereg.md``,
committed before this module existed.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from geometrodynamics.bulk.finite_mouth import (
    BULK_RADIUS,
    MOUTH_ANGLE,
    geometry,
    misner_sharp,
    throat_radius,
)

__all__ = [
    "ricci_null_contraction",
    "lanczos_null_contraction_closed_form",
    "lanczos_null_contraction",
    "gauss_bonnet_ratio",
    "matter_null_stress",
    "coupling_threshold",
    "relative_gauss_bonnet_size",
    "HETEROTIC_SIGN",
    "measure_the_lanczos_tensor_is_correct",
    "measure_gauss_bonnet_reinforces_einstein",
    "measure_the_required_coupling",
    "measure_the_expansion_breaks_down",
    "measure_the_gauss_bonnet_ledger",
]

#: In heterotic string theory ``alpha_GB = alpha'/8 > 0``. The sign matters.
HETEROTIC_SIGN = +1


# ── the two null contractions, in closed form ───────────────────────────────

def ricci_null_contraction(s, neck_radius: float,
                           lapse=None, lapse_derivative=None, step: float = 1e-5):
    """``R_kk = -3(N f'' - N' f')/(N f)`` on the scalar-flat profile.

    At a neck ``f' = 0`` kills the lapse term, so ``R_kk = -3 f''/f`` whatever
    ``N'`` does — the same structure that made the P4 price lapse-independent.
    """
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    fp = s / f
    fpp = neck_radius * neck_radius / f ** 3
    if lapse is None:
        return -3.0 * fpp / f
    n = np.asarray(lapse(s, neck_radius), dtype=float)
    if lapse_derivative is not None:
        dn = np.asarray(lapse_derivative(s, neck_radius), dtype=float)
    else:
        dn = (np.asarray(lapse(s + step, neck_radius), dtype=float)
              - np.asarray(lapse(s - step, neck_radius), dtype=float)) / (2.0 * step)
    return -3.0 * (n * fpp - dn * fp) / (n * f)


def lanczos_null_contraction_closed_form(s, neck_radius: float,
                                         lapse=None, lapse_derivative=None,
                                         step: float = 1e-5):
    """``H_kk = 12(f'^2 - 1)(N f'' - N' f')/(N f^3)``.

    Equivalently ``H_kk = (4 mu / f^4) R_kk`` with ``mu = f^2(1 - f'^2)``.
    """
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    fp = s / f
    ricci = ricci_null_contraction(s, neck_radius, lapse, lapse_derivative, step)
    return 4.0 * (1.0 - fp * fp) / (f * f) * ricci


def gauss_bonnet_ratio(s, neck_radius: float):
    """``H_kk/R_kk = 4(1 - f'^2)/f^2 = 4 mu/f^4``, the reinforcement factor.

    Positive wherever the Misner-Sharp parameter is positive, which is the
    whole throat here since ``mu = b^2``. Sign-positive means Gauss-Bonnet adds
    to the Einstein violation instead of cancelling it.
    """
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    return 4.0 * misner_sharp(s, neck_radius) / f ** 4


def matter_null_stress(s, neck_radius: float, coupling: float,
                       lapse=None, lapse_derivative=None):
    """``8 pi G_5 T_kk = R_kk + alpha_GB H_kk``, what the matter must supply."""
    ricci = ricci_null_contraction(s, neck_radius, lapse, lapse_derivative)
    lanczos = lanczos_null_contraction_closed_form(
        s, neck_radius, lapse, lapse_derivative)
    return ricci + coupling * lanczos


def coupling_threshold(neck_radius: float) -> float:
    """``alpha_GB = -f_0^2/4``: the largest coupling that lets ``T_kk = 0``."""
    return -0.25 * neck_radius * neck_radius


def relative_gauss_bonnet_size(coupling: float, neck_radius: float) -> float:
    """``alpha_GB H_kk / R_kk = 4 alpha_GB / f_0^2`` at the neck.

    Exactly ``-1`` at the threshold: the correction equals the leading term,
    which is where a derivative expansion truncated at Gauss-Bonnet stops
    meaning anything.
    """
    return 4.0 * coupling / (neck_radius * neck_radius)


# ── an independent numerical Lanczos tensor ─────────────────────────────────

def _christoffel(metric: Callable, point: np.ndarray, step: float) -> np.ndarray:
    dim = point.size
    g = metric(point)
    inverse = np.linalg.inv(g)
    dg = np.zeros((dim, dim, dim))
    for c in range(dim):
        forward, backward = point.copy(), point.copy()
        forward[c] += step
        backward[c] -= step
        dg[c] = (metric(forward) - metric(backward)) / (2.0 * step)
    gamma = np.zeros((dim, dim, dim))
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                gamma[a, b, c] = 0.5 * sum(
                    inverse[a, d] * (dg[c][d, b] + dg[b][d, c] - dg[d][b, c])
                    for d in range(dim))
    return gamma


def _riemann(metric: Callable, point: np.ndarray, step: float) -> np.ndarray:
    dim = point.size
    gamma = _christoffel(metric, point, step)
    dgamma = np.zeros((dim, dim, dim, dim))
    for c in range(dim):
        forward, backward = point.copy(), point.copy()
        forward[c] += step
        backward[c] -= step
        dgamma[c] = (_christoffel(metric, forward, step)
                     - _christoffel(metric, backward, step)) / (2.0 * step)
    riemann = np.zeros((dim, dim, dim, dim))
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    value = dgamma[c][a, b, d] - dgamma[d][a, b, c]
                    for e in range(dim):
                        value += (gamma[a, c, e] * gamma[e, b, d]
                                  - gamma[a, d, e] * gamma[e, b, c])
                    riemann[a, b, c, d] = value
    return riemann


def lanczos_null_contraction(metric: Callable, point: np.ndarray,
                             null_vector: np.ndarray,
                             step: float = 1e-4) -> Tuple[float, float]:
    """``(R_kk, H_kk)`` built from a numerically differentiated Riemann tensor.

    Shares no algebra with the closed forms above, so agreement is a regression
    rather than an identity. The ``-(1/2) g_ab L_GB`` term is omitted because
    it drops under a null contraction; ``measure_the_lanczos_tensor_is_correct``
    checks the null condition rather than assuming it.
    """
    g = metric(point)
    inverse = np.linalg.inv(g)
    dim = point.size
    riemann_up = _riemann(metric, point, step)
    riemann = np.einsum("am,mbcd->abcd", g, riemann_up)
    ricci = np.einsum("mamb->ab", riemann_up)
    scalar = float(np.einsum("ab,ab->", inverse, ricci))

    k = null_vector
    ricci_kk = float(k @ ricci @ k)
    # 2[ R R_kk - 2 (R_ac k^a)(R^c_b k^b) - 2 R^{cd} R_acbd k^a k^b
    #    + R_acde R_b^{cde} k^a k^b ]
    rk_low = ricci @ k
    rk_up = inverse @ rk_low
    term2 = float(rk_low @ rk_up)
    ricci_uu = inverse @ ricci @ inverse
    term3 = float(np.einsum("cd,acbd,a,b->", ricci_uu, riemann, k, k))
    riemann_up3 = np.einsum("abcd,bp,cq,dr->apqr", riemann, inverse, inverse,
                            inverse)
    term4 = float(np.einsum("acde,bcde,a,b->", riemann, riemann_up3, k, k))
    lanczos_kk = 2.0 * (scalar * ricci_kk - 2.0 * term2 - 2.0 * term3 + term4)
    return ricci_kk, lanczos_kk


def _throat_metric(neck_radius: float, lapse_value: float = 1.0) -> Callable:
    """``-N^2 dt^2 + ds^2 + f^2 dOmega_3^2`` as a callable, for the checker."""
    def metric(point: np.ndarray) -> np.ndarray:
        s, a1, a2 = point[1], point[2], point[3]
        f = math.sqrt(s * s + neck_radius * neck_radius)
        return np.diag([-lapse_value ** 2, 1.0, f * f,
                        f * f * math.sin(a1) ** 2,
                        f * f * math.sin(a1) ** 2 * math.sin(a2) ** 2])
    return metric


def _schwarzschild4(mass: float) -> Callable:
    def metric(point: np.ndarray) -> np.ndarray:
        r, theta = point[1], point[2]
        lapse = 1.0 - 2.0 * mass / r
        return np.diag([-lapse, 1.0 / lapse, r * r,
                        r * r * math.sin(theta) ** 2])
    return metric


def _general4(profile: Callable) -> Callable:
    """A general non-vacuum static 4D metric, so the ``D = 4`` zero is not
    an artefact of vacuum structure."""
    def metric(point: np.ndarray) -> np.ndarray:
        r, theta = point[1], point[2]
        value = profile(r)
        return np.diag([-value, 1.0 / value, r * r,
                        r * r * math.sin(theta) ** 2])
    return metric


# ── measurements ────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def measure_the_lanczos_tensor_is_correct() -> Dict[str, object]:
    """B0 — validate the implementation where the answer is already known.

    In ``D = 4`` the Gauss-Bonnet term is topological, so ``H_ab`` vanishes
    identically. Checked for Schwarzschild *and* for a general non-vacuum
    profile, so a zero cannot be blamed on ``R_ab = 0``.
    """
    rows = []
    point4 = np.array([0.0, 3.0, 1.1, 0.0])
    for name, metric in (("4D Schwarzschild (vacuum)", _schwarzschild4(1.0)),
                         ("4D general A(r) = 1 - 2/r + 0.3/r^2 (non-vacuum)",
                          _general4(lambda r: 1.0 - 2.0 / r + 0.3 / r ** 2)),
                         ("4D general A(r) = 1 + 0.2 r^2 (non-vacuum)",
                          _general4(lambda r: 1.0 + 0.2 * r ** 2))):
        g = metric(point4)
        lapse = -g[0, 0]
        k = np.array([1.0 / lapse, 1.0, 0.0, 0.0])
        null_residual = float(k @ g @ k)
        ricci_kk, lanczos_kk = lanczos_null_contraction(metric, point4, k,
                                                        step=2e-4)
        rows.append({"metric": name, "null_residual": null_residual,
                     "ricci_kk": ricci_kk, "lanczos_kk": lanczos_kk,
                     "vanishes": bool(abs(lanczos_kk) < 1e-4)})
    return {
        "rows": rows,
        "topological_in_four_dimensions": bool(all(r["vanishes"] for r in rows)),
        "worst_null_residual": max(abs(r["null_residual"]) for r in rows),
        "why": ("In D = 4 the Gauss-Bonnet invariant is a total derivative, so "
                "its variation H_ab vanishes identically. Getting zero for a "
                "general NON-VACUUM A(r) as well as for Schwarzschild shows the "
                "zero is a property of the implementation and the dimension, "
                "not of R_ab = 0."),
    }


@lru_cache(maxsize=8)
def measure_gauss_bonnet_reinforces_einstein(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """B1 — the decisive sign: ``H_kk/R_kk = 4 mu/f^4 > 0``."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    s = np.linspace(-g.half_length, g.half_length, 401)
    ricci = ricci_null_contraction(s, b)
    lanczos = lanczos_null_contraction_closed_form(s, b)
    ratio = gauss_bonnet_ratio(s, b)
    metric = _throat_metric(b)
    checks = []
    for value in (0.0, 0.3 * g.half_length, 0.7 * g.half_length):
        point = np.array([0.0, value, 1.1, 0.9, 0.0])
        k = np.array([1.0, 1.0, 0.0, 0.0, 0.0])
        numeric_r, numeric_h = lanczos_null_contraction(metric, point, k,
                                                        step=1e-3)
        closed_r = float(ricci_null_contraction(np.array([value]), b)[0])
        closed_h = float(lanczos_null_contraction_closed_form(
            np.array([value]), b)[0])
        checks.append({
            "s": value,
            "ricci_closed": closed_r, "ricci_numeric": numeric_r,
            "lanczos_closed": closed_h, "lanczos_numeric": numeric_h,
            "ricci_relative_error": abs(numeric_r - closed_r) / abs(closed_r),
            "lanczos_relative_error": abs(numeric_h - closed_h) / abs(closed_h)})
    mid = len(s) // 2
    return {
        "neck_radius": b,
        "ricci_at_neck": float(ricci[mid]),
        "lanczos_at_neck": float(lanczos[mid]),
        "ratio_at_neck": float(ratio[mid]),
        "expected_ratio_at_neck": 4.0 / (b * b),
        "ratio_is_everywhere_positive": bool(np.all(ratio > 0.0)),
        "misner_sharp_is_positive": bool(
            np.all(misner_sharp(s, b) > 0.0)),
        "independent_checks": checks,
        "worst_relative_error": max(
            max(c["ricci_relative_error"], c["lanczos_relative_error"])
            for c in checks),
        "reinforces": bool(float(ratio[mid]) > 0.0),
        "why_this_decides_the_branch": (
            "The branch was invoked so that alpha_GB H_kk could supply the "
            "negative part geometrically and leave T_kk >= 0. Instead H_kk has "
            "the SAME sign as R_kk at every neck, so for the heterotic sign "
            "alpha_GB > 0 the Gauss-Bonnet term makes the violation WORSE by "
            "the factor (1 + 4 alpha_GB/f_0^2)."),
        "the_ratio_is_the_misner_sharp_parameter": (
            "H_kk/R_kk = 4(1-f'^2)/f^2 = 4 mu/f^4 with mu = f^2(1-f'^2), the "
            "same quantity finite_mouth's P2 showed is continuous across the "
            "seam. Since mu = b^2 > 0 throughout, the reinforcement holds along "
            "the whole throat and not only at the neck."),
    }


@lru_cache(maxsize=8)
def measure_the_required_coupling(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """B2 — what coupling the matter NEC would demand, and its sign."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    threshold = coupling_threshold(b)
    rows = []
    for name, coupling in (("heterotic sign, alpha = +b^2/4", 0.25 * b * b),
                           ("alpha = 0 (pure Einstein)", 0.0),
                           ("alpha = -b^2/8 (half-way)", -0.125 * b * b),
                           ("threshold alpha = -b^2/4", threshold),
                           ("beyond, alpha = -b^2/2", -0.5 * b * b)):
        value = float(matter_null_stress(np.array([0.0]), b, coupling)[0])
        rows.append({"case": name, "coupling": coupling,
                     "matter_null_stress": value,
                     "nec_satisfied": bool(value >= -1e-9)})
    positive_ok = [r for r in rows if r["coupling"] > 0.0 and r["nec_satisfied"]]
    return {
        "neck_radius": b,
        "threshold": threshold,
        "threshold_formula": "alpha_GB <= -f_0^2/4 = -R^2 sin^4 a / 4",
        "rows": rows,
        "no_positive_coupling_works": bool(not positive_ok),
        "heterotic_sign": HETEROTIC_SIGN,
        "heterotic_makes_it_worse": bool(rows[0]["matter_null_stress"]
                                         < rows[1]["matter_null_stress"]),
        "why_the_sign_matters": (
            "Heterotic string theory gives alpha_GB = alpha'/8 > 0. That is "
            "exactly the sign that deepens the violation here, so the "
            "best-motivated value of the coupling is the one that fails "
            "hardest."),
    }


@lru_cache(maxsize=8)
def measure_the_expansion_breaks_down(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """B3 — even granting the wrong sign, the truncation is unjustified."""
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    threshold = coupling_threshold(b)
    rows = [{"coupling": coupling,
             "relative_size": relative_gauss_bonnet_size(coupling, b)}
            for coupling in (0.01 * threshold, 0.1 * threshold,
                             0.5 * threshold, threshold)]
    return {
        "threshold": threshold,
        "relative_size_at_threshold": relative_gauss_bonnet_size(threshold, b),
        "rows": rows,
        "it_equals_the_leading_term": bool(
            abs(relative_gauss_bonnet_size(threshold, b) + 1.0) < 1e-12),
        "gauss_bonnet_length": math.sqrt(abs(threshold)),
        "neck_radius": b,
        "length_ratio": math.sqrt(abs(threshold)) / b,
        "why": ("alpha_GB H_kk/R_kk = 4 alpha_GB/f_0^2, which is exactly -1 at "
                "the threshold: the 'correction' equals the term it corrects. "
                "A derivative expansion truncated at Gauss-Bonnet has stopped "
                "meaning anything there, and the whole Lovelock tower would "
                "contribute. The required Gauss-Bonnet length is f_0/2, i.e. "
                "tied to the throat radius rather than to a separate short "
                "scale."),
    }


@lru_cache(maxsize=4)
def measure_the_gauss_bonnet_ledger() -> Dict[str, object]:
    """B4 — the verdict, and what is left after it."""
    validation = measure_the_lanczos_tensor_is_correct()
    sign = measure_gauss_bonnet_reinforces_einstein()
    coupling = measure_the_required_coupling()
    expansion = measure_the_expansion_breaks_down()
    closed = bool(sign["reinforces"] and coupling["no_positive_coupling_works"]
                  and expansion["it_equals_the_leading_term"])
    return {
        "branch_is_closed": closed,
        "verdict": ("Gauss-Bonnet does NOT reopen the throat: it reinforces "
                    "the Einstein violation" if closed
                    else "the branch survives -- see the rows"),
        "entries": [
            {"claim": "the Lanczos implementation is trustworthy",
             "verdict": "VALIDATED IN D = 4",
             "evidence": "H_ab vanishes identically for Schwarzschild AND for "
                         "two general non-vacuum A(r), so the zero is not an "
                         "artefact of R_ab = 0"},
            {"claim": "Gauss-Bonnet can supply the negative null stress",
             "verdict": "NO -- IT HAS THE SAME SIGN",
             "evidence": f"H_kk/R_kk = 4 mu/f^4 = "
                         f"{sign['ratio_at_neck']:.4f} > 0 at the neck, "
                         "positive along the whole throat since mu = b^2 > 0"},
            {"claim": "the lapse can help in Gauss-Bonnet",
             "verdict": "NO -- IT DROPS OUT AT THE NECK",
             "evidence": "N' multiplies f', and f'(0) = 0 is what makes s = 0 "
                         "a neck; the same structure as P4"},
            {"claim": "the heterotic coupling alpha' /8 > 0 helps",
             "verdict": "NO -- IT MAKES IT WORSE",
             "evidence": "the best-motivated sign deepens the violation by "
                         "(1 + 4 alpha_GB/f_0^2); the matter NEC would need "
                         f"alpha_GB <= {coupling['threshold']:.9f}"},
            {"claim": "a negative coupling of the right size would work",
             "verdict": "ONLY OUTSIDE THE THEORY'S OWN REGIME",
             "evidence": "at threshold alpha_GB H_kk/R_kk = "
                         f"{expansion['relative_size_at_threshold']:.1f}, so "
                         "the correction equals the leading term and "
                         "truncating Lovelock at Gauss-Bonnet is unjustified"},
        ],
        "what_this_closes": (
            "The classical-geometry escape from the source audit. Gauss-Bonnet "
            "was the only branch keeping both a classical geometry and a "
            "traversable throat, and for the natural D = 5 higher-curvature "
            "term with its best-motivated sign it fails."),
        "what_remains": {
            "1 accept the horizon": "the Tangherlini branch N(0) = 0 as the "
                                    "particle, abandoning MTY traversability",
            "2 ghost": "a wrong-sign field, with its stability problem",
            "3 quantum stress": "Casimir-type support, so the geometry is no "
                                "longer classical",
            "4 reinterpret": "particle exchange needs no traversable throat",
        },
        "not_refuted_here": (
            "Dilatonic alpha(phi) L_GB, where the scalar's own stress enters "
            "and known 5D solutions exist; the full Lovelock tower; and f(R). "
            "Those are separate premises, and this round tested constant-"
            "coupling EGB only."),
    }
