"""Does Gauss–Bonnet reopen the finite-mouth throat? Narrowed, not closed.

THE BRANCH, AND WHY IT WAS WORTH COMPUTING
──────────────────────────────────────────
``source_audit`` closed with five open branches, each a *premise* of the
flare-out theorem rather than a consequence. Gauss-Bonnet is the one to compute
first: in ``D = 5`` the invariant is **dynamical** rather than topological, it
is the only escape that keeps both a classical geometry and a traversable
throat, and it is a calculation instead of a philosophical choice. In ``D = 5``
it is also the *last* Lovelock term — see ``lovelock_status``.

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

Satisfying the matter NEC **at the neck** therefore needs

    alpha_GB <= -f_0^2/4

— a *negative* coupling. That does not close the branch, and an earlier draft
of this module wrongly said it did. Two corrections, both from review:

**Lovelock terminates at Gauss-Bonnet in D = 5.** The k-th Lovelock density
antisymmetrises 2k indices, so it vanishes identically for 2k > D; cubic
Lovelock needs six and is identically zero here. There is therefore no "rest of
the tower" to invoke, and ``alpha_GB H_kk/R_kk = -1`` does **not** mean exact
constant-coupling EGB has left its own regime — EGB is a complete classical
theory with second-order equations, not a truncation of one. It means only that
the Gauss-Bonnet term is order-one relative to Einstein.

**The heterotic sign is a comparison, not a closure argument.** In heterotic
string theory the Gauss-Bonnet term is dilaton-coupled, which this module
explicitly excludes. Its positive sign is suggestive about where a
string-motivated coupling would sit; it does not constrain constant EGB.

WHAT ACTUALLY NARROWS THE BRANCH
────────────────────────────────
A neck-only cancellation is not a wormhole. Since
``8 pi G_5 T_kk = R_kk (1 + 4 alpha mu/f^4)`` with ``R_kk < 0``, the NEC needs

    alpha_GB <= -f^4/(4 mu) = -f^4/(4 b^2)      pointwise

which is *weakest* at the neck and strengthens outward as ``f`` grows. Imposing
it across the whole throat is controlled by the mouth, and there

    f_m^4/(4 b^2) = (R sin a)^4 / (4 R^2 sin^4 a) = R^2/4      exactly

independently of the mouth angle. So a **global** matter-NEC solution needs

    alpha_GB <= -R^2/4 ,   i.e.   sqrt|alpha_GB| >= R/2

the Gauss-Bonnet length must be half the **bulk radius** — the size of the
universe — rather than a short-distance scale. At ``a = 0.3`` the neck-only
value ``-b^2/4`` leaves ``T_kk`` as negative as ``-98.3`` elsewhere on the
throat, missing the global requirement by ``1/sin^4 a = 131``.

That is the honest verdict: **narrowed, not closed.** Positive/string-sign
``alpha`` reinforces the violation; a sufficiently negative constant coupling
can satisfy the matter NEC, but only at a Gauss-Bonnet length comparable to the
whole closed universe. Global existence and stability of such a solution are
untouched here.

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
Constant-coupling Einstein-Gauss-Bonnet only, and it is **narrowed rather than
refuted**: negative ``alpha_GB <= -R^2/4`` does satisfy the matter NEC along
the throat, with global existence and stability untouched. A *dilatonic*
``alpha(phi)L_GB`` — where the scalar's own stress enters, and where the
heterotic term actually lives — and ``f(R)`` are separate premises and are not
tested. The Lovelock tower is **not** a separate premise in ``D = 5``: it
terminates at Gauss-Bonnet.

Predictions G1-G5 are frozen in ``docs/gauss_bonnet_prereg.md``, committed
before this module existed. G6 was frozen too and is **withdrawn** — see
``measure_the_negative_coupling_requirement``.
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
    "measure_the_negative_coupling_requirement",
    "lovelock_status",
    "global_coupling_threshold",
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


def lanczos_mixed(metric: Callable, point: np.ndarray,
                  step: float = 1e-3) -> Tuple[np.ndarray, float]:
    """The full ``H^a_b`` and ``L_GB``, keeping the ``-(1/2) g_ab L_GB`` term.

    ``lanczos_null_contraction`` may drop that term because it dies under a
    null contraction. Here it does not: the spatial block is exactly where it
    survives, so it has to be carried.

    The index placement is the delicate part. ``R_abcd R^abcd`` needs **all
    four** indices raised; a draft of this calculation raised only three and
    contracted a down index against a down index. The tell was that the error
    did not shrink under refinement — see the note in
    ``docs/negative_egb_prereg.md``.
    """
    g = metric(point)
    inverse = np.linalg.inv(g)
    riemann_up = _riemann(metric, point, step)
    riemann = np.einsum("am,mbcd->abcd", g, riemann_up)          # R_abcd
    ricci = np.einsum("mamb->ab", riemann_up)                    # R_ab
    scalar = float(np.einsum("ab,ab->", inverse, ricci))
    ricci_uu = inverse @ ricci @ inverse                         # R^ab
    riemann_u3 = np.einsum("abcd,bp,cq,dr->apqr", riemann,
                           inverse, inverse, inverse)            # R_a^{pqr}
    riemann_u4 = np.einsum("apqr,am->mpqr", riemann_u3, inverse)  # R^{mpqr}

    l_gb = (scalar ** 2
            - 4.0 * float(np.einsum("ab,ab->", ricci_uu, ricci))
            + float(np.einsum("abcd,abcd->", riemann, riemann_u4)))

    term1 = scalar * ricci
    term2 = np.einsum("ac,cb->ab", ricci @ inverse, ricci)
    term3 = np.einsum("cd,acbd->ab", ricci_uu, riemann)
    term4 = np.einsum("acde,bcde->ab", riemann, riemann_u3)
    lanczos = (2.0 * (term1 - 2.0 * term2 - 2.0 * term3 + term4)
               - 0.5 * g * l_gb)
    return inverse @ lanczos, l_gb


def _ultrastatic_lumpy() -> Callable:
    """``-dt^2 + h_4`` with a deliberately generic, inhomogeneous, non-diagonal
    spatial metric — no Killing vector, no symmetry to lean on."""
    def metric(point: np.ndarray) -> np.ndarray:
        x, y, z, w = point[1], point[2], point[3], point[4]
        g = np.eye(5)
        g[0, 0] = -1.0
        g[1, 1] = 1.0 + 0.30 * math.sin(x) + 0.10 * math.cos(2.0 * y)
        g[2, 2] = 1.0 + 0.20 * math.cos(y) * math.sin(z)
        g[3, 3] = 1.0 + 0.25 * math.sin(z + 0.4 * w)
        g[4, 4] = 1.0 + 0.15 * math.cos(w) * math.cos(x)
        g[1, 2] = g[2, 1] = 0.12 * math.sin(x + y) * math.cos(z)
        g[3, 4] = g[4, 3] = 0.08 * math.sin(w - z)
        return g
    return metric


def _lumpy_with_lapse() -> Callable:
    """The control: the same spatial lumps, but a **nonconstant lapse**, so the
    product structure is broken. If the spatial block vanished here too, the
    check could not detect anything."""
    base = _ultrastatic_lumpy()

    def metric(point: np.ndarray) -> np.ndarray:
        g = base(point).copy()
        g[0, 0] = -(1.0 + 0.3 * math.sin(point[1])) ** 2
        return g
    return metric


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


@lru_cache(maxsize=4)
def measure_the_spatial_block_vanishes() -> Dict[str, object]:
    """B0b — ``H^i_j = 0`` for **any** ultrastatic product, not just a
    maximally symmetric slice.

    Three earlier docstrings in this programme attributed the vanishing spatial
    Lanczos block to maximal symmetry of the ``S^4_R`` slice. That is the right
    value for the wrong reason, and it understates the result. For any
    ``-dt^2 + h_4`` in ``D = 5`` the spatial block of ``H^a_b`` reduces to the
    **four-dimensional Gauss-Bonnet (Euler) tensor of ``h_4``**, which vanishes
    identically because Gauss-Bonnet is topological in ``D = 4``. No symmetry
    is used.

    The consequence is not cosmetic: Gauss-Bonnet cannot touch the pressures
    *anywhere* on this geometry — throat as well as exterior — so the whole
    correction lands in the density, and the throat's ``p_s``, ``p_Omega`` are
    the Einstein ones at every coupling.
    """
    rows = []
    steps = (2e-3, 1e-3, 5e-4)
    cases = [
        ("throat, s = 0.9 (ultrastatic, NOT max. symmetric)",
         _throat_metric(0.7), np.array([0.0, 0.9, 1.0, 0.8, 0.0]), True),
        ("throat, at the neck s = 0",
         _throat_metric(0.7), np.array([0.0, 0.0, 1.1, 0.7, 0.0]), True),
        ("generic lumpy 4-slice (no symmetry at all)",
         _ultrastatic_lumpy(), np.array([0.0, 0.3, -0.5, 0.9, 0.2]), True),
        ("CONTROL: same slice, nonconstant lapse",
         _lumpy_with_lapse(), np.array([0.0, 0.3, -0.5, 0.9, 0.2]), False),
    ]
    for name, metric, point, ultrastatic in cases:
        residuals = []
        for step in steps:
            mixed, _ = lanczos_mixed(metric, point, step)
            residuals.append(float(np.abs(mixed[1:, 1:]).max()))
        # A discretisation zero shrinks like step^2; a structural nonzero does
        # not shrink at all.
        shrinks = bool(residuals[0] > 4.0 * residuals[-1])
        rows.append({"metric": name, "ultrastatic": ultrastatic,
                     "residuals": residuals,
                     "refinement_ratio": float(residuals[0] / residuals[-1])
                     if residuals[-1] > 0.0 else float("inf"),
                     "converges_to_zero": shrinks})

    ultrastatic_rows = [r for r in rows if r["ultrastatic"]]
    control = next(r for r in rows if not r["ultrastatic"])
    return {
        "rows": rows,
        "vanishes_for_every_ultrastatic_case": bool(
            all(r["converges_to_zero"] for r in ultrastatic_rows)),
        "control_does_not_vanish": bool(not control["converges_to_zero"]),
        "worst_ultrastatic_residual": max(r["residuals"][-1]
                                          for r in ultrastatic_rows),
        "control_residual": control["residuals"][-1],
        "why": ("For -dt^2 + h_4 in D = 5 the spatial block of H^a_b is the 4D "
                "Gauss-Bonnet (Euler) tensor of h_4, which vanishes "
                "identically because Gauss-Bonnet is topological in D = 4. "
                "Maximal symmetry is not used, and the throat slice is not "
                "maximally symmetric."),
        "why_the_control_matters": (
            "Breaking the product with a nonconstant lapse leaves a spatial "
            f"block of {control['residuals'][-1]:.2e} that does NOT shrink "
            "under refinement, so the zeros above are a result rather than a "
            "property of the instrument."),
        "what_it_means": (
            "Gauss-Bonnet cannot touch the pressures anywhere on this "
            "geometry. The entire alpha-dependence lands in the density, in "
            "the throat as well as in the exterior."),
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


def lovelock_status(order: int, dimension: int) -> str:
    """Whether the ``k``-th Lovelock density is dynamical, topological, or zero.

    The density antisymmetrises ``2k`` indices, so it vanishes identically for
    ``2k > D`` and is a total derivative at ``2k = D``. In ``D = 5`` that makes
    Gauss-Bonnet (``k = 2``) the **last** dynamical term: cubic Lovelock needs
    six indices and is identically zero. There is no further tower to invoke.
    """
    if 2 * order > dimension:
        return "identically zero"
    if 2 * order == dimension:
        return "topological"
    return "dynamical"


def global_coupling_threshold(bulk_radius: float = BULK_RADIUS) -> float:
    """``alpha_GB <= -R^2/4``: the NEC along the *whole* throat, not just the neck.

    Pointwise the requirement is ``alpha <= -f^4/(4 mu)``, weakest at the neck
    and strengthening outward. At the mouth
    ``f_m^4/(4b^2) = (R sin a)^4/(4 R^2 sin^4 a) = R^2/4`` **exactly**,
    independently of the mouth angle — so the binding constraint is set by the
    bulk radius alone.
    """
    return -0.25 * bulk_radius * bulk_radius


@lru_cache(maxsize=8)
def measure_the_negative_coupling_requirement(
        bulk_radius: float = BULK_RADIUS,
        mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """B3 — a neck-only cancellation is not a wormhole.

    Replaces an earlier draft's claim that ``alpha H_kk/R_kk = -1`` meant the
    derivative expansion had broken down. It does not: in ``D = 5`` Lovelock
    terminates at Gauss-Bonnet, and exact EGB is a complete classical theory
    rather than a truncation. What actually constrains the branch is that the
    NEC must hold along the throat, not only at its narrowest point.
    """
    g = geometry(bulk_radius, mouth_angle)
    b = g.neck_radius
    neck = coupling_threshold(b)
    glob = global_coupling_threshold(bulk_radius)
    s = np.linspace(-g.half_length, g.half_length, 601)
    with_neck = matter_null_stress(s, b, neck)
    with_global = matter_null_stress(s, b, glob)
    lovelock = [{"order": k, "status": lovelock_status(k, 5)}
                for k in (1, 2, 3, 4)]
    return {
        "neck_threshold": neck,
        "global_threshold": glob,
        "global_threshold_is_minus_quarter_R_squared": bool(
            abs(glob + 0.25 * bulk_radius ** 2) < 1e-15),
        "ratio_global_to_neck": glob / neck,
        "expected_ratio": 1.0 / math.sin(mouth_angle) ** 4,
        "neck_only_min_over_throat": float(np.min(with_neck)),
        "neck_only_satisfies_nec_globally": bool(np.all(with_neck >= -1e-9)),
        "global_min_over_throat": float(np.min(with_global)),
        "global_satisfies_nec": bool(np.all(with_global >= -1e-9)),
        "gauss_bonnet_length": math.sqrt(abs(glob)),
        "bulk_radius": bulk_radius,
        "length_over_bulk_radius": math.sqrt(abs(glob)) / bulk_radius,
        "relative_size_at_neck_threshold": relative_gauss_bonnet_size(neck, b),
        "lovelock_in_five_dimensions": lovelock,
        "tower_terminates_at_gauss_bonnet": bool(
            lovelock_status(3, 5) == "identically zero"),
        "what_was_withdrawn": (
            "An earlier draft read alpha_GB H_kk/R_kk = -1 as 'the derivative "
            "expansion has broken down and the whole Lovelock tower "
            "contributes'. Both halves are wrong in D = 5: cubic Lovelock "
            "antisymmetrises six indices and is IDENTICALLY ZERO here, so there "
            "is no further tower; and exact Einstein-Gauss-Bonnet is a complete "
            "classical theory with second-order equations, not a truncation "
            "that can lose validity. Order-one is just order-one."),
        "what_replaces_it": (
            "The NEC has to hold along the throat, not only at the neck. "
            "Pointwise it demands alpha <= -f^4/(4 mu), which is WEAKEST at the "
            "neck; at the mouth f_m^4/(4b^2) = R^2/4 exactly, independent of "
            "the mouth angle. So a global solution needs alpha <= -R^2/4, i.e. "
            "a Gauss-Bonnet length of at least R/2 -- half the radius of the "
            "closed universe, not a short-distance scale."),
        "the_honest_verdict": (
            "NARROWED, NOT CLOSED. Positive/string-sign alpha reinforces the "
            "violation. A sufficiently negative constant coupling can satisfy "
            "the matter NEC, but only at a Gauss-Bonnet length comparable to "
            "the whole closed universe. Global existence and stability of such "
            "a solution are untouched here."),
    }


@lru_cache(maxsize=4)
def measure_the_gauss_bonnet_ledger() -> Dict[str, object]:
    """B4 — the verdict, and what is left after it."""
    validation = measure_the_lanczos_tensor_is_correct()
    sign = measure_gauss_bonnet_reinforces_einstein()
    coupling = measure_the_required_coupling()
    requirement = measure_the_negative_coupling_requirement()
    narrowed = bool(sign["reinforces"] and coupling["no_positive_coupling_works"]
                    and not requirement["neck_only_satisfies_nec_globally"])
    return {
        "branch_is_narrowed_not_closed": narrowed,
        "verdict": ("Gauss-Bonnet is NARROWED, not closed: positive/string-sign "
                    "alpha reinforces the violation, and a global matter-NEC "
                    "solution needs a Gauss-Bonnet length of at least half the "
                    "BULK radius" if narrowed
                    else "the branch is unconstrained -- see the rows"),
        "entries": [
            {"claim": "the Lanczos implementation is trustworthy",
             "verdict": "VALIDATED IN D = 4",
             "evidence": "H_ab vanishes identically for Schwarzschild AND for "
                         "two general non-vacuum A(r), so the zero is not an "
                         "artefact of R_ab = 0"},
            {"claim": "the vanishing spatial block H^i_j = 0 is due to the "
                      "slice being maximally symmetric",
             "verdict": "NO -- IT IS THE 4D EULER TENSOR",
             "evidence": "H^i_j vanishes (ratio 16 under a 4x step reduction) "
                         "for the throat slice and for a generic lumpy 4-slice "
                         "with no symmetry; a nonconstant-lapse control stays "
                         "at 2.3e-2 and does not shrink. So Gauss-Bonnet cannot "
                         "touch the pressures ANYWHERE on this geometry"},
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
            {"claim": "alpha H_kk/R_kk = -1 means the expansion has broken down",
             "verdict": "WITHDRAWN -- WRONG IN D = 5",
             "evidence": "cubic Lovelock antisymmetrises six indices and is "
                         "IDENTICALLY ZERO in five dimensions, so there is no "
                         "further tower; and exact EGB is a complete classical "
                         "theory, not a truncation that can lose validity"},
            {"claim": "a negative coupling that works at the neck suffices",
             "verdict": "NO -- THE NECK IS THE EASIEST POINT",
             "evidence": "the NEC needs alpha <= -f^4/(4 mu) pointwise, weakest "
                         f"at the neck; alpha = -b^2/4 leaves T_kk as low as "
                         f"{requirement['neck_only_min_over_throat']:.1f} "
                         "elsewhere on the throat"},
            {"claim": "a global matter-NEC solution is available cheaply",
             "verdict": "ONLY AT A COSMOLOGICAL GAUSS-BONNET LENGTH",
             "evidence": "the mouth sets f_m^4/(4b^2) = R^2/4 exactly, "
                         "independent of the mouth angle, so alpha <= -R^2/4 "
                         f"and sqrt|alpha| >= R/2 -- "
                         f"{requirement['ratio_global_to_neck']:.0f}x the "
                         "neck-only requirement here"},
        ],
        "what_this_narrows": (
            "The classical-geometry escape from the source audit. Gauss-Bonnet "
            "is the only branch keeping both a classical geometry and a "
            "traversable throat, and it is not closed: a sufficiently negative "
            "constant coupling does satisfy the matter NEC. What is closed is "
            "the cheap version -- the string-motivated sign fails, and a global "
            "solution needs the Gauss-Bonnet length to be half the radius of "
            "the closed universe rather than a short-distance scale."),
        "what_remains": {
            "0 negative-coupling EGB": "not refuted -- alpha <= -R^2/4 satisfies "
                                       "the matter NEC along the throat, at a "
                                       "Gauss-Bonnet length of half the bulk "
                                       "radius; global existence and stability "
                                       "are open",
            "1 accept the horizon": "the Tangherlini branch N(0) = 0 as the "
                                    "particle, abandoning MTY traversability",
            "2 ghost": "a wrong-sign field, with its stability problem",
            "3 quantum stress": "Casimir-type support, so the geometry is no "
                                "longer classical",
            "4 reinterpret": "particle exchange needs no traversable throat",
        },
        "not_refuted_here": (
            "Negative-coupling constant EGB itself, which this round NARROWS "
            "rather than closes. Also dilatonic alpha(phi) L_GB, where the "
            "scalar's own stress enters and known 5D solutions exist, and "
            "f(R). The Lovelock tower is not a separate premise in D = 5: it "
            "terminates at Gauss-Bonnet."),
    }
