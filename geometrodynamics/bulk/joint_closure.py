"""What data does joint closure retain? Two independently prepared triangles.

Pre-registered in ``docs/joint_closure_composition_prereg.md`` at ``b78157a``
(amendment A1), committed before this file. The question is whether the
existing BAM closure rules require the joint weight to depend on the two
triangle invariants separately, on their product alone, or on additional
joint invariants.

Scope, restated from the freeze. The reference preparation is product Haar
with separate phase windows in each triangle. Round 8 established that phase
conditioning is **chosen**, justified by ``history/closure.py:11``, not forced
by the zero set; conditioning on the numerators instead gives the uniform
measure. Nothing here promotes that choice to a geometric necessity.

``PRODUCT_STATISTIC_SUFFICIENT`` for the reference rule is a design
consequence stated in the freeze, not a discovery of this module. What is
open is whether the repository supplies any *other* physically required joint
rule, and whether the reference factorization transfers to other weights.
It does not: see :func:`level_set_controls`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


PUBLIC_PREREG = "b78157a"
SEED = 2026090710
SECTOR_SIGNS: Tuple[Tuple[int, int], ...] = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def _unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n == 0.0:
        raise ValueError("cannot normalise a zero or non-finite vector")
    return v / n


@dataclass(frozen=True)
class Triangle:
    """One prepared pair: settings ``a, b`` and sector signs ``sA, sB``.

    Round 9's singlet convention: ``u = sA a``, ``w = -sB b``. ``q = u x w``
    is the closure normal, ``t = 1 + u.w`` and ``s = u + w`` with
    ``|s|^2 = 2t``. On the closure circle ``D = t + sqrt(2t) cos psi``.
    """

    a: Tuple[float, float, float]
    b: Tuple[float, float, float]
    sA: int
    sB: int

    def __post_init__(self):
        for v in (self.a, self.b):
            if abs(float(np.linalg.norm(v)) - 1.0) > 1e-12:
                raise ValueError("settings must be unit vectors")
        if self.sA not in (-1, 1) or self.sB not in (-1, 1):
            raise ValueError("sector signs must be +-1")
        if np.linalg.norm(np.cross(self.a, self.b)) < 1e-9:
            raise ValueError("collinear settings are excluded from the regular locus")

    @property
    def u(self) -> np.ndarray:
        return self.sA * np.asarray(self.a, dtype=float)

    @property
    def w(self) -> np.ndarray:
        return -self.sB * np.asarray(self.b, dtype=float)

    @property
    def q(self) -> np.ndarray:
        return np.cross(self.u, self.w)

    @property
    def s(self) -> np.ndarray:
        return self.u + self.w

    @property
    def t(self) -> float:
        return float(1.0 + self.u @ self.w)

    def invariants(self, x: np.ndarray) -> Tuple[float, float, float]:
        """``(N, D, theta)`` at a point of the sphere."""
        x = np.asarray(x, dtype=float)
        N = float(x @ self.q)
        D = float(self.t + x @ self.s)
        return N, D, math.atan2(N, D)

    def circle_point(self, psi: float) -> np.ndarray:
        """Arclength parametrisation of the closure circle ``N = 0``.

        ``psi`` is measured from ``s/|s|``; ``D(psi) = t + sqrt(2t) cos psi``
        follows because ``s`` is orthogonal to ``q``.
        """
        shat = _unit(self.s)
        perp = np.cross(_unit(self.q), shat)
        return math.cos(psi) * shat + math.sin(psi) * perp

    def D_of_psi(self, psi):
        psi = np.asarray(psi, dtype=float)
        return self.t + math.sqrt(2.0 * self.t) * np.cos(psi)


def _tangent_frame(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    ref = np.array([0.0, 0.0, 1.0]) if abs(x[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = _unit(np.cross(x, ref))
    return e1, np.cross(x, e1)


def _exp_map(x: np.ndarray, h1: float, h2: float,
             e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    v = h1 * e1 + h2 * e2
    r = float(np.linalg.norm(v))
    if r == 0.0:
        return x.copy()
    return math.cos(r) * x + math.sin(r) * (v / r)


def _wrap_pi(d: float) -> float:
    """Fold a phase difference into ``(-pi/2, pi/2]``.

    ``theta = atan2(N, D)`` jumps by ``2 pi`` across the branch cut whenever
    ``D < 0``, and the closure locus ``N = 0`` lies exactly on that cut. Only
    the smooth local branch is physical: closure is a condition on
    ``dist(theta, pi Z)``, so differences are defined modulo ``pi``. Without
    this fold the measured gradient diverges as the step shrinks on every
    ``D < 0`` point, which is half the closure circle.
    """
    return ((d + math.pi / 2.0) % math.pi) - math.pi / 2.0


def _directional_derivative(f, x, e, step, wrap=False):
    """Central difference of ``f`` along a tangent direction, via ``exp_x``."""
    other = np.cross(x, e)
    d = (f(_exp_map(x, step, 0.0, e, other))
         - f(_exp_map(x, -step, 0.0, e, other)))
    return (_wrap_pi(d) if wrap else d) / (2.0 * step)


def phase_gradient(tri: Triangle, x: np.ndarray, step: float) -> np.ndarray:
    """Central-difference gradient of ``theta`` in an orthonormal tangent frame.

    Computed from the map itself, not from the closed form it is compared to.
    """
    x = np.asarray(x, dtype=float)
    return np.array([_directional_derivative(lambda y: tri.invariants(y)[2], x, e,
                                             step, wrap=True)
                     for e in _tangent_frame(x)])


def numerator_gradient(tri: Triangle, x: np.ndarray, step: float) -> np.ndarray:
    """Central-difference gradient of ``N``, the round-8 conditioning control."""
    x = np.asarray(x, dtype=float)
    return np.array([_directional_derivative(lambda y: tri.invariants(y)[0], x, e, step)
                     for e in _tangent_frame(x)])


def joint_normal_gram(tri1: Triangle, tri2: Triangle, x1: np.ndarray,
                      x2: np.ndarray, step: float) -> np.ndarray:
    """Gram matrix of ``dF dF^T`` for ``F = (theta1, theta2)`` on ``S2 x S2``.

    The two differentials occupy orthogonal tangent blocks because ``theta_i``
    depends only on ``x_i`` — that is what "independently prepared" means, so
    the vanishing off-diagonal is a structural consequence of the map, not a
    measurement. What is measured is ``sqrt(det)`` against the closed form.
    """
    g1 = phase_gradient(tri1, x1, step)
    g2 = phase_gradient(tri2, x2, step)
    rows = np.array([np.concatenate([g1, np.zeros(2)]),
                     np.concatenate([np.zeros(2), g2])])
    return rows @ rows.T


def joint_closure_jacobian(tri1: Triangle, tri2: Triangle,
                           x1: np.ndarray, x2: np.ndarray) -> float:
    """Closed form ``|q1||q2| / |D1 D2|`` of the joint coarea Jacobian."""
    D1 = tri1.invariants(x1)[1]
    D2 = tri2.invariants(x2)[1]
    return (float(np.linalg.norm(tri1.q)) * float(np.linalg.norm(tri2.q))
            / abs(D1 * D2))


def joint_coarea_density(tri1: Triangle, tri2: Triangle,
                         x1: np.ndarray, x2: np.ndarray) -> float:
    """Reference joint density ``|D1 D2| / (|q1||q2|)`` per product arclength.

    A function of the **absolute** product alone at fixed settings. The freeze
    records this as a design consequence, not a discovery.
    """
    return 1.0 / joint_closure_jacobian(tri1, tri2, x1, x2)


def W1_closed(t: float) -> float:
    """``int_0^{2pi} |t + sqrt(2t) cos psi| dpsi`` in closed form.

    For ``t >= 2`` the integrand has no sign change and the value is ``2 pi t``.
    """
    if t <= 0:
        raise ValueError("t must be positive")
    if t >= 2.0:
        return 2.0 * math.pi * t
    k = math.sqrt(2.0 / t)
    psi0 = math.acos(-1.0 / k)
    return t * (4.0 * psi0 + 4.0 * k * math.sin(psi0) - 2.0 * math.pi)


def W1_quadrature(t: float) -> float:
    """Independent quadrature for :func:`W1_closed`, split at the punctures."""
    root = math.sqrt(2.0 * t)
    f = lambda ps: abs(t + root * math.cos(ps))
    if t >= 2.0:
        return quad(f, 0.0, 2.0 * math.pi, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
    psi0 = math.acos(-t / root)
    cuts = [0.0, psi0, 2.0 * math.pi - psi0, 2.0 * math.pi]
    return sum(quad(f, lo, hi, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
               for lo, hi in zip(cuts[:-1], cuts[1:]))


def W1_grid(t: float, n: int) -> float:
    """Uniform arclength midpoint rule with ``n`` points, for grid refinement."""
    psi = (np.arange(n) + 0.5) * 2.0 * math.pi / n
    return float(np.sum(np.abs(t + math.sqrt(2.0 * t) * np.cos(psi)))) * 2.0 * math.pi / n


def sector_triangles(a: Sequence[float], b: Sequence[float]) -> List[Triangle]:
    """The four sectors of one prepared pair, in fixed order."""
    return [Triangle(tuple(a), tuple(b), sA, sB) for sA, sB in SECTOR_SIGNS]


def reference_sector_masses(a1, b1, a2, b2) -> Dict[str, object]:
    """Unnormalised and normalised masses over all 16 joint sectors.

    The coarea density is ``|D1 D2| / (|q1||q2|)`` per product arclength, so a
    joint sector mass is ``W1(t1) W1(t2) / (|q1||q2|)`` times the ``1/16``
    prior. Normalisation is taken once over all 16 sectors, as the freeze
    requires, not separately inside each sector.
    """
    tri1, tri2 = sector_triangles(a1, b1), sector_triangles(a2, b2)
    unnormalised, labels = [], []
    for i, T1 in enumerate(tri1):
        for j, T2 in enumerate(tri2):
            m = (W1_closed(T1.t) * W1_closed(T2.t)
                 / (float(np.linalg.norm(T1.q)) * float(np.linalg.norm(T2.q))))
            unnormalised.append(m / 16.0)
            labels.append(f"({T1.sA:+d},{T1.sB:+d})x({T2.sA:+d},{T2.sB:+d})")
    total = float(sum(unnormalised))
    probs = [m / total for m in unnormalised]
    single1 = [W1_closed(T.t) / float(np.linalg.norm(T.q)) for T in tri1]
    single2 = [W1_closed(T.t) / float(np.linalg.norm(T.q)) for T in tri2]
    s1 = [v / sum(single1) for v in single1]
    s2 = [v / sum(single2) for v in single2]
    marginal_error = max(
        abs(probs[4 * i + j] - s1[i] * s2[j]) for i in range(4) for j in range(4))
    return {"labels": labels, "unnormalised": unnormalised, "probabilities": probs,
            "total_unnormalised": total,
            "single_pair_probabilities": [s1, s2],
            "product_marginal_error": marginal_error}


def three_factor_associativity(a1, b1, a2, b2) -> float:
    """Regression control: grouping of three independent factors is immaterial.

    Structurally guaranteed by the product construction; recorded as a
    regression oracle, never as selection evidence.
    """
    tri = [sector_triangles(a1, b1), sector_triangles(a2, b2), sector_triangles(a1, b1)]
    def mass(T):
        return W1_closed(T.t) / float(np.linalg.norm(T.q))
    left = right = 0.0
    worst = 0.0
    for A in tri[0]:
        for B in tri[1]:
            for C in tri[2]:
                left = (mass(A) * mass(B)) * mass(C)
                right = mass(A) * (mass(B) * mass(C))
                worst = max(worst, abs(left - right) / abs(right))
    return worst


def _negative_set_measure(g, lo: float, hi: float, nodes: int = 257
                          ) -> Tuple[float, int]:
    """Measure of ``{g < 0}`` on ``[lo, hi]``, plus the number of components.

    Correction N25. The first implementation assumed that if ``g`` was
    negative at both ends the whole interval was accepted. The accepted set
    can be disconnected -- at ``gamma = 0.1``, sector ``(+,-)``, ``psi = pi``,
    ``epsilon = 0.1`` it is a neighbourhood of ``z = 0`` together with one of
    ``z = 1``, and the excluded middle was being counted.
    """
    zs = np.linspace(lo, hi, nodes)
    vals = np.array([g(float(z)) for z in zs])
    negative = vals < 0.0
    total, components = 0.0, 0
    i = 0
    while i < nodes:
        if not negative[i]:
            i += 1
            continue
        j = i
        while j + 1 < nodes and negative[j + 1]:
            j += 1
        left = zs[i] if i == 0 else brentq(g, zs[i - 1], zs[i],
                                           xtol=1e-14, rtol=8.9e-16)
        right = zs[j] if j == nodes - 1 else brentq(g, zs[j], zs[j + 1],
                                                    xtol=1e-14, rtol=8.9e-16)
        total += right - left
        components += 1
        i = j + 1
    return total, components


def window_slice(tri: Triangle, psi: float, epsilon: float) -> Tuple[float, int]:
    """Measure of the accepted ``z`` set at fixed ``psi``, and its components.

    Coordinates ``x = sqrt(1-z^2) r(psi) + z qhat`` with area element
    ``dpsi dz``, so ``N = z|q|`` and
    ``D = t + sqrt(1-z^2) sqrt(2t) cos psi``. ``dist(theta, pi Z) < epsilon``
    is exactly ``|N| < |D| tan(epsilon)``; the inequality is solved directly
    and the limiting coarea density is never inserted as an integrand.
    ``g`` is even in ``z``, so the measure on ``[0, 1]`` is doubled.
    """
    Aq = float(np.linalg.norm(tri.q))
    root = math.sqrt(2.0 * tri.t) * math.cos(psi)
    tan_e = math.tan(epsilon)
    g = lambda z: Aq * z - abs(tri.t + math.sqrt(max(1.0 - z * z, 0.0)) * root) * tan_e
    measure, components = _negative_set_measure(g, 0.0, 1.0)
    return 2.0 * measure, components


def finite_window_mass(tri: Triangle, epsilon: float) -> Tuple[float, int]:
    """``int dpsi dz`` over the finite window, and the worst component count."""
    worst = 0

    def integrand(ps):
        nonlocal worst
        measure, components = window_slice(tri, ps, epsilon)
        worst = max(worst, components)
        return measure

    if tri.t >= 2.0:
        cuts = [0.0, math.pi, 2.0 * math.pi]
    else:
        psi0 = math.acos(-tri.t / math.sqrt(2.0 * tri.t))
        cuts = [0.0, psi0, math.pi, 2.0 * math.pi - psi0, 2.0 * math.pi]
    total = sum(quad(integrand, lo, hi, limit=200, epsabs=1e-10, epsrel=1e-10)[0]
                for lo, hi in zip(cuts[:-1], cuts[1:]))
    return total, worst


def window_convergence(a1, b1, a2, b2,
                       widths: Sequence[float] = (0.04, 0.02, 0.01, 0.005)
                       ) -> Dict[str, object]:
    """Finite windows against the coarea limit, with asymmetric width pairs.

    The joint window region is a product set and the measure is a product
    measure, so joint sector probabilities factor exactly. That factorisation
    is a regression control. What is measured here is the limit: whether the
    normalised finite-window probabilities approach the coarea prediction.
    """
    tri1, tri2 = sector_triangles(a1, b1), sector_triangles(a2, b2)
    target = reference_sector_masses(a1, b1, a2, b2)["probabilities"]
    rows, factor_error, components = [], 0.0, 0
    for e in widths:
        for scale1, scale2 in ((1, 1), (1, 2), (2, 1)):
            sliced1 = [finite_window_mass(T, scale1 * e) for T in tri1]
            sliced2 = [finite_window_mass(T, scale2 * e) for T in tri2]
            m1 = [v for v, _ in sliced1]
            m2 = [v for v, _ in sliced2]
            components = max([components] + [c for _, c in sliced1 + sliced2])
            joint = [m1[i] * m2[j] for i in range(4) for j in range(4)]
            total = sum(joint)
            probs = [v / total for v in joint]
            discrepancy = max(abs(p - q) for p, q in zip(probs, target))
            p1 = [v / sum(m1) for v in m1]
            p2 = [v / sum(m2) for v in m2]
            factor_error = max(factor_error, max(
                abs(probs[4 * i + j] - p1[i] * p2[j]) for i in range(4) for j in range(4)))
            rows.append({"epsilon": (scale1 * e, scale2 * e),
                         "max_discrepancy_from_coarea": discrepancy})
    smallest = min(widths)
    final = max(r["max_discrepancy_from_coarea"] for r in rows
                if min(r["epsilon"]) <= smallest + 1e-15)
    return {"rows": rows, "final_window_discrepancy": final,
            "joint_factorisation_regression": factor_error,
            # the registered settings and widths keep the accepted z-set
            # connected; a value above 1 flags the regime N25 describes
            "max_accepted_components": components,
            "monotone": all(rows[i]["max_discrepancy_from_coarea"]
                            >= rows[i + 3]["max_discrepancy_from_coarea"] - 1e-12
                            for i in range(len(rows) - 3))}


def puncture_geometry(tri: Triangle) -> Dict[str, float]:
    """The two punctures of the closure circle and the exact excision law.

    ``D`` vanishes on ``Gamma`` exactly at ``x = -u`` and ``x = -w``. There
    ``|dD/dpsi| = sqrt(t(2-t)) = |q|`` identically. That derivative identity
    is exact. It makes the arc ``|D|/|q| < eta`` equal to ``|psi - psi0| < eta``
    only to leading order, so the excised mass is ``2 eta^2`` asymptotically
    and not exactly; see :func:`excision_two_term` for the quartic term
    (correction N23).

    Note ``t = 1 -+ cos gamma`` lies strictly inside ``(0, 2)`` for any
    non-collinear settings, so every regular closure circle has exactly two
    punctures; ``t >= 2`` is unreachable and is handled only defensively.
    """
    if tri.t >= 2.0:
        return {"has_punctures": 0.0, "slope_minus_q": 0.0}
    psi0 = math.acos(-tri.t / math.sqrt(2.0 * tri.t))
    slope = math.sqrt(2.0 * tri.t) * abs(math.sin(psi0))
    qn = float(np.linalg.norm(tri.q))
    at_u = min(float(np.linalg.norm(tri.circle_point(s * psi0) + tri.u))
               for s in (1.0, -1.0))
    at_w = min(float(np.linalg.norm(tri.circle_point(s * psi0) + tri.w))
               for s in (1.0, -1.0))
    return {"has_punctures": 1.0, "psi0": psi0, "slope": slope,
            "slope_minus_q": abs(slope - qn),
            "puncture_is_minus_u": at_u, "puncture_is_minus_w": at_w}


def puncture_arc_bounds(tri: Triangle, psi0: float, eta: float) -> Tuple[float, float]:
    """Solve ``|D|/|q| = eta`` around a puncture, the domain the freeze names.

    Correction N23. The first implementation integrated ``|psi - psi0| < eta``
    instead. Those domains agree only to leading order, because
    ``|dD/dpsi| = |q|`` holds *at* the puncture and not on the whole arc.
    """
    qn = float(np.linalg.norm(tri.q))
    f = lambda ps: abs(float(tri.D_of_psi(ps))) / qn - eta
    span = min(0.5, 0.9 * psi0)
    return (brentq(f, psi0 - span, psi0 - 1e-15, xtol=1e-15, rtol=8.9e-16),
            brentq(f, psi0 + 1e-15, psi0 + span, xtol=1e-15, rtol=8.9e-16))


def excision_two_term(t: float, eta: float) -> float:
    """``M_eta = 2 eta^2 + (1+t)/(2-t) eta^4 + O(eta^6)`` for the two punctures.

    Derived from ``D = -|q| phi + t phi^2/2 + |q| phi^3/6 + O(phi^4)`` at a
    puncture: inverting to ``phi(D/|q|)`` and integrating ``|D|/|q|`` over
    ``|D|/|q| < eta`` gives ``eta^2 + 3(2 alpha^2 + 1/6) eta^4 / 2`` per
    puncture with ``alpha = t / (2|q|)``, and ``3(2 alpha^2 + 1/6)``
    collapses to ``(1+t)/(2-t)``. The leading term is settings-independent;
    the quartic correction is not.
    """
    if t >= 2.0:
        return 0.0
    return 2.0 * eta ** 2 + (1.0 + t) / (2.0 - t) * eta ** 4


def excision_masses(tri: Triangle, etas: Sequence[float] = (0.02, 0.01, 0.005)
                    ) -> Dict[str, object]:
    """Excised coarea mass on ``|D|/|q| < eta`` against the two-term law.

    The mass vanishes as ``2 eta^2``; that leading behaviour is exact and
    settings-independent, but the mass itself is not ``2 eta^2`` exactly.
    The gate is the two-term expansion, whose residual must fall like
    ``eta^6``.
    """
    qn = float(np.linalg.norm(tri.q))
    total = W1_closed(tri.t) / qn
    rows = []
    for eta in etas:
        if tri.t >= 2.0:
            measured, lo, hi = 0.0, float("nan"), float("nan")
        else:
            psi0 = math.acos(-tri.t / math.sqrt(2.0 * tri.t))
            lo, hi = puncture_arc_bounds(tri, psi0, eta)
            # D has an elementary primitive, F(psi) = t psi + sqrt(2t) sin psi,
            # and a single sign on each side of the puncture. Using it instead
            # of adaptive quadrature removes a ~1e-10 noise floor that was
            # swamping the genuine O(eta^6) remainder.
            root = math.sqrt(2.0 * tri.t)
            F = lambda ps: tri.t * ps + root * math.sin(ps)
            measured = 2.0 * (abs(F(psi0) - F(lo)) + abs(F(hi) - F(psi0))) / qn
        leading = 0.0 if tri.t >= 2.0 else 2.0 * eta * eta
        two_term = excision_two_term(tri.t, eta)
        rows.append({"eta": eta, "excised": measured, "arc": [lo, hi],
                     "leading_2eta2": leading, "two_term": two_term,
                     "residual_vs_two_term": abs(measured - two_term),
                     "deviation_from_leading": (abs(measured - leading) / leading
                                                if leading > 0 else 0.0),
                     "fraction": measured / total})
    scaled = [r["residual_vs_two_term"] / r["eta"] ** 6
              for r in rows if r["two_term"] > 0]
    improvement = [r["deviation_from_leading"] * r["leading_2eta2"]
                   / max(r["residual_vs_two_term"], 1e-300)
                   for r in rows if r["two_term"] > 0]
    return {"rows": rows, "total_mass": total,
            "max_relative_residual_vs_two_term": max(
                (r["residual_vs_two_term"] / r["two_term"]
                 for r in rows if r["two_term"] > 0), default=0.0),
            "residual_over_eta6": scaled,
            # the substantive claim: the remainder after the quartic term is
            # genuinely sixth order, so this ratio is constant in eta (its
            # value depends on t and diverges as t -> 2)
            "residual_scales_as_eta6": bool(
                len(scaled) < 2
                or max(scaled) / max(min(scaled), 1e-300) < 2.0),
            "quartic_improvement_factor": min(improvement, default=math.inf)}


def joint_excision_bound(tri1: Triangle, tri2: Triangle,
                         etas: Sequence[float] = (0.02, 0.01, 0.005)
                         ) -> List[Dict[str, float]]:
    """Inclusion-exclusion bound for the two excised sets on the product."""
    e1 = {r["eta"]: r["fraction"] for r in excision_masses(tri1, etas)["rows"]}
    e2 = {r["eta"]: r["fraction"] for r in excision_masses(tri2, etas)["rows"]}
    return [{"eta": eta, "factor1": e1[eta], "factor2": e2[eta],
             "joint_excluded_fraction": e1[eta] + e2[eta] - e1[eta] * e2[eta]}
            for eta in etas]


def _psi_for_D(t: float, D: float) -> float:
    """Arclength angle realising a target ``D`` on the closure circle."""
    c = (D - t) / math.sqrt(2.0 * t)
    if abs(c) > 1.0:
        raise ValueError(f"D={D} is outside the range attainable at t={t}")
    return math.acos(c)


def cubic_weight(D: float) -> float:
    """Round 9's nonnegative witness ``Phi(D) = D^2 (1 - D/5)``, for controls."""
    return float(D) ** 2 * (1.0 - float(D) / 5.0)


def level_set_controls(gamma: float = math.pi / 2) -> Dict[str, object]:
    """Do configurations sharing a candidate statistic share the weight?

    Three exact controls at ``t1 = t2 = 1``, all regular and inside the
    attainable range ``[1-sqrt2, 1+sqrt2]``:

    * **reflection** ``psi_i -> -psi_i`` preserves ``(D1, D2)`` exactly, hence
      also the reference density. This is the pair-statistic level set.
    * **same signed product** ``(1,2)`` against ``(sqrt2, sqrt2)``: the
      reference density agrees, the factorwise cubic does not.
    * **same absolute product, opposite sign** ``(1/4,1)`` against
      ``(-1/4,1)``: again the reference density agrees and the cubic does not,
      because ``Phi(d) - Phi(-d) = -2 d^3 / 5``.

    The reference density is a function of the absolute product, so it is a
    function of the signed product too. The separations below show that this
    does **not** transfer to another weight; it is not a universal sufficiency
    theorem.
    """
    a = (0.0, 0.0, 1.0)
    b = (math.sin(gamma), 0.0, math.cos(gamma))
    tri = Triangle(a, b, 1, 1)
    if abs(tri.t - 1.0) > 1e-12:
        raise ValueError("these controls are stated at t = 1")
    qn = float(np.linalg.norm(tri.q))

    def config(D1, D2, sign1=1.0, sign2=1.0):
        p1, p2 = sign1 * _psi_for_D(1.0, D1), sign2 * _psi_for_D(1.0, D2)
        x1, x2 = tri.circle_point(p1), tri.circle_point(p2)
        d1, d2 = tri.invariants(x1)[1], tri.invariants(x2)[1]
        return {"psi": (p1, p2), "D": (d1, d2),
                "reference_density": abs(d1 * d2) / (qn * qn),
                "cubic_product": cubic_weight(d1) * cubic_weight(d2),
                "realisation_error": max(abs(d1 - D1), abs(d2 - D2))}

    reflect_a, reflect_b = config(0.7, 1.9), config(0.7, 1.9, -1.0, -1.0)
    prod_a, prod_b = config(1.0, 2.0), config(math.sqrt(2.0), math.sqrt(2.0))
    sign_a, sign_b = config(0.25, 1.0), config(-0.25, 1.0)

    def compare(p, q, statistic):
        vp, vq = statistic_value(p, statistic), statistic_value(q, statistic)
        match = (max(abs(i - j) for i, j in zip(vp, vq))
                 if isinstance(vp, tuple) else abs(vp - vq))
        # Distance in the product S2 x S2: the joint configurations differ if
        # either factor does. (A per-factor minimum would read zero whenever
        # one factor is deliberately held fixed, as in the sign control.)
        separation = math.sqrt(sum(
            float(np.linalg.norm(tri.circle_point(p["psi"][i])
                                 - tri.circle_point(q["psi"][i]))) ** 2
            for i in (0, 1)))
        return {"statistic": statistic, "statistic_match": match,
                "reference_density_gap": abs(p["reference_density"] - q["reference_density"]),
                "cubic_gap": abs(p["cubic_product"] - q["cubic_product"]),
                "values": (p["D"], q["D"]),
                # A level set is only informative if it contains distinct
                # histories: the freeze notes that different points on one
                # level set can remain distinct field histories.
                "chord_separation": separation,
                "realisation_error": max(p["realisation_error"], q["realisation_error"])}

    return {"reflection": compare(reflect_a, reflect_b, "pair"),
            "signed_product": compare(prod_a, prod_b, "product"),
            "absolute_product": compare(sign_a, sign_b, "absolute_product")}


def statistic_value(cfg: Dict[str, object], statistic: str) -> float:
    D1, D2 = cfg["D"]
    if statistic == "pair":
        return (float(D1), float(D2))
    if statistic == "product":
        return float(D1 * D2)
    if statistic == "absolute_product":
        return abs(float(D1 * D2))
    raise ValueError(f"unknown statistic {statistic!r}")


def generic_closure_rule_is_rank_one() -> Dict[str, object]:
    """Q2: what does the repository's own closure checker do to a union?

    ``history/closure.py`` sums every event and worldline phase into one
    ``total_phase`` and accepts the history when that single number is within
    the branch tolerance of ``pi Z``. Applied to a history that is the union
    of two disconnected preparations, that is the rank-one condition
    ``theta1 + theta2 in pi Z``, which the freeze explicitly forbids as a
    substitute for the two independent conditions: it accepts cancellation
    between two subsystems that do not close separately.

    Demonstrated, not asserted: two sub-loops with phases ``+delta`` and
    ``-delta`` are each far from closure, yet the union is accepted.
    """
    from geometrodynamics.history.closure import (
        Event, EventType, History, Worldline)

    # pi/2 is the largest attainable distance to pi Z, and at the module's
    # default sigma = 0.6 even that scores weight 0.0324 > 0.01, so no phase
    # mismatch whatsoever can fail the gate. Both facts are reported.
    delta = math.pi / 2.0
    sigma = 0.3
    hist = History()
    for i, phase in enumerate((delta, -delta)):
        base = 10 * i
        for k in range(2):
            hist.add_event(Event(event_id=base + k, event_type=EventType.CREATION,
                                 p4=np.array([1.0, 0.0, 0.0, 0.0]), t=0.0,
                                 orientation=+1 if k == 0 else -1))
        hist.add_worldline(Worldline(start_event_id=base, end_event_id=base + 1,
                                     phase_accumulated=phase))
    union = hist.check_closure(sigma=sigma)

    singles = []
    for phase in (delta, -delta):
        h = History()
        for k in range(2):
            h.add_event(Event(event_id=k, event_type=EventType.CREATION,
                              p4=np.array([1.0, 0.0, 0.0, 0.0]), t=0.0,
                              orientation=+1 if k == 0 else -1))
        h.add_worldline(Worldline(start_event_id=0, end_event_id=1,
                                  phase_accumulated=phase))
        singles.append(h.check_closure(sigma=sigma))

    default_worst = History()
    for k in range(2):
        default_worst.add_event(Event(event_id=k, event_type=EventType.CREATION,
                                      p4=np.array([1.0, 0.0, 0.0, 0.0]), t=0.0,
                                      orientation=+1 if k == 0 else -1))
    default_worst.add_worldline(Worldline(start_event_id=0, end_event_id=1,
                                          phase_accumulated=delta))
    worst = default_worst.check_closure()

    return {"subsystem_phase": delta, "sigma": sigma,
            "union_total_phase": hist.total_phase(),
            "union_is_closed": bool(union.is_closed),
            "union_mismatch": union.phase_mismatch,
            "subsystem_is_closed": [bool(r.is_closed) for r in singles],
            "subsystem_mismatch": [r.phase_mismatch for r in singles],
            "default_sigma_worst_case_weight": worst.weight,
            "default_sigma_cannot_reject": bool(worst.is_closed),
            "rank_one_cancellation_demonstrated": bool(
                union.is_closed and not any(r.is_closed for r in singles))}


def _quaternion_multiply(g: np.ndarray, h: np.ndarray) -> np.ndarray:
    return np.concatenate(([g[0] * h[0] - g[1:] @ h[1:]],
                           g[0] * h[1:] + h[0] * g[1:] + np.cross(g[1:], h[1:])))


def based_loop_composition_scope() -> Dict[str, object]:
    """Q2: what the inherited composition theorem gives on the closure locus.

    Correction N24. The first implementation sampled arbitrary holonomy angles
    and reported their generic ``SU(2)`` non-commutativity as the obstruction.
    That is an **off-closure** statement and does not apply to conditioned
    triangles. On ``Gamma_1 x Gamma_2`` we have ``theta_i in pi Z``, so each
    reduced holonomy ``cos theta + sin theta x`` is ``+-1`` — central — and the
    two commute exactly. Inside finite windows of half-width ``epsilon_i`` the
    commutator norm is bounded by ``2 sin(epsilon_1) sin(epsilon_2)``, which
    vanishes with the windows.

    The conclusion is therefore **stronger and different**: composition is
    available on the closure locus, and what it delivers is precisely
    ``theta_1 + theta_2 in pi Z`` — the rank-one condition the freeze forbids
    as a substitute for the two independent conditions. The repository's one
    composition rule is the wrong one, not an unavailable one.

    The generic sampling is retained only as a labelled off-closure control.
    It also assumes a common quaternion frame for the two triangles; the
    freeze forbids treating that as physical without deriving the transport
    identification, and none is derived here.
    """
    rng = np.random.default_rng(SEED)

    def su2(axis, angle):
        return np.concatenate(([math.cos(angle)], math.sin(angle) * axis))

    # (a) same base point: the inherited additivity theorem, reproduced.
    same_base = 0.0
    for _ in range(400):
        x = _unit(rng.normal(size=3))
        t1, t2 = rng.uniform(-2.0, 2.0, size=2)
        same_base = max(same_base, abs(
            _quaternion_multiply(su2(x, t1), su2(x, t2))[0] - math.cos(t1 + t2)))

    # (b) ON the closure locus: theta in pi Z, so the holonomies are central.
    on_closure_commutator = on_closure_central = 0.0
    for _ in range(400):
        x1, x2 = _unit(rng.normal(size=3)), _unit(rng.normal(size=3))
        k1, k2 = rng.integers(-3, 4, size=2)
        G1, G2 = su2(x1, math.pi * k1), su2(x2, math.pi * k2)
        on_closure_commutator = max(on_closure_commutator, float(np.linalg.norm(
            _quaternion_multiply(G1, G2) - _quaternion_multiply(G2, G1))))
        for G in (G1, G2):
            on_closure_central = max(on_closure_central,
                                     float(np.linalg.norm(G[1:])),
                                     abs(abs(G[0]) - 1.0))

    # (c) finite windows: the commutator obeys 2 sin(e1) sin(e2).
    window_rows = []
    for e1, e2 in ((0.04, 0.04), (0.04, 0.08), (0.08, 0.04), (0.01, 0.01)):
        worst = 0.0
        for _ in range(2000):
            d1, d2 = rng.uniform(-e1, e1), rng.uniform(-e2, e2)
            k1, k2 = rng.integers(-2, 3, size=2)
            G1 = su2(_unit(rng.normal(size=3)), math.pi * k1 + d1)
            G2 = su2(_unit(rng.normal(size=3)), math.pi * k2 + d2)
            worst = max(worst, float(np.linalg.norm(
                _quaternion_multiply(G1, G2) - _quaternion_multiply(G2, G1))))
        bound = 2.0 * math.sin(e1) * math.sin(e2)
        window_rows.append({"epsilon": (e1, e2), "max_commutator": worst,
                            "bound_2_sin_sin": bound, "respects_bound": worst <= bound})

    # (d) OFF-CLOSURE CONTROL ONLY. Generic angles, and a common quaternion
    #     frame that is assumed rather than derived.
    off_closure = 0.0
    for _ in range(400):
        t1, t2 = rng.uniform(-2.0, 2.0, size=2)
        off_closure = max(off_closure, float(np.linalg.norm(
            _quaternion_multiply(su2(_unit(rng.normal(size=3)), t1),
                                 su2(_unit(rng.normal(size=3)), t2))
            - _quaternion_multiply(su2(_unit(rng.normal(size=3)), t2),
                                   su2(_unit(rng.normal(size=3)), t1)))))

    return {"same_base_additivity_residual": same_base,
            "on_closure_commutator": on_closure_commutator,
            "on_closure_holonomy_is_central": on_closure_central,
            "window_rows": window_rows,
            "windows_respect_bound": all(r["respects_bound"] for r in window_rows),
            "off_closure_control_commutator": off_closure,
            "off_closure_control_is_not_evidence": True,
            "common_frame_transport_derived": False,
            "composition_on_closure_is_rank_one": True}


def repository_joint_rule_audit() -> Dict[str, object]:
    """Q2: does any inherited module supply a joint weight for two preparations?

    Each entry records the module's actual variables and whether it applies to
    two *disconnected* preparations or only after a connection is added. This
    is a search over the inspected machinery; absence here is a repository gap,
    not a theorem that no BAM completion can supply such a rule.
    """
    rank_one = generic_closure_rule_is_rank_one()
    based = based_loop_composition_scope()
    entries = [
        {"module": "geometrodynamics/history/closure.py",
         "variables": "event phases, worldline phases, orientation-weighted charges",
         "rule": "single summed total_phase within tolerance of pi Z; Gaussian weight",
         "applies_to_disconnected_pairs": True,
         "supplies_joint_weight_rule": False,
         "reason": "the summed rule is rank one on a union and accepts "
                   "cancellation between subsystems that do not close; the "
                   "freeze forbids it as a substitute for the two conditions"},
        {"module": "geometrodynamics/bulk/history_action.py",
         "variables": "SU(2) closure holonomy G of one based loop, theta, S_H",
         "rule": "theta additive for loops based at a common x",
         "applies_to_disconnected_pairs": True,
         "supplies_joint_weight_rule": False,
         "reason": "on the closure locus theta_i lie in pi Z, so the reduced "
                   "holonomies are central and compose exactly -- but what "
                   "that composition delivers is theta_1 + theta_2 in pi Z, "
                   "the rank-one condition the freeze forbids as a substitute "
                   "for the two independent conditions"},
        {"module": "geometrodynamics/bulk/closure_current.py",
         "variables": "one triangle's D, |u x v|, Pin branch label",
         "rule": "positive or holonomy-weighted coarea on a single closure set",
         "applies_to_disconnected_pairs": False,
         "supplies_joint_weight_rule": False,
         "reason": "single-pair measures; no two-preparation rule is defined"},
        {"module": "geometrodynamics/transaction/network.py",
         "variables": "mouth clocks, throat transfer t_AB, loop eigenvalue",
         "rule": "self-consistent field on a throat-connected loop",
         "applies_to_disconnected_pairs": False,
         "supplies_joint_weight_rule": False,
         "reason": "presupposes a throat connecting the mouths; two "
                   "independent preparations have no such channel"},
        {"module": "geometrodynamics/transaction/derived_network.py",
         "variables": "Lambda_l(omega, Delta), G_eff, eta_topo",
         "rule": "derived-geometry loop eigenvalue and effective Green function",
         "applies_to_disconnected_pairs": False,
         "supplies_joint_weight_rule": False,
         "reason": "wires a derived throat into the same connected loop"},
    ]
    return {"entries": entries,
            "rank_one_demonstration": rank_one,
            "based_loop_scope": based,
            "any_module_supplies_joint_weight_rule": any(
                e["supplies_joint_weight_rule"] for e in entries),
            "search_scope": "the five modules named in the freeze's Q2, at the "
                            "pinned baseline; not an exhaustive search of BAM"}


def _rotation(axis: str, angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)
    raise ValueError(axis)


ROTATION_PAIR2 = _rotation("x", 0.61) @ _rotation("z", 0.37)


def gram_and_metric_checks(tri1: Triangle, tri2: Triangle,
                           samples: int = 24,
                           steps: Sequence[float] = (1e-4, 5e-5, 2.5e-5),
                           regular_floor: float = 0.05) -> Dict[str, object]:
    """Measured joint Gram matrix against ``|q1||q2| / |D1 D2|``.

    Only regular points with ``|D_i| / |q_i| >= regular_floor`` are used, as
    the freeze requires. The step schedule is fixed in advance; the finest step
    must improve on the coarsest unless the coarsest is already below 1e-10.
    """
    rng = np.random.default_rng(SEED)
    per_step = {step: 0.0 for step in steps}
    offdiag = 0.0
    used = attempts = 0
    while used < samples:
        attempts += 1
        if attempts > 200 * samples:
            raise RuntimeError("could not find enough regular joint points")
        p1, p2 = rng.uniform(0.0, 2.0 * math.pi, size=2)
        x1, x2 = tri1.circle_point(p1), tri2.circle_point(p2)
        D1, D2 = tri1.invariants(x1)[1], tri2.invariants(x2)[1]
        if (abs(D1) / float(np.linalg.norm(tri1.q)) < regular_floor
                or abs(D2) / float(np.linalg.norm(tri2.q)) < regular_floor):
            continue
        used += 1
        closed = joint_closure_jacobian(tri1, tri2, x1, x2)
        for step in steps:
            G = joint_normal_gram(tri1, tri2, x1, x2, step)
            offdiag = max(offdiag, abs(G[0, 1]), abs(G[1, 0]))
            per_step[step] = max(per_step[step],
                                 abs(math.sqrt(abs(np.linalg.det(G))) - closed) / closed)
    coarse, fine = per_step[steps[0]], per_step[steps[-1]]
    return {"samples": used, "per_step_relative_error": dict(per_step),
            "max_offdiagonal_gram": offdiag,
            "finest_relative_error": fine,
            "improves": bool(fine <= coarse or coarse < 1e-10)}


def analytic_gradient(tri: Triangle, x: np.ndarray) -> np.ndarray:
    """Independent analytic route to ``grad theta``, valid off the closure set.

    ``grad theta = (D grad N - N grad D) / (N^2 + D^2)`` with the spherical
    gradient ``grad_{S2}(x.v) = v - (x.v) x``. This never uses ``N = 0``, so
    comparing it to ``|q| / |D|`` on the closure circle is an independent
    check of the closed form rather than a restatement of it.
    """
    x = np.asarray(x, dtype=float)
    N, D, _ = tri.invariants(x)
    gN = tri.q - (x @ tri.q) * x
    gD = tri.s - (x @ tri.s) * x
    return (D * gN - N * gD) / (N * N + D * D)


def analytic_gram_route(tri1: Triangle, tri2: Triangle, samples: int = 24,
                        regular_floor: float = 0.05) -> Dict[str, float]:
    """Analytic Gram against the closed form, to the freeze's 1e-10 target."""
    rng = np.random.default_rng(SEED + 2)
    worst_single = worst_joint = 0.0
    used = attempts = 0
    while used < samples:
        attempts += 1
        if attempts > 200 * samples:
            raise RuntimeError("could not find enough regular joint points")
        p1, p2 = rng.uniform(0.0, 2.0 * math.pi, size=2)
        x1, x2 = tri1.circle_point(p1), tri2.circle_point(p2)
        D1, D2 = tri1.invariants(x1)[1], tri2.invariants(x2)[1]
        n1, n2 = float(np.linalg.norm(tri1.q)), float(np.linalg.norm(tri2.q))
        if abs(D1) / n1 < regular_floor or abs(D2) / n2 < regular_floor:
            continue
        used += 1
        for tri, x, D, n in ((tri1, x1, D1, n1), (tri2, x2, D2, n2)):
            worst_single = max(worst_single, abs(
                float(np.linalg.norm(analytic_gradient(tri, x))) - n / abs(D)) / (n / abs(D)))
        joint = (float(np.linalg.norm(analytic_gradient(tri1, x1)))
                 * float(np.linalg.norm(analytic_gradient(tri2, x2))))
        closed = joint_closure_jacobian(tri1, tri2, x1, x2)
        worst_joint = max(worst_joint, abs(joint - closed) / closed)
    return {"single_factor_relative_error": worst_single,
            "joint_jacobian_relative_error": worst_joint, "samples": used}


def sector_probability_grids(a1, b1, a2, b2,
                             sizes: Sequence[int] = (512, 1024, 2048)
                             ) -> Dict[str, object]:
    """Normalised joint sector probabilities from uniform grids per factor."""
    target = reference_sector_masses(a1, b1, a2, b2)["probabilities"]
    tri1, tri2 = sector_triangles(a1, b1), sector_triangles(a2, b2)
    rows = []
    for n in sizes:
        m1 = [W1_grid(T.t, n) / float(np.linalg.norm(T.q)) for T in tri1]
        m2 = [W1_grid(T.t, n) / float(np.linalg.norm(T.q)) for T in tri2]
        joint = [m1[i] * m2[j] for i in range(4) for j in range(4)]
        total = sum(joint)
        probs = [v / total for v in joint]
        rows.append({"n": n, "max_error": max(abs(p - q) for p, q in zip(probs, target))})
    return {"rows": rows, "final_grid_error": rows[-1]["max_error"]}


def covariance_and_exchange(a1, b1, a2, b2, samples: int = 24) -> Dict[str, float]:
    """A common SO(3) frame change and a copy exchange leave the density fixed."""
    rng = np.random.default_rng(SEED + 1)
    R = ROTATION_PAIR2
    tri1, tri2 = Triangle(tuple(a1), tuple(b1), 1, 1), Triangle(tuple(a2), tuple(b2), 1, -1)
    rot1 = Triangle(tuple(R @ np.asarray(a1)), tuple(R @ np.asarray(b1)), 1, 1)
    rot2 = Triangle(tuple(R @ np.asarray(a2)), tuple(R @ np.asarray(b2)), 1, -1)
    cov = exch = 0.0
    for _ in range(samples):
        p1, p2 = rng.uniform(0.0, 2.0 * math.pi, size=2)
        x1, x2 = tri1.circle_point(p1), tri2.circle_point(p2)
        base = joint_coarea_density(tri1, tri2, x1, x2)
        rotated = joint_coarea_density(rot1, rot2, R @ x1, R @ x2)
        cov = max(cov, abs(rotated - base) / base)
        exch = max(exch, abs(joint_coarea_density(tri2, tri1, x2, x1) - base) / base)
    # copy_exchange is symmetric by construction of the density expression and
    # is reported as a structural regression, not as independent evidence.
    return {"common_frame_covariance": cov, "copy_exchange": exch}


def grid_refinement(t: float, sizes: Sequence[int] = (512, 1024, 2048)) -> Dict[str, object]:
    """Uniform arclength grids against the split independent quadrature."""
    reference = W1_quadrature(t)
    errors = [abs(W1_grid(t, n) - reference) / reference for n in sizes]
    return {"t": t, "sizes": list(sizes), "relative_errors": errors,
            "closed_form_error": abs(W1_closed(t) - reference) / reference,
            "final_error": errors[-1],
            # |D| has a corner at each puncture, so the unsplit uniform grid
            # converges without being monotone. The freeze gates the final
            # grid, not the sequence; monotonicity is reported, not required.
            "decreasing": all(b <= a for a, b in zip(errors[:-1], errors[1:]))}


def verdict(checks: Dict[str, bool], audit: Dict[str, object],
            controls: Dict[str, object]) -> Dict[str, str]:
    """Three separate fields, per the freeze. Never one conflated label.

    A failed mandatory check yields ``UNRESOLVED`` for the affected claim; the
    existence of a passing baseline must not hide a failed factorisation.
    """
    fields = ("reference_composition", "reduction_of_specified_rules",
              "reduction_scope", "additional_physical_rule",
              "consequence_for_selection", "P3_hypotheses_established")
    if not checks or not all(checks.values()):
        # Correction N26: the failure branch must carry every field the
        # renderer and the archive read, or reporting crashes before the
        # UNRESOLVED verdict is ever written.
        failed = sorted(k for k, v in checks.items() if not v) or ["no checks supplied"]
        return dict.fromkeys(fields, "UNRESOLVED") | {
            "failed_checks": failed,
            "reduction_scope": "UNRESOLVED -- required checks failed: "
                               + "; ".join(failed)}
    reference = "INDEPENDENT_PHASE_PRODUCT_VERIFIED"

    # Field 2 has two parts: the reference rule, and any further rule the
    # repository actually specifies. Product sufficiency for the reference
    # rule needs the factorisation identity AND the level-set controls
    # confirming the density is constant on those level sets.
    factorises = (controls["signed_product"]["reference_density_gap"] < 1e-10
                  and controls["absolute_product"]["reference_density_gap"] < 1e-10
                  and controls["reflection"]["reference_density_gap"] < 1e-10)
    reference_reduction = ("PRODUCT_STATISTIC_SUFFICIENT" if factorises
                           else "UNRESOLVED")
    additional = ("JOINT_WEIGHT_RULE_UNSPECIFIED"
                  if not audit["any_module_supplies_joint_weight_rule"]
                  else "UNRESOLVED")

    # Independence alone selects nothing. The conditional obstruction is NOT
    # invoked: P3's hypotheses are not established here, and the factorwise
    # cubic remains an admissible extension.
    consequence = ("NO_SELECTION_FROM_INDEPENDENCE"
                   if reference_reduction == "PRODUCT_STATISTIC_SUFFICIENT"
                   else "UNRESOLVED")
    return {"reference_composition": reference,
            "reduction_of_specified_rules": reference_reduction,
            "reduction_scope": "the reference rule only, through the ABSOLUTE "
                               "product; the separations in the level-set "
                               "controls show the factorisation does not "
                               "transfer to the cubic weight",
            "additional_physical_rule": additional,
            "consequence_for_selection": consequence,
            "failed_checks": [],
            "P3_hypotheses_established": "NO -- neither the physical "
                                         "completeness of the allowed rule "
                                         "family, nor the composite-scalar "
                                         "identification, nor the domain "
                                         "hypothesis is established here"}
