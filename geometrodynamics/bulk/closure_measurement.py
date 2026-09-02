"""Does sharp phase closure of a closed history induce measurement dependence at the source?

Pre-registered in ``docs/closure_measurement_dependence_prereg.md`` (``1b0144e``),
committed before this file. The model, stated there and repeated here so the
code carries its own assumptions:

* the created pair's frame direction ``x in S^2`` has the Haar prior;
* a closed history ``source -> A -> B -> source`` is a loop of frame
  directions ``x -> u -> v -> x`` with ``u = s_A a``, ``v = s_B b`` (outcome
  signs are boundary data of the history) and geodesic legs;
* its spin-1/2 geometric phase is half the solid angle
  ``Omega = 2 atan2(N, D)``, ``N = x.(u x v)``, ``D = 1 + x.u + u.v + v.x``;
* the repository's closure axiom ``Omega/2 = 0 or pi (mod 2 pi)`` is
  ``N = 0``: ``x`` on the great circle through ``a`` and ``b``;
* Haar conditioned on ``N = 0`` is the coarea measure, density ``|D|/(2|u x v|)``,
  the ``eps -> 0`` limit of a uniform window ``|Omega mod 2 pi| < eps``.

No amplitude, no Gaussian, no width, no projector, no Born rule: the quantum
correlation ``cos(gamma)`` appears only as a comparison target and, in P5, as
the identity that the *signed* density ``D/4`` is the real part of the
Bargmann invariant ``Tr(P_x P_u P_v)`` -- a quasi-probability, not a measure
this module adopts.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "solid_angle", "great_circle", "closure_weights", "closed_form_weights",
    "correlation", "closed_form_correlation", "marginals", "chsh", "chsh_max",
    "window_monte_carlo", "bargmann_identity", "measurement_dependence",
    "two_leg_loop_control", "local_detector_control", "gaussian_width_control",
    "verdict",
]

_STD = (0.0, 0.5 * math.pi, 0.25 * math.pi, -0.25 * math.pi)


def _unit(theta: float) -> np.ndarray:
    return np.array([math.sin(theta), 0.0, math.cos(theta)])


def solid_angle(x: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """``Omega = 2 atan2(x.(u x v), 1 + x.u + u.v + v.x)`` for ``x`` of shape (n, 3)."""
    x = np.atleast_2d(x)
    N = x @ np.cross(u, v)
    D = 1.0 + x @ u + float(u @ v) + x @ v
    return 2.0 * np.arctan2(N, D)


def great_circle(a: np.ndarray, b: np.ndarray, n: int = 100000) -> np.ndarray:
    """Points of the great circle through ``a`` and ``b`` (the closure set ``N = 0``)."""
    e1 = a / np.linalg.norm(a)
    w = b - (b @ e1) * e1
    e2 = w / np.linalg.norm(w)
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    return np.outer(np.cos(t), e1) + np.outer(np.sin(t), e2)


def closure_weights(gamma: float, variant: str = "abs", n: int = 100000) -> Dict[Tuple[int, int], float]:
    """``W(s_A, s_B) = (1/2|u x v|) mean_Gamma w(D)`` with ``w = |D|`` (coarea),
    ``D`` (signed, P5) or ``max(D, 0)`` (strict zero, P6); normalised over the
    four outcome pairs."""
    a, b = _unit(0.0), _unit(gamma)
    x = great_circle(a, b, n)
    out = {}
    for sA in (+1, -1):
        for sB in (+1, -1):
            u, v = sA * a, sB * b
            D = 1.0 + x @ u + float(u @ v) + x @ v
            w = {"abs": np.abs(D), "signed": D, "pos": np.maximum(D, 0.0)}[variant]
            out[(sA, sB)] = float(np.mean(w)) / (2.0 * np.linalg.norm(np.cross(u, v)))
    Z = sum(out.values())
    return {k: val / Z for k, val in out.items()}


def closed_form_weights(gamma: float) -> Dict[Tuple[int, int], float]:
    """The hand-derived closed forms for the coarea variant:
    ``W(+,+) = W(-,-) = 2 + c(pi-gamma)/s``, ``W(+,-) = W(-,+) = 2 + s gamma/c``."""
    c, s = math.cos(0.5 * gamma), math.sin(0.5 * gamma)
    like = 2.0 + c * (math.pi - gamma) / s
    unlike = 2.0 + s * gamma / c
    Z = 2.0 * (like + unlike)
    return {(1, 1): like / Z, (-1, -1): like / Z, (1, -1): unlike / Z, (-1, 1): unlike / Z}


def correlation(gamma: float, variant: str = "abs", n: int = 100000) -> float:
    P = closure_weights(gamma, variant, n)
    return sum(sa * sb * P[(sa, sb)] for (sa, sb) in P)


def closed_form_correlation(gamma) -> np.ndarray:
    """``E = [c^2 (pi-gamma) - s^2 gamma] / [2 sin gamma + c^2 (pi-gamma) + s^2 gamma]``."""
    g = np.asarray(gamma, dtype=float)
    c2, s2 = np.cos(0.5 * g) ** 2, np.sin(0.5 * g) ** 2
    return (c2 * (math.pi - g) - s2 * g) / (2.0 * np.sin(g) + c2 * (math.pi - g) + s2 * g)


def marginals(gamma: float, variant: str = "abs") -> Dict[str, float]:
    P = closure_weights(gamma, variant)
    return {"P(A=+)": P[(1, 1)] + P[(1, -1)], "P(B=+)": P[(1, 1)] + P[(-1, 1)]}


def _E_signed_angle(delta: float, variant: str, cache: Dict[float, float]) -> float:
    g = abs(delta) % (2.0 * math.pi)
    g = min(g, 2.0 * math.pi - g)
    key = round(g, 9)
    if key not in cache:
        cache[key] = correlation(max(g, 1e-9), variant, n=20000) if g < math.pi - 1e-9 else \
            -correlation(1e-9, variant, n=20000)
    return cache[key]


def chsh(settings: Sequence[float] = _STD, variant: str = "abs") -> float:
    a, ap, b, bp = settings
    cache: Dict[float, float] = {}
    E = lambda d: _E_signed_angle(d, variant, cache)
    return abs(E(a - b) + E(a - bp) + E(ap - b) - E(ap - bp))


def chsh_max(variant: str = "abs", n_grid: int = 49) -> Dict[str, object]:
    """Maximum over ``(a', b, b')`` on a grid with ``a = 0`` (rotational
    invariance), using the closed form for the coarea variant."""
    if variant == "abs":
        E = lambda d: float(closed_form_correlation(min(abs(d) % (2 * math.pi),
                                                        2 * math.pi - abs(d) % (2 * math.pi)) or 1e-9))
    else:
        cache: Dict[float, float] = {}
        E = lambda d: _E_signed_angle(d, variant, cache)
    grid = np.linspace(-math.pi, math.pi, n_grid)
    best, arg = 0.0, None
    for ap in grid:
        for b in grid:
            for bp in grid:
                S = abs(E(0 - b) + E(0 - bp) + E(ap - b) - E(ap - bp))
                if S > best:
                    best, arg = S, (0.0, float(ap), float(b), float(bp))
    return {"S_max": best, "settings": arg, "below_tsirelson": bool(best < 2.0 * math.sqrt(2.0) - 1e-9)}


def window_monte_carlo(gamma: float, eps: float, n: int = 2000000, seed: int = 0) -> float:
    """Uniform ``x ~ Haar(S^2)``, keep ``|Omega mod 2 pi| < eps``, count outcome
    pairs. The ``eps -> 0`` limit is the coarea measure: the independent
    construction of P7."""
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    a, b = _unit(0.0), _unit(gamma)
    cnt = {}
    for sA in (+1, -1):
        for sB in (+1, -1):
            Om = solid_angle(x, sA * a, sB * b)
            r = np.abs((Om + math.pi) % (2.0 * math.pi) - math.pi)
            cnt[(sA, sB)] = int(np.sum(r < eps))
    Z = sum(cnt.values())
    return sum(sa * sb * cnt[(sa, sb)] for (sa, sb) in cnt) / Z


def bargmann_identity(n: int = 300, seed: int = 1) -> Dict[str, object]:
    """P5 — ``D/4 = Re Tr(P_x P_u P_v)`` for spin-1/2 projectors
    ``P_n = (1 + n.sigma)/2``. The signed closure density is the real part of
    the Bargmann invariant: a quasi-probability."""
    rng = np.random.default_rng(seed)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    def proj(v):
        return 0.5 * (np.eye(2) + v[0] * sx + v[1] * sy + v[2] * sz)
    worst = 0.0
    for _ in range(n):
        x, u, v = [q / np.linalg.norm(q) for q in rng.standard_normal((3, 3))]
        D = 1.0 + x @ u + u @ v + v @ x
        worst = max(worst, abs(D / 4.0 - np.trace(proj(x) @ proj(u) @ proj(v)).real))
    # the signed variant is cos(gamma): check at a few angles
    dev = max(abs(correlation(g, "signed", n=40000) - math.cos(g)) for g in (0.3, 1.0, 2.0))
    negative_arc = {g: float(np.mean((1.0 + great_circle(_unit(0), _unit(g), 20000) @ (_unit(0) + _unit(g))
                                       + math.cos(g)) < 0)) for g in (0.5, 1.0, 2.0)}
    return {"max_residual": worst, "signed_variant_is_cos": bool(dev < 1e-4),
            "signed_variant_max_deviation_from_cos": dev,
            "negative_weight_fraction_of_circle_like_outcomes": negative_arc,
            "expected_fraction": {g: g / (2 * math.pi) for g in (0.5, 1.0, 2.0)}}


def _unit3(theta: float, phi: float = 0.0) -> np.ndarray:
    return np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])


def source_density_on_circle(a: np.ndarray, b: np.ndarray, x: np.ndarray) -> np.ndarray:
    """The outcome-summed conditioned source density on the closure circle
    ``Gamma(a, b)``: ``rho(x|a,b) ∝ sum_{s_A, s_B} |D_s(x)| / (2|u x v|)``,
    evaluated at the points ``x`` (which must lie on ``Gamma``)."""
    dens = np.zeros(len(x))
    for sA in (+1, -1):
        for sB in (+1, -1):
            u, v = sA * a, sB * b
            dens += np.abs(1.0 + x @ u + float(u @ v) + x @ v) / (2.0 * np.linalg.norm(np.cross(u, v)))
    return dens / np.mean(dens)


def measurement_dependence(a: Tuple[float, float], b: Tuple[float, float],
                           a2: Tuple[float, float], b2: Tuple[float, float],
                           n: int = 100000) -> Dict[str, object]:
    """P1 — total-variation distance between the conditioned source measures
    ``rho(x|a,b)`` and ``rho(x|a2,b2)``, settings given as ``(theta, phi)``.

    *Correction note (post-implementation).* The pre-registration said the
    supports are distinct great circles meeting in two points, so ``TV = 1``.
    That holds for **non-coplanar** settings. In the standard Bell
    configuration all analyzers lie in one plane, every pair shares the same
    great circle, and the dependence is in the **density** on it, not the
    support. Both cases are computed here.
    """
    va, vb, va2, vb2 = _unit3(*a), _unit3(*b), _unit3(*a2), _unit3(*b2)
    n1, n2 = np.cross(va, vb), np.cross(va2, vb2)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    coplanar = bool(np.allclose(np.abs(n1 @ n2), 1.0, atol=1e-12))
    if not coplanar:
        return {"coplanar": False, "total_variation": 1.0, "measurement_independence_holds": False,
                "why": "distinct closure circles; supports meet in two points"}
    x = great_circle(va, vb, n)
    d1, d2 = source_density_on_circle(va, vb, x), source_density_on_circle(va2, vb2, x)
    tv = 0.5 * float(np.mean(np.abs(d1 - d2)))
    return {"coplanar": True, "total_variation": tv,
            "measurement_independence_holds": bool(tv < 1e-9),
            "max_density_ratio": float(np.max(d1 / d2)),
            "why": "same closure circle; the conditioned density on it depends on the settings"}


def two_leg_loop_control(n: int = 20000, seed: int = 2) -> Dict[str, object]:
    """Loop-topology control: ``x -> u -> x`` has zero solid angle for every
    ``x``, closure is automatic, and the source measure is untouched."""
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    u = _unit(0.7)
    Om = solid_angle(x, u, x[0])       # degenerate: v = x itself, per point
    Om_pointwise = np.array([solid_angle(xi[None, :], u, xi)[0] for xi in x[:2000]])
    return {"max_abs_solid_angle": float(np.max(np.abs(Om_pointwise))),
            "closure_automatic": bool(np.max(np.abs(Om_pointwise)) < 1e-12),
            "source_measure_conditioned": False}


def local_detector_control(n: int = 400000, seed: int = 3) -> Dict[str, object]:
    """Local-detector control: the previous round's C5 detectors with
    independent Haar ``y_A, y_B`` on the same Haar source ``x`` (partner ``-x``).
    A local model: ``S <= 2`` by Bell's theorem."""
    rng = np.random.default_rng(seed)
    x, yA, yB = [q / np.linalg.norm(q, axis=1, keepdims=True) for q in rng.standard_normal((3, n, 3))]

    def E(a, b):
        A = np.sign(x @ _unit(a) + yA @ _unit(a))
        B = np.sign(-x @ _unit(b) + yB @ _unit(b))
        return float(np.mean(A * B))
    a, ap, b, bp = _STD
    S = abs(E(a, b) + E(a, bp) + E(ap, b) - E(ap, bp))
    return {"E_at_gamma": {g: E(0.0, g) for g in (0.3, 1.0, 2.0)}, "CHSH_standard": S,
            "local_bound_respected": bool(S <= 2.0 + 0.01)}


def gaussian_width_control(gamma: float = 1.0, sigmas: Sequence[float] = (0.6, 0.3, 0.1, 0.03),
                           n: int = 400000, seed: int = 5) -> Dict[str, object]:
    """Width control: replace the sharp condition by the repository's Gaussian
    ``exp(-mismatch^2 / 2 sigma^2)`` in the phase mismatch ``Omega/2 mod pi``.
    ``E`` then depends on ``sigma``; the sharp limit is the coarea value."""
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    a, b = _unit(0.0), _unit(gamma)
    rows = []
    for sig in sigmas:
        w = {}
        for sA in (+1, -1):
            for sB in (+1, -1):
                Om = solid_angle(x, sA * a, sB * b)
                mism = np.abs((0.5 * Om + 0.5 * math.pi) % math.pi - 0.5 * math.pi)
                w[(sA, sB)] = float(np.sum(np.exp(-mism ** 2 / (2 * sig ** 2))))
        Z = sum(w.values())
        rows.append({"sigma": sig, "E": sum(sa * sb * w[(sa, sb)] for (sa, sb) in w) / Z})
    return {"rows": rows, "sharp_limit": float(closed_form_correlation(gamma)),
            "depends_on_sigma": bool(abs(rows[0]["E"] - rows[-1]["E"]) > 1e-2),
            "repository_sigma": 0.6}


def verdict(E_dev_from_cos: float, S_max: float, marginal_dev: float, tv: float) -> str:
    if marginal_dev > 1e-9:
        return "CLOSURE_SIGNALS"
    if tv < 1e-9:
        return "CLOSURE_INDUCES_NO_MEASUREMENT_DEPENDENCE"
    if E_dev_from_cos < 1e-4:
        return "CLOSURE_REPRODUCES_QUANTUM_CORRELATIONS"
    if S_max < 2.0 * math.sqrt(2.0):
        return "CLOSURE_INDUCES_SETTING_DEPENDENT_SOURCE_MEASURE_NO_SIGNALLING_NOT_BORN"
    return "UNDETERMINED"
