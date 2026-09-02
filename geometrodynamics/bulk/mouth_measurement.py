"""Can a deterministic classical detector on the mouth spin-frame bundle produce Born frequencies?

Pre-registered in ``docs/classical_born_prereg.md`` (``7ff6e41``), committed
before this file. The classical state is a point ``q`` of the mouth spin-frame
bundle ``S^3 = P_Spin(S^2)`` (``mouth_spin_frame``): a Bloch direction
``x = q^{-1} i q`` fixed by the preparation and a fibre phase ``phi`` left
unresolved, with the tangent frame rotating by ``2 phi``. A detector with
analyzer axis ``a`` is a deterministic map ``D_a(q) in {+1, -1}``. The ensemble
measure is Haar on the fibre, or Haar on ``S^2`` where a candidate says so.

*On the status of fibre Haar (post-review).* The fibre is the Spin(2)
spin-frame fibre. Whether its coordinate is gauge or a physical microscopic
coordinate is a fork the repository has not decided. If gauge, averaging over
it is gauge averaging and cannot produce stochastic outcomes at all. If
physical, Haar is the natural invariant measure but its emergence from
preparation dynamics (phase mixing) is not shown anywhere. "Unresolved" is
epistemic. The ledger therefore carries fibre Haar as *natural invariant
measure; physical preparation derivation open* -- and the negative result
holds even granting it.

Nothing here uses a projector, ``|amplitude|^2`` as a probability, a random
number drawn with quantum weights, a tensor-product state, CHSH or QED. The
Born function ``cos^2(theta/2)`` appears only as the comparison target.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from geometrodynamics.bulk.mouth_spin_frame import spin_frame, qmul

__all__ = [
    "born", "relative_configuration", "arc_measure",
    "linear_threshold", "intensity_detector", "two_harmonic", "tuned_born_basin",
    "induced_probability", "linear_threshold_closed_form",
    "classification_theorem", "linear_family_best_fit", "two_harmonic_natural_weightings",
    "symmetric_basin", "detector_symmetry_check", "detector_mouth_pushforward",
    "archimedes_uniformity", "archimedes_probability", "archimedes_monte_carlo",
    "repository_winding_detector", "measure_control", "reversal_control", "verdict",
    "narrow_verdict",
]

_I = np.array([0, 1.0, 0, 0])


def born(theta) -> np.ndarray:
    """The comparison target ``cos^2(theta/2)``. Never used to build a detector."""
    return np.cos(np.asarray(theta, dtype=float) / 2.0) ** 2


# ── the relative configuration of analyzer and frame ───────────────────────

def relative_configuration(a: np.ndarray, q: np.ndarray) -> Dict[str, float]:
    """``theta`` = angle between ``a`` and the Bloch direction ``x``; ``psi`` =
    azimuth of ``a`` in the tangent frame ``(e_1, e_2)`` at ``x``. Under the
    fibre ``q -> e^{i phi} q`` the frame turns by ``2 phi``, so ``psi`` shifts
    by ``-2 phi`` and ``x`` is fixed (checked in ``mouth_spin_frame``)."""
    x, e1, e2 = spin_frame(q)
    a = np.asarray(a, dtype=float)
    a = a / np.linalg.norm(a)
    theta = math.acos(float(np.clip(a @ x, -1.0, 1.0)))
    psi = math.atan2(float(a @ e2), float(a @ e1))
    return {"theta": theta, "psi": psi}


def arc_measure(detector: Callable[[float, np.ndarray], np.ndarray], theta: float,
                n: int = 20000, weight: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> float:
    """Measure of ``{phi : D(theta, phi) = +1}`` on the fibre circle, under Haar
    (uniform) or an optional non-invariant weight for the measure control."""
    phi = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    plus = detector(theta, phi) > 0
    if weight is None:
        return float(np.mean(plus))
    w = weight(phi)
    return float(np.sum(w * plus) / np.sum(w))


# ── detector families, each an actual coupling, each deterministic ─────────

def linear_threshold(alpha_over_rho: float, psi0: float = 0.0):
    """C1 — sign of a linear functional of the frame: a torque or gradient
    coupling ``a . (alpha x + rho (cos psi0 e_1 + sin psi0 e_2))``. In terms
    of the fibre phase, ``a . e_1 = sin(theta) cos(psi)`` with
    ``psi = psi_a - 2 phi``."""
    r = alpha_over_rho

    def D(theta, phi):
        psi = -2.0 * np.asarray(phi, dtype=float)
        return r * math.cos(theta) + math.sin(theta) * np.cos(psi - psi0)
    return D


def linear_threshold_closed_form(alpha_over_rho: float, theta) -> np.ndarray:
    """``f(theta) = arccos(-(alpha/rho) cot theta)/pi`` clipped to ``[0, 1]``."""
    theta = np.asarray(theta, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = -alpha_over_rho / np.tan(theta)
    out = np.where(c <= -1.0, 1.0, np.where(c >= 1.0, 0.0, np.arccos(np.clip(c, -1.0, 1.0)) / math.pi))
    out = np.where(theta == 0.0, 1.0 if alpha_over_rho > 0 else 0.5, out)
    out = np.where(theta == math.pi, 0.0 if alpha_over_rho > 0 else 0.5, out)
    return out


def intensity_detector():
    """C2 — the classical Malus intensities ``|<a,±|q>|^2 = cos^2, sin^2 (theta/2)``
    are phase-independent; a detector comparing them is ``sign(cos theta)``."""
    def D(theta, phi):
        return np.full_like(np.asarray(phi, dtype=float), math.cos(theta))
    return D


def two_harmonic(A: Callable[[float], float], B: Callable[[float], float],
                 d1: float = 0.0, d2: float = 0.0):
    """C3 — a coupling that sees the spinor (``e^{i phi}``) and the frame
    (``e^{2 i phi}``): ``sign(A(theta) cos(phi + d1) + B(theta) cos(2 phi + d2))``."""
    def D(theta, phi):
        phi = np.asarray(phi, dtype=float)
        return A(theta) * np.cos(phi + d1) + B(theta) * np.cos(2.0 * phi + d2)
    return D


def _circular_distance_to(psi: np.ndarray, centre: float) -> np.ndarray:
    return np.abs((psi - centre + math.pi) % (2.0 * math.pi) - math.pi)


def symmetric_basin(f: Callable[[float], float]):
    """The constructive basin of H1b: ``D = +1 iff |psi - pi/2|_circ < pi f(theta)``.

    Centred at ``psi = pi/2`` so that it is invariant under the complementarity
    constraint ``psi -> pi - psi`` and, its complement being centred at
    ``-pi/2``, satisfies reversal ``D(pi - theta, psi + pi) = -D(theta, psi)``
    whenever ``f(pi - theta) = 1 - f(theta)``. *Correction note:* the first
    version used ``|psi| < pi f``, centred at ``0``, which is invariant under
    ``psi -> -psi`` but not under ``psi -> pi - psi``; the arc-measure
    conclusion was unaffected, the constructive proof was not what the
    pre-registration claimed. ``detector_symmetry_check`` now tests the
    detector-level conditions directly."""
    def D(theta, phi):
        psi = -2.0 * np.asarray(phi, dtype=float)
        return math.pi * f(theta) - _circular_distance_to(psi, 0.5 * math.pi)
    return D


def tuned_born_basin():
    """Basin control — ``symmetric_basin(cos^2(theta/2))``. Reproduces Born by
    construction; classified TUNED and never counted as derived."""
    return symmetric_basin(lambda t: math.cos(t / 2.0) ** 2)


def detector_symmetry_check(f: Callable[[float], float], n_theta: int = 61,
                            n_phi: int = 720) -> Dict[str, object]:
    """Detector-level (not probability-level) checks of the two constraints on
    ``D(theta, psi)`` for the constructive basin of ``f``: complementarity
    ``D(theta, pi - psi) = D(theta, psi)`` and reversal
    ``D(pi - theta, psi + pi) = -D(theta, psi)``."""
    D = symmetric_basin(f)
    th = np.linspace(0.05, math.pi - 0.05, n_theta)
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    psi = -2.0 * phi
    worst_rev = worst_comp = 0
    for t in th:
        d = np.sign(D(t, phi))
        d_comp = np.sign(D(t, -(math.pi - psi) / 2.0))        # psi -> pi - psi
        d_rev = np.sign(D(math.pi - t, phi - 0.5 * math.pi))   # psi -> psi + pi
        worst_comp = max(worst_comp, int(np.sum(d_comp != d)))
        worst_rev = max(worst_rev, int(np.sum(d_rev != -d)))
    return {"complementarity_violations": worst_comp, "reversal_violations": worst_rev,
            "both_hold": worst_comp <= 2 and worst_rev <= 2}   # boundary grid points


def induced_probability(detector, thetas: Sequence[float], n: int = 20000,
                        weight=None) -> np.ndarray:
    return np.array([arc_measure(detector, float(t), n, weight) for t in thetas])


# ── H1: the classification ─────────────────────────────────────────────────

def classification_theorem(n_theta: int = 181) -> Dict[str, object]:
    """H1a/H1b — reversal ``f(pi - theta) = 1 - f(theta)`` is satisfied by Born,
    by the straight line and by the step; every such ``f`` is realised by a
    symmetric basin; so symmetry + fibre Haar do not select Born."""
    th = np.linspace(0.0, math.pi, n_theta)
    candidates = {"born": born(th), "line": 1.0 - th / math.pi,
                  "step": (th < math.pi / 2).astype(float) + 0.5 * (th == math.pi / 2)}
    reversal = {k: float(np.max(np.abs(v[::-1] - (1.0 - v)))) for k, v in candidates.items()}
    # H1b: the basin |psi| < pi f(theta) realises f under fibre Haar
    realised = {}
    for k, v in candidates.items():
        fvals = dict(zip(np.round(th, 12), v))
        D = symmetric_basin(lambda theta, fv=fvals: fv[round(theta, 12)])
        got = induced_probability(D, th[1:-1], n=20000)
        realised[k] = float(np.max(np.abs(got - v[1:-1])))
    return {"reversal_residuals": reversal,
            "all_satisfy_reversal": all(r < 1e-12 for r in reversal.values()),
            "basin_realisation_residuals": realised,
            "every_f_is_realisable": all(r < 2e-4 for r in realised.values()),
            "conclusion": ("symmetry and fibre Haar constrain f only by f(pi-theta) = 1-f(theta); "
                           "Born, the line and the step all pass; the basin shape decides.")}


def linear_family_best_fit(ratios: Sequence[float] = tuple(np.linspace(0.05, 20.0, 400)),
                           n_theta: int = 721) -> Dict[str, object]:
    """C1 — the linear-threshold family never reaches Born: best member and its miss."""
    th = np.linspace(0.0, math.pi, n_theta)
    b = born(th)
    misses = [(float(np.max(np.abs(linear_threshold_closed_form(r, th) - b))), float(r)) for r in ratios]
    best_miss, best_r = min(misses)
    # regression: closed form against the arc measure of the actual sign functional
    check_th = np.array([0.3, 0.8, 1.2, 2.0, 2.7])
    numeric = induced_probability(linear_threshold(best_r), check_th, n=40000)
    closed = linear_threshold_closed_form(best_r, check_th)
    return {"best_alpha_over_rho": best_r, "best_max_miss": best_miss,
            "closed_form_vs_arc_measure": float(np.max(np.abs(numeric - closed))),
            "plateaus": "f = 1 for cot theta > rho/alpha and 0 for cot theta < -rho/alpha",
            "reaches_born": bool(best_miss < 1e-3)}


def two_harmonic_natural_weightings(n_theta: int = 181) -> List[Dict[str, object]]:
    """C3 — the four natural geometric weightings, none near Born."""
    th = np.linspace(0.0, math.pi, n_theta)
    b = born(th)
    rows = []
    for name, A, B, d1, d2 in (
            ("(cos th/2, sin th/2), 0, 0", lambda t: math.cos(t / 2), lambda t: math.sin(t / 2), 0.0, 0.0),
            ("(cos th/2, sin th/2), 0, pi/2", lambda t: math.cos(t / 2), lambda t: math.sin(t / 2), 0.0, math.pi / 2),
            ("(sin th, cos th), 0, 0", math.sin, math.cos, 0.0, 0.0),
            ("(cos th, sin th), pi/2, 0", math.cos, math.sin, math.pi / 2, 0.0)):
        f = induced_probability(two_harmonic(A, B, d1, d2), th, n=4000)
        rows.append({"weighting": name, "max_miss": float(np.max(np.abs(f - b)))})
    return rows


# ── C5: the Archimedes route ───────────────────────────────────────────────

def detector_mouth_pushforward(n: int = 20000, seed: int = 4) -> Dict[str, object]:
    """C5's measure, stated precisely: Haar on ``S^2`` is the base marginal
    (pushforward under the Hopf map) of Haar on ``SU(2) = S^3``, not the same
    object. An identical, unprepared detector mouth ``q_D ~ Haar(S^3)`` has
    ``y = h(q_D) ~ Haar(S^2)`` by isotropy: checked by the Kolmogorov distance
    of ``a . y`` from uniform on ``[-1, 1]``. What this does NOT give is the
    coupling weight ``kappa = 1``; that would need a symmetric polarisation
    coupling between two identical mouths, which is not built here."""
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    y = np.array([spin_frame(qi)[0] for qi in q])
    u = np.sort(y[:, 2])
    ks = float(np.max(np.abs(np.arange(1, n + 1) / n - (u + 1.0) / 2.0)))
    return {"kolmogorov_distance_of_a_dot_y": ks, "base_marginal_is_haar_S2": bool(ks < 1.5e-2),
            "kappa_derived": False,
            "open_route": ("identical unprepared detector mouth [isotropy gives y Haar on S^2] + "
                           "symmetric classical polarisation coupling [would give kappa = 1]")}


def archimedes_uniformity(n: int = 400000, seed: int = 0) -> Dict[str, object]:
    """``a . y`` is uniform on ``[-1, 1]`` for ``y`` Haar on ``S^2`` (hat-box theorem):
    checked by the Kolmogorov distance of the empirical CDF."""
    rng = np.random.default_rng(seed)
    y = rng.standard_normal((n, 3))
    y /= np.linalg.norm(y, axis=1, keepdims=True)
    u = np.sort(y[:, 2])
    emp = np.arange(1, n + 1) / n
    ks = float(np.max(np.abs(emp - (u + 1.0) / 2.0)))
    return {"kolmogorov_distance": ks, "uniform": bool(ks < 5e-3)}


def archimedes_probability(kappa: float, theta) -> np.ndarray:
    """``D = sign(a.(x + kappa y))`` with ``y`` Haar on ``S^2``:
    ``f = clip((1 + cos theta / kappa)/2, 0, 1)``."""
    return np.clip((1.0 + np.cos(np.asarray(theta, dtype=float)) / kappa) / 2.0, 0.0, 1.0)


def archimedes_monte_carlo(kappa: float, theta: float, n: int = 400000, seed: int = 1) -> float:
    rng = np.random.default_rng(seed)
    y = rng.standard_normal((n, 3))
    y /= np.linalg.norm(y, axis=1, keepdims=True)
    a = np.array([0.0, 0.0, 1.0])
    x = np.array([math.sin(theta), 0.0, math.cos(theta)])
    return float(np.mean(a @ x + kappa * (y @ a) > 0))


# ── C4: the repository's own detector ──────────────────────────────────────

def repository_winding_detector() -> Dict[str, object]:
    """PR #238's derived pairing: a fibre-``U(1)`` setting with a winding-diagonal
    (``sigma_z``) detector. Winding-diagonal means phase-independent: the
    induced ``f`` is the step, and the repository's own measurement round
    (``measurement_sector_probe``) supplies outcome statistics for a
    superposed preparation from an explicitly *Born* throat ensemble."""
    th = np.linspace(0.0, math.pi, 181)
    f = induced_probability(intensity_detector(), th, n=2000)
    return {"phase_independent": True, "induced_f_is_step": bool(
        np.all(f[th < math.pi / 2] == 1.0) and np.all(f[th > math.pi / 2] == 0.0)),
            "max_miss_from_born": float(np.max(np.abs(f - born(th)))),
            "repository_source_of_born_statistics":
                "measurement_sector_probe: 'The winding Stern-Gerlach with a Born throat ensemble'"}


# ── controls ───────────────────────────────────────────────────────────────

def measure_control(alpha_over_rho: float = 0.8) -> Dict[str, object]:
    """Replace fibre Haar by the non-invariant weight ``1 + 0.8 cos(2 phi + 1)``:
    the C1 curve changes, so any result that needed such a weight fails
    Theorem 2.

    *Correction note.* A first version used ``1 + cos(phi)``. Every
    frame-coupled detector has period ``pi`` in ``phi`` (the frame turns by
    ``2 phi``), so a period-``2 pi`` tilt integrates to zero against it and
    the control could not fire (change ``9e-5``). The weight now has the
    detector's period. This changes the control, not any criterion.
    """
    th = np.linspace(0.1, math.pi - 0.1, 61)
    D = linear_threshold(alpha_over_rho)
    haar = induced_probability(D, th, n=8000)
    tilted = induced_probability(D, th, n=8000,
                                 weight=lambda p: 1.0 + 0.8 * np.cos(2.0 * p + 1.0))
    return {"max_change": float(np.max(np.abs(haar - tilted))),
            "measure_matters": bool(np.max(np.abs(haar - tilted)) > 1e-3),
            "weight": "1 + 0.8 cos(2 phi + 1)"}


def reversal_control(curves: Dict[str, np.ndarray]) -> Dict[str, float]:
    return {k: float(np.max(np.abs(v[::-1] - (1.0 - v)))) for k, v in curves.items()}


def verdict(linear_best_miss: float, two_harmonic_misses: Sequence[float],
            archimedes_kappa1_miss: float, archimedes_off_misses: Sequence[float]) -> str:
    repo_coupling_reaches_born = linear_best_miss < 1e-3 or min(two_harmonic_misses) < 1e-3
    if repo_coupling_reaches_born:
        return "BORN_RULE_DERIVED_FROM_CLASSICAL_BAM_MEASURE"
    if archimedes_kappa1_miss < 1e-3 and min(archimedes_off_misses) > 1e-2:
        return "BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW"
    return "CLASSICAL_INTENSITY_ONLY_NO_OUTCOME_PROBABILITY"


def narrow_verdict(pre_registered: str) -> str:
    """*Correction note (post-review).* The pre-registered label reads as a
    no-go against every classical BAM detector. It is not: C5 is an analytic
    Born reproduction whose two inputs (``y`` Haar on ``S^2``, ``kappa = 1``)
    are underived rather than impossible. The scope-correct statement is that
    the preparation and detector dynamics currently in the repository do not
    derive Born."""
    if pre_registered == "BORN_REQUIRES_AN_IMPORTED_MEASURE_OR_DETECTOR_LAW":
        return "CURRENT_BAM_PREPARATION_AND_DETECTOR_DYNAMICS_DO_NOT_DERIVE_BORN"
    return pre_registered
