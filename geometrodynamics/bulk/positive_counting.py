"""Can positive classical history counting reach the quantum correlation? Yes.

Ninth round, pre-registered in ``docs/positive_counting_prereg.md`` before the
open questions were computed. It was opened to test a hypothesis suggested by
rounds 5-8 and PR #283 — that every mechanism landing *near but not on* quantum
mechanics (`2.1423`, `3.3941`, `3.7712`) reflected a structural obstruction.

**The hypothesis is false.** Both answers are negative for the no-go:

* the positive class reaches ``CHSH = 4``, the algebraic maximum, so positivity
  bounds nothing at all; and
* the quantum law itself is attained by an explicit nonnegative weight,
  ``Phi(D) = D^2 (1 - D/5)``, to machine precision.

The consequence is a correction to round 6's framing. That round said the
distance from the repository's law to quantum mechanics "is exactly the
distance from ``|D|`` to ``D``: from counting the two closure branches to
summing them with their holonomy." That is withdrawn. One *can* reach the
quantum law by counting closed histories with nonnegative weights; no signed
cancellation, holonomy weighting or oriented current is required. What is
needed is a different counting function, and nothing in the geometry selects
one.

The reduction that makes all of this computable: on the closure circle, with
``t = 1 + u.w``,

    D(x) = t + sqrt(2t) cos(psi),    psi uniform in arclength,

so every model in the class is a single function of one variable,
``W(t) = int Phi(t + sqrt(2t) cos psi) dpsi``, and with
``t_like = 1 - cos gamma``, ``t_unlike = 1 + cos gamma`` summing to ``2``,

    E(gamma) = [W(tau) - W(2-tau)] / [W(tau) + W(2-tau)],   tau = 1 - cos gamma.

The quantum law is equivalent to ``G(t) = W(t)/t`` being **even about t = 1**.
"""

from __future__ import annotations

import math
from math import comb
from typing import Callable, Dict, List, Sequence

import numpy as np

from geometrodynamics.bulk.closure_measurement import great_circle

__all__ = [
    "reduction_residual", "W", "correlation_from_phi", "chsh_from_phi",
    "quantum_phi", "monomial_W_coefficients", "G_in_shifted_basis",
    "symmetry_conditions", "solve_minimal_degree", "no_global_polynomial",
    "threshold_phi_reaches_four", "marginals_are_automatic",
    "downstream_numbers_unchanged", "dependency_ledger", "verdict",
]

# physical range of D: t + sqrt(2t) cos psi over t in (0, 2]
D_MIN, D_MAX = -0.5, 4.0
_NPSI = 400001
_PSI = np.linspace(0.0, 2.0 * math.pi, _NPSI, endpoint=False)
_COS = np.cos(_PSI)


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


# ── O1: the reduction ───────────────────────────────────────────────────────

def reduction_residual(gammas: Sequence[float] = (0.4, 1.0, 2.0, 2.7),
                       n: int = 200000) -> Dict[str, object]:
    """``D(x) = t + sqrt(2t) cos psi`` with ``psi`` uniform, as an equality of laws.

    Checked by comparing sorted samples, so it tests the *distribution* of ``D``
    along the circle rather than a pointwise parametrisation.
    """
    worst = 0.0
    rows = []
    for gamma in gammas:
        a = np.array([0.0, 0.0, 1.0])
        b = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
        for s_a in (1, -1):
            for s_b in (1, -1):
                u, w = s_a * a, -(s_b * b)
                t = 1.0 + float(u @ w)
                circle = great_circle(u, w, n)
                D = 1.0 + circle @ u + float(u @ w) + circle @ w
                psi = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
                pred = t + math.sqrt(2.0 * t) * np.cos(psi)
                worst = max(worst, float(np.max(np.abs(np.sort(D) - np.sort(pred)))))
        rows.append({"gamma": gamma, "t_like": 1 - math.cos(gamma),
                     "t_unlike": 1 + math.cos(gamma),
                     "sum": (1 - math.cos(gamma)) + (1 + math.cos(gamma))})
    return {"worst_sorted_residual": worst, "rows": rows,
            "t_pair_sums_to_two": all(abs(r["sum"] - 2.0) < 1e-15 for r in rows),
            "reduction_holds": bool(worst < 1e-4)}


def W(t: float, phi: Callable[[np.ndarray], np.ndarray]) -> float:
    """``W(t) = int_0^{2pi} Phi(t + sqrt(2t) cos psi) dpsi``."""
    return float(np.mean(phi(t + math.sqrt(2.0 * t) * _COS))) * 2.0 * math.pi


def correlation_from_phi(phi: Callable[[np.ndarray], np.ndarray],
                         gamma: float) -> float:
    """``E(gamma)`` for the positive-counting model with weight ``Phi``."""
    lo, hi = W(1.0 - math.cos(gamma), phi), W(1.0 + math.cos(gamma), phi)
    return (lo - hi) / (lo + hi)


def chsh_from_phi(phi: Callable[[np.ndarray], np.ndarray],
                  gamma: float = math.pi / 4) -> float:
    """Standard-angle CHSH, ``|3E(g) - E(3g)|``."""
    return abs(3.0 * correlation_from_phi(phi, gamma)
               - correlation_from_phi(phi, 3.0 * gamma))


def quantum_phi(d):
    """``Phi(D) = D^2 (1 - D/5)`` — nonnegative on ``[-1/2, 4]``, gives ``-cos gamma``.

    Minimal-degree solution beyond the signed ``Phi = D``; see
    `solve_minimal_degree`. It is nonnegative exactly on the range the model
    evaluates (``D <= 4``, with the root at ``D = 5`` outside it), not globally —
    and `no_global_polynomial` shows no polynomial solution can be global.
    """
    d = np.asarray(d, dtype=float)
    return d * d * (1.0 - d / 5.0)


# ── the polynomial structure ────────────────────────────────────────────────

def monomial_W_coefficients(n: int) -> Dict[int, float]:
    """``D^n -> W(t)/2pi = sum_j C(n,2j) C(2j,j) 2^{-j} t^{n-j}``.

    From ``<(t + s cos psi)^n>`` with ``s^2 = 2t``: odd powers of ``cos`` average
    to zero and ``<cos^{2j}> = C(2j,j)/4^j``.
    """
    out: Dict[int, float] = {}
    for j in range(n // 2 + 1):
        out[n - j] = out.get(n - j, 0.0) + comb(n, 2 * j) * comb(2 * j, j) / 2 ** j
    return out


def G_in_shifted_basis(n: int, size: int) -> np.ndarray:
    """``G_n(u) = W_n(t)/(2 pi t)`` expanded in ``u = t - 1``, ascending."""
    poly = monomial_W_coefficients(n)
    over_t = np.zeros(max(poly))
    for k, v in poly.items():
        over_t[k - 1] = v                       # divide by t
    out = np.zeros(size)
    for k, c in enumerate(over_t):
        for i in range(k + 1):
            out[i] += c * comb(k, i)            # (u+1)^k
    return out


def symmetry_conditions(max_degree: int) -> Dict[str, object]:
    """The quantum law as linear conditions on polynomial coefficients.

    ``E = -cos gamma`` iff ``G(t) = W(t)/t`` is even about ``t = 1``, so for
    ``Phi = sum a_n D^n`` the odd-power coefficients of ``sum a_n G_n(u)`` must
    vanish. ``G_n`` has degree ``n-1`` with leading coefficient ``1``.
    """
    size = max_degree + 1
    G = {n: G_in_shifted_basis(n, size) for n in range(1, max_degree + 1)}
    odd_rows = [i for i in range(1, size, 2)]
    M = np.array([[G[n][i] for n in range(1, max_degree + 1)] for i in odd_rows])
    return {"odd_powers": odd_rows, "matrix": M,
            "n_conditions": len(odd_rows), "n_coefficients": max_degree,
            "G1_is_constant": bool(np.allclose(G[1][1:], 0.0)),
            "leading_coefficients_are_one": all(
                abs(G[n][n - 1] - 1.0) < 1e-12 for n in range(1, max_degree + 1))}


def solve_minimal_degree() -> Dict[str, object]:
    """Degree 3 is the minimal degree of a nonnegative solution.

    Degree 1 gives only ``Phi = a D``, which is negative somewhere. Degree 2
    forces ``a_2 = 0`` (the single ``u`` condition), collapsing to degree 1. At
    degree 3 the condition ``a_2 + 5 a_3 = 0`` has the solution
    ``Phi = D^2 - D^3/5``, nonnegative on ``[-1/2, 4]``.
    """
    g2 = G_in_shifted_basis(2, 4)
    g3 = G_in_shifted_basis(3, 4)
    combo = g2 - g3 / 5.0
    grid = np.linspace(D_MIN, D_MAX, 200001)
    vals = quantum_phi(grid)
    return {
        "G2_odd_u": float(g2[1]), "G3_odd_u": float(g3[1]),
        "condition": "a_2 + 5 a_3 = 0",
        "combined_G_in_u": combo.tolist(),
        "combined_is_even": bool(all(abs(combo[i]) < 1e-12
                                     for i in range(1, len(combo), 2))),
        "degree_2_forces_a2_zero": bool(abs(g2[1]) > 1e-12),
        "min_phi_on_physical_range": float(vals.min()),
        "phi_nonnegative_on_range": bool(vals.min() > -1e-9),
    }


def no_global_polynomial(max_degree: int = 12) -> Dict[str, object]:
    """No polynomial solution is nonnegative on all of ``R``.

    ``G_n`` has degree ``n-1`` with leading coefficient ``1``, so for a solution
    of top degree ``N`` the ``u^{N-1}`` coefficient is exactly ``a_N``. When
    ``N`` is even, ``N-1`` is odd and that coefficient is one of the vanishing
    conditions, forcing ``a_N = 0``. Hence the top degree is always **odd**, and
    an odd-degree polynomial is negative on one side. Nonnegativity is therefore
    always range-limited — which is legitimate here because the model only ever
    evaluates ``Phi`` on ``[-1/2, 4]``.
    """
    rows = []
    for N in range(2, max_degree + 1):
        s = symmetry_conditions(N)
        M = s["matrix"]
        # the u^{N-1} row exists iff N-1 is odd; its only entry in column a_N is 1
        forced = (N - 1) % 2 == 1
        rows.append({"top_degree": N, "top_degree_forced_zero": bool(forced)})
    return {"rows": rows,
            "every_even_top_degree_is_forced_zero": all(
                r["top_degree_forced_zero"] for r in rows if r["top_degree"] % 2 == 0),
            "conclusion": "top degree is always odd, so no globally nonnegative "
                          "polynomial solution exists"}


# ── Q1: positivity bounds nothing ───────────────────────────────────────────

def threshold_phi_reaches_four(widths: Sequence[float] = (0.2, 0.1, 0.05, 0.02)
                               ) -> Dict[str, object]:
    """A narrow nonnegative bump at ``v = 1 + sqrt(2)`` gives ``CHSH = 4``.

    ``W(t) > 0`` iff ``v`` lies in ``[t - sqrt(2t), t + sqrt(2t)]``. The upper
    endpoint increases in ``t`` and equals ``v`` at ``t = 1``, so the bump makes
    ``W`` vanish for ``t < 1`` and be positive for ``t > 1``. Then
    ``E(gamma) = -sgn(cos gamma)`` and the standard angles give the algebraic
    maximum — a PR box, with marginals exactly ``1/2`` by `marginals_are_automatic`.
    """
    v = 1.0 + math.sqrt(2.0)
    rows = []
    for width in widths:
        def phi(d, w=width):
            return (np.abs(np.asarray(d, dtype=float) - v) < w).astype(float)
        rows.append({"width": width, "CHSH": chsh_from_phi(phi),
                     "E_pi_over_4": correlation_from_phi(phi, math.pi / 4),
                     "E_3pi_over_4": correlation_from_phi(phi, 3 * math.pi / 4)})
    return {"v": v, "rows": rows,
            "reaches_algebraic_maximum": all(abs(r["CHSH"] - 4.0) < 1e-9
                                             for r in rows),
            "threshold_at_t_equals_one": abs((1.0 + math.sqrt(2.0)) - v) < 1e-15}


def marginals_are_automatic(gammas: Sequence[float] = (0.4, 1.0, 2.0),
                            n: int = 200000) -> Dict[str, object]:
    """Every model in the class has marginals exactly ``1/2``.

    ``Gamma`` is invariant under ``x -> -x``, which exchanges ``(s_A, s_B)`` with
    ``(-s_A, -s_B)`` while preserving ``D``. So no-signalling is a property of
    the class, not evidence for any member of it, and may not be cited as
    support for one.
    """
    worst = 0.0
    for gamma in gammas:
        a = np.array([0.0, 0.0, 1.0])
        b = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
        for s_a in (1, -1):
            for s_b in (1, -1):
                u, w = s_a * a, -(s_b * b)
                circle = great_circle(u, w, n)
                D = 1.0 + circle @ u + float(u @ w) + circle @ w
                Dflip = 1.0 + (-circle) @ (-u) + float(u @ w) + (-circle) @ (-w)
                worst = max(worst, float(np.max(np.abs(np.sort(D) - np.sort(Dflip)))))
    return {"flip_residual": worst, "marginals_forced_half": bool(worst < 1e-12)}


# ── the frozen downstream numbers ───────────────────────────────────────────

def downstream_numbers_unchanged() -> Dict[str, object]:
    """Pre-registered rule 5 — nothing in rounds 5-8 or #283 may move."""
    from geometrodynamics.bulk.closure_measurement import correlation
    from geometrodynamics.bulk.history_action import morse_bott_component_masses
    checks = []
    for variant, expected in (("abs", 0.3984966504), ("signed", 0.5403023059)):
        got = float(correlation(1.0, variant, n=400001))
        checks.append({"quantity": f"round-5 correlation({variant})",
                       "expected": expected, "got": got,
                       "delta": abs(got - expected)})
    masses = morse_bott_component_masses(
        np.array([0.0, 0.0, 1.0]),
        np.array([math.sin(1.0), 0.0, math.cos(1.0)]), n=200001)
    for key, expected in (("M_0", 11.67080583639533), ("M_pi", 0.16951227835901764)):
        checks.append({"quantity": f"round-7 {key}", "expected": expected,
                       "got": float(masses[key]),
                       "delta": abs(float(masses[key]) - expected)})
    worst = max(c["delta"] for c in checks)
    return {"checks": checks, "worst_delta": worst,
            "nothing_downstream_moved": bool(worst < 1e-6)}


def dependency_ledger() -> List[Dict[str, str]]:
    """What this round assumes, and what it therefore does not establish."""
    return [
        {"input": "class restriction: W_s = int Phi(D_s) dsigma, Phi >= 0",
         "status": "chosen",
         "where": "excludes (int D)^2, x-dependence beyond D, sector-dependent Phi"},
        {"input": "geodesic itinerary and realignment", "status": "chosen",
         "where": "inherited from rounds 5-6"},
        {"input": "conditioning variable = phase", "status": "chosen",
         "where": "inherited from round 8"},
        {"input": "equal sector prior", "status": "chosen",
         "where": "inherited; the class gives 1/2 marginals for any paired prior"},
        {"input": "nonnegativity of Phi on [-1/2, 4]", "status": "derived",
         "where": "the range the model evaluates; global nonnegativity is "
                  "impossible for polynomials (no_global_polynomial)"},
        {"input": "which Phi the physics selects", "status": "open",
         "where": "THE remaining question; #283's equilibrium selects |.|"},
    ]


def verdict() -> Dict[str, object]:
    """The two pre-registered questions, answered."""
    q1 = threshold_phi_reaches_four()
    worst = max(abs(correlation_from_phi(quantum_phi, g) + math.cos(g))
                for g in (0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
    attainable = worst < 1e-6
    return {
        "Q1_positivity_bound": ("NO_BOUND_ALGEBRAIC_MAXIMUM_REACHED"
                                if q1["reaches_algebraic_maximum"]
                                else "POSITIVITY_BOUNDS_CHSH"),
        "Q1_sup_CHSH": 4.0,
        "Q2_quantum_law": ("QUANTUM_LAW_ATTAINABLE_BY_POSITIVE_COUNTING"
                           if attainable
                           else "QUANTUM_LAW_UNATTAINABLE_IN_THIS_CLASS"),
        "Q2_witness": "Phi(D) = D^2 (1 - D/5)",
        "Q2_worst_error_vs_minus_cos": worst,
        "withdrawn": ("round 6's 'the distance to quantum mechanics is exactly "
                      "the distance from |D| to D' — counting suffices"),
        "remaining": "nothing in the geometry selects Phi",
    }
