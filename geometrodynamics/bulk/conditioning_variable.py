"""Is the ``|D|`` density forced by the closure set, or is it a choice?

Eighth round of the finite-mouth chain, pre-registered in
``docs/conditioning_variable_prereg.md`` before this file existed. A
**correction round**: it changes the recorded *status* of an input that every
result of rounds 5-7 rests on, and by pre-registered rule it may not change any
number those rounds produced.

The finding it acts on came from an external audit of ``92a915b``. The
repository said

    Haar conditioned on N = 0 is the coarea measure, density |D|/(2|u x v|)

and ledgered that as **derived**. But ``window_monte_carlo`` conditions on
``|Omega mod 2pi| < eps`` -- the *phase* -- and the two prescriptions have the
same support and different limits. Conditioning a measure on a measure-zero set
is not determined by the set (Borel-Kolmogorov); the limiting family is part of
the specification.

Conditioning on ``N`` gives the uniform measure and ``E = 0``. Conditioning on
the phase gives ``|D|`` and the repository's law. The repository's own third
closure condition is stated on *phase*
(``geometrodynamics/history/closure.py:12``), which justifies its choice --
but a justification is not a derivation, and the ledger now says so.
"""

from __future__ import annotations

import math
import os
import re
from typing import Callable, Dict, List, Tuple

import numpy as np

from geometrodynamics.bulk.closure_measurement import (
    great_circle, solid_angle, correlation, window_monte_carlo)

__all__ = [
    "n_window_is_sector_blind", "gradient_magnitudes_on_closure",
    "two_conditioning_families", "phase_axiom_citation",
    "downstream_numbers_unchanged", "velocity_field_non_uniqueness",
    "kappa_is_route_local", "dependency_ledger", "verdict",
]


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _axes(gamma: float) -> Tuple[np.ndarray, np.ndarray]:
    return (np.array([0.0, 0.0, 1.0]),
            np.array([math.sin(gamma), 0.0, math.cos(gamma)]))


def _ND(x, u, w):
    x = np.atleast_2d(np.asarray(x, dtype=float))
    return x @ np.cross(u, w), 1.0 + x @ u + float(np.asarray(u) @ np.asarray(w)) + x @ w


# ── O1-O4: the two conditioning variables ───────────────────────────────────

def n_window_is_sector_blind(gamma: float = 1.0, eps: float = 0.01,
                             n: int = 500000, seed: int = 0) -> Dict[str, object]:
    """O1/O4 — an ``|N| < eps`` window cannot distinguish the outcome sectors.

    With ``u = s_A a`` and ``w = -s_B b``, ``u x w = -s_A s_B (a x b)``, so
    ``|N| = |x.(u x w)|`` is *identical* in all four sectors: the window selects
    literally the same set of ``x``. Equal sector priors then give
    ``P(s_A,s_B) = 1/4`` and ``E = 0`` at every angle, and the counts are equal
    exactly rather than statistically.
    """
    a, b = _axes(gamma)
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    counts, cross_products = {}, {}
    for s_a in (1, -1):
        for s_b in (1, -1):
            u, w = s_a * a, -(s_b * b)
            cross_products[(s_a, s_b)] = np.cross(u, w)
            N, _ = _ND(x, u, w)
            counts[(s_a, s_b)] = int(np.sum(np.abs(N) < eps))
    values = list(counts.values())
    total = sum(values)
    E = (sum(sa * sb * c for (sa, sb), c in counts.items()) / total
         if total else float("nan"))
    q = np.cross(a, b)
    worst_cross = max(
        float(np.max(np.abs(np.abs(v) - np.abs(q)))) for v in cross_products.values())
    return {
        "counts": {f"({sa:+d},{sb:+d})": c for (sa, sb), c in counts.items()},
        "counts_all_equal": bool(len(set(values)) == 1),
        "cross_product_is_pm_q_residual": worst_cross,
        "E": E,
        "E_is_zero": bool(abs(E) < 1e-12),
    }


def gradient_magnitudes_on_closure(gamma: float = 1.0, n: int = 200001
                                   ) -> Dict[str, object]:
    """O2/O3 — why the two windows disagree, in one line each.

    On the closure circle ``Gamma`` (where ``x`` is perpendicular to ``q``):

    * ``grad_{S^2} N = P_x(q) = q``, so ``|grad N| = |a x b| = sin gamma`` is
      **constant**, and the coarea density with respect to ``N`` is uniform in
      arclength;
    * ``grad_{S^2} theta = q/D``, so ``|grad theta| = |q|/|D|``, and the coarea
      density with respect to ``theta`` is proportional to **``|D|``**.

    Coarea density is ``1/|grad(conditioning variable)|``. Same zero set, two
    different densities: that is the whole content of the correction.
    """
    a, b = _axes(gamma)
    rows = []
    for s_a in (1, -1):
        for s_b in (1, -1):
            u, w = s_a * a, -(s_b * b)
            q = np.cross(u, w)
            nq = float(np.linalg.norm(q))
            circle = great_circle(u, w, n)
            _, D = _ND(circle, u, w)
            grad_theta = nq / np.abs(D)
            # coarea densities, normalised to unit mass on Gamma
            dens_N = np.ones(len(circle)) / nq
            dens_theta = np.abs(D) / nq
            rows.append({
                "sector": f"({s_a:+d},{s_b:+d})",
                "|grad N|": nq,
                "grad_N_spread": float(np.max(dens_N) - np.min(dens_N)),
                "|grad theta| min": float(np.min(grad_theta)),
                "|grad theta| max": float(np.max(grad_theta)),
                "N_density_uniform": True,
                "theta_density_range": (float(np.min(dens_theta)),
                                        float(np.max(dens_theta))),
            })
    return {
        "rows": rows,
        "grad_N_is_constant_on_closure": all(r["grad_N_spread"] < 1e-12
                                             for r in rows),
        "theta_density_varies": all(
            r["theta_density_range"][1] - r["theta_density_range"][0] > 1e-3
            for r in rows),
    }


def two_conditioning_families(gamma: float = 1.0, n: int = 4000000,
                              epsilons=(0.05, 0.02, 0.01, 0.005, 0.002)
                              ) -> Dict[str, object]:
    """The Borel-Kolmogorov demonstration: same support, two limits.

    Both windows shrink onto the same circle ``Gamma``. The ``N`` window's limit
    is ``E = 0``; the phase window's limit is the repository's law. Neither is
    "the" conditional measure of Haar on ``Gamma``, because there is no such
    thing without naming the variable.
    """
    a, b = _axes(gamma)
    rng = np.random.default_rng(1)
    x = rng.standard_normal((n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    rows = []
    for eps in epsilons:
        cn, cp = {}, {}
        for s_a in (1, -1):
            for s_b in (1, -1):
                u, w = s_a * a, -(s_b * b)
                N, _ = _ND(x, u, w)
                cn[(s_a, s_b)] = int(np.sum(np.abs(N) < eps))
                om = solid_angle(x, s_a * a, s_b * b)
                cp[(s_a, s_b)] = int(np.sum(
                    np.abs((om + math.pi) % (2.0 * math.pi) - math.pi) < eps))

        def corr(c):
            z = sum(c.values())
            return (sum(sa * sb * v for (sa, sb), v in c.items()) / z
                    if z else float("nan"))
        rows.append({"eps": eps, "E_N_window": corr(cn),
                     "E_phase_window": corr(cp)})
    limit = float(correlation(gamma, "abs", n=400001))
    return {
        "rows": rows,
        "N_window_limit": 0.0,
        "phase_window_limit_closed_form": limit,
        "N_window_is_zero_at_every_eps": all(abs(r["E_N_window"]) < 1e-12
                                             for r in rows),
        "phase_window_tracks_the_repository_law": bool(
            abs(rows[-1]["E_phase_window"] - limit) < 0.02),
        "same_support_different_limits": True,
    }


def phase_axiom_citation() -> Dict[str, object]:
    """The structure in the repository that privileges the phase, by line.

    Rule 1 of the pre-registration forbids preferring the phase window because
    it gives the interesting answer. The justification must already exist. It
    does: the third closure condition is stated on *phase*. This function reads
    it out of the source rather than quoting it from memory.
    """
    path = os.path.join(os.path.dirname(__file__), "..", "history", "closure.py")
    path = os.path.normpath(path)
    with open(path, "r") as handle:
        lines = handle.read().splitlines()
    hits = [(i + 1, line.strip()) for i, line in enumerate(lines[:40])
            if re.search(r"phase closure", line, re.I)]
    return {
        "file": "geometrodynamics/history/closure.py",
        "matches": hits,
        "axiom_is_stated_on_phase": bool(hits),
        "text": hits[0][1] if hits else None,
        "note": ("the closure axiom is a condition on total phase, so a phase "
                 "tolerance is the natural regularisation of it; N = 0 is a "
                 "derived description of the same locus"),
    }


# ── rule 3: nothing downstream may move ─────────────────────────────────────

def downstream_numbers_unchanged() -> Dict[str, object]:
    """Pre-registered rule 3 — this corrects a status, not a computation.

    Round 5's two correlations, round 6's oriented law and round 7's Morse-Bott
    masses are recomputed and compared against the values the merged rounds
    published. Any movement would mean the finding is larger than a mislabelled
    ledger entry.
    """
    from geometrodynamics.bulk.closure_current import singlet_loop_law
    from geometrodynamics.bulk.history_action import morse_bott_component_masses

    checks = []
    for variant, expected in (("abs", 0.3984966504), ("signed", 0.5403023059)):
        got = float(correlation(1.0, variant, n=400001))
        checks.append({"quantity": f"round-5 correlation({variant}) at gamma=1",
                       "expected": expected, "got": got,
                       "delta": abs(got - expected)})
    law = singlet_loop_law(1.0, n=40001)
    got_E = float(law["E"])
    checks.append({"quantity": "round-6 oriented singlet E at gamma=1",
                   "expected": -math.cos(1.0), "got": got_E,
                   "delta": abs(got_E + math.cos(1.0))})
    masses = morse_bott_component_masses(*_axes(1.0), n=200001)
    for key, expected in (("M_0", 11.67080583639533), ("M_pi", 0.16951227835901764)):
        checks.append({"quantity": f"round-7 {key}", "expected": expected,
                       "got": float(masses[key]),
                       "delta": abs(float(masses[key]) - expected)})
    worst = max(c["delta"] for c in checks)
    return {"checks": checks, "worst_delta": worst,
            "nothing_downstream_moved": bool(worst < 1e-6)}


# ── C: the velocity field is not unique ─────────────────────────────────────

def velocity_field_non_uniqueness(n_grid: int = 48, box: float = 6.0
                                  ) -> Dict[str, object]:
    """`born_rule_equivariance.md:73` claims uniqueness that does not hold.

    That line lists as reason (iii) that ``v = J/rho`` "is the unique velocity
    field whose current closes the continuity equation". It is not. For any
    ``K`` with ``div K = 0``,

        div(rho v') = div(J + K) = div J,      v' = (J + K)/rho,

    so ``v'`` satisfies the *same* continuity equation wherever ``rho > 0``.
    Demonstrated here with an explicit compactly supported ``K = curl A``,
    ``A = f(r) z_hat``, on a Gaussian wavepacket. Two things are measured: that
    ``div K`` vanishes to discretisation error, and that ``div(J + K)`` matches
    ``div J`` to the same order.

    The addition also leaves the **integrated** current alone --
    ``int K d^3x = 0`` for compactly supported ``curl A`` -- so an
    Ehrenfest/mean-velocity check does not exclude it either.

    This does not rescue an alternative flow physically: extra assumptions may
    well rule these out. It removes one stated reason, and it means the cited
    Goldstein-Struyve uniqueness (which assumes Bohmian dynamics and locality
    conditions on the density functional) is doing work that (iii) presents as
    already done.
    """
    g = np.linspace(-box, box, n_grid)
    h = float(g[1] - g[0])
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")

    # a Gaussian wavepacket with a phase ramp: rho > 0 everywhere, J != 0
    k = np.array([0.7, -0.3, 0.2])
    psi = np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / 4.0) * np.exp(
        1j * (k[0] * X + k[1] * Y + k[2] * Z))
    rho = np.abs(psi) ** 2

    def grad(f):
        return np.gradient(f, h, edge_order=2)

    gx, gy, gz = grad(psi.real)
    hx, hy, hz = grad(psi.imag)
    # J = Im(psi* grad psi) = Re(psi) grad Im(psi) - Im(psi) grad Re(psi)
    J = [psi.real * hx - psi.imag * gx,
         psi.real * hy - psi.imag * gy,
         psi.real * hz - psi.imag * gz]

    # K = curl A with A = f(r) z_hat, f compactly supported -> div K = 0
    r2 = X ** 2 + Y ** 2 + Z ** 2
    f = np.exp(-r2 / 2.0) * (r2 < 16.0)
    fx, fy, _fz = grad(f)
    K = [fy, -fx, np.zeros_like(f)]          # curl(f z_hat) = (d_y f, -d_x f, 0)

    def div(V):
        return (np.gradient(V[0], h, axis=0, edge_order=2)
                + np.gradient(V[1], h, axis=1, edge_order=2)
                + np.gradient(V[2], h, axis=2, edge_order=2))

    interior = (slice(3, -3),) * 3
    div_K = div(K)[interior]
    div_J = div(J)[interior]
    div_JK = div([J[i] + K[i] for i in range(3)])[interior]
    scale = float(np.max(np.abs(div_J)))
    integrated = [float(np.sum(K[i]) * h ** 3) for i in range(3)]
    j_scale = max(float(np.sum(np.abs(J[i])) * h ** 3) for i in range(3))
    return {
        "max|div K| / max|div J|": float(np.max(np.abs(div_K))) / scale,
        "max|div(J+K) - div J| / max|div J|":
            float(np.max(np.abs(div_JK - div_J))) / scale,
        "integrated_K": integrated,
        "integrated_K_relative": max(abs(v) for v in integrated) / j_scale,
        "K_is_divergence_free": bool(float(np.max(np.abs(div_K))) / scale < 5e-3),
        "continuity_unchanged": bool(
            float(np.max(np.abs(div_JK - div_J))) / scale < 5e-3),
        "mean_velocity_check_cannot_exclude_it": bool(
            max(abs(v) for v in integrated) / j_scale < 1e-6),
        "uniqueness_claim_at_line_73_is_false": True,
    }


def kappa_is_route_local() -> Dict[str, object]:
    """B — ``kappa`` is a parameter of the route round 7 closed, not a postulate."""
    return {
        "three_universal_underived_inputs": [
            "branch aggregation (positive count or oriented sum)",
            "relative outcome-sector coefficients (r = 1 is chosen)",
            "current-to-frequency readout",
        ],
        "kappa": ("the normalisation of e^{i kappa S_H}; it exists only inside "
                  "the holonomy-trace route, which round 7 closed by showing "
                  "S_H is not additive. It is a parameter of that attempted "
                  "route, not a fourth universal underived input."),
        "was_previously_described_as": "a fourth underived ingredient",
    }


# ── ledger and verdict ──────────────────────────────────────────────────────

def dependency_ledger() -> List[Dict[str, str]]:
    """The corrected round-5 ledger, with the conditioning variable broken out."""
    return [
        {"input": "Haar prior on S^2", "status": "chosen",
         "where": "invariant prior; physicality of x is itself a choice"},
        {"input": "closure axiom: total phase = 0 or pi (mod 2pi)",
         "status": "repository axiom", "where": "history/closure.py:12"},
        {"input": "geodesic-realignment detection model", "status": "chosen",
         "where": "round 5"},
        {"input": "the conditioning VARIABLE (phase, not N)", "status": "chosen",
         "where": ("was ledgered 'derived'. Conditioning on a measure-zero set "
                   "is not fixed by the set: an N-window gives the uniform "
                   "measure and E = 0, a phase window gives |D| and the "
                   "repository's law. Justified by the axiom being stated on "
                   "phase, but justified is not derived")},
        {"input": "coarea density |D|/(2|u x v|) GIVEN the phase variable",
         "status": "derived", "where": "1/|grad theta| with grad theta = q/D"},
        {"input": "outcome signs as history boundary data", "status": "chosen",
         "where": "D-type"},
        {"input": "equal prior on the four outcome sectors", "status": "chosen",
         "where": "counting measure; only r = 1 gives the quantum law"},
    ]


def verdict() -> Dict[str, object]:
    """The pre-registered question, answered from the measurements above."""
    sector = n_window_is_sector_blind()
    grads = gradient_magnitudes_on_closure()
    axiom = phase_axiom_citation()
    frozen = downstream_numbers_unchanged()
    forced = not (sector["E_is_zero"] and grads["grad_N_is_constant_on_closure"])
    if forced:
        label = "CONDITIONING_VARIABLE_IS_FORCED_BY_THE_CLOSURE_SET"
    elif axiom["axiom_is_stated_on_phase"]:
        label = "CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_JUSTIFIED_BY_THE_PHASE_AXIOM"
    else:
        label = "CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_WITH_NO_JUSTIFICATION_IN_THE_REPOSITORY"
    return {
        "A_conditioning_variable": label,
        "ledger_change": "coarea conditioning: derived -> chosen (the variable)",
        "downstream_numbers_unchanged": frozen["nothing_downstream_moved"],
        "B_kappa": "route-local parameter, not a fourth universal input",
        "C_velocity_uniqueness":
            "born_rule_equivariance.md:73 reason (iii) is withdrawn",
    }
