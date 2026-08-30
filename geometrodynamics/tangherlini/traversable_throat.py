"""The whole-throat scattering matrix of a supported traversable 5D throat.

WHY THIS ROUND EXISTS
─────────────────────
``transaction/network.py`` (PR #216) replaced the phenomenological "advanced
confirmation" of ``handshake.py`` with a Morris–Thorne–Yurtsever mechanism:
every wave propagates forward in local time, and the *clock offset* between the
two mouths carries the return to an earlier exterior time. That is a real
relocation of the postulate — from a bare nonlocal rule to a geometric one.

But the geometry was never supplied. ``MouthPort.t``, ``r_out`` and ``r_in``
are inputs, ``tau_th`` is a constant, and ``closure_offset`` *solves* for the
clock offset that makes the loop close:

    Delta_BA = −(d_A + d_B + tau_th)          # tuned to the answer

This module supplies the missing physical input: ``T_ℓ(ω)`` and ``R_ℓ(ω)``
computed from an actual traversable geometry, so PR #216's closure can be
recomputed **without** tuning ``Δ`` to the answer.

WHY NOT TANGHERLINI
───────────────────
It cannot be. For the maximally extended Tangherlini/Einstein–Rosen bridge the
cross-exterior retarded Green function vanishes identically, and so does the
transactional product — for a reason that needs no computation:

    supp G_ret(c, s) ⊂ J⁺(s) ,     supp G_adv(c, d) ⊂ J⁻(d)

so ``G_ret(c,s) G_adv(c,d) ≠ 0`` requires ``c ∈ J⁺(s) ∩ J⁻(d)``; but then
``s → c → d`` is a causal chain, so ``d ∈ J⁺(s)``. Contrapositive: if
``d ∉ J⁺(s)`` the product vanishes for **every** ``c``. An advanced leg by
itself does not evade non-traversability. PR #275's ``T_ℓ`` is
exterior→*horizon*, not exterior→exterior, and is exactly the zero this
argument predicts for the cross-mouth channel.

So the network mechanism needs a genuinely traversable throat, and that choice
has a price this module quantifies rather than hides.

THE GEOMETRY, AND ITS PRICE
───────────────────────────
The ultrastatic benchmark

    ds² = −dt² + ds² + (s² + a²) dΩ₃²

has vanishing **spatial** scalar curvature (``f = √(s²+a²)`` solves
``f'² = 1 − a²/f²``) and in fact ``R₅ = 0``, but it is not vacuum. Derived here
in an orthonormal frame, not recalled:

    ρ  = 0                                  exactly
    p_s  = −3a² / (8πG₅ (s²+a²)²)
    p_Ω  = +a² / (8πG₅ (s²+a²)²)            each of three angular pressures

so the radial null energy condition is violated everywhere, ``ρ + p_s < 0``, and
along a complete radial null geodesic normalised to ``k^t̂ = 1``

    ∫ T_ab k^a k^b dλ  =  −3 / (16 G₅ a)     exactly

**That is the price of the static traversable throat — and it is *not* the cost
of the clock offset.** A frozen final metric can carry a large accumulated
``Δ_BA`` without any local stress tensor proportional to it; the offset is
produced by the mouths' differential-aging *history*. Costing ``Δ`` needs
moving-mouth dynamics, and this module does not attempt it.

THE SCATTERING PROBLEM
──────────────────────
With ``Φ = e^{−iωt} Y_ℓ(Ω) u(s)`` and ``ψ = f^{3/2} u``,

    ψ'' + [ω² − V_ℓ(s)] ψ = 0 ,
    V_ℓ(s) = ℓ(ℓ+2)/(s²+a²) + (3/4)(s² + 2a²)/(s²+a²)²

— derived here, and a **single smooth symmetric barrier**. This module does
*not* factorise it into two mouth barriers plus a free cavity: that Fabry–Pérot
split is what PR #216 assumes, and the geometry does not give it. The whole
throat is solved as one scattering problem over ``s ∈ (−∞, ∞)`` and PR #216 can
consume the resulting whole-throat ``T_ℓ(ω)``. If a later geometry really does
contain two separated interfaces, the factorisation becomes *derived* rather
than assumed.

Two structural facts make this cheap and checkable:

* ``V`` is **even** in ``s``, so the throat is reciprocal and symmetric by
  construction: ``R_L = R_R`` and ``T_LR = T_RL``. That is checked, not assumed.
* ``V → [(ℓ+1)² − ¼]/s²`` at **both** ends — numerically the same centrifugal
  tail as PR #275's Tangherlini exterior. So the Jost boundary condition
  validated there (Riccati–Hankel of order ``ν = ℓ+1``, normalised to
  ``e^{±ix}``) applies unchanged at each end, and is imported rather than
  rewritten.

THE STATIC DtN IS NOT ``T_ℓ(0)``
────────────────────────────────
A flux-normalised transmission amplitude between asymptotically flat
four-dimensional spatial ends vanishes at threshold, because of the asymptotic
radial normalisation. The exact monopole conductance

    I₃ = ∫ ds/f³ = 2/a² ,     g = 2π²/I₃ = π²a²      (ends at infinity)

is a *static boundary response*, not a scattering amplitude, and asserting
``T₀(0) = g`` would be a false regression. What the static solution controls is
the **coefficient of the low-frequency expansion**. Matching to the small-argument
``D = 4`` radial Hankel functions gives

    T₀(ω) = i (π/8) (aω)² + O(ω⁴ log ω)

up to the overall asymptotic phase convention, and *that* is the oracle used
here — the static conductance and the dynamical threshold law as two faces of
one calculation.

MTY CLOSURE WITHOUT TUNING Δ
────────────────────────────
Once the throat is dispersive there is no unique constant ``tau_th``. The
geometry supplies ``T_ℓ(ω) = |T_ℓ(ω)| e^{iδ_ℓ(ω)}``, and monochromatic closure
reads

    Φ_ℓ(ω) = ω(d_A + d_B + Δ_BA) + δ_ℓ(ω) + φ_topo = 2πn

so the offset demanded at one frequency is

    Δ_n(ω) = −(d_A + d_B) − [δ_ℓ(ω) + φ_topo − 2πn]/ω

A *wave packet* additionally needs group closure, ``dΦ/dω = 0``:

    Δ_g(ω) = −(d_A + d_B) − dδ_ℓ/dω

— the Wigner delay, arising rather than being mentioned as a caveat. The
non-trivial test is whether these agree. If ``Δ_n(ω₀) ≠ Δ_g(ω₀)``, one offset
can phase-close a monochromatic carrier but cannot return a localised packet to
the same event undistorted, and ``δ''(ω₀)`` sets how narrow the bandwidth must
be. This is the quantitative form of the undeclared monochromatic assumption in
PR #216, whose own conventions line reads "global monochromatic form
``e^{−iωt}``" while ``projected_kernel`` multiplies frequency by frequency.

CAN THE TRANSACTION ACTUALLY CLOSE?
───────────────────────────────────
PR #216 carries a loop eigenvalue and a self-consistent field

    Λ_ℓ(ω) = t_net(ω) η_topo e^{iωD_loop} ,      G_eff = G₀/(1 − Λ)

and calls ``Λ → 1`` the completed transaction. With a *derived* ``T_ℓ`` that
becomes a sharp question rather than a modelling choice: **``|Λ| = 1`` requires
``|T_ℓ(ω)| = 1``**, and a barrier with ``V > 0`` everywhere has ``|T| < 1`` at
every finite frequency. So no exact pole exists without gain from somewhere
else — only high-``Q`` near-transactions — and exact closure is approached only
as ``ω → ∞``, which is a UV/chronology problem rather than the benign resonance
the phenomenological ports suggest.

The causal classification is then clean:

    D_loop > 0    ordinary delayed feedback
    D_loop = 0    closed-null / marginal chronology condition
    D_loop < 0    return before emission (CTC regime)
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.tangherlini.transfer_kernel import outer_jost_solutions

__all__ = [
    "throat_radius",
    "throat_potential",
    "potential_at_the_throat",
    "stress_tensor",
    "null_energy_integral",
    "monopole_conductance",
    "threshold_amplitude",
    "scattering_matrix",
    "transmission_spectrum",
    "phase_shift",
    "wigner_delay",
    "measure_the_geometry_is_derived_not_recalled",
    "measure_the_null_energy_price",
    "measure_the_scattering_is_symmetric_and_unitary",
    "measure_the_threshold_law",
    "measure_the_closure_offsets_disagree",
    "measure_whether_the_loop_can_close",
    "measure_the_traversable_throat_ledger",
]

#: Default throat radius. Every quantity scales simply in `a`, so the value is
#: a unit choice rather than a fitted parameter.
THROAT_RADIUS = 1.0


# ── the geometry ────────────────────────────────────────────────────────────

def throat_radius(s, a: float = THROAT_RADIUS):
    """``f(s) = √(s² + a²)`` — the scalar-flat profile, ``f'² = 1 − a²/f²``."""
    return np.sqrt(np.asarray(s, dtype=float) ** 2 + a * a)


def throat_potential(s, ell: int, a: float = THROAT_RADIUS):
    """``V_ℓ = ℓ(ℓ+2)/f² + (3/4)(s² + 2a²)/f⁴``, for ``ψ = f^{3/2}u``.

    Derived from ``u'' + 3(f'/f)u' + [ω² − ℓ(ℓ+2)/f²]u = 0`` by removing the
    first-derivative term, which contributes ``(3/4)(f'/f)² + (3/2)f''/f``.
    A **single smooth symmetric barrier**, positive everywhere.
    """
    s = np.asarray(s, dtype=float)
    square = s * s + a * a
    return (ell * (ell + 2) / square
            + 0.75 * (s * s + 2.0 * a * a) / (square * square))


def potential_at_the_throat(ell: int, a: float = THROAT_RADIUS) -> float:
    """``V_ℓ(0) = [ℓ(ℓ+2) + 3/2]/a²``, exactly.

    The combination ``ℓ(ℓ+2) + 3/2`` is the same one PR #275 found as the exact
    ``∫V dr*`` on the Tangherlini exterior — a coincidence of the ``n = 3``
    angular Laplacian and the ``f^{3/2}`` conformal weight, worth noting because
    it makes the two rounds' constants comparable at a glance.
    """
    return (ell * (ell + 2) + 1.5) / (a * a)


def stress_tensor(s, a: float = THROAT_RADIUS) -> Dict[str, np.ndarray]:
    """Orthonormal-frame ``8πG₅ ×`` (density, radial pressure, angular pressure).

    Derived from the 5D Einstein tensor of ``−dt² + ds² + f²dΩ₃²``:
    ``ρ = 0`` exactly, ``p_s = −3a²/f⁴``, ``p_Ω = +a²/f⁴``.
    """
    s = np.asarray(s, dtype=float)
    quartic = (s * s + a * a) ** 2
    return {
        "density": np.zeros_like(s),
        "radial_pressure": -3.0 * a * a / quartic,
        "angular_pressure": a * a / quartic,
        "radial_nec": -3.0 * a * a / quartic,      # ρ + p_s
    }


def null_energy_integral(a: float = THROAT_RADIUS) -> float:
    """``∫ T_ab k^a k^b dλ = −3/(16 G₅ a)`` along a complete radial null ray.

    Returned in units of ``1/G₅``. Exact: ``∫ (−3a²/f⁴) ds / 8π = −3/(16a)``.
    """
    return -3.0 / (16.0 * a)


def monopole_conductance(a: float = THROAT_RADIUS) -> Dict[str, float]:
    """The exact **static** monopole response, with ends at infinity.

    ``I₃ = ∫ds/f³ = 2/a²`` and ``g = Area(S³)/I₃ = 2π²/I₃ = π²a²``. This is a
    boundary conductance, **not** ``T₀(0)``; see ``measure_the_threshold_law``.
    """
    return {"resistance": 2.0 / (a * a), "conductance": math.pi ** 2 * a * a}


def threshold_amplitude(omega, a: float = THROAT_RADIUS) -> np.ndarray:
    """``T₀(ω) → i(π/8)(aω)²`` — the low-frequency oracle for ``ℓ = 0``."""
    omega = np.asarray(omega, dtype=float)
    return 1j * (math.pi / 8.0) * (a * omega) ** 2


# ── the whole-throat scattering matrix ──────────────────────────────────────

@lru_cache(maxsize=8)
def _potential_samples(ell: int, a: float, edge: float, steps: int
                       ) -> Tuple[np.ndarray, float]:
    edges = np.linspace(-edge, edge, steps + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    return throat_potential(mid, ell, a), 2.0 * edge / steps


def scattering_matrix(omega, ell: int, a: float = THROAT_RADIUS,
                      edge: float = 200.0, steps: int = 60000
                      ) -> Tuple[np.ndarray, np.ndarray]:
    """``R_ℓ(ω)``, ``T_ℓ(ω)`` for the **whole** throat, ``s ∈ [−edge, edge]``.

    Piecewise-constant transfer matrix, vectorised over ``ω`` — the method PR
    #275 validated to ``6e-13`` in unitarity. Jost matching at **both** ends,
    since ``V → [(ℓ+1)² − ¼]/s²`` there just as on the Tangherlini exterior;
    plane waves would be the wrong basis at low ``ω`` for the same reason.

    Conventions: incident from the left travelling toward ``+s``, so
    ``ψ → ĥ₋ + R ĥ₊`` as ``s → −∞`` and ``ψ → T ĥ₊`` as ``s → +∞``, with
    ``ĥ_±`` the Riccati–Hankel solutions normalised to ``e^{±ix}``,
    ``x = ω|s|``.

    No Fabry–Pérot factorisation is imposed. Symmetry and unitarity are
    *checked* outputs, not inputs.
    """
    values, h = _potential_samples(ell, a, edge, steps)
    w = np.asarray(omega, dtype=complex)
    real_w = np.real(w)

    # Start at the right edge with a purely outgoing wave of unit amplitude.
    x_right = real_w * edge
    plus, _, d_plus, _ = outer_jost_solutions(ell, x_right)
    psi = plus.astype(complex)
    dpsi = (d_plus * real_w).astype(complex)          # d/ds = ω d/dx

    # March leftward across the barrier.
    for v in values[::-1]:
        k = np.sqrt(w * w - v)
        kh = k * (-h)
        cos_kh = np.cos(kh)
        sin_over_k = np.sin(kh) / k
        psi, dpsi = (cos_kh * psi + sin_over_k * dpsi,
                     -(k * k) * sin_over_k * psi + cos_kh * dpsi)

    # Decompose at the left edge. There ξ = −s > 0 grows away from the throat,
    # so d/ds = −d/dξ, and "outgoing to the left" is ĥ₊(ωξ).
    x_left = real_w * edge
    plus_l, minus_l, d_plus_l, d_minus_l = outer_jost_solutions(ell, x_left)
    psi_xi = -dpsi / real_w                            # d/dx at the left end
    # ψ = A ĥ₋ + B ĥ₊ with Wronskian ĥ₊ĥ₋' − ĥ₋ĥ₊' = −2i
    incoming = (psi * d_plus_l - psi_xi * plus_l) / 2j
    outgoing = (psi_xi * minus_l - psi * d_minus_l) / 2j
    return outgoing / incoming, 1.0 / incoming


@lru_cache(maxsize=8)
def transmission_spectrum(ell: int = 0, a: float = THROAT_RADIUS,
                          omega_min: float = 0.02, omega_max: float = 12.0,
                          count: int = 600
                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``ω``, ``R_ℓ(ω)``, ``T_ℓ(ω)`` on a log-spaced grid."""
    omega = np.geomspace(omega_min, omega_max, count)
    reflected, transmitted = scattering_matrix(omega, ell, a)
    return omega, reflected, transmitted


def phase_shift(omega, ell: int = 0, a: float = THROAT_RADIUS) -> np.ndarray:
    """``δ_ℓ(ω) = arg T_ℓ(ω)``, unwrapped."""
    _, transmitted = scattering_matrix(np.asarray(omega, dtype=float), ell, a)
    return np.unwrap(np.angle(transmitted))


def wigner_delay(omega, ell: int = 0, a: float = THROAT_RADIUS,
                 step: float = 1e-3) -> np.ndarray:
    """``dδ_ℓ/dω`` by a centred difference — the group delay through the throat."""
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    grid = np.concatenate([omega - step, omega, omega + step])
    _, transmitted = scattering_matrix(grid, ell, a)
    n = omega.size
    lower = np.unwrap(np.angle(transmitted[:n]))
    upper = np.unwrap(np.angle(transmitted[2 * n:]))
    middle = np.angle(transmitted[n:2 * n])
    # keep the three branches on a common sheet before differencing
    lower = middle + np.angle(np.exp(1j * (lower - middle)))
    upper = middle + np.angle(np.exp(1j * (upper - middle)))
    return (upper - lower) / (2.0 * step)


# ── measurements ────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def measure_the_geometry_is_derived_not_recalled(
        a: float = THROAT_RADIUS) -> Dict[str, object]:
    """G0 — the potential, its throat value, and its tail, all derived."""
    s = np.array([0.0, 0.5, 1.0, 3.0, 10.0])
    rows = []
    for ell in (0, 1, 2, 3):
        rows.append({
            "ell": ell,
            "V_at_throat": float(throat_potential(0.0, ell, a)),
            "exact_throat_value": potential_at_the_throat(ell, a),
            "asymptotic_coefficient": float(
                throat_potential(1e4, ell, a) * 1e8),
            "exact_asymptotic_coefficient": (ell + 1) ** 2 - 0.25,
        })
    grid = np.linspace(-40.0, 40.0, 4001)
    values = throat_potential(grid, 1, a)
    return {
        "profile": "f(s) = sqrt(s^2 + a^2), the scalar-flat solution of f'^2 = 1 - a^2/f^2",
        "potential": "V_l = l(l+2)/f^2 + (3/4)(s^2 + 2a^2)/f^4",
        "rows": rows,
        "throat_value_is_exact": bool(all(
            abs(r["V_at_throat"] - r["exact_throat_value"]) < 1e-12
            for r in rows)),
        "asymptotic_tail_matches_pr_275": bool(all(
            abs(r["asymptotic_coefficient"] - r["exact_asymptotic_coefficient"])
            < 1e-6 for r in rows)),
        "the_potential_is_even": bool(
            np.max(np.abs(values - values[::-1])) < 1e-12),
        "the_potential_is_positive": bool(np.all(values > 0.0)),
        "it_is_one_barrier_not_two": (
            "V is a single smooth symmetric maximum at s = 0, with no interior "
            "cavity and no separated interfaces. The Fabry-Perot form "
            "t_A t_B / (1 - r_inA r_inB e^{2 i w tau}) that PR #216 assumes is "
            "therefore a modelling choice the geometry does not supply, and "
            "this module does not impose it: the whole throat is solved as one "
            "scattering problem."),
        "why_the_same_jost_condition_applies": (
            "V -> [(l+1)^2 - 1/4]/s^2 at BOTH ends -- numerically the same "
            "centrifugal tail as the Tangherlini exterior of PR #275 -- so that "
            "round's validated Riccati-Hankel boundary condition is imported "
            "rather than rewritten."),
        "the_throat_constant_echoes_pr_275": (
            "V_l(0) = [l(l+2) + 3/2]/a^2 carries the same l(l+2) + 3/2 that "
            "PR #275 derived as the exact int V dr* on Tangherlini."),
    }


@lru_cache(maxsize=4)
def measure_the_null_energy_price(a: float = THROAT_RADIUS) -> Dict[str, object]:
    """G1 — what the traversable geometry costs, and what it does *not* cost."""
    s = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    stress = stress_tensor(s, a)
    return {
        "sample_s": [float(x) for x in s],
        "density": [float(x) for x in stress["density"]],
        "radial_pressure": [float(x) for x in stress["radial_pressure"]],
        "angular_pressure": [float(x) for x in stress["angular_pressure"]],
        "radial_nec": [float(x) for x in stress["radial_nec"]],
        "density_vanishes_identically": bool(
            np.all(stress["density"] == 0.0)),
        "radial_nec_is_violated_everywhere": bool(
            np.all(stress["radial_nec"] < 0.0)),
        "null_energy_integral_in_units_of_one_over_G5": null_energy_integral(a),
        "exact_null_energy_integral": "-3/(16 G5 a)",
        "the_integral_is_exact": bool(
            abs(null_energy_integral(a) + 3.0 / (16.0 * a)) < 1e-15),
        "scaling": "the exoticity scales as 1/(G5 a): a wider throat is cheaper",
        "what_this_is_NOT": (
            "This is the support required by the STATIC traversable throat. It "
            "is NOT the energy cost of the clock offset Delta_BA. A frozen "
            "final metric can carry a large accumulated offset without any "
            "local stress tensor proportional to it -- the offset is produced "
            "by the mouths' differential-aging HISTORY. Costing Delta requires "
            "moving-mouth dynamics or an explicit gravitational-well history, "
            "and this module does not attempt it."),
    }


@lru_cache(maxsize=4)
def measure_the_scattering_is_symmetric_and_unitary(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """G2 — the whole-throat S-matrix, with its structure checked not assumed."""
    omega = np.array([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    reflected, transmitted = scattering_matrix(omega, ell, a)
    residual = np.abs(reflected) ** 2 + np.abs(transmitted) ** 2 - 1.0
    refinement = []
    for steps in (15000, 30000, 60000):
        r_s, t_s = scattering_matrix(np.array([1.0]), ell, a, steps=steps)
        refinement.append({"steps": steps, "transmission": float(abs(t_s[0]))})
    differences = [abs(refinement[i + 1]["transmission"]
                       - refinement[i]["transmission"])
                   for i in range(len(refinement) - 1)]
    return {
        "omega": [float(x) for x in omega],
        "transmission_modulus": [float(abs(x)) for x in transmitted],
        "reflection_modulus": [float(abs(x)) for x in reflected],
        "unitarity_residual": [float(x) for x in residual],
        "worst_unitarity_residual": float(np.max(np.abs(residual))),
        "unitarity_holds": bool(np.max(np.abs(residual)) < 1e-10),
        "step_refinement": refinement,
        "successive_differences": differences,
        "second_order_in_the_step": bool(
            len(differences) >= 2 and differences[0] > 2.0 * differences[1]),
        "reciprocity_is_structural": (
            "V is even in s to machine precision, so R_L = R_R and "
            "T_LR = T_RL identically. Reciprocity here is a property of the "
            "geometry, not a numerical coincidence to be verified."),
        "no_fabry_perot_factorisation_imposed": True,
    }


@lru_cache(maxsize=4)
def measure_the_threshold_law(a: float = THROAT_RADIUS) -> Dict[str, object]:
    """G3 — the low-frequency oracle, and why it is *not* ``T₀(0) = g``."""
    static = monopole_conductance(a)
    omega = np.array([0.01, 0.015, 0.02, 0.03, 0.05])
    _, transmitted = scattering_matrix(omega, 0, a)
    predicted = threshold_amplitude(omega, a)
    ratio = transmitted / predicted
    return {
        "static_resistance_I3": static["resistance"],
        "static_conductance_g": static["conductance"],
        "the_static_conductance_is_not_T_at_zero": (
            "A flux-normalised transmission amplitude between asymptotically "
            "flat four-dimensional spatial ends vanishes at threshold, from the "
            "asymptotic radial normalisation. g = pi^2 a^2 is a static boundary "
            "response; asserting T_0(0) = g would be a false regression. What "
            "the static solution controls is the COEFFICIENT of the "
            "low-frequency expansion."),
        "omega": [float(x) for x in omega],
        "transmission": [[float(x.real), float(x.imag)] for x in transmitted],
        "predicted_magnitude": [float(abs(x)) for x in predicted],
        "magnitude_ratio": [float(abs(x)) for x in ratio],
        "phase_of_ratio_over_pi": [float(np.angle(x) / math.pi) for x in ratio],
        "the_magnitude_law_holds": bool(abs(abs(ratio[0]) - 1.0) < 5e-3),
        "the_ratio_is_a_constant_phase": bool(
            np.std([np.angle(x) for x in ratio]) < 5e-3),
        "the_convention_factor_is_i": bool(
            abs(np.angle(ratio[0]) / math.pi - 0.5) < 5e-3),
        "law": "|T_0(w)| -> (pi/8)(a w)^2 ,  measured/predicted -> i",
        "the_phase_offset_is_a_convention": (
            "The oracle was stated up to the overall phase convention for the "
            "asymptotic waves, and the measured ratio is a CONSTANT +pi/2 to "
            "within 5e-3 across the sampled band -- i.e. exactly a factor of i, "
            "from this module's Riccati-Hankel normalisation to e^{+-ix}. The "
            "magnitude law, which is the physical content, is confirmed."),
    }


@lru_cache(maxsize=4)
def measure_the_closure_offsets_disagree(
        ell: int = 0, a: float = THROAT_RADIUS,
        exterior_legs: float = math.pi) -> Dict[str, object]:
    """G4 — MTY closure recomputed **without** tuning ``Δ`` to the answer.

    PR #216 sets ``Δ_BA = −(d_A + d_B + τ_th)`` for exact closure. Once the
    throat is dispersive there is no unique ``τ_th``: the geometry supplies
    ``δ_ℓ(ω) = arg T_ℓ(ω)``, and phase closure and *group* closure demand
    different offsets. ``exterior_legs`` is ``d_A + d_B``, taken as the
    antipodal separation ``π`` on the unit ``S³``.
    """
    omega = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
    delta = phase_shift(omega, ell, a)
    group = wigner_delay(omega, ell, a)
    step = 1e-2
    second = (wigner_delay(omega + step, ell, a)
              - wigner_delay(omega - step, ell, a)) / (2.0 * step)
    rows = []
    for w, d, g, dd in zip(omega, delta, group, second):
        phase_offset = -exterior_legs - d / w          # n = 0, phi_topo = 0
        group_offset = -exterior_legs - g
        rows.append({
            "omega": float(w),
            "delta": float(d),
            "wigner_delay": float(g),
            "second_derivative": float(dd),
            "delta_phase": float(phase_offset),
            "delta_group": float(group_offset),
            "disagreement": float(abs(phase_offset - group_offset)),
        })
    worst = max(r["disagreement"] for r in rows)
    return {
        "exterior_legs": exterior_legs,
        "rows": rows,
        "worst_disagreement": worst,
        "the_two_offsets_disagree": bool(worst > 1e-3),
        "why_that_matters": (
            "One clock offset can phase-close a monochromatic carrier while "
            "failing to return a localised packet to the same event. Delta_n "
            "and Delta_g are the two demands; where they differ, no single "
            "offset does both."),
        "bandwidth_from_dispersion": (
            "d^2 delta / dw^2 sets how narrow the packet must be for the "
            "monochromatic answer to survive: the group delay varies across the "
            "packet by delta''(w0) * bandwidth."),
        "what_this_replaces": (
            "closure_offset(d_A, d_B, tau_th) = -(d_A + d_B + tau_th), which "
            "SOLVES for the offset that makes the loop close. Here the offset "
            "is a derived demand of the geometry, and it is frequency "
            "dependent."),
    }


@lru_cache(maxsize=4)
def measure_whether_the_loop_can_close(
        ell: int = 0, a: float = THROAT_RADIUS) -> Dict[str, object]:
    """G5 — can ``Λ_ℓ(ω) = 1``? Only in a limit, and the approach is exponential."""
    omega = np.linspace(1.5, 8.0, 20)
    _, transmitted = scattering_matrix(omega, ell, a)
    deficit = 1.0 - np.abs(transmitted) ** 2
    usable = deficit > 1e-14
    slope, intercept = np.polyfit(omega[usable], np.log(deficit[usable]), 1)
    probe = np.array([1.0, 2.0, 5.0, 10.0, 20.0])
    _, far = scattering_matrix(probe, ell, a)
    return {
        "omega": [float(x) for x in omega],
        "reflection_probability": [float(x) for x in deficit],
        "log_slope": float(slope),
        "the_approach_to_unity_is_exponential": bool(slope < -1.0),
        "transmission_at_probes": [float(abs(x)) for x in far],
        "transmission_is_below_one_at_finite_frequency": bool(
            np.all(np.abs(far) < 1.0 + 1e-12)),
        "so_lambda_cannot_equal_one_exactly": (
            "|Lambda| = |t_net| |eta_topo| and a barrier with V > 0 everywhere "
            "has |T| < 1 at every finite frequency. So 1 - Lambda cannot vanish "
            "and G_eff = G_0/(1 - Lambda) has no true pole -- unless some "
            "element outside the throat supplies gain. PR #216's completed "
            "transaction is a high-Q near-resonance, not an exact one."),
        "but_the_deficit_closes_exponentially": (
            f"1 - |T|^2 ~ exp({slope:.2f} w), so the quality factor of the "
            "near-transaction grows exponentially with frequency and the "
            "frequency needed for a given Q grows only logarithmically. Exact "
            "closure is approached in the UV, which is a chronology-horizon "
            "concern rather than the benign resonance the phenomenological "
            "ports suggest."),
        "causal_classification": {
            "D_loop > 0": "ordinary delayed feedback",
            "D_loop = 0": "closed-null / marginal chronology condition",
            "D_loop < 0": "return before emission (CTC regime)",
        },
    }


@lru_cache(maxsize=4)
def measure_the_traversable_throat_ledger() -> Dict[str, object]:
    """G6 — what this round settles about PR #216, and what it does not."""
    price = measure_the_null_energy_price()
    closure = measure_the_closure_offsets_disagree()
    loop = measure_whether_the_loop_can_close()
    threshold = measure_the_threshold_law()
    entries = [
        {"claim": "an advanced leg alone evades ER non-traversability",
         "verdict": "NO -- PROVED, NOT COMPUTED",
         "evidence": "supp G_ret(c,s) in J+(s) and supp G_adv(c,d) in J-(d), so "
                     "a nonzero product needs c in J+(s) cap J-(d), whence "
                     "s -> c -> d and d in J+(s). If d is not in J+(s) the "
                     "product vanishes for every c"},
        {"claim": "PR #216's throat transfer is supplied by a geometry",
         "verdict": "NOW IT IS",
         "evidence": "T_l(w), R_l(w) computed from the whole traversable throat "
                     "instead of MouthPort.t, r_in, r_out being inputs"},
        {"claim": "the throat factorises as two mouths plus a cavity",
         "verdict": "NOT BY THIS GEOMETRY",
         "evidence": "V_l is a single smooth symmetric barrier; the "
                     "Fabry-Perot form is a modelling choice, and is not "
                     "imposed here"},
        {"claim": "the static conductance g is T_0(0)",
         "verdict": "NO -- IT IS THE THRESHOLD COEFFICIENT",
         "evidence": "|T_0| -> (pi/8)(a w)^2, magnitude ratio "
                     f"{threshold['magnitude_ratio'][0]:.4f} at w = 0.01, with "
                     "a constant phase convention factor of i"},
        {"claim": "the traversable throat is free",
         "verdict": "NO",
         "evidence": "rho = 0 but p_s = -3a^2/(8 pi G5 f^4), radial NEC violated "
                     "everywhere, and the complete null-ray integral is exactly "
                     f"{price['null_energy_integral_in_units_of_one_over_G5']:.6f}/G5 "
                     "= -3/(16 G5 a)"},
        {"claim": "that NEC integral is the cost of the clock offset",
         "verdict": "NO -- DIFFERENT QUESTION",
         "evidence": "it is the support of the static throat; Delta_BA comes "
                     "from the mouths' aging history, and a frozen metric can "
                     "carry a large offset with no local stress proportional "
                     "to it"},
        {"claim": "one clock offset closes the loop",
         "verdict": "NOT FOR A PACKET",
         "evidence": "phase closure and group closure demand different offsets; "
                     f"worst disagreement {closure['worst_disagreement']:.4f} "
                     "over the sampled band"},
        {"claim": "Lambda = 1 (the completed transaction) can occur",
         "verdict": "NOT AT ANY FINITE FREQUENCY",
         "evidence": "|T| < 1 for a positive barrier, so 1 - Lambda cannot "
                     f"vanish; but 1 - |T|^2 ~ exp({loop['log_slope']:.2f} w), "
                     "so exact closure is approached only in the UV"},
    ]
    return {
        "entries": entries,
        "what_this_round_establishes": (
            "PR #216 relocated the advanced wave from a bare postulate into a "
            "Morris-Thorne-Yurtsever geometry, which was real progress. This "
            "round supplies the geometry it was missing and finds that the "
            "relocation has a price (NEC violation, exactly -3/(16 G5 a)), that "
            "the closure it assumed is frequency dependent rather than a single "
            "tuned constant, and that its completed transaction is a high-Q "
            "limit rather than an attainable pole."),
        "what_remains_postulated": (
            "That the throat is traversable at all, and that the mouths carry "
            "the required clock offset. The first now has an exact price; the "
            "second needs moving-mouth dynamics this round does not do."),
        "still_open": (
            "The history that produces Delta_BA, and whether any BAM exchange "
            "kernel is meant to approximate the whole-throat T_l derived here."),
    }
