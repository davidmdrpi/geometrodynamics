"""The retarded outer→inner transfer kernel on the benchmarked D = 5 background.

WHAT THIS IS
────────────
PR #270 deferred the retarded transfer function because it is a ratio of two
signals it could not trust. PR #274 settled which signal was right, against a
published spectrum. This round builds the object itself.

For the master equation ``∂_t²ψ − ∂_{r*}²ψ + V_ℓ ψ = 0``, a wave sent in from
the far region arrives at the horizon filtered. In frequency space that filter
is the transmission amplitude ``T_ℓ(ω)``; in time it is a convolution kernel

    ψ_transmitted(v) = (K_ℓ ⋆ ψ_incident)(v) ,        v = t + r*

and ``K_ℓ`` is the object. Because ``T_ℓ(ω) → 1`` as ``|ω| → ∞`` (high
frequencies fly over the barrier), the kernel splits as

    K_ℓ(t) = δ(t) + K_ℓ^reg(t)

— an instantaneous part plus **memory**.

THE ANSWER: THE TRANSFER IS NOT RIGID
─────────────────────────────────────
A rigid exchange kernel is ``δ(t)`` alone, possibly with a delay: whatever goes
in comes out, undistorted. The measured kernel is not that, and not marginally.

    ∫ K_ℓ^reg dt  =  −1      exactly
    ∫ |K_ℓ^reg| dt ≈ 2.02     against the δ's weight of 1

The first is a **sum rule**, not a fit: ``∫K dt = T(0)``, and ``T(0) = 0``
because the centrifugal barrier reflects zero-frequency waves completely. So
the memory term exactly cancels the instantaneous term at DC. A rigid kernel
transmits a static signal perfectly; the real geometry blocks it completely,
and does so *entirely* through the memory. The memory is not a correction to
rigid exchange — it is the same size as the thing it corrects.

WHY THE FREQUENCY DOMAIN IS THE RIGHT TOOL HERE, HAVING FAILED BEFORE
─────────────────────────────────────────────────────────────────────
PR #270 and #274 both failed to get quasinormal frequencies by shooting in real
``r``, because for ``Im ω < 0`` one solution grows like ``e^{|Im ω|R}`` and
swamps the coefficient being zeroed. **That objection does not apply here.**
For *real* ``ω`` both ``e^{±iωr*}`` have unit modulus and the problem is
perfectly conditioned. The scattering computation used here is not a repaired
version of the shoot that failed; it is a different, well-posed problem.

The evidence is unitarity, which is not imposed anywhere:

    | |R|² + |T|² − 1 |  ≲  8e-13

THREE EXACT ANCHORS, DERIVED NOT RECALLED
─────────────────────────────────────────
``D = 5`` supplies three closed forms that pin this calculation.

**1. The barrier peak at ℓ = 1 is rational.** With ``x = 1/r²`` the potential is
a cubic in ``x``; at ``ℓ = 1`` its maximum sits at ``x = 5/9`` exactly, so

    r² = 9/5 ,        V_max = 100/81         (exactly)

and ℓ = 1 is the *only* non-negative integer ℓ for which this is rational —
which is a **theorem**, not an inference from spot checks. With ``m = ℓ+1`` the
discriminant is ``(16m⁴ + 28m² + 73)/4``, and for ``m ≥ 4``

    (4m² + 3)²  <  16m⁴ + 28m² + 73  <  (4m² + 4)²

(left: ``4m² + 64 > 0``; right: ``4m² > 57``), so it is strictly between
consecutive squares and never square. ``m = 1, 2, 3`` give ``117``,
``441 = 21²``, ``1621`` — only ``m = 2``. ℓ = 1 is also the mode #270's two
codes disagreed on.

**2. The integrated potential is exact.** Since ``dr* = dr/f`` cancels the ``f``
in ``V``,

    ∫ V_ℓ dr*  =  ∫₁^∞ (L/r² + (9/4)/r⁴) dr  =  ℓ(ℓ+2) + 3/2

with ``L = (ℓ+1)² − ¼``. At ℓ = 1 this is ``9/2``.

**3. Hence the kernel's high-frequency phase is known in closed form.** The
leading eikonal phase through the barrier is ``−(1/2ω)∫V dr*``, so

    T_ℓ(ω) → exp(−i c_ℓ/ω) ,   c_ℓ = (ℓ(ℓ+2) + 3/2)/2 ,   c₁ = 9/4 exactly

This is what makes the kernel computable: ``T − 1 ~ −i c/ω`` decays too slowly
to transform numerically, and the exact ``c`` lets it be removed analytically
instead of windowed away. The **exact** value is what gets subtracted — a fitted
one would leave ``−i(c_exact − c_fit)/ω`` behind, still ``1/ω``, defeating the
purpose. The fit is kept as a *measurement against* the exact value, and agrees
to ``0.047 %`` uniformly in the outer edge.

THE OUTER BOUNDARY CONDITION IS NOT PLANE WAVES
───────────────────────────────────────────────
Matching to free ``e^{±iωr*}`` at a finite outer edge assumes ``ωr* ≫ ν``. At
the low-frequency end of a kernel grid that is badly false: for
``ω ≈ 0.005`` the outer turning scale is ``√L/ω ≈ 397``, and ``V`` at the edge
(``1.7e-4``) still exceeds ``ω²`` (``2.4e-5``) — the first bin sits *inside* the
centrifugal tail. That bin sets the DC end of the transform and hence the
numerical realisation of the ``−1`` sum rule, so it is the worst place to be
sloppy.

The fix is the exact solutions of the centrifugal tail, ``√x H^{(1,2)}_ν(x)``
with ``x = ωr*``, normalised to ``e^{±ix}`` so the high-frequency convention is
untouched. With them the low-``ω`` spectrum is **independent of the outer edge**
(relative spread ``5.6e-5`` across ``r*_out = 150, 300, 600``, where plane waves
drifted), and ``T(ω→0)`` fell from ``1.73e-06`` to ``4.10e-07`` — four times
closer to its exact value of zero.

THE THREE GATES, ALL MEASURED
─────────────────────────────
**Causal support.** ``K(t) = 0`` for ``t < 0``, to ``~1e-6`` away from the
front. This is not decoration: the exact zero for ``t < 0`` is a free read-out
of the **acausal** transform artefacts, and it caught two of them that would
otherwise have been read as physics — see below. It does *not* bound the total
numerical error, since a causal error can live entirely at ``t > 0``.

**Flux conservation.** ``|R|² + |T|² = 1`` to ``8e-13``, nowhere imposed.

**Late-time ringdown.** Fitting the kernel's own ringdown gives, against the
*external* published value ``1.01601691149 − 0.36232802385i``, real part within
``0.062 %–1.17 %`` and damping within ``0.11 %–3.80 %`` across nine
sample-spacing/window combinations. The band is the honest statement — as PR
#274 established for the time-domain solver, extraction choices move the fit, so
the spread is reported rather than the best row. Checked against the published
number, never against this repository's own solver.

**Independent cross-check.** The kernel is validated against the characteristic
null-grid evolution of PR #274 — a different propagation algorithm on the same
operator. Convolving ``K`` with the incident profile reproduces the transmitted
wave to **0.73 % peak, 0.13 % rms** at launch radius ``r* = 400``, with the
residual halving as the launch radius doubles.

WHAT THE CAUSALITY GATE CAUGHT
──────────────────────────────
Both were invisible without it, and both would have been mistaken for signal.

*A missing DC cell.* Sampling ``ω`` at right endpoints left the interval
``[0, dω]`` uncovered. Since ``T(0) = 0``, ``S(0) = −1``, and the omission put a
constant ``≈ −dω/π ≈ −1.9e-3`` under the whole kernel — including ``t < 0``,
where the answer is exactly zero. Midpoint sampling fixed it.

*Gibbs ringing from the ``1/ω`` tail.* Before the analytic subtraction, the
kernel sat on a ``1.6e-3`` plateau at large ``|t|`` — on both sides. Reading the
``t > 0`` plateau as a late-time tail would have been an error of interpretation
about numerics; the identical plateau at ``t < 0`` said it was an artefact.

SCOPE, AND ONE THING NOT DELIVERED
──────────────────────────────────
Fixed exact Tangherlini background, test scalar, no backreaction. One angular
channel at a time; no mode mixing (the background is spherically symmetric, so
there is none to have).

**The late-time power-law tail is NOT measured.** By ``t ≈ 40`` the ringdown has
decayed into the ``~1e-6`` noise floor established by the causality gate, and a
tail would be orders of magnitude below that. No exponent is quoted. The
frequency-domain route has excellent dynamic range near the barrier and poor
dynamic range in the tail, which is the opposite of what a time-domain
evolution offers; extracting the tail is a separate calculation with a separate
method, not a refinement of this one.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.special import hankel1, hankel2

from geometrodynamics.tangherlini.ringdown import (
    HORIZON,
    PUBLISHED_FUNDAMENTAL,
    PUBLISHED_SOURCE,
    characteristic_evolution,
    potential,
    radius_of_tortoise,
)

__all__ = [
    "integrated_potential",
    "high_frequency_phase_constant",
    "barrier_peak",
    "potential_samples",
    "outer_jost_solutions",
    "scattering_amplitudes",
    "transfer_spectrum",
    "transfer_kernel",
    "measure_the_exact_background_anchors",
    "measure_only_ell_one_has_a_rational_peak",
    "measure_the_high_frequency_tail_is_the_exact_one",
    "measure_the_low_frequency_outer_matching_is_converged",
    "measure_the_scattering_is_well_conditioned",
    "measure_the_kernel_is_causal",
    "measure_the_kernel_reproduces_the_published_ringdown",
    "measure_the_kernel_against_the_time_domain_evolution",
    "kernel_integrals",
    "measure_the_transfer_is_not_rigid",
    "measure_the_kernel_integrals_are_converged",
    "measure_the_subtraction_leaves_an_inverse_square_tail",
    "measure_the_late_time_tail_is_not_resolved",
    "measure_what_the_causality_gate_caught",
    "measure_the_transfer_kernel_ledger",
]

#: Default extent of the scattering domain in `r*`. Inside `-40` the potential
#: is below double precision (`V ~ e^{2 r*}`); outside, it falls only as
#: `L/r*^2`, which is why the outer edge is the one that has to be argued for.
INNER_EDGE = -50.0
OUTER_EDGE = 150.0


# ── exact closed forms ──────────────────────────────────────────────────────

def integrated_potential(ell: int) -> float:
    """``∫V_ℓ dr* = ℓ(ℓ+2) + 3/2``, exactly.

    ``dr* = dr/f`` cancels the overall ``f`` in ``V``, leaving an elementary
    integral of ``L/r² + (9/4)/r⁴`` from the horizon to infinity.
    """
    return float(ell * (ell + 2)) + 1.5


def high_frequency_phase_constant(ell: int) -> float:
    """``c_ℓ = ½∫V dr*``, so that ``T_ℓ(ω) → exp(−i c_ℓ/ω)``. ``c₁ = 9/4``."""
    return 0.5 * integrated_potential(ell)


def barrier_peak(ell: int) -> Dict[str, float]:
    """Peak of ``V_ℓ``. At ``ℓ = 1`` it is ``100/81`` at ``r² = 9/5``, exactly.

    In ``x = 1/r²`` the potential is the cubic
    ``V = Lx + (9/4 − L)x² − (9/4)x³``, so the peak solves a quadratic. At
    ``ℓ = 1`` the discriminant is a perfect square and everything is rational;
    at other ``ℓ`` it is not.
    """
    limit = (ell + 1) ** 2 - 0.25
    # 27/4 x^2 - 2(9/4 - L) x - L = 0
    a, b, c = 6.75, -2.0 * (2.25 - limit), -limit
    x = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
    r_squared = 1.0 / x
    return {
        "ell": ell,
        "x_peak": x,
        "r_squared": r_squared,
        "r": math.sqrt(r_squared),
        "r_star": float(_tortoise(math.sqrt(r_squared))),
        "V_max": float(potential(math.sqrt(r_squared), ell)),
    }


def _tortoise(r: float) -> float:
    return r + 0.5 * math.log((r - HORIZON) / (r + HORIZON))


# ── the scattering problem, well conditioned for real omega ─────────────────

@lru_cache(maxsize=8)
def potential_samples(ell: int, inner: float = INNER_EDGE,
                      outer: float = OUTER_EDGE,
                      steps: int = 40000) -> Tuple[np.ndarray, float]:
    """``V_ℓ`` at the midpoint of each step of a uniform ``r*`` grid."""
    edges = np.linspace(inner, outer, steps + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    values = np.array([0.0 if x < -40.0
                       else float(potential(radius_of_tortoise(x), ell))
                       for x in mid])
    return values, (outer - inner) / steps


def outer_jost_solutions(ell: int, x):
    """Exact outer solutions of the centrifugal tail, normalised to ``e^{±ix}``.

    Matching to *free* plane waves at a finite outer edge is wrong at low
    frequency, and this is the dominant error there rather than a refinement.
    Asymptotically the ``r*`` equation is

        ψ'' + (ω² − L/r*²) ψ = 0 ,      L = ν² − ¼ ,   ν = ℓ + 1

    whose **exact** solutions are ``√x H^{(1,2)}_ν(x)`` with ``x = ω r*``. The
    plane-wave form is only their ``x → ∞`` limit, so a plane-wave match
    silently assumes ``ω r* ≫ ν``. At the lowest frequency of a kernel grid that
    is badly false: for ``ω ≈ 0.005`` the outer turning scale is ``√L/ω ≈ 400``,
    far beyond any affordable edge, and ``V(r*_out) ≈ 1.7e-4`` still exceeds
    ``ω² ≈ 2.4e-5``. That first bin is what fixes the DC end of the transform and
    hence the numerical realisation of the ``−1`` sum rule.

    Normalisation carries the phase ``e^{±i(νπ/2 + π/4)}`` so these reduce to
    ``e^{±ix}`` exactly, leaving the high-frequency convention untouched — only
    the low-frequency matching moves. Their Wronskian is exactly ``−2i``, which
    the tests check rather than assume.
    """
    order = ell + 1
    x = np.asarray(x, dtype=float)
    prefactor = np.sqrt(np.pi * x / 2.0)
    d_prefactor = np.sqrt(np.pi / 2.0) / (2.0 * np.sqrt(x))
    phase = np.exp(1j * (order * np.pi / 2.0 + np.pi / 4.0))
    first, second = hankel1(order, x), hankel2(order, x)
    d_first = 0.5 * (hankel1(order - 1, x) - hankel1(order + 1, x))
    d_second = 0.5 * (hankel2(order - 1, x) - hankel2(order + 1, x))
    plus = prefactor * first * phase
    minus = prefactor * second * np.conj(phase)
    d_plus = (d_prefactor * first + prefactor * d_first) * phase
    d_minus = (d_prefactor * second + prefactor * d_second) * np.conj(phase)
    return plus, minus, d_plus, d_minus


def scattering_amplitudes(omega, ell: int, inner: float = INNER_EDGE,
                          outer: float = OUTER_EDGE, steps: int = 40000
                          ) -> Tuple[np.ndarray, np.ndarray]:
    """``R(ω)``, ``T(ω)`` by a piecewise-constant transfer matrix.

    Boundary conditions: purely ingoing at the horizon, ``e^{−iωr*} + R e^{iωr*}``
    at infinity, so ``T`` is the outer→inner transmission amplitude.

    On each step ``V ≈ V_mid`` and the propagator is exact:

        ψ(x+h)  =  cos(kh) ψ + sin(kh)/k ψ' ,   k = √(ω² − V_mid)
        ψ'(x+h) = −k sin(kh) ψ + cos(kh) ψ'

    with complex ``k`` covering the classically forbidden region.

    **Why this succeeds where the quasinormal shoots failed.** For real ``ω``
    the *asymptotic* in and out waves are unit-flux, so the coefficient being
    extracted is not swamped by an exponentially larger companion — unlike
    ``Im ω < 0``, where the outgoing solution grows like ``e^{|Im ω|R}`` over
    the whole matching range. That is a claim about the asymptotic
    normalisation, **not** about the propagation everywhere: under the barrier
    ``k = √(ω² − V)`` is imaginary and the local propagator is hyperbolic, so
    one solution does grow relative to the other across the forbidden region.
    The method is therefore well conditioned *on the tested real-frequency
    range* — evidenced by unitarity and step refinement — rather than
    structurally immune to conditioning problems as ``ω → 0`` or ``ℓ`` grows.

    Vectorised over ``ω``: the spatial loop is shared across all frequencies.
    """
    values, h = potential_samples(ell, inner, outer, steps)
    w = np.asarray(omega, dtype=complex)
    psi = np.exp(-1j * w * inner)
    dpsi = -1j * w * psi
    for v in values:
        k = np.sqrt(w * w - v)
        kh = k * h
        cos_kh = np.cos(kh)
        sin_over_k = np.sin(kh) / k
        psi, dpsi = (cos_kh * psi + sin_over_k * dpsi,
                     -(k * k) * sin_over_k * psi + cos_kh * dpsi)
    x = np.real(w) * outer
    plus, minus, d_plus, d_minus = outer_jost_solutions(ell, x)
    psi_x = dpsi / np.real(w)                  # d/dx = (1/ω) d/dr*
    incident = (psi * d_plus - psi_x * plus) / 2j
    reflected = (psi_x * minus - psi * d_minus) / 2j
    return reflected / incident, 1.0 / incident


@lru_cache(maxsize=8)
def transfer_spectrum(ell: int, omega_max: float = 40.0, count: int = 4096
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """``T_ℓ(ω)`` on a **midpoint** grid covering ``[0, ω_max]``.

    Midpoint sampling is not a detail. Right-endpoint sampling leaves ``[0, dω]``
    uncovered, and since ``T(0) = 0`` that omission puts a constant ``≈ −dω/π``
    under the entire kernel — including ``t < 0``, where the answer is exactly
    zero. The causality gate is what exposed it.
    """
    spacing = omega_max / count
    omega = (np.arange(count) + 0.5) * spacing
    reflected, transmitted = scattering_amplitudes(omega, ell)
    return omega, reflected, transmitted, spacing


def transfer_kernel(times, ell: int = 1, omega_max: float = 40.0,
                    count: int = 4096, decay: float = 1.0) -> np.ndarray:
    """``K_ℓ^reg(t)``, the kernel with its ``δ(t)`` removed.

    ``T → 1`` at large ``ω``, so ``K = δ(t) + K^reg``. What is left,
    ``S = T − 1``, decays only as ``−i c_ℓ/ω`` — too slowly to transform
    numerically. Since ``c_ℓ`` is known exactly, subtract a causal function with
    the same tail and add its transform back in closed form:

        A(ω) = −i c /(ω + i a)        →        −c e^{−at} θ(t)

    whose only pole is at ``ω = −ia``, in the lower half plane, so the
    subtraction cannot itself introduce an acausal piece.

    The coefficient used is the **exact** ``c_ℓ = ½∫V dr*``, not a fit to the
    computed spectrum. That distinction is not cosmetic: subtracting a fitted
    ``c_fit ≠ c_exact`` leaves ``−i(c_exact − c_fit)/ω`` in the remainder, which
    still decays only as ``1/ω``, so the whole purpose of the subtraction — a
    remainder falling like ``1/ω²`` — would be silently forfeited. The fitted
    value is retained, but as a *measurement against* the exact one; see
    ``measure_the_high_frequency_tail_is_the_exact_one``.
    """
    omega, _, transmitted, spacing = transfer_spectrum(ell, omega_max, count)
    residual = transmitted - 1.0
    exact = high_frequency_phase_constant(ell)
    analytic = -1j * exact / (omega + 1j * decay)
    remainder = residual - analytic

    stamps = np.atleast_1d(np.asarray(times, dtype=float))
    out = np.empty(stamps.size)
    for start in range(0, stamps.size, 4000):        # bound the outer product
        block = stamps[start:start + 4000]
        out[start:start + 4000] = 2.0 * np.real(
            np.exp(-1j * np.outer(block, omega)) @ remainder)
    out *= spacing / (2.0 * np.pi)
    out += np.where(stamps > 0.0, -exact * np.exp(-decay * stamps), 0.0)
    return out


# ── measurements ────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def measure_the_exact_background_anchors() -> Dict[str, object]:
    """R0 — three closed forms this round derives and then checks numerically."""
    peaks = []
    for ell in (0, 1, 2, 3):
        row = barrier_peak(ell)
        peaks.append(row)
    one = barrier_peak(1)

    integrals = []
    for ell in (0, 1, 2, 3):
        values, h = potential_samples(ell)
        integrals.append({
            "ell": ell,
            "exact": integrated_potential(ell),
            "quadrature_on_the_truncated_domain": float(np.sum(values) * h),
            "missing_tail_beyond_the_outer_edge": (
                ((ell + 1) ** 2 - 0.25) / OUTER_EDGE),
        })

    return {
        "barrier_peaks": peaks,
        "ell_1_peak_is_exactly_100_over_81": bool(
            abs(one["V_max"] - 100.0 / 81.0) < 1e-12),
        "ell_1_peak_radius_squared_is_exactly_9_over_5": bool(
            abs(one["r_squared"] - 1.8) < 1e-12),
        "only_ell_1_is_rational": (
            "In x = 1/r^2 the peak solves 27/4 x^2 - 2(9/4 - L) x - L = 0. At "
            "l = 1 the discriminant is a perfect square (110.25 = 10.5^2) and "
            "x = 5/9 exactly; l = 0, 2, 3 give sqrt(13), sqrt(1621), sqrt(57). "
            "l = 1 is also the mode PR #270's two codes disagreed on."),
        "integrated_potential": integrals,
        "the_integral_is_elementary": (
            "dr* = dr/f cancels the f in V, leaving int_1^inf (L/r^2 + "
            "(9/4)/r^4) dr = L + 3/4 = l(l+2) + 3/2."),
        "high_frequency_phase_constants": {
            str(ell): high_frequency_phase_constant(ell) for ell in (0, 1, 2, 3)},
        "c_of_ell_1_is_exactly_nine_quarters": bool(
            abs(high_frequency_phase_constant(1) - 2.25) < 1e-15),
    }


def measure_only_ell_one_has_a_rational_peak(
        limit: int = 400) -> Dict[str, object]:
    """R0b — a **theorem**, not three spot checks.

    With ``m = ℓ + 1`` the peak condition's discriminant is
    ``(16m⁴ + 28m² + 73)/4``, so the peak is rational exactly when
    ``16m⁴ + 28m² + 73`` is a perfect square. For ``m ≥ 4``,

        (4m² + 3)²  <  16m⁴ + 28m² + 73  <  (4m² + 4)²

    — the left inequality is ``4m² + 64 > 0`` and the right is ``4m² > 57``,
    i.e. ``m ≥ 4``. Strictly between consecutive squares, so never square.
    That leaves ``m = 1, 2, 3``, giving ``117``, ``441 = 21²``, ``1621``. Only
    ``m = 2``, i.e. ``ℓ = 1``, survives.

    An earlier draft asserted this from checking ``ℓ = 0, 2, 3``. It is provable,
    and the bracketing step is verified here over a range rather than trusted.
    """
    def form(m: int) -> int:
        return 16 * m ** 4 + 28 * m ** 2 + 73

    small = {m: form(m) for m in (1, 2, 3)}
    squares = {m: (v, int(round(math.isqrt(v))) ** 2 == v)
               for m, v in small.items()}
    bracket_holds = all(
        (4 * m * m + 3) ** 2 < form(m) < (4 * m * m + 4) ** 2
        for m in range(4, limit + 1))
    # and no accidental square anywhere in range, as a belt-and-braces sweep
    exhaustive = [m for m in range(1, limit + 1)
                  if math.isqrt(form(m)) ** 2 == form(m)]
    return {
        "discriminant_numerator": "16 m^4 + 28 m^2 + 73,  m = l + 1",
        "small_cases": {str(m): {"value": v, "is_square": squares[m][1]}
                        for m, v in small.items()},
        "bracketing_holds_for_m_at_least_four": bool(bracket_holds),
        "bracketing_argument": (
            "(4m^2+3)^2 < 16m^4+28m^2+73 < (4m^2+4)^2 for m >= 4: the left is "
            "4m^2 + 64 > 0, the right is 4m^2 > 57. Strictly between "
            "consecutive squares, hence never a square."),
        "squares_found_in_range": exhaustive,
        "only_m_equals_two": bool(exhaustive == [2]),
        "therefore": (
            "l = 1 is the ONLY non-negative integer l with a rational barrier "
            "peak. This is a theorem, not an inference from spot checks."),
    }


@lru_cache(maxsize=4)
def measure_the_high_frequency_tail_is_the_exact_one(
        ell: int = 1) -> Dict[str, object]:
    """R1b — the fitted ``c`` as a *measurement against* the exact one.

    ``transfer_kernel`` subtracts the exact ``c_ℓ = ½∫V dr*``. Subtracting a
    fitted coefficient instead would leave ``−i(c_exact − c_fit)/ω`` in the
    remainder — still ``1/ω``, forfeiting the point of the subtraction.

    **The estimator matters, and an earlier draft used a biased one.**
    ``Re[iω(T−1)]`` is not ``c`` at finite ``ω`` even for the exact asymptote
    ``T = e^{−ic/ω}``:

        Re[iω(e^{−ic/ω} − 1)] = ω sin(c/ω) = c − c³/(6ω²) + O(ω⁻⁴)

    At ``c = 2.25`` over the sampled band that deterministic bias is ``−1.3e-3``
    — the *entire* size of the "deficit" the earlier draft reported, and which it
    attributed to the ``1/r⁴`` and ``1/r⁶`` parts of ``V``. That attribution was
    wrong. The unbiased estimator ``c(ω) = −ω arg T(ω)`` is exact for the toy
    asymptote, and against the real spectrum it gives a residual of ``+2.6e-4``
    — four times smaller and of the **opposite sign**.

    That the biased estimator's deficit is *itself* independent of the outer edge
    is the tell: a truncation effect would move with the edge, a fixed-form bias
    would not.
    """
    exact = high_frequency_phase_constant(ell)
    rows = []
    for outer in (150.0, 300.0, 600.0):
        steps = int((outer + 50.0) / 0.005)
        omega = (np.arange(1024) + 0.5) * (40.0 / 1024)
        _, transmitted = scattering_amplitudes(
            omega, ell, outer=outer, steps=steps)
        phase = float(np.mean(-omega[-100:] * np.angle(transmitted[-100:])))
        biased = float(np.mean((1j * omega * (transmitted - 1.0)).real[-100:]))
        rows.append({"outer_edge": outer,
                     "phase_estimator": phase,
                     "biased_estimator": biased,
                     "exact": exact,
                     "phase_deviation": phase - exact,
                     "biased_deviation": biased - exact,
                     "relative_error": abs(phase - exact) / exact})
    spread = (max(r["phase_estimator"] for r in rows)
              - min(r["phase_estimator"] for r in rows))
    biased_spread = (max(r["biased_estimator"] for r in rows)
                     - min(r["biased_estimator"] for r in rows))
    predicted_bias = -exact ** 3 / (6.0 * 38.0 ** 2)
    return {
        "exact": exact,
        "rows": rows,
        "spread_across_outer_edges": spread,
        "independent_of_the_outer_edge": bool(spread < 1e-6),
        "worst_relative_error": max(r["relative_error"] for r in rows),
        "agrees_with_the_exact_value": bool(
            max(r["relative_error"] for r in rows) < 1e-3),
        "predicted_bias_of_the_naive_estimator": predicted_bias,
        "observed_bias_of_the_naive_estimator": float(
            np.mean([r["biased_deviation"] for r in rows])),
        "the_biased_estimator_is_also_edge_independent": bool(
            biased_spread < 1e-6),
        "why_the_exact_value_is_the_one_subtracted": (
            "If T - 1 = -i c_exact/w + O(w^-2) and a fitted c_fit != c_exact is "
            "subtracted, the remainder retains -i(c_exact - c_fit)/w and still "
            "falls only as 1/w. The subtraction exists to make the remainder "
            "O(1/w^2); using a fit would quietly defeat it."),
        "a_correction_to_an_earlier_draft": (
            "That draft used Re[i w (T-1)], which carries a deterministic "
            "-c^3/(6 w^2) bias, and read the resulting 1.06e-3 shortfall as the "
            "1/r^4 and 1/r^6 part of V that a centrifugal Jost condition cannot "
            "capture. It was the estimator. With -w arg T the residual is "
            "+2.6e-4, four times smaller and opposite in sign, and no physical "
            "attribution for it is offered here."),
    }


@lru_cache(maxsize=4)
def measure_the_subtraction_leaves_an_inverse_square_tail(
        ell: int = 1) -> Dict[str, object]:
    """R1d — the property the inverse transform actually needs.

    Subtracting ``A(ω) = −i c/(ω + ia)`` is only useful if the remainder falls
    like ``1/ω²``. Checking the coefficient ``c`` is necessary but not
    sufficient, so check the remainder directly: ``ω²|S − A|`` must stay
    bounded rather than growing.
    """
    omega, _, transmitted, _ = transfer_spectrum(ell)
    residual = transmitted - 1.0
    exact = high_frequency_phase_constant(ell)
    rows = []
    for decay in (0.5, 1.0, 2.0):
        remainder = residual - (-1j * exact / (omega + 1j * decay))
        picks = [len(omega) // 8, len(omega) // 4, len(omega) // 2, -1]
        scaled = [float(omega[i] ** 2 * abs(remainder[i])) for i in picks]
        rows.append({"decay": decay,
                     "omega": [float(omega[i]) for i in picks],
                     "omega_squared_times_remainder": scaled,
                     "bounded": bool(max(scaled) < 2.0 * min(scaled))})
    return {
        "rows": rows,
        "the_remainder_falls_like_one_over_omega_squared": bool(
            all(r["bounded"] for r in rows)),
        "the_plateau_depends_on_the_subtraction_parameter": (
            "w^2|S - A| settles to a different constant for each a, because a "
            "different A redistributes the O(1/w^2) tail between the analytic "
            "and numerical pieces. The KERNEL must not depend on a at all -- "
            "that is checked separately in the convergence study."),
    }


@lru_cache(maxsize=4)
def measure_the_low_frequency_outer_matching_is_converged(
        ell: int = 1) -> Dict[str, object]:
    """R1c — the low-``ω`` end, where plane-wave matching is simply wrong.

    At ``ω ≈ 0.005`` the outer turning scale is ``√L/ω ≈ 400`` and
    ``V(r*=150) ≈ 1.7e-4`` still exceeds ``ω² ≈ 2.4e-5``: the first bin of the
    kernel grid is *inside* the centrifugal tail. Matching it to free plane
    waves there assumes the opposite. Since that bin sets the DC end of the
    transform, it is the bin the ``−1`` sum rule depends on most.

    With the Jost condition the low-``ω`` spectrum is independent of the outer
    edge; with plane waves it drifts. That is the check, run at the frequency
    where it matters rather than at a convenient one.
    """
    omega = np.array([40.0 / (2 * 4096), 0.02, 0.1])
    turning = [math.sqrt((ell + 1) ** 2 - 0.25) / w for w in omega]
    v_edge = float(potential(radius_of_tortoise(OUTER_EDGE), ell))
    rows = []
    for outer in (150.0, 300.0, 600.0):
        steps = int((outer + 50.0) / 0.005)
        _, transmitted = scattering_amplitudes(
            omega, ell, outer=outer, steps=steps)
        rows.append({"outer_edge": outer,
                     "transmission_modulus": [float(abs(t)) for t in transmitted]})
    spreads = [max(r["transmission_modulus"][i] for r in rows)
               - min(r["transmission_modulus"][i] for r in rows)
               for i in range(len(omega))]
    relative = [s / max(r["transmission_modulus"][i] for r in rows)
                for i, s in enumerate(spreads)]
    return {
        "omega": [float(w) for w in omega],
        "outer_turning_scale": turning,
        "potential_at_the_outer_edge": v_edge,
        "omega_squared_at_the_lowest_bin": float(omega[0] ** 2),
        "the_lowest_bin_is_inside_the_centrifugal_tail": bool(
            v_edge > omega[0] ** 2),
        "rows": rows,
        "relative_spread_across_outer_edges": relative,
        "converged_in_the_outer_edge": bool(max(relative) < 1e-3),
        "why_this_bin_matters": (
            "It sets the DC end of the inverse transform, and therefore the "
            "numerical realisation of the int K_reg dt = -1 sum rule."),
        "what_changed": (
            "Replacing plane-wave matching with the exact centrifugal Jost "
            "solutions moved T(w -> 0) from 1.73e-06 to 4.10e-07 -- a factor of "
            "four closer to its exact value of zero -- and made the low-w "
            "spectrum independent of the outer edge instead of drifting with "
            "it."),
    }


@lru_cache(maxsize=4)
def measure_the_scattering_is_well_conditioned(
        ell: int = 1,
        probes: Sequence[float] = (0.1, 0.3, 0.7, 1.0, 1.5, 3.0, 10.0, 30.0)
        ) -> Dict[str, object]:
    """R1 — unitarity, nowhere imposed, and why this is not the failed shoot."""
    omega = np.array(probes, dtype=float)
    reflected, transmitted = scattering_amplitudes(omega, ell)
    residual = np.abs(reflected) ** 2 + np.abs(transmitted) ** 2 - 1.0

    refinement = []
    for steps in (10000, 20000, 40000):
        r_s, t_s = scattering_amplitudes(np.array([1.0]), ell, steps=steps)
        refinement.append({"steps": steps,
                           "transmission_modulus": float(abs(t_s[0])),
                           "unitarity_residual": float(
                               abs(r_s[0]) ** 2 + abs(t_s[0]) ** 2 - 1.0)})
    differences = [abs(refinement[i + 1]["transmission_modulus"]
                       - refinement[i]["transmission_modulus"])
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
        "transmission_goes_to_zero_at_dc": bool(abs(transmitted[0]) < 1e-2),
        "transmission_goes_to_one_at_high_frequency": bool(
            abs(abs(transmitted[-1]) - 1.0) < 1e-6),
        "why_this_is_not_the_failed_shoot": (
            "PR #270 and #274 could not converge a quasinormal frequency by "
            "shooting in real r, because for Im w < 0 the outgoing solution "
            "grows like e^{|Im w| R} across the whole matching range and swamps "
            "the coefficient being zeroed. Here w is REAL, so the ASYMPTOTIC in "
            "and out waves are unit-flux and the extracted coefficient is not "
            "swamped. This is a different, well-posed problem -- not a repair "
            "of the one that failed."),
        "the_conditioning_claim_is_scoped": (
            "Unit modulus is a property of the asymptotic normalisation, not of "
            "the propagation everywhere: under the barrier k = sqrt(w^2 - V) is "
            "imaginary and the propagator is hyperbolic, so one solution does "
            "grow relative to the other. The correct statement is WELL "
            "CONDITIONED ON THE TESTED REAL-FREQUENCY RANGE -- evidenced by "
            "unitarity to ~1e-13 and second-order step refinement -- not "
            "structurally immune as w -> 0 or l grows."),
    }


@lru_cache(maxsize=4)
def measure_the_kernel_is_causal(ell: int = 1) -> Dict[str, object]:
    """R2 — ``K(t < 0) = 0``, and the residual as a measured noise floor."""
    early = np.array([-300.0, -200.0, -100.0, -50.0, -20.0, -10.0, -5.0,
                      -2.0, -1.0, -0.5])
    before = transfer_kernel(early, ell)
    far = np.array([-300.0, -200.0, -100.0, -50.0])
    floor = float(np.max(np.abs(transfer_kernel(far, ell))))

    late = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, 20.0,
                     30.0, 40.0, 60.0, 100.0])
    after = transfer_kernel(late, ell)

    return {
        "times_before_zero": [float(x) for x in early],
        "kernel_before_zero": [float(x) for x in before],
        "worst_acausal_value": float(np.max(np.abs(before))),
        "noise_floor_far_from_the_front": floor,
        "the_kernel_is_causal": bool(np.max(np.abs(before)) < 1e-3),
        "times_after_zero": [float(x) for x in late],
        "kernel_after_zero": [float(x) for x in after],
        "the_exact_zero_is_an_acausal_artifact_monitor": (
            "K(t) vanishes identically for t < 0, so whatever the computation "
            "returns there is a direct read-out of its ACAUSAL transform "
            "artifacts -- no reference value needed. That is what caught the "
            "missing DC cell and the symmetric truncation ringing, both of "
            "which contaminate t < 0 and t > 0 alike."),
        "it_is_not_a_total_error_bar": (
            "A numerical error can be perfectly causal and live only at t > 0: "
            "outer-boundary error, finite w_max, the arbitrary subtraction "
            "parameter, quadrature bias. This round demonstrated exactly that "
            "-- replacing plane-wave outer matching with the Jost condition "
            "moved positive-time results with no matching negative-time "
            "signature. So the negative-time residual bounds one FAMILY of "
            "error, not the total, and the late-time tail is judged out of "
            "reach on positive-time parameter variation as well."),
    }


@lru_cache(maxsize=4)
def measure_the_kernel_reproduces_the_published_ringdown(
        ell: int = 1) -> Dict[str, object]:
    """R3 — the kernel's own ringdown, against the EXTERNAL published value."""
    published = PUBLISHED_FUNDAMENTAL[ell]
    rows = []
    for spacing in (0.2, 0.3, 0.5):
        for window in ((3.0, 14.0), (4.0, 16.0), (5.0, 18.0)):
            stamps = np.arange(window[0], window[1], spacing)
            values = transfer_kernel(stamps, ell)
            found = _prony(stamps, values, 2)
            candidates = [w for w in found
                          if 0.5 < w.real < 2.0 and -1.0 < w.imag < -0.05]
            if not candidates:
                continue
            best = min(candidates, key=lambda w: abs(w.imag))
            rows.append({
                "sample_spacing": spacing,
                "window": list(window),
                "omega": [best.real, best.imag],
                "real_relative_error": abs(best.real - published.real) / abs(published.real),
                "damping_relative_error": abs(best.imag - published.imag) / abs(published.imag),
            })
    real_errors = [r["real_relative_error"] for r in rows]
    damping_errors = [r["damping_relative_error"] for r in rows]
    return {
        "published": [published.real, published.imag],
        "source": PUBLISHED_SOURCE,
        "rows": rows,
        "real_part_band": {"min": min(real_errors), "max": max(real_errors)},
        "damping_band": {"min": min(damping_errors), "max": max(damping_errors)},
        # Stated at the level the measurement supports: every extraction lands
        # on the fundamental mode, and the spread between them is the honest
        # uncertainty -- not a claim that every row is sub-percent.
        "every_fit_lands_on_the_fundamental_mode": bool(
            max(real_errors) < 2e-2),
        "the_best_fit_is_sub_tenth_percent": bool(min(real_errors) < 1e-3),
        "the_ringdown_is_the_published_one": bool(
            max(real_errors) < 2e-2 and min(real_errors) < 1e-3),
        "scored_against_an_external_value": (
            "The comparison is with the published continued-fraction value, "
            "never with this repository's own time-domain solver. Scoring a "
            "kernel against a frequency extracted from the same machinery that "
            "produced it would not be a check."),
        "the_band_is_reported_not_the_best_row": (
            "Extraction choices move the fit, exactly as PR #274 measured for "
            "the time-domain solver. The spread is the honest statement."),
    }


def _prony(times: np.ndarray, values: np.ndarray, modes: int) -> np.ndarray:
    spacing = float(times[1] - times[0])
    width = 2 * modes
    design = np.array([values[i:i + width] for i in range(len(values) - width)])
    target = -values[width:]
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    roots = np.roots(np.concatenate([[1.0], coefficients[::-1]]))
    roots = roots[np.abs(roots) > 1e-12]
    return 1j * np.log(roots) / spacing


@lru_cache(maxsize=4)
def measure_the_kernel_against_the_time_domain_evolution(
        ell: int = 1, width: float = 3.0) -> Dict[str, object]:
    """R4 — the independent check: ``K ⋆ g`` against a characteristic evolution.

    Deep inside, the transmitted wave as a function of ``v = t + r*`` is exactly
    the convolution of the incident profile with the kernel. PR #274's
    time-domain evolution uses a completely different propagation algorithm --
    a characteristic null-grid march rather than a frequency-domain transfer
    matrix — so agreement is real cross-validation. It is *not* independent in
    every sense: both rest on the same ``potential`` and tortoise definitions,
    and this module imports the characteristic solver. The accurate claim is
    **independent numerical propagation methods on the same operator**, which
    tests the propagation and not the operator.

    **The launch radius is the whole subtlety, and it is now quantified.** The
    incident amplitude is only defined where the wave is free. On the ``u = 0``
    line the pulse sits at ``r* = v_c/2``, and the phase it has *not yet*
    accumulated is set by the potential remaining beyond that point, which is
    known exactly:

        ∫_{r*_launch}^∞ V dr*  ≈  L / r*_launch

    The residual disagreement tracks that quantity — halving as the launch
    radius doubles — which is what identifies it as placement rather than a
    disagreement between methods.

    *An earlier draft reported ``0.92 %`` at ``r* = 100`` while matching the
    outer boundary to plane waves. That number was flattered by cancellation:
    the plane-wave outer condition carried an error of its own, in the opposite
    direction. With the Jost condition both errors are attributed correctly, the
    same launch gives ``2.73 %``, and the series converges cleanly.*
    """
    omega, _, transmitted, spacing = transfer_spectrum(ell)
    rows = []
    for centre, step, horizon in ((200.0, 0.05, 500.0),
                                  (400.0, 0.10, 900.0),
                                  (800.0, 0.10, 1500.0)):
        launch = 0.5 * centre
        times, signal = characteristic_evolution(
            ell, step=step, t_max=horizon, observer_rstar=-30.0,
            pulse_centre=centre, pulse_width=width)
        v = times - 30.0
        mask = (v > centre - 40.0) & (v < centre + 50.0)
        if mask.sum() < 10:
            continue
        profile = (width * math.sqrt(math.pi)
                   * np.exp(-(omega * width / 2.0) ** 2)
                   * np.exp(1j * omega * centre))
        predicted = (2.0 * np.real(
            np.exp(-1j * np.outer(v[mask], omega)) @ (transmitted * profile))
            ) * spacing / (2.0 * np.pi)
        measured = signal[mask]
        peak = float(np.max(np.abs(predicted)))
        rows.append({
            "pulse_centre": centre,
            "launch_r_star": launch,
            "potential_beyond_the_launch_point": (
                ((ell + 1) ** 2 - 0.25) / launch),
            "peak_relative_max_difference": float(
                np.max(np.abs(measured - predicted)) / peak),
            "peak_relative_rms_difference": float(
                math.sqrt(float(np.mean((measured - predicted) ** 2))) / peak),
        })
    best = min(rows, key=lambda r: r["peak_relative_max_difference"])
    ratios = [a["peak_relative_max_difference"] / b["peak_relative_max_difference"]
              for a, b in zip(rows[:-1], rows[1:])]
    return {
        "rows": rows,
        "best": best,
        "successive_ratios": ratios,
        "the_two_methods_agree": bool(
            best["peak_relative_max_difference"] < 0.01),
        "the_residual_halves_as_the_launch_radius_doubles": bool(
            all(1.7 < r < 2.3 for r in ratios)),
        "what_this_exposed": (
            "The incident amplitude is only defined where the wave is free, and "
            "the residual tracks the exactly-known potential remaining beyond "
            "the launch point. PR #274 launched at r* = 6, where V ~ 0.1 -- "
            "harmless for a quasinormal frequency, since a ringdown does not "
            "care how it was excited, and fatal for a transmission ratio."),
        "an_earlier_number_was_flattered_by_cancellation": (
            "With plane-wave outer matching this check read 0.92% at r* = 100. "
            "That was two errors partly cancelling: the plane-wave outer "
            "condition carried its own error in the opposite direction. Under "
            "the correct Jost condition the same launch reads 2.73%, and the "
            "series converges as 1/r*_launch. The larger number is the honest "
            "one."),
    }


def kernel_integrals(ell: int = 1, spacing: float = 0.005,
                     horizon: float = 300.0, **kwargs) -> Tuple[float, float]:
    """``∫K_reg dt`` and ``∫|K_reg| dt`` on a **midpoint** time grid.

    Midpoint sampling matters here for the same reason it did on the frequency
    grid. ``K_reg`` jumps at ``t = 0⁺`` — the analytic piece contributes
    ``−c`` there — so a grid starting at ``t₀ > 0`` silently omits
    ``∫₀^{t₀} ≈ −c t₀``. An earlier draft integrated from ``t = 0.001`` and
    reported ``−0.997757`` for a quantity whose exact value is ``−1``; that gap
    was almost entirely the missing first cell, and it shrank as ``O(dt)``
    rather than converging. A grid of ``(j + ½)dt`` covering ``[0, t_max]``
    exactly gives ``−0.999996``.
    """
    stamps = (np.arange(int(horizon / spacing)) + 0.5) * spacing
    values = transfer_kernel(stamps, ell, **kwargs)
    return (float(np.sum(values) * spacing),
            float(np.sum(np.abs(values)) * spacing))


@lru_cache(maxsize=4)
def measure_the_transfer_is_not_rigid(ell: int = 1) -> Dict[str, object]:
    """R5 — **the deliverable**: how far the kernel is from instantaneous.

    The durable statement is not the memory mass but the pair of limits. With
    the chosen asymptotic normalisation ``T(∞) = 1`` while ``T(0) = 0``;
    therefore the transfer function cannot be a rigid phase or delay, and the
    regular part must carry signed weight ``−1``. That is analytic. Everything
    below is a numerical check of it.
    """
    omega, _, transmitted, _ = transfer_spectrum(ell)
    signed, absolute = kernel_integrals(ell)
    return {
        "transmission_at_lowest_bin": float(abs(transmitted[0])),
        "lowest_bin_frequency": float(omega[0]),
        "transmission_at_exact_dc_is_zero_analytically": (
            "T(w = 0) = 0 exactly: a zero-frequency wave is completely "
            "reflected by the centrifugal barrier. The tabulated value above is "
            "the lowest sampled bin, not w = 0, and is reported as such."),
        "transmission_at_high_frequency": float(abs(transmitted[-1])),
        "sum_rule_exact_value": -1.0,
        "sum_rule_measured": signed,
        "sum_rule_relative_error": abs(signed + 1.0),
        "the_sum_rule_holds": bool(abs(signed + 1.0) < 1e-4),
        "memory_absolute_mass": absolute,
        "instantaneous_weight": 1.0,
        "the_kernel_is_not_rigid": bool(absolute > 1.0),
        "the_durable_statement": (
            "With the chosen asymptotic normalisation T(inf) = 1 while "
            "T(0) = 0. Therefore the transfer function cannot be a rigid phase "
            "or delay, and the regular memory must carry signed weight -1. "
            "This is analytic and does not depend on any measured number."),
        "why_the_sum_rule_is_exact": (
            "int K dt = T(w = 0), and T(0) = 0 because the centrifugal barrier "
            "reflects zero-frequency waves completely. With K = delta(t) + "
            "K_reg, that forces int K_reg dt = -1 EXACTLY. It is a constraint, "
            "not a fit, and the measured value is checked against it."),
        "what_this_says_about_rigid_exchange": (
            "A rigid exchange kernel is delta(t), possibly delayed: whatever "
            "enters leaves undistorted, and a static signal passes perfectly. "
            "The real geometry blocks a static signal COMPLETELY, and does so "
            "entirely through the memory term, which exactly cancels the "
            "instantaneous one at DC. In absolute mass the memory is about "
            "twice the delta."),
        "scope_of_that_claim": (
            "This is a statement about the transfer kernel of a test scalar on "
            "a fixed D = 5 Tangherlini background, per angular channel. It says "
            "what the causal geometry does. Whether any particular BAM exchange "
            "kernel is meant to approximate THIS object is a separate question "
            "this round does not settle. Note also that T maps the asymptotic "
            "exterior to flux crossing the FUTURE HORIZON: it is a one-way "
            "exterior-to-interior channel, NOT mouth-to-mouth transmission."),
    }


@lru_cache(maxsize=4)
def measure_the_kernel_integrals_are_converged(
        ell: int = 1) -> Dict[str, object]:
    """R5b — the memory mass earns its digits, or it does not get quoted.

    Three axes. ``decay`` is the cleanest: the kernel cannot depend on the
    arbitrary subtraction parameter *at all*, since ``A`` is removed
    numerically and added back in closed form, so any dependence measures the
    inadequacy of ``ω_max``. The time quadrature is the axis that was actually
    broken before midpoint sampling.
    """
    by_decay = [{"decay": a, "sum": kernel_integrals(ell, decay=a)[0],
                 "mass": kernel_integrals(ell, decay=a)[1]}
                for a in (0.5, 1.0, 2.0, 4.0)]
    by_spectrum = []
    for omega_max, count in ((20.0, 2048), (40.0, 4096), (40.0, 8192)):
        signed, mass = kernel_integrals(ell, omega_max=omega_max, count=count)
        by_spectrum.append({"omega_max": omega_max, "count": count,
                            "sum": signed, "mass": mass})
    by_time = []
    for spacing, horizon in ((0.02, 300.0), (0.01, 300.0), (0.005, 300.0),
                             (0.005, 150.0), (0.005, 600.0)):
        signed, mass = kernel_integrals(ell, spacing=spacing, horizon=horizon)
        by_time.append({"spacing": spacing, "horizon": horizon,
                        "sum": signed, "mass": mass})

    def spread(rows, key):
        return max(r[key] for r in rows) - min(r[key] for r in rows)

    spreads = {"decay": spread(by_decay, "mass"),
               "spectrum": spread(by_spectrum, "mass"),
               "time": spread(by_time, "mass")}
    worst = max(spreads.values())
    sum_spread = max(spread(rows, "sum")
                     for rows in (by_decay, by_spectrum, by_time))
    return {
        "by_decay": by_decay,
        "by_spectrum": by_spectrum,
        "by_time_quadrature": by_time,
        "mass_spreads": spreads,
        "worst_mass_spread": worst,
        "worst_sum_spread": sum_spread,
        # The spread sets the precision, so the quoted value carries three
        # significant figures and no more. The widest axis is the subtraction
        # parameter at a = 4, where the analytic piece decays fastest and the
        # most weight falls on the numerical transform at finite w_max.
        "the_mass_is_converged": bool(worst < 5e-3),
        "significant_figures_earned": 3,
        "the_mass_to_quoted_precision": round(
            0.5 * (max(r["mass"] for r in by_decay)
                   + min(r["mass"] for r in by_decay)), 2),
        "the_sum_rule_holds_across_every_knob": bool(
            all(abs(r["sum"] + 1.0) < 1e-3
                for rows in (by_decay, by_spectrum, by_time) for r in rows)),
        "independent_of_the_subtraction_parameter": bool(
            spreads["decay"] < 5e-3),
        "the_precision_the_spread_supports": (
            "Worst spread ~2e-3, on the subtraction parameter. So the memory "
            "mass is 2.03 to the digits this computation earns; quoting "
            "2.0286 or 2.0309 would claim precision the knobs do not support."),
        "why_decay_independence_is_the_clean_check": (
            "K cannot depend on a at all: A(w) is subtracted numerically and "
            "its exact transform added back. Any residual dependence is "
            "measuring the finite w_max, not physics."),
        "what_was_broken_before": (
            "The time quadrature. A left-endpoint grid starting at t = 0.001 "
            "omitted int_0^t0 ~ -c t0 across the jump at t = 0+, giving "
            "-0.997757 for a quantity whose exact value is -1, and shrinking "
            "only as O(dt) rather than converging. Midpoint sampling over "
            "[0, t_max] gives -0.999996, a thousandfold tighter check of the "
            "same analytic constraint."),
    }


@lru_cache(maxsize=4)
def measure_the_late_time_tail_is_not_resolved(
        ell: int = 1) -> Dict[str, object]:
    """R6b — the tail is out of reach, argued from **positive** time.

    The negative-time floor says a transform artefact is small; it does not say
    a positive-time feature is real, because a causal error lives only at
    ``t > 0``. So the claim that the late-time tail is unmeasured is made the
    way it should be: vary the parameters and see whether the late-time values
    agree with each other. They do not — the spread across settings exceeds the
    values, and they change sign.
    """
    stamps = np.array([40.0, 60.0, 100.0, 150.0, 200.0])
    settings = [
        ("base", {}),
        ("omega_max=20, count=2048", {"omega_max": 20.0, "count": 2048}),
        ("count=8192", {"count": 8192}),
        ("decay=0.5", {"decay": 0.5}),
        ("decay=2.0", {"decay": 2.0}),
    ]
    rows = []
    for label, kwargs in settings:
        rows.append({"setting": label,
                     "values": [float(x) for x in
                                transfer_kernel(stamps, ell, **kwargs)]})
    stacked = np.array([r["values"] for r in rows])
    spread = (stacked.max(axis=0) - stacked.min(axis=0))
    magnitude = np.abs(stacked).max(axis=0)
    signs_disagree = bool(np.any(stacked.max(axis=0) * stacked.min(axis=0) < 0))
    return {
        "times": [float(t) for t in stamps],
        "rows": rows,
        "spread_across_settings": [float(x) for x in spread],
        "largest_magnitude": [float(x) for x in magnitude],
        "spread_exceeds_the_values": bool(np.all(spread > 0.5 * magnitude)),
        "the_sign_is_not_even_stable": signs_disagree,
        "the_tail_is_not_resolved": bool(
            np.all(spread > 0.5 * magnitude) and signs_disagree),
        "why_this_is_the_right_argument": (
            "The negative-time floor bounds acausal transform artefacts only. "
            "Whether a positive-time feature is real has to be settled at "
            "positive times, by varying the parameters it would be independent "
            "of. Here it is not: the spread exceeds the values and the sign "
            "flips, so no exponent is quoted."),
        "what_would_be_needed": (
            "A method with dynamic range where this one has none -- a long "
            "time-domain evolution in extended precision. Not a refinement of "
            "the frequency-domain route, whose range is excellent near the "
            "barrier and absent in the tail."),
    }


def measure_what_the_causality_gate_caught() -> Dict[str, object]:
    """R6 — two artefacts that would have been read as physics."""
    return {
        "missing_dc_cell": {
            "symptom": ("a constant -1.9e-3 under the whole kernel, at t < 0 "
                        "as well as t > 0"),
            "cause": ("omega sampled at right endpoints left [0, dw] "
                      "uncovered; since T(0) = 0, S(0) = -1 and the omitted "
                      "cell contributes about -dw/pi everywhere"),
            "fix": "midpoint sampling, which covers [0, w_max] exactly",
            "would_have_been_read_as": (
                "a constant offset in the kernel, i.e. a spurious permanent "
                "memory"),
        },
        "gibbs_ringing_from_the_slow_tail": {
            "symptom": ("a 1.6e-3 plateau at large |t|, on BOTH sides of the "
                        "origin"),
            "cause": ("S = T - 1 decays only as -i c/w, so a truncated "
                      "transform rings"),
            "fix": ("subtract A(w) = -i c/(w + i a), whose only pole is at "
                    "w = -ia in the lower half plane so it cannot introduce an "
                    "acausal piece, and add back its closed-form transform "
                    "-c e^{-at} theta(t); the remainder decays like 1/w^2"),
            "would_have_been_read_as": (
                "a late-time tail -- the very quantity this round would most "
                "like to measure, which is exactly why the gate mattered"),
        },
        "the_general_point": (
            "A quantity that is exactly zero on part of its domain is worth "
            "more than an accuracy claim: it calibrates the noise floor for "
            "free, on the same run, with no reference value. Both artefacts "
            "above sat at the level of the physics being looked for, and "
            "neither was visible from the t > 0 side alone."),
    }


def measure_the_transfer_kernel_ledger() -> Dict[str, object]:
    """R7 — what is settled and what is not.

    The numbers here are **derived from the measurements**, not embedded as
    literals. An earlier draft hard-coded the cross-check percentages, and when
    the Jost patch changed them the ledger silently kept the old ones -- the
    exact stale-number failure this arc has been trying to eliminate.
    """
    cross = measure_the_kernel_against_the_time_domain_evolution()
    rigid = measure_the_transfer_is_not_rigid()
    causal = measure_the_kernel_is_causal()
    ringdown = measure_the_kernel_reproduces_the_published_ringdown()
    best = cross["best"]
    entries = [
        {"claim": "the retarded outer-to-inner transfer kernel exists as a "
                  "computable object",
         "verdict": "YES, DELIVERED",
         "evidence": "K = delta(t) + K_reg from T(w) on a well-conditioned "
                     "real-frequency scattering problem"},
        {"claim": "the frequency domain can be used here despite PR #270's "
                  "and #274's shooting failures",
         "verdict": "YES -- DIFFERENT PROBLEM",
         "evidence": "those failed for Im w < 0, where one solution grows "
                     "exponentially; for real w nothing dominates, and "
                     "unitarity holds to ~1e-13 without being imposed"},
        {"claim": "the kernel is causal",
         "verdict": f"YES, TO ~{causal['noise_floor_far_from_the_front']:.0e}",
         "evidence": "K(t < 0) measured directly; the residual bounds ACAUSAL "
                     "transform artifacts, not the total numerical error"},
        {"claim": "the kernel carries the published ringdown",
         "verdict": "YES",
         "evidence": "fitted against the EXTERNAL continued-fraction value: "
                     f"real part within "
                     f"{100*ringdown['real_part_band']['min']:.3f}%-"
                     f"{100*ringdown['real_part_band']['max']:.2f}%, damping "
                     f"{100*ringdown['damping_band']['min']:.2f}%-"
                     f"{100*ringdown['damping_band']['max']:.2f}% over "
                     f"{len(ringdown['rows'])} extraction settings; band "
                     "reported, not best row"},
        {"claim": "an independent propagation method reproduces the kernel",
         "verdict": (f"YES, TO {100*best['peak_relative_max_difference']:.2f}% "
                     f"PEAK / {100*best['peak_relative_rms_difference']:.2f}% "
                     f"RMS at launch r* = {best['launch_r_star']:g}"),
         "evidence": "convolution against PR #274's characteristic null-grid "
                     "evolution -- a different propagation algorithm on the "
                     "same operator, not an independent operator; the residual "
                     "halves as the launch radius doubles"},
        {"claim": "the transfer is rigid / instantaneous",
         "verdict": "NO -- AND ANALYTICALLY SO",
         "evidence": "T(inf) = 1 while T(0) = 0, so int K_reg dt = -1 is "
                     "FORCED; measured "
                     f"{rigid['sum_rule_measured']:.6f}. The memory mass "
                     f"{rigid['memory_absolute_mass']:.4f} against the delta's "
                     "1 is a quantitative extra, not the argument"},
        {"claim": "the late-time power-law tail is measured",
         "verdict": "NO",
         "evidence": "argued at POSITIVE times, not from the causality floor: "
                     "varying (w_max, count) and the subtraction parameter "
                     "gives late-time values whose spread exceeds the values "
                     "and whose sign is not stable. No exponent is quoted"},
        {"claim": "PR #274's pulse placement was adequate for this object",
         "verdict": "NO -- AND IT WAS FOR ITS OWN",
         "evidence": "the incident amplitude is only defined where the wave is "
                     "free; launching inside the barrier's reach was harmless "
                     "for a quasinormal frequency and fatal for a transmission "
                     "ratio"},
        {"claim": "the causality gate bounds the total numerical error",
         "verdict": "NO -- ACAUSAL ARTIFACTS ONLY",
         "evidence": "a causal error can live entirely at t > 0; this round's "
                     "own Jost patch moved positive-time results with no "
                     "matching negative-time signature"},
        {"claim": "the high-frequency deficit was the 1/r^4 and 1/r^6 tail",
         "verdict": "NO -- IT WAS THE ESTIMATOR",
         "evidence": "Re[i w (T-1)] carries a deterministic -c^3/(6 w^2) bias "
                     "of -1.3e-3, the whole size of the reported deficit; the "
                     "unbiased -w arg T gives +2.6e-4, opposite in sign"},
        {"claim": "the memory mass 2.03 is a converged number",
         "verdict": "YES, NOW",
         "evidence": "stable across the subtraction parameter, (w_max, count) "
                     "and the time quadrature; an earlier draft's -0.997757 "
                     "sum rule was a left-endpoint grid missing the [0, dt] "
                     "cell across the jump at t = 0+"},
    ]
    return {
        "entries": entries,
        "the_lesson_this_round_adds": (
            "An exactly-zero region is a free monitor for the artefacts that "
            "violate it. Causality gave this round a stretch of the domain "
            "where the answer is known to be zero, and two separate ACAUSAL "
            "artefacts -- a missing DC cell and Gibbs ringing -- were caught "
            "there at exactly the amplitude of the physics being sought. It "
            "bounds that family of error and no other: a causal error lives "
            "only at t > 0, as this round's own Jost patch demonstrated."),
        "the_second_lesson": (
            "Check what an estimator measures before attributing what it "
            "returns. Re[i w (T-1)] carries a deterministic -c^3/(6 w^2) bias, "
            "and an earlier draft read the resulting shortfall as a physical "
            "property of the potential's sub-leading tail. The tell was there: "
            "the deficit was independent of the outer edge, which a truncation "
            "effect would not be and a fixed-form bias would."),
        "the_next_object": (
            "The late-time tail, which needs a method with dynamic range where "
            "this one has none -- a long time-domain evolution in extended "
            "precision rather than a refinement of the frequency-domain route. "
            "Separately: whether any BAM exchange kernel is intended to "
            "approximate this object, which is a modelling question and not a "
            "numerical one."),
        "still_blocked": (
            "The nonlinear backreaction oracle, by the C(v) / inner-flux issue, "
            "unchanged and unrelated to anything here."),
    }
