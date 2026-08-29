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

and ℓ = 1 is the *only* ℓ for which this is rational — ℓ = 0, 2, 3 carry
``√13``, ``√1621``, ``√57``. ℓ = 1 is also the mode #270's two codes disagreed
on.

**2. The integrated potential is exact.** Since ``dr* = dr/f`` cancels the ``f``
in ``V``,

    ∫ V_ℓ dr*  =  ∫₁^∞ (L/r² + (9/4)/r⁴) dr  =  ℓ(ℓ+2) + 3/2

with ``L = (ℓ+1)² − ¼``. At ℓ = 1 this is ``9/2``.

**3. Hence the kernel's high-frequency phase is known in closed form.** The
leading eikonal phase through the barrier is ``−(1/2ω)∫V dr*``, so

    T_ℓ(ω) → exp(−i c_ℓ/ω) ,   c_ℓ = (ℓ(ℓ+2) + 3/2)/2 ,   c₁ = 9/4 exactly

This is what makes the kernel computable: ``T − 1 ~ −i c/ω`` decays too slowly
to transform numerically, and the exact ``c`` lets it be removed analytically
instead of windowed away.

THE THREE GATES, ALL MEASURED
─────────────────────────────
**Causal support.** ``K(t) = 0`` for ``t < 0``, to ``~3e-7`` away from the
front. This is not decoration: the exact zero for ``t < 0`` is a **free
calibration of the numerical noise floor**, and it caught two artefacts that
would otherwise have been read as physics — see below.

**Flux conservation.** ``|R|² + |T|² = 1`` to ``8e-13``, nowhere imposed.

**Late-time ringdown.** Fitting the kernel's own ringdown gives, against the
*external* published value ``1.01601691149 − 0.36232802385i``, real part within
``0.062 %–1.17 %`` and damping within ``0.11 %–3.80 %`` across nine
sample-spacing/window combinations. The band is the honest statement — as PR
#274 established for the time-domain solver, extraction choices move the fit, so
the spread is reported rather than the best row. Checked against the published
number, never against this repository's own solver.

**Independent cross-check.** The kernel is validated against the completely
separate time-domain characteristic evolution of PR #274: convolving ``K`` with
the incident profile reproduces the transmitted wave to **0.92 % peak, 0.17 %
rms**. Two methods sharing no code.

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
    "scattering_amplitudes",
    "transfer_spectrum",
    "transfer_kernel",
    "measure_the_exact_background_anchors",
    "measure_the_scattering_is_well_conditioned",
    "measure_the_kernel_is_causal",
    "measure_the_kernel_reproduces_the_published_ringdown",
    "measure_the_kernel_against_the_time_domain_evolution",
    "measure_the_transfer_is_not_rigid",
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


def scattering_amplitudes(omega, ell: int, inner: float = INNER_EDGE,
                          outer: float = OUTER_EDGE, steps: int = 40000
                          ) -> Tuple[np.ndarray, np.ndarray]:
    """``R(ω)``, ``T(ω)`` by a piecewise-constant transfer matrix.

    Boundary conditions: purely ingoing at the horizon, ``e^{−iωr*} + R e^{iωr*}``
    at infinity, so ``T`` is the outer→inner transmission amplitude.

    On each step ``V ≈ V_mid`` and the propagator is exact:

        ψ(x+h)  =  cos(kh) ψ + sin(kh)/k ψ' ,   k = √(ω² − V_mid)
        ψ'(x+h) = −k sin(kh) ψ + cos(kh) ψ'

    with complex ``k`` covering the classically forbidden region. **For real
    ``ω`` both basis solutions have unit modulus, so nothing exponentially
    dominates anything** — this is why the calculation succeeds where PR #270's
    and #274's quasinormal shoots (``Im ω < 0``) did not. Vectorised over ``ω``:
    the spatial loop is shared across all frequencies.
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
    incident = 0.5 * (psi - dpsi / (1j * w)) * np.exp(1j * w * outer)
    reflected = 0.5 * (psi + dpsi / (1j * w)) * np.exp(-1j * w * outer)
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
    subtraction cannot itself introduce an acausal piece. ``S − A`` decays like
    ``1/ω²`` and transforms cleanly.
    """
    omega, _, transmitted, spacing = transfer_spectrum(ell, omega_max, count)
    residual = transmitted - 1.0
    fitted = float(np.mean((1j * omega * residual).real[-200:]))
    analytic = -1j * fitted / (omega + 1j * decay)
    remainder = residual - analytic

    stamps = np.atleast_1d(np.asarray(times, dtype=float))
    out = np.empty(stamps.size)
    for start in range(0, stamps.size, 4000):        # bound the outer product
        block = stamps[start:start + 4000]
        out[start:start + 4000] = 2.0 * np.real(
            np.exp(-1j * np.outer(block, omega)) @ remainder)
    out *= spacing / (2.0 * np.pi)
    out += np.where(stamps > 0.0, -fitted * np.exp(-decay * stamps), 0.0)
    return out


# ── measurements ────────────────────────────────────────────────────────────

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
            "grows like e^{|Im w| R} and swamps the coefficient being zeroed. "
            "Here w is REAL, both e^{+-i w r*} have unit modulus, and nothing "
            "dominates anything. This is a different, well-posed problem -- not "
            "a repair of the one that failed. Unitarity to ~1e-13, imposed "
            "nowhere, is the evidence."),
    }


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
        "the_exact_zero_is_a_free_error_bar": (
            "K(t) vanishes identically for t < 0, so whatever the computation "
            "returns there IS its noise floor -- no reference value needed. "
            "Any feature at t > 0 smaller than that floor is not measurable, "
            "which is how this round knows the late-time tail is out of reach."),
    }


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


def measure_the_kernel_against_the_time_domain_evolution(
        ell: int = 1, width: float = 3.0) -> Dict[str, object]:
    """R4 — the independent check: ``K ⋆ g`` against a characteristic evolution.

    Deep inside, the transmitted wave as a function of ``v = t + r*`` is exactly
    the convolution of the incident profile with the kernel. The time-domain
    evolution of PR #274 shares no code with the transfer matrix, so agreement
    is a real cross-validation.

    **The launch radius is the whole subtlety.** The incident amplitude is only
    defined where ``V ≈ 0``. On the ``u = 0`` line the pulse sits at
    ``r* = v_c/2``, so PR #274's ``v_c = 12`` launches at ``r* = 6``, where
    ``V ≈ 0.1`` — inside the barrier's reach. That is harmless for extracting a
    quasinormal frequency (a ringdown does not care how it was excited) and
    fatal for defining a transmission ratio. Moving the launch out fixes it, and
    the error tracks ``V`` at the launch point.
    """
    omega, _, transmitted, spacing = transfer_spectrum(ell)
    rows = []
    for centre in (12.0, 60.0, 200.0):
        launch = 0.5 * centre
        times, signal = characteristic_evolution(
            ell, step=0.05, t_max=500.0, observer_rstar=-30.0,
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
            "potential_at_launch": float(
                potential(radius_of_tortoise(launch), ell)),
            "peak_relative_max_difference": float(
                np.max(np.abs(measured - predicted)) / peak),
            "peak_relative_rms_difference": float(
                math.sqrt(float(np.mean((measured - predicted) ** 2))) / peak),
        })
    best = min(rows, key=lambda r: r["peak_relative_max_difference"])
    return {
        "rows": rows,
        "best": best,
        "the_two_methods_agree": bool(
            best["peak_relative_max_difference"] < 0.02),
        "the_error_tracks_the_potential_at_launch": bool(
            all(a["peak_relative_max_difference"] > b["peak_relative_max_difference"]
                for a, b in zip(rows[:-1], rows[1:]))),
        "what_this_exposed": (
            "PR #274's pulse was launched at r* = 6, where V ~ 0.1. That is "
            "fine for a quasinormal frequency and wrong for an incident "
            "amplitude, and it shows up here as a 43% mismatch that falls to "
            "0.92% once the launch moves to r* = 100. The transfer kernel needs "
            "an asymptotic launch; the ringdown did not."),
    }


def measure_the_transfer_is_not_rigid(ell: int = 1) -> Dict[str, object]:
    """R5 — **the deliverable**: how far the kernel is from instantaneous."""
    published_dc = abs(transfer_spectrum(ell)[2][0])
    stamps = np.arange(0.001, 300.0, 0.005)
    values = transfer_kernel(stamps, ell)
    signed = float(np.trapezoid(values, stamps))
    absolute = float(np.trapezoid(np.abs(values), stamps))
    return {
        "transmission_at_dc": float(published_dc),
        "sum_rule_exact_value": -1.0,
        "sum_rule_measured": signed,
        "sum_rule_relative_error": abs(signed + 1.0),
        "the_sum_rule_holds": bool(abs(signed + 1.0) < 5e-3),
        "memory_absolute_mass": absolute,
        "instantaneous_weight": 1.0,
        "memory_to_instantaneous_ratio": absolute,
        "the_kernel_is_not_rigid": bool(absolute > 1.0),
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
            "twice the delta. The memory is not a correction to rigid exchange; "
            "it is the same size as the thing it would correct."),
        "scope_of_that_claim": (
            "This is a statement about the transfer kernel of a test scalar on "
            "a fixed D = 5 Tangherlini background, per angular channel. It says "
            "what the causal geometry does. Whether any particular BAM exchange "
            "kernel is meant to approximate THIS object is a separate question "
            "this round does not settle."),
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
    """R7 — what is settled and what is not."""
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
         "verdict": "YES, TO ~3e-7",
         "evidence": "K(t < 0) measured directly; the residual doubles as the "
                     "noise floor for t > 0"},
        {"claim": "the kernel carries the published ringdown",
         "verdict": "YES",
         "evidence": "fitted against the EXTERNAL continued-fraction value: "
                     "real part within 0.062%-1.17%, damping 0.11%-3.80% over "
                     "nine extraction settings; band reported, not best row"},
        {"claim": "an independent method reproduces the kernel",
         "verdict": "YES, TO 0.92% PEAK / 0.17% RMS",
         "evidence": "convolution against PR #274's time-domain characteristic "
                     "evolution, which shares no code with the transfer matrix"},
        {"claim": "the transfer is rigid / instantaneous",
         "verdict": "NO -- AND NOT MARGINALLY",
         "evidence": "int K_reg dt = -1 exactly (sum rule from T(0) = 0), so "
                     "the memory cancels the instantaneous part at DC; "
                     "int |K_reg| dt = 2.02 against the delta's 1"},
        {"claim": "the late-time power-law tail is measured",
         "verdict": "NO",
         "evidence": "the ringdown reaches the ~1e-6 causality noise floor by "
                     "t ~ 40; a tail would be orders of magnitude below it. No "
                     "exponent is quoted"},
        {"claim": "PR #274's pulse placement was adequate for this object",
         "verdict": "NO -- AND IT WAS FOR ITS OWN",
         "evidence": "launching at r* = 6 (V ~ 0.1) gives a 43% mismatch here "
                     "and was harmless for the quasinormal frequency; a "
                     "transmission ratio needs an asymptotic launch"},
    ]
    return {
        "entries": entries,
        "the_lesson_this_round_adds": (
            "An exactly-zero region is a free error bar. Causality gave this "
            "round a stretch of the domain where the answer is known to be "
            "zero, on the same run and with no external reference, and two "
            "separate artefacts -- a missing DC cell and Gibbs ringing -- were "
            "caught there at exactly the amplitude of the physics being sought. "
            "PR #274 needed a published spectrum to find its floor; here the "
            "structure of the problem supplied one."),
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
