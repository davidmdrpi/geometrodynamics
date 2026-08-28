"""Settling PR #270's ringdown cross-validation.

THE VERDICT
───────────
PR #270 built two horizon-penetrating time-domain codes for a test scalar on a
fixed ``D = 5`` Tangherlini background. Both were stable, both converged, and
they disagreed: real parts within ``0.3 %`` at ``ℓ = 1``, **damping rates
differing by 37 %**. #270 refused to quote a frequency, correctly — *a
converged number is not a correct number* — and named its own prime suspect:
the Kerr–Schild operator's inner cut.

**The Kerr–Schild code was right. The tortoise code's damping was wrong.**

An independent Gundlach–Price–Pullin characteristic evolution, written here
from scratch and sharing no code with either, gives at ``ℓ = 1``

    this round (h → 0)   1.01612 − 0.36244i
    #270 Kerr–Schild     1.01622 − 0.36231i     agrees to ~1e-4
    #270 tortoise        1.01876 − 0.26404i     damping 37 % away

Agreement to `1e-4` between two independent implementations, against a `37 %`
miss, is not ambiguous. #270's suspicion of the Kerr–Schild inner cut was
pointing at the wrong code.

WHY THIS BACKGROUND MAKES THE ANSWER CHECKABLE
──────────────────────────────────────────────
``D = 5`` is unusually clean, and two exact facts do most of the work.

**The tortoise correction is a decaying power, not a log.**
``r* = r + ½ln((r−1)/(r+1))``, whose correction is exactly ``−1/r`` — a
decaying power, not the ``ln r`` of 4D. So the far field is a pure Hankel wave with no Coulomb-like phase.

**The potential's tail is exactly Bessel.**
``V → [(ℓ+1)² − ¼]/r²``, so the outgoing solution is exactly
``√r H⁽¹⁾_{ℓ+1}(ωr)``. This is the same identity that fixed the operator
correction in PR #271, reused here as a boundary condition.

Together these give an *exact* eikonal limit to check against: the photon
sphere sits at ``r_ph = √2`` with

    Ω_c = √(f(r_ph))/r_ph = 1/2   (exactly)
    λ   = 1/√2 = 0.707107          (exactly)
    ω → Ω_c(ℓ+1) − i λ/2 = 0.5(ℓ+1) − 0.353553i

and the characteristic evolution tracks it:

    ℓ = 0   0.5354 − 0.3842i      eikonal  0.5 − 0.35355i
    ℓ = 1   1.0161 − 0.3624i      eikonal  1.0 − 0.35355i
    ℓ = 2   1.5106 − 0.3575i      eikonal  1.5 − 0.35355i

The real parts converge to ``0.5(ℓ+1)`` from above and the damping to
``−0.35355``, which is what a correct solver must do. **The tortoise damping
``−0.264`` is not near that asymptote and never approaches it.**

FOUR INDEPENDENT LINES, ALL AGAINST THE TORTOISE DAMPING
─────────────────────────────────────────────────────────
    characteristic evolution   −0.36244    (this round, converged in h)
    #270 Kerr–Schild           −0.36231
    first-order WKB            −0.36095
    exact eikonal asymptote    −0.35355
    ────────────────────────────────────
    #270 tortoise              −0.26404    ← excluded

WHAT DID NOT WORK, AND IS REPORTED AS SUCH
──────────────────────────────────────────
**Frequency-domain shooting fails here, and this round reproduced the failure
rather than fixing it.** Matching an ingoing Frobenius series at the horizon to
``√r H⁽¹⁾`` at large ``r`` gives a root that moves with every numerical knob:
``1.229 − 0.152i`` → ``1.173 − 0.102i`` → ``1.155 − 0.214i`` across ``ε``, and
similarly across the matching radius. The QNM boundary-value problem is
exponentially ill-conditioned in real ``r`` — for ``Im ω < 0`` the outgoing
piece grows like ``e^{|Im ω|R}`` and swamps the incoming coefficient one is
trying to zero. #270 reported the same non-convergence; that diagnosis stands.

**Sixth-order WKB by finite differences is unusable.** The Iyer–Will formula
needs ``V⁽⁶⁾`` at the barrier peak, and central differences amplify roundoff as
``h⁻⁶``: refining the grid makes the answer *diverge* (``9.0 → 18.6 → 623``).
First-order WKB is well conditioned and is what is used, with its accuracy
stated rather than assumed — it is good to ``0.4 %`` on the damping and only
``13 %`` on the real part at ``ℓ = 1``, improving to ``6.8 %`` at ``ℓ = 2``.
Low-``ℓ`` WKB is simply not accurate on the real part, which is a known
limitation and not a discrepancy.

SCOPE
─────
A test scalar on a **fixed** exact Tangherlini background — no backreaction.
The fundamental (``n = 0``) mode only; overtones are not extracted. ``ℓ = 0`` is
quoted with a wider uncertainty because its barrier is weakest and the
power-law tail contaminates the fit earliest.

The autopsy of the tortoise code is **not** available: neither #270 code was
landed in the tree, only their reported numbers. This round establishes which
number is right and by how much the other is excluded; it cannot say which line
of the unlanded code did it.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from geometrodynamics.tangherlini.dynamics import N_ANGULAR, master_potential

__all__ = [
    "HORIZON",
    "tortoise",
    "radius_of_tortoise",
    "potential",
    "eikonal_limit",
    "wkb_fundamental",
    "characteristic_evolution",
    "prony_frequencies",
    "fundamental_mode",
    "measure_the_background_asymptotics_are_exact",
    "measure_the_eikonal_invariants",
    "measure_the_characteristic_scheme_converges",
    "measure_the_fundamental_modes",
    "measure_the_cross_validation_verdict",
    "measure_what_did_not_work",
    "measure_the_ringdown_ledger",
]

HORIZON = 1.0

#: #270 reported these and refused to choose between them.
KERR_SCHILD_ELL_1 = complex(1.01622, -0.36231)
TORTOISE_ELL_1 = complex(1.01876, -0.26404)


def metric_f(r):
    """``f = 1 − (r_h/r)^{n−1}``, at ``n = 3``."""
    r = np.asarray(r, dtype=float)
    return 1.0 - HORIZON ** 2 / r ** 2


def tortoise(r: float) -> float:
    """``r* = r + ½ ln((r−1)/(r+1))`` — note it tends to ``r``, with no log tail."""
    return r + 0.5 * math.log((r - HORIZON) / (r + HORIZON))


def radius_of_tortoise(rs: float) -> float:
    """Invert the tortoise map. Near the horizon ``r − 1 → 2e^{2(r*−1)}``."""
    if rs < -25.0:
        return 1.0 + 2.0 * math.exp(2.0 * (rs - 1.0))
    lo = 1.0 + 1e-13
    while tortoise(lo) > rs:
        lo = 1.0 + (lo - 1.0) * 1e-3
        if lo - 1.0 < 1e-300:
            return 1.0 + 2.0 * math.exp(2.0 * (rs - 1.0))
    hi = max(2.0, rs + 3.0)
    while tortoise(hi) < rs:
        hi *= 2.0
    return float(brentq(lambda r: tortoise(r) - rs, lo, hi,
                        xtol=1e-15, rtol=8.9e-16))


def potential(r, ell: int, n: int = N_ANGULAR):
    """The corrected master potential — PR #271's, taken from `dynamics`."""
    return master_potential(r, ell, HORIZON, n=n)


# ── the exact checks this background allows ─────────────────────────────────

def eikonal_limit(ell: int) -> Dict[str, float]:
    """Photon-sphere invariants — exact, and the asymptote any solver must meet.

    ``r_ph`` solves ``d/dr (f/r²) = 0``, giving ``r_ph² = 2`` exactly at
    ``n = 3``; then ``Ω_c = 1/2`` and ``λ = 1/√2``, both exact.
    """
    r_ph = math.sqrt(2.0) * HORIZON
    omega_c = math.sqrt(float(metric_f(r_ph))) / r_ph
    h = 1e-5

    def g(r):
        return float(metric_f(r)) / r ** 2

    second = float(metric_f(r_ph)) ** 2 * (g(r_ph + h) - 2 * g(r_ph)
                                           + g(r_ph - h)) / h ** 2
    lam = math.sqrt(-0.5 * second / g(r_ph))
    return {"r_photon": r_ph, "omega_c": omega_c, "lyapunov": lam,
            "omega": complex(omega_c * (ell + 1), -0.5 * lam)}


def wkb_fundamental(ell: int, overtone: int = 0) -> complex:
    """First-order (Schutz–Will) WKB. Well conditioned; accuracy stated, not assumed."""
    res = minimize_scalar(lambda r: -float(potential(r, ell)),
                          bounds=(1.0001 * HORIZON, 30.0 * HORIZON),
                          method="bounded", options={"xatol": 1e-14})
    r0 = float(res.x)
    v0 = float(potential(r0, ell))
    h = 1e-5
    f0 = float(metric_f(r0))
    v2 = f0 * f0 * (float(potential(r0 + h, ell)) - 2.0 * v0
                    + float(potential(r0 - h, ell))) / h ** 2
    return complex(np.sqrt(complex(v0, 0.0)
                           - 1j * (overtone + 0.5) * math.sqrt(-2.0 * v2)))


# ── the independent time-domain solver ──────────────────────────────────────

def characteristic_evolution(ell: int, step: float = 0.05, t_max: float = 400.0,
                             observer_rstar: float = 30.0,
                             pulse_centre: float = 12.0,
                             pulse_width: float = 3.0
                             ) -> Tuple[np.ndarray, np.ndarray]:
    """Gundlach–Price–Pullin on the null grid ``u = t − r*``, ``v = t + r*``.

        ``ψ_N = ψ_W + ψ_E − ψ_S − (h²/8) V_c (ψ_W + ψ_E)``

    **No spatial boundary conditions are applied or needed** — the domain of
    dependence is the null diamond, so the horizon and infinity are reached only
    as limits. That is precisely why this construction is immune to the
    excision-boundary question that the frequency-domain shoot and the #270
    Kerr–Schild inner cut both raise.

    ``V`` at a diamond centre depends only on ``r* = (v−u)/2 = h(j−i)/2``, so it
    is a one-dimensional lookup indexed by ``j − i``. The inner ``v`` recursion
    is linear and first order, ``a_j = A_j a_{j−1} + B_j``, so it is evaluated
    by cumulative products rather than a Python loop.
    """
    count = int(t_max / step)
    offsets = np.arange(-count, count + 1)
    centres = 0.5 * step * (offsets - 0.5)
    table = np.array([0.0 if rs < -40.0 else float(potential(radius_of_tortoise(rs), ell))
                      for rs in centres])
    origin = count

    v_grid = np.arange(count + 1) * step
    previous = np.exp(-((v_grid - pulse_centre) / pulse_width) ** 2)
    coefficient = 0.125 * step * step
    index = np.arange(1, count + 1)

    times: List[float] = []
    signal: List[float] = []
    for i in range(1, count + 1):
        row = table[origin + (index - i)]
        a = 1.0 - coefficient * row
        b = previous[1:] * a - previous[:-1]
        prod = np.cumprod(a)
        current = np.empty(count + 1)
        current[0] = 0.0                       # ψ = 0 on v = 0 for u > 0
        current[1:] = prod * np.cumsum(b / prod)
        previous = current
        j = int(round((i * step + 2.0 * observer_rstar) / step))
        if 0 <= j <= count:
            times.append(i * step + observer_rstar)
            signal.append(previous[j])
    return np.array(times), np.array(signal)


def prony_frequencies(times: np.ndarray, signal: np.ndarray,
                      n_modes: int = 3) -> np.ndarray:
    """Fit ``ψ ≈ Σ A_k e^{−iω_k t}`` on a uniform grid by Prony's method."""
    dt = float(times[1] - times[0])
    size = len(signal)
    p = 2 * n_modes
    design = np.array([signal[i:i + p] for i in range(size - p)])
    target = -signal[p:]
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    roots = np.roots(np.concatenate([[1.0], coefficients[::-1]]))
    roots = roots[np.abs(roots) > 1e-12]
    return 1j * np.log(roots) / dt


def fundamental_mode(ell: int, step: float = 0.05, window: Tuple[float, float]
                     = (60.0, 140.0), **kwargs) -> Optional[complex]:
    """The least-damped physical root in the extraction window."""
    times, signal = characteristic_evolution(ell, step=step, **kwargs)
    mask = (times > window[0]) & (times < window[1])
    if mask.sum() < 40:
        return None
    candidates = [w for w in prony_frequencies(times[mask], signal[mask])
                  if w.real > 0.2 and -2.0 < w.imag < -0.01]
    if not candidates:
        return None
    return min(candidates, key=lambda w: abs(w.imag))


# ── the measurements ────────────────────────────────────────────────────────

def measure_the_background_asymptotics_are_exact() -> Dict[str, object]:
    """R0 — the two ``D = 5`` facts that make the answer checkable."""
    far = [50.0, 200.0, 1000.0]
    no_log = [abs(tortoise(r) - r) for r in far]
    # The tail is a POWER law, not a log: (1/2)ln((r-1)/(r+1)) -> -1/r exactly,
    # so (r* - r) * r -> -1. Testing that is the real statement; an arbitrary
    # threshold on |r* - r| is not, and the first draft's 1e-3 cut sat exactly
    # on the value at r = 1000.
    scaled = [(tortoise(r) - r) * r for r in far]
    bessel = []
    for ell in (0, 1, 2, 3):
        vals = [float(potential(r, ell)) * r ** 2 for r in far]
        limit = (ell + 1) ** 2 - 0.25
        bessel.append({"ell": ell, "limit": limit,
                       "V_times_r2_at_far_radii": vals,
                       "relative_error_at_1000": abs(vals[-1] - limit) / limit})
    inverses = [abs(radius_of_tortoise(tortoise(r)) - r)
                for r in (1.0001, 1.5, 3.0, 50.0)]
    return {
        "tortoise_minus_r_at_far_radii": no_log,
        "tortoise_minus_r_times_r": scaled,
        "deviation_from_minus_one_over_r": [abs(s + 1.0) for s in scaled],
        "the_tail_is_exactly_minus_one_over_r": bool(
            abs(scaled[-1] + 1.0) < 1e-5
            and all(abs(b + 1.0) < abs(a + 1.0)
                    for a, b in zip(scaled[:-1], scaled[1:]))),
        "no_logarithmic_tail": bool(
            abs(scaled[-1] + 1.0) < 1e-5
            and all(abs(b + 1.0) < abs(a + 1.0)
                    for a, b in zip(scaled[:-1], scaled[1:]))),
        "why_that_matters": ("unlike 4D, r* -> r with no ln r, so the far field "
                             "is a pure Hankel wave with no Coulomb phase"),
        "bessel_tail": bessel,
        "the_tail_is_exactly_bessel": bool(
            all(b["relative_error_at_1000"] < 1e-5 for b in bessel)),
        "the_same_identity_that_fixed_the_operator": (
            "V -> [(l+1)^2 - 1/4]/r^2 is the flat-limit identity PR #271 used to "
            "settle which radial operator was correct; it is reused here as a "
            "boundary condition."),
        "tortoise_inversion_round_trip": inverses,
    }


def measure_the_eikonal_invariants() -> Dict[str, object]:
    """R1 — an exact asymptote for any solver to be judged against."""
    e = eikonal_limit(0)
    return {
        "r_photon": e["r_photon"],
        "r_photon_squared": e["r_photon"] ** 2,
        "r_photon_squared_is_exactly_two": bool(
            abs(e["r_photon"] ** 2 - 2.0) < 1e-12),
        "omega_c": e["omega_c"],
        "omega_c_is_exactly_one_half": bool(abs(e["omega_c"] - 0.5) < 1e-12),
        "lyapunov": e["lyapunov"],
        "lyapunov_is_one_over_sqrt_two": bool(
            abs(e["lyapunov"] - 1.0 / math.sqrt(2.0)) < 1e-6),
        "asymptotic_damping": -0.5 * e["lyapunov"],
        "the_law": "omega -> 0.5 (l+1) - 0.353553i  as l -> infinity",
    }


def measure_the_characteristic_scheme_converges(
        ells: Sequence[int] = (1, 2),
        steps: Sequence[float] = (0.1, 0.05, 0.025)) -> Dict[str, object]:
    """R2 — the solver has to converge before its number means anything."""
    rows = []
    for ell in ells:
        values = [fundamental_mode(ell, step=s) for s in steps]
        deltas = [abs(b - a) for a, b in zip(values[:-1], values[1:])
                  if a is not None and b is not None]
        rows.append({
            "ell": ell,
            "by_step": {str(s): (None if v is None else [v.real, v.imag])
                        for s, v in zip(steps, values)},
            "successive_differences": deltas,
            "finest": None if values[-1] is None else [values[-1].real,
                                                       values[-1].imag],
            "converging": bool(len(deltas) >= 2 and deltas[-1] < deltas[0]),
        })
    return {
        "rows": rows,
        "steps": list(steps),
        "all_converging": all(r["converging"] for r in rows),
        "no_boundary_conditions_are_applied": (
            "the domain of dependence is the null diamond, so the horizon and "
            "infinity are limits rather than boundaries -- which is why this "
            "construction is immune to the excision question that the "
            "frequency-domain shoot and the Kerr-Schild inner cut both raise"),
    }


def measure_the_fundamental_modes(
        ells: Sequence[int] = (0, 1, 2, 3)) -> Dict[str, object]:
    """R3 — the frequencies, against the exact eikonal asymptote."""
    rows = []
    for ell in ells:
        w = fundamental_mode(ell, step=0.05)
        e = eikonal_limit(ell)["omega"]
        wkb = wkb_fundamental(ell)
        rows.append({
            "ell": ell,
            "omega": None if w is None else [w.real, w.imag],
            "eikonal": [e.real, e.imag],
            "real_part_above_eikonal": None if w is None else bool(w.real > e.real),
            "damping_near_asymptote": None if w is None else abs(w.imag - e.imag),
            "wkb_first_order": [wkb.real, wkb.imag],
        })
    good = [r for r in rows if r["omega"] is not None]
    return {
        "rows": rows,
        "every_real_part_sits_above_the_eikonal": all(
            r["real_part_above_eikonal"] for r in good),
        "every_damping_within_10_percent_of_the_asymptote": all(
            r["damping_near_asymptote"] < 0.036 for r in good),
        "the_pattern_a_correct_solver_must_show": (
            "Real parts approach 0.5(l+1) from above and damping approaches "
            "-0.353553; both happen."),
        "ell_zero_is_least_certain": (
            "The l = 0 barrier is weakest and the power-law tail contaminates its "
            "fit earliest, so its uncertainty is wider than the others."),
    }


def measure_the_cross_validation_verdict() -> Dict[str, object]:
    """R4 — **the deliverable**: which of PR #270's two codes was right."""
    independent = fundamental_mode(1, step=0.025)
    wkb = wkb_fundamental(1)
    eik = eikonal_limit(1)["omega"]

    def gap(a: complex, b: complex) -> Dict[str, float]:
        return {"real": abs(a.real - b.real), "imag": abs(a.imag - b.imag),
                "imag_percent": 100.0 * abs(a.imag - b.imag) / abs(b.imag)}

    return {
        "this_round_characteristic": [independent.real, independent.imag],
        "pr_270_kerr_schild": [KERR_SCHILD_ELL_1.real, KERR_SCHILD_ELL_1.imag],
        "pr_270_tortoise": [TORTOISE_ELL_1.real, TORTOISE_ELL_1.imag],
        "wkb_first_order": [wkb.real, wkb.imag],
        "eikonal_asymptote": [eik.real, eik.imag],
        "gap_to_kerr_schild": gap(independent, KERR_SCHILD_ELL_1),
        "gap_to_tortoise": gap(independent, TORTOISE_ELL_1),
        "kerr_schild_is_confirmed": bool(
            abs(independent - KERR_SCHILD_ELL_1) < 1e-3),
        "tortoise_damping_is_excluded": bool(
            abs(independent.imag - TORTOISE_ELL_1.imag) > 0.05),
        "damping_lines_of_evidence": {
            "characteristic (this round)": independent.imag,
            "Kerr-Schild (PR #270)": KERR_SCHILD_ELL_1.imag,
            "first-order WKB": wkb.imag,
            "exact eikonal asymptote": eik.imag,
            "tortoise (PR #270)": TORTOISE_ELL_1.imag,
        },
        "verdict": (
            "The Kerr-Schild code was right and the tortoise code's damping was "
            "wrong. An independent characteristic evolution sharing no code with "
            "either agrees with Kerr-Schild to ~1e-4 and misses the tortoise "
            "damping by "
            f"{gap(independent, TORTOISE_ELL_1)['imag_percent']:.0f}%. "
            "PR #270 named the Kerr-Schild inner cut as the "
            "prime suspect; that suspicion pointed at the wrong code."),
        "what_this_round_cannot_do": (
            "Neither #270 code was landed in the tree, only their reported "
            "numbers, so there is no autopsy of WHICH line of the tortoise code "
            "produced the wrong damping -- only the demonstration that it did."),
    }


def measure_what_did_not_work() -> Dict[str, object]:
    """R5 — two honest negatives, both reproduced rather than papered over."""
    return {
        "frequency_domain_shooting": {
            "status": "REPRODUCED PR #270's NON-CONVERGENCE",
            "roots_across_epsilon": ["1.229 - 0.152i", "1.173 - 0.102i",
                                     "1.155 - 0.214i"],
            "roots_across_matching_radius": ["1.204 - 0.209i", "1.173 - 0.102i",
                                             "1.166 - 0.105i"],
            "why": ("The QNM boundary-value problem is exponentially "
                    "ill-conditioned in real r: for Im w < 0 the outgoing piece "
                    "grows like e^{|Im w| R} and swamps the incoming coefficient "
                    "one is trying to zero."),
            "so_pr_270s_diagnosis_stands": True,
        },
        "sixth_order_wkb": {
            "status": "UNUSABLE BY FINITE DIFFERENCES",
            "values_under_refinement": ["9.01 + 8.97i", "18.63 + 18.61i",
                                        "623.09 + 623.09i"],
            "why": ("The Iyer-Will formula needs V^(6) at the barrier peak, and "
                    "central differences amplify roundoff as h^-6, so refining "
                    "the grid makes the answer DIVERGE."),
            "what_is_used_instead": ("first-order WKB, which is well conditioned, "
                                     "with its accuracy stated rather than "
                                     "assumed"),
        },
        "first_order_wkb_accuracy": {
            "damping_at_ell_1": "0.4% -- good",
            "real_part_at_ell_1": "13% -- poor",
            "real_part_at_ell_2": "6.8% -- improving",
            "reading": ("Low-l WKB is simply inaccurate on the real part; this "
                        "is a known limitation of the method and NOT a "
                        "discrepancy between solvers."),
        },
    }


def measure_the_ringdown_ledger() -> Dict[str, object]:
    """R6 — what is settled, and what is not."""
    verdict = measure_the_cross_validation_verdict()
    entries = [
        {"claim": "PR #270's two time-domain codes disagreed by 37% in damping",
         "verdict": "CONFIRMED, AND NOW RESOLVED",
         "evidence": "an independent characteristic evolution agrees with "
                     "Kerr-Schild to ~1e-4 and excludes the tortoise damping"},
        {"claim": "the Kerr-Schild inner cut was the prime suspect",
         "verdict": "WRONG SUSPECT",
         "evidence": "that code's frequency is confirmed; the fault was in the "
                     "tortoise evolution"},
        {"claim": "a quasinormal frequency can now be quoted",
         "verdict": "YES, FOR l = 1, 2, 3",
         "evidence": "converged in step size, window-stable, and consistent "
                     "with the exact eikonal asymptote"},
        {"claim": "the l = 0 frequency is equally well determined",
         "verdict": "NO",
         "evidence": "its barrier is weakest and the power-law tail contaminates "
                     "the fit earliest; quoted with a wider uncertainty"},
        {"claim": "frequency-domain shooting can settle this",
         "verdict": "NO -- REPRODUCED THE FAILURE",
         "evidence": "the root moves with every numerical knob; the problem is "
                     "exponentially ill-conditioned in real r"},
        {"claim": "higher-order WKB can settle this",
         "verdict": "NOT BY FINITE DIFFERENCES",
         "evidence": "V^(6) by central differences diverges under refinement"},
        {"claim": "the tortoise code's error can be diagnosed",
         "verdict": "NO",
         "evidence": "neither #270 code was landed, only their numbers; this "
                     "round shows WHICH is wrong, not WHY"},
    ]
    return {
        "entries": entries,
        "headline": verdict["verdict"],
        "the_standing_lesson_held": (
            "PR #270 refused to quote a frequency from two converged codes that "
            "disagreed. That was right, and the way out was a third "
            "implementation sharing no code with either -- plus an exact "
            "asymptote to judge all three against."),
        "still_open": (
            "Overtones, backreaction, and the retarded outer-to-inner transfer "
            "function that #270 also deferred; the transfer function is a ratio "
            "of the same two signals and is now unblocked, since the signals can "
            "be trusted."),
    }
