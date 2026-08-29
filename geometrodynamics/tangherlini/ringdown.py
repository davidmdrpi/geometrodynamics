"""Settling PR #270's ringdown cross-validation, against a published spectrum.

THE VERDICT
───────────
PR #270 built two horizon-penetrating time-domain codes for a test scalar on a
fixed ``D = 5`` Tangherlini background. Both were stable, both converged, and
they disagreed: real parts within ``0.3 %`` at ``ℓ = 1``, **damping rates apart
by 37 % of the smaller value** — equivalently, the wrong one was 27 % low; see
STATE THE DENOMINATOR below. #270 refused to quote a frequency, correctly — *a
converged number is not a correct number* — and named its own prime suspect:
the Kerr–Schild operator's inner cut.

**The Kerr–Schild code was right. The tortoise code's damping was wrong.**

Two things establish it. An independent Gundlach–Price–Pullin characteristic
evolution, written here from scratch and sharing no code with either; and a
**published high-precision spectrum** external to this repository entirely
(Matyjasek 2021, continued fractions cross-checked against Hill determinants,
agreeing to 11 digits). At ``ℓ = 1``, with ``r_h = 1``:

    published (external)   1.01601691 − 0.36232802i    the reference
    #270 Kerr–Schild       1.01622    − 0.36231i       damping  0.005 % off
    this round (h = 0.025) 1.01612    − 0.36244i       damping  0.031 % off
    this round (h = 0.05)  1.01618    − 0.36240i       damping  0.019 % off
    #270 tortoise          1.01876    − 0.26404i       damping 27.1   % off

**No ``h → 0`` value is quoted, deliberately.** The two rows above are raw
values at two step sizes, not an extrapolated limit. The sequence is not
monotone in its distance to the published value, so Richardson extrapolation is
not justified here and none is performed; an earlier draft labelled the
``h = 0.025`` row "``h → 0``", which it is not.

So this round is not merely a tie-breaker between two internal codes. It is an
independent implementation reproducing a known spectrum, which is a much
stronger check on the corrected radial operator and on the GPP machinery than
internal arbitration could ever be.

**Note the ordering.** #270's Kerr–Schild code is roughly six times *more*
accurate than this round's characteristic evolution. The characteristic scheme
arbitrated correctly, but it is not the most accurate of the three, and saying
otherwise would misread what it was for.

WHAT THE EXTERNAL REFERENCE EXPOSED ABOUT THIS SOLVER
─────────────────────────────────────────────────────
Refining the step gave ``−0.3621571 → −0.3623949 → −0.3624352``, a last
successive difference of ``4.0e-5``. The actual distance from the finest value
to the published one is ``1.07e-4`` — **2.7× larger**. The ``h = 0.05`` value is
in fact *closer* to the truth than the ``h = 0.025`` one.

What that establishes, precisely: **self-convergence estimates the error
component associated with the refinement parameter, not the total error.** There
is a component of order ``1e-4`` that step refinement does not control.

Where that component lives is a separate question, and it is now *measured*
rather than asserted — see ``measure_the_extraction_systematics``, which varies
the extraction window, observer radius and ``t_max``. The window dominates by
orders of magnitude (late windows admit the power-law tail); ``t_max`` is
bit-irrelevant; and the band over reasonable choices is comparable to the gap to
the published value. An earlier draft asserted this without varying anything,
and also claimed "nothing internal to a solver can reveal this" — which is too
strong, since that function is internal and reveals it. What the external
reference uniquely supplies is the **anchor**: a spread says the answer is
uncertain, but only an external value says which point in the spread is right.

WHY THIS BACKGROUND MAKES THE ANSWER CHECKABLE
──────────────────────────────────────────────
``D = 5`` is unusually clean, and two exact facts do most of the work.

**The tortoise correction is a decaying power, not a log.** Exactly,

    r* − r = ½ln((r−1)/(r+1)) = −artanh(1/r)
           = −1/r − 1/(3r³) − 1/(5r⁵) − ⋯

so ``−1/r`` is the leading asymptotic behaviour, *not* an exact equality — the
next term is ``−1/(3r³)`` and the tests check its coefficient. What matters
physically is that every term decays: unlike 4D's growing ``2M ln r`` there is
no Coulomb-like phase, and the far field is a pure Hankel wave.

**The potential's tail is asymptotically Bessel.** Exactly,

    V = L/r² + (9/4 − L)/r⁴ − (9/4)/r⁶ ,     L = (ℓ+1)² − ¼

(verified symbolically), so ``V → L/r²`` only *asymptotically* and
``√r H⁽¹⁾_{ℓ+1}(ωr)`` is the **leading** outgoing solution, not the exact one at
finite ``r`` — the exact solution carries a further radial series. An earlier
draft said "exactly Bessel", which the ``1/r⁴`` and ``1/r⁶`` terms contradict.
This is the same identity that fixed the operator correction in PR #271, reused
here as an asymptotic boundary condition.

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

FIVE INDEPENDENT LINES, ALL AGAINST THE TORTOISE DAMPING
─────────────────────────────────────────────────────────
    published (external)       −0.36233    ← the reference
    characteristic evolution   −0.36244    (this round, h = 0.025)
    #270 Kerr–Schild           −0.36231
    first-order WKB            −0.36095
    exact eikonal asymptote    −0.35355
    ────────────────────────────────────
    #270 tortoise              −0.26404    ← excluded

STATE THE DENOMINATOR
─────────────────────
"The tortoise damping is X % off" is ambiguous until the denominator is named,
and this module previously quoted the two conventions in different places. Both,
once, explicitly:

    |Δ| / |published|  = 27.1 %   the tortoise result's relative error
    |Δ| / |tortoise|   = 37.3 %   how much larger the correct damping is

``27.1 %`` is the conventional relative error against truth and is the number to
quote. ``37.3 %`` is what #270 measured, because it was comparing its two codes
to each other with no reference available — which is exactly the situation this
round removed. Both are reported by ``measure_the_cross_validation_verdict``.

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
**But it is not shown to be the sole cause.** The outer condition was truncated
to *pure* Hankel at finite ``R``, and since ``V`` carries ``1/r⁴`` and ``1/r⁶``
terms the exact outgoing solution has a further series — so some of the
matching-radius drift is boundary truncation, not conditioning. This round did
not separate them. Imposing the asymptotic series to several orders and
re-scanning ``R`` would: conditioning-driven drift persists, truncation-driven
drift shrinks with series order.

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
    "PUBLISHED_FUNDAMENTAL",
    "PUBLISHED_SOURCE",
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
    "measure_against_the_published_spectrum",
    "measure_the_extraction_systematics",
    "measure_what_did_not_work",
    "measure_the_ringdown_ledger",
]

HORIZON = 1.0

#: #270 reported these and refused to choose between them.
KERR_SCHILD_ELL_1 = complex(1.01622, -0.36231)
TORTOISE_ELL_1 = complex(1.01876, -0.26404)

#: Published fundamental scalar QNMs of the D = 5 Schwarzschild-Tangherlini hole,
#: **external to this repository**: J. Matyjasek, "Accurate quasinormal modes of
#: the five-dimensional Schwarzschild-Tangherlini black holes", Phys. Rev. D 104,
#: 084066 (2021), arXiv:2107.04815. Computed there by continued fractions and
#: cross-checked against Hill determinants, agreeing to 11 digits.
#:
#: The paper tabulates the scaled frequency `w~ = w / T_H` with `T_H = 1/(2 pi)`,
#: so these are `w~ / (2 pi)` at `r_h = 1`. The l = 1 entry is quoted there to
#: full precision (`6.38382253011 - 2.27657411582i`); l = 0 and l = 2 are given
#: to fewer digits, which is why the tolerances below differ by l.
PUBLISHED_FUNDAMENTAL = {
    0: complex(0.53383557, -0.38337537),
    1: complex(1.01601691149, -0.36232802385),
    2: complex(1.51056745177, -0.35753725529),
}

PUBLISHED_SOURCE = ("Matyjasek, Phys. Rev. D 104, 084066 (2021), "
                    "arXiv:2107.04815 -- continued fractions cross-checked "
                    "against Hill determinants, agreeing to 11 digits")


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

    ``V`` at a diamond centre depends only on ``r*``, so it is a
    one-dimensional lookup indexed by ``j − i``. The inner ``v`` recursion is
    linear and first order, ``a_j = A_j a_{j−1} + B_j``, so it is evaluated by
    cumulative products rather than a Python loop.

    **The diamond centre, derived.** For the update ``N = (i, j)`` with
    ``S = (i−1, j−1)``, the centre sits at ``u_c = (i − ½)h`` and
    ``v_c = (j − ½)h``, so

        ``r*_c = (v_c − u_c)/2 = [(j − ½)h − (i − ½)h] / 2 = (j − i)h/2``

    — **the two half-steps cancel**. An earlier version of this function
    sampled ``0.5·step·((j−i) − 0.5)``, i.e. ``r*_c − h/4``, contradicting its
    own docstring. That is fixed here and locked by a test. Measured effect of
    the correction: ``~1e-6`` in the extracted frequency, three orders of
    magnitude below this solver's error floor — a real derivation bug whose
    numerical consequence at these step sizes is negligible, and **not** the
    cause of the floor (see ``measure_the_extraction_systematics``).
    """
    count = int(t_max / step)
    offsets = np.arange(-count, count + 1)
    centres = 0.5 * step * offsets
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
    # The tail is a POWER law, not a log. Exactly,
    #     r* - r = (1/2)ln((r-1)/(r+1)) = -artanh(1/r)
    #            = -1/r - 1/(3r^3) - 1/(5r^5) - ...
    # so -1/r is the LEADING BEHAVIOUR, not an equality. Testing (r* - r)*r -> -1
    # is the real statement; an arbitrary threshold on |r* - r| is not, and the
    # first draft's 1e-3 cut sat exactly on the value at r = 1000. Better still,
    # the series predicts the next coefficient, so the deviation is checked
    # against -1/(3r^2) rather than merely being required to shrink.
    scaled = [(tortoise(r) - r) * r for r in far]
    predicted = [-1.0 / (3.0 * r ** 2) for r in far]
    coefficient_errors = [abs((s + 1.0) - p) / abs(p)
                          for s, p in zip(scaled, predicted)]
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
        "predicted_next_term_minus_one_over_three_r_squared": predicted,
        "next_coefficient_relative_errors": coefficient_errors,
        "the_leading_tail_is_minus_one_over_r": bool(
            abs(scaled[-1] + 1.0) < 1e-5
            and all(abs(b + 1.0) < abs(a + 1.0)
                    for a, b in zip(scaled[:-1], scaled[1:]))),
        "the_next_series_coefficient_is_confirmed": bool(
            all(e < 1e-3 for e in coefficient_errors)),
        "the_exact_closed_form": (
            "r* - r = -artanh(1/r) = -1/r - 1/(3r^3) - 1/(5r^5) - ... , so -1/r "
            "is the leading asymptotic behaviour and NOT an exact equality. The "
            "next term is checked against its predicted coefficient."),
        "no_logarithmic_tail": bool(
            abs(scaled[-1] + 1.0) < 1e-5
            and all(abs(b + 1.0) < abs(a + 1.0)
                    for a, b in zip(scaled[:-1], scaled[1:]))),
        "why_that_matters": ("every term in the series decays, so unlike 4D's "
                             "growing 2M ln r there is no Coulomb-like phase and "
                             "the far field is a pure Hankel wave"),
        "bessel_tail": bessel,
        "the_tail_is_asymptotically_bessel": bool(
            all(b["relative_error_at_1000"] < 1e-5 for b in bessel)),
        "the_exact_potential_is_not_pure_inverse_square": (
            "V = L/r^2 + (9/4 - L)/r^4 - (9/4)/r^6 with L = (l+1)^2 - 1/4, "
            "verified symbolically. So V -> L/r^2 only ASYMPTOTICALLY, and "
            "sqrt(r) H^(1)_{l+1}(w r) is the leading outgoing solution rather "
            "than the exact one at finite r -- the exact outgoing solution "
            "carries a further radial series. An earlier draft said 'exactly "
            "Bessel', which the 1/r^4 and 1/r^6 terms contradict."),
        "the_same_identity_that_fixed_the_operator": (
            "V -> [(l+1)^2 - 1/4]/r^2 is the flat-limit identity PR #271 used to "
            "settle which radial operator was correct; it is reused here as an "
            "asymptotic boundary condition."),
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
    step = 0.025
    independent = fundamental_mode(1, step=step)
    wkb = wkb_fundamental(1)
    eik = eikonal_limit(1)["omega"]

    reference = PUBLISHED_FUNDAMENTAL[1]

    def gap(a: complex, b: complex) -> Dict[str, float]:
        """Both denominators, named. `b` is the thing `a` is being compared to.

        A bare "X % off" is ambiguous, and this module used to quote the two
        conventions in different places. `percent_of_reference` divides by the
        published value and is the conventional relative error; the other
        divides by `b` and is what #270 could compute with no reference in hand.
        """
        delta = abs(a.imag - b.imag)
        return {
            "real": abs(a.real - b.real),
            "imag": delta,
            "imag_percent_of_reference": 100.0 * delta / abs(reference.imag),
            "imag_percent_of_comparison": 100.0 * delta / abs(b.imag),
        }

    return {
        "this_round_characteristic": [independent.real, independent.imag],
        "this_round_step_size": step,
        "this_round_is_a_raw_value_not_an_extrapolated_limit": (
            "Quoted at h = 0.025. This is the finest raw point, NOT an h -> 0 "
            "limit: the sequence is non-monotone in its distance to the "
            "published value, so Richardson extrapolation is not justified and "
            "none is performed. R7 reports h = 0.05 separately; the two differ "
            "because they are different step sizes, not because the solver has "
            "two answers."),
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
        "published_reference": [reference.real, reference.imag],
        "published_source": PUBLISHED_SOURCE,
        "gap_to_published": gap(independent, reference),
        "damping_lines_of_evidence": {
            "published (external reference)": reference.imag,
            f"characteristic (this round, h = {step})": independent.imag,
            "Kerr-Schild (PR #270)": KERR_SCHILD_ELL_1.imag,
            "first-order WKB": wkb.imag,
            "exact eikonal asymptote": eik.imag,
            "tortoise (PR #270)": TORTOISE_ELL_1.imag,
        },
        "the_denominator_is_named": {
            "tortoise_relative_error_against_published": (
                100.0 * abs(TORTOISE_ELL_1.imag - reference.imag)
                / abs(reference.imag)),
            "correct_damping_is_larger_than_tortoise_by": (
                100.0 * abs(TORTOISE_ELL_1.imag - reference.imag)
                / abs(TORTOISE_ELL_1.imag)),
            "which_to_quote": (
                "The relative error against the published value (~27.1%) is the "
                "conventional statement. The ~37.3% figure divides by the "
                "tortoise value instead -- it says the correct damping is that "
                "much LARGER -- and is what PR #270 measured because it had two "
                "codes and no reference. Both are true; neither should be "
                "quoted without naming its denominator."),
        },
        "verdict": (
            "The Kerr-Schild code was right and the tortoise code's damping was "
            "wrong. A published high-precision spectrum external to this "
            "repository puts the answer at "
            f"{reference.real:.8f}{reference.imag:+.8f}i, which confirms "
            "Kerr-Schild to "
            f"{gap(KERR_SCHILD_ELL_1, reference)['imag_percent_of_reference']:.3f}% "
            "and this round's independent characteristic evolution to "
            f"{gap(independent, reference)['imag_percent_of_reference']:.3f}%, "
            "while the tortoise damping is off by "
            f"{gap(TORTOISE_ELL_1, reference)['imag_percent_of_reference']:.1f}%. "
            "PR #270 named the Kerr-Schild inner cut as the prime suspect; that "
            "suspicion pointed at the wrong code."),
        "what_this_round_cannot_do": (
            "Neither #270 code was landed in the tree, only their reported "
            "numbers, so there is no autopsy of WHICH line of the tortoise code "
            "produced the wrong damping -- only the demonstration that it did."),
    }


def measure_against_the_published_spectrum(
        steps: Sequence[float] = (0.1, 0.05, 0.025),
        step: float = 0.05) -> Dict[str, object]:
    """R7 — the external check, and what it says about this solver's error bar.

    Two separate things come out of comparing against a published spectrum,
    and only the first is the one that was expected.

    **The spectrum is reproduced.** Three modes, independently computed here,
    land on values obtained by continued fractions and Hill determinants.

    **The step-size study was over-optimistic.** Refining `h` produced a last
    successive difference of ~`4e-5`, while the finest value actually sits
    ~`1.1e-4` from the published one — and the middle step is *closer* than the
    finest. So discretization is not the limiting error; extraction systematics
    are. Self-convergence measures only the error being refined away. It is a
    consistency check, **not an error bar** — which is exactly the failure mode
    #270 warned about, in a form only an external reference could reveal.
    """
    rows = []
    for ell, published in sorted(PUBLISHED_FUNDAMENTAL.items()):
        got = fundamental_mode(ell, step=step)
        if got is None:
            rows.append({"ell": ell, "characteristic": None})
            continue
        rows.append({
            "ell": ell,
            "characteristic": [got.real, got.imag],
            "published": [published.real, published.imag],
            "real_relative_error": abs(got.real - published.real) / abs(published.real),
            "damping_relative_error": abs(got.imag - published.imag) / abs(published.imag),
        })
    good = [r for r in rows if r["characteristic"] is not None]

    # What refinement claimed, against what the reference says is true.
    published_1 = PUBLISHED_FUNDAMENTAL[1].imag
    sequence = [fundamental_mode(1, step=s) for s in steps]
    damping = [None if w is None else w.imag for w in sequence]
    deltas = [abs(b - a) for a, b in zip(damping[:-1], damping[1:])
              if a is not None and b is not None]
    distances = [None if d is None else abs(d - published_1) for d in damping]
    finest_error = distances[-1]
    last_delta = deltas[-1] if deltas else None

    return {
        "source": PUBLISHED_SOURCE,
        "step_size_of_the_rows": step,
        "rows": rows,
        "every_mode_within_0p3_percent": bool(
            all(r["damping_relative_error"] < 3e-3
                and r["real_relative_error"] < 3e-3 for r in good)),
        "ell_1_and_2_within_0p05_percent": bool(
            all(r["damping_relative_error"] < 5e-4
                and r["real_relative_error"] < 5e-4
                for r in good if r["ell"] in (1, 2))),
        "ell_0_is_the_loosest": bool(
            max(good, key=lambda r: r["damping_relative_error"])["ell"] == 0),
        "refinement_versus_truth": {
            "steps": list(steps),
            "damping_by_step": damping,
            "distance_to_published_by_step": distances,
            "last_successive_difference": last_delta,
            "distance_from_finest_to_published": finest_error,
            "understatement_factor": (None if not last_delta else
                                      finest_error / last_delta),
            "the_finest_step_is_not_the_closest": bool(
                any(d is not None and d < finest_error for d in distances[:-1])),
        },
        "the_lesson": (
            "Self-convergence measures only the error it is refining away. The "
            "step-size study's last successive difference was ~2.7x smaller "
            "than the finest value's actual distance to the published one, and "
            "the middle step lands closer than the finest -- so discretization "
            "is not what limits this solver, extraction systematics are. A "
            "convergence study is a consistency check, not an error bar."),
        "the_reframing": (
            "With an external reference in hand this is no longer only a "
            "tie-breaker between two internal codes. It is an independent "
            "implementation reproducing a known high-precision spectrum, which "
            "is a considerably stronger check on PR #271's corrected radial "
            "operator and on the GPP machinery than internal arbitration."),
        "who_is_most_accurate": (
            "PR #270's Kerr-Schild code, at 0.005% in damping against 0.031% "
            "for this round's characteristic evolution -- about 6x better. The "
            "characteristic scheme arbitrated correctly; it is not the most "
            "accurate of the three, and should not be described as though it "
            "were."),
    }


def measure_the_extraction_systematics() -> Dict[str, object]:
    """R8 — where this solver's error floor actually lives, by varying knobs.

    An earlier draft asserted that the residual error was extraction
    systematics. That was a hypothesis with nothing varied behind it. This
    measurement varies them, and the assertion survives — but it is now a
    measurement, and two of its corollaries did not survive.

    **The extraction window dominates**, by orders of magnitude: late windows
    are contaminated by the power-law tail and degrade without limit. Prony
    order and observer radius matter at the `1e-3` level. `t_max` does not
    matter at all (bit-identical), because the extraction window sits well
    inside it.

    **What this corrects.** An internal study *can* expose a floor that step
    refinement hides — this function is that study, and it needed no external
    value. What the published spectrum uniquely supplies is the *anchor*: the
    spread tells you the answer is uncertain, but only an external reference
    tells you which point in the spread is right.
    """
    published = PUBLISHED_FUNDAMENTAL[1]

    def error_of(**kwargs) -> Optional[float]:
        w = fundamental_mode(1, step=0.05, **kwargs)
        return None if w is None else abs(w.imag - published.imag) / abs(published.imag)

    windows = [(50.0, 130.0), (60.0, 140.0), (70.0, 150.0),
               (80.0, 160.0), (90.0, 180.0)]
    by_window = [{"window": list(w), "damping_relative_error": error_of(window=w)}
                 for w in windows]
    by_observer = [{"observer_rstar": o,
                    "damping_relative_error": error_of(observer_rstar=o)}
                   for o in (20.0, 30.0, 40.0)]
    by_t_max = [{"t_max": t, "damping_relative_error": error_of(t_max=t)}
                for t in (300.0, 400.0, 500.0)]

    live = [d["damping_relative_error"] for d in by_window + by_observer
            if d["damping_relative_error"] is not None]
    # "Reasonable" = the settings a careful person would actually pick: a window
    # that holds ringdown without running into the tail, and an observer the
    # signal has cleanly reached. Not the best-case cherry-pick.
    reasonable = [error_of(window=(60.0, 140.0)), error_of(window=(70.0, 150.0)),
                  error_of(observer_rstar=20.0)]
    reasonable = [x for x in reasonable if x is not None]

    t_max_values = [d["damping_relative_error"] for d in by_t_max
                    if d["damping_relative_error"] is not None]

    # The step-refinement difference, in the same units, to compare against.
    damping_by_step = [fundamental_mode(1, step=st) for st in (0.05, 0.025)]
    step_difference = (
        abs(damping_by_step[1].imag - damping_by_step[0].imag) / abs(published.imag)
        if all(w is not None for w in damping_by_step) else None)

    return {
        "by_extraction_window": by_window,
        "by_observer_radius": by_observer,
        "by_t_max": by_t_max,
        "worst_over_all_variations": max(live) if live else None,
        "band_over_reasonable_choices": (
            {"min": min(reasonable), "max": max(reasonable)} if reasonable else None),
        "the_window_dominates": bool(
            live and max(live) > 100.0 * min(live)),
        "t_max_is_irrelevant": bool(
            len(set(f"{v:.12e}" for v in t_max_values)) == 1),
        # The meaningful comparison, not an arbitrary cut: does varying the
        # extraction move the answer by more than refining the step does?
        "step_refinement_last_relative_difference": step_difference,
        "extraction_band_exceeds_step_refinement": bool(
            reasonable and step_difference is not None
            and max(reasonable) > step_difference),
        "how_many_times_larger": (
            None if not reasonable or not step_difference
            else max(reasonable) / step_difference),
        "what_this_establishes": (
            "The residual is dominated by extraction choices -- above all the "
            "window, which degrades without limit once the power-law tail "
            "enters it. Over reasonable choices the extraction band is several "
            "times the step-refinement difference, and comparable to the gap "
            "between this solver and the published value -- which is why step "
            "refinement alone was the wrong error bar."),
        "what_this_corrects": (
            "An internal study CAN expose this floor -- this function is one, "
            "and it used no external value. The earlier claim that 'nothing "
            "internal to a solver can reveal this' was too strong. What the "
            "published spectrum uniquely provides is the anchor: a systematic "
            "spread says the answer is uncertain, but only an external "
            "reference says which point in the spread is correct."),
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
            "a_second_confounder_not_separated": (
                "The outer condition was truncated to PURE Hankel at finite R. "
                "Since V carries 1/r^4 and 1/r^6 terms, the exact outgoing "
                "solution has a further radial series, so part of the observed "
                "matching-radius drift is boundary truncation rather than "
                "conditioning. This round did not separate the two, so the "
                "ill-conditioning above is a CONTRIBUTING cause and not a "
                "demonstrated sole cause."),
            "what_would_separate_them": (
                "Impose the asymptotic series to several orders and re-scan the "
                "matching radius: conditioning-driven drift would persist, "
                "truncation-driven drift would shrink with series order."),
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
        {"claim": "PR #270's two time-domain codes disagreed in damping "
                  "(37% of the smaller value; the wrong one 27% low)",
         "verdict": "CONFIRMED, AND NOW RESOLVED",
         "evidence": "an independent characteristic evolution agrees with "
                     "Kerr-Schild to ~1e-4 and excludes the tortoise damping"},
        {"claim": "the Kerr-Schild inner cut was the prime suspect",
         "verdict": "WRONG SUSPECT",
         "evidence": "that code's frequency is confirmed; the fault was in the "
                     "tortoise evolution"},
        {"claim": "the verdict rests only on internal arbitration",
         "verdict": "NO -- CONFIRMED EXTERNALLY",
         "evidence": "a published high-precision spectrum (continued fractions "
                     "+ Hill determinants) puts l = 1 at 1.01601691-0.36232802i, "
                     "confirming Kerr-Schild to 0.005% and excluding the "
                     "tortoise damping at 27.1%"},
        {"claim": "the characteristic evolution is the most accurate of the three",
         "verdict": "NO",
         "evidence": "PR #270's Kerr-Schild code is ~6x closer to the published "
                     "value (0.005% against 0.031%); this round arbitrated, it "
                     "did not out-resolve"},
        {"claim": "the step-size study bounded this solver's error",
         "verdict": "NO -- IT UNDERSTATED IT 2.7x",
         "evidence": "last successive difference 4.0e-5, actual distance to the "
                     "published value 1.1e-4, and h = 0.05 lands closer than "
                     "h = 0.025; self-convergence estimates only the error "
                     "component tied to the refinement parameter"},
        {"claim": "the residual is extraction systematics",
         "verdict": "MEASURED, NOT ASSERTED",
         "evidence": "varying the extraction window, observer radius and t_max "
                     "shows the window dominating by orders of magnitude and "
                     "t_max being bit-irrelevant; the band over reasonable "
                     "choices is comparable to the gap to the published value"},
        {"claim": "only an external reference could expose that floor",
         "verdict": "NO -- TOO STRONG",
         "evidence": "the internal systematics scan exposes it using no "
                     "external value; what the published spectrum uniquely "
                     "supplies is the anchor, i.e. WHICH point in the spread "
                     "is right"},
        {"claim": "the GPP potential was sampled at the diamond centre",
         "verdict": "NO -- IT WAS OFF BY h/4",
         "evidence": "the centre is r* = (j-i)h/2, the two half-steps "
                     "cancelling; the code sampled r* - h/4, contradicting its "
                     "own docstring. Fixed and locked by a test; measured "
                     "effect ~1e-6, so a real bug that does NOT explain the "
                     "error floor"},
        {"claim": "the far-field solution is exactly Hankel",
         "verdict": "NO -- ASYMPTOTICALLY",
         "evidence": "V = L/r^2 + (9/4-L)/r^4 - (9/4)/r^6 exactly, so the "
                     "outgoing solution carries a further radial series; the "
                     "failed shoot truncated to pure Hankel, making boundary "
                     "truncation a second unseparated confounder"},
        {"claim": "a quasinormal frequency can now be quoted",
         "verdict": "YES, FOR l = 1, 2, 3",
         "evidence": "converged in step size, window-stable, consistent with the "
                     "exact eikonal asymptote, and matching a published spectrum "
                     "to <0.05% at l = 1 and 2"},
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
        "the_lesson_this_round_adds": (
            "Self-convergence estimates the error component associated with the "
            "refinement parameter, not the total numerical, model or extraction "
            "error. This round's step-size study would have quoted ~4e-5 when "
            "the true error was ~1.1e-4. An internal scan over extraction "
            "choices does expose the missing component -- what the external "
            "reference uniquely supplies is the anchor that says which point in "
            "the spread is correct. Look for a published benchmark BEFORE "
            "building a third implementation to break a tie."),
        "still_open": (
            "Overtones, backreaction, and the retarded outer-to-inner transfer "
            "function that #270 also deferred; the transfer function is a ratio "
            "of the same two signals and is now unblocked, since the signals can "
            "be trusted."),
        "the_next_object": (
            "The retarded transfer function G_l(t; r_obs, r_src) from the same "
            "characteristic evolution with a compact ingoing excitation, gated "
            "on three checks before any physical reading: causal support "
            "G(t < t_null) = 0; flux conservation |R_l|^2 + |T_l|^2 = 1 on the "
            "fixed background; and late-time ringdown consistent with the "
            "EXTERNAL value 1.01601691149-0.36232802385i rather than with this "
            "solver's own fitted number."),
    }
