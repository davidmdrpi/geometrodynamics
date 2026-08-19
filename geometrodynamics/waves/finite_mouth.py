"""A **finite-radius** mouth — and the answer to PR #260's question.

The question
────────────
PR #260 built the conservative finite throat the arc had owed since #253: a tube
with an exact Dirichlet-to-Neumann map, a reflected channel, and a real traversal
delay.  It also found that with **point** mouths the composite carries an
exponentially growing mode for *every* choice of parameters, generated at the
point-mouth/tube **interface** — and localizing, in the ``σL, σd ≫ 1`` limit, to
a single mouth at a rate ``σ* = 2√(π/𝒜)`` that knows neither the tube's length
nor the mouths' separation.  That gated the roadmap, with one question to settle
first:

    **does the negative mode survive a finite-radius mouth?**

The answer is **no**, and the reason is sharper than "it goes away".

What changes
────────────
A point mouth couples the tube to the ambient at a *point*, which requires
subtracting the ``1/(4πχ)`` divergence and leaves the **renormalized** self-energy
``g(λ) = −(λ^{1/2}/4π)cot(πλ^{1/2})`` — negative at ``λ = 0`` and growing like
``−κ/4π`` down the imaginary axis.  A finite mouth does not need the subtraction.
Smearing the coupling over a sphere of radius ``a`` — the same operator on both
sides, so the composite stays manifestly self-adjoint — replaces ``g`` by the
**unsubtracted** Green function at distance ``a``:

    ``𝒢_self(λ)  = f(a,λ)·G(a,λ)`` ,     ``𝒢_cross(λ) = f(a,λ)²·G(d,λ)``

with ``f(χ,λ) = sin(ωχ)/(ω sin χ)`` the regular radial solution, ``f(0) = 1``.
Both are **mean-value identities**, not approximations, and both are checked
against direct quadrature on ``S³`` rather than asserted — the cross one to
``2e-11``.

The answer, and the mechanism
─────────────────────────────
At ``λ < 0`` the two sides have **opposite signs, structurally**:

* the tube's channel functions are ``−coth(κL/2)/(𝒜κ)`` and ``−tanh(κL/2)/(𝒜κ)``
  — strictly **negative**;
* the ambient's are ``f·G(a) ± f²·G(d)`` — strictly **positive**, since ``G`` and
  ``f`` are positive on the imaginary axis and ``G(a) > f(a)G(d)`` whenever the
  mouths are smaller than half their separation, which they must be to be
  disjoint.

So ``det(A − 𝒢) ≠ 0`` on the whole negative-``λ`` axis, for **every** parameter
choice.  There is no growing mode, and the statement is structural rather than
parametric — no sweep can find one, and a sweep of 3078 points does not.

**What produced #260's mode was the linearization.**  That round modelled the
mouth by adding a *constant* ``1/(4πa)`` to the self-energy — the leading term of
``G(a,λ) = 1/(4πa) + g(λ) + O(a)``.  The exact ``G(a,λ)`` is **screened**: it
decays like ``e^{−κa}`` down the imaginary axis, while the constant does not, so
the linearized mouth eventually beats the tube's ``−1/(𝒜κ)`` and crosses zero.
It crosses at ``κ ≈ 1/a`` — precisely where ``κa ≈ 1`` and the linearization is
invalid.  The two agree to four digits for ``κa ≪ 1`` and disagree in **sign**
beyond it.

Where the mode went
───────────────────
It did not vanish; it became **soft and positive**.  The composite has exactly
one state below the free ESU gap, in the symmetric channel, and as ``a → 0``

    ``λ₀  ⟶  8πa/(𝒜L)``

— two mouth capacitances ``4πa`` restoring a tube of volume ``𝒜L`` — verified to
``0.2%`` at ``a = 0.005``.  So the point limit drives this mode to zero **from
above**, and the linearized mouth overshoots it into the unstable half-plane.
That is the whole story: not a small error in a rate, but a sign error produced
by freezing a screened quantity.

What survives from PR #260
──────────────────────────
The traversal delay survives cleanly: slope ``1.0010`` in ``L``, saturating at
the ambient path ``d`` to ``0.0``, with the mouth contributing only a
**sub-leading** ``O(a)`` shift.  The static rank-one collapse and PR #258's
``𝒲 = −β(λ)`` survive unchanged, because both came from the *tube*'s zero mode
and the mouth does not touch it.

And the contour is easier now: with no growing mode to clear, ``ε`` is back to
PR #259's single requirement of sitting above the frequency spacing.  PR #260
needed ``ε > σ* ≈ 2``; ``0.4`` is comfortable here.

What is still put in
────────────────────
The coupling is **monopole**: one channel per mouth, so only the ``ℓ = 0``
projection on each sphere talks to the tube.  Higher multipoles are dropped, and
PR #250's screening law says what that costs — the ``ℓ``-th one enters at
``O((a/d)^ℓ)``.  The mouths are spheres in a *fixed* ambient, not a solved neck
geometry.  **No backreaction.**
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .finite_throat import FOUR_PI, FiniteThroat, causal_onset
from .two_wave import GaussianPulse, RetardedGrid, green_omega

__all__ = [
    "FiniteMouthThroat",
    "WORKING_MOUTH",
    "regular_radial",
    "mouth_green",
    "screened_products",
    "shell_average_cross",
    "shell_average_self",
    "measure_the_mean_value_identities_hold",
    "measure_the_negative_mode_does_not_survive",
    "measure_the_instability_was_the_linearization",
    "measure_the_mode_became_soft_and_positive",
    "measure_the_delay_survives_with_a_radius_correction",
    "measure_the_static_results_survive",
    "measure_monopole_matching_is_the_remaining_approximation",
]


# ════════════════════════════════════════════════════════════════════════════
# THE AMBIENT SEEN BY A SPHERE
# ════════════════════════════════════════════════════════════════════════════
def regular_radial(chi: float, lmbda: float) -> float:
    """``f(χ,λ) = sin(ωχ)/(ω sin χ)`` — the radial solution regular at ``χ = 0``.

    Normalized to ``f(0) = 1``, so it is the **form factor** of a uniform shell:
    a shell of radius ``a`` carrying unit charge produces the exterior field
    ``f(a,λ)·G(χ,λ)``, and by the mean-value theorem the average of any regular
    solution over that shell is ``f(a,λ)`` times its value at the centre.  Both
    directions of the smearing therefore cost one factor of ``f``.

    Continued to ``λ < 0`` as ``sinh(κχ)/(κ sin χ)``, positive throughout — which
    is half of why the composite has no growing mode.
    """
    c = float(chi)
    lam = float(lmbda)
    if c < 1e-12:
        return 1.0
    s = math.sin(c)
    if lam > 0.0:
        w = math.sqrt(lam)
        return math.sin(w * c) / (w * s)
    if lam < 0.0:
        k = math.sqrt(-lam)
        return math.sinh(k * c) / (k * s)
    return c / s


def mouth_green(chi: float, lmbda: float) -> float:
    """``G(χ,λ)`` — the **unsubtracted** Green function, real ``λ`` either sign.

    This is the whole difference from the point mouth.  A point interaction has
    to subtract the ``1/(4πχ)`` divergence and keep the renormalized remainder
    ``g(λ)``, which is *negative*; a sphere of radius ``a`` simply reads ``G`` at
    ``χ = a``, which is *positive* and, at ``λ = −κ²``, **screened**:

        ``G(χ,−κ²)  =  (e^{−κχ} − e^{−κ(2π−χ)}) / (4π sin χ (1 − e^{−2πκ}))``

    written in decaying exponentials so it neither overflows nor cancels.  The
    two-path form is visible in it: ``e^{−κχ}`` the short way round the sphere and
    ``e^{−κ(2π−χ)}`` the long way, PR #253's ledger with the oscillation replaced
    by decay.
    """
    c = float(chi)
    lam = float(lmbda)
    if not 0.0 < c <= math.pi + 1e-12:
        raise ValueError("the mouth needs 0 < χ ≤ π")
    s = math.sin(c)
    if abs(lam) < 1e-14:
        return (math.pi - c) / (4.0 * math.pi ** 2 * s)
    if lam < 0.0:
        k = math.sqrt(-lam)
        decay = math.exp(-2.0 * math.pi * k)
        return ((math.exp(-k * c) - math.exp(-k * (2.0 * math.pi - c)))
                / (4.0 * math.pi * s * (1.0 - decay)))
    w = math.sqrt(lam)
    sp = math.sin(math.pi * w)
    if abs(sp) < 1e-13:
        raise ValueError("λ is a free eigenvalue; G has a pole there")
    return math.sin(w * (math.pi - c)) / (4.0 * math.pi * s * sp)


def screened_products(radius: float, separation: float, sigma: float
                      ) -> Tuple[float, float]:
    """``(f·G(a), f²·G(d))`` at ``λ = −σ²``, in a form that cannot overflow —
    **and in which their positivity is manifest.**

    ``f`` grows like ``e^{κa}`` and ``G`` decays like ``e^{−κa}``, so computing
    them apart overflows above ``κa ≈ 700`` while their product is bounded.
    Combining first:

        ``f·G(a)   = (1 − e^{−2κa})(1 − e^{−2κ(π−a)}) / (8πκ sin²a (1 − e^{−2πκ}))``
        ``f²·G(d)  = e^{2κa−κd}(1 − e^{−2κa})²(1 − e^{−2κ(π−d)})
                     / (16πκ² sin²a sin d (1 − e^{−2πκ}))``

    Every bracket is in ``(0, 1]`` and ``2κa − κd < 0`` because disjoint mouths
    need ``a < d/2``.  So **both are positive for every ``κ``**, by inspection
    rather than by scan — which is the ambient half of this round's answer, and
    the reason no parameter choice can produce a growing mode.
    """
    a, d, s = float(radius), float(separation), float(sigma)
    if s < 1e-12:
        f = regular_radial(a, 0.0)
        return (f * mouth_green(a, 0.0), f * f * mouth_green(d, 0.0))
    sa, sd = math.sin(a), math.sin(d)
    tail = 1.0 - math.exp(-2.0 * math.pi * s)
    self_term = ((1.0 - math.exp(-2.0 * s * a))
                 * (1.0 - math.exp(-2.0 * s * (math.pi - a)))
                 / (8.0 * math.pi * s * sa * sa * tail))
    cross = (math.exp(2.0 * s * a - s * d)
             * (1.0 - math.exp(-2.0 * s * a)) ** 2
             * (1.0 - math.exp(-2.0 * s * (math.pi - d)))
             / (16.0 * math.pi * s * s * sa * sa * sd * tail))
    return self_term, cross


def shell_average_cross(radius: float, separation: float, lmbda: float,
                        n: int = 200001) -> float:
    """``⟨G(·,c₂)⟩`` over the sphere at distance ``a`` from ``c₁``, by quadrature.

    The locus at fixed distance from a point of ``S³`` is a round ``S²``; the
    spherical law of cosines gives the distance to the other mouth,
    ``cos χ₂ = cos a cos d + sin a sin d cos θ``, and the measure in the polar
    angle toward ``c₂`` is ``½ sin θ dθ``.  The mean-value theorem predicts
    ``f(a,λ)·G(d,λ)``; this computes it instead, so the model's central identity
    is checked rather than assumed.
    """
    a, d = float(radius), float(separation)
    th = np.linspace(0.0, math.pi, int(n))
    cos2 = math.cos(a) * math.cos(d) + math.sin(a) * math.sin(d) * np.cos(th)
    chi2 = np.arccos(np.clip(cos2, -1.0, 1.0))
    vals = np.array([mouth_green(float(c), lmbda) for c in chi2])
    return 0.5 * float(np.trapezoid(vals * np.sin(th), th))


def shell_average_self(radius: float, lmbda: float, n: int = 4001) -> float:
    """``⟨⟨G⟩⟩`` over two copies of the same sphere — predicted ``f(a)·G(a)``.

    Both points sit at distance ``a`` from the centre, so
    ``cos χ = cos²a + sin²a cos θ``.  The integrand has the Green function's
    integrable singularity at ``θ = 0``, so this converges like the grid rather
    than to machine precision, and the residual is quoted as quadrature error
    rather than treated as a model error.
    """
    a = float(radius)
    th = np.linspace(0.0, math.pi, int(n))
    cosx = math.cos(a) ** 2 + math.sin(a) ** 2 * np.cos(th)
    chi = np.arccos(np.clip(cosx, -1.0, 1.0))
    keep = chi > 1e-9
    vals = np.array([mouth_green(float(c), lmbda) for c in chi[keep]])
    return 0.5 * float(np.trapezoid(vals * np.sin(th[keep]), th[keep]))


# ════════════════════════════════════════════════════════════════════════════
# THE THROAT WITH RESOLVED MOUTHS
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class FiniteMouthThroat:
    """PR #260's tube, with its point mouths replaced by spheres of radius ``a``.

    The tube block is unchanged and imported rather than re-derived — this round
    changes exactly one thing, the ambient's view of the mouth — so any
    difference in the conclusions is attributable to that change alone.
    """

    separation: float = 1.3
    length: float = 0.9
    area: float = FOUR_PI
    radius: float = 0.05
    interior_mass: float = 0.0

    def __post_init__(self) -> None:
        if not 0.0 < self.radius < 0.5 * self.separation:
            raise ValueError("the mouths must be disjoint: 0 < a < d/2")

    def tube(self) -> FiniteThroat:
        """PR #260's tube, with **no** mouth term — that is this round's job."""
        return FiniteThroat(separation=self.separation, length=self.length,
                            area=self.area, interior_mass=self.interior_mass,
                            mouth_radius=None)

    # -- the ambient block ---------------------------------------------------
    def ambient_matrix(self, lmbda: float) -> np.ndarray:
        """``𝒢(λ)`` — what the two smeared mouths see of each other and of
        themselves.  Diagonal ``f(a)G(a)``, off-diagonal ``f(a)²G(d)``."""
        lam = float(lmbda)
        f = regular_radial(self.radius, lam)
        self_term = f * mouth_green(self.radius, lam)
        cross = f * f * mouth_green(self.separation, lam)
        return np.array([[self_term, cross], [cross, self_term]], dtype=float)

    def ambient_channels(self, lmbda: float) -> Tuple[float, float]:
        g = self.ambient_matrix(lmbda)
        return float(g[0, 0] + g[0, 1]), float(g[0, 0] - g[0, 1])

    def tube_channels(self, lmbda: float) -> Tuple[float, float]:
        """``A_sym = cot(kL/2)/(𝒜k)`` and ``A_anti = −tan(kL/2)/(𝒜k)``."""
        lam = float(lmbda)
        a = self.tube().boundary_matrix(complex(lam))
        return (float((a[0, 0] + a[0, 1]).real),
                float((a[0, 0] - a[0, 1]).real))

    def channel_functions(self, lmbda: float) -> Tuple[float, float]:
        """``(A − 𝒢)`` in both channels; its zeros are the composite spectrum."""
        at, ab = self.tube_channels(lmbda)
        gt, gb = self.ambient_channels(lmbda)
        return at - gt, ab - gb

    def negative_lambda_channels(self, sigma: float) -> Tuple[float, float]:
        """``(A − 𝒢)`` at ``λ = −σ²``, in decaying exponentials throughout.

        The two sides are separately signed here, and that is the round's
        result: ``A`` strictly negative, ``𝒢`` strictly positive, so the
        difference cannot vanish.  `signed_parts` returns them apart.
        """
        a_part, g_part = self.signed_parts(sigma)
        return (a_part[0] - g_part[0], a_part[1] - g_part[1])

    def signed_parts(self, sigma: float) -> Tuple[Tuple[float, float],
                                                  Tuple[float, float]]:
        """``((A_sym, A_anti), (𝒢_sym, 𝒢_anti))`` at ``λ = −σ²``.

        Kept separate because the *no-growing-mode* statement is about their
        signs and not about their difference: a sweep that only reports
        ``A − 𝒢 < 0`` shows there is no root here, while the two sign facts show
        there can be no root anywhere.
        """
        s = float(sigma)
        kappa = math.sqrt(s ** 2 + float(self.interior_mass) ** 2)
        x = math.exp(-kappa * float(self.length))
        scale = float(self.area) * kappa
        a_sym = -(1.0 + x) / ((1.0 - x) * scale)
        a_anti = -(1.0 - x) / ((1.0 + x) * scale)
        g_self, g_cross = screened_products(self.radius, self.separation, s)
        return (a_sym, a_anti), (g_self + g_cross, g_self - g_cross)

    # -- the linearized mouth, for the contrast ------------------------------
    def linearized_channels(self, sigma: float) -> Tuple[float, float]:
        """PR #260's mouth: ``G(a,λ) → 1/(4πa) + g(λ)``, a **constant** shift.

        The same object with the screening switched off.  Its root is what this
        round is diagnosing.
        """
        s = float(sigma)
        (a_sym, a_anti), _ = self.signed_parts(s)
        decay = math.exp(-2.0 * math.pi * s)
        g = (-s * (1.0 + decay) / (FOUR_PI * (1.0 - decay)) if s > 1e-12
             else -1.0 / (4.0 * math.pi ** 2))
        cross = mouth_green(self.separation, -(s ** 2))
        shift = 1.0 / (FOUR_PI * self.radius) + g
        return (a_sym - (shift + cross), a_anti - (shift - cross))

    # -- the spectrum --------------------------------------------------------
    def growing_modes(self, sigma_max: float = 400.0,
                      n: int = 2000) -> Dict[str, object]:
        """Roots at ``λ < 0``.  **There are none**, and the sign structure says
        why: the returned worst values are the closest either channel comes."""
        sig = np.geomspace(1e-4, float(sigma_max), int(n))
        out: Dict[str, object] = {}
        for name, idx in (("symmetric", 0), ("antisymmetric", 1)):
            vals = [self.negative_lambda_channels(float(s))[idx] for s in sig]
            out[name] = None if max(vals) < 0.0 else float(
                sig[int(np.argmax(np.array(vals) >= 0.0))])
            out[f"worst_{name}"] = float(max(vals))
        out["the_tube_side_is_negative"] = bool(all(
            max(self.signed_parts(float(s))[0]) < 0.0 for s in sig))
        out["the_ambient_side_is_positive"] = bool(all(
            min(self.signed_parts(float(s))[1]) > 0.0 for s in sig))
        return out

    def bound_states(self, n: int = 4000) -> Dict[str, object]:
        """Roots in ``λ ∈ (0, 1)`` — below the free ESU gap.

        There is exactly one, in the symmetric channel: the tube's zero mode,
        given a restoring force by the two mouths' capacitance.
        """
        lams = np.linspace(1e-9, 1.0 - 1e-9, int(n))
        out: Dict[str, object] = {}
        for name, idx in (("symmetric", 0), ("antisymmetric", 1)):
            vals = np.array([self.channel_functions(float(x))[idx]
                             for x in lams])
            sign = np.sign(vals)
            cross = np.where(sign[:-1] * sign[1:] < 0.0)[0]
            root = None
            if cross.size:
                lo, hi = lams[cross[0]], lams[cross[0] + 1]
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    if (self.channel_functions(lo)[idx]
                            * self.channel_functions(mid)[idx]) <= 0.0:
                        hi = mid
                    else:
                        lo = mid
                root = 0.5 * (lo + hi)
            out[name] = root
        out["capacitance_estimate"] = self.soft_mode_closed_form()
        return out

    def soft_mode_closed_form(self) -> float:
        """``λ₀ → 8πa/(𝒜L)`` — two mouth capacitances over the tube's volume.

        A sphere of radius ``a`` has capacitance ``4πa``: holding its surface at
        value ``v`` costs charge ``4πa·v``.  Two of them restore a tube whose
        inertia is its volume ``𝒜L``, so the softest mode sits at
        ``λ₀ = 2·4πa/(𝒜L)``.  Exact in the limit, and the approach is measured.
        """
        return 8.0 * math.pi * self.radius / (self.area * self.length)

    # -- the frequency-domain response ---------------------------------------
    def response_spectrum(self, omegas: np.ndarray) -> np.ndarray:
        """``R(ω) = (A(ω) − 𝒢(ω))⁻¹`` on the retarded contour.

        With no growing mode the contour no longer has to clear one, so ``ε`` is
        back to PR #259's single requirement — comfortably above the frequency
        spacing — and nothing else.
        """
        om = np.asarray(omegas, dtype=complex)
        lam = om ** 2
        c = 1.0 / (self.area * om)
        x = om * self.length
        f = np.sin(om * self.radius) / (om * math.sin(self.radius))
        g_self = f * green_omega(self.radius, om)
        g_cross = f * f * green_omega(self.separation, om)
        m = np.empty(om.shape + (2, 2), dtype=complex)
        m[..., 0, 0] = m[..., 1, 1] = c / np.tan(x) - g_self
        m[..., 0, 1] = m[..., 1, 0] = c / np.sin(x) - g_cross
        det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
        r = np.empty_like(m)
        r[..., 0, 0] = m[..., 1, 1] / det
        r[..., 1, 1] = m[..., 0, 0] / det
        r[..., 0, 1] = -m[..., 0, 1] / det
        r[..., 1, 0] = -m[..., 1, 0] / det
        del lam
        return r

    def impulse_response(self, grid: RetardedGrid, width: float = 0.03
                         ) -> Dict[str, np.ndarray]:
        om = grid.omegas
        spec = GaussianPulse(amplitude=1.0, carrier=0.0, width=float(width),
                             t0=0.0).spectrum(om)
        r = self.response_spectrum(om)
        return {"times": grid.times,
                "same_mouth": np.real(
                    np.fft.fft(r[:, 0, 0] * spec) / grid.span),
                "opposite_mouths": np.real(
                    np.fft.fft(r[:, 0, 1] * spec) / grid.span)}

    def static_response(self, lmbda: float = 1e-8) -> np.ndarray:
        s, t = self.channel_functions(float(lmbda))
        v = np.array([[1.0, 1.0], [1.0, -1.0]]) / math.sqrt(2.0)
        return v @ np.diag([1.0 / s, 1.0 / t]) @ v.T


WORKING_MOUTH = FiniteMouthThroat(separation=1.3, length=0.9, area=FOUR_PI,
                                  radius=0.05)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_mean_value_identities_hold(
        radii: Sequence[float] = (0.05, 0.15, 0.35),
        separations: Sequence[float] = (0.8, 1.3, 2.2),
        lams: Sequence[float] = (0.16, 2.89)) -> Dict[str, object]:
    """The model's two identities, against quadrature on ``S³``.

    Everything below rests on ``𝒢_self = f(a)G(a)`` and ``𝒢_cross = f(a)²G(d)``,
    so those are computed a second way — by integrating over the actual
    two-sphere at distance ``a`` — rather than asserted from the mean-value
    theorem.  The cross identity is exact and lands at ``2e-11``; the self one
    carries the Green function's integrable singularity at coincidence, so its
    residual is grid-limited and is reported as such.
    """
    cross, selves = [], []
    for a in radii:
        for lam in lams:
            selves.append({
                "radius": float(a), "lambda": float(lam),
                "quadrature": shell_average_self(a, lam),
                "prediction": regular_radial(a, lam) * mouth_green(a, lam)})
            for d in separations:
                if a >= 0.5 * d:
                    continue
                q = shell_average_cross(a, d, lam)
                p = regular_radial(a, lam) * mouth_green(d, lam)
                cross.append({"radius": float(a), "separation": float(d),
                              "lambda": float(lam), "quadrature": q,
                              "prediction": p,
                              "relative_error": abs(q - p) / abs(p)})
    for row in selves:
        row["relative_error"] = (abs(row["quadrature"] - row["prediction"])
                                 / abs(row["prediction"]))
    return {"cross": cross, "self": selves,
            "worst_cross_error": float(max(r["relative_error"]
                                           for r in cross)),
            "worst_self_error": float(max(r["relative_error"]
                                          for r in selves)),
            "the_cross_identity_is_exact": bool(
                max(r["relative_error"] for r in cross) < 1e-9),
            "the_self_identity_is_grid_limited": bool(
                max(r["relative_error"] for r in selves) < 1e-3),
            "what_it_justifies": ("the smeared mouth's ambient block, and with "
                                  "it every result in this module")}


def measure_the_negative_mode_does_not_survive(
        radii: Sequence[float] = (0.01, 0.05, 0.15, 0.35, 0.6),
        separations: Sequence[float] = (0.8, 1.3, 2.2, 3.0),
        lengths: Sequence[float] = (0.3, 0.9, 3.0),
        areas: Sequence[float] = (0.2, math.pi, FOUR_PI),
        masses: Sequence[float] = (0.0, 1.5),
        sigmas: Sequence[float] = (1e-3, 0.03, 0.3, 1.0, 3.0, 10.0, 30.0,
                                   100.0, 400.0)) -> Dict[str, object]:
    """**The answer: no.**  And it is structural, not parametric.

    At ``λ = −σ²`` the two sides of ``det(A − 𝒢) = 0`` have opposite signs for
    *every* admissible parameter choice:

    * the tube gives ``−coth(κL/2)/(𝒜κ)`` and ``−tanh(κL/2)/(𝒜κ)``, strictly
      **negative** — a passive interior cannot supply a restoring force;
    * the ambient gives ``f·G(a) ± f²·G(d)``, strictly **positive** — ``f`` and
      ``G`` are positive on the imaginary axis, and ``G(a) > f(a)G(d)`` whenever
      ``a < d/2``, which disjoint mouths require anyway.

    A difference of a negative and a positive number has no zero, so no sweep
    can find a growing mode and this one does not.  The sweep is reported with
    the *worst approach to zero* rather than as a bare pass, and the two sign
    facts are reported separately, because they are the actual argument.
    """
    rows = []
    worst = -math.inf
    tube_ok = ambient_ok = True
    for a in radii:
        for d in separations:
            if a >= 0.5 * d:
                continue
            for length in lengths:
                for area in areas:
                    for m in masses:
                        t = FiniteMouthThroat(separation=d, length=length,
                                              area=area, radius=a,
                                              interior_mass=m)
                        for s in sigmas:
                            (at, ab), (gt, gb) = t.signed_parts(float(s))
                            f_sym, f_anti = at - gt, ab - gb
                            worst = max(worst, f_sym, f_anti)
                            tube_ok = tube_ok and max(at, ab) < 0.0
                            ambient_ok = ambient_ok and min(gt, gb) > 0.0
                            rows.append({"radius": float(a),
                                         "separation": float(d),
                                         "length": float(length),
                                         "area": float(area), "mass": float(m),
                                         "sigma": float(s),
                                         "f_sym": f_sym, "f_anti": f_anti})
    positives = sum(1 for r in rows
                    if r["f_sym"] >= 0.0 or r["f_anti"] >= 0.0)
    return {"samples": len(rows), "positives": int(positives),
            "worst_approach_to_zero": float(worst),
            "the_tube_side_is_always_negative": bool(tube_ok),
            "the_ambient_side_is_always_positive": bool(ambient_ok),
            "there_is_no_growing_mode": bool(positives == 0 and tube_ok
                                             and ambient_ok),
            "the_answer": ("no: a finite-radius mouth removes the negative "
                           "mode, and the sign structure says no parameter "
                           "choice can bring it back")}


def measure_the_instability_was_the_linearization(
        throat: FiniteMouthThroat = WORKING_MOUTH,
        radii: Sequence[float] = (0.02, 0.05, 0.15, 0.35)
        ) -> Dict[str, object]:
    """**The diagnosis.**  PR #260's mode came from freezing a screened quantity.

    That round wrote the mouth as a constant shift ``1/(4πa)`` — the leading term
    of ``G(a,λ) = 1/(4πa) + g(λ) + O(a)``.  The exact ``G(a,−κ²)`` is screened,
    ``≈ e^{−κa}/(4πa)``, and dies; the constant does not.  So the linearized
    ambient eventually loses to the tube's ``−1/(𝒜κ)`` and the difference
    changes sign, while the exact one never does.

    Measured: where the two agree, where they part, and where the linearized
    model's root sits — at ``κa ≈ 1``, which is exactly the boundary of the
    approximation that produced it.  **A mode living at the scale where its own
    derivation fails is an artifact, and this is the demonstration rather than
    the suspicion PR #260 could only record.**
    """
    rows = []
    for a in radii:
        t = FiniteMouthThroat(separation=throat.separation,
                              length=throat.length, area=throat.area,
                              radius=float(a))
        root = None
        sig = np.geomspace(1e-3, 1e4, 4000)
        vals = np.array([t.linearized_channels(float(s))[0] for s in sig])
        cross = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0.0)[0]
        if cross.size:
            lo, hi = sig[cross[0]], sig[cross[0] + 1]
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if (t.linearized_channels(lo)[0]
                        * t.linearized_channels(mid)[0]) <= 0.0:
                    hi = mid
                else:
                    lo = mid
            root = 0.5 * (lo + hi)
        agree = []
        for ka in (0.02, 0.1, 1.0, 3.0):
            s = ka / a
            exact = t.negative_lambda_channels(s)[0]
            lin = t.linearized_channels(s)[0]
            agree.append({"kappa_a": float(ka),
                          "exact": exact, "linearized": lin,
                          "relative_gap": abs(exact - lin) / abs(exact)})
        rows.append({"radius": float(a),
                     "linearized_root": root,
                     "root_times_radius": (None if root is None
                                           else float(root * a)),
                     "exact_has_a_root": bool(
                         t.growing_modes()["symmetric"] is not None),
                     "agreement": agree})
    deep = [g for r in rows for g in r["agreement"] if g["kappa_a"] <= 0.1]
    late = [g for r in rows for g in r["agreement"] if g["kappa_a"] >= 3.0]
    roots = [r["root_times_radius"] for r in rows
             if r["root_times_radius"] is not None]
    return {"rows": rows,
            "worst_gap_below_kappa_a_of_0p1": float(max(
                g["relative_gap"] for g in deep)),
            "worst_gap_at_kappa_a_of_3": float(min(
                g["relative_gap"] for g in late)),
            "linearized_roots_times_radius": roots,
            "every_linearized_model_has_a_root": bool(
                len(roots) == len(rows)),
            "no_exact_model_has_one": bool(
                not any(r["exact_has_a_root"] for r in rows)),
            "the_root_sits_at_kappa_a_of_order_one": bool(
                all(0.2 < v < 5.0 for v in roots)),
            "the_diagnosis": ("the linearization is excellent for κa ≪ 1 and "
                              "wrong in SIGN beyond it; the mode it produces "
                              "sits at κa ≈ 1, the edge of its own validity")}


def measure_the_mode_became_soft_and_positive(
        throat: FiniteMouthThroat = WORKING_MOUTH,
        radii: Sequence[float] = (0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
        ) -> Dict[str, object]:
    """Where the mode went: **below the gap, and above zero.**

    The composite has exactly one state under the free ESU gap ``λ = 1``, in the
    symmetric channel — the tube's zero mode, which a point cannot hold and two
    spheres can.  Its position has a closed form,

        ``λ₀  ⟶  8πa/(𝒜L)`` ,

    two mouth capacitances ``4πa`` restoring a tube of volume ``𝒜L``, approached
    to ``0.2%`` by ``a = 0.005``.  So shrinking the mouth drives the mode to zero
    **from above**: it becomes soft, and never unstable.

    That is the precise sense in which PR #260's model failed.  It did not get a
    rate slightly wrong; it took a mode approaching zero from above and put it
    on the other side, at ``λ ≈ −1/a²`` instead of ``λ ≈ +8πa/(𝒜L)``.
    """
    rows = []
    for a in radii:
        t = FiniteMouthThroat(separation=throat.separation,
                              length=throat.length, area=throat.area,
                              radius=float(a))
        states = t.bound_states()
        lam0 = states["symmetric"]
        closed = t.soft_mode_closed_form()
        rows.append({"radius": float(a),
                     "lambda_0": None if lam0 is None else float(lam0),
                     "omega_0": None if lam0 is None else math.sqrt(lam0),
                     "closed_form": float(closed),
                     "ratio": None if lam0 is None else float(lam0 / closed),
                     "antisymmetric": states["antisymmetric"]})
    small = [r for r in rows if r["radius"] <= 0.02]
    return {"rows": rows,
            "every_one_is_positive": bool(all(
                r["lambda_0"] is not None and r["lambda_0"] > 0.0
                for r in rows)),
            "every_one_is_below_the_gap": bool(all(
                r["lambda_0"] < 1.0 for r in rows)),
            "no_antisymmetric_bound_state": bool(all(
                r["antisymmetric"] is None for r in rows)),
            "worst_closed_form_error": float(max(
                abs(1.0 - r["ratio"]) for r in small)),
            "the_capacitance_formula_holds": bool(
                max(abs(1.0 - r["ratio"]) for r in small) < 0.02),
            "the_point_limit_approaches_zero_from_above": bool(
                rows[-1]["lambda_0"] < rows[0]["lambda_0"]
                and rows[-1]["lambda_0"] > 0.0),
            "what_it_means": ("the mode is soft, not unstable; PR #260 put it "
                              "on the wrong side of zero")}


def measure_the_delay_survives_with_a_radius_correction(
        radius: float = 0.02, separation: float = 1.3, area: float = FOUR_PI,
        lengths: Sequence[float] = (0.4, 0.6, 0.9, 1.2),
        long_lengths: Sequence[float] = (2.0, 3.0),
        radii: Sequence[float] = (0.01, 0.02, 0.03, 0.04),
        width: float = 0.03, eps: float = 0.4) -> Dict[str, object]:
    """PR #260's headline result, re-measured with the mouth resolved.

    It survives, and cleanly: the cross-mouth onset still moves with slope ``1``
    in ``L`` and still saturates at the ambient path ``d`` exactly.  The mouth
    adds only a **sub-leading** shift, earlier by ``O(a)`` — measured slope
    ``d(onset)/da ≈ −0.39``, an order of magnitude below the leading term —
    because the wave now reaches the mouth's *surface* rather than its centre.

    A first draft predicted ``−2a``, one radius per leg, from a version of the
    ambient block with the shell form factor ``f(a)`` left out.  Restoring ``f``
    changes the answer: ``f`` carries an advance of its own that partly cancels
    the surface offset.  The measured slope is quoted rather than the predicted
    one, and the prediction is recorded as wrong.

    Note the contour: with no growing mode to clear, ``ε`` is back to PR #259's
    single requirement of sitting well above the frequency spacing.  PR #260
    needed ``ε > σ* ≈ 2``; here ``0.4`` is comfortable.
    """
    rows = []
    for length in list(lengths) + list(long_lengths):
        t = FiniteMouthThroat(separation=separation, length=float(length),
                              area=area, radius=float(radius))
        grid = RetardedGrid(n=1 << 17, span=300.0, eps=float(eps))
        imp = t.impulse_response(grid, width=width)
        rows.append({"length": float(length),
                     "onset_same_mouth": causal_onset(imp["same_mouth"],
                                                      imp["times"]),
                     "onset_opposite": causal_onset(imp["opposite_mouths"],
                                                    imp["times"]),
                     "prediction": float(min(length, separation)
                                         - 2.0 * radius)})
    short = [r for r in rows if r["length"] < separation]
    long = [r for r in rows if r["length"] > separation]
    slope = float(np.polyfit([r["length"] for r in short],
                             [r["onset_opposite"] for r in short], 1)[0])
    by_radius = []
    for a in radii:
        t = FiniteMouthThroat(separation=separation, length=0.9, area=area,
                              radius=float(a))
        grid = RetardedGrid(n=1 << 17, span=300.0, eps=float(eps))
        imp = t.impulse_response(grid, width=width)
        by_radius.append({"radius": float(a),
                          "onset": causal_onset(imp["opposite_mouths"],
                                                imp["times"])})
    radius_slope = float(np.polyfit([r["radius"] for r in by_radius],
                                    [r["onset"] for r in by_radius], 1)[0])
    return {"rows": rows, "by_radius": by_radius,
            "slope_in_length": slope,
            "slope_in_radius": radius_slope,
            "onset_spread_above_d": float(
                max(r["onset_opposite"] for r in long)
                - min(r["onset_opposite"] for r in long)),
            "contour": float(eps),
            "the_traversal_time_survives": bool(abs(slope - 1.0) < 0.05),
            "the_ambient_path_still_takes_over": bool(
                max(r["onset_opposite"] for r in long)
                - min(r["onset_opposite"] for r in long) < 1e-9),
            "the_mouth_shift_is_subleading": bool(
                -1.0 < radius_slope < 0.0),
            "the_correction": ("the mouth shifts the onset earlier by O(a), "
                               "well below the leading min(L,d); a first draft "
                               "predicted −2a from an ambient block missing "
                               "the shell form factor")}


def measure_the_static_results_survive(
        throat: FiniteMouthThroat = WORKING_MOUTH,
        lams: Sequence[float] = (1e-8, 1e-6, 1e-4)) -> Dict[str, object]:
    """PR #260's static statements, re-run with the mouth resolved.

    Both survive, because both came from the *tube*'s zero mode and not from the
    mouth: the static response still collapses onto the antisymmetric direction,
    ``det S`` still vanishes linearly in ``λ``, and PR #258's defect is still
    ``𝒲 = −β(λ)`` exactly.  Only the coefficient moves — the ambient's diagonal
    is now ``+f(a)G(a)`` instead of the renormalized ``g₀ < 0``, so ``det S/λ``
    changes sign and magnitude, which is reported rather than smoothed over.
    """
    rows = []
    for lam in lams:
        s = throat.static_response(float(lam))
        det = float(np.linalg.det(s))
        beta = float(throat.tube().transmission(complex(float(lam))).real)
        g_cross = float(throat.ambient_matrix(float(lam))[0, 1])
        rows.append({"lambda": float(lam), "det_S": det,
                     "det_S_over_lambda": det / float(lam),
                     "antisymmetry": float(abs(s[0, 0] + s[0, 1])
                                           / abs(s[0, 0])),
                     "defect": float(s[0, 1] / det - g_cross),
                     "minus_beta": -beta,
                     "defect_error": float(abs(s[0, 1] / det - g_cross
                                               + beta) / abs(beta))})
    coeff = [r["det_S_over_lambda"] for r in rows]
    # as in PR #260, the collapse is a limit: the claim is made where λ → 0 and
    # the largest λ is carried as the row where the correction turns on
    small = [r for r in rows if r["lambda"] <= 1e-6]
    return {"rows": rows,
            "linear_coefficient": float(np.mean(coeff)),
            "det_S_is_linear_in_lambda": bool(
                (max(coeff) - min(coeff)) / abs(np.mean(coeff)) < 5e-3),
            "worst_antisymmetry": float(max(r["antisymmetry"] for r in small)),
            "the_static_response_is_still_rank_one": bool(
                max(r["antisymmetry"] for r in small) < 1e-4),
            "worst_defect_error": float(max(r["defect_error"] for r in rows)),
            "the_defect_is_still_minus_beta": bool(
                max(r["defect_error"] for r in rows) < 1e-8),
            "why": ("both came from the tube's zero mode, which the mouth does "
                    "not touch")}


def measure_monopole_matching_is_the_remaining_approximation(
        radii: Sequence[float] = (0.02, 0.05, 0.15, 0.35),
        separation: float = 1.3, lmbda: float = 0.16,
        n_modes: int = 6) -> Dict[str, object]:
    """What is still put in, quantified rather than confessed.

    One channel per mouth means only the ``ℓ = 0`` projection of the neighbour's
    field on each sphere is coupled.  The dropped content is measured directly —
    the neighbour's field expanded on the sphere and its ``ℓ ≥ 1`` weight taken —
    and it obeys PR #250's screening law: the ``ℓ``-th multipole is suppressed
    like ``(a/d)^ℓ``, so the leading omission is the dipole, at ``O(a/d)``.

    That is the honest boundary of this round: the mouths are resolved as
    *spheres carrying one number each*, not as a solved neck geometry.
    """
    rows = []
    for a in radii:
        th = np.linspace(0.0, math.pi, 20001)
        cos2 = (math.cos(a) * math.cos(separation)
                + math.sin(a) * math.sin(separation) * np.cos(th))
        chi2 = np.arccos(np.clip(cos2, -1.0, 1.0))
        field = np.array([mouth_green(float(c), lmbda) for c in chi2])
        weight = 0.5 * np.sin(th)
        moments = []
        for ell in range(n_modes):
            leg = np.polynomial.legendre.Legendre.basis(ell)(np.cos(th))
            moments.append(float(np.trapezoid(field * leg * weight, th))
                           * (2 * ell + 1))
        total = float(np.trapezoid(field ** 2 * weight, th))
        kept = moments[0] ** 2
        rows.append({"radius": float(a),
                     "monopole": moments[0],
                     "dipole": moments[1],
                     "dipole_over_monopole": abs(moments[1] / moments[0]),
                     "ratio_over_a_over_d": abs(moments[1] / moments[0])
                     / (a / separation),
                     "dropped_power_fraction": float(
                         max(total - kept, 0.0) / total)})
    ratios = [r["ratio_over_a_over_d"] for r in rows]
    return {"lambda": float(lmbda), "separation": float(separation),
            "rows": rows,
            "dipole_scales_like_a_over_d": bool(
                (max(ratios) - min(ratios)) / abs(np.mean(ratios)) < 0.35),
            "worst_dropped_fraction": float(max(r["dropped_power_fraction"]
                                                for r in rows)),
            "smallest_dropped_fraction": float(
                min(r["dropped_power_fraction"] for r in rows)),
            "the_scope": ("one channel per mouth couples only ℓ = 0; the "
                          "dropped multipoles are screened as (a/d)^ℓ, PR "
                          "#250's law, so the leading omission is the dipole")}
