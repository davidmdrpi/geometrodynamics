"""The **neck**: the balls actually removed, and the answer made exact.

Where this sits
───────────────
PR #260 found a growing mode and gated the roadmap on whether it survives a
finite-radius mouth.  PR #261 answered **no** — but with a stated limitation, and
it is the load-bearing one: those mouths were *spheres in a fixed ambient*.  The
balls were never removed.  The model smeared the coupling over ``∂B_a`` while
still using the **whole sphere's** Green function, and coupled through **one
channel** per mouth, so only ``ℓ = 0`` talked to the tube.

This round removes the balls.  The ambient is ``Ω = S³ ∖ (B_a(c₁) ∪ B_a(c₂))``,
the tube is glued along the boundary spheres, and the question is re-asked.

**The answer is still no — and here it is a theorem rather than a sweep.**

The theorem
───────────
With the balls removed there is no subtraction anywhere, so the composite energy
is a sum of manifestly non-negative terms:

    ``E[φ,u] = ∫_Ω (|∇φ|² + φ²) dV  +  𝒜∫₀^L (|u'|² + m²|u|²) ds``

Every term is ``≥ 0``, and ``E = 0`` forces ``φ ≡ 0`` on ``Ω`` — whereupon the
matching gives ``u(0) = u(L) = 0`` and Poincaré gives
``𝒜∫|u'|² ≥ (π/L)² 𝒜∫|u|²``, so ``u ≡ 0`` too.  Hence

    ``λ  =  E/‖·‖²  >  0``   for every configuration, **all multipoles, no
    truncation, no parameter sweep.**

That is the whole argument, and it is why removing the balls is not a refinement
of PR #261 but a change of footing: that round had to establish a *sign* on a
reduced ``2×2``; this one has positivity of the form itself.  The renormalized
``g(λ) < 0`` that PR #260's mode fed on has nowhere to enter, because nothing is
renormalized.

What is measured anyway
───────────────────────
A theorem with no numbers attached is a claim, so:

* the **exterior Dirichlet-to-Neumann map** ``N_ℓ(λ)`` is computed by shooting
  the radial equation ``v'' + [λ − ℓ(ℓ+2)/sin²χ]v = 0`` from the far pole, and
  agrees with the closed form in the ``ℓ = 0`` channel to ``1e-13``.  It is
  **positive for every ``ℓ`` and every ``λ < 0``**, and *increasing* in ``ℓ`` —
  higher multipoles are stiffer, so they cannot be the ones to go soft;
* the **Rayleigh quotient** of explicit trial configurations is positive, with
  its infimum matching the computed lowest mode;
* the ``ℓ ≥ 1`` sectors **decouple exactly** from a one-channel tube.  They are
  spectators with their own positive spectrum, which means PR #261's monopole
  truncation was never a limitation *for stability* — only for amplitudes;
* and what PR #261's fixed ambient actually cost is quantified: its
  ``f(a)G(a)`` against this round's ``1/N₀(a)``, agreeing to ``3.8e-03`` at the
  working radius ``a = 0.05`` (and to ``1.3e-04`` at ``a = 0.02``, ``λ = 0``),
  departing to ``11%`` at ``a = 0.35``, ``λ = −4``.

What survives
─────────────
The soft mode.  With the ball removed ``λ₀`` moves by ``1.2e-04`` relative at
``a = 0.02`` and still tends to ``8πa/(𝒜L)`` — two mouth capacitances restoring
a tube of volume ``𝒜L``.  PR #261's answer was right, and this round says by how
much its approximation was right.

And it is now **forced** rather than found.  The same style of argument that
kills the growing mode produces this one, from the two ends of the gap:
``F_sym → +∞`` as ``λ → 0⁺`` (the tube's ``2/(𝒜λL)``) and ``F_sym → −∞`` as
``λ → 1⁻``, because the exterior stiffness vanishes there —
``N₀ → 2π(π − a + sin a cos a)(1 − λ)``, reproduced to ``3e-05``.  A continuous
function running from ``+∞`` to ``−∞`` has a zero.

One correction to PR #261
─────────────────────────
That round reported "exactly one state below the gap".  That is a statement
about ``L < π``, not a structural one.  The channel functions have **poles** at
``λ = (2πj/L)²`` and ``λ = (π(2j−1)/L)²``; above ``L = π`` these enter the gap
and each brings a genuine extra state just above it — three states at ``L = 8``.
A pole is a sign change with no zero, so crossing-counting alone miscounts;
classifying by the residual separates roots from poles by fifteen orders of
magnitude.  None of it touches the stability answer: every one of those states
is positive, as the quadratic form requires.

What is still put in
────────────────────
The tube has **one transverse channel**, so ``𝒜`` is a coupling and the neck is
a quantum-graph edge rather than a solved geometry.  The ambient is a **fixed**
round ``S³`` with two balls cut out — the metric is not solved for.  **No
backreaction.**
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .finite_mouth import FOUR_PI, FiniteMouthThroat, mouth_green, regular_radial

__all__ = [
    "NeckThroat",
    "WORKING_NECK",
    "exterior_log_derivative",
    "exterior_dtn",
    "exterior_dtn_monopole",
    "rayleigh_quotient",
    "measure_the_exterior_dtn_is_positive_in_every_channel",
    "measure_the_quadratic_form_is_positive",
    "measure_the_negative_mode_does_not_survive_the_neck",
    "measure_the_higher_multipoles_decouple",
    "measure_what_the_fixed_ambient_cost",
    "measure_the_soft_mode_is_forced_by_the_two_ends",
    "measure_the_soft_mode_survives_the_removal",
]


# ════════════════════════════════════════════════════════════════════════════
# THE EXTERIOR OF A REMOVED BALL
# ════════════════════════════════════════════════════════════════════════════
def exterior_log_derivative(radius: float, lmbda: float, ell: int,
                            eps: float = 1e-7) -> float:
    """``ψ'/ψ`` at ``χ = a`` for the solution regular at the far pole.

    Substituting ``ψ = v/sin χ`` turns the ``S³`` radial equation for
    ``−∇² + 1`` into a one-dimensional Schrödinger problem with **no first
    derivative**,

        ``v'' + [λ − ℓ(ℓ+2)/sin²χ] v = 0`` ,

    which is both the cleanest thing to integrate and the clearest statement of
    the physics: at ``λ < 0`` the bracket is negative everywhere, so ``v`` cannot
    oscillate and the exterior map has a sign.

    Regularity at the antipode ``χ = π`` fixes the branch: with ``e = π − χ`` the
    two behaviours are ``e^{ℓ+1}`` and ``e^{−ℓ}``, and only the first leaves
    ``ψ = v/sin χ`` finite.  The integration starts from it.
    """
    a, lam, l = float(radius), float(lmbda), int(ell)
    if not 0.0 < a < math.pi:
        raise ValueError("the removed ball needs 0 < a < π")

    def rhs(chi, y):
        return [y[1], -(lam - l * (l + 2) / math.sin(chi) ** 2) * y[0]]

    start = math.pi - eps
    y0 = [eps ** (l + 1), -(l + 1) * eps ** l]
    sol = solve_ivp(rhs, (start, a), y0, rtol=1e-11, atol=1e-300,
                    method="DOP853")
    v, dv = float(sol.y[0, -1]), float(sol.y[1, -1])
    return dv / v - 1.0 / math.tan(a)


def exterior_dtn_monopole(radius: float, lmbda: float) -> float:
    """``N₀(λ)`` in closed form — the ``ℓ = 0`` exterior map.

    The monopole solution regular at the antipode is the free Green function
    itself, so this needs no integration:
    ``ψ ∝ sin(ω(π−χ))/sin χ`` and its hyperbolic continuation.  Kept as an
    independent construction, so the shooting integrator is checked against
    something rather than against a finer version of itself.
    """
    a, lam = float(radius), float(lmbda)
    area = FOUR_PI * math.sin(a) ** 2
    cot_a = 1.0 / math.tan(a)
    if abs(lam) < 1e-14:
        return -area * (-1.0 / (math.pi - a) - cot_a)
    if lam < 0.0:
        k = math.sqrt(-lam)
        return -area * (-k / math.tanh(k * (math.pi - a)) - cot_a)
    w = math.sqrt(lam)
    tan_arg = math.tan(w * (math.pi - a))
    if abs(tan_arg) < 1e-14:
        raise ValueError("λ sits on an exterior Dirichlet resonance")
    return -area * (-w / tan_arg - cot_a)


def exterior_dtn(radius: float, lmbda: float, ell: int = 0) -> float:
    """``N_ℓ(λ) = −4π sin²a · ψ'(a)/ψ(a)`` — the flux the exterior returns.

    Positive at ``λ < 0`` in every channel, and **increasing in ``ℓ``**: a
    stiffer angular pattern costs more, so the monopole is the softest channel
    and the higher ones cannot be the first to go unstable.  In the small-ball
    limit ``N₀ → 4πa``, the capacitance of a sphere, which is the sanity check
    that the normalization is the physical one.
    """
    if int(ell) == 0:
        return exterior_dtn_monopole(radius, lmbda)
    a = float(radius)
    return -(FOUR_PI * math.sin(a) ** 2) * exterior_log_derivative(
        a, lmbda, int(ell))


# ════════════════════════════════════════════════════════════════════════════
# THE NECK
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class NeckThroat:
    """PR #261's throat with the balls **removed** from the ambient.

    One line changes: the ambient's diagonal block is ``1/N₀(a,λ)`` — the
    inverse exterior DtN — where PR #261 used ``f(a)G(a)``, the double-shell
    average taken over an ambient that still contained the balls.

    The off-diagonal keeps the smeared form ``f(a)²G(d)``.  That is the one
    approximation left in the reduced ``2×2``, and it is a **single-scattering**
    one: the exact two-ball exterior would add the series in which a monopole on
    one boundary sphere drives the other and is driven back.  Its expansion
    parameter is the ratio the model already computes, ``cross/self ≈ 0.8·(a/d)``
    — ``0.031`` at the working point — so the neglected terms enter at
    ``~1e-03`` there and at worst ``0.16`` anywhere sampled.  They cannot reach
    the sign, and `measure_what_the_fixed_ambient_cost` reports the bound rather
    than leaving it as a word.

    None of this touches the theorem, which is about the exact problem on
    ``Ω``, not about this reduction.
    """

    separation: float = 1.3
    length: float = 0.9
    area: float = FOUR_PI
    radius: float = 0.05
    interior_mass: float = 0.0

    def __post_init__(self) -> None:
        if not 0.0 < self.radius < 0.5 * self.separation:
            raise ValueError("the removed balls must be disjoint: 0 < a < d/2")

    def fixed_ambient(self) -> FiniteMouthThroat:
        """PR #261's object with the same parameters, for the comparison."""
        return FiniteMouthThroat(separation=self.separation,
                                 length=self.length, area=self.area,
                                 radius=self.radius,
                                 interior_mass=self.interior_mass)

    def ambient_channels(self, lmbda: float) -> Tuple[float, float]:
        """``(𝒢_sym, 𝒢_anti)`` with the balls removed."""
        lam = float(lmbda)
        self_term = 1.0 / exterior_dtn_monopole(self.radius, lam)
        cross = (regular_radial(self.radius, lam) ** 2
                 * mouth_green(self.separation, lam))
        return self_term + cross, self_term - cross

    def tube_channels(self, lmbda: float) -> Tuple[float, float]:
        a = self.fixed_ambient().tube().boundary_matrix(complex(float(lmbda)))
        return (float((a[0, 0] + a[0, 1]).real),
                float((a[0, 0] - a[0, 1]).real))

    def channel_functions(self, lmbda: float) -> Tuple[float, float]:
        at, ab = self.tube_channels(lmbda)
        gt, gb = self.ambient_channels(lmbda)
        return at - gt, ab - gb

    def signed_parts(self, sigma: float) -> Tuple[Tuple[float, float],
                                                  Tuple[float, float]]:
        """The two sides at ``λ = −σ²``, kept apart — the tube's strictly
        negative, the exterior's strictly positive."""
        s = float(sigma)
        kappa = math.sqrt(s ** 2 + float(self.interior_mass) ** 2)
        x = math.exp(-kappa * float(self.length))
        scale = float(self.area) * kappa
        tube = (-(1.0 + x) / ((1.0 - x) * scale),
                -(1.0 - x) / ((1.0 + x) * scale))
        return tube, self.ambient_channels(-(s ** 2))

    def growing_modes(self, sigma_max: float = 200.0,
                      n: int = 400) -> Dict[str, object]:
        """Scan ``λ < 0`` for a zero.  There is none — the signs forbid it.

        Routed through `signed_parts` rather than `channel_functions`: the
        tube's ``cot``/``csc`` form overflows above ``σ ≈ 350`` while the
        decaying-exponential form does not, and a sweep meant to establish an
        absence should not be the place where an overflow decides the answer.
        """
        sig = np.geomspace(1e-3, float(sigma_max), int(n))
        parts = [self.signed_parts(float(s)) for s in sig]
        out: Dict[str, object] = {}
        for name, idx in (("symmetric", 0), ("antisymmetric", 1)):
            vals = [t[idx] - g[idx] for t, g in parts]
            out[name] = None if max(vals) < 0.0 else float(
                sig[int(np.argmax(np.array(vals) >= 0.0))])
            out[f"worst_{name}"] = float(max(vals))
        return out

    def tube_poles_in_the_gap(self) -> Dict[str, List[float]]:
        """Where the tube's own longitudinal harmonics fall below ``λ = 1``.

        ``A_sym = cot(kL/2)/(𝒜k)`` blows up at ``kL/2 = jπ`` and ``A_anti =
        −tan(kL/2)/(𝒜k)`` at ``kL/2 = (j−½)π``, so the channel functions have
        **poles** at ``λ = (2πj/L)²`` and ``λ = (π(2j−1)/L)²``.  The first of
        them enters the gap at ``L = π``, and it matters because a pole is a
        sign change without a zero: counting crossings alone reports a bound
        state that is not there.
        """
        out: Dict[str, List[float]] = {"symmetric": [], "antisymmetric": []}
        j = 1
        while (2.0 * math.pi * j / self.length) ** 2 < 1.0:
            out["symmetric"].append((2.0 * math.pi * j / self.length) ** 2)
            j += 1
        j = 1
        while (math.pi * (2 * j - 1) / self.length) ** 2 < 1.0:
            out["antisymmetric"].append(
                (math.pi * (2 * j - 1) / self.length) ** 2)
            j += 1
        return out

    def bound_states(self, n: int = 1200,
                     pole_cut: float = 1.0) -> Dict[str, object]:
        """Every state below the free ESU gap, with poles thrown out.

        Each bracketed sign change is bisected and then **classified** by the
        size of the residual: a genuine root leaves ``|f| ~ 1e-15``, a pole
        leaves ``|f| ~ 1e+15``.  The separation is fifteen orders of magnitude
        wide, so ``pole_cut`` is not a tuned threshold.

        ``symmetric`` and ``antisymmetric`` keep returning the **lowest** root
        in each channel, which is what the earlier rounds meant by them; the
        full lists are alongside, because for ``L > π`` there is more than one.
        """
        lams = np.linspace(1e-9, 1.0 - 1e-9, int(n))
        out: Dict[str, object] = {}
        for name, idx in (("symmetric", 0), ("antisymmetric", 1)):
            vals = np.array([self.channel_functions(float(x))[idx]
                             for x in lams])
            roots, poles = [], []
            for i in np.where(np.sign(vals[:-1])
                              * np.sign(vals[1:]) < 0.0)[0]:
                lo, hi = float(lams[i]), float(lams[i + 1])
                for _ in range(120):
                    mid = 0.5 * (lo + hi)
                    if (self.channel_functions(lo)[idx]
                            * self.channel_functions(mid)[idx]) <= 0.0:
                        hi = mid
                    else:
                        lo = mid
                mid = 0.5 * (lo + hi)
                (roots if abs(self.channel_functions(mid)[idx]) < pole_cut
                 else poles).append(mid)
            out[name] = roots[0] if roots else None
            out[f"{name}_roots"] = roots
            out[f"{name}_poles"] = poles
        out["states_below_the_gap"] = (len(out["symmetric_roots"])
                                       + len(out["antisymmetric_roots"]))
        out["capacitance_estimate"] = (8.0 * math.pi * self.radius
                                       / (self.area * self.length))
        return out


WORKING_NECK = NeckThroat(separation=1.3, length=0.9, area=FOUR_PI,
                          radius=0.05)


# ════════════════════════════════════════════════════════════════════════════
# THE QUADRATIC FORM
# ════════════════════════════════════════════════════════════════════════════
def rayleigh_quotient(throat: NeckThroat, profile, tube_profile,
                      n: int = 4001) -> Dict[str, float]:
    """``E/‖·‖²`` for an explicit configuration, in the separable geometry.

    ``profile(χ)`` is a radial ambient field on ``(a, π)`` and
    ``tube_profile(s)`` the interior one, matched at the mouth.  Both integrals
    use the honest measures — ``4π sin²χ dχ`` outside and ``𝒜 ds`` inside — so
    the number returned is the physical Rayleigh quotient and not a proxy.

    The point of computing it is that the theorem is a two-line argument, and a
    two-line argument deserves a check that the object it is about was built the
    way the argument assumes.
    """
    a = float(throat.radius)
    chi = np.linspace(a, math.pi - 1e-9, int(n))
    phi = np.asarray([profile(float(c)) for c in chi], dtype=float)
    dphi = np.gradient(phi, chi, edge_order=2)
    weight = FOUR_PI * np.sin(chi) ** 2
    e_amb = float(np.trapezoid((dphi ** 2 + phi ** 2) * weight, chi))
    n_amb = float(np.trapezoid(phi ** 2 * weight, chi))

    s = np.linspace(0.0, float(throat.length), int(n))
    u = np.asarray([tube_profile(float(x)) for x in s], dtype=float)
    du = np.gradient(u, s, edge_order=2)
    m2 = float(throat.interior_mass) ** 2
    e_tube = float(throat.area * np.trapezoid(du ** 2 + m2 * u ** 2, s))
    n_tube = float(throat.area * np.trapezoid(u ** 2, s))

    total_e, total_n = e_amb + e_tube, n_amb + n_tube
    return {"energy": total_e, "norm": total_n,
            "quotient": total_e / total_n,
            "ambient_energy": e_amb, "tube_energy": e_tube,
            "ambient_norm": n_amb, "tube_norm": n_tube,
            "poincare_floor": (math.pi / throat.length) ** 2}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_exterior_dtn_is_positive_in_every_channel(
        radii: Sequence[float] = (0.02, 0.2, 0.5),
        sigmas: Sequence[float] = (0.05, 1.0, 5.0, 30.0, 150.0),
        ells: Sequence[int] = (0, 1, 2, 5, 10)) -> Dict[str, object]:
    """The ambient half of the answer, with the ball genuinely gone.

    ``N_ℓ(λ) > 0`` for every channel and every ``λ < 0`` — computed by shooting
    the radial equation from the far pole, and checked in the ``ℓ = 0`` channel
    against an independent closed form.  Two further facts fall out and both
    matter:

    * ``N_ℓ`` **increases with ``ℓ``**, so the monopole is the softest channel.
      Whatever goes soft first is the one the tube couples to, which is the one
      PR #261 kept;
    * ``N₀ → 4πa`` as the ball shrinks — the capacitance of a sphere — which is
      what fixes the normalization to the physical one rather than a convention.
    """
    rows, worst_closed = [], 0.0
    monotone = True
    for a in radii:
        for s in sigmas:
            lam = -(s ** 2)
            values = [exterior_dtn(a, lam, int(l)) for l in ells]
            if int(ells[0]) == 0:
                shot = -(FOUR_PI * math.sin(a) ** 2) * \
                    exterior_log_derivative(a, lam, 0)
                worst_closed = max(worst_closed,
                                   abs(shot - values[0]) / abs(values[0]))
            monotone = monotone and all(
                b > c for c, b in zip(values, values[1:]))
            rows.append({"radius": float(a), "sigma": float(s),
                         "by_ell": [float(v) for v in values],
                         "min_over_ell": float(min(values))})
    caps = []
    for a in (0.02, 0.01, 0.005):
        caps.append({"radius": float(a),
                     "N0_at_zero": exterior_dtn_monopole(a, 0.0),
                     "four_pi_a": FOUR_PI * a,
                     "ratio": exterior_dtn_monopole(a, 0.0) / (FOUR_PI * a)})
    return {"ells": [int(x) for x in ells], "rows": rows,
            "capacitance": caps,
            "worst_shooting_vs_closed_form": float(worst_closed),
            "smallest_dtn_seen": float(min(r["min_over_ell"] for r in rows)),
            "it_is_positive_everywhere": bool(
                min(r["min_over_ell"] for r in rows) > 0.0),
            "it_increases_with_ell": bool(monotone),
            "the_capacitance_normalization_holds": bool(
                all(abs(c["ratio"] - 1.0) < 0.05 for c in caps)),
            "what_it_replaces": ("PR #261's f(a)G(a), which averaged over a "
                                 "sphere while keeping the ball in the ambient")}


def measure_the_quadratic_form_is_positive(
        throat: NeckThroat = WORKING_NECK, trials: int = 40,
        seed: int = 20260819) -> Dict[str, object]:
    """**The theorem**, and a check that the object obeys it.

    ``E = ∫_Ω(|∇φ|² + φ²) + 𝒜∫(|u'|² + m²|u|²)`` is a sum of non-negative terms
    with the balls removed — there is no subtraction anywhere to make one of
    them negative, which is the entire difference from a point interaction.  So
    ``λ = E/‖·‖² > 0`` for **every** configuration, with all multipoles and no
    truncation.

    The degenerate case is the one worth naming: a configuration with ``φ ≡ 0``
    outside has ``u(0) = u(L) = 0`` by matching, and Poincaré then gives
    ``𝒜∫|u'|² ≥ (π/L)²𝒜∫|u|²``.  So even a purely interior configuration is
    bounded below by ``π²/L²``, and nothing is left that could reach zero.

    Random trial configurations are evaluated with the honest measures.  Every
    quotient is positive, and the smallest sits above the lowest computed mode —
    which is what a variational bound is supposed to do.
    """
    rng = np.random.default_rng(int(seed))
    a = float(throat.radius)
    lowest = throat.bound_states()["symmetric"]
    rows = []
    for _ in range(int(trials)):
        coeffs = rng.normal(size=4)
        decay = float(abs(rng.normal()) + 0.4)

        def profile(chi, c=coeffs, k=decay, a0=a):
            x = (chi - a0) / (math.pi - a0)
            return float(math.exp(-k * x)
                         * (c[0] + c[1] * x + c[2] * x ** 2 + c[3] * x ** 3))

        end = profile(a)
        amp = float(rng.normal())

        def tube_profile(s, e=end, b=amp, L=float(throat.length)):
            return float(e + b * math.sin(math.pi * s / L))

        r = rayleigh_quotient(throat, profile, tube_profile, n=2001)
        rows.append({"quotient": r["quotient"],
                     "ambient_fraction": r["ambient_norm"] / r["norm"]})
    quotients = [r["quotient"] for r in rows]
    # the purely interior configuration, where Poincaré is the only floor
    interior = rayleigh_quotient(
        throat, lambda chi: 0.0,
        lambda s: math.sin(math.pi * s / float(throat.length)), n=4001)
    return {"trials": int(trials), "rows": rows[:8],
            "smallest_quotient": float(min(quotients)),
            "largest_quotient": float(max(quotients)),
            "lowest_computed_mode": float(lowest),
            "interior_only_quotient": float(interior["quotient"]),
            "poincare_floor": float(interior["poincare_floor"]),
            "every_quotient_is_positive": bool(min(quotients) > 0.0),
            "no_trial_beats_the_lowest_mode": bool(
                min(quotients) > float(lowest)),
            "the_interior_only_case_hits_the_poincare_floor": bool(
                abs(interior["quotient"] / interior["poincare_floor"] - 1.0)
                < 1e-6),
            "the_theorem": ("every term is non-negative and the degenerate "
                            "case is floored by π²/L², so λ > 0 for all "
                            "multipoles with no truncation and no sweep")}


def measure_the_negative_mode_does_not_survive_the_neck(
        radii: Sequence[float] = (0.01, 0.05, 0.15, 0.35, 0.6),
        separations: Sequence[float] = (0.8, 1.3, 2.2, 3.0),
        lengths: Sequence[float] = (0.3, 0.9, 3.0),
        areas: Sequence[float] = (0.2, math.pi, FOUR_PI),
        sigmas: Sequence[float] = (1e-3, 0.05, 0.5, 2.0, 10.0, 60.0, 200.0)
        ) -> Dict[str, object]:
    """The binary question, re-asked with the balls removed.

    Still no.  The sweep is the same shape as PR #261's and reaches the same
    place from a different construction — the diagonal is now ``1/N₀``, the
    inverse *exterior* map, rather than a smeared average over an ambient that
    still contained the balls.  Two independent models, one answer.
    """
    rows, worst, tube_ok, amb_ok = [], -math.inf, True, True
    for a in radii:
        for d in separations:
            if a >= 0.5 * d:
                continue
            for length in lengths:
                for area in areas:
                    t = NeckThroat(separation=d, length=length, area=area,
                                   radius=a)
                    for s in sigmas:
                        (at, ab), (gt, gb) = t.signed_parts(float(s))
                        worst = max(worst, at - gt, ab - gb)
                        tube_ok = tube_ok and max(at, ab) < 0.0
                        amb_ok = amb_ok and min(gt, gb) > 0.0
                        rows.append({"radius": float(a), "separation": float(d),
                                     "length": float(length),
                                     "area": float(area), "sigma": float(s),
                                     "f_sym": at - gt, "f_anti": ab - gb})
    positives = sum(1 for r in rows
                    if r["f_sym"] >= 0.0 or r["f_anti"] >= 0.0)
    return {"samples": len(rows), "positives": int(positives),
            "worst_approach_to_zero": float(worst),
            "the_tube_side_is_always_negative": bool(tube_ok),
            "the_exterior_side_is_always_positive": bool(amb_ok),
            "there_is_no_growing_mode": bool(positives == 0 and tube_ok
                                             and amb_ok),
            "the_answer": ("no, with the balls removed as well — and here the "
                           "quadratic form settles it without any sweep")}


def measure_the_higher_multipoles_decouple(
        throat: NeckThroat = WORKING_NECK,
        ells: Sequence[int] = (1, 2, 3, 5),
        sigmas: Sequence[float] = (0.1, 1.0, 10.0)) -> Dict[str, object]:
    """PR #261's monopole truncation was never a *stability* limitation.

    A tube with one transverse channel presents a single number at each mouth,
    so it drives the ``ℓ = 0`` projection of the boundary data and nothing else.
    The ``ℓ ≥ 1`` sectors do not see the tube at all: they are the exterior's
    own modes, with its own boundary condition, and their DtN is positive.  So
    they can neither be driven nor go unstable, and dropping them costs
    amplitudes rather than the answer.

    That is worth stating precisely because PR #261 listed monopole matching as
    its main open approximation, quantified by PR #250's ``(a/d)^ℓ`` screening.
    The screening bounds how much *response* is missed.  For the question this
    arc gated the roadmap on, the higher channels were never in play.
    """
    rows = []
    for l in ells:
        for s in sigmas:
            lam = -(s ** 2)
            rows.append({"ell": int(l), "sigma": float(s),
                         "dtn": exterior_dtn(throat.radius, lam, int(l)),
                         "monopole_dtn": exterior_dtn(throat.radius, lam, 0)})
    for r in rows:
        r["stiffer_than_monopole"] = bool(r["dtn"] > r["monopole_dtn"])
    return {"rows": rows,
            "every_channel_is_positive": bool(all(r["dtn"] > 0.0
                                                  for r in rows)),
            "every_channel_is_stiffer_than_the_monopole": bool(
                all(r["stiffer_than_monopole"] for r in rows)),
            "smallest_ratio_to_the_monopole": float(min(
                r["dtn"] / r["monopole_dtn"] for r in rows)),
            "why_it_matters": ("a one-channel tube drives only ℓ = 0, so the "
                               "higher sectors are spectators; the (a/d)^ℓ "
                               "screening bounds missed amplitude, not the "
                               "stability answer")}


def measure_what_the_fixed_ambient_cost(
        radii: Sequence[float] = (0.02, 0.05, 0.15, 0.35),
        lams: Sequence[float] = (0.0, -4.0)) -> Dict[str, object]:
    """How wrong PR #261's "spheres in a fixed ambient" actually was.

    That round used ``f(a)G(a)`` — a double-shell average taken with the balls
    still in place.  This one uses ``1/N₀(a,λ)``, the inverse exterior map with
    them gone.  The two agree to ``3.8e-03`` at the working radius ``a = 0.05``
    (``1.3e-04`` at ``a = 0.02``, ``λ = 0``, where the ball occupies less of the
    sphere) and part company at large radius and deep ``λ``, reaching ``11%`` at
    ``a = 0.35``, ``λ = −4``.

    The error tracks ``a`` and ``|λ|`` together, which is the expected shape: it
    is the fraction of the sphere wrongly left in, weighted by how much of the
    field's support sits there.

    So the approximation was good where it was used.  Quoting the number is the
    point: PR #261 flagged this as its main limitation, and a limitation with a
    measured size is a different object from one with a caveat.

    The same measurement prices the approximation this round has **not**
    removed.  The reduced ``2×2`` keeps a **single-scattering** cross term
    ``f(a)²G(d)``; the exact two-ball exterior would add the series in which a
    monopole on one boundary sphere drives the other and is driven back.  Its
    expansion parameter is ``cross/self``, which the model computes anyway and
    which comes out as ``0.8·(a/d)`` — ``0.031`` at the working point, so the
    neglected terms enter at ``~1e-03`` there and never exceed ``0.16`` anywhere
    sampled.  A correction that small cannot flip a sign, which is the only
    thing the reduced model is being asked to decide — and the theorem does not
    go through the reduced model at all.
    """
    rows = []
    for a in radii:
        for lam in lams:
            kept = regular_radial(a, lam) * mouth_green(a, lam)
            removed = 1.0 / exterior_dtn_monopole(a, lam)
            rows.append({"radius": float(a), "lambda": float(lam),
                         "fixed_ambient": float(kept),
                         "balls_removed": float(removed),
                         "ratio": float(kept / removed),
                         "relative_error": float(abs(kept - removed)
                                                 / abs(removed))})
    # the one approximation left in the reduced 2x2: the cross term is
    # single-scattering, and its expansion parameter is cross/self
    scatter = []
    for a, d in ((0.02, 1.3), (0.05, 1.3), (0.15, 1.3), (0.35, 1.3),
                 (0.35, 0.8)):
        for lam in lams:
            t = NeckThroat(separation=d, length=0.9, area=FOUR_PI, radius=a)
            self_term = 1.0 / exterior_dtn_monopole(a, float(lam))
            cross = (regular_radial(a, float(lam)) ** 2
                     * mouth_green(d, float(lam)))
            scatter.append({"radius": float(a), "separation": float(d),
                            "lambda": float(lam),
                            "cross_over_self": float(cross / self_term),
                            "over_a_over_d": float(cross / self_term
                                                   / (a / d)),
                            "neglected_order": float((cross / self_term) ** 2)})
    small = [r for r in rows if r["radius"] <= 0.05]
    working = [r for r in scatter
               if r["radius"] == 0.05 and r["separation"] == 1.3]
    return {"rows": rows,
            "single_scattering": scatter,
            "worst_error_at_the_working_radius": float(
                max(r["relative_error"] for r in small)),
            "worst_error_overall": float(max(r["relative_error"]
                                             for r in rows)),
            "the_error_grows_with_the_radius": bool(
                rows[-1]["relative_error"] > rows[0]["relative_error"]),
            "pr_261_was_right_where_it_looked": bool(
                max(r["relative_error"] for r in small) < 1e-2),
            "cross_over_self_scales_like_a_over_d": bool(
                max(r["over_a_over_d"] for r in scatter if r["lambda"] == 0.0)
                / min(r["over_a_over_d"] for r in scatter
                      if r["lambda"] == 0.0) < 1.3),
            "neglected_scattering_at_the_working_point": float(
                max(r["neglected_order"] for r in working)),
            "neglected_scattering_worst": float(
                max(r["neglected_order"] for r in scatter)),
            "the_neglected_series_cannot_reach_the_sign": bool(
                max(r["cross_over_self"] for r in scatter) < 0.5)}


def measure_the_soft_mode_is_forced_by_the_two_ends(
        lengths: Sequence[float] = (0.3, 0.9, 2.0, 3.0, 4.0, 6.0, 8.0),
        throat: NeckThroat = WORKING_NECK) -> Dict[str, object]:
    """Existence, by the same kind of argument that gave non-existence.

    The growing mode was ruled out by a **sign**: negative tube against positive
    exterior, no zero.  The soft mode is ruled *in* by the **two ends of the
    gap**, and no scan is needed for that either:

    * at ``λ → 0⁺`` the tube's symmetric channel is ``cot(kL/2)/(𝒜k) ≈
      2/(𝒜λL) → +∞`` while the exterior block stays finite, so ``F_sym > 0``;
    * at ``λ → 1⁻`` the exterior map ``N₀ → 0⁺`` — the free ESU threshold is
      exactly where the regular exterior solution goes flat, ``ψ ∝
      sin(ω(π−χ))/sin χ → 1`` — so the ambient diagonal ``1/N₀ → +∞`` and
      ``F_sym → −∞``.  The rate is closed-form: writing ``ω = 1 − δ``, the
      first-order expansion gives

          ``N₀ → 2π(π − a + sin a cos a)·(1 − λ)``

      which the integrator reproduces to ``7e-07`` at ``a = 1`` and ``3e-05``
      at ``a = 0.05`` (the residual is the ``O(δ²)`` truncation, and it shrinks
      as the edge is approached).  At ``a → 0`` the slope is ``2π²``.

    A continuous function running from ``+∞`` to ``−∞`` has a zero.  **The soft
    mode is not a numerical find; it is forced**, and it is forced by the
    exterior stiffness vanishing at the gap edge rather than by anything about
    the tube.

    The antisymmetric channel starts at ``−L/(2𝒜) < 0`` instead, which is why it
    has no state — *while the tube is short enough*.

    A correction to PR #261
    ───────────────────────
    That round reported "exactly one state below the gap, in the symmetric
    channel".  That is true for the lengths it used, and **not** structural.
    ``A_sym`` has poles at ``λ = (2πj/L)²`` and ``A_anti`` at
    ``λ = (π(2j−1)/L)²``; the first enters the gap at ``L = π``, and each tube
    harmonic that falls in brings a genuine extra state just above it.  At
    ``L = 8`` there are three.

    The bad reading is easy to get and worth naming: a pole is a **sign change
    with no zero**, so crossing-counting alone reports states that are not
    there.  Classifying by the residual separates them by fifteen orders of
    magnitude.  None of this touches the stability answer — every one of these
    states is positive, as the quadratic form says they must be.
    """
    rows = []
    for length in lengths:
        t = NeckThroat(separation=throat.separation, length=float(length),
                       area=throat.area, radius=throat.radius,
                       interior_mass=throat.interior_mass)
        b = t.bound_states(n=6000)
        poles = t.tube_poles_in_the_gap()
        rows.append({
            "length": float(length),
            "f_sym_at_zero": t.channel_functions(1e-9)[0],
            "f_sym_at_the_gap": t.channel_functions(1.0 - 1e-9)[0],
            "f_anti_at_zero": t.channel_functions(1e-9)[1],
            "symmetric_roots": [float(x) for x in b["symmetric_roots"]],
            "antisymmetric_roots": [float(x)
                                    for x in b["antisymmetric_roots"]],
            "spurious_pole_crossings": (len(b["symmetric_poles"])
                                        + len(b["antisymmetric_poles"])),
            "tube_harmonics_in_the_gap": (len(poles["symmetric"])
                                          + len(poles["antisymmetric"])),
            "states": int(b["states_below_the_gap"]),
        })
    a = float(throat.radius)
    slope = 2.0 * math.pi * (math.pi - a + math.sin(a) * math.cos(a))
    edge = []
    for x in (0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999):
        n0 = exterior_dtn_monopole(a, float(x))
        edge.append({"lmbda": float(x), "N0": float(n0),
                     "N0_over_one_minus_lambda": float(n0 / (1.0 - x)),
                     "closed_form_slope": float(slope),
                     "relative": float(abs(n0 / (1.0 - x) - slope) / slope)})
    return {
        "rows": rows,
        "gap_edge": edge,
        "the_symmetric_channel_starts_positive": bool(
            all(r["f_sym_at_zero"] > 0.0 for r in rows)),
        "and_ends_negative": bool(all(r["f_sym_at_the_gap"] < 0.0
                                      for r in rows)),
        "so_a_symmetric_state_is_forced": bool(
            all(len(r["symmetric_roots"]) >= 1 for r in rows)),
        "the_antisymmetric_channel_starts_negative": bool(
            all(r["f_anti_at_zero"] < 0.0 for r in rows)),
        "gap_edge_slope_error": float(edge[-1]["relative"]),
        "gap_edge_error_ratio_per_decade": float(edge[-1]["relative"]
                                                 / edge[-2]["relative"]),
        "the_exterior_stiffness_vanishes_linearly_at_the_gap_edge": bool(
            edge[-1]["N0"] > 0.0
            and abs(edge[-1]["relative"] / edge[-2]["relative"] - 0.1) < 0.01),
        "short_tubes_hold_exactly_one_state": bool(
            all(r["states"] == 1 for r in rows if r["length"] < math.pi)),
        "long_tubes_hold_more": bool(
            any(r["states"] > 1 for r in rows if r["length"] > math.pi)),
        "every_extra_state_follows_a_tube_harmonic": bool(
            all(r["states"] == 1 + r["tube_harmonics_in_the_gap"]
                for r in rows)),
        "pole_crossings_that_would_have_been_miscounted": int(
            sum(r["spurious_pole_crossings"] for r in rows)),
        "what_pr_261_got_wrong": ("'exactly one state below the gap' is a "
                                  "statement about L < π, not a structural "
                                  "one; it holds for every length that round "
                                  "used and fails above L = π"),
    }


def measure_the_soft_mode_survives_the_removal(
        radii: Sequence[float] = (0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
        ) -> Dict[str, object]:
    """The one positive result of PR #261, re-measured with the balls gone.

    ``λ₀`` moves by ``1.2e-04`` relative at ``a = 0.02`` and still tends to
    ``8πa/(𝒜L)``: two mouth capacitances restoring a tube of volume ``𝒜L``.  The
    capacitance reading is if anything cleaner here, since ``N₀ → 4πa`` is now
    the exterior map's own small-ball limit rather than an identification made
    after the fact.
    """
    rows = []
    for a in radii:
        neck = NeckThroat(separation=1.3, length=0.9, area=FOUR_PI,
                          radius=float(a))
        lam_neck = neck.bound_states()["symmetric"]
        lam_fixed = neck.fixed_ambient().bound_states()["symmetric"]
        closed = 8.0 * math.pi * a / (FOUR_PI * 0.9)
        rows.append({"radius": float(a),
                     "lambda_0_neck": float(lam_neck),
                     "lambda_0_fixed_ambient": float(lam_fixed),
                     "closed_form": float(closed),
                     "shift_from_removal": float(
                         abs(lam_neck - lam_fixed) / lam_fixed),
                     "ratio_to_closed_form": float(lam_neck / closed)})
    small = [r for r in rows if r["radius"] <= 0.02]
    return {"rows": rows,
            "every_one_is_positive_and_below_the_gap": bool(all(
                0.0 < r["lambda_0_neck"] < 1.0 for r in rows)),
            "worst_shift_from_removal": float(max(r["shift_from_removal"]
                                                  for r in rows)),
            "shift_at_the_working_radius": float(
                [r for r in rows if r["radius"] == 0.05][0]
                ["shift_from_removal"]),
            "worst_closed_form_error": float(max(
                abs(1.0 - r["ratio_to_closed_form"]) for r in small)),
            "the_capacitance_formula_survives": bool(max(
                abs(1.0 - r["ratio_to_closed_form"]) for r in small) < 0.02),
            "removing_the_balls_barely_moves_it": bool(
                max(r["shift_from_removal"] for r in rows) < 0.02)}
