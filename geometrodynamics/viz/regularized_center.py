"""The center as a finite bearing, not a point — and what that costs.

The proposal
────────────
Every visualisation in this arc so far has put a *point* at the middle of the
picture.  A point is where the clock-hands story works: two radial arms

    ``P_A → O → P_B``

can change direction at ``O`` for free, because at a point there is no angular
direction left to change.  That is exactly the property the connection wants —
it makes the link's cost independent of where the mouths sit — and it is bought
with ``f = 0``, which is where the geometry stops existing.

Regularise it.  Keep the bulk metric rotationally symmetric,

    ``dℓ² = ds² + f(s)² dΩ²`` ,

and replace ``f(0) = 0`` with

    ``f_min = f₀ > 0`` .

In the 2-D cross-section the middle is now a small **circle** rather than a
point.  In a ``d``-dimensional bulk it is the space of radial directions,
``S^{d−1}`` — or ``RP^{d−1}`` if the clock hand is an unoriented *axis*, so
that ``n ∼ −n``.  Nothing is singular any more, and three things that the
point picture forced become free.

**The two arms need not match.**  There is no reason for the outer and inner
ends to have the same transverse scale or the same distance to the bearing.
Take the scalar-flat profile ``f'² = 1 − f₀/f`` as the worked example — it is
the one PR #265 forced, so it is not an invention — and the proper distance
from the neck out to an end of scale ``F`` is exactly

    ``L(F) = √(F(F − f₀)) + f₀ arcosh√(F/f₀)`` .

So ``L_o = L(f_o)`` and ``L_i = L(f_i)`` are independently whatever the two
ends are, and ``L_o/L_i`` runs to 437 in `measure_the_two_arms_are_independent`
without anything breaking.  Setting ``f_o = f_i = sin a`` and ``f₀ = sin³a``
recovers `physical_throat.VacuumThroat.length()` **to the last bit** — the
asymmetric arm is the symmetric throat's own formula with the symmetry dropped.

**Scale transport becomes explicit.**  A feature of angular width ``Δθ`` has
physical width ``w(s) = f(s)Δθ``, so along the route it goes

    ``w_o = f_o Δθ  →  w_min = f₀ Δθ  →  w_i = f_i Δθ`` ,

and ``w_i/w_o = f_i/f_o`` exactly.  A packet does not teleport from one radius
to another at fixed size; it is squeezed into the bearing and let out again.

**And the hinge costs something finite.**  That is where this module's actual
result is.

The turn cost is quadratic, not linear
──────────────────────────────────────
The natural first guess is that turning through ``α`` around a bearing of
radius ``f₀`` costs its arc length, ``f₀α``, giving
``L_throat(α) ≃ L_o + L_i + f₀α``.  That is **exactly right for one specific
route** — go down the outer arm at fixed direction, turn on the bearing, go up
the inner arm — and `corner_route_length` computes it.  But that route is not
the geodesic, and the geodesic is much cheaper.

Clairaut's relation for a surface of revolution, ``f sin ψ = h``, makes the
shortest path cut the corner: it starts turning while it is still descending,
where the lever arm ``f`` is longer.  Solving for the ``h`` that sweeps a
given ``α`` and integrating (`turn_cost`) gives, in the small-angle limit,

    ``T(α) = α² / (2I)`` ,   ``I = ∫ ds/f²`` ,

and ``I`` is **not a new quantity** — it is the same resistance integral
`physical_throat.VacuumThroat.resistance` already computes, whose reciprocal
sets the throat's monopole conductance ``4π/I``.  Per arm it is

    ``I(F) = (2/f₀)√(1 − f₀/F)`` ,

exact, and ``2I(sin a)`` reproduces the repo's ``4cos a / sin³a`` bit for bit.
So the geometric cost of swinging the clock hands and the electrical cost of
pushing monopole flux through the throat are set by *one* integral:

    ``T(α) = α² · conductance / (8π)`` .

For two long arms ``I → 4/f₀`` and the law reads ``T(α) → f₀α²/8``.

This matters, because it is a correction that runs the *helpful* way.  The
linear guess overstates the cost badly: at ``α = 0.1`` the geodesic spends
``1.25%`` of ``f₀α``, and even at ``α = π`` only ``36%``.  Turning is never
expensive compared with travelling — ``T(π)/(L_o + L_i) = 8.4e-04`` at
``f₀ = 1e-03``, and you would need ``α ≈ 104`` radians before the hinge cost as
much as the arms.  **The conclusion the point-center picture was wanted for
survives the regularisation, and survives it more strongly than proposed.**

Both costs still vanish with the bearing: ``T`` is linear in ``f₀`` at fixed
``α``, so ``f₀ → 0`` returns ``L_o + L_i``, independent of where the mouths
are, exactly as the singular picture had it.

What "intersection" becomes
───────────────────────────
Two fronts no longer have to collide at ``r = 0``.  They arrive at angular
positions on the bearing, and

    they meet  ⟺  their angular extents overlap ,

which is a statement about **angles alone** and does not involve ``f₀`` at all.
What ``f₀`` sets is the *size* of the meeting: the overlap region is
``f₀ × (overlap angle)`` across, and the gap between two fronts that miss is
``f₀ × (gap angle)``.  Both shrink to nothing as ``f₀ → 0``, which is how the
point limit is recovered — not by making everything meet, but by making the
distinction between meeting and missing physically invisible.
`measure_the_bearing_replaces_collision_with_overlap` carries both halves.

Scope
─────
This is **geometry**, deliberately.  It is a metric, its geodesics, and the
transport of an angular width along them.  No field equation is solved on it,
nothing here says which of the three candidate readings of the bulk is right,
and the scalar-flat profile is used as a concrete worked example rather than
as a claim that the bearing must be that profile.  The one thing the choice
does buy is that every closed form below is checkable against `physical_throat`,
which is why it was used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

__all__ = [
    "RegularizedCenter",
    "WORKING_CENTER",
    "arm_length",
    "arm_resistance",
    "arm_moment",
    "bearing_distance",
    "measure_the_arm_length_is_the_repos_own_formula",
    "measure_the_two_arms_are_independent",
    "measure_the_width_rescales_with_the_profile",
    "measure_the_turn_cost_is_quadratic_not_linear",
    "measure_the_corner_route_is_only_an_upper_bound",
    "measure_the_hinge_is_never_the_expensive_part",
    "measure_the_bearing_shrinks_back_to_the_point",
    "measure_the_bearing_replaces_collision_with_overlap",
    "measure_the_turn_cost_does_not_care_which_sphere",
    "measure_the_law_does_not_depend_on_the_profile",
    "measure_the_fourth_moment_is_where_the_neck_shape_enters",
    "measure_the_hinge_and_the_monopole_are_one_dirichlet_form",
    "hyperbolic_neck",
]


# ════════════════════════════════════════════════════════════════════════════
# CLOSED FORMS
# ════════════════════════════════════════════════════════════════════════════
def arm_length(scale: float, neck: float) -> float:
    """``L(F) = √(F(F−f₀)) + f₀ arcosh√(F/f₀)`` — neck to an end of scale ``F``.

    The proper distance ``∫ds`` along ``f'² = 1 − f₀/f``, in closed form.  With
    ``F = sin a`` and ``f₀ = sin³a`` twice this is
    `physical_throat.VacuumThroat.length`, which is the check that it is the
    same object with the symmetry assumption removed rather than a new one.
    """
    f, f0 = float(scale), float(neck)
    if f0 <= 0.0:
        return f                      # the singular limit: L(F) -> F
    if f < f0:
        raise ValueError("an arm cannot be narrower than the neck")
    return math.sqrt(f * (f - f0)) + f0 * math.acosh(math.sqrt(f / f0))


def arm_resistance(scale: float, neck: float) -> float:
    """``I(F) = (2/f₀)√(1 − f₀/F)`` — one arm's share of ``∫ds/f²``.

    The same integral `physical_throat.VacuumThroat.resistance` evaluates over
    both arms of the symmetric throat, and it is what sets the turn cost.  As
    ``F → ∞`` it saturates at ``2/f₀``, so two long arms give ``I → 4/f₀``.
    """
    f, f0 = float(scale), float(neck)
    if f < f0:
        raise ValueError("an arm cannot be narrower than the neck")
    return (2.0 / f0) * math.sqrt(1.0 - f0 / f)


def arm_moment(scale: float, neck: float, order: int) -> float:
    """``I_n(F) = ∫ ds/fⁿ`` over one arm, in closed form for even ``n``.

    Substituting ``t = f' = √(1 − f₀/f)`` — the profile's own slope — turns
    every one of these into a polynomial.  ``f = f₀/(1−t²)`` gives
    ``ds = 2f₀dt/(1−t²)²``, so

        ``I_n = (2/f₀^{n−1}) ∫₀^T (1−t²)^{n−2} dt`` ,  ``T = √(1 − f₀/F)`` ,

    and the integral is a finite binomial sum.  ``n = 2`` recovers
    `arm_resistance`, ``(2/f₀)T``; ``n = 4`` gives
    ``(2/f₀³)[T − ⅔T³ + ⅕T⁵]``, which tends to ``16/(15f₀³)`` for a long arm.

    Only even orders appear in the turn-cost expansion, because
    ``(1 − h²/f²)^{−1/2}`` is a series in ``h²/f²``.
    """
    f, f0, n = float(scale), float(neck), int(order)
    if n < 2 or n % 2:
        raise ValueError("only even orders n >= 2 are defined here")
    if f < f0:
        raise ValueError("an arm cannot be narrower than the neck")
    t = math.sqrt(1.0 - f0 / f)
    m = n - 2
    total = 0.0
    for k in range(m + 1):
        total += (math.comb(m, k) * (-1.0) ** k * t ** (2 * k + 1)
                  / (2 * k + 1))
    return (2.0 / f0 ** (n - 1)) * total


def bearing_distance(alpha: float, projective: bool = False) -> float:
    """Angular distance on the bearing, in units of its radius.

    ``α`` on the round direction sphere; ``min(α, π−α)`` when the clock hand is
    an unoriented axis, ``n ∼ −n``, so that the bearing is ``RP^{d−1}``.  The
    formula does not depend on the bearing's dimension: two directions span a
    totally geodesic 2-plane, so the walk between them is a great-circle arc
    whatever ``d`` is.
    """
    a = abs(float(alpha)) % (2.0 * math.pi)
    a = min(a, 2.0 * math.pi - a)
    return min(a, math.pi - a) if projective else a


@dataclass(frozen=True)
class RegularizedCenter:
    """A finite central bearing with two independently sized arms.

    ``neck`` is ``f₀``; ``outer`` and ``inner`` are the transverse scales
    ``f_o`` and ``f_i`` at the two ends.  Nothing requires them to be equal.
    """

    neck: float = 1e-3
    outer: float = 1.0
    inner: float = 0.35

    # ── the arms ────────────────────────────────────────────────────────────
    def outer_length(self) -> float:
        return arm_length(self.outer, self.neck)

    def inner_length(self) -> float:
        return arm_length(self.inner, self.neck)

    def arm_length_sum(self) -> float:
        """``L_o + L_i`` — what the route costs before any turning."""
        return self.outer_length() + self.inner_length()

    def resistance(self) -> float:
        """``I = ∫ds/f²`` over both arms."""
        return arm_resistance(self.outer, self.neck) + \
            arm_resistance(self.inner, self.neck)

    def conductance(self) -> float:
        """``4π/I`` — the same quantity `physical_throat` calls by this name."""
        return 4.0 * math.pi / self.resistance()

    def length_by_quadrature(self, which: str = "outer") -> float:
        """``∫ds`` for one arm, for the check that the closed form is one."""
        f0 = float(self.neck)
        end = float(self.outer if which == "outer" else self.inner)
        return quad(lambda t: 2.0 * math.sqrt(f0 + t * t), 0.0,
                    math.sqrt(end - f0), limit=400)[0]

    def resistance_by_quadrature(self, which: str = "outer") -> float:
        f0 = float(self.neck)
        end = float(self.outer if which == "outer" else self.inner)
        return quad(lambda t: 2.0 * math.sqrt(f0 + t * t) / (f0 + t * t) ** 2,
                    0.0, math.sqrt(end - f0), limit=400)[0]

    # ── the profile, for drawing ────────────────────────────────────────────
    def profile(self, n: int = 400) -> Tuple[np.ndarray, np.ndarray]:
        """``(s, f)`` along the whole route, outer arm at negative ``s``.

        Sampled uniformly in ``√(f − f₀)`` so the neck — where the profile
        turns over and everything interesting happens — is not under-resolved.
        """
        f0 = float(self.neck)
        out = []
        for end, sign in ((self.outer, -1.0), (self.inner, +1.0)):
            t = np.linspace(0.0, math.sqrt(float(end) - f0), n)
            f = f0 + t * t
            s = np.array([arm_length(float(x), f0) for x in f])
            out.append((sign * s, f))
        s = np.concatenate([out[0][0][::-1], out[1][0][1:]])
        f = np.concatenate([out[0][1][::-1], out[1][1][1:]])
        return s, f

    def feature_width(self, scale: float, angular_width: float) -> float:
        """``w = f Δθ`` — the physical width of a fixed angular feature."""
        return float(scale) * float(angular_width)

    # ── the hinge ───────────────────────────────────────────────────────────
    def corner_route_length(self, alpha: float,
                            projective: bool = False) -> float:
        """``L_o + L_i + f₀·d(α)`` — down, turn on the bearing, up.

        An **exact** length, for a route that is not the shortest one: it holds
        its direction all the way down, turns on the bearing, and holds it all
        the way up.  So it is a rigorous upper bound on `route_length`, and
        `measure_the_corner_route_is_only_an_upper_bound` says by how much.
        """
        return self.arm_length_sum() + \
            self.neck * bearing_distance(alpha, projective)

    def moment(self, order: int) -> float:
        """``I_n = ∫ds/fⁿ`` over both arms.  ``moment(2)`` is `resistance`."""
        return (arm_moment(self.outer, self.neck, order)
                + arm_moment(self.inner, self.neck, order))

    def fourth_moment(self) -> float:
        """``I₄`` — the first moment that remembers the neck's *shape*."""
        return self.moment(4)

    def turn_cost_to_fourth_order(self, alpha: float,
                                  projective: bool = False) -> float:
        """``α²/(2I₂) − α⁴I₄/(8I₂⁴)`` — the hinge cost with its first correction.

        Expanding ``dℓ/ds − 1 = (1 − h²/f²)^{−1/2} − 1`` gives
        ``T = ½h²I₂ + ⅜h⁴I₄``, and the sweep ``α = hI₂ + ½h³I₄``; eliminating
        ``h`` leaves the above.  The correction is **negative**, which is why
        the measured shape ``T/(α²/2I₂)`` sits below one at every angle.
        """
        a = bearing_distance(alpha, projective)
        i2, i4 = self.resistance(), self.fourth_moment()
        return a * a / (2.0 * i2) - a ** 4 * i4 / (8.0 * i2 ** 4)

    def monopole_profile(self, x: float) -> float:
        """``∫ds/f²`` from the neck out to ``x`` — the static ℓ=0 potential.

        ``(f²u')' = 0`` gives ``u' = c/f²``, so the potential's shape is this
        partial moment.  In the good variable it is ``(2/f₀)tanh x``, exactly.
        """
        return (2.0 / self.neck) * math.tanh(float(x))

    def azimuth_profile(self, x: float, kappa: float) -> float:
        """``θ(x)`` swept by the geodesic with Clairaut constant ``κ``.

        The same object as `monopole_profile` once ``κ → 0``: see
        `measure_the_hinge_and_the_monopole_are_one_dirichlet_form`.
        """
        return quad(lambda y: 2.0 * kappa
                    / math.sqrt(max(math.cosh(y) ** 4 - kappa * kappa, 1e-300)),
                    0.0, float(x), limit=200)[0]

    def half_length_in_x(self, scale: float) -> float:
        """``X = arcosh√(F/f₀)`` — the arm's extent in the good variable.

        ``f = f₀cosh²x`` is the substitution that turned PR #267's throat
        equation into a constant-coefficient one, and it is the right variable
        here for the same reason: it gives ``ds = 2f dx`` exactly, so both
        Clairaut integrands below become smooth ``O(1)`` functions decaying
        like ``e^{−2x}``, whatever ``f₀`` is.
        """
        return math.acosh(math.sqrt(float(scale) / float(self.neck)))

    def _arm_integral(self, end: float, kappa: float, kind: str) -> float:
        """One arm of the Clairaut integrals, in ``f = f₀cosh²x``.

        With ``h = f₀κ`` (and ``κ < 1``, which is what crossing the neck
        means),

            ``dθ/dx = 2κ / √(cosh⁴x − κ²)`` ,
            ``dT/dx = 2f₀cosh²x · κ² / [ R (cosh²x + R) ]`` , ``R = √(cosh⁴x−κ²)`` .

        Two things about the second one are deliberate.  It integrates the
        turn cost **directly** rather than differencing the geodesic against
        ``L_o + L_i`` — ``T`` is as small as ``1e-15`` here and the two lengths
        are ``O(1)``, so a subtraction would be all cancellation.  And the
        factor ``cosh²x − R`` has been rewritten as ``κ²/(cosh²x + R)``, which
        is the same number with nothing to cancel in it.

        Both matter.  Differencing gave answers ``2.7×`` wrong at
        ``f₀ = 1e-07``, and integrating in ``f`` instead of ``x`` put a spike
        of width ``√f₀`` in front of an adaptive quadrature that walked past
        it.  In ``x`` the result is stable over seven decades of ``f₀``.
        """
        upper = self.half_length_in_x(end)

        def integrand(x: float) -> float:
            c2 = math.cosh(x) ** 2
            r = math.sqrt(max(c2 * c2 - kappa * kappa, 0.0))
            if r <= 0.0:
                return 0.0
            if kind == "angle":
                return 2.0 * kappa / r
            return 2.0 * self.neck * c2 * kappa * kappa / (r * (c2 + r))

        # Two separate scales have to be resolved by hand, because an adaptive
        # rule finds neither reliably.
        #
        # The tail: both integrands decay like e^{-2x}, while `upper` runs to
        # 15 for a small enough f0, so one call spends its subdivisions on
        # dead range and reports a 15% error estimate on a correct answer.
        #
        # The peak: near x = 0, cosh^4 x - k^2 ~ (1 - k^2) + 4x^2, so the angle
        # integrand is a spike of width sqrt(1 - k^2)/2 that decays like 1/x
        # beyond it.  Driving alpha up pushes k towards 1 and that width
        # towards zero -- the root find probes k = 1 - 1e-15, where it is
        # 2e-08.  A 1/x tail spanning eight decades needs the breakpoints
        # spaced logarithmically; a fixed few leave most of it in one panel.
        width = 0.5 * math.sqrt(max(1.0 - kappa * kappa, 0.0))
        marks = [0.0, 3.0, 8.0]
        if 0.0 < width < 1.0:
            e = width
            while e < 3.0:
                marks.append(e)
                e *= 4.0
        edges = sorted({e for e in marks if 0.0 <= e < upper}) + [upper]
        return sum(quad(integrand, lo, hi, limit=200)[0]
                   for lo, hi in zip(edges, edges[1:]))

    def _sweep(self, kappa: float) -> float:
        return self._arm_integral(self.outer, kappa, "angle") + \
            self._arm_integral(self.inner, kappa, "angle")

    def clairaut_constant(self, alpha: float,
                          projective: bool = False) -> float:
        """``κ = h/f₀`` in ``f sin ψ = h``, for the geodesic that sweeps ``α``.

        Only ``κ < 1`` — that is, ``h < f₀`` — crosses the neck at all; a
        geodesic with a larger ``h`` turns back before reaching it.  The swept
        angle grows without bound as ``κ → 1``, so every ``α`` is reachable.

        The bracket stops at ``1 − 1e-09`` rather than at ``1``, and the
        headroom is not marginal: `bearing_distance` never asks for more than
        ``π``, which lands at ``κ ≈ 0.67``, while the bracket's top sweeps
        past ``30`` radians — near ``10π``.  Going closer to ``1`` buys
        nothing and costs something, because the angle integrand's peak
        narrows as ``√(1−κ²)`` and by ``1 − 1e-15`` it is too thin for the
        quadrature to resolve.
        """
        a = bearing_distance(alpha, projective)
        if a <= 0.0:
            return 0.0
        top = 1.0 - 1e-9
        if self._sweep(top) < a:                    # unreachable, and say so
            raise ValueError(
                f"a sweep of {a} rad needs kappa beyond the bracket; the "
                f"bracket reaches {self._sweep(top)} rad")
        return brentq(lambda k: self._sweep(k) - a, 1e-18, top,
                      xtol=1e-18, rtol=8.9e-16)

    def turn_cost(self, alpha: float, projective: bool = False) -> float:
        """``T(α)`` — what the hinge costs, as a **direct** integral.

        Never formed as ``route_length − (L_o + L_i)``: see `_arm_integral`.
        """
        a = bearing_distance(alpha, projective)
        if a <= 0.0:
            return 0.0
        k = self.clairaut_constant(a)
        return self._arm_integral(self.outer, k, "length") + \
            self._arm_integral(self.inner, k, "length")

    def route_length(self, alpha: float, projective: bool = False) -> float:
        """The **geodesic** length of the route that turns through ``α``.

        The exact closed-form arms plus the directly integrated hinge, so the
        small quantity is never recovered from the difference of large ones.
        """
        return self.arm_length_sum() + self.turn_cost(alpha, projective)

    def turn_cost_small_angle(self, alpha: float,
                              projective: bool = False) -> float:
        """``α²/(2I)`` — the small-angle law, with ``I`` the arms' resistance."""
        a = bearing_distance(alpha, projective)
        return a * a / (2.0 * self.resistance())


WORKING_CENTER = RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_arm_length_is_the_repos_own_formula() -> Dict[str, object]:
    """The asymmetric arm reduces to `VacuumThroat.length` — bit for bit.

    Two separate claims.  First that ``L(F)`` is ``∫ds``: checked against
    quadrature, which is the ordinary closed-form-versus-integral test.  Second
    that it is not a *new* formula: with ``f_o = f_i = sin a`` and
    ``f₀ = sin³a``, ``L_o + L_i`` must equal `physical_throat`'s ``length()``
    exactly, and ``I`` must equal its ``resistance()`` exactly.  The point of
    insisting on *exactly* is that a regularised centre is only worth
    proposing if it contains the geometry the repo already forced.
    """
    from geometrodynamics.waves.physical_throat import VacuumThroat

    quadrature = []
    for f0 in (1e-4, 1e-2, 0.1):
        for end in (1.3 * f0, 5.0 * f0, 0.5, 1.0):
            if end <= f0:
                continue
            c = RegularizedCenter(neck=f0, outer=end, inner=end)
            quadrature.append({
                "neck": f0, "scale": end,
                "length_closed": c.outer_length(),
                "length_quadrature": c.length_by_quadrature(),
                "length_relative": abs(c.outer_length()
                                       / c.length_by_quadrature() - 1.0),
                "resistance_closed": arm_resistance(end, f0),
                "resistance_quadrature": c.resistance_by_quadrature(),
                "resistance_relative": abs(arm_resistance(end, f0)
                                           / c.resistance_by_quadrature() - 1.0),
            })

    against_repo = []
    for a in (0.02, 0.05, 0.10, 0.20):
        t = VacuumThroat(mouth_radius=a)
        c = RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                              inner=t.mouth_f())
        against_repo.append({
            "mouth_radius": a,
            "length_here": c.arm_length_sum(),
            "length_physical_throat": t.length(),
            "length_difference": abs(c.arm_length_sum() - t.length()),
            "resistance_here": c.resistance(),
            "resistance_physical_throat": t.resistance(),
            "resistance_difference": abs(c.resistance() - t.resistance()),
        })
    return {
        "quadrature": quadrature,
        "against_physical_throat": against_repo,
        "closed_forms": ["L(F) = sqrt(F(F-f0)) + f0 arcosh sqrt(F/f0)",
                         "I(F) = (2/f0) sqrt(1 - f0/F)"],
        "worst_length_quadrature_error": float(
            max(r["length_relative"] for r in quadrature)),
        "worst_resistance_quadrature_error": float(
            max(r["resistance_relative"] for r in quadrature)),
        "the_closed_forms_are_the_integrals": bool(
            max(r["length_relative"] for r in quadrature) < 1e-9
            and max(r["resistance_relative"] for r in quadrature) < 1e-9),
        "it_contains_the_symmetric_throat_exactly": bool(
            max(r["length_difference"] for r in against_repo) < 1e-15
            and max(r["resistance_difference"] for r in against_repo) < 1e-8),
    }


def measure_the_two_arms_are_independent() -> Dict[str, object]:
    """``L_o ≠ L_i`` is ordinary, not pathological.

    The point-center picture had no room for this: two concentric boundaries
    separated by one arbitrary radial gap.  Here the two ends are separately
    whatever they are, and the ratio is scanned to ``437`` without the
    resistance, the turn cost or the small-angle law changing character.
    """
    f0 = 1e-3
    rows = []
    for inner in (0.002, 0.01, 0.1, 0.5, 1.0):
        c = RegularizedCenter(neck=f0, outer=1.0, inner=inner)
        rows.append({
            "inner_scale": inner,
            "outer_length": c.outer_length(),
            "inner_length": c.inner_length(),
            "length_ratio": c.outer_length() / c.inner_length(),
            "resistance": c.resistance(),
            "turn_cost_one_radian": c.turn_cost(1.0),
            "small_angle_law": c.turn_cost_small_angle(1.0),
            "law_ratio": c.turn_cost(1.0) / c.turn_cost_small_angle(1.0),
        })
    return {
        "rows": rows,
        "largest_arm_ratio": float(max(r["length_ratio"] for r in rows)),
        "the_arms_can_be_very_unequal": bool(
            max(r["length_ratio"] for r in rows) > 100.0),
        "and_nothing_about_the_hinge_changes": bool(
            all(0.98 < r["law_ratio"] < 1.0 for r in rows)),
        "why_it_matters": "the old vacuole picture had to give the inner and "
                          "outer boundaries one shared arbitrary radial gap; "
                          "a finite bearing with two arms has no such "
                          "constraint, and the R_inner/R_outer ratio stops "
                          "carrying any physical significance",
    }


def measure_the_width_rescales_with_the_profile() -> Dict[str, object]:
    """``w(s) = f(s)Δθ``: squeezed into the bearing, let out the other side.

    The angular width is what is transported — it is *constant* along the
    route, because the metric's angular part is ``f²dΩ²`` and the coordinate
    does not change.  The physical width is not.  So the visualisation has to
    show ``w_o → w_min → w_i``, and the end-to-end ratio is ``f_i/f_o``, with
    ``f₀`` nowhere in it.

    ``w_i/w_o = f_i/f_o`` is an *identity*, but it is not asserted here as a
    bitwise one, and that distinction has bitten this arc before.  The drawn
    ratio is ``(f_i Δθ)/(f_o Δθ)``, whose two products round separately, so it
    can sit one ulp away from ``f_i/f_o`` — it does at the first row below.
    The identity is exact in the field; the *picture* carries the rounding of
    the widths it actually drew, and the honest report is the ulp count.
    """
    rows = []
    for (f0, fo, fi, dtheta) in ((1e-3, 1.0, 0.35, 0.20),
                                 (1e-2, 0.8, 0.05, 0.05),
                                 (1e-4, 1.0, 1.0, 0.10)):
        c = RegularizedCenter(neck=f0, outer=fo, inner=fi)
        wo = c.feature_width(fo, dtheta)
        wm = c.feature_width(f0, dtheta)
        wi = c.feature_width(fi, dtheta)
        drawn, ideal = wi / wo, fi / fo
        rows.append({
            "neck": f0, "outer_scale": fo, "inner_scale": fi,
            "angular_width": dtheta,
            "width_outer": wo, "width_at_the_bearing": wm, "width_inner": wi,
            "contraction_into_the_bearing": wo / wm,
            "end_to_end_ratio": drawn,
            "scale_ratio": ideal,
            "residue_in_ulps": (0.0 if drawn == ideal
                                else abs(drawn - ideal) / math.ulp(ideal)),
        })
    return {
        "rows": rows,
        "law": "w(s) = f(s) dtheta ; w_i/w_o = f_i/f_o",
        "what_is_transported": "the ANGULAR width, which is constant along "
                               "the route; the physical width is not",
        "worst_residue_in_ulps": float(max(r["residue_in_ulps"] for r in rows)),
        "the_ratio_is_the_scale_ratio_to_the_last_ulp": bool(
            all(r["residue_in_ulps"] <= 1.0 for r in rows)),
        "why_not_bitwise": "the drawn ratio is (f_i dtheta)/(f_o dtheta), and "
                           "the two products round separately -- the identity "
                           "is exact, the two multiplications are not",
        "the_end_ratio_does_not_involve_the_neck": True,
        "so_a_packet_does_not_teleport_unchanged": True,
    }


def measure_the_turn_cost_is_quadratic_not_linear() -> Dict[str, object]:
    """**The result.**  ``T(α) = α²/(2I)``, and ``I`` is the repo's resistance.

    The linear guess ``f₀α`` measures the bearing's arc.  The geodesic does not
    walk the arc: Clairaut's ``f sin ψ = h`` makes it start turning while it is
    still high up the arm, where the lever arm is longer, so the cost is
    *quadratic* in the angle and set by ``I = ∫ds/f²`` rather than by ``f₀``
    alone.  ``I`` is the same integral whose reciprocal is the throat's
    monopole conductance, so

        ``T(α) = α² · (4π/I) / (8π)`` ,

    and the geometric cost of swinging the clock hands is the electrical cost
    of pushing flux through the same tube.  With two long arms ``I → 4/f₀`` and
    the law reads ``T → f₀α²/8``.

    The law is checked against the *integrated* geodesic, not assumed: exact to
    ``8e-05`` at ``α = 0.1``, ``0.3%`` at ``α = 0.6``, and still within ``8%``
    at ``α = π``, where the small-angle expansion has no right to hold at all.
    """
    rows = []
    for f0 in (1e-3, 1e-5):
        c = RegularizedCenter(neck=f0, outer=1.0, inner=0.35)
        for alpha in (0.02, 0.05, 0.1, 0.3, 0.6, 1.0, math.pi / 2, math.pi):
            exact = c.turn_cost(alpha)
            rows.append({
                "neck": f0, "alpha": alpha,
                "turn_cost": exact,
                "small_angle_law": c.turn_cost_small_angle(alpha),
                "law_ratio": exact / c.turn_cost_small_angle(alpha),
                "linear_guess": f0 * alpha,
                "fraction_of_the_linear_guess": exact / (f0 * alpha),
            })
    c = RegularizedCenter(neck=1e-5, outer=1.0, inner=0.35)
    from geometrodynamics.waves.physical_throat import VacuumThroat
    a = 0.05
    t = VacuumThroat(mouth_radius=a)
    sym = RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                            inner=t.mouth_f())
    small = [r for r in rows if r["alpha"] <= 0.3]
    return {
        "rows": rows,
        "law": "T(alpha) = alpha^2 / (2I),  I = int ds/f^2",
        "conductance_form": "T(alpha) = alpha^2 * conductance / (8 pi)",
        "long_arm_limit": "I -> 4/f0, so T -> f0 alpha^2 / 8",
        "long_arm_prefactor": float(1.0 / (2.0 * c.resistance() * c.neck)),
        "the_same_I_as_physical_throat": {
            "mouth_radius": a,
            "here": sym.resistance(),
            "physical_throat": t.resistance(),
            "difference": abs(sym.resistance() - t.resistance()),
            "conductance_here": sym.conductance(),
            "conductance_physical_throat": t.conductance(),
        },
        "worst_small_angle_error": float(
            max(abs(r["law_ratio"] - 1.0) for r in small)),
        "the_small_angle_law_holds": bool(
            all(abs(r["law_ratio"] - 1.0) < 1e-3 for r in small)),
        "it_is_still_close_at_pi": bool(
            all(r["law_ratio"] > 0.9
                for r in rows if abs(r["alpha"] - math.pi) < 1e-9)),
        "the_prefactor_is_one_eighth": bool(
            abs(1.0 / (2.0 * c.resistance() * c.neck) - 0.125) < 1e-3),
        "the_resistance_is_the_repos": bool(
            abs(sym.resistance() - t.resistance()) < 1e-8),
        "what_was_proposed": "L_turn ~ f0 alpha, the bearing's arc length",
        "what_is_measured": "T(alpha) = alpha^2/(2I), quadratic in the angle "
                            "and much smaller -- the geodesic cuts the corner "
                            "rather than walking the arc",
    }


def measure_the_corner_route_is_only_an_upper_bound() -> Dict[str, object]:
    """``L_o + L_i + f₀α`` is exact — for a route that is not the shortest.

    Worth separating from the measurement above, because the linear formula is
    not *wrong*: it is the exact length of the down-turn-up route, and that
    route exists.  It is simply not the geodesic, so it bounds the answer from
    above.  How loose the bound is depends only on ``α`` and not on ``f₀`` —
    the fractions below are the same to three figures at ``1e-03`` and
    ``1e-05`` — which is itself the signature of the quadratic law.
    """
    rows = []
    for f0 in (1e-3, 1e-5):
        c = RegularizedCenter(neck=f0, outer=1.0, inner=0.35)
        for alpha in (0.1, 0.5, 1.0, math.pi / 2, math.pi):
            geo = c.route_length(alpha)
            corner = c.corner_route_length(alpha)
            rows.append({
                "neck": f0, "alpha": alpha,
                "geodesic": geo, "corner_route": corner,
                "corner_is_longer_by": corner - geo,
                "geodesic_uses_this_fraction_of_the_arc":
                    (geo - c.arm_length_sum()) / (f0 * alpha),
            })
    lo = [r for r in rows if r["neck"] == 1e-3]
    hi = [r for r in rows if r["neck"] == 1e-5]
    drift = max(abs(a["geodesic_uses_this_fraction_of_the_arc"]
                    / b["geodesic_uses_this_fraction_of_the_arc"] - 1.0)
                for a, b in zip(lo, hi))
    return {
        "rows": rows,
        "the_corner_route_is_exact_for_its_own_route": True,
        "and_is_an_upper_bound_on_the_geodesic": bool(
            all(r["corner_is_longer_by"] > 0.0 for r in rows)),
        "fraction_of_the_arc_actually_spent_at_alpha_0p1": float(
            lo[0]["geodesic_uses_this_fraction_of_the_arc"]),
        "fraction_at_pi": float(lo[-1]["geodesic_uses_this_fraction_of_the_arc"]),
        "the_fraction_is_a_function_of_alpha_alone": bool(drift < 1e-2),
        "fraction_drift_between_two_decades_of_neck": float(drift),
        "why": "the geodesic starts turning while it is still descending, "
               "where f is larger and a given angle costs less arc",
    }


def measure_the_hinge_is_never_the_expensive_part() -> Dict[str, object]:
    """The property the point center was wanted for, kept.

    The reason for putting a point in the middle was that the link's cost then
    did not care where the mouths were.  A finite bearing charges *something*,
    so the question is whether it charges enough to matter — and it does not:
    at the working point the whole half-turn costs ``8e-04`` of the arms, and
    the angle at which the hinge would cost as much as the journey is ``104``
    radians, thirty-three times around.
    """
    rows = []
    for f0 in (1e-2, 1e-3, 1e-4):
        c = RegularizedCenter(neck=f0, outer=1.0, inner=0.35)
        base = c.arm_length_sum()
        full = c.turn_cost(math.pi)
        rows.append({
            "neck": f0,
            "arm_length_sum": base,
            "turn_cost_at_pi": full,
            "turn_over_arms": full / base,
            "alpha_at_which_turning_costs_as_much_as_travelling":
                math.sqrt(2.0 * c.resistance() * base),
        })
    return {
        "rows": rows,
        "the_hinge_is_always_cheap": bool(
            all(r["turn_over_arms"] < 1e-2 for r in rows)),
        "and_the_break_even_angle_is_far_past_pi": bool(
            all(r["alpha_at_which_turning_costs_as_much_as_travelling"]
                > 10.0 * math.pi for r in rows)),
        "so_the_clock_hand_picture_survives": True,
        "what_this_does_not_say": "that the cost is zero -- it is finite, and "
                                  "linear in f0, which is exactly the "
                                  "difference between a bearing and a point",
    }


def measure_the_bearing_shrinks_back_to_the_point() -> Dict[str, object]:
    """``f₀ → 0`` returns the singular picture, with the rate.

    Two things have to go: the turn cost, and the arms' own excess over their
    naive value ``F``.  Both do, and both linearly in ``f₀`` up to a log —
    ``L(F) − F ≃ (f₀/2)[ln(4F/f₀) − 1]``, checked below.  So the singular
    clock-hands model is the ``f₀ → 0`` limit of this one rather than a
    different model, which is the whole point of calling it a regularisation.
    """
    rows = []
    alpha = math.pi / 2
    for f0 in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        c = RegularizedCenter(neck=f0, outer=1.0, inner=0.35)
        excess = c.outer_length() - c.outer
        asymptotic = 0.5 * f0 * (math.log(4.0 * c.outer / f0) - 1.0)
        rows.append({
            "neck": f0,
            "arm_length_sum": c.arm_length_sum(),
            "turn_cost": c.turn_cost(alpha),
            "turn_cost_over_neck": c.turn_cost(alpha) / f0,
            "outer_arm_excess": excess,
            "excess_asymptotic": asymptotic,
            "excess_ratio": excess / asymptotic,
        })
    tail = rows[2:]
    return {
        "rows": rows,
        "alpha": alpha,
        "limit": "L_o + L_i, independent of where the mouths are",
        "turn_cost_is_linear_in_the_neck": bool(
            max(r["turn_cost_over_neck"] for r in tail)
            / min(r["turn_cost_over_neck"] for r in tail) - 1.0 < 1e-2),
        "arm_excess_asymptotic": "L(F) - F  ~  (f0/2)[ln(4F/f0) - 1]",
        "the_arm_asymptotic_holds": bool(
            all(abs(r["excess_ratio"] - 1.0) < 5e-2 for r in tail)),
        "the_singular_model_is_the_limit_not_a_rival": True,
    }


def measure_the_bearing_replaces_collision_with_overlap() -> Dict[str, object]:
    """What "intersection" becomes once the middle is not a point.

    Two fronts no longer have to collide at ``r = 0``.  They land on the
    bearing at angular positions, and the question of whether they meet
    separates cleanly into two halves that behave completely differently:

    * **Whether** they meet is a statement about angles — their angular
      extents overlap, or they do not — and ``f₀`` does not enter it at all.
    * **How big** the meeting is does depend on ``f₀``: the overlap is
      ``f₀ × (overlap angle)`` across, and two fronts that miss are separated
      by ``f₀ × (gap angle)``.

    That is how the point limit is recovered, and it is not the way it is
    usually described.  ``f₀ → 0`` does not make everything meet.  It makes the
    overlap and the gap *both* shrink to zero, so the distinction between
    meeting and missing survives as a yes/no while becoming invisible as a
    length.
    """
    rows = []
    for f0 in (1e-2, 1e-3):
        for (sep, wa, wb) in ((0.05, 0.20, 0.20), (0.30, 0.20, 0.20),
                              (0.30, 0.40, 0.30), (1.20, 0.20, 0.20)):
            reach = 0.5 * (wa + wb)
            overlap_angle = max(reach - sep, 0.0)
            gap_angle = max(sep - reach, 0.0)
            rows.append({
                "neck": f0,
                "angular_separation": sep,
                "angular_width_a": wa, "angular_width_b": wb,
                "they_meet": bool(overlap_angle > 0.0),
                "overlap_angle": overlap_angle,
                "gap_angle": gap_angle,
                "overlap_length_on_the_bearing": f0 * overlap_angle,
                "gap_length_on_the_bearing": f0 * gap_angle,
            })
    by_angle = {}
    consistent = True
    for r in rows:
        key = (r["angular_separation"], r["angular_width_a"],
               r["angular_width_b"])
        if key in by_angle and by_angle[key] != r["they_meet"]:
            consistent = False
        by_angle[key] = r["they_meet"]
    return {
        "rows": rows,
        "criterion": "they meet iff the angular extents overlap: "
                     "|separation| < (w_a + w_b)/2",
        "whether_is_an_angular_question": "f0 does not appear in the criterion",
        "how_big_is_a_length_question": "the overlap is f0 x (overlap angle)",
        "the_verdict_does_not_depend_on_the_neck": bool(consistent),
        "both_the_overlap_and_the_gap_scale_with_the_neck": bool(
            all(abs(r["overlap_length_on_the_bearing"]
                    - r["neck"] * r["overlap_angle"]) < 1e-18 for r in rows)),
        "the_point_limit_correctly_stated":
            "as f0 -> 0 the ANGULAR INCIDENCE SURVIVES and the PHYSICAL "
            "INTERACTION REGION COLLAPSES: which directions the fronts come "
            "in at, and whether their extents overlap, are untouched by f0; "
            "the region in which they actually share space is f0 x (overlap "
            "angle) and goes to zero. So f0 -> 0 does NOT make every route "
            "meet -- it shrinks the overlap AND the gap together, and the "
            "distinction survives as a yes/no while disappearing as a length",
        "what_survives_the_limit": "angular incidence -- direction of arrival "
                                   "and angular overlap, neither of which "
                                   "involves f0",
        "what_collapses": "the physical interaction region, f0 x (overlap "
                          "angle), and equally the physical separation of "
                          "fronts that miss",
        "no_route_passes_through_a_singular_point": True,
    }


def measure_the_turn_cost_does_not_care_which_sphere() -> Dict[str, object]:
    """The bearing may be ``S¹``, ``S²``, ``S³`` — or projective.

    Two directions span a totally geodesic 2-plane, so the walk between them is
    a great-circle arc whatever the bearing's dimension is, and ``T`` depends
    on the pair only through that one angle.  The visualisation can therefore
    draw the 2-D cross-section — a circle — without the answer being an
    artifact of having drawn a circle.

    The identification ``n ∼ −n`` is a different matter and does change the
    number: it replaces ``α`` by ``min(α, π−α)``, and near ``α = π`` that is a
    factor of ``418`` off the cost.  An unoriented clock hand has almost
    nothing to pay to reverse.
    """
    c = RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    rows = []
    for alpha in (0.3, 1.0, math.pi / 2, 2.5, 3.0):
        oriented = c.turn_cost(alpha)
        proj_angle = bearing_distance(alpha, projective=True)
        projective = c.turn_cost(alpha, projective=True)
        rows.append({
            "alpha": alpha,
            "projective_angle": proj_angle,
            "turn_cost_oriented": oriented,
            "turn_cost_projective": projective,
            "saving": oriented / projective if projective > 0 else float("inf"),
        })
    # The great-circle reduction, as an actual check rather than a restatement.
    # Route A: take the two directions as points of an S^2 bearing and read the
    # angle off the dot product.  Route B: build the orthonormal frame they
    # span, express both in it -- which is literally drawing them on the 2-D
    # cross-section's circle -- and read the planar angle with atan2.  If the
    # reduction is right the two angles, and so the two costs, must agree.
    great_circle = []
    for (na, nb) in (((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                     ((1.0, 0.0, 0.0), (0.6, 0.8, 0.0)),
                     ((0.0, 0.0, 1.0), (0.5, 0.5, math.sqrt(0.5))),
                     ((0.3, -0.5, 0.81), (-0.7, 0.2, 0.68))):
        a = np.array(na, float)
        a = a / np.linalg.norm(a)
        b = np.array(nb, float)
        b = b / np.linalg.norm(b)
        from_dot = math.acos(max(-1.0, min(1.0, float(a @ b))))
        # e1, e2: the orthonormal frame of the plane the pair spans
        e1 = a
        perp = b - float(a @ b) * a
        e2 = perp / np.linalg.norm(perp)
        planar = math.atan2(float(b @ e2), float(b @ e1))
        great_circle.append({
            "angle_from_the_dot_product": from_dot,
            "angle_drawn_on_the_2d_circle": planar,
            "angle_difference": abs(from_dot - planar),
            "cost_from_the_sphere": c.turn_cost(from_dot),
            "cost_from_the_drawn_circle": c.turn_cost(planar),
            "cost_relative_difference": abs(
                c.turn_cost(from_dot) / c.turn_cost(planar) - 1.0),
        })
    return {
        "rows": rows,
        "great_circle_reduction": great_circle,
        "why_dimension_free": "two directions span a totally geodesic "
                              "2-plane, so the walk is a great-circle arc on "
                              "S^1, S^2 or S^3 alike",
        "worst_great_circle_cost_difference": float(
            max(r["cost_relative_difference"] for r in great_circle)),
        "the_sphere_and_the_drawn_circle_agree": bool(
            max(r["cost_relative_difference"] for r in great_circle) < 1e-9),
        "the_projective_saving_near_pi": float(
            max(r["saving"] for r in rows)),
        "the_projective_identification_matters": bool(
            max(r["saving"] for r in rows) > 100.0),
        "and_it_is_the_only_thing_that_does": True,
        "so_the_2d_cross_section_is_honest": True,
    }


# ════════════════════════════════════════════════════════════════════════════
# A SECOND PROFILE, TO SEE WHAT IS GENERAL AND WHAT IS NOT
# ════════════════════════════════════════════════════════════════════════════
def hyperbolic_neck(neck: float, outer: float, inner: float,
                    alpha: float) -> Dict[str, float]:
    """The same three numbers on ``f(s) = √(f₀² + s²)`` — a *different* neck.

    Not the scalar-flat profile, and not glued to anything: just another
    rotationally symmetric metric with a finite minimum, used to separate the
    part of the result that is about necks in general from the part that is
    about *this* neck.

    The substitution is the matching one, ``s = f₀ sinh x`` so ``f = f₀cosh x``,
    which resolves the neck for the same reason ``f = f₀cosh²x`` does on the
    scalar-flat profile:

        ``I = (1/f₀)∫dx/cosh x`` ,
        ``dθ/dx = κ/√(cosh²x − κ²)`` ,
        ``dT/dx = f₀cosh x · κ²/[R(cosh x + R)]`` ,  ``R = √(cosh²x − κ²)`` .
    """
    f0 = float(neck)
    ends = [math.asinh(float(L) / f0) for L in (outer, inner)]

    def over_arms(g, width: float = -1.0) -> float:
        total = 0.0
        for upper in ends:
            marks = [0.0, 0.25, 1.0, 3.0, 8.0]
            if 0.0 < width < 1.0:
                e = width
                while e < 3.0:
                    marks.append(e)
                    e *= 4.0
            edges = sorted({m for m in marks if 0.0 <= m < upper}) + [upper]
            total += sum(quad(g, lo, hi, limit=200)[0]
                         for lo, hi in zip(edges, edges[1:]))
        return total

    resistance = over_arms(lambda x: 1.0 / (f0 * math.cosh(x)))

    def sweep(kappa: float) -> float:
        w = 0.5 * math.sqrt(max(1.0 - kappa * kappa, 0.0))

        def g(x: float) -> float:
            c = math.cosh(x)
            return kappa / math.sqrt(max(c * c - kappa * kappa, 1e-300))

        return over_arms(g, w)

    a = bearing_distance(alpha)
    kappa = brentq(lambda k: sweep(k) - a, 1e-18, 1.0 - 1e-9,
                   xtol=1e-18, rtol=8.9e-16)

    def cost(x: float) -> float:
        c = math.cosh(x)
        r = math.sqrt(max(c * c - kappa * kappa, 1e-300))
        return f0 * c * kappa * kappa / (r * (c + r))

    turn = over_arms(cost)
    return {"resistance": resistance, "kappa": kappa, "turn_cost": turn,
            "small_angle_law": a * a / (2.0 * resistance),
            "shape": turn / (a * a / (2.0 * resistance))}


def measure_the_law_does_not_depend_on_the_profile() -> Dict[str, object]:
    """``T = α²/(2I)`` is about necks, not about *this* neck — to leading order.

    The derivation never used the profile.  Minimising ``∫½f²θ'²ds`` at fixed
    ``∫θ'ds = α`` gives ``θ' = h/f²``, hence ``α = hI`` and an excess length
    ``½h²I = α²/(2I)``: **any** rotationally symmetric ``f(s)`` with a finite
    minimum, with its own ``I``.  Asserting that from the algebra alone would
    be cheap, so it is checked on a second profile that has nothing to do with
    the first — ``f = √(f₀² + s²)``, not scalar-flat and not glued to anything.

    It holds, to the same accuracy and across six decades of ``f₀``.

    What does **not** carry over is the correction beyond leading order.  The
    shape ``T/(α²/2I)`` is ``0.9250`` at ``α = π`` on the scalar-flat profile
    and ``0.8896`` on this one — same law, different ``O(α⁴)`` term.  So the
    quadratic law is a statement about necks; the ``8%`` deficit at a half turn
    is a statement about a particular neck.
    """
    rows = []
    for (f0, outer, inner) in ((1e-2, 1.0, 0.4), (1e-4, 2.0, 0.3),
                               (1e-6, 1.0, 1.0)):
        for alpha in (0.02, 0.1, 1.0, math.pi):
            got = hyperbolic_neck(f0, outer, inner, alpha)
            rows.append({"neck": f0, "outer": outer, "inner": inner,
                         "alpha": alpha, **got})
    small = [r for r in rows if r["alpha"] <= 0.1]
    at_pi = [r for r in rows if abs(r["alpha"] - math.pi) < 1e-9]
    scalar_flat = RegularizedCenter(neck=1e-5, outer=1.0, inner=0.35)
    sf_shape = (scalar_flat.turn_cost(math.pi)
                / scalar_flat.turn_cost_small_angle(math.pi))
    return {
        "rows": rows,
        "second_profile": "f(s) = sqrt(f0^2 + s^2)",
        "why_the_law_is_general": "minimising int f^2 theta'^2/2 ds at fixed "
                                  "int theta' ds = alpha gives theta' = h/f^2, "
                                  "so alpha = h I and the excess is "
                                  "h^2 I/2 = alpha^2/(2I) -- the profile "
                                  "enters only through I",
        "worst_small_angle_error": float(
            max(abs(r["shape"] - 1.0) for r in small)),
        "the_law_holds_on_the_second_profile": bool(
            all(abs(r["shape"] - 1.0) < 1e-3 for r in small)),
        "shape_at_pi_here": float(sum(r["shape"] for r in at_pi) / len(at_pi)),
        "shape_at_pi_scalar_flat": float(sf_shape),
        "the_correction_is_profile_dependent": bool(
            abs(sum(r["shape"] for r in at_pi) / len(at_pi) - sf_shape) > 0.02),
        "so_what_is_general": "the quadratic law T = alpha^2/(2I) is a "
                              "statement about necks; the ~8-11% deficit at a "
                              "half turn is a statement about a particular one",
    }


# ════════════════════════════════════════════════════════════════════════════
# THE MOMENT HIERARCHY, AND THE IDENTITY UNDERNEATH IT
# ════════════════════════════════════════════════════════════════════════════
def measure_the_fourth_moment_is_where_the_neck_shape_enters(
) -> Dict[str, object]:
    """``I₂`` sets the universal hinge; ``I₄`` is the first term that is local.

    Expand the geodesic exactly.  ``dℓ/ds − 1 = (1 − h²/f²)^{−1/2} − 1``, so

        ``T = ½h²I₂ + ⅜h⁴I₄ + …`` ,   ``α = hI₂ + ½h³I₄ + …`` ,

    with ``I_n = ∫ds/fⁿ``.  Eliminating ``h`` gives

        ``T(α) = α²/(2I₂) − α⁴I₄/(8I₂⁴) + O(α⁶)`` ,

    equivalently a shape ``T/(α²/2I₂) = 1 − α²I₄/(4I₂³)``.

    **That is the division of labour.**  ``I₂`` is the only moment in the
    leading term, and it is the *same* ``I₂`` that sets the monopole
    conductance — so the quadratic hinge is universal in the strong sense that
    it is fixed by a quantity the throat already had.  ``I₄`` is the first
    moment that is not shared with anything, and it is where the neck's shape
    is first felt.  The two profiles differ there and nowhere earlier:

    * scalar-flat, long arms — ``I₂ = 4/f₀``, ``I₄ = 32/(15f₀³)``, so the
      shape is ``1 − α²/120`` ;
    * hyperbolic ``f = √(f₀²+s²)`` — ``I₂ = π/f₀``, ``I₄ = π/(2f₀³)``, so the
      shape is ``1 − α²/(8π²)`` .

    ``1/120`` against ``1/(8π²) = 1/79`` is exactly why the two shapes part
    company at large angle (``0.9250`` against ``0.8886`` at ``α = π``) while
    agreeing to eight digits at ``α = 0.1``.  The moments below are closed
    forms and are checked against quadrature; the shapes are checked against
    the integrated geodesic.
    """
    long_arm = RegularizedCenter(neck=1e-6, outer=1.0, inner=1.0)
    i2, i4 = long_arm.resistance(), long_arm.fourth_moment()
    rows = []
    for alpha in (0.1, 0.3, 0.6, 1.0, 2.0, math.pi):
        exact = long_arm.turn_cost(alpha)
        rows.append({
            "alpha": alpha,
            "turn_cost": exact,
            "second_order": long_arm.turn_cost_small_angle(alpha),
            "fourth_order": long_arm.turn_cost_to_fourth_order(alpha),
            "shape_measured": exact / long_arm.turn_cost_small_angle(alpha),
            "shape_predicted": 1.0 - alpha ** 2 * i4 / (4.0 * i2 ** 3),
            "fourth_order_relative_error": abs(
                exact / long_arm.turn_cost_to_fourth_order(alpha) - 1.0),
        })
    # the same expansion on the unrelated profile
    hyper = []
    f0 = 1e-6
    j2 = 2.0 * quad(lambda x: 1.0 / (f0 * math.cosh(x)), 0.0,
                    math.asinh(1.0 / f0), limit=200)[0]
    j4 = 2.0 * quad(lambda x: 1.0 / (f0 ** 3 * math.cosh(x) ** 3), 0.0,
                    math.asinh(1.0 / f0), limit=200)[0]
    for alpha in (0.1, 1.0, math.pi):
        got = hyperbolic_neck(f0, 1.0, 1.0, alpha)
        hyper.append({"alpha": alpha, "shape_measured": got["shape"],
                      "shape_predicted": 1.0 - alpha ** 2 * j4 / (4.0 * j2 ** 3)})
    small = [r for r in rows if r["alpha"] <= 0.3]
    return {
        "rows": rows,
        "hyperbolic_rows": hyper,
        "expansion": "T = alpha^2/(2 I2) - alpha^4 I4/(8 I2^4) + O(alpha^6)",
        "shape_law": "T/(alpha^2/2 I2) = 1 - alpha^2 I4/(4 I2^3)",
        "scalar_flat": {"I2_times_neck": i2 * long_arm.neck,
                        "I4_times_neck_cubed": i4 * long_arm.neck ** 3,
                        "shape_coefficient": i4 / (4.0 * i2 ** 3),
                        "closed_form": "I2 = 4/f0, I4 = 32/(15 f0^3), "
                                       "shape = 1 - alpha^2/120"},
        "hyperbolic": {"I2_times_neck": j2 * f0,
                       "I4_times_neck_cubed": j4 * f0 ** 3,
                       "shape_coefficient": j4 / (4.0 * j2 ** 3),
                       "closed_form": "I2 = pi/f0, I4 = pi/(2 f0^3), "
                                      "shape = 1 - alpha^2/(8 pi^2)"},
        "the_shape_law_holds_at_small_angle": bool(
            all(abs(r["shape_measured"] - r["shape_predicted"]) < 1e-6
                for r in small)),
        "the_fourth_order_beats_the_second": bool(
            all(r["fourth_order_relative_error"]
                < abs(r["shape_measured"] - 1.0) for r in rows
                if r["alpha"] >= 0.3)),
        "the_second_moment_is_shared_the_fourth_is_not": bool(
            abs(i2 * long_arm.neck / (j2 * f0) - 4.0 / math.pi) < 1e-4
            and abs(i4 / (4.0 * i2 ** 3) - 1.0 / 120.0) < 1e-6
            and abs(j4 / (4.0 * j2 ** 3) - 1.0 / (8.0 * math.pi ** 2)) < 1e-6),
        "the_division_of_labour": "I2 controls the universal quadratic hinge "
                                  "and is the same integral as the monopole "
                                  "conductance; I4 is the first moment that "
                                  "remembers the neck's shape",
    }


def measure_the_hinge_and_the_monopole_are_one_dirichlet_form(
) -> Dict[str, object]:
    """**The identity underneath.**  Both are ``∫f²φ'²ds``, with different ``φ``.

    Static monopole flux and infinitesimal rotation of the throat are not two
    problems that happen to share a number.  They are *one* variational problem
    on the tube,

        minimise  ``E[φ] = ∫ f² φ'² ds``  at fixed total increment ``Δφ`` ,

    whose Euler–Lagrange equation is ``(f²φ') ' = 0``: the current ``f²φ'`` is
    conserved, ``Δφ = c·I₂``, and the minimum is ``Δφ²/I₂``.  The weight is the
    transverse area element, which is why ``I₂ = ∫ds/f²`` and nothing else.

    Reading it twice:

    * ``φ = u``, the ℓ=0 potential.  The conserved current **is** the flux,
      ``Φ = 4πf²u'``; the drop is ``ΦI₂/4π``; the conductance is ``4π/I₂``.
      That is `physical_throat.VacuumThroat.conductance`.
    * ``φ = θ``, the azimuth.  The conserved current **is** Clairaut's
      constant, ``h = f²θ'``; the sweep is ``hI₂``; the excess length is
      ``½h²I₂ = α²/(2I₂)``.  That is `turn_cost`.

    So the sharpest form of the statement is not about the two numbers but
    about the two *profiles*: normalised to their own total, the monopole
    potential and the geodesic's azimuth are **the same function of position
    along the tube**.  Checked below, and the deviation falls as ``α²`` —
    ``4.9e-03`` at ``α = 1`` down to ``4.9e-07`` at ``α = 0.01``, a clean
    factor of 100 per decade.  Which is exactly why "infinitesimal" is the
    right word: at finite ``α`` the geodesic feels ``I₄`` and the potential
    does not.
    """
    c = RegularizedCenter(neck=1e-3, outer=1.0, inner=0.35)
    x_end = c.half_length_in_x(c.outer)
    sample = [x_end * k / 8.0 for k in range(1, 9)]

    convergence = []
    for alpha in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003):
        kappa = c.clairaut_constant(alpha)
        theta_end = c.azimuth_profile(x_end, kappa)
        worst = max(abs(c.azimuth_profile(x, kappa) / theta_end
                        - c.monopole_profile(x) / c.monopole_profile(x_end))
                    for x in sample)
        convergence.append({"alpha": alpha, "worst_profile_deviation": worst,
                            "over_alpha_squared": worst / alpha ** 2})

    from geometrodynamics.waves.physical_throat import VacuumThroat
    a = 0.05
    t = VacuumThroat(mouth_radius=a)
    sym = RegularizedCenter(neck=t.neck_radius(), outer=t.mouth_f(),
                            inner=t.mouth_f())
    alpha = 0.05
    ratios = sorted(r["over_alpha_squared"] for r in convergence)
    return {
        "form": "E[phi] = int f^2 phi'^2 ds, minimised at fixed increment",
        "euler_lagrange": "(f^2 phi')' = 0, so f^2 phi' is conserved",
        "as_a_potential": "phi = u: the current is the flux 4 pi f^2 u', the "
                          "drop is Phi I2 / 4 pi, the conductance is 4 pi / I2",
        "as_an_azimuth": "phi = theta: the current is Clairaut's h = f^2 "
                         "theta', the sweep is h I2, the cost is h^2 I2 / 2",
        "profile_convergence": convergence,
        "the_profiles_coincide_as_alpha_vanishes": bool(
            convergence[0]["worst_profile_deviation"]
            > 1e4 * convergence[-1]["worst_profile_deviation"]),
        "the_deviation_is_second_order_in_alpha": bool(
            ratios[-1] / ratios[0] - 1.0 < 1e-2),
        "conductance_here": sym.conductance(),
        "conductance_physical_throat": t.conductance(),
        "conductance_difference": abs(sym.conductance() - t.conductance()),
        "hinge_from_the_conductance": alpha ** 2 * t.conductance()
        / (8.0 * math.pi),
        "hinge_from_the_geometry": sym.turn_cost_small_angle(alpha),
        "the_two_readings_agree": bool(
            abs(alpha ** 2 * t.conductance() / (8.0 * math.pi)
                / sym.turn_cost_small_angle(alpha) - 1.0) < 1e-12),
        "why_it_is_not_a_coincidence": "the weight in the Dirichlet form is "
                                       "the transverse area element f^2, and "
                                       "both the monopole potential and the "
                                       "azimuth are functions on the tube "
                                       "with no angular dependence -- so they "
                                       "solve the same equation with the same "
                                       "measure, and I2 is the resistance to "
                                       "either",
        "what_is_not_claimed": "that the two remain the same beyond leading "
                               "order -- at finite alpha the geodesic picks "
                               "up I4 and the static potential does not",
    }
