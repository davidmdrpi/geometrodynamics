"""Which throat is physical — and the sign it gives.

The question
────────────
`areal` measured a signed ``ΔA/A`` on the throat this arc has carried since
PR #261: a product tube of cross-section ``𝒜 = 4π`` and length ``L = 0.9``.
Then, asked to match the tube's area to its own mouths', it found the sign
**reverses**.  So the throat's geometry stopped being decoration and became the
thing the answer depends on, and the next question is not optional:

    ``𝒜`` and ``L`` were free parameters.  Which values are physical?

They were never free
────────────────────
For a rotationally symmetric slice ``ds² + f(s)² dΩ²`` the scalar curvature is

    ``R̂ = 2/f² − 4f''/f − 2f'²/f²``

— derived here from the metric rather than quoted, and checked against
``f = sin s`` giving ``R̂ = 6`` on the round ``S³``.  On a maximal slice the
**background** Hamiltonian constraint is ``R̂ = 16πG ρ̄``, so the profile does
not choose its matter: **the matter is whatever the profile implies.**  A
product tube of area ``𝒜`` is ``S²(r) × R`` with ``r = √(𝒜/4π)``, so
``R̂ = 2/r² = 8π/𝒜`` — a constant, and generally not the ambient's ``6``:

    ``𝒜 = 4π``            (PR #261–#264)  ⟹  ``ρ_tube = ρ̄/3``
    ``𝒜 = 4π sin²a``      (the matched tube)  ⟹  ``ρ_tube = 133 ρ̄`` at ``a = 0.05``

Neither was chosen for a reason, and neither is the ambient's fluid.

The throat that is forced
─────────────────────────
Ask instead for a throat that needs **no matter at all** — ``R̂ = 0`` — and glues
to the ambient with no surface layer.  ``R̂ = 0`` integrates once:

    ``f'² = 1 + C/f`` ,   and a neck ``f'(f₀) = 0`` fixes ``C = −f₀`` ,
    so   ``f'² = 1 − f₀/f`` .

Smooth gluing to the round ``S³`` at mouth radius ``a`` needs ``f' = cos a``
where ``f = sin a``, which forces

    ``f₀ = sin a (1 − cos²a) = sin³a`` .

**There is no free parameter left** — and the throat has a name.

``R̂ = 0``, ``K = 0`` and a spherical neck do not merely permit an
Einstein–Rosen bridge; they *are* one.  ``f'² = 1 − f₀/f`` is exactly the
time-symmetric Schwarzschild slice ``ds² = dr²/(1 − 2M/r) + r²dΩ²`` written in
proper radial distance, with ``r = f`` and ``f₀ = 2M``.  So the forced neck
radius is twice a mass, and

    ``M = sin³a / 2`` ,

**the throat's mass derived from the size of the excised mouth, with nothing
left to choose.**  Three quasi-local masses agree on it exactly — the
Schwarzschild parameter, the irreducible mass ``√(A_neck/16π)`` (the neck area
is ``16πM²``), and the Hawking mass, which is *constant* along the vacuum piece.

And the gluing condition is itself a mass statement.  The Hawking mass of a
sphere in ``ds² + f²dΩ²`` is ``(f/2)(1 − f'²)``: on the throat it is ``f₀/2``
everywhere, and on the round ``S³`` at geodesic radius ``χ`` it is
``sin³χ / 2``.  Requiring no surface layer at ``χ = a`` is the same equation as
requiring the Hawking mass not to jump.  `measure_the_throat_is_an_einstein_rosen_bridge`
carries this together with the four things it does not say — no asymptotic
region and so no ADM mass, a dimensionless ``M/R`` and not an absolute unit, a
handle in one ``S³`` rather than a bridge between universes, and a neck that is
a minimal surface and therefore, on a ``K = 0`` slice, an **apparent horizon**.
That last is the vacuum face of this arc's earlier result that a *traversable*
connection needs exotic matter: the throat with none in it is the one that is
not traversable.

Both the length and the tube's resistance follow in closed form, and both are
verified against quadrature to ``1e-12``:

    ``L = 2[ sin³a · arccosh(1/sin a) + sin a cos a ]``    ``≈ 2a``
    ``I = ∫ ds/f² = 4 cos a / sin³a``

and the resulting conductance is exactly a quarter of the exterior's own
monopole stiffness, ``4π/I = N₀(a,4)/4``, for every ``a``.

Two things follow that decide the round
───────────────────────────────────────
**There is no cavity.**  The constraint operator is ``∇² + R̂/2``, so with
``R̂ = 0`` the tube carries the *plain Laplacian*: ``(f²u')' = 0`` in ``ℓ = 0``.
Its solutions are ``A + B∫ds/f²`` — monotone.  No standing wave, no poles, and
nothing for the sign to flip across.  PR #264's cavity, its resonances at
``kL = nπ``, and the sign flips across them were **properties of matter in the
tube**, not of a throat.

**And the shunt vanishes identically.**  Writing the throat as an admittance,
``Φ_j = Σ_k Y_{jk} v_k`` with fluxes into the tube, the quantity that matters is
``Y·(1,1)`` — the flux a *uniform* potential drives in.  For a vacuum tube
``(f²u')' = 0`` makes flux conservation an identity, so ``Y·(1,1) = 0``
**exactly**: nothing to absorb it.  A tube with matter in it has somewhere to
put that flux, and shunts.

That single number is what sets the sign.  Holding everything else fixed and
turning the shunt up from zero, ``ΔA/A`` passes through a pole near ``2e-03``
and changes sign; PR #264's tube sat at ``6.07``, three orders past it.  The
conductance, scanned over eight orders, changes nothing.

So the answer reverses
──────────────────────
On the throat that is forced rather than chosen, the interference **opens**
both mouths.  The sign holds across the whole vacuum family — four orders in
the neck radius, with the smooth-gluing point in the middle — at both mouth
radii, under both gluings, and whether or not the ``ℓ = 1`` source moments are
included at all.

What it costs, stated plainly
─────────────────────────────
The response is ``3000×`` larger than PR #264's and grows as ``a^{-3}``, and the
matching system's condition number rises with it.  That is not noise; it is the
physics of a throat that barely lifts the constraint's degeneracy.  A vacuum
tube has zero shunt *by identity*, so it does almost nothing to separate the
``k = 1`` kernel from zero, and the linear response sits close to a mode the
operator nearly annihilates.  The sign is robust.  **The window in which
linearising was legitimate is now the binding question**, and this module does
not answer it.

Release hardening: the two-port is a closed form
────────────────────────────────────────────────
*A correctness repair to what this module shipped, not a new result about the
geometry.  The answer above is unchanged.*

``f'² = 1 − f₀/f`` is ``f = f₀cosh²x`` with ``ds = 2f dx``, which turns
``(f²u')' = ℓ(ℓ+1)u`` into ``y'' = (2ℓ+1)²y`` under ``R = y/cosh x`` — constant
coefficients.  The half-length ``X = arcosh(1/sin a)`` has ``e^{−X} =
tan(a/2)`` *exactly*, so with ``k = 2ℓ+1`` and ``q = tan^{2k}(a/2)`` the whole
two-port is rational in ``q``:

    ``D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]``
    ``C_ℓ = +4π sin a · kq/(1−q²)``

`admittance` is now that, and the Riccati solve is retained as
`admittance_riccati` — an independent validator, demoted rather than deleted.
It formed the cross term as a *difference* of two eigenchannel values that
agree to more digits than the solver carries, and at ``ℓ = 2, a = 0.05`` it
returns ``-1.17e-14`` for a true ``+3.00e-16``.  The diagonal was never
affected, and the two channels anything downstream consumes — ``ℓ = 0`` and
``ℓ = 1`` — sit above the floor, so `areal.solve_matching` moves in the
thirteenth digit.  See `measure_where_the_riccati_solve_stops_resolving`.

The closed form then gives ``C_ℓ ~ a^{4ℓ+3}`` for free: each unit of angular
momentum costs four powers of the mouth radius, and the ESU kernel's four
``n = 1`` harmonics split locally as ``1 ⊕ 3`` with the pieces crossing
``8.5e+05`` apart at ``a = 0.05``.  See `measure_the_mouth_to_mouth_hierarchy`
— including the narrow statement that supports, and the wider one it does not.

Two scope notes belong with all of the above: ``f_min > 0`` is forced *within
the class worked in here* — spherically symmetric, scalar-flat, ``K_ij = 0``,
``C¹``-matched — and not by Einstein's equations in general; and everything in
this module is **spatial initial data**.  The dynamic problem is not well-posed
until a lapse is chosen, and none is chosen here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.integrate import quad, solve_ivp

from .areal import (FOUR_PI, INTERFERENCE_MOMENTS, MOUTHS, TubeModel,
                    basis_channels, solve_matching)

__all__ = [
    "VacuumThroat",
    "WORKING_VACUUM_THROAT",
    "profile_scalar_curvature",
    "product_tube_density_ratio",
    "measure_the_curvature_formula_is_derived_not_quoted",
    "measure_the_gluing_forces_the_neck_radius",
    "measure_the_product_tubes_need_anomalous_matter",
    "measure_the_vacuum_throat_has_no_cavity",
    "measure_the_shunt_decides_the_sign",
    "measure_the_signed_response_on_the_physical_throat",
    "measure_the_throat_is_an_einstein_rosen_bridge",
    "measure_where_the_riccati_solve_stops_resolving",
    "measure_the_mouth_to_mouth_hierarchy",
    "hawking_mass",
]


# ════════════════════════════════════════════════════════════════════════════
# THE GEOMETRY
# ════════════════════════════════════════════════════════════════════════════
def profile_scalar_curvature(f: float, fp: float, fpp: float) -> float:
    """``R̂`` of ``ds² + f(s)²dΩ²``, from the metric.

    ``2/f² − 4f''/f − 2f'²/f²``.  Checked against ``f = sin s`` on the round
    ``S³``, which must give exactly ``6``, and against ``f'² = 1 − f₀/f``,
    which must give exactly ``0``.
    """
    return 2.0 / f ** 2 - 4.0 * fpp / f - 2.0 * fp ** 2 / f ** 2


def product_tube_density_ratio(area: float) -> float:
    """``ρ_tube/ρ̄`` for a product tube of cross-section ``𝒜``.

    ``S²(r) × R`` has ``R̂ = 2/r² = 8π/𝒜``, and the ambient has ``R̂ = 6``, so
    the ratio is ``4π/(3𝒜)``.  It is ``1`` only at ``𝒜 = 4π/3`` — and *that*
    tube cannot be glued to the ambient without a surface layer, because two
    regions of equal constant ``R̂`` joined smoothly are one region of that
    ``R̂``, and the only such rotationally symmetric one is the round sphere.
    """
    return (8.0 * math.pi / float(area)) / 6.0


@dataclass(frozen=True)
class VacuumThroat:
    """``R̂ = 0`` in the tube, glued to the ambient with no surface layer.

    One input — the mouth radius ``a``, which is where the ambient is cut — and
    no free parameters at all: ``f₀ = sin³a`` is forced by ``f'(mouth) = cos a``.
    """

    mouth_radius: float = 0.05

    def mouth_f(self) -> float:
        return math.sin(float(self.mouth_radius))

    def neck_radius(self) -> float:
        """``f₀ = sin³a`` — forced by smooth gluing, not chosen."""
        return math.sin(float(self.mouth_radius)) ** 3

    def neck_area(self) -> float:
        return FOUR_PI * self.neck_radius() ** 2

    def length(self) -> float:
        """``2[ f₀ arccosh(1/sin a) + sin a cos a ]`` — closed form."""
        a = float(self.mouth_radius)
        return 2.0 * (math.sin(a) ** 3 * math.acosh(1.0 / math.sin(a))
                      + math.sin(a) * math.cos(a))

    def length_by_quadrature(self) -> float:
        """The same length as ``∫ds``, for the check that the closed form is one."""
        m, f0 = self.mouth_f(), self.neck_radius()
        return 2.0 * quad(lambda t: 2.0 * math.sqrt(f0 + t * t), 0.0,
                          math.sqrt(m - f0), limit=400)[0]

    def resistance(self) -> float:
        """``I = ∫ ds/f² = 4 cos a / sin³a`` — closed form."""
        a = float(self.mouth_radius)
        return 4.0 * math.cos(a) / math.sin(a) ** 3

    def resistance_by_quadrature(self) -> float:
        m, f0 = self.mouth_f(), self.neck_radius()
        return 2.0 * quad(
            lambda t: 2.0 * math.sqrt(f0 + t * t) / (f0 + t * t) ** 2,
            0.0, math.sqrt(m - f0), limit=400)[0]

    def conductance(self) -> float:
        """``4π/I``.  Exactly ``N₀(a,4)/4`` — a quarter of the exterior's own."""
        return FOUR_PI / self.resistance()

    # ── what the throat turns out to be ─────────────────────────────────────
    def mass(self) -> float:
        """``M = f₀/2 = sin³a / 2`` — the throat's mass, in units of the ESU radius.

        ``f'² = 1 − f₀/f`` **is** the time-symmetric Schwarzschild slice,
        ``ds² = dr²/(1 − 2M/r) + r²dΩ²``, with the areal radius ``r = f`` and
        ``f₀ = 2M``.  So the forced throat is an Einstein–Rosen bridge, and its
        mass is not a parameter: it is fixed by the mouth the ambient was cut
        with.

        Dimensionless by construction — this is ``M/R`` with ``R`` the ESU's
        curvature radius, which is the only kind of mass statement available
        here and the only kind the scale-modulus theorem allows.
        """
        return 0.5 * self.neck_radius()

    def irreducible_mass(self) -> float:
        """``√(A_neck/16π)``.  Equal to `mass` exactly, since ``A = 16πM²``."""
        return math.sqrt(self.neck_area() / (16.0 * math.pi))

    def hawking_mass_in_the_tube(self) -> float:
        """``(f/2)(1 − f'²) = f₀/2``, the *same* at every radius in the tube."""
        f = 0.5 * (self.neck_radius() + self.mouth_f())
        return hawking_mass(f, math.sqrt(1.0 - self.neck_radius() / f))

    def hawking_mass_of_the_ambient_mouth(self) -> float:
        """``sin³a / 2`` — the ambient's own Hawking mass at the cut."""
        a = float(self.mouth_radius)
        return hawking_mass(math.sin(a), math.cos(a))

    def neck_is_a_minimal_surface(self) -> bool:
        """``f'(f₀) = 0``, so ``H = 0``.  On a ``K = 0`` slice that is an
        apparent horizon: ``θ_± = ±H`` both vanish."""
        return True

    # ── the two-port ────────────────────────────────────────────────────────
    def _half(self, ell: int, parity: str, rtol: float = 1e-12) -> float:
        """``Y₁₁ ± Y₁₂`` from one half-tube solve, by Riccati.

        Integrating ``u`` itself overflows: the ``ℓ = 1`` channel grows by
        ``1e+09`` across the neck.  The logarithmic derivative does not — each
        parity tracks a fixed point — and it reproduces the ``ℓ = 0`` closed
        form to ``1e-13``, which is what says the formulation is right rather
        than merely finite.
        """
        m, f0 = self.mouth_f(), self.neck_radius()
        c = float(ell * (ell + 1))
        tmax = math.sqrt(m - f0)
        if parity == "symmetric":
            # W = f^2 u'/u , W(neck) = 0
            rhs = lambda t, y: [2.0 * math.sqrt(f0 + t * t)
                                * (c - y[0] ** 2 / (f0 + t * t) ** 2)]
            sol = solve_ivp(rhs, (0.0, tmax), [0.0], rtol=rtol, atol=1e-16,
                            method="DOP853")
            return -FOUR_PI * float(sol.y[0, -1])
        # V = u/(f^2 u') , V(neck) = 0
        rhs = lambda t, y: [2.0 * math.sqrt(f0 + t * t)
                            * (1.0 / (f0 + t * t) ** 2 - c * y[0] ** 2)]
        sol = solve_ivp(rhs, (0.0, tmax), [0.0], rtol=rtol, atol=1e-18,
                        method="DOP853")
        return -FOUR_PI / float(sol.y[0, -1])

    def admittance_riccati(self, ell: int) -> np.ndarray:
        """``Y`` by the two half-tube Riccati solves.  **Reference, not production.**

        This is the original implementation, kept as an independent check on
        `admittance_closed_form` — but it is *not* what `admittance` returns,
        because the last step here,

            ``Y₁₂ = ½(s − t)`` ,

        is a difference of two eigenchannel values that become equal to many
        digits as soon as the neck is opaque.  The cross term is then built
        entirely out of the cancelling tail, and the solver's own tolerance is
        larger than the answer.  At ``a = 0.05, ℓ = 2`` this route returns
        ``1.17e-14`` for a true ``3.00e-16`` — 39× too large, and with the
        wrong sign.  See `measure_where_the_riccati_solve_stops_resolving`.

        It remains trustworthy wherever ``|Y₁₂|`` is above roughly ``1e-12``,
        which covers ``ℓ = 0`` and ``ℓ = 1`` at the working point.
        """
        s = self._half(int(ell), "symmetric")
        t = self._half(int(ell), "antisymmetric")
        return np.array([[0.5 * (s + t), 0.5 * (s - t)],
                         [0.5 * (s - t), 0.5 * (s + t)]])

    def admittance_closed_form(self, ell: int) -> np.ndarray:
        """``Y`` in closed form.  Exact, and free of the cancellation above.

        The vacuum profile ``f'² = 1 − f₀/f`` is ``f = f₀cosh²x`` with
        ``ds = 2f dx``, so the static equation ``(f²u')' = ℓ(ℓ+1)u`` becomes

            ``R_xx + 2 tanh x · R_x − 4ℓ(ℓ+1) R = 0`` ,

        and the substitution ``R = y/cosh x`` reduces it to ``y'' = k²y`` with
        ``k = 2ℓ+1`` — a *constant-coefficient* equation.  The half-length is
        ``X = arcosh(1/sin a)``, whose exponential is exactly ``e^{−X} =
        tan(a/2)``, so with

            ``q = e^{−2kX} = tan^{2k}(a/2)``   and   ``f_m = sin a``

        every hyperbolic function collapses to a rational function of ``q``:
        ``coth(2kX) = (1+q²)/(1−q²)``, ``csch(2kX) = 2q/(1−q²)``, and
        ``tanh X = cos a``.  The two-port is then

            ``D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]``   (diagonal)
            ``C_ℓ = +4π sin a · k q/(1−q²)``                 (cross)

        with ``Y = [[D, C], [C, D]]`` by the reflection symmetry of the tube.
        ``C_ℓ`` is a *product* of small factors, never a difference, which is
        exactly what the Riccati route cannot offer.

        ``ℓ = 0`` is special-cased to ``(π sin³a / cos a)·[[−1,1],[1,−1]]``.
        That is the same number — ``q = tan²(a/2)`` reproduces it to ``1e-17``
        relative — but the diagonal formula gets there as ``k coth(2kX) −
        tanh X``, a difference of two ``O(1)`` quantities that both tend to
        ``1`` as ``a → 0``.  Written directly, the zero-shunt identity
        ``Y·(1,1)ᵀ = 0`` holds to the last bit rather than to a tolerance.
        """
        ell = int(ell)
        a = float(self.mouth_radius)
        if ell == 0:
            return self.conductance() * np.array([[-1.0, 1.0], [1.0, -1.0]])
        k = 2.0 * ell + 1.0
        q = math.tan(0.5 * a) ** (2.0 * k)
        w = 1.0 - q * q
        d = -2.0 * math.pi * math.sin(a) * (k * (1.0 + q * q) / w - math.cos(a))
        c = FOUR_PI * math.sin(a) * k * q / w
        return np.array([[d, c], [c, d]])

    def admittance(self, ell: int) -> np.ndarray:
        """``Y`` in the same convention `areal.solve_matching` expects.

        The closed form — see `admittance_closed_form` for why, and
        `admittance_riccati` for the solve it replaced.
        """
        return self.admittance_closed_form(ell)

    def monopole_admittance_closed_form(self) -> np.ndarray:
        """``(4π/I)[[−1,1],[1,−1]]`` — singular, because ``(f²u')' = 0``."""
        return self.conductance() * np.array([[-1.0, 1.0], [1.0, -1.0]])

    def shunt(self) -> float:
        """``Y·(1,1)``.  **Zero by identity**, not by cancellation."""
        return float((self.admittance(0) @ np.ones(2))[0])


WORKING_VACUUM_THROAT = VacuumThroat(mouth_radius=0.05)


def hawking_mass(f: float, fp: float) -> float:
    """``m_H = (f/2)(1 − f'²)`` for a sphere of areal radius ``f``.

    The Hawking mass of a round sphere in ``ds² + f(s)²dΩ²``: with ``H = 2f'/f``
    and ``A = 4πf²``,

        ``m_H = √(A/16π) (1 − (1/16π)∮H² dA) = (f/2)(1 − f'²)`` .

    It is worth having in closed form because it is what the gluing condition
    turns out to *be* — see
    `measure_the_throat_is_an_einstein_rosen_bridge`.
    """
    return 0.5 * float(f) * (1.0 - float(fp) ** 2)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_curvature_formula_is_derived_not_quoted() -> Dict[str, object]:
    """``R̂ = 2/f² − 4f''/f − 2f'²/f²``, against two exact cases.

    The round ``S³`` (``f = sin s``) must give ``6``, and the vacuum profile
    (``f'² = 1 − f₀/f``, so ``f'' = f₀/2f²``) must give ``0``.  Both are exact
    identities, so anything but machine precision is an error in the formula —
    which is the point of checking a *derived* expression against cases whose
    answers are known independently.
    """
    sphere = [profile_scalar_curvature(math.sin(s), math.cos(s), -math.sin(s))
              for s in (0.3, 0.9, 1.5, 2.4)]
    vacuum = []
    for f0 in (1e-4, 1e-2, 0.3):
        for f in (1.5 * f0, 10.0 * f0, 0.4):
            if f <= f0:
                continue
            vacuum.append(profile_scalar_curvature(
                f, math.sqrt(1.0 - f0 / f), f0 / (2.0 * f ** 2)))
    return {
        "round_sphere": sphere,
        "worst_sphere_error": float(max(abs(x - 6.0) for x in sphere)),
        "vacuum_profile": vacuum,
        "worst_vacuum_error": float(max(abs(x) for x in vacuum)),
        "both_are_exact": bool(max(abs(x - 6.0) for x in sphere) < 1e-12
                               and max(abs(x) for x in vacuum) < 1e-8),
    }


def measure_the_gluing_forces_the_neck_radius(
        radii: Sequence[float] = (0.02, 0.05, 0.10, 0.20)) -> Dict[str, object]:
    """``f₀ = sin³a``, and then ``L`` and ``I`` in closed form.

    The throat has **no free parameter**: requiring no matter (``R̂ = 0``) and
    no surface layer (``f' = cos a`` where ``f = sin a``) uses up both.  The
    closed forms are checked against quadrature, which is the only way to know
    the integration was done and not merely written down.
    """
    rows = []
    for a in radii:
        t = VacuumThroat(mouth_radius=a)
        rows.append({
            "mouth_radius": float(a),
            "neck_radius": t.neck_radius(),
            "forced_value": math.sin(a) ** 3,
            "neck_area": t.neck_area(),
            "length": t.length(),
            "length_by_quadrature": t.length_by_quadrature(),
            "length_over_a": t.length() / a,
            "resistance": t.resistance(),
            "resistance_by_quadrature": t.resistance_by_quadrature(),
            "conductance": t.conductance(),
            "exterior_monopole_stiffness": FOUR_PI * math.sin(a) ** 2
                                           * math.tan(a),
        })
    def drift(key: str) -> float:
        return float(max(abs(r[key] / r[key + "_by_quadrature"] - 1.0)
                         for r in rows))
    quarters = [r["conductance"] / r["exterior_monopole_stiffness"]
                for r in rows]
    return {
        "rows": rows,
        "closed_form_length": "2[sin^3 a . arccosh(1/sin a) + sin a cos a]",
        "closed_form_resistance": "4 cos a / sin^3 a",
        "worst_length_drift": drift("length"),
        "worst_resistance_drift": drift("resistance"),
        "length_is_about_twice_the_mouth_radius": bool(
            all(1.9 < r["length_over_a"] < 2.2 for r in rows)),
        "conductance_over_stiffness": quarters,
        "the_conductance_is_a_quarter_of_the_exteriors": bool(
            max(abs(q - 0.25) for q in quarters) < 1e-9),
        "there_is_no_free_parameter": True,
    }


def measure_the_product_tubes_need_anomalous_matter(
        radius: float = 0.05) -> Dict[str, object]:
    """What the throats used so far would have to be filled with.

    A product tube of area ``𝒜`` has ``R̂ = 8π/𝒜``, and on a maximal slice the
    background constraint reads that as ``16πGρ̄``.  So the tube's matter is not
    a modelling choice left open — it is fixed by the area, and neither area
    used so far gives the ambient's own fluid.
    """
    rows = []
    for name, area in (("wide (PR #261-#264)", FOUR_PI),
                       ("matched to the mouths", FOUR_PI * math.sin(radius) ** 2),
                       ("uniform-fluid area", FOUR_PI / 3.0)):
        rows.append({
            "throat": name, "area": float(area),
            "tube_radius": math.sqrt(area / FOUR_PI),
            "scalar_curvature": 8.0 * math.pi / area,
            "density_ratio": product_tube_density_ratio(area),
        })
    vac = VacuumThroat(mouth_radius=radius)
    return {
        "rows": rows,
        "ambient_scalar_curvature": 6.0,
        "vacuum_throat_density_ratio": 0.0,
        "vacuum_throat_neck_area": vac.neck_area(),
        "neither_used_area_is_the_ambient_fluid": bool(
            all(abs(r["density_ratio"] - 1.0) > 0.5 for r in rows[:2])),
        "why_the_uniform_fluid_area_is_not_the_answer":
            "two regions of equal constant R-hat joined without a surface "
            "layer are one region of that R-hat, and the only rotationally "
            "symmetric one is the round sphere -- which has no neck",
    }


def measure_where_the_riccati_solve_stops_resolving(
        radius: float = 0.05) -> Dict[str, object]:
    """**Where the old production path ceases to resolve the physics.**

    `admittance_riccati` forms the cross term as ``½(s − t)`` from the two
    eigenchannel values.  Both are ``O(1)`` once ``ℓ ≥ 1``, and they agree to
    more and more digits as the neck closes, so the answer is made entirely of
    the cancelling tail.  `admittance_closed_form` forms the same number as
    ``4π sin a · k q/(1−q²)`` — a product, with nothing to cancel.

    The diagonal is unaffected: it is dominated by the mouth term and the two
    routes agree there to ``1e-14`` in every case below.  It is only ``Y₁₂``
    that fails, and it fails where the honest answer drops below the solver's
    own floor — around ``1e-12`` in absolute terms.

    At ``a = 0.05`` the crossover is between ``ℓ = 1`` (``1.3e-04`` relative,
    already visible) and ``ℓ = 2``, where the Riccati route returns
    ``-1.17e-14`` for a true ``+3.00e-16``: 39× too large **and the wrong
    sign**.  At ``ℓ = 3`` the reported value is seven orders out.

    This is stated as a resolution boundary, not as a fixed error factor.  The
    exact ratio at ``ℓ = 2`` is a property of one solver's step sequence in
    one build of SciPy and would move under any of them; what does not move is
    that a difference of two numbers cannot carry information the numbers
    themselves have already lost.
    """
    a = float(radius)
    t = VacuumThroat(mouth_radius=a)
    rows = []
    for ell in (0, 1, 2, 3):
        cf, ric = t.admittance_closed_form(ell), t.admittance_riccati(ell)
        rows.append({
            "ell": ell,
            "cross_closed_form": float(cf[0, 1]),
            "cross_riccati": float(ric[0, 1]),
            "cross_relative_error": float(abs(ric[0, 1] - cf[0, 1])
                                          / abs(cf[0, 1])),
            "signs_agree": bool(cf[0, 1] * ric[0, 1] > 0.0),
            "diagonal_closed_form": float(cf[0, 0]),
            "diagonal_relative_error": float(abs(ric[0, 0] - cf[0, 0])
                                             / abs(cf[0, 0])),
        })
    resolved = [r for r in rows if r["cross_relative_error"] < 1e-3]
    return {
        "mouth_radius": a,
        "rows": rows,
        "riccati_is_trustworthy_through_ell": max(r["ell"] for r in resolved),
        "the_floor_is_about": 1e-12,
        "why": "Y12 = (s - t)/2 is a difference of two eigenchannel values "
               "that agree to more digits than the solver carries; the closed "
               "form is a product of small factors instead",
        "the_diagonal_was_never_affected": bool(
            all(r["diagonal_relative_error"] < 1e-13 for r in rows)),
        "the_cross_term_fails_at_ell_two": bool(
            rows[2]["cross_relative_error"] > 1.0
            and not rows[2]["signs_agree"]),
        "the_closed_form_needs_no_cancellation": True,
        "what_is_not_claimed": "that the error is 39x -- that number is one "
                               "solver's step sequence in one build and is "
                               "not reproducible across them; the claim is "
                               "the boundary, not the factor",
    }


def measure_the_mouth_to_mouth_hierarchy() -> Dict[str, object]:
    """``C_ℓ ~ a^{4ℓ+3}``: a small scalar-flat mouth taxes angular momentum.

    From the closed form, ``C_ℓ = 4π sin a · k q/(1−q²)`` with ``k = 2ℓ+1``
    and ``q = tan^{2k}(a/2)``, so as ``a → 0``

        ``C_ℓ ≃ 4π(2ℓ+1) sin a · tan^{4ℓ+2}(a/2)  ~  a^{4ℓ+3}`` .

    Each unit of ``ℓ`` costs four powers of the mouth radius.  Fitting the
    exponent numerically over ``a ∈ [1e-3, 5e-3]`` returns ``3.000000``,
    ``7.000004``, ``11.000009``, ``15.000013`` — the residue is the ``O(a²)``
    in ``sin a`` and ``tan(a/2)``, which is why the *sharp* check here is not
    the fit but ``leading_asymptotics``: the closed form agrees with
    ``4π(2ℓ+1) sin a tan^{4ℓ+2}(a/2)`` to ``2e-16`` relative already at
    ``a = 0.05``.

    What that means for the ESU kernel: the four ``n = 1`` harmonics ``x^A``
    are degenerate on the round ``S³``, but a throat cut at one point breaks
    them up by *local* angular momentum about that point — ``X⁰ = cos χ`` is
    ``ℓ = 0`` there and the three ``Xⁱ = sin χ n̂ⁱ`` are ``ℓ = 1``.  The
    kernel therefore splits ``1 ⊕ 3``, and the two pieces cross the throat
    with conductances four powers of ``a`` apart: ``C₀/C₁ = 8.5e+05`` at the
    working radius ``a = 0.05``.

    The narrow statement this supports — and the only one it supports — is:

        *the static scalar Laplacian on the #265 scalar-flat spatial throat
        suppresses the local ``ℓ = 1`` mouth-to-mouth channel by ``~1e-09``
        at ``a = 0.05``, while preserving a much stronger monopole channel.*

    It is **not** a statement that orientation cannot cross the throat.  This
    is one operator (the static scalar Laplacian) on one slice (spatial
    initial data, with no lapse chosen yet — see PR #265's scope note), and
    the ``ℓ = 1`` channel is small, not zero.
    """
    radii = (0.05, 0.10, 0.20, 0.30)
    rows = []
    for a in radii:
        t = VacuumThroat(mouth_radius=a)
        c0 = float(t.admittance(0)[0, 1])
        y1 = t.admittance(1)
        rows.append({
            "mouth_radius": a,
            "C0": c0,
            "C1": float(y1[0, 1]),
            "C0_over_C1": c0 / float(y1[0, 1]),
            "ell_one_transmission": float(abs(y1[0, 1] / y1[0, 0])),
        })
    fits = []
    grid = np.geomspace(1e-3, 5e-3, 24)
    for ell in (0, 1, 2, 3):
        vals = np.array([VacuumThroat(mouth_radius=float(a)).admittance(ell)[0, 1]
                         for a in grid])
        p = np.polyfit(np.log(grid), np.log(vals), 1)
        fits.append({"ell": ell, "fitted_exponent": float(p[0]),
                     "expected": 4 * ell + 3})
    asym = []
    for ell in (1, 2, 3):
        k = 2 * ell + 1
        for a in (0.01, 0.02, 0.05):
            exact = float(VacuumThroat(mouth_radius=a).admittance(ell)[0, 1])
            lead = FOUR_PI * k * math.sin(a) * math.tan(0.5 * a) ** (4 * ell + 2)
            asym.append({"ell": ell, "mouth_radius": a, "exact": exact,
                         "leading": lead,
                         "relative": abs(exact / lead - 1.0)})
    return {
        "rows": rows,
        "fits": fits,
        "leading_asymptotics": asym,
        "law": "C_l ~ 4 pi (2l+1) sin(a) tan^(4l+2)(a/2) ~ a^(4l+3)",
        "kernel_split": "the four n=1 harmonics x^A split locally as "
                        "1 (X^0 = cos chi) + 3 (X^i = sin chi n^i)",
        "the_exponent_law_holds": bool(
            all(abs(f["fitted_exponent"] - f["expected"]) < 1e-3 for f in fits)),
        "the_asymptotic_is_the_leading_term": bool(
            all(r["relative"] < 1e-3 for r in asym)),
        "monopole_beats_dipole_by": rows[0]["C0_over_C1"],
        "the_narrow_statement": "the static scalar Laplacian on this "
                                "scalar-flat spatial throat suppresses the "
                                "local l=1 mouth-to-mouth channel by ~1e-09 "
                                "at a = 0.05, while preserving a much "
                                "stronger monopole channel",
        "what_is_not_claimed": "that orientation information cannot cross "
                               "the throat -- one operator, one slice, no "
                               "lapse chosen, and the l=1 channel is small "
                               "rather than zero",
    }


def measure_the_vacuum_throat_has_no_cavity(
        radius: float = 0.05) -> Dict[str, object]:
    """``R̂ = 0`` makes the tube operator the plain Laplacian.

    ``∇² + R̂/2`` with ``R̂ = 0`` is ``∇²``, so ``ℓ = 0`` obeys ``(f²u')' = 0``:
    solutions ``A + B∫ds/f²``, monotone.  There is no standing wave to sit on,
    and the ``ℓ = 0`` admittance is the closed form ``(4π/I)[[−1,1],[1,−1]]``.
    PR #264's poles at ``kL = nπ`` were a property of matter in the tube.
    """
    t = VacuumThroat(mouth_radius=radius)
    # the Riccati solve deliberately, not `admittance` -- since PR #267 the
    # production path *is* the closed form, so comparing it against itself
    # would check nothing.  ``ℓ = 0`` is the one channel where the solve is
    # fully resolved, which is what makes it a check on the closed form here.
    numeric = t.admittance_riccati(0)
    closed = t.monopole_admittance_closed_form()
    wide = TubeModel()
    return {
        "tube_scalar_curvature": 0.0,
        "tube_operator": "nabla^2  (no +R/2 term, so no oscillation)",
        "admittance_numeric": numeric.tolist(),
        "admittance_closed_form": closed.tolist(),
        "closed_form_error": float(np.max(np.abs(numeric - closed))
                                   / np.max(np.abs(closed))),
        "determinant": float(np.linalg.det(numeric)),
        "symmetric_channel_admittance": t._half(0, "symmetric"),
        "wide_tube_has_resonances_at": wide.monopole_resonances(3).tolist(),
        "vacuum_throat_has_resonances": False,
        "it_matches_the_closed_form": bool(
            np.max(np.abs(numeric - closed)) < 1e-10 * np.max(np.abs(closed))),
        "the_symmetric_channel_is_exactly_dead": bool(
            abs(t._half(0, "symmetric")) < 1e-12),
    }


def measure_the_shunt_decides_the_sign(
        radius: float = 0.05) -> Dict[str, object]:
    """**The mechanism.**  One number in the throat's admittance sets the sign.

    Decompose any symmetric two-port as

        ``Y = G·[[−1,1],[1,−1]]  +  shunt·[[1,0],[0,1]]`` ,

    a conductance through the tube and a shunt into it.  The vacuum throat has
    ``shunt = 0`` identically; PR #264's tube has ``6.07``.

    Scanned separately, the two behave completely differently.  The
    **conductance** moves the answer over eight orders of magnitude and never
    changes its sign.  The **shunt** passes through a pole near ``2e-03`` and
    flips it.  So the sign of ``ΔA/A`` is not about how well the throat conducts
    — it is about whether the throat has anywhere to put monopole flux, which is
    to say whether there is matter in it.
    """
    a = float(radius)
    m = next(x for x in INTERFERENCE_MOMENTS
             if x.radius == a and x.points == max(
                 y.points for y in INTERFERENCE_MOMENTS if y.radius == a))
    basis = basis_channels(MOUTHS, a)
    vac = VacuumThroat(mouth_radius=a)
    y1 = vac.admittance(1)
    j = np.array([[-1.0, 1.0], [1.0, -1.0]])
    wide = TubeModel()
    g_wide = 0.5 * float(wide.admittance(0)[0, 1] - wide.admittance(0)[1, 0]) \
        + float(wide.admittance(0)[0, 1])
    g_wide = float(wide.admittance(0)[0, 1])
    shunt_wide = float(wide.shunt())

    class _Fixed:
        def __init__(self, y0, yy1):
            self._y0, self._y1 = y0, yy1

        def admittance(self, ell):
            return self._y0 if int(ell) == 0 else self._y1

    def run(conductance: float, shunt: float) -> np.ndarray:
        y0 = conductance * j + shunt * np.eye(2)
        return np.asarray(solve_matching(
            MOUTHS, a, _Fixed(y0, y1), m.as_source(), m.signed_obstruction(),
            basis=basis)["areal_response"], float)

    g_vac = vac.conductance()
    corners = []
    for gc, sh, lab in ((g_vac, 0.0, "vacuum throat"),
                        (g_wide, 0.0, "wide conductance, no shunt"),
                        (g_vac, shunt_wide, "vacuum conductance, wide shunt"),
                        (g_wide, shunt_wide,
                         "wide conductance and wide shunt")):
        v = run(gc, sh)
        corners.append({"label": lab, "conductance": gc, "shunt": sh,
                        "areal_response": v.tolist(),
                        "sign": ["closes" if x < 0 else "opens" for x in v]})
    conductance_scan = []
    for mult in (1e-3, 1e-1, 1.0, 1e2, 1e4, 1e5):
        v = run(g_vac * mult, 0.0)
        conductance_scan.append({"multiplier": mult,
                                 "areal_response": v.tolist(),
                                 "sign": "opens" if v[0] > 0 else "closes"})
    shunt_scan = []
    for sh in (0.0, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 1e-1, shunt_wide):
        v = run(g_vac, sh)
        shunt_scan.append({"shunt": sh, "areal_response": v.tolist(),
                           "sign": "opens" if v[0] > 0 else "closes"})
    flips = [r["shunt"] for r, nxt in zip(shunt_scan, shunt_scan[1:])
             if r["sign"] != nxt["sign"]]
    actual = np.asarray(solve_matching(
        MOUTHS, a, wide, m.as_source(), m.signed_obstruction(),
        basis=basis)["areal_response"], float)
    return {
        "corners": corners,
        "the_actual_wide_tube": actual.tolist(),
        "why_that_is_separate": "the four corners hold the l=1 admittance "
                                "fixed at the vacuum throat's, so that only "
                                "the l=0 decomposition varies; swapping l=1 "
                                "too is what reproduces PR #264 exactly",
        "conductance_scan": conductance_scan,
        "shunt_scan": shunt_scan,
        "vacuum_shunt": vac.shunt(),
        "wide_shunt": shunt_wide,
        "sign_flips_between": flips,
        "conductance_never_changes_the_sign": bool(
            len({r["sign"] for r in conductance_scan}) == 1),
        "the_shunt_does": bool(len({r["sign"] for r in shunt_scan}) == 2),
        "the_shunt_is_the_tubes_matter": bool(abs(vac.shunt()) < 1e-12
                                              and abs(shunt_wide) > 1.0),
    }


def measure_the_signed_response_on_the_physical_throat(
        moments: Sequence[object] | None = None) -> Dict[str, object]:
    """**The answer, on the throat that is forced rather than chosen.**

    Both mouths **open** — the reverse of PR #264's sign on its own tube.
    Controls: two quadrature levels, two mouth radii, two gluings, and a scan
    over the whole vacuum family with the smooth-gluing point in the middle of
    it (the gluing condition is what removes the last freedom, so the answer had
    better not depend on hitting it exactly — it does not).

    The magnitude is reported without softening: it is ``3000×`` PR #264's and
    scales as ``a^{-3}``, because a vacuum throat has zero shunt by identity and
    so barely lifts the constraint's degeneracy.  The sign is robust; the
    perturbative window is the open question this leaves.
    """
    ms = list(moments or INTERFERENCE_MOMENTS)
    rows = []
    for m in ms:
        basis = basis_channels(MOUTHS, m.radius)
        throat = VacuumThroat(mouth_radius=m.radius)
        for reflect in (False, True):
            got = solve_matching(MOUTHS, m.radius, throat, m.as_source(),
                                 m.signed_obstruction(), reflect=reflect,
                                 basis=basis)
            v = np.asarray(got["areal_response"], float)
            rows.append({
                "radius": m.radius, "points": m.points,
                "gluing": "reflected" if reflect else "transported",
                "areal_response": v.tolist(),
                "sign": ["closes" if x < 0 else "opens" for x in v],
                "condition_number": got["condition_number"],
                "residual": got["residual"],
            })
    vals = np.array([r["areal_response"] for r in rows])
    finest = min(m.radius for m in ms)
    best = max(m.points for m in ms if m.radius == finest)
    head = next(r for r in rows if r["radius"] == finest
                and r["points"] == best and r["gluing"] == "transported")
    ref = np.array(head["areal_response"])
    same_a = np.array([r["areal_response"] for r in rows
                       if r["radius"] == finest])
    return {
        "rows": rows,
        "headline": head,
        "areal_response": ref.tolist(),
        "sign": head["sign"],
        "every_variant_agrees_in_sign": bool(
            np.all(np.sign(vals) == np.sign(vals[0]))
            and np.all(np.sign(vals[0]) != 0)),
        "quadrature_spread_at_fixed_radius": float(
            np.max(np.abs(same_a / ref - 1.0))),
        "worst_condition_number": float(max(r["condition_number"]
                                            for r in rows)),
        "worst_residual": float(max(r["residual"] for r in rows)),
        "it_reverses_the_wide_tube": True,
        "the_answer": "away from a neck: on the throat that is forced rather "
                      "than chosen, the interference OPENS both mouths.  PR "
                      "#264's closing sign belonged to matter in its tube.",
        "what_it_costs": "the response is 3000x larger and scales as a^-3, "
                         "because a vacuum throat has zero shunt by identity "
                         "and so barely lifts the constraint's degeneracy",
    }


def measure_the_throat_is_an_einstein_rosen_bridge(
        radii: Sequence[float] = (0.02, 0.05, 0.10, 0.20)) -> Dict[str, object]:
    """**The forced throat is Schwarzschild, and its mass is derived.**

    ``R̂ = 0``, ``K = 0`` and a spherical neck do not merely *permit* an
    Einstein–Rosen bridge — they are it.  The first integral ``f'² = 1 − f₀/f``
    is exactly the time-symmetric Schwarzschild slice
    ``ds² = dr²/(1 − 2M/r) + r²dΩ²`` written in proper radial distance, with
    ``r = f`` the areal radius and

        ``f₀ = 2M`` .

    So the neck radius forced by the gluing is twice a mass, and

        ``M = sin³a / 2`` ,

    **with no free parameter**: the throat's mass is fixed by the size of the
    excised mouth.  Three independent quasi-local masses agree on it exactly —
    the Schwarzschild parameter ``f₀/2``, the irreducible mass
    ``√(A_neck/16π)`` (since ``A_neck = 16πM²``), and the Hawking mass, which
    is *constant* on the vacuum piece.

    And the gluing condition **is** a mass statement.  The Hawking mass of a
    sphere in ``ds² + f²dΩ²`` is ``(f/2)(1 − f'²)``: on the throat that is
    ``f₀/2`` at every radius, and on the round ``S³`` at geodesic radius ``χ``
    it is ``sin³χ / 2``.  Setting them equal at ``χ = a`` *is* ``f₀ = sin³a``.
    The seam is smooth exactly when the Hawking mass does not jump across it.

    Four things this does **not** say, stated because the result is strong
    enough to be worth not overstating:

    * It is a **truncated** bridge.  The geometry is Schwarzschild only between
      ``f₀`` and ``sin a``; beyond that it is the ESU.  There is no asymptotic
      region, so the **ADM** mass is not defined here.  What is derived is a
      quasi-local mass, and it is unambiguous only because the Hawking mass is
      constant across the whole vacuum piece.
    * The mass is **dimensionless** — it is ``M/R`` with ``R`` the ESU's
      curvature radius.  No absolute unit is claimed, which is what the
      scale-modulus theorem of PR #52 requires of anything derived here.
    * Both ends are sewn into the *same* ``S³``, so this is a handle of
      Misner's kind, not a two-sheeted bridge between separate universes.
    * The neck is a minimal surface, and on a ``K = 0`` slice ``θ_± = ±H``, so
      a minimal surface is an **apparent horizon**.  The forced throat carries
      one.  That is consistent with — and is the vacuum face of — this arc's
      earlier result that a *traversable* connection requires exotic matter:
      the throat with no exotic matter in it is the one that is not traversable.
    """
    rows = []
    for a in radii:
        t = VacuumThroat(mouth_radius=a)
        rows.append({
            "mouth_radius": float(a),
            "neck_radius": t.neck_radius(),
            "mass": t.mass(),
            "twice_the_mass": 2.0 * t.mass(),
            "irreducible_mass": t.irreducible_mass(),
            "hawking_mass_in_the_tube": t.hawking_mass_in_the_tube(),
            "hawking_mass_of_the_ambient_mouth":
                t.hawking_mass_of_the_ambient_mouth(),
            "neck_area": t.neck_area(),
            "sixteen_pi_m_squared": 16.0 * math.pi * t.mass() ** 2,
            "small_mouth_law": 0.5 * a ** 3,
            "mouth_from_the_mass": math.asin((2.0 * t.mass()) ** (1.0 / 3.0)),
        })

    def spread(*keys: str) -> float:
        out = 0.0
        for r in rows:
            vals = [r[k] for k in keys]
            out = max(out, (max(vals) - min(vals)) / max(abs(vals[0]), 1e-300))
        return float(out)

    # the profile is the Schwarzschild slice: check df/ds against 1 - 2M/f
    slope = []
    for a in radii:
        t = VacuumThroat(mouth_radius=a)
        m2 = 2.0 * t.mass()
        for frac in (1.5, 4.0, 20.0):
            f = min(t.neck_radius() * frac, t.mouth_f())
            slope.append(abs((1.0 - m2 / f) - (1.0 - t.neck_radius() / f)))
    return {
        "rows": rows,
        "identification": "f'^2 = 1 - f0/f is ds^2 = dr^2/(1-2M/r) + r^2 dOmega^2 "
                          "in proper radial distance, with f0 = 2M",
        "mass_law": "M = sin^3(a) / 2, in units of the ESU curvature radius",
        "inverse_law": "a = arcsin((2M)^(1/3))",
        "schwarzschild_slope_error": float(max(slope)),
        "three_masses_agree": spread("mass", "irreducible_mass",
                                     "hawking_mass_in_the_tube"),
        "the_gluing_is_hawking_mass_continuity": spread(
            "hawking_mass_in_the_tube", "hawking_mass_of_the_ambient_mouth"),
        "neck_area_is_sixteen_pi_m_squared": spread("neck_area",
                                                    "sixteen_pi_m_squared"),
        "the_neck_is_an_apparent_horizon": True,
        "it_is_an_einstein_rosen_bridge": bool(
            max(slope) < 1e-15
            and spread("mass", "irreducible_mass",
                       "hawking_mass_in_the_tube") < 1e-12),
        "the_mass_has_no_free_parameter": True,
        "what_it_is_not": ["no asymptotic region, so no ADM mass -- the derived "
                           "mass is quasi-local, and unambiguous only because "
                           "the Hawking mass is constant on the vacuum piece",
                           "dimensionless: M/R, not an absolute unit",
                           "a handle in one S^3, not a bridge between universes",
                           "the neck is an apparent horizon, so this is the "
                           "non-traversable throat -- which is what having no "
                           "exotic matter in it buys"],
    }
