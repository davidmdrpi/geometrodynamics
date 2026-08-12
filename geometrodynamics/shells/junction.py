"""
Darmois–Israel junctions: what a detached shell's gluing forces about its stress.

SCOPE OF EVERY CLAIM BELOW
──────────────────────────
**Einstein gravity, Darmois–Israel thin shells, spherical symmetry, vacuum
bulk.**  Nothing here bounds thick shells, non-spherical configurations,
modified gravity, or shells with non-vacuum interiors.  ``G = 1``.  The spatial
dimension is a parameter: ``D = 4`` is the regression case and ``D = 5``
(Tangherlini) is the one this program actually cares about.

THE QUESTION
────────────
A wormhole throat needs negative surface energy.  The hope is that a *detached*
closed shell in the bulk, glued with the opposite orientation, might supply the
exotic-looking restoring stress while itself being ordinary matter.

That hope decomposes into **three observables that are logically independent**:

1. **Does the shell itself require exotic matter?**  ``σ = −S^τ_τ`` from
   ``Sⁱⱼ = −(1/8πG)([Kⁱⱼ] − δⁱⱼ[K])``.
2. **Does it support the throat?**  The shell's contribution to the throat's
   potential gradient.
3. **Is the configuration stable?**  The stiffness ``V''(b₀)``.

A positive stiffness alone would mean "restoring" and would *not* establish that
the shell supplied it, so all three are reported separately.

THE ORIENTATION IS DERIVED FROM THE GLUING, NOT ASSUMED
───────────────────────────────────────────────────────
Earlier drafts carried ``ε = ±1`` as a free flag, which made the central result
conditional on a hand-set sign.  Here each side of the surface is a **branch** —
which radial half of its region is retained — and ``ε`` follows:

* the ``−`` side is approached from within along the ``−→+`` normal, so
  ``ε₋ = +1`` if it retains ``r ≤ R`` (INNER) and ``−1`` if it retains ``r ≥ R``
  (OUTER);
* the ``+`` side is entered leaving the surface, so ``ε₊ = +1`` for OUTER and
  ``−1`` for INNER.

There are therefore **four** gluings, not two, and with
``σ = −(D−2)(ε₊β₊ − ε₋β₋)/(8πG R)``:

============  ==========  ============  ===================================
``−`` branch  ``+`` branch  ``η = ε₊ε₋``  what it is, and what ``σ`` must be
============  ==========  ============  ===================================
INNER         OUTER       ``+1``        ordinary bubble — either sign
OUTER         OUTER       ``−1``        **minimal surface** — ``σ < 0`` always
INNER         INNER       ``−1``        **maximal surface** — ``σ > 0`` always
OUTER         INNER       ``+1``        anti-bubble — either sign
============  ==========  ============  ===================================

So ``η = −1`` alone decides nothing: it covers two gluings whose signs are
*opposite and both forced*.  The sign is a property of the branch pair, not of
``η``.

WHAT IS ACTUALLY FORCED
───────────────────────
For a **minimal surface** — ``r`` increasing away on both sides, which is what a
throat is — the two terms add:

    ``σ = −(D−2)(β₊ + β₋)/(8πG R) < 0``

with ``β± = √(f± + Ṙ²) ≥ 0`` for any timelike shell.  For a **maximal surface**
they add with the other sign and ``σ > 0`` always.  Both are identities, in
every ``D``.

The dichotomy this produces is the real result, and it is sharper than "the
oppositely-glued shell is exotic":

    a detached surface that **connects** to the throat's asymptotic region
    through a minimal surface is necessarily exotic; a detached surface that is
    non-exotic by its gluing is a **maximal** surface, which caps off on both
    sides and is therefore a closed region sharing no bulk with the throat — so
    it cannot support anything.

Within Einstein–Israel spherical thin shells, exotic matter is relocated, never
removed.

WHAT AN ORDINARY BUBBLE CAN DO
──────────────────────────────
An aligned bubble can be ordinary and does shift the throat's potential by
screening mass.  That shift is reported as a **potential-gradient contribution**
rather than a force: it is taken at fixed throat rest mass and omits the
equation-of-state response, so it is not an equilibrium-consistent force.  The
radial acceleration contribution is ``−½ ∂ΔV/∂b`` from ``ḃ² + V = 0``, and both
are reported.

STABILITY
─────────
``ḃ² + V(b) = 0``; a static solution needs ``V(b₀) = 0`` and ``V'(b₀) = 0`` and
is stable when ``V''(b₀) > 0``.  A fixed global barotropic index admits no
stable static throat at all — ``V''(b₀) = 2GM(n−1)/b₀³`` with ``n = 2 + 4w`` in
``D = 4``, so static needs ``w < −1/2`` and stability ``w > −1/4`` — so
``β² ≡ (dp/dσ)|₀`` is left free at the equilibrium, as usual.  ``β²`` is an
**equation-of-state derivative parameter**; for exotic matter its sign carries
no sound-speed reading and none is claimed.

Because Birkhoff makes the stiffness matrix diagonal in this model, its
eigenvalues are just the two stiffnesses.  They are **not** normal-mode
frequencies: no kinetic metric has been derived, so nothing here is a
generalised eigenproblem, and the quantities are named stiffnesses throughout.

BIRKHOFF
────────
In this vacuum spherical model the region between the two surfaces is
Schwarzschild/Tangherlini with a constant mass parameter, so the throat's data
cannot depend on where the shell sits.  That establishes **no
separation-dependent coupling in this model** — not that every spherical trapped
resonator is impossible.  The vanishing ``∂²V/∂a∂b`` is structural, since
Birkhoff is imported the moment the intervening region is written that way;
what is *measured* is that a genuinely different family of shells (surface
density varying by hundreds) leaves the throat bit-for-bit unchanged.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "INNER",
    "OUTER",
    "Gluing",
    "MINIMAL_SURFACE",
    "MAXIMAL_SURFACE",
    "ORDINARY_BUBBLE",
    "ANTI_BUBBLE",
    "GLUINGS",
    "DetachedShell",
    "Region",
    "ThroatShellSystem",
    "Z2Throat",
    "surface_stress",
    "measure_the_junction_reproduces_known_shells",
    "measure_the_gluing_fixes_the_sign",
    "measure_the_forced_signs_hold_in_any_dimension",
    "measure_the_detached_shell_can_be_ordinary",
    "measure_the_shell_potential_gradient",
    "measure_the_throat_and_shell_are_decoupled",
    "measure_the_stability_window",
    "measure_the_three_observables_are_independent",
    "measure_the_stiffnesses_scale_dimensionally",
]

G = 1.0
INNER = "INNER"     # the branch r ≤ R is retained on this side
OUTER = "OUTER"     # the branch r ≥ R is retained on this side


class Gluing:
    """Which radial branch each side of the surface retains.

    ``ε`` is *derived* from this, not supplied: the ``−`` side is approached
    from within along the ``−→+`` normal and the ``+`` side is entered leaving
    the surface, which fixes both signs with no freedom left.
    """

    def __init__(self, minus: str, plus: str) -> None:
        for b in (minus, plus):
            if b not in (INNER, OUTER):
                raise ValueError(f"branch must be INNER or OUTER, got {b!r}")
        self.minus, self.plus = minus, plus

    @property
    def eps_minus(self) -> int:
        return +1 if self.minus == INNER else -1

    @property
    def eps_plus(self) -> int:
        return +1 if self.plus == OUTER else -1

    @property
    def eta(self) -> int:
        return self.eps_minus * self.eps_plus

    @property
    def name(self) -> str:
        return {(OUTER, OUTER): "minimal surface",
                (INNER, INNER): "maximal surface",
                (INNER, OUTER): "ordinary bubble",
                (OUTER, INNER): "anti-bubble"}[(self.minus, self.plus)]

    @property
    def forced_sign(self) -> Optional[int]:
        """``−1``/``+1`` where the gluing alone fixes ``sign(σ)``, else ``None``."""
        if (self.minus, self.plus) == (OUTER, OUTER):
            return -1
        if (self.minus, self.plus) == (INNER, INNER):
            return +1
        return None

    @property
    def connects_to_infinity_on_both_sides(self) -> bool:
        return self.minus == OUTER and self.plus == OUTER

    def __repr__(self) -> str:
        return f"Gluing({self.minus}->{self.plus}, {self.name})"


MINIMAL_SURFACE = Gluing(OUTER, OUTER)
MAXIMAL_SURFACE = Gluing(INNER, INNER)
ORDINARY_BUBBLE = Gluing(INNER, OUTER)
ANTI_BUBBLE = Gluing(OUTER, INNER)
GLUINGS = [ORDINARY_BUBBLE, MINIMAL_SURFACE, MAXIMAL_SURFACE, ANTI_BUBBLE]


# ════════════════════════════════════════════════════════════════════════════
# THE BULK
# ════════════════════════════════════════════════════════════════════════════
class Region:
    """A static spherically symmetric vacuum region in ``D`` dimensions.

    ``f(r) = 1 − μ/r^{D−3} + q²/r^{2(D−3)} − 2Λr²/((D−1)(D−2))`` — Tangherlini,
    reducing to Schwarzschild / Reissner–Nordström / Schwarzschild–de Sitter at
    ``D = 4``.  Parametrised by the **mass parameter** ``μ`` rather than an ADM
    mass, so no ``Ω_{D−2}`` convention has to be adopted; at ``D = 4``,
    ``μ = 2GM``.
    """

    def __init__(self, mu: float = 0.0, charge: float = 0.0,
                 lambda_: float = 0.0, dim: int = 4) -> None:
        if dim < 4:
            raise ValueError("dim must be at least 4")
        self.mu = float(mu)
        self.charge = float(charge)
        self.lambda_ = float(lambda_)
        self.dim = int(dim)

    @classmethod
    def from_mass(cls, mass: float, dim: int = 4, **kw) -> "Region":
        """``D = 4`` convenience: ``μ = 2GM``."""
        if dim != 4:
            raise ValueError("from_mass is a D=4 convenience; pass mu for D>4")
        return cls(mu=2.0 * G * mass, dim=4, **kw)

    @property
    def n(self) -> int:
        """``D − 3``, the radial falloff exponent."""
        return self.dim - 3

    def f(self, r: float) -> float:
        d = self.dim
        return (1.0 - self.mu / r ** self.n
                + self.charge ** 2 / r ** (2 * self.n)
                - 2.0 * self.lambda_ * r * r / ((d - 1) * (d - 2)))

    def df(self, r: float) -> float:
        d, n = self.dim, self.n
        return (n * self.mu / r ** (n + 1)
                - 2 * n * self.charge ** 2 / r ** (2 * n + 1)
                - 4.0 * self.lambda_ * r / ((d - 1) * (d - 2)))

    def d2f(self, r: float) -> float:
        d, n = self.dim, self.n
        return (-n * (n + 1) * self.mu / r ** (n + 2)
                + 2 * n * (2 * n + 1) * self.charge ** 2 / r ** (2 * n + 2)
                - 4.0 * self.lambda_ / ((d - 1) * (d - 2)))

    def beta(self, r: float, rdot: float = 0.0) -> float:
        v = self.f(r) + rdot * rdot
        if v < 0.0:
            raise ValueError(f"shell is not timelike at r={r}: f + Ṙ² = {v}")
        return math.sqrt(v)


# ════════════════════════════════════════════════════════════════════════════
# THE JUNCTION
# ════════════════════════════════════════════════════════════════════════════
def surface_stress(minus: Region, plus: Region, gluing: Gluing, r: float,
                   rdot: float = 0.0) -> Dict[str, object]:
    """``σ = −(D−2)[K^θ_θ]/(8πG)``, with the orientation taken from ``gluing``.

    ``K^θ_θ = (ε/R)√(f + Ṙ²)`` in every ``D``; only the ``(D−2)`` trace weight
    carries the dimension.
    """
    if minus.dim != plus.dim:
        raise ValueError("both sides must have the same dimension")
    d = minus.dim
    b_m, b_p = minus.beta(r, rdot), plus.beta(r, rdot)
    jump = (gluing.eps_plus * b_p - gluing.eps_minus * b_m) / r
    sigma = -(d - 2) * jump / (8.0 * math.pi * G)
    return {"sigma": sigma, "jump_K_theta": jump, "dim": d,
            "beta_minus": b_m, "beta_plus": b_p,
            "eps": (gluing.eps_minus, gluing.eps_plus),
            "eta": gluing.eta, "gluing": gluing.name,
            "forced_sign": gluing.forced_sign,
            "is_exotic": bool(sigma < 0.0)}


def _potential(minus: Region, plus: Region, gluing: Gluing, r: float,
               sigma: float) -> float:
    """``V(R) = f₊ − A²`` with ``A = ε₊β₊`` solved from the junction.

    ``A − B = −K`` and ``A² − B² = f₊ − f₋`` (since ``ε² = 1``) separate the two
    roots exactly, so no square root survives to be differentiated and every
    gluing is handled by the same expression.
    """
    d = minus.dim
    k = 8.0 * math.pi * G * r * sigma / (d - 2)
    if abs(k) < 1e-300:
        raise ValueError("massless shell: the junction does not determine Ṙ")
    df = plus.f(r) - minus.f(r)
    a = 0.5 * (-k - df / k)
    return plus.f(r) - a * a


def _sigma_series(r0: float, sigma0: float, p0: float, beta2: float,
                  d: int) -> Tuple[float, float]:
    """``σ'`` and ``σ''`` at ``r₀`` from conservation and ``β² = dp/dσ``.

    ``σ' = −(D−2)(σ + p)/R``, differentiated once more with ``p' = β²σ'``.
    """
    sp = -(d - 2) * (sigma0 + p0) / r0
    spp = -(d - 2) * sp * (1.0 + beta2) / r0 + (d - 2) * (sigma0 + p0) / r0 ** 2
    return sp, spp


def _stiffness(minus: Region, plus: Region, gluing: Gluing, r0: float,
               sigma0: float, p0: float, beta2: float) -> float:
    """``V''(r₀)`` by Richardson-extrapolated differences of the exact model.

    ``σ`` is known at ``r₀`` to second order, which is exactly what ``V''``
    needs, so the local quadratic is not an approximation to the answer — only
    the differencing is, and Richardson removes its leading error.
    """
    d = minus.dim
    sp, spp = _sigma_series(r0, sigma0, p0, beta2, d)

    def v(dr: float) -> float:
        s = sigma0 + sp * dr + 0.5 * spp * dr * dr
        return _potential(minus, plus, gluing, r0 + dr, s)

    def second(h: float) -> float:
        return (v(h) - 2.0 * v(0.0) + v(-h)) / (h * h)

    h = 1e-3 * r0
    return (4.0 * second(0.5 * h) - second(h)) / 3.0


def _static_data(minus: Region, plus: Region, gluing: Gluing, r0: float
                 ) -> Dict[str, float]:
    """``σ₀`` from ``V(r₀) = 0`` and ``p₀`` from ``V'(r₀) = 0``.

    ``σ₀`` comes straight from the junction; ``p₀`` is found by solving the
    linear condition ``∂V/∂σ' = 0`` numerically, which keeps one code path for
    all four gluings instead of four hand-derived expressions.
    """
    d = minus.dim
    b_m, b_p = minus.beta(r0), plus.beta(r0)
    jump = (gluing.eps_plus * b_p - gluing.eps_minus * b_m) / r0
    sigma0 = -(d - 2) * jump / (8.0 * math.pi * G)
    if abs(sigma0) < 1e-300:
        raise ValueError("degenerate junction: σ₀ = 0")

    def vprime(p: float) -> float:
        sp = -(d - 2) * (sigma0 + p) / r0
        h = 1e-6 * r0

        def v(dr: float) -> float:
            return _potential(minus, plus, gluing, r0 + dr, sigma0 + sp * dr)
        return (v(h) - v(-h)) / (2.0 * h)

    # V' is linear in p only in the h → 0 limit; at finite h it keeps an
    # O(h·p²) term, so probing at p = 1 — thousands of times the physical scale
    # σ₀ — lets that contamination dominate the extrapolation and returns a
    # root that is wrong by a factor of three.  Probe at the scale of σ₀ and
    # then iterate a secant, which converges on the true root in a few steps.
    scale = abs(sigma0)
    a, b = 0.0, scale
    fa, fb = vprime(a), vprime(b)
    if abs(fb - fa) < 1e-300:
        raise ValueError("V' does not depend on p: degenerate configuration")
    for _ in range(60):
        if abs(fb - fa) < 1e-300:
            break
        c = b - fb * (b - a) / (fb - fa)
        a, fa, b = b, fb, c
        fb = vprime(c)
        if abs(fb) < 1e-18 * max(1.0, abs(fa)):
            break
    p0 = b
    return {"sigma": sigma0, "p": p0, "residual_Vprime": vprime(p0),
            "is_exotic": bool(sigma0 < 0.0)}


# ════════════════════════════════════════════════════════════════════════════
# SURFACES
# ════════════════════════════════════════════════════════════════════════════
class Z2Throat:
    """A ``Z2``-symmetric throat: a **minimal surface**, so ``σ < 0`` is forced.

    The gluing is not configurable — a throat retains the OUTER branch on both
    sides by definition, and that is what makes it exotic.
    """

    def __init__(self, mu: float = 2.0, b0: float = 5.0, dim: int = 4) -> None:
        self.region = Region(mu=mu, dim=dim)
        self.b0 = float(b0)
        self.dim = int(dim)
        self.gluing = MINIMAL_SURFACE
        if self.region.f(self.b0) <= 0.0:
            raise ValueError("throat must sit outside the horizon")

    def static(self) -> Dict[str, float]:
        return _static_data(self.region, self.region, self.gluing, self.b0)

    def stiffness(self, beta2: float) -> float:
        s = self.static()
        return _stiffness(self.region, self.region, self.gluing, self.b0,
                          s["sigma"], s["p"], beta2)

    def stability_window(self) -> Dict[str, float]:
        v0, v1 = self.stiffness(0.0), self.stiffness(1.0)
        slope = v1 - v0
        crit = (-v0 / slope) if abs(slope) > 1e-14 else float("nan")
        return {"stiffness_at_zero": v0, "slope": slope,
                "beta2_critical": crit,
                "stable_below_critical": bool(slope < 0.0),
                "window_needs_negative_beta2": bool(
                    not math.isnan(crit) and crit < 0.0)}

    def potential(self, b: float, mu: Optional[float] = None) -> float:
        """``V(b)`` continued at **fixed rest mass**, for the screening shift.

        The fixed-rest-mass continuation omits the equation-of-state response,
        so this is a potential shift and not an equilibrium-consistent
        trajectory; it is used only to compare two bulk mass parameters.
        """
        reg = self.region if mu is None else Region(mu=mu, dim=self.dim)
        s = self.static()
        u0 = -4.0 * math.pi * G * self.b0 * s["sigma"] / (self.dim - 2)
        u = u0 * self.b0 / b
        return reg.f(b) - u * u


class DetachedShell:
    """A closed shell with an explicit gluing between two vacuum regions."""

    def __init__(self, minus: Region, plus: Region, a: float,
                 gluing: Gluing = ORDINARY_BUBBLE) -> None:
        if minus.dim != plus.dim:
            raise ValueError("both sides must have the same dimension")
        self.minus, self.plus = minus, plus
        self.a = float(a)
        self.gluing = gluing

    def stress(self, adot: float = 0.0) -> Dict[str, object]:
        return surface_stress(self.minus, self.plus, self.gluing, self.a, adot)

    def static(self) -> Dict[str, float]:
        """Valid for **every** gluing — the junction is solved generically."""
        return _static_data(self.minus, self.plus, self.gluing, self.a)

    def stiffness(self, beta2: float) -> float:
        s = self.static()
        return _stiffness(self.minus, self.plus, self.gluing, self.a,
                          s["sigma"], s["p"], beta2)

    def stability_window(self) -> Dict[str, float]:
        v0, v1 = self.stiffness(0.0), self.stiffness(1.0)
        slope = v1 - v0
        crit = (-v0 / slope) if abs(slope) > 1e-14 else float("nan")
        return {"stiffness_at_zero": v0, "slope": slope,
                "beta2_critical": crit}

    def screened_mu(self) -> float:
        return self.plus.mu - self.minus.mu


# ════════════════════════════════════════════════════════════════════════════
# THE PAIR
# ════════════════════════════════════════════════════════════════════════════
class ThroatShellSystem:
    """Throat at ``b`` and detached shell at ``a > b`` in a common bulk."""

    def __init__(self, mu_outer: float = 2.0, screened: float = 0.6,
                 b0: float = 5.0, a0: float = 20.0, dim: int = 4,
                 gluing: Gluing = ORDINARY_BUBBLE) -> None:
        self.dim = int(dim)
        self.mu_outer = float(mu_outer)
        self.mu_inner = float(mu_outer) - float(screened)
        self.b0, self.a0 = float(b0), float(a0)
        self.inner = Region(mu=self.mu_inner, dim=dim)
        self.outer = Region(mu=self.mu_outer, dim=dim)
        self.throat = Z2Throat(mu=self.mu_inner, b0=b0, dim=dim)
        self.shell = DetachedShell(self.inner, self.outer, a0, gluing)

    def shell_potential_shift(self, b: float) -> float:
        return (self.throat.potential(b, self.mu_inner)
                - self.throat.potential(b, self.mu_outer))

    def shell_potential_gradient(self, b: Optional[float] = None,
                                 h: float = 1e-5) -> Dict[str, float]:
        """``−∂ΔV/∂b``, and the acceleration contribution ``−½∂ΔV/∂b``.

        Not an equilibrium-consistent force: the continuation holds the
        throat's rest mass fixed and omits its equation-of-state response.  It
        is the gradient of the potential *shift* the shell's screening produces,
        which is what the sign statement is about.
        """
        b = self.b0 if b is None else b
        grad = (self.shell_potential_shift(b + h)
                - self.shell_potential_shift(b - h)) / (2.0 * h)
        return {"minus_dV_db": -grad,
                "acceleration_contribution": -0.5 * grad,
                "opposes_closure": bool(-grad > 0.0)}

    def stiffness_matrix(self, beta2_throat: float = -1.0,
                         beta2_shell: float = 0.5) -> np.ndarray:
        """``diag`` of the two stiffnesses, with the off-diagonal computed.

        **Not** a normal-mode matrix: no kinetic metric has been derived, so
        these are stiffnesses and their eigenvalues are not frequencies.
        """
        h = 1e-5 * self.a0

        def throat_stiff_at(a: float) -> float:
            probe = ThroatShellSystem(
                mu_outer=self.mu_outer,
                screened=self.mu_outer - self.mu_inner,
                b0=self.b0, a0=a, dim=self.dim, gluing=self.shell.gluing)
            return probe.throat.stiffness(beta2_throat)

        v_ab = (throat_stiff_at(self.a0 + h)
                - throat_stiff_at(self.a0 - h)) / (2.0 * h)
        return np.array([[self.throat.stiffness(beta2_throat), v_ab],
                         [v_ab, self.shell.stiffness(beta2_shell)]])

    def stiffnesses(self, beta2_throat: float = -1.0,
                    beta2_shell: float = 0.5) -> Dict[str, object]:
        m = self.stiffness_matrix(beta2_throat, beta2_shell)
        return {"matrix": m.tolist(),
                "throat_stiffness": float(m[0, 0]),
                "shell_stiffness": float(m[1, 1]),
                "off_diagonal": float(m[0, 1]),
                "both_stiffnesses_positive": bool(m[0, 0] > 0 and m[1, 1] > 0),
                "is_diagonal": bool(abs(m[0, 1])
                                    < 1e-8 * max(abs(m[0, 0]), abs(m[1, 1]))),
                "note": ("stiffnesses, not normal-mode frequencies — no "
                         "kinetic metric has been derived")}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_junction_reproduces_known_shells() -> Dict[str, object]:
    """Two ``D = 4`` controls with published answers, before anything new."""
    rows = []
    for mass in (0.001, 0.01, 0.1):
        r = 10.0
        s = surface_stress(Region.from_mass(0.0), Region.from_mass(mass),
                           ORDINARY_BUBBLE, r)
        rows.append({"case": "bubble", "mass": mass, "sigma": s["sigma"],
                     "rest_mass": 4.0 * math.pi * r * r * s["sigma"],
                     "rest_mass_error": abs(4.0 * math.pi * r * r * s["sigma"]
                                            - mass) / mass})
    worst = 0.0
    for mass in (0.1, 0.5):
        for r in (3.0, 5.0):
            reg = Region.from_mass(mass)
            s = surface_stress(reg, reg, MINIMAL_SURFACE, r)
            ref = -math.sqrt(reg.f(r)) / (2.0 * math.pi * G * r)
            worst = max(worst, abs(s["sigma"] - ref))
            rows.append({"case": "visser_throat", "mass": mass, "radius": r,
                         "sigma": s["sigma"], "reference": ref})
    bub = [r for r in rows if r["case"] == "bubble"]
    return {"rows": rows, "worst_visser_error": worst,
            "the_bubble_is_ordinary": bool(all(r["sigma"] > 0 for r in bub)),
            "its_rest_mass_is_the_bulk_mass": bool(
                bub[0]["rest_mass_error"] < 1e-3),
            "visser_is_reproduced": bool(worst < 1e-12)}


def measure_the_gluing_fixes_the_sign(dim: int = 4) -> Dict[str, object]:
    """``ε`` derived from the branches, and which gluings force a sign.

    The point of the table: ``η = −1`` covers **two** gluings whose forced signs
    are *opposite*, so ``η`` alone decides nothing and the earlier framing in
    terms of it was too coarse.
    """
    rows = []
    minus = Region(mu=0.4, dim=dim)
    plus = Region(mu=0.8, dim=dim)
    for g in GLUINGS:
        s = surface_stress(minus, plus, g, 10.0)
        rows.append({"gluing": g.name, "minus_branch": g.minus,
                     "plus_branch": g.plus,
                     "eps": [g.eps_minus, g.eps_plus], "eta": g.eta,
                     "forced_sign": g.forced_sign, "sigma": s["sigma"],
                     "sign_agrees": (g.forced_sign is None
                                     or int(np.sign(s["sigma"]))
                                     == g.forced_sign)})
    eta_minus = [r for r in rows if r["eta"] == -1]
    return {
        "rows": rows,
        "eta_minus_one_covers": [r["gluing"] for r in eta_minus],
        "their_forced_signs": [r["forced_sign"] for r in eta_minus],
        "eta_alone_does_not_decide": bool(
            len({r["forced_sign"] for r in eta_minus}) == 2),
        "every_forced_sign_is_realised": bool(all(r["sign_agrees"]
                                                  for r in rows)),
        "the_minimal_surface_is_exotic": bool(
            [r for r in rows if r["gluing"] == "minimal surface"][0]["sigma"]
            < 0.0),
        "the_maximal_surface_is_not": bool(
            [r for r in rows if r["gluing"] == "maximal surface"][0]["sigma"]
            > 0.0),
        "and_the_maximal_one_is_closed_off": (
            "a maximal surface retains r ≤ R on both sides, so the manifold "
            "caps off either way and shares no bulk with the throat — it is "
            "non-exotic precisely because it is disconnected"),
    }


def measure_the_forced_signs_hold_in_any_dimension(
        dims: Sequence[int] = (4, 5, 6), samples: int = 40_000,
        seed: int = 11) -> Dict[str, object]:
    """The two identities, swept over random bulks in ``D = 4, 5, 6``.

    ``D = 5`` is the Tangherlini case this program actually cares about; ``D =
    4`` is kept as the regression.  The sweep checks the implementation — the
    claims are identities, since ``σ`` for a minimal surface is
    ``−(D−2)(β₊+β₋)/8πGR`` and for a maximal one the same with the other sign.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for d in dims:
        n = min_bad = max_bad = 0
        worst_min = -math.inf
        worst_max = math.inf
        per = samples // len(dims)
        for _ in range(per):
            r = float(rng.uniform(1.0, 40.0))

            def draw() -> Region:
                kind = int(rng.integers(0, 3))
                if kind == 0:
                    return Region(mu=float(rng.uniform(0, 0.4 * r ** (d - 3))),
                                  dim=d)
                if kind == 1:
                    big = float(rng.uniform(2.0 * r, 20.0 * r))
                    lam = (d - 1) * (d - 2) / (2.0 * big * big)
                    return Region(lambda_=lam, dim=d)
                return Region(mu=float(rng.uniform(0, 0.3 * r ** (d - 3))),
                              charge=float(rng.uniform(0, 0.2 * r ** (d - 3))),
                              dim=d)

            m_, p_ = draw(), draw()
            rdot = float(rng.normal(0.0, 0.5))
            if (m_.f(r) + rdot ** 2 <= 0.0) or (p_.f(r) + rdot ** 2 <= 0.0):
                continue
            n += 1
            s_min = surface_stress(m_, p_, MINIMAL_SURFACE, r, rdot)["sigma"]
            s_max = surface_stress(m_, p_, MAXIMAL_SURFACE, r, rdot)["sigma"]
            worst_min = max(worst_min, s_min)
            worst_max = min(worst_max, s_max)
            if s_min >= 0.0:
                min_bad += 1
            if s_max <= 0.0:
                max_bad += 1
        rows.append({"dim": d, "samples": n,
                     "minimal_surface_violations": min_bad,
                     "maximal_surface_violations": max_bad,
                     "worst_minimal_sigma": worst_min,
                     "worst_maximal_sigma": worst_max})
    return {
        "rows": rows,
        "dims": list(dims),
        "no_violations_in_any_dimension": bool(
            all(r["minimal_surface_violations"] == 0
                and r["maximal_surface_violations"] == 0 for r in rows)),
        "d5_is_the_bam_case": True,
        "it_is_an_identity_not_a_statistic": (
            "σ_minimal = −(D−2)(β₊+β₋)/8πGR and σ_maximal = +the same, with "
            "β± ≥ 0 for any timelike shell"),
        "scope": ("Einstein gravity, Darmois–Israel thin shells, spherical "
                  "symmetry, vacuum bulk"),
    }


def measure_the_detached_shell_can_be_ordinary(
        mu_outer: float = 2.0, screened: float = 0.6,
        radii: Sequence[float] = (8.0, 20.0, 60.0), dim: int = 4
        ) -> Dict[str, object]:
    """The same two regions, four gluings, at several radii."""
    rows = []
    for a in radii:
        m_ = Region(mu=mu_outer - screened, dim=dim)
        p_ = Region(mu=mu_outer, dim=dim)
        entry = {"radius": a}
        for g in GLUINGS:
            entry[g.name] = DetachedShell(m_, p_, a, g).stress()["sigma"]
        rows.append(entry)
    return {
        "rows": rows,
        "the_bubble_is_ordinary": bool(
            all(r["ordinary bubble"] > 0.0 for r in rows)),
        "the_minimal_surface_is_exotic": bool(
            all(r["minimal surface"] < 0.0 for r in rows)),
        "the_maximal_surface_is_ordinary": bool(
            all(r["maximal surface"] > 0.0 for r in rows)),
        "but_it_is_disconnected": True,
    }


def measure_the_shell_potential_gradient(
        mu_outer: float = 2.0, b0: float = 5.0,
        screens: Sequence[float] = (0.0, 0.2, 0.6, 1.2), dim: int = 4
        ) -> Dict[str, object]:
    """Observable 2, reported as a potential gradient rather than a force.

    Taken at fixed throat rest mass with no equation-of-state response, so it
    is the gradient of the shift that screening produces and not an
    equilibrium-consistent force.  The acceleration contribution is half it,
    from ``ḃ² + V = 0``.
    """
    rows = []
    for dm in screens:
        sysm = ThroatShellSystem(mu_outer=mu_outer, screened=dm, b0=b0,
                                 dim=dim)
        g = sysm.shell_potential_gradient()
        rows.append({"screened_mu": dm,
                     "potential_shift_at_b0": sysm.shell_potential_shift(b0),
                     "minus_dV_db": g["minus_dV_db"],
                     "acceleration_contribution":
                         g["acceleration_contribution"]})
    nz = [r for r in rows if r["screened_mu"] > 0.0]
    return {
        "rows": rows,
        "the_gradient_opposes_closure": bool(all(r["minus_dV_db"] > 0.0
                                                 for r in nz)),
        "it_grows_with_the_screened_mass": bool(
            all(b["minus_dV_db"] > a["minus_dV_db"]
                for a, b in zip(nz, nz[1:]))),
        "zero_shell_gives_zero_gradient": bool(
            abs(rows[0]["minus_dV_db"]) < 1e-12),
        "it_is_not_an_equilibrium_consistent_force": (
            "fixed throat rest mass, no equation-of-state response; the "
            "acceleration contribution is −½∂ΔV/∂b"),
    }


def measure_the_throat_and_shell_are_decoupled(
        mu_outer: float = 2.0, screened: float = 0.6, b0: float = 5.0,
        radii: Sequence[float] = (8.0, 20.0, 60.0, 200.0), dim: int = 4
        ) -> Dict[str, object]:
    """A genuinely different family of shells that the throat cannot tell apart.

    The vanishing off-diagonal is **structural** — Birkhoff is imported the
    moment the intervening region is written with a constant mass parameter —
    and establishes no separation-dependent coupling *in this model*, not that
    every spherical trapped resonator is impossible.
    """
    rows = []
    for a in radii:
        sysm = ThroatShellSystem(mu_outer=mu_outer, screened=screened,
                                 b0=b0, a0=a, dim=dim)
        st = sysm.shell.static()
        rows.append({"radius": a, "shell_sigma": st["sigma"],
                     "shell_pressure": st["p"],
                     "throat_sigma": sysm.throat.static()["sigma"],
                     "off_diagonal": sysm.stiffnesses()["off_diagonal"]})
    sig = [r["shell_sigma"] for r in rows]
    th = [r["throat_sigma"] for r in rows]
    return {
        "rows": rows,
        "shell_sigma_varies_by": max(sig) / min(sig),
        "throat_sigma_spread": max(th) - min(th),
        "worst_off_diagonal": max(abs(r["off_diagonal"]) for r in rows),
        "the_shells_are_genuinely_different": bool(max(sig) / min(sig) > 2.0),
        "the_throat_never_notices": bool(max(th) - min(th) == 0.0),
        "the_off_diagonal_vanishes": bool(
            max(abs(r["off_diagonal"]) for r in rows) < 1e-12),
        "but_that_is_structural_not_measured": (
            "Birkhoff is assumed when the region between is written with a "
            "constant mass parameter; the zero confirms the implementation is "
            "consistent with it and proves nothing more"),
        "what_it_establishes": (
            "no separation-dependent coupling in this vacuum spherical model "
            "— not that every spherical trapped resonator is impossible"),
    }


def measure_the_stability_window(
        b0: float = 5.0, mus: Sequence[float] = (2.0, 1.8, 1.6, 1.4, 1.0),
        dim: int = 4) -> Dict[str, object]:
    """Observable 3: the stiffness, and what screening does to its window.

    ``β²`` is an equation-of-state derivative parameter at the equilibrium.  No
    sound-speed reading is attached to its sign, which for exotic matter would
    not be meaningful.
    """
    rows = []
    for mu in mus:
        t = Z2Throat(mu=mu, b0=b0, dim=dim)
        w, s = t.stability_window(), t.static()
        rows.append({"mu_interior": mu, "sigma": s["sigma"],
                     "is_exotic": s["is_exotic"],
                     "beta2_critical": w["beta2_critical"],
                     "stable_below_critical": w["stable_below_critical"]})
    crits = [r["beta2_critical"] for r in rows]
    return {
        "rows": rows, "throat_radius": b0, "dim": dim,
        "screening_raises_the_critical_beta2": bool(
            all(b > a for a, b in zip(crits, crits[1:]))),
        "the_window_always_needs_negative_beta2": bool(all(c < 0.0
                                                           for c in crits)),
        "the_throat_is_always_exotic": bool(all(r["is_exotic"] for r in rows)),
        "beta2_is_an_eos_derivative_not_a_sound_speed": True,
    }


def measure_the_three_observables_are_independent(
        b0: float = 5.0, mu_outer: float = 2.0, screened: float = 0.6,
        dim: int = 4) -> Dict[str, object]:
    """One system, three questions, three different signs."""
    sysm = ThroatShellSystem(mu_outer=mu_outer, screened=screened, b0=b0,
                             dim=dim)
    shell = sysm.shell.stress()
    throat = sysm.throat.static()
    window = sysm.throat.stability_window()
    grad = sysm.shell_potential_gradient()
    return {
        "dim": dim,
        "observable_1_shell_sigma": shell["sigma"],
        "observable_1_shell_is_exotic": shell["is_exotic"],
        "observable_2_minus_dV_db": grad["minus_dV_db"],
        "observable_2_opposes_closure": grad["opposes_closure"],
        "observable_3_beta2_critical": window["beta2_critical"],
        "observable_3_needs_negative_beta2": window[
            "window_needs_negative_beta2"],
        "throat_sigma": throat["sigma"],
        "the_throat_is_still_exotic": throat["is_exotic"],
        "they_do_not_agree": bool((not shell["is_exotic"])
                                  and grad["opposes_closure"]
                                  and throat["is_exotic"]),
    }


def measure_the_stiffnesses_scale_dimensionally(
        b0: float = 5.0, a0: float = 20.0, beta2: float = -2.0,
        scales: Sequence[float] = (1.0, 2.0, 4.0), dim: int = 4
        ) -> Dict[str, object]:
    """A dimensional check on the stiffnesses, and nothing stronger.

    Rescaling every length and mass parameter together must send the
    stiffnesses as ``1/L²``.  This is a units check on the implementation; it
    is **not** a demonstration that a fixed system has no dilation mode, which
    would need the kinetic term this module does not derive.
    """
    rows = []
    for s in scales:
        p = float(dim - 3)
        sysm = ThroatShellSystem(mu_outer=2.0 * s ** p,
                                 screened=0.6 * s ** p,
                                 b0=b0 * s, a0=a0 * s, dim=dim)
        st = sysm.stiffnesses(beta2, 0.5)
        vals = [st["throat_stiffness"], st["shell_stiffness"]]
        rows.append({"scale": s, "stiffnesses": vals,
                     "scaled": [v * s * s for v in vals]})
    ref = rows[0]["scaled"]
    drift = max(max(abs(r["scaled"][i] - ref[i]) / max(abs(ref[i]), 1e-30)
                    for i in range(2)) for r in rows)
    return {
        "rows": rows, "dim": dim, "worst_scaling_drift": drift,
        "stiffnesses_scale_as_inverse_length_squared": bool(drift < 1e-6),
        "this_is_a_units_check_only": (
            "it does not show a fixed system has no dilation mode — that "
            "needs a kinetic term, which is not derived here"),
    }
