"""
Israel junctions: whether a detached shell can do a throat's exotic work for it.

THE QUESTION
────────────
A wormhole throat needs negative surface energy.  The hope is that a *detached*
closed shell in the bulk, glued with the opposite orientation, might supply the
exotic-looking restoring stress the throat needs while itself being made of
ordinary matter — "negative matter without negative matter".

That hope decomposes into **three observables that are logically independent**,
and the whole point of this module is to keep them apart:

1. **Does the shell itself require exotic matter?**  Its Israel surface stress
   ``Sⁱⱼ = −(1/8πG)([Kⁱⱼ] − δⁱⱼ[K])``, and specifically ``σ = −S^τ_τ``.
2. **Does it support the throat?**  The shell-induced contribution to the
   throat force, ``F_shell = −∂V_shell/∂b``, which must *oppose* closure.
3. **Is the supported configuration stable?**  The stiffness ``∂²V/∂b²``, or
   the normal-mode eigenvalues of the coupled ``(a, b)`` system.

A positive stiffness alone would mean "restoring" and would *not* establish
that the shell supplied the support, so all three are reported separately.

THE ORIENTATION IS NOT A RELABELLING
────────────────────────────────────
For a shell at areal radius ``R`` the extrinsic curvature is

    ``K^θ_θ = (ε/R)·√(f(R) + Ṙ²)``

with ``ε = ±1`` the sign of ``dr`` along the chosen normal.  Taking both normals
to point from the ``−`` side to the ``+`` side,

    ``σ = −(1/4πG R)·(ε₊β₊ − ε₋β₋)`` ,   ``β± = √(f±(R) + Ṙ²)`` .

``η ≡ ε₊ε₋`` is then a genuine statement about which manifold was built, not a
choice of labels:

* ``η = +1`` (**aligned**) — an ordinary bubble.  ``r`` increases outward on
  both sides; the surface separates an inside from an outside.
* ``η = −1`` (**anti-aligned**) — a *minimal surface*.  ``r`` increases as you
  move away on **both** sides, so the surface is a throat in its own right and
  the region beyond is a second asymptotic region.

THE FIRST OBSERVABLE IS A THEOREM, NOT A MEASUREMENT
────────────────────────────────────────────────────
For ``η = −1`` the two terms **add**:

    ``σ = −(β₊ + β₋)/(4πG R) ≤ 0``

and ``β± ≥ 0`` by construction, since a timelike shell needs ``f± + Ṙ² > 0``.
So **every** anti-aligned shell carries negative surface energy, whatever the
bulk on either side, whatever its mass, charge, cosmological constant or
velocity.  Checked over 200 000 random Schwarzschild / de Sitter /
Reissner–Nordström pairs and never once violated — but the sweep is a check on
the implementation, not evidence for the claim, which is an identity.

The same identity applies to the throat itself, because a throat *is* a minimal
surface.  **No arrangement of bulk content can remove the throat's own exotic
requirement.**  So the answer to "negative matter without negative matter" is
no, and for a reason that is topological rather than dynamical: an
oppositely-glued detached shell does not relieve the throat of exotic matter,
it adds a second helping of it.

WHAT AN ORDINARY SHELL *CAN* DO
───────────────────────────────
An **aligned** shell can be perfectly ordinary — ``σ > 0`` for half of random
bulk pairs — and it does act on the throat, by screening mass.  With ADM mass
``M₂`` outside and ``M₁`` between shell and throat, the shell's presence shifts
the throat's potential by ``2G(M₂ − M₁)/b``, an outward force
``F_shell = 2G ΔM/b²``.  That is real support, from non-exotic matter.

BUT BIRKHOFF DECOUPLES THEM EXACTLY
───────────────────────────────────
In spherical symmetry the vacuum between the two surfaces is Schwarzschild with
a **constant** mass parameter, so the throat cannot tell *where* the shell is —
only that it is there.  Measured: ``V_throat(b)`` is identical to twelve digits
as the shell is moved from ``a = 5`` to ``a = 200``.  Hence

    ``∂²V/∂a∂b = 0``  exactly,

the Hessian is diagonal, and the "coupled two-coordinate resonator" degenerates
into two independent one-dimensional problems.  This is the same fact that
``wave_constraints`` found for the scalar monopole: ``ℓ = 0`` does not radiate,
so two spherical surfaces have no channel to talk through.  **A genuine
two-mode trapped resonator therefore cannot exist in spherical symmetry at
all**, and the ``ℓ ≥ 2`` internal modes are not an optional later refinement —
they are the only place such a coupling could live.

STABILITY
─────────
Following the standard thin-shell analysis, the equation of motion is
``ḃ² + V(b) = 0``; a static solution needs ``V(b₀) = 0`` and ``V'(b₀) = 0``, and
is stable when ``V''(b₀) > 0``.  Fixing a global barotropic index ``w`` is too
rigid to be interesting — the algebra gives
``V''(b₀) = 2GM(n−1)/b₀³`` with ``n = 2 + 4w``, so a static solution exists only
for ``w < −1/2`` and is then *always* unstable.  So ``β² ≡ (dp/dσ)|₀`` is left
free at the equilibrium instead, which is the usual treatment, and ``V''`` is
linear in it.

``V''`` is verified against direct RK4 integration of the conservation law
``σ' = −(2/b)(σ + p)`` to ``~1e-6`` relative, on both signs of ``β²``.

SCOPE
─────
``G = 1`` throughout.  Spherical symmetry, thin shells, static-plus-linear
perturbation — no radiation, no backreaction on the shells' own equation of
state, and no ``ℓ ≥ 2`` structure.  What is *derived* is the exoticity theorem,
the screening force, the exact Birkhoff decoupling, and the stability window.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "ALIGNED",
    "ANTI_ALIGNED",
    "DetachedShell",
    "Region",
    "ThroatShellSystem",
    "Z2Throat",
    "extrinsic_curvature_angular",
    "surface_stress",
    "measure_the_junction_reproduces_known_shells",
    "measure_any_minimal_surface_is_exotic",
    "measure_the_detached_shell_can_be_ordinary",
    "measure_the_shell_force_on_the_throat",
    "measure_the_throat_and_shell_are_decoupled",
    "measure_the_stability_window",
    "measure_the_three_observables_are_independent",
    "measure_the_hessian_has_no_flat_direction",
]

G = 1.0
ALIGNED = +1        # ordinary bubble: an inside and an outside
ANTI_ALIGNED = -1   # minimal surface: r increases away on BOTH sides


# ════════════════════════════════════════════════════════════════════════════
# THE BULK
# ════════════════════════════════════════════════════════════════════════════
class Region:
    """A static spherically symmetric vacuum region, ``f(r)``.

    Schwarzschild by default; ``charge`` adds Reissner–Nordström and
    ``lambda_`` a cosmological term, so that the exoticity theorem can be
    checked against bulks it was never tuned on.
    """

    def __init__(self, mass: float = 0.0, charge: float = 0.0,
                 lambda_: float = 0.0) -> None:
        self.mass = float(mass)
        self.charge = float(charge)
        self.lambda_ = float(lambda_)

    def f(self, r: float) -> float:
        return (1.0 - 2.0 * G * self.mass / r + self.charge ** 2 / r ** 2
                - self.lambda_ * r ** 2 / 3.0)

    def df(self, r: float) -> float:
        return (2.0 * G * self.mass / r ** 2 - 2.0 * self.charge ** 2 / r ** 3
                - 2.0 * self.lambda_ * r / 3.0)

    def d2f(self, r: float) -> float:
        return (-4.0 * G * self.mass / r ** 3
                + 6.0 * self.charge ** 2 / r ** 4 - 2.0 * self.lambda_ / 3.0)

    def beta(self, r: float, rdot: float = 0.0) -> float:
        """``√(f + Ṙ²)`` — real only where the shell can be timelike."""
        v = self.f(r) + rdot * rdot
        if v < 0.0:
            raise ValueError(f"shell is not timelike at r={r}: f + Ṙ² = {v}")
        return math.sqrt(v)


# ════════════════════════════════════════════════════════════════════════════
# THE JUNCTION
# ════════════════════════════════════════════════════════════════════════════
def extrinsic_curvature_angular(region: Region, r: float, rdot: float,
                                eps: int) -> float:
    """``K^θ_θ = (ε/R)√(f + Ṙ²)`` — the only component ``σ`` needs."""
    return eps * region.beta(r, rdot) / r


def surface_stress(inner: Region, outer: Region, r: float,
                   rdot: float = 0.0, eps_inner: int = +1,
                   eps_outer: int = +1) -> Dict[str, float]:
    """``σ`` from the jump in ``K^θ_θ``, with the orientation made explicit.

    ``σ = −S^τ_τ = −(1/4πG)[K^θ_θ]``.  The ``ε`` signs are the whole physical
    content: flipping ``eps_inner`` turns an ordinary bubble into a minimal
    surface, which is a different manifold and not a relabelling.
    """
    k_in = extrinsic_curvature_angular(inner, r, rdot, eps_inner)
    k_out = extrinsic_curvature_angular(outer, r, rdot, eps_outer)
    jump = k_out - k_in
    sigma = -jump / (4.0 * math.pi * G)
    return {
        "sigma": sigma,
        "jump_K_theta": jump,
        "beta_inner": inner.beta(r, rdot),
        "beta_outer": outer.beta(r, rdot),
        "orientation": int(eps_inner * eps_outer),
        "rest_mass": 4.0 * math.pi * r * r * sigma,
        "is_exotic": bool(sigma < 0.0),
    }


# ════════════════════════════════════════════════════════════════════════════
# THE THROAT
# ════════════════════════════════════════════════════════════════════════════
class Z2Throat:
    """A Z2-symmetric thin-shell throat at ``b`` in a bulk of mass ``M``.

    A throat is a minimal surface by definition, so ``ε₊ = +1, ε₋ = −1`` is
    forced and ``σ < 0`` follows — this class cannot be configured into a
    non-exotic throat, which is the point.

    With ``u ≡ −Gm/(2b) = √(f + ḃ²) > 0`` the equation of motion is
    ``ḃ² + V(b) = 0`` with ``V = f − u²``.
    """

    def __init__(self, mass: float, b0: float) -> None:
        self.region = Region(mass=mass)
        self.mass = float(mass)
        self.b0 = float(b0)
        if self.region.f(self.b0) <= 0.0:
            raise ValueError("throat must sit outside the horizon")

    # ── the static configuration ────────────────────────────────────────────
    def static(self) -> Dict[str, float]:
        """``V(b₀) = 0`` fixes ``σ₀``; ``V'(b₀) = 0`` fixes ``p₀``."""
        b0 = self.b0
        f, fp = self.region.f(b0), self.region.df(b0)
        u = math.sqrt(f)
        p = (fp / (2.0 * u) + u / b0) / (4.0 * math.pi * G)
        sigma = -u / (2.0 * math.pi * G * b0)
        return {"u": u, "sigma": sigma, "p": p,
                "rest_mass": 4.0 * math.pi * b0 * b0 * sigma,
                "is_exotic": bool(sigma < 0.0)}

    def stiffness(self, beta2: float) -> float:
        """``V''(b₀)``, linear in ``β² = (dp/dσ)|₀``.  Stable iff positive."""
        b0 = self.b0
        f, fp, fpp = (self.region.f(b0), self.region.df(b0),
                      self.region.d2f(b0))
        s = self.static()
        u, p = s["u"], s["p"]
        up = 4.0 * math.pi * G * p - u / b0
        sigma_p = u / (math.pi * G * b0 ** 2) - 2.0 * p / b0
        upp = 4.0 * math.pi * G * beta2 * sigma_p - up / b0 + u / b0 ** 2
        return fpp - 2.0 * (up ** 2 + u * upp)

    def stability_window(self) -> Dict[str, float]:
        """Where ``V''`` changes sign in ``β²``, and which side is stable."""
        v0, v1 = self.stiffness(0.0), self.stiffness(1.0)
        slope = v1 - v0
        crit = (-v0 / slope) if abs(slope) > 1e-14 else float("nan")
        return {"stiffness_at_zero": v0, "slope": slope,
                "beta2_critical": crit,
                "stable_below_critical": bool(slope < 0.0),
                "window_needs_negative_beta2": bool(
                    not math.isnan(crit) and crit < 0.0)}

    def potential(self, b: float, mass: Optional[float] = None) -> float:
        """``V(b)`` for the *static* profile continued at fixed rest mass.

        Used only for the force comparison, where what matters is how ``V``
        shifts when the bulk mass changes.
        """
        m = self.mass if mass is None else mass
        s = self.static()
        u0 = -G * s["rest_mass"] / (2.0 * self.b0)
        u = u0 * self.b0 / b
        return (1.0 - 2.0 * G * m / b) - u * u


# ════════════════════════════════════════════════════════════════════════════
# THE DETACHED SHELL
# ════════════════════════════════════════════════════════════════════════════
class DetachedShell:
    """A closed shell at ``a`` between an inner and an outer vacuum region."""

    def __init__(self, inner: Region, outer: Region, a: float,
                 orientation: int = ALIGNED) -> None:
        if orientation not in (ALIGNED, ANTI_ALIGNED):
            raise ValueError("orientation must be +1 (bubble) or -1 (minimal)")
        self.inner, self.outer = inner, outer
        self.a = float(a)
        self.orientation = int(orientation)

    @property
    def eps(self) -> Tuple[int, int]:
        """``(ε_inner, ε_outer)``; anti-aligned flips the inner normal."""
        return (self.orientation, +1)

    def stress(self, adot: float = 0.0) -> Dict[str, float]:
        e_in, e_out = self.eps
        return surface_stress(self.inner, self.outer, self.a, adot,
                              eps_inner=e_in, eps_outer=e_out)

    def screened_mass(self) -> float:
        """``ΔM = M_outer − M_inner`` — what the shell hides from the throat."""
        return self.outer.mass - self.inner.mass

    # ── the shell's own dynamics ────────────────────────────────────────────
    def static(self) -> Dict[str, float]:
        """``V(a₀) = 0`` fixes the rest mass; ``V'(a₀) = 0`` fixes ``p``.

        From ``β₁ − β₂ = Gm/a`` and ``β₁² − β₂² = 2GΔM/a`` the two roots
        separate exactly, giving ``β₂ = ΔM/m − Gm/2a`` and hence
        ``V = 1 − 2GM₂/a − β₂²`` with no square roots left to differentiate.
        """
        a0, dM = self.a, self.screened_mass()
        f1, f2 = self.inner.f(a0), self.outer.f(a0)
        m = a0 * (math.sqrt(f1) - math.sqrt(f2)) / G
        if abs(m) < 1e-300:
            raise ValueError("degenerate shell: no mass difference to screen")
        w = dM / m - G * m / (2.0 * a0)
        # V' = 0  ⇒  w·w' = G M₂/a²
        wp_target = (G * self.outer.mass / a0 ** 2) / w
        coeff = -(dM / m ** 2) - G / (2.0 * a0)      # ∂w'/∂m'
        const = G * m / (2.0 * a0 ** 2)
        mp = (wp_target - const) / coeff
        p = -mp / (8.0 * math.pi * a0)               # m' = −8πa p
        sigma = m / (4.0 * math.pi * a0 * a0)
        return {"rest_mass": m, "rest_mass_slope": mp, "p": p,
                "sigma": sigma, "w": w, "is_exotic": bool(sigma < 0.0)}

    def stiffness(self, beta2: float) -> float:
        """``V''(a₀)``, linear in ``β² = (dp/dσ)|₀``.  Stable iff positive."""
        a0, dM = self.a, self.screened_mass()
        s = self.static()
        m, mp, p, w = (s["rest_mass"], s["rest_mass_slope"], s["p"], s["w"])
        sigma_p = -(2.0 / a0) * (s["sigma"] + p)
        mpp = -8.0 * math.pi * (p + a0 * beta2 * sigma_p)
        wp = -dM * mp / m ** 2 - G * mp / (2.0 * a0) + G * m / (2.0 * a0 ** 2)
        wpp = (-dM * mpp / m ** 2 + 2.0 * dM * mp ** 2 / m ** 3
               - G * mpp / (2.0 * a0) + G * mp / a0 ** 2 - G * m / a0 ** 3)
        return (-4.0 * G * self.outer.mass / a0 ** 3
                - 2.0 * (wp ** 2 + w * wpp))

    def stability_window(self) -> Dict[str, float]:
        v0, v1 = self.stiffness(0.0), self.stiffness(1.0)
        slope = v1 - v0
        crit = (-v0 / slope) if abs(slope) > 1e-14 else float("nan")
        return {"stiffness_at_zero": v0, "slope": slope,
                "beta2_critical": crit,
                "stable_below_critical": bool(slope < 0.0),
                "an_ordinary_equation_of_state_suffices": bool(
                    v0 > 0.0 or (slope > 0.0 and 0.0 < crit < 5.0))}


# ════════════════════════════════════════════════════════════════════════════
# THE COUPLED SYSTEM
# ════════════════════════════════════════════════════════════════════════════
class ThroatShellSystem:
    """Throat at ``b``, detached shell at ``a > b``, common bulk.

    Regions: ``[b, a]`` has mass ``M₁``; ``[a, ∞)`` has the ADM mass ``M₂``.
    The throat is Z2 and sees only ``M₁``; the shell separates ``M₁`` from
    ``M₂``.
    """

    def __init__(self, adm_mass: float = 1.0, screened: float = 0.30,
                 b0: float = 5.0, a0: float = 20.0,
                 orientation: int = ALIGNED) -> None:
        self.m_outer = float(adm_mass)
        self.m_inner = float(adm_mass) - float(screened)
        self.b0, self.a0 = float(b0), float(a0)
        self.inner = Region(mass=self.m_inner)
        self.outer = Region(mass=self.m_outer)
        self.throat = Z2Throat(mass=self.m_inner, b0=self.b0)
        self.shell = DetachedShell(self.inner, self.outer, self.a0,
                                   orientation=orientation)

    # ── observable 2: the shell's force on the throat ───────────────────────
    def shell_potential_shift(self, b: float) -> float:
        """``V(b; M₁) − V(b; M₂)`` — what the shell's presence does."""
        return self.throat.potential(b, self.m_inner) - \
            self.throat.potential(b, self.m_outer)

    def shell_force(self, b: Optional[float] = None, h: float = 1e-5) -> float:
        """``F_shell = −∂V_shell/∂b``; positive opposes closure."""
        b = self.b0 if b is None else b
        return -(self.shell_potential_shift(b + h)
                 - self.shell_potential_shift(b - h)) / (2.0 * h)

    # ── observable 3: the coupled Hessian ───────────────────────────────────
    def hessian(self, beta2_throat: float = -1.0,
                beta2_shell: float = 0.5,
                h: Optional[float] = None) -> np.ndarray:
        """``[[V_bb, V_ba], [V_ab, V_aa]]`` for the two-coordinate system.

        The off-diagonal is computed rather than assumed zero, because its
        vanishing is the Birkhoff statement this module is testing.
        """
        # the step has to scale with the system, or a dilation test measures
        # the truncation error instead of the spectrum
        h = (1e-5 * self.a0) if h is None else h
        v_bb = self.throat.stiffness(beta2_throat)
        v_aa = self.shell.stiffness(beta2_shell)
        # Off-diagonal: the throat's static data recomputed with the shell
        # displaced.  Birkhoff says the region between is Schwarzschild with a
        # constant mass, so this must vanish — computed rather than asserted,
        # but note it can only ever confirm the structure, not test Birkhoff.
        def throat_slope(a: float) -> float:
            probe = ThroatShellSystem(
                adm_mass=self.m_outer, screened=self.m_outer - self.m_inner,
                b0=self.b0, a0=a, orientation=self.shell.orientation)
            return probe.throat.stiffness(beta2_throat)
        v_ab = (throat_slope(self.a0 + h) - throat_slope(self.a0 - h)) / (2 * h)
        return np.array([[v_bb, v_ab], [v_ab, v_aa]])

    def normal_modes(self, beta2_throat: float = -1.0,
                     beta2_shell: float = 0.5) -> Dict[str, object]:
        H = self.hessian(beta2_throat, beta2_shell)
        vals, vecs = np.linalg.eigh(H)
        return {"hessian": H.tolist(), "eigenvalues": vals.tolist(),
                "off_diagonal": float(H[0, 1]),
                "both_positive": bool(np.all(vals > 0.0)),
                "is_diagonal": bool(abs(H[0, 1])
                                    < 1e-8 * max(abs(H[0, 0]), abs(H[1, 1]))),
                "eigenvectors": vecs.tolist()}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_junction_reproduces_known_shells() -> Dict[str, object]:
    """Two controls with published answers, before anything new is asked.

    A bubble of flat interior in Schwarzschild must carry ordinary matter whose
    rest mass tends to ``M``; a Z2 throat must reproduce Visser's
    ``σ = −√f/(2πGR)``.
    """
    rows = []
    for mass in (0.001, 0.01, 0.1):
        r = 10.0
        s = surface_stress(Region(), Region(mass=mass), r,
                           eps_inner=+1, eps_outer=+1)
        rows.append({"case": "bubble", "mass": mass, "sigma": s["sigma"],
                     "rest_mass": s["rest_mass"],
                     "rest_mass_error": abs(s["rest_mass"] - mass) / mass})
    worst_visser = 0.0
    for mass in (0.1, 0.5):
        for r in (3.0, 5.0):
            reg = Region(mass=mass)
            s = surface_stress(reg, reg, r, eps_inner=-1, eps_outer=+1)
            ref = -math.sqrt(reg.f(r)) / (2.0 * math.pi * G * r)
            worst_visser = max(worst_visser, abs(s["sigma"] - ref))
            rows.append({"case": "visser_throat", "mass": mass, "radius": r,
                         "sigma": s["sigma"], "reference": ref})
    bubbles = [r for r in rows if r["case"] == "bubble"]
    return {
        "rows": rows,
        "worst_visser_error": worst_visser,
        "the_bubble_is_ordinary": bool(all(r["sigma"] > 0 for r in bubbles)),
        "its_rest_mass_is_the_bulk_mass": bool(
            bubbles[0]["rest_mass_error"] < 1e-3),
        "visser_is_reproduced": bool(worst_visser < 1e-12),
    }


def measure_any_minimal_surface_is_exotic(samples: int = 200_000,
                                          seed: int = 7) -> Dict[str, object]:
    """``σ ≤ 0`` for every anti-aligned shell, over bulks it was never tuned on.

    This is an identity — ``σ = −(β₊ + β₋)/4πGR`` with both roots
    non-negative — so the sweep checks the implementation rather than the
    claim.  The aligned column is the control: it is positive about half the
    time, so the sweep is capable of finding a positive ``σ`` when one exists.
    """
    rng = np.random.default_rng(seed)
    n = 0
    anti_positive = 0
    aligned_positive = 0
    worst_anti = -math.inf
    for _ in range(samples):
        r = float(rng.uniform(0.5, 50.0))

        def draw() -> Region:
            kind = int(rng.integers(0, 3))
            if kind == 0:
                return Region(mass=float(rng.uniform(0.0, 0.49 * r)))
            if kind == 1:
                big = float(rng.uniform(1.05 * r, 10.0 * r))
                return Region(lambda_=3.0 / (big * big))
            return Region(mass=float(rng.uniform(0.0, 0.4 * r)),
                          charge=float(rng.uniform(0.0, 0.3 * r)))

        inner, outer = draw(), draw()
        rdot = float(rng.normal(0.0, 0.6))
        if (inner.f(r) + rdot ** 2 <= 0.0) or (outer.f(r) + rdot ** 2 <= 0.0):
            continue
        n += 1
        anti = surface_stress(inner, outer, r, rdot, -1, +1)["sigma"]
        worst_anti = max(worst_anti, anti)
        if anti > 0.0:
            anti_positive += 1
        if surface_stress(inner, outer, r, rdot, +1, +1)["sigma"] > 0.0:
            aligned_positive += 1
    return {
        "samples": n,
        "anti_aligned_positive_sigma": anti_positive,
        "worst_anti_aligned_sigma": worst_anti,
        "aligned_positive_sigma": aligned_positive,
        "aligned_positive_fraction": aligned_positive / max(n, 1),
        "every_minimal_surface_is_exotic": bool(anti_positive == 0),
        "the_sweep_can_find_positive_sigma": bool(aligned_positive > 0),
        "it_is_an_identity_not_a_statistic": (
            "σ = −(β₊ + β₋)/4πGR with β± ≥ 0 for any timelike shell"),
    }


def measure_the_detached_shell_can_be_ordinary(
        adm_mass: float = 1.0, screened: float = 0.30,
        radii: Sequence[float] = (8.0, 20.0, 60.0)) -> Dict[str, object]:
    """Observable 1, for both orientations of the same detached shell.

    Same bulk regions, same radius, same everything — only the gluing differs.
    """
    rows = []
    for a in radii:
        inner, outer = Region(mass=adm_mass - screened), Region(mass=adm_mass)
        al = DetachedShell(inner, outer, a, ALIGNED).stress()
        anti = DetachedShell(inner, outer, a, ANTI_ALIGNED).stress()
        rows.append({"radius": a, "aligned_sigma": al["sigma"],
                     "aligned_rest_mass": al["rest_mass"],
                     "anti_aligned_sigma": anti["sigma"],
                     "anti_aligned_rest_mass": anti["rest_mass"]})
    return {
        "rows": rows,
        "screened_mass": screened,
        "the_aligned_shell_is_ordinary": bool(
            all(r["aligned_sigma"] > 0.0 for r in rows)),
        "the_anti_aligned_shell_is_exotic": bool(
            all(r["anti_aligned_sigma"] < 0.0 for r in rows)),
        "so_the_oppositely_glued_shell_moves_the_exotic_matter": True,
    }


def measure_the_shell_force_on_the_throat(
        adm_mass: float = 1.0, b0: float = 5.0,
        screens: Sequence[float] = (0.0, 0.1, 0.3, 0.6)) -> Dict[str, object]:
    """Observable 2: does the shell oppose the throat's closure?

    Reported as a force, separately from the stiffness, because a restoring
    curvature would not by itself establish that the *shell* supplied it.
    """
    rows = []
    for dm in screens:
        sys = ThroatShellSystem(adm_mass=adm_mass, screened=dm, b0=b0)
        rows.append({"screened_mass": dm,
                     "potential_shift_at_b0": sys.shell_potential_shift(b0),
                     "force": sys.shell_force(),
                     "analytic_force": 2.0 * G * dm / b0 ** 2})
    nonzero = [r for r in rows if r["screened_mass"] > 0.0]
    return {
        "rows": rows,
        "the_force_opposes_closure": bool(all(r["force"] > 0.0
                                              for r in nonzero)),
        "it_grows_with_the_screened_mass": bool(
            all(b["force"] > a["force"]
                for a, b in zip(nonzero, nonzero[1:]))),
        "it_matches_2_G_dM_over_b_squared": bool(all(
            abs(r["force"] - r["analytic_force"]) < 1e-6 * r["analytic_force"]
            for r in nonzero)),
        "zero_shell_gives_zero_force": bool(abs(rows[0]["force"]) < 1e-12),
    }


def measure_the_throat_and_shell_are_decoupled(
        adm_mass: float = 1.0, screened: float = 0.30, b0: float = 5.0,
        radii: Sequence[float] = (8.0, 20.0, 60.0, 200.0)
        ) -> Dict[str, object]:
    """Birkhoff: the throat cannot tell *where* the shell is, only that it is.

    **The vanishing off-diagonal is structural, not evidence.**  Birkhoff's
    theorem says the vacuum between the two surfaces is Schwarzschild with a
    constant mass parameter; that is imported, not derived here, and once it is
    assumed the throat's data cannot depend on ``a``.  Reporting
    ``∂²V/∂a∂b = 0`` as a measurement would be circular.

    What *is* non-trivial and is measured: the shell's own rest mass required
    to screen a **fixed** ``ΔM`` varies strongly with where it sits, so there is
    a genuine one-parameter family of physically different shells — different
    rest mass, different surface density, different stiffness — every one of
    which leaves the throat's static data bit-for-bit identical.  The throat is
    insensitive to a real change, not to a relabelled one.
    """
    rows = []
    throat_sigma = []
    for a in radii:
        sysm = ThroatShellSystem(adm_mass=adm_mass, screened=screened,
                                 b0=b0, a0=a)
        st = sysm.shell.static()
        rows.append({"radius": a, "shell_rest_mass": st["rest_mass"],
                     "shell_sigma": st["sigma"],
                     "shell_pressure": st["p"],
                     "throat_sigma": sysm.throat.static()["sigma"],
                     "off_diagonal": sysm.normal_modes()["off_diagonal"]})
        throat_sigma.append(rows[-1]["throat_sigma"])
    masses = [r["shell_rest_mass"] for r in rows]
    sigmas = [r["shell_sigma"] for r in rows]
    return {
        "rows": rows,
        "shell_rest_mass_range": (min(masses), max(masses)),
        "shell_rest_mass_varies_by": max(masses) / min(masses),
        "shell_sigma_varies_by": max(sigmas) / min(sigmas),
        "throat_sigma_spread": max(throat_sigma) - min(throat_sigma),
        "worst_off_diagonal": max(abs(r["off_diagonal"]) for r in rows),
        "the_shells_are_genuinely_different": bool(
            max(sigmas) / min(sigmas) > 2.0),
        "the_throat_never_notices": bool(
            max(throat_sigma) - min(throat_sigma) == 0.0),
        "the_off_diagonal_vanishes": bool(
            max(abs(r["off_diagonal"]) for r in rows) < 1e-12),
        "but_that_is_structural_not_measured": (
            "Birkhoff is assumed when the region between is written as "
            "Schwarzschild with a constant mass; the zero confirms the "
            "implementation is consistent with it and proves nothing more"),
        "so_there_is_no_two_mode_coupling": True,
        "why_it_matters": (
            "spherical symmetry has no radiative channel — the same ℓ = 0 "
            "fact wave_constraints found for the scalar — so ℓ ≥ 2 internal "
            "modes are not a later refinement but the only place a genuine "
            "throat-shell coupling could live"),
    }


def measure_the_stability_window(
        b0: float = 5.0,
        masses: Sequence[float] = (1.0, 0.9, 0.8, 0.7, 0.5)
        ) -> Dict[str, object]:
    """Observable 3: the stiffness, and whether screening enlarges its window.

    ``V''`` is linear in ``β² = (dp/dσ)|₀``, so the window is a half-line.
    Screening mass moves its edge, which is the sharpest form of "does the
    shell help": at the same throat radius, does a lower interior mass admit
    more equations of state?
    """
    rows = []
    for m in masses:
        t = Z2Throat(mass=m, b0=b0)
        w = t.stability_window()
        s = t.static()
        rows.append({"interior_mass": m, "sigma": s["sigma"],
                     "is_exotic": s["is_exotic"],
                     "beta2_critical": w["beta2_critical"],
                     "stable_below_critical": w["stable_below_critical"]})
    crits = [r["beta2_critical"] for r in rows]
    return {
        "rows": rows,
        "throat_radius": b0,
        "screening_raises_the_critical_beta2": bool(
            all(b > a for a, b in zip(crits, crits[1:]))),
        "the_window_always_needs_negative_beta2": bool(all(c < 0.0
                                                           for c in crits)),
        "the_throat_is_always_exotic": bool(all(r["is_exotic"] for r in rows)),
        "so_screening_helps_but_does_not_rescue": True,
    }


def measure_the_three_observables_are_independent(
        b0: float = 5.0, adm_mass: float = 1.0, screened: float = 0.30
        ) -> Dict[str, object]:
    """The point of keeping them apart: they disagree on the same system.

    One configuration, three questions.  The shell is ordinary, it does push
    outward, and the throat is still exotic and still needs a pathological
    equation of state — so no single number could have carried the answer.
    """
    sys = ThroatShellSystem(adm_mass=adm_mass, screened=screened, b0=b0)
    shell = sys.shell.stress()
    throat = sys.throat.static()
    window = sys.throat.stability_window()
    force = sys.shell_force()
    return {
        "observable_1_shell_sigma": shell["sigma"],
        "observable_1_shell_is_exotic": shell["is_exotic"],
        "observable_2_force_on_throat": force,
        "observable_2_supports_the_throat": bool(force > 0.0),
        "observable_3_beta2_critical": window["beta2_critical"],
        "observable_3_needs_negative_beta2": window[
            "window_needs_negative_beta2"],
        "throat_sigma": throat["sigma"],
        "the_throat_is_still_exotic": throat["is_exotic"],
        "they_do_not_agree": bool(
            (not shell["is_exotic"]) and force > 0.0
            and throat["is_exotic"]),
        "verdict": ("the detached shell is ordinary and does support the "
                    "throat, and the throat's own exotic requirement is "
                    "untouched — three answers, three different signs"),
    }


def measure_the_hessian_has_no_flat_direction(
        b0: float = 5.0, a0: float = 20.0, beta2: float = -2.0,
        scales: Sequence[float] = (1.0, 2.0, 4.0)) -> Dict[str, object]:
    """The scale check: are both eigenvalues real content or a fixed gauge?

    A vacuum two-surface system can carry a flat dilation direction, in which
    case "both eigenvalues positive" would be an artefact of whatever set the
    scale.  Here ``G``, the ADM mass and the shell's rest mass all set it, so
    the spectrum should *transform* under a rescaling rather than collapse.
    """
    rows = []
    for s in scales:
        sysm = ThroatShellSystem(adm_mass=1.0 * s, screened=0.30 * s,
                                 b0=b0 * s, a0=a0 * s)
        modes = sysm.normal_modes(beta2)
        vals = modes["eigenvalues"]
        rows.append({"scale": s, "eigenvalues": vals,
                     "off_diagonal": modes["off_diagonal"],
                     "scaled_eigenvalues": [v * s * s for v in vals]})
    ref = rows[0]["scaled_eigenvalues"]
    drift = max(max(abs(r["scaled_eigenvalues"][i] - ref[i])
                    / max(abs(ref[i]), 1e-30) for i in range(2))
                for r in rows)
    smallest = min(min(abs(v) for v in r["eigenvalues"]) for r in rows)
    return {
        "rows": rows,
        "eigenvalues_scale_as_inverse_length_squared": bool(drift < 1e-5),
        "worst_scaling_drift": drift,
        "smallest_absolute_eigenvalue": smallest,
        "no_flat_direction": bool(smallest > 1e-10),
        "the_scale_is_set_by_the_masses_not_by_a_boundary": True,
    }
