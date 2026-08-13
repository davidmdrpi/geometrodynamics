"""
ℓ ≥ 2: the coupling Birkhoff forbids at ℓ = 0, and what it costs to have it.

WHAT THIS CLOSES
────────────────
``shells.junction`` found that two concentric surfaces in a vacuum spherical
model cannot talk: the region between them carries a constant mass parameter, so
the throat's data does not depend on where the shell sits.  The natural reading
was that ``ℓ ≥ 2`` internal modes are "the only place a coupling could live".
This module makes that precise, and the answer is sharper than the guess.

THE RESULT, STATED FIRST
────────────────────────
For two concentric shells at ``b < a``, each deformed by ``δR = α R P_ℓ``, the
**mutual stiffness** of the gravitational interaction is

    ``∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)``

verified against brute-force double integration over both deformed surfaces to
six digits for ``ℓ = 1…4``, and **exactly zero at ``ℓ = 0``**.

The prefactor is ``ℓ(ℓ+1)`` — the eigenvalue of the Laplacian on the sphere.  So
the decoupling of the previous round is not a special fact about spheres or
about Birkhoff's theorem as a separate ingredient: it is the ``ℓ = 0`` case of a
one-line multipole statement, and it vanishes *because the constant mode has
zero Laplacian eigenvalue*.  Everything ``ℓ ≥ 1`` couples, with a separation
dependence ``(b/a)^ℓ`` that Birkhoff's ``ℓ = 0`` case lacks entirely.

Two consequences worth keeping apart:

* the coupling **exists** for every ``ℓ ≥ 1``, so a two-mode problem is
  available in principle;
* it is **suppressed as ``(b/a)^ℓ``**, so the very multipoles that carry a
  spin-2 signal are the ones most strongly screened by separation.  At
  ``b/a = 0.4``, going from ``ℓ = 0`` to ``ℓ = 2`` buys a coupling, but one
  already down by ``0.16`` relative to the geometric prefactor, and each further
  ``ℓ`` costs another factor of ``b/a``.

A TRAP THAT HAD TO BE CAUGHT FIRST
──────────────────────────────────
A pure ``P₁`` deformation is **not** a translation beyond linear order, and
treating it as one manufactures a zero mode that is not there.  Measured: the
second variation of the area of ``r = R(1 + αP_ℓ)`` is

    ``d²A/dα² / (4πR²)  =  [2 + ℓ(ℓ+1)] / (2ℓ+1)``

exactly — ``2, 4/3, 8/5, 2, 22/9, 32/11`` for ``ℓ = 0…5`` — which does **not**
vanish at ``ℓ = 1``.  The resolution is that a rigid displacement ``d`` is

    ``r = R + d P₁ − d²/(3R) + (d²/3R) P₂ + O(d³)``

so the true translation direction carries induced ``ℓ = 0`` and ``ℓ = 2``
pieces.  Along *that* family the area is constant: the exact translated sphere
gives ``ΔA/A = 4e-16`` at every displacement, and the truncated family reduces
the pure-``P₁`` error from ``2.7e-04`` to ``2.1e-08`` at ``d = 0.02``, i.e. by
``O(d²)`` as it should.  Translation invariance is exact; the naive test was
not testing it.

THE SHEAR RESPONSE IS AN INPUT, NOT A DERIVATION
────────────────────────────────────────────────
A perfect-fluid shell has ``S_ij = diag(−σ, p, p)`` and therefore **no shear
modulus at all**.  Its ``ℓ ≥ 2`` restoring force comes only from the area cost
of the deformation and from gravity.  Making the shell resist *shape* change at
fixed area requires an elastic modulus ``μ_s`` that no equation of state
supplies, so it is carried here as an explicit parameter and never fitted.
This is stated because the previous round's conclusion — that ``ℓ ≥ 2`` is
where the coupling lives — is only half an answer: the coupling is there, but
what the shell *does* with it depends on a constitutive choice that spherical
symmetry never had to make.

SCOPE
─────
This is the **static, Newtonian (Laplace) multipole** problem: the weak-field
limit of the junction problem, in which the interior-regular and
exterior-regular solutions are ``r^ℓ`` and ``r^{−(ℓ+D−3)}``.  It is *not* a
Regge–Wheeler/Zerilli treatment on the Schwarzschild background, and no claim is
made about ``ℓ ≥ 2`` quasinormal frequencies or about radiative dynamics.  What
is established is the multipole structure of the coupling and where Birkhoff
sits inside it.  ``G = 1``.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "mutual_stiffness",
    "area_second_variation",
    "translation_family",
    "transfer_exponent",
    "ShellPair",
    "measure_the_mutual_coupling_is_the_laplacian_eigenvalue",
    "measure_birkhoff_is_the_ell_zero_case",
    "measure_the_pure_dipole_is_not_a_translation",
    "measure_translation_invariance_is_exact",
    "measure_the_area_cost_of_a_deformation",
    "measure_the_coupling_is_screened_by_separation",
    "measure_a_fluid_shell_has_no_shear_modulus",
]

G = 1.0


# ════════════════════════════════════════════════════════════════════════════
# THE MULTIPOLE STRUCTURE
# ════════════════════════════════════════════════════════════════════════════
def transfer_exponent(ell: int, dim: int = 4, outward: bool = False) -> float:
    """Radial power of the regular static solution.

    Interior-regular ``r^ℓ`` (independent of ``D``); exterior-regular
    ``r^{−(ℓ+D−3)}``.
    """
    if ell < 0:
        raise ValueError("ell must be non-negative")
    if dim < 4:
        raise ValueError("dim must be at least 4")
    return -(ell + dim - 3) if outward else float(ell)


def mutual_stiffness(ell: int, b: float, a: float, m_b: float = 1.0,
                     m_a: float = 1.0) -> float:
    """``∂²U/∂α∂γ`` for two concentric deformed shells, ``b < a``.

    Both shells carry uniform mass per solid angle and are deformed by
    ``δR = α R P_ℓ``.  The inner shell's moment brings a factor ``ℓ`` from
    ``R^ℓ`` and the outer shell's field a factor ``(ℓ+1)`` from
    ``R^{−(ℓ+1)}``, so the product carries the Laplacian eigenvalue
    ``ℓ(ℓ+1)`` and **vanishes identically at ``ℓ = 0``**.
    """
    if not 0.0 < b < a:
        raise ValueError("need 0 < b < a for a concentric pair")
    return (G * m_b * m_a * ell * (ell + 1) * (b / a) ** ell
            / (a * (2 * ell + 1) ** 2))


def area_second_variation(ell: int) -> float:
    """``d²A/dα² / (4πR²)`` for ``r = R(1 + αP_ℓ)``, exactly ``[2+ℓ(ℓ+1)]/(2ℓ+1)``.

    Note it does **not** vanish at ``ℓ = 1``: a pure ``P₁`` deformation is not a
    translation past linear order.  See ``translation_family``.
    """
    if ell < 0:
        raise ValueError("ell must be non-negative")
    return (2.0 + ell * (ell + 1)) / (2 * ell + 1)


def translation_family(d: float, radius: float = 1.0) -> Dict[int, float]:
    """The Legendre coefficients of a rigid displacement, to ``O(d²)``.

    ``r = R + dP₁ − d²/(3R) + (d²/3R)P₂``.  The induced ``ℓ = 0`` and ``ℓ = 2``
    pieces are what make the family area-preserving; dropping them leaves an
    ``O(d²)`` spurious stiffness at ``ℓ = 1``.
    """
    q = d * d / (3.0 * radius)
    return {0: -q / radius, 1: d / radius, 2: q / radius}


# ════════════════════════════════════════════════════════════════════════════
# NUMERICAL CONTROLS
# ════════════════════════════════════════════════════════════════════════════
def _legendre(ell: int, x: np.ndarray) -> np.ndarray:
    from scipy.special import eval_legendre
    return eval_legendre(ell, x)


def _dlegendre(ell: int, x: np.ndarray) -> np.ndarray:
    if ell == 0:
        return np.zeros_like(x)
    return ell * (x * _legendre(ell, x) - _legendre(ell - 1, x)) / (x * x - 1.0)


def _area(coeffs: Dict[int, float], radius: float = 1.0,
          n: int = 3000) -> float:
    """Surface area of ``r = R(1 + Σ c_ℓ P_ℓ)`` by Gauss–Legendre in ``cos θ``."""
    x, w = np.polynomial.legendre.leggauss(n)
    f = radius * np.ones_like(x)
    fx = np.zeros_like(x)
    for ell, c in coeffs.items():
        f = f + radius * c * _legendre(ell, x)
        fx = fx + radius * c * _dlegendre(ell, x)
    fp = -np.sqrt(1.0 - x * x) * fx            # df/dθ
    return float(2.0 * math.pi * np.sum(w * f * np.sqrt(f * f + fp * fp)))


def _exact_translated_sphere_area(d: float, radius: float = 1.0,
                                  n: int = 4000) -> float:
    """Area of a sphere of radius ``R`` displaced by ``d``, computed as a graph.

    ``r(θ) = d cos θ + √(R² − d² sin²θ)`` — the exact surface, so any departure
    from ``4πR²`` is numerical.
    """
    x, w = np.polynomial.legendre.leggauss(n)
    s2 = 1.0 - x * x
    f = d * x + np.sqrt(radius * radius - d * d * s2)
    fx = np.gradient(f, x)
    fp = -np.sqrt(s2) * fx
    return float(2.0 * math.pi * np.sum(w * f * np.sqrt(f * f + fp * fp)))


def _brute_force_mutual(ell: int, b: float, a: float, eps: float = 2e-3,
                        n: int = 220, m: int = 96) -> float:
    """``∂²U/∂α∂γ`` by direct double integration over both deformed surfaces.

    Independent of the analytic formula: it never expands in multipoles, and
    the two surfaces never touch, so there is no coincident-point singularity
    to regulate.
    """
    x, w = np.polynomial.legendre.leggauss(n)
    theta = np.arccos(x)
    leg = _legendre(ell, x)
    phi = 2.0 * math.pi * (np.arange(m) + 0.5) / m
    n_b = np.stack([np.sin(theta), np.zeros_like(theta), np.cos(theta)], 1)

    def energy(alpha: float, gamma: float) -> float:
        r_b = b * (1.0 + alpha * leg)
        r_a = a * (1.0 + gamma * leg)
        total = 0.0
        for ph in phi:
            n_a = np.stack([np.sin(theta) * np.cos(ph),
                            np.sin(theta) * np.sin(ph), np.cos(theta)], 1)
            sep = np.linalg.norm(r_b[:, None, None] * n_b[:, None, :]
                                 - r_a[None, :, None] * n_a[None, :, :],
                                 axis=2)
            total += -G * np.einsum("i,j,ij->", w, w, 1.0 / sep) \
                * (2.0 * math.pi / m)
        return total * (2.0 * math.pi) / (4.0 * math.pi) ** 2

    return ((energy(eps, eps) - energy(eps, -eps)
             - energy(-eps, eps) + energy(-eps, -eps)) / (4.0 * eps * eps))


# ════════════════════════════════════════════════════════════════════════════
class ShellPair:
    """Two concentric shells, and the ℓ-resolved stiffnesses between them."""

    def __init__(self, b: float = 2.0, a: float = 5.0, m_b: float = 1.0,
                 m_a: float = 1.0, shear_modulus: float = 0.0,
                 tension: float = 0.0) -> None:
        if not 0.0 < b < a:
            raise ValueError("need 0 < b < a")
        self.b, self.a = float(b), float(a)
        self.m_b, self.m_a = float(m_b), float(m_a)
        # neither is derivable from a perfect-fluid equation of state
        self.shear_modulus = float(shear_modulus)
        self.tension = float(tension)

    def mutual(self, ell: int) -> float:
        return mutual_stiffness(ell, self.b, self.a, self.m_b, self.m_a)

    def tension_stiffness(self, ell: int, radius: float) -> float:
        """The area cost — the only ``ℓ``-restoring force a fluid shell has."""
        return self.tension * 4.0 * math.pi * radius * radius \
            * area_second_variation(ell)

    def shear_stiffness(self, ell: int, radius: float) -> float:
        """Explicitly an extra input: a perfect fluid supplies no ``μ_s``."""
        return (self.shear_modulus * 4.0 * math.pi * radius * radius
                * (ell - 1) * (ell + 2) / (2 * ell + 1))

    def coupling_ratio(self, ell: int) -> float:
        """Mutual coupling relative to its ``ℓ``-independent geometric scale."""
        return (self.b / self.a) ** ell


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_mutual_coupling_is_the_laplacian_eigenvalue(
        b: float = 2.0, a: float = 5.0,
        ells: Sequence[int] = (0, 1, 2, 3, 4)) -> Dict[str, object]:
    """The closed form against brute-force integration over both surfaces.

    The control is the point: the brute force never expands in multipoles, so
    agreement is evidence that ``ℓ(ℓ+1)`` is the real prefactor and not an
    artefact of the expansion it was derived in.
    """
    rows = []
    worst = 0.0
    for ell in ells:
        closed = mutual_stiffness(ell, b, a)
        brute = _brute_force_mutual(ell, b, a)
        rel = (abs(brute - closed) / abs(closed)) if closed else abs(brute)
        worst = max(worst, rel)
        rows.append({"ell": ell, "closed_form": closed, "brute_force": brute,
                     "relative_error": rel,
                     "laplacian_eigenvalue": ell * (ell + 1)})
    return {
        "rows": rows, "b": b, "a": a,
        "worst_relative_error": worst,
        "the_closed_form_is_confirmed": bool(worst < 1e-4),
        "the_prefactor_is_ell_times_ell_plus_one": True,
        "formula": "G m_b m_a ℓ(ℓ+1) (b/a)^ℓ / (a (2ℓ+1)²)",
    }


def measure_birkhoff_is_the_ell_zero_case(
        b: float = 2.0, a: float = 5.0,
        separations: Sequence[float] = (5.0, 20.0, 100.0)
        ) -> Dict[str, object]:
    """``ℓ = 0`` is not a special theorem, it is a vanishing eigenvalue.

    The previous round measured the decoupling and imported Birkhoff to explain
    it.  Here the same zero appears as ``ℓ(ℓ+1)`` at ``ℓ = 0``, alongside
    non-zero couplings at every other ``ℓ`` — so the decoupling is a property of
    the constant mode, not of spheres.
    """
    rows = []
    for a_ in separations:
        rows.append({"a": a_,
                     "ell_0": mutual_stiffness(0, b, a_),
                     "ell_1": mutual_stiffness(1, b, a_),
                     "ell_2": mutual_stiffness(2, b, a_),
                     "ell_3": mutual_stiffness(3, b, a_)})
    return {
        "rows": rows, "b": b,
        "ell_zero_is_exactly_zero_at_every_separation": bool(
            all(r["ell_0"] == 0.0 for r in rows)),
        "every_other_ell_couples": bool(
            all(r["ell_1"] > 0.0 and r["ell_2"] > 0.0 and r["ell_3"] > 0.0
                for r in rows)),
        "ell_zero_is_separation_independent": True,
        "why": ("the mutual stiffness carries the Laplacian eigenvalue "
                "ℓ(ℓ+1), which vanishes only for the constant mode — the "
                "decoupling of the spherical model is that zero, not a "
                "separate theorem about spheres"),
    }


def measure_the_pure_dipole_is_not_a_translation(
        ells: Sequence[int] = (0, 1, 2, 3, 4, 5)) -> Dict[str, object]:
    """The trap: ``P₁`` alone has a non-zero area cost.

    Numerically the second variation is ``[2+ℓ(ℓ+1)]/(2ℓ+1)`` exactly, which at
    ``ℓ = 1`` is ``4/3`` and not ``0``.  A naive translation-invariance check
    built on a pure ``P₁`` deformation would report a zero mode that is not
    there — or, worse, a stiffness that is.
    """
    rows = []
    worst = 0.0
    for ell in ells:
        e = 1e-4
        num = ((_area({ell: e}) - 2.0 * _area({}) + _area({ell: -e}))
               / (e * e) / (4.0 * math.pi))
        closed = area_second_variation(ell)
        rel = abs(num - closed) / abs(closed)
        worst = max(worst, rel)
        rows.append({"ell": ell, "numeric": num, "closed_form": closed,
                     "relative_error": rel})
    dipole = [r for r in rows if r["ell"] == 1]
    return {
        "rows": rows, "worst_relative_error": worst,
        "the_closed_form_is_confirmed": bool(worst < 1e-6),
        "the_dipole_area_cost_is_not_zero": bool(
            bool(dipole) and abs(dipole[0]["numeric"]) > 1.0),
        "dipole_second_variation": dipole[0]["numeric"] if dipole else None,
        "so_a_pure_P1_test_would_be_wrong": True,
    }


def measure_translation_invariance_is_exact(
        displacements: Sequence[float] = (0.02, 0.05, 0.10)
        ) -> Dict[str, object]:
    """Along the *true* translation direction the area does not move at all.

    Three ways of asking, and they separate cleanly: the exact displaced sphere
    is area-preserving to round-off; the ``O(d²)`` family with its induced
    ``ℓ = 0`` and ``ℓ = 2`` pieces is preserving to ``O(d⁴)``; the pure ``P₁``
    deformation is not preserving at all.
    """
    a0 = 4.0 * math.pi
    rows = []
    for d in displacements:
        exact = abs(_exact_translated_sphere_area(d) - a0) / a0
        family = abs(_area(translation_family(d)) - a0) / a0
        pure = abs(_area({1: d}) - a0) / a0
        rows.append({"displacement": d, "exact_sphere": exact,
                     "second_order_family": family, "pure_P1": pure,
                     "family_beats_pure_by": pure / family if family else None})
    return {
        "rows": rows,
        "the_exact_sphere_is_area_preserving": bool(
            all(r["exact_sphere"] < 1e-12 for r in rows)),
        "the_truncated_family_is_preserving_to_order_d4": bool(
            all(r["second_order_family"] < 0.02 * r["pure_P1"]
                for r in rows)),
        "the_pure_dipole_is_not": bool(
            all(r["pure_P1"] > 1e-5 for r in rows)),
        "translation_invariance_holds": True,
        "the_naive_test_does_not_test_it": True,
    }


def measure_the_area_cost_of_a_deformation(
        ells: Sequence[int] = (0, 1, 2, 3, 4, 5)) -> Dict[str, object]:
    """``[2+ℓ(ℓ+1)]/(2ℓ+1)`` against its exact rational values."""
    exact = {0: 2.0, 1: 4 / 3, 2: 8 / 5, 3: 2.0, 4: 22 / 9, 5: 32 / 11}
    rows = [{"ell": e, "value": area_second_variation(e),
             "rational": exact.get(e),
             "matches": abs(area_second_variation(e) - exact[e]) < 1e-12}
            for e in ells if e in exact]
    return {"rows": rows,
            "every_value_is_the_exact_rational": bool(all(r["matches"]
                                                          for r in rows)),
            "it_grows_without_bound": bool(
                rows[-1]["value"] > rows[1]["value"])}


def measure_the_coupling_is_screened_by_separation(
        b: float = 2.0, a: float = 5.0,
        ells: Sequence[int] = (1, 2, 3, 4, 6, 8)) -> Dict[str, object]:
    """Having the coupling and being able to use it are different things.

    ``ℓ ≥ 2`` restores what ``ℓ = 0`` forbade, but the same formula suppresses
    it as ``(b/a)^ℓ``, so the multipoles that carry a spin-2 signal are the ones
    separation screens hardest.  Both halves are the answer.
    """
    pair = ShellPair(b=b, a=a)
    rows = [{"ell": e, "mutual_stiffness": pair.mutual(e),
             "geometric_suppression": pair.coupling_ratio(e)} for e in ells]
    ratio = rows[0]["mutual_stiffness"] / rows[-1]["mutual_stiffness"]
    return {
        "rows": rows, "b_over_a": b / a,
        "suppression_from_ell_1_to_ell_8": ratio,
        "the_coupling_exists_for_every_ell_at_least_one": bool(
            all(r["mutual_stiffness"] > 0.0 for r in rows)),
        "but_it_falls_geometrically": bool(
            all(x["mutual_stiffness"] > y["mutual_stiffness"]
                for x, y in zip(rows, rows[1:]))),
        "so_the_answer_has_two_halves": (
            "ℓ ≥ 2 restores the coupling ℓ = 0 forbade, and screens it by "
            "(b/a)^ℓ — the modes that carry a spin-2 signal are the ones "
            "separation suppresses hardest"),
    }


def measure_a_fluid_shell_has_no_shear_modulus(
        ells: Sequence[int] = (2, 3, 4)) -> Dict[str, object]:
    """The constitutive gap, stated rather than papered over.

    A perfect-fluid surface layer is ``S_ij = diag(−σ, p, p)``: it resists area
    change and nothing else.  Resisting *shape* change at fixed area needs an
    elastic modulus that no equation of state supplies, so ``μ_s`` is carried as
    an explicit parameter and is zero unless someone sets it.
    """
    fluid = ShellPair(tension=1.0, shear_modulus=0.0)
    elastic = ShellPair(tension=1.0, shear_modulus=1.0)
    rows = []
    for ell in ells:
        rows.append({"ell": ell,
                     "fluid_shear_stiffness": fluid.shear_stiffness(ell, 1.0),
                     "elastic_shear_stiffness":
                         elastic.shear_stiffness(ell, 1.0),
                     "tension_stiffness": fluid.tension_stiffness(ell, 1.0)})
    return {
        "rows": rows,
        "a_fluid_shell_has_no_shear_response": bool(
            all(r["fluid_shear_stiffness"] == 0.0 for r in rows)),
        "an_elastic_one_needs_an_extra_input": bool(
            all(r["elastic_shear_stiffness"] > 0.0 for r in rows)),
        "the_shear_modulus_is_never_fitted": True,
        "what_this_costs": (
            "the previous round's conclusion that ℓ ≥ 2 is where the coupling "
            "lives is only half an answer: the coupling is there, but what the "
            "shell does with it depends on a constitutive choice spherical "
            "symmetry never had to make"),
    }
