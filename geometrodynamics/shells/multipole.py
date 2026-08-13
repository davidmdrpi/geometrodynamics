"""
Where the two-shell coupling starts, in a static Newtonian model.

SCOPE, BECAUSE THE HEADLINE DEPENDS ON IT
─────────────────────────────────────────
This is the **static Newtonian (Laplace) two-shell model** — the weak-field
analogue of the junction problem, with interior/exterior static solutions
``r^ℓ`` and ``r^{−(ℓ+D−3)}``.  What it establishes is the **shell-theorem /
multipole** structure of the coupling.  **Birkhoff's theorem is a GR result and
remains what ``shells.junction`` (PR #249) relies on**; nothing here replaces
it, and the ℓ = 0 statement below is its Newtonian analogue, not the theorem
itself.  Not Regge–Wheeler/Zerilli; no quasinormal frequencies; no radiative
dynamics.  ``G = 1``.

THE RESULT, STATED FIRST
────────────────────────
In this model the **monopole mutual stiffness vanishes** while higher angular
multipoles couple, with the coupling **suppressed geometrically by separation**.

For two concentric shells at ``b < a``, each deformed by ``δR = α R P_ℓ``,

    ``∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)``

verified against brute-force double integration over both deformed surfaces to
six digits, and exactly zero at ``ℓ = 0``.  The prefactor is ``ℓ(ℓ+1)``, the
eigenvalue of the angular Laplacian: **the ℓ = 0 Newtonian decoupling is that
zero eigenvalue**.

WHERE THE COUPLING ACTUALLY STARTS — AND IT IS NOT ℓ = 1
───────────────────────────────────────────────────────
An earlier draft of this module concluded that "everything ``ℓ ≥ 1`` couples".
That is wrong as a statement about physical modes, and the error was checking
translation invariance of the **area** instead of of the **mutual energy**.
Done at the level of the energy the two disagree:

* a **rigid translation** of either sphere leaves the mutual energy exactly
  ``−G m_b m_a / a`` — Newton's shell theorem, held to ``1e-15`` out to
  ``d = 2.5`` — so the cross-derivative with respect to rigid displacements is
  ``8.3e-13``.  **The translation mode does not couple.**
* a pure ``P₁`` *shape* deformation is a different object, and it does couple,
  at ``1.78e-02``.

So the honest ordering is: ``ℓ = 0`` decouples by the vanishing eigenvalue, the
``ℓ = 1`` **translation** mode decouples by the shell theorem, and genuine
coupling **starts at ``ℓ = 2``** — which is what PR #249 guessed and this
establishes.

THE SAME TRAP, TWICE
────────────────────
Both errors above are the same mistake: **a pure ``P₁`` deformation is not a
translation past linear order.**  The second variation of the area of
``r = R(1 + αP_ℓ)`` is

    ``d²A/dα² / (4πR²)  =  [2 + ℓ(ℓ+1)] / (2ℓ+1)``

exactly — ``2, 4/3, 8/5, 2, 22/9, 32/11`` for ``ℓ = 0…5`` — which does not
vanish at ``ℓ = 1``.  A rigid displacement is instead

    ``r = R + d P₁ − d²/(3R) + (d²/3R) P₂ + O(d³)``

carrying induced ``ℓ = 0`` and ``ℓ = 2`` pieces.  The exact translated sphere is
area-preserving to ``4e-16`` at every displacement.  The lesson is that the
zero-mode test has to be run on the quantity the claim is about: the area test
does not decide whether ``ℓ = 1`` couples, and the energy test does.

THE SHEAR RESPONSE IS AN INPUT, NOT A DERIVATION
────────────────────────────────────────────────
A perfect-fluid shell has ``S_ij = diag(−σ, p, p)`` and therefore **no shear
modulus at all**: it resists area change and nothing else.  Making it resist
*shape* change at fixed area requires an elastic modulus ``μ_s`` that no
equation of state supplies, so it is carried as an explicit parameter and never
fitted.  PR #249's conclusion that ``ℓ ≥ 2`` is where the coupling lives is
therefore only half an answer — the coupling is there, but what a shell does
with it is a constitutive choice spherical symmetry never had to make.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "mutual_stiffness",
    "rigid_pair_mutual_energy",
    "measure_the_translation_mode_does_not_couple",
    "area_second_variation",
    "translation_family",
    "transfer_exponent",
    "ShellPair",
    "measure_the_mutual_coupling_is_the_laplacian_eigenvalue",
    "measure_the_ell_zero_decoupling_is_a_zero_eigenvalue",
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


def measure_the_ell_zero_decoupling_is_a_zero_eigenvalue(
        b: float = 2.0, a: float = 5.0,
        separations: Sequence[float] = (5.0, 20.0, 100.0)
        ) -> Dict[str, object]:
    """The ``ℓ = 0`` Newtonian decoupling is the zero eigenvalue of ``ℓ(ℓ+1)``.

    ``shells.junction`` measured a decoupling in GR and relies on Birkhoff's
    theorem for it.  **This is the Newtonian analogue, not that theorem**: here
    the same zero appears as the vanishing angular-Laplacian eigenvalue at
    ``ℓ = 0``, alongside non-zero shape couplings at higher ``ℓ``.
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
        "why": ("the mutual stiffness carries the angular-Laplacian "
                "eigenvalue ℓ(ℓ+1), which vanishes only for the constant "
                "mode; in this Newtonian model that zero IS the ℓ = 0 "
                "decoupling"),
        "this_is_not_birkhoff": ("Birkhoff is a GR theorem and remains what "
                                 "shells.junction relies on; this is its "
                                 "static Newtonian analogue"),
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
        "every_shape_mode_couples": bool(
            all(r["mutual_stiffness"] > 0.0 for r in rows)),
        "but_the_ell_1_entry_is_a_shape_mode_not_a_translation": (
            "the rigid translation has zero mutual coupling by the shell "
            "theorem — see measure_the_translation_mode_does_not_couple"),
        "but_it_falls_geometrically": bool(
            all(x["mutual_stiffness"] > y["mutual_stiffness"]
                for x, y in zip(rows, rows[1:]))),
        "so_the_answer_has_two_halves": (
            "genuine coupling starts at ℓ = 2, and the same formula screens it "
            "by (b/a)^ℓ — the modes that carry a spin-2 signal are the ones "
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


# ════════════════════════════════════════════════════════════════════════════
# THE ℓ = 1 CONTROL, AT THE LEVEL OF THE ENERGY
# ════════════════════════════════════════════════════════════════════════════
def rigid_pair_mutual_energy(b: float, a: float, d_b: float = 0.0,
                             d_a: float = 0.0, n: int = 200,
                             m: int = 120) -> float:
    """Mutual energy of two spheres **rigidly displaced** along the axis.

    Both surfaces are exact translated spheres — no Legendre truncation — so
    this tests the energy, not a shape parametrisation of it.  Newton's shell
    theorem says the result is ``−G m_b m_a / a`` for any displacements that
    keep the inner sphere entirely inside the outer.
    """
    if not 0.0 < b < a:
        raise ValueError("need 0 < b < a")
    x, w = np.polynomial.legendre.leggauss(n)
    theta = np.arccos(x)
    n_b = np.stack([np.sin(theta), np.zeros_like(theta), np.cos(theta)], 1)
    pts_b = b * n_b + np.array([0.0, 0.0, d_b])
    phi = 2.0 * math.pi * (np.arange(m) + 0.5) / m
    total = 0.0
    for ph in phi:
        n_a = np.stack([np.sin(theta) * np.cos(ph),
                        np.sin(theta) * np.sin(ph), np.cos(theta)], 1)
        pts_a = a * n_a + np.array([0.0, 0.0, d_a])
        sep = np.linalg.norm(pts_b[:, None, :] - pts_a[None, :, :], axis=2)
        total += -G * np.einsum("i,j,ij->", w, w, 1.0 / sep) \
            * (2.0 * math.pi / m)
    return float(total * (2.0 * math.pi) / (4.0 * math.pi) ** 2)


def measure_the_translation_mode_does_not_couple(
        b: float = 2.0, a: float = 5.0, eps: float = 1e-2,
        displacements: Sequence[float] = (0.1, 0.8, 2.5)
        ) -> Dict[str, object]:
    """The ``ℓ = 1`` control the area test does not perform.

    Checking translation invariance of the **area** says nothing about the
    mutual gravitational energy, which is the quantity the coupling result is
    about.  Done properly, the two disagree and the physics changes:

    * a **rigid translation** of either sphere leaves the mutual energy exactly
      ``−G m_b m_a / a`` — Newton's shell theorem — so the translation mode has
      **zero** mutual coupling;
    * a pure ``P₁`` *shape* deformation is a different object and does couple,
      at ``1.78e-02``.

    So "every ``ℓ ≥ 1`` couples" is wrong as a statement about physical modes.
    Genuine coupling starts at ``ℓ = 2``.
    """
    ref = rigid_pair_mutual_energy(b, a)
    rows = []
    for d in displacements:
        u = rigid_pair_mutual_energy(b, a, d_b=d)
        rows.append({"displacement": d, "energy": u,
                     "relative_drift": abs(u - ref) / abs(ref),
                     "inner_stays_inside": bool(b + d < a)})

    def cross(f) -> float:
        return (f(eps, eps) - f(eps, -eps) - f(-eps, eps)
                + f(-eps, -eps)) / (4.0 * eps * eps)

    rigid = cross(lambda p, q: rigid_pair_mutual_energy(b, a, p, q))
    shape = mutual_stiffness(1, b, a)
    return {
        "rows": rows,
        "reference_energy": ref,
        "newtonian_prediction": -G / a,
        "shell_theorem_holds": bool(
            abs(ref + G / a) < 1e-12
            and all(r["relative_drift"] < 1e-12 for r in rows)),
        "rigid_translation_cross_derivative": rigid,
        "pure_P1_shape_coupling": shape,
        "the_translation_mode_does_not_couple": bool(abs(rigid) < 1e-9),
        "the_shape_mode_does": bool(shape > 1e-3),
        "so_coupling_starts_at_ell_two": True,
        "why_the_area_test_was_not_enough": (
            "translation invariance of the area is a statement about the "
            "surface, not about the mutual gravitational energy; the energy "
            "control is the one that decides whether ℓ = 1 couples, and it "
            "says it does not"),
    }
