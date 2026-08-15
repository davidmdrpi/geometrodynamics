"""
A flux-conserving two-mouth throat operator, and where its family is stable.

WHAT PR #255 OWED
─────────────────
`waves.branch_coupling` solved a mouth relation self-consistently and was
explicit that the relation was not a boundary operator: it carried field
**values** through the free Green function, with no normal-derivative matching,
no reflected channel, a ``1×1`` mouth object where a conserving junction needs
``2×2``, and ``κ²`` power throughput.

This module builds the real object.  A point-supported throat is a
**self-adjoint extension** of the Laplacian on ``S³ ∖ {M⁺, M⁻}``, and von
Neumann's theorem says the extensions are parametrized by a unitary between the
deficiency spaces — here ``U(2)`` — equivalently, by Krein's formula, a
Hermitian ``2×2`` matrix.

THE FREE GREEN FUNCTION, IN CLOSED FORM
───────────────────────────────────────
Fourier-transforming PR #254's image sum (or summing PR #255's branch series)
gives, on the conformally coupled ESU,

    ``G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))``

**real** for real ``ω``, with poles exactly at the free spectrum ``ω = n+1`` and
— because the numerator's zero cancels ``sin χ`` — *finite* at the antipode.
Its short-distance expansion is ``1/(4πχ) + g(ω) + O(χ)`` with

    ``g(ω) = −(ω/4π) · cot(πω)``

the **regularized** value at coincidence: the divergence is the universal
Coulomb one, so the subtraction a point interaction needs is forced rather than
chosen.

THE OPERATOR
────────────
Write the field as free plus two point sources,
``φ(x) = φ_in(x) + Σ_j G(χ_j, ω) q_j``.  Near mouth ``j`` this is
``q_j/(4πχ_j) + φ_j^reg + O(χ_j)`` with ``φ^reg = φ_in|_mouths + Γ q``, and a
**linear boundary condition** is a pair of matrices,

    ``B φ^reg = C q``      ⟹      ``(C − BΓ) q = B φ_in|_mouths``

so the mouth-active spectrum is ``det(C − BΓ(ω)) = 0``.  The extension is
self-adjoint iff ``rank[B|C] = 2`` and ``BC†`` is Hermitian.  The familiar case
``B = I``, ``C = A`` recovers ``φ^reg = A q`` with ``A`` Hermitian, and the
general form is needed because **PR #255's relation is not of that shape** — see
below.

WHAT SELF-ADJOINTNESS BUYS, AND WHAT IT DOES NOT
────────────────────────────────────────────────
* **It buys flux conservation, exactly.**  The current through a small sphere at
  mouth ``j`` is ``Im(q_j* φ_j^reg)``, independent of the sphere, so the total
  absorbed is ``Im(q† A q)`` — zero for *every* ``q`` iff ``A = A†``.  Measured.
* **It buys a real spectrum in ``λ = ω²``**, which is the eigenvalue of the
  spatial operator.  ``Γ`` is real symmetric for real ``λ`` of either sign, so
  the secular function is real and its roots are real ``λ``.
* **It does not buy ``λ ≥ 0``, and therefore does not buy stability.**  A
  negative eigenvalue means ``ω = ±i√|λ|`` and one member of the pair grows like
  ``e^{√|λ| t}``.  Positivity is a **separate** condition on the boundary data,
  and for the exchange-symmetric family it is exactly

      ``α + β ≥ g₀ + G₀``  and  ``α − β ≥ g₀ − G₀`` ,
      ``g₀ = −1/(4π²)`` ,  ``G₀ = (π−d)/(4π² sin d)``

  — the threshold values of ``g ± G_d`` at ``λ = 0``, both monotone decreasing
  along the imaginary axis.  Outside that wedge the throat has a growing mode.
  This module maps the region rather than asserting it.

THE SECTORS OF THE SPECTRUM
───────────────────────────
``det(C − BΓ) = 0`` is the **mouth-active** sector only.  Level ``n`` on ``S³``
has degeneracy ``(n+1)²``, and only the two combinations with support at the
mouths can move: ``(n+1)² − 2`` modes stay exactly at the free eigenvalue.  Any
full-spectrum statement has to carry them.

Within the mouth-active sector there are three regimes, and the first two are
easy to miss by scanning ``ω`` on an interval that starts above ``1``:

* ``λ < 0`` — the growing modes above;
* ``0 ≤ λ < 1`` — stable modes **below the free ground state**;
* ``λ ∈ (m², (m+1)²)`` — the interlacing regime, where each channel function
  climbs across the gap.  ``g + G_d`` runs the full ``−∞ → +∞`` there, but
  ``g − G_d`` does **not** at the first gap: the ``m = 1`` pole cancels because
  the constant mode has the same value at both mouths, so the antisymmetric
  channel's endpoint is finite and a root there is conditional on ``α − β``.
  Measured, not assumed.

WHAT IS STILL PUT IN
────────────────────
A **linear** field on a **fixed** round background.  The boundary data is a
choice — four real parameters — not a derivation: `shells.junction` (PR #249) is
what would fix it from a matter model, and nothing here computes the
exotic-matter bill.  The throat is **point-supported**, so it has no interior and
no proper length; the "delay" of PRs #251–#255 is not a parameter of a point
extension and does not survive into it.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "TWO_PI",
    "free_green",
    "regularized_green",
    "esu_free_spectrum",
    "untouched_multiplicity",
    "gamma_at",
    "BoundaryCondition",
    "MouthPair",
    "DirectionalThroat",
    "boundary_mixing",
    "boundary_from_scattering",
    "stability_thresholds",
    "negative_lambda_modes",
    "is_stable",
    "stability_map",
    "mouth_active_spectrum",
    "spectrum_by_channel",
    "channel_endpoints",
    "mode_charges",
    "mouth_flux",
    "measure_the_green_function_has_a_finite_part",
    "measure_hermiticity_is_flux_conservation",
    "measure_self_adjointness_makes_lambda_real_not_omega",
    "measure_the_stability_region_in_the_boundary_family",
    "measure_the_mouth_active_sector_is_rank_two",
    "measure_the_pr255_boundary_condition_is_not_self_adjoint",
    "measure_the_spectrum_is_conjugation_symmetric_in_beta",
]

TWO_PI = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE FREE GREEN FUNCTION
# ════════════════════════════════════════════════════════════════════════════
def free_green(chi: float, omega: float) -> float:
    """``G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))``.

    The frequency-domain retarded Green function of the conformally coupled
    massless scalar on the ESU, in closed form — the Fourier transform of PR
    #254's image sum and the ``γ → 0`` limit of PR #255's branch series, checked
    against the latter rather than asserted.  Real for real ``ω``.
    """
    s = math.sin(chi)
    if abs(s) < 1e-13:
        raise ValueError("the Green function is singular at χ = 0 and χ = π")
    sp = math.sin(math.pi * omega)
    if abs(sp) < 1e-13:
        raise ValueError("ω is a free eigenfrequency; G has a pole there")
    return math.sin(omega * (math.pi - chi)) / (4.0 * math.pi * s * sp)


def regularized_green(omega: float) -> float:
    """``g(ω) = −(ω/4π) cot(πω)`` — the finite part at coincidence.

    ``G(χ,ω) = 1/(4πχ) + g(ω) + O(χ)``.  The ``1/(4πχ)`` is the universal
    Coulomb singularity every point interaction subtracts; what is left is what
    the throat couples to.  Poles at ``ω ∈ ℤ`` with residue ``−m/(4π²)``.
    """
    sp = math.sin(math.pi * omega)
    if abs(sp) < 1e-13:
        raise ValueError("ω is a free eigenfrequency; g has a pole there")
    return -omega * math.cos(math.pi * omega) / (4.0 * math.pi * sp)


def esu_free_spectrum(n_max: int = 8) -> List[float]:
    """``ω_n = n + 1``, i.e. ``λ_n = (n+1)²`` — the free ESU spectrum."""
    return [float(n + 1) for n in range(int(n_max))]


def untouched_multiplicity(n: int) -> int:
    """``(n+1)² − 2`` — how many modes of level ``n`` the throat cannot move.

    Level ``n`` on ``S³`` has degeneracy ``(n+1)²``.  A two-point interaction is
    a **rank-two** perturbation: only the two combinations with support at the
    mouths shift, and every other eigenfunction — every one that vanishes at
    both mouths — stays exactly at ``λ = (n+1)²``.  ``det(C − BΓ) = 0`` never
    sees them, so no statement derived from it is a full-spectrum statement.
    """
    return max(0, (int(n) + 1) ** 2 - 2)


# ════════════════════════════════════════════════════════════════════════════
# Γ, ON BOTH SIDES OF λ = 0
# ════════════════════════════════════════════════════════════════════════════
def gamma_at(lmbda: float, separation: float) -> np.ndarray:
    """``Γ(λ)`` — the free propagator between the mouths, for real ``λ`` of
    either sign.

    For ``λ > 0`` this is ``ω = √λ`` and the trigonometric forms.  For ``λ < 0``
    it is ``ω = iσ``, where

        ``g = −(σ/4π) coth(πσ)`` ,
        ``G_d = sinh(σ(π−d)) / (4π sin d · sinh(πσ))``

    both **real**, both monotone decreasing in ``σ``, and both continuous
    through ``λ = 0`` where they take the threshold values ``g₀ = −1/(4π²)`` and
    ``G₀ = (π−d)/(4π² sin d)``.  Working in ``λ`` rather than ``ω`` is the whole
    point: the negative-``λ`` sector is where the growing modes live, and an
    ``ω``-scan over an interval of the positive real axis cannot see it.
    """
    d = float(separation)
    sd = math.sin(d)
    if abs(sd) < 1e-13:
        raise ValueError("the mouths must be a nondegenerate distance apart")
    lam = float(lmbda)
    if abs(lam) < 1e-14:
        g = -1.0 / (4.0 * math.pi ** 2)
        gd = (math.pi - d) / (4.0 * math.pi ** 2 * sd)
    elif lam < 0.0:
        s = math.sqrt(-lam)
        g = -s / math.tanh(math.pi * s) / (4.0 * math.pi)
        gd = (math.sinh(s * (math.pi - d))
              / (4.0 * math.pi * sd * math.sinh(math.pi * s)))
    else:
        w = math.sqrt(lam)
        sp = math.sin(math.pi * w)
        if abs(sp) < 1e-13:
            raise ValueError("λ is a free eigenvalue; Γ has a pole there")
        g = -w * math.cos(math.pi * w) / (4.0 * math.pi * sp)
        gd = math.sin(w * (math.pi - d)) / (4.0 * math.pi * sd * sp)
    return np.array([[g, gd], [gd, g]], dtype=complex)


def stability_thresholds(separation: float) -> Dict[str, float]:
    """The ``λ = 0`` values of ``g ± G_d`` — the edges of the stable wedge.

    Both channel functions are monotone decreasing along the imaginary axis from
    these values to ``−∞``, so ``g ± G_d = α ± β`` has a negative-``λ`` root
    **iff** the right-hand side falls below the threshold.  That makes stability
    of an exchange-symmetric throat a closed-form condition rather than a scan.
    """
    d = float(separation)
    g0 = -1.0 / (4.0 * math.pi ** 2)
    gd0 = (math.pi - d) / (4.0 * math.pi ** 2 * math.sin(d))
    return {"g_at_zero": g0, "G_d_at_zero": gd0,
            "symmetric_threshold": g0 + gd0,
            "antisymmetric_threshold": g0 - gd0,
            "the_rule": ("a negative-λ mode exists iff α+β < symmetric "
                         "threshold or α−β < antisymmetric threshold")}


# ════════════════════════════════════════════════════════════════════════════
# BOUNDARY CONDITIONS
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class BoundaryCondition:
    """``B φ^reg = C q`` — a linear two-mouth boundary condition.

    General enough to hold both the self-adjoint family (``B = I``, ``C = A``
    Hermitian) *and* PR #255's relation, which is **not** of that shape and
    needs a singular ``B``.  The two criteria are separate and both measured:

    * **maximal** — ``rank[B | C] = 2``, so the condition is a boundary
      condition rather than an under- or over-determination;
    * **self-adjoint** — ``B C†`` Hermitian.  This is the statement that no flux
      is created or destroyed at the mouths.
    """

    B: np.ndarray
    C: np.ndarray

    def is_maximal(self, tol: float = 1e-12) -> bool:
        return bool(np.linalg.matrix_rank(np.hstack([self.B, self.C]),
                                          tol=tol) == 2)

    def is_self_adjoint(self, tol: float = 1e-12) -> bool:
        m = self.B @ self.C.conjugate().T
        return bool(np.abs(m - m.conjugate().T).max() <= tol)

    def self_adjointness_defect(self) -> float:
        m = self.B @ self.C.conjugate().T
        return float(np.abs(m - m.conjugate().T).max())

    def secular_at(self, gamma: np.ndarray) -> complex:
        """``det(C − BΓ)`` — zero on the mouth-active spectrum."""
        return complex(np.linalg.det(self.C - self.B @ gamma))


@dataclass(frozen=True)
class MouthPair:
    """A **self-adjoint** two-mouth throat: a Hermitian ``2×2`` matrix ``A``.

    ``alpha1``, ``alpha2`` are the mouths' own real self-energies — inverse
    scattering lengths, so ``α → ∞`` *decouples* and ``α = 0`` is a strongly
    coupled throat.  ``beta`` is the complex mouth-to-mouth amplitude, the only
    part of the operator with no local realization on ``S³``: that non-locality
    is the wormhole.

    Hermiticity gives flux conservation and a real spectrum **in ``λ``**.  It
    does *not* give ``λ ≥ 0``; see `is_stable`.
    """

    separation: float
    alpha1: float = 0.05
    alpha2: float = 0.05
    beta: complex = 0.03

    def boundary_matrix(self) -> np.ndarray:
        b = complex(self.beta)
        return np.array([[self.alpha1, b], [b.conjugate(), self.alpha2]],
                        dtype=complex)

    def boundary_condition(self) -> BoundaryCondition:
        """``B = I``, ``C = A`` — the shape ``φ^reg = A q``."""
        return BoundaryCondition(np.eye(2, dtype=complex),
                                 self.boundary_matrix())

    def is_self_adjoint(self, tol: float = 1e-15) -> bool:
        a = self.boundary_matrix()
        return bool(np.abs(a - a.conjugate().T).max() <= tol)

    def is_exchange_symmetric(self, tol: float = 1e-15) -> bool:
        return bool(abs(self.alpha1 - self.alpha2) <= tol
                    and abs(complex(self.beta).imag) <= tol)

    def gamma(self, lmbda: float) -> np.ndarray:
        return gamma_at(lmbda, self.separation)

    def krein_matrix(self, lmbda: float) -> np.ndarray:
        """``M(λ) = A − Γ(λ)``.  Hermitian for real ``λ``, either sign."""
        return self.boundary_matrix() - self.gamma(lmbda)

    def secular(self, lmbda: float) -> float:
        """``det M(λ)`` — **real** for real ``λ``; its zeros are the sector."""
        return float(np.linalg.det(self.krein_matrix(lmbda)).real)

    def channel_functions(self, lmbda: float) -> Tuple[float, float]:
        """``g ± G_d`` at ``λ`` — the symmetric and antisymmetric channels."""
        gam = self.gamma(lmbda)
        g, gd = float(gam[0, 0].real), float(gam[0, 1].real)
        return g + gd, g - gd


@dataclass(frozen=True)
class DirectionalThroat:
    """PR #255's relation, embedded **exactly** — and it is not self-adjoint.

    That round drove ``M⁻`` from what arrived at ``M⁺`` and not the reverse:
    ``q₁ = 0`` (``M⁺`` absorbs but does not radiate) and
    ``q₂ = gain · φ₁^reg``.  In the ``B φ^reg = C q`` form that is

        ``B = [[0, 0], [gain, 0]]`` ,  ``C = I``

    for which ``det(C − BΓ) = 1 − gain·G_d`` — **exactly** PR #255's pole
    condition ``1 − L = 0``.  It is a maximal boundary condition (``rank[B|C] =
    2``) but ``BC† = B`` is not Hermitian, so it is not self-adjoint, and no
    finite Hermitian ``A`` reproduces it: it needs the singular ``B``.

    A previous version of this module wrote the control as
    ``A = [[0,0],[1/gain,0]]`` with ``B = I``.  That gives
    ``det(A − Γ) = g² − G_d² + G_d/gain``, which is a *different function* — so
    the old control was not PR #255's model at all.  The correction is recorded
    here rather than quietly dropped.
    """

    separation: float
    delay: float = 1.0
    eta: int = 1
    kappa: float = 0.3

    def gain(self, omega: complex) -> complex:
        return self.eta * self.kappa * cmath.exp(-1j * complex(omega)
                                                 * self.delay)

    def boundary_condition(self, omega: complex) -> BoundaryCondition:
        g = self.gain(omega)
        return BoundaryCondition(np.array([[0.0, 0.0], [g, 0.0]],
                                          dtype=complex),
                                 np.eye(2, dtype=complex))

    def secular(self, omega: complex) -> complex:
        """``det(C − BΓ(ω))``, continued to complex ``ω``."""
        w = complex(omega)
        sp = cmath.sin(math.pi * w)
        d = self.separation
        g = -w * cmath.cos(math.pi * w) / (4.0 * math.pi * sp)
        gd = cmath.sin(w * (math.pi - d)) / (4.0 * math.pi * math.sin(d) * sp)
        gam = np.array([[g, gd], [gd, g]], dtype=complex)
        return self.boundary_condition(w).secular_at(gam)

    def pr255_pole_condition(self, omega: complex) -> complex:
        """``1 − L(ω)`` as PR #255 wrote it — the independent construction."""
        w = complex(omega)
        d = self.separation
        gd = (cmath.sin(w * (math.pi - d))
              / (4.0 * math.pi * math.sin(d) * cmath.sin(math.pi * w)))
        return 1.0 - self.gain(w) * gd


# ════════════════════════════════════════════════════════════════════════════
# THE CAYLEY FORM — AND WHAT IT IS NOT
# ════════════════════════════════════════════════════════════════════════════
def boundary_mixing(boundary: np.ndarray,
                    scale: float = 0.1) -> Dict[str, object]:
    """``S = (A − ic)(A + ic)^{-1}`` — unitary for every Hermitian ``A``.

    This is von Neumann's parametrization made concrete: the self-adjoint
    two-mouth conditions **are** ``U(2)``, four real parameters.

    **Its entries are not reflection and transmission amplitudes.**  Their
    magnitudes depend on the arbitrary reference scale ``c`` — the same ``A``
    gives ``|S₁₁| = 0.654`` at ``c = 0.05`` and ``0.941`` at ``c = 0.2`` — so
    they are *boundary-mixing coefficients*, not a physical scattering matrix.
    A physical ``r`` and ``t`` would need a flux normalization, which a closed
    universe with no asymptotic region does not supply.  The physical
    conservation statement is `mouth_flux`: Hermitian ``A`` ⇔ zero net boundary
    flux.
    """
    a = np.asarray(boundary, dtype=complex)
    eye = np.eye(2, dtype=complex)
    s = (a - 1j * scale * eye) @ np.linalg.inv(a + 1j * scale * eye)
    return {
        "S": s, "scale": float(scale),
        "diagonal_mixing": [abs(complex(s[0, 0])), abs(complex(s[1, 1]))],
        "off_diagonal_mixing": [abs(complex(s[0, 1])), abs(complex(s[1, 0]))],
        "unitarity_defect": float(
            np.abs(s.conjugate().T @ s - np.eye(2)).max()),
        "column_norms": [float(abs(s[0, 0]) ** 2 + abs(s[1, 0]) ** 2),
                         float(abs(s[0, 1]) ** 2 + abs(s[1, 1]) ** 2)],
        "caveat": ("magnitudes depend on the reference scale c; these are "
                   "boundary-mixing coefficients, not physical r and t"),
    }


def boundary_from_scattering(s: np.ndarray, scale: float = 0.1) -> np.ndarray:
    """The inverse Cayley, ``A = ic (I + S)(I − S)^{-1}`` — Hermitian again."""
    eye = np.eye(2, dtype=complex)
    return 1j * scale * (eye + s) @ np.linalg.inv(eye - s)


# ════════════════════════════════════════════════════════════════════════════
# THE SPECTRUM, SECTOR BY SECTOR
# ════════════════════════════════════════════════════════════════════════════
def _bisect(f, lo: float, hi: float, n: int = 300,
            tol: float = 1e-14) -> Optional[float]:
    a, b = float(lo), float(hi)
    fa, fb = f(a), f(b)
    if not (math.isfinite(fa) and math.isfinite(fb)) or fa * fb > 0.0:
        return None
    for _ in range(int(n)):
        m = 0.5 * (a + b)
        fm = f(m)
        if not math.isfinite(fm):
            return None
        if fa * fm <= 0.0:
            b = m
        else:
            a, fa = m, fm
        if b - a < tol * max(1.0, abs(b)):
            break
    return 0.5 * (a + b)


def negative_lambda_modes(pair: MouthPair, lambda_min: float = -4000.0,
                          n_grid: int = 4000) -> List[Dict[str, float]]:
    """Roots of the secular function with ``λ < 0`` — the **growing** modes.

    A negative eigenvalue of the spatial operator means ``ω = ±i√|λ|`` and one
    member of that pair goes like ``e^{√|λ| t}``.  Self-adjointness does not
    forbid it: Hermiticity makes ``λ`` real, positivity is a separate condition.
    """
    out: List[Dict[str, float]] = []
    grid = -np.geomspace(1e-8, abs(float(lambda_min)), int(n_grid))[::-1]
    grid = np.concatenate([grid, [0.0]])
    vals = [pair.secular(float(x)) for x in grid]
    for i in range(len(grid) - 1):
        a, b = vals[i], vals[i + 1]
        if not (math.isfinite(a) and math.isfinite(b)) or a * b > 0.0:
            continue
        r = _bisect(pair.secular, float(grid[i]), float(grid[i + 1]))
        if r is None or r >= 0.0:
            continue
        out.append({"lmbda": float(r), "sigma": math.sqrt(-r),
                    "growth_rate": math.sqrt(-r),
                    "secular": float(pair.secular(r))})
    return sorted(out, key=lambda r: r["lmbda"])


def is_stable(pair: MouthPair) -> Dict[str, object]:
    """Is the boundary data in the positive-``λ`` region?

    For an exchange-symmetric pair the answer is closed-form — both channel
    functions are monotone decreasing along the imaginary axis from their
    ``λ = 0`` values, so a growing mode exists iff ``α ± β`` falls below the
    corresponding threshold.  Otherwise the negative-``λ`` axis is scanned.  For
    symmetric data both are computed and compared, which is how the closed form
    is checked rather than trusted.
    """
    th = stability_thresholds(pair.separation)
    found = negative_lambda_modes(pair)
    out: Dict[str, object] = {
        "n_growing_modes": len(found), "modes": found,
        "stable": bool(not found),
        "worst_growth_rate": (max(m["growth_rate"] for m in found)
                              if found else 0.0),
    }
    if pair.is_exchange_symmetric():
        a, b = pair.alpha1, float(np.real(pair.beta))
        sym_bad = (a + b) < th["symmetric_threshold"]
        anti_bad = (a - b) < th["antisymmetric_threshold"]
        out.update({
            "closed_form_available": True,
            "symmetric_channel_unstable": bool(sym_bad),
            "antisymmetric_channel_unstable": bool(anti_bad),
            "closed_form_stable": bool(not (sym_bad or anti_bad)),
            "closed_form_agrees_with_the_scan": bool(
                (not (sym_bad or anti_bad)) == (not found)),
            **th})
    else:
        out["closed_form_available"] = False
    return out


def stability_map(separation: float = 1.3,
                  alphas: Sequence[float] = tuple(np.linspace(-0.15, 0.15, 13)),
                  betas: Sequence[float] = tuple(np.linspace(-0.2, 0.2, 17))
                  ) -> Dict[str, object]:
    """The stable region of the exchange-symmetric family, as a grid.

    The predicted region is the wedge ``α + β ≥ g₀ + G₀``, ``α − β ≥ g₀ − G₀``.
    Every grid point is also scanned for a negative-``λ`` root, so the wedge is
    verified point by point rather than drawn from the formula.
    """
    th = stability_thresholds(separation)
    grid = []
    mismatches = 0
    for a in alphas:
        row = []
        for b in betas:
            p = MouthPair(separation, float(a), float(a), float(b))
            pred = ((a + b) >= th["symmetric_threshold"]
                    and (a - b) >= th["antisymmetric_threshold"])
            got = not negative_lambda_modes(p)
            mismatches += int(pred != got)
            row.append(bool(got))
        grid.append(row)
    n = len(alphas) * len(betas)
    return {"alphas": [float(a) for a in alphas],
            "betas": [float(b) for b in betas],
            "stable": grid, "n_points": n,
            "n_stable": sum(sum(1 for x in r if x) for r in grid),
            "n_mismatches": mismatches,
            "the_closed_form_matches_everywhere": bool(mismatches == 0),
            **th}


def channel_endpoints(pair: MouthPair, gap: int) -> Dict[str, object]:
    """What each channel function actually does at the ends of a ``λ`` gap.

    The convenient statement "both run from ``−∞`` to ``+∞`` across every gap"
    is **false at the first one**.  The constant ``n = 0`` mode has the same
    value at both mouths, so it couples only to the symmetric combination and
    the antisymmetric channel's pole at ``ω = 1`` cancels: its endpoints there
    are finite, and whether a root exists in that gap depends on ``α − β``.
    """
    lo, hi = float(gap) ** 2, float(gap + 1) ** 2
    eps = 1e-9
    a = pair.channel_functions(lo + eps * max(1.0, lo))
    b = pair.channel_functions(hi - eps * max(1.0, hi))
    return {"gap": gap, "lambda_range": (lo, hi),
            "symmetric_at_lower": a[0], "symmetric_at_upper": b[0],
            "antisymmetric_at_lower": a[1], "antisymmetric_at_upper": b[1],
            "symmetric_spans_the_whole_line": bool(abs(a[0]) > 1e3
                                                   and abs(b[0]) > 1e3),
            "antisymmetric_spans_the_whole_line": bool(abs(a[1]) > 1e3
                                                       and abs(b[1]) > 1e3)}


def mouth_active_spectrum(pair: MouthPair, n_gaps: int = 6,
                          n_grid: int = 600) -> List[Dict[str, object]]:
    """Every mouth-active eigenvalue up to ``λ = (n_gaps+1)²``, by sector.

    Three sectors are scanned, and the first two are the ones an ``ω``-scan
    starting above ``1`` misses entirely: ``λ < 0`` (growing), ``0 ≤ λ < 1``
    (stable, **below the free ground state**), and each gap
    ``(m², (m+1)²)``.  The grid inside a gap is Chebyshev-clustered at both
    edges, because a strongly decoupled throat pushes both roots to within
    ``~1/‖A‖`` of an endpoint and a uniform grid then reports no root at all.
    """
    out: List[Dict[str, object]] = []
    for m in negative_lambda_modes(pair):
        out.append({"lmbda": m["lmbda"], "omega": None, "sector": "growing",
                    "growth_rate": m["growth_rate"]})

    def scan(lo: float, hi: float, sector: str, gap: Optional[int]) -> None:
        t = np.linspace(0.0, 1.0, int(n_grid))
        xs = lo + (hi - lo) * 0.5 * (1.0 - np.cos(math.pi * t))
        vals = []
        for x in xs:
            try:
                vals.append(pair.secular(float(x)))
            except ValueError:
                vals.append(math.nan)
        for i in range(len(xs) - 1):
            a, b = vals[i], vals[i + 1]
            if not (math.isfinite(a) and math.isfinite(b)) or a * b > 0.0:
                continue
            r = _bisect(pair.secular, float(xs[i]), float(xs[i + 1]))
            if r is None:
                continue
            out.append({"lmbda": float(r), "omega": math.sqrt(max(r, 0.0)),
                        "sector": sector, "gap": gap,
                        "secular": float(pair.secular(r))})

    scan(1e-9, 1.0 - 1e-9, "below the free ground state", 0)
    for m in range(1, int(n_gaps) + 1):
        scan(float(m) ** 2 + 1e-9, float(m + 1) ** 2 - 1e-9, "interlacing", m)
    return sorted(out, key=lambda r: r["lmbda"])


def spectrum_by_channel(pair: MouthPair, n_gaps: int = 8) -> Dict[str, object]:
    """For an exchange-symmetric pair, the two channels solved separately.

    ``g + G_d = α + β`` and ``g − G_d = α − β``.  Reported per gap **with the
    count**, because the antisymmetric channel does not always have a root: see
    `channel_endpoints`.
    """
    if not pair.is_exchange_symmetric():
        raise ValueError("channel splitting needs α₁ = α₂ and real β")
    a, b = pair.alpha1, float(np.real(pair.beta))
    rows = []
    for m in range(1, int(n_gaps) + 1):
        lo, hi = float(m) ** 2 + 1e-9, float(m + 1) ** 2 - 1e-9
        sym = _bisect(lambda x: pair.channel_functions(x)[0] - (a + b), lo, hi)
        anti = _bisect(lambda x: pair.channel_functions(x)[1] - (a - b), lo, hi)
        rows.append({"gap": m, "symmetric": sym, "antisymmetric": anti,
                     "n_roots": int(sym is not None) + int(anti is not None)})
    return {"rows": rows,
            "n_symmetric": sum(1 for r in rows if r["symmetric"] is not None),
            "n_antisymmetric": sum(1 for r in rows
                                   if r["antisymmetric"] is not None)}


def mode_charges(lmbda: float, pair: MouthPair) -> np.ndarray:
    """The null vector ``q`` of ``M(λ)`` — how strongly each mouth radiates."""
    _, _, vh = np.linalg.svd(pair.krein_matrix(lmbda))
    q = vh[-1].conjugate()
    return q / np.linalg.norm(q)


def mouth_flux(q: np.ndarray, boundary: np.ndarray) -> Dict[str, object]:
    """Flux absorbed at each mouth, ``Im(q_j* (A q)_j)``, and the net.

    Around mouth ``j`` the field is ``q_j/(4πχ) + φ_j^reg``, and the radial
    current integrated over a small sphere is ``Im(q_j* φ_j^reg)``, independent
    of the sphere's radius.  With ``φ^reg = A q`` the total is ``Im(q† A q)``,
    zero for Hermitian ``A`` **identically**.  This — not the Cayley entries —
    is the physical conservation statement.
    """
    q = np.asarray(q, dtype=complex)
    aq = np.asarray(boundary, dtype=complex) @ q
    per = [float(np.imag(q[j].conjugate() * aq[j])) for j in range(len(q))]
    return {"per_mouth": per, "net": float(sum(per)),
            "scale": float(np.abs(q).max() * np.abs(aq).max() + 1e-300)}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_green_function_has_a_finite_part(
        omegas: Sequence[float] = (0.37, 1.63, 2.4, 5.21),
        chis: Sequence[float] = (0.4, 1.3, 2.6),
        radii: Sequence[float] = (1e-2, 1e-3, 1e-4)) -> Dict[str, object]:
    """The closed form, checked against PR #255, and its short-distance split.

    ``G`` is the ``γ → 0`` limit of PR #255's branch series — an independent
    construction — and ``G(χ,ω) − 1/(4πχ) → g(ω)`` **linearly in χ**, which is
    what makes a point interaction definable: the divergence is the universal
    Coulomb one, so the subtraction is forced rather than chosen.
    """
    from geometrodynamics.waves.branch_coupling import mouth_transfer

    rows = []
    worst = 0.0
    for w in omegas:
        for chi in chis:
            got = mouth_transfer(w, chi, 1e-6)
            want = free_green(chi, w)
            err = abs(got.real - want)
            worst = max(worst, err)
            rows.append({"omega": w, "chi": chi, "closed_form": want,
                         "branch_series_limit": got.real, "abs_error": err})

    fin, ratios = [], []
    for w in omegas[:3]:
        errs = []
        for r in radii:
            fp = free_green(float(r), w) - 1.0 / (4.0 * math.pi * r)
            errs.append(abs(fp - regularized_green(w)))
            fin.append({"omega": w, "chi": r, "finite_part": fp,
                        "g": regularized_green(w), "error": errs[-1]})
        ratios += [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]

    # and the antipode, where the numerator zero cancels sin χ
    anti = []
    for w in (0.7, 2.3, 4.4):
        want = w / (4.0 * math.pi * math.sin(math.pi * w))
        anti.append({"omega": w, "limit": want,
                     "at_pi_minus_1e-5": free_green(math.pi - 1e-5, w),
                     "relative_error": abs(free_green(math.pi - 1e-5, w)
                                           - want) / abs(want)})
    return {
        "rows": rows, "worst_abs_error": worst,
        "the_closed_form_is_the_branch_series": bool(worst < 1e-10),
        "finite_part": fin, "convergence_ratios_per_decade": ratios,
        "the_remainder_is_first_order_in_chi": bool(
            ratios and all(abs(r - 10.0) < 1.0 for r in ratios)),
        "antipode": anti,
        "the_antipodal_focus_is_finite": bool(
            all(a["relative_error"] < 1e-8 for a in anti)),
        "what_this_shows": ("the coincidence divergence is the universal "
                            "1/(4πχ), so g(ω) exists and a point interaction "
                            "is definable"),
    }


def measure_hermiticity_is_flux_conservation(
        separation: float = 1.3, n_draws: int = 200, seed: int = 20260815,
        scales: Sequence[float] = (0.05, 0.1, 0.2),
        omega: float = 2.3, kappa: float = 0.3) -> Dict[str, object]:
    """``Im(q† A q) = 0`` ⟺ ``A = A†`` — the physical conservation statement.

    The radial current through a small sphere at mouth ``j`` is
    ``Im(q_j* φ_j^reg)``, independent of the sphere, so the total absorbed is
    ``Im(q† A q)``.  For Hermitian ``A`` that vanishes for **every** ``q``, and
    a purely off-diagonal throat moves flux from one mouth to the other exactly.

    Reported alongside — and this is a **correction to an earlier version of
    this module** — the Cayley entries are *not* reflection and transmission
    amplitudes.  Their magnitudes depend on the arbitrary reference scale ``c``:
    the same ``A`` gives quite different numbers at different ``c``, tabulated
    here.  Unitarity of the Cayley transform is a real fact about the
    parametrization; reading its entries as a scattering matrix is not, because
    a closed universe supplies no asymptotic flux normalization.
    """
    rng = np.random.default_rng(seed)
    worst, worst_pair = 0.0, 0.0
    for _ in range(int(n_draws)):
        a1, a2 = rng.normal(0, 0.3, 2)
        b = complex(*rng.normal(0, 0.3, 2))
        p = MouthPair(separation, float(a1), float(a2), b)
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        f = mouth_flux(q, p.boundary_matrix())
        worst = max(worst, abs(f["net"]) / f["scale"])
        f0 = mouth_flux(q, MouthPair(separation, 0.0, 0.0, b).boundary_matrix())
        worst_pair = max(worst_pair,
                         abs(f0["per_mouth"][0] + f0["per_mouth"][1])
                         / f0["scale"])

    # the non-self-adjoint control, in the same units
    ctl = DirectionalThroat(separation, 1.0, +1, kappa)
    a_ctl = np.array([[0.0, 0.0], [1.0 / ctl.gain(omega), 0.0]], dtype=complex)
    nets = []
    for _ in range(32):
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        nets.append(abs(mouth_flux(q, a_ctl)["net"])
                    / mouth_flux(q, a_ctl)["scale"])

    # the Cayley caveat, as numbers
    p = MouthPair(separation, 0.2, -0.13, 0.15 + 0.07j)
    cay = []
    for c in scales:
        m = boundary_mixing(p.boundary_matrix(), float(c))
        cay.append({"scale": c, "diagonal": m["diagonal_mixing"][0],
                    "off_diagonal": m["off_diagonal_mixing"][0],
                    "unitarity_defect": m["unitarity_defect"],
                    "column_norm": m["column_norms"][0]})
    spread = (max(r["diagonal"] for r in cay)
              - min(r["diagonal"] for r in cay))
    return {
        "n_draws": int(n_draws), "worst_relative_net_flux": worst,
        "flux_is_conserved_identically": bool(worst < 1e-14),
        "worst_pairwise_imbalance": worst_pair,
        "what_one_mouth_absorbs_the_other_emits": bool(worst_pair < 1e-14),
        "control_median_net_flux": float(np.median(nets)),
        "the_control_does_not_conserve": bool(np.median(nets) > 1e-3),
        "cayley": cay,
        "the_cayley_transform_is_unitary": bool(
            all(r["unitarity_defect"] < 1e-13 for r in cay)),
        "cayley_diagonal_spread_over_the_reference_scale": spread,
        "the_cayley_entries_are_not_physical_amplitudes": bool(spread > 0.1),
        "the_correction": ("an earlier version of this module read the Cayley "
                           "entries as reflection and transmission; their "
                           "magnitudes depend on the arbitrary scale c, so "
                           "they are boundary-mixing coefficients and the "
                           "physical conservation result is the flux identity"),
        "the_identity": "Im(q† A q) = 0 for all q  ⟺  A = A†",
    }


def measure_self_adjointness_makes_lambda_real_not_omega(
        separation: float = 1.3,
        params: Sequence[Tuple[float, float, complex]] = (
            (0.05, 0.05, 0.03), (0.2, -0.13, 0.15 + 0.07j),
            (-0.4, 0.07, -0.09 + 0.31j)),
        n_gaps: int = 6) -> Dict[str, object]:
    """**The correction this round needed.**  Hermiticity buys ``λ`` real, and
    that is all it buys.

    ``Γ(λ)`` is real symmetric for real ``λ`` of *either sign*, so ``M = A − Γ``
    is Hermitian there and the secular function is real: the eigenvalues of the
    spatial operator are real.  They are not thereby **positive**, and a
    negative ``λ`` means ``ω = ±i√|λ|`` with one member of the pair growing like
    ``e^{√|λ| t}``.

    Two of the three boundary matrices this module previously advertised have
    exactly that.  They were missed because the earlier search seeded only
    ``Re ω ∈ [1.1, 6.9]`` and discarded roots leaving that window — a search
    that structurally could not see a root on the imaginary axis.  The claim
    "real spectrum for every coupling, so a conserving throat cannot ring up" is
    **withdrawn**; what replaces it is the stability region.
    """
    rows = []
    worst_im = 0.0
    for (a1, a2, b) in params:
        p = MouthPair(separation, a1, a2, b)
        im = 0.0
        for lam in np.linspace(-40.0, float(n_gaps + 1) ** 2 - 0.3, 500):
            try:
                det = complex(np.linalg.det(p.krein_matrix(float(lam))))
            except ValueError:
                continue
            im = max(im, abs(det.imag) / max(abs(det), 1e-300))
        st = is_stable(p)
        worst_im = max(worst_im, im)
        rows.append({
            "alpha1": a1, "alpha2": a2, "beta": b,
            "hermitian": p.is_self_adjoint(),
            "worst_relative_imaginary_part_of_det": im,
            "n_growing_modes": st["n_growing_modes"],
            "lambda_of_the_growing_modes": [m["lmbda"] for m in st["modes"]],
            "growth_rates": [m["growth_rate"] for m in st["modes"]],
            "stable": st["stable"]})

    unstable = [r for r in rows if not r["stable"]]
    return {
        "rows": rows,
        "worst_relative_imaginary_part": worst_im,
        "the_secular_function_is_real_in_lambda": bool(worst_im < 1e-12),
        "hermiticity_gives_real_lambda": True,
        "n_unstable_examples": len(unstable),
        "hermiticity_does_not_give_positivity": bool(len(unstable) > 0),
        "the_withdrawn_claim": ("'the spectrum is real for every coupling, so "
                                "a conserving throat cannot ring up' — false: "
                                "λ is real, ω need not be"),
        "why_it_was_missed": ("the earlier complex-root search seeded only "
                              "Re ω in [1.1, 6.9] and discarded roots leaving "
                              "that window, so a root on the imaginary axis "
                              "was outside its reach by construction"),
        "what_replaces_it": ("a positivity audit: solve in λ, and map the "
                             "region of boundary data with λ_min ≥ 0"),
    }


def measure_the_stability_region_in_the_boundary_family(
        separation: float = 1.3,
        probes: Sequence[Tuple[float, float]] = ((0.05, 0.03), (0.0, 0.0),
                                                 (-0.05, 0.0), (0.05, 0.30),
                                                 (0.2, 0.0)),
        n_alpha: int = 13, n_beta: int = 17) -> Dict[str, object]:
    """Where the family **is** stable — in closed form, and verified pointwise.

    Along the imaginary axis both channel functions are monotone decreasing from
    their ``λ = 0`` values ``g₀ ± G₀`` to ``−∞``, so ``g ± G_d = α ± β`` has a
    negative-``λ`` root iff the right-hand side falls below the threshold.  The
    stable set is therefore the wedge

        ``α + β ≥ g₀ + G₀``   and   ``α − β ≥ g₀ − G₀``

    with ``g₀ = −1/(4π²)`` and ``G₀ = (π−d)/(4π² sin d)``.  Every point of a
    grid is *also* scanned for a negative-``λ`` root, so the wedge is checked
    rather than asserted — and the monotonicity it rests on is checked too.
    """
    th = stability_thresholds(separation)

    # ordered by increasing σ, i.e. decreasing λ, both channels must fall
    mono = True
    lams = -np.geomspace(1e-6, 1600.0, 4000)
    p0 = MouthPair(separation, 0.0, 0.0, 0.0)
    ends = {}
    for k in (0, 1):
        v = np.array([p0.channel_functions(float(x))[k] for x in lams])
        mono = mono and bool((np.diff(v) < 0).all())
        ends[k] = (float(v[0]), float(v[-1]))

    rows = []
    for (a, b) in probes:
        p = MouthPair(separation, a, a, b)
        st = is_stable(p)
        rows.append({"alpha": a, "beta": b, "stable": st["stable"],
                     "closed_form_stable": st["closed_form_stable"],
                     "agrees": st["closed_form_agrees_with_the_scan"],
                     "n_growing_modes": st["n_growing_modes"],
                     "worst_growth_rate": st["worst_growth_rate"]})

    grid = stability_map(separation,
                         tuple(np.linspace(-0.15, 0.15, int(n_alpha))),
                         tuple(np.linspace(-0.2, 0.2, int(n_beta))))
    return {
        "thresholds": th, "probes": rows,
        "the_channel_functions_are_monotone_on_the_imaginary_axis": mono,
        "channel_run_along_the_imaginary_axis": {
            "symmetric": ends[0], "antisymmetric": ends[1]},
        "the_closed_form_agrees_with_every_probe": bool(
            all(r["agrees"] for r in rows)),
        "grid_points": grid["n_points"], "grid_stable": grid["n_stable"],
        "grid_mismatches": grid["n_mismatches"],
        "the_closed_form_matches_everywhere": grid[
            "the_closed_form_matches_everywhere"],
        "both_signs_are_represented": bool(
            0 < grid["n_stable"] < grid["n_points"]),
        "what_this_shows": ("positivity is a condition on the boundary data, "
                            "separate from self-adjointness, and it has a "
                            "closed form for the exchange-symmetric family"),
    }


def measure_the_mouth_active_sector_is_rank_two(
        separation: float = 1.3, alpha: float = 0.05, beta: float = 0.03,
        n_gaps: int = 6, levels: Sequence[int] = (0, 1, 2, 3, 4)
        ) -> Dict[str, object]:
    """What ``det(C − BΓ) = 0`` describes — and what it leaves out.

    A two-point interaction is a **rank-two** perturbation.  Level ``n`` has
    degeneracy ``(n+1)²`` and only two combinations can move, so ``(n+1)² − 2``
    modes stay exactly at the free eigenvalue and the secular equation never
    sees them.  Scoping matters: "two roots per gap" is a statement about the
    mouth-active sector, not about the spectrum.

    Also measured, because the convenient version is false: the claim that both
    channel functions run ``−∞ → +∞`` across **every** gap.  The ``n = 0``
    constant mode has the same value at both mouths, so it couples only to the
    symmetric channel; the antisymmetric channel's pole at ``ω = 1`` cancels,
    its endpoints in the first gap are finite, and a root there is conditional
    on ``α − β``.
    """
    pair = MouthPair(separation, alpha, alpha, beta)
    spec = mouth_active_spectrum(pair, n_gaps)
    by_ch = spectrum_by_channel(pair, n_gaps)

    ends = [channel_endpoints(pair, m) for m in range(1, n_gaps + 1)]
    first = ends[0]

    counts = {}
    for r in spec:
        counts[r["sector"]] = counts.get(r["sector"], 0) + 1
    per_gap = {}
    for r in spec:
        if r["sector"] == "interlacing":
            per_gap[r["gap"]] = per_gap.get(r["gap"], 0) + 1

    book = [{"level": n, "degeneracy": (n + 1) ** 2,
             "mouth_active": min(2, (n + 1) ** 2),
             "untouched": untouched_multiplicity(n)} for n in levels]

    # does the antisymmetric channel have a root in gap 1, and is that a choice?
    scan = []
    for ab in (0.02, -0.05, -0.09):
        p = MouthPair(separation, alpha, alpha, alpha - ab)
        r = spectrum_by_channel(p, 1)["rows"][0]
        scan.append({"alpha_minus_beta": ab,
                     "antisymmetric_root_in_gap_1": r["antisymmetric"]
                     is not None})

    return {
        "sector_counts": counts, "roots_per_interlacing_gap": per_gap,
        "two_per_interlacing_gap": bool(per_gap
                                        and set(per_gap.values()) == {2}),
        "n_below_the_free_ground_state": counts.get(
            "below the free ground state", 0),
        "there_is_a_sector_below_the_ground_state": bool(
            counts.get("below the free ground state", 0) > 0),
        "level_bookkeeping": book,
        "at_most_two_modes_move_per_level": bool(
            all(r["mouth_active"] <= 2 for r in book)),
        "untouched_modes_at_level_4": untouched_multiplicity(4),
        "channel_endpoints": ends,
        "the_first_gap_antisymmetric_endpoints_are_finite": bool(
            not first["antisymmetric_spans_the_whole_line"]),
        "the_symmetric_channel_does_span_it": bool(
            first["symmetric_spans_the_whole_line"]),
        "first_gap_root_depends_on_alpha_minus_beta": scan,
        "existence_in_the_first_gap_is_conditional": bool(
            len({r["antisymmetric_root_in_gap_1"] for r in scan}) > 1),
        "by_channel": by_ch["rows"],
        "the_scope": ("det(C − BΓ) = 0 is the rank-two mouth-active sector; "
                      "(n+1)² − 2 modes per level are untouched and are not in "
                      "it"),
    }


def measure_the_pr255_boundary_condition_is_not_self_adjoint(
        separation: float = 1.3, delay: float = 1.0,
        kappas: Sequence[float] = (0.3, 1.0),
        omegas: Sequence[complex] = (1.3 + 0.2j, 2.77 - 0.4j, 4.11 + 0.05j)
        ) -> Dict[str, object]:
    """PR #255's relation, embedded **exactly** — and a correction on the way.

    That round set ``q₁ = 0`` and ``q₂ = gain · φ₁^reg``: mouth ``M⁺`` absorbs
    but does not radiate.  In the ``B φ^reg = C q`` form that is
    ``B = [[0,0],[gain,0]]``, ``C = I``, for which ``det(C − BΓ) = 1 − gain·G_d``
    — **exactly** PR #255's own ``1 − L``, checked here against that round's
    expression rather than re-derived twice.  It is maximal (``rank[B|C] = 2``)
    but ``BC† = B`` is not Hermitian, so it is not self-adjoint; and no finite
    Hermitian ``A`` reproduces it, because it needs the singular ``B``.

    **The correction:** an earlier version of this module used
    ``A = [[0,0],[1/gain,0]]`` with ``B = I`` as the control.  That gives
    ``det(A − Γ) = g² − G_d² + G_d/gain``, a different function — so the old
    control was not PR #255's model, and any conclusion drawn by comparing
    against it was unsupported.  Both are computed here so the difference is a
    number rather than a claim.
    """
    rows = []
    worst_embed = 0.0
    worst_old = 0.0
    for kap in kappas:
        ctl = DirectionalThroat(separation, delay, +1, kap)
        for w in omegas:
            embed = ctl.secular(w)
            own = ctl.pr255_pole_condition(w)
            worst_embed = max(worst_embed, abs(embed - own))
            # the old, wrong control
            sp = cmath.sin(math.pi * w)
            g = -w * cmath.cos(math.pi * w) / (4.0 * math.pi * sp)
            gd = (cmath.sin(w * (math.pi - separation))
                  / (4.0 * math.pi * math.sin(separation) * sp))
            a_old = np.array([[0.0, 0.0], [1.0 / ctl.gain(w), 0.0]],
                             dtype=complex)
            old = complex(np.linalg.det(
                a_old - np.array([[g, gd], [gd, g]], dtype=complex)))
            worst_old = max(worst_old, abs(old - own))
            bc = ctl.boundary_condition(w)
            rows.append({"kappa": kap, "omega": w,
                         "det_of_the_embedding": embed,
                         "pr255_one_minus_L": own,
                         "embedding_error": abs(embed - own),
                         "old_wrong_control": old,
                         "old_control_error": abs(old - own),
                         "maximal": bc.is_maximal(),
                         "self_adjointness_defect":
                             bc.self_adjointness_defect()})
    return {
        "rows": rows,
        "worst_embedding_error": worst_embed,
        "the_embedding_is_exact": bool(worst_embed < 1e-12),
        "worst_old_control_error": worst_old,
        "the_old_control_was_a_different_model": bool(worst_old > 1e-3),
        "every_embedding_is_maximal": bool(all(r["maximal"] for r in rows)),
        "none_is_self_adjoint": bool(
            all(r["self_adjointness_defect"] > 1e-9 for r in rows)),
        "it_needs_a_singular_B": ("q₁ = 0 is not of the form φ^reg = A q for "
                                  "any finite Hermitian A"),
        "what_is_not_concluded": ("that PR #255's off-axis poles were *caused* "
                                  "by non-self-adjointness — a self-adjoint "
                                  "throat can also have growing modes, as the "
                                  "stability measurement shows, so this is a "
                                  "classification, not a diagnosis"),
    }


def measure_the_spectrum_is_conjugation_symmetric_in_beta(
        separation: float = 1.3, alpha1: float = 0.2, alpha2: float = -0.13,
        modulus: float = 0.1655, n_phase: int = 24,
        lmbda: float = 5.3) -> Dict[str, object]:
    """What the phase of ``β`` does, and what it does not.

    The secular function is ``(α₁−g)(α₂−g) − |β − G_d|²`` with ``G_d`` real, so
    it depends on ``Re β`` as well as ``|β|`` — the mouths are connected through
    the bulk too, and that fixes the relative phase.  So ``arg β`` **is**
    physical here.

    What is *not* physical is a preferred direction: the secular function is
    invariant under ``β → β*``, which is time reversal.  An earlier version of
    this module called a complex ``β`` "non-reciprocal", reading the asymmetry
    of the Cayley entries as a physical asymmetry.  It is not: those entries
    depend on the reference scale, and the spectrum itself is conjugation
    symmetric.
    """
    gam = gamma_at(lmbda, separation)
    gd = float(gam[0, 1].real)
    g = float(gam[0, 0].real)

    def det_of(b: complex) -> float:
        return (alpha1 - g) * (alpha2 - g) - abs(b - gd) ** 2

    rows = []
    for th in np.linspace(0.0, math.pi, int(n_phase)):
        b = modulus * cmath.exp(1j * float(th))
        rows.append({"arg_beta": float(th), "det": det_of(b),
                     "det_conjugate": det_of(b.conjugate()),
                     "conjugation_defect": abs(det_of(b)
                                               - det_of(b.conjugate()))})
    spread = max(r["det"] for r in rows) - min(r["det"] for r in rows)
    return {
        "rows": rows, "modulus": modulus,
        "spread_over_the_phase": spread,
        "the_phase_of_beta_is_physical": bool(spread > 1e-6),
        "worst_conjugation_defect": max(r["conjugation_defect"]
                                        for r in rows),
        "the_spectrum_is_conjugation_symmetric": bool(
            max(r["conjugation_defect"] for r in rows) < 1e-15),
        "why_the_phase_survives": ("Γ is not diagonal — the mouths are joined "
                                   "through the bulk as well as through the "
                                   "throat, and that fixes the relative phase"),
        "the_withdrawn_claim": ("that a complex β makes the throat "
                                "non-reciprocal; that read the Cayley entries, "
                                "which depend on an arbitrary scale, as "
                                "physical amplitudes"),
    }
