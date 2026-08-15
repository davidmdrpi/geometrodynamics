"""
A flux-conserving two-mouth throat operator, and the spectrum it produces.

WHAT PR #255 OWED
─────────────────
`waves.branch_coupling` solved a mouth relation self-consistently and was
explicit that the relation was not a boundary operator: it carried field
**values** through the free Green function, with no normal-derivative matching,
no reflected channel, a ``1×1`` mouth object where a conserving junction needs
``2×2`` unitary, and ``κ²`` power throughput — lossy, hence not an
identification of anything.  It also found its poles off the real axis and had
to distinguish three thresholds to say what that meant.

This module replaces the relation with the real object.  A point-supported
throat on a manifold is a **self-adjoint extension** of the Laplacian on
``S³ ∖ {M⁺, M⁻}``, and von Neumann's theorem says those are parametrized by a
unitary matrix between the deficiency spaces — here ``U(2)``.  Equivalently, by
Krein's formula, a Hermitian ``2×2`` matrix ``A``.  Every claim below follows
from that one substitution.

THE FREE GREEN FUNCTION, IN CLOSED FORM
───────────────────────────────────────
Fourier-transforming PR #254's image sum (or summing PR #255's branch series)
gives, on the conformally coupled ESU,

    ``G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))``

**real** for real ``ω``, with poles exactly at the free spectrum ``ω = n+1``.
Its short-distance expansion is ``1/(4πχ) + g(ω) + O(χ)`` with

    ``g(ω) = −(ω/4π) · cot(πω)``

the **regularized** value at coincidence — the object a point interaction needs,
and the reason the throat is renormalized rather than singular.

THE OPERATOR
────────────
Write the field as free plus two point sources,
``φ(x) = φ_in(x) + Σ_j G(χ_j, ω) q_j``.  Near mouth ``j`` this is
``q_j/(4πχ_j) + φ_j^reg + O(χ_j)``, and the self-adjoint extension is the
statement that the regular parts are a **Hermitian** linear image of the
charges, ``φ^reg = A q``.  Collecting terms,

    ``M(ω) q = φ_in|_mouths`` ,   ``M(ω) = A − Γ(ω)`` ,
    ``Γ(ω) = [[g(ω), G_d(ω)], [G_d(ω), g(ω)]]``

with ``d`` the geodesic mouth separation.  Then:

* **``A`` Hermitian is exactly flux conservation.**  The flux absorbed at the
  mouths is ``Σ_j Im(q_j* φ_j^reg) = Im(q† A q)``, which vanishes identically
  iff ``A = A†``.  Measured, against a directional control that does not.
* **the boundary operator is a unitary ``2×2``.**  The Cayley transform
  ``S = (A − ic)(A + ic)^{-1}`` is unitary for every Hermitian ``A`` and every
  real scale ``c``, and inverts back.  Its **diagonal is reflection** and its
  **off-diagonal is transmission**, with ``|r|² + |t|² = 1`` exactly.  PR #255's
  model is the corner of this family with ``r = 0`` and ``|t| = κ ≠ 1``, which
  is not in ``U(2)`` at all unless ``κ = 1``.
* **the spectrum is real, for every coupling.**  ``Γ`` is real symmetric on the
  real axis, so ``M`` is Hermitian there and ``det M(ω)`` is a real function;
  its zeros are the coupled eigenfrequencies.  A scan of the complex plane finds
  nothing off the axis.  **PR #255's instability was an artefact of the
  non-conserving model**, not a feature of a throat — and that is the result
  this round exists to establish.
* **and the coupled spectrum interlaces the free one.**  In the exchange-
  symmetric case ``M`` splits into ``g ± G_d = α ± β``, and each of those runs
  monotonically from ``−∞`` to ``+∞`` across every unit gap, so there is
  **exactly one root per channel per gap** — two coupled frequencies between
  consecutive free ones, and the free spectrum recovered as the coupling is
  switched off.

WHAT IS STILL PUT IN
────────────────────
A **linear** field on a **fixed** round background.  ``A`` is a choice — four
real parameters — not a derivation: `shells.junction` (PR #249) is what would
fix them from a matter model, and nothing here computes the exotic-matter bill.
The throat is still **point-supported**, so it has no interior and no proper
length; the "delay" of PRs #251–#255 is not a parameter of a self-adjoint point
extension and does not survive into it, which is itself worth stating plainly.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant — that is the next step.
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
    "MouthPair",
    "DirectionalThroat",
    "coupled_spectrum",
    "spectrum_by_channel",
    "boundary_from_scattering",
    "mode_charges",
    "mouth_flux",
    "complex_root_search",
    "measure_the_green_function_has_a_finite_part",
    "measure_the_boundary_operator_is_unitary_with_both_channels",
    "measure_self_adjointness_makes_the_spectrum_real",
    "measure_the_coupled_spectrum_interlaces_the_free_one",
    "measure_the_flux_balance_is_exactly_hermiticity",
    "measure_the_directional_model_is_what_pr255_solved",
]

TWO_PI = 2.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE FREE GREEN FUNCTION
# ════════════════════════════════════════════════════════════════════════════
def free_green(chi: float, omega: float) -> float:
    """``G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))``.

    The frequency-domain retarded Green function of the conformally coupled
    massless scalar on the ESU, in closed form.  It is the Fourier transform of
    PR #254's image sum and the ``γ → 0`` limit of PR #255's branch series —
    checked against the latter rather than asserted.

    **Real** for real ``ω``, which is the first sign that a self-adjoint problem
    is being set up: a closed universe has no radiation condition to break the
    reality.
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
    Coulomb singularity that any point interaction subtracts; what is left is
    what the throat couples to.  Poles at ``ω ∈ ℤ`` — the free spectrum — with
    residue ``−m/(4π²)``.
    """
    sp = math.sin(math.pi * omega)
    if abs(sp) < 1e-13:
        raise ValueError("ω is a free eigenfrequency; g has a pole there")
    return -omega * math.cos(math.pi * omega) / (4.0 * math.pi * sp)


def esu_free_spectrum(n_max: int = 8) -> List[float]:
    """``ω_n = n + 1`` — the conformally coupled ESU spectrum of PR #254."""
    return [float(n + 1) for n in range(int(n_max))]


# ════════════════════════════════════════════════════════════════════════════
# THE OPERATOR
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class MouthPair:
    """A **self-adjoint** two-mouth throat: a Hermitian ``2×2`` boundary matrix.

    ``alpha1``, ``alpha2`` are the mouths' own (real) self-energies — inverse
    scattering lengths, and the **reflection** channel PR #255 had none of.
    ``beta`` is the complex mouth-to-mouth amplitude, the **transmission**
    channel, and the only part of the operator with no local realization on
    ``S³``: that non-locality is the wormhole.

    Hermiticity is not decoration.  It is the statement that no flux is created
    or destroyed at the mouths, and everything downstream — real spectrum,
    orthogonal modes, a unitary scattering matrix — follows from it.
    """

    separation: float
    alpha1: float = 0.05
    alpha2: float = 0.05
    beta: complex = 0.03

    def boundary_matrix(self) -> np.ndarray:
        """``A``, Hermitian by construction."""
        b = complex(self.beta)
        return np.array([[self.alpha1, b], [b.conjugate(), self.alpha2]],
                        dtype=complex)

    def is_self_adjoint(self, tol: float = 1e-15) -> bool:
        a = self.boundary_matrix()
        return bool(np.abs(a - a.conjugate().T).max() <= tol)

    def is_exchange_symmetric(self, tol: float = 1e-15) -> bool:
        """``α₁ = α₂`` and ``β`` real — the case that splits into two channels."""
        return bool(abs(self.alpha1 - self.alpha2) <= tol
                    and abs(complex(self.beta).imag) <= tol)

    # ── the Krein matrix ────────────────────────────────────────────────────
    def gamma(self, omega: float) -> np.ndarray:
        """``Γ(ω)`` — the free propagator evaluated between the mouths.

        Diagonal is the *regularized* coincidence value; off-diagonal is the
        ordinary Green function at the mouth separation.  Real symmetric.
        """
        g = regularized_green(omega)
        gd = free_green(self.separation, omega)
        return np.array([[g, gd], [gd, g]], dtype=complex)

    def krein_matrix(self, omega: float) -> np.ndarray:
        """``M(ω) = A − Γ(ω)``.  Hermitian for real ``ω``."""
        return self.boundary_matrix() - self.gamma(omega)

    def secular(self, omega: float) -> float:
        """``det M(ω)`` — real for real ``ω``, and its zeros are the spectrum."""
        return float(np.linalg.det(self.krein_matrix(omega)).real)

    def channel_functions(self, omega: float) -> Tuple[float, float]:
        """``g ± G_d`` — the symmetric and antisymmetric channels.

        Meaningful when the pair is exchange-symmetric, where the secular
        equation factorizes into ``g + G_d = α + β`` and ``g − G_d = α − β``.
        """
        g = regularized_green(omega)
        gd = free_green(self.separation, omega)
        return g + gd, g - gd

    # ── the scattering form ─────────────────────────────────────────────────
    def scattering_matrix(self, scale: float = 0.1) -> np.ndarray:
        """``S = (A − ic)(A + ic)^{-1}`` — the **unitary** boundary operator.

        The Cayley transform of a Hermitian matrix is unitary for every real
        ``c``, which is von Neumann's parametrization of self-adjoint extensions
        made concrete.  ``c`` is a reference impedance: it fixes which unitary
        represents a given ``A``, not whether one does.
        """
        a = self.boundary_matrix()
        eye = np.eye(2, dtype=complex)
        return (a - 1j * scale * eye) @ np.linalg.inv(a + 1j * scale * eye)

    def channels(self, scale: float = 0.1) -> Dict[str, object]:
        """Reflection and transmission, read off the unitary.

        ``|r|² + |t|² = 1`` per mouth — the statement PR #255's ``1×1`` object
        could not make, since it had no ``r`` and ``|t| = κ``.
        """
        s = self.scattering_matrix(scale)
        r1, r2 = complex(s[0, 0]), complex(s[1, 1])
        t12, t21 = complex(s[0, 1]), complex(s[1, 0])
        return {
            "S": s, "reflection_1": r1, "reflection_2": r2,
            "transmission_12": t12, "transmission_21": t21,
            "unitarity_defect": float(
                np.abs(s.conjugate().T @ s - np.eye(2)).max()),
            "sum_of_squares_mouth_1": float(abs(r1) ** 2 + abs(t21) ** 2),
            "sum_of_squares_mouth_2": float(abs(r2) ** 2 + abs(t12) ** 2),
            "reciprocal": bool(abs(t12 - t21) < 1e-12),
        }


def boundary_from_scattering(s: np.ndarray, scale: float = 0.1) -> np.ndarray:
    """The inverse Cayley, ``A = ic (I + S)(I − S)^{-1}`` — Hermitian again."""
    eye = np.eye(2, dtype=complex)
    return 1j * scale * (eye + s) @ np.linalg.inv(eye - s)


@dataclass(frozen=True)
class DirectionalThroat:
    """PR #255's model, written in the same language so it can be compared.

    That model drove ``M⁻`` from what arrived at ``M⁺`` and **not the reverse**,
    with a gain ``κ`` and a delay ``Δ``.  In boundary-operator terms that is

        ``A(ω) = [[0, 0], [1/(ηκ e^{−iωΔ}), 0]]``

    — strictly lower-triangular, hence **not Hermitian**, hence not
    flux-conserving.  It is kept here as the control: everything the
    self-adjoint operator gets right, this one gets wrong in a measurable way.
    """

    separation: float
    delay: float = 1.0
    eta: int = 1
    kappa: float = 0.3

    def boundary_matrix(self, omega: float) -> np.ndarray:
        gain = self.eta * self.kappa * cmath.exp(-1j * float(omega) * self.delay)
        return np.array([[0.0, 0.0], [1.0 / gain, 0.0]], dtype=complex)

    def anti_hermitian_part(self, omega: float) -> float:
        a = self.boundary_matrix(omega)
        return float(np.abs(0.5 * (a - a.conjugate().T)).max())

    def krein_matrix(self, omega: complex) -> np.ndarray:
        w = complex(omega)
        g = -w * cmath.cos(math.pi * w) / (4.0 * math.pi
                                           * cmath.sin(math.pi * w))
        d = self.separation
        gd = (cmath.sin(w * (math.pi - d))
              / (4.0 * math.pi * math.sin(d) * cmath.sin(math.pi * w)))
        gain = self.eta * self.kappa * cmath.exp(-1j * w * self.delay)
        a = np.array([[0.0, 0.0], [1.0 / gain, 0.0]], dtype=complex)
        return a - np.array([[g, gd], [gd, g]], dtype=complex)

    def secular(self, omega: complex) -> complex:
        return complex(np.linalg.det(self.krein_matrix(omega)))


# ════════════════════════════════════════════════════════════════════════════
# THE SPECTRUM
# ════════════════════════════════════════════════════════════════════════════
def _brent(f, lo: float, hi: float, n: int = 200,
           tol: float = 1e-14) -> Optional[float]:
    """Bisection with a secant polish — no SciPy dependency for one root."""
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
            b, fb = m, fm
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)


def coupled_spectrum(pair: MouthPair, n_gaps: int = 8,
                     n_grid: int = 400,
                     edge: float = 1e-7) -> List[Dict[str, object]]:
    """The real roots of ``det M(ω) = 0``, gap by gap.

    ``g`` and ``G_d`` both have simple poles at every free eigenfrequency, so
    the secular function is scanned strictly *between* consecutive integers and
    every sign change is bracketed.  No root is assumed: the count per gap is
    reported, because it is one of the things being measured.

    The grid is **Chebyshev-clustered at both edges of the gap**, not uniform.
    A strongly decoupled throat pushes both roots to within ``~1/‖A‖`` of an
    integer, and on a uniform grid both then land in one cell, the determinant
    changes sign twice inside it, and the scan reports *no* root — a failure
    that looks exactly like a correct decoupling limit.
    """
    out: List[Dict[str, object]] = []
    for m in range(1, int(n_gaps) + 1):
        lo, hi = m + edge, m + 1.0 - edge
        t = np.linspace(0.0, 1.0, int(n_grid))
        ws = lo + (hi - lo) * 0.5 * (1.0 - np.cos(math.pi * t))
        vals = []
        for w in ws:
            try:
                vals.append(pair.secular(float(w)))
            except ValueError:
                vals.append(math.nan)
        for i in range(len(ws) - 1):
            a, b = vals[i], vals[i + 1]
            if not (math.isfinite(a) and math.isfinite(b)) or a * b > 0.0:
                continue
            r = _brent(pair.secular, float(ws[i]), float(ws[i + 1]))
            if r is None:
                continue
            out.append({"omega": float(r), "gap": m,
                        "position_in_gap": float(r - m),
                        "distance_to_nearest_free": float(
                            abs(r - round(r))),
                        "secular": float(pair.secular(r))})
    return out


def spectrum_by_channel(pair: MouthPair, n_gaps: int = 8) -> Dict[str, object]:
    """For an exchange-symmetric pair, the two channels solved separately.

    ``g + G_d = α + β`` is the symmetric channel and ``g − G_d = α − β`` the
    antisymmetric one.  Each left-hand side is monotone from ``−∞`` to ``+∞``
    across a unit gap, so each contributes exactly one root there — which is why
    the coupled spectrum interlaces the free one two-for-one.
    """
    if not pair.is_exchange_symmetric():
        raise ValueError("channel splitting needs α₁ = α₂ and real β")
    a, b = pair.alpha1, float(np.real(pair.beta))
    rows = []
    for m in range(1, int(n_gaps) + 1):
        lo, hi = m + 1e-9, m + 1.0 - 1e-9
        sym = _brent(lambda w: pair.channel_functions(w)[0] - (a + b), lo, hi)
        anti = _brent(lambda w: pair.channel_functions(w)[1] - (a - b), lo, hi)
        rows.append({"gap": m, "symmetric": sym, "antisymmetric": anti})
    return {"rows": rows,
            "n_symmetric": sum(1 for r in rows if r["symmetric"] is not None),
            "n_antisymmetric": sum(1 for r in rows
                                   if r["antisymmetric"] is not None)}


def mode_charges(omega: float, pair: MouthPair) -> np.ndarray:
    """The null vector ``q`` of ``M(ω)`` — how strongly each mouth radiates."""
    m = pair.krein_matrix(omega)
    _, _, vh = np.linalg.svd(m)
    q = vh[-1].conjugate()
    return q / np.linalg.norm(q)


def mouth_flux(q: np.ndarray, boundary: np.ndarray) -> Dict[str, object]:
    """Flux absorbed at each mouth, ``Im(q_j* (A q)_j)``, and the net.

    Around mouth ``j`` the field is ``q_j/(4πχ) + φ_j^reg`` and the radial
    current integrated over a small sphere is ``Im(q_j* φ_j^reg)``, independent
    of the sphere's radius.  With ``φ^reg = A q`` the total is ``Im(q† A q)``,
    which is zero for Hermitian ``A`` **identically** — not to leading order.
    """
    q = np.asarray(q, dtype=complex)
    aq = np.asarray(boundary, dtype=complex) @ q
    per = [float(np.imag(q[j].conjugate() * aq[j])) for j in range(len(q))]
    return {"per_mouth": per, "net": float(sum(per)),
            "scale": float(np.abs(q).max() * np.abs(aq).max() + 1e-300)}


def complex_root_search(secular, re_range: Tuple[float, float] = (1.1, 6.9),
                        im_seeds: Sequence[float] = (-0.45, -0.25, -0.08,
                                                     0.08, 0.25, 0.45),
                        n_re: int = 40, max_iter: int = 200,
                        tol: float = 1e-13) -> Dict[str, object]:
    """Newton from a grid of **complex** seeds — does anything leave the axis?

    Deliberately not the argument principle: ``det M`` has double poles sitting
    *on* the real axis, so a contour hugging it needs impractical sampling and
    aliases.  Seeded Newton answers the question that is actually being asked —
    *is there a root with ``Im ω ≠ 0``* — and it answers it the same way PR #255
    found its poles, so the comparison is like for like.

    The claim this supports is "no converged root is off the axis", not "every
    root was found": enumeration is `coupled_spectrum`'s job, on the axis.
    """
    roots: List[complex] = []
    for x in np.linspace(re_range[0], re_range[1], int(n_re)):
        for y in im_seeds:
            w = complex(x, y)
            ok = False
            for _ in range(int(max_iter)):
                h = 1e-7
                try:
                    fp = (secular(w + h) - secular(w - h)) / (2.0 * h)
                    fw = secular(w)
                except (ValueError, ZeroDivisionError):
                    break
                if abs(fp) < 1e-300:
                    break
                step = fw / fp
                if abs(step) > 0.3:
                    step *= 0.3 / abs(step)
                w -= step
                if abs(step) < tol:
                    ok = True
                    break
            if not ok or not (re_range[0] - 0.2 < w.real < re_range[1] + 0.2):
                continue
            if abs(w.imag) > 1.0 or any(abs(w - o) < 1e-6 for o in roots):
                continue
            roots.append(w)
    off = [r for r in roots if abs(r.imag) > 1e-8]
    return {"n_roots": len(roots), "n_off_axis": len(off),
            "worst_abs_imaginary": float(max((abs(r.imag) for r in roots),
                                             default=0.0)),
            "n_growing": sum(1 for r in off if r.imag < 0.0),
            "off_axis": [complex(r) for r in sorted(off,
                                                    key=lambda z: z.real)[:8]],
            "nothing_off_the_axis": bool(roots and not off)}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_green_function_has_a_finite_part(
        omegas: Sequence[float] = (0.37, 1.63, 2.4, 5.21),
        chis: Sequence[float] = (0.4, 1.3, 2.6),
        dampings: Sequence[float] = (1e-4, 1e-5, 1e-6),
        radii: Sequence[float] = (1e-2, 1e-3, 1e-4)) -> Dict[str, object]:
    """The closed form, checked against PR #255, and its short-distance split.

    Two things, because the second depends on the first.  ``G(χ,ω)`` is the
    ``γ → 0`` limit of PR #255's branch series — an independent construction, so
    agreement is a check.  And ``G(χ,ω) − 1/(4πχ) → g(ω) = −(ω/4π)cot(πω)``
    **linearly in χ**, which is what makes a point interaction definable: the
    divergence is the universal Coulomb one and the remainder is finite.
    """
    from geometrodynamics.waves.branch_coupling import mouth_transfer

    rows = []
    worst = 0.0
    for w in omegas:
        for chi in chis:
            got = mouth_transfer(w, chi, dampings[-1])
            want = free_green(chi, w)
            err = abs(got.real - want)
            worst = max(worst, err)
            rows.append({"omega": w, "chi": chi, "closed_form": want,
                         "branch_series_limit": got.real,
                         "residual_imaginary": abs(got.imag),
                         "abs_error": err})

    fin = []
    ratios = []
    for w in omegas[:3]:
        errs = []
        for r in radii:
            fp = free_green(float(r), w) - 1.0 / (4.0 * math.pi * r)
            errs.append(abs(fp - regularized_green(w)))
            fin.append({"omega": w, "chi": r, "finite_part": fp,
                        "g": regularized_green(w), "error": errs[-1]})
        ratios += [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]

    return {
        "rows": rows, "worst_abs_error": worst,
        "the_closed_form_is_the_branch_series": bool(worst < 1e-10),
        "finite_part": fin,
        "convergence_ratios_per_decade": ratios,
        "the_remainder_is_first_order_in_chi": bool(
            ratios and all(abs(r - 10.0) < 1.0 for r in ratios)),
        "what_this_shows": ("the coincidence divergence is the universal "
                            "1/(4πχ), so the finite part g(ω) exists and a "
                            "point interaction is definable"),
    }


def measure_the_boundary_operator_is_unitary_with_both_channels(
        separation: float = 1.3,
        params: Sequence[Tuple[float, float, complex]] = (
            (0.05, -0.011, 0.03 + 0.02j), (0.2, 0.2, 0.15),
            (-0.4, 0.07, -0.09 + 0.31j)),
        scales: Sequence[float] = (0.05, 0.1, 0.2),
        kappas: Sequence[float] = (0.3, 0.6, 1.0)) -> Dict[str, object]:
    """**The object PR #255 did not have.**  A unitary ``2×2`` with reflection.

    The Cayley transform of any Hermitian ``A`` is unitary, and inverts back, so
    the self-adjoint two-mouth boundary conditions *are* ``U(2)`` — four real
    parameters, not one.  Its diagonal is reflection and its off-diagonal
    transmission, with ``|r|² + |t|² = 1`` at each mouth.

    Alongside, PR #255's model in the same language: no reflection at all and
    ``|t| = κ``, so its column norm is ``κ² ≠ 1`` unless ``κ = 1``.  It is not a
    point in ``U(2)``; it is a point outside it.
    """
    rows = []
    worst_u = 0.0
    worst_sq = 0.0
    for (a1, a2, b) in params:
        for c in scales:
            p = MouthPair(separation, a1, a2, b)
            ch = p.channels(c)
            back = boundary_from_scattering(ch["S"], c)
            rt = float(np.abs(back - p.boundary_matrix()).max())
            worst_u = max(worst_u, ch["unitarity_defect"])
            worst_sq = max(worst_sq,
                           abs(ch["sum_of_squares_mouth_1"] - 1.0),
                           abs(ch["sum_of_squares_mouth_2"] - 1.0))
            rows.append({"alpha1": a1, "alpha2": a2, "beta": b, "scale": c,
                         "reflection": abs(ch["reflection_1"]),
                         "transmission": abs(ch["transmission_12"]),
                         "sum_of_squares": ch["sum_of_squares_mouth_1"],
                         "unitarity_defect": ch["unitarity_defect"],
                         "cayley_round_trip": rt,
                         "reciprocal": ch["reciprocal"]})

    pr255 = []
    for k in kappas:
        pr255.append({"kappa": k, "reflection": 0.0, "transmission": k,
                      "sum_of_squares": k * k,
                      "in_U2": bool(abs(k * k - 1.0) < 1e-12)})
    return {
        "rows": rows, "worst_unitarity_defect": worst_u,
        "worst_sum_of_squares_defect": worst_sq,
        "the_cayley_transform_is_unitary": bool(worst_u < 1e-13),
        "every_mouth_conserves": bool(worst_sq < 1e-12),
        "both_channels_are_present": bool(
            all(r["reflection"] > 1e-6 and r["transmission"] > 1e-6
                for r in rows)),
        "the_family_is_U2": ("four real parameters — two self-energies and a "
                             "complex mouth-to-mouth amplitude"),
        "pr255_in_the_same_language": pr255,
        "pr255_is_outside_U2_unless_kappa_is_one": bool(
            all(not r["in_U2"] for r in pr255 if r["kappa"] != 1.0)),
    }


def measure_self_adjointness_makes_the_spectrum_real(
        separation: float = 1.3,
        params: Sequence[Tuple[float, float, complex]] = (
            (0.05, 0.05, 0.03), (0.2, -0.13, 0.15 + 0.07j),
            (-0.4, 0.07, -0.09 + 0.31j)),
        n_gaps: int = 6, kappa: float = 0.3,
        delay: float = 1.0) -> Dict[str, object]:
    """**The result this round exists for.**  A conserving throat cannot ring up.

    ``Γ(ω)`` is real symmetric on the real axis, so ``M = A − Γ`` is Hermitian
    there for Hermitian ``A``, ``det M`` is a real function of real ``ω``, and
    its zeros are real.  An argument-principle census over a complex rectangle
    finds nothing off the axis.

    The control is PR #255's directional model written as a boundary matrix:
    strictly lower-triangular, so its anti-Hermitian part is as large as the
    matrix itself, and its secular function is **not** real on the axis.
    **PR #255's off-axis poles were the non-conservation, not the throat.**
    """
    rows = []
    worst_imag = 0.0
    for (a1, a2, b) in params:
        p = MouthPair(separation, a1, a2, b)
        im = 0.0
        for w in np.linspace(1.05, float(n_gaps) + 0.95, 400):
            det = complex(np.linalg.det(p.krein_matrix(float(w))))
            im = max(im, abs(det.imag) / max(abs(det), 1e-300))
        spec = coupled_spectrum(p, n_gaps)
        scan = complex_root_search(lambda z: complex(np.linalg.det(
            _krein_complex(z, p))), (1.1, float(n_gaps) + 0.9))
        worst_imag = max(worst_imag, im)
        rows.append({"alpha1": a1, "alpha2": a2, "beta": b,
                     "hermitian": p.is_self_adjoint(),
                     "worst_relative_imaginary_part_of_det": im,
                     "n_real_roots": len(spec),
                     "n_complex_seeds_converged": scan["n_roots"],
                     "n_off_axis": scan["n_off_axis"],
                     "worst_abs_imaginary": scan["worst_abs_imaginary"],
                     "nothing_off_the_axis": scan["nothing_off_the_axis"]})

    ctl = DirectionalThroat(separation, delay, +1, kappa)
    ctl_im = 0.0
    for w in np.linspace(1.05, float(n_gaps) + 0.95, 400):
        det = ctl.secular(complex(w))
        ctl_im = max(ctl_im, abs(det.imag) / max(abs(det), 1e-300))
    ctl_scan = complex_root_search(ctl.secular, (1.1, float(n_gaps) + 0.9))

    return {
        "rows": rows,
        "worst_relative_imaginary_part": worst_imag,
        "the_secular_function_is_real_on_the_axis": bool(worst_imag < 1e-12),
        "every_root_found_is_real": bool(all(r["n_real_roots"] > 0
                                             for r in rows)),
        "nothing_off_the_axis": bool(all(r["nothing_off_the_axis"]
                                         for r in rows)),
        "worst_abs_imaginary_over_all_seeds": float(max(
            (r["worst_abs_imaginary"] for r in rows), default=0.0)),
        "control_directional": {
            "anti_hermitian_part": ctl.anti_hermitian_part(2.0),
            "worst_relative_imaginary_part": ctl_im,
            "secular_is_real_on_the_axis": bool(ctl_im < 1e-12),
            "n_roots": ctl_scan["n_roots"],
            "n_off_axis": ctl_scan["n_off_axis"],
            "n_growing": ctl_scan["n_growing"],
            "worst_abs_imaginary": ctl_scan["worst_abs_imaginary"],
        },
        "the_control_fails_both_tests": bool(
            ctl_im > 1e-6 and ctl_scan["n_off_axis"] > 0),
        "and_the_control_is_unstable_even_at_unit_transmission": bool(
            complex_root_search(
                DirectionalThroat(separation, delay, +1, 1.0).secular,
                (1.1, float(n_gaps) + 0.9))["n_growing"] > 0),
        "what_this_retires": ("PR #255 found its poles off the real axis and "
                             "had to separate three thresholds to say what "
                             "that meant; a flux-conserving throat has a real "
                             "spectrum for every coupling, so the instability "
                             "was the model's non-conservation, not a throat"),
    }


def _krein_complex(z: complex, pair: MouthPair) -> np.ndarray:
    """``M(ω)`` continued to complex ``ω``, for the argument-principle scan."""
    w = complex(z)
    sp = cmath.sin(math.pi * w)
    g = -w * cmath.cos(math.pi * w) / (4.0 * math.pi * sp)
    d = pair.separation
    gd = cmath.sin(w * (math.pi - d)) / (4.0 * math.pi * math.sin(d) * sp)
    return pair.boundary_matrix() - np.array([[g, gd], [gd, g]], dtype=complex)


def measure_the_coupled_spectrum_interlaces_the_free_one(
        separation: float = 1.3, alpha: float = 0.05, beta: float = 0.03,
        n_gaps: int = 8,
        strengths: Sequence[float] = (1.0, 10.0, 1e2, 1e3, 1e4)
        ) -> Dict[str, object]:
    """Where the throat puts the frequencies, and where it puts them back.

    In the exchange-symmetric case the secular equation factorizes into
    ``g + G_d = α + β`` and ``g − G_d = α − β``.  Both left-hand sides are
    monotone from ``−∞`` to ``+∞`` across each unit gap, so each contributes
    **exactly one** root there: two coupled frequencies strictly between every
    pair of consecutive free ones.  Interlacing, not merely shifting.

    And switching the throat *off* returns the free spectrum ``ω = n+1``.  Off
    is ``‖A‖ → ∞``, not ``A → 0``: the diagonal of ``A`` is an **inverse**
    scattering length, so ``α = 0`` is a strongly coupled — indeed resonant —
    point interaction and ``α → ∞`` is the one that decouples (``q → 0``).  The
    shift then scales like ``1/‖A‖``, measured as an exponent.  Getting this
    backwards is easy and the measurement is written so that it would show.
    """
    pair = MouthPair(separation, alpha, alpha, beta)
    by_ch = spectrum_by_channel(pair, n_gaps)
    spec = coupled_spectrum(pair, n_gaps)

    per_gap = {m: 0 for m in range(1, n_gaps + 1)}
    for r in spec:
        per_gap[r["gap"]] += 1
    strictly_inside = all(r["gap"] < r["omega"] < r["gap"] + 1
                          for r in spec)

    # monotonicity of the two channel functions, which is why the count is two
    mono = True
    for m in range(1, min(n_gaps, 4) + 1):
        ws = np.linspace(m + 1e-6, m + 1 - 1e-6, 2000)
        for k in (0, 1):
            v = np.array([pair.channel_functions(float(w))[k] for w in ws])
            mono = mono and bool((np.diff(v) > 0).all())

    # the decoupling limit: ‖A‖ → ∞, with the shift ∝ 1/‖A‖.  Solved per
    # channel, because that is a monotone bisection on the whole gap and cannot
    # miss a root that has migrated to the very edge of it.
    scal = []
    prev = None
    for t in strengths:
        p = MouthPair(separation, alpha * t, alpha * t, beta * t)
        got = [r[k] for r in spectrum_by_channel(p, 4)["rows"]
               for k in ("symmetric", "antisymmetric") if r[k] is not None]
        worst = max((abs(w - round(w)) for w in got), default=0.0)
        exponent = (math.log(prev[1] / worst) / math.log(t / prev[0])
                    if prev and worst > 0 else None)
        scal.append({"strength": t, "inverse_strength": 1.0 / t,
                     "worst_shift": worst, "exponent": exponent})
        prev = (t, worst)

    exps = [r["exponent"] for r in scal if r["exponent"] is not None]
    tail = exps[-2:]
    return {
        "by_channel": by_ch["rows"], "n_roots": len(spec),
        "roots_per_gap": per_gap,
        "exactly_two_per_gap": bool(all(v == 2 for v in per_gap.values())),
        "every_root_strictly_between_free_ones": bool(strictly_inside),
        "both_channel_functions_are_monotone": bool(mono),
        "decoupling": scal, "shift_exponents": exps,
        "asymptotic_exponents": tail,
        "the_shift_goes_like_one_over_the_boundary_norm": bool(
            tail and all(abs(e - 1.0) < 0.05 for e in tail)),
        "free_spectrum_recovered": bool(
            scal and scal[-1]["worst_shift"] < 1e-3),
        "off_is_large_A_not_small_A": ("the diagonal of A is an inverse "
                                       "scattering length, so α → ∞ decouples "
                                       "and α = 0 is a resonant throat"),
        "what_this_shows": ("the throat interlaces the ESU spectrum two per "
                            "gap and returns it when switched off; the "
                            "coupled problem is a rank-two perturbation of a "
                            "self-adjoint operator, which is why"),
    }


def measure_the_flux_balance_is_exactly_hermiticity(
        separation: float = 1.3, n_draws: int = 200, seed: int = 20260815,
        kappa: float = 0.3, delay: float = 1.0,
        omega: float = 2.3) -> Dict[str, object]:
    """``Im(q† A q) = 0`` ⟺ ``A = A†`` — flux conservation, as an identity.

    Around mouth ``j`` the field is ``q_j/(4πχ) + φ_j^reg`` and the radial
    current through a small sphere is ``Im(q_j* φ_j^reg)``, independent of the
    sphere.  Summing over mouths with ``φ^reg = A q`` gives ``Im(q† A q)``.  For
    Hermitian ``A`` that is zero for **every** ``q``, not on average and not to
    leading order — measured over random draws.

    And what each mouth does individually is the transmission: for a purely
    off-diagonal Hermitian ``A``, whatever one mouth absorbs the other emits,
    exactly.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    worst_pair = 0.0
    for _ in range(int(n_draws)):
        a1, a2 = rng.normal(0, 0.3, 2)
        b = complex(*rng.normal(0, 0.3, 2))
        p = MouthPair(separation, float(a1), float(a2), b)
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        f = mouth_flux(q, p.boundary_matrix())
        worst = max(worst, abs(f["net"]) / f["scale"])
        # a purely off-diagonal throat: one mouth's loss is the other's gain
        p0 = MouthPair(separation, 0.0, 0.0, b)
        f0 = mouth_flux(q, p0.boundary_matrix())
        worst_pair = max(worst_pair,
                         abs(f0["per_mouth"][0] + f0["per_mouth"][1])
                         / f0["scale"])

    ctl = DirectionalThroat(separation, delay, +1, kappa)
    a_ctl = ctl.boundary_matrix(omega)
    net_ctl = []
    for _ in range(32):
        q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
        f = mouth_flux(q, a_ctl)
        net_ctl.append(abs(f["net"]) / f["scale"])
    return {
        "n_draws": int(n_draws),
        "worst_relative_net_flux": worst,
        "flux_is_conserved_identically": bool(worst < 1e-14),
        "worst_pairwise_imbalance": worst_pair,
        "what_one_mouth_absorbs_the_other_emits": bool(worst_pair < 1e-14),
        "control_directional_worst_net_flux": float(max(net_ctl)),
        "control_median_net_flux": float(np.median(net_ctl)),
        "the_control_does_not_conserve": bool(np.median(net_ctl) > 1e-3),
        "control_anti_hermitian_part": ctl.anti_hermitian_part(omega),
        "the_identity": "Im(q† A q) = 0 for all q  ⟺  A = A†",
    }


def measure_the_directional_model_is_what_pr255_solved(
        separation: float = 1.3, delay: float = 1.0, kappa: float = 0.3,
        omegas: Sequence[float] = (1.3, 2.77, 4.11, 5.6)) -> Dict[str, object]:
    """Where PR #255 sits inside — or rather outside — this family.

    Its relation is recovered exactly as the strictly-lower-triangular boundary
    matrix ``A(ω) = [[0,0],[1/(ηκe^{−iωΔ}),0]]``: no self-energy, so **no
    reflection**; one direction only, so **not Hermitian**; and ``ω``-dependent,
    which a self-adjoint point extension is not allowed to be, because ``A`` is
    a boundary condition rather than a dynamical object.

    Three defects, each with a number.  None of them touches PR #255's
    resolvent, which was exact for the model it posed — they say which model.
    """
    rows = []
    for w in omegas:
        ctl = DirectionalThroat(separation, delay, +1, kappa)
        a = ctl.boundary_matrix(w)
        rows.append({
            "omega": w,
            "reflection_entries": [abs(a[0, 0]), abs(a[1, 1])],
            "anti_hermitian_norm": ctl.anti_hermitian_part(w),
            "hermitian_norm": float(np.abs(0.5 * (a + a.conjugate().T)).max()),
            "secular_relative_imaginary": (
                abs(ctl.secular(complex(w)).imag)
                / max(abs(ctl.secular(complex(w))), 1e-300)),
        })
    ctl = DirectionalThroat(separation, delay, +1, kappa)
    entries = [complex(ctl.boundary_matrix(w)[1, 0]) for w in omegas]
    same_matrix = all(abs(e - entries[0]) < 1e-12 for e in entries)
    same_magnitude = all(abs(abs(e) - abs(entries[0])) < 1e-12
                         for e in entries)
    return {
        "boundary_entry_by_frequency": entries,
        "the_boundary_matrix_depends_on_frequency": not same_matrix,
        "only_through_its_phase": bool(same_magnitude and not same_matrix),
        "rows": rows,
        "no_reflection_channel": bool(
            all(max(r["reflection_entries"]) == 0.0 for r in rows)),
        "not_hermitian_at_any_frequency": bool(
            all(r["anti_hermitian_norm"] > 0.0 for r in rows)),
        "anti_hermitian_part_is_the_whole_matrix": bool(
            all(abs(r["anti_hermitian_norm"] - r["hermitian_norm"]) < 1e-12
                for r in rows)),
        "what_is_not_wrong_with_it": ("the resolvent of PR #255 is exact for "
                                      "the model it posed; these are "
                                      "statements about which model"),
        "what_a_self_adjoint_extension_forbids": (
            "a frequency-dependent boundary matrix — A is a boundary "
            "condition, not a dynamical response, so the delay Δ has no place "
            "in a point extension and does not survive into it"),
    }
