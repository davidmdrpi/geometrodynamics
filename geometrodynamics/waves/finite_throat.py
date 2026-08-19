"""A **finite** conservative throat: a Dirichlet-to-Neumann map, a point limit,
and a traversal delay.

The debt this pays
──────────────────
Every round from PR #253 to PR #259 carried the same disclaimer: *the throat is
point-supported — no interior, no proper length, no delay*, and what is solved
is a rank-one **mouth-transfer** model, field values only, ``1×1`` where a
conserving junction needs ``2×2``, and lossy for ``κ < 1``.  This module
replaces it with an object that has an interior.

A tube of proper length ``L``, cross-section ``𝒜`` and interior mass ``m``
joins the two mouths.  Its **Dirichlet-to-Neumann map** is exact,

    ``N(λ)  =  𝒜k · [[cot kL, −csc kL], [−csc kL, cot kL]]`` ,  ``k² = λ − m²``

and the matching to the ambient ``S³`` is value continuity and flux continuity
at the mouths.  Everything below follows from those two lines; nothing else is
put in.

Where the self-adjointness lives
───────────────────────────────
The conservative object is the **enlarged system** — ambient ``⊕`` tube — with
the *fixed*, ``λ``-independent matching above.  That is one self-adjoint
operator on ``L²(S³) ⊕ L²([0,L])``, and it is what conserves flux.

Eliminating the tube leaves an ambient-only problem with a **``λ``-dependent**
boundary condition, ``A(λ) = −N(λ)⁻¹``: the Weyl (``M``-) function of that
elimination.  ``A(λ)`` is *not* itself a self-adjoint operator on the ambient
space, and this module does not claim it is — an energy-dependent boundary
condition never is.  What it is, is a matrix **Nevanlinna** function, monotone
in ``λ`` between its poles, and that monotonicity is the enlarged system's
self-adjointness showing through after the elimination.  It is measured, not
assumed.

That ``λ``-dependence *is* the interior.  Four consequences, each measured:

1. **A traversal time.**  Expanding on the retarded contour,

       ``cot x = −i − 2i Σ_{k≥1} e^{2ikx}`` ,   ``csc x = −2i Σ_{k≥0} e^{i(2k+1)x}``

   so the same-mouth entry carries delays ``0, 2L, 4L, …`` — an instantaneous
   reflection and its echoes — and the cross-mouth entry carries
   ``L, 3L, 5L, …``.  **The throat transmits at ``t = L``, not at ``t = 0``.**
   A frozen ``A`` transmits instantaneously.

2. **The short-tube limit exists — and it leaves the chart.**  As ``L → 0`` the
   boundary pair does *not* approach a generic finite ``A``; rescaled, it
   approaches the projector pair ``(B, C) → (P_a, −P_s)``, i.e. the stratum

       ``q_sym = 0``   and   ``Φ_anti = 0``

   — a **mixed Dirichlet–Neumann** point extension, maximal and self-adjoint
   (``rank[B|C] = 2`` throughout), and reachable by no finite Hermitian ``A``
   because both ``B`` and ``C`` are singular there.  It is exactly the kind of
   stratum PR #257's review insisted the finite-``A`` chart does not cover.  So
   the honest statement is *not* "there is no point limit" but **"no finite-``A``
   point limit"**.

3. **The static response is rank one — which falsifies the finite-``A`` chart,
   not point-ness.**  A massless tube holds a zero mode: a constant field on it
   costs nothing, carries no current, and the symmetric channel decouples.  So
   ``det S → 0`` linearly in ``λ − m²`` and PR #258's defect
   ``𝒲 = S₁₂/det S − G₀`` diverges.  But the limiting stratum in (2) is rank one
   *too* — ``R → diag(0, −1/Γ_a)`` in the channel basis — so what rank one
   distinguishes is this tube **from the generic rank-two finite-``A`` family**,
   and not a finite throat from a point one.  Off the collapse ``𝒲 = −β(λ)``
   exactly, with ``β(λ) = csc(kL)/(𝒜k)``: PR #258's theorem survives the
   generalization and returns the interior's amplitude rather than a constant.

4. **The model fails the stability gate, and the failure is the point mouth's.**
   ``A(λ)`` decreases and ``Γ(λ)`` increases (PR #257's Gram identity), so
   ``A − Γ`` is strictly monotone between poles and each channel has at most one
   root — a count, not a scan.  The symmetric channel always has exactly one,
   and it is at ``λ < 0``: **an exponentially growing mode, for every choice of
   parameters.**  Its rate,

       ``σ*  =  ½[ 1/a + √(1/a² + 16π/𝒜) ]``   →  ``2√(π/𝒜)``  for point mouths

   contains **neither ``L`` nor the mouth separation ``d``**, and the two
   channels split by ``1.04·e^{−σ*d}`` — the Euclidean propagator between the
   mouths.  A mode blind to the tube's length and to the separation, and
   degenerate between the channels, is a **single-mouth object**: the growing
   mode is the *point-mouth matching's*, not the tube's.

   This is the round's falsification result, and it is not cured anywhere below.
   Placing the retarded contour above ``σ*`` evaluates the correct retarded
   solution **of an unstable system**; it does not stabilize it.  Whether a
   finite-radius mouth or neck geometry removes the mode is open, and is the
   thing to settle before any stationary-action or backreaction work.

Which frequencies each result uses
──────────────────────────────────
Worth stating plainly, because the band is not uniform across the round.  The
delay and the bounce ledger are statements about the **analytic structure of the
exact model at all frequencies**: a causal onset is a UV object, and the probe
pulse used to resolve it (width ``0.03``) carries content out to ``ω ∼ 30``, far
above ``σ* ∼ 1.4``.  They are exact results *about this model*, not predictions
about a resolved physical mouth.  The static rank collapse, the defect and the
Weyl-function monotonicity are low-frequency statements and sit inside the band.

Relatedly, ``𝒜`` is a **one-dimensional coupling strength**, not an area with a
mouth radius attached: the interior has a single transverse channel and no
higher modes.  Reading ``√(𝒜/4π)`` as a mouth radius at the working point would
give a radius of order the whole unit ``S³``, which is another way of saying the
same thing.

What is put in
──────────────
The background, the mouth positions, and three numbers with dimensions —
``L``, ``𝒜``, ``m`` — where the *real*-field point sector of PRs #257–#259 has
three real numbers without them (``α₁, α₂, β`` with ``β`` real, because a real
scalar forces ``A = A*``).  **No backreaction**: the throat is still a fixed
background, and nothing here solves for the geometry that would produce it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .throat_operator import DirectionalThroat, MouthPair, gamma_at
from .throat_positivity import positivity_defect
from .two_wave import GaussianPulse, RetardedGrid, gamma_omega

__all__ = [
    "FiniteThroat",
    "WORKING_THROAT",
    "dtn_matrix",
    "interior_profile",
    "green_identity_residual",
    "bounce_delays",
    "response_spectrum",
    "impulse_response",
    "causal_onset",
    "channel_basis",
    "short_tube_stratum",
    "measure_the_enlarged_system_is_conservative",
    "measure_the_throat_transmits_at_the_traversal_time",
    "measure_the_delay_ledger_is_the_bounce_series",
    "measure_the_short_tube_limit_is_a_mixed_stratum",
    "measure_the_static_limit_is_rank_one_and_the_defect_diverges",
    "measure_the_interior_mass_is_a_transmission_cutoff",
    "measure_the_growing_mode_belongs_to_the_mouth",
    "measure_the_contour_must_clear_the_growing_mode",
]

FOUR_PI = 4.0 * math.pi


# ════════════════════════════════════════════════════════════════════════════
# THE TUBE'S DIRICHLET-TO-NEUMANN MAP
# ════════════════════════════════════════════════════════════════════════════
def dtn_matrix(k, length: float, area: float) -> np.ndarray:
    """``N = 𝒜k[[cot kL, −csc kL], [−csc kL, cot kL]]`` — array-aware in ``k``.

    The interior problem is ``u'' + k²u = 0`` on ``[0, L]`` with the mouth
    values as Dirichlet data; ``N`` returns the outward flux at each end.  Both
    entries are **even** in ``k``, so ``k = √(λ − m²)`` needs no branch choice —
    ``N`` is a meromorphic function of ``λ`` with no cut, and the square root is
    cosmetic.

    Poles at ``kL = nπ`` are the tube's own Dirichlet resonances: at those
    frequencies the interior supports a mode with both ends held at zero, and
    the coupling to the ambient diverges.  They are real, they are the
    interior's spectrum, and a point throat has none of them.
    """
    z = np.asarray(k, dtype=complex)
    x = z * float(length)
    ct = float(area) * z / np.tan(x)
    cs = float(area) * z / np.sin(x)
    if z.ndim == 0:
        return np.array([[ct, -cs], [-cs, ct]], dtype=complex)
    out = np.empty(z.shape + (2, 2), dtype=complex)
    out[..., 0, 0] = out[..., 1, 1] = ct
    out[..., 0, 1] = out[..., 1, 0] = -cs
    return out


def interior_profile(k: complex, length: float, ends: Sequence[complex],
                     n: int = 4001) -> Tuple[np.ndarray, np.ndarray,
                                             np.ndarray]:
    """``(s, u(s), u'(s))`` — the interior solution with the given end values.

    Used to check the DtN map against the interior it claims to summarize,
    rather than against itself.  The derivative is analytic: PR #248's lesson
    was that `np.gradient`'s one-sided end differences are exactly what a
    boundary-value identity is most sensitive to.
    """
    s = np.linspace(0.0, float(length), int(n))
    a, b = complex(ends[0]), complex(ends[1])
    sk = np.sin(k * length)
    c = (b - a * np.cos(k * length)) / sk
    u = a * np.cos(k * s) + c * np.sin(k * s)
    du = k * (-a * np.sin(k * s) + c * np.cos(k * s))
    return s, u, du


def green_identity_residual(k: float, length: float, area: float,
                            ends: Sequence[complex], n: int = 40001) -> float:
    """``|Φ†NΦ − 𝒜∫(|u'|² − k²|u|²)ds|`` — the DtN's defining identity.

    Green's identity for ``u'' + k²u = 0`` says the boundary form built from the
    outward fluxes equals the interior energy.  Checking it by quadrature is the
    one test that the matrix above really is the tube's map and not a plausible
    matrix: an error in a sign or in a ``cot``/``csc`` fails it immediately.

    The form is **sesquilinear** — ``Φ†NΦ`` and ``|u'|²``, not ``ΦᵀNΦ`` and
    ``u'²``.  For real ``k`` and real end values the two agree, which is exactly
    why the bilinear version passed while saying something weaker: it is the
    sesquilinear identity that expresses *energy*, and only that one fails when
    the data is complex.  ``k`` is required real here for the same reason —
    below the cutoff the identity is the one for ``κ² = −k²``.
    """
    kk = float(np.real(k))
    if abs(np.imag(k)) > 1e-14:
        raise ValueError("the energy identity is stated for real k; below the "
                         "cutoff use kappa² = −k²")
    s, u, du = interior_profile(complex(kk), length, ends, n)
    bulk = float(area) * float(np.trapezoid(
        np.abs(du) ** 2 - (kk ** 2) * np.abs(u) ** 2, s).real)
    phi = np.array([complex(ends[0]), complex(ends[1])], dtype=complex)
    form = complex(phi.conjugate() @ dtn_matrix(kk, length, area) @ phi)
    return float(abs(form - bulk))


_CHANNEL_BASIS = np.array([[1.0, 1.0], [1.0, -1.0]]) / math.sqrt(2.0)
_P_SYM = np.diag([1.0, 0.0])
_P_ANTI = np.diag([0.0, 1.0])


def channel_basis() -> np.ndarray:
    """``V = [v_sym | v_anti]`` — the basis both ``N`` and ``Γ`` diagonalize."""
    return _CHANNEL_BASIS.copy()


def short_tube_stratum() -> Tuple[np.ndarray, np.ndarray]:
    """The boundary pair the tube approaches as ``L → 0``, in the channel basis.

    ``(B, C) = (P_anti, −P_sym)``, i.e. ``Φ_anti = 0`` and ``q_sym = 0``: a
    **mixed Dirichlet–Neumann** point extension.  Maximal and self-adjoint —
    ``rank[B|C] = 2``, ``BC† = 0`` — and reached by no finite Hermitian ``A``,
    since both blocks are singular.

    A tube shorter than every other scale is therefore not "no throat" and not a
    weakly coupled one: it is a *specific* self-adjoint extension, the one that
    forbids net charge in the symmetric channel and forces the antisymmetric
    mouth values to agree.  Which is what a very short pipe does — it shorts the
    two mouths together and stores nothing.
    """
    return _P_ANTI.copy(), -_P_SYM.copy()


def bounce_delays(order: int = 3) -> Dict[str, List[float]]:
    """The delays the two entries of ``A(λ)`` carry, in units of ``L``.

    Same mouth in and out is an even number of traversals — ``0`` is the
    instantaneous reflection off the mouth itself — and opposite mouths is an
    odd number.  This is the ledger PR #253's arrival list was missing, and it
    is a derivation from ``cot`` and ``csc``, not a fit.
    """
    return {"same_mouth": [2.0 * i for i in range(int(order))],
            "opposite_mouths": [2.0 * i + 1.0 for i in range(int(order))]}


# ════════════════════════════════════════════════════════════════════════════
# THE THROAT
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class FiniteThroat:
    """A tube of length ``L``, area ``𝒜`` and interior mass ``m`` between two
    mouths a geodesic distance ``d`` apart.

    ``mouth_radius`` is optional.  ``None`` is the point-mouth matching, where
    the tube's end value is the ambient's *regularized* value; a finite radius
    ``a`` says the tube sees the field at radius ``a``, which exceeds the
    regular part by ``q/(4πa)``.  It changes ``σ*`` and nothing qualitative,
    and it is carried because the growing mode's closed form is the clearest
    statement of what the point-mouth idealization costs.
    """

    separation: float = 1.3
    length: float = 0.9
    area: float = FOUR_PI
    interior_mass: float = 0.0
    mouth_radius: Optional[float] = None

    # -- the interior --------------------------------------------------------
    def wavenumber(self, lmbda):
        """``k = √(λ − m²)``.  Only ``k²`` is ever used; see `dtn_matrix`."""
        return np.sqrt(np.asarray(lmbda, dtype=complex)
                       - float(self.interior_mass) ** 2)

    def dtn(self, lmbda) -> np.ndarray:
        return dtn_matrix(self.wavenumber(lmbda), self.length, self.area)

    def mouth_term(self) -> float:
        """``1/(4πa)``, or zero for point mouths."""
        if self.mouth_radius is None:
            return 0.0
        return 1.0 / (FOUR_PI * float(self.mouth_radius))

    # -- the boundary condition ----------------------------------------------
    def normalized_pair(self, lmbda: float) -> Tuple[np.ndarray, np.ndarray]:
        """``(B, C)`` in the channel basis, row-scaled so ``L → 0`` is finite.

        A boundary pair is defined only up to ``(B, C) → (MB, MC)`` for
        invertible ``M``, so a limit that looks singular in one normalization
        can be perfectly finite in another — which is exactly what happens here.
        Scaling the symmetric row by 1 and the antisymmetric row by ``1/N_anti``
        gives ``B = diag(N_sym, 1)`` and ``C = −diag(1, 1/N_anti)``, both of
        which converge as ``L → 0``, to `short_tube_stratum`.
        """
        lam = float(lmbda)
        n = self.dtn(complex(lam))
        v = _CHANNEL_BASIS
        diag = np.real(v.T @ n @ v)
        n_sym, n_anti = float(diag[0, 0]), float(diag[1, 1])
        return (np.diag([n_sym, 1.0]),
                -np.diag([1.0, 1.0 / n_anti]))

    def boundary_pair(self, lmbda) -> Tuple[np.ndarray, np.ndarray]:
        """``(B, C)`` with ``B φ^reg = C q`` — valid where ``A`` is not.

        Flux continuity is ``q = −N Φ`` and value continuity is
        ``Φ = φ^reg + q/(4πa)``, so ``B = N`` and ``C = −(I + N/(4πa))``.  The
        pair form is the honest one: at ``λ = m²`` the tube's symmetric channel
        decouples, ``N`` is singular, and no finite ``A`` exists — exactly the
        Dirichlet stratum PR #257's review insisted the chart does not reach.
        """
        n = self.dtn(lmbda)
        eye = np.eye(2, dtype=complex)
        return n, -(eye + n * self.mouth_term())

    def boundary_matrix(self, lmbda) -> np.ndarray:
        """``A(λ) = −N(λ)⁻¹ − I/(4πa)`` — the chart, where it exists.

        In closed form ``−N⁻¹ = (1/𝒜k)[[cot kL, csc kL], [csc kL, cot kL]]``,
        since ``det N = −(𝒜k)²`` and ``cot² − csc² = −1``.  So the throat's
        **transmission amplitude** is ``β(λ) = csc(kL)/(𝒜k)`` and its
        **self-energy** is ``α(λ) = cot(kL)/(𝒜k) − 1/(4πa)``: one function of
        the interior, evaluated at two arguments.
        """
        z = np.asarray(lmbda, dtype=complex)
        k = self.wavenumber(z)
        x = k * float(self.length)
        c = 1.0 / (float(self.area) * k)
        diag = c / np.tan(x) - self.mouth_term()
        off = c / np.sin(x)
        if z.ndim == 0:
            return np.array([[diag, off], [off, diag]], dtype=complex)
        out = np.empty(z.shape + (2, 2), dtype=complex)
        out[..., 0, 0] = out[..., 1, 1] = diag
        out[..., 0, 1] = out[..., 1, 0] = off
        return out

    def transmission(self, lmbda):
        """``β(λ) = csc(kL)/(𝒜k)`` — what PR #258's defect returns."""
        k = self.wavenumber(lmbda)
        return 1.0 / (float(self.area) * k * np.sin(k * float(self.length)))

    def self_energy(self, lmbda):
        """``α(λ) = cot(kL)/(𝒜k) − 1/(4πa)``."""
        k = self.wavenumber(lmbda)
        return (1.0 / (float(self.area) * k * np.tan(k * float(self.length)))
                - self.mouth_term())

    def krein_matrix(self, lmbda) -> np.ndarray:
        """``M(λ) = A(λ) − Γ(λ)``, the object whose zeros are the spectrum."""
        z = np.asarray(lmbda, dtype=complex)
        if z.ndim == 0:
            return (self.boundary_matrix(z)
                    - gamma_at(float(z.real), self.separation))
        return self.boundary_matrix(z) - gamma_omega(np.sqrt(z),
                                                     self.separation)

    def channel_functions(self, lmbda: float) -> Tuple[float, float]:
        """``(A−Γ)`` in the symmetric and antisymmetric channels, real ``λ``.

        Both matrices are ``pI + q·swap``, so they commute and the ``2×2``
        problem is two scalar ones.  ``A_sym = −tan(kL/2)/(𝒜k)`` and
        ``A_anti = cot(kL/2)/(𝒜k)`` — the half-angle is the tube folded about
        its middle, which is what the symmetric and antisymmetric modes are.
        """
        lam = float(lmbda)
        a = self.boundary_matrix(complex(lam))
        g = gamma_at(lam, self.separation)
        s = float((a[0, 0] + a[0, 1] - g[0, 0] - g[0, 1]).real)
        t = float((a[0, 0] - a[0, 1] - g[0, 0] + g[0, 1]).real)
        return s, t

    # -- the growing mode ----------------------------------------------------
    def negative_lambda_channels(self, sigma: float) -> Tuple[float, float]:
        """``(A−Γ)`` in both channels at ``λ = −σ²``, without overflowing.

        Written in decaying exponentials throughout.  ``throat_operator``'s
        ``gamma_at`` forms ``sinh(σ(π−d))/sinh(πσ)`` directly and overflows
        above ``σ ≈ 350``, which is inside the range this scan needs; the same
        ratio in the form

            ``G_d(σ)  =  (e^{−σd} − e^{−σ(2π−d)}) / (4π sin d (1 − e^{−2πσ}))``

        is stable everywhere — and it is also the statement that the Euclidean
        mouth-to-mouth propagator is a sum over the **two ways round the
        sphere**, PR #253's ledger with the oscillation replaced by decay.  The
        two constructions are checked against each other where both are finite.

        The tube's side is ``A_sym = −coth(κL/2)/(𝒜κ)`` and
        ``A_anti = −tanh(κL/2)/(𝒜κ)`` with ``κ = √(σ² + m²)`` — the half-angle
        is the tube folded about its middle.
        """
        s = float(sigma)
        d = float(self.separation)
        e = math.pi - d
        two_pi_s = math.exp(-2.0 * math.pi * s)
        g = -s * (1.0 + two_pi_s) / (FOUR_PI * (1.0 - two_pi_s))
        if abs(math.sin(e)) < 1e-12:                     # the antipode, by limit
            gd = 2.0 * s * math.exp(-math.pi * s) / (FOUR_PI * (1.0 - two_pi_s))
        else:
            gd = ((math.exp(-s * d) - math.exp(-s * (2.0 * math.pi - d)))
                  / (FOUR_PI * math.sin(d) * (1.0 - two_pi_s)))
        kappa = math.sqrt(s ** 2 + float(self.interior_mass) ** 2)
        x = math.exp(-kappa * float(self.length))
        scale = float(self.area) * kappa
        a_sym = -(1.0 + x) / ((1.0 - x) * scale) - self.mouth_term()
        a_anti = -(1.0 - x) / ((1.0 + x) * scale) - self.mouth_term()
        return a_sym - (g + gd), a_anti - (g - gd)

    def sigma_star_closed_form(self) -> float:
        """``½[1/a + √(1/a² + 16π/𝒜)]`` — where the growing mode sits.

        Valid once ``σL ≫ 1`` and ``σd ≫ 1``, where ``coth(σL/2) → 1``,
        ``coth(πσ) → 1`` and the mouth-to-mouth term has died: the crossing
        condition ``−1/(𝒜σ) − 1/(4πa) = −σ/4π`` is then a quadratic.  **``L``
        and ``d`` have both dropped out**, which is the whole point — the mode
        knows the mouth's scale and nothing else.
        """
        inv_a = FOUR_PI * self.mouth_term()          # 1/a, or 0
        return 0.5 * (inv_a + math.sqrt(inv_a ** 2
                                        + 16.0 * math.pi / float(self.area)))

    def growing_modes(self, sigma_max: float = 200.0,
                      n: int = 4000) -> Dict[str, object]:
        """Roots of the channel functions at ``λ < 0``, found by monotonicity.

        ``A(λ)`` decreases and ``Γ(λ)`` increases, so each channel function is
        strictly monotone between poles and a sign change is a root — bisection
        is exact here rather than a hopeful scan.  The symmetric channel always
        has one; the antisymmetric one does only for a thin enough tube.
        """
        sig = np.geomspace(1e-3, float(sigma_max), int(n))
        out: Dict[str, object] = {}
        for name, idx in (("symmetric", 0), ("antisymmetric", 1)):
            vals = np.array([self.negative_lambda_channels(s)[idx]
                             for s in sig])
            root = None
            sign = np.sign(vals)
            cross = np.where(sign[:-1] * sign[1:] < 0.0)[0]
            if cross.size:
                lo, hi = sig[cross[0]], sig[cross[0] + 1]
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    if (self.negative_lambda_channels(lo)[idx]
                            * self.negative_lambda_channels(mid)[idx]) <= 0.0:
                        hi = mid
                    else:
                        lo = mid
                root = 0.5 * (lo + hi)
            out[name] = root
        s_sym = out["symmetric"]
        out["closed_form"] = self.sigma_star_closed_form()
        out["channel_split"] = (
            abs(float(s_sym) - float(out["antisymmetric"]))
            if out["antisymmetric"] is not None and s_sym is not None
            else None)
        out["the_symmetric_channel_always_has_one"] = bool(s_sym is not None)
        return out

    # -- comparisons ---------------------------------------------------------
    def matched_point_throat(self, lmbda: float) -> MouthPair:
        """The constant-``A`` `MouthPair` that agrees with this throat at one
        ``λ``.  Everything PRs #257–#259 did applies to it exactly — and only
        at that one point, which is the content of `the_point_throat_is_a
        _single_frequency_match`."""
        a = self.boundary_matrix(complex(float(lmbda)))
        return MouthPair(self.separation, float(a[0, 0].real),
                         float(a[1, 1].real), complex(a[0, 1].real))

    def loewner_margin(self, lmbda: float) -> float:
        """``λ_min(A(λ) − Γ(0))`` — PR #257's margin, read at one frequency."""
        return float(positivity_defect(
            self.boundary_matrix(complex(float(lmbda))).real.astype(complex),
            self.separation)["min_eigenvalue"])

    def resonances(self, n: int = 4) -> List[float]:
        """``λ = m² + (nπ/L)²`` — the interior's Dirichlet spectrum."""
        return [float(self.interior_mass ** 2
                      + (i * math.pi / self.length) ** 2)
                for i in range(1, int(n) + 1)]


WORKING_THROAT = FiniteThroat(separation=1.3, length=0.9, area=FOUR_PI)


# ════════════════════════════════════════════════════════════════════════════
# THE THROAT'S OWN IMPULSE RESPONSE
# ════════════════════════════════════════════════════════════════════════════
def response_spectrum(throat: FiniteThroat, omegas: np.ndarray,
                      constant: Optional[np.ndarray] = None) -> np.ndarray:
    """``R(ω) = (A(ω) − Γ(ω))⁻¹`` on the contour, inverted in closed form.

    Passing ``constant`` freezes ``A`` at a matrix, which is the point-throat
    control: same ambient, same mouths, same everything, one function replaced
    by one number.
    """
    om = np.asarray(omegas, dtype=complex)
    if constant is None:
        a = throat.boundary_matrix(om ** 2)
    else:
        a = np.broadcast_to(np.asarray(constant, dtype=complex),
                            om.shape + (2, 2))
    m = a - gamma_omega(om, throat.separation)
    det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
    r = np.empty_like(m)
    r[..., 0, 0] = m[..., 1, 1] / det
    r[..., 1, 1] = m[..., 0, 0] / det
    r[..., 0, 1] = -m[..., 0, 1] / det
    r[..., 1, 0] = -m[..., 1, 0] / det
    return r


def impulse_response(throat: FiniteThroat, grid: RetardedGrid,
                     width: float = 0.03,
                     constant: Optional[np.ndarray] = None
                     ) -> Dict[str, np.ndarray]:
    """``r_ij(t)`` — the **coupled** mouth-to-mouth response to a short pulse.

    ``R = (A − Γ)⁻¹`` is the ambient *and* tube together: ``Γ`` carries the
    ambient's own mouth-to-mouth propagation, so this is not the throat in
    isolation.  That is exactly why the cross-mouth onset below is ``min(L, d)``
    and not ``L``.  What the construction does remove is the *source and
    observer* geometry — the legs from a source to a mouth and from a mouth to
    an observer — leaving the two-mouth block on its own.

    Probed with a narrow Gaussian rather than a delta, because ``R(ω)`` decays
    only like ``1/ω`` and an unsmoothed inversion is all Gibbs.  The pulse width
    is the time resolution, and it is what the onset measurement is quoted
    against.

    The returned series is the **contour amplitude** ``φ(t)e^{−εt}``, not the
    field: the growing mode makes ``e^{εt}`` enormous at late times, so a
    threshold relative to the field's peak would sit far above the causal onset.
    Causality is a statement about the contour amplitude, and this is it.
    """
    om = grid.omegas
    spec = GaussianPulse(amplitude=1.0, carrier=0.0, width=float(width),
                         t0=0.0).spectrum(om)
    r = response_spectrum(throat, om, constant)
    out = {}
    for name, (i, j) in (("same_mouth", (0, 0)), ("opposite_mouths", (0, 1))):
        out[name] = np.real(np.fft.fft(r[:, i, j] * spec) / grid.span)
    out["times"] = grid.times
    return out


def causal_onset(series: np.ndarray, times: np.ndarray,
                 floor: float = 1e-9) -> float:
    """The first time the series leaves its own numerical zero.

    ``floor`` is relative to the series' peak, so the answer is early by the
    probe pulse's tail — a fixed offset, identical across a sweep, which is why
    the sweep's **slope** rather than its intercept is what gets quoted.
    """
    s = np.asarray(series, dtype=float)
    peak = float(np.abs(s).max())
    idx = int(np.argmax(np.abs(s) > float(floor) * peak))
    return float(np.asarray(times)[idx])


def _grid_for(throat: FiniteThroat, clearance: float = 0.8,
              n: int = 1 << 17, span: float = 300.0) -> RetardedGrid:
    """A contour placed above the growing mode, by ``clearance``."""
    modes = throat.growing_modes()
    sig = max([float(v) for v in (modes["symmetric"], modes["antisymmetric"])
               if v is not None] or [0.0])
    return RetardedGrid(n=int(n), span=float(span), eps=sig + float(clearance))


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_enlarged_system_is_conservative(
        throat: FiniteThroat = WORKING_THROAT) -> Dict[str, object]:
    """Conservation, stated where it actually lives.

    The self-adjoint object is the **enlarged system**, ambient ``⊕`` tube, with
    the ``λ``-independent matching ``q = −NΦ``, ``Φ = φ^reg + q/(4πa)``.  What is
    checked here is that eliminating the tube leaves, at each ``λ``, a *maximal
    self-adjoint boundary condition* for the ambient problem at that ``λ``:
    ``rank[B|C] = 2`` with ``BC†`` Hermitian, both exactly, since ``N`` is real
    symmetric.  That is a statement about the elimination being faithful — it is
    **not** the claim that ``A(λ)`` is one self-adjoint operator on the ambient
    space, which an energy-dependent boundary condition never is.

    The fingerprint of the enlarged system's self-adjointness, on the ambient
    side, is that ``A(λ)`` is a matrix **Nevanlinna function** — monotone in
    ``λ`` between its poles.  That is measured in
    `the_boundary_matrix_decreases_in_lambda` and is what makes the mode count
    below a count rather than a scan.

    The control is PR #255's `DirectionalThroat`, the rank-one mouth-transfer
    model this module replaces: it drives ``M⁻`` from ``M⁺`` and not the
    reverse, so ``BC† = B`` is **not** Hermitian, and the defect is the size of
    the coupling itself.  That is the ``κ < 1`` loss, stated as an operator
    property rather than as a bookkeeping slip.

    The DtN map is also checked against the interior it summarizes — by the
    **sesquilinear** Green's identity under quadrature, the one that expresses
    energy — rather than against itself.
    """
    lams = [-4.0, -0.7, 0.15, 0.9, 2.5, 7.0, 19.0]
    rows = []
    for lam in lams:
        b, c = throat.boundary_pair(complex(lam))
        prod = b @ c.conjugate().T
        rows.append({
            "lambda": float(lam),
            "hermiticity_defect": float(np.abs(
                prod - prod.conjugate().T).max()),
            "rank": int(np.linalg.matrix_rank(np.hstack([b, c]), tol=1e-10)),
            "imaginary_part": float(np.abs(np.hstack([b, c]).imag).max())})
    control = DirectionalThroat(separation=throat.separation, kappa=0.3)
    control_bc = control.boundary_condition(1.0)
    cprod = control_bc.B @ control_bc.C.conjugate().T
    green = [green_identity_residual(complex(k), throat.length, throat.area,
                                     (1.0, -0.4))
             for k in (0.7, 1.9, 3.3)]
    return {"rows": rows,
            "worst_hermiticity_defect": float(
                max(r["hermiticity_defect"] for r in rows)),
            "every_condition_is_maximal": bool(
                all(r["rank"] == 2 for r in rows)),
            "worst_imaginary_part": float(
                max(r["imaginary_part"] for r in rows)),
            "worst_green_identity_residual": float(max(green)),
            "the_rank_one_control_defect": float(np.abs(
                cprod - cprod.conjugate().T).max()),
            "the_finite_throat_is_conservative": bool(
                max(r["hermiticity_defect"] for r in rows) < 1e-12
                and all(r["rank"] == 2 for r in rows)),
            "the_control_is_not": bool(np.abs(
                cprod - cprod.conjugate().T).max() > 1e-6),
            "what_it_checks": ("maximality and BC† = CB† at every λ, the DtN "
                               "against its own interior, and the rank-one "
                               "transfer model as the failing control")}


def measure_the_throat_transmits_at_the_traversal_time(
        lengths: Sequence[float] = (0.4, 0.6, 0.9, 1.2),
        long_lengths: Sequence[float] = (2.0, 3.0),
        area: float = FOUR_PI, separation: float = 1.3,
        width: float = 0.03) -> Dict[str, object]:
    """**The result.**  The throat's transmission is delayed by ``L``.

    The measured object is the **two-mouth block's** impulse response, not a
    field at some observer: ``r_ij(t)``, the inverse transform of
    ``R(ω) = (A(ω) − Γ(ω))⁻¹``.  That removes the source and observer geometry
    from the answer entirely.  It does *not* isolate the throat: ``Γ`` is the
    ambient's own mouth-to-mouth propagator and stays in, which is the whole
    reason the second prediction below reads ``min(L, d)`` rather than ``L``.

    Two predictions, and they are different:

    * ``r₁₁`` — same mouth in and out — starts at ``t = 0``.  A wave that
      reaches a mouth is partly **reflected instantaneously**; the tube's echoes
      arrive later, at ``2L, 4L, …``.  The point throat has the reflection and
      not the echoes.
    * ``r₁₂`` — opposite mouths — starts at ``min(L, d)``.  The tube's own path
      takes the traversal time ``L``; the *ambient* also connects the two mouths,
      along a geodesic of length ``d``, and that path is there whether or not the
      mouths are joined.  This is PR #258's cross-mouth channel and PR #259's
      ``β = 0`` control, now **separated in time** rather than by rank counting:
      for ``L < d`` the throat arrives first, for ``L > d`` the ambient does and
      the onset stops moving.

    The probe pulse's tail puts every onset early by the same fixed amount, so
    the quoted number is the **slope** ``d(onset)/dL``: ``1`` below ``d``,
    ``0`` above it.
    """
    rows = []
    for lam_len in list(lengths) + list(long_lengths):
        throat = FiniteThroat(separation=separation, length=float(lam_len),
                              area=area)
        grid = _grid_for(throat)
        imp = impulse_response(throat, grid, width=width)
        rows.append({
            "length": float(lam_len),
            "sigma_star": float(throat.growing_modes()["symmetric"]),
            "contour": float(grid.eps),
            "onset_same_mouth": causal_onset(imp["same_mouth"], imp["times"]),
            "onset_opposite": causal_onset(imp["opposite_mouths"],
                                           imp["times"]),
            "prediction": float(min(lam_len, separation))})
    short = [r for r in rows if r["length"] < separation]
    long = [r for r in rows if r["length"] > separation]
    slope = float(np.polyfit([r["length"] for r in short],
                             [r["onset_opposite"] for r in short], 1)[0])
    long_spread = float(max(r["onset_opposite"] for r in long)
                        - min(r["onset_opposite"] for r in long))
    # the point-throat control: matched at λ = 1, transmission is instantaneous
    ref = FiniteThroat(separation=separation, length=0.9, area=area)
    grid = _grid_for(ref)
    frozen = ref.boundary_matrix(complex(1.0)).real.astype(complex)
    pt = impulse_response(ref, grid, width=width, constant=frozen)
    return {"rows": rows,
            "slope_below_the_ambient_path": slope,
            "onset_spread_above_it": long_spread,
            "point_throat_onset_same_mouth": causal_onset(pt["same_mouth"],
                                                          pt["times"]),
            "point_throat_onset_opposite": causal_onset(pt["opposite_mouths"],
                                                        pt["times"]),
            "reflection_is_instantaneous": bool(
                max(r["onset_same_mouth"] for r in rows) < 4.0 * width),
            "the_delay_is_the_traversal_time": bool(abs(slope - 1.0) < 0.05),
            "the_ambient_path_takes_over": bool(long_spread < 1e-9),
            "the_point_throat_transmits_instantly": bool(
                causal_onset(pt["opposite_mouths"], pt["times"])
                < 4.0 * width),
            "the_two_channels": ("the tube arrives at L, the ambient at d; a "
                                 "point throat arrives at 0, which is what a "
                                 "point throat is")}


def measure_the_delay_ledger_is_the_bounce_series(
        throat: FiniteThroat = WORKING_THROAT, order: int = 400,
        eps: float = 1.6) -> Dict[str, object]:
    """The delays are a derivation, not a fit.

    On the retarded contour ``Im x > 0``, so the geometric series converge:

        ``cot x = −i − 2i Σ_{k≥1} e^{2ikx}`` ,  ``csc x = −2i Σ_{k≥0} e^{i(2k+1)x}``

    with ``x = kL``.  Reading the delays off: the same-mouth entry carries
    ``0, 2L, 4L, …`` and the cross-mouth entry ``L, 3L, 5L, …``.  **The
    parities are the physics** — an even number of traversals returns to the
    mouth it came from — and this is the reflected channel the rank-one
    transfer model does not have at all.

    Checked against the closed forms on the actual contour, not asserted.
    """
    grid = RetardedGrid(n=1 << 12, span=60.0, eps=float(eps))
    om = grid.omegas
    x = om * throat.length
    cot_s = -1j - 2j * sum(np.exp(2j * i * x) for i in range(1, int(order)))
    csc_s = -2j * sum(np.exp(1j * (2 * i + 1) * x) for i in range(int(order)))
    return {"contour": float(eps), "order": int(order),
            "cot_series_error": float(np.abs(1.0 / np.tan(x) - cot_s).max()),
            "csc_series_error": float(np.abs(1.0 / np.sin(x) - csc_s).max()),
            "delays": bounce_delays(4),
            "the_series_converge_on_the_contour": bool(
                max(float(np.abs(1.0 / np.tan(x) - cot_s).max()),
                    float(np.abs(1.0 / np.sin(x) - csc_s).max())) < 1e-11),
            "the_parity_rule": ("even multiples of L return to the mouth they "
                                "entered; odd multiples cross")}


def measure_the_short_tube_limit_is_a_mixed_stratum(
        throat: FiniteThroat = WORKING_THROAT, lam0: float = 1.0,
        factors: Sequence[float] = (1.0, 1.05, 1.2, 1.5, 2.0, 3.0)
        ) -> Dict[str, object]:
    """There **is** a point limit.  It is not a finite ``A``.

    A first draft of this round called the limit non-existent, on the evidence
    below: freezing ``A`` at ``A(λ₀)`` is exact at ``λ₀`` and nowhere else, and
    as ``L → 0`` the antisymmetric channel converges to ``−L/(2𝒜)`` while the
    symmetric one *diverges* like ``2/(𝒜λL)``.  Both facts are true and both are
    re-measured here.  The conclusion drawn from them was not.

    A boundary pair is defined up to ``(B, C) → (MB, MC)``, so a chart matrix
    running off to infinity means the limit has **left the chart**, not that it
    is absent.  Row-scaled — `normalized_pair` — the pair converges cleanly:

        ``(B, C)  ⟶  (P_anti, −P_sym)`` ,   i.e.   ``Φ_anti = 0`` , ``q_sym = 0``

    a **mixed Dirichlet–Neumann** stratum, maximal (``rank[B|C] = 2`` all the
    way) and self-adjoint, and reachable by no finite Hermitian ``A`` because
    both blocks are singular.  Exactly the stratum PR #257's review said the
    chart does not reach.  Convergence is linear in ``L``, measured.

    So the correct statement is **"no finite-``A`` point limit"**, and the
    constant-``A`` family of PRs #257–#259 is this throat read at one frequency
    — a band whose width in ``ω`` is ``∼ 1/L`` (in ``λ``, ``∼ 2√λ/L``; the two
    are not the same and the first draft mixed them).
    """
    frozen = throat.boundary_matrix(complex(lam0))
    band = []
    for f in factors:
        lam = float(lam0) * float(f)
        exact = throat.boundary_matrix(complex(lam))
        band.append({
            "lambda": lam,
            "beta_exact": float(exact[0, 1].real),
            "beta_frozen": float(frozen[0, 1].real),
            "relative_error": float(abs(exact[0, 1] - frozen[0, 1])
                                    / abs(exact[0, 1]))})
    short = []
    b_star, c_star = short_tube_stratum()
    for length in (0.4, 0.2, 0.1, 0.05, 0.02):
        t = FiniteThroat(separation=throat.separation, length=length,
                         area=throat.area)
        a = t.boundary_matrix(complex(lam0))
        anti = float((a[0, 0] - a[0, 1]).real)
        sym = float((a[0, 0] + a[0, 1]).real)
        limit = -length / (2.0 * t.area)
        b, c = t.normalized_pair(lam0)
        short.append({"length": float(length), "anti": anti,
                      "anti_limit": limit,
                      "anti_error": abs(anti - limit),
                      "anti_error_over_L_squared": abs(anti - limit)
                      / length ** 2,
                      "sym": sym,
                      "sym_prediction": 2.0 / (t.area * lam0 * length),
                      "distance_to_the_stratum": float(
                          max(np.abs(b - b_star).max(),
                              np.abs(c - c_star).max())),
                      "rank": int(np.linalg.matrix_rank(np.hstack([b, c]),
                                                        tol=1e-12))})
    for row in short:
        row["distance_over_L"] = row["distance_to_the_stratum"] / row["length"]
    rates = [r["distance_over_L"] for r in short]
    return {"lambda_0": float(lam0), "band": band, "short_tubes": short,
            "the_stratum": {"B": b_star.tolist(), "C": c_star.tolist(),
                            "condition": "Phi_anti = 0 and q_sym = 0"},
            "the_band_error_reaches_one": bool(
                band[-1]["relative_error"] > 0.5),
            "the_antisymmetric_channel_has_a_limit": bool(
                max(r["anti_error"] for r in short) < 1e-2
                and max(r["anti_error_over_L_squared"] for r in short) < 1.0),
            "the_chart_matrix_diverges": bool(
                short[-1]["sym"] > 4.0 * short[0]["sym"]),
            "worst_symmetric_prediction_error": float(max(
                abs(r["sym"] - r["sym_prediction"]) / abs(r["sym"])
                for r in short)),
            "distance_to_the_stratum": float(
                short[-1]["distance_to_the_stratum"]),
            "convergence_is_linear_in_L": bool(
                (max(rates) - min(rates)) / abs(np.mean(rates)) < 0.05),
            "every_pair_is_maximal": bool(all(r["rank"] == 2 for r in short)),
            "the_limit_exists_and_is_not_a_finite_A": bool(
                short[-1]["distance_to_the_stratum"] < 0.2
                and all(r["rank"] == 2 for r in short)
                and short[-1]["sym"] > 4.0 * short[0]["sym"]),
            "the_scope": ("the limit is a mixed Dirichlet-Neumann stratum, "
                          "not a finite A; the chart matrix diverges because "
                          "the limit leaves the chart, not because it is "
                          "absent")}


def measure_the_static_limit_is_rank_one_and_the_defect_diverges(
        throat: FiniteThroat = WORKING_THROAT,
        lams: Sequence[float] = (1e-8, 1e-6, 1e-4, 1e-2),
        masses: Sequence[float] = (0.05, 0.1, 0.2, 0.4)) -> Dict[str, object]:
    """PR #258's tomography, run on a throat with an interior — and it breaks.

    A massless tube's symmetric channel decouples at ``λ = 0``: a constant field
    on it costs nothing and carries no current, so ``A_sym → ∞`` and the static
    response ``S = Re R`` collapses onto the antisymmetric direction.  Measured:
    ``S ∝ [[1,−1],[−1,1]]`` and ``det S → 0`` **linearly in ``λ``**, so PR
    #258's defect ``𝒲 = S₁₂/det S − G₀`` diverges like ``1/λ``.

    **What that does and does not falsify.**  It is a falsifiable statement
    against the **generic finite-``A`` family** — every one of those has
    ``rank S = 2`` — and it is *not* a way to tell a finite throat from a point
    one.  The tube's own short-tube limit, the mixed stratum
    ``(P_anti, −P_sym)``, gives ``R → diag(0, −1/Γ_anti)`` in the channel basis:
    **rank one as well**, and the tube converges to it.  Both facts are measured
    here, because the first draft of this round claimed the stronger and wrong
    version.

    Give the tube an interior mass and the rank comes back, with
    ``det S ∝ (λ − m²)`` and the same coefficient.  And then the closure: off
    the collapse, ``𝒲 = −β(λ)`` **exactly**, with ``β`` the tube's own
    transmission amplitude.  PR #258's theorem survives the generalization; what
    it returns is no longer a constant but the interior's amplitude.
    """
    g0 = float(gamma_at(0.0, throat.separation)[0, 1].real)
    rows = []
    for lam in lams:
        s = np.linalg.inv(throat.krein_matrix(complex(lam))).real
        det = float(np.linalg.det(s))
        rows.append({"lambda": float(lam), "det_S": det,
                     "det_S_over_lambda": det / float(lam),
                     "antisymmetry": float(abs(s[0, 0] + s[0, 1])
                                           / abs(s[0, 0])),
                     "defect": float(s[0, 1] / det - g0),
                     "minus_beta": -float(throat.transmission(complex(lam)).real
                                          )})
    massive = []
    for m in masses:
        t = FiniteThroat(separation=throat.separation, length=throat.length,
                         area=throat.area, interior_mass=float(m))
        s = np.linalg.inv(t.krein_matrix(complex(0.0))).real
        det = float(np.linalg.det(s))
        massive.append({"mass": float(m), "det_S": det,
                        "det_S_over_mass_squared": det / m ** 2,
                        "defect": float(s[0, 1] / det - g0),
                        "minus_beta": -float(
                            t.transmission(complex(0.0)).real),
                        "defect_error": float(abs(
                            s[0, 1] / det - g0
                            + float(t.transmission(complex(0.0)).real)))})
    # the collapse is a *limit*: the claims are made where λ → 0, and the
    # largest λ is carried as the row that shows the correction turning on
    small = [r for r in rows if r["lambda"] <= 1e-4]
    coeff = [r["det_S_over_lambda"] for r in small]
    # the limiting stratum's own static response, for the scope statement
    gam = gamma_at(0.0, throat.separation).real
    gamma_anti = float(gam[0, 0] - gam[0, 1])
    stratum_r = _CHANNEL_BASIS @ np.diag([0.0, -1.0 / gamma_anti]) \
        @ _CHANNEL_BASIS.T
    approach = []
    for length in (0.4, 0.1, 0.02):
        t = FiniteThroat(separation=throat.separation, length=length,
                         area=throat.area, interior_mass=throat.interior_mass)
        s = np.linalg.inv(t.krein_matrix(complex(1e-9))).real
        approach.append({"length": float(length),
                         "distance_to_the_stratum_response": float(
                             np.abs(s - stratum_r).max()),
                         "det_S": float(np.linalg.det(s))})
    return {"rows": rows, "massive": massive,
            "stratum_response": stratum_r.tolist(),
            "stratum_rank": int(np.linalg.matrix_rank(stratum_r, tol=1e-12)),
            "approach_to_the_stratum": approach,
            "the_stratum_is_rank_one_too": bool(
                np.linalg.matrix_rank(stratum_r, tol=1e-12) == 1),
            "it_falsifies_the_finite_A_family_not_point_ness": bool(
                np.linalg.matrix_rank(stratum_r, tol=1e-12) == 1
                and approach[-1]["distance_to_the_stratum_response"]
                < approach[0]["distance_to_the_stratum_response"]),
            "det_S_is_linear_in_lambda": bool(
                (max(coeff) - min(coeff)) / abs(np.mean(coeff)) < 1e-3),
            "linear_coefficient": float(np.mean(coeff)),
            "departure_at_the_largest_lambda": float(
                abs(rows[-1]["det_S_over_lambda"] - np.mean(coeff))
                / abs(np.mean(coeff))),
            "worst_antisymmetry": float(max(r["antisymmetry"] for r in small)),
            "worst_defect_error_over_lambda": float(max(
                abs(r["defect"] - r["minus_beta"]) / abs(r["minus_beta"])
                for r in rows)),
            "the_static_response_is_rank_one": bool(
                max(r["antisymmetry"] for r in small) < 1e-4),
            "the_defect_diverges": bool(
                abs(rows[0]["defect"]) > 100.0 * abs(rows[-1]["defect"])),
            "worst_defect_error": float(max(r["defect_error"]
                                            for r in massive)),
            "the_defect_is_still_minus_beta": bool(
                max(r["defect_error"] for r in massive) < 1e-6),
            "the_falsifiable_difference": ("a generic finite-A throat is "
                                           "statically rank two; this tube is "
                                           "rank one and 𝒲 diverges — but so "
                                           "is its own point limit, so rank "
                                           "one falsifies the chart, not "
                                           "point-ness")}


def measure_the_interior_mass_is_a_transmission_cutoff(
        separation: float = 1.3, length: float = 3.0, area: float = FOUR_PI,
        mass: float = 2.0,
        lams: Sequence[float] = (0.0, 0.5, 1.5, 3.0, 3.9, 4.1, 6.0, 12.0)
        ) -> Dict[str, object]:
    """Below ``λ = m²`` the tube is evanescent, and the throat stops
    transmitting.

    ``k² = λ − m²`` turns imaginary below the cutoff, so

        ``β(λ)  =  csc(kL)/(𝒜k)  →  −csch(κL)/(𝒜κ)  ≈  −2e^{−κL}/(𝒜κ)``

    — **negative**, and exponentially suppressed rather than oscillating.  The
    interior has a **mass gap**, and below it the channel is *evanescent*: the
    two mouths look like two ordinary scatterers with a tunnelling correction,
    while the same throat above the gap is fully connected.  (An earlier draft
    called this "low-pass", which is backwards — it is the *low* frequencies
    that are suppressed.  It is an evanescent cutoff, not a filter passband.)  The sign flip is not
    cosmetic — it is the difference between a propagating and a tunnelling
    channel, and PR #258's defect reads it as ``𝒲 = −β`` either way.

    The cutoff is also where the **rank collapses**: ``k → 0`` makes
    ``β → 1/(𝒜k²L)`` diverge, which is the same zero-mode that empties the
    symmetric channel at ``λ = 0`` for a massless tube.  The two statements are
    one statement, read at ``λ = m²``.

    Both regimes stay real and self-adjoint: the square root is even in ``k``,
    so there is no branch cut and nothing to choose.
    """
    t = FiniteThroat(separation=separation, length=length, area=area,
                     interior_mass=float(mass))
    rows = []
    for lam in lams:
        beta = complex(t.transmission(complex(lam)))
        kappa = math.sqrt(max(mass ** 2 - lam, 0.0))
        asymptote = (-2.0 * math.exp(-kappa * length) / (area * kappa)
                     if kappa > 1e-9 else float("nan"))
        rows.append({"lambda": float(lam), "propagating": bool(lam > mass ** 2),
                     "beta": float(beta.real),
                     "imaginary_part": float(abs(beta.imag)),
                     "kappa": float(kappa),
                     "exponential_asymptote": float(asymptote)})
    deep = [r for r in rows if r["kappa"] * length > 4.0]
    high = [r for r in rows if r["propagating"] and r["lambda"] > 5.0]
    suppression = float(abs(deep[0]["beta"]) / abs(high[0]["beta"]))
    # the discriminator is monotone decay against oscillation, not the sign
    # itself: above the cutoff β changes sign every time kL passes a multiple
    # of π, and a single high-λ sample can land on either branch
    below = np.array([float(t.transmission(complex(x)).real)
                      for x in np.linspace(0.0, mass ** 2 - 1e-3, 400)])
    above = np.array([float(t.transmission(complex(x)).real)
                      for x in np.linspace(mass ** 2 + 1e-3, 60.0, 400)])
    return {"mass": float(mass), "cutoff": float(mass ** 2),
            "length": float(length), "rows": rows,
            "worst_imaginary_part": float(max(r["imaginary_part"]
                                              for r in rows)),
            "asymptote_error": float(max(
                abs(r["beta"] - r["exponential_asymptote"]) / abs(r["beta"])
                for r in deep)),
            "deepest_kappa_L": float(max(r["kappa"] for r in rows) * length),
            "suppression_ratio": suppression,
            "sign_changes_below": int(np.sum(np.diff(np.sign(below)) != 0)),
            "sign_changes_above": int(np.sum(np.diff(np.sign(above)) != 0)),
            "below_is_monotone": bool(np.all(np.diff(np.abs(below)) > 0.0)),
            "the_evanescent_side_does_not_oscillate": bool(
                int(np.sum(np.diff(np.sign(below)) != 0)) == 0
                and int(np.sum(np.diff(np.sign(above)) != 0)) > 1),
            "the_transmission_is_suppressed_below_cutoff": bool(
                suppression < 0.05),
            "everything_stays_real": bool(
                max(r["imaginary_part"] for r in rows) < 1e-14),
            "what_it_means": ("a massive interior gives the channel a mass "
                              "gap; below m² transmission is evanescent and "
                              "the mouths look like two ordinary scatterers "
                              "with a tunnelling term")}


def measure_the_growing_mode_belongs_to_the_mouth(
        areas: Sequence[float] = (0.2, 0.5, 1.0),
        lengths: Sequence[float] = (1.5, 3.0, 6.0),
        separations: Sequence[float] = (0.8, 1.3, 2.4, 3.0)
        ) -> Dict[str, object]:
    """**The model fails the stability gate** — and the failure is the
    *mouth's*, not the tube's.

    ``A(λ)`` decreases in ``λ`` and ``Γ(λ)`` increases (PR #257's Gram
    identity), so ``A − Γ`` is strictly monotone and each channel has at most
    one root per pole-interval.  In the symmetric channel there is always
    exactly one, at ``λ < 0``: the tube's zero mode, pushed below zero by
    coupling to a point.

    Three facts identify what it is, and all three are **limits** — they hold as
    ``σ*L`` and ``σ*d`` grow, and the convergence is measured rather than the
    conclusion asserted:

    * its rate tends to ``σ* = 2√(π/𝒜)`` — closed form, **no ``L``**;
    * its dependence on the mouth separation ``d`` dies exponentially: at
      ``σ*d ≈ 19`` two separations agree to ``4e-9``;
    * the symmetric and antisymmetric channels become **degenerate**, and the
      splitting is not merely small but *is* the Euclidean mouth-to-mouth
      propagator: ``σ_sym − σ_anti = 1.04·e^{−σ*d}`` across a decade, which is
      the mechanism rather than a bound on it.

    A mode that ignores the tube's length and the mouths' separation, and does
    not distinguish the two channels, is a single-mouth object.  Its scale is
    ``1/σ* = √(𝒜/4π)`` — the mouth's own radius, which is exactly where "point
    mouth" stops being an approximation.  Every statement in this module is made
    at ``|λ| ≪ σ*²``, and that band is reported here rather than assumed.

    Both approximations fail in the same direction and for the same reason: at
    ``σ*L ≲ 2`` the tube is shorter than the mode, so folding it about its
    middle is not yet ``coth → 1``, and the closed form is ``15%`` out.  That
    row is kept in the table.
    """
    rows = []
    for area in areas:
        for length in lengths:
            t = FiniteThroat(separation=1.3, length=length, area=area)
            m = t.growing_modes()
            rows.append({"area": float(area), "length": float(length),
                         "sigma_star": float(m["symmetric"]),
                         "closed_form": float(m["closed_form"]),
                         "relative_error": abs(float(m["symmetric"])
                                               - float(m["closed_form"]))
                         / float(m["closed_form"]),
                         "channel_split": (float(m["channel_split"])
                                           if m["channel_split"] is not None
                                           else None),
                         "split_over_exponential": (
                             float(m["channel_split"])
                             / math.exp(-float(m["symmetric"]) * 1.3)
                             if m["channel_split"] is not None else None),
                         "sigma_times_L": float(m["symmetric"]) * length})
    by_sep = []
    for d in separations:
        t = FiniteThroat(separation=d, length=3.0, area=min(areas))
        s = float(t.growing_modes()["symmetric"])
        by_sep.append({"separation": float(d), "sigma_star": s,
                       "sigma_times_d": s * float(d)})
    spread = float(max(r["sigma_star"] for r in by_sep)
                   - min(r["sigma_star"] for r in by_sep))
    far = [r for r in by_sep if r["sigma_times_d"] > 15.0]
    far_spread = float(max(r["sigma_star"] for r in far)
                       - min(r["sigma_star"] for r in far))
    resolved = [r for r in rows if r["sigma_times_L"] > 6.0]
    splits = [r["channel_split"] for r in resolved
              if r["channel_split"] is not None]
    # the splitting has two exponentially small sources — the tube, through
    # coth − tanh = 2/sinh(κL), and the ambient, through G_d — and only past
    # σ*L ≈ 14 is the tube's gone and the ambient's law clean
    deep = [r for r in rows if r["sigma_times_L"] > 14.0]
    ratios = [r["split_over_exponential"] for r in deep
              if r["split_over_exponential"] is not None]
    work = WORKING_THROAT.growing_modes()
    return {"rows": rows, "by_separation": by_sep,
            "separation_spread": spread,
            "separation_spread_far": far_spread,
            "worst_closed_form_error": float(max(r["relative_error"]
                                                 for r in resolved)),
            "worst_channel_split": float(max(splits)) if splits else None,
            "split_over_exponential": ratios,
            "every_throat_has_one": bool(all(r["sigma_star"] > 0.0
                                             for r in rows)),
            "it_stops_knowing_the_separation": bool(far_spread < 1e-7),
            "the_closed_form_holds_once_sigma_L_is_large": bool(
                max(r["relative_error"] for r in resolved) < 2e-3),
            "the_split_is_the_euclidean_propagator": bool(
                max(ratios) / min(ratios) < 1.1 if ratios else False),
            "the_working_band": float(work["symmetric"]) ** 2,
            "the_diagnosis": ("a mode that ignores L and d and does not split "
                              "the channels is a single-mouth object, so the "
                              "instability belongs to the point-mouth matching "
                              "and not to the interior"),
            "the_open_question": ("whether a finite-radius mouth or neck "
                                  "geometry removes it — unresolved here, and "
                                  "the thing to settle before stationary "
                                  "action or backreaction")}


def measure_the_contour_must_clear_the_growing_mode(
        throat: FiniteThroat = FiniteThroat(separation=1.3, length=0.6,
                                            area=FOUR_PI),
        clearances: Sequence[float] = (-0.03, 0.02, 0.3, 0.8),
        width: float = 0.03) -> Dict[str, object]:
    """The numerical edge, with its failure mode shown rather than described.

    The retarded contour ``Im ω = ε`` must lie **above every singularity of the
    response**, and a finite throat puts one at ``ω = iσ*``.  Place ``ε`` just
    below it and the inversion returns a field with support *before its own
    light cone*: a constant pedestal, at the level of the pole's residue divided
    by the frequency span, arriving at ``t = 0`` for an event that cannot begin
    until ``t = L``.

    It is the same species of error as PR #259's under-resolved contour — a
    plausible-looking number produced by a contour in the wrong place — and it
    is reported the same way: both values, and the rule.  ``σ*`` has a closed
    form, so the rule is checkable before the solve rather than after it.

    **What clearing the contour does not do is stabilize anything.**  Above
    ``σ*`` the inversion returns the correct retarded solution *of an unstable
    system*; the growing mode is still there, and every time series in this
    module is dominated by it at late times.  The delay is read from the causal
    **onset**, which is a statement about the analytic structure and is immune
    to what the solution does afterwards — but the instability is a property of
    the model, not of the contour, and moving the contour never touches it.
    """
    sig = float(throat.growing_modes()["symmetric"])
    onset_true = min(throat.length, throat.separation)
    rows = []
    for c in clearances:
        grid = RetardedGrid(n=1 << 17, span=300.0, eps=sig + float(c))
        imp = impulse_response(throat, grid, width=width)
        series = imp["opposite_mouths"]
        early = series[imp["times"] < 0.5 * onset_true]
        rows.append({
            "clearance": float(c), "contour": float(grid.eps),
            "above_the_mode": bool(c > 0.0),
            "onset": causal_onset(series, imp["times"]),
            "pedestal": float(np.abs(early).max() / np.abs(series).max())})
    below = [r for r in rows if not r["above_the_mode"]]
    above = [r for r in rows if r["above_the_mode"] and r["clearance"] > 0.1]
    return {"sigma_star": sig, "true_onset": float(onset_true), "rows": rows,
            "pedestal_below": float(max(r["pedestal"] for r in below)),
            "pedestal_above": float(max(r["pedestal"] for r in above)),
            "onset_below": float(min(r["onset"] for r in below)),
            "onset_above": float(min(r["onset"] for r in above)),
            "a_contour_below_the_mode_breaks_causality": bool(
                max(r["pedestal"] for r in below)
                > 1e4 * max(r["pedestal"] for r in above)),
            "the_rule": ("ε must exceed σ*, which has a closed form, so the "
                         "contour can be placed before the solve rather than "
                         "diagnosed after it"),
            "what_it_does_not_do": ("clearing the contour evaluates the "
                                    "correct retarded solution of an unstable "
                                    "system; it does not cure the instability")}
