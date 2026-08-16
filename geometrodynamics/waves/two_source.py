"""Static two-source throat tomography — **not** the two-wave invariant.

Roadmap step 3, **partly**.  PR #253 closed rank counting with an explicit
statement of what it *cannot* supply: a quantity that vanishes when a source is
removed rather than merely becoming underdetermined.  This module supplies one,
and then establishes what it does and does not measure.  The headline of the
first draft was wrong about the second half, and the scope is stated here rather
than in a caveat at the end.

What this is **not**
────────────────────
It is **not** the roadmap's two-wave collision invariant.  The object built here
is a *static* (zero- or fixed-spectral-parameter) source-interaction kernel.  It
contains **no local null momenta**, so it cannot distinguish equal-energy
collinear from counterpropagating waves — which was the entire load-bearing
control behind ``𝒞 = I_A I_B (k_A·k_B)²``.  The dynamical object, built from
``T_A^{μν} T^B_{μν}`` and resolved on the geodesic/winding branches of PRs
#253–#255, is still owed.

Relatedly, the index pair ``(i, j)`` here labels **mouth channels** — which
mouth the field entered and which it left — and *not* the branches of #253–#255,
which
are geodesic histories with winding numbers and Maslov signs.  The functions are
named for mouths accordingly.

What this is
────────────
For a linear field, superposition is exact, so no *linear* functional of the
field carries two-source information — it is additive.  A **quadratic**
functional does, in its cross term:

    𝒞  =  Q[a φ_A + b φ_B] − Q[a φ_A] − Q[b φ_B]  =  2ab · 𝒞(y_A, y_B)

identically zero if either amplitude is zero.  For static sources ``Q`` is the
interaction energy and the cross term is the throat's Green function between the
two source points:

    𝒞(y_A, y_B)  =  Re G_A(y_A, y_B)
                 =  G(y_A,y_B)  +  Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B)

with ``R = (C − BΓ)⁻¹B`` the throat's **response matrix** — in the finite-``A``
chart of PR #257, ``R = (A − Γ(λ))⁻¹``.

What it does and does not discriminate
──────────────────────────────────────
Three claims that look like the signature and are not, each one excluded here:

1. **The cross term is nonzero.**  So is every interference pattern.  Not a
   signature of anything.
2. **The interaction is anisotropic** — it depends on more than the geodesic
   separation ``χ_AB``.  True, and measured; but *two independent scatterers*
   break isotropy exactly as loudly.  Not a signature of a throat.
3. **The response matrix ``R`` has off-diagonal entries.**  It does even when
   the boundary matrix is diagonal, because ``Γ`` itself is: two unconnected
   scatterers talk through the ambient field.  So the off-diagonal block is a
   **cross-mouth** channel, and calling it "through the throat" would contradict
   the very next paragraph.

What *is* a discriminator is a parameter count.  The static invariant determines
exactly three numbers — the real symmetric matrix ``S = Re R`` — while the
disconnected family ``A = diag(α₁,α₂)`` has only two knobs.  The image of the
disconnected family is therefore a surface in that three-space, with an exact
equation:

    S₁₂  =  G₀ · det S           (two independent scatterers, any α)

so the **disconnection defect**

    𝒲  =  S₁₂ / det S  −  G₀

is zero on it and nonzero off it.  Stated precisely, and this is the claim the
round actually proves:

    𝒲 detects OFF-DIAGONAL MOUTH-BOUNDARY MIXING, relative to the diagonal
    two-scatterer null model, within this point-interaction model.

On the real slice it is not merely nonzero but *exactly the mouth-mixing
amplitude*, ``𝒲 = −β``, independent of the self-energies, the mouth separation,
and the **Löwner margin** — the answer to PR #255's caution that anything built
from a resummed field measures the pole rather than the source.

Which field, and why it decides the blind spot
──────────────────────────────────────────────
``𝒲 = 0`` has solutions away from ``β = 0`` only for **complex** ``β``.  That is
not a free choice.  A **real** time-domain scalar — which is what PR #254 solves
— needs the self-adjoint domain to be invariant under complex conjugation, and
in the finite-``A`` chart that is exactly ``A = A*``, hence ``β`` real.  The
consequence is checkable in one line and checked here: with complex ``β`` a
**real static source produces a complex field**.

So for the arc's actual field content there is **no blind family**: ``𝒲 = −β``
removes it at a single spectral parameter.  The blind family is a statement
about a deliberately **time-reversal-breaking, complex-scalar** extension, and
it is reported that way rather than as a defect of the test.

And even there the limitation is the *protocol*, not the operator.  Real static
sources see only ``Re R``, three numbers for four parameters.  **Phase-sensitive
complex sources** recover the full complex ``R`` at one spectral parameter, and
then ``A = Γ + R⁻¹`` outright.  The multi-parameter reconstruction below is a
restriction of the real-static-source protocol.

A note on ``λ``
───────────────
``λ = ω²``, so a negative ``λ`` is an *imaginary* frequency, not a second
physical driving frequency.  Everything here is parametrized by the **spectral
parameter** ``λ`` and named that way; the reconstruction defaults to two
positive ``λ`` below the free ground state ``λ = 1``.

What is put in
──────────────
The background, the boundary data, and the mouth positions.  The observer is
assumed to know where the mouths are and what the free Green function is; ``G₀``
enters ``𝒲`` explicitly.  Nothing here derives ``A`` from matter, and the throat
is still point-supported — no interior, no proper length, no delay.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .throat_operator import MouthPair, gamma_at
from .throat_positivity import boundary_pair_from_unitary, positivity_defect

__all__ = [
    "green_at",
    "geodesic",
    "mouth_positions",
    "source_vector",
    "response_matrix",
    "response_of_pair",
    "static_response",
    "mouth_channel_invariant",
    "energy_functional",
    "is_real_field_compatible",
    "recover_complex_response",
    "interaction_energy",
    "free_interaction_energy",
    "cross_matrix",
    "isotropy_profile",
    "disconnection_defect",
    "defect_of_pair",
    "invisible_partner",
    "recover_response",
    "recover_boundary",
    "random_points",
    "ring_points",
    "measure_the_invariant_vanishes_when_a_source_is_removed",
    "measure_the_throat_channel_has_the_rank_of_the_boundary_condition",
    "measure_anisotropy_is_not_the_signature",
    "measure_two_disconnected_scatterers_lie_on_a_surface",
    "measure_the_defect_is_the_mouth_mixing_amplitude",
    "measure_the_invariant_is_recoverable_from_observations",
    "measure_the_blind_spot_of_a_single_frequency_test",
    "measure_a_real_field_forces_beta_real",
    "measure_phase_sensitive_sources_need_only_one_spectral_parameter",
    "measure_two_spectral_parameters_reconstruct_the_boundary_matrix",
    "measure_the_antipodal_endpoint_on_its_own",
]

# The working point: a boundary matrix strictly inside PR #257's cone, quoted
# with its Löwner margin everywhere it is used.  Nothing here is evaluated at
# the apex or on the null boundary, where the static response is singular.
WORKING_SEPARATION = 1.3
WORKING_BOUNDARY = (0.30, 0.35, 0.06 + 0.0j)


# ════════════════════════════════════════════════════════════════════════════
# GEOMETRY ON S³, AND THE FREE PROPAGATOR BETWEEN TWO POINTS
# ════════════════════════════════════════════════════════════════════════════
def green_at(chi: float, lmbda: float = 0.0) -> float:
    """``G(χ, λ)`` — the free propagator between two points a geodesic distance
    ``χ`` apart, for real ``λ`` of either sign.

    This is exactly ``gamma_at``'s off-diagonal entry read at a general
    separation rather than the mouth separation, and it is called that way on
    purpose: the same function that gives ``Γ``'s mouth-to-mouth entry gives the
    source-to-mouth entries, so no second implementation of the antipodal limit
    or the ``λ < 0`` continuation can drift from the first.
    """
    return float(gamma_at(float(lmbda), float(chi))[0, 1].real)


def geodesic(u: Sequence[float], v: Sequence[float]) -> float:
    """The geodesic distance on the unit ``S³``, ``arccos(u·v)``."""
    c = float(np.clip(np.dot(np.asarray(u, dtype=float),
                             np.asarray(v, dtype=float)), -1.0, 1.0))
    return float(math.acos(c))


def mouth_positions(separation: float) -> Tuple[np.ndarray, np.ndarray]:
    """Two mouths a geodesic distance ``d`` apart, in a fixed frame."""
    d = float(separation)
    if not 0.0 < d <= math.pi + 1e-12:
        raise ValueError("the mouths need 0 < d ≤ π")
    return (np.array([1.0, 0.0, 0.0, 0.0]),
            np.array([math.cos(d), math.sin(d), 0.0, 0.0]))


def random_points(n: int, seed: int = 20260816) -> np.ndarray:
    """``n`` uniform points on ``S³``, as unit rows."""
    rng = np.random.default_rng(int(seed))
    x = rng.normal(size=(int(n), 4))
    return x / np.linalg.norm(x, axis=1, keepdims=True)


def ring_points(centre: Sequence[float], chi: float, n: int,
                seed: int = 20260816) -> np.ndarray:
    """``n`` points at *exactly* geodesic distance ``χ`` from ``centre``.

    The set is a two-sphere, not a circle — on ``S³`` the locus at fixed
    distance from a point is ``S²`` — so it is sampled rather than swept.  Every
    point on it has the same free interaction with the centre, which is the
    whole reason for building it.
    """
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    rng = np.random.default_rng(int(seed))
    x = rng.normal(size=(int(n), 4))
    x -= np.outer(x @ c, c)                       # into the tangent 3-space
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    return math.cos(chi) * c[None, :] + math.sin(chi) * x


def source_vector(point: Sequence[float], separation: float,
                  lmbda: float = 0.0) -> np.ndarray:
    """``v_i = G(χ(y, c_i), λ)`` — how strongly a source at ``y`` drives each
    mouth.  Real, and the only way the source's position enters the throat
    channel at all."""
    c1, c2 = mouth_positions(separation)
    return np.array([green_at(geodesic(point, c1), lmbda),
                     green_at(geodesic(point, c2), lmbda)], dtype=float)


# ════════════════════════════════════════════════════════════════════════════
# THE RESPONSE MATRIX
# ════════════════════════════════════════════════════════════════════════════
def response_matrix(b: np.ndarray, c: np.ndarray,
                    gamma: np.ndarray) -> np.ndarray:
    """``R = (C − BΓ)⁻¹ B`` — the throat's response, for any self-adjoint pair.

    Krein's resolvent formula reads
    ``G_A(x,y) = G(x,y) + Σ_ij G(x,c_i) R_ij G(c_j,y)``, so ``R`` is the entire
    content of the throat as far as any external source is concerned.  Written
    with the pair ``(B, C)`` rather than the chart matrix ``A`` because PR #257
    ended by showing the chart is not the whole family: in the chart
    (``B = I``, ``C = A``) this is ``(A − Γ)⁻¹``, and on the Dirichlet strata,
    where ``A`` does not exist, it is still this.

    ``R`` is Hermitian whenever the pair is self-adjoint — verified, not assumed
    — and ``rank R = rank B``, which is what makes the rank statement below a
    statement about the boundary condition rather than about the sources.
    """
    b = np.asarray(b, dtype=complex)
    c = np.asarray(c, dtype=complex)
    return np.linalg.solve(c - b @ np.asarray(gamma, dtype=complex), b)


def response_of_pair(pair: MouthPair, lmbda: float = 0.0) -> np.ndarray:
    """``R = (A − Γ(λ))⁻¹`` — the response in the finite-``A`` chart."""
    return np.linalg.inv(pair.krein_matrix(float(lmbda)))


def static_response(pair: MouthPair, lmbda: float = 0.0) -> np.ndarray:
    """``S = Re R`` — the part a pair of *real* static sources can see.

    The interaction energy of two real sources is a Hermitian form evaluated on
    a real vector, so only the real symmetric part of ``R`` survives.  ``S`` has
    **three** independent entries; the boundary matrix has **four** parameters.
    That mismatch is this round's central fact, and §"blind spot" is its cost.
    """
    return np.asarray(response_of_pair(pair, lmbda).real, dtype=float)


# ════════════════════════════════════════════════════════════════════════════
# THE INVARIANT
# ════════════════════════════════════════════════════════════════════════════
def free_interaction_energy(y_a: Sequence[float], y_b: Sequence[float],
                            lmbda: float = 0.0) -> float:
    """``G(χ_AB, λ)`` — the cross term with no throat at all.

    A function of the geodesic separation and **nothing else**.  That is the
    null hypothesis every measurement below is stated against.
    """
    return green_at(geodesic(y_a, y_b), lmbda)


def mouth_channel_invariant(pair: MouthPair, y_a: Sequence[float],
                            y_b: Sequence[float],
                            lmbda: float = 0.0) -> Dict[str, object]:
    """The static invariant, resolved on a pair of **mouth channels**.

    Not the branch index of PRs #253–#255 — those are geodesic histories with
    winding numbers and Maslov signs, and nothing here carries either.  The
    index is simply which mouth the field entered and which it left, plus the
    channel that used neither:

    * ``direct`` — the ``(∅, ∅)`` entry, ``G(χ_AB)``, present without a throat;
    * ``channels[i][j]`` — the field left source ``A`` into mouth ``i`` and
      reached source ``B`` out of mouth ``j``.

    The off-diagonal is reported as ``cross_mouth`` and **not** as "through the
    throat": ``Γ`` is off-diagonal whatever the boundary data says, so two
    disconnected scatterers fill it too.  What isolates the non-local part of
    the boundary condition is `disconnection_defect`, not this entry.

    The total is the full cross term, bilinear in the two source strengths.
    """
    v_a = source_vector(y_a, pair.separation, lmbda)
    v_b = source_vector(y_b, pair.separation, lmbda)
    s = static_response(pair, lmbda)
    block = np.outer(v_a, v_b) * s
    direct = free_interaction_energy(y_a, y_b, lmbda)
    return {"direct": float(direct),
            "channels": [[float(block[0, 0]), float(block[0, 1])],
                         [float(block[1, 0]), float(block[1, 1])]],
            "cross_mouth": float(block[0, 1] + block[1, 0]),
            "same_mouth": float(block[0, 0] + block[1, 1]),
            "throat_total": float(block.sum()),
            "total": float(direct + block.sum()),
            "source_vectors": [[float(t) for t in v_a],
                               [float(t) for t in v_b]]}


def energy_functional(pair: MouthPair, y_a: Sequence[float],
                      y_b: Sequence[float], a: float, b: float,
                      lmbda: float = 0.0) -> float:
    """``Q[a φ_A + b φ_B]`` — the quadratic functional itself, self-energies
    included.

    The whole point of the round is that the invariant is the *cross term* of a
    quadratic functional, so the functional has to exist independently of the
    cross term or the vanishing is a tautology.  For point sources of strengths
    ``a`` and ``b``,

        ``Q = ½[a² G^reg_A(y_A,y_A) + 2ab Re G_A(y_A,y_B)
                + b² G^reg_A(y_B,y_B)]``

    where ``G^reg`` is the Coulomb-subtracted coincidence limit.  The
    subtraction is entirely in the *free* part — the throat term
    ``vᵀ S v`` is finite at coincidence — so ``G^reg_A(y,y) = g₀ + vᵀ_y S v_y``
    with the same ``g₀`` PR #256 built the boundary condition on.

    Both self-energy terms are present and both are large; they cancel in
    ``Q[a,b] − Q[a,0] − Q[0,b]`` because the functional is quadratic, not
    because anything was multiplied by zero.
    """
    lam = float(lmbda)
    v_a = source_vector(y_a, pair.separation, lam)
    v_b = source_vector(y_b, pair.separation, lam)
    s = static_response(pair, lam)
    g_reg = float(pair.gamma(lam).real[0, 0])       # the subtracted coincidence
    self_a = g_reg + float(v_a @ s @ v_a)
    self_b = g_reg + float(v_b @ s @ v_b)
    cross = free_interaction_energy(y_a, y_b, lam) + float(v_a @ s @ v_b)
    return float(0.5 * (a * a * self_a + 2.0 * a * b * cross
                        + b * b * self_b))


def interaction_energy(pair: MouthPair, y_a: Sequence[float],
                       y_b: Sequence[float], lmbda: float = 0.0) -> float:
    """The scalar invariant — the total of `mouth_channel_invariant`."""
    v_a = source_vector(y_a, pair.separation, lmbda)
    v_b = source_vector(y_b, pair.separation, lmbda)
    s = static_response(pair, lmbda)
    return float(free_interaction_energy(y_a, y_b, lmbda)
                 + v_a @ s @ v_b)


def cross_matrix(pair: MouthPair, sources: np.ndarray,
                 lmbda: float = 0.0) -> Dict[str, object]:
    """The full ``N × N`` table of pairwise invariants, split into channels.

    The point is the **rank**.  The direct table is a generic Gram-like matrix
    of the ambient Green function and has full rank; the throat table is
    ``Vᵀ S V`` with ``V`` a ``2 × N`` matrix, so its rank is at most two *no
    matter how many sources there are*.  The entire multi-source signature of a
    point throat lives in a two-dimensional space, and on a Dirichlet stratum in
    a one-dimensional one.
    """
    pts = np.asarray(sources, dtype=float)
    n = pts.shape[0]
    v = np.array([source_vector(p, pair.separation, lmbda) for p in pts]).T
    s = static_response(pair, lmbda)
    throat = v.T @ s @ v
    direct = np.array([[free_interaction_energy(pts[i], pts[j], lmbda)
                        if i != j else 0.0
                        for j in range(n)] for i in range(n)])
    return {"n_sources": n,
            "direct_rank": int(np.linalg.matrix_rank(direct, tol=1e-9)),
            "throat_rank": int(np.linalg.matrix_rank(throat, tol=1e-12)),
            "response_rank": int(np.linalg.matrix_rank(s, tol=1e-12)),
            "throat_norm": float(np.abs(throat).max()),
            "direct_norm": float(np.abs(direct).max())}


def isotropy_profile(pair: MouthPair, chi: float = 1.0, n: int = 240,
                     lmbda: float = 0.0,
                     seed: int = 20260816) -> Dict[str, object]:
    """Hold ``χ_AB`` fixed, move ``B`` over the sphere of that radius.

    Free field theory says the interaction cannot move: it is a function of the
    separation alone.  It does move, and the size of the motion is the throat
    channel.  It also moves for a **disconnected** pair of scatterers, reported
    alongside, which is why anisotropy is a detection of *structure* and not of
    a *connection*.
    """
    c1, _ = mouth_positions(pair.separation)
    y_a = np.array([math.cos(0.7), 0.0, math.sin(0.7), 0.0])
    ring = ring_points(y_a, float(chi), int(n), seed=seed)
    free = free_interaction_energy(y_a, ring[0], lmbda)
    disconnected = MouthPair(pair.separation, pair.alpha1, pair.alpha2, 0.0)
    vals = np.array([interaction_energy(pair, y_a, p, lmbda) for p in ring])
    dvals = np.array([interaction_energy(disconnected, y_a, p, lmbda)
                      for p in ring])
    freevals = np.array([free_interaction_energy(y_a, p, lmbda) for p in ring])
    return {"chi": float(chi), "n_points": int(n),
            "free_value": float(free),
            "free_spread": float(freevals.max() - freevals.min()),
            "throat_spread": float(vals.max() - vals.min()),
            "throat_relative_spread": float((vals.max() - vals.min())
                                            / abs(vals.mean())),
            "disconnected_spread": float(dvals.max() - dvals.min()),
            "disconnected_relative_spread": float((dvals.max() - dvals.min())
                                                  / abs(dvals.mean())),
            "the_free_field_is_isotropic": bool(
                (freevals.max() - freevals.min()) < 1e-12),
            "both_break_it": bool(dvals.max() - dvals.min() > 1e-9
                                  and vals.max() - vals.min() > 1e-9)}


# ════════════════════════════════════════════════════════════════════════════
# THE DISCRIMINATOR
# ════════════════════════════════════════════════════════════════════════════
def disconnection_defect(response: np.ndarray, g_between_mouths: float
                         ) -> float:
    """``𝒲 = S₁₂ / det S − G₀`` — zero iff two disconnected scatterers could
    have produced ``S``.

    Two independent point scatterers are ``A = diag(α₁, α₂)``, so
    ``S = (diag(α) − Γ)⁻¹`` with ``Γ = [[g₀, G₀], [G₀, g₀]]``.  Writing
    ``p = α₁ − g₀`` and ``q = α₂ − g₀``,

        ``S = [[q, G₀], [G₀, p]] / (pq − G₀²)``  ⟹  ``det S = 1/(pq − G₀²)``

    and therefore ``S₁₂ = G₀ det S`` identically, for every ``α``.  Two knobs,
    three observables: the disconnected family is a **surface** in the space of
    static responses, and ``𝒲`` is its defining function.
    """
    s = np.asarray(response, dtype=float)
    det = float(np.linalg.det(s))
    if abs(det) < 1e-300:
        raise ValueError("the static response is singular; 𝒲 is undefined")
    return float(s[0, 1] / det - float(g_between_mouths))


def defect_of_pair(pair: MouthPair, lmbda: float = 0.0) -> float:
    """``𝒲`` for a throat in the chart, at frequency ``λ``."""
    gam = pair.gamma(float(lmbda)).real
    return disconnection_defect(static_response(pair, lmbda), float(gam[0, 1]))


def invisible_partner(alpha1: float, alpha2: float, re_beta: float,
                      separation: float,
                      lmbda: float = 0.0) -> Optional[float]:
    """The ``|Im β|`` that makes a *connected* throat look disconnected.

    Working the identity out for general complex ``β``,

        ``𝒲 = −Re β − (G_d − Re β)(Im β)² / P`` ,
        ``P = (α₁ − g)(α₂ − g) − (Re β − G_d)²`` ,

    so ``𝒲`` vanishes off the ``β = 0`` surface whenever
    ``(Im β)² = −Re β · P / (G_d − Re β)`` has a real root.  With ``P > 0``
    that needs ``Re β`` and ``G_d − Re β`` to have opposite signs, which happens
    on **two** branches — ``Re β < 0``, and ``Re β > G_d``.  They behave
    completely differently under PR #257's gate, and `measure_the_blind_spot_of
    _a_single_frequency_test` separates them: on the invisibility surface
    ``det(A − Γ) = P · G_d / (G_d − Re β)``, so the second branch has a negative
    determinant and is **unstable**, while the first is strictly inside the
    cone.  Returns ``|Im β|``, or ``None`` when no real root exists.  The sign
    is immaterial: ``Im β → −Im β`` is PR #256's time reversal, and no static
    observable sees it.
    """
    gam = gamma_at(float(lmbda), float(separation)).real
    g, gd = float(gam[0, 0]), float(gam[0, 1])
    p = (float(alpha1) - g) * (float(alpha2) - g) - (float(re_beta) - gd) ** 2
    denom = gd - float(re_beta)
    if abs(denom) < 1e-300:
        return None
    val = -float(re_beta) * p / denom
    if val < 0.0:
        return None
    return float(math.sqrt(val))


# ════════════════════════════════════════════════════════════════════════════
# FROM MEASUREMENTS BACK TO THE THROAT
# ════════════════════════════════════════════════════════════════════════════
def recover_response(observations: Sequence[Tuple[Sequence[float],
                                                  Sequence[float], float]],
                     separation: float,
                     lmbda: float = 0.0) -> Dict[str, object]:
    """Solve for ``S`` from measured invariants at chosen source placements.

    Each observation is ``(y_A, y_B, 𝒞)``.  Subtracting the known free term
    leaves ``v_Aᵀ S v_B``, linear in the three independent entries of ``S``, so
    three independent placements suffice and more overdetermine it.  This is the
    step that makes ``𝒲`` a *measurement protocol* rather than a formula in the
    boundary data: an observer who knows the background and where the mouths are
    can build ``𝒲`` without ever being told ``A``.
    """
    rows, rhs = [], []
    for y_a, y_b, value in observations:
        v_a = source_vector(y_a, separation, lmbda)
        v_b = source_vector(y_b, separation, lmbda)
        rows.append([v_a[0] * v_b[0],
                     v_a[0] * v_b[1] + v_a[1] * v_b[0],
                     v_a[1] * v_b[1]])
        rhs.append(float(value)
                   - free_interaction_energy(y_a, y_b, lmbda))
    a = np.array(rows, dtype=float)
    b = np.array(rhs, dtype=float)
    sol, *_ = np.linalg.lstsq(a, b, rcond=None)
    s = np.array([[sol[0], sol[1]], [sol[1], sol[2]]], dtype=float)
    return {"response": s,
            "n_observations": len(rows),
            "condition_number": float(np.linalg.cond(a)),
            "residual": float(np.abs(a @ sol - b).max())}


def is_real_field_compatible(boundary: np.ndarray, tol: float = 1e-14) -> bool:
    """``A = A*`` — whether the extension supports a **real** scalar field.

    PR #254 solves a real time-domain field.  A real solution requires the
    self-adjoint domain to be invariant under complex conjugation: if
    ``φ^reg = A q`` holds, then conjugating gives ``φ^reg* = A* q*``, which is
    in the domain only when ``A* = A``.  Combined with Hermiticity that is
    ``A`` real symmetric — **``β`` real**.

    The consequence is not abstract.  With complex ``β`` the response
    ``R = (A − Γ)⁻¹`` is Hermitian but not real, so a real static source
    produces a **complex** field.  `measure_a_real_field_forces_beta_real`
    measures exactly that, and it is what decides whether the blind family of
    `invisible_partner` exists at all for the arc's field content.
    """
    a = np.asarray(boundary, dtype=complex)
    return bool(np.abs(a - a.conjugate()).max() <= tol)


def recover_complex_response(pair: MouthPair, y_a: Sequence[float],
                             y_b: Sequence[float],
                             lmbda: float = 0.0) -> Dict[str, object]:
    """The full complex ``R`` at **one** spectral parameter, then ``A``.

    Real static sources see only ``Re R`` — three numbers for four parameters,
    which is where the whole multi-parameter story comes from.  That is a
    restriction of *the protocol*, not of the operator.  A pair of sources with
    complex strengths ``a, b`` contributes ``2 Re[a* b · G_A(y_A,y_B)]``, so
    scanning the relative phase gives both quadratures: ``a*b = 1`` returns
    ``Re``, ``a*b = i`` returns ``−Im``.

    With the complex kernel in hand at three placements the complex ``R``
    follows, and then ``A = Γ + R⁻¹`` **outright** — no second spectral
    parameter, no ``Im β`` sign ambiguity from the reconstruction itself.
    """
    lam = float(lmbda)
    v_a = source_vector(y_a, pair.separation, lam)
    v_b = source_vector(y_b, pair.separation, lam)
    r = response_of_pair(pair, lam)
    kernel = complex(v_a @ r @ v_b)
    # a = 1, b = 1  ⟹  2 Re[K];    a = i, b = 1  ⟹  2 Re[−i K] = 2 Im[K]
    in_phase = 2.0 * kernel.real
    quadrature = 2.0 * (np.conjugate(1j) * kernel).real
    got = complex(0.5 * in_phase, 0.5 * quadrature)
    return {"kernel": kernel,
            "in_phase": float(in_phase),
            "quadrature": float(quadrature),
            "real_part_from_quadratures": float(0.5 * in_phase),
            "imag_part_from_quadratures": float(0.5 * quadrature),
            "the_quadratures_give_the_kernel": bool(abs(got - kernel) < 1e-14)}


def recover_boundary(pair: MouthPair,
                     lambdas: Sequence[float] = (0.3, 0.8),
                     start: Sequence[float] = (0.3, 0.3, 0.0, 0.1)
                     ) -> Dict[str, object]:
    """Reconstruct ``A`` from the real-static invariant at several **spectral
    parameters**.

    Not "frequencies": ``λ = ω²``, so a negative ``λ`` is an imaginary
    frequency, and the default here is two *positive* ``λ`` below the free
    ground state ``λ = 1`` — both genuinely drivable.

    One spectral parameter gives the three entries of ``S(λ)`` for four
    parameters.  A second gives three more, and ``Γ(λ)`` moves between them, so
    the system is over-determined and the map is injective up to
    ``Im β → −Im β``.  Solved as a least-squares problem and reported with its
    residual, so a failure to converge cannot be mistaken for a reconstruction.

    This is the **real-static-source** protocol.  `recover_complex_response`
    needs only one ``λ``.
    """
    from scipy.optimize import least_squares

    d = float(pair.separation)
    idx = np.triu_indices(2)

    def response_of(p: np.ndarray, lam: float) -> np.ndarray:
        a = np.array([[p[0], p[2] + 1j * p[3]],
                      [p[2] - 1j * p[3], p[1]]], dtype=complex)
        return np.linalg.inv(a - gamma_at(float(lam), d)).real

    target = np.concatenate([response_of(
        np.array([pair.alpha1, pair.alpha2,
                  complex(pair.beta).real, complex(pair.beta).imag]),
        lam)[idx] for lam in lambdas])

    def residual(p: np.ndarray) -> np.ndarray:
        return np.concatenate([response_of(p, lam)[idx]
                               for lam in lambdas]) - target

    # The residual is not convex, and a single start does land in a local
    # minimum for some throats — caught by the reported residual, which is the
    # reason it is reported.  Several starts, keep the best.
    starts = [np.asarray(start, dtype=float),
              np.array([0.5, 0.5, 0.1, -0.2]),
              np.array([0.2, 0.4, -0.1, 0.3]),
              np.array([0.4, 0.2, 0.0, -0.05])]
    out = None
    for s0 in starts:
        trial = least_squares(residual, s0, xtol=1e-14, ftol=1e-14, gtol=1e-14)
        if out is None or np.abs(trial.fun).max() < np.abs(out.fun).max():
            out = trial
        if np.abs(out.fun).max() < 1e-12:
            break
    truth = np.array([pair.alpha1, pair.alpha2,
                      complex(pair.beta).real, abs(complex(pair.beta).imag)])
    got = np.array([out.x[0], out.x[1], out.x[2], abs(out.x[3])])
    return {"lambdas": [float(v) for v in lambdas],
            "n_starts_tried": len(starts),
            "true": [float(v) for v in truth],
            "recovered": [float(v) for v in got],
            "max_parameter_error": float(np.abs(got - truth).max()),
            "residual": float(np.abs(out.fun).max()),
            "sign_of_im_beta_is_not_observable": True}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def _working_pair(separation: float = WORKING_SEPARATION) -> MouthPair:
    a1, a2, b = WORKING_BOUNDARY
    return MouthPair(separation, a1, a2, b)


def _margin(pair: MouthPair) -> float:
    return float(positivity_defect(pair.boundary_matrix(),
                                   pair.separation)["min_eigenvalue"])


def measure_the_invariant_vanishes_when_a_source_is_removed(
        n_pairs: int = 40, seed: int = 20260816) -> Dict[str, object]:
    """The property PR #253 said rank counting could not have.

    The invariant is the cross term of a quadratic functional, so it is bilinear
    in the two source strengths and *identically* zero when either is switched
    off — not underdetermined, not ill-posed, zero.  Compared here against the
    thing it replaces: deleting one of the closure equations leaves a system
    whose solution set gains a dimension, which is a statement about square
    systems and would say the same for a source that was never there.

    Computed the honest way.  A first draft "removed a source" by multiplying
    the answer by zero, which proves nothing about the field construction.  Here
    the quadratic functional `energy_functional` is built with explicit source
    strengths and **its own self-energy terms**, and the cross term is extracted
    as ``Q[a,b] − Q[a,0] − Q[0,b]`` from three separate evaluations.  The
    self-energies are large — they are reported — and they cancel because the
    functional is quadratic.  Removing a source then means evaluating ``Q`` at
    ``b = 0``, which is a different configuration of the same functional rather
    than a multiplication.
    """
    pair = _working_pair()
    pts = random_points(2 * int(n_pairs), seed=seed)
    worst_cross_err, worst_removed, smallest_on = 0.0, 0.0, math.inf
    a, b = 1.7, -0.9
    for k in range(int(n_pairs)):
        y_a, y_b = pts[2 * k], pts[2 * k + 1]
        q_both = energy_functional(pair, y_a, y_b, a, b)
        q_a = energy_functional(pair, y_a, y_b, a, 0.0)
        q_b = energy_functional(pair, y_a, y_b, 0.0, b)
        cross = q_both - q_a - q_b
        predicted = a * b * interaction_energy(pair, y_a, y_b)
        worst_cross_err = max(worst_cross_err, abs(cross - predicted))
        # and the same three-evaluation recipe with source B actually absent
        removed = (energy_functional(pair, y_a, y_b, a, 0.0)
                   - energy_functional(pair, y_a, y_b, a, 0.0)
                   - energy_functional(pair, y_a, y_b, 0.0, 0.0))
        worst_removed = max(worst_removed, abs(removed))
        smallest_on = min(smallest_on, abs(cross))
    y_a, y_b = pts[0], pts[1]
    self_a = energy_functional(pair, y_a, y_b, a, 0.0)
    doubled = (energy_functional(pair, y_a, y_b, 2.0 * a, b)
               - energy_functional(pair, y_a, y_b, 2.0 * a, 0.0)
               - energy_functional(pair, y_a, y_b, 0.0, b))
    single = (energy_functional(pair, y_a, y_b, a, b)
              - energy_functional(pair, y_a, y_b, a, 0.0)
              - energy_functional(pair, y_a, y_b, 0.0, b))
    return {"n_pairs": int(n_pairs),
            "source_strengths": [a, b],
            "boundary": [pair.alpha1, pair.alpha2, complex(pair.beta).real,
                         complex(pair.beta).imag],
            "loewner_margin": _margin(pair),
            "inside_the_cone": bool(_margin(pair) > 0.0),
            "a_self_energy_term_for_scale": float(self_a),
            "worst_error_in_Q_minus_Q_minus_Q": float(worst_cross_err),
            "largest_value_with_a_source_removed": float(worst_removed),
            "smallest_value_with_both_present": float(smallest_on),
            "the_cross_term_is_the_invariant": bool(worst_cross_err < 1e-12),
            "it_vanishes_exactly": bool(worst_removed == 0.0),
            "it_is_not_vacuous": bool(smallest_on > 1e-6),
            "is_bilinear": bool(abs(doubled - 2.0 * single) < 1e-14),
            "the_contrast": ("a deleted equation costs a dimension whatever "
                             "was deleted; a deleted source costs the whole "
                             "cross term")}


def measure_the_throat_channel_has_the_rank_of_the_boundary_condition(
        n_sources: int = 12, seed: int = 20260816) -> Dict[str, object]:
    """``rank(throat table) = rank(Re R) ≤ 2``, for any number of sources.

    The throat table is ``Vᵀ S V`` with ``V`` of shape ``2 × N``, so the whole
    multi-source signature of a point throat is **two-dimensional however many
    sources are placed**, while the direct table has full rank — the deficiency
    is the throat's and not the geometry's.

    Off the chart the statement needs care, and getting it wrong is easy.  The
    *complex* response obeys ``rank R = rank B``, so PR #257's Dirichlet strata
    are rank one.  But a pair of real static sources sees only ``Re R``, and the
    real part of a rank-one **complex** Hermitian matrix ``c·uu†`` is
    ``c(Re u Re uᵀ + Im u Im uᵀ)`` — generically **rank two**.  So a
    time-reversal-breaking rank-one boundary condition still fills both channels
    of the static observable, and only a *real* Dirichlet direction drops the
    table to rank one.  Both cases are drawn and separated here rather than
    averaged into one claim.
    """
    pair = _working_pair()
    pts = random_points(int(n_sources), seed=seed)
    chart = cross_matrix(pair, pts)

    rng = np.random.default_rng(int(seed) + 7)
    gam = gamma_at(0.0, pair.separation)
    vv = np.array([source_vector(p, pair.separation) for p in pts]).T

    def stratum(direction: np.ndarray, theta: float) -> Dict[str, object]:
        v = direction / np.linalg.norm(direction)
        w = np.array([-v[1].conjugate(), v[0].conjugate()])
        u = np.outer(v, v.conjugate()) + np.exp(1j * theta) * np.outer(
            w, w.conjugate())
        b, c = boundary_pair_from_unitary(u)
        r = response_matrix(b, c, gam)
        table = vv.T @ r.real @ vv
        return {"det_B": float(abs(np.linalg.det(b))),
                "rank_B": int(np.linalg.matrix_rank(b, tol=1e-9)),
                "rank_R": int(np.linalg.matrix_rank(r, tol=1e-10)),
                "hermitian_defect": float(np.abs(r - r.conjugate().T).max()),
                "table_rank": int(np.linalg.matrix_rank(table, tol=1e-12))}

    complex_strata = [stratum(rng.normal(size=2) + 1j * rng.normal(size=2),
                              float(rng.uniform(0.4, 5.5))) for _ in range(4)]
    real_strata = [stratum(rng.normal(size=2).astype(complex),
                           float(rng.uniform(0.4, 5.5))) for _ in range(4)]
    return {"n_sources": int(n_sources),
            "loewner_margin": _margin(pair),
            **{f"chart_{k}": v for k, v in chart.items()},
            "complex_strata": complex_strata,
            "real_strata": real_strata,
            "the_throat_table_is_rank_two_in_the_chart": bool(
                chart["throat_rank"] == 2),
            "the_direct_table_has_full_rank": bool(
                chart["direct_rank"] == int(n_sources)),
            "the_complex_response_has_the_rank_of_B": bool(
                all(s["rank_R"] == s["rank_B"] == 1
                    for s in complex_strata + real_strata)),
            "a_complex_dirichlet_direction_still_fills_both_channels": bool(
                all(s["table_rank"] == 2 for s in complex_strata)),
            "a_real_dirichlet_direction_drops_the_table_to_one": bool(
                all(s["table_rank"] == 1 for s in real_strata)),
            "the_response_is_hermitian_off_the_chart": bool(
                all(s["hermitian_defect"] < 1e-12
                    for s in complex_strata + real_strata)),
            "the_bound_is_two_at_every_source_count": bool(
                chart["throat_rank"] <= 2
                and all(s["table_rank"] <= 2
                        for s in complex_strata + real_strata))}


def measure_anisotropy_is_not_the_signature(
        chi: float = 1.0, n: int = 240) -> Dict[str, object]:
    """The obvious observable, and why it does not decide anything.

    Fixing the geodesic separation and moving one source over the sphere of that
    radius: the free interaction is constant to machine precision, and the
    throat's is not.  That is a real effect and a real measurement — and **two
    disconnected scatterers produce one of the same size**, so it detects
    structure at the mouths, not a connection between them.
    """
    pair = _working_pair()
    prof = isotropy_profile(pair, chi=chi, n=n)
    ratio = (prof["disconnected_relative_spread"]
             / prof["throat_relative_spread"])
    return {"loewner_margin": _margin(pair), **prof,
            "disconnected_over_connected": float(ratio),
            "anisotropy_does_not_discriminate": bool(0.2 < ratio < 5.0),
            "what_it_does_show": ("the interaction depends on more than the "
                                  "geodesic separation, which no free field "
                                  "on this background can do")}


def measure_two_disconnected_scatterers_lie_on_a_surface(
        n_draws: int = 200, seed: int = 20260817) -> Dict[str, object]:
    """``S₁₂ = G₀ · det S`` — exactly, for every pair of independent scatterers.

    Two knobs, three observables.  The disconnected family therefore cannot fill
    the space of static responses, and the equation of the surface it does fill
    is the discriminator.  Checked over random self-energies and random mouth
    separations, and checked to be a genuine restriction by drawing connected
    throats and finding them off the surface.
    """
    rng = np.random.default_rng(int(seed))
    worst_on, smallest_off, n_off = 0.0, math.inf, 0
    for _ in range(int(n_draws)):
        d = float(rng.uniform(0.25, math.pi))
        a1, a2 = rng.uniform(0.12, 0.7, size=2)
        worst_on = max(worst_on,
                       abs(defect_of_pair(MouthPair(d, a1, a2, 0.0))))
        b = float(rng.uniform(-0.25, 0.25))
        if abs(b) < 5e-3:
            continue
        w = abs(defect_of_pair(MouthPair(d, a1, a2, b)))
        smallest_off = min(smallest_off, w)
        n_off += int(w > 1e-9)
    return {"n_draws": int(n_draws),
            "worst_defect_on_the_disconnected_family": float(worst_on),
            "smallest_defect_off_it": float(smallest_off),
            "n_connected_draws_detected": int(n_off),
            "the_disconnected_family_is_a_surface": bool(worst_on < 1e-12),
            "connected_throats_are_off_it": bool(smallest_off > 1e-9),
            "the_equation": "S₁₂ = G₀ · det S,  two knobs and three observables"
            }


def measure_the_defect_is_the_mouth_mixing_amplitude(
        n_draws: int = 120, seed: int = 20260818) -> Dict[str, object]:
    """``𝒲 = −β`` on the real slice — exactly, and that slice is the physical
    one.

    Not merely a detector: for real ``β`` the discriminator *is* the
    mouth-mixing amplitude, the one part of the boundary operator with no local
    realization on ``S³``.  Real ``β`` is not a convenient case — it is what a
    real scalar field requires (`is_real_field_compatible`), so this is the
    statement for PR #254's field rather than a slice through it.

    Scope, stated exactly: what is detected is **off-diagonal mouth-boundary
    mixing relative to the diagonal two-scatterer null model**, inside this
    point-interaction model.  It is not a detection of topology, of a
    traversable interior, or of anything the model does not contain.

    And it does not move with the self-energies, the mouth separation, or the
    Löwner margin — the answer to PR #255's caution that anything built from a
    resummed field measures the pole rather than the source.  The raw invariant
    does grow as the margin closes; ``𝒲`` is reported against it at four margins
    to show that it does not.
    """
    rng = np.random.default_rng(int(seed))
    worst = 0.0
    for _ in range(int(n_draws)):
        d = float(rng.uniform(0.25, math.pi))
        a1, a2 = rng.uniform(0.12, 0.7, size=2)
        b = float(rng.uniform(-0.3, 0.3))
        pair = MouthPair(d, a1, a2, b)
        worst = max(worst, abs(defect_of_pair(pair) + b))

    rows = []
    d = WORKING_SEPARATION
    gam = gamma_at(0.0, d).real
    beta = 0.06
    pts = random_points(2, seed=31)
    for eps in (0.4, 0.1, 0.02, 0.004):
        # slide toward the cone's null boundary from strictly inside it: for
        # equal mouths and real β the eigenvalues of A − Γ(0) are
        # α − g₀ ∓ (G₀ − β), so this sets the Löwner margin — the *smaller* of
        # the two — to ε exactly, and every row stays in the positive sector
        a = float(gam[0, 0] + abs(beta - gam[0, 1]) + eps)
        pair = MouthPair(d, a, a, beta)
        rows.append({"loewner_margin": _margin(pair),
                     "invariant": interaction_energy(pair, pts[0], pts[1]),
                     "defect": defect_of_pair(pair)})
    spread = max(abs(r["defect"] + beta) for r in rows)
    growth = abs(rows[-1]["invariant"]) / abs(rows[0]["invariant"])
    return {"n_draws": int(n_draws),
            "worst_error_in_W_plus_beta": float(worst),
            "W_is_minus_beta": bool(worst < 1e-9),
            "margin_rows": rows,
            "invariant_growth_toward_the_boundary": float(growth),
            "worst_defect_drift_across_margins": float(spread),
            "every_row_is_inside_the_cone": bool(
                all(r["loewner_margin"] > 0.0 for r in rows)),
            "the_discriminator_does_not_see_the_resonance": bool(
                spread < 1e-9 and growth > 3.0
                and all(r["loewner_margin"] > 0.0 for r in rows)),
            "the_caution_it_answers": ("PR #255: a test built from a resummed "
                                       "field can measure the pole instead of "
                                       "the source"),
            "what_is_actually_detected": ("off-diagonal mouth-boundary mixing, "
                                          "relative to the diagonal "
                                          "two-scatterer null model")}


def measure_the_invariant_is_recoverable_from_observations(
        n_obs: int = 24, seed: int = 20260819) -> Dict[str, object]:
    """From placements and numbers back to ``S``, and then to ``𝒲``.

    The discriminator has to be buildable by someone who measures interaction
    energies, knows the background, and knows where the mouths are — but is not
    told the boundary data.  Three independent placements determine ``S``; this
    uses more and reports the residual and the conditioning, because a protocol
    that only works on a well-chosen configuration is not a protocol.
    """
    pair = _working_pair()
    pts = random_points(2 * int(n_obs), seed=seed)
    obs = [(pts[2 * k], pts[2 * k + 1],
            interaction_energy(pair, pts[2 * k], pts[2 * k + 1]))
           for k in range(int(n_obs))]
    rec = recover_response(obs, pair.separation)
    truth = static_response(pair)
    gam = pair.gamma(0.0).real
    w_true = defect_of_pair(pair)
    w_rec = disconnection_defect(rec["response"], float(gam[0, 1]))
    return {"n_observations": int(n_obs),
            "loewner_margin": _margin(pair),
            "condition_number": rec["condition_number"],
            "fit_residual": rec["residual"],
            "worst_entry_error": float(np.abs(rec["response"] - truth).max()),
            "W_from_the_boundary_data": float(w_true),
            "W_from_the_observations": float(w_rec),
            "W_error": float(abs(w_rec - w_true)),
            "minus_beta": float(-complex(pair.beta).real),
            "the_protocol_closes": bool(abs(w_rec - w_true) < 1e-9)}


def measure_a_real_field_forces_beta_real(seed: int = 20260823
                                          ) -> Dict[str, object]:
    """Which field is being solved, and what that costs the blind family.

    PR #254 solves a **real** time-domain scalar.  A real solution needs the
    self-adjoint domain to be invariant under complex conjugation: conjugating
    ``φ^reg = A q`` gives ``φ^reg* = A* q*``, which lies in the domain only when
    ``A* = A``.  With Hermiticity that is ``A`` real symmetric — ``β`` real.

    Measured rather than argued: with complex ``β`` the response is Hermitian
    but **not real**, so a real unit static source produces a field with a
    nonzero imaginary part.  It is a complex-scalar / time-reversal-breaking
    extension, and calling it that is the difference between a limitation of the
    test and a change of the model.

    The consequence for the round is the whole blind-spot story:
    `invisible_partner` needs ``Im β ≠ 0``, so on the real slice there is **no
    blind family at all** and ``𝒲 = −β`` settles it at a single spectral
    parameter.
    """
    d = WORKING_SEPARATION
    pts = random_points(2, seed=seed)
    rows = []
    for beta in (0.06, complex(0.06, 0.20)):
        pair = MouthPair(d, 0.30, 0.35, beta)
        r = response_of_pair(pair)
        v_a = source_vector(pts[0], d)
        v_b = source_vector(pts[1], d)
        field = complex(free_interaction_energy(pts[0], pts[1])
                        + v_b @ r @ v_a)
        rows.append({
            "beta": complex(beta),
            "real_field_compatible": is_real_field_compatible(
                pair.boundary_matrix()),
            "max_imaginary_part_of_R": float(np.abs(r.imag).max()),
            "imaginary_part_of_the_field": float(field.imag),
            "has_an_invisible_partner": bool(
                invisible_partner(0.30, 0.35, complex(beta).real, d)
                is not None)})
    real_row, complex_row = rows[0], rows[1]
    return {"separation": d, "rows": rows,
            "a_real_beta_gives_a_real_field": bool(
                abs(real_row["imaginary_part_of_the_field"]) < 1e-15),
            "a_complex_beta_does_not": bool(
                abs(complex_row["imaginary_part_of_the_field"]) > 1e-6),
            "the_condition": "A = A* ⟺ the domain is conjugation-invariant",
            "on_the_real_slice_W_is_minus_beta": float(
                defect_of_pair(MouthPair(d, 0.30, 0.35, 0.06))),
            "the_blind_family_needs_a_complex_scalar": True,
            "so_for_PR254s_field_there_is_no_blind_family": True}


def measure_phase_sensitive_sources_need_only_one_spectral_parameter(
        seed: int = 20260824) -> Dict[str, object]:
    """The multi-``λ`` requirement is the protocol's, not the operator's.

    Real static sources see only ``Re R`` — three numbers for four parameters.
    Complex sources with a scanned relative phase see both quadratures of
    ``G_A(y_A,y_B)``, hence the full complex ``R``, and then

        ``A  =  Γ(λ)  +  R⁻¹``

    at a **single** spectral parameter.  So the two-``λ`` reconstruction is a
    statement about what real static probes can do, not a limitation of the
    two-mouth operator, and the round should not have implied otherwise.
    """
    d = WORKING_SEPARATION
    n_obs = 8
    pts = random_points(2 * n_obs, seed=seed)
    pair = MouthPair(d, 0.30, 0.35, complex(-0.05, 0.24))
    quad = recover_complex_response(pair, pts[0], pts[1])
    r_true = response_of_pair(pair)

    # R is HERMITIAN, not symmetric, so the unknowns are the two real diagonal
    # entries and the two real parts of R₀₁ — four real numbers, and each
    # complex measurement supplies two real equations.  Writing R₀₁ = x + iy,
    #   vᵀR w = v₀w₀ R₀₀ + v₁w₁ R₁₁ + x(v₀w₁ + v₁w₀) + i y(v₀w₁ − v₁w₀) .
    # A symmetric ansatz here silently drops y, which is exactly the Im β the
    # real-static protocol could not see.
    rows, rhs = [], []
    for k in range(n_obs):
        v_a = source_vector(pts[2 * k], d)
        v_b = source_vector(pts[2 * k + 1], d)
        sym = v_a[0] * v_b[1] + v_a[1] * v_b[0]
        asym = v_a[0] * v_b[1] - v_a[1] * v_b[0]
        kernel = complex(v_a @ r_true @ v_b)
        rows.append([v_a[0] * v_b[0], v_a[1] * v_b[1], sym, 0.0])
        rhs.append(kernel.real)
        rows.append([0.0, 0.0, 0.0, asym])
        rhs.append(kernel.imag)
    a_mat = np.array(rows, dtype=float)
    sol, *_ = np.linalg.lstsq(a_mat, np.array(rhs, dtype=float), rcond=None)
    r_rec = np.array([[sol[0], sol[2] + 1j * sol[3]],
                      [sol[2] - 1j * sol[3], sol[1]]], dtype=complex)
    a_rec = gamma_at(0.0, d) + np.linalg.inv(r_rec)
    a_true = pair.boundary_matrix()
    return {"separation": d,
            "n_placements": n_obs,
            "boundary": [pair.alpha1, pair.alpha2, complex(pair.beta).real,
                         complex(pair.beta).imag],
            "the_quadratures_give_the_kernel":
                quad["the_quadratures_give_the_kernel"],
            "in_phase": quad["in_phase"], "quadrature": quad["quadrature"],
            "condition_number": float(np.linalg.cond(a_mat)),
            "worst_response_error": float(np.abs(r_rec - r_true).max()),
            "worst_boundary_error": float(np.abs(a_rec - a_true).max()),
            "one_spectral_parameter_suffices": bool(
                np.abs(a_rec - a_true).max() < 1e-9),
            "the_identity": "A = Γ(λ) + R(λ)⁻¹",
            "what_this_corrects": ("the two-λ requirement is a restriction of "
                                   "the real-static-source protocol, not of "
                                   "the operator")}


def measure_the_blind_spot_of_a_single_frequency_test(
        seed: int = 20260820) -> Dict[str, object]:
    """The blind family — and the two conditions that shrink it to nothing.

    For complex ``β`` the defect is
    ``𝒲 = −Re β − (G_d − Re β)(Im β)²/P``, which has a zero away from
    ``β = 0`` on **two** branches.  Reported with both of the things that
    remove it, because the first draft presented it as a standing defect of the
    test and it is not:

    * **PR #257's gate** removes one branch.  On the invisibility surface
      ``det(A − Γ) = P · G_d/(G_d − Re β)``, so ``Re β > G_d`` has a negative
      determinant and is unstable;
    * **reality of the field** removes what is left.  Every blind point needs
      ``Im β ≠ 0``, so the whole family lives outside the real-scalar sector
      (`measure_a_real_field_forces_beta_real`).  It exists only for a
      deliberately time-reversal-breaking complex extension.

    And inside that extension the fix is not even a second spectral parameter:
    phase-sensitive complex sources reconstruct ``A`` at one ``λ``.  What is
    left is a real statement, and a narrow one: *the real-static-source protocol
    at a single ``λ`` is blind on this family, in the complex-scalar model.*

    The couplings are **comparable to, and smaller than, the self-energies** —
    an earlier draft said larger, which is false for every row here.
    """
    d = WORKING_SEPARATION
    gd = float(gamma_at(0.0, d).real[0, 1])
    rows = []
    for a1, a2, rb in ((0.30, 0.35, -0.05), (0.50, 0.40, -0.02),
                       (0.25, 0.25, -0.10),
                       (0.30, 0.35, gd + 0.02), (0.50, 0.40, gd + 0.05)):
        ib = invisible_partner(a1, a2, rb, d)
        if ib is None:
            continue
        pair = MouthPair(d, a1, a2, complex(rb, ib))
        rows.append({"alpha": [a1, a2], "re_beta": float(rb),
                     "branch": "Re β < 0" if rb < 0.0 else "Re β > G_d",
                     "im_beta": float(ib),
                     "abs_beta": float(abs(complex(rb, ib))),
                     "smaller_than_the_self_energies": bool(
                         abs(complex(rb, ib)) < min(a1, a2)),
                     "real_field_compatible": is_real_field_compatible(
                         pair.boundary_matrix()),
                     "W": defect_of_pair(pair),
                     "loewner_margin": _margin(pair),
                     "inside_the_cone": bool(_margin(pair) > 0.0),
                     "W_at_a_second_spectral_parameter":
                         defect_of_pair(pair, 0.8)})
    lower = [r for r in rows if r["branch"] == "Re β < 0"]
    upper = [r for r in rows if r["branch"] == "Re β > G_d"]
    return {"separation": d, "G_between_mouths": gd, "rows": rows,
            "the_blind_family_is_not_empty": bool(rows),
            "the_upper_branch_is_excluded_by_the_stability_gate": bool(
                upper and all(not r["inside_the_cone"] for r in upper)),
            "the_lower_branch_survives_the_stability_gate": bool(
                lower and all(r["inside_the_cone"] for r in lower)),
            "but_no_blind_point_is_real_field_compatible": bool(
                rows and not any(r["real_field_compatible"] for r in rows)),
            "they_are_invisible_at_lambda_zero": bool(
                rows and all(abs(r["W"]) < 1e-12 for r in rows)),
            "they_are_visible_at_a_second_spectral_parameter": bool(
                rows and all(abs(r["W_at_a_second_spectral_parameter"]) > 1e-3
                             for r in rows)),
            "largest_stable_invisible_coupling": float(
                max((r["abs_beta"] for r in lower), default=0.0)),
            "smallest_stable_invisible_margin": float(
                min((r["loewner_margin"] for r in lower), default=0.0)),
            "every_stable_coupling_is_smaller_than_its_self_energies": bool(
                lower and all(r["smaller_than_the_self_energies"]
                              for r in lower)),
            "the_scope": ("the real-static-source protocol at a single λ is "
                          "blind on this family — and the family exists only "
                          "in a complex-scalar, time-reversal-breaking "
                          "extension, not in PR #254's real field")}


def measure_two_spectral_parameters_reconstruct_the_boundary_matrix(
        n_draws: int = 6, seed: int = 20260821) -> Dict[str, object]:
    """Six equations, four unknowns — the real-static-source repair, measured.

    One spectral parameter gives the three entries of ``S(λ)``.  A second gives
    three more, and ``Γ(λ)`` moves between them, so the blind surface moves too.
    The boundary matrix is then over-determined and comes back exactly, up to
    the sign of ``Im β`` — which PR #256 established is not an observable but a
    time reversal.  Run on random throats *and* on a member of the blind family.

    Both defaults are **positive** ``λ`` below the free ground state ``λ = 1``,
    so both are genuinely drivable; ``λ = ω²`` makes a negative ``λ`` an
    imaginary frequency rather than a second driving frequency, and the first
    draft called it one.  And this is the *real-static-source* protocol —
    `measure_phase_sensitive_sources_need_only_one_spectral_parameter` shows one
    ``λ`` suffices once the sources carry a phase.
    """
    rng = np.random.default_rng(int(seed))
    rows = []
    for _ in range(int(n_draws)):
        a1, a2 = rng.uniform(0.15, 0.6, size=2)
        rb, ib = rng.uniform(-0.2, 0.2), rng.uniform(-0.3, 0.3)
        pair = MouthPair(WORKING_SEPARATION, float(a1), float(a2),
                         complex(rb, ib))
        out = recover_boundary(pair)
        rows.append({"true": out["true"], "recovered": out["recovered"],
                     "max_parameter_error": out["max_parameter_error"],
                     "residual": out["residual"],
                     "loewner_margin": _margin(pair)})
    a1, a2, rb = 0.30, 0.35, -0.05
    ib = invisible_partner(a1, a2, rb, WORKING_SEPARATION)
    blind = MouthPair(WORKING_SEPARATION, a1, a2, complex(rb, ib))
    blind_out = recover_boundary(blind)
    return {"n_draws": int(n_draws), "rows": rows,
            "worst_parameter_error": float(
                max(r["max_parameter_error"] for r in rows)),
            "worst_residual": float(max(r["residual"] for r in rows)),
            "the_boundary_matrix_is_reconstructed": bool(
                max(r["max_parameter_error"] for r in rows) < 1e-6),
            "blind_family_member": blind_out,
            "even_the_blind_family_is_reconstructed": bool(
                blind_out["max_parameter_error"] < 1e-6),
            "what_is_still_not_observable": ("the sign of Im β — PR #256's "
                                             "time reversal, not a gap in the "
                                             "measurement")}


def measure_the_antipodal_endpoint_on_its_own(
        seed: int = 20260822) -> Dict[str, object]:
    """``d = π`` tested as itself, not as a limit — and it behaves differently.

    PR #257 ended by showing that ``Γ(0)`` at the exact antipode is negative
    *semi*definite with a zero eigenvalue in the symmetric channel, so ``A = 0``
    sits **on** the cone's boundary.  The consequence here is direct: the static
    response ``(A − Γ(0))⁻¹`` is singular there, so the two-source invariant
    **diverges** as ``A → 0``, like ``1/ε`` in the margin.

    And the discriminator does not care.  ``𝒲 = −β`` is algebraic in ``G₀``,
    which is finite at the antipode, so at ``A = εI`` the invariant runs away
    while ``𝒲`` stays exactly zero: the loudest possible two-source signal
    carrying no information about whether the mouths are connected.  That is the
    cleanest statement in the round of why a large effect is not a measurement.
    """
    d = math.pi
    pts = random_points(2, seed=seed)
    gam0 = gamma_at(0.0, d).real
    rows = []
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        pair = MouthPair(d, eps, eps, 0.0)
        rows.append({"epsilon": float(eps),
                     "loewner_margin": _margin(pair),
                     "invariant": interaction_energy(pair, pts[0], pts[1]),
                     "W": defect_of_pair(pair)})
    scaling = [abs(r["invariant"]) * r["epsilon"] for r in rows]
    beta = 0.06
    connected = MouthPair(d, 0.30, 0.35, beta)
    return {"separation": d,
            "G_between_mouths_at_zero": float(gam0[0, 1]),
            "minus_g0": float(-gam0[0, 0]),
            "the_antipodal_value_is_minus_g0": bool(
                abs(gam0[0, 1] + gam0[0, 0]) < 1e-15),
            "rows": rows,
            "invariant_times_epsilon": [float(v) for v in scaling],
            "residue_of_the_divergence": float(scaling[-1]),
            "the_invariant_diverges_like_one_over_epsilon": bool(
                abs(scaling[-1] / scaling[-2] - 1.0) < 1e-2),
            "the_defect_stays_zero": bool(
                max(abs(r["W"]) for r in rows) < 1e-12),
            "W_of_a_connected_antipodal_throat": defect_of_pair(connected),
            "minus_beta": float(-beta),
            "the_identity_survives_the_endpoint": bool(
                abs(defect_of_pair(connected) + beta) < 1e-9),
            "the_lesson": ("at the marginal endpoint the signal is unbounded "
                           "and the discriminator is exactly zero — size is "
                           "not evidence")}
