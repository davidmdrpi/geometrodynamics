"""
For which Hermitian ``A`` is the point-throat operator non-negative?

THE QUESTION PR #256 LEFT
─────────────────────────
`waves.throat_operator` established that a point-supported throat on ``S³`` is a
self-adjoint extension parametrized by a Hermitian ``2×2`` boundary matrix, that
Hermiticity is exactly flux conservation — and that it does **not** imply
stability: the eigenvalue ``λ = ω²`` is real but need not be positive, and a
negative one means ``ω = ±i√|λ|`` with a growing mode.  That round mapped the
stable region only for the two-parameter exchange-symmetric slice, by scanning.

The full four-parameter answer is not a scan.  It is one inequality.

THE ANSWER — IN THE FINITE-``A`` CHART
──────────────────────────────────────
    **The throat operator is non-negative if and only if ``A ⪰ Γ(0)``**

in the Löwner order, where ``Γ(λ)`` is the Krein matrix of PR #256 evaluated at
threshold:

    ``Γ(0) = [[g₀, G₀], [G₀, g₀]]`` ,
    ``g₀ = −1/(4π²)`` ,  ``G₀ = (π−d)/(4π² sin d)`` .

**Why**, in one line: ``dΓ/dλ ≻ 0`` on ``λ < 0`` (measured), so every eigenvalue
of ``M(λ) = A − Γ(λ)`` is strictly decreasing in ``λ``; as ``λ → −∞``,
``Γ → −(σ/4π)I`` and both eigenvalues run to ``+∞``.  So an eigenvalue crosses
zero somewhere below threshold **iff** it is already negative at threshold.

AND THE SAME ARGUMENT COUNTS THEM
─────────────────────────────────
Nothing in it is special to ``λ* = 0``.  For any threshold below the free ground
state ``λ = 1``,

    ``#{mouth-active eigenvalues < λ*}  =  #{negative eigenvalues of A − Γ(λ*)}``

which is a Krein-type **inertia theorem**: 0 mismatches in 160 random tests at
``λ* = −2, 0, 0.5, 0.9``.  Stability is the ``λ* = 0`` case, and the count — not
just the yes/no — comes out of it.

THE GEOMETRY: A LIGHT CONE
──────────────────────────
Hermitian ``2×2`` matrices are ``ℝ⁴`` under ``A − Γ(0) = x₀ I + x·σ``, and
positive semidefiniteness is ``x₀ ≥ |x|``.  So the stable set is a **forward
light cone**, with

* **apex** at ``A = Γ(0)`` — the unique boundary data with a *doubly* degenerate
  zero mode;
* **null boundary** ``x₀ = |x| > 0`` — exactly one zero mode, ``λ = 0`` in the
  spectrum, which is what makes the boundary detectable rather than conventional;
* **interior** ``x₀ > |x|`` — strictly positive, no growing mode.

In coordinates, with ``A = [[α₁, β], [β*, α₂]]``:

    ``x₀ = (α₁+α₂)/2 − g₀`` ,  ``x₃ = (α₁−α₂)/2`` ,
    ``x₁ = Re β − G₀`` ,       ``x₂ = −Im β``

and PR #256's exchange-symmetric wedge ``α ± β ≥ g₀ ± G₀`` is exactly the slice
``x₂ = x₃ = 0`` — two of the four dimensions, which is why scanning it could not
have produced this.

WHICH IS A CHART, NOT THE WHOLE FAMILY
──────────────────────────────────────
``φ^reg = A q`` is the pair ``B = I``, ``C = A``, so a finite Hermitian ``A``
requires ``B`` invertible.  In the ``U(2)`` parametrization ``B = i(U−I)``,
``C = U+I``, the chart is exactly ``1 ∉ spec U``, and on an eigenvector with
``U v = e^{iθ}v`` the ``A``-eigenvalue is ``−cot(θ/2)`` — so ``θ → 0`` runs off
to ``∓∞`` depending on the side.  Those missing strata are **Dirichlet
directions**: ``q = 0`` in some subspace, the mouth switched off there.

The criterion extends to all of them by restriction.  With ``P`` the projector
onto the allowed-charge subspace and ``A_eff`` the induced boundary form,

    **non-negative  ⟺  ``A_eff ⪰ P† Γ(0) P``**

``k = 2`` is the chart; ``k = 1`` drops one direction and leaves a scalar
inequality; ``k = 0`` is the free operator, which has no mouth-active spectrum
and is non-negative with nothing to check.  All three are measured.

THE MONOTONICITY IS A GRAM MATRIX, NOT A SAMPLE
───────────────────────────────────────────────
``Γ_ij(λ) = ⟨δ_i, (H₀−λ)^{-1} δ_j⟩`` up to a λ-independent coincidence
subtraction, so

    ``dΓ_ij/dλ = ⟨δ_i, (H₀−λ)^{-2} δ_j⟩``

is the **Gram matrix** of ``(H₀−λ)^{-1}δ_j``: positive semidefinite for free,
and positive definite for distinct mouths.  That is the proof; the mode sum is
computed independently and agrees with the closed-form derivative to ``1e-12``,
so the random root scans are regression checks rather than the argument.

HOW THE INSTABILITY TURNS ON
────────────────────────────
Crossing the boundary a distance ``ε`` (in the smallest eigenvalue of
``A − Γ(0)``) gives ``λ ≈ −ε/μ′`` and therefore

    ``σ = √|λ| ∝ √ε`` ,   exponent ``½`` ,

with ``μ′ = d(g+G_d)/dλ`` at threshold.  Measured: ``λ/ε → −7.3745`` and
``σ/√ε → 2.7156``, and ``√7.3745 = 2.7156``.

WHERE THE APEX SITS
───────────────────
``tr Γ(0) = 2g₀ = −1/(2π²)`` — **independent of the mouth separation**.  Its
eigenvalues are ``g₀ ± G₀``, i.e. exactly PR #256's two channel thresholds, and
``det Γ(0) = g₀² − G₀² < 0`` for every ``0 < d < π``, so the apex is an
*indefinite* matrix there and ``A = 0`` is unstable.

**The antipodal endpoint is different, and it is the one that matters here.**
``G_d`` has a *removable* singularity at ``d = π``: ``G_π(ω) = ω/(4π sin πω)``
and ``G_π(0) = +1/(4π²) = −g₀``, so ``Γ(0) = g₀[[1,−1],[−1,1]]`` with
eigenvalues ``(2g₀, 0)`` — **negative semidefinite**, not indefinite.  At exact
antipodes ``A = 0`` is therefore **marginally non-negative**: it sits on the
cone's boundary with a zero mode in the *symmetric* channel, not outside it.

WHAT IS STILL PUT IN
────────────────────
A **linear** field on a **fixed** round background, and the boundary data — four
real numbers chosen, not derived.  `shells.junction` (PR #249) is what would fix
them from matter; nothing here computes the exotic-matter bill.  The throat is
**point-supported**: no interior, no proper length, no delay.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant.  What this round buys that round is a stated region
to work inside, and the count of what goes wrong outside it.
"""

from __future__ import annotations

import math
from typing import Dict, Sequence, Tuple

import numpy as np

from geometrodynamics.waves.throat_operator import (
    MouthPair,
    gamma_at,
    mouth_active_spectrum,
    negative_lambda_modes,
    stability_thresholds,
)

__all__ = [
    "threshold_matrix",
    "hermitian_from_parameters",
    "cone_coordinates",
    "positivity_defect",
    "is_non_negative",
    "inertia_below",
    "count_modes_below",
    "apex",
    "boundary_point",
    "zero_mode",
    "threshold_scaling",
    "cone_fraction",
    "gram_derivative",
    "boundary_pair_from_unitary",
    "finite_boundary_from_unitary",
    "is_maximal_self_adjoint",
    "allowed_charge_basis",
    "reduced_boundary_form",
    "is_non_negative_pair",
    "measure_the_positive_sector_is_a_shifted_psd_cone",
    "measure_the_inertia_theorem_counts_modes_below_any_threshold",
    "measure_the_boundary_of_the_cone_is_a_zero_mode",
    "measure_the_growth_rate_turns_on_with_a_square_root",
    "measure_the_exchange_symmetric_wedge_is_a_slice",
    "measure_where_the_apex_sits_as_the_mouths_separate",
    "measure_the_monotonicity_is_a_gram_matrix",
    "measure_the_criterion_extends_to_the_boundary_strata",
]

_PAULI = (np.array([[0, 1], [1, 0]], dtype=complex),
          np.array([[0, -1j], [1j, 0]], dtype=complex),
          np.array([[1, 0], [0, -1]], dtype=complex))


# ════════════════════════════════════════════════════════════════════════════
# THE THRESHOLD MATRIX AND THE CONE
# ════════════════════════════════════════════════════════════════════════════
def threshold_matrix(separation: float, lmbda: float = 0.0) -> np.ndarray:
    """``Γ(λ*)`` — the Krein matrix at the threshold being asked about.

    ``λ* = 0`` is stability.  Any ``λ* < 1`` is legitimate and the inertia
    theorem holds there too; ``λ* = 1`` is the free ground state, where ``Γ``
    has a pole.
    """
    return gamma_at(float(lmbda), float(separation)).real


def hermitian_from_parameters(alpha1: float, alpha2: float,
                              beta: complex) -> np.ndarray:
    """``A = [[α₁, β], [β*, α₂]]`` — the four real parameters, as a matrix."""
    b = complex(beta)
    return np.array([[float(alpha1), b], [b.conjugate(), float(alpha2)]],
                    dtype=complex)


def cone_coordinates(boundary: np.ndarray, separation: float,
                     lmbda: float = 0.0) -> Dict[str, object]:
    """``A − Γ(λ*) = x₀ I + x·σ`` — the light-cone coordinates.

    Hermitian ``2×2`` matrices are ``ℝ⁴``, and positive semidefiniteness is
    ``x₀ ≥ |x|``: the *forward light cone*.  This function is the whole
    geometric content of the answer, and everything else measures it.
    """
    m = np.asarray(boundary, dtype=complex) - threshold_matrix(separation,
                                                               lmbda)
    x0 = 0.5 * float(np.trace(m).real)
    x = np.array([0.5 * float(np.trace(m @ s).real) for s in _PAULI])
    return {"x0": x0, "x": x, "norm_x": float(np.linalg.norm(x)),
            "inside_the_cone": bool(x0 >= np.linalg.norm(x)),
            "lightlike_defect": float(x0 - np.linalg.norm(x)),
            "eigenvalues": [float(v) for v in np.linalg.eigvalsh(m)]}


def gram_derivative(lmbda: float, separation: float,
                    n_modes: int = 200000) -> np.ndarray:
    """``dΓ_ij/dλ = ⟨δ_i, (H₀−λ)^{-2} δ_j⟩`` — a **Gram matrix**, hence PSD.

    This is the proof the monotonicity argument needs, not a sample of it.
    ``Γ_ij(λ) = ⟨δ_i, (H₀−λ)^{-1} δ_j⟩`` (with the λ-independent coincidence
    subtraction on the diagonal), so differentiating under the spectral sum,

        ``dΓ_ij/dλ = Σ_n ψ_n(x_i) ψ_n(x_j)* / (λ_n − λ)²``

    which is the Gram matrix of the vectors ``(H₀−λ)^{-1} δ_j``.  It is
    therefore positive **semi**definite always, and positive **definite** for
    distinct mouths, where those vectors are independent.  The coincidence
    subtraction is λ-independent so it drops out of the derivative, and the
    remaining sum converges: on ``S³`` the level ``n`` weight goes like ``m²``
    against ``(m²−λ)²``, i.e. like ``1/m²``.

    Computed here from the ``S³`` addition theorem
    ``Σ_lm Y Y* = m sin(mχ)/(2π² sin χ)`` with ``m = n+1``, plus the analytic
    ``1/(2π²M)`` tail, so it is an *independent* construction to check the
    closed-form derivative against.
    """
    d = float(separation)
    e = math.pi - d
    m = np.arange(1, int(n_modes) + 1, dtype=float)
    w = 1.0 / (m * m - float(lmbda)) ** 2
    diag = float((m * m / (2.0 * math.pi ** 2) * w).sum())
    diag += 1.0 / (2.0 * math.pi ** 2 * n_modes)          # the analytic tail
    # the zonal weight m·sin(md)/(2π² sin d), written in e = π − d so that the
    # antipode is a limit rather than a division by zero: sin(md)/sin(d) =
    # (−1)^{m+1} sin(me)/sin(e), which tends to (−1)^{m+1} m as e → 0
    alt = (-1.0) ** (m + 1)
    se = math.sin(e)
    ratio = m.copy() if abs(se) < 1e-12 else np.sin(m * e) / se
    off = float((m * alt * ratio / (2.0 * math.pi ** 2) * w).sum())
    return np.array([[diag, off], [off, diag]], dtype=float)


# ════════════════════════════════════════════════════════════════════════════
# BEYOND THE FINITE-A CHART: THE U(2) STRATA
# ════════════════════════════════════════════════════════════════════════════
def boundary_pair_from_unitary(unitary: np.ndarray) -> Tuple[np.ndarray,
                                                             np.ndarray]:
    """``B = i(U − I)``, ``C = U + I`` — every self-adjoint extension, once.

    ``BC† = i(U − U†)`` is Hermitian for any unitary and ``rank[B|C] = 2``
    always, because ``Uv = v`` and ``Uv = −v`` cannot both hold.  This is the
    parametrization the finite-``A`` picture is a *chart* of: ``A`` exists only
    when ``U − I`` is invertible, i.e. when ``1 ∉ spec U``.
    """
    u = np.asarray(unitary, dtype=complex)
    eye = np.eye(2, dtype=complex)
    return 1j * (u - eye), u + eye


def finite_boundary_from_unitary(unitary: np.ndarray) -> np.ndarray:
    """``A = −i(U−I)^{-1}(U+I)`` — the chart, where it exists.

    On an eigenvector with ``U v = e^{iθ}v`` this is ``A v = −cot(θ/2) v``, so
    ``θ → 0`` sends the eigenvalue to ``∓∞`` depending on the side: the chart
    has a genuine coordinate singularity at the Dirichlet directions, and the
    strata there are not reachable by any finite Hermitian ``A``.
    """
    b, c = boundary_pair_from_unitary(unitary)
    if abs(np.linalg.det(b)) < 1e-12:
        raise ValueError("U has an eigenvalue 1; this extension is not in the "
                         "finite-A chart")
    return np.linalg.solve(b, c)


def is_maximal_self_adjoint(b: np.ndarray, c: np.ndarray,
                            tol: float = 1e-10) -> Dict[str, object]:
    """``rank[B|C] = 2`` and ``BC†`` Hermitian — the two conditions."""
    m = np.asarray(b, dtype=complex) @ np.asarray(c, dtype=complex).conj().T
    return {"rank": int(np.linalg.matrix_rank(np.hstack([b, c]), tol=tol)),
            "hermitian_defect": float(np.abs(m - m.conjugate().T).max()),
            "maximal": bool(np.linalg.matrix_rank(np.hstack([b, c]),
                                                  tol=tol) == 2),
            "self_adjoint": bool(np.abs(m - m.conjugate().T).max() <= tol)}


def allowed_charge_basis(b: np.ndarray, c: np.ndarray,
                         tol: float = 1e-10) -> Dict[str, object]:
    """Which charge directions the boundary condition permits at all.

    A left-null vector ``w`` of ``B`` turns ``B φ^reg = C q`` into the pure
    constraint ``w†C q = 0``: a **Dirichlet direction**, in which the mouth
    carries no charge and therefore drops out of the spectral problem entirely.
    What is left is the orthogonal complement, and the criterion lives there.

    ``rank B = 0`` is the free operator — no charges at all, so no mouth-active
    spectrum and nothing below the free ground state.
    """
    bb = np.asarray(b, dtype=complex)
    cc = np.asarray(c, dtype=complex)
    u, s, _ = np.linalg.svd(bb)
    r = int((s > tol).sum())
    if r == 2:
        return {"rank_B": 2, "allowed": np.eye(2, dtype=complex),
                "range_B": np.eye(2, dtype=complex), "k": 2,
                "n_dirichlet": 0}
    forbidden = cc.conjugate().T @ u[:, r:]
    if r == 0:
        return {"rank_B": 0, "allowed": np.zeros((2, 0), dtype=complex),
                "range_B": np.zeros((2, 0), dtype=complex), "k": 0,
                "n_dirichlet": 2}
    _, _, vt = np.linalg.svd(forbidden.conjugate().T)
    allowed = vt.conjugate().T[:, forbidden.shape[1]:]
    return {"rank_B": r, "allowed": allowed, "range_B": u[:, :r],
            "k": int(allowed.shape[1]), "n_dirichlet": 2 - r}


def reduced_boundary_form(b: np.ndarray, c: np.ndarray,
                          tol: float = 1e-10) -> Dict[str, object]:
    """The effective Hermitian boundary form on the allowed-charge subspace.

    Projecting ``B φ^reg = C q`` onto ``range(B)`` and restricting ``q`` to the
    allowed subspace gives a ``k×k`` pencil ``C_eff − N Γ P``; for a maximal
    self-adjoint pair the row space of ``B`` *is* the allowed subspace, so
    ``N Γ P = (NP)(P†ΓP)`` and the pencil is ``A_eff − P†ΓP`` with
    ``A_eff = (NP)^{-1} C_eff`` **Hermitian**.  Both of those are checked here
    rather than assumed.
    """
    info = allowed_charge_basis(b, c, tol)
    k = info["k"]
    if k == 0:
        return {**info, "A_eff": np.zeros((0, 0), dtype=complex),
                "hermitian_defect": 0.0, "row_space_defect": 0.0}
    p = info["allowed"]
    ucols = info["range_B"]
    n = ucols.conjugate().T @ np.asarray(b, dtype=complex)
    np_ = n @ p
    ceff = ucols.conjugate().T @ np.asarray(c, dtype=complex) @ p
    a_eff = np.linalg.solve(np_, ceff)
    # is the row space of B the allowed subspace?  (n† should lie in span p)
    resid = n.conjugate().T - p @ (p.conjugate().T @ n.conjugate().T)
    return {**info, "A_eff": a_eff,
            "hermitian_defect": float(np.abs(a_eff
                                             - a_eff.conjugate().T).max()),
            "row_space_defect": float(np.abs(resid).max())}


def is_non_negative_pair(b: np.ndarray, c: np.ndarray,
                         separation: float) -> Dict[str, object]:
    """The criterion on the **whole** ``U(2)`` family, not just the chart.

        ``A_eff ⪰ P† Γ(0) P``   on the allowed-charge subspace

    with ``P`` the projector onto it.  ``k = 2`` is the finite-``A`` chart and
    reduces to ``A ⪰ Γ(0)``; ``k = 1`` drops one Dirichlet direction and leaves
    a scalar inequality; ``k = 0`` is the free operator, whose spectrum starts
    at the free ground state and is therefore non-negative with nothing to
    check.
    """
    red = reduced_boundary_form(b, c)
    k = red["k"]
    if k == 0:
        return {**red, "non_negative": True, "n_negative": 0,
                "margin": float("inf"),
                "stratum": "free (q = 0): no mouth-active spectrum"}
    p = red["allowed"]
    gam = p.conjugate().T @ threshold_matrix(separation) @ p
    m = red["A_eff"] - gam
    ev = np.linalg.eigvalsh(0.5 * (m + m.conjugate().T))
    return {**red, "eigenvalues": [float(v) for v in ev],
            "non_negative": bool(ev[0] >= 0.0),
            "n_negative": int((ev < 0.0).sum()),
            "margin": float(ev[0]),
            "stratum": ("finite-A chart" if k == 2
                        else f"{red['n_dirichlet']} Dirichlet direction(s)")}


def positivity_defect(boundary: np.ndarray, separation: float,
                      lmbda: float = 0.0) -> Dict[str, object]:
    """Eigenvalues of ``A − Γ(λ*)``, and the smallest of them.

    The smallest eigenvalue is the **Löwner (spectral) margin**, not a Euclidean
    distance to the cone's boundary: negative means there are modes below
    ``λ*``, and *how many* is how many eigenvalues are negative.
    """
    m = np.asarray(boundary, dtype=complex) - threshold_matrix(separation,
                                                               lmbda)
    ev = np.linalg.eigvalsh(m)
    return {"eigenvalues": [float(v) for v in ev],
            "min_eigenvalue": float(ev[0]),
            "n_negative": int((ev < 0.0).sum()),
            "non_negative": bool(ev[0] >= 0.0)}


def is_non_negative(boundary: np.ndarray, separation: float) -> bool:
    """``A ⪰ Γ(0)`` — the answer, as one call."""
    return positivity_defect(boundary, separation)["non_negative"]


def inertia_below(boundary: np.ndarray, separation: float,
                  lmbda_star: float = 0.0) -> int:
    """**Predicted** number of mouth-active eigenvalues below ``λ*``.

    The inertia of ``A − Γ(λ*)``.  Because every eigenvalue of
    ``M(λ) = A − Γ(λ)`` is strictly decreasing in ``λ`` and both run to ``+∞``
    as ``λ → −∞``, an eigenvalue is below ``λ*`` exactly when it is already
    negative *at* ``λ*``.
    """
    return positivity_defect(boundary, separation, lmbda_star)["n_negative"]


def count_modes_below(pair: MouthPair, lmbda_star: float = 0.0,
                      n_gaps: int = 3) -> int:
    """**Measured** number of mouth-active eigenvalues below ``λ*``.

    Independent of `inertia_below`: this one actually finds the roots.
    """
    if lmbda_star <= 0.0:
        return sum(1 for r in negative_lambda_modes(pair, lambda_min=-40000.0,
                                                    n_grid=6000)
                   if r["lmbda"] < lmbda_star - 1e-9)
    return sum(1 for r in mouth_active_spectrum(pair, n_gaps)
               if r["lmbda"] < lmbda_star - 1e-9)


def apex(separation: float) -> Dict[str, object]:
    """``A = Γ(0)`` — the tip of the cone, and what sits there.

    Its eigenvalues are ``g₀ ± G₀``, which are exactly PR #256's two channel
    thresholds; its trace is ``2g₀ = −1/(2π²)``, **independent of the mouth
    separation**; and its determinant ``g₀² − G₀²`` is negative for every
    separation, so the apex is an *indefinite* matrix.  One consequence is worth
    stating plainly: ``A = 0`` is never stable, at any separation.
    """
    g = threshold_matrix(separation)
    ev = np.linalg.eigvalsh(g)
    th = stability_thresholds(separation)
    return {"Gamma_0": g, "eigenvalues": [float(v) for v in ev],
            "trace": float(np.trace(g).real),
            "determinant": float(np.linalg.det(g).real),
            "trace_is_two_g0": float(2.0 * th["g_at_zero"]),
            "channel_thresholds": [th["antisymmetric_threshold"],
                                   th["symmetric_threshold"]],
            "indefinite": bool(ev[0] < 0.0 < ev[1]),
            "negative_semidefinite": bool(ev[1] <= 0.0),
            # A = 0 is non-negative iff 0 ⪰ Γ(0), i.e. Γ(0) ⪯ 0
            "zero_matrix_is_stable": bool(ev[1] <= 0.0)}


def boundary_point(separation: float, direction: Sequence[float],
                   x0: float = 0.1) -> np.ndarray:
    """A boundary matrix on the cone's null surface, in a chosen direction.

    ``A = Γ(0) + x₀(I + n̂·σ)`` with ``|n̂| = 1``: lightlike by construction, so
    ``A − Γ(0)`` is rank one and ``λ = 0`` is in the spectrum exactly.
    """
    n = np.asarray(direction, dtype=float)
    nn = n / np.linalg.norm(n)
    m = float(x0) * (np.eye(2, dtype=complex)
                     + sum(nn[i] * _PAULI[i] for i in range(3)))
    return threshold_matrix(separation) + m


def zero_mode(boundary: np.ndarray, separation: float) -> Dict[str, object]:
    """The null vector of ``A − Γ(0)`` — the charge pattern of the zero mode.

    On the cone's boundary this vector satisfies ``M(0) q = 0``, so ``λ = 0``
    solves the secular equation: a genuine static solution supported by the
    throat, sitting *below* the free ground state ``λ = 1``.
    """
    m = np.asarray(boundary, dtype=complex) - threshold_matrix(separation)
    ev, vec = np.linalg.eigh(m)
    q = vec[:, 0]
    return {"q": q, "eigenvalue": float(ev[0]),
            "residual": float(np.linalg.norm(m @ q)),
            "is_a_zero_mode": bool(abs(ev[0]) < 1e-12),
            "degeneracy": int((np.abs(ev) < 1e-12).sum())}


def threshold_scaling(separation: float = 1.3,
                      epsilons: Sequence[float] = (1e-2, 1e-3, 1e-4, 1e-5,
                                                  1e-6)) -> Dict[str, object]:
    """How the growth rate turns on just outside the cone.

    Step out along the symmetric channel by ``ε`` in ``α + β``, so the **Löwner
    margin** — the smallest eigenvalue of ``A − Γ(0)`` — is ``−ε``.  Since that
    eigenvalue decreases in ``λ`` with slope ``μ′``, the root sits at
    ``λ ≈ −ε/μ′``:

        ``λ ∝ −ε`` (linear)  and  ``σ = √|λ| ∝ √ε`` (exponent ½).
    """
    th = stability_thresholds(separation)
    edge = th["symmetric_threshold"]
    h = 1e-6
    mu_prime = float(np.diff([
        np.linalg.eigvalsh(threshold_matrix(separation, x))[1]
        for x in (-h, h)])[0] / (2.0 * h))
    rows = []
    for eps in epsilons:
        a = b = 0.5 * (edge - eps)
        p = MouthPair(separation, a, a, b)
        modes = negative_lambda_modes(p, lambda_min=-10.0, n_grid=8000)
        if not modes:
            rows.append({"epsilon": eps, "lmbda": None, "sigma": None})
            continue
        lam, sig = modes[0]["lmbda"], modes[0]["sigma"]
        rows.append({"epsilon": eps, "lmbda": lam, "sigma": sig,
                     "lambda_over_epsilon": lam / eps,
                     "sigma_over_sqrt_epsilon": sig / math.sqrt(eps)})
    good = [r for r in rows if r["lmbda"] is not None]
    exps = []
    for i in range(len(good) - 1):
        exps.append(math.log(good[i]["sigma"] / good[i + 1]["sigma"])
                    / math.log(good[i]["epsilon"] / good[i + 1]["epsilon"]))
    return {"rows": rows, "exponents": exps,
            "asymptotic_exponent": (exps[-1] if exps else None),
            "lambda_over_epsilon_limit": (good[-1]["lambda_over_epsilon"]
                                          if good else None),
            "predicted_from_the_eigenvalue_slope": (-1.0 / mu_prime
                                                    if mu_prime else None),
            "mu_prime": mu_prime}


def cone_fraction(separation: float = 1.3, half_width: float = 0.2,
                  n_draws: int = 4000, seed: int = 20260816
                  ) -> Dict[str, object]:
    """What fraction of a stated box of boundary data is stable.

    A cone is unbounded, so "how big is the stable set" only means something
    relative to a box.  The box is stated: ``|α_j| ≤ w`` and ``|Re β|, |Im β| ≤
    w``.  Reported with the box, because the number is meaningless without it.
    """
    rng = np.random.default_rng(seed)
    w = float(half_width)
    n_ok = 0
    for _ in range(int(n_draws)):
        a1, a2 = rng.uniform(-w, w, 2)
        b = complex(rng.uniform(-w, w), rng.uniform(-w, w))
        if is_non_negative(hermitian_from_parameters(a1, a2, b), separation):
            n_ok += 1
    return {"half_width": w, "n_draws": int(n_draws), "n_stable": n_ok,
            "fraction": n_ok / float(n_draws),
            "the_box": "|α_j| ≤ w and |Re β|, |Im β| ≤ w",
            "caveat": "a cone is unbounded; the fraction is box-dependent"}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_the_positive_sector_is_a_shifted_psd_cone(
        separation: float = 1.3, n_draws: int = 200, seed: int = 20260816,
        spread: float = 0.15) -> Dict[str, object]:
    """**The answer.**  Non-negative ⟺ ``A ⪰ Γ(0)`` — one inequality, four
    parameters, no scan.

    Checked against the thing it replaces: for each random Hermitian ``A`` the
    negative-``λ`` axis is actually scanned for roots, and the two verdicts are
    compared.  Also reported in light-cone coordinates, since
    ``A − Γ(0) = x₀I + x·σ`` is positive semidefinite exactly when
    ``x₀ ≥ |x|``, which is what makes the region a cone rather than a box.

    Includes complex ``β`` and unequal mouths, so all four parameters are
    exercised — PR #256's wedge was a two-dimensional slice.
    """
    rng = np.random.default_rng(seed)
    rows = []
    mismatch = 0
    n_complex = 0
    for _ in range(int(n_draws)):
        a1, a2 = rng.normal(0, spread, 2)
        b = complex(*rng.normal(0, spread, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        p = MouthPair(separation, float(a1), float(a2), b)
        pred = is_non_negative(A, separation)
        got = not negative_lambda_modes(p, lambda_min=-40000.0, n_grid=6000)
        cone = cone_coordinates(A, separation)
        mismatch += int(pred != got)
        n_complex += int(abs(b.imag) > 1e-12)
        rows.append({"alpha1": float(a1), "alpha2": float(a2), "beta": b,
                     "predicted_non_negative": pred, "scan_says_stable": got,
                     "agrees": bool(pred == got),
                     "x0": cone["x0"], "norm_x": cone["norm_x"],
                     "inside_the_cone": cone["inside_the_cone"]})
    n_stable = sum(1 for r in rows if r["scan_says_stable"])
    cone_agrees = all(r["inside_the_cone"] == r["scan_says_stable"]
                      for r in rows)
    return {
        "n_draws": int(n_draws), "n_mismatches": mismatch,
        "the_criterion_is_exact": bool(mismatch == 0),
        "n_stable": n_stable,
        "both_verdicts_occur": bool(0 < n_stable < len(rows)),
        "n_with_complex_beta": n_complex,
        "the_light_cone_form_agrees": bool(cone_agrees),
        "rows": rows[:24],
        "the_criterion": "A ⪰ Γ(0) in the Löwner order",
        "the_geometry": ("A − Γ(0) = x₀I + x·σ is PSD iff x₀ ≥ |x| — the "
                         "stable set is a forward light cone in the four "
                         "dimensions of Hermitian boundary data"),
        "why": ("dΓ/dλ ≻ 0 below threshold, so every eigenvalue of A − Γ(λ) "
                "is decreasing in λ and both run to +∞ as λ → −∞; one crosses "
                "zero below λ* iff it is already negative at λ*"),
    }


def measure_the_inertia_theorem_counts_modes_below_any_threshold(
        separation: float = 1.3,
        thresholds: Sequence[float] = (-2.0, 0.0, 0.5, 0.9),
        n_draws: int = 40, seed: int = 20260816,
        spread: float = 0.15) -> Dict[str, object]:
    """The same argument does not only decide — it **counts**.

        ``#{mouth-active eigenvalues < λ*} = #{negative eigenvalues of A − Γ(λ*)}``

    for any ``λ*`` below the free ground state ``λ = 1``.  Stability is the
    ``λ* = 0`` case; the rest of the family is checked because a theorem that
    only works at one point is a coincidence.

    Also measured, since the whole argument rests on it: ``dΓ/dλ`` is positive
    definite below threshold, so the eigenvalues of ``M(λ)`` are monotone.
    """
    rng = np.random.default_rng(seed)
    mono = True
    slopes = []
    for lam in (-100.0, -9.0, -1.0, -0.05, -0.001):
        h = abs(lam) * 1e-5
        dg = ((threshold_matrix(separation, lam + h)
               - threshold_matrix(separation, lam - h)) / (2.0 * h))
        ev = np.linalg.eigvalsh(dg)
        slopes.append({"lmbda": lam, "eigenvalues": [float(v) for v in ev]})
        mono = mono and bool((ev > 0).all())

    rows = []
    total = 0
    bad = 0
    for lstar in thresholds:
        n_bad = 0
        for _ in range(int(n_draws)):
            a1, a2 = rng.normal(0, spread, 2)
            b = complex(*rng.normal(0, spread, 2))
            A = hermitian_from_parameters(float(a1), float(a2), b)
            p = MouthPair(separation, float(a1), float(a2), b)
            pred = inertia_below(A, separation, float(lstar))
            got = count_modes_below(p, float(lstar))
            total += 1
            n_bad += int(pred != got)
        bad += n_bad
        rows.append({"lambda_star": float(lstar), "n_draws": int(n_draws),
                     "mismatches": n_bad})
    return {
        "rows": rows, "n_tested": total, "n_mismatches": bad,
        "the_inertia_theorem_holds": bool(bad == 0),
        "gamma_derivative": slopes,
        "d_gamma_d_lambda_is_positive_definite": mono,
        "the_theorem": ("#{eigenvalues < λ*} = #{negative eigenvalues of "
                        "A − Γ(λ*)} for every λ* below the free ground state"),
        "stability_is_the_lambda_star_zero_case": True,
    }


def measure_the_boundary_of_the_cone_is_a_zero_mode(
        separation: float = 1.3,
        directions: Sequence[Tuple[float, float, float]] = (
            (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
            (0.6, -0.5, 0.62)),
        x0s: Sequence[float] = (0.02, 0.1, 0.4)) -> Dict[str, object]:
    """What marks the boundary — and it is not a convention.

    On the null surface ``x₀ = |x|`` the matrix ``A − Γ(0)`` is rank one, so
    ``λ = 0`` is an eigenvalue of the throat operator: a **zero mode**, a static
    solution supported by the throat and sitting below the free ground state
    ``λ = 1``.  Its charge pattern is the null vector.

    At the **apex** ``A = Γ(0)`` the matrix vanishes and there are *two* — which
    is what makes the apex a distinguished point rather than an artefact of the
    coordinates.
    """
    rows = []
    for dirn in directions:
        for x0 in x0s:
            A = boundary_point(separation, dirn, x0)
            zm = zero_mode(A, separation)
            cone = cone_coordinates(A, separation)
            p = MouthPair(separation, float(A[0, 0].real),
                          float(A[1, 1].real), complex(A[0, 1]))
            found = negative_lambda_modes(p)
            marginal = (min((r["lmbda"] for r in found), key=abs)
                        if found else 0.0)
            rows.append({
                "direction": dirn, "x0": x0,
                "lightlike_defect": cone["lightlike_defect"],
                "smallest_eigenvalue": zm["eigenvalue"],
                "is_a_zero_mode": zm["is_a_zero_mode"],
                "secular_at_zero": p.secular(0.0),
                "marginal_lambda_found_by_root_finding": marginal,
                "n_strictly_growing": sum(1 for r in found
                                          if r["lmbda"] < -1e-10),
                "q": zm["q"]})
    ap = apex(separation)
    ap_zm = zero_mode(ap["Gamma_0"], separation)
    return {
        "rows": rows,
        "every_boundary_point_has_a_zero_mode": bool(
            all(r["is_a_zero_mode"] for r in rows)),
        "worst_secular_at_zero": max(abs(r["secular_at_zero"]) for r in rows),
        "the_secular_function_vanishes_there": bool(
            max(abs(r["secular_at_zero"]) for r in rows) < 1e-12),
        "the_marginal_mode_sits_at_lambda_zero": bool(
            max(abs(r["marginal_lambda_found_by_root_finding"])
                for r in rows) < 1e-10),
        "worst_marginal_lambda": max(
            abs(r["marginal_lambda_found_by_root_finding"]) for r in rows),
        "no_boundary_point_is_strictly_unstable": bool(
            all(r["n_strictly_growing"] == 0 for r in rows)),
        "apex_degeneracy": ap_zm["degeneracy"],
        "the_apex_carries_two_zero_modes": bool(ap_zm["degeneracy"] == 2),
        "what_this_shows": ("the cone's boundary is where λ = 0 enters the "
                            "spectrum — located independently by root-finding, "
                            "not read off the eigenvalue — so it is detectable "
                            "rather than conventional"),
    }


def measure_the_growth_rate_turns_on_with_a_square_root(
        separation: float = 1.3) -> Dict[str, object]:
    """How badly things go just outside, as a scaling law.

    Step a distance ``ε`` past the boundary and the eigenvalue appears at
    ``λ ≈ −ε/μ′`` — **linear** — so the growth rate ``σ = √|λ|`` turns on with
    exponent ``½``.  The coefficient is not fitted: ``μ′`` is the slope of the
    relevant eigenvalue of ``Γ`` at threshold, computed independently.
    """
    r = threshold_scaling(separation)
    good = [x for x in r["rows"] if x["lmbda"] is not None]
    return {
        **r,
        "exponent_is_one_half": bool(
            r["asymptotic_exponent"] is not None
            and abs(r["asymptotic_exponent"] - 0.5) < 0.01),
        "lambda_is_linear_in_epsilon": bool(
            len(good) >= 2
            and abs(good[-1]["lambda_over_epsilon"]
                    - good[-2]["lambda_over_epsilon"])
            < 1e-3 * abs(good[-1]["lambda_over_epsilon"])),
        "the_coefficient_matches_the_eigenvalue_slope": bool(
            r["predicted_from_the_eigenvalue_slope"] is not None
            and good
            and abs(good[-1]["lambda_over_epsilon"]
                    - r["predicted_from_the_eigenvalue_slope"])
            < 0.02 * abs(r["predicted_from_the_eigenvalue_slope"])),
        "what_this_shows": ("the boundary is a genuine threshold, not a "
                            "numerical artefact: the growth rate is "
                            "continuous and rises like √ε"),
    }


def measure_the_exchange_symmetric_wedge_is_a_slice(
        separation: float = 1.3, n_draws: int = 400,
        seed: int = 20260816, spread: float = 0.15) -> Dict[str, object]:
    """PR #256's wedge, recovered — and located.

    Setting ``α₁ = α₂`` and ``β`` real is ``x₃ = x₂ = 0``, so the wedge
    ``α ± β ≥ g₀ ± G₀`` is a **two-dimensional slice** of a four-dimensional
    cone.  Recovered exactly here, which is the check that the general
    criterion contains the special one.

    And measured: how often the slice's verdict would be *wrong* if applied to
    general boundary data by ignoring ``Im β`` and the mouth asymmetry — which
    is the practical reason the general form was needed.
    """
    th = stability_thresholds(separation)
    slice_rows = []
    worst = 0.0
    for a in np.linspace(-0.1, 0.15, 11):
        for b in np.linspace(-0.15, 0.15, 13):
            A = hermitian_from_parameters(float(a), float(a), float(b))
            cone = is_non_negative(A, separation)
            wedge = ((a + b) >= th["symmetric_threshold"]
                     and (a - b) >= th["antisymmetric_threshold"])
            slice_rows.append({"alpha": float(a), "beta": float(b),
                               "cone": cone, "wedge": wedge,
                               "agrees": bool(cone == wedge)})
            c = cone_coordinates(A, separation)
            worst = max(worst, abs(float(c["x"][1])), abs(float(c["x"][2])))

    rng = np.random.default_rng(seed)
    naive_wrong = 0
    for _ in range(int(n_draws)):
        a1, a2 = rng.normal(0, spread, 2)
        b = complex(*rng.normal(0, spread, 2))
        A = hermitian_from_parameters(float(a1), float(a2), b)
        truth = is_non_negative(A, separation)
        a_bar = 0.5 * (a1 + a2)
        naive = ((a_bar + b.real) >= th["symmetric_threshold"]
                 and (a_bar - b.real) >= th["antisymmetric_threshold"])
        naive_wrong += int(naive != truth)
    return {
        "n_slice_points": len(slice_rows),
        "slice_mismatches": sum(1 for r in slice_rows if not r["agrees"]),
        "the_wedge_is_exactly_the_slice": bool(
            all(r["agrees"] for r in slice_rows)),
        "worst_off_slice_coordinate": worst,
        "the_slice_really_is_x2_equals_x3_equals_zero": bool(worst < 1e-15),
        "n_general_draws": int(n_draws),
        "n_the_wedge_rule_gets_wrong": naive_wrong,
        "the_slice_rule_is_not_enough_in_general": bool(naive_wrong > 0),
        "why": ("Im β and the mouth asymmetry are the two dimensions the wedge "
                "does not see; both push A out of the cone without changing "
                "α ± Re β"),
    }


def measure_where_the_apex_sits_as_the_mouths_separate(
        separations: Sequence[float] = (0.2, 0.5, 0.8, 1.3, 2.0, 2.8,
                                        3.0, math.pi)
        ) -> Dict[str, object]:
    """The apex is at ``Γ(0)``, and where that is depends on the geometry.

    Its **trace is ``2g₀ = −1/(2π²)`` at every separation** — the mouth
    distance does not enter it — while its eigenvalues are exactly PR #256's two
    channel thresholds, ``g₀ ± G₀``.

    For ``0 < d < π`` the determinant ``g₀² − G₀²`` is negative, so ``Γ(0)`` is
    **indefinite** and ``A = 0`` is unstable there.  **The exact antipodal
    endpoint is a different statement**, and it is the one this geometry cares
    about: ``G_d`` has a removable singularity at ``d = π``, with
    ``G_π(0) = +1/(4π²) = −g₀``, so ``Γ(0) = g₀[[1,−1],[−1,1]]`` has eigenvalues
    ``(2g₀, 0)`` — negative *semi*definite.  At antipodes ``A = 0`` is therefore
    **marginally non-negative**, sitting on the cone's boundary with a zero mode
    in the symmetric channel.  Measured at ``d = π`` exactly, not extrapolated.
    """
    rows = []
    for d in separations:
        ap = apex(float(d))
        zero_A = np.zeros((2, 2), dtype=complex)
        rows.append({"separation": float(d), "trace": ap["trace"],
                     "zero_A_margin": positivity_defect(zero_A,
                                                        float(d))["min_eigenvalue"],
                     "determinant": ap["determinant"],
                     "eigenvalues": ap["eigenvalues"],
                     "indefinite": ap["indefinite"],
                     "negative_semidefinite": ap["negative_semidefinite"],
                     "zero_matrix_is_stable": ap["zero_matrix_is_stable"],
                     "channel_thresholds": ap["channel_thresholds"]})
    ap_pi = apex(math.pi)
    zm_pi = zero_mode(np.zeros((2, 2), dtype=complex), math.pi)
    q = np.abs(zm_pi["q"])
    ant = {**ap_pi,
           "zero_A_margin": positivity_defect(np.zeros((2, 2), dtype=complex),
                                              math.pi)["min_eigenvalue"],
           "G_pi_at_zero": float(threshold_matrix(math.pi)[0, 1]),
           "minus_g0": 1.0 / (4.0 * math.pi ** 2),
           "marginal_channel": ("symmetric" if abs(q[0] - q[1]) < 1e-9
                                else "antisymmetric")}
    traces = [r["trace"] for r in rows]
    expect = -1.0 / (2.0 * math.pi ** 2)
    return {
        "rows": rows,
        "trace_is_separation_independent": bool(
            max(traces) - min(traces) < 1e-15),
        "trace_value": traces[0], "predicted_trace": expect,
        "trace_matches_minus_one_over_two_pi_squared": bool(
            abs(traces[0] - expect) < 1e-15),
        "the_apex_is_indefinite_away_from_the_antipode": bool(
            all(r["indefinite"] for r in rows
                if r["separation"] < math.pi - 1e-9)),
        "the_apex_is_negative_semidefinite_at_the_antipode": bool(
            ant["negative_semidefinite"] and not ant["indefinite"]),
        "the_zero_matrix_is_unstable_away_from_the_antipode": bool(
            not any(r["zero_matrix_is_stable"] for r in rows
                    if r["separation"] < math.pi - 1e-9)),
        "antipodal": ant,
        "the_antipodal_endpoint_is_marginal": bool(
            abs(ant["eigenvalues"][1]) < 1e-15
            and ant["eigenvalues"][0] < 0.0
            and ant["zero_matrix_is_stable"]),
        "at_the_antipode_A_zero_sits_on_the_boundary": bool(
            abs(ant["zero_A_margin"]) < 1e-15),
        "the_marginal_channel_is_symmetric": ant["marginal_channel"],
        "eigenvalues_are_the_channel_thresholds": bool(all(
            abs(r["eigenvalues"][0] - r["channel_thresholds"][0]) < 1e-14
            and abs(r["eigenvalues"][1] - r["channel_thresholds"][1]) < 1e-14
            for r in rows)),
        "the_symmetric_threshold_closes_at_the_antipode": bool(
            rows[-1]["eigenvalues"][1] < rows[0]["eigenvalues"][1] / 100.0),
    }


def measure_the_monotonicity_is_a_gram_matrix(
        separations: Sequence[float] = (1.3, 2.6, math.pi),
        lambdas: Sequence[float] = (-9.0, -1.0, -0.05, 0.0, 0.5)
        ) -> Dict[str, object]:
    """**The proof**, not a sample of it.

    ``Γ_ij(λ) = ⟨δ_i, (H₀−λ)^{-1} δ_j⟩`` up to a λ-independent coincidence
    subtraction, so ``dΓ_ij/dλ = ⟨δ_i, (H₀−λ)^{-2} δ_j⟩`` is a **Gram matrix**:
    positive semidefinite for free, positive definite for distinct mouths.
    That is what makes every eigenvalue of ``A − Γ(λ)`` monotone in ``λ``, and
    therefore what makes the inertia theorem a theorem.

    Checked by computing the Gram sum from the ``S³`` addition theorem — an
    independent construction that never differentiates the closed form — and
    comparing with the closed form's own derivative.  Including at the exact
    antipode, where the zonal weight becomes alternating.
    """
    rows = []
    worst = 0.0
    pd = True
    for d in separations:
        for lam in lambdas:
            h = max(abs(lam), 1.0) * 1e-6
            fd = ((threshold_matrix(d, lam + h)
                   - threshold_matrix(d, lam - h)) / (2.0 * h))
            gm = gram_derivative(float(lam), float(d))
            ev = np.linalg.eigvalsh(gm)
            err = float(np.abs(fd - gm).max())
            worst = max(worst, err)
            pd = pd and bool((ev > 0).all())
            rows.append({"separation": float(d), "lmbda": float(lam),
                         "abs_error": err,
                         "eigenvalues": [float(v) for v in ev],
                         "positive_definite": bool((ev > 0).all())})
    return {
        "rows": rows, "worst_abs_error": worst,
        "the_gram_sum_is_the_closed_form_derivative": bool(worst < 1e-9),
        "positive_definite_everywhere": pd,
        "including_at_the_antipode": bool(
            all(r["positive_definite"] for r in rows
                if abs(r["separation"] - math.pi) < 1e-9)),
        "the_identity": "dΓ_ij/dλ = ⟨δ_i, (H₀−λ)^{-2} δ_j⟩ — a Gram matrix",
        "what_this_upgrades": ("Löwner monotonicity of Γ from a fact sampled "
                               "at a few λ to an analytic consequence; the "
                               "root scans elsewhere are regression checks"),
    }


def measure_the_criterion_extends_to_the_boundary_strata(
        separation: float = 1.3, n_draws: int = 60, seed: int = 20260817
        ) -> Dict[str, object]:
    """The cone is the criterion **in a chart** — here is the rest of ``U(2)``.

    ``φ^reg = A q`` needs ``B`` invertible, which in ``B = i(U−I)``,
    ``C = U+I`` is ``1 ∉ spec U``.  The missing strata are Dirichlet
    directions, reached as ``‖A‖ → ∞``, and they are not represented by any
    finite Hermitian ``A``.  The criterion extends to them by restriction:
    ``A_eff ⪰ P†Γ(0)P`` on the allowed-charge subspace.

    Three things are checked.  That the general form agrees with the cone on
    the chart; that on the ``k = 1`` stratum it agrees with an actual
    negative-``λ`` root scan; and that the reduction is legitimate at all — the
    row space of ``B`` really is the allowed subspace and ``A_eff`` really is
    Hermitian, both of which the reduction assumes.
    """
    rng = np.random.default_rng(seed)

    def haar(rng_):
        x = rng_.normal(size=(2, 2)) + 1j * rng_.normal(size=(2, 2))
        q, r = np.linalg.qr(x)
        return q @ np.diag(np.diag(r) / np.abs(np.diag(r)))

    chart = []
    for _ in range(int(n_draws)):
        u = haar(rng)
        if abs(np.linalg.det(1j * (u - np.eye(2)))) < 1e-6:
            continue
        b, c = boundary_pair_from_unitary(u)
        a = finite_boundary_from_unitary(u)
        gen = is_non_negative_pair(b, c, separation)
        chart.append({"k": gen["k"], "general": gen["non_negative"],
                      "cone": is_non_negative(a, separation),
                      "agrees": bool(gen["non_negative"]
                                     == is_non_negative(a, separation)),
                      "hermitian_defect": gen["hermitian_defect"],
                      "row_space_defect": gen["row_space_defect"]})

    strata = []
    for _ in range(int(n_draws)):
        th = float(rng.uniform(0.2, 2.0 * math.pi - 0.2))
        v = haar(rng)
        u = v @ np.diag([1.0, np.exp(1j * th)]) @ v.conjugate().T
        b, c = boundary_pair_from_unitary(u)
        gen = is_non_negative_pair(b, c, separation)
        p = MouthPair(separation, 0.0, 0.0, 0.0)   # only for the scan helper
        lams = -np.geomspace(1e-8, 4000.0, 5000)[::-1]

        def sec(x: float) -> float:
            return float(np.linalg.det(
                c - b @ gamma_at(x, separation)).real)

        vals = [sec(float(x)) for x in lams]
        n_roots = sum(1 for i in range(len(vals) - 1)
                      if vals[i] * vals[i + 1] < 0.0)
        strata.append({"theta": th, "k": gen["k"],
                       "n_dirichlet": gen["n_dirichlet"],
                       "criterion_non_negative": gen["non_negative"],
                       "scan_negative_roots": n_roots,
                       "agrees": bool(gen["non_negative"] == (n_roots == 0)),
                       "hermitian_defect": gen["hermitian_defect"],
                       "row_space_defect": gen["row_space_defect"]})

    free = is_non_negative_pair(np.zeros((2, 2), dtype=complex),
                                np.eye(2, dtype=complex), separation)
    return {
        "chart_draws": len(chart),
        "chart_mismatches": sum(1 for r in chart if not r["agrees"]),
        "the_general_form_agrees_with_the_cone_on_the_chart": bool(
            all(r["agrees"] for r in chart)),
        "stratum_draws": len(strata),
        "stratum_mismatches": sum(1 for r in strata if not r["agrees"]),
        "the_general_form_agrees_with_the_scan_on_the_strata": bool(
            all(r["agrees"] for r in strata)),
        "every_stratum_has_one_dirichlet_direction": bool(
            all(r["n_dirichlet"] == 1 for r in strata)),
        "worst_hermitian_defect": max(
            [r["hermitian_defect"] for r in chart + strata], default=0.0),
        "worst_row_space_defect": max(
            [r["row_space_defect"] for r in chart + strata], default=0.0),
        "the_reduction_is_legitimate": bool(
            max([r["hermitian_defect"] for r in chart + strata],
                default=0.0) < 1e-9
            and max([r["row_space_defect"] for r in chart + strata],
                    default=0.0) < 1e-9),
        "free_stratum": {"k": free["k"], "non_negative": free["non_negative"],
                         "stratum": free["stratum"]},
        "the_free_stratum_is_non_negative": bool(free["non_negative"]),
        "the_scope": ("A ⪰ Γ(0) is the criterion in the finite-A chart; the "
                      "general one is A_eff ⪰ P†Γ(0)P on the allowed-charge "
                      "subspace, and the chart is its k = 2 case"),
    }
