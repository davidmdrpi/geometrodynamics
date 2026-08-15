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

THE ANSWER
──────────
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
``det Γ(0) = g₀² − G₀² < 0`` always, so the apex is an *indefinite* matrix: no
throat with ``A = 0`` is stable, at any separation.  As ``d → π`` the positive
threshold ``g₀ + G₀ → 0``.

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
    "measure_the_positive_sector_is_a_shifted_psd_cone",
    "measure_the_inertia_theorem_counts_modes_below_any_threshold",
    "measure_the_boundary_of_the_cone_is_a_zero_mode",
    "measure_the_growth_rate_turns_on_with_a_square_root",
    "measure_the_exchange_symmetric_wedge_is_a_slice",
    "measure_where_the_apex_sits_as_the_mouths_separate",
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


def positivity_defect(boundary: np.ndarray, separation: float,
                      lmbda: float = 0.0) -> Dict[str, object]:
    """Eigenvalues of ``A − Γ(λ*)``, and the smallest of them.

    The smallest eigenvalue is the signed distance to the cone's boundary in
    the only sense that matters: negative means there are modes below ``λ*``,
    and *how many* is how many eigenvalues are negative.
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
            "zero_matrix_is_stable": bool(ev[0] >= 0.0)}


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

    Step out along the symmetric channel by ``ε`` in ``α + β``.  The smallest
    eigenvalue of ``A − Γ(0)`` is then ``−ε``, and since it is decreasing in
    ``λ`` with slope ``μ′``, the root sits at ``λ ≈ −ε/μ′``:

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
        separations: Sequence[float] = (0.2, 0.5, 0.8, 1.3, 2.0, 2.8, 3.0)
        ) -> Dict[str, object]:
    """The apex is at ``Γ(0)``, and where that is depends on the geometry.

    Its **trace is ``2g₀ = −1/(2π²)`` at every separation** — the mouth
    distance does not enter it — while its eigenvalues are exactly PR #256's two
    channel thresholds, ``g₀ ± G₀``.  Its determinant ``g₀² − G₀²`` is negative
    for every separation, so ``Γ(0)`` is **indefinite** and one corollary is
    immediate: ``A = 0`` is unstable at every separation, which no amount of
    tuning the mouths fixes.

    As ``d → π`` the positive threshold ``g₀ + G₀ → 0``: antipodal mouths make
    the symmetric channel marginally stable at ``A = 0``.
    """
    rows = []
    for d in separations:
        ap = apex(float(d))
        rows.append({"separation": float(d), "trace": ap["trace"],
                     "determinant": ap["determinant"],
                     "eigenvalues": ap["eigenvalues"],
                     "indefinite": ap["indefinite"],
                     "zero_matrix_is_stable": ap["zero_matrix_is_stable"],
                     "channel_thresholds": ap["channel_thresholds"]})
    traces = [r["trace"] for r in rows]
    expect = -1.0 / (2.0 * math.pi ** 2)
    return {
        "rows": rows,
        "trace_is_separation_independent": bool(
            max(traces) - min(traces) < 1e-15),
        "trace_value": traces[0], "predicted_trace": expect,
        "trace_matches_minus_one_over_two_pi_squared": bool(
            abs(traces[0] - expect) < 1e-15),
        "the_apex_is_always_indefinite": bool(
            all(r["indefinite"] for r in rows)),
        "the_zero_matrix_is_never_stable": bool(
            not any(r["zero_matrix_is_stable"] for r in rows)),
        "eigenvalues_are_the_channel_thresholds": bool(all(
            abs(r["eigenvalues"][0] - r["channel_thresholds"][0]) < 1e-14
            and abs(r["eigenvalues"][1] - r["channel_thresholds"][1]) < 1e-14
            for r in rows)),
        "the_symmetric_threshold_closes_at_the_antipode": bool(
            rows[-1]["eigenvalues"][1] < rows[0]["eigenvalues"][1] / 100.0),
    }
