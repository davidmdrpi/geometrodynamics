"""Does the Hamiltonian constraint's solvability force antipodal parity?

Public freeze ``495f1f1`` (``docs/parity_solvability_prereg.md``), committed
and published before this file. The linearized constraint
``(Delta + 3/a^2) u = -kappa delta_rho/4`` on the round ``S^3`` has a genuine
four-dimensional ``l=1`` kernel, so a solution exists only when the source is
orthogonal to the four ambient coordinates ``x^A``.

The question is whether the two antipodal parity eigenspaces are the maximal
linear subspaces on which that obstruction vanishes identically. They are
not: the obstruction is carried by a degree-1 triple overlap, which is
nonzero only for adjacent degrees.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Iterable, List, Sequence, Tuple
import math

import numpy as np

from geometrodynamics.waves.reciprocal_scalar_tt import (
    ReciprocalModel, harmonic_multiplet, monomials)


PUBLIC_PREREG = "495f1f185411a80cfd688b0b33291d1986f9702a"
SEED = 2026090713
MAX_DEGREE = 6


def lambda_n(n: int, radius: float = 1.0) -> float:
    """Scalar Laplacian eigenvalue ``n(n+2)/a^2`` on the round ``S^3``."""
    return n * (n + 2) / radius ** 2


def dipole_coefficient(n: int, m: int, radius: float = 1.0) -> float:
    """Frozen P2 coefficient ``(n(n+2) + m(m+2) + 1) / (2 a^2)``."""
    return (n * (n + 2) + m * (m + 2) + 1) / (2.0 * radius ** 2)


def _degree_one_basis() -> np.ndarray:
    """The four ``l=1`` kernel modes, orthonormal for unit-``S^3`` Haar."""
    return harmonic_multiplet(1).B


@lru_cache(maxsize=None)
def triple_overlap(n: int, m: int) -> np.ndarray:
    """``T[A, i, j] = int Y_1^A Y_n^i Y_m^j dV`` over normalized Haar.

    Evaluated from the exact monomial moments of the inherited multiplet
    machinery, so no quadrature rule enters. The result is the object the
    freeze predicts to vanish unless ``|n - m| = 1``.
    """
    mult_n, mult_m, mult_1 = (harmonic_multiplet(k) for k in (n, m, 1))
    from geometrodynamics.waves.reciprocal_scalar_tt import sphere_moment
    e1 = np.array(mult_1.exponents)
    en = np.array(mult_n.exponents)
    em = np.array(mult_m.exponents)
    # raw monomial-moment tensor, then rotate into the orthonormal mode bases
    raw = np.zeros((len(e1), len(en), len(em)))
    for a in range(len(e1)):
        for b in range(len(en)):
            base = e1[a] + en[b]
            for c in range(len(em)):
                raw[a, b, c] = sphere_moment(tuple(base + em[c]))
    overlap = np.einsum("pa,qb,rc,abc->pqr", mult_1.B.T, mult_n.B.T, mult_m.B.T, raw)
    overlap.flags.writeable = False        # cached; never mutate in place
    return overlap


def selection_rule_table(max_degree: int = MAX_DEGREE) -> Dict[str, object]:
    """Full ``|n - m|`` table of triple-overlap norms, adjacent and not.

    Reports every pair, not only the vanishing ones: a table of zeros alone
    would not establish a selection rule.
    """
    rows = []
    for n in range(1, max_degree + 1):
        for m in range(1, max_degree + 1):
            norm = float(np.linalg.norm(triple_overlap(n, m)))
            rows.append({"n": n, "m": m, "gap": abs(n - m),
                         "adjacent": abs(n - m) == 1,
                         "overlap_norm": norm})
    adjacent = [r["overlap_norm"] for r in rows if r["adjacent"]]
    other = [r["overlap_norm"] for r in rows if not r["adjacent"]]
    return {"rows": rows,
            "min_adjacent_norm": min(adjacent) if adjacent else 0.0,
            "max_non_adjacent_norm": max(other) if other else 0.0,
            "selection_rule_holds": bool(
                adjacent and min(adjacent) > 1e-3 and max(other) < 1e-12)}


def dipole_from_overlap(field: Dict[int, np.ndarray],
                        momentum: Dict[int, np.ndarray] | None = None,
                        radius: float = 1.0) -> np.ndarray:
    """``P^A`` from the frozen P2 reduction: overlaps and the closed coefficient."""
    momentum = momentum or {}
    degrees = sorted(set(field) | set(momentum))
    total = np.zeros(4)
    for n in degrees:
        for m in degrees:
            T = triple_overlap(n, m)
            if n in field and m in field:
                total += dipole_coefficient(n, m, radius) * np.einsum(
                    "pqr,q,r->p", T, field[n], field[m])
            if n in momentum and m in momentum:
                total += np.einsum("pqr,q,r->p", T, momentum[n], momentum[m])
    # Correction C1, found by the frozen P2 cross-check against the inherited
    # improved stress. The P2 formula is already the TOTAL cross contribution,
    # so summing over ordered pairs (n,m) and (m,n) counts it twice; the
    # diagonal n=m term is likewise half the bilinear form. One overall factor
    # of 1/2 fixes both. The two routes disagreed by exactly 2.000 before this.
    return 0.5 * total


def _jets_on_points(degree: int, coefficients: np.ndarray, points: np.ndarray):
    return harmonic_multiplet(degree).polynomial_jets(points, coefficients)


def dipole_from_improved_stress(field: Dict[int, np.ndarray],
                                momentum: Dict[int, np.ndarray] | None = None,
                                radial_order: int = 12,
                                angular_order: int = 24) -> np.ndarray:
    """``P^A`` from the inherited improved stress, integrated independently.

    This never uses :func:`triple_overlap`; it evaluates the pointwise
    energy density of `reciprocal_scalar_tt` and integrates against the four
    kernel modes. It is the independent check of the P2 reduction.
    """
    from geometrodynamics.waves.reciprocal_scalar_tt import (
        sphere_quadrature, scalar_jets, ReciprocalModel as RM)
    from geometrodynamics.waves.backreaction import stress_series
    momentum = momentum or {}
    points, weights = sphere_quadrature(radial_order, angular_order)
    degrees = sorted(set(field) | set(momentum))
    model = RM(degree=degrees[0])
    volume = model.volume

    total_phi = np.zeros(len(points))
    total_dt = np.zeros(len(points))
    total_grad = np.zeros((len(points), 3))
    total_dtgrad = np.zeros((len(points), 3))
    total_hess = np.zeros((len(points), 3, 3))
    for n in degrees:
        q = field.get(n, np.zeros(harmonic_multiplet(n).dimension))
        p = momentum.get(n, np.zeros(harmonic_multiplet(n).dimension))
        sub = RM(degree=n)
        jets, _ = scalar_jets(sub, q, p, -sub.omega_scalar2 * q, points)
        total_phi += jets["phi"][:, 0]
        total_dt += jets["dt"][:, 0]
        total_grad += jets["grad"][:, 0, :]
        total_dtgrad += jets["dtgrad"][:, 0, :]
        total_hess += jets["hess"][:, 0, :, :]
    jets = {"phi": total_phi[:, None], "dt": total_dt[:, None],
            "dtt": np.zeros_like(total_phi)[:, None],
            "grad": total_grad[:, None, :], "dtgrad": total_dtgrad[:, None, :],
            "hess": total_hess[:, None, :, :],
            "laplacian": np.trace(total_hess, axis1=1, axis2=2)[:, None]}
    rho = stress_series(jets)[:, 0, 0, 0]
    kernel = monomials(points, np.array(harmonic_multiplet(1).exponents)) @ _degree_one_basis()
    del volume
    return np.einsum("p,p,pa->a", weights, rho, kernel)


def parity_split(degrees: Iterable[int]) -> Tuple[List[int], List[int]]:
    ds = sorted(degrees)
    return [d for d in ds if d % 2 == 0], [d for d in ds if d % 2 == 1]


def has_adjacent_pair(degrees: Iterable[int]) -> bool:
    ds = sorted(set(degrees))
    return any(b - a == 1 for a, b in zip(ds[:-1], ds[1:]))


def random_state(degrees: Sequence[int], rng, amplitude: float = 1.0
                 ) -> Dict[int, np.ndarray]:
    out = {}
    for n in degrees:
        v = rng.normal(size=harmonic_multiplet(n).dimension)
        out[n] = amplitude * v / np.linalg.norm(v)
    return out


def scan_degree_sets(sets: Sequence[Sequence[int]], samples: int = 200,
                     amplitudes: Sequence[float] = (1.0, 0.5, 0.25)
                     ) -> List[Dict[str, object]]:
    """Sample each degree set and record the largest dipole found."""
    rng = np.random.default_rng(SEED)
    rows = []
    for degrees in sets:
        worst = 0.0
        for _ in range(samples):
            for amp in amplitudes:
                state = random_state(degrees, rng, amp)
                worst = max(worst, float(np.max(np.abs(dipole_from_overlap(state)))))
        even, odd = parity_split(degrees)
        rows.append({"degrees": list(degrees), "mixed_parity": bool(even and odd),
                     "has_adjacent_pair": has_adjacent_pair(degrees),
                     "max_abs_dipole": worst})
    return rows


def quadrature_exactness(degrees: Sequence[int] = (2, 3)) -> Dict[str, float]:
    """The sphere rules are exact for these polynomial integrands.

    The improved-stress dipole integrand is a polynomial, so the
    Gauss-Legendre times trapezoid rule reproduces it exactly once the order
    is high enough. Demonstrating that licenses the coarse grid for the
    high-sample bilinearity scan without loss of fidelity.
    """
    rng = np.random.default_rng(SEED + 6)
    state = random_state(degrees, rng)
    fine = dipole_from_improved_stress(state, None, 12, 24)
    coarse = dipole_from_improved_stress(state, None, 6, 12)
    scale = max(float(np.max(np.abs(fine))), 1e-30)
    return {"fine_vs_coarse_absolute": float(np.max(np.abs(fine - coarse))),
            "fine_vs_coarse_relative": float(np.max(np.abs(fine - coarse))) / scale,
            "exact": bool(np.max(np.abs(fine - coarse)) / scale < 1e-12)}


def bilinearity_in_amplitude(degrees: Sequence[int], samples: int = 40,
                             radial_order: int = 6, angular_order: int = 12
                             ) -> Dict[str, float]:
    """Halving one parity component must halve the dipole, not merely shrink it.

    Guards against reporting a generic nonzero number as the obstruction.
    """
    rng = np.random.default_rng(SEED + 1)
    even, odd = parity_split(degrees)
    worst = 0.0
    for _ in range(samples):
        state = random_state(degrees, rng)
        halved = dict(state)
        for n in odd:
            halved[n] = state[n] / 2.0
        # Measured through the INDEPENDENT improved-stress route. Using the
        # overlap route here would be bilinear by construction and could not
        # fail, so it would be a structural regression rather than evidence.
        full = float(np.max(np.abs(dipole_from_improved_stress(
            state, None, radial_order, angular_order))))
        half = float(np.max(np.abs(dipole_from_improved_stress(
            halved, None, radial_order, angular_order))))
        if full > 1e-8:
            worst = max(worst, abs(half / full - 0.5))
    return {"max_halving_error": worst, "route": "improved_stress",
            "grid": [radial_order, angular_order]}


def subspace_maximality(degree_a: int = 1, degree_b: int = 2) -> Dict[str, object]:
    """Can any nonzero subspace of ``V_b`` evade the obstruction against all of ``V_a``?

    For each ``u`` in ``V_b`` the map ``v -> P(v, u)`` is linear on ``V_a``.
    Its matrix is linear in ``u``, so the set of evading ``u`` is the kernel
    of one linear map. A trivial kernel means no subspace ``V_a + U`` with
    ``U`` nonzero avoids the dipole.
    """
    T = triple_overlap(degree_a, degree_b)
    coeff = dipole_coefficient(degree_a, degree_b)
    dim_a = harmonic_multiplet(degree_a).dimension
    dim_b = harmonic_multiplet(degree_b).dimension
    # Correction C2: this retained the factor of two that C1 removed from the
    # dipole normalization. It rescales every singular value and cannot change
    # the kernel, but the reported magnitude was twice the consistent one.
    flat = (coeff * T).transpose(2, 0, 1).reshape(dim_b, -1).T
    singular = np.linalg.svd(flat, compute_uv=False)
    return {"degrees": [degree_a, degree_b], "dim_a": dim_a, "dim_b": dim_b,
            "singular_values": [float(s) for s in singular],
            "smallest_singular_value": float(singular.min()),
            "kernel_is_trivial": bool(singular.min() > 1e-9)}


def graph_subspace_search(degree_a: int = 1, degree_b: int = 2) -> Dict[str, object]:
    """Is there a nonzero ``T: V_a -> V_b`` whose graph has vanishing dipole?

    On the graph ``{v + Tv}`` the dipole is the quadratic form
    ``v -> 2 c P(v, Tv)``; it vanishes identically exactly when the
    symmetric part vanishes. That is a linear condition on ``T``, so the
    answer is the nullity of an explicit matrix rather than a search.
    """
    T = triple_overlap(degree_a, degree_b)
    dim_a = harmonic_multiplet(degree_a).dimension
    dim_b = harmonic_multiplet(degree_b).dimension
    rows = []
    for i in range(dim_a):
        for j in range(dim_a):
            if j < i:
                continue
            # symmetric part of v_i v_j coefficient, as a functional of T
            selector = np.zeros((4, dim_a, dim_b))
            if i == j:
                selector[:, i, :] = T[:, i, :]
            else:
                selector[:, j, :] = 0.5 * T[:, i, :]
                selector[:, i, :] = selector[:, i, :] + 0.5 * T[:, j, :]
            rows.append(selector.reshape(4, dim_a * dim_b))
    matrix = np.concatenate(rows, axis=0)
    singular = np.linalg.svd(matrix, compute_uv=False)
    nullity = int(np.sum(singular < 1e-9)) + max(0, dim_a * dim_b - len(singular))
    return {"degrees": [degree_a, degree_b],
            "matrix_shape": list(matrix.shape),
            "smallest_singular_value": float(singular.min()),
            "nullity": nullity,
            "only_trivial_graph": bool(nullity == 0)}


VERDICT_FIELDS = ("dipole_obstruction_structure", "antipodal_parity_status",
                  "f6_consequence", "momentum_sector", "triangle_map", "readout")

#: Correction C3. The first version gated only on ``all(checks.values())``, so
#: ``verdict({"unrelated": True}, ...)`` returned the full affirmative result
#: and a probe that silently dropped every real check would still have passed.
#: Presence is now required as well as truth.
REQUIRED_CHECKS = (
    "selection rule holds for every pair up to degree 6",
    "predicted free mixed-parity sets have no obstruction",
    "predicted obstructed sets do obstruct",
    "parity-pure sets have no obstruction",
    "overlap reduction matches the inherited improved stress",
    "the obstruction is bilinear, measured independently",
    "no nonzero subspace evades an adjacent partner",
    "no nonzero graph subspace evades an adjacent pair",
    "a small-projection mixed-parity subspace evades an adjacent pair",
    "the complete six Killing and four gradient charges are audited",
    "momentum charges are not a parity condition",
)


def verdict(checks: Dict[str, bool], momentum_note: str) -> Dict[str, object]:
    """Stable schema on every path; a failure names itself and blocks all fields."""
    missing = [k for k in REQUIRED_CHECKS if k not in checks]
    failed = sorted(k for k, v in checks.items() if not v)
    if missing or failed or not checks:
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"),
                "failed_checks": failed or ["no checks supplied"],
                "missing_checks": missing}
    return {"dipole_obstruction_structure":
                "BILINEAR_CROSS_PARITY_ADJACENT_DEGREE_ONLY",
            "antipodal_parity_status": "SUFFICIENT_NOT_NECESSARY",
            "f6_consequence":
                "CONSTRAINT_SOLVABILITY_DOES_NOT_DERIVE_THE_ANTIPODAL_CONDITION",
            "momentum_sector": momentum_note,
            "triangle_map": "NOT_DERIVED", "readout": "NOT_DERIVED",
            "failed_checks": [], "missing_checks": []}


def route_agreement(field: Dict[int, np.ndarray],
                    momentum: Dict[int, np.ndarray] | None = None,
                    floor: float = 1e-9) -> Dict[str, float]:
    """Compare the two routes with a scale floor, not a bare ratio.

    A pure relative measure is meaningless when both routes return numbers at
    rounding level, which is exactly the predicted-zero case. Below ``floor``
    the absolute difference is the meaningful quantity.
    """
    overlap = dipole_from_overlap(field, momentum)
    stress = dipole_from_improved_stress(field, momentum)
    absolute = float(np.max(np.abs(overlap - stress)))
    scale = max(float(np.max(np.abs(overlap))), float(np.max(np.abs(stress))))
    return {"absolute": absolute, "scale": scale,
            "relative": absolute / scale if scale > floor else 0.0,
            "below_floor": bool(scale <= floor)}


def killing_charges(field: Dict[int, np.ndarray],
                    momentum: Dict[int, np.ndarray]) -> np.ndarray:
    """The six ``SO(4)`` Killing charges ``-p^T D_i q`` of the momentum sector.

    The momentum constraint's own solvability conditions come from the Killing
    fields, not from the ``l=1`` scalar kernel. The generators act within a
    single multiplet, so these charges are diagonal in degree.
    """
    charges = []
    for n in sorted(set(field) & set(momentum)):
        mult = harmonic_multiplet(n)
        q, p = field[n], momentum[n]
        charges.append(np.array([-p @ (D @ q) for D in mult.D]))
    if not charges:
        return np.zeros(3)
    return np.sum(charges, axis=0)


def momentum_sector_report() -> Dict[str, object]:
    """Is the momentum obstruction a parity condition? Measured, not predicted.

    The freeze deliberately fixes no prediction here. What is reported is
    whether a *parity-pure* multiplet — which always kills the Hamiltonian
    dipole — can still carry a nonzero Killing charge.
    """
    rng = np.random.default_rng(SEED + 2)
    rows = []
    for n in (1, 2, 3, 4):
        mult = harmonic_multiplet(n)
        q = rng.normal(size=mult.dimension)
        q /= np.linalg.norm(q)
        aligned = mult.D[0] @ q                     # p proportional to D_1 q
        random_p = rng.normal(size=mult.dimension)
        random_p /= np.linalg.norm(random_p)
        rows.append({
            "degree": n, "parity": "even" if n % 2 == 0 else "odd",
            "hamiltonian_dipole": float(np.max(np.abs(
                dipole_from_overlap({n: q}, {n: aligned})))),
            "killing_charge_aligned": float(np.max(np.abs(
                killing_charges({n: q}, {n: aligned})))),
            "killing_charge_random": float(np.max(np.abs(
                killing_charges({n: q}, {n: random_p})))),
        })
    return {"rows": rows,
            "parity_pure_can_carry_charge": bool(
                max(r["killing_charge_aligned"] for r in rows) > 1e-3),
            "parity_pure_kills_dipole": bool(
                max(r["hamiltonian_dipole"] for r in rows) < 1e-12),
            "note": "NOT_PREDICTED_IN_ADVANCE"}


def _monomial_mode(degree: int, exponent: Tuple[int, ...]) -> np.ndarray:
    """Modal coefficients of a single ambient monomial inside one multiplet."""
    from geometrodynamics.waves.reciprocal_scalar_tt import sphere_moment
    mult = harmonic_multiplet(degree)
    exps = np.array(mult.exponents)
    target = np.zeros(len(exps))
    for i, e in enumerate(exps):
        if tuple(int(v) for v in e) == exponent:
            target[i] = 1.0
    moments = np.array([[sphere_moment(tuple(a + b)) for b in exps] for a in exps])
    return mult.B.T @ (moments @ target)


def small_projection_counterexample() -> Dict[str, object]:
    """A mixed-parity subspace inside ADJACENT degrees with no obstruction.

    Correction C4. The first version claimed the parity answer "fails only
    globally", on the strength of two maximality tests. Both enlarge or map
    *all* of ``V_n``, so neither family contains a subspace whose projection
    into ``V_n`` is smaller. Such subspaces exist:

        S = span{ x_0, x_1 x_2 }  in  V_1 + V_2,

    for which ``int x^A x_0 x_1 x_2 dV = 0`` for every ambient coordinate, so
    the dipole vanishes on all of ``S`` at every amplitude. The parity
    eigenspaces are therefore not maximal even within an adjacent truncation.
    This strengthens the round's negative conclusion rather than weakening it.
    """
    rng = np.random.default_rng(SEED + 4)
    q1 = _monomial_mode(1, (1, 0, 0, 0))          # x_0
    q2 = _monomial_mode(2, (0, 1, 1, 0))          # x_1 x_2
    worst_overlap = worst_stress = 0.0
    for _ in range(200):
        a, b = rng.normal(size=2) * rng.choice([0.2, 1.0, 5.0])
        state = {1: a * q1, 2: b * q2}
        worst_overlap = max(worst_overlap,
                            float(np.max(np.abs(dipole_from_overlap(state)))))
        worst_stress = max(worst_stress,
                           float(np.max(np.abs(dipole_from_improved_stress(state)))))
    generic = 0.0
    for _ in range(100):
        generic = max(generic, float(np.max(np.abs(
            dipole_from_overlap(random_state([1, 2], rng))))))
    return {"subspace": "span{x_0, x_1 x_2} in V_1 + V_2",
            "mixed_parity": True, "degrees_adjacent": True,
            "max_dipole_overlap_route": worst_overlap,
            "max_dipole_stress_route": worst_stress,
            "generic_adjacent_dipole": generic,
            "evades": bool(max(worst_overlap, worst_stress) < 1e-12
                           and generic > 1e-3)}


def killing_generators(degree: int) -> np.ndarray:
    """All six ``so(4)`` generators acting on one multiplet.

    Correction C5. The first version used only ``mult.D``, three of the six.
    ``D`` spans one ``su(2)`` factor and ``rotation_generators`` the diagonal
    ``so(3)``; the two intersect trivially, so together they span ``so(4)``.
    """
    mult = harmonic_multiplet(degree)
    return np.concatenate([np.asarray(mult.D),
                           np.asarray(mult.rotation_generators)], axis=0)


def killing_charges(field: Dict[int, np.ndarray],
                    momentum: Dict[int, np.ndarray]) -> np.ndarray:
    """All six ``SO(4)`` Killing charges ``-p^T G q``.

    A Killing field is divergence free, so the improved-stress correction
    ``-(1/6) d_i d_t(phi^2)`` integrates away against it and the charge is the
    canonical one. That is not true of the gradient conformal fields; see
    :func:`gradient_conformal_charges`.
    """
    charges = np.zeros(6)
    for n in sorted(set(field) & set(momentum)):
        q, p = field[n], momentum[n]
        charges += np.array([-p @ (G @ q) for G in killing_generators(n)])
    return charges


def gradient_conformal_charges(field: Dict[int, np.ndarray],
                               momentum: Dict[int, np.ndarray],
                               radial_order: int = 12,
                               angular_order: int = 24) -> np.ndarray:
    """The four gradient conformal-Killing charges ``int j . grad(x^A) dV``.

    These need the full improved momentum density, because ``grad(x^A)`` is not
    divergence free: ``div grad(x^A) = -lambda_1 x^A``. Evaluated pointwise on
    the inherited sphere quadrature rather than through a modal shortcut.
    """
    from geometrodynamics.waves.reciprocal_scalar_tt import (
        sphere_quadrature, scalar_jets, ReciprocalModel as RM)
    points, weights = sphere_quadrature(radial_order, angular_order)
    degrees = sorted(set(field) | set(momentum))
    phi = np.zeros(len(points))
    dt = np.zeros(len(points))
    grad = np.zeros((len(points), 3))
    dtgrad = np.zeros((len(points), 3))
    for n in degrees:
        mult = harmonic_multiplet(n)
        q = field.get(n, np.zeros(mult.dimension))
        p = momentum.get(n, np.zeros(mult.dimension))
        sub = RM(degree=n)
        jets, _ = scalar_jets(sub, q, p, -sub.omega_scalar2 * q, points)
        phi += jets["phi"][:, 0]
        dt += jets["dt"][:, 0]
        grad += jets["grad"][:, 0, :]
        dtgrad += jets["dtgrad"][:, 0, :]
    # improved momentum density j_i = -T_0i = -(2/3) phidot d_i phi
    #                                        + (1/3) phi d_i phidot
    current = -(2.0 / 3.0) * dt[:, None] * grad + (1.0 / 3.0) * phi[:, None] * dtgrad
    one = RM(degree=1)
    charges = []
    for basis in np.eye(harmonic_multiplet(1).dimension):
        jets, _ = scalar_jets(one, basis, np.zeros_like(basis),
                              -one.omega_scalar2 * basis, points)
        charges.append(float(np.einsum("p,pi,pi->", weights, current,
                                       jets["grad"][:, 0, :])))
    return np.array(charges)


def complete_momentum_audit() -> Dict[str, object]:
    """Frozen check 6, delivered in full: six Killing and four gradient charges.

    Correction C5. The first version reported three of the six Killing charges
    and none of the four gradient conformal charges, so it could report "zero"
    while an actual charge was large.

    Completing it exposes the structure. The Killing charges are **diagonal in
    degree** — a Killing field is divergence free, the improvement term drops,
    and the canonical charge ``-p^T G q`` needs field and momentum in the same
    multiplet. The gradient conformal charges are not: ``grad(x^A)`` has
    divergence ``-lambda_1 x^A``, the improvement survives, and the charge is a
    field-momentum pairing obeying the same **adjacency** rule as the
    Hamiltonian dipole. Neither is a parity condition.
    """
    rng = np.random.default_rng(SEED + 5)
    mult = harmonic_multiplet(3)
    q = rng.normal(size=mult.dimension); q /= np.linalg.norm(q)

    # (a) data annihilating the three originally reported generators
    rows = np.array([D @ q for D in np.asarray(mult.D)])
    p = rng.normal(size=mult.dimension)
    p = p - rows.T @ np.linalg.pinv(rows @ rows.T) @ rows @ p
    p /= np.linalg.norm(p)
    partial = np.array([-p @ (D @ q) for D in np.asarray(mult.D)])
    full = killing_charges({3: q}, {3: p})

    # (b) field and momentum in disjoint ADJACENT degrees: every Killing charge
    #     and the Hamiltonian dipole vanish, yet a gradient charge survives.
    q2 = rng.normal(size=harmonic_multiplet(2).dimension); q2 /= np.linalg.norm(q2)
    p3 = rng.normal(size=harmonic_multiplet(3).dimension); p3 /= np.linalg.norm(p3)
    split_field, split_momentum = {2: q2}, {3: p3}
    split_killing = killing_charges(split_field, split_momentum)
    split_gradient = gradient_conformal_charges(split_field, split_momentum)
    split_dipole = dipole_from_overlap(split_field, split_momentum)

    # (c) the gradient charge obeys adjacency, not parity
    far_gradient = gradient_conformal_charges({2: q2}, {5: rng.normal(
        size=harmonic_multiplet(5).dimension)})
    same_gradient = gradient_conformal_charges({3: q}, {3: p})

    return {
        "originally_reported_three": float(np.max(np.abs(partial))),
        "complete_six_on_the_same_data": float(np.max(np.abs(full))),
        "incomplete_audit_would_report_zero": bool(np.max(np.abs(partial)) < 1e-10
                                                   and np.max(np.abs(full)) > 1e-3),
        "split_degree_killing": float(np.max(np.abs(split_killing))),
        "split_degree_hamiltonian_dipole": float(np.max(np.abs(split_dipole))),
        "split_degree_gradient_charge": float(np.max(np.abs(split_gradient))),
        "gradient_sector_is_independent": bool(
            np.max(np.abs(split_killing)) < 1e-12
            and np.max(np.abs(split_dipole)) < 1e-12
            and np.max(np.abs(split_gradient)) > 1e-3),
        "gradient_same_degree": float(np.max(np.abs(same_gradient))),
        "gradient_non_adjacent": float(np.max(np.abs(far_gradient))),
        "gradient_obeys_adjacency": bool(np.max(np.abs(same_gradient)) < 1e-12
                                         and np.max(np.abs(far_gradient)) < 1e-12
                                         and np.max(np.abs(split_gradient)) > 1e-3),
        "killing_is_diagonal_in_degree": bool(np.max(np.abs(split_killing)) < 1e-12
                                              and np.max(np.abs(full)) > 1e-3),
    }
