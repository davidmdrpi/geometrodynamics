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


def bilinearity_in_amplitude(degrees: Sequence[int], samples: int = 40
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
        full = float(np.max(np.abs(dipole_from_improved_stress(state))))
        half = float(np.max(np.abs(dipole_from_improved_stress(halved))))
        if full > 1e-8:
            worst = max(worst, abs(half / full - 0.5))
    return {"max_halving_error": worst, "route": "improved_stress"}


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
    # P^A is the symmetrized cross term; for a != b the cross term appears twice
    flat = (2.0 * coeff * T).transpose(2, 0, 1).reshape(dim_b, -1).T
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


def verdict(checks: Dict[str, bool], momentum_note: str) -> Dict[str, object]:
    """Stable schema on every path; a failure names itself and blocks all fields."""
    if not checks or not all(checks.values()):
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"),
                "failed_checks": sorted(k for k, v in checks.items() if not v)
                                 or ["no checks supplied"]}
    return {"dipole_obstruction_structure":
                "BILINEAR_CROSS_PARITY_ADJACENT_DEGREE_ONLY",
            "antipodal_parity_status": "SUFFICIENT_NOT_NECESSARY",
            "f6_consequence":
                "CONSTRAINT_SOLVABILITY_DOES_NOT_DERIVE_THE_ANTIPODAL_CONDITION",
            "momentum_sector": momentum_note,
            "triangle_map": "NOT_DERIVED", "readout": "NOT_DERIVED",
            "failed_checks": []}


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
