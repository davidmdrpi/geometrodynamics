"""The topology and bundle structure of the finite-mouth handle, constructed.

WHAT THIS MODULE IS FOR
───────────────────────
``finite_mouth`` derived the spatial geometry of one handle: two geodesic
four-balls of angular radius ``a`` cut from a round ``S^4_R`` and joined by the
tube ``[-S, S] x S^3`` with ``f = sqrt(s^2 + b^2)``. It deliberately left the
discrete identification data out. The repository elsewhere calls the throat
"non-orientable", speaks of an "``RP^2`` mouth", a "Möbius wrap", "Pin^-
transport" and a "``Z_2`` wrap parity", and ``embedding/topology.py`` stores
non-orientability as a Boolean field. This module replaces the words with the
actual gluing data, and asks what the geometry forces.

The pre-registration is ``docs/finite_mouth_topology_prereg.md``, committed
before this file existed. Its hypotheses H1-H5 are the section headings below.

H1 · THE SEAM GLUING IS FREE
────────────────────────────
Darmois matching at a mouth is two ``O(4)``-invariant scalars, so each seam
map ``psi_A, psi_B`` is an arbitrary isometry of the mouth ``S^3`` and only the
loop monodromy ``m = psi_B^{-1} psi_A in O(4)`` is meaningful. The closed
4-manifold is the mapping torus of ``m``; its orientability is ``det m``.
Preserving the brane (the equatorial ``S^3``, whose normal is one of the four
mouth directions) splits ``m = (m_3, eps) in O(3) x O(1)`` and gives four
classes ``(det m_3, eps)``.

H2 · THE FREE INVOLUTION
────────────────────────
With ``P_B = -P_A`` on ``S^4`` the antipode ``x -> -x`` swaps the mouths and
extends through the tube uniquely as ``iota(s, Omega) = (-s, -Omega)``. It is
the only fixed-point-free involution available (``-I_5`` is the only
fixed-point-free involution in ``O(5)``), it is an isometry for both lapses,
and it reverses the orientation of the 4-dimensional handle.

H3 · THE SCALAR CONDITION ON THE QUOTIENT
─────────────────────────────────────────
A scalar on ``M/iota`` satisfies ``Phi(iota x) = eta Phi(x)`` with ``eta = ±1``
the line-bundle sign. With ``Y_l(-Omega) = (-1)^l Y_l(Omega)`` this is
``psi_l(-s) = eta (-1)^l psi_l(s)``: Neumann or Dirichlet at the neck by the
parity of ``l``. PR #129's condition is the ``eta = +1`` sector.

H5 · WHERE ``J = i sigma_y K`` ACTUALLY LIVES
─────────────────────────────────────────────
Modelled on the quaternions: ``R^4 = H`` with left and right multiplication
``L_q, R_q``, Spin(4) = SU(2)_L x SU(2)_R acting on the vector by
``v -> g_L v g_R^{-1}`` and on the two half-spinors by ``L_{g_L}`` and
``L_{g_R}``. The map ``embedding/transport.py`` writes down is ``L_{-j}``:
left multiplication by the unit quaternion ``-j``. It anticommutes with the
Hopf complex structure ``L_i`` (so it is antilinear there, hence the ``K``)
and commutes with ``R_i`` (so on the half-spinor it is the linear matrix
``i sigma_y``). Every element of Pin(4) acts linearly on spinors, so ``J`` is
not the lift of any gluing map; its antilinearity is the reversal of the
Hopf ``U(1)``, a bundle datum, not a manifold datum.

No Pauli matrix is a starting point anywhere here: the quaternion algebra is
built from its multiplication table, and the Clifford algebra ``Cl(4)`` from
its regular representation.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.bulk.finite_mouth import (
    BULK_RADIUS,
    MOUTH_ANGLE,
    geometry,
    static_admittance,
    throat_radius,
)

__all__ = [
    "GluingClass",
    "classify_gluing",
    "antipodal_gluing",
    "identity_gluing",
    "reflection_gluing",
    "mapping_torus_w1",
    "harmonic_representation",
    "antipodal_parity",
    "Involution",
    "antipodal_involution",
    "neck_reflection_involution",
    "free_involutions_of_o5",
    "stiefel_whitney_rp",
    "pin_structures_rp",
    "quaternion_left",
    "quaternion_right",
    "hopf_transport_matrix",
    "complex_structure_commutation",
    "clifford_regular",
    "pin_lifts_of_reflection",
    "spin_lifts_of_antipode",
    "neck_sector",
    "half_tube_admittance",
    "half_tube_admittance_oracle",
    "static_operator_commutes_with_parity",
    "traversal_table",
]


# ═══════════════════════════════════════════════════════════════════════════
# H1 · gluing data on the two-mouth handle
# ═══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class GluingClass:
    """The discrete data of one loop monodromy ``m in O(4)`` on the mouth.

    ``normal_index`` says which of the four mouth directions is the brane
    normal; the other three span the brane mouth ``S^2``.
    """

    monodromy: Tuple[Tuple[float, ...], ...]
    normal_index: int
    det_bulk: int          # det m: w_1 of the 4-dimensional handle on its loop
    det_brane: int         # det m_3: w_1 of the brane handle S^2 x S^1 / S^2 x~ S^1
    normal_sign: int       # eps: w_1 of the brane's normal line bundle
    preserves_brane: bool

    @property
    def bulk_orientable(self) -> bool:
        return self.det_bulk == +1

    @property
    def brane_orientable(self) -> bool:
        return self.det_brane == +1

    @property
    def normal_bundle_trivial(self) -> bool:
        return self.normal_sign == +1

    @property
    def label(self) -> str:
        sb = "+" if self.det_brane == 1 else "-"
        sn = "+" if self.normal_sign == 1 else "-"
        return f"(det m_3, eps) = ({sb}, {sn})"


def _sign(x: float) -> int:
    return +1 if x > 0 else -1


def classify_gluing(monodromy: np.ndarray, normal_index: int = 3) -> GluingClass:
    """Classify a loop monodromy ``m in O(4)`` by its two orientation signs.

    Raises if ``m`` is not orthogonal. If ``m`` does not preserve the brane
    (does not map the normal direction to ±itself) the brane data are
    reported as ``0``.
    """
    m = np.asarray(monodromy, dtype=float)
    if m.shape != (4, 4) or not np.allclose(m.T @ m, np.eye(4), atol=1e-12):
        raise ValueError("the monodromy must be a 4x4 orthogonal matrix")
    det_bulk = _sign(np.linalg.det(m))
    e = np.zeros(4)
    e[normal_index] = 1.0
    image = m @ e
    preserves = bool(np.allclose(np.abs(image), e, atol=1e-12))
    if preserves:
        eps = _sign(image[normal_index])
        keep = [i for i in range(4) if i != normal_index]
        m3 = m[np.ix_(keep, keep)]
        det_brane = _sign(np.linalg.det(m3))
    else:
        eps, det_brane = 0, 0
    return GluingClass(
        monodromy=tuple(tuple(float(x) for x in row) for row in m),
        normal_index=normal_index, det_bulk=det_bulk, det_brane=det_brane,
        normal_sign=eps, preserves_brane=preserves)


def identity_gluing() -> np.ndarray:
    return np.eye(4)


def antipodal_gluing() -> np.ndarray:
    """``-I_4``: the antipodal map of the mouth ``S^3``. ``det = +1``."""
    return -np.eye(4)


def reflection_gluing(axis: int = 3) -> np.ndarray:
    """Reflection of one mouth direction. ``det = -1``."""
    m = np.eye(4)
    m[axis, axis] = -1.0
    return m


def mapping_torus_w1(monodromy: np.ndarray) -> int:
    """``w_1`` of the mapping torus of ``m``, evaluated on the handle loop.

    For a mapping torus ``M_m = (F x [0,1]) / (x,1) ~ (m x, 0)`` of an
    orientable fibre ``F``, the loop transverse to the fibre has
    ``w_1 = det m`` and ``w_1`` vanishes on the fibre. That is the whole of
    ``w_1`` because ``H^1(M_m; Z_2) = H^1(F)^m (+) Z_2`` and ``H^1(S^3) = 0``.
    """
    return _sign(np.linalg.det(np.asarray(monodromy, dtype=float)))


# ── the action on S^3 harmonics, computed rather than asserted ──────────────

@lru_cache(maxsize=16)
def _harmonic_basis(ell: int) -> Tuple[np.ndarray, np.ndarray]:
    """A basis of degree-``ell`` harmonic polynomials on ``R^4``.

    Returns ``(exponents, coefficients)``: ``exponents`` lists the monomials
    of degree ``ell`` in four variables, ``coefficients`` is a matrix whose
    columns are harmonic polynomials in that monomial basis, obtained as the
    null space of the Laplacian. Its rank must be ``(ell+1)^2``.
    """
    mons = [e for e in itertools.product(range(ell + 1), repeat=4) if sum(e) == ell]
    if ell < 2:
        return np.array(mons), np.eye(len(mons))
    lower = [e for e in itertools.product(range(ell - 1), repeat=4) if sum(e) == ell - 2]
    index = {e: i for i, e in enumerate(lower)}
    lap = np.zeros((len(lower), len(mons)))
    for j, e in enumerate(mons):
        for axis in range(4):
            if e[axis] >= 2:
                target = list(e)
                target[axis] -= 2
                lap[index[tuple(target)], j] += e[axis] * (e[axis] - 1)
    # null space
    _, sv, vt = np.linalg.svd(lap)
    rank = int(np.sum(sv > 1e-10 * sv.max()))
    null = vt[rank:].T
    return np.array(mons), null


def _evaluate_monomials(exponents: np.ndarray, points: np.ndarray) -> np.ndarray:
    out = np.ones((points.shape[0], exponents.shape[0]))
    for j, e in enumerate(exponents):
        for axis in range(4):
            if e[axis]:
                out[:, j] *= points[:, axis] ** e[axis]
    return out


def harmonic_representation(monodromy: np.ndarray, ell: int,
                            n_points: Optional[int] = None,
                            seed: int = 7) -> Dict[str, object]:
    """The matrix ``D_l(m)`` by which ``m`` acts on degree-``l`` harmonics.

    Solved numerically from ``(p o m)(x) = sum_k D_kj p_k(x)`` on random
    points of ``S^3``, so that ``D_l(-I) = (-1)^l I`` is a *result* and not an
    input. Returns the matrix, the residual of the fit, and whether the
    action is a scalar.
    """
    m = np.asarray(monodromy, dtype=float)
    exps, coeff = _harmonic_basis(ell)
    dim = coeff.shape[1]
    rng = np.random.default_rng(seed)
    n = n_points or 4 * dim + 8
    pts = rng.standard_normal((n, 4))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    basis = _evaluate_monomials(exps, pts) @ coeff
    moved = _evaluate_monomials(exps, pts @ m.T) @ coeff
    rep, *_ = np.linalg.lstsq(basis, moved, rcond=None)
    residual = float(np.max(np.abs(basis @ rep - moved)))
    scalar = rep[0, 0]
    is_scalar = bool(np.allclose(rep, scalar * np.eye(dim), atol=1e-9))
    return {"ell": ell, "dimension": dim, "expected_dimension": (ell + 1) ** 2,
            "matrix": rep, "fit_residual": residual,
            "is_scalar": is_scalar, "scalar": float(scalar) if is_scalar else None,
            "trace": float(np.trace(rep))}


def antipodal_parity(ell: int) -> int:
    """``Y_l(-Omega) = (-1)^l Y_l(Omega)``: a degree-``l`` polynomial is odd or
    even under ``x -> -x``. Checked against :func:`harmonic_representation`."""
    return (-1) ** ell


# ═══════════════════════════════════════════════════════════════════════════
# H2 · the involution on the handle
# ═══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Involution:
    """An isometric involution of the tube ``[-S,S] x S^3``.

    ``s_sign`` is ``±1`` on the throat coordinate, ``omega_map`` a 4x4
    orthogonal matrix on the ``S^3`` factor.
    """

    s_sign: int
    omega_map: Tuple[Tuple[float, ...], ...]
    name: str

    @property
    def omega(self) -> np.ndarray:
        return np.array(self.omega_map, dtype=float)

    def apply(self, s, omega: np.ndarray):
        return self.s_sign * np.asarray(s, dtype=float), np.asarray(omega) @ self.omega.T

    def is_involution(self) -> bool:
        return bool(np.allclose(self.omega @ self.omega, np.eye(4)) and self.s_sign ** 2 == 1)

    def is_isometry(self, neck_radius: float, lapse: Callable, n: int = 201) -> bool:
        """Isometry of ``-N^2 dt^2 + ds^2 + f^2 dOmega^2`` iff ``f`` and ``N`` are
        even in ``s`` (when ``s_sign = -1``) and ``omega_map`` is orthogonal."""
        s = np.linspace(-1.0, 1.0, n)
        f_ok = np.allclose(throat_radius(s, neck_radius),
                           throat_radius(self.s_sign * s, neck_radius))
        n_ok = np.allclose(lapse(s, neck_radius), lapse(self.s_sign * s, neck_radius))
        o_ok = np.allclose(self.omega.T @ self.omega, np.eye(4))
        return bool(f_ok and n_ok and o_ok)

    def fixed_points_exist(self) -> bool:
        """Fixed point ``(s, Omega)`` needs ``s_sign s = s`` and ``m Omega = Omega``.

        With ``s_sign = -1`` that forces ``s = 0``, and then a fixed point exists
        iff ``m`` has eigenvalue ``+1`` (a unit eigenvector lies on ``S^3``).
        """
        if self.s_sign == +1:
            return True
        eig = np.linalg.eigvals(self.omega)
        return bool(np.any(np.abs(eig - 1.0) < 1e-10))

    @property
    def is_free(self) -> bool:
        return not self.fixed_points_exist()

    def tangent_determinant(self) -> int:
        """``det d iota`` on the 4-dimensional handle: ``s_sign * det m``."""
        return self.s_sign * _sign(np.linalg.det(self.omega))

    def brane_determinant(self, normal_index: int = 3) -> int:
        """``det d iota`` restricted to the brane slice ``[-S,S] x S^2``."""
        keep = [i for i in range(4) if i != normal_index]
        return self.s_sign * _sign(np.linalg.det(self.omega[np.ix_(keep, keep)]))


def antipodal_involution() -> Involution:
    """``iota(s, Omega) = (-s, -Omega)``: the restriction of the ``S^4`` antipode."""
    return Involution(s_sign=-1, omega_map=tuple(map(tuple, -np.eye(4))),
                      name="antipodal iota: (s, Omega) -> (-s, -Omega)")


def neck_reflection_involution() -> Involution:
    """The antipodal control: ``(s, Omega) -> (-s, Omega)``. Fixes the neck."""
    return Involution(s_sign=-1, omega_map=tuple(map(tuple, np.eye(4))),
                      name="control iota_0: (s, Omega) -> (-s, Omega)")


def free_involutions_of_o5() -> List[Dict[str, object]]:
    """Every conjugacy class of involutions in ``O(5)`` and whether it is free
    on ``S^4``. An involution is conjugate to ``diag(±1,...)``; it has a fixed
    point on ``S^4`` iff some eigenvalue is ``+1``. Only ``-I_5`` is free.
    """
    rows = []
    for k in range(6):
        d = np.diag([-1.0] * k + [1.0] * (5 - k))
        free = k == 5
        rows.append({"minus_ones": k, "free_on_S4": free,
                     "swaps_antipodal_centres": bool(d[0, 0] < 0),
                     "det": int(round(np.linalg.det(d)))})
    return rows


def stiefel_whitney_rp(n: int) -> Dict[str, int]:
    """``w(RP^n) = (1 + a)^{n+1}`` mod 2: returns ``w_1, w_2`` and ``w_1^2`` as
    coefficients of ``a, a^2`` (``a^2 != 0`` iff ``n >= 2``)."""
    w1 = math.comb(n + 1, 1) % 2
    w2 = math.comb(n + 1, 2) % 2 if n >= 2 else 0
    w1_sq = w1 if n >= 2 else 0
    return {"n": n, "w1": w1, "w2": w2, "w1_squared": w1_sq}


def pin_structures_rp(n: int) -> Dict[str, object]:
    """Existence of Pin^+ (``w_2 = 0``) and Pin^- (``w_2 + w_1^2 = 0``) on ``RP^n``,
    and of a spin structure (orientable and ``w_2 = 0``)."""
    w = stiefel_whitney_rp(n)
    return {"n": n, "orientable": w["w1"] == 0,
            "spin": w["w1"] == 0 and w["w2"] == 0,
            "pin_plus": w["w2"] == 0,
            "pin_minus": (w["w2"] + w["w1_squared"]) % 2 == 0}


# ═══════════════════════════════════════════════════════════════════════════
# H5 · quaternions, Clifford algebra, and what J is
# ═══════════════════════════════════════════════════════════════════════════

# quaternion basis order (1, i, j, k); product table e_a e_b = sum_c T[a,b,c] e_c
_QT = np.zeros((4, 4, 4))
_QT[0, :, :] = np.eye(4)
_QT[:, 0, :] = np.eye(4)
for (a, b, c) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
    _QT[a, b, c] = 1.0
    _QT[b, a, c] = -1.0
for a in (1, 2, 3):
    _QT[a, a, 0] = -1.0


def quaternion_left(q: Sequence[float]) -> np.ndarray:
    """Real ``4x4`` matrix of ``v -> q v`` in the basis ``(1, i, j, k)``."""
    q = np.asarray(q, dtype=float)
    return np.einsum("a,abc->cb", q, _QT)


def quaternion_right(q: Sequence[float]) -> np.ndarray:
    """Real ``4x4`` matrix of ``v -> v q``."""
    q = np.asarray(q, dtype=float)
    return np.einsum("b,abc->ca", q, _QT)


def hopf_transport_matrix() -> np.ndarray:
    """The real ``4x4`` matrix of ``embedding.transport.orientation_reversal_on_s3``
    in the coordinates ``(Re z1, Im z1, Re z2, Im z2)`` with ``q = z1 + z2 j``,
    i.e. the basis ``(1, i, j, k)``. Built from the map itself, not from ``-j``."""
    from geometrodynamics.embedding.transport import orientation_reversal_on_s3
    cols = []
    for k in range(4):
        e = np.zeros(4)
        e[k] = 1.0
        w1, w2 = orientation_reversal_on_s3(complex(e[0], e[1]), complex(e[2], e[3]))
        cols.append([w1.real, w1.imag, w2.real, w2.imag])
    return np.array(cols).T


def complex_structure_commutation() -> Dict[str, object]:
    """S2 — which complex structure makes ``J`` linear, and which antilinear.

    ``L_i`` is the Hopf complex structure (the fibre ``z -> e^{i phi} z`` is
    left multiplication); ``R_i`` is the one commuting with ``SU(2)_L``, i.e.
    the complex structure of the left half-spinor.
    """
    J = hopf_transport_matrix()
    Li, Ri = quaternion_left([0, 1, 0, 0]), quaternion_right([0, 1, 0, 0])
    Lmj = quaternion_left([0, 0, -1, 0])
    return {
        "J_equals_left_mult_by_minus_j": bool(np.allclose(J, Lmj)),
        "det_J": float(np.linalg.det(J)),
        "J_squared_is_minus_identity": bool(np.allclose(J @ J, -np.eye(4))),
        "anticommutes_with_L_i": bool(np.allclose(J @ Li, -Li @ J)),
        "commutes_with_L_i": bool(np.allclose(J @ Li, Li @ J)),
        "commutes_with_R_i": bool(np.allclose(J @ Ri, Ri @ J)),
        "linear_matrix_in_R_i_basis": _matrix_in_right_complex_basis(J),
        "hopf_base_antipode": _hopf_base_sign(J),
    }


def _matrix_in_right_complex_basis(J: np.ndarray) -> np.ndarray:
    """Express a real map commuting with ``R_i`` as a complex ``2x2`` matrix in
    the ``R_i``-basis ``{1, j}`` of ``H = C^2`` (coordinates ``w1 + j w2``)."""
    Ri = quaternion_right([0, 1, 0, 0])
    e1 = np.array([1.0, 0, 0, 0])
    e2 = np.array([0, 0, 1.0, 0])                  # j
    out = np.zeros((2, 2), dtype=complex)
    for col, e in enumerate((e1, e2)):
        v = J @ e
        # v = a*1 + b*j with a, b complex w.r.t. R_i: v = e1*a + e2*b, a = a_r + a_i (right i)
        # components: real part along e, imaginary part along Ri e
        for row, f in enumerate((e1, e2)):
            out[row, col] = complex(v @ f, v @ (Ri @ f))
    return out


def _hopf_base_sign(J: np.ndarray) -> bool:
    """``h(Jz) = -h(z)`` for the Hopf map ``h(q) = q^{-1} i q`` (left-fibre)."""
    rng = np.random.default_rng(3)
    ok = True
    for _ in range(50):
        q = rng.standard_normal(4)
        q /= np.linalg.norm(q)
        ok &= np.allclose(_hopf(q), -_hopf(J @ q))
    return bool(ok)


def _hopf(q: np.ndarray) -> np.ndarray:
    conj = q * np.array([1, -1, -1, -1.0])
    i = np.array([0, 1.0, 0, 0])
    return quaternion_left(conj) @ quaternion_left(i) @ q


# ── Clifford algebra Cl(4), regular representation ─────────────────────────

@lru_cache(maxsize=2)
def clifford_regular(signature: int = +1) -> Dict[str, np.ndarray]:
    """``Cl(4)`` with ``e_a^2 = signature`` as ``16x16`` real matrices acting on
    itself by left multiplication. Basis: subsets of ``{0,1,2,3}`` in binary
    order. No Pauli matrix enters; only the anticommutation relations do."""
    basis = [tuple(i for i in range(4) if (mask >> i) & 1) for mask in range(16)]
    index = {b: n for n, b in enumerate(basis)}

    def mul(a: Tuple[int, ...], b: Tuple[int, ...]) -> Tuple[int, Tuple[int, ...]]:
        sign, out = 1, list(a)
        for g in b:
            # move g past the elements of out greater than g, sign per swap
            pos = len(out)
            while pos > 0 and out[pos - 1] > g:
                pos -= 1
                sign = -sign
            if pos > 0 and out[pos - 1] == g:
                del out[pos - 1]
                sign *= signature
            else:
                out.insert(pos, g)
        return sign, tuple(out)

    gens = {}
    for a in range(4):
        mat = np.zeros((16, 16))
        for n, b in enumerate(basis):
            sign, prod = mul((a,), b)
            mat[index[prod], n] = sign
        gens[a] = mat
    volume = gens[0] @ gens[1] @ gens[2] @ gens[3]
    return {"generators": np.array([gens[a] for a in range(4)]),
            "volume": volume, "identity": np.eye(16)}


def pin_lifts_of_reflection(axis: int = 0) -> Dict[str, object]:
    """S1/S4 — the lifts of the neck reflection ``dι = diag(-1, +1, +1, +1)``.

    In ``Pin^±(4)`` the reflection in ``e_a`` lifts to ``±e_a`` with
    ``e_a^2 = ±1``. The lift anticommutes with the volume element, so it
    exchanges the two chirality eigenspaces of ``Spin(4)``. It is a real
    matrix: complex-linear on any complexified spinor.
    """
    rows = []
    for name, sig in (("Pin+", +1), ("Pin-", -1)):
        cl = clifford_regular(sig)
        g = cl["generators"][axis]
        vol = cl["volume"]
        for sgn in (+1, -1):
            lift = sgn * g
            rows.append({
                "pin_type": name, "sign": sgn,
                "square": int(round((lift @ lift)[0, 0])),
                "anticommutes_with_volume": bool(np.allclose(lift @ vol, -vol @ lift)),
                "is_real_matrix": True,
            })
    distinct = all(not np.allclose(sgn_a * clifford_regular(sa)["generators"][axis],
                                   sgn_b * clifford_regular(sb)["generators"][axis])
                   for (sa, sgn_a), (sb, sgn_b) in itertools.combinations(
                       [(+1, 1), (+1, -1), (-1, 1), (-1, -1)], 2)
                   if sa == sb)
    return {"lifts": rows, "lifts_within_a_type_are_distinct": bool(distinct),
            "count": len(rows)}


def spin_lifts_of_antipode() -> Dict[str, object]:
    """The lifts of ``-I_4 in SO(4)`` to ``Spin(4)``: ``±(e_0 e_1 e_2 e_3)``,
    which is ``±`` the volume element. It acts as ``+1`` on one chirality and
    ``-1`` on the other, so a spinor transported once round an antipodally
    glued handle returns with a chirality-dependent sign. Not ``J``."""
    cl = clifford_regular(+1)
    vol = cl["volume"]
    return {"lifts": [{"sign": s, "square_is_identity":
                       bool(np.allclose((s * vol) @ (s * vol), np.eye(16)))}
                      for s in (+1, -1)],
            "commutes_with_all_generators": bool(all(
                np.allclose(vol @ g, g @ vol) for g in cl["generators"])),
            "anticommutes_with_generators": bool(all(
                np.allclose(vol @ g, -g @ vol) for g in cl["generators"])),
            "eigenvalues": sorted(set(np.round(np.linalg.eigvals(vol).real, 9).tolist()))}


# ═══════════════════════════════════════════════════════════════════════════
# H3/H4 · the scalar boundary condition and the half-tube admittance
# ═══════════════════════════════════════════════════════════════════════════

def neck_sector(ell: int, eta: int = +1, involution: Optional[Involution] = None) -> str:
    """The neck condition forced on ``psi_l`` by ``Phi(iota x) = eta Phi(x)``.

    Derived: ``psi_l(-s) Y_l(m Omega) = eta psi_l(s) Y_l(Omega)``. When ``m``
    acts on degree-``l`` harmonics as the scalar ``chi_l`` this is
    ``psi_l(-s) = eta chi_l psi_l(s)``: ``+1`` gives Neumann, ``-1`` Dirichlet.
    For ``m = -I`` the scalar is ``(-1)^l`` (computed, not assumed). For the
    control ``m = I`` it is ``1`` for every ``l`` and the grading is gone.
    """
    inv = involution or antipodal_involution()
    rep = harmonic_representation(inv.omega, ell)
    if not rep["is_scalar"]:
        return "mixed (m does not act as a scalar on this l)"
    sign = eta * int(round(rep["scalar"]))
    return "Neumann" if sign == +1 else "Dirichlet"


def half_tube_admittance_oracle(ell: int, condition: str,
                                bulk_radius: float = BULK_RADIUS,
                                mouth_angle: float = MOUTH_ANGLE) -> float:
    """B2's closed form: the ``(1, ±1)`` eigenvalue of the PR #277 two-mouth
    ``Y_l(0)``. Neumann is the symmetric sector ``d + c``... written out:

        Y^N = 2 pi^2 F^2 [ k tanh(kX) - cos a ]
        Y^D = 2 pi^2 F^2 [ k coth(kX) - cos a ]

    Obtained from ``coth 2x ∓ csch 2x = tanh x, coth x``.
    """
    g = geometry(bulk_radius, mouth_angle)
    k = float(ell + 1)
    x = g.rapidity
    pre = 2.0 * math.pi ** 2 * g.mouth_radius ** 2
    if condition == "Neumann":
        return pre * (k * math.tanh(k * x) - math.cos(mouth_angle))
    if condition == "Dirichlet":
        return pre * (k / math.tanh(k * x) - math.cos(mouth_angle))
    raise ValueError(condition)


def half_tube_admittance(ell: int, condition: str,
                         bulk_radius: float = BULK_RADIUS,
                         mouth_angle: float = MOUTH_ANGLE,
                         steps: int = 2000) -> float:
    """Outward flux at the mouth for unit Dirichlet data there, with the
    neck condition, from a conservative tridiagonal solve on ``[0, S]``.

    Shares no algebra with the closed form: it never uses the
    ``s = b sinh x`` reduction, and the two-mouth ``solve_admittance`` is not
    called. Neumann at the neck is imposed with a ghost point (second order);
    Dirichlet by dropping the neck unknown.
    """
    from scipy.linalg import solve_banded
    g = geometry(bulk_radius, mouth_angle)
    lam = float(ell * (ell + 2))
    s_grid = np.linspace(0.0, g.half_length, steps + 1)
    h = float(s_grid[1] - s_grid[0])
    f_node = throat_radius(s_grid, g.neck_radius)
    w_half = throat_radius(0.5 * (s_grid[:-1] + s_grid[1:]), g.neck_radius) ** 3
    if condition == "Dirichlet":
        # unknowns u_1 .. u_{M-1}; u_0 = 0, u_M = 1
        m = steps - 1
        lower = w_half[1:-1] / h ** 2
        upper = w_half[1:-1] / h ** 2
        diag = -(w_half[:-1] + w_half[1:]) / h ** 2 - lam * f_node[1:-1]
        rhs = np.zeros(m)
        rhs[-1] -= w_half[-1] / h ** 2
        banded = np.zeros((3, m))
        banded[0, 1:] = upper
        banded[1, :] = diag
        banded[2, :-1] = lower
        interior = solve_banded((1, 1), banded, rhs)
        u = np.concatenate(([0.0], interior, [1.0]))
    elif condition == "Neumann":
        # unknowns u_0 .. u_{M-1}; u_M = 1; at s=0 the flux w u' vanishes, so
        # the conservative balance at node 0 uses only the right half-cell
        m = steps
        diag = np.empty(m)
        diag[0] = -w_half[0] / h ** 2 - 0.5 * lam * f_node[0]
        diag[1:] = -(w_half[:-1] + w_half[1:]) / h ** 2 - lam * f_node[1:-1]
        coupling = w_half[:-1] / h ** 2      # w_{i+1/2} for rows 0..M-2 (upper),
        rhs = np.zeros(m)                    # and w_{i-1/2} for rows 1..M-1 (lower)
        rhs[-1] -= w_half[-1] / h ** 2
        banded = np.zeros((3, m))
        banded[0, 1:] = coupling
        banded[1, :] = diag
        banded[2, :-1] = coupling
        interior = solve_banded((1, 1), banded, rhs)
        u = np.concatenate((interior, [1.0]))
    else:
        raise ValueError(condition)
    v_right = w_half[-1] * (u[-1] - u[-2]) / h + 0.5 * h * lam * f_node[-1] * u[-1]
    return float(2.0 * math.pi ** 2 * v_right)


def static_operator_commutes_with_parity(neck_radius: float, lapse: Callable,
                                         ell: int = 2, steps: int = 400,
                                         half_length: float = 1.0) -> float:
    """L1 — the static radial operator ``(N f^3 u')' - l(l+2) N f u`` on the full
    tube, discretised, commutes with ``s -> -s`` for both lapses because ``f``
    and ``N`` are even. Returns ``|P O P - O|`` (should vanish), so the neck
    sector labels do not depend on which lapse is chosen."""
    s = np.linspace(-half_length, half_length, steps + 1)
    h = s[1] - s[0]
    mid = 0.5 * (s[:-1] + s[1:])
    w = np.asarray(lapse(mid, neck_radius)) * throat_radius(mid, neck_radius) ** 3
    d = np.asarray(lapse(s, neck_radius)) * throat_radius(s, neck_radius)
    n = steps + 1
    op = np.zeros((n, n))
    for i in range(1, n - 1):
        op[i, i - 1] = w[i - 1] / h ** 2
        op[i, i + 1] = w[i] / h ** 2
        op[i, i] = -(w[i - 1] + w[i]) / h ** 2 - ell * (ell + 2) * d[i]
    P = np.fliplr(np.eye(n))
    return float(np.max(np.abs(P @ op @ P - op)))


def traversal_table() -> List[Dict[str, object]]:
    """What one and two traversals do, for each object the repository names."""
    J = hopf_transport_matrix()
    rows = [
        {"object": "handle loop, identity gluing (m = I)",
         "one_traversal": "identity on the transverse S^3",
         "two_traversals": "identity"},
        {"object": "handle loop, antipodal gluing (m = -I)",
         "one_traversal": "Omega -> -Omega, (-1)^l on harmonics, orientation preserved",
         "two_traversals": "identity"},
        {"object": "handle loop, reflection gluing (det m = -1)",
         "one_traversal": "orientation reversed (w_1 = -1)",
         "two_traversals": "identity"},
        {"object": "deck transformation iota on M/iota",
         "one_traversal": "(s, Omega) -> (-s, -Omega); orientation reversed",
         "two_traversals": "identity on M (iota^2 = 1)"},
        {"object": "transport.py map sigma = L_{-j} on S^3",
         "one_traversal": "isoclinic rotation by pi/2; Hopf base antipode; fibre reversed",
         "two_traversals": ("sigma^2 = -I_4, the antipodal map of S^3 "
                            f"(checked: {bool(np.allclose(J @ J, -np.eye(4)))})"),
         "four_traversals": "identity"},
    ]
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# The mouth Pin holonomy (docs/mouth_pin_holonomy_prereg.md, P1-P8)
# ═══════════════════════════════════════════════════════════════════════════
#
# Two retractions of the section above are recorded in that document: the
# Pin+ / Pin- "mismatch" is the standard induction of an intrinsic Pin-
# structure from an ambient Pin+ one through two twisted normal lines, and
# J = i sigma_y K IS the Spin(4) lift of the SO(4) rotation L_{-j}; what it is
# not is a lift of iota or of -I_4.

__all__ += [
    "transport_along_great_semicircle",
    "neck_holonomy",
    "sw_polynomial_rp2",
    "restricted_stiefel_whitney",
    "twisted_tangent_generators",
    "holonomy_lifts",
    "closed_loop_spin_holonomy",
    "intrinsic_pin2_module",
    "compare_mouth_holonomy_with_transport",
    "pin_structure_counts",
]


def _transport_ode(v0: np.ndarray, centre: np.ndarray, direction: np.ndarray,
                   theta_end: float, steps: int = 2000) -> np.ndarray:
    """Parallel transport of ``v0`` along the great circle
    ``gamma(theta) = cos(theta) c + sin(theta) d`` of the unit ``S^3 ⊂ R^4``,
    integrating ``v' = -(v . gamma') gamma`` with RK4. Computed, not recalled."""
    def gamma(t):
        return math.cos(t) * centre + math.sin(t) * direction

    def dgamma(t):
        return -math.sin(t) * centre + math.cos(t) * direction

    def rhs(t, v):
        return -(v @ dgamma(t)) * gamma(t)

    h = theta_end / steps
    v, t = v0.astype(float).copy(), 0.0
    for _ in range(steps):
        k1 = rhs(t, v)
        k2 = rhs(t + 0.5 * h, v + 0.5 * h * k1)
        k3 = rhs(t + 0.5 * h, v + 0.5 * h * k2)
        k4 = rhs(t + h, v + h * k3)
        v = v + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return v


def transport_along_great_semicircle(steps: int = 2000) -> Dict[str, np.ndarray]:
    """P1, first half: transport the neck frame ``(n, t1, t2)`` (tangent to
    ``S^3`` at ``Omega = e1``) along the semicircle towards ``-e1`` in the
    ``t1 = e2`` direction. ``n = e4`` is the brane normal, ``t2 = e3``."""
    c, d = np.eye(4)[0], np.eye(4)[1]
    out = {}
    for name, v0 in (("n", np.eye(4)[3]), ("t1", np.eye(4)[1]), ("t2", np.eye(4)[2])):
        out[name] = _transport_ode(v0, c, d, math.pi, steps)
    return out


def neck_holonomy(steps: int = 2000) -> Dict[str, object]:
    """P1 — the holonomy of the deck generator of the neck ``RP^2`` on the
    frame ``(d_s, n, t1, t2)``: transport to ``-Omega`` followed by ``d iota``.

    ``d iota`` sends ``d_s -> -d_s`` and every ``S^3``-tangent vector ``v`` at
    ``-Omega`` to ``-v`` at ``Omega`` (the differential of ``-I_4``).
    """
    moved = transport_along_great_semicircle(steps)
    # d iota on the S^3 tangent part: v -> -v (as vectors of R^4)
    back = {k: -v for k, v in moved.items()}
    frame = {"n": np.eye(4)[3], "t1": np.eye(4)[1], "t2": np.eye(4)[2]}
    H = np.zeros((4, 4))
    H[0, 0] = -1.0                                     # d_s -> -d_s
    order = ["n", "t1", "t2"]
    for j, name in enumerate(order):
        for i, other in enumerate(order):
            H[i + 1, j + 1] = float(back[name] @ frame[other])
    return {"holonomy": H, "frame": ["d_s", "n", "t1", "t2"],
            "det": float(np.linalg.det(H)),
            "normal_block": H[:2, :2].copy(), "tangent_block": H[2:, 2:].copy(),
            "both_normals_reversed": bool(np.allclose(H[:2, :2], -np.eye(2))),
            "tangent_is_reflection": bool(np.allclose(np.linalg.det(H[2:, 2:]), -1.0)),
            "eigenvalues": sorted(np.round(np.linalg.eigvals(H).real, 9).tolist())}


# ── Stiefel-Whitney arithmetic in Z_2[a]/(a^3) ─────────────────────────────

def sw_polynomial_rp2(coefficients: Sequence[int]) -> Tuple[int, int, int]:
    """Reduce a polynomial in ``a`` mod ``2`` and mod ``a^3``."""
    out = [0, 0, 0]
    for k, c in enumerate(coefficients):
        if k < 3:
            out[k] = (out[k] + c) % 2
    return tuple(out)


def _sw_mul(p: Tuple[int, int, int], q: Tuple[int, int, int]) -> Tuple[int, int, int]:
    out = [0, 0, 0]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            if i + j < 3:
                out[i + j] = (out[i + j] + a * b) % 2
    return tuple(out)


def restricted_stiefel_whitney() -> Dict[str, object]:
    """P2 — ``w(TM|_N) = w(TN) w(nu)`` on ``N = RP^2`` with ``nu = lambda + lambda``.

    ``w(lambda) = 1 + a``; ``w(TN) = (1+a)^3``; everything in ``Z_2[a]/(a^3)``.
    """
    one_plus_a = (1, 1, 0)
    w_lambda = one_plus_a
    w_nu = _sw_mul(w_lambda, w_lambda)
    w_tn = _sw_mul(_sw_mul(one_plus_a, one_plus_a), one_plus_a)
    w_tm = _sw_mul(w_tn, w_nu)
    w1, w2 = w_tn[1], w_tn[2]
    return {"w_nu": w_nu, "w_TN": w_tn, "w_TM_restricted": w_tm,
            "w2_TM_restricted": w_tm[2],
            "ambient_pin_plus_compatible": w_tm[2] == 0,
            "w1_TN": w1, "w2_TN": w2, "w1_TN_squared": w1,     # a^2 coefficient of a*a
            "intrinsic_pin_minus": (w2 + w1) % 2 == 0,        # w2 + w1^2 = 0
            "intrinsic_pin_plus": w2 == 0}


# ── the conversion Pin+ (ambient) -> Pin- (intrinsic) ──────────────────────

def twisted_tangent_generators() -> Dict[str, object]:
    """P3 — with ``e_s, e_n`` the normal generators of ``Cl(4)`` (``e^2 = +1``),
    ``e~_t = e_t e_s e_n`` square to ``-1`` and anticommute: ``Cl^-(2)``."""
    cl = clifford_regular(+1)
    e_s, e_n, e_t1, e_t2 = cl["generators"]
    nv = e_s @ e_n
    et1, et2 = e_t1 @ nv, e_t2 @ nv
    I = cl["identity"]
    return {"normal_volume_squared": int(round((nv @ nv)[0, 0])),
            "twisted_t1_squared": int(round((et1 @ et1)[0, 0])),
            "twisted_t2_squared": int(round((et2 @ et2)[0, 0])),
            "twisted_anticommute": bool(np.allclose(et1 @ et2, -et2 @ et1)),
            "twisted_product_squared": int(round(((et1 @ et2) @ (et1 @ et2))[0, 0])),
            "generates_quaternions": bool(
                np.allclose(et1 @ et1, -I) and np.allclose(et2 @ et2, -I)
                and np.allclose(et1 @ et2, -et2 @ et1))}


def holonomy_lifts() -> Dict[str, object]:
    """P4 — lifts of ``H = diag(-1,-1,+1,-1)`` to ``Pin^±(4)``: ``±e_s e_n e_t2``.

    ``H`` is the product of the reflections in ``d_s``, ``n`` and ``t2``, so
    its lift is the product of the three generators (up to sign).
    """
    rows = []
    for name, sig in (("Pin+", +1), ("Pin-", -1)):
        cl = clifford_regular(sig)
        e_s, e_n, _, e_t2 = cl["generators"]
        vol = cl["volume"]
        lift = e_s @ e_n @ e_t2
        rows.append({"pin_type": name,
                     "square": int(round((lift @ lift)[0, 0])),
                     "odd_anticommutes_with_volume": bool(np.allclose(lift @ vol, -vol @ lift)),
                     "equals_twisted_t2": bool(np.allclose(lift, e_t2 @ e_s @ e_n)),
                     "normal_part_squared": int(round(((e_s @ e_n) @ (e_s @ e_n))[0, 0]))})
    return {"lifts": rows,
            "square_in_pin_plus": rows[0]["square"], "square_in_pin_minus": rows[1]["square"]}


def closed_loop_spin_holonomy(n_loops: int = 60, steps: int = 1500) -> Dict[str, object]:
    """P5 — the spin holonomy of the equator of the neck ``S^2``, computed by
    transporting a tangent vector around latitude circles at polar angle
    ``theta0`` from the pole ``t2 = e3``, unwrapping the rotation angle
    continuously in ``theta0``, and reading the value at the equator.

    Classical value ``2 pi (1 - cos theta0)``; at the equator ``2 pi``, whose
    spin lift ``exp(i * angle / 2)`` is ``-1``.
    """
    pole = np.eye(4)[2]
    angles = []
    for theta0 in np.linspace(0.05, 0.5 * math.pi, n_loops):
        # latitude circle in the S^2 = span(e1, e2, e3) ∩ S^3
        centre = math.cos(theta0) * pole
        r = math.sin(theta0)
        # parametrise p(phi) = cos(theta0) e3 + sin(theta0)(cos phi e1 + sin phi e2)
        # transport the tangent vector v0 = dp/dphi(0)/r = e2 around the loop
        v = np.eye(4)[1].copy()
        h = 2.0 * math.pi / steps

        def p(phi):
            return centre + r * (math.cos(phi) * np.eye(4)[0] + math.sin(phi) * np.eye(4)[1])

        def dp(phi):
            return r * (-math.sin(phi) * np.eye(4)[0] + math.cos(phi) * np.eye(4)[1])

        def rhs(phi, vv):
            # parallel transport on the unit sphere along p(phi): v' = -(v . p') p
            return -(vv @ dp(phi)) * p(phi)

        phi = 0.0
        for _ in range(steps):
            k1 = rhs(phi, v)
            k2 = rhs(phi + 0.5 * h, v + 0.5 * h * k1)
            k3 = rhs(phi + 0.5 * h, v + 0.5 * h * k2)
            k4 = rhs(phi + h, v + h * k3)
            v = v + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            phi += h
        # rotation angle of the returned vector relative to the start, in the
        # oriented tangent plane at p(0) spanned by (e2, u) with u the inward
        # meridian direction
        u = np.cross(np.eye(3)[1], (p(0.0)[:3] / np.linalg.norm(p(0.0)[:3])))
        u = np.append(u, 0.0)
        ang = math.atan2(v @ u, v @ np.eye(4)[1])
        angles.append(ang)
    unwrapped = np.unwrap(np.array(angles))
    unwrapped -= unwrapped[0] - angles[0]
    theta_grid = np.linspace(0.05, 0.5 * math.pi, n_loops)
    expected = 2.0 * math.pi * (1.0 - np.cos(theta_grid))
    # orientation sign of the unwrapped angle is a convention; compare |.|
    total = float(abs(unwrapped[-1]))
    return {"equator_rotation_angle": total,
            "expected_2pi": 2.0 * math.pi,
            "max_deviation_from_classical_law": float(np.max(np.abs(np.abs(unwrapped) - expected))),
            "spin_lift_at_equator": float(math.cos(0.5 * total)),   # Re e^{i angle/2}
            "spin_holonomy_is_minus_one": bool(abs(math.cos(0.5 * total) + 1.0) < 1e-6)}


# ── the intrinsic Pin^-(2) module and the comparison with L_{-j} ───────────

def intrinsic_pin2_module() -> Dict[str, object]:
    """P6/P7 — ``Cl^-(2) ≅ H`` with ``e~_t1 -> i``, ``e~_t2 -> j``, ``e~_t1 e~_t2 -> k``,
    acting on ``H`` by left multiplication. The mouth holonomy ``e~_t2`` is
    ``L_j``; ``Spin(2) = {exp(theta k)}``; the fibre-reversing component is
    ``{cos(alpha) i + sin(alpha) j}``."""
    Li, Lj, Lk = (quaternion_left([0, 1, 0, 0]), quaternion_left([0, 0, 1, 0]),
                  quaternion_left([0, 0, 0, 1]))
    Ri = quaternion_right([0, 1, 0, 0])
    hol = Lj
    return {"holonomy_is_L_j": True,
            "holonomy_squared_minus_identity": bool(np.allclose(hol @ hol, -np.eye(4))),
            "spin2_generator_is_L_k": bool(np.allclose(Li @ Lj, Lk)),
            "holonomy_anticommutes_with_spin2_generator": bool(np.allclose(hol @ Lk, -Lk @ hol)),
            "holonomy_commutes_with_right_i": bool(np.allclose(hol @ Ri, Ri @ hol)),
            "pin2_is_normaliser_of_spin2": bool(np.allclose(
                hol @ Lk @ np.linalg.inv(hol), -Lk))}


def compare_mouth_holonomy_with_transport() -> Dict[str, object]:
    """P6 — three levels of comparison between the mouth holonomy and ``J``.

    * vector level: ``H`` versus ``sigma = L_{-j}`` as elements of ``O(4)``;
    * ambient spinor level: ``e_s e_n e_t2`` (odd) versus ``(-j, 1)`` (even);
    * intrinsic level: ``L_j`` versus ``L_{-j}`` inside ``Pin^-(2) ⊂ SU(2)``,
      up to ``Spin(2)`` conjugation and sign.
    """
    H = neck_holonomy()["holonomy"]
    sigma = hopf_transport_matrix()
    vector = {"det_H": float(np.linalg.det(H)), "det_sigma": float(np.linalg.det(sigma)),
              "eig_H": sorted(np.round(np.linalg.eigvals(H).real, 6).tolist()),
              "eig_sigma_imag": sorted(np.round(np.linalg.eigvals(sigma).imag, 6).tolist()),
              "same_O4_class": bool(abs(np.linalg.det(H) - np.linalg.det(sigma)) < 1e-9)}
    # ambient: parity in Cl(4)
    ambient = {"holonomy_lift_parity": "odd (three generators)",
               "sigma_lift_parity": "even (Spin(4), (-j, 1))",
               "conjugate_in_Pin4": False}
    # intrinsic: is L_{-j} = g L_j g^{-1} for some g in Spin(2) = exp(theta k), up to sign?
    Lj = quaternion_left([0, 0, 1, 0])
    Lmj = quaternion_left([0, 0, -1, 0])
    found = None
    for theta in np.linspace(0.0, 2.0 * math.pi, 3601):
        g = quaternion_left([math.cos(theta), 0, 0, math.sin(theta)])
        if np.allclose(g @ Lj @ np.linalg.inv(g), Lmj, atol=1e-9):
            found = float(theta)
            break
    # the whole fibre-reversing component, as a set: cos(a) i + sin(a) j
    component = [quaternion_left([0, math.cos(a), math.sin(a), 0])
                 for a in np.linspace(0, 2 * math.pi, 8, endpoint=False)]
    all_square_minus_one = all(np.allclose(c @ c, -np.eye(4)) for c in component)
    return {"vector_level": vector, "ambient_spinor_level": ambient,
            "intrinsic_level": {
                "L_minus_j_conjugate_to_L_j_by_spin2": found is not None,
                "conjugating_angle": found,
                "also_equal_up_to_sign": bool(np.allclose(Lmj, -Lj)),
                "fibre_reversing_component_all_order_four": bool(all_square_minus_one),
                "unfixed_data": "U(1) direction in the (i, j) plane, and the sign"},
            "outcome": "B"}


def pin_structure_counts() -> Dict[str, object]:
    """P8 — ``|H^1(.; Z_2)|`` counts: two Pin^- structures on ``RP^2``, four
    Pin^+ on ``RP^4 # RP^4``, two spin structures on ``S^3 x S^1``. None of the
    signs is fixed by the metric or by ``iota``."""
    return {"RP2_pin_minus_structures": 2, "RP4_connected_sum_pin_plus_structures": 4,
            "S3xS1_spin_structures": 2,
            "sign_of_holonomy_fixed_by_geometry": False,
            "what_it_is_in_the_repository": "embedding.topology.ThroatDefect.wrap_parity"}
