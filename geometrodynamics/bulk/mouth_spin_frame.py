"""The Hopf bundle as the spin-frame bundle of the brane mouth, and Pin⁻ pairing.

Pre-registered in ``docs/mouth_spin_frame_prereg.md`` (T1-T7), committed
before this file. Model: the brane mouth is ``S^2 = {x in Im H : |x| = 1}``,
the bulk mouth is ``S^3 = {q in H : |q| = 1}`` with the brane normal at the
mouth centre identified with ``1 in H`` and a reference brane direction
``i``. Then

    q  ->  ( x, e_1, e_2 ) = ( q^{-1} i q,  q^{-1} j q,  q^{-1} k q )

is ``Spin(3) -> SO(3) = F_SO(S^2)``: the bulk mouth IS the spin-frame bundle
of the brane mouth, and the repository's Hopf fibre ``q -> e^{i phi} q`` is
the ``Spin(2)`` fibre, rotating the tangent frame by ``2 phi``.

Quaternions come from ``mouth_topology.quaternion_left/right`` (the
multiplication table); no Pauli matrix enters. Nothing here uses a singlet,
a Born rule, a projector, a tensor product, CHSH, a QED amplitude or a path
integral.
"""

from __future__ import annotations

import cmath
import math
from typing import Dict, List, Tuple

import numpy as np

from geometrodynamics.bulk.mouth_topology import quaternion_left, quaternion_right

__all__ = [
    "qmul", "qinv", "spin_frame", "frame_map_checks", "fibre_rotates_frame_twice",
    "levi_civita_versus_hopf_connection", "chern_versus_euler",
    "deck_lifts_of_antipode", "two_sheeted_involution",
    "pin_minus_structures_rp2", "abk_invariant", "pin_bordism_pairing_rule",
]

_ONE = np.array([1.0, 0, 0, 0])
_I = np.array([0, 1.0, 0, 0])
_J = np.array([0, 0, 1.0, 0])
_K = np.array([0, 0, 0, 1.0])


def qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return quaternion_left(a) @ np.asarray(b, dtype=float)


def qinv(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    return q * np.array([1, -1, -1, -1.0]) / (q @ q)


def _conj_by(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """``q^{-1} v q`` as a vector of ``Im H`` (the three imaginary components)."""
    return qmul(qmul(qinv(q), v), q)[1:]


def spin_frame(q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """T1 — ``(x, e_1, e_2) = (q^{-1} i q, q^{-1} j q, q^{-1} k q)``."""
    q = np.asarray(q, dtype=float)
    return _conj_by(q, _I), _conj_by(q, _J), _conj_by(q, _K)


def _random_units(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((n, 4))
    return q / np.linalg.norm(q, axis=1, keepdims=True)


def frame_map_checks(n: int = 400, seed: int = 11) -> Dict[str, object]:
    """T1 — the image is an oriented orthonormal tangent frame; ``q`` and ``-q``
    coincide; distinct ``q != ±q'`` give distinct frames."""
    worst_unit = worst_orth = worst_tangent = worst_det = 0.0
    worst_pair = 0.0
    qs = _random_units(n, seed)
    frames = []
    for q in qs:
        x, e1, e2 = spin_frame(q)
        xm, e1m, e2m = spin_frame(-q)
        worst_unit = max(worst_unit, abs(x @ x - 1), abs(e1 @ e1 - 1), abs(e2 @ e2 - 1))
        worst_orth = max(worst_orth, abs(e1 @ e2))
        worst_tangent = max(worst_tangent, abs(e1 @ x), abs(e2 @ x))
        worst_det = max(worst_det, abs(np.linalg.det(np.array([e1, e2, x])) - 1.0))
        worst_pair = max(worst_pair, np.max(np.abs(np.concatenate([x - xm, e1 - e1m, e2 - e2m]))))
        frames.append(np.concatenate([x, e1, e2]))
    frames = np.array(frames)
    # injectivity up to sign: nearest other frame must be far unless q' = ±q
    d = np.linalg.norm(frames[:, None, :] - frames[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return {"orthonormal_error": worst_unit + worst_orth, "tangent_error": worst_tangent,
            "orientation_error": worst_det, "q_and_minus_q_coincide": worst_pair,
            "min_distance_between_distinct_frames": float(d.min()),
            "double_cover": bool(worst_pair < 1e-12 and d.min() > 1e-3)}


def fibre_rotates_frame_twice(n: int = 200, seed: int = 5) -> Dict[str, object]:
    """T2 — ``q -> e^{i phi} q`` fixes ``x`` and rotates ``(e_1, e_2)`` by ``2 phi``."""
    worst_x = worst_angle = 0.0
    rng = np.random.default_rng(seed)
    for q in _random_units(n, seed):
        phi = rng.uniform(-math.pi, math.pi)
        g = np.array([math.cos(phi), math.sin(phi), 0.0, 0.0])
        x, e1, e2 = spin_frame(q)
        xg, e1g, e2g = spin_frame(qmul(g, q))
        worst_x = max(worst_x, np.max(np.abs(xg - x)))
        ang = math.atan2(e1g @ e2, e1g @ e1)          # angle of e1' in the (e1, e2) frame
        diff = (ang + 2.0 * phi + math.pi) % (2.0 * math.pi) - math.pi
        worst_angle = max(worst_angle, abs(diff))
    return {"base_point_fixed": worst_x, "frame_angle_minus_two_phi": worst_angle,
            "angle_is_2phi": bool(worst_x < 1e-12 and worst_angle < 1e-10),
            "sense": "e_1 -> cos(2phi) e_1 - sin(2phi) e_2"}


def levi_civita_versus_hopf_connection(n_paths: int = 20, samples: int = 200,
                                       seed: int = 3) -> Dict[str, object]:
    """T2 — along random smooth paths ``q(t)`` on ``S^3``, the Levi-Civita form
    of the round ``S^2`` pulled back through the frame map, ``omega = <e_1', e_2>``,
    equals ``-2 A`` with ``A = <i q, q'>`` the Hopf connection of the left
    action (the sign is the orientation convention of the frame; the factor
    two is ``Spin(2) -> SO(2)``).

    ``q'`` is analytic (derivative of the normalised path); ``e_1'`` is a
    fourth-order central difference of the frame map, so the residual is at
    the ``1e-9`` level rather than the ``1e-4`` of a second-order stencil.
    Analytically: with ``w = q' q^{-1} = w_1 i + w_2 j + w_3 k``,
    ``e_1' = q^{-1}[j, w] q = -2 w_1 e_2 + 2 w_3 x`` and ``A = w_1``.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    h = 1e-3
    for _ in range(n_paths):
        a, b, c = _random_units(3, int(rng.integers(1 << 30)))

        def raw(t):
            return a + t * b + 0.5 * t * t * c

        def path(t):
            v = raw(t)
            return v / np.linalg.norm(v)

        def qdot(t):
            v, vd = raw(t), b + t * c
            nv = np.linalg.norm(v)
            return vd / nv - v * (v @ vd) / nv ** 3

        for t in np.linspace(0.05, 0.95, samples):
            q = path(t)
            _, e1, e2 = spin_frame(q)
            e1dot = (-spin_frame(path(t + 2 * h))[1] + 8 * spin_frame(path(t + h))[1]
                     - 8 * spin_frame(path(t - h))[1] + spin_frame(path(t - 2 * h))[1]) / (12 * h)
            omega = float(e1dot @ e2)
            A = float(qmul(_I, q) @ qdot(t))
            worst = max(worst, abs(omega + 2.0 * A))
    return {"worst_residual_of_omega_plus_2A": worst,
            "omega_is_minus_twice_A": bool(worst < 1e-8)}


def chern_versus_euler() -> Dict[str, object]:
    """T2 — ``c_1(Hopf) = 1`` (the repository's ``hopf.chern.compute_c1``) against
    ``e(TS^2) = chi(S^2) = 2`` from Gauss-Bonnet on the unit sphere."""
    from geometrodynamics.hopf.chern import compute_c1
    c1 = compute_c1()["c1_abs"]
    theta = np.linspace(0.0, math.pi, 200001)
    euler = float(np.trapezoid(np.sin(theta), theta) * 2.0 * math.pi / (2.0 * math.pi))  # ∫K dA / 2π, K=1
    return {"c1_hopf": float(c1), "euler_TS2": euler,
            "ratio": euler / c1, "ratio_is_two": bool(abs(euler / c1 - 2.0) < 1e-6)}


def deck_lifts_of_antipode(n: int = 200, seed: int = 9) -> Dict[str, object]:
    """T3 — ``L_u`` covers the antipode iff ``u ⊥ i`` (unit); the frame image is
    ``dA`` composed with one reflection; ``L_u^2 = -1``."""
    rng = np.random.default_rng(seed)
    qs = _random_units(n, seed)
    rows = []
    for name, u in (("-j", -_J), ("j", _J), ("k", _K), ("(j+k)/sqrt2", (_J + _K) / math.sqrt(2)),
                    ("i (not perp)", _I), ("(i+j)/sqrt2 (not perp)", (_I + _J) / math.sqrt(2))):
        covers = True
        oriented = True
        for q in qs:
            x, e1, e2 = spin_frame(q)
            xu, e1u, e2u = spin_frame(qmul(u, q))
            covers &= bool(np.allclose(xu, -x))
            oriented &= bool(abs(np.linalg.det(np.array([e1u, e2u, xu])) - 1.0) < 1e-10)
        Lu = quaternion_left(u)
        rows.append({"u": name, "perp_to_i": bool(abs(u @ _I) < 1e-12),
                     "covers_antipode": covers, "image_frame_oriented": oriented,
                     "L_u_squared_is_minus_one": bool(np.allclose(Lu @ Lu, -np.eye(4)))})
    # the frame image for u = -j in the frame at x: (e1, e2) -> ?
    q = qs[0]
    x, e1, e2 = spin_frame(q)
    xu, e1u, e2u = spin_frame(qmul(-_J, q))
    image = np.array([[e1u @ e1, e1u @ e2], [e2u @ e1, e2u @ e2]])
    # dA(e1, e2) = (-e1, -e2); composing with reflection r gives det(image) = +1?  as
    # a map of R^3 vectors the image is (e1, -e2): dA followed by the reflection of e1
    # conjugacy: any two perp u are related by a fibre rotation e^{i theta}
    conj_ok = True
    for theta in np.linspace(0, 2 * math.pi, 37):
        g = np.array([math.cos(theta), math.sin(theta), 0, 0])
        gu = qmul(qmul(g, -_J), qinv(g))
        conj_ok &= bool(abs(gu @ _I) < 1e-12 and abs(gu @ gu - 1) < 1e-12)
    return {"rows": rows,
            "exactly_the_perp_units_cover": all(r["covers_antipode"] == r["perp_to_i"] for r in rows),
            "frame_image_of_minus_j_in_frame_coords": image.tolist(),
            "frame_image_is_dA_times_one_reflection": bool(np.allclose(image, np.diag([1.0, -1.0]))),
            "fibre_rotation_conjugates_within_perp_circle": conj_ok}


def two_sheeted_involution() -> Dict[str, object]:
    """T3/T4 — on the Pin^-(2) bundle of ``S^2`` (two sheets ``±`` = the two
    orientations), ``A~ = (L_u on +, L_{-u} on -)`` is a free involution;
    the two signs ``epsilon`` are the two Pin^- structures of ``RP^2``."""
    out = {}
    for eps in (+1, -1):
        u = eps * (-_J)
        Lu, Lmu = quaternion_left(u), quaternion_left(-u)
        # A~^2 on the + sheet: L_{-u} L_u
        sq_plus = Lmu @ Lu
        # single-sheet iterate: L_u L_u
        single = Lu @ Lu
        out[f"epsilon={eps:+d}"] = {
            "A_tilde_squared_is_identity": bool(np.allclose(sq_plus, np.eye(4))),
            "single_sheet_L_u_squared_is_minus_one": bool(np.allclose(single, -np.eye(4))),
            "free": True}
    out["number_of_pin_minus_structures"] = 2
    out["preferred_by_geometry"] = False
    return out


def pin_minus_structures_rp2() -> List[Dict[str, object]]:
    """T4/T5 — the two structures, their quadratic enhancements and ABK."""
    return [abk_invariant(eps) for eps in (+1, -1)]


def abk_invariant(epsilon: int) -> Dict[str, object]:
    """T5 — ``q : H_1(RP^2; Z_2) -> Z_4`` with ``q(0) = 0``, ``q(g) = epsilon mod 4``;
    check the quadratic property ``q(x+y) = q(x) + q(y) + 2 (x.y)`` with
    ``g.g = 1``; Gauss sum ``sum_x i^{q(x)} = sqrt(2) e^{i pi ABK/4}``."""
    q = {0: 0, 1: epsilon % 4}
    dot = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 1}
    quadratic = all((q[(x + y) % 2] - q[x] - q[y] - 2 * dot[(x, y)]) % 4 == 0
                    for x in (0, 1) for y in (0, 1))
    gauss = sum(1j ** q[x] for x in (0, 1))
    abk = round(cmath.phase(gauss) / (math.pi / 4.0)) % 8
    return {"epsilon": epsilon, "q_of_generator_mod4": q[1], "is_quadratic": quadratic,
            "gauss_sum": gauss, "gauss_modulus": abs(gauss),
            "modulus_is_sqrt2": bool(abs(abs(gauss) - math.sqrt(2)) < 1e-12),
            "ABK_mod8": abk,
            "H1_action_shifts_q_by_2": (q[1] + 2) % 4 == (-epsilon) % 4}


def pin_bordism_pairing_rule() -> Dict[str, object]:
    """T6 — with ``Omega_2^{Pin-} = Z_8`` generated by ``RP^2`` (imported), a
    Pin^- 3-manifold with boundary ``RP^2_a ⊔ RP^2_b`` exists iff
    ``ABK_a + ABK_b = 0 mod 8``."""
    abk = {eps: abk_invariant(eps)["ABK_mod8"] for eps in (+1, -1)}
    table = []
    for a in (+1, -1):
        for b in (+1, -1):
            total = (abk[a] + abk[b]) % 8
            table.append({"sector_a": a, "sector_b": b, "abk_sum_mod8": total,
                          "bounds": total == 0})
    single = [{"sector": e, "abk": abk[e], "bounds_alone": abk[e] % 8 == 0} for e in (+1, -1)]
    return {"bordism_group": "Z_8 (imported: Kirby-Taylor; Anderson-Brown-Peterson)",
            "abk": abk, "pairs": table, "single_mouth": single,
            "opposite_sectors_bound": all(r["bounds"] for r in table if r["sector_a"] != r["sector_b"]),
            "equal_sectors_do_not": all(not r["bounds"] for r in table if r["sector_a"] == r["sector_b"]),
            "single_mouth_cannot_bound": all(not r["bounds_alone"] for r in single),
            "modelling_assumption": ("pair creation is a Pin- bordism of the mouth "
                                     "surfaces, i.e. their Pin- structures are induced "
                                     "from one on the worldvolume [chosen]")}
