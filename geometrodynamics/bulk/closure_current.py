"""Positive coarea count or holonomy-weighted coarea current?

Pre-registered in ``docs/closure_current_prereg.md`` (``f954e3d``), committed
before this file. It carries four corrections to the previous round and asks
whether anything in the repository's classical structures selects between
the two parameter-free measures on the closure set:

* ``POSITIVE_COAREA``  --  ``|D| / (2|u x v|)``,           S_max = 2.1423
* ``HOLONOMY_WEIGHTED_COAREA``  --  ``e^{i Omega/2} |D| / (2|u x v|) = D / (2|u x v|)``,  S_max = 2 sqrt 2

The rule fixed there: the oriented current may be adopted only if the Pin
structure or a stationary classical equation requires it, never because it
gives ``2 sqrt 2``.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

from geometrodynamics.bulk.mouth_spin_frame import spin_frame, qmul, qinv
from geometrodynamics.bulk.closure_measurement import (
    solid_angle, great_circle, closed_form_weights, closure_weights, chsh,
)

__all__ = [
    "minimal_rotation_lift", "pin_loop_reduction", "branch_holonomy_is_sign_D",
    "holonomy_weighted_law", "sector_prior_control", "stationary_set_of_closure_phase",
    "stationarity_audit", "pin_label_versus_weight", "verdict",
]

_J = np.array([0.0, 0.0, -1.0, 0.0])          # the mouth transport L_{-j}


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def minimal_rotation_lift(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """``g = cos(theta/2) + sin(theta/2) n``, ``n = (p x q)/|p x q|``, ``theta = angle(p, q)``:
    the SU(2) lift of the geodesic (parallel) transport from ``p`` to ``q``."""
    p, q = _unit(p), _unit(q)
    th = math.acos(float(np.clip(p @ q, -1.0, 1.0)))
    n = np.cross(p, q)
    if np.linalg.norm(n) < 1e-12:
        n = np.cross(p, [1.0, 0, 0]) if abs(p[0]) < 0.9 else np.cross(p, [0, 1.0, 0])
    n = _unit(n)
    return np.array([math.cos(0.5 * th), *(math.sin(0.5 * th) * n)])


def _ad(g: np.ndarray, v: np.ndarray) -> np.ndarray:
    return qmul(qmul(g, np.array([0.0, *v])), qinv(g))[1:]


def _angle_about(g: np.ndarray, x: np.ndarray) -> float:
    ref = np.array([0.3, 0.5, 0.7])
    t = _unit(np.cross(x, ref) if abs(x @ ref) < 0.95 else np.cross(x, [1.0, 0, 0]))
    r = _ad(g, t)
    return math.atan2(float(np.cross(t, r) @ x), float(t @ r))


def pin_loop_reduction(gamma: float = 1.0, n: int = 400, seed: int = 2) -> Dict[str, object]:
    """R1 — the history ``source -> A -> J -> B -> source`` in the quaternion model.

    Frame rotation by ``R = Ad_g`` is ``q -> q g^{-1}`` (right); ``J`` is left
    multiplication; they commute, so ``q_5 = J [J q_0 g_1^{-1}] g_2^{-1} g_3^{-1}
    = -q_0 (g_3 g_2 g_1)^{-1}``. The frame holonomy is ``Ad_G``, a rotation about
    ``x`` by ``Omega(x; u, -v)``, and ``G = cos(Omega/2) + sin(Omega/2) x``.
    """
    rng = np.random.default_rng(seed)
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    worst_q5 = worst_phi = worst_angle = worst_lift = 0.0
    for _ in range(n):
        q0 = _unit(rng.standard_normal(4))
        x = spin_frame(q0)[0]
        u, v = rng.choice([1, -1]) * a, rng.choice([1, -1]) * b
        g1, g2, g3 = minimal_rotation_lift(x, u), minimal_rotation_lift(-u, v), minimal_rotation_lift(v, -x)
        q1 = qmul(q0, qinv(g1)); q2 = qmul(_J, q1); q3 = qmul(q2, qinv(g2))
        q4 = qmul(q3, qinv(g3)); q5 = qmul(_J, q4)
        G = qmul(qmul(g3, g2), g1)
        worst_q5 = max(worst_q5, float(np.max(np.abs(q5 + qmul(q0, qinv(G))))))
        # fibre independence of the SO(3) holonomy
        phi = rng.uniform(0, 2 * math.pi)
        q0p = qmul(np.array([math.cos(phi), math.sin(phi), 0, 0]), q0)
        q5p = qmul(_J, qmul(qmul(qmul(_J, qmul(q0p, qinv(g1))), qinv(g2)), qinv(g3)))
        H, Hp = qmul(q5, qinv(q0)), qmul(q5p, qinv(q0p))
        worst_phi = max(worst_phi, float(np.max(np.abs(_ad(Hp, [1, 0, 0]) - _ad(H, [1, 0, 0])))))
        Om = float(solid_angle(x[None, :], u, -v)[0])
        ang = _angle_about(G, x)
        worst_angle = max(worst_angle, abs((ang - Om + math.pi) % (2 * math.pi) - math.pi))
        worst_lift = max(worst_lift, float(np.max(np.abs(G - np.array([math.cos(0.5 * Om), *(math.sin(0.5 * Om) * x)])))))
    return {"q5_equals_minus_q0_Ginv": worst_q5, "holonomy_fibre_independent": worst_phi,
            "frame_holonomy_is_Omega_of_triangle_x_u_minus_v": worst_angle,
            "lift_is_cos_half_Omega_plus_sin_half_Omega_x": worst_lift,
            "reduces_to_triangle": bool(max(worst_q5, worst_phi, worst_angle, worst_lift) < 1e-12),
            "partner_sign": "v -> -v: the singlet sign", "extra_from_J_squared": -1}


def branch_holonomy_is_sign_D(gamma: float = 1.0, n: int = 4000) -> Dict[str, object]:
    """C1 — on the closure circle ``e^{i Omega/2} = sgn D`` (away from ``D = 0``),
    and the derived loop's spinor holonomy is ``-sgn D(x, u, -v)``."""
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    x = great_circle(a, b, n)
    worst = worst_lift = 0.0
    frac_pi_branch = {}
    for sA in (1, -1):
        for sB in (1, -1):
            u, v = sA * a, sB * b
            Om = solid_angle(x, u, v)
            D = 1.0 + x @ u + float(u @ v) + x @ v
            keep = np.abs(D) > 1e-6
            hol = np.exp(0.5j * Om[keep])
            worst = max(worst, float(np.max(np.abs(hol - np.sign(D[keep])))))
            frac_pi_branch[(sA, sB)] = float(np.mean(D < 0))
            # the derived loop: G on the circle equals sgn D(x, u, -v) as a quaternion
            for xi in x[::400]:
                G = qmul(qmul(minimal_rotation_lift(v, -xi), minimal_rotation_lift(-u, v)),
                         minimal_rotation_lift(xi, u))
                Dm = 1.0 + xi @ u + float(u @ (-v)) + (-v) @ xi
                if abs(Dm) > 1e-3:
                    worst_lift = max(worst_lift, float(np.max(np.abs(G - np.array([np.sign(Dm), 0, 0, 0])))))
    return {"exp_i_half_Omega_equals_sign_D": worst, "loop_lift_equals_sign_D_on_circle": worst_lift,
            "pi_branch_fraction_by_sector": frac_pi_branch,
            "signed_density_is_holonomy_times_coarea": bool(worst < 1e-12)}


def holonomy_weighted_law(gamma: float) -> Dict[str, object]:
    """R2 — ``W~(s) ∝ ∫_Γ D_s dσ / (2|u x v|)`` with ``∫_Γ x dσ = 0`` gives
    ``P_like = (1 + cos γ)/4``, ``P_unlike = (1 − cos γ)/4``, ``E = cos γ``
    (triplet sign; ``−cos γ`` with the partner sign of R1)."""
    P = closure_weights(gamma, "signed")
    analytic = {(1, 1): (1 + math.cos(gamma)) / 4, (-1, -1): (1 + math.cos(gamma)) / 4,
                (1, -1): (1 - math.cos(gamma)) / 4, (-1, 1): (1 - math.cos(gamma)) / 4}
    worst = max(abs(P[k] - analytic[k]) for k in P)
    E = sum(sa * sb * P[(sa, sb)] for (sa, sb) in P)
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    x = great_circle(a, b, 20000)
    return {"P": P, "analytic": analytic, "max_deviation": worst, "E": E,
            "E_is_cos": bool(abs(E - math.cos(gamma)) < 1e-9),
            "circle_mean_of_x": float(np.max(np.abs(x.mean(axis=0)))),
            "no_projectors_used": True}


def sector_prior_control(gamma: float = 1.0, ratios=(0.5, 1.0, 2.0)) -> Dict[str, object]:
    """R3 — under ``pi_like / pi_unlike = r`` the marginals stay ``1/2`` and
    ``E = (r W_l − W_u)/(r W_l + W_u)`` moves. The equal prior is the counting
    measure on sectors: chosen."""
    c, s = math.cos(0.5 * gamma), math.sin(0.5 * gamma)
    Wl, Wu = 2.0 + c * (math.pi - gamma) / s, 2.0 + s * gamma / c
    rows = []
    for r in ratios:
        Z = 2 * (r * Wl + Wu)
        P = {(1, 1): r * Wl / Z, (-1, -1): r * Wl / Z, (1, -1): Wu / Z, (-1, 1): Wu / Z}
        rows.append({"ratio": r, "E": (r * Wl - Wu) / (r * Wl + Wu),
                     "P(A=+)": P[(1, 1)] + P[(1, -1)]})
    return {"rows": rows, "marginals_stay_half": all(abs(r["P(A=+)"] - 0.5) < 1e-12 for r in rows),
            "E_moves": bool(abs(rows[0]["E"] - rows[-1]["E"]) > 0.1),
            "symmetry_fixing_ratio": None,
            "status": "equal prior = counting measure on sectors: chosen"}


def stationary_set_of_closure_phase(gamma: float = 1.0, n_theta: int = 721, n_phi: int = 1441,
                                    sA: int = 1, sB: int = 1) -> Dict[str, object]:
    """R4 — the stationary set ``∇_{S²} Ω = 0`` of the closure phase, found on a
    grid, against the closure set ``N = 0``. ``∇Ω = 2 (D ∇N − N ∇D)/(N² + D²)``."""
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    u, v = sA * a, sB * b
    th = np.linspace(1e-3, math.pi - 1e-3, n_theta)
    ph = np.linspace(0.0, 2 * math.pi, n_phi, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    x = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
    w = np.cross(u, v)
    N = x @ w
    D = 1.0 + x @ u + float(u @ v) + x @ v
    gradN = w[None, None, :] - N[..., None] * x
    gradD = (u + v)[None, None, :] - (x @ (u + v))[..., None] * x
    grad = 2.0 * (D[..., None] * gradN - N[..., None] * gradD) / (N * N + D * D)[..., None]
    gnorm = np.linalg.norm(grad, axis=-1)
    # local minima of |grad| below tolerance
    tol = 2e-2
    cand = np.argwhere(gnorm < tol)
    pts = x[gnorm < tol]
    # cluster
    clusters: List[np.ndarray] = []
    for p in pts:
        if all(np.linalg.norm(p - c) > 0.1 for c in clusters):
            clusters.append(p)
    dist_to_circle = [float(abs(c @ _unit(w))) for c in clusters]     # |x . n_Gamma| = sin of angular distance to Γ
    # Lexell: the level set Omega = const is a circle through -u and -v (the two
    # singular points D = N = 0); check coplanarity of a sampled level set
    level = 0.7
    sel = np.abs(2.0 * np.arctan2(N, D) - level) < 2e-3          # the signed level set
    pts = x[sel]
    if len(pts) > 10:
        centroid = pts.mean(axis=0)
        centred = pts - centroid
        U_, sv, Vt = np.linalg.svd(centred, full_matrices=False)
        coplanarity = float(sv[-1] / sv[0])
        normal = Vt[-1]
        # the fitted plane (hence the Lexell circle) must contain -u and -v;
        # grid samples cannot resolve the level set near those singular points
        singular_on_it = [float(abs((-u - centroid) @ normal)), float(abs((-v - centroid) @ normal))]
    else:
        coplanarity, singular_on_it = float("nan"), []
    singular_D_N = [float(abs(1.0 + (-u) @ u + u @ v + v @ (-u))), float(abs(1.0 + (-v) @ u + u @ v + v @ (-v)))]
    return {"n_stationary_points": len(clusters), "points": [c.tolist() for c in clusters],
            "sin_distance_to_closure_circle": dist_to_circle,
            "no_stationary_points": len(clusters) == 0,
            "coincides_with_closure_set": bool(len(clusters) > 0 and max(dist_to_circle) < 1e-3),
            "min_grad_on_grid": float(gnorm.min()),
            "lexell": {"level_set_is_a_circle (plane residual)": coplanarity,
                       "level_set_plane_contains_minus_u_minus_v (plane distances)": singular_on_it,
                       "D_at_minus_u_minus_v": singular_D_N}}


def stationarity_audit() -> Dict[str, object]:
    """R4 — is stationarity implemented anywhere in ``history/closure.py``?"""
    import inspect
    from geometrodynamics.history import closure as hc
    src = inspect.getsource(hc)
    # identifiers (not prose) that would implement it: any name containing
    # 'stationar', 'extremal' or 'action' (other than 'transaction'), collected
    # from the syntax tree so that docstrings and comments do not count
    import ast
    tree = ast.parse(src)
    idents = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            idents.add(node.name)
        elif isinstance(node, ast.Name):
            idents.add(node.id)
        elif isinstance(node, ast.Attribute):
            idents.add(node.attr)
    implemented = any(("stationar" in n.lower() or "extremal" in n.lower()
                       or ("action" in n.lower() and "transaction" not in n.lower())) for n in idents)
    in_docstring = "Stationarity" in src.split('"""')[1]
    stat = stationary_set_of_closure_phase()
    return {"named_in_module_docstring": in_docstring, "implemented": implemented,
            "branch_sign_used_in_weight": False,
            "weight_form": "exp(-mismatch^2 / (2 sigma^2)), sigma = 0.6, positive for both branches",
            "stationary_set": stat,
            "stationarity_decides_the_fork": False}


def pin_label_versus_weight() -> Dict[str, object]:
    """R5 — what the Pin structure supplies: the branch holonomy as a label on
    each closed history. Whether the label enters the measure is not a bundle
    property; no classical structure in the repository selects it."""
    return {"supplied_by_pin_structure": "closure holonomy -sgn D(x,u,-v) on every closed history",
            "selected_by_pin_structure": "nothing about the measure",
            "selected_by_implemented_classical_dynamics": "nothing (stationarity unimplemented; field equations linear, no history measure)",
            "the_open_step": "positive count versus oriented (holonomy-weighted) sum over the 0/pi closure branches",
            "candidates": {"POSITIVE_COAREA": {"E(1)": float(closed_form_weights(1.0)[(1, 1)] * 2 - closed_form_weights(1.0)[(1, -1)] * 2), "S_max": chsh()},
                           "HOLONOMY_WEIGHTED_COAREA": {"E(1)": math.cos(1.0), "S_max": chsh(variant="signed")}}}


def verdict(reduces: bool, holonomy_is_sign_D: bool, born_from_current: bool,
            stationarity_decides: bool, pin_selects_weight: bool) -> str:
    if not (reduces and holonomy_is_sign_D and born_from_current):
        return "NEITHER"
    if pin_selects_weight or stationarity_decides:
        return "HOLONOMY_WEIGHTED_COAREA"
    return "FORK_UNRESOLVED_BY_CURRENT_STRUCTURES"
