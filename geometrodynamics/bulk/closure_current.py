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

*Scope of "derived" (review finding 3).* What R1 derives is the **quaternionic
reduction** of the *chosen* ``source -> A -> J -> B -> source`` itinerary,
under the *chosen* geodesic realignment at the detectors, to the triangle
``x -> u -> -v -> x``, together with its holonomy. The itinerary and the
realignment remain model choices; only the reduction and the holonomy are
derived.
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
    "holonomy_weighted_law", "sector_prior_control", "closure_phase_gradient_on_closure_set",
    "singlet_loop_law", "oriented_sector_prior_control",
    "integrated_sector_weights", "oriented_current_audit", "underived_inputs",
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


def singlet_loop_law(gamma: float, n: int = 200001) -> Dict[str, object]:
    """R2b — the **actual** Pin-derived loop, computed directly.

    R1 shows the derived history is ``x -> u -> -v -> x`` with full holonomy
    ``-sgn D(x, u, -v)``, so its oriented current is

        -D(x, u, -v) / (2 |u x (-v)|)

    integrated over the closure circle of that triangle (the same great
    circle, since ``x . (u x -v) = -x . (u x v)``). With
    ``int_Gamma x dsigma = 0`` this gives ``2 pi (1 - u.v)`` per sector; the
    common factor ``-1`` from ``J^2`` and the common ``1/(2|u x v|)`` cancel in
    the normalisation, leaving

        P(s_A, s_B) = (1 - s_A s_B cos gamma) / 4 ,      E = -cos gamma

    the singlet. This is tested here on the derived object rather than
    inferred from the triplet ``D(x, u, v)`` by a verbal sign substitution --
    the gap that review finding 2 identified.
    """
    a = np.array([0.0, 0.0, 1.0])
    b = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    weights = {}
    for sA in (1, -1):
        for sB in (1, -1):
            u, w = sA * a, -(sB * b)                 # the derived third vertex is -v
            x = great_circle(u, w, n)
            Dp = 1.0 + x @ u + float(u @ w) + x @ w
            weights[(sA, sB)] = -float(np.mean(Dp)) * 2.0 * math.pi / (
                2.0 * float(np.linalg.norm(np.cross(u, w))))
    Z = sum(weights.values())
    P = {k: val / Z for k, val in weights.items()}
    analytic = {k: (1.0 - k[0] * k[1] * math.cos(gamma)) / 4.0 for k in P}
    E = sum(sa * sb * P[(sa, sb)] for (sa, sb) in P)
    return {"P": P, "analytic": analytic,
            "max_deviation": max(abs(P[k] - analytic[k]) for k in P),
            "E": E, "E_is_minus_cos": bool(abs(E + math.cos(gamma)) < 1e-9),
            "all_weights_positive_after_normalisation": all(v > 0 for v in P.values()),
            "computed_on_the_derived_loop": True}


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


def oriented_sector_prior_control(gamma: float = 1.0,
                                  ratios=(0.5, 1.0, 2.0)) -> Dict[str, object]:
    """R3b — the sector prior is a gap for the **oriented** branch as well.

    The oriented sector integrals are proportional to ``1 + u.v``, i.e.
    ``1 + s_A s_B cos gamma`` (triplet loop) or ``1 - s_A s_B cos gamma``
    (the derived singlet loop). Under a prior ratio
    ``r = pi_like / pi_unlike``,

        E_r^triplet = [r(1+cos g) - (1-cos g)] / [r(1+cos g) + (1-cos g)]
        E_r^singlet = [r(1-cos g) - (1+cos g)] / [r(1-cos g) + (1+cos g)]

    which equal ``+-cos gamma`` **only at r = 1**. So the equal prior is
    load-bearing on both sides of the fork: the quantum joint law is not
    recovered from the holonomy weighting alone. Marginals stay ``1/2`` for
    every ``r`` (the prior respects the ``(x,u,v) -> (-x,-u,-v)`` symmetry),
    so no-signalling does not constrain it.
    """
    c = math.cos(gamma)
    rows = []
    for r in ratios:
        trip = (r * (1 + c) - (1 - c)) / (r * (1 + c) + (1 - c))
        sing = (r * (1 - c) - (1 + c)) / (r * (1 - c) + (1 + c))
        rows.append({"ratio": r, "E_triplet": trip, "E_singlet": sing,
                     "equals_cos": bool(abs(trip - c) < 1e-12),
                     "equals_minus_cos": bool(abs(sing + c) < 1e-12),
                     "P(A=+)": 0.5})
    return {"rows": rows,
            "quantum_law_only_at_ratio_one": all(
                (abs(r["ratio"] - 1.0) < 1e-12) == r["equals_cos"] for r in rows),
            "marginals_stay_half": True,
            "status": ("the equal sector prior is chosen and load-bearing on BOTH branches: "
                       "the holonomy weighting alone does not give the quantum joint law")}


def closure_phase_gradient_on_closure_set(gamma: float = 1.0, n: int = 2001,
                                          sA: int = 1, sB: int = 1,
                                          step: float = 1e-6) -> Dict[str, object]:
    """R4 (analytic) — the phase-stationarity **proxy**, and why it is disjoint
    from sharp closure.

    With ``N = x.q``, ``q = u x v``, and ``D = A + x.p``, ``A = 1 + u.v``,
    ``p = u + v``, the closure phase is ``Omega = 2 arg(D + iN)`` and

        grad Omega = 2 (D grad N - N grad D) / (D^2 + N^2).

    On the closure set ``N = 0`` the tangential gradient of ``N`` is ``q``
    itself (``q`` is already tangent there, since ``x.q = 0``), so

        grad_{S^2} Omega |_{N=0}  =  2 q / D  =  2 (u x v) / D ,

    of norm ``2|u x v| / |D| > 0`` for non-collinear analyzers. Hence sharp
    closure and stationarity of the phase are **disjoint** away from the two
    points where ``D = 0``, and those points are ``x = -u`` and ``x = -v``,
    where ``D = N = 0``: the ``arg`` chart is singular there, which is not
    stationarity.

    The finite-difference check wraps the increment of ``arg`` into
    ``(-pi, pi]``: ``Omega`` itself jumps by ``4 pi`` across the closure
    circle wherever ``D < 0`` (the ``pi`` branch), while the holonomy
    ``e^{i Omega/2} = -1`` is continuous there.

    **This is a proxy and is labelled as one.** The repository's fourth
    closure condition is *stationarity of an action*, which is not
    implemented anywhere (see :func:`stationarity_audit`); nothing here
    tests it.
    """
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    u, v = sA * a, sB * b
    q, p, A = np.cross(u, v), u + v, 1.0 + float(u @ v)
    x = great_circle(u, v, n)
    D = A + x @ p
    keep = np.abs(D) > 1e-2
    xk, Dk = x[keep], D[keep]
    analytic = 2.0 * q[None, :] / Dk[:, None]
    analytic = analytic - (np.sum(analytic * xk, axis=1)[:, None]) * xk

    def _arg(y):
        y = y / np.linalg.norm(y, axis=-1, keepdims=True)
        return np.arctan2(y @ q, A + y @ p)

    fd = np.zeros_like(xk)
    for k in range(3):
        d = np.zeros(3)
        d[k] = 1.0
        tangent = d[None, :] - (xk @ d)[:, None] * xk
        raw = _arg(xk + step * tangent) - _arg(xk - step * tangent)
        fd[:, k] = 2.0 * ((raw + math.pi) % (2.0 * math.pi) - math.pi) / (2.0 * step)
    fd = fd - (np.sum(fd * xk, axis=1)[:, None]) * xk
    norms = np.linalg.norm(analytic, axis=1)
    singular = [float(abs(A + (-u) @ p)) + float(abs((-u) @ q)),
                float(abs(A + (-v) @ p)) + float(abs((-v) @ q))]
    return {"analytic_gradient": "2 (u x v) / D",
            "finite_difference_residual": float(np.max(np.abs(fd - analytic))),
            "min_gradient_norm_on_closure_set": float(norms.min()),
            "cross_product_norm": float(np.linalg.norm(q)),
            "gradient_never_vanishes_on_closure_set": bool(norms.min() > 1e-6),
            "singular_points_are_minus_u_and_minus_v": singular,
            "singular_points_are_chart_not_stationary": bool(max(singular) < 1e-12),
            "closure_and_phase_stationarity_are_disjoint": bool(norms.min() > 1e-6),
            "classification": ("phase-stationarity proxy: incompatible with sharp phase closure "
                               "for generic (non-collinear) settings. It does NOT test the "
                               "repository's unimplemented extremal-action condition.")}


def integrated_sector_weights(gamma: float = 1.0, n: int = 200001) -> Dict[str, object]:
    """The oriented current does not produce negative event probabilities.

    Although ``D`` changes sign pointwise on the closure circle, with
    ``int_Gamma x dsigma = 0`` the integrated sector weight is

        int_Gamma D_s dsigma = 2 pi (1 + u.v)  >=  0

    for every outcome pair, vanishing only at ``u = -v``. So the
    holonomy-weighted construction uses destructive cancellation
    *internally* and yields non-negative normalised sector weights, which
    removes the naive "negative probabilities" objection.

    **The analogy stays qualified** (review finding 4). Non-negativity of the
    integrated current is not yet an event-frequency law. If ``D`` is an
    amplitude or a current, a classical detector normally responds to a
    quadratic (energy) functional of it; if it is an event measure, the signed
    cancellation needs a reason. Even if local coefficients eventually force
    the sign, one still has to derive why observed frequencies are *linear* in
    the integrated current rather than quadratic in it, or given by some other
    readout functional. That is a third open item, not a corollary of the
    second.
    """
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    rows, worst = [], 0.0
    for sA in (1, -1):
        for sB in (1, -1):
            u, v = sA * a, sB * b
            x = great_circle(u, v, n)
            D = 1.0 + x @ u + float(u @ v) + x @ v
            quad = float(np.mean(D) * 2.0 * math.pi)
            exact = 2.0 * math.pi * (1.0 + float(u @ v))
            worst = max(worst, abs(quad - exact))
            rows.append({"sector": (sA, sB), "integral": quad, "exact_2pi_1_plus_udotv": exact,
                         "nonnegative": bool(exact >= -1e-12),
                         "negative_arc_fraction": float(np.mean(D < 0))})
    return {"rows": rows, "max_quadrature_residual": worst,
            "all_sector_integrals_nonnegative": all(r["nonnegative"] for r in rows),
            "cancellation_is_internal": True,
            "reading": ("destructive cancellation within a sector, non-negative normalised "
                        "sector weights: this removes the naive negative-probability "
                        "objection, and the wave analogy stays qualified"),
            "still_open_linear_versus_quadratic_readout": (
                "why observed frequencies would be LINEAR in the integrated current rather "
                "than quadratic in it (the usual classical detector response to an amplitude) "
                "or another readout functional -- a third open item, not a corollary")}


def oriented_current_audit() -> Dict[str, object]:
    """An open mathematical direction, recorded as an audit item and **not** a
    success criterion of the pre-registration.

    If the physical observable is the integral of a section of the local
    system that the Pin/Hopf data define over the closure locus, the sign
    ``e^{i Omega/2}`` is geometrically mandatory and the oriented current is
    forced. If the physical object is instead a measure on histories,
    positivity forces ``|D|``. That is a sharper statement of the fork than
    "quasi-probability versus probability", and neither side is established
    here.
    """
    return {"question": ("do the Pin/Hopf data make the closure locus naturally an oriented "
                         "current with local coefficients, whose observables are integrals of "
                         "sections, or is the physical object a measure on histories?"),
            "if_local_system": "the sign is geometrically mandatory: HOLONOMY_WEIGHTED_COAREA",
            "if_measure_on_histories": "positivity forces |D|: POSITIVE_COAREA",
            "status": "open; recorded as an audit item, not a success criterion",
            "established_here": False,
            "would_not_by_itself_complete_the_derivation": (
                "even if the sign were forced, the sector prior (R3b) and the "
                "current-to-frequency readout (linear vs quadratic) would remain open")}


def underived_inputs() -> List[str]:
    """The headline, narrowed after review: **three** things remain underived,
    not one binary choice."""
    return [
        "branch aggregation: positive count |D| versus oriented sum e^{i Omega/2}|D| = D",
        "the relative sector coefficients (the equal outcome-sector prior), which move the "
        "correlation on BOTH branches and are fixed by no symmetry of the model",
        "the readout: why observed event frequencies would be linear in the integrated "
        "current rather than quadratic in it or another functional",
    ]


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
    proxy = closure_phase_gradient_on_closure_set()
    return {"named_in_module_docstring": in_docstring,
            "repository_condition": "stationarity: the history has extremal action",
            "implemented": implemented,
            "branch_sign_used_in_weight": False,
            "weight_form": "exp(-mismatch^2 / (2 sigma^2)), sigma = 0.6, positive for both branches",
            "phase_stationarity_proxy": proxy,
            "proxy_tests_the_repository_condition": False,
            "stationarity_decides_the_fork": False,
            "why": ("The repository's fourth condition is extremal action and there is no action "
                    "functional in the module. Substituting stationarity of the geometric phase "
                    "would be a new assumption, and analytically it is incompatible with sharp "
                    "closure anyway: grad Omega = 2(u x v)/D never vanishes on the closure set. "
                    "So no variational principle available in the repository can choose between "
                    "positive and holonomy-weighted coarea.")}


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
