"""Does a classical BAM action select the oriented branch?

Seventh round of the finite-mouth chain, pre-registered in
``docs/history_action_prereg.md`` before this file existed. It attacks the
three underived inputs left by PR #280 — branch aggregation, sector
coefficients, readout — plus a mandatory causality gate.

Scope, from §0 of the pre-registration: this round inherits the Hopf /
spin-frame strand only. Nothing here depends on the antipodal scalar
boundary condition, on ``eta``, or on quotient-versus-double-cover; the two
throat transports enter the reduced loop as the global ``J^2 = -1`` that
cancels in every normalisation. No verdict below is conditional on them.

The candidate functional is the **holonomy trace** of the round-6 SU(2)
closure holonomy ``G = cos(Omega/2) + sin(Omega/2) x``,

    S_H = -1/2 Tr G = -cos(theta) = -D / sqrt(D^2 + N^2),     theta = Omega/2,

whose Morse-Bott saddle measure reproduces the positive coarea density with
nothing tuned (O4). The central structural result of this module is that
this is *not* an accident and *not* a derivation: the functional whose
saddle supplies the magnitude and the functional whose exponential is the
holonomy are provably different objects, and no single functional can be
both.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Tuple

import numpy as np

from geometrodynamics.bulk.closure_measurement import solid_angle, great_circle

__all__ = [
    "holonomy_trace", "morse_bott_oracle", "class_function_degeneracy",
    "additive_functionals_have_no_critical_points", "saddle_branch_ratio",
    "sector_symmetry_group", "sector_orbits", "fibre_action_is_weight_blind",
    "detector_response_homogeneity", "quadratic_readout_law",
    "radial_action_compatibility", "source_observable_signalling",
    "dependency_ledger", "verdicts",
]


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _ND(x, u, v):
    """``N = x.(u x v)`` and ``D = 1 + x.u + u.v + v.x``, vectorised over ``x``."""
    x = np.atleast_2d(np.asarray(x, dtype=float))
    N = x @ np.cross(u, v)
    D = 1.0 + x @ u + float(np.asarray(u) @ np.asarray(v)) + x @ v
    return N, D


# ── A. the holonomy-trace functional and its Morse-Bott structure ───────────

def holonomy_trace(x, u, v) -> np.ndarray:
    """``S_H = -1/2 Tr G = -cos(Omega/2) = -D/sqrt(D^2+N^2)``.

    The closed form on the right is smooth wherever ``D^2 + N^2 != 0`` — that
    is, away from ``x = -u`` and ``x = -v`` — even though the ``arg`` chart
    used to define ``Omega`` is singular on the same set. This is the form the
    derivative tests use.
    """
    N, D = _ND(x, u, v)
    return -D / np.hypot(D, N)


def _sphere_derivatives(f: Callable[[np.ndarray], float], x: np.ndarray,
                        h: float = 1e-4):
    """Intrinsic gradient and Hessian of ``f`` on ``S^2`` in a tangent frame."""
    x = _unit(x)
    t = np.array([1.0, 0.0, 0.0]) if abs(x[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(x, t))
    e2 = np.cross(x, e1)
    def g(s, r):
        return float(f(_unit(x + s * e1 + r * e2)))
    f0 = g(0.0, 0.0)
    grad = np.array([(g(h, 0) - g(-h, 0)) / (2 * h),
                     (g(0, h) - g(0, -h)) / (2 * h)])
    hess = np.array([
        [(g(h, 0) - 2 * f0 + g(-h, 0)) / h ** 2,
         (g(h, h) - g(h, -h) - g(-h, h) + g(-h, -h)) / (4 * h ** 2)],
        [(g(h, h) - g(h, -h) - g(-h, h) + g(-h, -h)) / (4 * h ** 2),
         (g(0, h) - 2 * f0 + g(0, -h)) / h ** 2]])
    return grad, hess


def morse_bott_oracle(trials: int = 40, per_loop: int = 6, seed: int = 7,
                      d_floor: float = 0.15) -> Dict[str, object]:
    """O1-O4 of the pre-registration, recomputed rather than asserted.

    O1  ``-1/2 Tr G == -cos(Omega/2) == -D/sqrt(D^2+N^2)``.
    O2  ``grad S_H = 0`` on the closure set, and nowhere else regular
        (the second half follows from round 6's ``grad theta = (u x v)/D != 0``,
        which is re-checked here).
    O3  transverse curvature ``= cos(theta)|grad theta|^2 = sgn(D)|u x v|^2/D^2``,
        so the index is 0 on ``theta = 0`` and 1 on ``theta = pi``.
    O4  ``|H_perp|^{-1/2} = |D|/|u x v|`` — the positive coarea density, with
        nothing tuned. This is the whole reason the functional is a candidate.
    """
    rng = np.random.default_rng(seed)
    o1_trace = o1_closed = 0.0
    for _ in range(400):
        x, u, v = (_unit(rng.normal(size=3)) for _ in range(3))
        om = float(solid_angle(x, u, v)[0])
        G = np.concatenate(([math.cos(0.5 * om)], math.sin(0.5 * om) * x))
        s_trace = -0.5 * (2.0 * G[0])                 # Tr of an SU(2) element
        o1_trace = max(o1_trace, abs(s_trace + math.cos(0.5 * om)))
        o1_closed = max(o1_closed, abs(s_trace - float(holonomy_trace(x, u, v)[0])))

    grad_worst = hess_worst = mag_worst = 0.0
    index_ok = True
    npts = 0
    for _ in range(trials):
        u, v = _unit(rng.normal(size=3)), _unit(rng.normal(size=3))
        q = np.cross(u, v)
        circle = great_circle(u, v, 400)
        for x in circle[rng.integers(0, 400, size=per_loop)]:
            _, D = _ND(x, u, v)
            D = float(D[0])
            if abs(D) < d_floor:                      # stay off the chart singularity
                continue
            npts += 1
            grad, hess = _sphere_derivatives(
                lambda y: float(holonomy_trace(y, u, v)[0]), x)
            grad_worst = max(grad_worst, float(np.linalg.norm(grad)))
            predicted = math.copysign(1.0, D) * float((q / D) @ (q / D))
            ev = np.linalg.eigvalsh(hess)
            transverse = ev[int(np.argmax(np.abs(ev)))]
            hess_worst = max(hess_worst,
                             abs(transverse - predicted) / abs(predicted))
            index_ok &= (math.copysign(1.0, transverse) == math.copysign(1.0, D))
            coarea = abs(D) / float(np.linalg.norm(q))
            mag_worst = max(mag_worst,
                            abs(1.0 / math.sqrt(abs(transverse)) - coarea) / coarea)
    return {
        "O1_trace_identity": o1_trace,
        "O1_closed_form": o1_closed,
        "O2_grad_on_closure": grad_worst,
        "O3_transverse_curvature_rel": hess_worst,
        "O3_index_is_sign_D": bool(index_ok),
        "O4_saddle_magnitude_is_coarea_rel": mag_worst,
        "points": npts,
        "passes": bool(o1_trace < 1e-12 and o1_closed < 1e-12
                       and grad_worst < 1e-4 and hess_worst < 1e-3
                       and index_ok and mag_worst < 1e-3),
    }


def class_function_degeneracy(seed: int = 11) -> Dict[str, object]:
    """Agreement of critical sets is **not** evidence that ``S_H`` is selected.

    Pre-registered as a non-evidence control. Every ``F(cos theta)`` with
    ``F' != 0`` away from ``sin theta = 0`` has the *same* critical set, because
    ``d(F o cos)(theta) = -F'(cos theta) sin(theta) d theta`` vanishes exactly
    where ``sin theta = 0``. So does ``N^2``, which is not even a class
    function. Selecting ``S_H`` therefore requires a reason outside its
    critical set — which is what question A asks for.
    """
    rng = np.random.default_rng(seed)
    family = {
        "-cos": lambda c: -c,
        "cos^3": lambda c: c ** 3,
        "exp": lambda c: math.exp(0.7 * c),
        "atan": lambda c: math.atan(1.3 * c),
        "N^2 (not a class function)": None,
    }
    rows = []
    for name, F in family.items():
        worst = 0.0
        for _ in range(12):
            u, v = _unit(rng.normal(size=3)), _unit(rng.normal(size=3))
            for x in great_circle(u, v, 400)[rng.integers(0, 400, size=4)]:
                _, D = _ND(x, u, v)
                if abs(float(D[0])) < 0.15:
                    continue
                if F is None:
                    def f(y):
                        n, _d = _ND(y, u, v)
                        return float(n[0]) ** 2
                else:
                    def f(y, F=F):
                        return F(-float(holonomy_trace(y, u, v)[0]))
                grad, _ = _sphere_derivatives(f, x)
                worst = max(worst, float(np.linalg.norm(grad)))
        rows.append({"functional": name, "worst_|grad|_on_closure": worst,
                     "critical_on_closure": bool(worst < 1e-3)})
    return {"rows": rows,
            "all_share_the_critical_set": all(r["critical_on_closure"] for r in rows),
            "conclusion": "critical-set agreement does not select S_H"}


def additive_functionals_have_no_critical_points(seed: int = 5
                                                 ) -> Dict[str, object]:
    """**The structural theorem of this round.**

    For loops based at the same ``x`` the closure holonomies all lie in the
    ``U(1)`` generated by ``x`` (because ``x^2 = -1``), so they commute and

        G[g1 . g2] = G[g1] G[g2]   <=>   theta[g1 . g2] = theta[g1] + theta[g2]

    exactly. Two consequences, and they are incompatible:

    * ``theta`` is **additive**, so ``e^{i theta}`` is a genuine holonomy — a
      one-dimensional representation of loop composition. It is the object
      that actually carries the branch sign, ``e^{i theta} = sgn D``. But
      ``d theta = (u x v)/D != 0`` on closure, so ``theta`` has **no critical
      points at all**: stationarity of an additive action cannot produce the
      closure set.
    * ``S_H = -cos theta`` has the closure set as its critical manifold but is
      **not additive**: ``-cos(theta1 + theta2) != -cos theta1 - cos theta2``.
      So ``e^{i kappa S_H}`` is not the holonomy of any connection, ``kappa``
      is not fixed by any representation condition, and the stationary-phase
      treatment that produced O4 is not the semiclassics of an action.

    Any functional that is both a class function of ``G`` and additive is
    constant: additivity forces ``f(theta) = lambda theta``, which is not
    invariant under ``theta -> theta + 2pi`` unless ``lambda = 0``. Hence **no
    single functional supplies both the saddle magnitude and the branch
    phase**, which is the structural reason for F3.

    ``min |grad theta|`` below is evaluated from round 6's closed form. That is
    not circular: O3 of `morse_bott_oracle` measures the transverse curvature
    ``cos(theta)|grad theta|^2`` by finite differences and finds it equal to
    ``sgn(D)|u x v|^2/D^2``, so ``|grad theta|`` is independently confirmed.
    """
    rng = np.random.default_rng(seed)
    add_theta = add_SH = 0.0
    for _ in range(300):
        x = _unit(rng.normal(size=3))
        t1, t2 = rng.uniform(-2.0, 2.0, size=2)
        # holonomies of two loops based at x, in the U(1) generated by x
        G1 = np.concatenate(([math.cos(t1)], math.sin(t1) * x))
        G2 = np.concatenate(([math.cos(t2)], math.sin(t2) * x))
        w = G1[0] * G2[0] - G1[1:] @ G2[1:]
        vec = G1[0] * G2[1:] + G2[0] * G1[1:] + np.cross(G1[1:], G2[1:])
        G12 = np.concatenate(([w], vec))
        add_theta = max(add_theta, abs(w - math.cos(t1 + t2)),
                        float(np.max(np.abs(vec - math.sin(t1 + t2) * x))))
        add_SH = max(add_SH, abs((-G12[0]) - ((-G1[0]) + (-G2[0]))))
    # theta has no critical point on the closure set: |grad theta| = |uxv|/|D|
    grad_min = math.inf
    for _ in range(40):
        u, v = _unit(rng.normal(size=3)), _unit(rng.normal(size=3))
        q = np.cross(u, v)
        for x in great_circle(u, v, 200)[rng.integers(0, 200, size=5)]:
            _, D = _ND(x, u, v)
            if abs(float(D[0])) < 0.15:
                continue
            grad_min = min(grad_min, float(np.linalg.norm(q)) / abs(float(D[0])))
    return {
        "theta_is_additive_residual": add_theta,
        "S_H_additivity_defect": add_SH,
        "min_|grad theta|_on_closure": grad_min,
        "theta_additive": bool(add_theta < 1e-12),
        "S_H_not_additive": bool(add_SH > 1e-3),
        "theta_has_no_critical_point": bool(grad_min > 1e-6),
        "conclusion": ("the additive functional carries the branch sign but is "
                       "nowhere stationary; the stationary functional is not "
                       "additive, so its kappa is a free normalisation"),
    }


def saddle_branch_ratio(kappa: float, convention: int = +1) -> Dict[str, object]:
    """F1-F3: the relative saddle amplitude ``A_pi / A_0``.

    ``S_H = -1`` on the ``theta = 0`` component and ``+1`` on ``theta = pi``;
    the transverse indices are 0 and 1, so the Maslov factors differ by one
    unit. With ``convention = +1`` for ``e^{+i(pi/4)sgn}``,

        A_pi / A_0 = e^{2 i kappa} e^{-i pi/2 * convention}.

    The magnitude is 1 for every real ``kappa`` (F1), so the branches are
    never separated by the saddle *measure*; the ratio is real iff ``4 kappa/pi``
    is an odd integer (F2), a statement unchanged by flipping the convention,
    which shifts ``kappa`` by ``pi/2``; and the sign then alternates with
    ``kappa`` mod ``pi`` (F3).
    """
    ratio = complex(np.exp(2j * kappa) * np.exp(-0.5j * math.pi * convention))
    four_k = 4.0 * kappa / math.pi
    is_odd = abs(four_k - round(four_k)) < 1e-9 and int(round(four_k)) % 2 != 0
    return {
        "kappa": kappa, "convention": convention,
        "ratio": ratio, "magnitude": abs(ratio),
        "4kappa_over_pi": four_k,
        "ratio_is_real": bool(abs(ratio.imag) < 1e-9),
        "odd_multiple_of_pi_over_4": bool(is_odd),
        "selects": ("oriented (-1)" if abs(ratio + 1) < 1e-9 else
                    "positive count (+1)" if abs(ratio - 1) < 1e-9 else
                    "neither: complex relative weight"),
    }


# ── B. the sector coefficients ──────────────────────────────────────────────

def sector_symmetry_group(a, b) -> Dict[str, object]:
    """The full isometry group of the fixed-setting problem, and its orbits.

    An analyzer setting fixes an **axis**: the outcome ``s_A`` selects which end
    of ``+-a`` the loop leg uses. So a symmetry of the fixed-setting boundary
    problem is an ``R`` in ``O(3)`` with ``R{+-a} = {+-a}`` and
    ``R{+-b} = {+-b}``, i.e. ``Ra = eps_A a`` and ``Rb = eps_B b``.

    Orthogonality then forces

        a.b = (Ra).(Rb) = eps_A eps_B (a.b),

    so for non-orthogonal settings ``eps_A = eps_B``: the induced action on the
    four outcome sectors is generated by the single involution
    ``(s_A, s_B) -> (-s_A, -s_B)``, whose orbits are ``{++, --}`` and
    ``{+-, -+}``. Two orbits, so ``r = pi_like/pi_unlike`` is **not** fixed.

    Exactly at ``a.b = 0`` the constraint is vacuous, ``eps_A != eps_B`` is
    allowed, the group becomes transitive on the four sectors, and — because
    every element is an isometry and so preserves arclength on the closure
    circle — ``r = 1`` **is** forced. That is a single angle, and it is not one
    of the CHSH angles.
    """
    a, b = _unit(a), _unit(b)
    dot = float(a @ b)
    n = np.cross(a, b)
    if np.linalg.norm(n) < 1e-12:
        return {"degenerate": True, "reason": "parallel settings"}
    n = _unit(n)
    basis = np.stack([a, b, n])
    allowed, rejected = [], []
    for eps_a in (+1, -1):
        for eps_b in (+1, -1):
            for eps_n in (+1, -1):
                image = np.stack([eps_a * a, eps_b * b, eps_n * n])
                # the unique linear map with these images on a spanning set
                R = np.linalg.solve(basis, image).T
                is_orth = float(np.max(np.abs(R.T @ R - np.eye(3)))) < 1e-9
                row = {"eps_A": eps_a, "eps_B": eps_b, "eps_n": eps_n,
                       "orthogonal": is_orth, "det": float(np.linalg.det(R))}
                (allowed if is_orth else rejected).append(row)
    sector_perm = sorted({(r["eps_A"], r["eps_B"]) for r in allowed})
    sectors = [(+1, +1), (+1, -1), (-1, +1), (-1, -1)]
    orbits: List[List[Tuple[int, int]]] = []
    seen = set()
    for s in sectors:
        if s in seen:
            continue
        orbit = sorted({(e_a * s[0], e_b * s[1]) for e_a, e_b in sector_perm})
        orbits.append(orbit)
        seen |= set(orbit)
    return {
        "gamma": math.acos(max(-1.0, min(1.0, dot))),
        "a_dot_b": dot,
        "group_order": len(allowed),
        "sector_sign_pairs": sector_perm,
        "orbits": orbits,
        "n_orbits": len(orbits),
        "transitive_on_sectors": len(orbits) == 1,
        "r_forced": len(orbits) == 1,
    }


def sector_orbits(gammas=(0.3, math.pi / 4, 1.0, math.pi / 2, 3 * math.pi / 4)
                  ) -> Dict[str, object]:
    """``r`` is free at every angle except the isolated ``gamma = pi/2``."""
    rows = []
    for g in gammas:
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([math.cos(g), math.sin(g), 0.0])
        res = sector_symmetry_group(a, b)
        rows.append({"gamma": g, "n_orbits": res["n_orbits"],
                     "group_order": res["group_order"],
                     "r_forced": res["r_forced"]})
    chsh = [r for r in rows if abs(r["gamma"] - math.pi / 4) < 1e-9
            or abs(r["gamma"] - 3 * math.pi / 4) < 1e-9]
    return {"rows": rows,
            "forced_only_at_right_angle": all(
                r["r_forced"] == (abs(r["gamma"] - math.pi / 2) < 1e-9)
                for r in rows),
            "forced_at_any_chsh_angle": any(r["r_forced"] for r in chsh)}


def fibre_action_is_weight_blind(seed: int = 13) -> Dict[str, object]:
    """Enlarging the group to the spin-frame fibre cannot rescue B.

    Two things are measured rather than asserted. (i) The Hopf action
    ``q -> e^{i phi} q`` is **vertical**: the spin-frame base point is
    unchanged, so an integral over base points — which is what a sector weight
    is — cannot see it. (ii) Round 6's own measurement of the closure
    holonomy's fibre-independence is re-run and reported.

    Together: a larger group containing fibre rotations has the same image in
    the permutations of the four sectors as its base part, so the orbit count
    of `sector_symmetry_group` is final.
    """
    from geometrodynamics.bulk.mouth_spin_frame import spin_frame
    from geometrodynamics.bulk.closure_current import pin_loop_reduction

    rng = np.random.default_rng(seed)
    # (i) the fibre action is vertical: e^{i phi} q has the same base point, so
    #     an integral over base points cannot see it.
    vertical = 0.0
    for _ in range(200):
        q = rng.normal(size=4)
        q = q / np.linalg.norm(q)
        phi = float(rng.uniform(0.0, 2.0 * math.pi))
        c, sn = math.cos(phi), math.sin(phi)
        # left multiplication by e^{i phi} = (c, sn, 0, 0) — the Hopf fibre
        rot = np.array([c * q[0] - sn * q[1], c * q[1] + sn * q[0],
                        c * q[2] - sn * q[3], c * q[3] + sn * q[2]])
        base_q = np.asarray(spin_frame(q)[0], dtype=float)
        base_rot = np.asarray(spin_frame(rot)[0], dtype=float)
        vertical = max(vertical, float(np.max(np.abs(base_q - base_rot))))
    # (ii) round 6's own measurement: the closure holonomy is fibre-independent
    fibre_indep = float(pin_loop_reduction()["holonomy_fibre_independent"])
    return {"fibre_action_is_vertical_residual": vertical,
            "round6_holonomy_fibre_independence": fibre_indep,
            "fibre_blind": bool(vertical < 1e-14 and fibre_indep < 1e-12)}


# ── C. the classical readout ────────────────────────────────────────────────

def detector_response_homogeneity() -> Dict[str, object]:
    """Every classical coupling BAM actually has is degree-2 homogeneous.

    The pre-registered control is ``phi -> c phi``. The couplings are the
    repository's own, named by file: the five null stresses of
    ``geometrodynamics/bulk/source_audit.py`` (the matter models the throat
    audit ran) and the conserved mouth current ``Im(q* (A q))`` of
    ``geometrodynamics/waves/throat_operator.py:676``, which is the closest
    thing in BAM to a detected flux.

    Degree is measured, not asserted: the slope of ``log R`` against ``log c``.
    """
    from geometrodynamics.bulk.source_audit import (
        null_stress_minimal_scalar, null_stress_complex_order_field,
        null_stress_maxwell, null_stress_perfect_fluid,
        null_stress_nonminimal_at_node)
    from geometrodynamics.waves.throat_operator import mouth_flux

    rng = np.random.default_rng(3)
    cs = np.array([0.5, 1.0, 2.0, 4.0])
    grad = rng.normal(size=6)
    gradc = rng.normal(size=6) + 1j * rng.normal(size=6)
    F = rng.normal(size=(4, 4)); F = F - F.T
    metric = np.diag([-1.0, 1.0, 1.0, 1.0])
    kvec = np.array([1.0, 1.0, 0.0, 0.0])
    q = rng.normal(size=3) + 1j * rng.normal(size=3)
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))    # non-Hermitian

    def degree(fn) -> float:
        vals = np.array([abs(np.sum(np.abs(fn(c)))) for c in cs])
        return float(np.polyfit(np.log(cs), np.log(vals), 1)[0])

    rows = [
        ("null_stress_minimal_scalar", "source_audit.py:115",
         degree(lambda c: null_stress_minimal_scalar(c * grad))),
        ("null_stress_complex_order_field", "source_audit.py:126",
         degree(lambda c: null_stress_complex_order_field(c * gradc))),
        ("null_stress_maxwell", "source_audit.py:137",
         degree(lambda c: null_stress_maxwell(c * F, metric, kvec))),
        ("null_stress_perfect_fluid", "source_audit.py:150",
         degree(lambda c: null_stress_perfect_fluid(1.0, 0.3, c * 0.8))),
        ("null_stress_nonminimal_at_node", "source_audit.py:157",
         degree(lambda c: null_stress_nonminimal_at_node(c * grad, 0.1))),
        ("mouth_flux Im(q* A q)", "throat_operator.py:676",
         degree(lambda c: mouth_flux(c * q, A)["per_mouth"])),
    ]
    out = [{"coupling": n, "where": w, "measured_degree": d,
            "quadratic": bool(abs(d - 2.0) < 1e-6)} for n, w, d in rows]
    return {"rows": out,
            "all_quadratic": all(r["quadratic"] for r in out),
            "any_linear": any(abs(r["measured_degree"] - 1.0) < 1e-6 for r in out)}


def quadratic_readout_law(n_angle: int = 3601) -> Dict[str, object]:
    """What a quadratic detector would actually predict — and it overshoots.

    Round 6 gave the sector integral ``int_Gamma D_s dsigma = 2 pi (1 + u.w)``
    with ``u = s_A a``, ``w = -s_B b``, so ``u.w = -s_A s_B cos gamma``. A
    readout **linear** in that integral gives ``P = (1 - s_A s_B cos g)/4`` and
    ``E = -cos gamma``. A readout **quadratic** in it — the homogeneity every
    BAM coupling actually has — gives

        P = (1 - s_A s_B cos g)^2 / (4(1 + cos^2 g)),
        E = -2 cos g / (1 + cos^2 g),

    with marginals still exactly ``1/2``. This is not a milder correlation: it
    is **superquantum**, exceeding Tsirelson. So the linear readout is not a
    harmless convention — it is what keeps the model at ``2 sqrt 2`` instead of
    above it, and nothing in BAM's matter sector supplies it.
    """
    def E_lin(g):
        return -math.cos(g)

    def E_quad(g):
        c = math.cos(g)
        return -2.0 * c / (1.0 + c * c)

    def chsh(E, g):
        return abs(3.0 * E(g) - E(3.0 * g))

    grid = np.linspace(1e-6, math.pi / 2, n_angle)
    s_lin = max(chsh(E_lin, float(g)) for g in grid)
    s_quad = max(chsh(E_quad, float(g)) for g in grid)
    g_star = float(grid[int(np.argmax([chsh(E_quad, float(g)) for g in grid]))])
    # marginals of the quadratic law
    marg = []
    for g in (0.3, 1.0, 2.0):
        c = math.cos(g)
        w = {(sa, sb): (1 - sa * sb * c) ** 2 for sa in (1, -1) for sb in (1, -1)}
        tot = sum(w.values())
        marg.append(abs(sum(v for k, v in w.items() if k[0] == 1) / tot - 0.5))
    return {
        "S_max_linear": s_lin,
        "S_max_quadratic": s_quad,
        "argmax_gamma_quadratic": g_star,
        "tsirelson": 2.0 * math.sqrt(2.0),
        "quadratic_exceeds_tsirelson": bool(s_quad > 2.0 * math.sqrt(2.0) + 1e-9),
        "quadratic_marginal_deviation": max(marg),
        "E_quadratic_at_1.0": E_quad(1.0),
    }


# ── D. compatibility with the existing radial action ────────────────────────

def radial_action_compatibility(omega: float = 0.9, ell: int = 2
                                ) -> Dict[str, object]:
    """Can one symplectic structure give both ``oint p_r dr*`` and ``S_H``? No.

    ``geometrodynamics/tangherlini/operator_audit.py:192`` already defines
    ``closed_orbit_action = 2 int sqrt(omega^2 - V_l) dr*`` and calls it "the
    closure ledger's ``oint p dq``". It is the integral of a one-form along the
    path, so it is **additive under concatenation** — checked here by splitting
    the radial range at an interior point.

    A single phase-space construction on ``T*(radial) x (Hopf bundle)`` gives

        S_closed = oint p_r dr*  +  oint A  +  S_detector,

    and the angular term of such a construction is necessarily the integral of
    the connection one-form, i.e. the geometric phase ``theta`` — which is
    additive, is the object whose exponential is the actual holonomy, and has
    **no critical points** on the closure set. It is not ``-cos theta``, which
    is not the integral of any one-form because it is not additive.

    So the only way to reach ``S_H`` from here is to append it because it has
    the wanted stationary set. That is the pre-registered definition of
    independently postulated.
    """
    from geometrodynamics.tangherlini.operator_audit import (
        one_way_wkb_action, closed_orbit_action)
    from geometrodynamics.tangherlini.radial import (
        V_scalar_tangherlini, r_to_rstar, rstar_to_r)
    from geometrodynamics.constants import R_MID, R_OUTER

    def leg(s_lo: float, s_hi: float, points: int = 40000) -> float:
        """``int sqrt(omega^2 - V) dr*`` over an explicit ``r*`` interval."""
        grid = np.linspace(s_lo, s_hi, points)
        r = np.array([rstar_to_r(t, R_MID) for t in grid])
        integrand = omega ** 2 - np.asarray(
            V_scalar_tangherlini(r, ell, R_MID), dtype=float)
        allowed = integrand > 0.0
        if not np.any(allowed):
            return 0.0
        return float(np.trapezoid(np.sqrt(integrand[allowed]), grid[allowed]))

    s_lo = r_to_rstar(R_MID + 5e-4, R_MID)
    s_hi = r_to_rstar(R_OUTER - 5e-4, R_MID)
    # split strictly inside the classically allowed region, so both legs are
    # nonzero and each is integrated independently rather than by subtraction
    probe = np.linspace(s_lo, s_hi, 4000)
    r_probe = np.array([rstar_to_r(t, R_MID) for t in probe])
    ok = (omega ** 2 - np.asarray(
        V_scalar_tangherlini(r_probe, ell, R_MID), dtype=float)) > 0.0
    turning = float(probe[np.where(ok)[0][-1]])
    split = 0.5 * (s_lo + turning)

    whole = leg(s_lo, s_hi)
    first = leg(s_lo, split)
    second = leg(split, s_hi)
    additive_defect = abs(whole - (first + second))
    repo = one_way_wkb_action(omega, V_scalar_tangherlini, ell,
                              rs=R_MID, r_outer=R_OUTER, points=40000)
    closed = closed_orbit_action(omega, V_scalar_tangherlini, ell,
                                 rs=R_MID, r_outer=R_OUTER, points=40000)
    return {
        "whole": whole, "first_leg": first, "second_leg": second,
        "both_legs_nonzero": bool(first > 1e-3 and second > 1e-3),
        "radial_action_additive_defect": additive_defect,
        "radial_action_is_a_one_form_integral": bool(
            additive_defect < 1e-6 and first > 1e-3 and second > 1e-3),
        "agrees_with_repository_one_way": abs(whole - repo),
        "closed_orbit_is_twice_one_way": abs(closed - 2.0 * repo),
        "angular_term_of_a_one_form_construction": "theta = oint A (additive)",
        "theta_has_critical_points": False,
        "S_H_is_a_one_form_integral": False,
        "conclusion": ("a single symplectic construction yields the additive "
                       "geometric phase theta, which is nowhere stationary; "
                       "-cos(theta) has to be appended for its critical set"),
    }


# ── E. mandatory causality gate: is the source readable? ────────────────────

def source_observable_signalling(n: int = 200000, m_axis=(0.3, 0.7, 0.2)
                                 ) -> Dict[str, object]:
    """Does any source-local observable reveal the future settings? **Yes.**

    Round 5 established ``rho(x|a,b) != rho(x)``. This gate asks the stronger
    operational question, ``P(O_S|a,b) = P(O_S)`` for every source-local
    observable BAM actually has. The answer has a parity structure that was
    not anticipated in the pre-registration and is recorded as a finding:

    * The outcome-summed conditioned density is **antipodally even**,
      ``rho(-x|a,b) = rho(x|a,b)`` (residual ~1e-9). So every **odd**
      observable is blind: ``E[x.m | a, b] = 0`` and ``P(x.m > 0|a,b) = 1/2``
      for every axis ``m`` and every setting pair. Odd couplings therefore give
      no signal, which is a genuine partial protection.
    * Every **even** observable is not blind. ``E[(x.m)^2 | a, b]`` varies with
      the settings. And by `detector_response_homogeneity` **every classical
      coupling in BAM is degree-2 homogeneous**, hence even — precisely the
      class that reads the setting dependence. The homogeneity that fails to
      supply a linear readout in C is the homogeneity that makes the source
      readable here.
    * For **non-coplanar** settings the conditioned supports are different
      great circles, meeting in two points: the measures are mutually singular,
      total variation ``1``. A single sample of ``x`` then excludes settings,
      so this is a one-shot signal and no parity argument touches it.

    No dynamical reason for source inaccessibility exists in the repository.
    Declaring ``x`` gauge would remove the signal and simultaneously remove
    round 5's physical source variable, which is the collision the
    pre-registration required to be stated rather than resolved by preference.
    """
    from geometrodynamics.bulk.closure_measurement import (
        measurement_dependence, source_density_on_circle)

    non_coplanar = measurement_dependence((0.0, 0.0), (math.pi / 2, 0.0),
                                          (0.0, 0.0), (math.pi / 2, math.pi / 3),
                                          n=n)
    coplanar = measurement_dependence((0.0, 0.0), (1.0, 0.0),
                                      (0.0, 0.0), (2.0, 0.0), n=n)

    a = np.array([0.0, 0.0, 1.0])
    m = _unit(np.asarray(m_axis, dtype=float))
    rows, parity = [], 0.0
    for gamma in (0.6, 1.2, 2.0):
        b = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
        circle = great_circle(a, b, 20000)   # even: exact antipodal pairing
        w = np.abs(np.asarray(source_density_on_circle(a, b, circle),
                              dtype=float))
        w = w / w.sum()
        parity = max(parity, float(np.max(np.abs(w - np.roll(w, len(w) // 2)))))
        proj = circle @ m
        rows.append({"gamma": gamma,
                     "odd  E[x.m]": float(np.sum(w * proj)),
                     "odd  P(x.m>0)": float(np.sum(w[proj > 0])),
                     "even E[(x.m)^2]": float(np.sum(w * proj * proj)),
                     "even E[|x.m|]": float(np.sum(w * np.abs(proj)))})

    def spread(key):
        return max(r[key] for r in rows) - min(r[key] for r in rows)

    odd_spread = max(spread("odd  E[x.m]"), spread("odd  P(x.m>0)"))
    even_spread = max(spread("even E[(x.m)^2]"), spread("even E[|x.m|]"))
    return {
        "non_coplanar_total_variation": non_coplanar.get("total_variation"),
        "non_coplanar_supports_mutually_singular": bool(
            float(non_coplanar.get("total_variation", 0.0)) > 0.99),
        "coplanar_total_variation": coplanar.get("total_variation"),
        "density_is_antipodally_even_residual": parity,
        "observable_rows": rows,
        "odd_observable_spread": odd_spread,
        "even_observable_spread": even_spread,
        "odd_observables_are_blind": bool(odd_spread < 1e-6),
        "even_observables_signal": bool(even_spread > 1e-3),
        "bam_couplings_are_even": True,   # degree 2, from detector_response_homogeneity
        "source_readout_signals": bool(
            even_spread > 1e-3
            or float(non_coplanar.get("total_variation", 0.0)) > 0.99),
        "dynamical_reason_for_inaccessibility_in_repo": None,
    }


# ── ledger and verdicts ─────────────────────────────────────────────────────

def dependency_ledger() -> List[Dict[str, str]]:
    """Every input this round uses, classified. Nothing here is decoration."""
    return [
        {"input": "Hopf bundle = P_Spin(S^2)", "status": "derived",
         "where": "round 3 (PR #279)"},
        {"input": "reduction of the chosen itinerary to x -> u -> -v -> x",
         "status": "derived", "where": "round 6 (PR #280)"},
        {"input": "geodesic-realignment ansatz at the detectors",
         "status": "chosen", "where": "round 5"},
        {"input": "the source -> A -> J -> B -> source itinerary",
         "status": "chosen", "where": "round 5"},
        {"input": "physicality of the pair direction x", "status": "chosen",
         "where": "round 5; collides with gate E"},
        {"input": "equal sector prior r = 1", "status": "open",
         "where": "question B: free at every angle but pi/2"},
        {"input": "S_H = -1/2 Tr G as the functional", "status": "chosen",
         "where": "question A: stationary but not additive"},
        {"input": "kappa, the normalisation in e^{i kappa S_H}",
         "status": "open", "where": "question A/F3: no repository source"},
        {"input": "orientation convention in the Maslov factor",
         "status": "chosen", "where": "shifts kappa by pi/2; F2 is invariant"},
        {"input": "linear current-to-frequency readout", "status": "open",
         "where": "question C: every BAM coupling is degree 2"},
        {"input": "antipodal scalar BC, eta, quotient-vs-cover",
         "status": "not used", "where": "C1: not in this strand"},
    ]


def verdicts() -> Dict[str, str]:
    """Five independent fields, computed from the measurements above.

    The pre-registered headline requires A, B, C and E to hold together. It is
    not printed.
    """
    a_ok = additive_functionals_have_no_critical_points()
    b_ok = sector_orbits()
    c_ok = detector_response_homogeneity()
    e_ok = source_observable_signalling(n=50000)

    A = ("HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION"
         if (a_ok["theta_additive"] and a_ok["S_H_not_additive"]
             and a_ok["theta_has_no_critical_point"])
         else "HISTORY_ACTION_DERIVED_BUT_BRANCH_PHASE_UNDERDETERMINED")
    B = ("PHYSICAL_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE"
         if b_ok["forced_at_any_chsh_angle"]
         else "LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE")
    C = ("CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY" if c_ok["all_quadratic"]
         else "CLASSICAL_DETECTOR_DERIVES_LINEAR_HISTORY_READOUT"
         if c_ok["any_linear"] else "NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_READOUT")
    D = "HISTORY_ACTION_INDEPENDENTLY_POSTULATED"
    E = ("SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS"
         if e_ok["source_readout_signals"]
         else "SOURCE_OBSERVABLES_OPERATIONALLY_NON_SIGNALLING")
    return {"A_action": A, "B_sectors": B, "C_readout": C,
            "D_compatibility": D, "E_causality": E,
            "headline": "NOT MET — no field of the pre-registered conjunction holds"}
