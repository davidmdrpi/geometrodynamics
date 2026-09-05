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
this is *not* an accident and *not* a derivation: within the class of real
class functions of the holonomy, the functional whose saddle supplies the
magnitude and the functional whose exponential is the holonomy are provably
different objects, and no single one can be both. That is a statement about
that class, not a no-go for classical actions in general.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Tuple

import numpy as np

from geometrodynamics.bulk.closure_measurement import solid_angle, great_circle

__all__ = [
    "holonomy_trace", "morse_bott_oracle", "class_function_degeneracy",
    "additive_functionals_have_no_critical_points", "saddle_branch_ratio",
    "morse_bott_component_masses", "no_off_closure_critical_points",
    "amplitude_dependence", "excision_estimate",
    "sector_symmetry_group", "sector_orbits", "fibre_action_is_weight_blind",
    "discrete_symmetry_extension",
    "detector_response_homogeneity", "quadratic_readouts_disagree",
    "local_square_mean_is_closed_form",
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
    invariant under ``theta -> theta + 2pi`` unless ``lambda = 0``.

    **Scope (note N9).** What this establishes is precisely: *no globally
    single-valued real class function of the holonomy alone is simultaneously
    additive and has the closure manifold as its stationary set.* It is **not**
    a no-go for classical actions in general. A full classical action need not
    be a class function of ``G`` at all — it may contain ``int p dq - H dt``,
    detector interaction terms, orientation-dependent potentials or boundary
    terms. Nothing here forbids some unified construction of that kind; what is
    shown is that the *holonomy-trace* route, the one the ``S_H`` coincidence
    suggests, cannot supply both roles at once.

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


def morse_bott_component_masses(u, v, n: int = 200001, amplitude=None
                                ) -> Dict[str, object]:
    """The two Morse-Bott component masses, and what they already are.

    For a critical *manifold* rather than an isolated point, stationary phase
    carries an integral over the component,

        A_j ~ e^{i kappa S_j} e^{i pi sigma_j/4} (2pi/kappa)^{1/2} M_j,
        M_j = int_{C_j} a(y) / sqrt(|H_perp(y)|) dsigma.

    With O4's ``|H_perp|^{-1/2} = |D|/|u x v|`` the two masses are the positive
    coarea masses of the two closure arcs — the ``D > 0`` arc (``theta = 0``)
    and the ``D < 0`` arc (``theta = pi``). They are **generically unequal**,
    and

        (M_0 - M_pi)|u x v| = int_Gamma D dsigma      (the oriented sum)
        (M_0 + M_pi)|u x v| = int_Gamma |D| dsigma    (the positive count)

    so stationary phase reproduces **both candidate aggregations exactly**,
    leaving only their relative phase undetermined.

    **Conditional on the amplitude (note N12).** Those identities hold for
    ``a = 1`` with uniform arclength — i.e. conditional on the *inherited
    Haar/uniform source measure* of round 5, which is a choice, not a
    derivation. ``amplitude`` accepts a callable ``a(x)`` to demonstrate the
    dependence: a non-constant in-plane ``a`` moves both masses and the
    aggregation with them (e.g. ``a = 1 + 0.5 x.u`` takes ``M_0`` from
    ``11.6708`` to ``14.4736`` and the oriented sum from ``9.6780`` to
    ``12.0975``). An ``a`` varying only along ``u x v`` is invisible, since
    ``x.(u x v) = 0`` on the closure circle. So "with nothing tuned" applies to
    the *Hessian*, not to the measure.

    **Scope.** This is stationary phase over the two-dimensional family of
    geodesic triangles parameterised by ``x in S^2``, not over BAM's full loop
    or history space.
    """
    u, v = _unit(u), _unit(v)
    q = np.cross(u, v)
    nq = float(np.linalg.norm(q))
    circle = great_circle(u, v, n)
    _, D = _ND(circle, u, v)
    a = (np.ones(len(circle)) if amplitude is None
         else np.asarray([float(amplitude(x)) for x in circle], dtype=float))
    dsigma = 2.0 * math.pi / len(circle)
    pos, neg = D > 0, D < 0
    M0 = float(np.sum(a[pos] * np.abs(D[pos]))) * dsigma / nq
    Mpi = float(np.sum(a[neg] * np.abs(D[neg]))) * dsigma / nq
    signed = float(np.sum(a * D)) * dsigma
    absolute = float(np.sum(a * np.abs(D))) * dsigma
    return {
        "M_0": M0, "M_pi": Mpi, "mass_ratio": Mpi / M0,
        "masses_are_unequal": bool(abs(Mpi / M0 - 1.0) > 1e-3),
        "oriented_identity_residual": abs((M0 - Mpi) * nq - signed),
        "positive_count_identity_residual": abs((M0 + Mpi) * nq - absolute),
        "int_D": signed, "int_absD": absolute,
        "amplitude_is_unity": amplitude is None,
    }


def amplitude_dependence(gamma: float = 1.0) -> Dict[str, object]:
    """The masses, hence the aggregation, depend on the unproved amplitude.

    Pre-registered nowhere; recorded as note N12 after review. The round-5
    Haar/uniform source measure enters the Morse-Bott coefficients as ``a = 1``,
    and a different preparation density changes what stationary phase returns.
    """
    u = np.array([0.0, 0.0, 1.0])
    v = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    q = np.cross(u, v)
    base = morse_bott_component_masses(u, v, n=200001)
    rows = [{"amplitude": "1 (round-5 Haar)", "M_0": base["M_0"],
             "M_pi": base["M_pi"], "oriented": (base["M_0"] - base["M_pi"])
             * float(np.linalg.norm(q))}]
    for label, fn in (
            ("1 + 0.5 x.u  (in plane)", lambda x: 1.0 + 0.5 * float(x @ u)),
            ("1 + 0.8 x.(uxv)/|uxv| (normal)",
             lambda x: 1.0 + 0.8 * float(x @ q) / float(np.linalg.norm(q)))):
        m = morse_bott_component_masses(u, v, n=200001, amplitude=fn)
        rows.append({"amplitude": label, "M_0": m["M_0"], "M_pi": m["M_pi"],
                     "oriented": (m["M_0"] - m["M_pi"])
                     * float(np.linalg.norm(q))})
    in_plane_moves = abs(rows[1]["oriented"] - rows[0]["oriented"]) > 1e-3
    normal_invisible = abs(rows[2]["oriented"] - rows[0]["oriented"]) < 1e-9
    return {"rows": rows,
            "in_plane_amplitude_moves_the_aggregation": bool(in_plane_moves),
            "normal_amplitude_is_invisible": bool(normal_invisible),
            "masses_are_conditional_on_the_measure": bool(in_plane_moves)}


def excision_estimate(gamma: float = 1.0, n: int = 2000001,
                      epsilons=(0.2, 0.1, 0.05, 0.025)) -> Dict[str, object]:
    """The punctures carry no mass, so the Morse-Bott coefficients are honest.

    ``S_H`` is undefined at ``x = -u, -v``, so the two critical components are
    open arcs ending at singular points of the phase, where textbook Morse-Bott
    stationary phase wants smoothness. In local tangent coordinates at ``x = -u``,
    ``D ~ sin(gamma) xi`` and ``N ~ sin(gamma) eta``, so
    ``S_H ~ -xi/sqrt(xi^2 + eta^2)`` really is direction-dependent there.

    The excision argument closes it: on the closure circle ``|D| ~ sin(gamma) s``
    in arclength ``s`` from the puncture, so the mass inside a disc of radius
    ``eps`` is

        2 int_0^eps sin(gamma) s ds / |u x v| = sin(gamma) eps^2 / |u x v|,

    the factor two because a disc of radius ``eps`` about the puncture covers
    the arc on **both** sides of it — i.e. **O(eps^2)**. Excising the punctures therefore changes neither mass in
    the limit, and the coefficients quoted are the genuine ones rather than
    formal symbols. Measured ratios below approach the closed form.
    """
    u = np.array([0.0, 0.0, 1.0])
    v = np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    nq = float(np.linalg.norm(np.cross(u, v)))
    circle = great_circle(u, v, n)
    _, D = _ND(circle, u, v)
    dsigma = 2.0 * math.pi / n
    arclen = np.arccos(np.clip(circle @ (-u), -1.0, 1.0))
    rows, worst = [], 0.0
    for eps in epsilons:
        mass = float(np.sum(np.abs(D[arclen < eps]))) * dsigma / nq
        predicted = math.sin(gamma) * eps ** 2 / nq          # two-sided
        rows.append({"eps": eps, "excised_mass": mass,
                     "predicted_sin_g_eps2_over_q": predicted,
                     "ratio": mass / predicted})
        worst = max(worst, abs(mass / predicted - 1.0))
    return {"rows": rows, "worst_relative_error": worst,
            "mass_is_O_eps_squared": bool(worst < 0.01),
            "excision_is_safe": bool(worst < 0.01)}


def saddle_branch_ratio(kappa: float, u=(0.0, 0.0, 1.0),
                        v=(math.sin(1.0), 0.0, math.cos(1.0)),
                        convention: int = +1) -> Dict[str, object]:
    """``A_pi / A_0 = (M_pi/M_0) e^{2 i kappa} e^{-i pi/2 * convention}``.

    ``S_H = -1`` on the ``theta = 0`` component and ``+1`` on ``theta = pi``,
    and the transverse indices are ``0`` and ``1``, so the Maslov factors
    differ by one unit. The **magnitude** of the ratio is the mass ratio
    ``M_pi/M_0``, which is generically far from ``1`` — it is *not* unity, and
    an earlier version of this module wrongly said so by treating the local
    phase prefactor as the whole component amplitude.

    What survives, and is the actual content:

        arg(A_pi/A_0) = 2 kappa - pi/2 * convention.

    The masses are the already-derived positive coarea masses (see
    `morse_bott_component_masses`), so stationary phase supplies the magnitudes
    correctly and leaves **only the relative phase** open. The ratio is real —
    i.e. the aggregation is one of the two candidates rather than something
    complex — iff ``4 kappa/pi`` is an odd integer, a statement unchanged by
    flipping the convention, which shifts ``kappa`` by ``pi/2``; the sign then
    alternates with ``kappa`` mod ``pi``. A relative factor of ``-1`` gives
    ``M_0 - M_pi`` (the oriented sum), ``+1`` gives ``M_0 + M_pi`` (the positive
    count).
    """
    masses = morse_bott_component_masses(u, v, n=40001)
    phase = complex(np.exp(2j * kappa) * np.exp(-0.5j * math.pi * convention))
    ratio = (masses["M_pi"] / masses["M_0"]) * phase
    four_k = 4.0 * kappa / math.pi
    is_odd = abs(four_k - round(four_k)) < 1e-9 and int(round(four_k)) % 2 != 0
    return {
        "kappa": kappa, "convention": convention,
        "phase_factor": phase, "mass_ratio": masses["M_pi"] / masses["M_0"],
        "ratio": ratio, "magnitude": abs(ratio),
        "magnitude_is_the_mass_ratio": True,
        "4kappa_over_pi": four_k,
        "arg_phase_factor": float(np.angle(phase)),
        "phase_factor_is_real": bool(abs(phase.imag) < 1e-9),
        "odd_multiple_of_pi_over_4": bool(is_odd),
        "selects": ("oriented sum  M_0 - M_pi" if abs(phase + 1) < 1e-9 else
                    "positive count  M_0 + M_pi" if abs(phase - 1) < 1e-9 else
                    "neither: complex relative weight"),
    }


def no_off_closure_critical_points(trials: int = 12, n_theta: int = 400,
                                   seed: int = 21) -> Dict[str, object]:
    """``Crit(S_H) = Gamma_closure`` exactly — the off-closure half, proved.

    Checking ``grad S_H = 0`` *on* the closure set does not exclude critical
    points elsewhere, and round 6's ``grad theta|_{N=0} != 0`` concerns only
    ``N = 0``. The missing argument is short. With ``p = u + v``, ``q = u x v``,
    ``A = 1 + u.v``, so ``N = x.q`` and ``D = A + x.p``,

        grad_{S^2} theta = P_x (D q - N p) / (D^2 + N^2),

    so a regular critical point needs ``D q - N p`` parallel to ``x``. Since
    ``p . q = 0``, write ``x`` in the orthonormal basis ``(p_hat, q_hat, r_hat)``
    with ``r_hat = p_hat x q_hat``. The vector ``D q - N p`` has no ``r_hat``
    component, so either

    * ``x`` has a nonzero ``r_hat`` component, forcing ``D q - N p = 0``, hence
      ``D = N = 0`` — the excluded chart singularity ``x = -u, -v``; or
    * ``x`` lies in ``span(p, q)``, where the parallelism condition reduces to
      ``|p|(x_p^2 + x_q^2) + A x_p = 0``, i.e. ``x_p = -|p|/A``.

    And ``|p|/A = 2cos(gamma/2) / (2cos^2(gamma/2)) = sec(gamma/2) > 1`` for
    ``0 < gamma < pi``, while ``|x_p| <= 1``. No solution. Both the closed form
    and a global grid search are reported.
    """
    rng = np.random.default_rng(seed)
    worst_off = math.inf
    sec_worst = 0.0
    rows = []
    for _ in range(trials):
        u, v = _unit(rng.normal(size=3)), _unit(rng.normal(size=3))
        p, q = u + v, np.cross(u, v)
        A = 1.0 + float(u @ v)
        gamma = math.acos(float(np.clip(u @ v, -1.0, 1.0)))
        required = -float(np.linalg.norm(p)) / A
        sec_form = -1.0 / math.cos(0.5 * gamma)
        sec_worst = max(sec_worst, abs(required - sec_form))
        th = np.arccos(np.linspace(-1.0, 1.0, n_theta))
        ph = np.linspace(0.0, 2.0 * math.pi, 2 * n_theta, endpoint=False)
        T, P = np.meshgrid(th, ph, indexing="ij")
        X = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P),
                      np.cos(T)], axis=-1).reshape(-1, 3)
        N = X @ q
        D = A + X @ p
        amb = D[:, None] * q[None, :] - N[:, None] * p[None, :]
        tang = amb - np.sum(amb * X, axis=1)[:, None] * X
        grad = np.linalg.norm(tang, axis=1) / (D ** 2 + N ** 2)
        off = np.abs(N) > 1e-3
        worst_off = min(worst_off, float(np.min(grad[off])))
        rows.append({"gamma": gamma, "required_x_p": required,
                     "impossible": bool(abs(required) > 1.0)})
    return {
        "min_grad_theta_off_closure": worst_off,
        "required_x_p_equals_minus_sec_half_gamma": sec_worst,
        "every_required_x_p_exceeds_the_sphere": all(r["impossible"] for r in rows),
        "rows": rows,
        "no_off_closure_critical_points": bool(
            worst_off > 1e-3 and sec_worst < 1e-12
            and all(r["impossible"] for r in rows)),
    }


# ── B. the sector coefficients ──────────────────────────────────────────────

def sector_symmetry_group(a, b) -> Dict[str, object]:
    """The ``O(3)`` base-isometry subgroup of the fixed-setting problem.

    This is the group of *base isometries* fixing both setting axes — not a
    classification of every discrete symmetry of the boundary-value problem.
    `discrete_symmetry_extension` covers detector exchange, history reversal
    and the Pin deck; `fibre_action_is_weight_blind` covers the vertical
    ``Spin(2)``. See note N11.

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


def discrete_symmetry_extension(gamma: float = 1.0, n: int = 200001,
                                seed: int = 29) -> Dict[str, object]:
    """Beyond base isometries: do the other discrete operations mix like/unlike?

    `sector_symmetry_group` enumerates the ``O(3)`` operations fixing both
    setting axes, and `fibre_action_is_weight_blind` closes the vertical
    ``Spin(2)`` extension. Neither is a classification of *every* discrete
    symmetry of the fixed-setting boundary problem (note N11), so the three
    further operations the model actually has are checked here directly. Each
    is applied to the sector data and tested for whether it can send a *like*
    sector to an *unlike* one — the only thing that could force ``r = 1``.

    * **detector-label exchange** ``A <-> B``: ``(s_A, s_B) -> (s_B, s_A)``.
      Preserves ``s_A s_B``, so like stays like.
    * **history reversal** of the derived loop ``x -> u -> w -> x``: swaps
      ``u <-> w``. ``D = 1 + x.u + u.w + x.w`` is symmetric in ``u, w`` so every
      weight is invariant; ``N -> -N``, which is zero on the closure set.
    * **Pin deck / sector flip**: acts on the frame, not on ``x``, and the
      weights are fibre-blind, so it is the identity on sector weights.

    All three preserve ``s_A s_B``. Together with the base result this is an
    enumeration over the operations the model supplies, **not** a proof that no
    further symmetry exists — which is why the verdict is phrased as *no
    identified symmetry* rather than as a classification.
    """
    a = np.array([0.0, 0.0, 1.0])
    b = np.array([math.sin(gamma), 0.0, math.cos(gamma)])

    def weight(s_a, s_b):
        u, w = s_a * a, -(s_b * b)
        circle = great_circle(u, w, n)
        _, D = _ND(circle, u, w)
        return float(np.mean(D))

    sectors = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    base = {t: weight(*t) for t in sectors}

    def like(t):
        return t[0] * t[1] == 1

    rows = []
    # 1. detector-label exchange
    exch = {t: base[(t[1], t[0])] for t in sectors}
    rows.append({"operation": "detector-label exchange A<->B",
                 "preserves_s_A_s_B": all(like(t) == like((t[1], t[0]))
                                          for t in sectors),
                 "max_weight_change": max(abs(exch[t] - base[t])
                                          for t in sectors if like(t)
                                          == like((t[1], t[0])))})
    # 2. history reversal: u <-> w
    rev_worst = 0.0
    for s_a, s_b in sectors:
        u, w = s_a * a, -(s_b * b)
        circle = great_circle(u, w, n)
        _, D1 = _ND(circle, u, w)
        _, D2 = _ND(circle, w, u)
        rev_worst = max(rev_worst, abs(float(np.mean(D1)) - float(np.mean(D2))))
    rows.append({"operation": "history reversal (u <-> w)",
                 "preserves_s_A_s_B": True, "max_weight_change": rev_worst})
    # 3. Pin deck / fibre sector flip
    rows.append({"operation": "Pin deck / fibre sector flip",
                 "preserves_s_A_s_B": True,
                 "max_weight_change": 0.0})

    return {
        "rows": rows,
        "all_preserve_like_unlike": all(r["preserves_s_A_s_B"] for r in rows),
        "no_identified_operation_mixes_like_and_unlike": bool(
            all(r["preserves_s_A_s_B"] for r in rows)),
        "is_a_classification_of_all_symmetries": False,
    }


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
    """The six audited stress/flux quantities are degree-2 homogeneous.

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


def quadratic_readouts_disagree(n_angle: int = 20001) -> Dict[str, object]:
    """Degree-2 homogeneity does not pick a readout: two quadratics disagree.

    Round 6 gave the sector integral ``int_Gamma D_s dsigma = 2 pi (1 + u.w)``
    with ``u = s_A a``, ``w = -s_B b``, so ``c_s = u.w = -s_A s_B cos gamma``. A
    readout **linear** in that integral gives ``E = -cos gamma`` and
    ``S_max = 2 sqrt 2``. But "quadratic" is ambiguous, and the ambiguity is
    not academic — two perfectly ordinary classical quadratic operations give
    different physics:

    * **square the integral** (a coherent, amplitude-like response):
      ``R_s ~ (int_Gamma D_s dsigma)^2``, giving
      ``E = -2 cos g/(1 + cos^2 g)`` and ``S_max = 8 sqrt2/3 = 3.7712``;
    * **square locally, then integrate** (an ordinary intensity/energy
      response): ``R_s ~ int_Gamma D_s^2 dsigma``. On the closure circle
      ``<D_s^2> = (1 + c_s)(2 + c_s)`` exactly, since ``<x> = 0`` and
      ``<(x.p)^2> = |p|^2/2 = 1 + c_s``; this gives
      ``E = -3 cos g/(2 + cos^2 g)`` and ``S_max = 3.3941``.

    Both keep the marginals at exactly ``1/2``. Both are quadratic. They
    disagree, and nothing measured in `detector_response_homogeneity` chooses
    between them — which is why C cannot be reported as a derived quadratic
    law. It is reported as *no derived readout at all*.
    """
    def E_linear(g):
        return -math.cos(g)

    def E_square_of_integral(g):
        c = math.cos(g)
        return -2.0 * c / (1.0 + c * c)

    def E_integral_of_square(g):
        c = math.cos(g)
        return -3.0 * c / (2.0 + c * c)

    def chsh_max(E):
        grid = np.linspace(1e-6, math.pi / 2, n_angle)
        vals = [abs(3.0 * E(float(g)) - E(3.0 * float(g))) for g in grid]
        k = int(np.argmax(vals))
        return float(vals[k]), float(grid[k])

    def marginal_dev(weight):
        worst = 0.0
        for g in (0.3, 1.0, 2.0):
            c = math.cos(g)
            w = {(sa, sb): weight(-sa * sb * c)
                 for sa in (1, -1) for sb in (1, -1)}
            tot = sum(w.values())
            worst = max(worst,
                        abs(sum(x for k, x in w.items() if k[0] == 1) / tot - 0.5))
        return worst

    s_lin, g_lin = chsh_max(E_linear)
    s_sq_int, g_sq_int = chsh_max(E_square_of_integral)
    s_int_sq, g_int_sq = chsh_max(E_integral_of_square)
    return {
        "linear": {"E_at_1": E_linear(1.0), "S_max": s_lin, "argmax": g_lin},
        "square_of_integral": {"E_at_1": E_square_of_integral(1.0),
                               "S_max": s_sq_int, "argmax": g_sq_int,
                               "marginal_dev": marginal_dev(lambda c: (1 + c) ** 2)},
        "integral_of_square": {"E_at_1": E_integral_of_square(1.0),
                               "S_max": s_int_sq, "argmax": g_int_sq,
                               "marginal_dev": marginal_dev(
                                   lambda c: (1 + c) * (2 + c))},
        "tsirelson": 2.0 * math.sqrt(2.0),
        "closed_form_square_of_integral": 8.0 * math.sqrt(2.0) / 3.0,
        "closed_form_integral_of_square": 12.0 * math.sqrt(2.0) / 5.0,
        "the_two_quadratics_disagree": bool(abs(s_sq_int - s_int_sq) > 1e-3),
        "both_exceed_tsirelson": bool(
            s_sq_int > 2.0 * math.sqrt(2.0) and s_int_sq > 2.0 * math.sqrt(2.0)),
        "homogeneity_does_not_select_a_readout": True,
    }


def local_square_mean_is_closed_form(seed: int = 17) -> Dict[str, object]:
    """``<D_s^2>_Gamma = (1 + c_s)(2 + c_s)``, checked against the circle."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(12):
        g = float(rng.uniform(0.2, 2.8))
        a = np.array([0.0, 0.0, 1.0])
        b = np.array([math.sin(g), 0.0, math.cos(g)])
        for s_a in (1, -1):
            for s_b in (1, -1):
                u, w = s_a * a, -(s_b * b)
                c = float(u @ w)
                circle = great_circle(u, w, 200001)
                _, D = _ND(circle, u, w)
                worst = max(worst,
                            abs(float(np.mean(D ** 2)) - (1 + c) * (2 + c)))
    return {"worst_residual": worst, "closed_form_holds": bool(worst < 1e-5)}


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

    **Scope (note N9).** This is a statement about the **repository's existing
    radial action**, not a general theorem. It shows that
    ``closed_orbit_action`` plus the natural connection term does not generate
    ``S_H``; it does not show that no unified classical action could.
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
    """Is future-setting information *present* at the source? Yes. Readable? Open.

    Round 5 established ``rho(x|a,b) != rho(x)``. This gate was pre-registered
    as ``P(O_S|a,b) = P(O_S)`` for every source-local observable. What the
    computation actually supports is narrower than the first version of this
    module claimed, and the scope is corrected here (notes N7-N8):

    * **Established.** The outcome-summed conditioned density is *exactly*
      antipodally even, ``rho(-x|a,b) = rho(x|a,b)``, so every **odd** function
      of ``x`` has vanishing expectation. This does NOT make its full law
      setting-independent: its variance can change. The sign tests below
      are blind for these nondegenerate continuous projections. And
      **some even functions do** distinguish them: ``E[(x.m)^2|a,b]`` and
      ``E[|x.m| |a,b]`` both move with ``gamma``. For non-coplanar settings the
      conditioned supports are different great circles meeting in two points —
      mutually singular, total variation ``1``.
    * **Not established, and no longer claimed: "every even observable
      signals".** That is false, and two counterexamples are computed here:
      constants, and ``x.x = 1`` on the sphere, are even and perfectly blind.
      The correct statement is that *some* even functions separate the
      conditioned ensembles.
    * **Not established: that BAM's couplings are even in ``x``.** Degree-2
      homogeneity in a field ``phi`` (question C) is not the same operation as
      ``x -> -x`` on the mouth direction. Deducing one from the other needs an
      explicit map from the source variable ``x`` to scalar / Maxwell / GL
      field configurations, and the repository has none. The earlier version of
      this module made that inference; it is withdrawn.
    * **Not established: an operational channel.** The observables used here
      (``x.m``, ``(x.m)^2``, ``|x.m|``) are synthetic functions of ``x``, not
      couplings BAM possesses. Nothing shows a BAM apparatus can measure one
      without altering the global two-boundary solution. A prescribed final
      pointer value would add a boundary condition, but a final record can
      instead be an output. See ``docs/source_readout.md`` for an exact
      conditional Hamiltonian pointer preserving both source boundaries.
      Its physical frame-to-field map and intervention law remain missing.

    So the result is a **causality hazard**, not a demonstrated signalling
    channel: the conditioned source ensemble carries information about future
    settings that an ideal faithful readout of ``x`` would reveal. Closing it
    requires either a dynamical non-readability theorem or a reformulation in
    which ``x`` is not an operational degree of freedom — and the latter
    collides with round 5's use of ``x`` as a physical source variable.
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
        rows.append({
            "gamma": gamma,
            "odd  E[x.m]": float(np.sum(w * proj)),
            "odd  P(x.m>0)": float(np.sum(w[proj > 0])),
            "even E[(x.m)^2]": float(np.sum(w * proj * proj)),
            "even E[|x.m|]": float(np.sum(w * np.abs(proj))),
            "even E[1] (blind)": float(np.sum(w)),
            "even E[x.x] (blind)": float(np.sum(w * np.sum(circle ** 2, axis=1))),
        })

    def spread(key):
        return max(r[key] for r in rows) - min(r[key] for r in rows)

    odd_spread = max(spread("odd  E[x.m]"), spread("odd  P(x.m>0)"))
    some_even = max(spread("even E[(x.m)^2]"), spread("even E[|x.m|]"))
    blind_even = max(spread("even E[1] (blind)"), spread("even E[x.x] (blind)"))
    return {
        "non_coplanar_total_variation": non_coplanar.get("total_variation"),
        "non_coplanar_supports_mutually_singular": bool(
            float(non_coplanar.get("total_variation", 0.0)) > 0.99),
        "coplanar_total_variation": coplanar.get("total_variation"),
        "density_is_antipodally_even_residual": parity,
        "observable_rows": rows,
        "odd_observable_spread": odd_spread,
        "some_even_observables_separate": some_even,
        "blind_even_observable_spread": blind_even,
        "odd_means_and_signs_are_blind": bool(odd_spread < 1e-6),
        "odd_projection_variance_spread": spread("even E[(x.m)^2]"),
        "odd_full_distributions_are_blind": bool(spread("even E[(x.m)^2]") < 1e-10),
        "some_even_functions_separate_the_ensembles": bool(some_even > 1e-3),
        "not_every_even_observable_separates": bool(blind_even < 1e-12),
        # the claims deliberately NOT made:
        "bam_couplings_shown_even_in_x": False,
        "operational_readout_constructed": False,
        "setting_information_present_at_source": bool(
            some_even > 1e-3
            or float(non_coplanar.get("total_variation", 0.0)) > 0.99),
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
        {"input": "current-to-frequency readout", "status": "open",
         "where": ("question C: six audited stress/flux quantities are degree 2, but two "
                   "ordinary quadratics disagree and none is a derived detector "
                   "coupling")},
        {"input": "the round-5 Haar/uniform source measure (a = 1 in the "
                  "Morse-Bott coefficients)", "status": "chosen",
         "where": "note N12: a non-constant in-plane amplitude moves both masses"},
        {"input": "restriction to the two-parameter family of geodesic triangles",
         "status": "chosen",
         "where": "note N12: not stationary phase over BAM's full history space"},
        {"input": "a map from the source variable x to field configurations",
         "status": "open",
         "where": "gate E: needed before degree in phi says anything about parity in x"},
        {"input": "a source-local readout that respects the two-boundary problem",
         "status": "open", "where": "gate E: not constructed; measurement is itself a BC"},
        {"input": "antipodal scalar BC, eta, quotient-vs-cover",
         "status": "not used", "where": "C1: not in this strand"},
    ]


def verdicts() -> Dict[str, str]:
    """Five independent fields, computed from the measurements above.

    The pre-registered headline requires A, B, C and E together. It is not
    printed. C and E carry labels narrower than the pre-registered options
    because the pre-registered options overstated what the computation
    supports; both narrowings are recorded as correction notes N6 and N8.
    """
    a_ok = additive_functionals_have_no_critical_points()
    b_ok = sector_orbits()
    c_ok = detector_response_homogeneity()
    c_amb = quadratic_readouts_disagree()
    e_ok = source_observable_signalling(n=50000)

    A = ("HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION"
         if (a_ok["theta_additive"] and a_ok["S_H_not_additive"]
             and a_ok["theta_has_no_critical_point"])
         else "HISTORY_ACTION_DERIVED_BUT_BRANCH_PHASE_UNDERDETERMINED")
    B = ("PHYSICAL_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE"
         if b_ok["forced_at_any_chsh_angle"]
         else "NO_IDENTIFIED_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE")
    # Six audited stress/flux quantities are quadratic, but that does not name a
    # readout: two ordinary quadratic operations give different physics, and
    # none of them is derived from a detector coupling.
    C = ("NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_READOUT"
         if (c_ok["all_quadratic"] and c_amb["the_two_quadratics_disagree"])
         else "CLASSICAL_DETECTOR_DERIVES_LINEAR_HISTORY_READOUT"
         if c_ok["any_linear"] else "CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY")
    D = "HISTORY_ACTION_INDEPENDENTLY_POSTULATED"
    # information is present at the source; no operational channel is built.
    E = ("SETTING_INFORMATION_IS_PRESENT_AT_SOURCE_READOUT_DYNAMICS_OPEN"
         if (e_ok["setting_information_present_at_source"]
             and not e_ok["operational_readout_constructed"])
         else "SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS"
         if e_ok["operational_readout_constructed"]
         else "SOURCE_OBSERVABLES_OPERATIONALLY_NON_SIGNALLING")
    return {"A_action": A, "B_sectors": B, "C_readout": C,
            "D_compatibility": D, "E_causality": E,
            "headline": "NOT MET — no field of the pre-registered conjunction holds"}
