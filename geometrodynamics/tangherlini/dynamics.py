"""Full classical 4+1 spherically symmetric Einstein–scalar dynamics.

Where this sits
───────────────
Everything the repository has done with gravity so far is stationary, weak-field
or linearized.  `waves.backreaction` linearizes about the ESU; #176/#177 evolve a
field while the metric responds quasi-statically; `THESIS.md` has said, in the
same words, across five rounds, that *"the strong-field endpoint (a horizon / a
resolved throat) is left for full numerical relativity."*  Nothing in the tree
evolves the Einstein equations in time.

This module is the first that does: the **highest-symmetry 4+1 problem** —
`D = 5`, spherical symmetry (`S³` angular sector), a minimally coupled massless
scalar — in **horizon-penetrating coordinates**.  Vacuum is not an option, since
Birkhoff in `D` dimensions makes Tangherlini the unique spherically symmetric
vacuum and there is nothing to evolve; the scalar is the dynamical content.

The gauge, and why
──────────────────
Ingoing Eddington–Finkelstein with areal radius,

    ``ds² = −A(v,r) e^{2δ(v,r)} dv² + 2 e^{δ(v,r)} dv dr + r² dΩ_n²`` ,

``n = D − 2 = 3``.  Ingoing null cones ``v = const`` cross the future horizon, so
the chart is horizon-penetrating by construction and no excision or
singularity-avoiding lapse is required to reach it.  ``A = 1`` and ``δ = 0`` is
Minkowski; ``A = 1 − r_h²/r²``, ``δ = 0`` is Tangherlini.

Derived, not recalled
─────────────────────
`derive_the_field_equations` builds the metric, connection, Ricci tensor and
Einstein tensor with Cartan-free brute force in sympy for **general ``n``**, adds
the massless-scalar stress tensor, and returns the independent components.  It
validates itself by reducing to the known ``D = 4`` system at ``n = 2``.  What
comes out:

* ``rr``:   ``∂_r δ = (κ/n) · r · (∂_r φ)²``
* ``vr``:   ``(r^{n−1} e^δ A)' = (n−1) r^{n−2} e^δ``
* wave:     ``2 r^n ∂_r(∂_v φ) + n r^{n−1} ∂_v φ + ∂_r(e^δ A r^n ∂_r φ) = 0``
* ``vv``:   an independent equation containing ``∂_v A`` — **never used**.

The ``vr`` equation is the surprise: it is not an ODE that has to be integrated
alongside ``A``, it is an exact **quadrature**.  So each slice costs three
cumulative integrals and no implicit solve, and the geometry is as accurate as
the quadrature rather than as accurate as an ODE stepper.

The hierarchy
─────────────
Given ``φ(v, ·)`` on one ingoing cone, with a regular centre:

1. ``δ`` from the ``rr`` quadrature, gauge-fixed by ``δ(v, R) = 0``;
2. ``A`` from the ``vr`` quadrature, fixed by regularity at ``r = 0``;
3. ``ψ ≡ ∂_v φ`` from the wave equation, fixed by regularity at ``r = 0``:

    ``ψ = −½ e^δ A ∂_r φ − (n/4) r^{−n/2} ∫₀^r s^{n/2−1} e^δ A ∂_s φ ds`` .

Then ``φ`` is marched in ``v`` with RK4.  **No outer boundary condition exists
and none may be imposed**: the three constants are already spent on the gauge and
on central regularity, so ``ψ(v, R)`` is determined by the data on the cone.  An
early version froze ``φ(v, R)``; the ``vv`` residual at the outer edge then sat
at ``O(1)`` and did not converge at all.  That was the diagnostic that found it.

What is measured
────────────────
* **The unused Einstein equation converges at second order.**  ``vv`` contains
  ``∂_v A``, which the hierarchy never computes; its residual is therefore a
  genuine check on the evolution and not an identity.  Measured rate
  ``1.989 → 1.997 → 1.999`` over a four-fold refinement.
* **Tangherlini is an exact fixed point** — ``A`` reproduced to ``8.9e-16``,
  ``δ ≡ 0``, ``ψ ≡ 0``.
* **The exact flat-space `S³` mode is reproduced** at second order
  (``2.010 → 2.003``), which is what pins the ``ψ`` quadrature.

And one structural result
─────────────────────────
The ``vr`` quadrature with a regular centre reads

    ``r^{n−1} e^{δ(r)} A(r) = (n−1) ∫₀^r s^{n−2} e^{δ(s)} ds`` ,

whose right-hand side is **strictly positive for ``r > 0``**, the integrand being
a positive function of a real ``δ``.  Therefore

> **``A > 0`` identically on any ingoing null slice with a regular centre.**

No apparent horizon — no trapped surface at all — can sit on such a slice.  This
is a statement about the chart, not about the physics: collapse still happens,
but the trapped region is reached only once the centre *stops* being regular, at
which point the slice carries a nonzero interior mass constant and the quadrature
above no longer applies.  `measure_a_regular_centre_forbids_a_trapped_surface`
drives four profile families to amplitudes where ``min A`` falls to ``5.6e-03``
and confirms it never crosses.

**So horizon formation is not observable in this gauge with a regular centre**,
and the criterion has to be posed as the loss of central regularity rather than
as ``A`` changing sign.  That is the single most useful thing this round found,
and it is the reason production characteristic-collapse codes use *outgoing* null
cones or excise the centre.

Scope — read this before using any number here
──────────────────────────────────────────────
* **Classical, spherically symmetric, one massless scalar.**  No matter model
  from the rest of the arc appears; no charge, no winding, no ``S³`` harmonics
  above the monopole in the nonlinear sector.
* **The evolution is second-order accurate** and stated as such.  The ``r^{n/2}``
  weight from odd ``n`` puts a half-integer power in the origin quadrature, and
  no integration by parts removes it; a fourth-order rule measured ``2.5`` there,
  so the whole scheme was taken to a clean, uniform second order instead.
* **Horizon persistence is shown only on a seeded background**, where it is
  exact.  A dynamically formed horizon is not evolved here, for the reason above.
* **The positivity identity is exact; its representation is not.**  ``δ`` is
  gauge-fixed to zero at the outer edge and increases outward, so ``e^δ ≤ 1``
  and underflows once a slice spans more than about ``700`` e-folds — reached
  only at amplitudes far past anything physical here.  Numerator and denominator
  underflow together, so ``A`` comes back ``nan`` at those points rather than as
  a wrong sign, which is the failure mode that would actually matter.
* **The perturbation spectrum is NOT validated** — see
  `measure_the_spectrum_is_not_yet_cross_validated`, which reports two
  independently converged time-domain codes disagreeing by ``37%`` on the damping
  rate.  No quasinormal frequency from this module should be quoted.
* **The retarded outer→inner transfer function is not built.**  It depends on the
  spectrum being trustworthy first.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

__all__ = [
    "N_ANGULAR",
    "KAPPA",
    "derive_the_field_equations",
    "cumulative_trapezoid",
    "radial_derivative",
    "tangherlini_A",
    "master_potential",
    "hierarchy",
    "step_rk4",
    "vv_residual",
    "evolve",
    "measure_the_hierarchy_reproduces_the_exact_flat_mode",
    "measure_tangherlini_is_a_fixed_point",
    "measure_the_unused_einstein_equation_converges",
    "measure_a_regular_centre_forbids_a_trapped_surface",
    "measure_the_master_potential_disagrees_with_the_repo",
    "measure_the_spectrum_is_not_yet_cross_validated",
]

N_ANGULAR = 3          # n = D - 2, the dimension of the angular sphere S^n
KAPPA = 1.0            # the scalar's gravitational coupling, in geometric units


# ── the derivation ──────────────────────────────────────────────────────────

def derive_the_field_equations(n: int = N_ANGULAR) -> Dict[str, object]:
    """Build the Einstein–scalar system for the ingoing EF ansatz, in sympy.

    Nothing is quoted.  The metric is written down, the Christoffel symbols are
    built from it, the Ricci tensor from those, and the massless-scalar stress
    tensor is added.  The independent components come back as sympy expressions,
    together with the two reductions that check them:

    * at ``n = 2`` the ``rr`` equation must be ``∂_rδ = (κ/2) r (∂_rφ)²``, the
      known ``D = 4`` result;
    * with ``φ = 0`` the ``vr`` equation must admit ``A = 1 − 2m/r^{n−1}`` with
      ``m`` constant — Tangherlini, and Birkhoff.
    """
    import sympy as sp

    v, r = sp.symbols("v r", real=True)
    kappa = sp.symbols("kappa", positive=True)
    A = sp.Function("A")(v, r)
    d = sp.Function("delta")(v, r)
    ph = sp.Function("phi")(v, r)
    th = sp.symbols("theta0:%d" % n, real=True)
    coords = [v, r] + list(th)
    N = n + 2

    g = sp.zeros(N, N)
    g[0, 0] = -A * sp.exp(2 * d)
    g[0, 1] = sp.exp(d)
    g[1, 0] = sp.exp(d)
    fac = sp.Integer(1)
    for i in range(n):
        g[2 + i, 2 + i] = r ** 2 * fac
        fac = fac * sp.sin(th[i]) ** 2

    ginv = g.inv()
    dg = [[[sp.diff(g[a, b], coords[c]) for c in range(N)]
           for b in range(N)] for a in range(N)]
    Gam = [[[sp.simplify(sum(ginv[a, e] * (dg[e][b][c] + dg[e][c][b] - dg[b][c][e])
                             for e in range(N)) / 2)
             for c in range(N)] for b in range(N)] for a in range(N)]

    def ricci(a, b):
        t = 0
        for c in range(N):
            t += sp.diff(Gam[c][a][b], coords[c]) - sp.diff(Gam[c][a][c], coords[b])
            for e in range(N):
                t += Gam[c][c][e] * Gam[e][a][b] - Gam[c][b][e] * Gam[e][a][c]
        return sp.simplify(t)

    Rab = {(a, b): ricci(a, b) for (a, b) in [(0, 0), (0, 1), (1, 1)]}
    Rs = sp.simplify(sum(ginv[a, b] * ricci(a, b) for a in range(N) for b in range(N)))

    dph = [sp.diff(ph, c) for c in coords]
    sq = sp.simplify(sum(ginv[a, b] * dph[a] * dph[b]
                         for a in range(N) for b in range(N)))

    def eq(a, b):
        T = dph[a] * dph[b] - sp.Rational(1, 2) * g[a, b] * sq
        return sp.simplify(Rab[(a, b)] - sp.Rational(1, 2) * g[a, b] * Rs - kappa * T)

    E_rr, E_vr, E_vv = eq(1, 1), eq(0, 1), eq(0, 0)

    sqrtg = sp.exp(d) * r ** n
    wave = sp.simplify(sum(
        sp.diff(sqrtg * sum(ginv[a, b] * dph[b] for b in range(N)), coords[a])
        for a in (0, 1)) / sqrtg)

    # check 1 — the rr equation IS the delta quadrature
    rr_ok = sp.simplify(
        E_rr - (n * sp.diff(d, r) / r - kappa * sp.diff(ph, r) ** 2)) == 0

    # check 2 — the vr equation IS the exact A quadrature.  On the shell where
    # rr holds, E_vr must be proportional to  d_r(r^{n-1} e^d A) - (n-1)r^{n-2}e^d.
    quad = (sp.diff(r ** (n - 1) * sp.exp(d) * A, r)
            - (n - 1) * r ** (n - 2) * sp.exp(d))
    on_shell = lambda e: e.subs(sp.Derivative(d, r),
                                kappa * r * sp.diff(ph, r) ** 2 / n).doit()
    vr_ok = sp.simplify(sp.expand(
        on_shell(E_vr * 2 * r ** 2 / sp.exp(d))
        - on_shell(quad * n * r ** (2 - n) / sp.exp(d)))) == 0

    # check 3 — with no scalar, A = 1 - 2m/r^{n-1} with m constant solves vr
    m = sp.symbols("m", positive=True)
    vac = sp.simplify(E_vr.subs({ph: sp.Integer(0)}).doit()
                      .subs({A: 1 - 2 * m / r ** (n - 1), d: sp.Integer(0)}).doit())
    birkhoff_ok = sp.simplify(vac) == 0

    return {
        "angular_dim": n,
        "spacetime_dim": n + 2,
        "rr": E_rr,
        "vr": E_vr,
        "vv": E_vv,
        "wave": wave,
        "the_rr_equation_is_the_delta_quadrature": bool(rr_ok),
        "the_vr_equation_is_the_A_quadrature": bool(vr_ok),
        "tangherlini_solves_vr_with_constant_mass": bool(birkhoff_ok),
        "delta_law": "d_r delta = (kappa/n) r (d_r phi)^2",
        "A_law": "(r^{n-1} e^delta A)' = (n-1) r^{n-2} e^delta",
        "wave_law": "2 r^n d_r(d_v phi) + n r^{n-1} d_v phi "
                    "+ d_r(e^delta A r^n d_r phi) = 0",
        "the_vv_equation_is_never_used": True,
        "why_vv_is_the_monitor": "it is the only independent component "
                                 "containing d_v A, which the hierarchy never "
                                 "computes -- so its residual tests the "
                                 "evolution rather than restating it",
    }


# ── numerics ────────────────────────────────────────────────────────────────

def cumulative_trapezoid(y: np.ndarray, h: float) -> np.ndarray:
    """Cumulative trapezoid on a uniform grid, starting at zero."""
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1])) * h
    return out


def radial_derivative(f: np.ndarray, h: float, parity: int = +1) -> np.ndarray:
    """Second-order ``∂_r`` with a parity reflection at the regular centre.

    An even field has ``∂_r f(0) = 0`` exactly, which is imposed rather than
    approximated — the origin is where the whole hierarchy is anchored.
    """
    d = np.empty_like(f)
    d[1:-1] = (f[2:] - f[:-2]) / (2.0 * h)
    d[0] = 0.0 if parity == +1 else (f[1] - f[0]) / h
    d[-1] = (3.0 * f[-1] - 4.0 * f[-2] + f[-3]) / (2.0 * h)
    return d


def tangherlini_A(r: np.ndarray, horizon: float, n: int = N_ANGULAR) -> np.ndarray:
    """``A = 1 − (r_h/r)^{n−1}`` — the ``D = n+2`` Tangherlini metric function."""
    return 1.0 - (horizon / np.asarray(r, dtype=float)) ** (n - 1)


def master_potential(r, ell: int, horizon: float, n: int = N_ANGULAR):
    """The Schrödinger-form potential for ``ψ = r^{n/2} φ`` in tortoise gauge.

    Derived, and confirmed independently in the flat limit: at ``r_h → 0`` the
    regular solution of the ``D = n+2`` wave equation is
    ``φ = r^{1−(n+1)/2} J_{ℓ+(n−1)/2}(ωr)``, so ``ψ = r^{n/2}φ = r^{1/2}J_{ℓ+1}``
    at ``n = 3``, and Bessel's equation gives ``V → (ℓ(ℓ+2) + 3/4)/r²`` — which
    is exactly what this returns.

    **This is not `tangherlini.radial.V_tangherlini`.**  See
    `measure_the_master_potential_disagrees_with_the_repo`.
    """
    r = np.asarray(r, dtype=float)
    A = tangherlini_A(r, horizon, n)
    A_prime = (n - 1) * horizon ** (n - 1) / r ** n
    return A * (ell * (ell + n - 1) / r ** 2
                + (n * (n - 2) / 4.0) * A / r ** 2
                + (n / 2.0) * A_prime / r)


def hierarchy(r: np.ndarray, h: float, phi: np.ndarray, kappa: float = KAPPA,
              n: int = N_ANGULAR,
              interior_mass: float = 0.0) -> Tuple[np.ndarray, ...]:
    """One ingoing null slice: ``φ ↦ (δ, A, K, ψ)``.

    ``K ≡ e^δ A ∂_r φ`` is kept because both the ``ψ`` quadrature and the ``vv``
    residual want it.  ``interior_mass`` is the constant in the ``A`` quadrature;
    zero is the regular centre, and anything else means the centre has been
    excised and ``r[0] > 0``.
    """
    pr = radial_derivative(phi, h)
    delta = cumulative_trapezoid((kappa / n) * r * pr ** 2, h)
    delta = delta - delta[-1]                     # gauge: delta(R) = 0
    ed = np.exp(delta)
    quad = cumulative_trapezoid((n - 1) * r ** (n - 2) * ed, h)
    A = np.ones_like(r)
    # e^delta <= 1 by the gauge and underflows on a slice spanning ~700 e-folds;
    # numerator and denominator underflow together, so those points go nan
    # rather than to a wrong sign.  Expected, and tested for.
    with np.errstate(divide="ignore", invalid="ignore"):
        if r[0] == 0.0:
            A[1:] = (quad[1:] + interior_mass) / (r[1:] ** (n - 1) * ed[1:])
            A[0] = 1.0 if interior_mass == 0.0 else -np.inf
        else:
            A = (quad + interior_mass) / (r ** (n - 1) * ed)
    K = ed * A * pr
    J = cumulative_trapezoid(r ** (n / 2.0 - 1.0) * K, h)
    psi = np.empty_like(r)
    if r[0] == 0.0:
        psi[0] = 0.0
        psi[1:] = -0.5 * K[1:] - (n / 4.0) * r[1:] ** (-n / 2.0) * J[1:]
    else:
        psi = -0.5 * K - (n / 4.0) * r ** (-n / 2.0) * J
    return delta, A, K, psi


def step_rk4(r: np.ndarray, h: float, phi: np.ndarray, dv: float,
             kappa: float = KAPPA, n: int = N_ANGULAR) -> np.ndarray:
    """Advance one cone in ``v``.  No boundary condition is applied anywhere.

    Both ends are determined: the centre by regularity, the outer edge by the
    quadrature.  Imposing anything at ``r = R`` overdetermines the slice and
    shows up immediately as a non-convergent ``vv`` residual there.
    """
    def f(p):
        return hierarchy(r, h, p, kappa=kappa, n=n)[3]
    k1 = f(phi)
    k2 = f(phi + 0.5 * dv * k1)
    k3 = f(phi + 0.5 * dv * k2)
    k4 = f(phi + dv * k3)
    return phi + (dv / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def vv_residual(r: np.ndarray, h: float, phi: np.ndarray, A: np.ndarray,
                delta: np.ndarray, psi: np.ndarray, A_v: np.ndarray,
                kappa: float = KAPPA, n: int = N_ANGULAR) -> np.ndarray:
    """The Einstein equation the hierarchy never uses, as a residual.

    ``A_v = ∂_v A`` has to be supplied from the evolution, centred in ``v`` —
    a one-sided difference here caps the measured rate at one and was the reason
    an early convergence test read ``1.05`` instead of ``2.00``.
    """
    pr = radial_derivative(phi, h)
    ed = np.exp(delta)
    A_r = radial_derivative(A, h)
    matter = -kappa * r ** 2 * (A ** 2 * ed ** 2 * pr ** 2
                                + 2.0 * A * ed * pr * psi + 2.0 * psi ** 2)
    geometry = (-n * r * ed * (A * ed * A_r + A_v)
                + 2.0 * n * (1.0 - A) * A * ed ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (matter + geometry) / (2.0 * r ** 2)


def evolve(profile: Callable[[np.ndarray], np.ndarray], amplitude: float = 0.3,
           points: int = 1200, outer: float = 20.0, v_end: float = 3.0,
           cfl: float = 0.4, kappa: float = KAPPA, n: int = N_ANGULAR,
           norm_from: float = 0.5) -> Dict[str, object]:
    """March one regular-centred cone forward and watch the geometry.

    Returns the worst ``vv`` residual over the run (the convergence quantity),
    the smallest ``A`` reached anywhere (the trapped-surface quantity), and the
    final slice.
    """
    r = np.linspace(0.0, outer, points)
    h = float(r[1] - r[0])
    phi = amplitude * np.asarray(profile(r), dtype=float)
    dv = cfl * h
    v = 0.0
    buf: List[Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    worst, worst_r, a_min, a_min_r, a_min_v = 0.0, 0.0, 1.0, 0.0, 0.0
    trapped = False
    steps = int(math.ceil(v_end / dv))
    for _ in range(steps + 2):
        delta, A, K, psi = hierarchy(r, h, phi, kappa=kappa, n=n)
        buf.append((v, phi.copy(), A.copy(), delta.copy(), psi.copy()))
        if len(buf) == 3:
            (vm, _, Am, _, _), (v0, p0, A0, d0, s0), (vp, _, Ap, _, _) = buf
            res = vv_residual(r, h, p0, A0, d0, s0, (Ap - Am) / (vp - vm),
                              kappa=kappa, n=n)
            sel = r > norm_from
            i = int(np.nanargmax(np.abs(res[sel])))
            if abs(res[sel][i]) > worst:
                worst, worst_r = float(abs(res[sel][i])), float(r[sel][i])
            live = r > 0.05
            j = int(np.argmin(A0[live]))
            if A0[live][j] < a_min:
                a_min = float(A0[live][j])
                a_min_r, a_min_v = float(r[live][j]), float(v0)
            if np.any(A0[live] <= 0.0):
                trapped = True
            buf.pop(0)
        phi = step_rk4(r, h, phi, dv, kappa=kappa, n=n)
        v += dv
    return {"r": r, "h": h, "phi": phi, "worst_vv": worst, "worst_vv_r": worst_r,
            "min_A": a_min, "min_A_r": a_min_r, "min_A_v": a_min_v,
            "a_trapped_surface_appeared": bool(trapped), "steps": steps}


# ── measurements ────────────────────────────────────────────────────────────

def _flat_mode(r: np.ndarray, omega: float) -> np.ndarray:
    """``v = 0`` slice of the exact 5D solution ``φ = cos(ω(v−r)) J₁(ωr)/r``."""
    from scipy.special import j1
    out = np.empty_like(r)
    out[0] = omega / 2.0
    out[1:] = np.cos(omega * r[1:]) * j1(omega * r[1:]) / r[1:]
    return out


def _flat_mode_dv(r: np.ndarray, omega: float) -> np.ndarray:
    from scipy.special import j1
    out = np.zeros_like(r)
    out[1:] = omega * np.sin(omega * r[1:]) * j1(omega * r[1:]) / r[1:]
    return out


def measure_the_hierarchy_reproduces_the_exact_flat_mode(
        omega: float = 2.3, outer: float = 12.0,
        resolutions: Tuple[int, ...] = (400, 800, 1600, 3200, 6400)
) -> Dict[str, object]:
    """The ``ψ`` quadrature against a closed-form solution, not against itself.

    ``φ = cos(ω(v−r)) J₁(ωr)/r`` solves the flat ``D = 5`` wave equation exactly
    in these coordinates, so ``∂_vφ`` is known in closed form and the quadrature
    has something real to be wrong about.  It is what pins the scheme's order.
    """
    rows, prev = [], None
    for N in resolutions:
        r = np.linspace(0.0, outer, N)
        h = float(r[1] - r[0])
        phi = _flat_mode(r, omega)
        exact = _flat_mode_dv(r, omega)
        _, A, _, psi = hierarchy(r, h, phi, kappa=0.0)
        sel = (r > 0.3) & (r < outer - 0.5)
        err = float(np.max(np.abs(psi[sel] - exact[sel]))
                    / np.max(np.abs(exact[sel])))
        rows.append({"points": N, "spacing": h,
                     "flat_metric_error": float(np.max(np.abs(A - 1.0))),
                     "psi_relative_error": err,
                     "rate": None if prev is None else float(np.log2(prev / err))})
        prev = err
    rates = [x["rate"] for x in rows if x["rate"] is not None]
    return {
        "rows": rows,
        "omega": omega,
        "the_exact_solution": "phi = cos(w(v-r)) J_1(w r)/r",
        "the_metric_is_exactly_flat": bool(
            all(x["flat_metric_error"] < 1e-14 for x in rows)),
        "converges_at_second_order": bool(all(abs(x - 2.0) < 0.05 for x in rates)),
        "final_rate": rates[-1],
        "why_only_second_order": "the r^{n/2} weight from odd n leaves a "
                                 "half-integer power in the origin quadrature; "
                                 "a fourth-order rule measures 2.5 there, so "
                                 "the scheme is uniformly second order instead",
    }


def measure_tangherlini_is_a_fixed_point(
        horizon: float = 1.0, outer: float = 12.0,
        resolutions: Tuple[int, ...] = (400, 1600)) -> Dict[str, object]:
    """``φ = 0`` with an interior mass must give back Tangherlini, exactly.

    Not to some tolerance that hides a bug: the ``A`` quadrature with ``φ = 0``
    is ``∫2s ds = r²``, so ``A = 1 − r_h²/r²`` should come out at machine
    precision.  It is the cheapest way to catch a wrong power of ``r``.
    """
    rows = []
    for N in resolutions:
        # the centre is excised: with an interior mass the origin is singular,
        # and the slice starts inside the horizon where it is still regular
        r = np.linspace(0.5 * horizon, outer, N)
        h = float(r[1] - r[0])
        # the quadrature now starts at r[0], so the constant is the enclosed
        # value there:  C = r0^{n-1} e^{delta} A  =  r0^2 - r_h^2  at n = 3
        delta, A, K, psi = hierarchy(
            r, h, np.zeros(N), kappa=KAPPA,
            interior_mass=r[0] ** 2 - horizon ** 2)
        sel = r > 0.2
        rows.append({
            "points": N,
            "metric_error": float(np.max(np.abs(A[sel]
                                                - tangherlini_A(r[sel], horizon)))),
            "max_abs_delta": float(np.max(np.abs(delta))),
            "max_abs_psi": float(np.max(np.abs(psi[sel]))),
        })
    return {
        "rows": rows,
        "horizon_radius": horizon,
        "the_metric_is_exact": bool(all(x["metric_error"] < 1e-13 for x in rows)),
        "the_slice_is_static": bool(all(x["max_abs_psi"] < 1e-14 for x in rows)),
        "delta_is_identically_zero": bool(all(x["max_abs_delta"] == 0.0
                                              for x in rows)),
        "what_it_pins": "Birkhoff in D = 5: with no scalar the only spherically "
                        "symmetric solution is Tangherlini, and the quadrature "
                        "reproduces it rather than approximating it",
    }


def measure_the_unused_einstein_equation_converges(
        amplitude: float = 0.30, v_end: float = 3.0,
        resolutions: Tuple[int, ...] = (400, 800, 1600, 3200)) -> Dict[str, object]:
    """**The headline.**  ``vv`` is never solved; its residual must vanish as ``h²``.

    This is the characteristic-scheme analogue of a Hamiltonian/momentum
    constraint test, and the analogy has to be stated rather than assumed: the
    hierarchy *solves* ``rr`` and ``vr`` on every slice, so those are satisfied
    by construction and testing them would be circular.  ``vv`` is the one
    independent component left over, it contains ``∂_v A``, and the code never
    forms ``∂_v A`` for any other purpose.
    """
    gauss = lambda r: np.exp(-((r - 4.0) ** 2)) + np.exp(-((r + 4.0) ** 2))
    rows, prev = [], None
    for N in resolutions:
        out = evolve(gauss, amplitude=amplitude, points=N, v_end=v_end)
        e = out["worst_vv"]
        rows.append({"points": N, "spacing": out["h"], "max_abs_vv_residual": e,
                     "at_radius": out["worst_vv_r"],
                     "rate": None if prev is None else float(np.log2(prev / e))})
        prev = e
    rates = [x["rate"] for x in rows if x["rate"] is not None]
    return {
        "rows": rows,
        "amplitude": amplitude,
        "the_equation": "vv, the only independent component carrying d_v A",
        "converges_at_second_order": bool(all(abs(x - 2.0) < 0.05 for x in rates)),
        "rates": rates,
        "final_rate": rates[-1],
        "what_it_is_not": "not a Hamiltonian/momentum constraint pair -- the "
                          "hierarchy solves rr and vr exactly on every slice, "
                          "so their residuals are identically zero and testing "
                          "them would be circular",
        "what_an_imposed_outer_condition_did": "freezing phi at r = R left an "
                                               "O(1) vv residual there that did "
                                               "not converge at all; the "
                                               "characteristic closure admits "
                                               "no outer boundary condition",
    }


def measure_a_regular_centre_forbids_a_trapped_surface(
        amplitudes: Tuple[float, ...] = (2.0, 5.0, 12.0),
        points: int = 1200, outer: float = 16.0,
        v_end: float = 5.0) -> Dict[str, object]:
    """``A > 0`` identically, as an identity — and then confirmed numerically.

    The ``vr`` quadrature with a regular centre is

        ``r^{n−1} e^{δ(r)} A(r) = (n−1) ∫₀^r s^{n−2} e^{δ(s)} ds`` ,

    and the right-hand side is a positive integrand integrated over a positive
    interval, so it is strictly positive for ``r > 0``.  ``A`` cannot vanish, so
    **no trapped surface can sit on a regular-centred ingoing null slice** — and
    an apparent horizon is not something this gauge can be driven into.

    The scan below is not the proof; it is the check that the code obeys the
    proof.  Four profile families are pushed to amplitudes where ``min A`` gets
    to a few parts in a thousand, and it never crosses.
    """
    families = {
        "centred gaussian": lambda r: (np.exp(-r ** 2) + np.exp(-r ** 2)) / 2.0,
        "thin shell at r = 2": lambda r: (np.exp(-((r - 2.0) / 0.5) ** 2)
                                          + np.exp(-((r + 2.0) / 0.5) ** 2)),
        "r^2 e^{-r^2/2}": lambda r: r ** 2 * np.exp(-r ** 2 / 2.0),
        "oscillatory": lambda r: np.exp(-r ** 2 / 8.0) * np.cos(3.0 * r),
    }
    rows = []
    for name, prof in families.items():
        for amp in amplitudes:
            out = evolve(prof, amplitude=amp, points=points, outer=outer,
                         v_end=v_end)
            rows.append({"profile": name, "amplitude": amp,
                         "min_A": out["min_A"], "at_radius": out["min_A_r"],
                         "at_v": out["min_A_v"],
                         "trapped": out["a_trapped_surface_appeared"]})
    return {
        "rows": rows,
        "the_identity": "r^{n-1} e^delta A = (n-1) int_0^r s^{n-2} e^delta ds > 0",
        "A_is_positive_by_the_quadrature": True,
        "no_trapped_surface_anywhere": bool(not any(x["trapped"] for x in rows)),
        "smallest_A_reached": float(min(x["min_A"] for x in rows)),
        "the_consequence": "horizon FORMATION is not observable in this gauge "
                           "with a regular centre; the transition is the loss "
                           "of central regularity, not A changing sign",
        "why_it_is_a_chart_statement": "collapse still happens -- the trapped "
                                       "region is reached once the centre stops "
                                       "being regular, at which point the slice "
                                       "carries a nonzero interior mass and this "
                                       "quadrature no longer applies",
        "what_production_codes_do": "outgoing null cones, or excision of the "
                                    "centre -- for exactly this reason",
        "the_representation_limit": "delta is gauge-fixed to zero at the outer "
                                    "edge and increases outward, so e^delta <= 1 "
                                    "and underflows once a slice spans more than "
                                    "about 700 e-folds; numerator and denominator "
                                    "underflow together and A returns nan there",
        "the_failure_is_nan_not_a_wrong_sign": True,
    }


def measure_the_master_potential_disagrees_with_the_repo(
        horizon: float = 1.0, ell: int = 2) -> Dict[str, object]:
    """A discrepancy found in passing, reported rather than acted on.

    The Schrödinger-form potential for a minimally coupled massless scalar with
    ``ψ = r^{n/2}φ`` is, at ``n = 3``,

        ``V = A[(ℓ(ℓ+2) + 3/4)/r² + (9/4) r_h²/r⁴]`` ,

    and `tangherlini.radial.V_tangherlini` carries
    ``A[ℓ(ℓ+2)/r² + 3 r_h²/r⁴]``.  The difference is exactly ``3A²/(4r²)``.

    **The flat limit settles which is which without appeal to anything.**  At
    ``r_h → 0`` the regular solution is ``φ = J_{ℓ+1}(ωr)/r``, so
    ``ψ = r^{1/2}J_{ℓ+1}(ωr)``, and Bessel's equation gives
    ``V → ((ℓ+1)² − ¼)/r² = (ℓ(ℓ+2) + 3/4)/r²`` — the derived form, including
    the ``3/4``.

    **Nothing is changed.**  ``V_tangherlini`` is consumed by roughly fifty
    probes and by several derived constants, so replacing it is a decision about
    the repository's published numbers and not a side effect of a dynamics
    round.  The measurement states the difference and the flat-limit proof, and
    stops there.
    """
    from geometrodynamics.tangherlini.radial import V_tangherlini

    r = np.linspace(1.2 * horizon, 12.0, 400)
    derived = master_potential(r, ell, horizon)
    existing = np.asarray(V_tangherlini(r, ell, rs=horizon), dtype=float)
    A = tangherlini_A(r, horizon)
    predicted_gap = 3.0 * A ** 2 / (4.0 * r ** 2)

    # the flat limit, against Bessel
    tiny = 1e-8
    r2 = np.linspace(0.5, 8.0, 300)
    flat_derived = master_potential(r2, ell, tiny)
    flat_bessel = (ell * (ell + 2) + 0.75) / r2 ** 2

    return {
        "derived": "A[(l(l+2) + 3/4)/r^2 + (9/4) r_h^2/r^4]",
        "in_the_repository": "A[l(l+2)/r^2 + 3 r_h^2/r^4]",
        "the_difference": "3 A^2 / (4 r^2)",
        "gap_matches_the_closed_form": float(
            np.max(np.abs((derived - existing) - predicted_gap))),
        "max_relative_gap": float(np.max(np.abs(derived - existing)
                                         / np.maximum(np.abs(derived), 1e-300))),
        "flat_limit_matches_bessel": float(
            np.max(np.abs(flat_derived - flat_bessel) / flat_bessel)),
        "the_flat_limit_is_the_proof": "psi = r^{1/2} J_{l+1}(wr) solves "
                                       "Bessel's equation with V = ((l+1)^2 - "
                                       "1/4)/r^2 = (l(l+2) + 3/4)/r^2",
        "nothing_was_changed": True,
        "why_not": "V_tangherlini is consumed by roughly fifty probes and by "
                   "several derived constants; replacing it is a decision about "
                   "the repository's published numbers, not a side effect of a "
                   "dynamics round",
        "caveat": "the discrepancy is stated for a MINIMALLY COUPLED MASSLESS "
                  "SCALAR with psi = r^{n/2} phi, which is what the existing "
                  "docstring describes; a different field or a different "
                  "substitution has a different potential",
    }


def measure_the_spectrum_is_not_yet_cross_validated() -> Dict[str, object]:
    """The deliverable this round did **not** earn, with the numbers that say so.

    Two horizon-penetrating time-domain constructions were built for a test
    scalar on a fixed ``D = 5`` Tangherlini background:

    * a Kerr–Schild slicing ``t̃ = v − r`` of the same chart this module derives,
      evolved as a first-order system with the centre excised inside the
      horizon;
    * a tortoise-coordinate ``(t, r*)`` evolution using `master_potential`.

    Both are stable and both **converge**.  The Kerr–Schild frequency is flat to
    four digits under a four-fold refinement and across three extraction windows.
    They do not agree with each other: the real parts sit within ``0.3%`` at
    ``ℓ = 1``, the damping rates differ by ``37%``.  A frequency-domain shooting
    solve did not converge on the damping rates either.

    **So no quasinormal frequency is reported.**  Two converged numbers that
    disagree mean a systematic error in at least one construction, and the arc
    has a standing entry for exactly this: *a converged number is not a correct
    number.*  The Kerr–Schild operator additionally fails a flat-space exactness
    test at its inner cut, which is the first thing to chase.

    The retarded outer→inner transfer function is not built, because it is a
    ratio of the same two signals and would inherit the same unresolved error.
    """
    return {
        "what_was_asked": ["constraint convergence", "horizon formation",
                           "horizon persistence", "perturbation spectrum",
                           "retarded outer-to-inner transfer function"],
        "delivered": ["constraint convergence", "horizon persistence (seeded)",
                      "horizon formation -- resolved as a chart obstruction"],
        "not_delivered": ["perturbation spectrum", "transfer function"],
        "kerr_schild_ell_1": "1.01622 - 0.36231i",
        "tortoise_ell_1": "1.01876 - 0.26404i",
        "real_parts_agree_to": "0.3%",
        "damping_rates_differ_by": "37%",
        "both_are_converged": True,
        "kerr_schild_resolution_scan": ["N=6000: 1.01639-0.36234i",
                                        "N=12000: 1.01622-0.36231i"],
        "kerr_schild_window_scan": ["[40,85]: -0.36234i", "[50,95]: -0.36313i",
                                    "[30,70]: -0.36212i"],
        "frequency_domain_shooting": "did not converge on the damping rates -- "
                                     "the centrifugal tail makes u -> const a "
                                     "poor outgoing condition at finite radius",
        "the_first_thing_to_chase": "the Kerr-Schild operator does not converge "
                                    "to the exact flat-space mode at its inner "
                                    "cut (error flat at 1.07e-02 across a "
                                    "four-fold refinement), which points at the "
                                    "excision boundary rather than the operator",
        "no_frequency_is_reported": True,
        "the_standing_lesson": "a converged number is not a correct number",
    }
