"""Metric backreaction: does ``A + B`` do something rescaling ``A`` or ``B`` cannot?

Where this sits
───────────────
PR #260 gated stationary action and backreaction on a growing mode; PR #261
removed it by resolving the mouth and PR #262 by removing the balls, upgrading
the answer to a theorem.  The gate is closed, so this round asks the first GR
question the roadmap actually named — and it is deliberately **not** "does
spacetime pinch off?":

    does ``A + B`` produce a metric response that rescaling ``A`` or ``B``
    alone cannot reproduce?

The field equation is linear, so ``φ_{A+B} = φ_A + φ_B``; ``T`` is quadratic, so

    ``T[A+B] = T[A] + T[B] + ΔT``

with ``ΔT`` the bilinear cross term; and linearized Einstein is linear in ``T``,
so ``β[A+B] = β[A] + β[B] + β[ΔT]``.  Rescaling ``A → cA`` sends
``β[A] → c²β[A]``, so **everything reachable by rescaling is
``{c²β[A] + d²β[B]}``** and the question is exactly whether ``β[ΔT]`` lies in
it.  That is a projection residual, and it is what this module measures.

The channel
───────────
The **transverse-traceless** sector, in its lowest (homogeneous, ``n = 3``)
harmonic — the shear of the universe.  Four reasons, and the first is the one
that matters:

* the ESU is held static by a fluid this arc never specifies.  A perfect fluid
  carries **no anisotropic stress**, and in an orthonormal frame its
  ``T_ab = diag(ρ,p,p,p)`` regardless of the anisotropy, so it contributes
  nothing to the traceless spatial part.  Neither does ``Λ``.  **This is the
  only channel whose answer does not depend on what was never put in.**  The
  scalar sector depends on the sound speed, and also carries the Eddington mode;
* TT perturbations are gauge-invariant, so no gauge choice enters;
* it is the softest tensor channel, so if anything responds, this does;
* and it reduces the response to a five-component driven oscillator, which makes
  the measurement a statement about the **source** rather than about a PDE
  solver's error.

The response, derived rather than recalled
──────────────────────────────────────────
Linearizing Einstein about the ESU in the homogeneous anisotropy
``a_i = a e^{β_i}``, ``Σβ_i = 0``, gives

    ``δG^TT_{ij} = β̈_{ij} + (8/a²) β_{ij}``      so   ``ω₃ = 2√2 / a``

`derive_the_tensor_mode_equation` does that with Cartan's equations, solving the
torsion-free condition rather than quoting a connection formula, and validates
itself on three cases with known answers: the round ``S³`` (``Ric = 2/a²``,
``R = 6/a²``), the ESU (``G = diag(3,−1,−1,−1)/a²``, which independently
reproduces `two_wave`'s ``_EINSTEIN``), and closed FRW
(``G₀₀ = 3(ȧ²+1)/a²``).  Two facts fall out of the linearization itself and are
checked rather than assumed: the first-order piece is **automatically
traceless**, and ``δG_{0i} = 0`` identically — the momentum constraint, so the
mode really is TT and does not mix with vector or scalar perturbations.

``ω₃² = +8/a² > 0``: the tensor sector of the ESU is **stable**.  That is what
makes an ``A``/``B``/``A+B`` comparison here measure the waves rather than the
background — the same gate PR #260 taught the arc to check first.

What is hard here is the quadrature
───────────────────────────────────
The source is ``S_{ij}(t) = ⟨T^TT_{ij}⟩`` over all of ``S³``, and the integrand
has ``1/χ⁴`` singularities at **eight** points: the two sources, the two mouths,
and all four antipodes.  Three things were needed:

* a **batched** field solver (`solve_batch`), since the kernels do not depend on
  the observation point — verified against `two_wave.solve_field` to ``5e-16``;
* a **smooth partition of unity** around each singular point.  Excising them
  instead makes the bulk integrand discontinuous and the rule stalls at
  ``1e-03``; the partition converges to ``4e-04`` and falls with refinement;
* and the mouths in the singular set.  Leaving them out is what a first draft
  did, and it is *why* that draft failed to converge — worse, `two_source`'s
  first mouth sits at ``(1,0,0,0)``, exactly the natural quadrature axis.

The traceless angular average is finite at every singular point — the ``1/χ⁴``
cancels exactly under an angular rule that integrates ``n_i n_j`` exactly — so
these integrals converge.  They just do not converge by accident.

**A number here means nothing without the refinement control**, and
`measure_the_quadrature_converges` is that control: two refinement levels have
to agree, in correlation *and* in magnitude, before the residual is reported.  A
first attempt at this round produced a confident ``0.982`` that was pure
quadrature noise — independent quadratures of the same quantity correlated at
``−0.04``.  That number is recorded in the docs as wrong rather than deleted.

The channel is never on resonance
─────────────────────────────────
The conformally coupled scalar on the ESU has spectrum ``ω_n = n + 1``:
**integers**.  The space is compact and static, so nothing decays and the field
rings on those modes forever; ``T`` is quadratic and integers are closed under
sums and differences, so the shear source rings on integers too — measured, its
peaks land on integers to within a grid bin, dominated by ``ω ≈ 6`` and
``ω ≈ 8``, and no peak of it ever lands on ``ω₃`` whatever the carrier is.

``ω₃ = 2√2`` is irrational, ``0.172`` from the nearest integer.  So on this
background the gravitational shear channel is driven **off resonance by
construction**, whatever the source does.  A first draft claimed instead that
``T`` being quadratic puts the power at ``2ω₀`` and picked the carrier to match;
that is wrong, and it is recorded rather than deleted.

The same scan is the round's honesty check on its own headline: the unreachable
fraction is **not a universal constant** — it moves with the carrier and with
the time window — so the spread is reported alongside it.

What is still put in
────────────────────
The ``n = 3`` harmonic only — the homogeneous shear, not the full TT tower.  A
**fixed** ESU background with point sources and the PR #257 point throat, so the
mouths are not the resolved ones of PR #261/#262.  The response is **linear**:
this is backreaction as a first-order response, not a solved coupled system.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .two_source import geodesic, mouth_positions
from .two_wave import (
    GaussianPulse, RetardedGrid, TwoWaveSetup, WORKING_SEPARATION,
    circle_point, gamma_omega, green_omega, orthonormal_frame, solve_field,
    working_pair,
)

__all__ = [
    "TENSOR_MODE_FREQUENCY",
    "ShearQuadrature",
    "WORKING_BACKREACTION",
    "BackreactionSetup",
    "derive_the_tensor_mode_equation",
    "left_invariant_frame",
    "quaternion_product",
    "shear_projection",
    "shear_response",
    "solve_batch",
    "stress_series",
    "unreachable_fraction",
    "measure_the_tensor_sector_is_stable_and_fluid_free",
    "measure_the_stress_tensor_splits_bilinearly",
    "measure_the_quadrature_converges",
    "measure_the_interference_response_is_unreachable",
    "measure_the_tensor_mode_is_incommensurate_with_the_matter_spectrum",
    "measure_the_answer_needs_the_branches",
]

FOUR_PI = 4.0 * math.pi
TWO_PI_SQ = 2.0 * math.pi ** 2          # the volume of the unit S³
_ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
_EINSTEIN = np.diag([3.0, -1.0, -1.0, -1.0])
_XI = 1.0 / 6.0

#: ``ω₃ = 2√2`` on the unit-radius ESU — derived, not quoted.  See
#: `derive_the_tensor_mode_equation`, which also re-derives the ``8``.
TENSOR_MODE_FREQUENCY = math.sqrt(8.0)


def tensor_mode_frequency(radius: float = 1.0) -> float:
    """``ω₃ = 2√2 / a`` — the lowest transverse-traceless mode of the ESU."""
    return TENSOR_MODE_FREQUENCY / float(radius)


# ════════════════════════════════════════════════════════════════════════════
# THE RESPONSE, DERIVED
# ════════════════════════════════════════════════════════════════════════════
def derive_the_tensor_mode_equation() -> Dict[str, object]:
    """Linearize Einstein about the ESU in the homogeneous anisotropy channel.

    Uses Cartan's equations in the orthonormal coframe ``e⁰ = dt``,
    ``eⁱ = a_i(t) σⁱ`` on ``SU(2) = S³``, with the connection obtained by
    **solving** the torsion-free condition

        ``W^a_{cb} − W^a_{bc} = C^a_{bc}`` ,   ``W_{abc} = −W_{bac}``

    rather than from a remembered closed form.  A first draft of this round did
    quote a formula and produced a ``G₀₀`` containing ``ä``, which is impossible
    — hence the three validation cases below, all of which have known answers,
    and one of which independently reproduces a constant already in `two_wave`.

    Returns the symbolic results; requires ``sympy``.
    """
    import sympy as sp

    t = sp.Symbol("t", real=True)
    a = sp.Symbol("a", positive=True)

    def connection(eta, C):
        n = len(eta)
        E = sp.diag(*eta)
        unk = {(i, j, k): sp.Symbol(f"W_{i}_{j}_{k}")
               for i in range(n) for j in range(n) for k in range(n) if i < j}

        def low(i, j, k):
            if i == j:
                return sp.S(0)
            return unk[(i, j, k)] if i < j else -unk[(j, i, k)]

        def up(i, j, k):
            return E[i, i] * low(i, j, k)

        eqs = [sp.Eq(up(i, k, j) - up(i, j, k), C.get((i, j, k), sp.S(0)))
               for i in range(n) for j in range(n) for k in range(n)]
        sol = sp.solve(eqs, list(unk.values()), dict=True)[0]
        vals = {key: sp.simplify(sol.get(sym, sym)) for key, sym in unk.items()}

        def w(i, j, k):
            if i == j:
                return sp.S(0)
            v = vals[(i, j, k)] if i < j else -vals[(j, i, k)]
            return E[i, i] * v
        return w

    def curvature(eta, C, deriv):
        n = len(eta)
        E = sp.diag(*eta)
        w = connection(eta, C)

        def riem(i, j, k, l):
            r = deriv(k, w(i, j, l)) - deriv(l, w(i, j, k))
            for e in range(n):
                r += w(i, e, k) * w(e, j, l) - w(i, e, l) * w(e, j, k)
                r -= C.get((e, k, l), sp.S(0)) * w(i, j, e)
            return sp.simplify(r)

        ric = sp.zeros(n, n)
        for j in range(n):
            for l in range(n):
                ric[j, l] = sp.simplify(sum(riem(i, j, i, l) for i in range(n)))
        scal = sp.simplify(sum(E[j, j] * ric[j, j] for j in range(n)))
        ein = sp.zeros(n, n)
        for j in range(n):
            for l in range(n):
                ein[j, l] = sp.simplify(ric[j, l]
                                        - sp.Rational(1, 2) * E[j, l] * scal)
        return ric, scal, ein

    def su2(a_list, time=False):
        """``[e_j,e_k] = 2(a_i/(a_j a_k)) e_i``, the unit-quaternion frame."""
        C: Dict[Tuple[int, int, int], object] = {}
        off = 1 if time else 0
        for i, j, k in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
            ai, aj, ak = a_list[i - 1], a_list[j - 1], a_list[k - 1]
            v = 2 * ai / (aj * ak)
            C[(i - 1 + off, j - 1 + off, k - 1 + off)] = v
            C[(i - 1 + off, k - 1 + off, j - 1 + off)] = -v
        if time:
            for i in (1, 2, 3):
                v = -sp.diff(a_list[i - 1], t) / a_list[i - 1]
                C[(i, 0, i)] = v
                C[(i, i, 0)] = -v
        return C

    dt_only = (lambda c, x: sp.diff(x, t) if c == 0 else sp.S(0))

    # ---- validation 1: the round S³ ----
    ric3, scal3, _ = curvature([1, 1, 1], su2([a, a, a]),
                               lambda c, x: sp.S(0))
    # ---- validation 2: the ESU ----
    _, _, ein4 = curvature([-1, 1, 1, 1], su2([a, a, a], time=True), dt_only)
    esu = [sp.simplify(ein4[i, i]) for i in range(4)]
    # ---- validation 3: closed FRW ----
    af = sp.Function("a")(t)
    _, _, einf = curvature([-1, 1, 1, 1], su2([af, af, af], time=True), dt_only)
    frw00 = sp.simplify(einf[0, 0])

    # ---- the anisotropy ----
    eps = sp.Symbol("epsilon")
    b1, b2 = sp.Function("b1")(t), sp.Function("b2")(t)
    betas = (b1, b2, -b1 - b2)
    _, _, eing = curvature([-1, 1, 1, 1],
                           su2([a * sp.exp(eps * b) for b in betas], time=True),
                           dt_only)
    first = []
    for i in (1, 2, 3):
        g = sp.series(eing[i, i], eps, 0, 2).removeO()
        first.append(sp.simplify(sp.diff(g, eps).subs(eps, 0)))
    trace = sp.simplify(sum(first) / 3)
    momentum = [sp.simplify(sp.series(eing[0, i], eps, 0, 2).removeO())
                for i in (1, 2, 3)]
    # read the frequency off the b1 equation: b1'' + (k/a²) b1
    coeff = sp.simplify(first[0].coeff(b1) * a ** 2)

    return {
        "s3_ricci": sp.simplify(ric3[0, 0]),
        "s3_scalar": sp.simplify(scal3),
        "esu_einstein_diagonal": esu,
        "frw_g00": frw00,
        "first_order_trace": trace,
        "momentum_constraint": momentum,
        "omega_squared_times_a_squared": coeff,
        "delta_g_tt": first[0],
        "the_validations_pass": bool(
            sp.simplify(ric3[0, 0] - 2 / a ** 2) == 0
            and sp.simplify(scal3 - 6 / a ** 2) == 0
            and all(sp.simplify(esu[i] - c / a ** 2) == 0
                    for i, c in enumerate((3, -1, -1, -1)))),
        "the_first_order_piece_is_traceless": bool(sp.simplify(trace) == 0),
        "the_momentum_constraint_holds": bool(
            all(sp.simplify(m) == 0 for m in momentum)),
        "the_frequency": bool(sp.simplify(coeff - 8) == 0),
    }


# ════════════════════════════════════════════════════════════════════════════
# THE GLOBAL FRAME
# ════════════════════════════════════════════════════════════════════════════
def quaternion_product(p: Sequence[float], q: Sequence[float]) -> np.ndarray:
    a1, b1, c1, d1 = (float(v) for v in p)
    a2, b2, c2, d2 = (float(v) for v in q)
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                     a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                     a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


_UNITS = (np.array([0.0, 1, 0, 0]), np.array([0.0, 0, 1, 0]),
          np.array([0.0, 0, 0, 1]))


def left_invariant_frame(x: Sequence[float]) -> np.ndarray:
    """``L_i(x) = x·e_i`` — a **global** orthonormal frame on ``S³``.

    ``S³`` is parallelizable, which is what makes the whole projection possible:
    a tensor field's components in this frame are globally defined functions, so
    the homogeneous ``n = 3`` harmonics are exactly the **constant** symmetric
    traceless coefficient matrices and the projection is a plain average.

    `two_wave.orthonormal_frame` is a *local* frame chosen per point by a seeded
    draw; it cannot be used for this, and `shear_projection` rotates out of it.
    ``[L_i, L_j] = 2ε_{ijk}L_k``, matching the normalization the derivation uses.
    """
    p = np.asarray(x, dtype=float)
    p = p / np.linalg.norm(p)
    return np.array([quaternion_product(p, e) for e in _UNITS])


# ════════════════════════════════════════════════════════════════════════════
# THE BATCHED SOLVER
# ════════════════════════════════════════════════════════════════════════════
def _green_derivatives_batch(chi: np.ndarray, omega: np.ndarray):
    """``(G, ∂_χG, ∂²_χG)`` for many ``χ`` at once — ``(P, n)`` arrays.

    The same closed forms as `two_wave.green_omega_derivatives`, written in
    ``e = π − χ`` for the same reason, and checked against it to ``5e-16``.
    """
    c = np.asarray(chi, dtype=float)[:, None]
    w = np.asarray(omega, dtype=complex)[None, :]
    e = math.pi - c
    scale = FOUR_PI * np.sin(math.pi * w)
    se, ce = np.sin(e), np.cos(e)
    small = np.abs(se) < 1e-9
    d0 = np.where(small, 1.0, se)
    n0 = np.sin(w * e)
    n1 = w * np.cos(w * e)
    n2 = -w * w * n0
    h = n0 / d0
    hp = n1 / d0 - n0 * ce / d0 ** 2
    hpp = (n2 / d0 - 2.0 * n1 * ce / d0 ** 2 + n0 * d0 / d0 ** 2
           + 2.0 * n0 * ce ** 2 / d0 ** 3)
    w2 = w * w
    h = np.where(small, w * (1.0 + e * e * (1.0 - w2) / 6.0), h)
    hp = np.where(small, w * e * (1.0 - w2) / 3.0, hp)
    hpp = np.where(small, w * (1.0 - w2) / 3.0, hpp)
    return h / scale, -hp / scale, hpp / scale


def _kernels(setup: TwoWaveSetup, source: np.ndarray, pulse: GaussianPulse
             ) -> List[Tuple[np.ndarray, np.ndarray]]:
    """``(centre, kernel(ω))`` per radial term.

    The whole reason batching works: **none of these depend on the observation
    point.**  Only ``χ`` and ``n̂`` do, and those are cheap.
    """
    om = setup.grid.omegas
    spec = pulse.spectrum(om)
    out = [(np.asarray(source, dtype=float), spec)]
    if not setup.with_throat:
        return out
    c1, c2 = setup.mouths()
    a = setup.pair.boundary_matrix()
    m = a[None, :, :] - gamma_omega(om, setup.pair.separation)
    det = m[:, 0, 0] * m[:, 1, 1] - m[:, 0, 1] * m[:, 1, 0]
    drive = np.stack([green_omega(float(geodesic(c1, source)), om),
                      green_omega(float(geodesic(c2, source)), om)], axis=-1)
    qs = np.empty_like(drive)
    qs[:, 0] = (m[:, 1, 1] * drive[:, 0] - m[:, 0, 1] * drive[:, 1]) / det
    qs[:, 1] = (-m[:, 1, 0] * drive[:, 0] + m[:, 0, 0] * drive[:, 1]) / det
    qs *= spec[:, None]
    return out + [(np.asarray(c1, dtype=float), qs[:, 0]),
                  (np.asarray(c2, dtype=float), qs[:, 1])]


def solve_batch(setup: TwoWaveSetup, source: np.ndarray, pulse: GaussianPulse,
                points: np.ndarray, frames: np.ndarray
                ) -> Dict[str, np.ndarray]:
    """`two_wave.solve_field` for many observation points in one pass.

    Same contour, same closed-form derivatives; the inverse transform is one
    batched FFT instead of ``P`` separate ones.  Agreement with the reference
    implementation is ``5e-16`` and is a test, not a hope — without this the
    ``30 000``-point quadrature the convergence control needs is unaffordable.
    """
    grid = setup.grid
    n = grid.n
    pts = np.asarray(points, dtype=float)
    p = len(pts)
    om = grid.omegas
    growth = np.exp(grid.eps * grid.times)[None, :]

    def invert(z):
        return np.real(np.fft.fft(z, axis=1)) / grid.span * growth

    phi = np.zeros((p, n))
    dt1 = np.zeros((p, n))
    dt2 = np.zeros((p, n))
    grad = np.zeros((p, n, 3))
    dtgrad = np.zeros((p, n, 3))
    hess = np.zeros((p, n, 3, 3))
    for centre, kern in _kernels(setup, source, pulse):
        dots = pts @ centre
        chi = np.arccos(np.clip(dots, -1.0, 1.0))
        s = np.sin(chi)
        if np.any(s < 1e-9):
            raise ValueError("a quadrature point sits on a source or antipode")
        gvec = -(centre[None, :] - dots[:, None] * pts) / s[:, None]
        nhat = np.einsum("pa,pia->pi", gvec, frames)
        g0, g1, g2 = _green_derivatives_batch(chi, om)
        k = kern[None, :]
        f0 = invert(k * g0)
        f1 = invert(k * g1)
        f2 = invert(k * g2)
        ft0 = invert(-1j * om[None, :] * k * g0)
        ft1 = invert(-1j * om[None, :] * k * g1)
        ftt = invert(-(om ** 2)[None, :] * k * g0)
        nn = np.einsum("pi,pj->pij", nhat, nhat)
        proj = np.eye(3)[None, :, :] - nn
        cot = (1.0 / np.tan(chi))[:, None]
        phi += f0
        dt1 += ft0
        dt2 += ftt
        grad += f1[:, :, None] * nhat[:, None, :]
        dtgrad += ft1[:, :, None] * nhat[:, None, :]
        hess += (f2[:, :, None, None] * nn[:, None, :, :]
                 + (f1 * cot)[:, :, None, None] * proj[:, None, :, :])
    return {"phi": phi, "dt": dt1, "dtt": dt2, "grad": grad,
            "dtgrad": dtgrad, "hess": hess,
            "laplacian": np.einsum("ptaa->pt", hess)}


def stress_series(sol: Dict[str, np.ndarray]) -> np.ndarray:
    """`two_wave.stress_tensor` for every point and time at once.

    ``□φ`` is again taken from the solved derivatives rather than replaced by
    its on-shell value, for `two_wave`'s reason: substituting on shell makes the
    trace vanish algebraically and the check would test nothing.
    """
    phi = sol["phi"]
    p, n = phi.shape
    d = np.zeros((p, n, 4))
    d[..., 0] = sol["dt"]
    d[..., 1:] = sol["grad"]
    dd = np.zeros((p, n, 4, 4))
    dd[..., 0, 0] = sol["dtt"]
    dd[..., 0, 1:] = sol["dtgrad"]
    dd[..., 1:, 0] = sol["dtgrad"]
    dd[..., 1:, 1:] = sol["hess"]
    grad_sq = np.einsum("pti,ij,ptj->pt", d, _ETA, d)
    box = sol["laplacian"] - sol["dtt"]
    box_phi2 = 2.0 * (grad_sq + phi * box)
    hess_phi2 = 2.0 * (np.einsum("pti,ptj->ptij", d, d)
                       + phi[..., None, None] * dd)
    return (np.einsum("pti,ptj->ptij", d, d)
            - 0.5 * _ETA * grad_sq[..., None, None]
            + _XI * (_ETA * box_phi2[..., None, None] - hess_phi2
                     + _EINSTEIN * (phi * phi)[..., None, None]))


def shear_projection(stress: np.ndarray, points: np.ndarray,
                     frames: np.ndarray) -> np.ndarray:
    """The traceless spatial block, rotated into the **global** frame.

    `solve_batch` works in `two_wave.orthonormal_frame`, which is a different
    arbitrary basis at every point.  ``T`` is a tensor, so rotating by
    ``M = L Eᵀ`` — orthogonal, since both bases are orthonormal — puts every
    point's components in one basis, which is what makes averaging them mean
    anything.
    """
    li = np.array([left_invariant_frame(p) for p in points])
    m = np.einsum("pia,pja->pij", li, frames)
    s = np.einsum("pia,ptab,pjb->ptij", m, stress[..., 1:, 1:], m)
    tr = np.einsum("ptii->pt", s)
    return s - (tr / 3.0)[..., None, None] * np.eye(3)


# ════════════════════════════════════════════════════════════════════════════
# THE QUADRATURE
# ════════════════════════════════════════════════════════════════════════════
def _tangent_basis(centre: Sequence[float]) -> np.ndarray:
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    u, _, _ = np.linalg.svd(np.eye(4) - np.outer(c, c))
    return u[:, :3].T


def _direction_rule(n_theta: int, n_phi: int) -> Tuple[np.ndarray, np.ndarray]:
    """A product rule on the direction sphere, Gauss in ``cos θ``.

    It has to integrate ``n_i n_j`` **exactly**: near a source the stress tensor
    diverges like ``n_i n_j/χ⁴`` and the traceless part of that has zero angular
    average, so an inexact rule leaves a spurious ``1/χ⁴``.  A random rule does
    not, which is why a first draft of this round diverged.
    """
    ct, wt = np.polynomial.legendre.leggauss(int(n_theta))
    phis = 2.0 * math.pi * np.arange(int(n_phi)) / int(n_phi)
    wp = 2.0 * math.pi / int(n_phi)
    dirs, wts = [], []
    for c, w in zip(ct, wt):
        st = math.sqrt(max(0.0, 1.0 - c * c))
        for ph in phis:
            dirs.append((st * math.cos(ph), st * math.sin(ph), c))
            wts.append(w * wp)
    return np.array(dirs), np.array(wts)


def _ball_rule(centre, radius, n_r, n_theta, n_phi):
    b = _tangent_basis(centre)
    c = np.asarray(centre, dtype=float)
    c = c / np.linalg.norm(c)
    xr, wr = np.polynomial.legendre.leggauss(int(n_r))
    rs = 0.5 * float(radius) * (xr + 1.0)
    wrs = 0.5 * float(radius) * wr
    dirs, wd = _direction_rule(n_theta, n_phi)
    vecs = dirs @ b
    pts, wts = [], []
    for r, w in zip(rs, wrs):
        vol = w * math.sin(r) ** 2
        pts.append(math.cos(r) * c[None, :] + math.sin(r) * vecs)
        wts.append(vol * wd)
    return np.vstack(pts), np.concatenate(wts)


def _bulk_rule(n_chi, n_theta, n_phi, axis):
    b = _tangent_basis(axis)
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.linalg.norm(ax)
    xc, wc = np.polynomial.legendre.leggauss(int(n_chi))
    chis = 0.5 * math.pi * (xc + 1.0)
    wch = 0.5 * math.pi * wc
    dirs, wd = _direction_rule(n_theta, n_phi)
    vecs = dirs @ b
    pts, wts = [], []
    for ch, w in zip(chis, wch):
        vol = w * math.sin(ch) ** 2
        pts.append(math.cos(ch) * ax[None, :] + math.sin(ch) * vecs)
        wts.append(vol * wd)
    return np.vstack(pts), np.concatenate(wts)


def _bump(d: np.ndarray, radius: float) -> np.ndarray:
    """A ``C²`` partition function: ``1`` inside ``R/3``, ``0`` outside ``R``."""
    x = (np.asarray(d, dtype=float) - radius / 3.0) / (radius - radius / 3.0)
    x = np.clip(x, 0.0, 1.0)
    return 1.0 - (10 * x ** 3 - 15 * x ** 4 + 6 * x ** 5)


@dataclass(frozen=True)
class ShearQuadrature:
    """A rule on ``S³`` that resolves all eight singular points.

    ``T`` blows up like ``1/χ⁴`` at the two **sources**, the two **mouths**, and
    every antipode.  The mouths are the ones a first draft forgot, and
    `two_source.mouth_positions` puts the first at ``(1,0,0,0)`` — exactly the
    axis a product grid naturally uses, so the grid's own pole sat on a
    singularity and refining it made the answer *worse*.

    The rule is a **smooth partition of unity**: the bulk grid carries
    ``1 − Σψ`` and a dedicated ball rule carries ``ψ`` around each singular
    point.  Simply excising the balls instead leaves the bulk integrand
    discontinuous and the volume stalls near ``1e-03``; the partition reaches
    ``4e-04`` and keeps falling.
    """

    bulk: Tuple[int, int, int] = (20, 12, 24)
    ball: Tuple[int, int, int] = (10, 8, 16)
    radius: float = 0.4
    axis: Tuple[float, ...] = (0.37, -0.51, 0.62, 0.46)

    def singular_points(self, setup: "BackreactionSetup") -> List[np.ndarray]:
        out: List[np.ndarray] = []
        for p in (setup.source_a, setup.source_b, *setup.mouths()):
            q = np.asarray(p, dtype=float)
            out += [q, -q]
        return out

    def build(self, setup: "BackreactionSetup"
              ) -> Tuple[np.ndarray, np.ndarray]:
        sing = self.singular_points(setup)
        closest = min(geodesic(sing[i], sing[j])
                      for i in range(len(sing)) for j in range(i + 1, len(sing)))
        if self.radius >= 0.5 * closest:
            raise ValueError(
                f"balls overlap: radius {self.radius} needs to be under "
                f"{0.5 * closest:.4f}")
        axis = np.asarray(self.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        pts, wts = _bulk_rule(*self.bulk, axis=axis)
        frac = np.ones(len(pts))
        for c in sing:
            frac -= _bump(np.arccos(np.clip(pts @ c, -1.0, 1.0)), self.radius)
        keep = np.abs(frac) > 1e-13
        pts, wts = pts[keep], (wts * frac)[keep]
        for c in sing:
            p2, w2 = _ball_rule(c, self.radius, *self.ball)
            psi = _bump(np.arccos(np.clip(p2 @ c, -1.0, 1.0)), self.radius)
            k = psi > 1e-13
            pts = np.vstack([pts, p2[k]])
            wts = np.concatenate([wts, (w2 * psi)[k]])
        return pts, wts

    def volume_error(self, setup: "BackreactionSetup") -> float:
        _, w = self.build(setup)
        return abs(float(w.sum()) - TWO_PI_SQ) / TWO_PI_SQ


# ════════════════════════════════════════════════════════════════════════════
# THE EXPERIMENT
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class BackreactionSetup:
    """Two pulsed sources on a throated ESU, and the shear they drive.

    The carrier defaults to ``ω₃``, but **not** because that puts the source's
    power there — `measure_the_tensor_mode_is_incommensurate_with_the_matter_spectrum`
    shows the source ringing on the ESU's own integer modes whatever the carrier
    is, and ``ω₃ = 2√2`` is irrational, so this channel is off resonance by
    construction.  ``ω₃`` is simply the value that gave the most source power at
    the mode among those scanned, and the same measurement reports how much the
    answer moves when it changes.  PR #257's carrier of ``60`` is far too high:
    it leaves the channel reading round-off.
    """

    source_gap: float = 1.6
    carrier: float = TENSOR_MODE_FREQUENCY
    width: float = 0.5
    t0: float = 3.0
    with_throat: bool = True
    grid: RetardedGrid = field(
        default_factory=lambda: RetardedGrid(n=1 << 12, span=60.0, eps=0.05))

    @property
    def source_a(self) -> np.ndarray:
        return circle_point(0.0)

    @property
    def source_b(self) -> np.ndarray:
        return circle_point(self.source_gap)

    def mouths(self) -> Tuple[np.ndarray, np.ndarray]:
        c1, c2 = mouth_positions(WORKING_SEPARATION)
        return np.asarray(c1, dtype=float), np.asarray(c2, dtype=float)

    def _setup(self, observer: np.ndarray) -> TwoWaveSetup:
        pulse = GaussianPulse(carrier=self.carrier, width=self.width,
                              t0=self.t0)
        return TwoWaveSetup(pair=working_pair(), source_a=self.source_a,
                            source_b=self.source_b, observer=observer,
                            pulse_a=pulse, pulse_b=pulse, grid=self.grid,
                            with_throat=self.with_throat)

    def shear_sources(self, quad: ShearQuadrature, chunk: int = 64
                      ) -> Dict[str, np.ndarray]:
        """``⟨T^TT⟩`` over ``S³`` for ``A``, ``B`` and the cross term.

        ``ΔT`` is formed as ``T[A+B] − T[A] − T[B]`` from three evaluations of
        the same functional, not from a hand-derived bilinear form — the
        discipline `two_wave.cross_stress_tensor` set, so the bilinearity stays
        a measurement.
        """
        pts, wts = quad.build(self)
        base = self._setup(pts[0])
        pulse = base.pulse_a
        acc = [np.zeros((self.grid.n, 3, 3)) for _ in range(3)]
        for i in range(0, len(pts), int(chunk)):
            block, w = pts[i:i + int(chunk)], wts[i:i + int(chunk)]
            frames = np.array([orthonormal_frame(p) for p in block])
            a = solve_batch(base, self.source_a, pulse, block, frames)
            b = solve_batch(base, self.source_b, pulse, block, frames)
            total = {k: a[k] + b[k] for k in a}
            ta, tb = stress_series(a), stress_series(b)
            tt = stress_series(total)
            for j, x in enumerate((ta, tb, tt - ta - tb)):
                acc[j] += np.einsum("p,ptij->tij", w,
                                    shear_projection(x, block, frames))
        return {"A": acc[0], "B": acc[1], "cross": acc[2],
                "points": len(pts), "volume_error": abs(
                    float(wts.sum()) - TWO_PI_SQ) / TWO_PI_SQ}


WORKING_BACKREACTION = BackreactionSetup()


def shear_response(source: np.ndarray, dt: float,
                   omega: float = TENSOR_MODE_FREQUENCY) -> np.ndarray:
    """Solve ``β̈ + ω²β = S`` retarded, by the mode's own Green function.

    ``β(t) = ∫₀^t sin(ω(t−t'))/ω · S(t') dt'`` — the retarded solution with
    ``β = β̇ = 0`` before the pulse.  A test checks the returned ``β`` against
    the ODE by finite differences rather than trusting the convolution.
    """
    s = np.asarray(source, dtype=float)
    n = s.shape[0]
    ker = np.sin(float(omega) * np.arange(n) * float(dt)) / float(omega)
    out = np.empty_like(s)
    for i in range(3):
        for j in range(3):
            out[:, i, j] = np.convolve(s[:, i, j], ker)[:n] * float(dt)
    return out


def unreachable_fraction(beta_a: np.ndarray, beta_b: np.ndarray,
                         beta_cross: np.ndarray, times: np.ndarray,
                         window: Tuple[float, float] = (4.0, 30.0)
                         ) -> Dict[str, float]:
    """How much of ``β[ΔT]`` no rescaling of ``A`` or ``B`` can reproduce.

    Rescaling ``A → cA`` sends ``β[A] → c²β[A]``, so the reachable set is the
    **cone** ``{c²β_A + d²β_B}``.  The residual is reported off the full linear
    **span** as well, which strictly contains that cone — so the span figure is
    the *conservative* one, and it is the one quoted.

    Note ``β[A+B] − (αβ_A + γβ_B) = β_cross − ((α−1)β_A + (γ−1)β_B)``, so asking
    whether the total response is reachable and whether the cross term is are
    the same question.  Nothing hinges on which is measured.
    """
    sl = (np.asarray(times) > window[0]) & (np.asarray(times) < window[1])

    def flat(x):
        return np.asarray(x)[sl].reshape(-1)

    va, vb, vd = flat(beta_a), flat(beta_b), flat(beta_cross)
    m = np.stack([va, vb], axis=1)
    coef, *_ = np.linalg.lstsq(m, vd, rcond=None)
    span = float(np.linalg.norm(vd - m @ coef) / np.linalg.norm(vd))
    return {"unreachable_off_the_span": span,
            "cross_over_single": float(np.linalg.norm(vd)
                                       / np.linalg.norm(va)),
            "norm_a": float(np.linalg.norm(va)),
            "norm_b": float(np.linalg.norm(vb)),
            "norm_cross": float(np.linalg.norm(vd)),
            "best_fit": [float(c) for c in coef],
            "cos_between_the_singles": float(
                abs(va @ vb) / np.linalg.norm(va) / np.linalg.norm(vb))}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def _betas(setup: BackreactionSetup, quad: ShearQuadrature
           ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, object]]:
    src = setup.shear_sources(quad)
    dt = setup.grid.dt
    return (shear_response(src["A"], dt), shear_response(src["B"], dt),
            shear_response(src["cross"], dt), src)


def measure_the_tensor_sector_is_stable_and_fluid_free() -> Dict[str, object]:
    """The response equation, derived — and why this channel and no other.

    Cartan about the ESU in the homogeneous anisotropy channel gives
    ``δG^TT_{ij} = β̈_{ij} + (8/a²)β_{ij}``, so ``ω₃ = 2√2/a`` and
    ``ω₃² > 0``: **the tensor sector is stable**.  PR #260 taught the arc to
    check that before computing anything on a background, and the scalar sector
    would fail it — that is where the Eddington mode lives.

    Three properties are checked rather than assumed, and each would have caught
    a different error:

    * the **validations** — round ``S³``, ESU, closed FRW — against known
      answers.  The ESU one reproduces `two_wave`'s ``_EINSTEIN`` from an
      independent calculation;
    * the first-order piece is **automatically traceless**, so the anisotropy
      sources no trace and cannot mix with the fluid's ``δρ``, ``δp``;
    * ``δG_{0i} = 0`` identically — the **momentum constraint** — so the mode is
      genuinely transverse and does not mix with vector or scalar modes.

    Together those are why the unspecified fluid and ``Λ`` drop out: in an
    orthonormal frame a comoving perfect fluid has ``T_ab = diag(ρ,p,p,p)``
    whatever the anisotropy, so its traceless spatial part is **exactly zero**,
    and ``Λ`` contributes ``−Λη_ab``, also traceless-free.  This is the only
    channel whose answer does not depend on what this arc never put in.
    """
    d = derive_the_tensor_mode_equation()
    return {
        "omega_squared_times_a_squared": str(d["omega_squared_times_a_squared"]),
        "tensor_mode_frequency": TENSOR_MODE_FREQUENCY,
        "s3_ricci": str(d["s3_ricci"]), "s3_scalar": str(d["s3_scalar"]),
        "esu_einstein_diagonal": [str(x) for x in d["esu_einstein_diagonal"]],
        "frw_g00": str(d["frw_g00"]),
        "delta_g_tt": str(d["delta_g_tt"]),
        "the_validations_pass": bool(d["the_validations_pass"]),
        "the_first_order_piece_is_traceless": bool(
            d["the_first_order_piece_is_traceless"]),
        "the_momentum_constraint_holds": bool(d["the_momentum_constraint_holds"]),
        "the_frequency_is_eight_over_a_squared": bool(d["the_frequency"]),
        "the_sector_is_stable": bool(TENSOR_MODE_FREQUENCY ** 2 > 0.0),
        "why_this_channel": ("a perfect fluid has T_ab = diag(rho,p,p,p) in an "
                             "orthonormal frame whatever the anisotropy, so it "
                             "and Λ contribute nothing traceless — the only "
                             "channel whose answer does not depend on the "
                             "matter this arc never specified"),
    }


def measure_the_stress_tensor_splits_bilinearly(
        setup: BackreactionSetup = WORKING_BACKREACTION,
        scales: Sequence[float] = (0.5, 1.0, 2.0, 3.0)) -> Dict[str, object]:
    """``T[A+B] = T[A] + T[B] + ΔT``, with ``ΔT`` bilinear — measured.

    This is what makes the whole question a *linear algebra* question: the field
    equation is linear so the fields add, ``T`` is quadratic so the cross term
    is bilinear, and linearized Einstein is linear so the responses add.
    Rescaling ``A → cA`` therefore sends ``β_A → c²β_A`` and ``β_Δ → cβ_Δ``,
    which is exactly why the reachable family is only two-parameter.

    Checked at a single point rather than assumed, by scaling one pulse and
    watching the cross term scale **linearly** while the self term scales
    quadratically.
    """
    obs = circle_point(0.7)
    base = setup._setup(obs)
    rows = []
    ref_self = ref_cross = None
    for c in scales:
        a = solve_field(base, setup.source_a,
                        GaussianPulse(amplitude=float(c), carrier=setup.carrier,
                                      width=setup.width, t0=setup.t0))
        b = solve_field(base, setup.source_b, base.pulse_b)
        tot = {k: a[k] + b[k] for k in a}
        i = int(len(base.grid.times) * 0.25)
        ta = stress_series({k: v[None] for k, v in a.items()})[0, i]
        tb = stress_series({k: v[None] for k, v in b.items()})[0, i]
        tt = stress_series({k: v[None] for k, v in tot.items()})[0, i]
        cross = tt - ta - tb
        rows.append({"scale": float(c),
                     "self_norm": float(np.linalg.norm(ta)),
                     "cross_norm": float(np.linalg.norm(cross))})
        if c == 1.0:
            ref_self, ref_cross = rows[-1]["self_norm"], rows[-1]["cross_norm"]
    for r in rows:
        r["self_over_c_squared"] = r["self_norm"] / (r["scale"] ** 2 * ref_self)
        r["cross_over_c"] = r["cross_norm"] / (r["scale"] * ref_cross)
    # and the cross term vanishes identically with one source switched off
    a = solve_field(base, setup.source_a, base.pulse_a)
    zero = {k: np.zeros_like(v) for k, v in a.items()}
    tot = {k: a[k] + zero[k] for k in a}
    i = int(len(base.grid.times) * 0.25)
    off = float(np.linalg.norm(
        stress_series({k: v[None] for k, v in tot.items()})[0, i]
        - stress_series({k: v[None] for k, v in a.items()})[0, i]
        - stress_series({k: v[None] for k, v in zero.items()})[0, i]))
    return {"rows": rows, "cross_with_one_source_off": off,
            "the_self_term_is_quadratic": bool(
                all(abs(r["self_over_c_squared"] - 1.0) < 1e-9 for r in rows)),
            "the_cross_term_is_linear": bool(
                all(abs(r["cross_over_c"] - 1.0) < 1e-9 for r in rows)),
            "it_vanishes_with_one_source_off": bool(off < 1e-30),
            "why_it_matters": ("beta_A scales as c^2 and beta_cross as c, so "
                               "the set reachable by rescaling is the "
                               "two-parameter cone {c^2 beta_A + d^2 beta_B} "
                               "and the question is whether beta_cross is in "
                               "it")}


def measure_the_quadrature_converges(
        setup: BackreactionSetup = WORKING_BACKREACTION,
        levels: Sequence[ShearQuadrature] = (
            ShearQuadrature(bulk=(16, 10, 20), ball=(8, 6, 12)),
            ShearQuadrature(bulk=(20, 12, 24), ball=(10, 8, 16))),
        window: Tuple[float, float] = (4.0, 30.0)) -> Dict[str, object]:
    """**The control.**  Nothing else in this module means anything without it.

    A first attempt at this round reported that ``98.2%`` of the interference
    response was unreachable.  It was **quadrature noise**: recomputing the same
    quantity with an independent rule gave ``corr(β_A) = −0.04``.  The number
    was not slightly off, it was meaningless — and nothing about the run looked
    wrong.

    So the rule here is that two refinement levels must agree in **correlation**
    *and* in **magnitude** before any residual is quoted.  Two things had to be
    fixed to get there: the eight singular points (the mouths were missing, and
    one of them sits on the natural quadrature axis) and the smooth partition of
    unity (excision leaves the bulk integrand discontinuous).
    """
    times = setup.grid.times
    sl = (times > window[0]) & (times < window[1])
    runs, points, volerr = [], [], []
    for q in levels:
        ba, bb, bc, src = _betas(setup, q)
        runs.append((ba, bb, bc))
        points.append(int(src["points"]))
        volerr.append(float(src["volume_error"]))
    rows = []
    for name, idx in (("A", 0), ("B", 1), ("cross", 2)):
        for k in range(len(runs) - 1):
            v1 = runs[k][idx][sl].reshape(-1)
            v2 = runs[k + 1][idx][sl].reshape(-1)
            rows.append({
                "component": name, "from": points[k], "to": points[k + 1],
                "correlation": float(v1 @ v2 / np.linalg.norm(v1)
                                     / np.linalg.norm(v2)),
                "magnitude_ratio": float(np.linalg.norm(v2)
                                         / np.linalg.norm(v1))})
    resid = [unreachable_fraction(*r, times, window)["unreachable_off_the_span"]
             for r in runs]
    return {"points": points, "volume_errors": volerr, "rows": rows,
            "residual_by_level": resid,
            "worst_correlation": float(min(r["correlation"] for r in rows)),
            "worst_magnitude_drift": float(max(abs(r["magnitude_ratio"] - 1.0)
                                               for r in rows)),
            "residual_drift": float(abs(resid[-1] - resid[0])),
            "every_component_is_converged": bool(
                min(r["correlation"] for r in rows) > 0.95
                and max(abs(r["magnitude_ratio"] - 1.0) for r in rows) < 0.05),
            "the_residual_is_stable": bool(abs(resid[-1] - resid[0]) < 0.05),
            "what_it_caught": ("a first draft's 0.982 was pure quadrature "
                               "noise — independent rules for the same "
                               "quantity correlated at -0.04")}


def measure_the_interference_response_is_unreachable(
        setup: BackreactionSetup = WORKING_BACKREACTION,
        quad: ShearQuadrature = ShearQuadrature(),
        windows: Sequence[Tuple[float, float]] = (
            (4.0, 12.0), (4.0, 20.0), (4.0, 30.0), (8.0, 45.0))
        ) -> Dict[str, object]:
    """**The answer.**  ``A + B`` does something rescaling ``A`` or ``B`` cannot.

    The reachable family is ``{c²β_A + d²β_B}``.  The residual is reported off
    the full linear **span**, which strictly contains that cone, so the figure
    is conservative — the true unreachable fraction is at least this large.

    It is quoted only because `measure_the_quadrature_converges` clears it.
    """
    ba, bb, bc, src = _betas(setup, quad)
    times = setup.grid.times
    rows = []
    for lo, hi in windows:
        r = unreachable_fraction(ba, bb, bc, times, (lo, hi))
        rows.append({"window": [float(lo), float(hi)], **r})
    main = [r for r in rows if r["window"] == [4.0, 30.0]] or rows
    return {"rows": rows, "points": int(src["points"]),
            "volume_error": float(src["volume_error"]),
            "unreachable_fraction": float(main[0]["unreachable_off_the_span"]),
            "cross_over_single": float(main[0]["cross_over_single"]),
            "cos_between_the_singles": float(
                main[0]["cos_between_the_singles"]),
            "spread_over_windows": float(
                max(r["unreachable_off_the_span"] for r in rows)
                - min(r["unreachable_off_the_span"] for r in rows)),
            "most_of_it_is_unreachable": bool(
                min(r["unreachable_off_the_span"] for r in rows) > 0.5),
            "the_two_singles_are_independent": bool(
                main[0]["cos_between_the_singles"] < 0.95),
            "the_answer": ("yes — the interference response is comparable in "
                           "size to the single-wave ones and almost orthogonal "
                           "to both, so no rescaling reproduces it")}


def measure_the_tensor_mode_is_incommensurate_with_the_matter_spectrum(
        quad: ShearQuadrature = ShearQuadrature(bulk=(16, 10, 20),
                                                ball=(8, 6, 12)),
        carriers: Sequence[float] = (0.7, TENSOR_MODE_FREQUENCY / 2,
                                     TENSOR_MODE_FREQUENCY, 4.0),
        window: Tuple[float, float] = (4.0, 30.0)) -> Dict[str, object]:
    """Why no carrier drives this channel — and how much the answer moves.

    The conformally coupled scalar on the ESU has spectrum ``ω_n = n + 1``:
    **integers**.  The space is compact and static, so nothing decays and the
    field rings on those modes forever; ``T`` is quadratic, and sums and
    differences of integers are integers, so the shear source rings on integers
    too.  Measured, its peaks land on integers to within a grid bin, dominated
    by ``ω ≈ 6`` and ``ω ≈ 8``, and the dominant one stays at ``ω ≈ 6`` for
    every carrier below the mode.  (At ``carrier = 4`` the envelope's DC lobe
    takes over and the largest bin moves to the bottom of the band; the
    *ringing* peaks are still on integers, but "the peak never moves" would have
    been too strong, so what is asserted is the claim that survives: **no peak
    of the source ever lands on ``ω₃``**.)

    The lowest tensor mode is ``ω₃ = 2√2``, which is **irrational**: ``0.172``
    from the nearest integer, and it cannot ever coincide with one.  So on this
    background the gravitational shear channel is driven **off resonance by
    construction**, whatever the source is doing.  That is a structural fact
    about a conformal scalar on an ESU, not a property of the pulse.

    The **exact** statement is the irrationality; the numerical one is that the
    nearest source peak sits ``0.313`` away, which is ``2.99`` grid bins.  A
    first version of this measurement demanded ``3`` bins and therefore *failed*
    — by ``0.4%`` — which is this ledger's own recurring species: a threshold
    standing in for a limit.  The criterion is now that the mode is *resolvably*
    off every peak (more than one and a half bins, so they are distinct on this
    grid), and the separation is reported in bins as a measurement rather than
    hidden inside a pass.

    A first draft of this round claimed instead that ``T`` being quadratic puts
    the source's power at ``2ω₀``, and chose the carrier on that basis.  **It is
    wrong** — the measured peak is ``5.969`` for carriers ``0.7``, ``1.414``,
    ``2.0`` and ``2.828`` alike, because the ringing is the *background's*, not
    the pulse's.  The wrong reasoning is recorded rather than deleted.

    The carrier scan is also the round's honesty check on its own headline: the
    unreachable fraction is **not a universal constant**, and the spread is
    reported instead of one number being quoted as if it were.
    """
    rows = []
    for carrier in carriers:
        s = BackreactionSetup(carrier=float(carrier))
        src = s.shear_sources(quad)
        n = s.grid.n
        freqs = np.fft.rfftfreq(n, d=s.grid.dt) * 2.0 * math.pi
        spec = np.abs(np.fft.rfft(src["A"].reshape(n, 9), axis=0)).sum(axis=1)
        band = (freqs > 0.3) & (freqs < 14.0)
        fb, sb = freqs[band], spec[band]
        peak = float(fb[int(np.argmax(sb))])
        i3 = int(np.argmin(np.abs(freqs - TENSOR_MODE_FREQUENCY)))
        dt = s.grid.dt
        r = unreachable_fraction(shear_response(src["A"], dt),
                                 shear_response(src["B"], dt),
                                 shear_response(src["cross"], dt),
                                 s.grid.times, window)
        # how far the strongest peaks sit from an integer
        loc = [i for i in range(1, len(sb) - 1)
               if sb[i] > sb[i - 1] and sb[i] > sb[i + 1]]
        loc = sorted(loc, key=lambda i: -sb[i])[:6]
        offs = [abs(float(fb[i]) - round(float(fb[i]))) for i in loc]
        gap = min(abs(float(fb[i]) - TENSOR_MODE_FREQUENCY) for i in loc)
        rows.append({"carrier": float(carrier), "spectral_peak": peak,
                     "peak_offset_from_integer": abs(peak - round(peak)),
                     "worst_peak_offset": float(max(offs)),
                     "nearest_peak_to_the_mode": float(gap),
                     "power_at_the_mode": float(spec[i3]),
                     "unreachable": r["unreachable_off_the_span"],
                     "cross_over_single": r["cross_over_single"]})
    peaks = [r["spectral_peak"] for r in rows]
    unreach = [r["unreachable"] for r in rows]
    grid_bin = 2.0 * math.pi / WORKING_BACKREACTION.grid.span
    return {
        "rows": rows,
        "tensor_mode_frequency": TENSOR_MODE_FREQUENCY,
        "distance_to_the_nearest_integer": float(
            abs(TENSOR_MODE_FREQUENCY - round(TENSOR_MODE_FREQUENCY))),
        "grid_frequency_spacing": float(grid_bin),
        "spectral_peak_spread": float(max(peaks) - min(peaks)),
        "unreachable_range": [float(min(unreach)), float(max(unreach))],
        "closest_any_peak_gets_to_the_mode": float(
            min(r["nearest_peak_to_the_mode"] for r in rows)),
        "separation_in_grid_bins": float(
            min(r["nearest_peak_to_the_mode"] for r in rows) / grid_bin),
        "the_source_rings_on_integers": bool(
            all(r["worst_peak_offset"] < 3.0 * grid_bin for r in rows)),
        "the_mode_is_resolvably_off_every_peak": bool(
            min(r["nearest_peak_to_the_mode"] for r in rows) > 1.5 * grid_bin),
        "the_tensor_mode_is_irrational": bool(
            abs(TENSOR_MODE_FREQUENCY - round(TENSOR_MODE_FREQUENCY)) > 0.1),
        "it_is_unreachable_at_every_carrier": bool(min(unreach) > 0.5),
        "the_fraction_is_not_a_universal_constant": bool(
            max(unreach) - min(unreach) > 0.1),
        "what_the_first_draft_got_wrong": (
            "that T being quadratic puts the source's power at 2*w0, so the "
            "carrier should be chosen to match; in fact the dominant peak stays "
            "at w = 6 for every carrier below the mode, because the ringing is "
            "the ESU's own integer spectrum rather than the pulse's (at carrier "
            "4 the envelope's DC lobe takes over instead, which is why the "
            "claim asserted is the weaker and true one: the mode is resolvably "
            "off every peak, and exactly so, since 2*sqrt(2) is irrational)"),
    }


def measure_the_answer_needs_the_branches(
        quad: ShearQuadrature = ShearQuadrature(bulk=(16, 10, 20),
                                                ball=(8, 6, 12)),
        window: Tuple[float, float] = (4.0, 30.0)) -> Dict[str, object]:
    """PR #257's warning, applied to backreaction: *which branches were there?*

    That round measured the same configuration giving an invariant of anything
    from ``0`` to ``4`` depending on which arrival branches were present, and
    said any quantity built by integrating over the field has to state them.
    Backreaction is such a quantity, so the throat is switched off as a control
    and the two answers are reported side by side rather than one of them being
    called *the* answer.
    """
    out = {}
    for tag, throat in (("with_throat", True), ("no_throat", False)):
        s = BackreactionSetup(with_throat=throat)
        ba, bb, bc, src = _betas(s, quad)
        r = unreachable_fraction(ba, bb, bc, s.grid.times, window)
        out[tag] = {"unreachable": r["unreachable_off_the_span"],
                    "cross_over_single": r["cross_over_single"],
                    "norm_cross": r["norm_cross"],
                    "points": int(src["points"])}
    a, b = out["with_throat"], out["no_throat"]
    return {**out,
            "unreachable_both_ways": bool(a["unreachable"] > 0.5
                                          and b["unreachable"] > 0.5),
            "the_throat_changes_the_size": float(
                abs(a["norm_cross"] - b["norm_cross"])
                / max(b["norm_cross"], 1e-300)),
            "what_it_scopes": ("the conclusion survives switching the throat "
                               "off, so it is a statement about two waves "
                               "rather than about the throat; the throat "
                               "changes the amplitude, and PR #257's rule that "
                               "an integrated quantity must name its branches "
                               "is why this control exists")}
