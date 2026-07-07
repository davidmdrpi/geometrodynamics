"""
The nonlinear no-signaling audit - companion probe (PR #204).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE REGISTER ITEM, AND THE OBSTRUCTION
--------------------------------------
#198 derived the Born rule at dBB grade CONDITIONAL on (2) the linear
measurement regime, and the #200 register carries "nonlinear measurement
theory beyond the linear test-throat regime" as an open item.  The sharpest
edge of that item is NO-SIGNALING: by the Gisin/Polchinski theorems,
deterministic NONLINEAR modifications of quantum mechanics generically
allow superluminal signaling through entangled states - and the BAM pilot
equation is genuinely nonlinear (Phi[rho] and q[rho]).  If BAM signals
superluminally in quantum equilibrium in its causal completion, it is
refuted by relativity and by experiment.  The standard no-signaling
theorem does not protect it (that theorem is linear), and the standard
Gisin proof does not automatically convict it (that proof uses the
projection postulate, which BAM does not have - measurement is dynamical
transport, #198).  So the question must be settled BY CONSTRUCTION, on
the running dynamics.  That is this probe.  Deliverable:
``docs/nonlinear_no_signaling_audit.md``.

THE RESULT (measured, all on the live psi-Phi-q dynamics)
  * THE CHANNEL IS EXACTLY THE GRAVITATIONAL FIELD.  A local unitary kick
    at A produces a density response at B (separation 30) that is 4e6x
    the linear kinematic floor - but clamping Phi to the no-kick history
    collapses it BACK to the floor (suppression ~6e6), and the response
    is linear in G (ratio 2.14 under G -> G/2).  Nothing else carries it
    (the order field q is local).
  * RETARDATION CONFINES IT TO THE CONE.  Replacing the instantaneous
    Poisson Phi (the Newtonian c -> infinity limit) by the causal wave
    equation Box Phi = 4 pi G rho, the response at B stays at the
    MACHINE floor (~1e-15) until the front arrives: c*t_front = 25-28
    for c in {8,12,16,20} vs geometric distance 25, and the c -> infinity
    limit recovers the instantaneous amplitude.  The "signaling" of the
    weak-field model is the action-at-a-distance of its NEWTONIAN
    APPROXIMATION, not of the theory.
  * EQUIVARIANCE SURVIVES RETARDATION.  The retarded Phi is still a REAL
    potential, so the #198 continuity theorem is untouched: norm exact,
    continuity residual at integrator error, a 20000-throat Born ensemble
    stays at sampling noise through the kicked retarded evolution.
    No-signaling and the Born rule hold SIMULTANEOUSLY.
  * THE ENTANGLED SECTOR, LINEAR PART: on a two-throat entangled state
    psi(x1,x2), a local unitary on throat 1 leaves the x2-marginal
    invariant to MACHINE precision (the no-signaling theorem, verified on
    the discrete flow), an equilibrium (Born) trajectory ensemble shows
    no kick dependence at B (KS at sampling noise), and a NON-equilibrium
    ensemble DOES (KS ~ 3x noise) - Valentini's signal-locality boundary,
    reproduced on the BAM transport: equilibrium is load-bearing, exactly
    as in dBB.
  * THE ENTANGLED SECTOR, NONLINEAR PART: adding the BAM mean-field
    gravity Phi[rho1+rho2] to the entangled pair, the x2-marginal
    invariance IS violated - at O(G), immediately, when Phi is
    instantaneous - and the violation is CONFINED BEHIND THE FRONT when
    Phi is retarded (response at B at machine floor until t ~ d/c).
    The Gisin channel exists, is exactly the physical gravitational
    coupling, and is causal in the completion.

Tests:
  T1. Goal (the register item; the stakes; the audit must be constructive).
  T2. The edge is live: the BAM flow is GENUINELY nonlinear (superposition
      defect measured; zero for the stripped linear control).
  T3. The channel identified (1D field sector): B-response >> kinematic
      floor; G-linear; Phi-clamp control collapses it to the floor.
  T4. The causal front: retarded Phi confines the response to t >= d/c
      (front scaling measured at three speeds; c -> infinity limit
      recovers Newton).
  T5. Equivariance survives retardation: norm exact, continuity residual
      at integrator error, Born ensemble at noise on the kicked retarded
      dynamics.
  T6. The entangled sector (linear): marginal invariance at machine
      precision; equilibrium ensembles signal-local; NON-equilibrium
      ensembles signal (the Valentini boundary, located).
  T7. The entangled sector (nonlinear): the mean-field gravitational
      coupling violates marginal invariance at O(G) instantaneously -
      and retardation confines the violation behind the front.
  T8. Honest scope + assessment.

Verdict:
  NO_SIGNALING_SURVIVES_AUDIT_THE_ONLY_NONLINEAR_CHANNEL_IS_THE_RETARDED_
  GRAVITATIONAL_FIELD_EQUIVARIANCE_INTACT_EQUILIBRIUM_SIGNAL_LOCALITY_AT_
  DBB_GRADE
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

# ========================================================================
# SECTION A - the 1D field sector (the #198 psi-Phi-q structure, two far
# throat-packets, a local unitary kick at A, the response measured at B)
# ========================================================================

_N = 1024
_L = 80.0
_X = np.linspace(-_L / 2, _L / 2, _N, endpoint=False)
_DX = _L / _N
_K = 2.0 * np.pi * np.fft.fftfreq(_N, d=_DX)
_K2 = _K ** 2
_K2NZ = _K2.copy()
_K2NZ[0] = 1.0
_G = 0.02          # Newton constant (#198's value)
_GC = 1.0          # density-order coupling g
_A0 = 0.4          # order-field mass a0
_LAMQ = 1.0        # order-field quartic
_KAPPA = 0.5       # order-field stiffness (local, diffusive)
_MU = 0.3          # q gravitates with weight mu
_DT = 5e-4
_T_END = 6.0
_XA, _XB = -15.0, 15.0      # kick site / observation site (separation 30)
_SIG = 2.0
_EPS_K, _W_K, _TAU_K = 0.8, 1.2, 0.5    # the local unitary kick at A
_BMASK = np.abs(_X - _XB) <= 5.0        # Bob's region
_EXP_T = np.exp(-1j * 0.5 * _K2 * _DT)
_VKICK = _EPS_K * np.exp(-(_X - _XA) ** 2 / (2 * _W_K ** 2))
_SNAP = 100                              # record every _SNAP steps
_I3 = 59                                 # snapshot index of t = 3.0
_CS = (8.0, 12.0, 16.0)                  # retardation speeds (front scan)
_C_INF = 20.0                            # the c -> infinity check
_D_GEOM = _XB - 5.0 - _XA                # kick center -> B inner edge = 25

_CACHE: dict = {}


def _poisson(rho: np.ndarray, g: float) -> np.ndarray:
    """Mean-subtracted spectral Poisson Phi'' = 4 pi G rho (the exact
    c -> infinity limit of the wave equation below)."""
    ph = np.fft.fft(rho - rho.mean())
    ph *= 4.0 * np.pi * g / (-_K2NZ)
    ph[0] = 0.0
    return np.real(np.fft.ifft(ph))


def _lap(f: np.ndarray) -> np.ndarray:
    return np.real(np.fft.ifft(-_K2 * np.fft.fft(f)))


def _initial_1d():
    psi = (np.exp(-(_X - _XA) ** 2 / (2 * _SIG ** 2))
           + np.exp(-(_X - _XB) ** 2 / (2 * _SIG ** 2))).astype(complex)
    psi /= math.sqrt(float(np.sum(np.abs(psi) ** 2) * _DX))
    q = 0.3 * (np.exp(-(_X - _XA) ** 2 / 8.0)
               + np.exp(-(_X - _XB) ** 2 / 8.0))
    return psi, q


def _halfkick(psi, Phi, q, t, kicked):
    v = Phi + 0.5 * _GC * q * q
    if kicked and t <= _TAU_K:
        v = v + _VKICK
    return np.exp(-1j * v * _DT / 2.0) * psi


def _qstep(q, psi):
    return q + _DT * (_KAPPA * _lap(q) - (_A0 - _GC * np.abs(psi) ** 2) * q
                      - _LAMQ * q ** 3)


def _d_at_b(ra, rb):
    return math.sqrt(float(np.sum((ra[_BMASK] - rb[_BMASK]) ** 2) * _DX))


def _run_pair(mode: str, g: float, c: Optional[float] = None,
              clamp: bool = False):
    """Evolve base and kicked runs in lockstep; if clamp, also a third
    kicked run forced to use the BASE run's Phi at every half-step (the
    Phi-clamp control: the only difference reaching B must then travel
    through psi/q, not through gravity).  Returns (t_grid, D_kick,
    D_clamp-or-None)."""
    psi_b, q_b = _initial_1d()
    psi_k, q_k = psi_b.copy(), q_b.copy()
    psi_c = psi_b.copy() if clamp else None
    q_c = q_b.copy() if clamp else None
    Phi_b = _poisson(np.abs(psi_b) ** 2 + _MU * q_b ** 2, g)
    Phi_k = Phi_b.copy()
    Pi_b = np.zeros(_N)
    Pi_k = np.zeros(_N)
    nst = int(round(_T_END / _DT))
    ts, dk, dc = [], [], []
    for it in range(1, nst + 1):
        t = it * _DT
        psi_b = _halfkick(psi_b, Phi_b, q_b, t, False)
        psi_k = _halfkick(psi_k, Phi_k, q_k, t, True)
        if clamp:
            psi_c = _halfkick(psi_c, Phi_b, q_c, t, True)
        psi_b = np.fft.ifft(_EXP_T * np.fft.fft(psi_b))
        psi_k = np.fft.ifft(_EXP_T * np.fft.fft(psi_k))
        if clamp:
            psi_c = np.fft.ifft(_EXP_T * np.fft.fft(psi_c))
        q_b = _qstep(q_b, psi_b)
        q_k = _qstep(q_k, psi_k)
        if clamp:
            q_c = _qstep(q_c, psi_c)
        rho_b = np.abs(psi_b) ** 2 + _MU * q_b ** 2
        rho_k = np.abs(psi_k) ** 2 + _MU * q_k ** 2
        if mode == "wave":
            src_b = 4.0 * np.pi * g * (rho_b - rho_b.mean())
            src_k = 4.0 * np.pi * g * (rho_k - rho_k.mean())
            Pi_b = Pi_b + _DT * c * c * (_lap(Phi_b) - src_b)
            Phi_b = Phi_b + _DT * Pi_b
            Pi_k = Pi_k + _DT * c * c * (_lap(Phi_k) - src_k)
            Phi_k = Phi_k + _DT * Pi_k
        else:
            Phi_b = _poisson(rho_b, g)
            Phi_k = _poisson(rho_k, g)
        psi_b = _halfkick(psi_b, Phi_b, q_b, t, False)
        psi_k = _halfkick(psi_k, Phi_k, q_k, t, True)
        if clamp:
            psi_c = _halfkick(psi_c, Phi_b, q_c, t, True)
        if it % _SNAP == 0:
            ts.append(round(t, 3))
            dk.append(_d_at_b(np.abs(psi_k) ** 2, np.abs(psi_b) ** 2))
            if clamp:
                dc.append(_d_at_b(np.abs(psi_c) ** 2, np.abs(psi_b) ** 2))
    return ts, dk, (dc if clamp else None)


def field_audit() -> dict:
    """All 1D field-sector runs (memoized)."""
    if "field" in _CACHE:
        return _CACHE["field"]
    out: dict = {}
    ts, d0, _ = _run_pair("newton", 0.0)
    out["t"] = ts
    out["floor"] = d0
    _, dg, dcl = _run_pair("newton", _G, clamp=True)
    out["newton"] = dg
    out["clamp"] = dcl
    _, dh, _ = _run_pair("newton", _G / 2)
    out["newton_halfG"] = dh
    for c in _CS + (_C_INF,):
        _, dw, _ = _run_pair("wave", _G, c=c)
        out[f"wave_c{c:g}"] = dw
    _CACHE["field"] = out
    return out


def _front_time(ts, ds, theta: float = 1e-9) -> Optional[float]:
    for t, d in zip(ts, ds):
        if d > theta:
            return t
    return None


# ========================================================================
# SECTION B - equivariance on the kicked RETARDED dynamics (the #198
# theorem re-verified where it now matters: real potential + retardation)
# ========================================================================

_NPART = 20000


def _velocity_1d(psi):
    dpsi = np.fft.ifft(1j * _K * np.fft.fft(psi))
    rho = np.abs(psi) ** 2
    return np.imag(np.conj(psi) * dpsi) / (rho + 1e-12)


def _sample_1d(dens, n, rng):
    csum = np.cumsum(dens)
    csum /= csum[-1]
    return np.interp(rng.random(n), csum, _X)


def _ks_1d(sample, dens):
    csum = np.cumsum(dens)
    csum /= csum[-1]
    f = np.interp(np.sort(sample), _X, csum)
    emp = np.arange(1, len(sample) + 1) / len(sample)
    return float(np.max(np.abs(f - emp)))


def equivariance_run(c: float = 12.0) -> dict:
    """Kicked, retarded evolution with a Born throat ensemble transported
    by v = grad S; continuity residual measured on the live flow."""
    if "equi" in _CACHE:
        return _CACHE["equi"]
    rng = np.random.default_rng(7)
    psi, q = _initial_1d()
    Phi = _poisson(np.abs(psi) ** 2 + _MU * q ** 2, _G)
    Pi = np.zeros(_N)
    ens = _sample_1d(np.abs(psi) ** 2, _NPART, rng)
    norm0 = float(np.sum(np.abs(psi) ** 2) * _DX)
    nst = int(round(_T_END / _DT))
    checks = sorted({int(nst * f) for f in (0.25, 0.5, 0.75, 1.0)})
    ks_series = []
    res_worst = 0.0
    for it in range(1, nst + 1):
        t = it * _DT
        v = np.real(_velocity_1d(psi))
        if it % 3000 == 0:
            # forward continuity residual on the live kicked retarded flow
            rho_m = np.abs(psi) ** 2
            ps2, q2, Ph2, Pi2 = psi.copy(), q.copy(), Phi.copy(), Pi.copy()
            ps2 = _halfkick(ps2, Ph2, q2, t, True)
            ps2 = np.fft.ifft(_EXP_T * np.fft.fft(ps2))
            q2 = _qstep(q2, ps2)
            rho2 = np.abs(ps2) ** 2 + _MU * q2 ** 2
            src = 4.0 * np.pi * _G * (rho2 - rho2.mean())
            Pi2 = Pi2 + _DT * c * c * (_lap(Ph2) - src)
            Ph2 = Ph2 + _DT * Pi2
            ps2 = _halfkick(ps2, Ph2, q2, t, True)
            drho = (np.abs(ps2) ** 2 - rho_m) / _DT
            div_j = np.gradient(rho_m * v, _DX)
            r = (math.sqrt(float(np.mean((drho + div_j) ** 2)))
                 / math.sqrt(float(np.mean(rho_m ** 2))))
            res_worst = max(res_worst, r)
        v1 = np.interp(ens, _X, v)
        e1 = ens + _DT * v1
        ens = ens + _DT / 2.0 * (v1 + np.interp(e1, _X, v))
        psi = _halfkick(psi, Phi, q, t, True)
        psi = np.fft.ifft(_EXP_T * np.fft.fft(psi))
        q = _qstep(q, psi)
        rho_tot = np.abs(psi) ** 2 + _MU * q ** 2
        src = 4.0 * np.pi * _G * (rho_tot - rho_tot.mean())
        Pi = Pi + _DT * c * c * (_lap(Phi) - src)
        Phi = Phi + _DT * Pi
        psi = _halfkick(psi, Phi, q, t, True)
        if it in checks:
            ks_series.append(round(_ks_1d(ens, np.abs(psi) ** 2), 4))
    out = {
        "c": c,
        "norm_drift": float(abs(np.sum(np.abs(psi) ** 2) * _DX - norm0)),
        "continuity_residual": res_worst,
        "ks_series": ks_series,
        "noise": 1.0 / math.sqrt(_NPART),
    }
    _CACHE["equi"] = out
    return out


# ========================================================================
# SECTION C - the nonlinearity witness (the edge is live)
# ========================================================================

def nonlinearity_witness() -> dict:
    """Superposition defect ||T[psi_a + psi_b] - T[psi_a] - T[psi_b]|| of
    the flow map T (T = 1 time unit): nonzero for the full BAM dynamics,
    zero for the stripped linear control (G = 0, q = 0)."""
    if "witness" in _CACHE:
        return _CACHE["witness"]

    def evolve(psi, q, g, live_q, tend=1.0):
        Phi = _poisson(np.abs(psi) ** 2 + _MU * q ** 2, g)
        for it in range(1, int(round(tend / _DT)) + 1):
            psi = np.exp(-1j * (Phi + 0.5 * _GC * q * q) * _DT / 2.0) * psi
            psi = np.fft.ifft(_EXP_T * np.fft.fft(psi))
            if live_q:
                q = _qstep(q, psi)
            Phi = _poisson(np.abs(psi) ** 2 + _MU * q ** 2, g)
            psi = np.exp(-1j * (Phi + 0.5 * _GC * q * q) * _DT / 2.0) * psi
        return psi

    pa = np.exp(-(_X - _XA) ** 2 / (2 * _SIG ** 2)).astype(complex)
    pb = np.exp(-(_X - _XB) ** 2 / (2 * _SIG ** 2)).astype(complex)
    qa = 0.3 * np.exp(-(_X - _XA) ** 2 / 8.0)
    qb = 0.3 * np.exp(-(_X - _XB) ** 2 / 8.0)
    z = np.zeros(_N)

    def defect(g, live_q):
        ea = evolve(pa.copy(), qa if live_q else z, g, live_q)
        eb = evolve(pb.copy(), qb if live_q else z, g, live_q)
        eab = evolve(pa + pb, (qa + qb) if live_q else z, g, live_q)
        num = math.sqrt(float(np.sum(np.abs(eab - ea - eb) ** 2) * _DX))
        den = math.sqrt(float(np.sum(np.abs(eab) ** 2) * _DX))
        return num / den

    out = {"defect_bam": defect(_G, True), "defect_linear": defect(0.0, False)}
    _CACHE["witness"] = out
    return out


# ========================================================================
# SECTION D - the entangled sector: two-throat configuration-space state
# (the effective linear-regime description per #198 condition 2), with
# and without the BAM mean-field gravitational dressing
# ========================================================================

_N2 = 192
_L2 = 48.0
_XG = np.linspace(-_L2 / 2, _L2 / 2, _N2, endpoint=False)
_DX2 = _L2 / _N2
_K1D = 2.0 * np.pi * np.fft.fftfreq(_N2, d=_DX2)
_KX1 = _K1D[:, None]
_KX2 = _K1D[None, :]
_DT2 = 5e-4
_EXP_T2 = np.exp(-1j * 0.5 * (_KX1 ** 2 + _KX2 ** 2) * _DT2)
_K21D = _K1D ** 2
_K21DNZ = _K21D.copy()
_K21DNZ[0] = 1.0
_SIG2 = 1.5
_EPSK2, _WK2, _TAUK2 = 1.0, 1.5, 0.4
_VMAX = 60.0
_NP2 = 12000

# T6 (linear, crossing): packets at -+6 with approach speed 2 -> the two
# branches (L,R) and (R,L) overlap in configuration space around t = 3,
# where the dBB nonlocality can act on x2.
_A6, _V6, _T6 = 6.0, 2.0, 4.5
# T7 (nonlinear, static): packets at -+8, kick at -8, observe x2 in B.
_A7, _T7, _G2, _C2 = 8.0, 4.0, 0.05, 6.0
_B7 = (_XG >= 4.0) & (_XG <= 12.0)
_D7_GEOM = 4.0 - (-_A7)          # kick center -> B inner edge = 12


def _pack(x0, v0):
    return (np.exp(-(_XG - x0) ** 2 / (2 * _SIG2 ** 2))
            * np.exp(1j * v0 * _XG))


def _entangled(a, v0):
    fL, fR = _pack(-a, v0), _pack(a, -v0)
    psi = fL[:, None] * fR[None, :] + fR[:, None] * fL[None, :]
    psi /= math.sqrt(float(np.sum(np.abs(psi) ** 2) * _DX2 * _DX2))
    return psi


def _branch_density(a, v0):
    fL, fR = _pack(-a, v0), _pack(a, -v0)
    return np.abs(fL[:, None] * fR[None, :]) ** 2


def _vkick2(a):
    return _EPSK2 * np.exp(-(_XG + a) ** 2 / (2 * _WK2 ** 2))


def _poisson1d(rho, g):
    ph = np.fft.fft(rho - rho.mean())
    ph *= 4.0 * np.pi * g / (-_K21DNZ)
    ph[0] = 0.0
    return np.real(np.fft.ifft(ph))


def _lap1d(f):
    return np.real(np.fft.ifft(-_K21D * np.fft.fft(f)))


def _marginals(psi):
    rho = np.abs(psi) ** 2
    return rho.sum(axis=1) * _DX2, rho.sum(axis=0) * _DX2


def _bilin(f, pts):
    gx = (pts[:, 0] + _L2 / 2) / _DX2
    gy = (pts[:, 1] + _L2 / 2) / _DX2
    i0 = np.floor(gx).astype(int) % _N2
    j0 = np.floor(gy).astype(int) % _N2
    i1 = (i0 + 1) % _N2
    j1 = (j0 + 1) % _N2
    fx = gx - np.floor(gx)
    fy = gy - np.floor(gy)
    return (f[i0, j0] * (1 - fx) * (1 - fy) + f[i1, j0] * fx * (1 - fy)
            + f[i0, j1] * (1 - fx) * fy + f[i1, j1] * fx * fy)


def _sample_2d(dens, n, rng):
    p = (dens / dens.sum()).ravel()
    idx = rng.choice(dens.size, size=n, p=p)
    ii, jj = np.unravel_index(idx, dens.shape)
    pts = np.empty((n, 2))
    pts[:, 0] = _XG[ii] + (rng.random(n) - 0.5) * _DX2
    pts[:, 1] = _XG[jj] + (rng.random(n) - 0.5) * _DX2
    return pts


def _ks_two_sample(a, b):
    a = np.sort(a)
    b = np.sort(b)
    allv = np.concatenate([a, b])
    ca = np.searchsorted(a, allv, side="right") / len(a)
    cb = np.searchsorted(b, allv, side="right") / len(b)
    return float(np.max(np.abs(ca - cb)))


def _evolve_linear_2d(psi, kick, nst, vk, ens_list):
    """Linear 2D flow (separable Hamiltonian + local unitary on x1 only),
    transporting dBB ensembles by v_i = d_i S."""
    for it in range(1, nst + 1):
        t = it * _DT2
        if kick and t <= _TAUK2:
            psi = np.exp(-1j * vk * _DT2 / 2.0)[:, None] * psi
        fh = np.fft.fft2(psi)
        if ens_list:
            d1 = np.fft.ifft2(1j * _KX1 * fh)
            d2 = np.fft.ifft2(1j * _KX2 * fh)
            rho = np.abs(psi) ** 2 + 1e-30
            v1g = np.imag(np.conj(psi) * d1) / rho
            v2g = np.imag(np.conj(psi) * d2) / rho
            np.clip(v1g, -_VMAX, _VMAX, out=v1g)
            np.clip(v2g, -_VMAX, _VMAX, out=v2g)
            for e in ens_list:
                v1a = _bilin(v1g, e)
                v2a = _bilin(v2g, e)
                ep = np.empty_like(e)
                ep[:, 0] = e[:, 0] + _DT2 * v1a
                ep[:, 1] = e[:, 1] + _DT2 * v2a
                v1b = _bilin(v1g, ep)
                v2b = _bilin(v2g, ep)
                e[:, 0] += _DT2 * 0.5 * (v1a + v1b)
                e[:, 1] += _DT2 * 0.5 * (v2a + v2b)
        psi = np.fft.ifft2(_EXP_T2 * fh)
        if kick and t <= _TAUK2:
            psi = np.exp(-1j * vk * _DT2 / 2.0)[:, None] * psi
    return psi


def entangled_linear() -> dict:
    """T6: the linear entangled sector - marginal invariance + the
    equilibrium/non-equilibrium signal-locality boundary (memoized)."""
    if "ent_lin" in _CACHE:
        return _CACHE["ent_lin"]
    rng = np.random.default_rng(7)
    psi0 = _entangled(_A6, _V6)
    vk = _vkick2(_A6)
    eq_b = _sample_2d(np.abs(psi0) ** 2, _NP2, rng)
    eq_k = eq_b.copy()
    ne_b = _sample_2d(_branch_density(_A6, _V6), _NP2, rng)
    ne_k = ne_b.copy()
    nst = int(round(_T6 / _DT2))
    psiB = _evolve_linear_2d(psi0.copy(), False, nst, vk, [eq_b, ne_b])
    psiK = _evolve_linear_2d(psi0.copy(), True, nst, vk, [eq_k, ne_k])
    _, m2B = _marginals(psiB)
    _, m2K = _marginals(psiK)
    out = {
        "marginal_invariance": float(np.max(np.abs(m2K - m2B))),
        "ks_equilibrium": _ks_two_sample(eq_k[:, 1], eq_b[:, 1]),
        "ks_nonequilibrium": _ks_two_sample(ne_k[:, 1], ne_b[:, 1]),
        "noise_two_sample": 1.36 * math.sqrt(2.0 / _NP2),
        "n_particles": _NP2,
    }
    _CACHE["ent_lin"] = out
    return out


def _evolve_meanfield_2d(mode, g, kick, nst, c=None):
    """The BAM mean-field dressing of the entangled pair: Phi(x) sourced
    by the PHYSICAL-space density rho1(x) + rho2(x); V = Phi(x1)+Phi(x2);
    a local unitary kick on x1 only.  Newton (instantaneous Poisson) or
    retarded wave-equation Phi.  Returns snapshots of the x2-marginal."""
    psi = _entangled(_A7, 0.0)
    vk = _vkick2(_A7)
    r1, r2 = _marginals(psi)
    Phi = _poisson1d(r1 + r2, g)
    Pi = np.zeros(_N2)
    snaps = []
    for it in range(1, nst + 1):
        t = it * _DT2
        V = Phi[:, None] + Phi[None, :]
        if kick and t <= _TAUK2:
            V = V + vk[:, None]
        psi = np.exp(-1j * V * _DT2 / 2.0) * psi
        psi = np.fft.ifft2(_EXP_T2 * np.fft.fft2(psi))
        r1, r2 = _marginals(psi)
        if mode == "wave":
            src = r1 + r2
            src = 4.0 * np.pi * g * (src - src.mean())
            Pi = Pi + _DT2 * c * c * (_lap1d(Phi) - src)
            Phi = Phi + _DT2 * Pi
        else:
            Phi = _poisson1d(r1 + r2, g)
        V = Phi[:, None] + Phi[None, :]
        if kick and t <= _TAUK2:
            V = V + vk[:, None]
        psi = np.exp(-1j * V * _DT2 / 2.0) * psi
        if it % 200 == 0:
            snaps.append(_marginals(psi)[1])
    return snaps


def entangled_nonlinear() -> dict:
    """T7: the nonlinear entangled sector (memoized)."""
    if "ent_nl" in _CACHE:
        return _CACHE["ent_nl"]
    nst = int(round(_T7 / _DT2))
    tgrid = [round(0.1 * (i + 1), 2) for i in range(nst // 200)]
    out: dict = {"t": tgrid}
    for label, mode, g, c in (("linear", "newton", 0.0, None),
                              ("newton", "newton", _G2, None),
                              ("wave", "wave", _G2, _C2)):
        sb = _evolve_meanfield_2d(mode, g, False, nst, c=c)
        sk = _evolve_meanfield_2d(mode, g, True, nst, c=c)
        out[label] = [math.sqrt(float(np.sum((a[_B7] - b[_B7]) ** 2) * _DX2))
                      for a, b in zip(sk, sb)]
    _CACHE["ent_nl"] = out
    return out


# ========================================================================
# TESTS
# ========================================================================

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "The register item, faced at its sharpest edge. #198 derived "
            "the Born rule at dBB grade CONDITIONAL on the linear "
            "measurement regime; the #200 register carries the nonlinear "
            "measurement theory as open. The sharpest refutation vector "
            "inside that item is NO-SIGNALING: by Gisin/Polchinski, "
            "deterministic nonlinear quantum evolutions generically "
            "allow superluminal signaling through entangled states - "
            "and BAM's pilot equation is genuinely nonlinear (Phi[rho], "
            "q[rho]). If BAM signals superluminally in equilibrium in "
            "its causal completion, it is refuted by relativity and "
            "experiment. The linear no-signaling theorem does not "
            "protect it; the Gisin proof (which uses the projection "
            "postulate BAM does not have) does not automatically "
            "convict it. The question must be settled BY CONSTRUCTION "
            "on the running dynamics - this probe."
        ),
        "deliverable": "docs/nonlinear_no_signaling_audit.md",
        "executes": "#198 condition 2 / the #200 register item, no-signaling edge",
        "framing": "QFT on the fixed classical throat geometry - not quantum gravity",
        "pass": True,
    }


def test_T2_edge_is_live() -> dict:
    w = nonlinearity_witness()
    genuinely_nonlinear = w["defect_bam"] > 1e-4
    control_linear = w["defect_linear"] < 1e-10
    ok = genuinely_nonlinear and control_linear
    return {
        "name": "T2_edge_is_live",
        "description": (
            "THE EDGE IS LIVE - BAM is genuinely nonlinear, so the "
            "Gisin threat applies to it and the audit has teeth. The "
            "superposition defect of the 1-time-unit flow map, "
            "||T[psi_a+psi_b] - T[psi_a] - T[psi_b]||/||T[psi_a+psi_b]||, "
            f"is {w['defect_bam']:.2e} for the full psi-Phi-q dynamics "
            "(self-consistent gravity + live order field) versus "
            f"{w['defect_linear']:.1e} for the stripped linear control "
            "(G = 0, q = 0 - machine zero). BAM cannot appeal to the "
            "linear no-signaling theorem: whether its nonlinearity "
            "signals must be measured, channel by channel."
        ),
        "superposition_defect_bam": float(f"{w['defect_bam']:.3e}"),
        "superposition_defect_linear_control": float(f"{w['defect_linear']:.2e}"),
        "pass": ok,
    }


def test_T3_channel_identified() -> dict:
    f = field_audit()
    t3 = f["t"][_I3]
    floor3 = f["floor"][_I3]
    newt3 = f["newton"][_I3]
    half3 = f["newton_halfG"][_I3]
    clamp3 = f["clamp"][_I3]
    g_ratio = newt3 / half3
    above_floor = newt3 > 1e-7 and newt3 / max(floor3, 1e-16) > 1e4
    g_linear = 1.7 < g_ratio < 2.4
    clamp_kills = clamp3 < 1e-10 and clamp3 / newt3 < 1e-4
    ok = abs(t3 - 3.0) < 1e-9 and above_floor and g_linear and clamp_kills
    return {
        "name": "T3_channel_identified",
        "description": (
            "THE CHANNEL IS EXACTLY THE GRAVITATIONAL FIELD - measured "
            "three ways on the 1D field sector (two throat-packets "
            "separated by 30; a local unitary kick at A for t <= 0.5; "
            "the density response D(t) = L2 norm of rho_kick - rho_base "
            "over Bob's region, at t = 3). (a) THE RESPONSE IS REAL: "
            f"D = {newt3:.2e} with instantaneous (Poisson) gravity vs "
            f"the kinematic (G = 0) floor {floor3:.2e} - a factor "
            f"{newt3/floor3:.1e}: the nonlinearity DOES transmit a "
            "local operation to a distant region, instantaneously, in "
            "the Newtonian model. Stated plainly. (b) IT IS O(G): "
            f"halving G gives D = {half3:.2e}, ratio {g_ratio:.3f} "
            "(linear response in the coupling). (c) NOTHING BUT PHI "
            "CARRIES IT: clamping Phi to the no-kick history (same "
            "kick, same live q, gravity blinded) collapses the response "
            f"to D = {clamp3:.2e} - back AT the kinematic floor, a "
            f"suppression of {newt3/clamp3:.1e}. The order field q is "
            "local (diffusive, kappa = 0.5) and carries nothing to B. "
            "One channel, at gravitational strength - a physical "
            "interaction, not a measurement-theoretic pathology."
        ),
        "t_measure": t3,
        "D_floor": float(f"{floor3:.3e}"),
        "D_newton": float(f"{newt3:.3e}"),
        "D_newton_halfG": float(f"{half3:.3e}"),
        "D_clamp": float(f"{clamp3:.3e}"),
        "G_ratio": round(g_ratio, 3),
        "clamp_suppression": float(f"{newt3/clamp3:.2e}"),
        "pass": ok,
    }


def test_T4_causal_front() -> dict:
    f = field_audit()
    ts = f["t"]
    fronts = {}
    pre_quiet = {}
    for c in _CS:
        ds = f[f"wave_c{c:g}"]
        tf = _front_time(ts, ds)
        fronts[c] = tf
        if tf is not None:
            ipre = max(0, int(round(tf / (_SNAP * _DT))) - 1
                       - int(round(0.5 / (_SNAP * _DT))))
            pre_quiet[c] = ds[ipre]
    all_found = all(fronts[c] is not None for c in _CS)
    cts = {c: (c * fronts[c] if fronts[c] else None) for c in _CS}
    scaling_ok = all_found and all(20.0 < cts[c] < 32.0 for c in _CS)
    monotone = all_found and fronts[8.0] > fronts[12.0] > fronts[16.0]
    quiet_ok = all(v < 1e-11 for v in pre_quiet.values())
    cinf3 = f[f"wave_c{_C_INF:g}"][_I3]
    newt3 = f["newton"][_I3]
    conv = cinf3 / newt3
    conv_ok = 0.5 < conv < 2.0
    ok = scaling_ok and monotone and quiet_ok and conv_ok
    return {
        "name": "T4_causal_front",
        "description": (
            "RETARDATION CONFINES THE CHANNEL TO THE CONE. Replacing "
            "the instantaneous Poisson equation by the causal wave "
            "equation Box Phi = 4 pi G rho (propagation speed c - the "
            "Newtonian model is its c -> infinity limit, recovered "
            "exactly by the same mean-subtracted spectral solver), the "
            "kick response at B stays at the MACHINE floor (~1e-15, "
            f"pre-front values { {f'{c:g}': float(f'{v:.1e}') for c, v in pre_quiet.items()} }) "
            "until the front arrives, then rises to the instantaneous "
            f"amplitude. Front times: { {f'{c:g}': fronts[c] for c in _CS} } "
            f"giving c*t_front = { {f'{c:g}': round(cts[c], 1) for c in _CS} } "
            f"versus the geometric distance {_D_GEOM:.0f} (kick center "
            "to Bob's near edge) - the front scales as d/c across a "
            "factor-2 speed sweep, and is monotone in c. The "
            f"c -> infinity check: at t = 3, D(c = {_C_INF:g}) / "
            f"D(Newton) = {conv:.2f} (the retarded amplitude recovers "
            "the instantaneous one, modulo wave-transient ringing). The "
            "apparent superluminal signaling of T3 is the "
            "action-at-a-distance OF THE NEWTONIAN APPROXIMATION, not "
            "of the theory: in the causal completion the channel is an "
            "ordinary retarded interaction."
        ),
        "front_times": {f"{c:g}": fronts[c] for c in _CS},
        "c_times_t": {f"{c:g}": round(cts[c], 1) for c in _CS},
        "geometric_distance": _D_GEOM,
        "pre_front_response": {f"{c:g}": float(f"{v:.2e}")
                               for c, v in pre_quiet.items()},
        "c_inf_over_newton_at_t3": round(conv, 3),
        "pass": ok,
    }


def test_T5_equivariance_survives() -> dict:
    e = equivariance_run()
    norm_ok = e["norm_drift"] < 1e-9
    res_ok = e["continuity_residual"] < 2e-3
    ks_ok = all(k <= 3.0 * e["noise"] for k in e["ks_series"])
    ok = norm_ok and res_ok and ks_ok
    return {
        "name": "T5_equivariance_survives_retardation",
        "description": (
            "NO-SIGNALING AND THE BORN RULE HOLD SIMULTANEOUSLY. The "
            "#198 equivariance theorem needed only the REALITY of the "
            "potentials - and the retarded Phi is exactly as real as "
            "the instantaneous one, so making gravity causal costs "
            "nothing. Verified on the kicked, retarded (c = "
            f"{e['c']:g}) live dynamics: norm drift "
            f"{e['norm_drift']:.1e} (unitary - the kick and the "
            f"retarded Phi are real potentials); continuity residual "
            f"d(rho)/dt + div(rho grad S) = {e['continuity_residual']:.1e} "
            f"(integrator error); a {_NPART}-throat Born ensemble "
            "transported by v = grad S stays at sampling noise through "
            f"the kicked retarded evolution (KS series {e['ks_series']} "
            f"vs noise {e['noise']:.4f}). The causal completion that "
            "removes the signaling PRESERVES the Born rule - the two "
            "halves of the audit are compatible, structurally."
        ),
        "norm_drift": float(f"{e['norm_drift']:.2e}"),
        "continuity_residual": float(f"{e['continuity_residual']:.2e}"),
        "ks_series": e["ks_series"],
        "sampling_noise": round(e["noise"], 4),
        "pass": ok,
    }


def test_T6_entangled_linear() -> dict:
    r = entangled_linear()
    noise = r["noise_two_sample"]
    minv_ok = r["marginal_invariance"] < 1e-12
    eq_ok = r["ks_equilibrium"] <= 1.5 * noise
    neq_ok = (r["ks_nonequilibrium"] >= 2.2 * noise
              and r["ks_nonequilibrium"] >= 2.0 * r["ks_equilibrium"])
    ok = minv_ok and eq_ok and neq_ok
    return {
        "name": "T6_entangled_linear_sector",
        "description": (
            "THE ENTANGLED SECTOR, LINEAR PART - where Gisin's theorem "
            "actually lives. The two-throat state psi(x1,x2) = "
            "phi_L(x1)phi_R(x2) + phi_R(x1)phi_L(x2) (the #198 "
            "effective linear-regime description; packets approach and "
            "the branches overlap in configuration space around t = 3). "
            "(a) THE THEOREM, ON THE DISCRETE FLOW: a local unitary on "
            "throat 1 (real potential pulse on x1 only) leaves the "
            "x2-marginal invariant to max|delta rho_2| = "
            f"{r['marginal_invariance']:.1e} - machine precision, the "
            "exact factorization Tr_1[(K x 1) rho (K x 1)^dag] = rho_2. "
            "(b) EQUILIBRIUM ENSEMBLES ARE SIGNAL-LOCAL: "
            f"{r['n_particles']} throat pairs sampled from |psi|^2 and "
            "transported by the dBB flow show NO kick dependence in the "
            f"x2-marginal - two-sample KS {r['ks_equilibrium']:.4f} vs "
            f"noise {noise:.4f}. (c) NON-EQUILIBRIUM ENSEMBLES SIGNAL: "
            "the same trajectories with a branch-only (non-Born) "
            f"preparation give KS {r['ks_nonequilibrium']:.4f} - "
            f"{r['ks_nonequilibrium']/noise:.1f}x noise: Alice's local "
            "choice IS visible at B out of equilibrium. This is "
            "Valentini's signal-locality boundary reproduced on the BAM "
            "transport: no-signaling is an EQUILIBRIUM property, with "
            "the #198 H-theorem relaxation as the mechanism that "
            "attains it - exactly dBB's epistemic position, no weaker "
            "and no stronger."
        ),
        "marginal_invariance": float(f"{r['marginal_invariance']:.2e}"),
        "ks_equilibrium": round(r["ks_equilibrium"], 4),
        "ks_nonequilibrium": round(r["ks_nonequilibrium"], 4),
        "noise_two_sample": round(noise, 4),
        "pass": ok,
    }


def test_T7_entangled_nonlinear() -> dict:
    r = entangled_nonlinear()
    ts = r["t"]
    lin_max = max(r["linear"])
    i10 = ts.index(1.0)
    iend = len(ts) - 1
    newt_front = _front_time(ts, r["newton"])
    wave_front = _front_time(ts, r["wave"])
    lin_ok = lin_max < 1e-12
    newt_ok = (newt_front is not None and newt_front <= 0.8
               and r["newton"][iend] > 1e-7)
    wave_pre = r["wave"][i10]
    wave_ok = (wave_front is not None and 1.3 <= wave_front <= 2.8
               and wave_pre < 1e-12 and r["wave"][iend] > 1e-7)
    ok = lin_ok and newt_ok and wave_ok
    return {
        "name": "T7_entangled_nonlinear_channel",
        "description": (
            "THE GISIN CHANNEL, EXHIBITED AND CONFINED. Dressing the "
            "entangled pair with the BAM mean-field gravity - Phi(x) "
            "sourced by the physical-space density rho_1(x) + rho_2(x), "
            "V = Phi(x_1) + Phi(x_2) - and kicking throat 1 at x = "
            f"{-_A7:g} (a local unitary on x1 ONLY, so the linear "
            "theorem makes the background exactly zero). The "
            "x2-marginal response in Bob's region [4, 12]: LINEAR "
            f"CONTROL: {lin_max:.1e} at all times - machine zero, the "
            "theorem. INSTANTANEOUS (Newtonian) MEAN FIELD: nonzero "
            f"from t = {newt_front} (immediately, O(G)) rising to "
            f"{r['newton'][iend]:.1e} - the marginal-invariance theorem "
            "IS violated by the nonlinearity, exactly as Gisin says it "
            "generically must be. Stated plainly. RETARDED MEAN FIELD "
            f"(c = {_C2:g}): the response stays at machine floor "
            f"({wave_pre:.1e} at t = 1.0) until t = {wave_front} - "
            f"versus the geometric d/c = {_D7_GEOM/_C2:.1f} - then "
            f"rises to {r['wave'][iend]:.1e}. The violation of the "
            "linear no-signaling theorem is real, is exactly the "
            "gravitational coupling, and travels ON THE CONE in the "
            "causal completion: at the entangled level too, BAM's "
            "nonlinear 'signaling' is ordinary retarded interaction."
        ),
        "linear_floor_max": float(f"{lin_max:.2e}"),
        "newton_front": newt_front,
        "newton_final": float(f"{r['newton'][iend]:.2e}"),
        "wave_front": wave_front,
        "wave_pre_arrival_t1": float(f"{wave_pre:.2e}"),
        "wave_final": float(f"{r['wave'][iend]:.2e}"),
        "geometric_arrival": round(_D7_GEOM / _C2, 2),
        "pass": ok,
    }


def test_T8_scope_and_assessment() -> dict:
    return {
        "name": "T8_scope_and_assessment",
        "description": (
            "THE AUDIT'S VERDICT, WITH ITS SCOPE. ESTABLISHED: the BAM "
            "nonlinearity opens exactly ONE channel from a local "
            "operation to a distant region - the gravitational field - "
            "measured (O(G), Phi-clamp kills it, q is local); in the "
            "Newtonian model that channel is instantaneous (the "
            "refutation edge was live and FIRED at the approximation "
            "level, stated plainly); in the causal completion (Box Phi "
            "= 4 pi G rho) it is confined to the cone (front = d/c at "
            "three speeds, machine-floor quiet outside, c -> infinity "
            "recovers Newton); the completion costs nothing - the "
            "retarded potential is real, so #198 equivariance, "
            "unitarity, and the Born ensemble survive untouched; and "
            "the entangled sector obeys the same structure (linear "
            "part: exact marginal invariance + equilibrium "
            "signal-locality + the Valentini non-equilibrium signal; "
            "nonlinear part: the O(G) violation confined behind the "
            "front). SCOPE, stated: (1) the wave-equation Phi is the "
            "minimal causal completion - the gauge-fixed weak-field "
            "form of the 5D Einstein equations whose Bianchi structure "
            "#199 verified; the full GR constraint analysis (lapse/"
            "shift, gauge-invariant vs gauge potentials) is not run "
            "here. (2) The configuration-space pair is the #198 "
            "EFFECTIVE description of the linear measurement regime, "
            "now shown to admit a causal gravitational dressing - the "
            "derivation of configuration-space structure from the "
            "single 3-space wave remains the standing open item (the "
            "register item is NARROWED: its no-signaling edge is "
            "audited, its emergence question is not). (3) Equilibrium "
            "is a hypothesis with a mechanism (the #198 relaxation), "
            "not a theorem - non-equilibrium ensembles signal, exactly "
            "as in dBB (Valentini), and this is a PREDICTION shared "
            "with that program, not a defect specific to BAM. (4) 1D/"
            "2D reductions; c is a model parameter (the physical value "
            "is the light speed); the theorems invoked are "
            "dimension-blind. THE LABEL: nonlinear no-signaling "
            "SURVIVES THE AUDIT - superluminal signaling in BAM is an "
            "artifact of the instantaneous-gravity approximation, "
            "removable at zero cost to the Born rule, with "
            "signal-locality holding at dBB grade (in equilibrium) and "
            "the one genuine nonlinear channel identified as ordinary "
            "retarded gravity."
        ),
        "classification": (
            "NO_SIGNALING_SURVIVES_AUDIT_THE_ONLY_NONLINEAR_CHANNEL_IS_"
            "THE_RETARDED_GRAVITATIONAL_FIELD_EQUIVARIANCE_INTACT_"
            "EQUILIBRIUM_SIGNAL_LOCALITY_AT_DBB_GRADE"
        ),
        "pass": True,
    }


# ========================================================================
# RUNNER + VERDICT
# ========================================================================

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_edge_is_live(),
        test_T3_channel_identified(),
        test_T4_causal_front(),
        test_T5_equivariance_survives(),
        test_T6_entangled_linear(),
        test_T7_entangled_nonlinear(),
        test_T8_scope_and_assessment(),
    ]
    t3, t4, t5, t6, t7 = tests[2], tests[3], tests[4], tests[5], tests[6]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "NO_SIGNALING_SURVIVES_AUDIT_THE_ONLY_NONLINEAR_CHANNEL_IS_"
            "THE_RETARDED_GRAVITATIONAL_FIELD_EQUIVARIANCE_INTACT_"
            "EQUILIBRIUM_SIGNAL_LOCALITY_AT_DBB_GRADE"
        )
        verdict = (
            "THE AUDIT RAN, THE EDGE FIRED WHERE IT SHOULD, AND THE "
            "THEORY SURVIVES IT (the argument is in "
            "docs/nonlinear_no_signaling_audit.md; this probe measures "
            "everything on the live dynamics).\n\n"
            "THE CHANNEL. BAM is genuinely nonlinear (superposition "
            "defect measured), so Gisin's threat applies - and indeed a "
            "local kick at A reaches B instantaneously in the Newtonian "
            f"model, {t3['D_newton']/t3['D_floor']:.0e} above the "
            "kinematic floor. But the channel is EXACTLY the "
            f"gravitational field: O(G) (ratio {t3['G_ratio']}), and "
            "clamping Phi collapses it to the floor (suppression "
            f"{t3['clamp_suppression']:.0e}). A physical interaction, "
            "not a measurement pathology.\n\n"
            "THE CONE. With the causal wave-equation Phi, the response "
            "is machine-floor quiet outside the light cone and arrives "
            f"at c*t = {list(t4['c_times_t'].values())} vs geometric "
            f"distance {t4['geometric_distance']:.0f}, at three speeds; "
            "c -> infinity recovers Newton. The superluminality belongs "
            "to the APPROXIMATION, not the theory.\n\n"
            "THE COEXISTENCE. The retarded potential is still real: "
            f"norm exact, continuity residual "
            f"{t5['continuity_residual']:.0e}, Born ensemble at noise "
            "through the kicked retarded evolution - removing the "
            "signaling costs the Born rule nothing.\n\n"
            "THE ENTANGLED SECTOR. Linear part: exact marginal "
            f"invariance ({t6['marginal_invariance']:.0e}); equilibrium "
            f"trajectories signal-local (KS {t6['ks_equilibrium']} ~ "
            f"noise {t6['noise_two_sample']}); non-equilibrium "
            f"trajectories signal (KS {t6['ks_nonequilibrium']}) - the "
            "Valentini boundary on the BAM transport. Nonlinear part: "
            "the mean-field gravitational dressing violates marginal "
            "invariance at O(G) instantaneously (the Gisin channel, "
            "exhibited) - and retardation confines it behind the front "
            f"(arrival {t7['wave_front']} vs geometric "
            f"{t7['geometric_arrival']}; machine floor before). "
            "No-signaling survives at dBB grade: equilibrium + causal "
            "fields, both now measured properties of the BAM dynamics."
        )
    else:
        verdict_class = "NO_SIGNALING_AUDIT_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A channel identification, front, "
            "equivariance, or entangled-sector check failed; re-examine "
            "before quoting the no-signaling status."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The nonlinear no-signaling audit: the only channel opened "
            "by the BAM nonlinearity is the gravitational field (O(G), "
            "Phi-clamp verified), instantaneous only in the Newtonian "
            "approximation, confined to the light cone by the causal "
            "wave-equation completion at zero cost to equivariance; "
            "entangled-sector marginal invariance exact in the linear "
            "part, violated at O(G) by the mean field and confined "
            "behind the front; equilibrium signal-locality at dBB "
            "grade, with the non-equilibrium Valentini signal exhibited"
        ),
        "executes": "#198 condition 2 (the nonlinear measurement theory), no-signaling edge",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The nonlinear no-signaling audit - companion probe (PR #204)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/nonlinear_no_signaling_audit.md` - the "
        "no-signaling audit of the nonlinear BAM psi-Phi-q dynamics. This "
        "probe measures every claim on the live dynamics. *(QFT on the "
        "fixed classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the register item; the Gisin stakes; audit by construction",
        "T2": "the edge is live: BAM genuinely nonlinear (defect measured)",
        "T3": "the channel is gravity: O(G), Phi-clamp kills it, q local",
        "T4": "retardation confines it to the cone (front = d/c, 3 speeds)",
        "T5": "equivariance survives retardation (Born at noise, kicked)",
        "T6": "entangled linear: exact invariance; Valentini boundary",
        "T7": "entangled nonlinear: O(G) violation, confined behind front",
        "T8": "no-signaling survives at dBB grade; scope stated",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    f = field_audit()
    out.append("## The field-sector response at B (t = 3, separation 30)")
    out.append("")
    out.append("| configuration | D_B(3) |")
    out.append("|---|---:|")
    out.append(f"| kinematic floor (G = 0) | {f['floor'][_I3]:.2e} |")
    out.append(f"| Newtonian (instantaneous) gravity | {f['newton'][_I3]:.2e} |")
    out.append(f"| half G | {f['newton_halfG'][_I3]:.2e} |")
    out.append(f"| Phi-clamp control | {f['clamp'][_I3]:.2e} |")
    for c in _CS:
        out.append(f"| retarded, c = {c:g} | {f[f'wave_c{c:g}'][_I3]:.2e} |")
    out.append(f"| retarded, c = {_C_INF:g} (the c to infinity check) "
               f"| {f[f'wave_c{_C_INF:g}'][_I3]:.2e} |")
    out.append("")
    t4 = s["tests"][3]
    out.append(f"Front times {t4['front_times']} -> c*t = {t4['c_times_t']} "
               f"vs geometric distance {t4['geometric_distance']:.0f}.")
    out.append("")
    t6, t7 = s["tests"][5], s["tests"][6]
    out.append("## The entangled sector")
    out.append("")
    out.append(f"- linear marginal invariance: {t6['marginal_invariance']:.1e} "
               "(machine); equilibrium KS "
               f"{t6['ks_equilibrium']} vs noise {t6['noise_two_sample']}; "
               f"non-equilibrium KS {t6['ks_nonequilibrium']} (the "
               "Valentini signal)")
    out.append(f"- nonlinear mean field: linear floor {t7['linear_floor_max']:.1e}; "
               f"Newton front {t7['newton_front']}; retarded front "
               f"{t7['wave_front']} vs geometric {t7['geometric_arrival']} "
               f"(pre-arrival {t7['wave_pre_arrival_t1']:.1e})")
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append("")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return "<array>"
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_nonlinear_no_signaling_audit_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
