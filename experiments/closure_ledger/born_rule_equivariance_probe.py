"""
The Born-rule equivariance test - companion probe (PR #198).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THIS PROBE IS THE APPENDIX, NOT THE ARGUMENT
--------------------------------------------
The deliverable of PR #198 is ``docs/born_rule_equivariance.md`` - item 3
of the theorem-shaped program, the deepest refutation vector: is |psi|^2
the measure preserved by the BAM transport flow (as in Bohmian mechanics,
where it is the unique equivariant density), or does the psi-Phi-q
dynamics equivariantly transport some OTHER functional of the field?

THE TWO THEOREMS (proved in the document; verified here on the running
dynamics)
----------------------------------------------------------------------
1. EQUIVARIANCE.  The real-time Hamiltonian flow of the #179/#180 energy
   functional is  i dpsi/dt = -1/2 lap psi + [Phi +- g/2 q^2] psi  with
   Phi and q REAL - so the polar decomposition gives the EXACT continuity
   equation  d(rho)/dt + div(rho grad S) = 0:  an ensemble of throats
   transported by the BAM velocity v = grad S = J/rho stays
   |psi|^2-distributed forever.  NOT automatic: any imaginary/dissipative
   coupling (including the repo's own imaginary-time RELAXATION flow),
   derivative coupling, or non-potential gravity would add a source and
   destroy it.
2. UNIQUENESS.  For any other density functional h(rho):
   d(h)/dt + div(h v) = (h - rho h') div v,  which vanishes identically
   iff h ~ rho, given the flow is compressible (div v != 0 - verified,
   O(10) during the soliton collision).  Only |psi|^2.

Plus the Valentini relaxation: a uniform (maximally non-Born) ensemble
relaxes toward |psi|^2 (the coarse-grained H-function drops by several
nats through the collision dynamics) - the Born rule is an ATTRACTOR,
not just a fixed point.

THE DEMONSTRATION DYNAMICS
--------------------------
The 1D reduction with the same structure as #179/#180 (standard kinetic,
1D Poisson gravity Phi'' = 4 pi G (rho + mu q^2), real order-field
coupling +g/2 q^2, overdamped real-time q): a two-soliton collision -
strongly compressible, far from stationary.  Split-step Fourier (norm
conserved exactly: the potential is real), Heun transport of 5 parallel
20 000-particle ensembles.

Tests:
  T1. Goal (item 3; the stakes; the document is the argument).
  T2. The flow + the guidance identification (norm exact; Ehrenfest
      d<x>/dt = <grad S>_rho on the nonlinear dynamics).
  T3. THEOREM 1 verified: grid continuity residual ~ integrator error;
      the Born ensemble stays at sampling noise through the collision.
  T4. THEOREM 2 verified: compressibility; the (h - rho h') div v
      identity pointwise; sqrt(rho)- and rho^2-ensembles FAIL.
  T5. Wrong-flow controls: v/2 and frozen transport fail - the
      equivariance is specific to the BAM phase-gradient flow.
  T6. The teeth: a dissipative deformation i gamma W psi produces
      exactly the predicted residual 2 gamma W rho, breaks norm
      conservation and equivariance - the reality of the BAM
      potentials is load-bearing; the refutation vector was live.
  T7. Relaxation: the coarse-grained H-function of a uniform ensemble
      drops by > 2 nats - the Born rule as attractor.
  T8. Assessment + honest scope.

Verdict:
  PSI_SQUARED_IS_THE_UNIQUE_EQUIVARIANT_DENSITY_OF_THE_BAM_TRANSPORT_
  FLOW_BORN_RULE_AT_DBB_GRADE_CONDITIONAL_ON_GUIDANCE_AND_EQUILIBRIUM
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

# -- the 1D BAM analog (the #179/#180 structure) -------------------------
_N = 1024
_L = 40.0
_X = np.linspace(-_L / 2, _L / 2, _N, endpoint=False)
_DX = _L / _N
_K = 2.0 * np.pi * np.fft.fftfreq(_N, d=_DX)
_G = 0.02          # Newton constant (1D Poisson kernel |x-y|)
_GC = 1.0          # density-order coupling g
_A0 = 0.4          # order-field mass a0
_LAMQ = 1.0        # order-field quartic
_KAPPA = 0.5       # order-field stiffness
_MU = 0.3          # q gravitates with weight mu
_DT = 5e-4
_T_END = 10.0
_NPART = 20000
_KER = np.abs(_X[:, None] - _X[None, :])
_EXP_T = np.exp(-1j * 0.5 * _K ** 2 * _DT)


def _phi_of(rho_tot: np.ndarray) -> np.ndarray:
    return 2.0 * np.pi * _G * (_KER @ rho_tot) * _DX


def _lap(f: np.ndarray) -> np.ndarray:
    return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / _DX ** 2


def _veff(psi: np.ndarray, q: np.ndarray) -> np.ndarray:
    return _phi_of(np.abs(psi) ** 2 + _MU * q * q) + 0.5 * _GC * q * q


def _step(psi: np.ndarray, q: np.ndarray, gamma_w: Optional[np.ndarray] = None):
    """One Strang split step of the REAL-TIME Hamiltonian flow (plus an
    optional dissipative deformation i*gammaW*psi for the T6 control)."""
    v = _veff(psi, q)
    half = np.exp(-1j * v * _DT / 2.0)
    if gamma_w is not None:
        half = half * np.exp(-gamma_w * _DT / 2.0)
    psi = half * psi
    psi = np.fft.ifft(_EXP_T * np.fft.fft(psi))
    q = q + _DT * (_KAPPA * _lap(q) - (_A0 - _GC * np.abs(psi) ** 2) * q
                   - _LAMQ * q ** 3)
    v = _veff(psi, q)
    half = np.exp(-1j * v * _DT / 2.0)
    if gamma_w is not None:
        half = half * np.exp(-gamma_w * _DT / 2.0)
    psi = half * psi
    return psi, q


def _velocity(psi: np.ndarray) -> np.ndarray:
    dpsi = np.fft.ifft(1j * _K * np.fft.fft(psi))
    rho = np.abs(psi) ** 2
    return np.imag(np.conj(psi) * dpsi) / (rho + 1e-12)


def _initial_state():
    psi = (np.exp(-(_X + 6.0) ** 2 / 4.0) * np.exp(1j * 0.4 * _X)
           + np.exp(-(_X - 6.0) ** 2 / 4.0) * np.exp(-1j * 0.4 * _X))
    psi = psi.astype(complex)
    psi /= math.sqrt(float(np.sum(np.abs(psi) ** 2) * _DX))
    q = 0.3 * np.exp(-_X ** 2 / 8.0)
    return psi, q


def _sample(dens: np.ndarray, n: int, rng) -> np.ndarray:
    c = np.cumsum(dens)
    c /= c[-1]
    return np.interp(rng.random(n), c, _X)


def _ks(sample: np.ndarray, dens: np.ndarray) -> float:
    c = np.cumsum(dens)
    c /= c[-1]
    f = np.interp(np.sort(sample), _X, c)
    emp = np.arange(1, len(sample) + 1) / len(sample)
    return float(np.max(np.abs(f - emp)))


def _hbar(sample: np.ndarray, dens: np.ndarray, nb: int = 48,
          lo: float = -14.0, hi: float = 14.0) -> float:
    hist, _ = np.histogram(sample, bins=nb, range=(lo, hi), density=True)
    idx = np.clip(((_X - lo) / (hi - lo) * nb).astype(int), 0, nb - 1)
    rb = np.zeros(nb)
    for i, r in zip(idx, dens):
        rb[i] += r
    rb /= rb.sum() * (hi - lo) / nb
    m = (hist > 1e-12) & (rb > 1e-12)
    return float(np.sum(hist[m] * np.log(hist[m] / rb[m])) * (hi - lo) / nb)


_CACHE: dict = {}


def main_run() -> dict:
    """The single main evolution with all diagnostics (memoized)."""
    if "run" in _CACHE:
        return _CACHE["run"]
    rng = np.random.default_rng(7)
    psi, q = _initial_state()
    rho0 = np.abs(psi) ** 2
    ens = {
        "born": _sample(rho0, _NPART, rng),
        "sqrt": _sample(np.sqrt(rho0), _NPART, rng),
        "rho2": _sample(rho0 ** 2, _NPART, rng),
        "half_flow": _sample(rho0, _NPART, rng),
        "frozen": _sample(rho0, _NPART, rng),
        "uniform": rng.uniform(-12.0, 12.0, _NPART),
    }
    norm0 = float(np.sum(rho0) * _DX)
    hbar0 = _hbar(ens["uniform"], rho0)
    nst = int(_T_END / _DT)
    checks = sorted({int(nst * f) for f in (0.25, 0.5, 0.75, 1.0)})
    series = []
    res_worst = 0.0
    dv_max = 0.0
    ident_err = 0.0
    for it in range(1, nst + 1):
        v = _velocity(psi)
        dv_max = max(dv_max, float(np.max(np.abs(np.gradient(v, _DX)))))
        if it % 4000 == 0:
            # grid continuity residual (forward in time) + the uniqueness
            # identity for h = sqrt(rho)
            rho_m = np.abs(psi) ** 2
            psi2, q2 = _step(psi.copy(), q.copy())
            rho_p = np.abs(psi2) ** 2
            drho = (rho_p - rho_m) / _DT
            div_j = np.gradient(rho_m * v, _DX)
            r = (math.sqrt(float(np.mean((drho + div_j) ** 2)))
                 / math.sqrt(float(np.mean(rho_m ** 2))))
            res_worst = max(res_worst, r)
            h = np.sqrt(rho_m + 1e-30)
            lhs = ((np.sqrt(rho_p + 1e-30) - h) / _DT
                   + np.gradient(h * v, _DX))
            rhs = (h - rho_m * 0.5 / h) * np.gradient(v, _DX)
            core = rho_m > 1e-4 * rho_m.max()
            ident_err = max(ident_err, float(
                np.max(np.abs(lhs[core] - rhs[core]))
                / max(float(np.max(np.abs(rhs[core]))), 1e-12)))

        def _move(e, fac=1.0):
            v1 = np.interp(e, _X, v) * fac
            e1 = e + _DT * v1
            return e + _DT / 2.0 * (v1 + np.interp(e1, _X, v) * fac)

        ens["born"] = _move(ens["born"])
        ens["sqrt"] = _move(ens["sqrt"])
        ens["rho2"] = _move(ens["rho2"])
        ens["half_flow"] = _move(ens["half_flow"], 0.5)
        ens["uniform"] = _move(ens["uniform"])
        psi, q = _step(psi, q)
        if it in checks:
            rho = np.abs(psi) ** 2
            series.append({
                "t": round(it * _DT, 2),
                "ks_born": round(_ks(ens["born"], rho), 4),
                "ks_sqrt": round(_ks(ens["sqrt"], np.sqrt(rho)), 4),
                "ks_rho2": round(_ks(ens["rho2"], rho ** 2), 4),
                "ks_half_flow": round(_ks(ens["half_flow"], rho), 4),
                "ks_frozen": round(_ks(ens["frozen"], rho), 4),
                "hbar_uniform": round(_hbar(ens["uniform"], rho), 4),
            })
    out = {
        "series": series,
        "norm0": norm0,
        "norm_end": float(np.sum(np.abs(psi) ** 2) * _DX),
        "hbar0": hbar0,
        "continuity_residual_worst": res_worst,
        "uniqueness_identity_relerr": ident_err,
        "max_abs_div_v": dv_max,
        "sampling_noise": 1.0 / math.sqrt(_NPART),
        "psi_end": psi, "q_end": q,
    }
    _CACHE["run"] = out
    return out


# ========================================================================
# TESTS
# ========================================================================

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "Item 3 of the theorem-shaped program - the deepest "
            "refutation vector. For a classical-field-foundational "
            "theory the one honest route to the Born rule is "
            "EQUIVARIANCE (Duerr-Goldstein-Zanghi): is |psi|^2 the "
            "measure preserved by the BAM transport flow, or does the "
            "psi-Phi-q dynamics equivariantly transport some OTHER "
            "functional? The stakes: |psi|^2 -> the first genuine "
            "foothold on the Born rule (the item every audit flagged as "
            "the deepest import); anything else -> BAM is refuted at "
            "its foundation. The argument is in "
            "docs/born_rule_equivariance.md (two theorems + controls); "
            "this probe verifies every claim on the running dynamics."
        ),
        "deliverable": "docs/born_rule_equivariance.md",
        "framing": "QFT on the fixed classical throat geometry - not quantum gravity",
        "pass": True,
    }


def test_T2_flow_and_guidance() -> dict:
    """Norm conservation (real potentials) + Ehrenfest guidance."""
    r = main_run()
    norm_exact = abs(r["norm_end"] - r["norm0"]) < 1e-9
    # Ehrenfest on a dedicated ASYMMETRIC state (the main run is
    # parity-symmetric, so its centroid is pinned at zero — trivial):
    # a displaced soliton with a phase ramp, centroid genuinely moving.
    psi = (np.exp(-(_X + 3.0) ** 2 / 4.0) * np.exp(1j * 0.35 * _X)).astype(complex)
    psi /= math.sqrt(float(np.sum(np.abs(psi) ** 2) * _DX))
    q = 0.3 * np.exp(-(_X + 3.0) ** 2 / 8.0)
    rho = np.abs(psi) ** 2
    x1 = float(np.sum(_X * rho) * _DX) / float(np.sum(rho) * _DX)
    vmean = float(np.sum(_velocity(psi) * rho) * _DX) / float(np.sum(rho) * _DX)
    nst = 200
    vmid = vmean
    for i in range(nst):
        psi, q = _step(psi, q)
        if i == nst // 2:        # midpoint <v>_rho for the centered ratio
            rr = np.abs(psi) ** 2
            vmid = float(np.sum(_velocity(psi) * rr) * _DX) / float(np.sum(rr) * _DX)
    rho2 = np.abs(psi) ** 2
    x2 = float(np.sum(_X * rho2) * _DX) / float(np.sum(rho2) * _DX)
    dxdt = (x2 - x1) / (nst * _DT)
    ehrenfest = abs(dxdt - vmid) < 2e-3 * abs(vmid)
    ok = norm_exact and ehrenfest
    return {
        "name": "T2_flow_and_guidance",
        "description": (
            "The BAM transport flow. The real-time Hamiltonian flow of "
            "the #179/#180 functional - i dpsi/dt = -1/2 lap psi + "
            "[Phi + g/2 q^2] psi with self-consistent 1D Poisson gravity "
            "and the live real order field - conserves the norm EXACTLY "
            f"(|dN| = {abs(r['norm_end']-r['norm0']):.1e}: the potentials "
            "are real, the split-step is unitary). THE GUIDANCE "
            "IDENTIFICATION: the throat velocity is the local momentum "
            "field v = grad S = J/rho; verified by the exact Ehrenfest "
            "relation on the full NONLINEAR dynamics for a genuinely "
            f"moving (asymmetric) soliton: d<x>/dt = {dxdt:.6f} vs "
            f"<grad S>_rho = {vmid:.6f} - the soliton "
            "centroid follows the phase-gradient flow. Not an extra "
            "postulate: it is also the unique velocity whose current "
            "closes the continuity equation (T3)."
        ),
        "norm_drift": float(f"{abs(r['norm_end']-r['norm0']):.2e}"),
        "ehrenfest_dxdt": round(dxdt, 6),
        "ehrenfest_mean_v": round(vmid, 6),
        "pass": ok,
    }


def test_T3_theorem1_equivariance() -> dict:
    """Continuity residual ~ 0; the Born ensemble stays at noise."""
    r = main_run()
    ks_end = r["series"][-1]["ks_born"]
    ks_all = [s["ks_born"] for s in r["series"]]
    noise = r["sampling_noise"]
    stays = all(k < 3.0 * noise for k in ks_all)
    resid_ok = r["continuity_residual_worst"] < 1e-3
    ok = stays and resid_ok
    return {
        "name": "T3_theorem1_equivariance",
        "description": (
            "THEOREM 1 VERIFIED ON THE RUNNING DYNAMICS. (a) The grid "
            "continuity residual d(rho)/dt + div(rho grad S), evaluated "
            "on the live evolution with the self-consistent Phi and the "
            f"live q, is {r['continuity_residual_worst']:.1e} (relative "
            "L2 - integrator error only): the polar-decomposition "
            "identity holds exactly for the BAM structure because every "
            "coupling is a REAL potential. (b) The ensemble statement: "
            f"{_NPART} throats sampled from |psi(0)|^2 and transported "
            "by v = grad S remain |psi(t)|^2-distributed through 10 time "
            "units INCLUDING a two-soliton collision - KS distance "
            f"{ks_all} vs sampling noise 1/sqrt(N) = {noise:.4f}: never "
            "above 3x noise, ending at "
            f"{ks_end:.4f}. If the throat ensemble is Born-distributed "
            "once, the BAM transport keeps it Born-distributed."
        ),
        "ks_born_series": ks_all,
        "sampling_noise": round(noise, 4),
        "continuity_residual": float(f"{r['continuity_residual_worst']:.2e}"),
        "pass": ok,
    }


def test_T4_theorem2_uniqueness() -> dict:
    """Compressibility; the (h - rho h') div v identity; alternatives
    fail."""
    r = main_run()
    compressible = r["max_abs_div_v"] > 1.0
    ident_ok = r["uniqueness_identity_relerr"] < 0.05
    ks_sqrt_end = r["series"][-1]["ks_sqrt"]
    ks_rho2_end = r["series"][-1]["ks_rho2"]
    ks_born_end = r["series"][-1]["ks_born"]
    alternatives_fail = (ks_sqrt_end > 0.05 and ks_rho2_end > 0.05
                         and ks_sqrt_end > 5 * ks_born_end
                         and ks_rho2_end > 5 * ks_born_end)
    ok = compressible and ident_ok and alternatives_fail
    return {
        "name": "T4_theorem2_uniqueness",
        "description": (
            "THEOREM 2 VERIFIED: ONLY |psi|^2. For any h(rho) the "
            "transported residual is d(h)/dt + div(h v) = "
            "(h - rho h') div v - zero for all configurations iff h ~ "
            "rho, GIVEN compressibility. (a) The flow is strongly "
            f"compressible: max|div v| = {r['max_abs_div_v']:.1f} "
            "through the collision - the uniqueness hypothesis is "
            "non-vacuous. (b) The identity itself is verified pointwise "
            "on the live dynamics for h = sqrt(rho): the directly "
            "computed residual matches (h - rho h') div v to "
            f"{100*r['uniqueness_identity_relerr']:.1f}% (core region). "
            "(c) The ensemble statement: sqrt(rho)- and rho^2-prepared "
            "ensembles transported by the SAME flow depart from their "
            f"functionals immediately - final KS {ks_sqrt_end:.3f} and "
            f"{ks_rho2_end:.3f} vs the Born ensemble's {ks_born_end:.4f} "
            "(>5x separation). Among density functionals, |psi|^2 is "
            "the unique equivariant one (Goldstein-Struyve extend to "
            "all local functionals - cited)."
        ),
        "max_abs_div_v": round(r["max_abs_div_v"], 2),
        "identity_rel_err": round(r["uniqueness_identity_relerr"], 4),
        "ks_sqrt_end": ks_sqrt_end,
        "ks_rho2_end": ks_rho2_end,
        "pass": ok,
    }


def test_T5_wrong_flow_controls() -> dict:
    """The equivariance is specific to the BAM phase-gradient flow."""
    r = main_run()
    ks_half = r["series"][-1]["ks_half_flow"]
    ks_frozen = r["series"][-1]["ks_frozen"]
    ks_born = r["series"][-1]["ks_born"]
    ok = ks_half > 0.05 and ks_frozen > 0.05 and ks_half > 5 * ks_born
    return {
        "name": "T5_wrong_flow_controls",
        "description": (
            "THE FLOW MATTERS. The same Born-prepared ensemble "
            "transported by the WRONG flow loses the Born distribution: "
            f"half-speed transport (v/2) ends at KS = {ks_half:.3f}, a "
            f"frozen ensemble at KS = {ks_frozen:.3f} - versus "
            f"{ks_born:.4f} for the true BAM flow (>5x separation). "
            "Equivariance is a property of the pair (|psi|^2, grad S), "
            "not of |psi|^2 alone: the phase-gradient transport is the "
            "unique flow under which the Born measure is preserved "
            "(the current rho*grad S is the one that closes the "
            "continuity equation, T3a)."
        ),
        "ks_half_flow_end": ks_half,
        "ks_frozen_end": ks_frozen,
        "pass": ok,
    }


def test_T6_teeth_dissipative_deformation() -> dict:
    """The refutation vector was live: breaking the reality of the
    potentials breaks equivariance exactly as predicted."""
    gamma = 0.15
    w = np.exp(-(_X + 4.0) ** 2 / 8.0)   # sits on the left soliton's path
    rng = np.random.default_rng(11)
    psi, q = _initial_state()
    ens = _sample(np.abs(psi) ** 2, 4000, rng)
    # predicted residual: d(rho)/dt + div(rho v) = -2 gamma W rho
    rho_m = np.abs(psi) ** 2
    v0 = _velocity(psi)
    psi2, q2 = _step(psi.copy(), q.copy(), gamma_w=gamma * w)
    rho_p = np.abs(psi2) ** 2
    lhs = (rho_p - rho_m) / _DT + np.gradient(rho_m * v0, _DX)
    rhs = -2.0 * gamma * w * rho_m
    core = rho_m > 1e-4 * rho_m.max()
    resid_match = float(np.max(np.abs(lhs[core] - rhs[core]))
                        / np.max(np.abs(rhs[core])))
    # evolve 2 time units with the deformation: norm decays, KS grows
    nst = int(2.0 / _DT)
    for _ in range(nst):
        v = _velocity(psi)
        v1 = np.interp(ens, _X, v)
        e1 = ens + _DT * v1
        ens = ens + _DT / 2.0 * (v1 + np.interp(e1, _X, v))
        psi, q = _step(psi, q, gamma_w=gamma * w)
    rho = np.abs(psi) ** 2
    norm_lost = 1.0 - float(np.sum(rho) * _DX)
    ks_broken = _ks(ens, rho)
    ok = resid_match < 0.10 and norm_lost > 0.05 and ks_broken > 0.03
    return {
        "name": "T6_teeth_dissipative_deformation",
        "description": (
            "THE TEETH - the reality of the BAM potentials is "
            "load-bearing, not automatic. Deforming the pilot equation "
            "by a dissipative term i*gamma*W(x)*psi (gamma = 0.15, W on "
            "the soliton's path - the structure of an imaginary-time/"
            "absorbing coupling): (a) the "
            "continuity equation acquires EXACTLY the predicted source "
            "-2 gamma W rho (direct grid computation matches to "
            f"{100*resid_match:.1f}%); (b) the norm decays "
            f"({100*norm_lost:.1f}% over 2 time units) - no functional "
            "of rho can be equivariant; (c) the Born-prepared ensemble "
            f"drifts off |psi|^2 (KS = {ks_broken:.3f}). The repo's own "
            "imaginary-time RELAXATION flow is of exactly this "
            "non-equivariant type (its norm is renormalized by hand "
            "each step) - the solver tool and the physical Hamiltonian "
            "flow are distinct, and only the latter carries the Born "
            "rule. Had the BAM energy functional contained any complex, "
            "dissipative, or derivative coupling, this test would have "
            "refuted the Born rule at the foundation. It does not."
        ),
        "predicted_residual_match_relerr": round(resid_match, 4),
        "norm_lost_2tu": round(norm_lost, 4),
        "ks_born_broken": round(ks_broken, 4),
        "pass": ok,
    }


def test_T7_relaxation() -> dict:
    """The subquantum H-theorem: uniform -> Born."""
    r = main_run()
    h0 = r["hbar0"]
    hs = [s["hbar_uniform"] for s in r["series"]]
    drop = h0 - hs[-1]
    monotone_trend = hs[0] < h0 and hs[-1] <= min(hs) + 0.5
    ok = drop > 2.0 and monotone_trend
    return {
        "name": "T7_relaxation",
        "description": (
            "THE BORN RULE AS ATTRACTOR (Valentini's subquantum "
            "H-theorem, on the BAM dynamics). A UNIFORM ensemble - "
            "maximally non-Born - transported by the same flow through "
            "the soliton-collision dynamics relaxes toward |psi|^2: the "
            "coarse-grained H-function Hbar = integral fbar ln(fbar/"
            f"rhobar) falls from {h0:.2f} to {hs[-1]:.2f} (a drop of "
            f"{drop:.2f} nats; series {hs}). Equivariance makes Born a "
            "fixed point; relaxation makes it an attractor of "
            "coarse-grained non-equilibrium - the two halves of the "
            "dBB-grade Born rule, both now demonstrated on the BAM "
            "transport."
        ),
        "hbar_initial": round(h0, 3),
        "hbar_series": hs,
        "drop_nats": round(drop, 3),
        "pass": ok,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "THE FOOTHOLD, WITH ITS CONDITIONS. Established: |psi|^2 is "
            "EXACTLY equivariant under the real-time BAM psi-Phi-q "
            "transport (structural - the reality of Phi and the "
            "q-coupling; residual at integrator error on the live "
            "dynamics); it is the UNIQUE equivariant density among "
            "functionals of rho (compressible flow; alternatives fail "
            "demonstrably); non-equilibrium relaxes toward it "
            "(H-theorem); and the property is falsifiable - the "
            "dissipative deformation breaks it exactly as predicted, so "
            "the refutation vector was live and did not fire. This is "
            "the Born rule at dBB GRADE - the same epistemic status it "
            "has in Bohmian mechanics, no weaker and no stronger. The "
            "bell-module import is now labeled: derived at dBB grade, "
            "CONDITIONAL on (1) the guidance identification (throat "
            "velocity = grad S - motivated by Galilean structure, "
            "Ehrenfest, and uniqueness of the closing current; its "
            "derivation from 5D throat motion is its own program), and "
            "(2) the measurement regime (test throat in an external "
            "pilot wave, where the nonlinear pilot equation "
            "linearizes; the general nonlinear measurement theory is "
            "open). 1D demonstration; the theorems are dimension-blind "
            "and apply verbatim to the radial 3D #180 system."
        ),
        "classification": (
            "PSI_SQUARED_IS_THE_UNIQUE_EQUIVARIANT_DENSITY_OF_THE_BAM_"
            "TRANSPORT_FLOW_BORN_RULE_AT_DBB_GRADE_CONDITIONAL_ON_"
            "GUIDANCE_AND_EQUILIBRIUM"
        ),
        "pass": True,
    }


# ========================================================================
# RUNNER + VERDICT
# ========================================================================

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_flow_and_guidance(),
        test_T3_theorem1_equivariance(),
        test_T4_theorem2_uniqueness(),
        test_T5_wrong_flow_controls(),
        test_T6_teeth_dissipative_deformation(),
        test_T7_relaxation(),
        test_T8_assessment(),
    ]
    t3, t4, t7 = tests[2], tests[3], tests[6]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "PSI_SQUARED_IS_THE_UNIQUE_EQUIVARIANT_DENSITY_OF_THE_BAM_"
            "TRANSPORT_FLOW_BORN_RULE_AT_DBB_GRADE_CONDITIONAL_ON_"
            "GUIDANCE_AND_EQUILIBRIUM"
        )
        verdict = (
            "THE DEEPEST REFUTATION VECTOR DID NOT FIRE (the argument "
            "is in docs/born_rule_equivariance.md; this probe verifies "
            "it on the running dynamics).\n\n"
            "EQUIVARIANCE. The real-time Hamiltonian flow of the BAM "
            "psi-Phi-q functional preserves |psi|^2 EXACTLY: the "
            "continuity residual on the live evolution (self-consistent "
            f"gravity + live order field) is "
            f"{t3['continuity_residual']:.0e}, and a 20000-throat "
            "Born ensemble transported by v = grad S stays at sampling "
            "noise through a two-soliton collision (KS series "
            f"{t3['ks_born_series']}, noise {t3['sampling_noise']}). "
            "The reason is structural: every BAM coupling - the "
            "self-consistent Phi, the order field q - enters the pilot "
            "equation as a REAL potential.\n\n"
            "UNIQUENESS AND FALSIFIABILITY. Among density functionals "
            "only |psi|^2 is equivariant ((h - rho h') div v with "
            f"max|div v| = {t4['max_abs_div_v']}; sqrt(rho)- and "
            "rho^2-ensembles fail, wrong flows fail), and the property "
            "is falsifiable: a dissipative deformation produces exactly "
            "the predicted continuity source and breaks the transport - "
            "the repo's own imaginary-time relaxation flow is of that "
            "non-equivariant type; only the Hamiltonian flow carries "
            "the Born rule. RELAXATION: a uniform ensemble falls "
            f"{t7['drop_nats']} nats toward Born - fixed point AND "
            "attractor.\n\n"
            "THE LABEL. The Born rule enters BAM at dBB GRADE - "
            "equivariance + uniqueness + relaxation, the same standing "
            "it has in Bohmian mechanics - conditional on the guidance "
            "identification (throat velocity = phase gradient; "
            "Ehrenfest-verified, uniquely current-closing, but its 5D "
            "derivation is its own program) and on the linear "
            "measurement regime (test throat in an external pilot "
            "wave). The deepest import is replaced by a theorem with "
            "stated hypotheses."
        )
    else:
        verdict_class = "BORN_RULE_EQUIVARIANCE_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A verification failed; re-examine the "
            "document's theorems against the dynamics before quoting "
            "the Born-rule status."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The Born-rule equivariance test: |psi|^2 is the unique "
            "equivariant density of the BAM psi-Phi-q transport flow "
            "(exact continuity from the reality of the potentials; "
            "alternatives and wrong flows fail; dissipative deformation "
            "refutes as predicted; uniform ensembles relax toward Born) "
            "- the Born rule at dBB grade, conditions stated"
        ),
        "answers": "item 3 - the deepest refutation vector (the flagged Born-rule import)",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The Born-rule equivariance test - companion probe (PR #198)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/born_rule_equivariance.md` - the "
        "equivariance and uniqueness theorems for the BAM psi-Phi-q "
        "transport flow. This probe verifies every claim on the running "
        "dynamics. *(QFT on the fixed classical throat geometry, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "item 3: the deepest refutation vector; the stakes",
        "T2": "the flow: norm exact; guidance = grad S (Ehrenfest)",
        "T3": "THEOREM 1: |psi|^2 equivariant (residual ~ 0; KS at noise)",
        "T4": "THEOREM 2: unique among h(rho); alternatives fail",
        "T5": "wrong flows fail: equivariance is specific to grad S",
        "T6": "teeth: dissipative deformation breaks it as predicted",
        "T7": "relaxation: uniform -> Born (H-theorem on BAM dynamics)",
        "T8": "Born rule at dBB grade; conditions stated",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    t3 = s["tests"][2]
    r = main_run()
    out.append("## The ensemble transport through the collision")
    out.append("")
    out.append("| t | KS(Born) | KS(sqrt-rho) | KS(rho^2) | KS(v/2) | KS(frozen) | Hbar(uniform) |")
    out.append("|---:|---:|---:|---:|---:|---:|---:|")
    for srow in r["series"]:
        out.append(f"| {srow['t']} | {srow['ks_born']} | {srow['ks_sqrt']} | "
                   f"{srow['ks_rho2']} | {srow['ks_half_flow']} | "
                   f"{srow['ks_frozen']} | {srow['hbar_uniform']} |")
    out.append("")
    out.append(f"(sampling noise 1/sqrt(N) = {r['sampling_noise']:.4f}; "
               f"initial Hbar(uniform) = {r['hbar0']:.2f}; "
               f"continuity residual {r['continuity_residual_worst']:.1e}; "
               f"max|div v| = {r['max_abs_div_v']:.1f})")
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
    # strip the non-serializable field arrays before dumping
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_born_rule_equivariance_probe"
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
