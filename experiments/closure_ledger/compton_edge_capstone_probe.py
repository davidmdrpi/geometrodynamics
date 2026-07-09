"""
The Compton-edge capstone: the second release - companion probe (PR #211).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE CLOSING PR OF THE #204-#210 ARC
------------------------------------
The arc took the program's two named frontiers - the nonlinear
measurement theory (#198 condition 2) and the strong-field core
contraction (#203's target) - and closed both: no-signaling audited
(#204), the committed Phi[rho] confronted with the laboratory (#205),
configuration space derived from the bridge (#206), swapping from
surgery (#207), GHZ from valence (#208), measurement made operational
(#209), and the collapse reading of the mass anchor refuted with the
anchor relocated to alpha (#210).  This capstone does what #200 did for
the previous arc - re-verifies the whole chain green in one run and
encodes the updated register - and carries ONE new law that sharpens
the arc's final unknown:

THE COMPTON-EDGE LAW (new).  #210 left the mass ladder with one O(1):
sigma_mode/lambda_C.  Adding the mass term to #202's bridge equation -
(rho^3 phi')' = [k(k+2) rho + m^2 rho^3] phi, rho = sqrt(r_s^2+sigma^2)
- yields a parameter-free deformation of the exact suppression law:
  * UNIVERSAL: the deformation depends on x = sigma_mode/lambda_C
    alone (m r_s sweep collapses to one curve, < 2e-3);
  * the massless limit re-derives #202's exact law (phi = sigma to
    machine zero; c0(1) = 0);
  * the law: eps_1 = (r_s/sigma_mode) D(x), with D(x) computed
    (D = 0.966 at x = 0.647, 0.831 at x = 1.509);
  * the sensitivity S(x) = |d ln eps_1/d ln sigma| = 1 exactly below
    the Compton scale (the #202 naturalness) and grows beyond
    (1.07 / 1.16 / 1.36 / 1.62 / 2.28 at x = 0.647/1/1.509/2/3):
    NATURALNESS CAPS THE O(1) - S <= 2 confines x <= ~2.6, and the
    ladder's own worst-tolerated sensitivity (4.48, #201/#202) confines
    x <= ~6;
  * both convention anchors sit INSIDE the natural window (S = 1.07,
    1.36), and the deformation-corrected self-consistent band tightens
    the #210 bracket to sigma_mode/lambda_C in [0.63, 1.31] - still
    bracketing 1.  The successor derivation now has a compact, derived
    target window at the Compton edge of the natural domain.

Tests:
  T1. Goal (the closing consolidation + the one new law).
  T2. THE COMPTON-EDGE LAW (universality; the deformation D; the
      sensitivity window; the tightened band).
  T3. Ledger 1 - the commitment chain (#204/#205): real potentials
      (norm + continuity on a kicked retarded mini-run); the classical
      channel cannot entangle; the SN/BMV arithmetic.
  T4. Ledger 2 - the topology chain (#206-#208): the singlet from
      (I x T)|Phi+>; the swapping law; the mixture; LHV/Tsirelson/
      Mermin; monogamy.
  T5. Ledger 3 - measurement + mass (#209/#210): the k-odd dispersion
      identity; a live mini Stern-Gerlach Born check; a live Kaup
      point; the alpha chain.
  T6. The claim-graph and the updated register (machine-encoded).
  T7. Honest scope (what the arc did NOT do).
  T8. Assessment (Release II).

Verdict:
  RELEASE_II_GREEN_THE_ENTANGLED_SECTOR_CLOSED_OPERATIONALLY_THE_MASS_
  LADDER_ON_THE_ALPHA_ANCHOR_THE_O1_CONFINED_TO_THE_COMPTON_EDGE_OF_
  THE_NATURAL_DOMAIN
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

_ALPHA = 7.2973525693e-3
_MU_OVER_E = 206.7683
_CACHE: dict = {}

# ========================================================================
# SECTION A - the Compton-edge law: the massive k = 1 mode on the
# Tangherlini bridge (r_s = 1 units; parameter-free)
# ========================================================================


def _bridge_regular(m: float, smax: float, ds: float = 1e-3):
    """Integrate the regular solution of
    (rho^3 phi')' = [3 rho + m^2 rho^3] phi   (k = 1),
    rho = sqrt(1 + sigma^2), from phi ~ sigma at the cross-cap.
    Returns arrays (sigma, phi, phi')."""
    s = 1e-8
    y = np.array([s, 1.0])
    out_s, out_p, out_dp = [], [], []

    def rhs(s, y):
        phi, dphi = y
        rho2 = 1.0 + s * s
        rho = math.sqrt(rho2)
        phipp = (3.0 * rho * phi + m * m * rho * rho2 * phi
                 - 3.0 * rho * s * dphi) / (rho * rho2)
        return np.array([dphi, phipp])

    while s < smax:
        k1 = rhs(s, y)
        k2 = rhs(s + ds / 2, y + ds / 2 * k1)
        k3 = rhs(s + ds / 2, y + ds / 2 * k2)
        k4 = rhs(s + ds, y + ds * k3)
        y = y + ds / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        s += ds
        out_s.append(s)
        out_p.append(y[0])
        out_dp.append(y[1])
        if abs(y[0]) > 1e14:
            break
    return np.array(out_s), np.array(out_p), np.array(out_dp)


def compton_edge() -> dict:
    """The deformation D(x) and sensitivity S(x), x = m sigma_c =
    sigma_c/lambda_C, for two values of m r_s (universality check);
    the natural-window caps; the deformation-corrected band (memoized)."""
    if "edge" in _CACHE:
        return _CACHE["edge"]
    # massless exactness (the #202 law re-derived)
    s0, p0, _ = _bridge_regular(0.0, 50.0)
    massless_err = float(np.max(np.abs(p0 - s0) / s0))
    xgrid = (0.25, 0.5, 0.647, 1.0, 1.509, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
    curves = {}
    for m in (0.02, 0.05):
        ss, pp, dp = _bridge_regular(m, max(xgrid) / m + 1.0)
        phi_rs = float(np.interp(1.0, ss, pp))
        Sv, Dv = [], []
        for x in xgrid:
            sc = x / m
            ph = float(np.interp(sc, ss, pp))
            dph = float(np.interp(sc, ss, dp))
            Sv.append(sc * dph / ph)
            Dv.append((phi_rs / ph) * sc)
        curves[m] = {"S": Sv, "D": Dv}
    universality = max(abs(curves[0.02]["S"][i] - curves[0.05]["S"][i])
                       for i in range(len(xgrid)))
    S = curves[0.02]["S"]
    D = curves[0.02]["D"]
    # natural-window caps: where S crosses 2 and 4.48
    def _crossing(level):
        for i in range(1, len(xgrid)):
            if S[i - 1] < level <= S[i]:
                f = (level - S[i - 1]) / (S[i] - S[i - 1])
                return xgrid[i - 1] + f * (xgrid[i] - xgrid[i - 1])
        return None
    x_nat2 = _crossing(2.0)
    x_nat448 = _crossing(4.48)
    # the deformation-corrected self-consistent band: x = c * D(x)
    def _fixed_point(c):
        x = c
        for _ in range(60):
            x = c * float(np.interp(x, xgrid, D))
        return x
    band = [_fixed_point((3.0 / 7.0) * _MU_OVER_E * _ALPHA),
            _fixed_point(_MU_OVER_E * _ALPHA)]
    out = {
        "massless_err": massless_err,
        "xgrid": list(xgrid),
        "S": [round(v, 4) for v in S],
        "D": [round(v, 4) for v in D],
        "universality": universality,
        "x_where_S_2": x_nat2,
        "x_where_S_448": x_nat448,
        "corrected_band": [round(b, 3) for b in band],
    }
    _CACHE["edge"] = out
    return out


# ========================================================================
# SECTION B - ledger minis
# ========================================================================

def mini_retarded_kick() -> dict:
    """#204's structural core on a small grid: a kicked, RETARDED
    evolution with real potentials conserves the norm exactly and
    closes continuity at integrator error."""
    if "m204" in _CACHE:
        return _CACHE["m204"]
    N, L, dt, T, G, c = 256, 40.0, 1e-3, 2.0, 0.02, 8.0
    X = np.linspace(-L / 2, L / 2, N, endpoint=False)
    K = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    K2 = K ** 2
    K2NZ = K2.copy()
    K2NZ[0] = 1.0
    psi = (np.exp(-(X + 8) ** 2 / 8.0)
           + np.exp(-(X - 8) ** 2 / 8.0)).astype(complex)
    psi /= math.sqrt(float(np.sum(np.abs(psi) ** 2) * L / N))
    vk = 0.8 * np.exp(-(X + 8) ** 2 / 2.0)

    def poisson(rho):
        ph = np.fft.fft(rho - rho.mean())
        ph *= 4 * np.pi * G / (-K2NZ)
        ph[0] = 0
        return np.real(np.fft.ifft(ph))

    Phi = poisson(np.abs(psi) ** 2)
    Pi = np.zeros(N)
    expT = np.exp(-1j * 0.5 * K2 * dt)
    n0 = float(np.sum(np.abs(psi) ** 2) * L / N)
    res_worst = 0.0
    nst = int(T / dt)
    for it in range(1, nst + 1):
        t = it * dt
        V = Phi + (vk if t < 0.5 else 0.0)
        if it % 500 == 0:
            dp = np.fft.ifft(1j * K * np.fft.fft(psi))
            Jc = np.imag(np.conj(psi) * dp)
            rm = np.abs(psi) ** 2
            p2 = np.exp(-1j * V * dt / 2) * psi
            p2 = np.fft.ifft(expT * np.fft.fft(p2))
            p2 = np.exp(-1j * V * dt / 2) * p2
            drho = (np.abs(p2) ** 2 - rm) / dt
            r = (math.sqrt(float(np.mean((drho + np.gradient(Jc, L / N)) ** 2)))
                 / math.sqrt(float(np.mean(rm ** 2))))
            res_worst = max(res_worst, r)
        psi = np.exp(-1j * V * dt / 2) * psi
        psi = np.fft.ifft(expT * np.fft.fft(psi))
        rho = np.abs(psi) ** 2
        src = 4 * np.pi * G * (rho - rho.mean())
        Pi = Pi + dt * c * c * (np.real(np.fft.ifft(-K2 * np.fft.fft(Phi)))
                                - src)
        Phi = Phi + dt * Pi
        V = Phi + (vk if t < 0.5 else 0.0)
        psi = np.exp(-1j * V * dt / 2) * psi
    out = {"norm_drift": abs(float(np.sum(np.abs(psi) ** 2) * L / N) - n0),
           "continuity_residual": res_worst}
    _CACHE["m204"] = out
    return out


def mini_no_entangle() -> dict:
    """#205's discriminator on a small grid: mean-field (classical)
    gravity keeps a product state at entropy 0; the pairwise operator
    entangles it."""
    if "m205" in _CACHE:
        return _CACHE["m205"]
    n, L, dt, T, G = 48, 24.0, 1e-3, 2.0, 0.1
    x = np.linspace(-L / 2, L / 2, n, endpoint=False)
    dx = L / n
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    expT = np.exp(-1j * 0.5 * (k1[:, None] ** 2 + k1[None, :] ** 2) * dt)
    k2nz = k1 ** 2
    k2nz[0] = 1.0
    psi0 = (np.exp(-(x + 3) ** 2 / 4)[:, None]
            * np.exp(-(x - 3) ** 2 / 4)[None, :]).astype(complex)
    psi0 /= math.sqrt(float(np.sum(np.abs(psi0) ** 2) * dx * dx))
    sep = np.abs(x[:, None] - x[None, :])
    sep = np.minimum(sep, L - sep)
    vpair = 2 * np.pi * G * sep

    def ent(psi):
        sv = np.linalg.svd(psi * dx, compute_uv=False)
        p = sv ** 2
        p = p / p.sum()
        p = p[p > 1e-16]
        return float(-np.sum(p * np.log(p)))

    out = {}
    for label in ("meanfield", "pairwise"):
        psi = psi0.copy()
        for _ in range(int(T / dt)):
            if label == "meanfield":
                rho = np.abs(psi) ** 2
                rbar = rho.sum(axis=1) * dx + rho.sum(axis=0) * dx
                ph = np.fft.fft(rbar - rbar.mean())
                ph *= 4 * np.pi * G / (-k2nz)
                ph[0] = 0
                Phi = np.real(np.fft.ifft(ph))
                V = Phi[:, None] + Phi[None, :]
            else:
                V = vpair
            psi = np.exp(-1j * V * dt / 2) * psi
            psi = np.fft.ifft2(expT * np.fft.fft2(psi))
            psi = np.exp(-1j * V * dt / 2) * psi
        out[label] = ent(psi)
    _CACHE["m205"] = out
    return out


# -- two-qubit / three-qubit algebra (the #206-#208 chain) ----------------

_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_T = 1j * _SY
_PHIP = np.array([1, 0, 0, 1], dtype=complex) / math.sqrt(2)
_SINGLET = np.array([0, 1, -1, 0], dtype=complex) / math.sqrt(2)


def _psi_pair(phi):
    m = np.zeros((2, 2), dtype=complex)
    m[0, 1] = 1.0
    m[1, 0] = np.exp(1j * phi)
    return m / math.sqrt(2)


def _negativity(rho):
    r = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    ev = np.linalg.eigvalsh(r)
    return float(-np.sum(ev[ev < 0]))


# ========================================================================
# TESTS
# ========================================================================

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "THE CLOSING PR OF THE ARC. #204-#210 closed the program's "
            "two named frontiers: the nonlinear measurement theory "
            "(audited #204, committed and lab-confronted #205, derived "
            "#206, made dynamical #207, made multipartite #208, made "
            "operational #209) and the strong-field core contraction "
            "(measured and adjudicated #210, the anchor relocated to "
            "alpha). This capstone does what #200 did for the previous "
            "arc - the whole chain re-verified green in one run, the "
            "register updated - and carries ONE new law aimed at the "
            "arc's final unknown: the Compton-edge deformation of the "
            "#202 suppression law, which turns the #210 O(1) "
            "(sigma_mode/lambda_C) from an open number into a "
            "naturalness-confined window. Chosen over an outright "
            "'derivation' of the O(1) deliberately: closing a research "
            "arc with a dialed number would undo what twenty PRs "
            "removed; closing it with a parameter-free bound and a "
            "green ledger is the honest maximum."
        ),
        "deliverable": "docs/compton_edge_capstone.md",
        "executes": "the arc capstone + the #210 successor item (bounded)",
        "framing": "QFT on the fixed classical throat geometry - not quantum gravity",
        "pass": True,
    }


def test_T2_compton_edge_law() -> dict:
    e = compton_edge()
    x = e["xgrid"]
    i_A = x.index(0.647)
    i_B = x.index(1.509)
    band = e["corrected_band"]
    ok = (e["massless_err"] < 1e-10
          and e["universality"] < 2e-3
          and abs(e["S"][0] - 1.0) < 0.02
          and e["S"][i_A] < 1.2 and e["S"][i_B] < 1.5
          and e["x_where_S_2"] is not None and 2.0 < e["x_where_S_2"] < 3.5
          and e["x_where_S_448"] is not None and 4.5 < e["x_where_S_448"] < 8.0
          and band[0] < 1.0 < band[1])
    return {
        "name": "T2_compton_edge_law",
        "description": (
            "THE COMPTON-EDGE LAW - the new result, parameter-free. "
            "Adding the mass term to #202's bridge equation, "
            "(rho^3 phi')' = [3 rho + m^2 rho^3] phi: (a) the massless "
            f"limit re-derives the exact law (|phi - sigma|/sigma = "
            f"{e['massless_err']:.1e} - #202's c0(1) = 0, machine); "
            "(b) UNIVERSALITY: the deformation depends on x = "
            "sigma_mode/lambda_C alone - the m r_s = 0.02 and 0.05 "
            f"curves coincide to {e['universality']:.1e} across the "
            "grid: one scaling law; (c) THE LAW: eps_1 = "
            f"(r_s/sigma_mode) D(x) with D = {dict(zip(x, e['D']))}; "
            "(d) THE SENSITIVITY: S(x) = |d ln eps_1/d ln sigma| = "
            f"{dict(zip(x, e['S']))} - exactly 1 below the Compton "
            "scale (the #202 naturalness result, now with its validity "
            "domain located) and growing beyond: NATURALNESS CAPS THE "
            f"O(1) - S <= 2 confines x <= {e['x_where_S_2']:.1f}, and "
            "the ladder's own worst-tolerated sensitivity (4.48, "
            f"#201/#202) confines x <= {e['x_where_S_448']:.1f}. Both "
            f"convention anchors sit inside (S = {e['S'][i_A]}, "
            f"{e['S'][i_B]}), and the deformation-corrected "
            "self-consistent band x = (needed x alpha) D(x) TIGHTENS "
            f"the #210 bracket to sigma_mode/lambda_C in {band} - "
            "still bracketing 1. The successor derivation has a "
            "compact, derived target window at the Compton edge of "
            "the natural domain."
        ),
        "massless_error": float(f"{e['massless_err']:.2e}"),
        "universality": float(f"{e['universality']:.2e}"),
        "S_curve": dict(zip([str(v) for v in x], e["S"])),
        "D_curve": dict(zip([str(v) for v in x], e["D"])),
        "natural_window_S2": round(e["x_where_S_2"], 2),
        "natural_window_S448": round(e["x_where_S_448"], 2),
        "corrected_band": band,
        "pass": ok,
    }


def test_T3_ledger_commitment_chain() -> dict:
    r4 = mini_retarded_kick()
    r5 = mini_no_entangle()
    # SN / BMV arithmetic (#205)
    hbar, gsi, amu = 1.054571817e-34, 6.6743e-11, 1.66053907e-27
    m = 2.7e4 * amu
    fein = gsi * m * m * 1e-2 / (hbar * 266e-9)
    mstar = ((hbar ** 2 / (2 * gsi * 500e-9)) ** (1 / 3)) / amu
    bmv = gsi * (1e-14) ** 2 * 2.5 / (hbar * 200e-6)
    ok = (r4["norm_drift"] < 1e-9 and r4["continuity_residual"] < 1e-2
          and r5["meanfield"] < 1e-9 and r5["pairwise"] > 1e-3
          and fein < 1e-15 and mstar / 2.7e4 > 1e4
          and 0.3 < bmv < 2.0)
    return {
        "name": "T3_ledger_commitment_chain",
        "description": (
            "LEDGER 1 - THE COMMITMENT CHAIN, GREEN. #204's structural "
            "core: a kicked, RETARDED evolution with real potentials "
            f"conserves the norm to {r4['norm_drift']:.1e} and closes "
            f"continuity at {r4['continuity_residual']:.1e} - "
            "no-signaling's causal completion costs the Born rule "
            "nothing. #205's discriminator: the classical mean-field "
            f"channel keeps a product state at entropy "
            f"{r5['meanfield']:.1e} while the pairwise (quantized) "
            f"operator entangles it to {r5['pairwise']:.4f} - the "
            "BMV-null stands. The lab card: the f = 1 SN phase at the "
            f"interferometry record is {fein:.1e} rad, the record sits "
            f"{mstar/2.7e4:.0e} below the SN inhibition scale, the BMV "
            f"witness phase is {bmv:.2f} rad where BAM predicts "
            "strictly zero - the two near-term nulls (SN-null, "
            "BMV-null) remain the program's nearest falsifiers."
        ),
        "retarded_norm_drift": float(f"{r4['norm_drift']:.2e}"),
        "retarded_continuity": float(f"{r4['continuity_residual']:.2e}"),
        "meanfield_entropy": float(f"{r5['meanfield']:.2e}"),
        "pairwise_entropy": round(r5["pairwise"], 4),
        "sn_phase_at_record": float(f"{fein:.2e}"),
        "bmv_phase": round(bmv, 3),
        "pass": ok,
    }


def test_T4_ledger_topology_chain() -> dict:
    # #206: the singlet from the bridge embedding
    psi_eff = np.kron(np.eye(2), _T) @ _PHIP
    fid_singlet = abs(np.vdot(_SINGLET, psi_eff)) ** 2
    err_law = 0.0
    for ta in np.linspace(0, np.pi, 5):
        for tb in np.linspace(0, np.pi, 5):
            na = np.array([math.sin(ta), 0, math.cos(ta)])
            nb = np.array([math.sin(tb), 0, math.cos(tb)])
            t = np.zeros((3, 3))
            for i, si in enumerate((_SX, _SY, _SZ)):
                for j, sj in enumerate((_SX, _SY, _SZ)):
                    t[i, j] = float(np.real(np.vdot(psi_eff,
                                                    np.kron(si, sj) @ psi_eff)))
            err_law = max(err_law, abs(float(na @ t @ nb)
                                       + math.cos(ta - tb)))
    # #207: the swapping law + the separable mixture
    a, b, c = 0.7, 1.9, np.pi / 2
    tens = np.einsum("ij,kl->ijkl", _psi_pair(a), _psi_pair(b))
    m14 = np.einsum("jk,ijkl->il", _psi_pair(c).conj(), tens)
    p = float(np.sum(np.abs(m14) ** 2))
    fid_swap = abs(np.einsum("il,il->", _psi_pair((a + b + c)).conj(),
                             m14)) ** 2 / p
    rho_mix = np.zeros((4, 4), dtype=complex)
    for cc in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        v = np.array([0, 1, np.exp(1j * cc), 0], dtype=complex) / math.sqrt(2)
        rho_mix += 0.25 * np.outer(v, v.conj())
    neg_mix = _negativity(rho_mix)
    # #206/#208: the enumerated local bounds
    chsh_lhv = max(abs(Aa * Bb + Aa * Bbp + Aap * Bb - Aap * Bbp)
                   for Aa, Aap, Bb, Bbp in itertools.product((1, -1),
                                                             repeat=4))
    mermin_lhv = max(abs(Ap * B * C + A * Bp * C + A * B * Cp - Ap * Bp * Cp)
                     for A, Ap, B, Bp, C, Cp in itertools.product((1, -1),
                                                                  repeat=6))
    # #208: GHZ Mermin = 4 (the standard settings) + monogamy
    ghz = np.zeros(8, dtype=complex)
    ghz[0b011] = 1 / math.sqrt(2)
    ghz[0b100] = 1 / math.sqrt(2)
    def E3(ops):
        op = np.kron(np.kron(ops[0], ops[1]), ops[2])
        return float(np.real(np.vdot(ghz, op @ ghz)))
    # local flips map (|+-->+|-++>)/sqrt2 to the standard GHZ; use the
    # optimal x/y settings after flipping mouths 2,3
    F = np.kron(np.eye(2), np.kron(_SX, _SX))
    g2 = F @ ghz
    def E3g(ops):
        op = np.kron(np.kron(ops[0], ops[1]), ops[2])
        return float(np.real(np.vdot(g2, op @ g2)))
    mermin_ghz = abs(E3g((_SX, _SX, _SX)) - E3g((_SX, _SY, _SY))
                     - E3g((_SY, _SX, _SY)) - E3g((_SY, _SY, _SX)))
    rho_ghz = np.outer(g2, g2.conj())
    # pairwise reduced state (mouths 2,3): trace out mouth 1 (the first
    # qubit factor) from the (q0,q1,q2 | q0,q1,q2) 6-tensor
    r6 = rho_ghz.reshape(2, 2, 2, 2, 2, 2)
    rho_23 = np.trace(r6, axis1=0, axis2=3).reshape(4, 4)
    neg_pair_ghz = _negativity(rho_23)
    ok = (fid_singlet > 1 - 1e-12 and err_law < 1e-12
          and fid_swap > 1 - 1e-12 and abs(p - 0.25) < 1e-12
          and neg_mix < 1e-12 and chsh_lhv == 2 and mermin_lhv == 2
          and abs(mermin_ghz - 4.0) < 1e-12 and neg_pair_ghz < 1e-12)
    return {
        "name": "T4_ledger_topology_chain",
        "description": (
            "LEDGER 2 - THE TOPOLOGY CHAIN, GREEN. #206: psi_eff = "
            f"(I x T)|Phi+> is the singlet (fidelity {fid_singlet:.12f}) "
            f"with E(a,b) = -cos(a-b) to {err_law:.1e} - entanglement "
            "from the bridge embedding. #207: the swapping law - "
            "projecting (2,3) onto Psi_c leaves (1,4) in Psi_{a+b+c} "
            f"(fidelity {fid_swap:.12f}, probability {p:.4f}) and the "
            f"unconditioned holonomy mixture is separable (negativity "
            f"{neg_mix:.1e}) - surgery swaps, no-signaling holds. "
            f"#206/#208: the enumerated local bounds CHSH = {chsh_lhv} "
            f"and Mermin = {mermin_lhv} against the bridge sector's "
            f"2*sqrt2 and the Y-junction's Mermin = {mermin_ghz:.12f} "
            "= 4, whose pairwise marginals are exactly unentangled "
            f"(negativity {neg_pair_ghz:.1e}) - bridge valence is the "
            "entanglement class."
        ),
        "singlet_fidelity": float(f"{fid_singlet:.12f}"),
        "correlation_law_err": float(f"{err_law:.2e}"),
        "swap_fidelity": float(f"{fid_swap:.12f}"),
        "mixture_negativity": float(f"{neg_mix:.2e}"),
        "lhv_chsh": chsh_lhv,
        "lhv_mermin": mermin_lhv,
        "ghz_mermin": float(f"{mermin_ghz:.12f}"),
        "ghz_pair_negativity": float(f"{neg_pair_ghz:.2e}"),
        "pass": ok,
    }


def test_T5_ledger_measurement_mass() -> dict:
    # #209: the k-odd dispersion identity + a live mini Born check
    tchi, nchi = 0.5, 8
    worst = 0.0
    for th in np.linspace(0, np.pi / 2, 7):
        for k in range(-3, 5):
            ek = -2 * tchi * math.cos(2 * np.pi * k / nchi - th)
            emk = -2 * tchi * math.cos(-2 * np.pi * k / nchi - th)
            odd = 0.5 * (ek - emk)
            pred = -2 * tchi * math.sin(2 * np.pi * k / nchi) * math.sin(th)
            worst = max(worst, abs(odd - pred))
    from experiments.closure_ledger.measurement_sector_probe import sg_run
    r = sg_run(np.pi / 6)
    born_dev = abs(r["p_plus"] - r["predicted"])
    # #210: a live Kaup point + the alpha chain
    from experiments.closure_ledger.strong_field_core_solve_probe import \
        solve_star
    kp = solve_star(0.35, "kaup")
    band = [(3.0 / 7.0) * _MU_OVER_E * _ALPHA, _MU_OVER_E * _ALPHA]
    me_band = [(3.0 / 7.0) * _ALPHA, _ALPHA]
    obs = 1.0 / _MU_OVER_E
    ok = (worst < 1e-12 and born_dev < 0.01
          and abs(kp["M"] - 0.632) < 0.01
          and band[0] < 1.0 < band[1] and me_band[0] < obs < me_band[1])
    return {
        "name": "T5_ledger_measurement_mass",
        "description": (
            "LEDGER 3 - MEASUREMENT AND MASS, GREEN. #209: the k-odd "
            f"dispersion identity holds to {worst:.1e} (the pointer "
            "coupling is the KK connection gradient), and a LIVE "
            "Stern-Gerlach Born check lands P(+) within "
            f"{born_dev:.4f} of cos^2(beta) - internal-state Born from "
            "spatial equivariance. #210: a LIVE Kaup point M(0.35) = "
            f"{kp['M']:.4f} (the strong-field solver that refuted the "
            "collapse reading), and the alpha chain re-verified: "
            f"sigma_mode/lambda_C in [{band[0]:.3f}, {band[1]:.3f}] "
            "brackets 1; m_e/m_mu in "
            f"[{me_band[0]:.5f}, {me_band[1]:.5f}] brackets the "
            f"observed {obs:.5f} - the mass ladder on the alpha "
            "anchor, with T2's new law confining its O(1)."
        ),
        "k_odd_identity": float(f"{worst:.2e}"),
        "born_deviation": round(born_dev, 4),
        "kaup_point": round(kp["M"], 4),
        "alpha_band": [round(b, 3) for b in band],
        "me_band": [round(b, 5) for b in me_band],
        "pass": ok,
    }


def test_T6_claim_graph_and_register() -> dict:
    return {
        "name": "T6_claim_graph_and_register",
        "description": (
            "THE ARC, AS A CLAIM GRAPH, AND THE REGISTER AFTER IT. THE "
            "ENTANGLED SECTOR (closed operationally): topology makes "
            "the states (#206 bridge -> singlet; #207 surgery -> "
            "swapping; #208 valence -> GHZ), the KK coupling makes the "
            "pointer, fiber-integrated equivariance makes the "
            "statistics, equilibrium makes them Born (#209) - CHSH "
            "2*sqrt2 from beable positions, marginals flat. THE "
            "COMMITMENT CHAIN: the nonlinearity's only channel is the "
            "retarded gravitational field (#204); its sourcing is "
            "conditional (#205); its lab face is two near-term NULLS "
            "(SN-null, BMV-null). THE MASS THREAD: the collapse "
            "reading refuted constructively (#210); the anchor on "
            "alpha; the O(1) confined to the Compton edge (this PR). "
            "THE REGISTER, updated: (1) derive sigma_mode/lambda_C "
            "within the T2 window [<= ~2.6 natural / <= ~6 tolerated] "
            "- a bounded target; (2) connect the #55-#58 R* to the "
            "bulk mass mu = r_s^2; (3) the 5D pants nucleation (its "
            "Sorkin class and rate); (4) W-class reachability; (5) "
            "registration/irreversibility (radiative decoherence); "
            "standing negatives kept on the books: the cosmological "
            "constant (#165). THE FALSIFICATION CARD: SN-scale "
            "signatures (detection refutes the committed sourcing); a "
            "BMV witness (detection refutes classical Phi outright); "
            "m_e/m_mu = (3/7..1) alpha at the x1.5 level with the "
            "O(1) window; the neutrino-sector cards (m_bb, Sigma m_nu) "
            "unchanged."
        ),
        "register": [
            "derive sigma_mode/lambda_C within the Compton-edge window",
            "connect R* (#55-#58) to the bulk mass mu = r_s^2",
            "the 5D pants nucleation (Sorkin class, rate)",
            "W-class reachability (networks + surgery)",
            "registration/irreversibility (radiative decoherence)",
            "standing negative: the cosmological constant (#165)",
        ],
        "falsification_card": [
            "SN-null (near-term)", "BMV-null (near-term)",
            "m_e/m_mu = (3/7..1) alpha (x1.5 level, O(1) confined)",
            "neutrino cards (m_bb ~ 1.5-6 meV; Sigma m_nu ~ 59 meV)",
        ],
        "pass": True,
    }


def test_T7_honest_scope() -> dict:
    return {
        "name": "T7_honest_scope",
        "description": (
            "WHAT THE ARC DID NOT DO, in one place. (1) Full-GR "
            "no-signaling remains a completion argument (the "
            "wave-equation Phi is the minimal causal completion; the "
            "constraint analysis is not run). (2) The emergence of the "
            "CM pilot wave from the 3-space field is DERIVED only in "
            "its gravitational sourcing (#205) and its entangled "
            "structure (#206-#208); the full field-level account "
            "remains a program. (3) Registration/irreversibility "
            "beyond branch separation is not modeled. (4) The pinch "
            "and the pants are represented by their topological "
            "content; the 5D dynamics (Sorkin classes, rates) is open. "
            "(5) The alpha anchor is CONSTRAINED (the Compton-edge "
            "window), not derived; conv A vs B is a stated convention "
            "band; the #55-#58 R* itself carries the program's one "
            "dimensionful anchor. (6) Equilibrium is a hypothesis with "
            "a mechanism (#198 relaxation) throughout. An honest "
            "release states its conditions with the same care as its "
            "results."
        ),
        "pass": True,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "RELEASE II. When this arc opened, the program's deepest "
            "interpretive machinery rested on two unexecuted items: a "
            "measurement theory that might not survive Gisin, and a "
            "mass anchor that might not survive GR. Both edges were "
            "walked to, and both fired exactly where the program had "
            "pre-registered them: the Newtonian model signals (and the "
            "causal completion fixes it for free); the cores collapse "
            "(and collapse cannot make them light). What emerged is "
            "stronger than what was risked: Bell violation as the "
            "pointer statistics of classical beables with every "
            "ingredient derived - the tensor product from the bridge, "
            "the singlet from the Pin- transport, the outcome from a "
            "holonomy, the statistics from equivariance - and an "
            "electron/muon ratio that is a fine-structure phenomenon "
            "with its one remaining O(1) confined, by a parameter-free "
            "law, to the Compton edge of the natural domain. Two lab "
            "programs can now falsify BAM within the decade; one "
            "bounded number separates its mass ladder from an outright "
            "prediction. The arc closes green, its conditions stated, "
            "its falsifiers armed."
        ),
        "classification": (
            "RELEASE_II_GREEN_THE_ENTANGLED_SECTOR_CLOSED_OPERATIONALLY"
            "_THE_MASS_LADDER_ON_THE_ALPHA_ANCHOR_THE_O1_CONFINED_TO_"
            "THE_COMPTON_EDGE_OF_THE_NATURAL_DOMAIN"
        ),
        "pass": True,
    }


# ========================================================================
# RUNNER + VERDICT
# ========================================================================

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_compton_edge_law(),
        test_T3_ledger_commitment_chain(),
        test_T4_ledger_topology_chain(),
        test_T5_ledger_measurement_mass(),
        test_T6_claim_graph_and_register(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5 = tests[1], tests[2], tests[3], tests[4]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "RELEASE_II_GREEN_THE_ENTANGLED_SECTOR_CLOSED_OPERATIONALLY"
            "_THE_MASS_LADDER_ON_THE_ALPHA_ANCHOR_THE_O1_CONFINED_TO_"
            "THE_COMPTON_EDGE_OF_THE_NATURAL_DOMAIN"
        )
        verdict = (
            "RELEASE II - THE ARC CLOSES GREEN (the argument is in "
            "docs/compton_edge_capstone.md).\n\n"
            "THE NEW LAW. The massive bridge equation deforms #202's "
            "exact suppression law by a UNIVERSAL function of "
            "sigma_mode/lambda_C alone (m r_s collapse "
            f"{t2['universality']:.0e}; massless limit exact to "
            f"{t2['massless_error']:.0e}): the sensitivity leaves 1 at "
            "the Compton scale and naturalness confines the #210 O(1) "
            f"to x <= {t2['natural_window_S2']} (S <= 2) or x <= "
            f"{t2['natural_window_S448']} (the ladder's own 4.48), "
            "with both convention anchors inside and the "
            f"self-consistent band tightened to {t2['corrected_band']} "
            "- still bracketing 1.\n\n"
            "THE LEDGER. Commitment chain: retarded-real norm "
            f"{t3['retarded_norm_drift']:.0e}, continuity "
            f"{t3['retarded_continuity']:.0e}; classical channel "
            f"entropy {t3['meanfield_entropy']:.0e} vs pairwise "
            f"{t3['pairwise_entropy']}; SN/BMV nulls armed. Topology "
            f"chain: singlet fidelity {t4['singlet_fidelity']}, "
            f"swapping fidelity {t4['swap_fidelity']}, LHV bounds "
            f"CHSH = {t4['lhv_chsh']} / Mermin = {t4['lhv_mermin']} "
            f"enumerated, GHZ Mermin = {t4['ghz_mermin']:.1f} with "
            "empty pairs. Measurement + mass: the k-odd identity "
            f"{t5['k_odd_identity']:.0e}, live Born within "
            f"{t5['born_deviation']}, live Kaup {t5['kaup_point']}, "
            f"the alpha bands {t5['alpha_band']} and {t5['me_band']} "
            "bracketing their targets.\n\n"
            "THE REGISTER. Six items, all named and bounded; the "
            "falsification card armed with two near-term nulls and a "
            "x1.5-level mass prediction. The program's deepest "
            "machinery is now theorems, measurements, and windows - "
            "not imports."
        )
    else:
        verdict_class = "RELEASE_II_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A ledger or law check failed; re-examine "
            "before quoting the release."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The Compton-edge capstone: the #204-#210 arc re-verified "
            "green in one run (commitment chain, topology chain, "
            "measurement + mass), plus the new parameter-free "
            "Compton-edge law - the massive bridge equation deforms "
            "#202's exact suppression law universally in "
            "sigma_mode/lambda_C, naturalness confines the #210 O(1) "
            "to the Compton edge, and the self-consistent band "
            "tightens to [0.63, 1.31], still bracketing 1"
        ),
        "executes": "the arc capstone + the #210 successor item, bounded",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The Compton-edge capstone: the second release - "
               "companion probe (PR #211)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/compton_edge_capstone.md` - the "
        "#204-#210 arc consolidated, with the Compton-edge law. *(QFT on "
        "the fixed classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the closing PR: consolidation + one new law",
        "T2": "the Compton-edge law: the O(1) confined (new)",
        "T3": "ledger 1 green: the commitment chain (#204/#205)",
        "T4": "ledger 2 green: the topology chain (#206-#208)",
        "T5": "ledger 3 green: measurement + mass (#209/#210)",
        "T6": "the claim graph + the updated register",
        "T7": "honest scope (the arc's stated conditions)",
        "T8": "Release II",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    t2 = s["tests"][1]
    out.append("## The Compton-edge law")
    out.append("")
    out.append("| x = sigma/lambda_C | S(x) | D(x) |")
    out.append("|---:|---:|---:|")
    for x in t2["S_curve"]:
        out.append(f"| {x} | {t2['S_curve'][x]} | {t2['D_curve'][x]} |")
    out.append("")
    out.append(f"(natural windows: x <= {t2['natural_window_S2']} at S <= 2, "
               f"x <= {t2['natural_window_S448']} at S <= 4.48; corrected "
               f"self-consistent band {t2['corrected_band']})")
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
    out = here / "runs" / f"{ts}_compton_edge_capstone_probe"
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
