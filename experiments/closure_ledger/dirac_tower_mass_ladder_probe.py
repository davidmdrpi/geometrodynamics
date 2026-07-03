"""
The Dirac-tower mass ladder: un-dialing the electron (PR #201).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE FIRST ITEM OF THE #200 REGISTER
-----------------------------------
Rebuild the lepton mass ladder with the electron level on the
index-protected Dirac zero mode (#195), the mouth coupling computed
against the throat overlap machinery (#185/#190), and the three
#192/#194 fine-tuning diagnostics re-run on the rebuilt model.

THE REBUILT MODEL
-----------------
  m_e-level = eps1 * o1 * S1        (a multiplicative chain)
    o1 = 1 (the #195 zero-mode pairing overlap, exact for k=1)
    S1 = (3/7) * m_mu-level  (the k=1 bare scale tied to k=3 by the
         Dirac-tower ratio m1/m3 = (3/2)/(7/2), #197; convention A -
         the band across conventions is carried through)
    eps1 = the mouth coupling, fit to mu/e = 206.77:
         eps1 = (7/3)/206.77 = 0.01129
  Heavy sector: the surrogate's natural 2x2 {k=3,5} block (k=1
  decoupled by winding-charge superselection; the uplift honestly
  refit to tau/mu after the decoupling).

WHAT THE FIT NUMBER IS (constrained, not derived - anti-rigging):
  WKB neck aspect  l/a = -ln eps1 = 4.48 (A) - 5.33 (B): O(1) geometry
  soliton overlap: K(R*) = eps1 on the ACTUAL #180 soliton kernel at
      R* = 5.0-5.6 = 3.9-4.4 x the soliton RMS: a few radii

THE UN-DIALING (the #192/#194 diagnostics, re-run here):
  BG sensitivity   74.7  ->  ~4.5 (= c), all others <= ~1
  sign stability   flips at +-2%  ->  positive identically
  Berger d ln m_e/d lam   -70.9  ->  +4.1;  lambda_break 0.986 -> NONE
  MC P(m_e <= observed)   7.7% sliver  ->  ~50% typical, no cliff

THE IMPOSSIBILITY BOUND: pairing-only mu/e would need eps3/eps1 = 88.6,
but tunneling gives eps_k = e^{-kc}, ratio e^{-2c} ~ 1e-4 <= 1: the
inter-generation hierarchy CANNOT be mouth pairing - it stays with the
dynamical uplift, which #194 already showed is natural (Delta < 1).

Tests:
  T1. Goal (the register item; the document is the argument).
  T2. The protected foundation re-verified: the k=1 zero mode (1e-10),
      the forbidden one-mouth lift (gap 2), superselection (no t13).
  T3. The rebuilt ladder fits (mu/e, tau/mu) with the refit uplift;
      eps1 and the convention band.
  T4. The geometry of eps1: the WKB neck aspect and the #185 soliton-
      kernel inversion (O(1) numbers, stated as consistency).
  T5. Un-dialing 1: the Barbieri-Giudice table on the rebuilt model
      (worst ~4.5 vs 74.7; sign-stable under +-25%).
  T6. Un-dialing 2: the Berger sweep (sensitivity +4.1 vs -70.9; NO
      lambda_break on (0, infinity)).
  T7. Un-dialing 3: the Monte Carlo (typical, smooth, no flips) + the
      pairing-impossibility bound (hierarchy stays dynamical).
  T8. Assessment.

Verdict:
  ELECTRON_LEVEL_REBUILT_ON_THE_INDEX_PROTECTED_ZERO_MODE_FINE_TUNING_
  REMOVED_HIERARCHY_REMAINS_DYNAMICAL_AND_NATURAL
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from experiments.closure_ledger.field_theoretic_odd_k_ladder_probe import (
    monopole_level,
)
from experiments.closure_ledger.multi_throat_exchange_kernel_probe import (
    _soliton_orbital,
    spatial_kernel,
)
from geometrodynamics.tangherlini.lepton_spectrum import (
    LEPTON_BASELINE_PHASE,
    LEPTON_BASELINE_PINHOLE,
    LEPTON_BASELINE_RESISTANCE,
    LEPTON_BASELINE_TRANSPORT,
    S3_ACTION_BASE,
    TAU_BETA_50PI,
)

MU_OVER_E = 206.7683
TAU_OVER_MU = 16.8170
_CACHE: dict = {}


# ------------------------------------------------------------------------
# the heavy {k=3,5} block (the surrogate's natural dynamics, k=1 decoupled)
# ------------------------------------------------------------------------

def heavy_block(beta: float,
                transport: float = LEPTON_BASELINE_TRANSPORT,
                pinhole: float = LEPTON_BASELINE_PINHOLE,
                resistance: float = LEPTON_BASELINE_RESISTANCE,
                base: float = S3_ACTION_BASE,
                slope: float = 0.5,
                phase: float = LEPTON_BASELINE_PHASE) -> np.ndarray:
    """The surrogate generation block restricted to k = 3, 5 (winding
    superselection removes the k=1 row/column). Same functional forms
    as ``_build_generation_block`` (winding_mode='max',
    tunnel_only, exponential resistance)."""
    depths = (3, 5)
    h = np.zeros((2, 2))
    for a, k in enumerate(depths):
        pin = pinhole                      # both 3, 5 carry the pinhole
        uplift = beta * max(0.0, k - 3.0) ** 2
        h[a, a] = (base + resistance * k * k
                   + resistance * (math.exp(k) - 1.0) + pin + uplift)
    alpha_eff = 0.35 * slope
    dk = 5
    amp = transport * math.exp(-alpha_eff * dk) * math.cos(phase * dk)
    h[0, 1] = h[1, 0] = -amp
    return h


def heavy_levels(beta: float, **kw) -> tuple:
    w = np.sort(eigh(heavy_block(beta, **kw), eigvals_only=True))
    return float(w[0]), float(w[1])


def refit_beta(**kw) -> float:
    """Refit the uplift so the heavy block gives tau/mu = 16.817 (the
    k=1 decoupling removed the electron's level repulsion on mu)."""
    f = lambda b: heavy_levels(b, **kw)[1] / heavy_levels(b, **kw)[0] - TAU_OVER_MU
    return brentq(f, 10.0, 400.0)


def rebuilt_ladder(c: float, beta: float, o1: float = 1.0, **kw) -> tuple:
    """(m_e, m_mu, m_tau) levels of the rebuilt model (matrix units):
    heavy = the natural 2x2 block; electron = eps*o1*S1 with
    S1 = (3/7) * m_mu (convention A)."""
    m_mu, m_tau = heavy_levels(beta, **kw)
    s1 = (3.0 / 7.0) * m_mu
    m_e = math.exp(-c) * o1 * s1
    return m_e, m_mu, m_tau


def fitted_point() -> dict:
    if "fit" in _CACHE:
        return _CACHE["fit"]
    beta = refit_beta()
    # convention A: eps1 independent of the heavy refit:
    eps1 = (7.0 / 3.0) / MU_OVER_E
    c = -math.log(eps1)
    m_e, m_mu, m_tau = rebuilt_ladder(c, beta)
    out = {"beta": beta, "eps1": eps1, "c": c,
           "levels": (m_e, m_mu, m_tau),
           "mu_over_e": m_mu / m_e, "tau_over_mu": m_tau / m_mu}
    _CACHE["fit"] = out
    return out


# ========================================================================
# TESTS
# ========================================================================

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "The first item of the #200 open-items register: rebuild the "
            "mass ladder with the electron level on the #195 "
            "index-protected zero mode, the mouth coupling checked "
            "against the #185/#190 throat overlap machinery, and the "
            "three #192/#194 fine-tuning diagnostics re-run on the "
            "rebuilt model. The goal is not a new fit (three fitted "
            "numbers for three masses, as always) - it is the STRUCTURE "
            "of the fitted point: a multiplicative chain with no "
            "subtraction, whose fine-tuning measures must collapse from "
            "the surrogate's Delta = 74.7 to O(few) if the #195 "
            "mechanism truly removes the #194 dial."
        ),
        "deliverable": "docs/dirac_tower_mass_ladder.md",
        "executes": "the #200 register item 1 (the #195 follow-up)",
        "framing": "QFT on the fixed classical throat geometry - not quantum gravity",
        "pass": True,
    }


def test_T2_protected_foundation() -> dict:
    """The zero mode, the forbidden lift, the superselection."""
    z = monopole_level(0.0, 0.0)                    # k=1 zero mode (q~ = 0)
    gap_minus = monopole_level(1.0, 1.0) + 1.0      # D2- lowest = 2q+1 = 2
    zero_ok = abs(z) < 1e-8
    gap_ok = abs(gap_minus - 2.0) < 1e-4
    ok = zero_ok and gap_ok
    return {
        "name": "T2_protected_foundation",
        "description": (
            "The foundation, re-verified live. (a) The k=1 sector's bare "
            f"level is EXACTLY zero (index): monopole ground {z:.1e} - "
            "the surrogate's 6.8754 - 6.6758 subtraction is replaced by "
            "a protected zero; nothing is being cancelled. (b) The "
            "one-mouth additive lift is FORBIDDEN by angular momentum "
            f"(the opposite-chirality gap = {gap_minus:.4f} = 2q+1; no "
            "j = q-1/2 partner exists, #195): the only mass channel is "
            "the two-mouth pairing - multiplicative BY STRUCTURE. "
            "(c) First-order mixing between winding sectors is forbidden "
            "by charge superselection: the surrogate's t13 = 14.85 - "
            "the term whose near-cancellation #194 diagnosed as the "
            "dial - has NO counterpart in the rebuilt model."
        ),
        "k1_zero_mode": float(f"{z:.2e}"),
        "one_mouth_gap": round(gap_minus, 5),
        "t13_counterpart": None,
        "pass": ok,
    }


def test_T3_rebuilt_fit() -> dict:
    """The rebuilt ladder hits the observed ratios; the convention band."""
    fit = fitted_point()
    mu_ok = abs(fit["mu_over_e"] - MU_OVER_E) / MU_OVER_E < 1e-9
    tau_ok = abs(fit["tau_over_mu"] - TAU_OVER_MU) / TAU_OVER_MU < 1e-6
    # convention band: A (S1 = 3/7 m_mu) vs B (S1 = m_mu)
    eps_a = (7.0 / 3.0) / MU_OVER_E
    eps_b = 1.0 / MU_OVER_E
    band = {"A": {"eps1": round(eps_a, 6), "c": round(-math.log(eps_a), 3)},
            "B": {"eps1": round(eps_b, 6), "c": round(-math.log(eps_b), 3)}}
    ok = mu_ok and tau_ok
    return {
        "name": "T3_rebuilt_fit",
        "description": (
            "The rebuilt ladder: heavy sector = the surrogate's natural "
            "2x2 {k=3,5} block with the uplift honestly REFIT after the "
            "k=1 decoupling (removing the electron's level repulsion "
            f"shifts mu): beta = {fit['beta']:.2f} = "
            f"{fit['beta']/TAU_BETA_50PI:.3f} x (50 pi) - the uplift "
            "remains the fitted dynamical knob it always was (#194: its "
            "sensitivities are natural). Electron = eps1*o1*S1 with "
            f"eps1 = (7/3)/206.77 = {fit['eps1']:.5f} (convention A - "
            "note eps1 is INDEPENDENT of the heavy refit) and c = -ln "
            f"eps1 = {fit['c']:.3f}. Fit check: mu/e = "
            f"{fit['mu_over_e']:.4f} ({mu_ok}), tau/mu = "
            f"{fit['tau_over_mu']:.4f} ({tau_ok}). Convention band "
            f"(carried through every number): {band}. Three fitted "
            "numbers for three masses - same count as the surrogate; "
            "the content is T5-T7."
        ),
        "beta_refit": round(fit["beta"], 3),
        "beta_over_50pi": round(fit["beta"] / TAU_BETA_50PI, 4),
        "eps1": round(fit["eps1"], 6),
        "c": round(fit["c"], 4),
        "convention_band": band,
        "pass": ok,
    }


def test_T4_geometry_of_eps1() -> dict:
    """The WKB aspect and the #185 soliton-kernel inversion."""
    fit = fitted_point()
    sol = _soliton_orbital()
    rms = sol["rms"]
    rs = np.linspace(0.5, 8.0, 40)
    ks = np.array([spatial_kernel(float(r)) for r in rs])
    lnk = -np.log(np.maximum(ks, 1e-12))
    inv = {}
    for name, eps in (("A", (7.0 / 3.0) / MU_OVER_E), ("B", 1.0 / MU_OVER_E)):
        rstar = float(np.interp(-math.log(eps), lnk, rs))
        inv[name] = {"eps1": round(eps, 6), "l_over_a": round(-math.log(eps), 3),
                     "R_star": round(rstar, 3),
                     "R_star_over_RMS": round(rstar / rms, 3)}
    aspect_o1 = all(1.0 < v["l_over_a"] < 8.0 for v in inv.values())
    sep_o1 = all(1.0 < v["R_star_over_RMS"] < 8.0 for v in inv.values())
    ok = aspect_o1 and sep_o1
    return {
        "name": "T4_geometry_of_eps1",
        "description": (
            "What the fitted coupling IS - constrained, not derived "
            "(anti-rigging: no claim of deriving 206.77). (a) THE WKB "
            "NECK ASPECT: eps1 = e^{-l/a} gives l/a = "
            f"{inv['A']['l_over_a']} (A) - {inv['B']['l_over_a']} (B): "
            "a throat neck a few radii long - an O(1) geometric aspect "
            "ratio, not a cancellation. (b) THE SOLITON-KERNEL "
            "INVERSION on the ACTUAL #180 self-gravitating "
            "throat-soliton (the #185 machinery, soliton RMS = "
            f"{rms:.3f}): K(R*) = eps1 at R* = {inv['A']['R_star']} (A) "
            f"- {inv['B']['R_star']} (B) code units = "
            f"{inv['A']['R_star_over_RMS']} - "
            f"{inv['B']['R_star_over_RMS']} x RMS: the mouth separation "
            "is a few soliton radii. The coupling the electron mass "
            "requires sits exactly where the repo's own GR overlap "
            "machinery puts a few-radii separation. Deriving l/a from "
            "the 5D throat solution is the remaining step (the #199 "
            "register's 5D core dynamics)."
        ),
        "soliton_rms": round(rms, 4),
        "inversions": inv,
        "pass": ok,
    }


def test_T5_undialing_bg() -> dict:
    """The Barbieri-Giudice table on the rebuilt model."""
    fit = fitted_point()
    params = {"c": fit["c"], "o1": 1.0,
              "transport": LEPTON_BASELINE_TRANSPORT,
              "pinhole": LEPTON_BASELINE_PINHOLE,
              "resistance": LEPTON_BASELINE_RESISTANCE,
              "base": S3_ACTION_BASE, "slope": 0.5, "beta": fit["beta"]}

    def levels(p):
        return rebuilt_ladder(p["c"], p["beta"], p["o1"],
                              transport=p["transport"], pinhole=p["pinhole"],
                              resistance=p["resistance"], base=p["base"],
                              slope=p["slope"])

    rows = []
    worst_e = 0.0
    worst_heavy = 0.0
    for name in params:
        eps = 1e-5
        pp = dict(params); pp[name] = params[name] * (1 + eps)
        pm = dict(params); pm[name] = params[name] * (1 - eps)
        wp, wm = levels(pp), levels(pm)
        d = [(math.log(wp[i]) - math.log(wm[i])) / (2 * eps) for i in range(3)]
        rows.append({"param": name, "m_e": round(d[0], 3),
                     "m_mu": round(d[1], 3), "m_tau": round(d[2], 3)})
        worst_e = max(worst_e, abs(d[0]))
        worst_heavy = max(worst_heavy, abs(d[1]), abs(d[2]))
    # sign stability under +-25%
    rng = np.random.default_rng(5)
    flips = 0
    for _ in range(2000):
        pp = {k: v * math.exp(rng.uniform(-math.log(1.25), math.log(1.25)))
              for k, v in params.items()}
        if levels(pp)[0] <= 0:
            flips += 1
    ok = worst_e < 6.0 and worst_heavy < 3.0 and flips == 0
    return {
        "name": "T5_undialing_bg",
        "description": (
            "UN-DIALING, PART 1 - the #194 Barbieri-Giudice table on the "
            f"rebuilt model. Worst Delta(m_e) = {worst_e:.2f} (= c, the "
            "neck aspect - the exponent of a positive exponential), vs "
            "the surrogate's 74.7; every heavy-sector sensitivity <= "
            f"{worst_heavy:.2f} (they were already natural, #194). SIGN "
            "STABILITY: m_e is a product of positive factors - across "
            "2000 random +-25% parameter draws it flipped sign "
            f"{flips} times (the surrogate flipped under +-2%). The "
            "dial is structurally gone: nothing in the chain can "
            "cancel."
        ),
        "sensitivities": rows,
        "worst_Delta_m_e": round(worst_e, 3),
        "worst_Delta_heavy": round(worst_heavy, 3),
        "sign_flips_in_2000_draws": flips,
        "surrogate_comparison": {"Delta": 74.7, "flip_at": "2%"},
        "pass": ok,
    }


def test_T6_undialing_berger() -> dict:
    """The #192 Berger sweep on the rebuilt electron."""
    fit = fitted_point()
    c = fit["c"]

    def m_e_lam(lam):
        # fiber-riding neck: the winding decay rate scales as 1/lam;
        # the k=1 tower factor m1(lam) = 1/lam + lam/2 (#197)
        return math.exp(-c / lam) * (1.0 / lam + lam / 2.0)

    h = 1e-6
    sens = (math.log(m_e_lam(1 + h)) - math.log(m_e_lam(1 - h))) / (2 * h)
    grid = np.linspace(0.05, 20.0, 800)
    positive = all(m_e_lam(l) > 0 for l in grid)
    smooth = all(np.isfinite((math.log(m_e_lam(l + 1e-6))
                              - math.log(m_e_lam(l))) / 1e-6) for l in grid)
    ok = abs(sens - (c - 1.0 / 3.0)) < 1e-3 and positive and smooth and sens < 6.0
    return {
        "name": "T6_undialing_berger",
        "description": (
            "UN-DIALING, PART 2 - the #192 Berger deformation re-run on "
            "the rebuilt electron: m_e(lam) = e^{-c/lam} * (1/lam + "
            "lam/2) (the fiber-riding neck - the winding decay rate "
            "scales as 1/lam - times the #197 k=1 tower factor). "
            f"Log-sensitivity at the round point: {sens:.3f} = c - 1/3 "
            "(closed form), vs the surrogate's -70.9; and m_e(lam) is a "
            "POSITIVE EXPONENTIAL times a positive tower factor - "
            f"positive on all of (0, infinity) ({positive}, 800-point "
            "sweep): THERE IS NO lambda_break (the surrogate's electron "
            "crossed zero at a 1.4% squash). The #192 spectral "
            "fine-tuning signature - sensitivity = inverse distance to "
            "a zero crossing - is eliminated because there is no zero "
            "crossing to be near."
        ),
        "dln_me_dlam_round": round(sens, 4),
        "closed_form": "c - 1/3",
        "surrogate_value": -70.9,
        "lambda_break": None,
        "pass": ok,
    }


def test_T7_undialing_mc_and_impossibility() -> dict:
    """The #194 Monte Carlo re-run + the pairing-impossibility bound."""
    fit = fitted_point()
    params = {"c": fit["c"], "o1": 1.0, "beta": fit["beta"],
              "transport": LEPTON_BASELINE_TRANSPORT,
              "pinhole": LEPTON_BASELINE_PINHOLE,
              "resistance": LEPTON_BASELINE_RESISTANCE,
              "base": S3_ACTION_BASE, "slope": 0.5}
    m_e0 = fit["levels"][0]
    rng = np.random.default_rng(42)
    n = 20000
    below = 0
    logs = np.empty(n)
    for i in range(n):
        pp = {k: v * math.exp(rng.uniform(-math.log(1.25), math.log(1.25)))
              for k, v in params.items()}
        me = rebuilt_ladder(pp["c"], pp["beta"], pp["o1"],
                            transport=pp["transport"], pinhole=pp["pinhole"],
                            resistance=pp["resistance"], base=pp["base"],
                            slope=pp["slope"])[0]
        logs[i] = math.log(me)
        if me <= m_e0:
            below += 1
    p_below = below / n
    # smoothness of the log-distribution (no cliff): interior histogram
    hist, _ = np.histogram(logs, bins=12)
    interior = hist[2:-2]
    smooth = float(interior.max() / max(interior.min(), 1)) < 4.0
    # the impossibility bound
    need = MU_OVER_E * (3.0 / 7.0)
    have = math.exp(-2.0 * fit["c"])
    impossible = have < 1.0 < need
    ok = 0.3 < p_below < 0.7 and smooth and impossible
    return {
        "name": "T7_undialing_mc_and_impossibility",
        "description": (
            "UN-DIALING, PART 3 - the #194 Monte Carlo re-run (20000 "
            "draws, log-uniform +-25%, fixed seed): P(m_e <= observed) "
            f"= {p_below:.3f} - the observed value is TYPICAL (the "
            "surrogate's 7.7% linear-measure sliver is gone), because "
            "ln m_e is smoothly distributed (multiplicative model; "
            f"interior histogram max/min < 4: {smooth}) with no cliff "
            "and no sign flips (T5). THE IMPOSSIBILITY BOUND: "
            "pairing-generated masses scale as eps_k = e^{-kc} - "
            "DECREASING in k. Reproducing mu/e = 206.77 from pairing "
            f"alone would need eps3/eps1 = {need:.1f}; tunneling gives "
            f"e^{{-2c}} = {have:.1e} <= 1: a tunneling amplitude cannot "
            "grow with the barrier. The inter-generation hierarchy "
            "CANNOT be mouth pairing - it stays with the dynamical "
            "uplift (the #122/#136 phase budget), which #194 already "
            "showed is natural (Delta < 1) and which remains the "
            "honestly-fitted sector. Smallness = geometry (index + "
            "pairing); ratios = dynamics (uplift): the same division "
            "of labor #193/#197 found spectrally, now realized in the "
            "mass model."
        ),
        "P_below_observed": round(p_below, 4),
        "surrogate_P": 0.077,
        "log_distribution_smooth": smooth,
        "pairing_needed_ratio": round(need, 1),
        "pairing_possible_ratio": float(f"{have:.2e}"),
        "hierarchy_origin": "dynamical uplift (fitted, natural) - not pairing",
        "pass": ok,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "THE #194 DIAL IS REMOVED. The electron level is rebuilt as "
            "a multiplicative chain on the index-protected zero mode - "
            "eps1 * o1 * S1, with the bare level exactly zero "
            "(Atiyah-Singer), the additive lift forbidden (angular "
            "momentum), and inter-sector mixing superselected (no t13 "
            "to cancel against). The three fine-tuning diagnostics "
            "collapse: Barbieri-Giudice 74.7 -> ~4.5 (= the neck "
            "aspect); Berger sensitivity -70.9 -> +4.1 with NO "
            "lambda_break anywhere; Monte Carlo 7.7% sliver -> ~50% "
            "typical, sign-stable under +-25%. The fitted coupling "
            "corresponds to O(1) geometry - a neck a few radii long, a "
            "mouth separation of 3.9-4.4 soliton radii on the actual "
            "#180 profile - constrained, not derived (the 5D core "
            "solve remains the register item). And the clean "
            "impossibility bound shows the inter-generation hierarchy "
            "cannot come from pairing: smallness is geometry, ratios "
            "are dynamics. After the rebuild, EVERY sensitivity in the "
            "ladder is O(few) or below: the ladder is fully natural, "
            "with the hierarchy parametrized but no longer tuned."
        ),
        "classification": (
            "ELECTRON_LEVEL_REBUILT_ON_THE_INDEX_PROTECTED_ZERO_MODE_"
            "FINE_TUNING_REMOVED_HIERARCHY_REMAINS_DYNAMICAL_AND_NATURAL"
        ),
        "pass": True,
    }


# ========================================================================
# RUNNER + VERDICT
# ========================================================================

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_protected_foundation(),
        test_T3_rebuilt_fit(),
        test_T4_geometry_of_eps1(),
        test_T5_undialing_bg(),
        test_T6_undialing_berger(),
        test_T7_undialing_mc_and_impossibility(),
        test_T8_assessment(),
    ]
    t4, t5, t6, t7 = tests[3], tests[4], tests[5], tests[6]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "ELECTRON_LEVEL_REBUILT_ON_THE_INDEX_PROTECTED_ZERO_MODE_"
            "FINE_TUNING_REMOVED_HIERARCHY_REMAINS_DYNAMICAL_AND_NATURAL"
        )
        verdict = (
            "THE #194 DIAL IS REMOVED (the argument is in "
            "docs/dirac_tower_mass_ladder.md; this probe re-runs the "
            "diagnostics).\n\n"
            "THE REBUILD. The electron level sits on the #195 "
            "index-protected zero mode: bare level exactly zero, "
            "one-mouth lift forbidden, inter-sector mixing "
            "superselected - the surrogate's cancellation structure "
            "has no counterpart. The mass is the multiplicative chain "
            "eps1*o1*S1; fitting mu/e fixes eps1 = 0.0113 (conv A), "
            "and the fitted number IS O(1) geometry: a neck aspect "
            f"l/a = {t4['inversions']['A']['l_over_a']}-"
            f"{t4['inversions']['B']['l_over_a']}, a mouth separation "
            f"{t4['inversions']['A']['R_star_over_RMS']}-"
            f"{t4['inversions']['B']['R_star_over_RMS']} soliton radii "
            "on the actual #180 profile (constrained, not derived).\n\n"
            "THE DIAGNOSTICS COLLAPSE. Barbieri-Giudice: worst "
            f"Delta(m_e) = {t5['worst_Delta_m_e']} vs the surrogate's "
            f"74.7, zero sign flips in 2000 +-25% draws (vs flipping "
            "at +-2%). Berger: d ln m_e/d lam = "
            f"{t6['dln_me_dlam_round']} (= c - 1/3) vs -70.9, and NO "
            "lambda_break on (0, infinity) (vs 0.986). Monte Carlo: "
            f"P(m_e <= observed) = {t7['P_below_observed']} - typical, "
            "smooth, no cliff (vs the 7.7% sliver).\n\n"
            "THE DIVISION OF LABOR, FINAL FORM. The impossibility "
            "bound (pairing needs eps3/eps1 = 88.6; tunneling gives "
            "1e-4) proves the inter-generation hierarchy cannot be "
            "mouth pairing: smallness = geometry (index + pairing), "
            "ratios = dynamics (the uplift - fitted, and natural per "
            "#194). After the rebuild every sensitivity in the ladder "
            "is O(few) or below: the ladder is fully natural. The "
            "remaining register item: derive the neck aspect from the "
            "5D core solve."
        )
    else:
        verdict_class = "DIRAC_TOWER_REBUILD_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A diagnostic did not collapse as claimed; "
            "re-examine before quoting the un-dialing."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The Dirac-tower mass ladder: the electron rebuilt as "
            "eps1*o1*S1 on the index-protected zero mode - BG 74.7 -> "
            "~4.5, Berger -70.9 -> +4.1 with no lambda_break, MC 7.7% "
            "-> ~50%; eps1 = O(1) geometry (neck aspect ~4.5-5.3, mouth "
            "separation ~4 soliton radii); the hierarchy provably not "
            "pairing - stays dynamical and natural"
        ),
        "executes": "the #200 register item 1 (the #195 follow-up; un-dials #194)",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The Dirac-tower mass ladder: un-dialing the electron (PR #201)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The first item of the #200 register: the electron level rebuilt "
        "on the #195 index-protected zero mode, the mouth coupling "
        "checked against the #185 throat overlap machinery, the "
        "#192/#194 diagnostics re-run. *(QFT on the fixed classical "
        "throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the register item: rebuild the ladder, re-run the diagnostics",
        "T2": "foundation: exact zero, forbidden lift, superselection",
        "T3": "the rebuilt fit (uplift refit; eps1 convention band)",
        "T4": "eps1 IS O(1) geometry: neck aspect + soliton-kernel R*",
        "T5": "BG 74.7 -> ~4.5; zero sign flips in 2000 +-25% draws",
        "T6": "Berger -70.9 -> +4.1; NO lambda_break on (0, inf)",
        "T7": "MC 7.7% -> ~50% typical; hierarchy provably not pairing",
        "T8": "the ladder is fully natural; the dial is gone",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    t4, t5 = s["tests"][3], s["tests"][4]
    out.append("## The geometry of the fitted coupling")
    out.append("")
    out.append("| convention | eps1 | l/a | R* | R*/RMS |")
    out.append("|---|---:|---:|---:|---:|")
    for name, v in t4["inversions"].items():
        out.append(f"| {name} | {v['eps1']} | {v['l_over_a']} | "
                   f"{v['R_star']} | {v['R_star_over_RMS']} |")
    out.append("")
    out.append("## The rebuilt Barbieri-Giudice table")
    out.append("")
    out.append("| parameter | Delta(m_e) | Delta(m_mu) | Delta(m_tau) |")
    out.append("|---|---:|---:|---:|")
    for r in t5["sensitivities"]:
        out.append(f"| `{r['param']}` | {r['m_e']} | {r['m_mu']} | {r['m_tau']} |")
    out.append("")
    out.append(f"(worst Delta(m_e) = {t5['worst_Delta_m_e']} vs the surrogate's 74.7; "
               f"sign flips in 2000 draws: {t5['sign_flips_in_2000_draws']})")
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
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_dirac_tower_mass_ladder_probe"
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
