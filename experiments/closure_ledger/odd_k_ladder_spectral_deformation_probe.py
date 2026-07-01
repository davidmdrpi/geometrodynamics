"""
The spectral deformation test: upgrading #183 from algebra to spectrum
(PR #192).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
> not quantum gravity.

THE GAP THIS PROBE CLOSES
─────────────────────────
PR #183 proved the odd-k {1,3,5} generation sector is protected by
metric-independent ALGEBRA — the deck determinant, the Pin⁻ spin-closure
class ``½ tr T² = −1``, the odd-parity grading.  But algebra is not
spectrum: the protection claim the thesis actually needs is that the
LADDER — three positive, ordered levels with the observed mass ratios —
survives a metric deformation.  #183 cannot distinguish "the topology
protects the spectrum" from "the round metric was doing spectral work the
topology can't protect".

This probe distinguishes them.  It rebuilds the locked lepton
Hamiltonian's geometric ingredients on a Berger-squashed S³_λ (the fiber
squashed by λ, the base S² round — the machinery built by the #165
R-unification audit, whose genuine SU(2) spectrum is imported here) and
tracks the ladder as λ moves off 1.

  PASS: the {1,3,5} structure and mass ratios deform SMOOTHLY over a
        FINITE λ window.
  FAIL: the ladder breaks at INFINITESIMAL squash — the round metric was
        doing spectral work the topology can't protect.

THE INGREDIENT MAP (declared BEFORE the sweep; no post-hoc reassignment)
────────────────────────────────────────────────────────────────────────
The Berger squash multiplies the Hopf-fiber length by λ and leaves the
base S² and the Hopf connection untouched.  Each locked ingredient of the
generation block (``_build_generation_block``, the LEPTON_BASELINE_*
constants) is classified by what it rides on:

  FIBER-RIDING (× λ, a path-length/action along the squashed fiber):
    - action_base = 2π            (the fiber circumference per pass)
    - action_slope = 0.5          (the per-winding tunnel action)
    - resistance_scale = 0.2179   (transport resistance along the winding)
    - k_uplift_beta = 50π         (the τ uplift, counted in fiber 2π-quanta)
  METRIC-BLIND (λ-independent):
    - phase_per_pass = 0.001      (the Hopf HOLONOMY — connection-level,
                                   the very invariant #183 protects)
    - hard_pinhole_gamma = 22.5   (a localized barrier on the round base S²)
    - transport_strength = 25.1   (an attempt-frequency prefactor)

Ambiguous assignments (resistance, uplift) are FLIPPED in the robustness
control T7; the conclusion must survive every map or it is not reported.

RESULT (mixed — and both halves are real)
─────────────────────────────────────────
* The STRUCTURE passes: three positive ordered levels persist over a
  finite window λ ∈ (0.986, ≥3]; the response at λ = 1 is LINEAR (no
  jump at infinitesimal squash) — the #183 fail-mode is ruled out, the
  protection claim upgrades from algebra to spectrum.
* The HIERARCHY does not: the electron level is a NEAR-ZERO (0.1996 in
  action units vs μ 41.26, τ 694.98) that the round metric holds just
  1.4% of squash above a zero crossing (λ_break = 0.98598).  The μ/e
  log-sensitivity −70.9 EQUALS −1/(1−λ_break) = −71.3 — the steepness IS
  the proximity to the spectral boundary, in every ingredient map.  τ/μ
  is gentle (+0.82).  The round metric is doing real spectral work for
  the e–μ hierarchy, which the topology cannot protect.

Tests:
  T1. Goal + framing (algebra → spectrum; the pass/fail criteria).
  T2. The machinery: the genuine SU(2) Berger spectrum (#165 audit,
      imported) + the λ=1 Hamiltonian reproduces the locked masses.
  T3. The ingredient map, declared before the sweep; the #183 algebra
      layer is λ-independent by construction (complementarity).
  T4. The spectral deformation: the ladder deforms smoothly over a
      finite window (positivity, ordering, finite derivatives).
  T5. No infinitesimal breakdown: the response at λ=1 is linear
      (central-difference slope converges) — the FAIL-mode ruled out.
  T6. The finite boundary and the fine-tuning: λ_break = 0.986 (the
      electron level's zero crossing), sensitivity = 1/(1−λ_break).
  T7. Robustness: all four ingredient maps agree (finite boundary,
      sensitivity–boundary identity, gentle τ/μ, intact stretch side).
  T8. Assessment.

Verdict:
  ODD_K_LADDER_SPECTRUM_DEFORMS_SMOOTHLY_OVER_A_FINITE_BERGER_WINDOW
  _BUT_THE_MU_E_HIERARCHY_IS_METRIC_FINE_TUNED
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.linalg import eigh

from experiments.closure_ledger.berger_r_unification_audit_probe import (
    berger_laplacian_level,
    conformal_frequencies,
)
from geometrodynamics.tangherlini.lepton_spectrum import (
    LEPTON_BASELINE_DEPTHS,
    LEPTON_BASELINE_PHASE,
    LEPTON_BASELINE_PINHOLE,
    LEPTON_BASELINE_RESISTANCE,
    LEPTON_BASELINE_TRANSPORT,
    S3_ACTION_BASE,
    TAU_BETA_50PI,
    _build_generation_block,
    solved_lepton_masses_mev,
)

# ── the #183 algebra layer (metric-independent by construction) ──────────
_T_CLOSURE = np.array([[0.0, 1.0], [-1.0, 0.0]])       # T = iσ_y, T² = −I
_HALF_TR_T2 = 0.5 * float(np.trace(_T_CLOSURE @ _T_CLOSURE))
_DECK_DET_BRANE = (-1.0) ** 3                          # RP² deck: −I in 3D frame
_DECK_DET_BULK = (-1.0) ** 4                           # RP³ deck: −I in 4D frame


# ════════════════════════════════════════════════════════════════════════
# THE BERGER-DEFORMED GENERATION BLOCK
# ════════════════════════════════════════════════════════════════════════

def deformed_levels(
    lam: float,
    res_fiber: bool = True,
    up_fiber: bool = True,
    slope_fiber: bool = True,
) -> np.ndarray:
    """The locked generation-block Hamiltonian rebuilt on S³_λ.

    Fiber-riding ingredients (path-length/action along the Hopf fiber)
    scale by λ; metric-blind ingredients (the Hopf holonomy phase, the
    base-S² pinhole, the transport prefactor) do not.  The flip flags
    move the AMBIGUOUS assignments (resistance, tunnel slope, uplift) to
    the metric-blind side for the T7 robustness control.  Returns the
    sorted eigenvalues (action units; masses = electron-calibrated
    positive eigenvalues)."""
    h = _build_generation_block(
        depths=LEPTON_BASELINE_DEPTHS,
        phase_per_pass=LEPTON_BASELINE_PHASE,
        transport_strength=LEPTON_BASELINE_TRANSPORT,
        resistance_model="exponential",
        resistance_scale=LEPTON_BASELINE_RESISTANCE * (lam if res_fiber else 1.0),
        hard_pinhole_gamma=LEPTON_BASELINE_PINHOLE,
        action_base=S3_ACTION_BASE * lam,
        action_slope=0.5 * (lam if slope_fiber else 1.0),
        depth_cost_mode="tunnel_only",
        winding_mode="max",
        k_uplift_beta=TAU_BETA_50PI * (lam if up_fiber else 1.0),
    )
    return np.sort(eigh(h, eigvals_only=True))


def ladder_ratios(lam: float, **kw) -> Optional[tuple]:
    """(μ/e, τ/μ) from the deformed levels, or None if the ladder has
    lost a positive level (the genuine breakdown; the library's
    edge-padding is NOT used here — a lost level is reported, not
    masked)."""
    w = deformed_levels(lam, **kw)
    if not (w[0] > 0.0 and w[1] > w[0] and w[2] > w[1]):
        return None
    return float(w[1] / w[0]), float(w[2] / w[1])


def find_lambda_break(**kw) -> float:
    """Bisect the squash-side boundary where the lowest (electron) level
    crosses zero — the spectral edge of the {1,3,5} window."""
    lo, hi = 0.5, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if deformed_levels(mid, **kw)[0] > 0.0:
            hi = mid
        else:
            lo = mid
    return hi


def log_sensitivities(h: float = 1e-4, **kw) -> tuple:
    """(d ln(μ/e)/dλ, d ln(τ/μ)/dλ) at the round point, by central
    difference."""
    rp, rm = ladder_ratios(1.0 + h, **kw), ladder_ratios(1.0 - h, **kw)
    d_mue = (math.log(rp[0]) - math.log(rm[0])) / (2.0 * h)
    d_tm = (math.log(rp[1]) - math.log(rm[1])) / (2.0 * h)
    return d_mue, d_tm


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "Upgrade #183 from ALGEBRA to SPECTRUM. #183 proved the odd-k "
            "{1,3,5} sector is protected by metric-independent invariants "
            "(deck determinant, ½trT²=−1, odd parity) — but algebra is not "
            "spectrum. This probe rebuilds the locked lepton Hamiltonian's "
            "geometric ingredients on a Berger-squashed S³_λ (fiber × λ, "
            "base round — the #165 machinery) and tracks the ladder as λ "
            "moves off 1. PASS: the {1,3,5} structure and mass ratios "
            "deform smoothly over a finite λ window. FAIL: the ladder "
            "breaks at infinitesimal squash — the round metric was doing "
            "spectral work the topology can't protect. Either outcome is a "
            "real result; #183 cannot distinguish them."
        ),
        "upgrades": "PR #183 (algebraic protection) → spectral protection",
        "machinery": "PR #165 Berger-sphere SU(2) spectrum (imported)",
        "framing": "QFT on the fixed classical throat geometry — not quantum gravity",
        "pass": True,
    }


def test_T2_machinery() -> dict:
    """The genuine SU(2) Berger spectrum + the locked λ=1 baseline."""
    # the imported #165 spectrum: round tower at λ=1
    round_ok = True
    for n in range(6):
        w = conformal_frequencies(n, 0.0)
        if not np.allclose(w, n + 1.0):
            round_ok = False
    # fiber modes stiffen as λ⁻²: the top |m|=j eigenvalue at λ=0.8
    a = 0.8 ** -2 - 1.0
    eig_sq, _ = berger_laplacian_level(4, a)     # j = 2
    eig_rd, _ = berger_laplacian_level(4, 0.0)
    stiffened = float(np.max(eig_sq)) > float(np.max(eig_rd))
    # the λ=1 Hamiltonian IS the locked baseline: masses reproduce
    w1 = deformed_levels(1.0)
    scale = 0.51099895 / w1[0]
    ours = np.array([scale * w1[0], scale * w1[1], scale * w1[2]])
    locked = solved_lepton_masses_mev()
    baseline_ok = bool(np.allclose(ours, locked, rtol=1e-12))
    ok = round_ok and stiffened and baseline_ok
    return {
        "name": "T2_machinery",
        "description": (
            "The machinery is the genuine article, not rebuilt for the "
            "occasion. (a) The SU(2) Berger spectrum is IMPORTED from the "
            "#165 R-unification audit: at λ=1 it collapses to the round "
            f"conformal tower (ω=n+1, verified = {round_ok}); squashing "
            "stiffens the fiber modes (Δ(j,m)=4j(j+1)+4m²(λ⁻²−1): max "
            f"eigenvalue at j=2 rises {np.max(eig_rd):.0f} → "
            f"{np.max(eig_sq):.1f} at λ=0.8). (b) At λ=1 the deformed "
            "Hamiltonian IS the locked baseline: electron-calibrated masses "
            f"(e, μ, τ) = ({ours[0]:.4f}, {ours[1]:.2f}, {ours[2]:.1f}) MeV "
            f"match solved_lepton_masses_mev() to machine precision "
            f"({baseline_ok}) — μ/e = {ours[1]/ours[0]:.2f}, τ/μ = "
            f"{ours[2]/ours[1]:.2f}."
        ),
        "round_tower_recovered": round_ok,
        "fiber_modes_stiffen": stiffened,
        "lambda1_reproduces_locked_masses": baseline_ok,
        "masses_mev_at_round": [round(float(x), 6) for x in ours],
        "pass": ok,
    }


def test_T3_ingredient_map() -> dict:
    """The fiber/base classification, declared before the sweep; the
    #183 algebra layer is λ-independent by construction."""
    # the algebra layer never references λ: recompute the invariants
    algebra_ok = (
        abs(_HALF_TR_T2 - (-1.0)) < 1e-15
        and _DECK_DET_BRANE == -1.0
        and _DECK_DET_BULK == 1.0
    )
    fiber = {
        "action_base": "2π → 2πλ (the fiber circumference per pass)",
        "action_slope": "0.5 → 0.5λ (per-winding tunnel action = fiber path length)",
        "resistance_scale": "0.2179 → 0.2179λ (transport resistance along the winding)",
        "k_uplift_beta": "50π → 50πλ (the τ uplift counted in fiber 2π-quanta)",
    }
    blind = {
        "phase_per_pass": "0.001 (the Hopf HOLONOMY — connection-level; the #183 invariant)",
        "hard_pinhole_gamma": "22.5 (a localized barrier on the round base S²)",
        "transport_strength": "25.1 (attempt-frequency prefactor)",
    }
    return {
        "name": "T3_ingredient_map",
        "description": (
            "The map is declared BEFORE the sweep and never reassigned "
            "post-hoc. The Berger squash multiplies the Hopf-fiber length "
            "by λ and leaves the base S² and the Hopf CONNECTION untouched, "
            "so: fiber-riding (× λ) = the 2π base action (fiber "
            "circumference), the tunnel slope, the winding resistance, the "
            "τ uplift; metric-blind = the Hopf holonomy phase (the very "
            "invariant #183 protects), the base-S² pinhole, the transport "
            "prefactor. The two AMBIGUOUS assignments (resistance, uplift) "
            "are flipped in T7 — the conclusion must survive every map. "
            "COMPLEMENTARITY: the #183 algebra layer (½trT² = "
            f"{_HALF_TR_T2:.0f}, brane deck det {_DECK_DET_BRANE:.0f}, bulk "
            f"{_DECK_DET_BULK:.0f}) never references the metric, so it is "
            "λ-independent BY CONSTRUCTION — the algebra cannot see the "
            "squash; only the spectrum can, which is exactly why this test "
            "is needed."
        ),
        "fiber_riding": fiber,
        "metric_blind": blind,
        "algebra_layer_lambda_independent": algebra_ok,
        "pass": algebra_ok,
    }


def test_T4_spectral_deformation() -> dict:
    """The ladder deforms smoothly over a finite window."""
    grid = [0.99, 0.995, 1.0, 1.01, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0]
    rows = []
    intact = True
    for lam in grid:
        w = deformed_levels(lam)
        r = ladder_ratios(lam)
        if r is None:
            intact = False
            rows.append({"lambda": lam, "levels": [round(float(x), 4) for x in w],
                         "mu_over_e": None, "tau_over_mu": None})
            continue
        rows.append({
            "lambda": lam,
            "levels": [round(float(x), 4) for x in w],
            "mu_over_e": round(r[0], 3),
            "tau_over_mu": round(r[1], 3),
        })
    # smoothness: finite derivatives everywhere on the window (sampled)
    smooth = True
    for lam in [0.99, 1.0, 1.1, 1.5, 2.5]:
        h = 1e-5
        rp, rm = ladder_ratios(lam + h), ladder_ratios(lam - h)
        if rp is None or rm is None:
            smooth = False
            continue
        d = abs(rp[0] - rm[0]) / (2 * h)
        if not np.isfinite(d):
            smooth = False
    # dense ordering/positivity check across the stretch side
    stretch_ok = all(ladder_ratios(l) is not None
                     for l in np.linspace(1.0, 3.0, 81))
    ok = intact and smooth and stretch_ok
    return {
        "name": "T4_spectral_deformation",
        "description": (
            "The {1,3,5} ladder — three positive, ordered, non-degenerate "
            "levels — persists across the whole sampled window λ ∈ [0.99, "
            f"3.0] (dense check λ ∈ [1,3], 81 points: {stretch_ok}), and "
            "the mass ratios deform SMOOTHLY (finite central-difference "
            "derivatives at every sampled λ). The ratios MOVE — μ/e from "
            f"{rows[0]['mu_over_e']} (λ=0.99) through "
            f"{rows[2]['mu_over_e']} (round) to {rows[-1]['mu_over_e']} "
            f"(λ=3), τ/μ from {rows[0]['tau_over_mu']} to "
            f"{rows[-1]['tau_over_mu']} — so the test has dynamic range: "
            "the spectrum genuinely responds to the metric, it is not "
            "frozen by construction."
        ),
        "ladder_curve": rows,
        "structure_intact_on_window": intact,
        "ratios_smooth": smooth,
        "stretch_side_intact_to_3": stretch_ok,
        "pass": ok,
    }


def test_T5_no_infinitesimal_breakdown() -> dict:
    """The response at λ=1 is linear: the FAIL-mode is ruled out."""
    slopes = {}
    for h in [1e-2, 1e-3, 1e-4, 1e-5]:
        rp, rm = ladder_ratios(1.0 + h), ladder_ratios(1.0 - h)
        slopes[f"{h:.0e}"] = round((rp[0] - rm[0]) / (2.0 * h), 1)
    vals = list(slopes.values())
    converged = abs(vals[-1] - vals[-2]) / abs(vals[-1]) < 1e-3
    finite = all(np.isfinite(v) for v in vals)
    ok = converged and finite
    return {
        "name": "T5_no_infinitesimal_breakdown",
        "description": (
            "THE FAIL-MODE TEST. If the round metric were doing spectral "
            "work that breaks at infinitesimal squash, the ratios would "
            "jump or the derivative would diverge as h → 0. Instead the "
            "central-difference slope d(μ/e)/dλ converges to a FINITE "
            f"limit ({vals[-1]:.1f}, stable to <0.1% between h=1e-4 and "
            "h=1e-5): the response is LINEAR in (λ−1) — an ordinary, "
            "differentiable spectral flow, no discontinuity. The ladder "
            "does NOT break at infinitesimal squash; the #183 protection "
            "claim survives its first genuinely spectral test."
        ),
        "d_mu_over_e_dlambda_by_h": slopes,
        "converged_linear_response": converged,
        "pass": ok,
    }


def test_T6_finite_boundary_and_fine_tuning() -> dict:
    """λ_break, and the sensitivity–boundary identity: the μ/e hierarchy
    is metric-fine-tuned."""
    lb = find_lambda_break()
    w1 = deformed_levels(1.0)
    d_mue, d_tm = log_sensitivities()
    inv_dist = 1.0 / (1.0 - lb)
    identity = abs(abs(d_mue) - inv_dist) / inv_dist < 0.05
    finite_not_infinitesimal = (1.0 - lb) > 1e-3
    gentle_tau = abs(d_tm) < 5.0
    ok = finite_not_infinitesimal and identity and gentle_tau
    return {
        "name": "T6_finite_boundary_and_fine_tuning",
        "description": (
            "THE DISCOVERY (the part #183's algebra could not see). The "
            "squash-side window boundary is the ELECTRON level's zero "
            f"crossing at λ_break = {lb:.5f} — a finite 1.4% squash, not "
            "infinitesimal, but CLOSE. The electron eigenvalue at the round "
            f"point is a NEAR-ZERO: {w1[0]:.4f} in action units against μ "
            f"= {w1[1]:.2f} and τ = {w1[2]:.2f} — the celebrated μ/e ≈ 207 "
            "is the ratio of a fine-tuned near-cancellation to a normal "
            "level. The Berger deformation exposes it: the μ/e "
            f"log-sensitivity d ln(μ/e)/dλ = {d_mue:.1f} EQUALS "
            f"−1/(1−λ_break) = −{inv_dist:.1f} (to "
            f"{100*abs(abs(d_mue)-inv_dist)/inv_dist:.1f}%) — the steepness "
            "IS the proximity to the spectral boundary. By contrast τ/μ is "
            f"GENTLE (d ln(τ/μ)/dλ = {d_tm:.2f}, O(1)). So: the topology "
            "protects the {1,3,5} STRUCTURE, and τ/μ is metric-robust — "
            "but the e–μ hierarchy is doing round-metric spectral work the "
            "topology cannot protect. A 1.4% fiber squash makes the "
            "electron level cross zero."
        ),
        "lambda_break": round(lb, 6),
        "one_minus_lambda_break": round(1.0 - lb, 6),
        "electron_level_round": round(float(w1[0]), 5),
        "mu_level_round": round(float(w1[1]), 3),
        "tau_level_round": round(float(w1[2]), 2),
        "dln_mue_dlambda": round(d_mue, 2),
        "inverse_distance_to_boundary": round(inv_dist, 2),
        "sensitivity_equals_inverse_distance": identity,
        "dln_taumu_dlambda": round(d_tm, 3),
        "pass": ok,
    }


def test_T7_robustness() -> dict:
    """All four ingredient maps agree; the conclusion is not an artifact
    of one assignment choice."""
    maps = {
        "default_all_fiber": {},
        "flip_resistance_blind": {"res_fiber": False},
        "flip_uplift_blind": {"up_fiber": False},
        "minimal_only_2pi_base": {
            "res_fiber": False, "up_fiber": False, "slope_fiber": False},
    }
    rows = []
    all_ok = True
    for name, kw in maps.items():
        lb = find_lambda_break(**kw)
        d_mue, d_tm = log_sensitivities(**kw)
        inv = 1.0 / (1.0 - lb)
        ident = abs(abs(d_mue) - inv) / inv < 0.05
        finite = 1e-3 < (1.0 - lb) < 0.2
        gentle = abs(d_tm) < 5.0
        stretch = all(ladder_ratios(l, **kw) is not None
                      for l in np.linspace(1.0, 3.0, 41))
        rows.append({
            "map": name,
            "lambda_break": round(lb, 5),
            "dln_mue_dlambda": round(d_mue, 2),
            "inv_distance": round(inv, 2),
            "identity_holds": ident,
            "dln_taumu_dlambda": round(d_tm, 3),
            "stretch_intact": stretch,
        })
        if not (ident and finite and gentle and stretch):
            all_ok = False
    return {
        "name": "T7_robustness",
        "description": (
            "The two ambiguous assignments (resistance, uplift) are "
            "flipped to the metric-blind side, and a minimal map (ONLY the "
            "unambiguous 2π fiber circumference rides λ) is run. All four "
            "maps agree on every conclusion: a FINITE squash-side boundary "
            "(λ_break ∈ [0.968, 0.986], never infinitesimal), the "
            "sensitivity–boundary identity |d ln(μ/e)/dλ| = 1/(1−λ_break) "
            "(within 5% in each map), a GENTLE τ/μ (|d ln/dλ| ≤ 0.9), and "
            "an intact stretch side to λ=3. The fine-tuned electron "
            "near-zero is a property of the LOCKED HAMILTONIAN, not of the "
            "fiber/base bookkeeping."
        ),
        "maps": rows,
        "all_maps_agree": all_ok,
        "pass": all_ok,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "MIXED — and both halves are real. (1) THE PASS: the {1,3,5} "
            "structure and the mass ratios deform smoothly over a finite λ "
            "window (linear response at λ=1; three ordered positive levels "
            "from λ=0.986 to beyond λ=3) — the ladder does NOT break at "
            "infinitesimal squash, so the #183 protection claim upgrades "
            "from algebra to spectrum. (2) THE CAVEAT: the protection is "
            "for the STRUCTURE, not the numbers. The electron level is a "
            "metric-fine-tuned near-zero — 1.4% of fiber squash away from "
            "a zero crossing, with the μ/e log-sensitivity (−71) exactly "
            "the inverse distance to that boundary — while τ/μ is "
            "metric-robust (+0.8). The round metric IS doing spectral work "
            "for the e–μ hierarchy; the topology guarantees three "
            "generations, not the mass hierarchy. That sharpens, not "
            "weakens, the thesis claim: what is topologically protected "
            "and what is metrically tuned are now separated by "
            "measurement."
        ),
        "classification": (
            "ODD_K_LADDER_SPECTRUM_DEFORMS_SMOOTHLY_OVER_A_FINITE_BERGER_"
            "WINDOW_BUT_THE_MU_E_HIERARCHY_IS_METRIC_FINE_TUNED"
        ),
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════
# RUNNER + VERDICT
# ════════════════════════════════════════════════════════════════════════

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_machinery(),
        test_T3_ingredient_map(),
        test_T4_spectral_deformation(),
        test_T5_no_infinitesimal_breakdown(),
        test_T6_finite_boundary_and_fine_tuning(),
        test_T7_robustness(),
        test_T8_assessment(),
    ]
    t6 = tests[5]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "ODD_K_LADDER_SPECTRUM_DEFORMS_SMOOTHLY_OVER_A_FINITE_BERGER_"
            "WINDOW_BUT_THE_MU_E_HIERARCHY_IS_METRIC_FINE_TUNED"
        )
        verdict = (
            "SPECTRUM, NOT JUST ALGEBRA — with a sharp caveat the algebra "
            "could not see.\n\n"
            "THE PASS. The locked lepton Hamiltonian, rebuilt on the "
            "Berger-squashed S³_λ (fiber-riding ingredients × λ, "
            "connection/base ingredients fixed), keeps its {1,3,5} "
            "structure — three positive, ordered levels — over a FINITE "
            f"window λ ∈ ({t6['lambda_break']:.4f}, ≥3], and the mass "
            "ratios deform SMOOTHLY with a linear response at the round "
            "point. The ladder does not break at infinitesimal squash: the "
            "round metric is not smuggling in the sector structure, and "
            "the #183 protection claim is upgraded from algebraic "
            "invariants to the actual spectrum.\n\n"
            "THE DISCOVERY. The protection does not extend to the "
            "HIERARCHY. The electron level at the round point is a "
            f"near-zero ({t6['electron_level_round']:.4f} in action units "
            f"vs μ {t6['mu_level_round']:.1f}, τ {t6['tau_level_round']:.0f}"
            ") that crosses zero at a 1.4% fiber squash (λ_break = "
            f"{t6['lambda_break']:.5f}); the μ/e log-sensitivity "
            f"{t6['dln_mue_dlambda']:.1f} equals −1/(1−λ_break) = "
            f"−{t6['inverse_distance_to_boundary']:.1f} — the steepness IS "
            "the proximity to the spectral boundary — while τ/μ moves "
            f"gently ({t6['dln_taumu_dlambda']:+.2f}). Robust across all "
            "four ingredient maps (flipped and minimal assignments). So "
            "the topology guarantees THREE GENERATIONS; the round metric "
            "tunes the e–μ hierarchy — the two claims are now separated "
            "by measurement, which is what the thesis needed."
        )
    else:
        verdict_class = "ODD_K_SPECTRAL_DEFORMATION_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A machinery validation, smoothness check, or "
            "robustness map failed; re-examine before reading the window."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The spectral deformation test: the locked {1,3,5} lepton "
            "ladder on the Berger-squashed S³ — structure survives a "
            "finite window (no infinitesimal breakdown), but the electron "
            "level is a metric-fine-tuned near-zero (λ_break = 0.986, "
            "sensitivity = inverse distance to the boundary)"
        ),
        "upgrades": "PR #183 (algebra) → spectrum; machinery from PR #165",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The spectral deformation test: the {1,3,5} ladder on the Berger-squashed S³ (PR #192)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "Upgrades #183 from algebra to spectrum: the locked lepton "
        "Hamiltonian's geometric ingredients are rebuilt on the Berger "
        "sphere S³_λ (the #165 machinery) and the ladder is tracked as λ "
        "moves off 1. *(QFT on the fixed classical throat geometry, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "goal: algebra → spectrum; the pass/fail criteria",
        "T2": "machinery: #165 SU(2) Berger spectrum + locked λ=1 baseline",
        "T3": "the ingredient map (declared before the sweep)",
        "T4": "the ladder deforms smoothly over a finite window",
        "T5": "NO infinitesimal breakdown (linear response at λ=1)",
        "T6": "λ_break = 0.986; μ/e sensitivity = 1/(1−λ_break)",
        "T7": "robust across all four ingredient maps",
        "T8": "structure protected; e–μ hierarchy metric-fine-tuned",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'—')} | {p} |")
    out.append("")
    t4, t6, t7 = s["tests"][3], s["tests"][5], s["tests"][6]
    out.append("## The deformed ladder (action units; ratios electron-calibrated)")
    out.append("")
    out.append("| λ | e level | μ level | τ level | μ/e | τ/μ |")
    out.append("|---:|---:|---:|---:|---:|---:|")
    for r in t4["ladder_curve"]:
        lv = r["levels"]
        out.append(f"| {r['lambda']} | {lv[0]} | {lv[1]} | {lv[2]} | "
                   f"{r['mu_over_e']} | {r['tau_over_mu']} |")
    out.append("")
    out.append("## The boundary and the fine-tuning")
    out.append("")
    out.append(f"- λ_break (electron level → 0): **{t6['lambda_break']}** "
               f"(a {100*t6['one_minus_lambda_break']:.1f}% squash — finite, not infinitesimal)")
    out.append(f"- electron level at round point: **{t6['electron_level_round']}** "
               f"(vs μ {t6['mu_level_round']}, τ {t6['tau_level_round']}) — a near-zero")
    out.append(f"- d ln(μ/e)/dλ = **{t6['dln_mue_dlambda']}** vs −1/(1−λ_break) = "
               f"−{t6['inverse_distance_to_boundary']} — the sensitivity IS the distance to the boundary")
    out.append(f"- d ln(τ/μ)/dλ = **{t6['dln_taumu_dlambda']}** — gentle; τ/μ is metric-robust")
    out.append("")
    out.append("## Robustness (all ingredient maps)")
    out.append("")
    out.append("| map | λ_break | d ln(μ/e)/dλ | 1/(1−λ_break) | identity | d ln(τ/μ)/dλ |")
    out.append("|---|---:|---:|---:|---|---:|")
    for r in t7["maps"]:
        out.append(f"| {r['map']} | {r['lambda_break']} | {r['dln_mue_dlambda']} | "
                   f"{r['inv_distance']} | {r['identity_holds']} | {r['dln_taumu_dlambda']} |")
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
    out = here / "runs" / f"{ts}_odd_k_ladder_spectral_deformation_probe"
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
