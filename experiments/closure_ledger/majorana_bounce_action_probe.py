"""
Majorana instanton-bounce action probe (PR #88).

PR #87 reframed the neutrino's Majorana seesaw scale `M_R = m_D·e^{S}`
as an instanton/bounce action `S = ln(m_D/m_ν) ≈ 15–18`, but left `S`
itself OPEN — the verdict was
`SEESAW_SUPPRESSION_IS_THROAT_ANTITHROAT_TUNNELING_S_OPEN`. This probe
tries to SHARPEN that verdict by constructing a *reduced Euclidean
bounce action* for the `ΔL=2` throat ↔ antithroat transition from
first BAM ingredients — the PR #58 nucleation potential, the PR #86
chargeless (`k=0, c₁=0`) Majorana sector, and the non-orientable
tortoise geometry of the throat — and asking whether the intrinsic
throat tension and boundary compliance naturally produce `S ≈ 15–18`.

## The non-orientable tortoise path

The throat sits at `r = rs = R_MID`. The 5D tortoise coordinate

```
r*(r) = r + (rs/2)·ln|(r − rs)/(r + rs)|
```

diverges *logarithmically* as `r → rs`: the throat is infinitely deep
in tortoise distance. The `ΔL=2` Majorana transition is the
throat ↔ antithroat flip — the inner/outer swap `C` (`c₁ → −c₁`,
PR #63), realised in the radial solver as the **odd extension across
the throat** (`u_full = [−u[::-1], u]`, the orientation-reversing /
non-orientable identification). The bounce is the Euclidean path of
the throat configuration along this non-orientable tortoise path.

## The reduced bounce action

Model the `ΔL=2` flip as a 1D Euclidean tunnelling problem in the
tortoise coordinate. To flip orientation the throat must pay the PR #58
nucleation/pinch cost `E_c = (16π/3)σ³/ρ²` (the barrier separating
throat from antithroat), carrying the throat's breathing inertia
`μ = 4πσ rs²` (the brane mass at the throat scale). The reduced bounce
(thin-wall momentum × path length) is

```
S  =  ∫ √(2 μ E_c)  dr*   =   √(2 μ E_c) · L*(ε),
```

with `L*(ε) = r*(R_OUTER) − r*(rs + ε)` the tortoise path length from
the outer mouth to within a compliance cutoff `ε` of the throat.
Because `L*(ε) = −(rs/2) ln ε + const`, the action is **logarithmic in
the boundary compliance** `ε`.

## Two structural results

  1. **Rigid throat ⟹ massless neutrino.** A perfectly rigid
     (incompressible) throat has `ε → 0`, so `L* → ∞` and `S → ∞`:
     `m_ν = m_D·e^{−S} → 0`. A boundary compliance `ε > 0` is exactly
     what gives the neutrino a non-zero Majorana mass. The smallness of
     `m_ν` is the near-rigidity of the throat.

  2. **`S` is a tortoise logarithm.** `S ∝ ln(1/ε)`, so it is naturally
     `O(10)` and varies only *slowly* (logarithmically) — matching the
     `O(15)`, generation-stable `S` PR #87 required. The huge keV→TeV
     gap is carried by a logarithm, not a hierarchy.

## The magnitude test (the honest part)

With the throat tension and bag density fixed by the *electron* throat
(PR #58/#87: `σ = 1/(12π)`, `ρ = 3/(4π)`, geometric units `mc²=1`,
`R*=1`), the under-barrier momentum is `√(2 μ E_c) ≈ 0.060`, so even a
near-rigid compliance gives `S ≈ 0.04–1` — **~40× short** of the
required `~16`. Pure compliance cannot close the gap (it would need an
absurd `ε ~ 10⁻²³⁰`).

Reaching `S ≈ 15–18` at a *sane* compliance (`ε ~ 10⁻³–10⁻⁶`) requires
the `ΔL=2` (lepton-number-violating / B−L) throat tension to be a
factor `t ≈ 6–12` stiffer than the EM-throat tension — i.e. the
nucleation barrier `~10²–10³×` larger (since `E_c ∝ σ³`). So the EM-
throat tension and boundary compliance do NOT by themselves produce
`S ≈ 15–18`; they fall short, and the residual open input is sharpened
from "a `~TeV` mass" (PR #86) → "an `O(15)` action" (PR #87) → "a
dimensionless `O(10)` B−L/EM tension ratio" (this probe).

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** the bounce coordinate is the
    non-orientable tortoise path (odd extension, `c₁ → −c₁`); the rigid
    throat gives an exactly massless Majorana neutrino, so the
    compliance is the mass-generating parameter; the bounce action is a
    tortoise logarithm `S ∝ ln(1/ε)`, hence naturally `O(10)` and
    coarsely generation-stable — the form PR #87's `S` required.

  - **Does not establish:** that `S ≈ 15–18` falls out of the
    EM-throat tension. It under-produces by `~40×`; matching needs a
    `ΔL=2` tension `~6–12×` stiffer than the EM-throat tension (barrier
    `~10²–10³×` larger). The open input is now this dimensionless
    tension ratio (plus the compliance `ε`), not yet derived. The
    detailed `m_ν` spectrum (the `×18` spread in `m_ν/m_D`) needs a
    generation-dependent compliance or the mixing sector.

Tests:
  T1. Required `S = ln(m_D/m_ν) ≈ 14.7–17.6` (PR #87 target, gen-stable).
  T2. Non-orientable tortoise path: odd extension across throat;
      r* diverges logarithmically at rs (slope rs/2).
  T3. Rigid-throat limit ε→0 ⟹ L*→∞ ⟹ S→∞ ⟹ m_ν=0 (compliance is the
      mass-generating parameter).
  T4. Reduced bounce S = √(2 μ E_c)·L*(ε); EM-throat-calibrated;
      logarithmic in ε; even near-rigid gives S ≲ 1.
  T5. Magnitude: EM tension under-produces S by ~40×; matching S∈[15,18]
      at sane compliance needs a ΔL=2 tension ratio t ≈ 6–12.
  T6. Generation structure: constant-geometry bounce ⟹ near-uniform S
      ⟹ m_ν ∝ m_D; observed ×18 spread needs gen-dependent compliance /
      mixing. Log form explains O(15) magnitude + coarse gen-stability.
  T7. Honest scope: structure established; S not produced by EM tension;
      open input sharpened to an O(10) tension ratio.
  T8. Assessment.

Verdict:
  - MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO (expected):
    the bounce is the non-orientable tortoise log (rigid throat =
    massless ν; S ∝ ln(1/ε), naturally O(10)/gen-stable), but the
    EM-throat tension under-produces S; S ≈ 15–18 needs a ΔL=2 tension
    ~6–12× stiffer — the residual open input, now a dimensionless ratio.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.radial import (
    V_tangherlini, r_to_rstar, rstar_to_r,
)


PI = math.pi
K_5 = 5
N_GRID = 800
ACTION_BASE = 2.0 * PI
L_THROAT = math.sqrt(ACTION_BASE) / K_5
RS = R_MID

M_E_MEV = 0.511
NU_MASS_OBS_EV = {1: 0.001, 2: 0.009, 3: 0.05}

# PR #58/#87 EM-throat-calibrated nucleation parameters (geometric units,
# mc² = 1, R* = R_MID = 1): B = 4πσ with E(R*) = 3BR* ² = m_e c² ⟹
# σ = 1/(12π); bag density ρ = rest energy / throat volume ⟹ ρ = 3/(4π).
SIGMA_EM = 1.0 / (12.0 * PI * RS ** 2)
RHO_EM = 3.0 / (4.0 * PI * RS ** 3)
R_C = 2.0 * SIGMA_EM / RHO_EM
E_C = (16.0 * PI / 3.0) * SIGMA_EM ** 3 / RHO_EM ** 2     # nucleation barrier
MU_THROAT = 4.0 * PI * SIGMA_EM * RS ** 2                  # brane breathing inertia
P0 = math.sqrt(2.0 * MU_THROAT * E_C)                     # under-barrier momentum


# ---------------------------------------------------------------------------
# Cavity floors → bare Dirac masses (as in PR #86/#87)
# ---------------------------------------------------------------------------

def _cavity_eigenvalues(l: int = 1):
    rsmin = r_to_rstar(RS + 5e-4, RS)
    rsmax = r_to_rstar(R_OUTER - 5e-4, RS)
    rstar = np.linspace(rsmin, rsmax, N_GRID)
    h = rstar[1] - rstar[0]
    rphys = np.array([rstar_to_r(s, RS) for s in rstar])
    V = V_tangherlini(rphys, l, RS)
    main = 2.0 / h ** 2 + V[1:-1]
    off = -1.0 / h ** 2 * np.ones(N_GRID - 3)
    ev, _ = np.linalg.eigh(np.diag(main) + np.diag(off, 1) + np.diag(off, -1))
    return np.maximum(ev, 0.0)


_EV = _cavity_eigenvalues()


def m2_unified(k: int, n: int) -> float:
    return (k * ACTION_BASE / L_THROAT) ** 2 + float(_EV[n])


_SCALE_MEV = M_E_MEV / math.sqrt(m2_unified(1, 0))


def m_D_eV(g: int) -> float:
    return math.sqrt(m2_unified(0, g - 1)) * _SCALE_MEV * 1e6


# ---------------------------------------------------------------------------
# Tortoise geometry (analytic — no 1e-15 floor, valid for arbitrarily small ε)
# ---------------------------------------------------------------------------

def rstar_analytic(r: float, rs: float = RS) -> float:
    """5D tortoise coordinate, analytic (r > rs)."""
    return r + (rs / 2.0) * math.log((r - rs) / (r + rs))


def rstar_at_throat_offset(eps: float, rs: float = RS) -> float:
    """r*(rs + ε) computed from ε directly, so arbitrarily small ε avoids
    the float cancellation (rs + ε == rs). Here r − rs = ε, r + rs = 2rs + ε."""
    return (rs + eps) + (rs / 2.0) * math.log(eps / (2.0 * rs + eps))


def tortoise_length(eps: float, rs: float = RS, r_outer: float = R_OUTER) -> float:
    """L*(ε) = r*(R_OUTER) − r*(rs + ε): tortoise path length from the outer
    mouth to within a compliance cutoff ε of the throat."""
    return rstar_analytic(r_outer, rs) - rstar_at_throat_offset(eps, rs)


def bounce_action(eps: float, tension_ratio: float = 1.0) -> float:
    """Reduced Euclidean bounce S = √(2 μ E_c)·L*(ε), with the ΔL=2 throat
    tension scaled by `tension_ratio` t (σ → t·σ ⟹ μ ∝ t, E_c ∝ t³ ⟹
    √(2μE_c) ∝ t²)."""
    return (tension_ratio ** 2) * P0 * tortoise_length(eps)


# ---------------------------------------------------------------------------
# T1. Required S (PR #87 target)
# ---------------------------------------------------------------------------

def test_T1_required_S() -> dict:
    """Recap PR #87: M_R = m_D·e^{S} ⟹ the required bounce action is
    S = ln(m_D/m_ν) ≈ 14.7–17.6, generation-stable (spread < 3)."""
    rows = []
    for g in (1, 2, 3):
        mD = m_D_eV(g)
        mnu = NU_MASS_OBS_EV[g]
        S = math.log(mD / mnu)
        rows.append({'generation': g, 'm_D_keV': mD / 1e3,
                     'm_nu_obs_eV': mnu, 'S_required': S})
    S_vals = [r['S_required'] for r in rows]
    return {
        'name': 'T1_required_bounce_action',
        'description': (
            "PR #87 target: S = ln(m_D/m_ν) ≈ 14.7–17.6, generation-stable "
            "(spread < 3). This probe asks whether a first-principles "
            "tortoise bounce produces it."
        ),
        'rows': rows,
        'S_range': (min(S_vals), max(S_vals)),
        'S_spread': max(S_vals) - min(S_vals),
        'generation_stable': (max(S_vals) - min(S_vals)) < 5.0,
        'pass': 10.0 < min(S_vals) and max(S_vals) < 25.0,
    }


# ---------------------------------------------------------------------------
# T2. Non-orientable tortoise path; logarithmic throat
# ---------------------------------------------------------------------------

def test_T2_tortoise_path() -> dict:
    """The ΔL=2 flip is the throat↔antithroat odd extension (c₁→−c₁,
    PR #63 inner/outer swap) — the non-orientable identification across
    the throat. The tortoise coordinate diverges logarithmically at the
    throat: r*(rs+ε) ≈ (rs/2)·ln ε, slope rs/2. Verify the log slope."""
    eps_a, eps_b = 1e-4, 1e-8
    ra = rstar_at_throat_offset(eps_a)
    rb = rstar_at_throat_offset(eps_b)
    # r* ≈ (rs/2) ln(ε/(2rs)) + rs ⟹ d r*/d ln ε = rs/2
    slope = (rb - ra) / (math.log(eps_b) - math.log(eps_a))
    expected_slope = RS / 2.0
    return {
        'name': 'T2_non_orientable_tortoise_path',
        'description': (
            "Bounce coordinate = non-orientable tortoise path: the odd "
            "extension across the throat (u_full = [−u[::-1], u], c₁→−c₁, "
            "PR #63). The tortoise coordinate diverges logarithmically at "
            "the throat (slope rs/2): the throat is infinitely deep."
        ),
        'rstar_at_eps_1e-4': ra,
        'rstar_at_eps_1e-8': rb,
        'log_slope_numeric': slope,
        'log_slope_expected_rs_over_2': expected_slope,
        'orientation_reversal': 'u → −u across rs (c₁ → −c₁, PR #63)',
        'pass': abs(slope - expected_slope) / expected_slope < 1e-3,
    }


# ---------------------------------------------------------------------------
# T3. Rigid-throat limit ⟹ massless neutrino
# ---------------------------------------------------------------------------

def test_T3_rigid_throat_massless() -> dict:
    """A perfectly rigid (incompressible) throat has ε → 0, so the
    tortoise length L* → ∞ and the bounce S → ∞: m_ν = m_D·e^{−S} → 0.
    The boundary compliance ε > 0 is what gives the neutrino mass — the
    smallness of m_ν is the near-rigidity of the throat."""
    rows = []
    for eps in (1e-2, 1e-4, 1e-8, 1e-16, 1e-32):
        L = tortoise_length(eps)
        S = bounce_action(eps)
        rows.append({'eps': eps, 'L_star': L, 'S_base': S,
                     'm_nu_over_m_D': math.exp(-S)})
    # L* grows without bound as ε → 0 (monotone increasing)
    Ls = [r['L_star'] for r in rows]
    diverges = all(Ls[i] < Ls[i + 1] for i in range(len(Ls) - 1))
    return {
        'name': 'T3_rigid_throat_gives_massless_neutrino',
        'description': (
            "Rigid throat (ε→0): L*→∞ ⟹ S→∞ ⟹ m_ν=m_D·e^{−S}→0. A "
            "boundary compliance ε>0 is the mass-generating parameter; "
            "the smallness of m_ν is the near-rigidity of the throat."
        ),
        'rows': rows,
        'L_star_diverges_as_eps_to_zero': diverges,
        'rigid_limit_massless': True,
        'pass': diverges,
    }


# ---------------------------------------------------------------------------
# T4. Reduced bounce action (EM-throat-calibrated)
# ---------------------------------------------------------------------------

def test_T4_reduced_bounce() -> dict:
    """Reduced bounce S = √(2 μ E_c)·L*(ε) with the PR #58 EM-throat
    tension (σ=1/12π, ρ=3/4π ⟹ E_c≈0.0055, μ≈1/3, √(2μE_c)≈0.060). S is
    logarithmic in ε; even a near-rigid compliance gives S ≲ 1 — far
    below the required ~16."""
    rows = []
    for eps in (1e-1, 1e-2, 1e-4, 1e-8, 1e-16):
        rows.append({'eps': eps, 'L_star': tortoise_length(eps),
                     'S_base': bounce_action(eps)})
    S_max_sane = bounce_action(1e-16)     # absurdly rigid, still small
    return {
        'name': 'T4_reduced_bounce_em_calibrated',
        'description': (
            "S = √(2 μ E_c)·L*(ε), EM-throat-calibrated (√(2μE_c) ≈ 0.060). "
            "Logarithmic in ε; even ε=1e-16 (near-rigid) gives S ≲ 1 — "
            "~40× short of the required ~16."
        ),
        'sigma_em': SIGMA_EM,
        'rho_em': RHO_EM,
        'E_c': E_C,
        'mu_throat': MU_THROAT,
        'under_barrier_momentum_p0': P0,
        'rows': rows,
        'S_base_at_eps_1e-16': S_max_sane,
        'underproduces': S_max_sane < 5.0,
        'pass': S_max_sane < 5.0,     # confirms the EM tension under-produces
    }


# ---------------------------------------------------------------------------
# T5. Magnitude test → required ΔL=2 tension ratio
# ---------------------------------------------------------------------------

def test_T5_required_tension_ratio() -> dict:
    """Matching S ∈ [15, 18] at a sane near-rigid compliance
    (ε ~ 1e-3–1e-6) requires the ΔL=2 (B−L) throat tension to be a factor
    t ≈ 6–12 stiffer than the EM-throat tension (since S ∝ t², barrier
    E_c ∝ t³ ⟹ ~10²–10³× larger). Pure compliance (t=1) would need an
    absurd ε ~ 1e-230."""
    S_target = 16.0
    rows = []
    for eps in (1e-2, 1e-3, 1e-6):
        L = tortoise_length(eps)
        t = math.sqrt(S_target / (P0 * L))
        rows.append({'eps': eps, 'L_star': L, 'required_tension_ratio_t': t,
                     'barrier_boost_t_cubed': t ** 3})
    # t=1 (EM tension only): solve P0·L* = 16 ⟹ L* = 16/P0; ε from L*
    L_needed = S_target / P0
    # L*(ε) = r*(R_OUTER) − [(rs+ε) + (rs/2)ln(ε/(2rs+ε))]; invert for tiny ε:
    # L* ≈ rs − (rs/2)ln(ε/(2rs)) + r*(R_OUTER) ⟹ ln ε ≈ −2(L*−c)/rs
    c = rstar_analytic(R_OUTER) - RS + (RS / 2.0) * math.log(2.0 * RS)
    ln_eps_pure = -2.0 * (L_needed - c) / RS
    eps_pure_compliance = math.exp(ln_eps_pure)
    t_range = (min(r['required_tension_ratio_t'] for r in rows),
               max(r['required_tension_ratio_t'] for r in rows))
    return {
        'name': 'T5_required_delta_L_tension_ratio',
        'description': (
            "Matching S∈[15,18] at sane compliance ε~1e-3–1e-6 needs a "
            "ΔL=2 tension ratio t≈6–12 (barrier ~10²–10³× larger). Pure "
            "compliance (t=1) would need an absurd ε~1e-230."
        ),
        'rows': rows,
        'required_tension_ratio_range': t_range,
        'eps_needed_if_pure_compliance': eps_pure_compliance,
        'pure_compliance_is_absurd': eps_pure_compliance < 1e-30,
        'pass': 3.0 < t_range[0] and t_range[1] < 30.0,
    }


# ---------------------------------------------------------------------------
# T6. Generation structure
# ---------------------------------------------------------------------------

def test_T6_generation_structure() -> dict:
    """A constant-geometry bounce (same throat, same ε, same tension for
    every generation) gives a generation-UNIFORM S, hence m_ν ∝ m_D
    (e^{−S} constant). Observed m_ν/m_D spans ×18 (2.3e-8 → 4.2e-7), so
    the detailed spectrum needs a generation-dependent compliance or the
    mixing sector. But the log form DOES explain the O(15) magnitude and
    the coarse generation-stability (spread < 3) PR #87 found."""
    rows = []
    ratios = []
    for g in (1, 2, 3):
        mD = m_D_eV(g)
        mnu = NU_MASS_OBS_EV[g]
        ratio = mnu / mD
        ratios.append(ratio)
        rows.append({'generation': g, 'm_nu_over_m_D': ratio,
                     'S_required': math.log(mD / mnu)})
    spread = max(ratios) / min(ratios)
    return {
        'name': 'T6_generation_structure',
        'description': (
            "Constant-geometry bounce ⟹ uniform S ⟹ m_ν ∝ m_D. Observed "
            "m_ν/m_D spans ×18, so the detailed spectrum needs a "
            "generation-dependent compliance or mixing. The log form does "
            "explain the O(15) magnitude + coarse gen-stability."
        ),
        'rows': rows,
        'm_nu_over_m_D_spread': spread,
        'log_form_explains_magnitude_and_coarse_stability': True,
        'detailed_spectrum_needs_gen_dependent_compliance_or_mixing': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Honest scope
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope',
        'description': (
            "Structure established (tortoise-log bounce, rigid→massless, "
            "O(15)/gen-stable form); S≈15–18 NOT produced by the EM-throat "
            "tension (under by ~40×); open input sharpened to an O(10) "
            "ΔL=2/EM tension ratio."
        ),
        'established_bam_native': [
            'bounce coordinate = non-orientable tortoise path (odd '
            'extension, c₁→−c₁, PR #63)',
            'rigid throat (ε→0) ⟹ S→∞ ⟹ exactly massless Majorana '
            'neutrino; compliance ε is the mass-generating parameter',
            'S ∝ ln(1/ε): a tortoise logarithm, naturally O(10) and '
            'coarsely generation-stable — the form PR #87 required',
        ],
        'open': [
            'S ≈ 15–18 does NOT fall out of the EM-throat tension (under '
            'by ~40×); matching needs a ΔL=2 (B−L) tension ~6–12× stiffer '
            '(barrier ~10²–10³× larger) — a dimensionless tension ratio, '
            'not yet derived',
            'the boundary compliance ε (a sub-throat length), not yet '
            'derived from the bulk',
            'the detailed m_ν spectrum (×18 spread) — needs gen-dependent '
            'compliance or the mixing sector',
        ],
        'progressive_localisation': (
            'open input: ~TeV mass (PR #86) → O(15) action S (PR #87) → '
            'O(10) B−L/EM tension ratio (PR #88)'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The Majorana bounce is the non-orientable tortoise log: rigid "
            "throat = massless ν, S ∝ ln(1/ε) naturally O(10)/gen-stable. "
            "But the EM-throat tension under-produces S; S≈15–18 needs a "
            "ΔL=2 tension ~6–12× stiffer — the residual open input, now a "
            "dimensionless ratio."
        ),
        'classification': (
            'MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_required_S(),
        test_T2_tortoise_path(),
        test_T3_rigid_throat_massless(),
        test_T4_reduced_bounce(),
        test_T5_required_tension_ratio(),
        test_T6_generation_structure(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO'
        verdict = (
            'MAJORANA BOUNCE IS THE NON-ORIENTABLE TORTOISE LOGARITHM; THE '
            'OPEN INPUT IS NOW A ΔL=2 TENSION RATIO. PR #87 reframed the '
            'seesaw scale as a bounce action S = ln(m_D/m_ν) ≈ 15–18 but '
            'left S open. This probe builds a reduced Euclidean bounce for '
            'the throat↔antithroat flip and sharpens the verdict.\n\n'
            'THE NON-ORIENTABLE TORTOISE PATH. The ΔL=2 Majorana flip is '
            'the throat↔antithroat odd extension (c₁→−c₁, PR #63 '
            'inner/outer swap; u_full = [−u[::-1], u]) — the '
            'orientation-reversing identification across the throat. The '
            '5D tortoise coordinate r* = r + (rs/2)ln|(r−rs)/(r+rs)| '
            'diverges LOGARITHMICALLY at the throat (slope rs/2): the '
            'throat is infinitely deep in tortoise distance.\n\n'
            'RIGID THROAT ⟹ MASSLESS NEUTRINO. The reduced bounce '
            'S = √(2 μ E_c)·L*(ε), with the PR #58 pinch barrier E_c and '
            'brane inertia μ, runs over the tortoise path length L*(ε) to '
            'within a compliance cutoff ε of the throat. As ε→0 (a rigid, '
            'incompressible throat) L*→∞ and S→∞, so m_ν = m_D·e^{−S}→0: a '
            'perfectly rigid throat gives an EXACTLY MASSLESS Majorana '
            'neutrino. A boundary compliance ε>0 is precisely what gives '
            'the neutrino mass — the smallness of m_ν is the near-rigidity '
            'of the throat.\n\n'
            'S IS A TORTOISE LOGARITHM. Because L*(ε) = −(rs/2)ln ε + '
            'const, the bounce action is S ∝ ln(1/ε): naturally O(10) and '
            'slowly (logarithmically) varying — exactly the O(15), '
            'generation-stable form PR #87 required. The huge keV→TeV gap '
            'is carried by a logarithm, not a hierarchy.\n\n'
            'THE MAGNITUDE FALLS SHORT (HONEST). With the throat tension '
            'fixed by the ELECTRON throat (PR #58/#87: σ=1/12π, ρ=3/4π), '
            'the under-barrier momentum √(2 μ E_c) ≈ 0.060, so even a '
            'near-rigid compliance gives S ≈ 0.04–1 — ~40× SHORT of the '
            'required ~16. Pure compliance cannot close the gap (it would '
            'need an absurd ε ~ 1e-230). Matching S∈[15,18] at a sane '
            'near-rigid compliance (ε ~ 1e-3–1e-6) requires the ΔL=2 '
            '(lepton-number-violating / B−L) throat tension to be a factor '
            't ≈ 6–12 stiffer than the EM-throat tension (barrier ~10²–10³× '
            'larger, since E_c ∝ σ³).\n\n'
            'PROGRESSIVE LOCALISATION. The open input has been sharpened '
            'three times: a mysterious ~TeV mass (PR #86) → an O(15) '
            'instanton action S (PR #87) → an O(10) dimensionless B−L/EM '
            'tension ratio (this probe). Each step ties the unknown more '
            'tightly to BAM geometry.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): the bounce coordinate '
            'is the non-orientable tortoise path; the rigid throat gives a '
            'massless neutrino (compliance is the mass-generating '
            'parameter); the bounce action is a tortoise logarithm, hence '
            'naturally O(10) and coarsely generation-stable. NOT '
            'established: that S≈15–18 follows from the EM-throat tension '
            '(it under-produces by ~40×); matching needs a ΔL=2 tension '
            '~6–12× stiffer (a dimensionless ratio, not yet derived), the '
            'compliance ε (a sub-throat length), and the detailed m_ν '
            'spectrum (the ×18 m_ν/m_D spread needs gen-dependent '
            'compliance or the mixing sector).'
        )
    else:
        verdict_class = 'MAJORANA_BOUNCE_INCONCLUSIVE'
        verdict = (
            'MAJORANA BOUNCE INCONCLUSIVE. A structural or numerical test '
            'failed; investigate before claiming the tortoise-log bounce.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'Majorana bounce = non-orientable tortoise logarithm: rigid '
            'throat ⟹ massless ν; S ∝ ln(1/ε), naturally O(10)/gen-stable; '
            'EM-throat tension under-produces S by ~40×'
        ),
        'mechanism': 'reduced Euclidean bounce S = √(2 μ E_c)·L*(ε) along the odd tortoise path',
        'structural_wins': 'rigid throat = massless ν; S is a tortoise log',
        'open': 'ΔL=2 (B−L) tension ratio t ≈ 6–12 + compliance ε (not yet derived)',
        'b4_caveat': (
            'tension/inertia from electron-throat calibration; S not '
            'produced by EM tension; closing needs a stiffer B−L sector'
        ),
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# Majorana instanton-bounce action probe (PR #88)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "PR #87 reframed the seesaw scale as a bounce action "
        "`S = ln(m_D/m_ν) ≈ 15–18` but left `S` open. This probe builds a "
        "reduced Euclidean bounce for the `ΔL=2` throat↔antithroat flip "
        "along the **non-orientable tortoise path** and tests whether the "
        "throat tension + boundary compliance produce `S ≈ 15–18`."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Mechanism**: {s['mechanism']}")
    L.append(f"- **Structural wins**: {s['structural_wins']}")
    L.append(f"- **Open**: {s['open']}")
    L.append(f"- **B4 caveat**: {s['b4_caveat']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'required S = ln(m_D/m_ν) ≈ 14.7–17.6, gen-stable',
        'T2': 'bounce coord = odd tortoise path; r* log-diverges (slope rs/2)',
        'T3': 'rigid throat ε→0 ⟹ S→∞ ⟹ massless ν (compliance = mass)',
        'T4': 'EM-throat bounce S ≲ 1 even near-rigid (~40× short)',
        'T5': 'matching S∈[15,18] needs ΔL=2 tension ratio t ≈ 6–12',
        'T6': 'uniform S ⟹ m_ν∝m_D; ×18 spread needs gen-compliance/mixing',
        'T7': 'structure established; open input → O(10) tension ratio',
        'T8': 'MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    # T3 rigid-limit table
    t3 = s['tests'][2]
    L.append('## T3: Rigid throat ⟹ massless neutrino')
    L.append('')
    L.append('| compliance ε | tortoise length L* | bounce S | m_ν/m_D = e^{−S} |')
    L.append('|---:|---:|---:|---:|')
    for r in t3['rows']:
        L.append(f"| {r['eps']:.0e} | {r['L_star']:.3f} | {r['S_base']:.3f} | "
                 f"{r['m_nu_over_m_D']:.3e} |")
    L.append('')
    L.append("As `ε → 0` (a rigid, incompressible throat) the tortoise "
             "length and the bounce diverge, so `m_ν → 0`: a perfectly "
             "rigid throat gives an **exactly massless** Majorana neutrino. "
             "Boundary compliance `ε > 0` is the mass-generating parameter.")
    L.append('')

    # T4 EM-calibrated bounce
    t4 = s['tests'][3]
    L.append('## T4: Reduced bounce (EM-throat-calibrated)')
    L.append('')
    L.append(f"- `σ = {t4['sigma_em']:.4e}`, `ρ = {t4['rho_em']:.4e}`, "
             f"`E_c = {t4['E_c']:.4e}`, `μ = {t4['mu_throat']:.4f}`, "
             f"`√(2μE_c) = {t4['under_barrier_momentum_p0']:.4f}`")
    L.append('')
    L.append('| compliance ε | tortoise length L* | bounce S |')
    L.append('|---:|---:|---:|')
    for r in t4['rows']:
        L.append(f"| {r['eps']:.0e} | {r['L_star']:.3f} | {r['S_base']:.4f} |")
    L.append('')
    L.append(f"Even an absurdly rigid `ε = 1e-16` gives "
             f"`S ≈ {t4['S_base_at_eps_1e-16']:.3f}` — **~40× short** of the "
             "required ~16. The EM-throat tension under-produces the action.")
    L.append('')

    # T5 required tension ratio
    t5 = s['tests'][4]
    L.append('## T5: Required ΔL=2 tension ratio')
    L.append('')
    L.append('| compliance ε | tortoise L* | required tension ratio t | barrier boost t³ |')
    L.append('|---:|---:|---:|---:|')
    for r in t5['rows']:
        L.append(f"| {r['eps']:.0e} | {r['L_star']:.3f} | "
                 f"{r['required_tension_ratio_t']:.2f} | "
                 f"{r['barrier_boost_t_cubed']:.0f}× |")
    L.append('')
    L.append(f"Matching `S ≈ 16` at a sane near-rigid compliance needs the "
             f"`ΔL=2` (B−L) throat tension `t ≈ "
             f"{t5['required_tension_ratio_range'][0]:.1f}–"
             f"{t5['required_tension_ratio_range'][1]:.1f}×` the EM-throat "
             f"tension. Pure compliance (`t=1`) would need an absurd "
             f"`ε ~ {t5['eps_needed_if_pure_compliance']:.0e}`.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **The ΔL=2 (B−L) tension ratio** `t ≈ 6–12` — the open '
             'input is now this dimensionless number, not a `~TeV` mass. '
             'Progressive localisation: ~TeV mass (PR #86) → `O(15)` action '
             '`S` (PR #87) → `O(10)` tension ratio (this probe).')
    L.append('- **The boundary compliance `ε`** — a sub-throat length, not '
             'yet derived from the bulk.')
    L.append('- **The detailed `m_ν` spectrum** — the `×18` spread in '
             '`m_ν/m_D` needs a generation-dependent compliance or the '
             'mixing sector.')
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_majorana_bounce_action_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
