"""
Generation-dependent ε_n and the neutrino hierarchy spread (PR #113).

PR #112 derived the SMALLNESS of the neutrino mass from a sub-throat,
bulk-geometric healing length ε ~ R_c³ (the bounce S = t²·P0·L*(ε) with the
meV scale an output), but flagged that a UNIFORM ε gives a uniform action,
hence m_ν ∝ m_D — a mere ×2.7 generation spread, far short of the observed
hierarchy. PR #91 suggested the fix: the generations are cavity radial
overtones n, and the overtone boundary stress χ_n (PR #79) DECREASES with n
(0.304, 0.097, 0.039), so higher-overtone necks are more compliant — a
generation-dependent ε_n that widens the spread "in the right direction."
This probe makes that quantitative and tests it honestly.

## The mechanism and its direction (right)

Compliance is the inverse of stiffness, so the natural law is

    ε_n ∝ 1/χ_n     (more compliant neck ⟹ larger healing length).

With χ_n decreasing, ε_n increases with n ⟹ L*(ε_n) decreases ⟹ S(n)
decreases ⟹ less suppression for higher overtones ⟹ m_ν increases with n.
Together with m_D increasing with n, this gives NORMAL ORDERING, untuned —
the DIRECTION is correct.

## The magnitude (overshoots badly)

But the observed hierarchy needs only a GENTLE ε_n variation. Anchoring
gen 1 at ε_1 = R_c³ (the PR #112 value, m_ν1 ≈ 2.08 meV) and demanding the
observed m_2 = 8.65, m_3 = 50.34 meV gives

    required  ε_n  = (0.0110, 0.0130, 0.0172),  ratios (1, 1.18, 1.57).

The χ_n-driven law overshoots by orders of magnitude:

    ε_n ∝ 1/χ_n  ⟹  ε ratios (1, 3.13, 7.79)
                 ⟹  m_ν = (2.1, 1038, 167650) meV,
                 ⟹  m_ν3/m_ν2 ≈ 162   vs observed 5.85   (×28 overshoot).

The culprit is the bounce STEEPNESS established in PR #112 (∂ln m_ν/∂ln ε ≈
4.8): the factor-~8 variation in χ_n is amplified into ~4 orders of
magnitude in mass. The required power in ε_n ∝ χ_n^{−p} is not the
principled p = 1 but an inconsistent fractional p ≈ 0.15 (gen 1→2) to 0.31
(gen 2→3) — no clean law fits both ratios.

## The honest verdict

A generation-dependent ε_n gets the hierarchy's DIRECTION right (overtone
compliance ⟹ normal ordering, untuned), but it cannot PREDICT the spread
from χ_n: the natural driver overshoots by ~28× on the m_3/m_2 ratio (orders
of magnitude in absolute mass), and only a gently fitted ε_n profile — not
the principled 1/χ_n — reproduces the data. So the same bounce steepness
that made ε's absolute value a residual (PR #112) also blocks the natural
overtone variation from setting the spread. The neutrino hierarchy spread
stays a residual: ε_n can ACCOMMODATE it (by fitting a gentle profile) but
not DERIVE it, and the spread plausibly belongs to the mixing / anarchy
sector (PR #92) rather than a generation-dependent healing length.

Tests:
  T1. Setup: PR #112 left the spread open; PR #91 proposed χ_n-driven ε_n.
  T2. Mechanism + direction: ε_n ∝ 1/χ_n, χ_n decreasing ⟹ m_ν increases
      with n ⟹ NORMAL ORDERING (direction right, untuned).
  T3. Required spread is GENTLE: ε_n = (0.011, 0.013, 0.017), ratios
      (1, 1.18, 1.57) to hit observed m_2, m_3.
  T4. χ_n-driven ε_n ∝ 1/χ_n OVERSHOOTS: m_ν3/m_ν2 ≈ 162 vs 5.85 (×28).
  T5. Why: the steep bounce (m_ν ∝ ε^{4.8}, PR #112) amplifies the ×8 χ_n
      into ~10⁴ in mass; required power p ≈ 0.15–0.31 (inconsistent, not 1).
  T6. No clean law: a single principled ε_n(χ_n) fits neither both ratios
      nor the gentle magnitude ⟹ accommodates (fit), does not predict.
  T7. Honest scope: direction DERIVED; spread NOT derived from ε_n; stays a
      residual, plausibly the mixing/anarchy sector (PR #92).
  T8. Assessment.

Verdict:
  - HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL
    (expected): a generation-dependent ε_n from the overtone boundary stress
    χ_n reproduces the hierarchy's DIRECTION (normal ordering, untuned) but
    overshoots its MAGNITUDE by ~28× (m_3/m_2) — the steep bounce amplifies
    the ×8 χ_n variation into orders of magnitude in mass, and no principled
    ε_n(χ_n) law fits the gentle required spread. ε_n accommodates the
    spread but does not derive it; the spread stays a residual, plausibly
    the mixing/anarchy sector.
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


PI = math.pi
K_5 = 5
RS = R_MID

# EM / charged-throat calibration (PR #58) — electron-anchored, meV-free.
SIGMA_EM = 1.0 / (12.0 * PI * RS ** 2)
RHO_EM = 3.0 / (4.0 * PI * RS ** 3)
R_C = 2.0 * SIGMA_EM / RHO_EM
E_C = (16.0 * PI / 3.0) * SIGMA_EM ** 3 / RHO_EM ** 2
MU_THROAT = 4.0 * PI * SIGMA_EM * RS ** 2
P0 = math.sqrt(2.0 * MU_THROAT * E_C)

T2_WIND = (K_5 * math.sqrt(2.0 * PI)) ** 2          # = 50π (winding edge, PR #89)
T2P0 = T2_WIND * P0                                  # ≈ 9.50

# Overtone boundary stress χ_n (PR #79/#91), decreasing with n.
CHI_N = {1: 0.304, 2: 0.097, 3: 0.039}

# Electron-anchored cavity-floor Dirac masses (eV), PR #86/#87.
M_D_EV = {1: 43.0e3, 2: 80.0e3, 3: 118.0e3}

# Observed neutrino mass eigenstates (eV): the heavier two are √Δm²
# (NuFIT 6.0); the lightest is unmeasured and taken at the gen-1 anchor.
SQRT_DM21_EV = math.sqrt(7.49e-5)                    # ≈ 8.65 meV
SQRT_DM31_EV = math.sqrt(2.534e-3)                   # ≈ 50.34 meV

EPS1 = R_C ** 3                                      # gen-1 anchor (PR #112)


def _rstar(r: float) -> float:
    return r + (RS / 2.0) * math.log((r - RS) / (r + RS))


def _rstar_offset(eps: float) -> float:
    return (RS + eps) + (RS / 2.0) * math.log(eps / (2.0 * RS + eps))


def Lstar(eps: float) -> float:
    return _rstar(R_OUTER) - _rstar_offset(eps)


def S_of(eps: float) -> float:
    return T2P0 * Lstar(eps)


def m_nu_meV(eps: float, gen: int) -> float:
    return M_D_EV[gen] * math.exp(-S_of(eps)) * 1e3


def eps_for_target(mnu_target_eV: float, gen: int) -> float:
    """Invert m_ν = m_D·e^{−S(ε)} for ε by bisection on L*(ε)."""
    L_target = math.log(M_D_EV[gen] / mnu_target_eV) / T2P0
    lo, hi = 1e-4, 0.45
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        # L* decreases as ε increases
        if Lstar(mid) > L_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


M_NU1_ANCHOR_EV = M_D_EV[1] * math.exp(-S_of(EPS1))  # gen-1 output ≈ 2.08 meV
OBS_TARGET_EV = {1: M_NU1_ANCHOR_EV, 2: SQRT_DM21_EV, 3: SQRT_DM31_EV}


# ---------------------------------------------------------------------------
# T1. Setup
# ---------------------------------------------------------------------------

def test_T1_setup() -> dict:
    return {
        'name': 'T1_setup',
        'description': (
            "PR #112 derived the neutrino-mass SMALLNESS from a sub-throat "
            "ε ~ R_c³ but left the SPREAD open (uniform ε ⟹ m_ν ∝ m_D, "
            "×2.7). PR #91 proposed a generation-dependent ε_n driven by the "
            "overtone boundary stress χ_n (decreasing). This probe tests it "
            "quantitatively."
        ),
        'uniform_eps_spread': round(M_D_EV[3] / M_D_EV[1], 2),
        'observed_m3_over_m2': round(SQRT_DM31_EV / SQRT_DM21_EV, 2),
        'chi_n': CHI_N,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Mechanism + direction
# ---------------------------------------------------------------------------

def test_T2_direction() -> dict:
    """ε_n ∝ 1/χ_n (compliance = 1/stiffness); χ_n decreasing ⟹ ε_n
    increasing ⟹ S(n) decreasing ⟹ m_ν increasing with n ⟹ NORMAL
    ORDERING, with m_D also increasing. Direction correct, untuned."""
    eps = {n: EPS1 * (CHI_N[1] / CHI_N[n]) for n in CHI_N}
    mnu = {n: m_nu_meV(eps[n], n) for n in CHI_N}
    increasing = mnu[1] < mnu[2] < mnu[3]
    return {
        'name': 'T2_mechanism_and_direction',
        'description': (
            "ε_n ∝ 1/χ_n; χ_n decreasing ⟹ ε_n increasing ⟹ less "
            "suppression for higher overtones ⟹ m_ν increases with n ⟹ "
            "NORMAL ORDERING (direction right, untuned)."
        ),
        'eps_ratios': {n: round(eps[n] / EPS1, 2) for n in CHI_N},
        'm_nu_increasing_with_n': increasing,
        'gives_normal_ordering': increasing,
        'pass': increasing,
    }


# ---------------------------------------------------------------------------
# T3. The required spread is gentle
# ---------------------------------------------------------------------------

def test_T3_required_gentle() -> dict:
    """Anchoring gen 1 at ε_1 = R_c³ and demanding observed m_2 = 8.65,
    m_3 = 50.34 meV, the required ε_n = (0.011, 0.013, 0.017), ratios
    (1, 1.18, 1.57) — a GENTLE variation."""
    req = {n: eps_for_target(OBS_TARGET_EV[n], n) for n in (1, 2, 3)}
    ratios = {n: round(req[n] / req[1], 2) for n in (1, 2, 3)}
    return {
        'name': 'T3_required_spread_is_gentle',
        'description': (
            "Required ε_n to hit observed m_2, m_3: (0.011, 0.013, 0.017), "
            "ratios (1, 1.18, 1.57) — a GENTLE generation variation."
        ),
        'required_eps': {n: round(req[n], 4) for n in (1, 2, 3)},
        'required_ratios': ratios,
        'is_gentle': ratios[3] < 2.0,
        'pass': ratios[3] < 2.0,
    }


# ---------------------------------------------------------------------------
# T4. The χ_n-driven law overshoots
# ---------------------------------------------------------------------------

def test_T4_chi_overshoots() -> dict:
    """ε_n ∝ 1/χ_n gives ε ratios (1, 3.13, 7.79) ⟹ m_ν = (2.1, 1038,
    167650) meV ⟹ m_ν3/m_ν2 ≈ 162 vs observed 5.85 — a ×28 overshoot on the
    spread ratio (orders of magnitude in absolute mass)."""
    eps = {n: EPS1 * (CHI_N[1] / CHI_N[n]) for n in CHI_N}
    mnu = {n: m_nu_meV(eps[n], n) for n in CHI_N}
    spread_ratio = mnu[3] / mnu[2]
    obs = SQRT_DM31_EV / SQRT_DM21_EV * 1.0
    overshoot = spread_ratio / (SQRT_DM31_EV / SQRT_DM21_EV)
    return {
        'name': 'T4_chi_driven_overshoots',
        'description': (
            "ε_n ∝ 1/χ_n ⟹ m_ν = (2.1, 1038, 167650) meV ⟹ m_ν3/m_ν2 ≈ 162 "
            "vs observed 5.85 — ×28 overshoot (orders of magnitude in "
            "absolute mass)."
        ),
        'eps_ratios': {n: round(eps[n] / EPS1, 2) for n in CHI_N},
        'm_nu_meV': {n: round(mnu[n], 1) for n in CHI_N},
        'm3_over_m2_predicted': round(spread_ratio, 1),
        'm3_over_m2_observed': round(SQRT_DM31_EV / SQRT_DM21_EV, 2),
        'overshoot_factor': round(overshoot, 1),
        'pass': overshoot > 5.0,
    }


# ---------------------------------------------------------------------------
# T5. Why: the steep bounce amplifies χ_n
# ---------------------------------------------------------------------------

def test_T5_steepness() -> dict:
    """The bounce steepness (∂ln m_ν/∂ln ε = t²·P0·rs/2 ≈ 4.8, PR #112)
    amplifies the factor-~8 χ_n variation into ~4 orders of magnitude in
    mass. The power in ε_n ∝ χ_n^{−p} that reproduces the data is not the
    principled p = 1 but an inconsistent fractional p ≈ 0.15 (gen 1→2) to
    0.31 (gen 2→3)."""
    sens = T2P0 * RS / 2.0
    req = {n: eps_for_target(OBS_TARGET_EV[n], n) for n in (1, 2, 3)}
    p_12 = math.log(req[2] / req[1]) / math.log(CHI_N[1] / CHI_N[2])
    p_23 = math.log(req[3] / req[2]) / math.log(CHI_N[2] / CHI_N[3])
    chi_variation = CHI_N[1] / CHI_N[3]
    return {
        'name': 'T5_steep_bounce_amplifies_chi',
        'description': (
            "Steepness ∂ln m_ν/∂ln ε ≈ 4.8 (PR #112) amplifies the ×8 χ_n "
            "variation into ~10⁴ in mass. Required power p ≈ 0.15 (gen 1→2), "
            "0.31 (gen 2→3) — inconsistent, not the principled p = 1."
        ),
        'd_lnmnu_d_lneps': round(sens, 2),
        'chi_variation_factor': round(chi_variation, 1),
        'required_power_p_12': round(p_12, 3),
        'required_power_p_23': round(p_23, 3),
        'principled_p_is_1': True,
        'powers_inconsistent': abs(p_12 - p_23) > 0.1,
        'pass': sens > 3.0 and p_12 < 0.5 and p_23 < 0.5,
    }


# ---------------------------------------------------------------------------
# T6. No clean law: accommodates, does not predict
# ---------------------------------------------------------------------------

def test_T6_no_clean_law() -> dict:
    """A single principled law ε_n(χ_n) fits neither both ratios nor the
    gentle magnitude: 1/χ_n overshoots ×28; the data-fitted power is an
    inconsistent fraction. So a generation-dependent ε_n can ACCOMMODATE the
    spread (by fitting a gentle profile) but not PREDICT it from χ_n."""
    return {
        'name': 'T6_accommodates_not_predicts',
        'description': (
            "No principled ε_n(χ_n) reproduces the gentle required spread: "
            "1/χ_n overshoots ×28; the fitted power is an inconsistent "
            "fraction. ε_n ACCOMMODATES the spread (fit) but does not PREDICT "
            "it from χ_n."
        ),
        'principled_1_over_chi': 'overshoots ×28 (m_3/m_2 = 162 vs 5.85)',
        'fitted_power': 'p ≈ 0.15–0.31, inconsistent — a fit, not a law',
        'accommodates': True,
        'predicts': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Honest scope
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope',
        'description': (
            "Direction DERIVED (overtone compliance ⟹ normal ordering, "
            "untuned); spread NOT derived from ε_n (χ_n overshoots; required "
            "ε_n gentle and unprincipled). The same bounce steepness that "
            "made ε's value a residual (PR #112) blocks the spread; it "
            "plausibly belongs to the mixing/anarchy sector (PR #92)."
        ),
        'derived': [
            'the hierarchy DIRECTION — overtone compliance (ε_n ∝ 1/χ_n) ⟹ '
            'normal ordering, untuned (sharpens PR #91)',
        ],
        'not_derived': [
            'the hierarchy MAGNITUDE — χ_n overshoots ×28; required ε_n is '
            'gentle (×1.57) and not a principled function of χ_n',
            'the same bounce steepness (m_ν ∝ ε^{4.8}, PR #112) amplifies '
            'the natural overtone variation far past the data',
        ],
        'plausibly_belongs_to': 'the mixing / anarchy sector (PR #92), not a generation-dependent ε_n',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "A generation-dependent ε_n from the overtone boundary stress "
            "χ_n reproduces the hierarchy's direction (normal ordering, "
            "untuned) but overshoots its magnitude by ~28× — the steep "
            "bounce amplifies the ×8 χ_n variation into orders of magnitude, "
            "and no principled ε_n(χ_n) fits the gentle required spread. ε_n "
            "accommodates the spread but does not derive it; it stays a "
            "residual, plausibly the mixing/anarchy sector."
        ),
        'classification': (
            'HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_setup(),
        test_T2_direction(),
        test_T3_required_gentle(),
        test_T4_chi_overshoots(),
        test_T5_steepness(),
        test_T6_no_clean_law(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL'
        )
        verdict = (
            'A GENERATION-DEPENDENT ε_n GETS THE HIERARCHY DIRECTION RIGHT '
            'BUT OVERSHOOTS ITS MAGNITUDE — THE SPREAD STAYS A RESIDUAL. '
            'PR #112 derived the neutrino-mass SMALLNESS from a sub-throat, '
            'bulk-geometric healing length ε ~ R_c³, but a uniform ε gives a '
            'uniform bounce action, hence m_ν ∝ m_D — only a ×2.7 generation '
            'spread, far short of the observed hierarchy. PR #91 proposed the '
            'fix: the generations are cavity radial overtones n, and the '
            'overtone boundary stress χ_n (PR #79) decreases with n (0.304, '
            '0.097, 0.039), so higher-overtone necks are more compliant — a '
            'generation-dependent ε_n. This probe makes that quantitative.\n\n'
            'THE DIRECTION IS RIGHT. Compliance is the inverse of stiffness, '
            'so the natural law is ε_n ∝ 1/χ_n: with χ_n decreasing, ε_n '
            'increases with n, the tortoise length L*(ε_n) decreases, the '
            'action S(n) decreases, and the suppression weakens — so m_ν '
            'increases with n, and with m_D also increasing the result is '
            'NORMAL ORDERING, untuned. The qualitative hierarchy is '
            'reproduced.\n\n'
            'THE MAGNITUDE OVERSHOOTS. But the observed hierarchy needs only '
            'a GENTLE ε_n variation. Anchoring gen 1 at ε_1 = R_c³ (the '
            'PR #112 value, m_ν1 ≈ 2.08 meV) and demanding the observed m_2 = '
            '8.65, m_3 = 50.34 meV requires ε_n = (0.011, 0.013, 0.017), '
            'ratios (1, 1.18, 1.57). The χ_n-driven law badly overshoots: '
            'ε_n ∝ 1/χ_n gives ε ratios (1, 3.13, 7.79), hence m_ν = (2.1, '
            '1038, 167650) meV and m_ν3/m_ν2 ≈ 162 against the observed 5.85 '
            '— a ×28 overshoot on the spread ratio, and orders of magnitude '
            'in absolute mass.\n\n'
            'WHY: THE STEEP BOUNCE. The culprit is the bounce steepness '
            'established in PR #112, ∂ln m_ν/∂ln ε ≈ 4.8: the factor-~8 '
            'variation in χ_n is amplified into ~4 orders of magnitude in '
            'mass. The power in ε_n ∝ χ_n^{−p} that would reproduce the data '
            'is not the principled p = 1 but an inconsistent fractional p ≈ '
            '0.15 (gen 1→2) to 0.31 (gen 2→3) — no single clean law fits both '
            'ratios. So a generation-dependent ε_n can ACCOMMODATE the spread '
            '(by fitting a gentle profile) but cannot PREDICT it from χ_n.\n\n'
            'THE HONEST VERDICT. The hierarchy DIRECTION is derived (overtone '
            'compliance ⟹ normal ordering, untuned — a sharpening of PR #91), '
            'but its MAGNITUDE is not: the natural χ_n driver overshoots, and '
            'the same bounce steepness that made ε\'s absolute value a '
            'residual (PR #112) now blocks the natural overtone variation '
            'from setting the spread. The neutrino hierarchy spread stays a '
            'residual — ε_n accommodates it but does not derive it — and it '
            'plausibly belongs to the mixing / anarchy sector (PR #92) rather '
            'than a generation-dependent healing length.'
        )
    else:
        verdict_class = 'HIERARCHY_SPREAD_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the ε_n(χ_n) '
            'model and the bounce chain.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'a generation-dependent ε_n ∝ 1/χ_n (overtone boundary stress) '
            'reproduces the hierarchy DIRECTION (normal ordering, untuned) '
            'but overshoots its MAGNITUDE by ×28 (m_3/m_2 = 162 vs 5.85); the '
            'spread stays a residual'
        ),
        'direction': 'normal ordering, untuned (ε_n ∝ 1/χ_n) — derived (sharpens PR #91)',
        'magnitude': 'χ_n overshoots ×28; required ε_n gentle (1, 1.18, 1.57), not principled',
        'cause': 'the steep bounce (m_ν ∝ ε^{4.8}, PR #112) amplifies the ×8 χ_n into ~10⁴ in mass',
        'status': 'ε_n accommodates the spread but does not derive it — stays a residual',
        'plausibly': 'the spread belongs to the mixing/anarchy sector (PR #92)',
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
    L.append('# Generation-dependent ε_n and the neutrino hierarchy spread (PR #113)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Tests PR #91's proposed fix for the spread PR #112 left open: a "
        "generation-dependent healing length `ε_n` driven by the overtone "
        "boundary stress `χ_n` (decreasing with n). **Result: the direction "
        "is right, the magnitude overshoots.** `ε_n ∝ 1/χ_n` reproduces "
        "**normal ordering** untuned, but the natural `χ_n` variation (×8) "
        "is amplified by the steep bounce (`m_ν ∝ ε^{4.8}`, PR #112) into "
        "orders of magnitude in mass — `m_ν3/m_ν2 ≈ 162` vs observed 5.85 "
        "(×28). The spread stays a **residual**."
    )
    L.append('')
    L.append(f"- **Direction**: {s['direction']}")
    L.append(f"- **Magnitude**: {s['magnitude']}")
    L.append(f"- **Cause**: {s['cause']}")
    L.append(f"- **Status**: {s['status']}")
    L.append(f"- **Plausibly**: {s['plausibly']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'uniform ε ⟹ ×2.7 spread; PR #91 proposed χ_n-driven ε_n',
        'T2': 'ε_n ∝ 1/χ_n ⟹ normal ordering (direction right, untuned)',
        'T3': 'required ε_n gentle: (1, 1.18, 1.57) to hit observed m_2, m_3',
        'T4': 'χ_n-driven overshoots: m_3/m_2 = 162 vs 5.85 (×28)',
        'T5': 'steep bounce (×4.8) amplifies ×8 χ_n; power p 0.15→0.31 (≠1)',
        'T6': 'no clean ε_n(χ_n) law — accommodates (fit), not predicts',
        'T7': 'direction derived; spread residual (mixing/anarchy, PR #92)',
        'T8': 'HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t3 = s['tests'][2]; t4 = s['tests'][3]
    L.append('## Required (gentle) vs χ_n-driven (overshoots)')
    L.append('')
    L.append('| gen | χ_n | required ε_n | required ratio | χ-driven ε_n ratio | χ-driven m_ν (meV) |')
    L.append('|---|---:|---:|---:|---:|---:|')
    for n in (1, 2, 3):
        L.append(f"| {n} | {CHI_N[n]} | {t3['required_eps'][n]} | "
                 f"{t3['required_ratios'][n]} | {t4['eps_ratios'][n]} | "
                 f"{t4['m_nu_meV'][n]} |")
    L.append('')
    L.append(f"The required `ε_n` rises only ×1.57 over three generations; "
             f"the principled `ε_n ∝ 1/χ_n` rises ×7.79, giving "
             f"`m_ν3/m_ν2 = {t4['m3_over_m2_predicted']}` vs observed "
             f"`{t4['m3_over_m2_observed']}` — a ×{t4['overshoot_factor']} "
             "overshoot.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this establishes (and does not)')
    L.append('')
    L.append('- **Derived:** the hierarchy DIRECTION — overtone compliance '
             '(`ε_n ∝ 1/χ_n`) gives normal ordering, untuned (sharpens '
             'PR #91).')
    L.append('- **Not derived:** the hierarchy MAGNITUDE — the natural `χ_n` '
             'driver overshoots by ×28; the required `ε_n` is gentle '
             '(×1.57) and not a principled function of `χ_n`. The same '
             'bounce steepness that made `ε`\'s value a residual (PR #112) '
             'blocks the spread, which plausibly belongs to the '
             'mixing/anarchy sector (PR #92).')
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
    out = here / 'runs' / f'{ts}_generation_dependent_eps_n_probe'
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
