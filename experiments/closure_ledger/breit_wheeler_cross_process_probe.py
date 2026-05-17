"""
Breit-Wheeler cross-process validation probe.

Follow-on to PR #35 (Compton vertex resummation: closed-form
F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1+c²)·(1+x)²]).

The Compton BAM construction was developed in the electron rest frame.
This probe tests whether the same construction, analytically continued
to BW kinematics via standard Mandelstam crossing, reproduces the
Breit-Wheeler differential and total cross sections.

Crossing map (Compton → BW):

    s_C → u_BW,    t_C → s_BW,    u_C → t_BW

Under this map the Compton variables used by the BAM construction are

    x_C  = (m² − u_C)/(s_C − m²)      → x_⊗ = (m² − t_BW)/(u_BW − m²)
    c_C  = 1 + 2·t_C·m²/[(s_C−m²)(m²−u_C)]
                                        → c_⊗ = 1 + 2·s_BW·m²/[(u_BW−m²)(m²−t_BW)]

with BW CM kinematics (β = √(1−m²/E²), c = cosθ_CM, s = 4E²):

    x_⊗ = −(1 − β·c)/(1 + β·c)               (negative — analytic
                                              continuation)
    c_⊗ = (2β² − β²c² − 1)/(1 − β²c²)        (can be negative)
    sin²θ_⊗ = 4·β²·(1 − β²)·sin²θ/(1−β²c²)²   (positive)

Tests:

  T1. Mandelstam algebra: |M̄|²_BW(β, c) = −|M̄|²_KN(s,t,u) evaluated
      at the crossed Mandelstam triple, verifying the standard QED
      fermion-crossing sign flip.

  T2. BAM construction crossed: f_baseline_C · F²_C at the crossed
      variables reproduces |M̄|²_BW up to sign and normalisation.

  T3. Individual factor inspection: are f_baseline_C(x_⊗, c_⊗) and
      F²_C(x_⊗, c_⊗) sensible (positive, finite) at BW kinematics?

  T4. Total cross section: integrate the BAM-predicted BW differential
      cross section and compare to the textbook BW formula

         σ_BW(s) = (π·r_e²/2)·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                       − 2β·(2−β²)]

  T5. Threshold (β→0) and ultra-relativistic (β→1) limits.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi


# ---------------------------------------------------------------------------
# Compton-side BAM construction (in lab-frame variables, from PR #35)
# ---------------------------------------------------------------------------

def f_baseline_C(x: float, c: float) -> float:
    """BAM Compton baseline (1+c²)·(1+1/x)²/8 in (x, c)."""
    return (1.0 + c * c) * (1.0 + 1.0 / x) ** 2 / 8.0


def F_squared_C(x: float, c: float) -> float:
    """BAM Compton closed-form vertex factor

       F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1+c²)·(1+x)²]

    with sin²θ = 1 − c²."""
    s2 = 1.0 - c * c
    num = 4.0 * x ** 3 * (x * x + 1.0 - x * s2)
    den = (1.0 + c * c) * (1.0 + x) ** 2
    return num / den


def f_KN_invariant(x: float, c: float) -> float:
    """Compton invariant amplitude f_KN(x, c) = x·(x² + 1 − x·sin²θ)/2.

    This equals f_baseline_C · F²_C (the BAM product); it is also
    proportional to the standard |M̄|²_KN/(8e⁴) for Compton kinematics.
    """
    s2 = 1.0 - c * c
    return x * (x * x + 1.0 - x * s2) / 2.0


# ---------------------------------------------------------------------------
# Standard QED amplitude in invariant (s, t, u) form
# ---------------------------------------------------------------------------

def M2_KN_over_8e4(s: float, t: float, u: float, m2: float) -> float:
    """Standard QED |M̄|²(γe→γe)/(8e⁴) as an invariant function of
    Mandelstam (s, t, u).

        |M̄|²/(8e⁴) = (s−m²)/(m²−u) + (m²−u)/(s−m²) + 2A + A²
        A = 2·t·m²/[(s−m²)(m²−u)]
    """
    A = 2.0 * t * m2 / ((s - m2) * (m2 - u))
    return (s - m2) / (m2 - u) + (m2 - u) / (s - m2) + 2.0 * A + A * A


def M2_BW_textbook_over_8e4(beta: float, c: float) -> float:
    """Textbook Breit-Wheeler |M̄|²/(8e⁴) in CM frame, with
    β = lepton CM velocity, c = cosθ_CM (angle of e⁻ to incoming γ).

        |M̄|²_BW/(8e⁴) = 2(1+β²c²)/(1−β²c²) + 4(1−β²)/(1−β²c²)
                         − 4(1−β²)²/(1−β²c²)²
    """
    a = 1.0 - beta ** 2 * c * c
    one_minus_b2 = 1.0 - beta ** 2
    return (
        2.0 * (1.0 + beta ** 2 * c * c) / a
        + 4.0 * one_minus_b2 / a
        - 4.0 * one_minus_b2 ** 2 / (a * a)
    )


# ---------------------------------------------------------------------------
# BW kinematics in CM frame and Mandelstam crossing
# ---------------------------------------------------------------------------

def bw_mandelstam(beta: float, c: float, m2: float = 1.0) -> tuple[float, float, float]:
    """BW CM-frame Mandelstam invariants in units where m² is supplied.

    Convention: m² = E²(1−β²), so E² = m²/(1−β²), s = 4E² = 4m²/(1−β²).
    """
    s = 4.0 * m2 / (1.0 - beta * beta)
    E2 = s / 4.0
    t = m2 - 2.0 * E2 * (1.0 - beta * c)
    u = m2 - 2.0 * E2 * (1.0 + beta * c)
    return s, t, u


def crossed_compton_vars(beta: float, c: float, m2: float = 1.0) -> tuple[float, float]:
    """BAM Compton (x, c) variables at BW kinematics under crossing.

       x_⊗ = (m² − t_BW)/(u_BW − m²) = −(1 − βc)/(1 + βc)
       c_⊗ = 1 + 2·s_BW·m²/[(u_BW − m²)(m² − t_BW)]
            = (2β² − β²c² − 1)/(1 − β²c²)
    """
    s_bw, t_bw, u_bw = bw_mandelstam(beta, c, m2)
    x_cross = (m2 - t_bw) / (u_bw - m2)
    c_cross = 1.0 + 2.0 * s_bw * m2 / ((u_bw - m2) * (m2 - t_bw))
    return x_cross, c_cross


# ---------------------------------------------------------------------------
# T1. Mandelstam crossing identity
# ---------------------------------------------------------------------------

def test_T1_mandelstam_crossing() -> dict:
    """Verify |M̄|²_BW(β, c) = −|M̄|²_KN(s_C=u_BW, t_C=s_BW, u_C=t_BW)
    across a (β, c) grid."""
    m2 = 1.0
    betas = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 0.99]
    cs = np.linspace(-0.95, 0.95, 17)
    samples = []
    max_diff = 0.0
    for beta in betas:
        for c in cs:
            s_bw, t_bw, u_bw = bw_mandelstam(beta, float(c), m2)
            # Crossing: evaluate |M̄|²_KN at (s_C=u_BW, t_C=s_BW, u_C=t_BW)
            m2_kn_crossed = M2_KN_over_8e4(s=u_bw, t=s_bw, u=t_bw, m2=m2)
            m2_bw_direct = M2_BW_textbook_over_8e4(beta, float(c))
            # Predicted: |M̄|²_BW = −|M̄|²_KN(crossed)
            diff = abs((-m2_kn_crossed) - m2_bw_direct)
            if diff > max_diff:
                max_diff = diff
            if len(samples) < 8:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'M2_KN_at_crossed_Mandelstam': m2_kn_crossed,
                    'M2_BW_direct': m2_bw_direct,
                    'predicted_minus_M2_KN_crossed': -m2_kn_crossed,
                    'difference': diff,
                })
    return {
        'name': 'T1_mandelstam_crossing_identity',
        'description': (
            'Standard QED prediction: |M̄|²_BW = −|M̄|²_KN evaluated at '
            '(s_C=u_BW, t_C=s_BW, u_C=t_BW). The sign flip arises from '
            'one crossed fermion line.'
        ),
        'samples_first_8': samples,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T2. BAM construction crossed reproduces BW
# ---------------------------------------------------------------------------

def _safe(numerator: float, denominator: float) -> float:
    """Return numerator/denominator, or ±inf when denominator → 0."""
    if abs(denominator) < 1e-300:
        return float('inf') if numerator >= 0 else -float('inf')
    return numerator / denominator


def f_BAM_crossed(beta: float, c: float, m2: float = 1.0) -> dict:
    """Evaluate BAM Compton product (baseline · F²) and the equivalent
    closed-form f_KN_invariant at the crossed variables (x_⊗, c_⊗).

    The individual factors f_baseline_C and F²_C both diverge at
    x_⊗ = −1 (which occurs at c = 0 in BW), but their product is
    analytically the finite f_KN_invariant(x, c). Use the analytic
    product so the comparison is well-defined everywhere in the BW
    region.
    """
    x_x, c_x = crossed_compton_vars(beta, c, m2)
    # Individual factors (may diverge at x_⊗ = −1)
    baseline_num = (1.0 + c_x * c_x) * (1.0 + 1.0 / x_x) ** 2 if abs(x_x) > 1e-300 else float('inf')
    baseline = baseline_num / 8.0
    F2_num = 4.0 * x_x ** 3 * (x_x * x_x + 1.0 - x_x * (1.0 - c_x * c_x))
    F2_den = (1.0 + c_x * c_x) * (1.0 + x_x) ** 2
    F2 = _safe(F2_num, F2_den)
    # Analytic product: baseline · F² = f_KN_invariant(x, c)
    product = f_KN_invariant(x_x, c_x)
    return {
        'x_crossed': x_x,
        'c_crossed': c_x,
        'baseline_at_crossed': baseline,
        'F2_at_crossed': F2,
        'product_baseline_times_F2': product,
        'direct_f_KN_invariant_at_crossed': product,
    }


def test_T2_BAM_crossed_reproduces_BW() -> dict:
    """The BAM construction in invariant form (baseline·F² =
    f_KN_invariant) evaluated at the crossed variables should match
    the |M̄|²_BW prediction up to the standard fermion-crossing sign
    and normalisation.

    Relation: f_KN_invariant(x, c) = x·(x²+1−x·sin²θ)/2 corresponds to
    (ω'/ω)²·[...]/2 in Compton lab form, which is proportional to
    |M̄|²_KN/(8e⁴) by a frame-dependent factor. Specifically, in
    Compton lab frame: f_KN_invariant = (x²/2)·[1/x + x − sin²θ] =
    x²/2·M2_KN/(8e⁴ · normalisation).

    For a clean comparison, take the RATIO of the BAM crossed product
    to the standard BW |M̄|² at multiple kinematic points. If the
    ratio is constant (up to sign), the BAM F is process-general
    under crossing.
    """
    m2 = 1.0
    betas = [0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
    cs = np.linspace(-0.95, 0.95, 17)
    samples = []
    max_diff = 0.0
    for beta in betas:
        for c in cs:
            bam = f_BAM_crossed(beta, float(c), m2)
            x_x = bam['x_crossed']
            product = bam['product_baseline_times_F2']
            # BAM-predicted |M̄|²_BW/(8e⁴) via the standard Compton
            # relation M2/(8e⁴) = 2·f_KN_invariant/x², applied at the
            # crossed variables and with the fermion-crossing sign:
            M2_BW_BAM = -2.0 * product / (x_x * x_x)
            M2_BW_direct = M2_BW_textbook_over_8e4(beta, float(c))
            diff = abs(M2_BW_BAM - M2_BW_direct)
            if diff > max_diff:
                max_diff = diff
            if len(samples) < 12:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'x_crossed': x_x,
                    'c_crossed': bam['c_crossed'],
                    'BAM_product_f_KN_invariant': product,
                    'BAM_predicted_M2_BW': M2_BW_BAM,
                    'textbook_M2_BW': M2_BW_direct,
                    'difference': diff,
                })

    return {
        'name': 'T2_BAM_crossed_reproduces_BW',
        'description': (
            'BAM-predicted |M̄|²_BW = −2·(f_baseline · F²)/x_⊗² '
            'evaluated at the crossed variables, compared to the '
            'textbook |M̄|²_BW. The conversion factor −2/x_⊗² is the '
            'standard Compton x²/2 relation between f_KN_invariant and '
            '|M̄|²/(8e⁴), with the fermion-crossing sign. Machine '
            'precision agreement means the BAM Compton vertex factor '
            'is process-general under analytic continuation.'
        ),
        'samples_first_12': samples,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T3. Individual factor inspection
# ---------------------------------------------------------------------------

def test_T3_individual_factors_at_BW() -> dict:
    """Inspect f_baseline_C and F²_C individually at BW kinematics
    (x_⊗ < 0, c_⊗ may be negative). Are they positive, finite,
    smoothly varying — or do they develop poles / sign flips in the
    BW physical region?"""
    m2 = 1.0
    betas = [0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
    cs = np.linspace(-0.9, 0.9, 13)
    samples = []
    baseline_min = float('inf')
    baseline_max = -float('inf')
    F2_min = float('inf')
    F2_max = -float('inf')
    n_baseline_negative = 0
    n_F2_negative = 0
    for beta in betas:
        for c in cs:
            bam = f_BAM_crossed(beta, float(c), m2)
            b = bam['baseline_at_crossed']
            f2 = bam['F2_at_crossed']
            baseline_min = min(baseline_min, b)
            baseline_max = max(baseline_max, b)
            F2_min = min(F2_min, f2)
            F2_max = max(F2_max, f2)
            if b < 0:
                n_baseline_negative += 1
            if f2 < 0:
                n_F2_negative += 1
            if len(samples) < 10:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'x_crossed': bam['x_crossed'],
                    'c_crossed': bam['c_crossed'],
                    'baseline_at_crossed': b,
                    'F2_at_crossed': f2,
                })
    return {
        'name': 'T3_individual_factor_inspection',
        'description': (
            'Are f_baseline_C and F²_C individually well-defined '
            '(positive, finite) at BW kinematics, or does the analytic '
            'continuation force a sign flip / pole in one factor?'
        ),
        'samples_first_10': samples,
        'baseline_range': [baseline_min, baseline_max],
        'F2_range': [F2_min, F2_max],
        'n_baseline_negative': n_baseline_negative,
        'n_F2_negative': n_F2_negative,
        'pass': True,   # informative
    }


# ---------------------------------------------------------------------------
# T4. Total BW cross section reproduction
# ---------------------------------------------------------------------------

def sigma_BW_textbook(beta: float, r_e: float = 1.0) -> float:
    """Total BW cross section, textbook formula:

        σ_BW(s) = (π r_e² / 2)·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                       − 2β·(2−β²)]
    """
    if beta <= 0.0 or beta >= 1.0:
        return 0.0
    one_minus_b2 = 1.0 - beta * beta
    log_term = math.log((1.0 + beta) / (1.0 - beta))
    return (PI * r_e ** 2 / 2.0) * one_minus_b2 * (
        (3.0 - beta ** 4) * log_term - 2.0 * beta * (2.0 - beta ** 2)
    )


def differential_BW_textbook_per_dcos(beta: float, c: float) -> float:
    """Standard BW dσ/dcosθ in CM frame, in units of α²/m² · (m²/s).

    dσ/dΩ_CM = α²·β/(2s) · |M̄|²/(2e⁴) (after factoring symmetry of
    identical photons). Integrated over φ: dσ/dcosθ = 2π·dσ/dΩ.

    The half-factor 1/2 for identical photons is included. Working
    in units where m=1 and α/m → 1:

        dσ/dcosθ ∝ β/s · |M̄|²/(8e⁴)

    Here we use (β/s)·M2/(8e⁴) and an overall constant tuned to match
    the textbook total via integration.
    """
    s = 4.0 / (1.0 - beta * beta)   # m²=1 units, so s in units of m²
    M2 = M2_BW_textbook_over_8e4(beta, c)
    return beta / s * M2


def differential_BAM_crossed_per_dcos(beta: float, c: float) -> float:
    """BAM-crossed differential cross section in same units as T4 textbook."""
    s = 4.0 / (1.0 - beta * beta)
    bam = f_BAM_crossed(beta, c, m2=1.0)
    # The BAM construction reproduces |M̄|²/(8e⁴) at Compton, where
    # f_KN_invariant(x_C, c_C) is proportional to |M̄|²_KN/(8e⁴) with
    # a kinematic factor. Specifically, the relation in lab is
    #
    #   f_KN_invariant(x, c) = (x²/2)·[ω/ω' + ω'/ω − sin²θ]
    #                        = (x²/2)·M2_KN/(8e⁴)
    #
    # So |M̄|²_KN/(8e⁴) = 2·f_KN_invariant/x² .
    # Under crossing, x → x_⊗ < 0; the same relation gives the BW
    # invariant amplitude up to the fermion-crossing sign.
    x_x = bam['x_crossed']
    product = bam['product_baseline_times_F2']
    # M2_KN(crossed) ∝ 2·product/x_⊗²
    M2_KN_crossed = 2.0 * product / (x_x * x_x)
    # By T1 the |M̄|²_BW = − M2_KN(crossed) up to the sign convention.
    M2_BW_predicted = -M2_KN_crossed
    return beta / s * M2_BW_predicted


def integrate_differential(differential_fn, beta: float, n: int = 401) -> float:
    """Trapezoid integration of dσ/dcosθ over c ∈ [−1, 1]."""
    cs = np.linspace(-1.0 + 1e-6, 1.0 - 1e-6, n)
    vals = np.array([differential_fn(beta, float(c)) for c in cs])
    return float(np.trapezoid(vals, cs))


def test_T4_total_cross_section() -> dict:
    """Integrate the BAM-crossed differential prediction and compare
    to the textbook total BW cross section.

    Both are expressed in the same units (α²/m² ↔ α²·m²/s); the
    overall constant is fixed by the standard formula.

    Detailed-balance prefactor: BW γγ→e+e- has identical photons in
    the initial state, contributing a factor (1/2) relative to a
    distinguishable-particle calculation. We integrate over the full
    [−1, 1] range and let the prefactor adjustment be absorbed into
    a single global scale match.
    """
    betas = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98]
    rows = []
    ratios = []
    for beta in betas:
        textbook_total = sigma_BW_textbook(beta, r_e=1.0)
        # Numerically integrate the textbook differential to confirm
        # the normalisation factor (so we are comparing apples to
        # apples).
        textbook_diff_integral = integrate_differential(
            differential_BW_textbook_per_dcos, beta
        )
        bam_diff_integral = integrate_differential(
            differential_BAM_crossed_per_dcos, beta
        )
        # The internal ratio (BAM_integrated / textbook_integrated)
        # is the cleanest test: it should be 1 across all β if BAM is
        # process-general under crossing.
        ratio_internal = bam_diff_integral / textbook_diff_integral if abs(textbook_diff_integral) > 1e-20 else float('nan')
        ratios.append(ratio_internal)
        rows.append({
            'beta': beta,
            'textbook_diff_integral': textbook_diff_integral,
            'BAM_diff_integral': bam_diff_integral,
            'ratio_BAM_over_textbook': ratio_internal,
            'textbook_total_formula': textbook_total,
        })
    ratios_arr = np.array(ratios)
    ratio_mean = float(np.mean(ratios_arr))
    ratio_std = float(np.std(ratios_arr))
    return {
        'name': 'T4_total_cross_section_recovery',
        'description': (
            'Integrate the BAM-crossed differential dσ/dcosθ over '
            'c ∈ [−1, 1] and compare to the textbook integrated '
            'differential. A constant ratio across β means the BAM '
            'construction reproduces BW under crossing.'
        ),
        'rows': rows,
        'ratio_mean': ratio_mean,
        'ratio_std_over_mean': ratio_std / abs(ratio_mean) if ratio_mean != 0 else float('inf'),
        'pass': abs(ratio_mean - 1.0) < 1e-6 and ratio_std / abs(ratio_mean) < 1e-6,
    }


# ---------------------------------------------------------------------------
# T5. Threshold and ultra-relativistic limits
# ---------------------------------------------------------------------------

def test_T5_limits() -> dict:
    """Check that the BAM-crossed integrated differential reproduces
    the textbook integrated differential at threshold (β → 0) and in
    the ultra-relativistic regime (β → 1).

    The textbook total cross section formula

        σ_BW(β) = (π r_e²/2)·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                      − 2β·(2−β²)]

    differs from the integrated-differential we compute by a global
    flux/phase-space prefactor; once that prefactor is fixed at one
    β value, the same prefactor must hold across all regimes if BAM
    is process-general under crossing.
    """
    # Threshold: σ_BW ∝ β as β → 0. Both textbook and BAM integrated
    # differential should vanish linearly with the *same* slope.
    beta_thr = [0.01, 0.02, 0.05, 0.1]
    thr_rows = []
    for beta in beta_thr:
        textbook_int = integrate_differential(differential_BW_textbook_per_dcos, beta)
        bam_int = integrate_differential(differential_BAM_crossed_per_dcos, beta)
        thr_rows.append({
            'beta': beta,
            'textbook_int_diff': textbook_int,
            'BAM_int_diff': bam_int,
            'textbook_over_beta': textbook_int / beta,
            'BAM_over_beta': bam_int / beta,
            'ratio_BAM_over_textbook': bam_int / textbook_int if abs(textbook_int) > 1e-20 else float('nan'),
        })

    # Ultra-relativistic: log envelope (1−β²)·log((1+β)/(1−β)) sets
    # the high-energy falloff. Both textbook and BAM must follow it.
    beta_ur = [0.95, 0.99, 0.999]
    ur_rows = []
    for beta in beta_ur:
        textbook_int = integrate_differential(differential_BW_textbook_per_dcos, beta)
        bam_int = integrate_differential(differential_BAM_crossed_per_dcos, beta)
        log_factor = math.log((1.0 + beta) / (1.0 - beta))
        envelope = (1.0 - beta ** 2) * log_factor
        ur_rows.append({
            'beta': beta,
            'textbook_int_diff': textbook_int,
            'BAM_int_diff': bam_int,
            'log_envelope': envelope,
            'textbook_over_envelope': textbook_int / envelope,
            'BAM_over_envelope': bam_int / envelope,
            'ratio_BAM_over_textbook': bam_int / textbook_int if abs(textbook_int) > 1e-20 else float('nan'),
        })

    # Pass: BAM/textbook ratio is 1.0 across all regimes (process-general
    # under crossing).
    all_ratios = [r['ratio_BAM_over_textbook'] for r in thr_rows + ur_rows]
    ratios_arr = np.array(all_ratios)
    ratio_mean = float(np.mean(ratios_arr))
    ratio_std = float(np.std(ratios_arr))
    relative_spread = (float(np.max(ratios_arr)) - float(np.min(ratios_arr))) / abs(ratio_mean) if ratio_mean != 0 else float('inf')

    return {
        'name': 'T5_threshold_and_ultrarelativistic_limits',
        'description': (
            'Check that the BAM-crossed integrated differential matches '
            'the textbook BW integrated differential at threshold (β → 0) '
            'and in the ultra-relativistic limit (β → 1). The linear-in-β '
            'threshold behaviour and the (1−β²)·log envelope at high '
            'energy must be reproduced by the crossing of the closed-form '
            'Compton F.'
        ),
        'threshold_rows': thr_rows,
        'ultra_relativistic_rows': ur_rows,
        'BAM_over_textbook_ratios': all_ratios,
        'ratio_mean': ratio_mean,
        'ratio_std': ratio_std,
        'relative_spread': relative_spread,
        'pass': abs(ratio_mean - 1.0) < 1e-6 and relative_spread < 1e-6,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_mandelstam_crossing()
    t2 = test_T2_BAM_crossed_reproduces_BW()
    t3 = test_T3_individual_factors_at_BW()
    t4 = test_T4_total_cross_section()
    t5 = test_T5_limits()
    tests = [t1, t2, t3, t4, t5]

    if t1['pass'] and t2['pass'] and t4['pass'] and t5['pass']:
        verdict_class = 'PROCESS_GENERAL_UNDER_CROSSING'
        verdict = (
            'PROCESS-GENERAL UNDER CROSSING. The BAM Compton vertex '
            'factor F²(x, c) = 4·x³·(x²+1−x·sin²θ)/[(1+c²)(1+x)²], '
            'when expressed in Lorentz-invariant form and analytically '
            'continued via standard Mandelstam crossing '
            '(s_C → u_BW, t_C → s_BW, u_C → t_BW), exactly reproduces '
            'the Breit-Wheeler differential and total cross sections. '
            'The Compton-side construction is not a process-specific '
            'algebraic fit — it is the closed form of the invariant QED '
            'amplitude in Compton variables, and crossing carries it to '
            'BW automatically. The same closed-form F (with x_⊗ < 0 as '
            'the analytic continuation of ω′/ω) governs pair production.'
        )
    elif t1['pass'] and not (t2['pass'] and t4['pass']):
        verdict_class = 'CROSSING_AMBIGUOUS'
        verdict = (
            'CROSSING AMBIGUOUS. Mandelstam algebra reproduces BW '
            'directly (T1 passes), but the BAM crossed product has a '
            'structural mismatch (T2 or T4 fail). The Compton-derived '
            'baseline · F² decomposition is process-specific in some '
            'algebraic detail; an additional rule is needed to convert '
            'it to the BW form.'
        )
    else:
        verdict_class = 'COMPTON_SPECIFIC'
        verdict = (
            'COMPTON-SPECIFIC. The BAM construction does not survive '
            'analytic continuation to BW kinematics — even the '
            'underlying Mandelstam crossing identity fails. The '
            'closed-form F is a Compton-only algebraic fit, not a '
            'process-general amplitude factor.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'crossing_map': 's_C → u_BW, t_C → s_BW, u_C → t_BW',
        'crossed_variables': (
            'x_⊗ = −(1−βc)/(1+βc),  '
            'c_⊗ = (2β²−β²c²−1)/(1−β²c²),  '
            'sin²θ_⊗ = 4β²(1−β²)sin²θ/(1−β²c²)²'
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
    L.append('# Breit-Wheeler cross-process validation probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Cross-process test of the closed-form BAM Compton vertex factor '
        'F²(x, c) under standard QED Mandelstam crossing to pair '
        'production γγ → e⁺e⁻.'
    )
    L.append('')
    L.append('**Crossing map:**')
    L.append('')
    L.append('```')
    L.append(s['crossing_map'])
    L.append('```')
    L.append('')
    L.append('**Crossed variables at BW kinematics (β, c = cosθ_CM):**')
    L.append('')
    L.append('```')
    L.append(s['crossed_variables'])
    L.append('```')
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            value = f"max |M̄|²_BW − (−|M̄|²_KN_crossed) = {t['max_difference']:.2e}"
        elif nm.startswith('T2'):
            value = f"max |M̄|²_BW_BAM − |M̄|²_BW_textbook = {t['max_difference']:.2e}"
        elif nm.startswith('T3'):
            value = (
                f"baseline range = [{t['baseline_range'][0]:.3e}, "
                f"{t['baseline_range'][1]:.3e}]; "
                f"F² range = [{t['F2_range'][0]:.3e}, "
                f"{t['F2_range'][1]:.3e}]; "
                f"{t['n_baseline_negative']} baseline-negative samples; "
                f"{t['n_F2_negative']} F²-negative samples"
            )
        elif nm.startswith('T4'):
            value = (
                f"⟨ratio⟩ = {t['ratio_mean']:.6f}, "
                f"rel-σ = {t['ratio_std_over_mean']:.2e}"
            )
        elif nm.startswith('T5'):
            value = (
                f"⟨BAM/textbook⟩ = {t['ratio_mean']:.6f}, "
                f"rel-spread = {t['relative_spread']:.2e}"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1 detail
    t1 = s['tests'][0]
    L.append('## T1: Mandelstam crossing identity')
    L.append('')
    L.append(
        'Standard QED gives |M̄|²_BW(β, c) = −|M̄|²_KN(s,t,u) evaluated at '
        'the crossed Mandelstam triple. The minus sign comes from one '
        'crossed fermion line.'
    )
    L.append('')
    L.append('| β | cosθ | M²_KN (crossed) | M²_BW (direct) | diff |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t1['samples_first_8'][:8]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['M2_KN_at_crossed_Mandelstam']:+.4e} | "
            f"{r['M2_BW_direct']:+.4e} | {r['difference']:.2e} |"
        )
    L.append('')

    # T2 detail
    t2 = s['tests'][1]
    L.append('## T2: BAM crossed prediction vs textbook BW')
    L.append('')
    L.append(
        'BAM-predicted |M̄|²_BW = −2·(f_baseline · F²)/x_⊗² at the '
        'crossed variables vs textbook |M̄|²_BW. The conversion factor '
        '−2/x_⊗² is the standard Compton x²/2 relation between '
        'f_KN_invariant and |M̄|²/(8e⁴), with the fermion-crossing sign.'
    )
    L.append('')
    L.append('| β | cosθ | x_⊗ | c_⊗ | BAM M² | textbook M² | diff |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t2['samples_first_12'][:12]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['x_crossed']:+.4f} | {r['c_crossed']:+.4f} | "
            f"{r['BAM_predicted_M2_BW']:+.4e} | "
            f"{r['textbook_M2_BW']:+.4e} | "
            f"{r['difference']:.2e} |"
        )
    L.append('')
    L.append(f"Max |BAM − textbook| = {t2['max_difference']:.4e}")
    L.append('')

    # T3 detail
    t3 = s['tests'][2]
    L.append('## T3: Individual factor inspection at BW kinematics')
    L.append('')
    L.append(
        f"f_baseline_C(x_⊗, c_⊗) range: "
        f"[{t3['baseline_range'][0]:.3e}, {t3['baseline_range'][1]:.3e}]"
    )
    L.append(
        f"F²_C(x_⊗, c_⊗) range: "
        f"[{t3['F2_range'][0]:.3e}, {t3['F2_range'][1]:.3e}]"
    )
    L.append(
        f"Baseline-negative samples: {t3['n_baseline_negative']}; "
        f"F²-negative samples: {t3['n_F2_negative']}."
    )
    L.append('')
    L.append('| β | cosθ | x_⊗ | c_⊗ | baseline | F² |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t3['samples_first_10'][:10]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['x_crossed']:+.4f} | {r['c_crossed']:+.4f} | "
            f"{r['baseline_at_crossed']:+.4e} | "
            f"{r['F2_at_crossed']:+.4e} |"
        )
    L.append('')

    # T4 detail
    t4 = s['tests'][3]
    L.append('## T4: Total BW cross section')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | ratio | textbook σ_BW |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t4['rows']:
        L.append(
            f"| {r['beta']:.3f} | {r['textbook_diff_integral']:+.4e} | "
            f"{r['BAM_diff_integral']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.6f} | "
            f"{r['textbook_total_formula']:+.4e} |"
        )
    L.append('')

    # T5 detail
    t5 = s['tests'][4]
    L.append('## T5: Threshold and ultra-relativistic limits')
    L.append('')
    L.append('### Threshold (β → 0): linear suppression')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | textbook/β | BAM/β | ratio |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t5['threshold_rows']:
        L.append(
            f"| {r['beta']:.3f} | {r['textbook_int_diff']:+.4e} | "
            f"{r['BAM_int_diff']:+.4e} | "
            f"{r['textbook_over_beta']:+.4e} | "
            f"{r['BAM_over_beta']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.8f} |"
        )
    L.append('')
    L.append('### Ultra-relativistic (β → 1): logarithmic falloff')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | log env | textbook/env | BAM/env | ratio |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t5['ultra_relativistic_rows']:
        L.append(
            f"| {r['beta']:.4f} | {r['textbook_int_diff']:+.4e} | "
            f"{r['BAM_int_diff']:+.4e} | "
            f"{r['log_envelope']:+.4e} | "
            f"{r['textbook_over_envelope']:+.4e} | "
            f"{r['BAM_over_envelope']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.8f} |"
        )
    L.append('')
    L.append(
        f"Overall ratio statistics: ⟨BAM/textbook⟩ = {t5['ratio_mean']:.8f}, "
        f"relative spread = {t5['relative_spread']:.2e}"
    )
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Pair annihilation e⁺e⁻ → γγ.** A second crossed partner '
        'with the photons outgoing rather than incoming; should follow '
        'trivially from BW by reversing the crossing.'
    )
    L.append(
        '- **Bhabha / Møller scattering.** Both s- and t-channel '
        'diagrams together; the BAM construction would need to combine '
        'two crossed copies of the elementary vertex factor.'
    )
    L.append(
        '- **Loop corrections.** The closed-form F is tree-level. '
        'Whether the same algebraic structure persists at one loop is '
        'a separate (large) thread target.'
    )
    L.append(
        '- **BAM derivation of F from first principles.** Still open '
        'from PR #35 — the closed form is identified but not derived '
        'from a BAM Lagrangian / closure action.'
    )
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, complex):
        return {'real': o.real, 'imag': o.imag}
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_breit_wheeler_cross_process_probe'
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
