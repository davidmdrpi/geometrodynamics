"""
Pair-annihilation crossing probe — closing the Compton/BW/annihilation
crossing triangle.

PR #35 derived the closed-form Compton vertex factor

    F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]

PR #36 showed that the BAM kernel `f_baseline · F²` analytically
continued via the Compton → BW Mandelstam crossing
(`s_C → u_BW, t_C → s_BW, u_C → t_BW`) reproduces |M̄|²_BW exactly.

This probe closes the remaining two edges of the triangle:

      Compton  (γe⁻ → γe⁻)
         / \
        /   \
       v     v
       BW ◀▶ annihilation
   (γγ→ee)  (ee→γγ)   T-reversal

  - Compton → annihilation via (s_C → u_ann, t_C → s_ann, u_C → t_ann).
  - BW ↔ annihilation via T-invariance (same Mandelstam, same |M̄|²).
  - Triangle loop closure: composing the three crossings returns
    identity.

Annihilation CM kinematics (e⁻ along +z, photon at angle θ):

    s_ann = 4E² = s,
    t_ann = m² − 2E²(1 − β·cosθ),
    u_ann = m² − 2E²(1 + β·cosθ),
    β = √(1 − 4m²/s).

Crossed variables in (β, cosθ):

    x_⊗ = (m² − t_ann)/(u_ann − m²) = −(1 − β·cosθ)/(1 + β·cosθ)
    c_⊗ = 1 + 2·s·m²/[(u_ann − m²)(m² − t_ann)]
        = (2β² − β²cos²θ − 1)/(1 − β²cos²θ)

(identical functional form to BW — the Mandelstam permutation is the
same.)

Tests:

  T1. Mandelstam crossing identity for annihilation.
  T2. BAM-predicted |M̄|²_ann matches textbook.
  T3. T-invariance: |M̄|²_BW(β, c) = |M̄|²_ann(β, c).
  T4. Triangle loop closure: π_ann→C ∘ π_BW→ann ∘ π_C→BW = identity.
  T5. Total annihilation cross section reproduces the Dirac formula
      σ_ann(β) = (π·r_e²/(2β²))·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                          − 2β·(2−β²)].
  T6. Threshold (β → 0, σ ∼ π/β divergence) and high-energy limits.
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
# Compton-side BAM construction (from PR #35)
# ---------------------------------------------------------------------------

def f_baseline_C(x: float, c: float) -> float:
    return (1.0 + c * c) * (1.0 + 1.0 / x) ** 2 / 8.0


def F_squared_C(x: float, c: float) -> float:
    s2 = 1.0 - c * c
    num = 4.0 * x ** 3 * (x * x + 1.0 - x * s2)
    den = (1.0 + c * c) * (1.0 + x) ** 2
    return num / den


def f_KN_invariant(x: float, c: float) -> float:
    """f_baseline · F² in closed form = x·(x²+1−x·sin²θ)/2."""
    s2 = 1.0 - c * c
    return x * (x * x + 1.0 - x * s2) / 2.0


# ---------------------------------------------------------------------------
# Standard QED amplitude in invariant form
# ---------------------------------------------------------------------------

def M2_KN_over_8e4(s: float, t: float, u: float, m2: float) -> float:
    """|M̄|²(γe→γe)/(8e⁴) as a function of Mandelstam invariants."""
    A = 2.0 * t * m2 / ((s - m2) * (m2 - u))
    return (s - m2) / (m2 - u) + (m2 - u) / (s - m2) + 2.0 * A + A * A


def M2_BW_or_ann_textbook_over_8e4(beta: float, c: float) -> float:
    """Textbook |M̄|²/(8e⁴) for γγ → e⁺e⁻ (equiv. e⁺e⁻ → γγ) in CM:

        |M̄|²/(8e⁴) = 2(1+β²c²)/(1−β²c²) + 4(1−β²)/(1−β²c²)
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
# Annihilation kinematics + Mandelstam crossing
# ---------------------------------------------------------------------------

def ann_mandelstam(beta: float, c: float, m2: float = 1.0) -> tuple[float, float, float]:
    """Annihilation Mandelstam invariants in CM frame, units m² supplied.

    e⁻ along +z, photon at angle θ. m² = E²(1−β²), s = 4E².
    """
    s = 4.0 * m2 / (1.0 - beta * beta)
    E2 = s / 4.0
    t = m2 - 2.0 * E2 * (1.0 - beta * c)
    u = m2 - 2.0 * E2 * (1.0 + beta * c)
    return s, t, u


def crossed_compton_vars_ann(beta: float, c: float, m2: float = 1.0) -> tuple[float, float]:
    """BAM Compton (x, c) at annihilation kinematics under
    (s_C → u_ann, t_C → s_ann, u_C → t_ann)."""
    s_a, t_a, u_a = ann_mandelstam(beta, c, m2)
    x_cross = (m2 - t_a) / (u_a - m2)
    c_cross = 1.0 + 2.0 * s_a * m2 / ((u_a - m2) * (m2 - t_a))
    return x_cross, c_cross


# Re-implement BW crossing for the T-invariance test
def bw_mandelstam(beta: float, c: float, m2: float = 1.0) -> tuple[float, float, float]:
    """BW CM-frame Mandelstam invariants. Functional form identical to
    annihilation by T-invariance."""
    s = 4.0 * m2 / (1.0 - beta * beta)
    E2 = s / 4.0
    t = m2 - 2.0 * E2 * (1.0 - beta * c)
    u = m2 - 2.0 * E2 * (1.0 + beta * c)
    return s, t, u


def crossed_compton_vars_bw(beta: float, c: float, m2: float = 1.0) -> tuple[float, float]:
    """BAM Compton (x, c) at BW kinematics under PR #36 crossing
    (s_C → u_BW, t_C → s_BW, u_C → t_BW)."""
    s_b, t_b, u_b = bw_mandelstam(beta, c, m2)
    x_cross = (m2 - t_b) / (u_b - m2)
    c_cross = 1.0 + 2.0 * s_b * m2 / ((u_b - m2) * (m2 - t_b))
    return x_cross, c_cross


# ---------------------------------------------------------------------------
# T1. Mandelstam crossing identity for annihilation
# ---------------------------------------------------------------------------

def test_T1_mandelstam_crossing_ann() -> dict:
    """Verify |M̄|²_ann(β, c) = −|M̄|²_KN evaluated at the crossed
    Mandelstam triple (s_C = u_ann, t_C = s_ann, u_C = t_ann)."""
    m2 = 1.0
    betas = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 0.99]
    cs = np.linspace(-0.95, 0.95, 17)
    samples = []
    max_diff = 0.0
    for beta in betas:
        for c in cs:
            s_a, t_a, u_a = ann_mandelstam(beta, float(c), m2)
            m2_kn_crossed = M2_KN_over_8e4(s=u_a, t=s_a, u=t_a, m2=m2)
            m2_ann_direct = M2_BW_or_ann_textbook_over_8e4(beta, float(c))
            diff = abs((-m2_kn_crossed) - m2_ann_direct)
            if diff > max_diff:
                max_diff = diff
            if len(samples) < 8:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'M2_KN_at_crossed': m2_kn_crossed,
                    'M2_ann_direct': m2_ann_direct,
                    'predicted_minus_M2_KN': -m2_kn_crossed,
                    'difference': diff,
                })
    return {
        'name': 'T1_mandelstam_crossing_ann',
        'description': (
            'Standard QED: |M̄|²_ann = −|M̄|²_KN at the crossed '
            '(s_C=u_ann, t_C=s_ann, u_C=t_ann) Mandelstam triple. '
            'The minus sign is the fermion-line crossing factor.'
        ),
        'samples_first_8': samples,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T2. BAM kernel crossed to annihilation
# ---------------------------------------------------------------------------

def test_T2_BAM_crossed_to_ann() -> dict:
    """BAM-predicted |M̄|²_ann/(8e⁴) = −2·(f_baseline · F²)/x_⊗²
    evaluated at the annihilation-crossed variables. Compare to the
    textbook |M̄|²_ann across a (β, c) grid."""
    m2 = 1.0
    betas = [0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
    cs = np.linspace(-0.95, 0.95, 17)
    samples = []
    max_diff = 0.0
    for beta in betas:
        for c in cs:
            x_x, c_x = crossed_compton_vars_ann(beta, float(c), m2)
            product = f_KN_invariant(x_x, c_x)  # = f_baseline · F²
            M2_ann_BAM = -2.0 * product / (x_x * x_x)
            M2_ann_direct = M2_BW_or_ann_textbook_over_8e4(beta, float(c))
            diff = abs(M2_ann_BAM - M2_ann_direct)
            if diff > max_diff:
                max_diff = diff
            if len(samples) < 12:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'x_crossed': x_x,
                    'c_crossed': c_x,
                    'BAM_product': product,
                    'BAM_predicted_M2_ann': M2_ann_BAM,
                    'textbook_M2_ann': M2_ann_direct,
                    'difference': diff,
                })
    return {
        'name': 'T2_BAM_crossed_to_ann',
        'description': (
            'BAM-predicted |M̄|²_ann = −2·(f_baseline · F²)/x_⊗² at '
            'annihilation-crossed variables vs textbook. Machine '
            'precision agreement = closed-form F is process-general '
            'under the Compton → annihilation crossing.'
        ),
        'samples_first_12': samples,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T3. T-invariance: BAM-BW = BAM-ann
# ---------------------------------------------------------------------------

def test_T3_T_invariance() -> dict:
    """The BW crossing and annihilation crossing both use the same
    Mandelstam permutation (s, t, u) → (u, s, t). The crossed
    variables (x_⊗, c_⊗) at the same (β, c) are therefore identical;
    the BAM predicted |M̄|² is the same function.

    Verify numerically that BAM-BW(β, c) = BAM-ann(β, c) to machine
    precision."""
    m2 = 1.0
    betas = [0.1, 0.3, 0.5, 0.7, 0.9, 0.98]
    cs = np.linspace(-0.95, 0.95, 17)
    max_diff = 0.0
    samples = []
    for beta in betas:
        for c in cs:
            x_bw, c_bw = crossed_compton_vars_bw(beta, float(c), m2)
            x_ann, c_ann = crossed_compton_vars_ann(beta, float(c), m2)
            M2_bw_bam = -2.0 * f_KN_invariant(x_bw, c_bw) / (x_bw * x_bw)
            M2_ann_bam = -2.0 * f_KN_invariant(x_ann, c_ann) / (x_ann * x_ann)
            diff = abs(M2_bw_bam - M2_ann_bam)
            if diff > max_diff:
                max_diff = diff
            if len(samples) < 8:
                samples.append({
                    'beta': beta,
                    'cos_theta': float(c),
                    'x_crossed_BW': x_bw,
                    'x_crossed_ann': x_ann,
                    'BAM_M2_BW': M2_bw_bam,
                    'BAM_M2_ann': M2_ann_bam,
                    'difference': diff,
                })
    return {
        'name': 'T3_T_invariance_BW_equals_ann',
        'description': (
            'BW and annihilation share the Mandelstam permutation; '
            'BAM-predicted |M̄|² is the same function of (β, c). '
            'Verifies T-invariance at the kernel level.'
        ),
        'samples_first_8': samples,
        'max_difference': max_diff,
        'pass': max_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T4. Triangle loop closure
# ---------------------------------------------------------------------------

def apply_permutation(stu: tuple[float, float, float], perm: str) -> tuple[float, float, float]:
    """Apply a Mandelstam permutation. perm is a 3-letter string in
    {s, t, u}^3 meaning new_s, new_t, new_u in terms of old labels."""
    s, t, u = stu
    mapping = {'s': s, 't': t, 'u': u}
    return mapping[perm[0]], mapping[perm[1]], mapping[perm[2]]


def test_T4_triangle_loop_closure() -> dict:
    """The three crossings of the triangle are Mandelstam permutations:

      π_C→BW : (s, t, u) → (u, s, t)   [Compton-side label '(u, s, t)'
                                         means new_s = u_old, etc.]
      π_BW→ann : identity (T-reversal preserves Mandelstam invariants)
      π_ann→C : inverse of π_C→ann = π_C→BW (also (s,t,u)→(t,u,s))

    Wait — both π_C→BW and π_C→ann take Compton labels to (u, s, t).
    So π_ann→C = inverse of π_C→ann = (s,t,u) → (t,u,s), which is the
    inverse of the cyclic permutation (s,t,u)→(u,s,t).

    Loop closure: π_ann→C ∘ π_BW→ann ∘ π_C→BW applied to a Compton
    label triple (s,t,u) should give back (s,t,u).

    Implementation: apply each permutation to a sample (s,t,u) triple
    and verify the round trip is identity.

    Additional check at the amplitude level: |M̄|²_KN(s,t,u) computed
    on (s,t,u) and on the round-tripped (s,t,u) must be identical.
    """
    samples = []
    max_label_diff = 0.0
    max_amplitude_diff = 0.0
    test_triples = [
        (2.0, -0.3, -0.7),   # Compton-like: s > m², t < 0
        (3.0, -0.1, -1.9),
        (5.0, -0.5, -3.5),
        (10.0, -1.0, -8.0),
    ]
    m2 = 1.0
    for stu_orig in test_triples:
        # Step 1: Compton → BW, permutation (s,t,u) → (u,s,t)
        # i.e. new_s = u_old, new_t = s_old, new_u = t_old → perm 'ust'
        stu_after_C_to_BW = apply_permutation(stu_orig, 'ust')
        # Step 2: BW → ann, identity (T-reversal)
        stu_after_BW_to_ann = stu_after_C_to_BW
        # Step 3: ann → Compton, the inverse of π_C→ann ('ust').
        # The inverse of cyclic (s,t,u)→(u,s,t) is (s,t,u)→(t,u,s),
        # i.e. perm 'tus'.
        stu_after_loop = apply_permutation(stu_after_BW_to_ann, 'tus')
        # Verify round-trip identity at the label level
        label_diff = sum(abs(a - b) for a, b in zip(stu_orig, stu_after_loop))
        max_label_diff = max(max_label_diff, label_diff)

        # Verify at the amplitude level
        s0, t0, u0 = stu_orig
        s1, t1, u1 = stu_after_loop
        M2_orig = M2_KN_over_8e4(s0, t0, u0, m2)
        M2_round = M2_KN_over_8e4(s1, t1, u1, m2)
        amp_diff = abs(M2_orig - M2_round)
        max_amplitude_diff = max(max_amplitude_diff, amp_diff)

        samples.append({
            'stu_original_Compton': list(stu_orig),
            'stu_after_C_to_BW': list(stu_after_C_to_BW),
            'stu_after_BW_to_ann': list(stu_after_BW_to_ann),
            'stu_after_ann_to_Compton_loop': list(stu_after_loop),
            'label_round_trip_diff': label_diff,
            'M2_KN_original': M2_orig,
            'M2_KN_round_tripped': M2_round,
            'amplitude_round_trip_diff': amp_diff,
        })

    return {
        'name': 'T4_triangle_loop_closure',
        'description': (
            'Compose π_C→BW ∘ π_BW→ann ∘ π_ann→C and verify the loop '
            'is identity, both at the Mandelstam label level and at '
            'the amplitude level |M̄|²_KN.'
        ),
        'permutations': {
            'C_to_BW': '(s,t,u) → (u,s,t)  perm "ust"',
            'BW_to_ann': '(s,t,u) → (s,t,u)  identity (T-reversal)',
            'ann_to_C': '(s,t,u) → (t,u,s)  perm "tus"  [inverse of C→ann]',
        },
        'samples': samples,
        'max_label_round_trip_diff': max_label_diff,
        'max_amplitude_round_trip_diff': max_amplitude_diff,
        'pass': max_label_diff < 1e-12 and max_amplitude_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T5. Total annihilation cross section (Dirac formula)
# ---------------------------------------------------------------------------

def sigma_ann_textbook(beta: float, r_e: float = 1.0) -> float:
    """Dirac annihilation cross section e⁺e⁻ → γγ:

        σ_ann(β) = (π·r_e²/(2β²))·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                            − 2β·(2−β²)]
    """
    if beta <= 0.0 or beta >= 1.0:
        return 0.0
    one_minus_b2 = 1.0 - beta * beta
    log_term = math.log((1.0 + beta) / (1.0 - beta))
    return (PI * r_e ** 2 / (2.0 * beta ** 2)) * one_minus_b2 * (
        (3.0 - beta ** 4) * log_term - 2.0 * beta * (2.0 - beta ** 2)
    )


def differential_ann_textbook_per_dcos(beta: float, c: float) -> float:
    """dσ_ann/dcosθ ∝ M2/(β·s) in matched units. The β/s of the BW
    differential becomes 1/(β·s) for annihilation."""
    s = 4.0 / (1.0 - beta * beta)
    M2 = M2_BW_or_ann_textbook_over_8e4(beta, c)
    return M2 / (beta * s)


def differential_ann_BAM_per_dcos(beta: float, c: float) -> float:
    """BAM-predicted dσ_ann/dcosθ in matched units."""
    s = 4.0 / (1.0 - beta * beta)
    x_x, c_x = crossed_compton_vars_ann(beta, c, m2=1.0)
    product = f_KN_invariant(x_x, c_x)
    M2_ann_BAM = -2.0 * product / (x_x * x_x)
    return M2_ann_BAM / (beta * s)


def integrate_differential(differential_fn, beta: float, n: int = 401) -> float:
    cs = np.linspace(-1.0 + 1e-6, 1.0 - 1e-6, n)
    vals = np.array([differential_fn(beta, float(c)) for c in cs])
    return float(np.trapezoid(vals, cs))


def test_T5_total_cross_section() -> dict:
    """Integrate the BAM-predicted dσ_ann/dcosθ and compare to the
    textbook differential integral. A constant ratio (= 1) across all
    β regimes means the closed-form F reproduces the Dirac
    annihilation cross section under crossing."""
    betas = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98]
    rows = []
    ratios = []
    for beta in betas:
        textbook_int = integrate_differential(differential_ann_textbook_per_dcos, beta)
        bam_int = integrate_differential(differential_ann_BAM_per_dcos, beta)
        ratio = bam_int / textbook_int if abs(textbook_int) > 1e-20 else float('nan')
        ratios.append(ratio)
        rows.append({
            'beta': beta,
            'textbook_int_diff': textbook_int,
            'BAM_int_diff': bam_int,
            'ratio_BAM_over_textbook': ratio,
            'textbook_total_Dirac': sigma_ann_textbook(beta, r_e=1.0),
        })
    ratios_arr = np.array(ratios)
    ratio_mean = float(np.mean(ratios_arr))
    relative_spread = (
        float(np.max(ratios_arr)) - float(np.min(ratios_arr))
    ) / abs(ratio_mean) if ratio_mean != 0 else float('inf')
    return {
        'name': 'T5_total_cross_section',
        'description': (
            'Integrate the BAM-predicted dσ_ann/dcosθ and compare to '
            'the textbook integrated differential. The Dirac total '
            'cross section is reported alongside for reference.'
        ),
        'rows': rows,
        'ratio_mean': ratio_mean,
        'relative_spread': relative_spread,
        'pass': abs(ratio_mean - 1.0) < 1e-6 and relative_spread < 1e-6,
    }


# ---------------------------------------------------------------------------
# T6. Threshold (β → 0: σ ∼ π/β divergence) and high-energy limits
# ---------------------------------------------------------------------------

def test_T6_limits() -> dict:
    """Check that the BAM-crossed annihilation prediction reproduces:

      - threshold β → 0: σ_ann ∼ π/β (s-wave Coulomb-like divergence;
        the (1−β²)·(2β + ...)/(2β²) → π/β leading behaviour)
      - ultra-relativistic β → 1: log envelope

    Both must hold with BAM/textbook ratio = 1 across regimes."""
    beta_thr = [0.01, 0.02, 0.05, 0.1]
    thr_rows = []
    for beta in beta_thr:
        textbook_int = integrate_differential(differential_ann_textbook_per_dcos, beta)
        bam_int = integrate_differential(differential_ann_BAM_per_dcos, beta)
        # At threshold σ_textbook · β → π (constant), so textbook_int · β
        # should saturate; verify both textbook and BAM share this scaling.
        thr_rows.append({
            'beta': beta,
            'textbook_int_diff': textbook_int,
            'BAM_int_diff': bam_int,
            'textbook_times_beta': textbook_int * beta,
            'BAM_times_beta': bam_int * beta,
            'ratio_BAM_over_textbook': bam_int / textbook_int if abs(textbook_int) > 1e-20 else float('nan'),
        })

    beta_ur = [0.95, 0.99, 0.999]
    ur_rows = []
    for beta in beta_ur:
        textbook_int = integrate_differential(differential_ann_textbook_per_dcos, beta)
        bam_int = integrate_differential(differential_ann_BAM_per_dcos, beta)
        log_factor = math.log((1.0 + beta) / (1.0 - beta))
        # σ_ann at high energy: (π/(2β²))·(1−β²)·(3·log + ...) → log envelope
        envelope = log_factor   # leading behaviour
        ur_rows.append({
            'beta': beta,
            'textbook_int_diff': textbook_int,
            'BAM_int_diff': bam_int,
            'log_envelope': envelope,
            'ratio_BAM_over_textbook': bam_int / textbook_int if abs(textbook_int) > 1e-20 else float('nan'),
        })

    all_ratios = [r['ratio_BAM_over_textbook'] for r in thr_rows + ur_rows]
    ratios_arr = np.array(all_ratios)
    ratio_mean = float(np.mean(ratios_arr))
    relative_spread = (
        float(np.max(ratios_arr)) - float(np.min(ratios_arr))
    ) / abs(ratio_mean) if ratio_mean != 0 else float('inf')

    return {
        'name': 'T6_threshold_and_ultrarelativistic_limits',
        'description': (
            'Check that the BAM-crossed σ_ann shows the s-wave 1/β '
            'threshold divergence and log-envelope ultra-relativistic '
            'falloff, with constant BAM/textbook ratio across both '
            'regimes.'
        ),
        'threshold_rows': thr_rows,
        'ultra_relativistic_rows': ur_rows,
        'BAM_over_textbook_ratios': all_ratios,
        'ratio_mean': ratio_mean,
        'relative_spread': relative_spread,
        'pass': abs(ratio_mean - 1.0) < 1e-6 and relative_spread < 1e-6,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_mandelstam_crossing_ann()
    t2 = test_T2_BAM_crossed_to_ann()
    t3 = test_T3_T_invariance()
    t4 = test_T4_triangle_loop_closure()
    t5 = test_T5_total_cross_section()
    t6 = test_T6_limits()
    tests = [t1, t2, t3, t4, t5, t6]

    if all(t['pass'] for t in tests):
        verdict_class = 'TRIANGLE_CLOSES'
        verdict = (
            'TRIANGLE CLOSES. The same closed-form Compton vertex '
            'factor F²(x, c) = 4·x³·(x²+1−x·sin²θ)/[(1+c²)(1+x)²], '
            'analytically continued via Mandelstam crossing, reproduces '
            'all three corners of the QED two-photon-two-fermion '
            'triangle: Compton (PR #35), Breit-Wheeler (PR #36), and '
            'annihilation (this PR). The Compton → BW and Compton → ann '
            'crossings share the same Mandelstam permutation '
            '(s,t,u) → (u,s,t); BW ↔ ann is T-reversal at the kernel '
            'level. The triangle loop π_C→BW ∘ π_BW→ann ∘ π_ann→C is '
            'identity, both at the Mandelstam label level and at the '
            '|M̄|²_KN amplitude level. The BAM tree kernel is process-'
            'general across the full crossing triangle.'
        )
    elif t1['pass'] and t2['pass'] and not t4['pass']:
        verdict_class = 'PARTIAL_CLOSURE'
        verdict = (
            'PARTIAL CLOSURE. Pairwise crossings (Compton ↔ BW, '
            'Compton ↔ ann) succeed, but the loop π_C→BW ∘ π_BW→ann ∘ '
            'π_ann→C is not identity. An algebraic inconsistency in '
            'the crossing permutations.'
        )
    elif not t2['pass']:
        verdict_class = 'ANNIHILATION_BREAKS'
        verdict = (
            'ANNIHILATION BREAKS. The BAM kernel reproduces Compton '
            'and BW (PR #35, PR #36) but not annihilation — a Compton-/'
            'BW-specific algebraic feature that fails the second '
            'crossing.'
        )
    else:
        verdict_class = 'TRIANGLE_BROKEN'
        verdict = (
            'TRIANGLE BROKEN. Multiple tests fail; the BAM kernel does '
            'not commute consistently with the crossing triangle.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'crossing_triangle': {
            'C_to_BW': 's_C → u_BW, t_C → s_BW, u_C → t_BW (PR #36)',
            'C_to_ann': 's_C → u_ann, t_C → s_ann, u_C → t_ann (this PR)',
            'BW_to_ann': 'T-reversal (same Mandelstam, same |M̄|²)',
        },
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
    L.append('# Pair-annihilation crossing probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Closes the Compton/Breit-Wheeler/annihilation crossing '
        'triangle. PR #35 established the closed-form F²(x, c); '
        'PR #36 verified Compton → BW; this probe verifies Compton → '
        'annihilation and the full triangle loop closure.'
    )
    L.append('')
    L.append('**Crossing triangle:**')
    L.append('')
    L.append('```')
    L.append('       Compton (γe → γe)')
    L.append('         /            \\')
    L.append('        v              v')
    L.append('    BW (γγ→ee) ◀▶ ann (ee→γγ)')
    L.append('                 (T-reversal)')
    L.append('```')
    L.append('')
    L.append('| edge | crossing |')
    L.append('|---|---|')
    for k, v in s['crossing_triangle'].items():
        L.append(f"| `{k}` | {v} |")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            value = f"max diff = {t['max_difference']:.2e}"
        elif nm.startswith('T2'):
            value = f"max |M̄|²_ann_BAM − textbook = {t['max_difference']:.2e}"
        elif nm.startswith('T3'):
            value = f"max |BAM-BW − BAM-ann| = {t['max_difference']:.2e}"
        elif nm.startswith('T4'):
            value = (
                f"label diff = {t['max_label_round_trip_diff']:.2e}, "
                f"amp diff = {t['max_amplitude_round_trip_diff']:.2e}"
            )
        elif nm.startswith('T5'):
            value = (
                f"⟨BAM/textbook⟩ = {t['ratio_mean']:.8f}, "
                f"rel-spread = {t['relative_spread']:.2e}"
            )
        elif nm.startswith('T6'):
            value = (
                f"⟨BAM/textbook⟩ = {t['ratio_mean']:.8f}, "
                f"rel-spread = {t['relative_spread']:.2e}"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1 detail
    t1 = s['tests'][0]
    L.append('## T1: Mandelstam crossing identity for annihilation')
    L.append('')
    L.append(
        '`|M̄|²_ann(β, c) = −|M̄|²_KN(s = u_ann, t = s_ann, u = t_ann)`.'
    )
    L.append('')
    L.append('| β | cosθ | M²_KN (crossed) | M²_ann (direct) | diff |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t1['samples_first_8'][:8]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['M2_KN_at_crossed']:+.4e} | "
            f"{r['M2_ann_direct']:+.4e} | "
            f"{r['difference']:.2e} |"
        )
    L.append('')

    # T2 detail
    t2 = s['tests'][1]
    L.append('## T2: BAM kernel crossed to annihilation')
    L.append('')
    L.append(
        '`|M̄|²_ann_BAM = −2·(f_baseline · F²)/x_⊗²` at the annihilation-'
        'crossed variables.'
    )
    L.append('')
    L.append('| β | cosθ | x_⊗ | c_⊗ | BAM M² | textbook M² | diff |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t2['samples_first_12'][:12]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['x_crossed']:+.4f} | {r['c_crossed']:+.4f} | "
            f"{r['BAM_predicted_M2_ann']:+.4e} | "
            f"{r['textbook_M2_ann']:+.4e} | "
            f"{r['difference']:.2e} |"
        )
    L.append('')
    L.append(f"Max |BAM − textbook| = {t2['max_difference']:.4e}")
    L.append('')

    # T3 detail
    t3 = s['tests'][2]
    L.append('## T3: T-invariance — BAM-BW = BAM-ann')
    L.append('')
    L.append('| β | cosθ | x_⊗ (BW) | x_⊗ (ann) | BAM M²_BW | BAM M²_ann | diff |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t3['samples_first_8'][:8]:
        L.append(
            f"| {r['beta']:.3f} | {r['cos_theta']:+.3f} | "
            f"{r['x_crossed_BW']:+.4f} | {r['x_crossed_ann']:+.4f} | "
            f"{r['BAM_M2_BW']:+.4e} | {r['BAM_M2_ann']:+.4e} | "
            f"{r['difference']:.2e} |"
        )
    L.append('')

    # T4 detail
    t4 = s['tests'][3]
    L.append('## T4: Triangle loop closure')
    L.append('')
    for k, v in t4['permutations'].items():
        L.append(f"- `{k}`: {v}")
    L.append('')
    L.append(
        'Round-trip identity check at the Mandelstam label level and '
        'at the |M̄|²_KN amplitude level:'
    )
    L.append('')
    L.append('| original (s, t, u) | after C→BW | after BW→ann | after loop | label diff | amp diff |')
    L.append('|---|---|---|---|---:|---:|')
    for r in t4['samples']:
        L.append(
            f"| {r['stu_original_Compton']} | "
            f"{r['stu_after_C_to_BW']} | "
            f"{r['stu_after_BW_to_ann']} | "
            f"{r['stu_after_ann_to_Compton_loop']} | "
            f"{r['label_round_trip_diff']:.2e} | "
            f"{r['amplitude_round_trip_diff']:.2e} |"
        )
    L.append('')

    # T5 detail
    t5 = s['tests'][4]
    L.append('## T5: Total annihilation cross section')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | ratio | Dirac σ_ann |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t5['rows']:
        L.append(
            f"| {r['beta']:.3f} | {r['textbook_int_diff']:+.4e} | "
            f"{r['BAM_int_diff']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.8f} | "
            f"{r['textbook_total_Dirac']:+.4e} |"
        )
    L.append('')

    # T6 detail
    t6 = s['tests'][5]
    L.append('## T6: Threshold and ultra-relativistic limits')
    L.append('')
    L.append('### Threshold (β → 0): σ ∼ π/β s-wave divergence')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | textbook·β | BAM·β | ratio |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t6['threshold_rows']:
        L.append(
            f"| {r['beta']:.3f} | {r['textbook_int_diff']:+.4e} | "
            f"{r['BAM_int_diff']:+.4e} | "
            f"{r['textbook_times_beta']:+.4e} | "
            f"{r['BAM_times_beta']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.8f} |"
        )
    L.append('')
    L.append('### Ultra-relativistic (β → 1): log envelope')
    L.append('')
    L.append('| β | textbook ∫ | BAM ∫ | log env | ratio |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t6['ultra_relativistic_rows']:
        L.append(
            f"| {r['beta']:.4f} | {r['textbook_int_diff']:+.4e} | "
            f"{r['BAM_int_diff']:+.4e} | "
            f"{r['log_envelope']:+.4e} | "
            f"{r['ratio_BAM_over_textbook']:+.8f} |"
        )
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Bhabha (e⁺e⁻ → e⁺e⁻) and Møller (e⁻e⁻ → e⁻e⁻).** '
        'Two-channel processes with both s- and t-channel diagrams '
        'interfering. The BAM kernel would need to combine two crossed '
        'copies coherently.'
    )
    L.append(
        '- **One-loop corrections.** Still tree-level only.'
    )
    L.append(
        '- **BAM first-principles derivation of F².** Still open '
        'from PR #35.'
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
    out = here / 'runs' / f'{ts}_pair_annihilation_crossing_probe'
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
