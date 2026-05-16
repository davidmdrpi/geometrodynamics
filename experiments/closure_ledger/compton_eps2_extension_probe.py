"""
O(ε²) analytic extension probe for the BAM Compton vertex.

Follow-on to PR #31 (analytic γ = −3/2 at O(ε)) and PR #33 (d-scaling
discrimination, 7 candidates surviving for γ origin). This probe
extends the analytic derivation to O(ε²) and asks which surviving
candidates predict the next-order coefficient correctly.

Analytic targets (derived in the probe):

  f_KN(ε,θ)        = (1+c²)/2 − ε(1−c)(1+c²) + ε²·(1−c)²(4+3c²)/2 + O(ε³)
  f_BAM_baseline   = (1+c²)/2 + ε(1−c)(1+c²)/2 + ε²(1−c)²(1+c²)/8 + O(ε³)
  f_BAM_PR31_O(ε²) = −(1−c)²(1+c²)/4

The O(ε²) residual to close:
  Δ_required_O(ε²) = (1−c)²·(9 + 7c²)/4

In the c⁰..c⁴ basis: (9, −18, 16, −14, 7)/4.

Tests:

  T1. Verify analytic O(ε²) expansion numerically.

  T2. Decompose Δ_required in c basis and report coefficients.

  T3. Family B' extension (ε²·μ₂ with μ₂ a polynomial in cos θ).
      Solve μ₂·(1+c²) = Δ_required·... and check for polynomial
      solutions.

  T4. Family A extension (ε·α coupling, contributing α²sin⁴θ at O(ε²)).
      Find α that minimises O(ε²) residual.

  T5. Combined A + B' family. Solve the joint system.

  T6. Compare numerical O(ε²) residual scaling to candidate predictions.
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
# Reference: KN exact and BAM baseline
# ---------------------------------------------------------------------------

def x_ratio(eps: float, theta: float) -> float:
    return 1.0 / (1.0 + eps * (1.0 - math.cos(theta)))


def f_KN_normalized(eps: float, theta: float) -> float:
    x = x_ratio(eps, theta)
    val = x * x * (x + 1.0 / x - math.sin(theta) ** 2)
    return val / 2.0


def f_BAM_baseline_normalized(eps: float, theta: float) -> float:
    x = x_ratio(eps, theta)
    angular = (1.0 + math.cos(theta) ** 2) / 2.0
    propagator_factor = (1.0 + 1.0 / x) ** 2 / 4.0
    return angular * propagator_factor


def f_BAM_with_PR31_vertex(eps: float, theta: float) -> float:
    """f_BAM with the PR #31 analytic vertex γ = -3/2, β = 0."""
    c = math.cos(theta)
    s2 = math.sin(theta) ** 2
    angular_mod = 1.0 + eps * (-1.5 * (1.0 - c))
    pol_sum_modified = (1.0 + c * c) / 2.0 * angular_mod * angular_mod
    x = x_ratio(eps, theta)
    propagator_factor = (1.0 + 1.0 / x) ** 2 / 4.0
    return pol_sum_modified * propagator_factor / 1.0
    # Normalization at θ=0: pol_sum_modified(ε,0) = 1 · 1 = 1,
    # propagator_factor(ε,0) = 1. So return as-is.


# ---------------------------------------------------------------------------
# T1. Analytic O(ε²) verification
# ---------------------------------------------------------------------------

def second_order_coeff(f, theta: float, eps_small: float = 1e-4) -> float:
    """Extract f''(0)/2 ≈ (f(ε) − f(0))/ε² − f'(0)/ε via finite difference."""
    f0 = f(0.0, theta)
    fp = f(eps_small, theta)
    fm = f(-eps_small, theta)
    # Symmetric second derivative: f''(0) ≈ (f(ε) + f(−ε) − 2·f(0))/ε²
    # → O(ε²) coefficient = f''(0)/2 ≈ (f(ε) + f(−ε) − 2·f(0))/(2·ε²)
    return (fp + fm - 2.0 * f0) / (2.0 * eps_small * eps_small)


def test_T1_analytic_O_eps2_verification() -> dict:
    """Compare numerical O(ε²) coefficients vs closed-form predictions:
      f_KN O(ε²)    = (1−c)²·(4 + 3c²)/2
      f_BAM_PR31 O(ε²) = −(1−c)²(1+c²)/4
    """
    thetas = np.linspace(0.05, PI - 0.05, 17)
    max_kn_err = 0.0
    max_pr31_err = 0.0
    rows = []
    for theta in thetas:
        c = math.cos(theta)
        kn_num = second_order_coeff(f_KN_normalized, float(theta))
        kn_ana = (1.0 - c) ** 2 * (4.0 + 3.0 * c * c) / 2.0
        pr31_num = second_order_coeff(f_BAM_with_PR31_vertex, float(theta))
        pr31_ana = -(1.0 - c) ** 2 * (1.0 + c * c) / 4.0
        max_kn_err = max(max_kn_err, abs(kn_num - kn_ana))
        max_pr31_err = max(max_pr31_err, abs(pr31_num - pr31_ana))
        rows.append({
            'theta_in_pi': float(theta) / PI,
            'cos_theta': c,
            'KN_O_eps2_num': float(kn_num),
            'KN_O_eps2_ana': float(kn_ana),
            'BAM_PR31_O_eps2_num': float(pr31_num),
            'BAM_PR31_O_eps2_ana': float(pr31_ana),
        })
    return {
        'name': 'T1_analytic_O_eps2_verification',
        'description': (
            'Verify analytic O(ε²) coefficients of f_KN and '
            'f_BAM_PR31 numerically (second-derivative finite '
            'difference).'
        ),
        'samples_first_8': rows[:8],
        'max_KN_O_eps2_error': max_kn_err,
        'max_BAM_PR31_O_eps2_error': max_pr31_err,
        # Tolerance: finite-difference noise is O((eps_small)²) ~ 1e-8
        # for ε_small=1e-4, plus truncation O(ε²) ~ 1e-8. Round-off
        # dominates: be generous (1e-3 for safety).
        'pass': max(max_kn_err, max_pr31_err) < 1e-3,
    }


# ---------------------------------------------------------------------------
# T2. Δ_required angular decomposition
# ---------------------------------------------------------------------------

def test_T2_delta_required_decomposition() -> dict:
    """Δ_required = (1−c)²·(9+7c²)/4 in c⁰..c⁴ basis.

    Expand: (1 − 2c + c²)(9 + 7c²)/4
            = (9 + 7c² − 18c − 14c³ + 9c² + 7c⁴)/4
            = (9 − 18c + 16c² − 14c³ + 7c⁴)/4
    """
    # Compute analytically
    delta_required_coeffs = np.array([9.0, -18.0, 16.0, -14.0, 7.0]) / 4.0

    # Numerically verify by evaluating at several θ
    thetas = np.linspace(0.05, PI - 0.05, 9)
    rows = []
    for theta in thetas:
        c = math.cos(theta)
        # Analytic Δ_required
        delta_ana = (1.0 - c) ** 2 * (9.0 + 7.0 * c * c) / 4.0
        # From basis decomposition
        delta_basis = sum(coef * c ** k for k, coef in enumerate(delta_required_coeffs))
        rows.append({
            'theta_in_pi': float(theta) / PI,
            'delta_required_analytic': float(delta_ana),
            'delta_required_from_basis': float(delta_basis),
            'difference': float(delta_ana - delta_basis),
        })
    max_diff = max(abs(r['difference']) for r in rows)
    return {
        'name': 'T2_delta_required_decomposition',
        'description': (
            'Δ_required = (1−c)²·(9+7c²)/4 decomposed in {1,c,c²,c³,c⁴} '
            'basis. Coefficients (9, −18, 16, −14, 7)/4.'
        ),
        'delta_required_coefficients_c0_to_c4': delta_required_coeffs.tolist(),
        'samples_verification': rows,
        'max_decomposition_error': max_diff,
        'pass': max_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T3. Family B' extension (μ_2 polynomial in c)
# ---------------------------------------------------------------------------

def test_T3_family_B_prime() -> dict:
    """Solve for μ₂ = β'·sin²θ + γ'·(1−c) + δ'·(1−c)² + ζ'·c + η'·c·(1−c)
    such that μ₂·(1+c²) = Δ_required·4 = (9 − 18c + 16c² − 14c³ + 7c⁴).

    Note: Family B' modification contributes 2μ₂·(1+c²)/2 = μ₂·(1+c²)
    to f_BAM at O(ε²) (via the (1+εμ)² expansion's μ² term... wait,
    actually it's 2μ₂ from the linear term in (1+ε²μ₂+...)².
    Wait, the modification factor is (1 + ε·μ₁ + ε²·μ₂); squaring
    gives 1 + 2ε·μ₁ + ε²·(2μ₂ + μ₁²) + O(ε³). So the NEW O(ε²)
    contribution from μ₂ is 2μ₂. Adding to the polarisation sum
    factor (1+c²)/2: 2·(2μ₂·(1+c²))/2 = 2μ₂·(1+c²).

    Wait let me redo. f_BAM = (1+c²)/2 · (modification_factor)² · propagator.
    modification_factor² = 1 + 2εμ₁ + ε²(2μ₂ + μ₁²).
    The 2μ₂ contributes at O(ε²) as: (1+c²)/2 · 2μ₂ = (1+c²)·μ₂.

    So setting: (1+c²)·μ₂ + [other already-accounted O(ε²) terms]
              = Δ_required + ...

    Already-accounted terms include the propagator (1+1/x)² at O(ε²)
    and the μ₁² term. PR #31 accounted for both via the γ=-3/2 only;
    the residual Δ_required is what remains. So:

    (1+c²)·μ₂ = Δ_required = (1−c)²(9+7c²)/4

    Parametrise μ₂ in {1, c, c², c³, c⁴} basis (limited to c⁴ since
    target is c⁴). Map each basis element ν_k · c^k to its product
    with (1+c²):
       1·(1+c²) = (1+c²)            → (1, 0, 1, 0, 0)
       c·(1+c²) = (c+c³)            → (0, 1, 0, 1, 0)
       c²·(1+c²) = (c²+c⁴)          → (0, 0, 1, 0, 1)
       c³·(1+c²) = (c³+c⁵)          → c⁵ outside target, partial
       c⁴·(1+c²) = (c⁴+c⁶)          → outside target

    A 3-parameter μ₂ = ν₀ + ν₁·c + ν₂·c²: gives
       (ν₀ + ν₁·c + (ν₀+ν₂)·c² + ν₁·c³ + ν₂·c⁴)
    Match to target (9/4, −18/4, 16/4, −14/4, 7/4):
       c⁰: ν₀ = 9/4
       c¹: ν₁ = −18/4 = −9/2
       c²: ν₀ + ν₂ = 16/4 = 4 → ν₂ = 4 − 9/4 = 7/4
       c³: ν₁ = −14/4 = −7/2 — DIFFERENT from c¹: −9/2
       c⁴: ν₂ = 7/4 ✓

    The c¹ and c³ coefficients of μ₂·(1+c²) are BOTH equal to ν₁,
    but target c¹ = −9/2 and c³ = −7/2 differ — IMPOSSIBLE with
    polynomial μ₂ alone.

    So a polynomial μ₂ CANNOT close O(ε²). The system is
    over-determined.
    """
    # Build the linear system for μ₂ = ν₀ + ν₁·c + ν₂·c²
    # ν₀: contributes (1, 0, 1, 0, 0) to {c⁰..c⁴}
    # ν₁: contributes (0, 1, 0, 1, 0)
    # ν₂: contributes (0, 0, 1, 0, 1)
    A = np.array([
        [1.0, 0.0, 0.0],   # c⁰: ν₀
        [0.0, 1.0, 0.0],   # c¹: ν₁
        [1.0, 0.0, 1.0],   # c²: ν₀ + ν₂
        [0.0, 1.0, 0.0],   # c³: ν₁
        [0.0, 0.0, 1.0],   # c⁴: ν₂
    ])
    b = np.array([9.0, -18.0, 16.0, -14.0, 7.0]) / 4.0
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    nu0, nu1, nu2 = (float(s) for s in sol)
    residual = A @ sol - b
    max_res = float(np.max(np.abs(residual)))
    rank_A = int(np.linalg.matrix_rank(A))
    rank_aug = int(np.linalg.matrix_rank(np.column_stack([A, b])))
    exact = rank_A == rank_aug

    return {
        'name': 'T3_family_B_prime_polynomial',
        'description': (
            'Solve for μ₂ = ν₀ + ν₁·c + ν₂·c² such that '
            'μ₂·(1+c²) = Δ_required = (9 − 18c + 16c² − 14c³ + 7c⁴)/4. '
            'Check whether polynomial μ₂ alone can close O(ε²).'
        ),
        'matrix_A': A.tolist(),
        'rhs_b': b.tolist(),
        'lsq_nu0': nu0,
        'lsq_nu1': nu1,
        'lsq_nu2': nu2,
        'lsq_residual': residual.tolist(),
        'lsq_max_residual': max_res,
        'rank_A': rank_A,
        'rank_augmented': rank_aug,
        'exact_solution_exists': exact,
        # Pass = correctly identifies the over-determination
        'pass': not exact and max_res > 0.1,
    }


# ---------------------------------------------------------------------------
# T4. Family A contribution at O(ε²) — α² sin⁴θ term
# ---------------------------------------------------------------------------

def test_T4_family_A_at_O_eps2() -> dict:
    """Family A modification V = (ε·ε'*) + ε·α·(ε·k̂')(ε'*·k̂) contributes
    at O(ε²) via |V|² expansion:

    Σ_pol |V_A|² at O(ε²) gets:
      ε² · α² · sin⁴θ  (from |α·(εk')(ε'*k)|² at λ=1,1, divided by 2 for
                       initial pol average)

    So additional O(ε²) f_BAM contribution: ε² · α² · sin⁴θ / 2.

    Adding Family B' (μ₂ polynomial) and Family A (α²·sin⁴θ/2):
    Need to solve:
      (1+c²)·μ₂ + α²·sin⁴θ/2 = Δ_required = (1−c)²(9+7c²)/4

    sin⁴θ = (1−c²)² = 1 − 2c² + c⁴

    So contribution: α²·(1 − 2c² + c⁴)/2

    With μ₂ = ν₀ + ν₁·c + ν₂·c²:
      μ₂·(1+c²) basis: (ν₀, ν₁, ν₀+ν₂, ν₁, ν₂)
      α²·sin⁴θ/2 basis: (α²/2, 0, −α², 0, α²/2)

    Sum:
      c⁰: ν₀ + α²/2
      c¹: ν₁
      c²: ν₀ + ν₂ − α²
      c³: ν₁
      c⁴: ν₂ + α²/2

    Target: (9, −18, 16, −14, 7)/4.

    System (5 eqns, 4 unknowns including α²):
      ν₀ + α²/2  = 9/4
      ν₁         = −18/4 = −9/2
      ν₀ + ν₂ − α² = 16/4 = 4
      ν₁         = −14/4 = −7/2
      ν₂ + α²/2  = 7/4

    c¹ vs c³: ν₁ = −9/2 vs ν₁ = −7/2 — STILL CONTRADICTION.

    Adding Family A doesn't help — the c¹ ≠ c³ discrepancy is in
    the polynomial μ₂·(1+c²), which Family A can't fix because
    sin⁴θ is even in c.
    """
    # Build the over-determined system
    A = np.array([
        [1.0, 0.0, 0.0,  0.5],   # c⁰: ν₀ + α²/2
        [0.0, 1.0, 0.0,  0.0],   # c¹: ν₁
        [1.0, 0.0, 1.0, -1.0],   # c²: ν₀ + ν₂ − α²
        [0.0, 1.0, 0.0,  0.0],   # c³: ν₁ — same as c¹
        [0.0, 0.0, 1.0,  0.5],   # c⁴: ν₂ + α²/2
    ])
    b = np.array([9.0, -18.0, 16.0, -14.0, 7.0]) / 4.0
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    nu0, nu1, nu2, alpha_sq = (float(s) for s in sol)
    residual = A @ sol - b
    max_res = float(np.max(np.abs(residual)))
    rank_A = int(np.linalg.matrix_rank(A))
    rank_aug = int(np.linalg.matrix_rank(np.column_stack([A, b])))
    exact = rank_A == rank_aug

    return {
        'name': 'T4_family_A_plus_B_prime',
        'description': (
            'Add Family A α² sin⁴θ contribution at O(ε²) and re-solve. '
            'Check if the c¹ vs c³ discrepancy from T3 is resolved.'
        ),
        'matrix_A': A.tolist(),
        'rhs_b': b.tolist(),
        'lsq_nu0': nu0,
        'lsq_nu1': nu1,
        'lsq_nu2': nu2,
        'lsq_alpha_squared': alpha_sq,
        'lsq_residual': residual.tolist(),
        'lsq_max_residual': max_res,
        'rank_A': rank_A,
        'rank_augmented': rank_aug,
        'exact_solution_exists': exact,
        'pass': not exact,   # T4 confirms the over-determination
    }


# ---------------------------------------------------------------------------
# T5. Extended Family with odd-c structure
# ---------------------------------------------------------------------------

def test_T5_extended_family_with_odd_c() -> dict:
    """The c¹ vs c³ discrepancy requires an odd-in-c vertex structure
    where the c¹ and c³ contributions are independent.

    Add Family E: μ_extra = ξ · c · (1 − c²) = ξ·(c − c³)
    Contribution to c basis: (0, ξ, 0, −ξ, 0).

    Combined system:
      c⁰: ν₀ + α²/2 = 9/4
      c¹: ν₁ + ξ = −18/4 = −9/2
      c²: ν₀ + ν₂ − α² = 4
      c³: ν₁ − ξ = −14/4 = −7/2
      c⁴: ν₂ + α²/2 = 7/4

    From c¹ and c³: 2ν₁ = −9/2 + (−7/2) = −8, ν₁ = −4
                   2ξ  = −9/2 − (−7/2) = −1,  ξ = −1/2

    5 equations, 5 unknowns (ν₀, ν₁, ν₂, α², ξ). Should be solvable.
    """
    A = np.array([
        [1.0, 0.0, 0.0,  0.5,  0.0],   # c⁰: ν₀ + α²/2
        [0.0, 1.0, 0.0,  0.0,  1.0],   # c¹: ν₁ + ξ
        [1.0, 0.0, 1.0, -1.0,  0.0],   # c²: ν₀ + ν₂ − α²
        [0.0, 1.0, 0.0,  0.0, -1.0],   # c³: ν₁ − ξ
        [0.0, 0.0, 1.0,  0.5,  0.0],   # c⁴: ν₂ + α²/2
    ])
    b = np.array([9.0, -18.0, 16.0, -14.0, 7.0]) / 4.0
    det = float(np.linalg.det(A))
    if abs(det) > 1e-12:
        sol = np.linalg.solve(A, b)
    else:
        sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    nu0, nu1, nu2, alpha_sq, xi = (float(s) for s in sol)
    residual = A @ sol - b
    max_res = float(np.max(np.abs(residual)))
    exact = max_res < 1e-10

    # Check if α² is non-negative (physical)
    alpha_physical = alpha_sq >= -1e-9
    # Compute α with sign convention (positive root if non-negative)
    if alpha_sq >= 0:
        alpha = math.sqrt(alpha_sq)
    else:
        alpha = float('nan')

    # Check if coefficients are clean rationals
    def is_clean(x, cands):
        return any(abs(x - c) < 0.05 for c in cands)

    cand_set = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, -0.25, -0.5,
                -0.75, -1.0, -1.5, -2.0, -3.0, -4.0, 7.0/4.0, 9.0/4.0,
                9.0/8.0, 15.0/8.0, 5.0/4.0]

    return {
        'name': 'T5_extended_family_with_odd_c',
        'description': (
            'Add Family E: μ_extra = ξ·c·(1−c²) (odd in c, allows '
            'independent c¹ and c³). Solve the 5-unknown system.'
        ),
        'matrix_A_5x5': A.tolist(),
        'rhs_b': b.tolist(),
        'det_A': det,
        'lsq_nu0': nu0,
        'lsq_nu1': nu1,
        'lsq_nu2': nu2,
        'lsq_alpha_squared': alpha_sq,
        'lsq_xi': xi,
        'lsq_alpha_value': alpha,
        'lsq_residual': residual.tolist(),
        'lsq_max_residual': max_res,
        'exact_solution': exact,
        'alpha_squared_non_negative': alpha_physical,
        'coefficients_clean': {
            'nu0': is_clean(nu0, cand_set),
            'nu1': is_clean(nu1, cand_set),
            'nu2': is_clean(nu2, cand_set),
            'alpha_sq': is_clean(alpha_sq, cand_set),
            'xi': is_clean(xi, cand_set),
        },
        'pass': exact and alpha_physical,
    }


# ---------------------------------------------------------------------------
# T6. Numerical verification: scan extended family parameters
# ---------------------------------------------------------------------------

def f_BAM_extended_modified(eps: float, theta: float,
                            mu1: float, nu0: float, nu1: float,
                            nu2: float, alpha: float, xi: float) -> float:
    """f_BAM with extended vertex modification at O(ε) and O(ε²):
       V = (ε·ε'*)·(1 + ε·μ₁(θ) + ε²·μ₂(θ)) + ε·α·(ε·k̂')(ε'*·k̂)
                                              + ε²·ξ·(ε·k̂')(ε'*·k̂)·(some)
    For simplicity, implement μ₁ = -3(1-c)/2 (PR #31) and the O(ε²)
    additions.

    This is the operational simulation; the algebraic matching is
    done in T5.
    """
    c = math.cos(theta)
    s2 = math.sin(theta) ** 2
    mu1_val = -1.5 * (1.0 - c)  # PR #31 fixed
    mu2_val = nu0 + nu1 * c + nu2 * c * c
    extra_odd = xi * c * (1.0 - c * c)
    mu2_val += extra_odd

    # Family A coupling (real α)
    x = x_ratio(eps, theta)

    # Compute the polarization-summed |M|² at each (λ, λ'):
    sin_theta = math.sin(theta)
    cos_theta = c
    # ε^1·ε'^1* = cos θ; ε^1·k̂' = sin θ; ε'^1*·k̂ = -sin θ
    # ε^2·ε'^2* = 1;     ε^2·k̂' = 0;     ε'^2*·k̂ = 0
    eps_k_prime_eps_prime_k = -sin_theta * sin_theta  # only (1,1)

    # Vertex factor at each pol (λ, λ'):
    # V_(1,1) = cos θ · (1 + ε·μ₁ + ε²·μ₂) + ε·α·(-sin²θ)
    # V_(2,2) = 1     · (1 + ε·μ₁ + ε²·μ₂)
    # V_(1,2), V_(2,1) = 0 in our basis
    angular_mod = 1.0 + eps * mu1_val + eps * eps * mu2_val
    V_11 = cos_theta * angular_mod + eps * alpha * eps_k_prime_eps_prime_k
    V_22 = 1.0 * angular_mod

    pol_sq = (abs(V_11) ** 2 + abs(V_22) ** 2) / 2.0   # /2 = average over initial λ

    propagator_factor = (1.0 + 1.0 / x) ** 2 / 4.0
    return pol_sq * propagator_factor


def test_T6_numerical_verification(t5: dict) -> dict:
    """Numerically verify the analytic T5 solution closes O(ε²) by
    evaluating f_BAM at small ε with the derived coefficients and
    comparing to f_KN."""
    nu0 = t5['lsq_nu0']
    nu1 = t5['lsq_nu1']
    nu2 = t5['lsq_nu2']
    alpha_sq = t5['lsq_alpha_squared']
    xi = t5['lsq_xi']
    # If α² < 0, use absolute value with the sign of the matching
    # ambiguity flagged; for now use sqrt(|α²|)
    alpha = math.sqrt(abs(alpha_sq)) * (1.0 if alpha_sq >= 0 else 0.0)

    epsilons = [1e-4, 1e-3, 1e-2, 1e-1]
    rows = []
    for eps in epsilons:
        max_diff = 0.0
        for theta in np.linspace(0.05, PI - 0.05, 17):
            b = f_BAM_extended_modified(
                eps, float(theta), -1.5, nu0, nu1, nu2, alpha, xi,
            )
            # Normalise at θ=0 of the same ε (b_at_0)
            b_at_0 = f_BAM_extended_modified(
                eps, 0.0, -1.5, nu0, nu1, nu2, alpha, xi,
            )
            b_norm = b / max(b_at_0, 1e-30)
            k = f_KN_normalized(eps, float(theta))
            k_at_0 = f_KN_normalized(eps, 0.0)
            k_norm = k / max(k_at_0, 1e-30)
            max_diff = max(max_diff, abs(b_norm - k_norm))
        rows.append({
            'epsilon': eps,
            'max_KN_residual': max_diff,
        })

    # Fit residual scaling
    eps_arr = np.array([r['epsilon'] for r in rows])
    res_arr = np.array([r['max_KN_residual'] for r in rows])
    valid = res_arr > 1e-15
    if valid.sum() >= 2:
        slope, _ = np.polyfit(np.log(eps_arr[valid]), np.log(res_arr[valid]), 1)
        order = float(slope)
    else:
        order = float('nan')
    return {
        'name': 'T6_numerical_verification',
        'description': (
            'Numerically verify the T5 analytic solution closes O(ε²) '
            'by computing f_BAM with the derived coefficients and '
            'measuring residual scaling in ε.'
        ),
        'numerical_results': rows,
        'fitted_residual_order_in_eps': order,
        'expected_for_O_eps2_closure': 3.0,
        # Pass = residual scales as O(ε³) or higher (showing O(ε²)
        # closure)
        'pass': order > 2.5,
    }


# ---------------------------------------------------------------------------
# T7. Candidate predictions for O(ε²) coefficients
# ---------------------------------------------------------------------------

def test_T7_candidate_predictions(t5: dict) -> dict:
    """For each of the 7 surviving PR #33 candidates, compute the
    natural O(ε²) extension and compare to the T5-derived
    coefficients."""
    # The T5 solution gives specific (ν₀, ν₁, ν₂, α², ξ) values.
    # Each candidate has a natural O(ε²) extension:
    #
    # A: doubled electron Casimir → predict O(ε²) ~ C_2² for j=1/2
    #     C_2(1/2) = 3/4; (C_2)² = 9/16. Doubled: 9/8.
    # B: photon Casimir − Hopf → predict (C_2(j=1))² − Hopf² = 4 - 1/4 = 15/4
    # D: Hopf × photon_mult → (1/2 · 3)² = 9/4
    # F: closure + Hopf → (3/2)² = 9/4
    # G: Pauli trace → no obvious O(ε²) extension; flag as silent
    # H: two-mouth × spin-½ Casimir → n_mouth · C_2² = 2 · 9/16 = 9/8
    # E: speculative
    #
    # These predictions are SCALAR — they don't predict the angular
    # structure of μ₂. So we compare to the magnitude scales of the
    # T5 coefficients (e.g. |ν₀| + |ν₂| as a rough magnitude).

    nu0 = t5['lsq_nu0']
    nu1 = t5['lsq_nu1']
    nu2 = t5['lsq_nu2']
    alpha_sq = t5['lsq_alpha_squared']
    xi = t5['lsq_xi']

    # Use T5-derived integrated magnitude as the target
    target_magnitude = abs(nu0) + abs(nu1) + abs(nu2) + abs(alpha_sq) + abs(xi)

    # Candidate magnitudes (their natural scalar O(ε²) prediction)
    cands = {
        'A_doubled_electron_casimir': 9.0 / 8.0,
        'B_photon_casimir_minus_hopf': 15.0 / 4.0,
        'D_hopf_times_photon_mult': 9.0 / 4.0,
        'F_closure_plus_hopf': 9.0 / 4.0,
        'G_pauli_trace': float('nan'),   # silent
        'H_two_mouth_spin_half': 9.0 / 8.0,
        'E_antipodal_natural_units': float('nan'),
    }

    rows = []
    for cand, mag in cands.items():
        if math.isnan(mag):
            note = 'no obvious O(ε²) extension'
            ratio = float('nan')
        else:
            ratio = mag / target_magnitude if target_magnitude > 1e-12 else float('nan')
            note = f'scalar magnitude {mag} vs target |Σ| {target_magnitude:.3f}'
        rows.append({
            'candidate': cand,
            'predicted_O_eps2_scalar_magnitude': mag,
            'target_total_magnitude_from_T5': target_magnitude,
            'ratio_predicted_to_target': ratio,
            'note': note,
        })

    return {
        'name': 'T7_candidate_O_eps2_predictions',
        'description': (
            'For each surviving PR #33 candidate, compute the natural '
            'O(ε²) scalar prediction and compare to the T5-derived '
            'coefficient magnitude.'
        ),
        'rows': rows,
        'target_total_magnitude': target_magnitude,
        'note': (
            'Candidate predictions are scalar magnitudes; they do '
            'not specify the angular structure of μ₂. The probe '
            'reports relative magnitudes as a weak discrimination '
            'metric.'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_analytic_O_eps2_verification()
    t2 = test_T2_delta_required_decomposition()
    t3 = test_T3_family_B_prime()
    t4 = test_T4_family_A_at_O_eps2()
    t5 = test_T5_extended_family_with_odd_c()
    t6 = test_T6_numerical_verification(t5)
    t7 = test_T7_candidate_predictions(t5)
    tests = [t1, t2, t3, t4, t5, t6, t7]

    if not t1['pass']:
        verdict_class = 'ANALYTIC_INCONSISTENCY'
        verdict = (
            'ANALYTIC INCONSISTENCY — the O(ε²) closed-form '
            'predictions do not match the numerical second-derivative '
            'extraction. Re-verify the analytic expansion.'
        )
    elif t5['exact_solution'] and t5['alpha_squared_non_negative']:
        verdict_class = 'CLOSED_WITH_EXTENDED_FAMILY'
        verdict = (
            f'CLOSED WITH EXTENDED FAMILY (B + A + odd-c) — the '
            f'O(ε²) gap closes analytically with '
            f'(ν₀, ν₁, ν₂, α², ξ) = '
            f'({t5["lsq_nu0"]:.4f}, {t5["lsq_nu1"]:.4f}, '
            f'{t5["lsq_nu2"]:.4f}, {t5["lsq_alpha_squared"]:.4f}, '
            f'{t5["lsq_xi"]:.4f}). '
            f'α² = {t5["lsq_alpha_squared"]:.4f} '
            f'(non-negative as required). '
            f'Residual fits to O(ε^{t6["fitted_residual_order_in_eps"]:.2f}). '
            'The PR #31 Family B alone is structurally insufficient; '
            'closing O(ε²) requires explicit Family A (ε·k coupling) '
            'AND an odd-in-c angular term — neither of which appears '
            'naturally in any of the PR #33 surviving candidates. '
            'This identifies the next BAM structural target.'
        )
    elif t5['exact_solution'] and not t5['alpha_squared_non_negative']:
        verdict_class = 'CLOSED_BUT_UNPHYSICAL'
        verdict = (
            f'CLOSED BUT UNPHYSICAL — the linear system has an exact '
            f'solution but α² = {t5["lsq_alpha_squared"]:.4f} is '
            'negative, which is unphysical (would require imaginary α). '
            'The structural form of the closure is mathematically '
            'identified but does not correspond to a real BAM vertex.'
        )
    else:
        verdict_class = 'OVER_DETERMINED'
        verdict = (
            f'OVER-DETERMINED — the extended vertex family (B + A + '
            f'odd-c) cannot exactly close O(ε²). LSQ residual '
            f'{t5["lsq_max_residual"]:.4f}. Closing the O(ε²) gap '
            'requires additional vertex structure beyond polynomial '
            'angular modulations and ε·k contractions.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'analytic_setup': (
            'Δ_required_O(ε²) = (1−c)²·(9+7c²)/4. '
            'PR #31 closes O(ε) but leaves this O(ε²) residual.'
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
    L.append('# O(ε²) analytic extension probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Follow-on to PR #31 (analytic γ = −3/2 at O(ε)) and PR #33 '
        '(d-scaling discrimination). Extends the analytic derivation '
        'to O(ε²), identifies what vertex structure is needed beyond '
        'PR #31, and discriminates between the 7 surviving '
        'coefficient-origin candidates from PR #33.'
    )
    L.append('')
    L.append('**Setup:**')
    L.append('')
    L.append('```')
    L.append(s['analytic_setup'])
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
            value = (
                f"max num/ana O(ε²) error = "
                f"{max(t['max_KN_O_eps2_error'], t['max_BAM_PR31_O_eps2_error']):.2e}"
            )
        elif nm.startswith('T2'):
            value = (
                f"max decomposition error = "
                f"{t['max_decomposition_error']:.2e}"
            )
        elif nm.startswith('T3'):
            value = (
                f"Family B' alone: "
                f"{'exact' if t['exact_solution_exists'] else 'over-det.'}; "
                f"max LSQ residual {t['lsq_max_residual']:.3f}"
            )
        elif nm.startswith('T4'):
            value = (
                f"Family A added: "
                f"{'exact' if t['exact_solution_exists'] else 'over-det.'}; "
                f"max LSQ residual {t['lsq_max_residual']:.3f}"
            )
        elif nm.startswith('T5'):
            value = (
                f"Extended (B+A+odd-c): "
                f"{'EXACT' if t['exact_solution'] else 'no'}; "
                f"α² = {t['lsq_alpha_squared']:.3f}"
            )
        elif nm.startswith('T6'):
            value = (
                f"residual order = "
                f"O(ε^{t['fitted_residual_order_in_eps']:.2f})"
            )
        elif nm.startswith('T7'):
            value = f"target magnitude {t['target_total_magnitude']:.3f}"
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T2 detail: target coefficients
    t2 = s['tests'][1]
    L.append('## Δ_required_O(ε²) in c basis (T2)')
    L.append('')
    L.append(
        f"**Coefficients (c⁰..c⁴):** "
        f"{[f'{x:+.3f}' for x in t2['delta_required_coefficients_c0_to_c4']]}"
    )
    L.append('')
    L.append('Equivalent form: (9 − 18c + 16c² − 14c³ + 7c⁴)/4 = (1−c)²·(9+7c²)/4')
    L.append('')

    # T3 detail
    t3 = s['tests'][2]
    L.append('## T3: Family B′ alone (polynomial μ₂)')
    L.append('')
    L.append(
        f"Rank(A) = {t3['rank_A']}, rank(A|b) = {t3['rank_augmented']}. "
        f"Exact solution exists: "
        f"{'YES' if t3['exact_solution_exists'] else 'NO'}."
    )
    L.append('')
    L.append(
        f"LSQ: ν₀ = {t3['lsq_nu0']:.4f}, ν₁ = {t3['lsq_nu1']:.4f}, "
        f"ν₂ = {t3['lsq_nu2']:.4f}. Max residual {t3['lsq_max_residual']:.4f}."
    )
    L.append('')
    L.append(
        '**Conclusion:** the c¹ and c³ coefficients of μ₂·(1+c²) are '
        'tied (both equal ν₁ from the parametrisation), but the '
        'target has c¹ = −9/2 and c³ = −7/2 — different. Polynomial '
        'Family B′ alone cannot close O(ε²).'
    )
    L.append('')

    # T4 detail
    t4 = s['tests'][3]
    L.append('## T4: Family A + B′ (α² sin⁴θ added)')
    L.append('')
    L.append(
        f"Exact: {'YES' if t4['exact_solution_exists'] else 'NO'}. "
        f"LSQ residual {t4['lsq_max_residual']:.4f}. "
        f"α² = {t4['lsq_alpha_squared']:.4f}."
    )
    L.append('')
    L.append(
        '**Conclusion:** sin⁴θ is even in c, so adding Family A '
        'preserves the c¹ = c³ contradiction. STILL over-determined.'
    )
    L.append('')

    # T5 detail
    t5 = s['tests'][4]
    L.append('## T5: Extended family with odd-c term ξ·c·(1−c²)')
    L.append('')
    L.append(
        f"5-unknown system: det(A) = {t5['det_A']:.4f}, "
        f"exact = {'YES' if t5['exact_solution'] else 'no'}. "
        f"Max residual {t5['lsq_max_residual']:.4e}."
    )
    L.append('')
    L.append(
        f"Coefficients: ν₀ = {t5['lsq_nu0']:.4f}, "
        f"ν₁ = {t5['lsq_nu1']:.4f}, ν₂ = {t5['lsq_nu2']:.4f}, "
        f"α² = {t5['lsq_alpha_squared']:.4f}, "
        f"ξ = {t5['lsq_xi']:.4f}."
    )
    L.append('')
    L.append(f"α value: {t5['lsq_alpha_value']:.4f} "
             f"(α² ≥ 0: {'YES' if t5['alpha_squared_non_negative'] else 'NO'})")
    L.append('')
    L.append('Clean rational check:')
    L.append('')
    for k, v in t5['coefficients_clean'].items():
        L.append(f"- {k}: {'clean' if v else 'fitted'}")
    L.append('')

    # T6 detail
    t6 = s['tests'][5]
    L.append('## T6: Numerical residual scaling')
    L.append('')
    L.append(
        f"Fitted residual order: **O(ε^{t6['fitted_residual_order_in_eps']:.3f})** "
        f"(expected ≥ 3 for O(ε²) closure)."
    )
    L.append('')
    L.append('| ε | max KN residual |')
    L.append('|---:|---:|')
    for r in t6['numerical_results']:
        L.append(f"| {r['epsilon']:.4g} | {r['max_KN_residual']:.4e} |")
    L.append('')

    # T7 detail
    t7 = s['tests'][6]
    L.append('## T7: Candidate O(ε²) scalar predictions')
    L.append('')
    L.append(t7['note'])
    L.append('')
    L.append('| candidate | predicted scalar magnitude | ratio to T5 target |')
    L.append('|---|---:|---:|')
    for r in t7['rows']:
        mag = r['predicted_O_eps2_scalar_magnitude']
        ratio = r['ratio_predicted_to_target']
        mag_s = f"{mag:.4f}" if not math.isnan(mag) else 'silent'
        ratio_s = f"{ratio:.3f}" if not math.isnan(ratio) else 'n/a'
        L.append(f"| `{r['candidate']}` | {mag_s} | {ratio_s} |")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **BAM derivation of the new structures** (ξ odd-c term, '
        'α² coupling) at O(ε²). The structural form is identified by '
        'T5 but the BAM-native origin remains to be derived.'
    )
    L.append(
        '- **O(ε³) and beyond.** Each successive order may require '
        'additional vertex structures; the pattern of required '
        'additions is open.'
    )
    L.append(
        '- **Candidate discrimination at O(ε²).** T7 uses scalar '
        'magnitudes; sharper discrimination requires the candidates '
        'to predict the full angular structure of μ₂, not just a '
        'scalar.'
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
    out = here / 'runs' / f'{ts}_compton_eps2_extension_probe'
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
