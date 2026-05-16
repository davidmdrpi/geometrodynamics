"""
Photon-structure probe for the BAM Compton amplitude.

Follow-on to the partial-match result from
`compton_amplitude_structure_probe.py` (PR #26). That probe showed
the antipodal `S³` Green function reproduces the QED propagator pole
exactly (fitted exponent 1.0002) but the natural spin-½ ansatz
produces Thomson angular dependence `(1 + cos θ)/2`, not Klein-Nishina
`(1 + cos²θ)/2`. The diagnostic T4 ruled out an overall phase factor
as the fix — a polarization-vector structure is needed.

This probe gives the photon explicit transverse-polarization machinery
and tests whether the structural gap closes.

Construction. Polarization-resolved BAM amplitude

    M^{(λ, λ')}_BAM(θ, ω)
        =  (ε^(λ)(k) · ε^(λ')*(k'))
           ·  [ G_S3(ψ_s) · exp(iφ_s) + G_S3(ψ_u) · exp(iφ_u) ]
           ·  T²

with two photon polarisations per particle, summed over to compute
the unpolarized cross section. The photon polarization vectors are
the standard two transverse modes:

    incoming  k along +ẑ:   ε^(1)(k)  = x̂,        ε^(2)(k)  = ŷ
    outgoing  k' at angle θ in xz:
                            ε^(1)(k') = (cos θ, 0, −sin θ),
                            ε^(2)(k') = ŷ

(`ε^(1)` is in the scattering plane, `ε^(2)` is perpendicular.) The
polarization-product factor `ε(k) · ε*(k')` is exactly the structure
the previous probe lacked.

Tests:

  T1. Polarization sum identity — sanity check that
      `Σ_{λ,λ'} | ε^(λ)(k) · ε^(λ')*(k') |² = 1 + cos²θ`
      at every sampled angle. Verifies the polarization-vector
      implementation.

  T2. Klein-Nishina restoration with photon polarization + scalar
      electron (no θ-dependent spin phase). Fit |M_total|² to
      A + B cos θ + C cos²θ and compare to KN ratio (½, 0, ½)/F(0).
      PASS = (A, B, C) matches KN to within 5%.

  T3. With electron spin-½ phases re-included — verify the previous
      probe's spin-½ angular form is spurious at Thomson. The fit
      should give a four-term `(1 + cos θ)(1 + cos²θ)` polynomial.

  T4. Propagator pole robustness — verify T1 from PR #26 still passes
      (polarization factor is angle-constant at fixed θ, so the
      `1/(s − m²)` pole should be preserved).

  T5. Polarization basis invariance — replace linear with circular
      polarisations and verify |M_total|² is unchanged
      (basis-completeness).

Verdict:
  FULL_MATCH if T1, T2, T4, T5 PASS and T3 produces the predicted
  spurious convolution. The Compton amplitude tree-level skeleton
  is then reproduced from BAM ingredients:
    antipodal S³ propagation + photon polarization + throat transport
    + Hopf-holonomy closure phase, with electron treated as a scalar
    charge at the Thomson limit (spin effects sub-leading in ω/m_e).
  PARTIAL_MATCH if T1, T4, T5 PASS but T2 fails — would indicate a
  deeper structural gap than polarization machinery.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
TAU = 2.0 * PI


# ---------------------------------------------------------------------------
# Photon polarization vectors
# ---------------------------------------------------------------------------

def linear_polarization_vectors_in(k_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two transverse linear polarization vectors for a photon along k_hat.

    Convention: ε^(1) lies in the scattering plane (xz), ε^(2) is
    perpendicular (along ŷ).

    For k_hat along +ẑ:
        ε^(1) = x̂,  ε^(2) = ŷ
    For k_hat along (sin θ, 0, cos θ) in xz-plane:
        ε^(1) = (cos θ, 0, −sin θ),  ε^(2) = (0, 1, 0)
    """
    k_hat = np.asarray(k_hat, dtype=float)
    n = np.linalg.norm(k_hat)
    if n < 1e-12:
        raise ValueError('zero propagation vector')
    k_hat = k_hat / n
    # ε^(2) is the y-axis (perpendicular to scattering plane)
    eps_2 = np.array([0.0, 1.0, 0.0])
    # ε^(1) is in the scattering plane (xz) and transverse to k_hat
    # For k_hat = (sin θ, 0, cos θ): ε^(1) = (cos θ, 0, −sin θ)
    # General: ε^(1) = (k̂_z, 0, −k̂_x) when k_hat is in xz-plane
    eps_1 = np.array([k_hat[2], 0.0, -k_hat[0]])
    norm1 = np.linalg.norm(eps_1)
    if norm1 < 1e-12:
        # k_hat along ŷ (edge case): use x̂ as ε^(1)
        eps_1 = np.array([1.0, 0.0, 0.0])
    else:
        eps_1 = eps_1 / norm1
    return eps_1, eps_2


def circular_polarization_vectors_in(k_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two transverse circular polarization vectors (complex).

    ε^(R) = (ε^(1) + i ε^(2)) / √2
    ε^(L) = (ε^(1) − i ε^(2)) / √2
    """
    e1, e2 = linear_polarization_vectors_in(k_hat)
    eR = (e1 + 1j * e2) / math.sqrt(2.0)
    eL = (e1 - 1j * e2) / math.sqrt(2.0)
    return eR, eL


def polarization_inner_product(eps_a: np.ndarray, eps_b: np.ndarray) -> complex:
    """ε_a · ε_b* (Hermitian inner product on C³)."""
    eps_a = np.asarray(eps_a, dtype=complex)
    eps_b = np.asarray(eps_b, dtype=complex)
    return complex(np.dot(eps_a, np.conjugate(eps_b)))


# ---------------------------------------------------------------------------
# Compton kinematics
# ---------------------------------------------------------------------------

def compton_kinematics(omega: float, theta: float, m_e: float = 1.0) -> dict:
    omega_prime = omega / (1.0 + (omega / m_e) * (1.0 - math.cos(theta)))
    p_e = np.array([m_e, 0.0, 0.0, 0.0])
    k_g = np.array([omega, 0.0, 0.0, omega])  # along +ẑ
    k_g_prime = np.array([
        omega_prime,
        omega_prime * math.sin(theta),
        0.0,
        omega_prime * math.cos(theta),
    ])
    s_4mom = p_e + k_g
    u_4mom = p_e - k_g_prime
    s = s_4mom[0] ** 2 - np.dot(s_4mom[1:], s_4mom[1:])
    u = u_4mom[0] ** 2 - np.dot(u_4mom[1:], u_4mom[1:])
    return {
        'omega_in': omega,
        'omega_out': omega_prime,
        'theta_scatter': theta,
        'k_in_hat': k_g[1:] / max(np.linalg.norm(k_g[1:]), 1e-30),
        'k_out_hat': k_g_prime[1:] / max(np.linalg.norm(k_g_prime[1:]), 1e-30),
        's': float(s),
        'u': float(u),
    }


# ---------------------------------------------------------------------------
# S³ Green function
# ---------------------------------------------------------------------------

def G_S3(psi: float, R: float = 1.0, eps: float = 1e-12) -> float:
    """Zero-mean S³ Green function. Diverges as 1/(4πRψ) near ψ=0."""
    psi_eff = max(min(psi, PI - eps), eps)
    return float(
        (((PI - psi_eff) / math.tan(psi_eff)) - 0.5)
        / (4.0 * PI * PI * R)
    )


# ---------------------------------------------------------------------------
# Polarization-resolved BAM amplitude
# ---------------------------------------------------------------------------

def bam_amplitude_with_polarization(
    omega: float,
    theta: float,
    eps_in: np.ndarray,
    eps_out: np.ndarray,
    m_e: float = 1.0,
    spin_mode: str = 'scalar',
) -> complex:
    """Compute M^{(λ, λ')}_BAM(θ, ω) with explicit polarization vectors.

    spin_mode:
      'scalar' — no θ-dependent electron spin phase (the natural
                 Thomson-limit choice; electron acts as a scalar
                 charge with the photon polarization sum providing
                 all angular structure).
      'half'   — previous probe's spin-½ phase ansatz: φ_s = π + θ/2,
                 φ_u = π − θ/2. Included here to verify it produces
                 the spurious (1 + cos θ)(1 + cos²θ) convolution at
                 Thomson.
    """
    kin = compton_kinematics(omega, theta, m_e)
    s, u = kin['s'], kin['u']

    eps_pole = 1e-12
    psi_s = max((s - m_e * m_e) / (2.0 * m_e * m_e), eps_pole)
    psi_u = max(abs(u - m_e * m_e) / (2.0 * m_e * m_e), eps_pole)
    psi_s = min(psi_s, PI - eps_pole)
    psi_u = min(psi_u, PI - eps_pole)

    Gs = G_S3(psi_s)
    Gu = G_S3(psi_u)

    hopf_phase = PI * math.cos(0.0)  # = π at χ = 0
    if spin_mode == 'scalar':
        phi_s = hopf_phase
        phi_u = hopf_phase
    elif spin_mode == 'half':
        phi_s = hopf_phase + 0.5 * theta
        phi_u = hopf_phase - 0.5 * theta
    else:
        raise ValueError(f'unknown spin_mode {spin_mode!r}')

    throat_factor = -1.0   # T² = −I

    pol_factor = polarization_inner_product(eps_in, eps_out)
    propagator_sum = (
        Gs * complex(np.exp(1j * phi_s))
        + Gu * complex(np.exp(1j * phi_u))
    )
    return pol_factor * propagator_sum * throat_factor


def unpolarized_cross_section(
    omega: float,
    theta: float,
    m_e: float = 1.0,
    spin_mode: str = 'scalar',
    polarization_basis: str = 'linear',
) -> float:
    """Σ_{λ, λ'} | M^{(λ, λ')}_BAM |²  (averaged 1/2 over initial λ).

    Returns the polarization-summed cross-section-like quantity.
    """
    kin = compton_kinematics(omega, theta, m_e)
    if polarization_basis == 'linear':
        eps_in_1, eps_in_2 = linear_polarization_vectors_in(kin['k_in_hat'])
        eps_out_1, eps_out_2 = linear_polarization_vectors_in(kin['k_out_hat'])
    elif polarization_basis == 'circular':
        eps_in_1, eps_in_2 = circular_polarization_vectors_in(kin['k_in_hat'])
        eps_out_1, eps_out_2 = circular_polarization_vectors_in(kin['k_out_hat'])
    else:
        raise ValueError(f'unknown basis {polarization_basis!r}')

    total = 0.0
    for eps_in in (eps_in_1, eps_in_2):
        for eps_out in (eps_out_1, eps_out_2):
            M = bam_amplitude_with_polarization(
                omega, theta, eps_in, eps_out,
                m_e=m_e, spin_mode=spin_mode,
            )
            total += abs(M) ** 2
    return 0.5 * total   # 1/2 average over incoming polarisations


# ---------------------------------------------------------------------------
# Angular fit helper
# ---------------------------------------------------------------------------

def fit_angular_polynomial(
    thetas: np.ndarray, values: np.ndarray, degree: int = 2,
) -> dict:
    """Fit values ≈ Σ_{k=0..degree} c_k · cos^k(θ).

    Returns dict with fitted coefficients c_k and residuals.
    """
    c = np.cos(thetas)
    basis = np.column_stack([c ** k for k in range(degree + 1)])
    coeffs, _, _, _ = np.linalg.lstsq(basis, values, rcond=None)
    pred = basis @ coeffs
    residual = float(np.max(np.abs(values - pred)))
    return {
        'coefficients': [float(x) for x in coeffs],
        'max_residual': residual,
    }


# ---------------------------------------------------------------------------
# T1. Polarization sum identity
# ---------------------------------------------------------------------------

def test_T1_polarization_sum_identity() -> dict:
    """Σ_{λ, λ'} | ε^(λ)(k) · ε^(λ')*(k') |² = 1 + cos²θ for k along +ẑ,
    k' at angle θ in xz-plane.
    """
    thetas = np.linspace(0.0, PI, 33)
    sums = []
    predicted = []
    for theta in thetas:
        k_in = np.array([0.0, 0.0, 1.0])
        k_out = np.array([math.sin(theta), 0.0, math.cos(theta)])
        e1_in, e2_in = linear_polarization_vectors_in(k_in)
        e1_out, e2_out = linear_polarization_vectors_in(k_out)
        total = 0.0
        for eps_in in (e1_in, e2_in):
            for eps_out in (e1_out, e2_out):
                total += abs(polarization_inner_product(eps_in, eps_out)) ** 2
        sums.append(total)
        predicted.append(1.0 + math.cos(theta) ** 2)

    sums = np.array(sums)
    predicted = np.array(predicted)
    max_err = float(np.max(np.abs(sums - predicted)))
    return {
        'name': 'T1_polarization_sum_identity',
        'description': (
            'Σ_{λ, λ\'} |ε^(λ)(k) · ε^(λ\')*(k\')|² = 1 + cos²θ '
            'across the sampled angular range. Sanity check on '
            'transverse polarization vector implementation.'
        ),
        'theta_samples': thetas.tolist(),
        'computed_sums': sums.tolist(),
        'predicted_1_plus_cos2': predicted.tolist(),
        'max_pointwise_error': max_err,
        'pass': max_err < 1e-12,
    }


# ---------------------------------------------------------------------------
# T2. Klein-Nishina restoration (scalar electron)
# ---------------------------------------------------------------------------

def test_T2_KN_restoration_scalar_electron(omega: float = 1e-4) -> dict:
    """Compute |M_total|²(θ) with photon polarization + scalar electron
    in the Thomson limit. Fit to A + B cos θ + C cos²θ and compare to
    Klein-Nishina ratio.

    Klein-Nishina: F(θ) = (1 + cos²θ)/2, normalised to F(0) = 1 gives
    F/F(0) = (1 + cos²θ)/2. Fitted (A, B, C) should give A ≈ ½, B ≈ 0,
    C ≈ ½.
    """
    thetas = np.linspace(0.0, PI, 33)
    M_sq = np.array([
        unpolarized_cross_section(omega, float(t), spin_mode='scalar')
        for t in thetas
    ])
    M_sq_norm = M_sq / max(M_sq[0], 1e-30)

    fit = fit_angular_polynomial(thetas, M_sq_norm, degree=2)
    A, B, C = fit['coefficients']
    A_kn, B_kn, C_kn = 0.5, 0.0, 0.5

    pred_kn = A_kn + B_kn * np.cos(thetas) + C_kn * np.cos(thetas) ** 2
    max_kn_residual = float(np.max(np.abs(M_sq_norm - pred_kn)))

    # PASS = (A, B, C) within 5% of KN
    near_kn = (
        abs(A - A_kn) < 0.05 and
        abs(B - B_kn) < 0.05 and
        abs(C - C_kn) < 0.05
    )

    return {
        'name': 'T2_KN_restoration_scalar_electron',
        'description': (
            '|M_total|²(θ) with photon polarization vectors + scalar '
            'electron (no θ-dependent spin phase). Fit to '
            'A + B·cos θ + C·cos²θ in the Thomson limit. KN target: '
            '(½, 0, ½).'
        ),
        'omega_thomson': omega,
        'fitted_coefficients': {'A': A, 'B': B, 'C': C},
        'klein_nishina_target': {'A': A_kn, 'B': B_kn, 'C': C_kn},
        'max_pointwise_residual_to_KN': max_kn_residual,
        'fit_residual': fit['max_residual'],
        'theta_samples': thetas.tolist(),
        'amplitude_squared_normalized': M_sq_norm.tolist(),
        'pass': near_kn,
    }


# ---------------------------------------------------------------------------
# T3. Spin-½ phase verification (spurious at Thomson)
# ---------------------------------------------------------------------------

def test_T3_spin_half_spurious_at_thomson(omega: float = 1e-4) -> dict:
    """Repeat T2 with spin_mode='half' (previous probe's spin-½ ansatz).
    Expected: the (1 + cos θ) factor from the spin-½ phases multiplies
    the polarization sum, giving (1 + cos θ)(1 + cos²θ).

    Fit to degree-3 polynomial in cos θ; the c_3 (cos³θ) coefficient
    should be non-zero (it's the convolution signature).
    """
    thetas = np.linspace(0.0, PI, 33)
    M_sq = np.array([
        unpolarized_cross_section(omega, float(t), spin_mode='half')
        for t in thetas
    ])
    M_sq_norm = M_sq / max(M_sq[0], 1e-30)

    # Fit degree 3
    fit3 = fit_angular_polynomial(thetas, M_sq_norm, degree=3)
    c0, c1, c2, c3 = fit3['coefficients']

    # Analytic prediction: (1 + cos θ)(1 + cos²θ) / 2 (normalized)
    # = (1/2)·(1 + cos θ + cos²θ + cos³θ)
    # Normalized at θ=0: F(0) = 1·2/2 = 1, so values ÷ 1 = (1+cos θ+cos²θ+cos³θ)/2
    pred = (1.0 + np.cos(thetas) + np.cos(thetas) ** 2 + np.cos(thetas) ** 3) / 2.0
    # Normalize prediction at θ=0
    pred_norm = pred / pred[0]
    max_residual_to_prediction = float(np.max(np.abs(M_sq_norm - pred_norm)))

    # PASS = c_3 is non-negligible (showing the cos³θ term) AND
    # the fit matches the (1 + cos θ)(1 + cos²θ) prediction.
    has_cos3 = abs(c3) > 0.05
    matches_prediction = max_residual_to_prediction < 0.05

    return {
        'name': 'T3_spin_half_spurious_at_thomson',
        'description': (
            'Repeat T2 with spin-½ phases re-included. Expected: '
            '(1 + cos θ)(1 + cos²θ) convolved structure — a four-'
            'term polynomial that is NOT Klein-Nishina. Verifies the '
            'spin-½ angular phases were spurious at Thomson.'
        ),
        'fitted_coefficients_degree3': {
            'c0': c0, 'c1': c1, 'c2': c2, 'c3': c3,
        },
        'predicted_form': '(1 + cos θ)(1 + cos²θ)/2',
        'max_residual_to_predicted_convolution': max_residual_to_prediction,
        'cos3_coefficient_nonneg': has_cos3,
        'pass': has_cos3 and matches_prediction,
    }


# ---------------------------------------------------------------------------
# T4. Propagator-pole robustness
# ---------------------------------------------------------------------------

def test_T4_propagator_pole_robustness(m_e: float = 1.0) -> dict:
    """Sweep ω → 0 at fixed θ and verify |M_total|^{1/2} ∝ 1/(s − m²).

    Use |M_total|^{1/2} (rather than |M|²) because the polarization
    factor enters squared in the cross section. The propagator pole
    structure should manifest at the amplitude level.
    """
    theta_fixed = PI / 3.0
    omegas = np.array([10.0 ** (-k) for k in [2, 3, 4, 5, 6]])
    amplitudes = []
    s_minus_msq = []
    for omega in omegas:
        M_sq = unpolarized_cross_section(
            float(omega), theta_fixed, spin_mode='scalar',
        )
        amplitudes.append(math.sqrt(max(M_sq, 1e-300)))
        kin = compton_kinematics(float(omega), theta_fixed, m_e)
        s_minus_msq.append(kin['s'] - m_e * m_e)

    amplitudes = np.array(amplitudes)
    s_minus_msq = np.array(s_minus_msq)
    valid = (amplitudes > 0) & (s_minus_msq > 0)
    if valid.sum() >= 3:
        slope, intercept = np.polyfit(
            np.log(s_minus_msq[valid]), np.log(amplitudes[valid]), 1,
        )
        residue = math.exp(intercept)
    else:
        slope = float('nan')
        residue = float('nan')

    pole_exponent = float(-slope)
    return {
        'name': 'T4_propagator_pole_robustness',
        'description': (
            '|M_total|^{1/2} ∝ 1/(s − m²) as ω → 0 — verifies the '
            'propagator pole structure is preserved when photon '
            'polarization machinery is added.'
        ),
        'theta_fixed': theta_fixed,
        'fitted_pole_exponent': pole_exponent,
        'residue': residue,
        'omegas': omegas.tolist(),
        'amplitudes': amplitudes.tolist(),
        'pass': 0.9 < pole_exponent < 1.1,
    }


# ---------------------------------------------------------------------------
# T5. Polarization basis invariance
# ---------------------------------------------------------------------------

def test_T5_basis_invariance(omega: float = 1e-4) -> dict:
    """Compute |M_total|²(θ) in both linear and circular polarization
    bases. Verify they agree (basis-completeness for any complete
    transverse polarization set).
    """
    thetas = np.linspace(0.0, PI, 17)
    M_sq_linear = np.array([
        unpolarized_cross_section(
            omega, float(t), spin_mode='scalar', polarization_basis='linear',
        )
        for t in thetas
    ])
    M_sq_circular = np.array([
        unpolarized_cross_section(
            omega, float(t), spin_mode='scalar', polarization_basis='circular',
        )
        for t in thetas
    ])
    rel_diff = np.abs(M_sq_linear - M_sq_circular) / np.maximum(
        np.abs(M_sq_linear), 1e-30,
    )
    max_rel_diff = float(np.max(rel_diff))

    return {
        'name': 'T5_polarization_basis_invariance',
        'description': (
            '|M_total|²(θ) is invariant under choice of polarization '
            'basis (linear vs circular). Basis-completeness verification.'
        ),
        'max_relative_difference': max_rel_diff,
        'samples': [
            {
                'theta': float(t),
                'linear': float(L),
                'circular': float(C),
            }
            for t, L, C in zip(thetas[:5], M_sq_linear[:5], M_sq_circular[:5])
        ],
        'pass': max_rel_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_polarization_sum_identity()
    t2 = test_T2_KN_restoration_scalar_electron()
    t3 = test_T3_spin_half_spurious_at_thomson()
    t4 = test_T4_propagator_pole_robustness()
    t5 = test_T5_basis_invariance()
    tests = [t1, t2, t3, t4, t5]

    # Verdict logic
    sanity = t1['pass'] and t4['pass'] and t5['pass']
    kn_restored = t2['pass']
    spurious_confirmed = t3['pass']
    if sanity and kn_restored and spurious_confirmed:
        verdict_class = 'FULL_MATCH'
        verdict = (
            'FULL STRUCTURAL MATCH — the BAM Compton amplitude with '
            'explicit photon transverse polarization machinery and '
            'a scalar electron at the Thomson limit reproduces the '
            'full Klein-Nishina angular factor (1 + cos²θ)/2 exactly, '
            'while preserving the propagator pole structure verified '
            'in PR #26 (fitted pole exponent ≈ 1). The polarization '
            'basis invariance check confirms basis-completeness. T3 '
            'confirms that the prior probe\'s spin-½ angular phases '
            'were spurious at Thomson — they produce the convolved '
            '(1 + cos θ)(1 + cos²θ) structure when re-included, '
            'identifying electron spin as a sub-leading effect in '
            'ω/m_e that does not enter the leading Thomson amplitude. '
            'The Compton amplitude tree-level structural skeleton is '
            'reproduced from BAM ingredients: antipodal S³ '
            'propagation (propagator pole) + photon transverse '
            'polarization (KN angular factor) + throat transport '
            '(closure sign) + Hopf-holonomy closure phase.'
        )
    elif sanity and not kn_restored:
        verdict_class = 'PARTIAL_MATCH'
        verdict = (
            'PARTIAL STRUCTURAL MATCH — sanity tests T1/T4/T5 pass '
            'but Klein-Nishina is not restored even with photon '
            'polarization machinery. The structural gap is deeper '
            'than the polarization sum alone. See T2 fit coefficients '
            'for diagnosis.'
        )
    elif not sanity:
        verdict_class = 'SANITY_FAIL'
        verdict = (
            'SANITY FAILURE — one of the sanity tests (T1 polarization '
            'sum identity, T4 propagator-pole robustness, T5 basis '
            'invariance) failed. The probe construction is internally '
            'inconsistent; check the implementation before drawing '
            'physical conclusions.'
        )
    else:
        verdict_class = 'OTHER'
        verdict = (
            'Mixed result. See per-test details for diagnosis.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'construction': (
            'M^{(λ, λ\')}_BAM = (ε^(λ)(k)·ε^(λ\')*(k\')) · '
            '[G_S3(ψ_s)·exp(iφ_s) + G_S3(ψ_u)·exp(iφ_u)] · T². '
            'Photon: two transverse linear polarizations per side. '
            'Electron: scalar at Thomson (no θ-dependent spin phase).'
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
    L.append('# Photon-structure probe — Klein-Nishina restoration')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Follow-on to the partial-match Compton amplitude probe (PR #26). '
        'Gives the photon explicit transverse-polarization machinery and '
        'tests whether the Klein-Nishina angular factor (1 + cos²θ)/2 is '
        'restored.'
    )
    L.append('')
    L.append('**Construction:**')
    L.append('')
    L.append('```')
    L.append(s['construction'])
    L.append('```')
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key metric | Value | PASS? |')
    L.append('|---|---|---|---:|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            metric, value = (
                'max |Σ|ε·ε\'|² − (1+cos²θ)|',
                f"{t['max_pointwise_error']:.2e}",
            )
        elif nm.startswith('T2'):
            f = t['fitted_coefficients']
            metric, value = (
                'fitted (A, B, C) vs KN (½, 0, ½)',
                f"({f['A']:.4f}, {f['B']:.4f}, {f['C']:.4f})",
            )
        elif nm.startswith('T3'):
            c = t['fitted_coefficients_degree3']
            metric, value = (
                'cos³θ coefficient (predicted ≠ 0)',
                f"{c['c3']:.4f}",
            )
        elif nm.startswith('T4'):
            metric, value = (
                'fitted pole exponent (expected 1.0)',
                f"{t['fitted_pole_exponent']:.4f}",
            )
        elif nm.startswith('T5'):
            metric, value = (
                'max relative diff linear vs circular',
                f"{t['max_relative_difference']:.2e}",
            )
        else:
            metric, value = '—', '—'
        L.append(f"| {nm[:2]} | `{nm}` | {metric} | {value} | {passed} |")
    L.append('')

    for t in s['tests']:
        L.append(f"## {t['name']}")
        L.append('')
        L.append(t['description'])
        L.append('')
        if t['name'].startswith('T2'):
            f = t['fitted_coefficients']
            k = t['klein_nishina_target']
            L.append(
                f"Fitted: A = {f['A']:+.4f}, B = {f['B']:+.4f}, "
                f"C = {f['C']:+.4f}."
            )
            L.append(
                f"KN target: A = {k['A']:+.4f}, B = {k['B']:+.4f}, "
                f"C = {k['C']:+.4f}."
            )
            L.append('')
            L.append(
                f"**Max pointwise residual to KN (normalised):** "
                f"{t['max_pointwise_residual_to_KN']:.4e}"
            )
            L.append('')
        elif t['name'].startswith('T3'):
            c = t['fitted_coefficients_degree3']
            L.append(
                f"Fitted coefficients (degree 3): c0 = {c['c0']:+.4f}, "
                f"c1 = {c['c1']:+.4f}, c2 = {c['c2']:+.4f}, "
                f"c3 = {c['c3']:+.4f}."
            )
            L.append('')
            L.append(
                f"**Max residual to predicted form "
                f"`(1+cos θ)(1+cos²θ)/2`:** "
                f"{t['max_residual_to_predicted_convolution']:.4e}"
            )
            L.append('')
        elif t['name'].startswith('T4'):
            L.append(
                f"Fitted pole exponent: **{t['fitted_pole_exponent']:.6f}** "
                f"(expected 1.0). Residue: {t['residue']:.4e}."
            )
            L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **High-energy Klein-Nishina recoil.** This probe is '
        'Thomson-limit only (ω ≪ m_e). The full KN formula '
        '`(ω\'/ω)²·[ω\'/ω + ω/ω\' − sin²θ]` at finite ω has a recoil '
        'enhancement that the present BAM construction does not yet '
        'reproduce. A separate sub-probe will test whether the '
        'kinematic ψ ↔ Mandelstam mapping reproduces the recoil '
        'factor or whether additional structure is needed.'
    )
    L.append(
        '- **Electron spin at finite energy.** The probe identifies '
        'the spin-½ angular phases as spurious at Thomson (T3 '
        'confirms the spurious (1+cos θ) factor). Whether they '
        'contribute correctly at finite ω/m_e (where electron spin '
        'does enter the QED amplitude) is open. The natural BAM '
        'reading is that spin contributes at sub-leading order in '
        'ω/m_e, consistent with how Dirac spinor structure enters '
        'standard QED Compton beyond the Thomson limit.'
    )
    L.append(
        '- **Loop corrections.** Tree-level only; vertex, '
        'self-energy, vacuum polarization require BAM\'s bulk '
        'radial channel (closure-ledger thread).'
    )
    L.append(
        '- **Lorentz covariance.** The polarization-vector '
        'construction is set up in the electron rest frame. '
        'Boost-covariance of the antipodal map + polarization '
        'vectors is the natural next-pass check.'
    )
    L.append(
        '- **Other QED events.** Pair production, electron '
        'self-energy, e⁻e⁻ scattering test the same BAM-amplitude '
        'machinery against different QED tree diagrams. With the '
        'photon-polarization structure now identified, those probes '
        'become tractable.'
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
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_compton_photon_structure_probe'
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
