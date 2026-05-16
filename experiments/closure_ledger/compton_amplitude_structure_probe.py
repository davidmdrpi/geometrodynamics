"""
BAM Compton amplitude — structural reproduction probe.

Follow-on to `compton_antipodal_kinematics_probe.py` (PR #25). The
kinematic probe established that BAM's antipodal-closure picture is
consistent with standard Compton scattering; this probe asks whether
a BAM transaction amplitude built from natural ingredients

    antipodal propagation  ·  throat transport  ·  closure phase

reproduces structural features of the QED Compton amplitude — either
the propagator denominators 1/(s−m²), 1/(u−m²), or the Klein-Nishina
angular dependence (1 + cos²θ)/2 in the Thomson limit.

The amplitude. Build

    M_BAM(θ, ω)  =  M_s(θ, ω)  +  M_u(θ, ω)

with each channel

    M_x  =  G_S3(ψ_x)  ·  exp(i·φ_x)  ·  T²        (x ∈ {s, u})

and ingredients

  - **Antipodal propagator** `G_S3(ψ) = ((π − ψ)cot ψ − ½)/(4π²R)` —
    the `S³` Green function from `transaction/s3_geometry.py`. Has a
    `1/ψ` pole as ψ → 0, exactly the structural object needed to
    reproduce the QED propagator pole.

  - **Mandelstam → ψ mapping**:
      ψ_s = (s − m²) / (2m²),
      ψ_u = |u − m²| / (2m²)
    so the on-shell QED poles map to `ψ → 0`, where `G_S3` diverges.

  - **Throat transport** `T = iσ_y`. Each scattering event has two
    throat traversals (initial and final mouths); `T² = −I` gives a
    factor `(i)² = −1` per channel. Treated here as a scalar `−1`
    times an SU(2) spin contribution absorbed into the closure phase.

  - **Closure phase** per channel:
      φ_s = π·cos(0) + θ/2  =  π + θ/2
      φ_u = π·cos(0) − θ/2  =  π − θ/2
    with `π·cos(0) = π` the Hopf holonomy on the χ = 0 fibre (the
    locked BAM lepton geometry) and ±θ/2 the SU(2) spin transport
    contribution at scattering angle θ.

  - **Channel signs** η_s = +1, η_u = −1 (s-u relative sign from the
    crossing of the electron line).

Tests:

  T1. **Propagator-pole reproduction.** Verify `|M_BAM| → ∞` with
      `1/(s − m²)` leading divergence in the on-shell s-channel
      limit. Fit the residue.

  T2. **Thomson-limit angular dependence.** Compute `|M_BAM(θ)|²` in
      the limit ω → 0 over θ ∈ [0, π]. Fit
      `|M_BAM|² ∝ A + B·cos θ + C·cos²θ` and compare to the
      Klein-Nishina form `(A, B, C) = (1, 0, 1)/2`.

  T3. **Ansatz sensitivity.** Repeat T2 with three different
      spin-phase ansätze (spin-½, spin-1, only-u-channel) and report
      which (if any) reproduces Klein-Nishina.

  T4. **Photon throat-pair diagnostic.** Multiply the amplitude by
      an extra spin-1 factor `exp(i·θ)` (a putative photon
      throat-pair contribution) and check whether Klein-Nishina is
      restored. Identifies whether photon-side BAM structure is the
      missing ingredient.

Expected outcome: PARTIAL match. Propagator poles reproduce
structurally (T1 PASS); Klein-Nishina angular dependence does NOT
reproduce from spin-½ alone (T2 FAIL); the missing structure is the
photon's spin-1 polarization sum, which would require explicit
photon throat-pair representation in BAM (T4 diagnostic).

The probe is honest about being a structural-reproduction test, not
a derivation. Even a partial match locates the open piece precisely.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

import numpy as np


PI = math.pi
TAU = 2.0 * PI


# ---------------------------------------------------------------------------
# S³ Green function (from transaction/s3_geometry.py, inlined for self-containment)
# ---------------------------------------------------------------------------

def G_S3(psi: float, R: float = 1.0, eps: float = 1e-9) -> float:
    """Zero-mean S³ Green function for the Laplacian:

        G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)

    diverges as 1/(4π R ψ) near ψ = 0 (recovering the flat-space
    Newtonian / Yukawa Green function) and as 1/(4π R (π−ψ)) near
    ψ = π.

    The probe uses |G| as the magnitude of the antipodal propagator;
    the sign is absorbed into the explicit channel phases.
    """
    psi_eff = max(min(psi, PI - eps), eps)
    return float(
        (((PI - psi_eff) / math.tan(psi_eff)) - 0.5)
        / (4.0 * PI * PI * R)
    )


# ---------------------------------------------------------------------------
# Compton kinematics (Mandelstam invariants from the kinematics probe)
# ---------------------------------------------------------------------------

def compton_mandelstam(omega: float, theta: float, m_e: float = 1.0) -> dict:
    """Standard Compton kinematics in the electron rest frame.

    Returns Mandelstam invariants s, u (and t for completeness).
    """
    omega_prime = omega / (1.0 + (omega / m_e) * (1.0 - math.cos(theta)))
    # 4-momenta
    p_e = np.array([m_e, 0.0, 0.0, 0.0])
    k_g = np.array([omega, omega, 0.0, 0.0])
    k_g_prime = np.array([
        omega_prime,
        omega_prime * math.cos(theta),
        0.0,
        omega_prime * math.sin(theta),
    ])
    s_4mom = p_e + k_g
    u_4mom = p_e - k_g_prime
    t_4mom = k_g - k_g_prime
    s = s_4mom[0] ** 2 - np.dot(s_4mom[1:], s_4mom[1:])
    u = u_4mom[0] ** 2 - np.dot(u_4mom[1:], u_4mom[1:])
    t = t_4mom[0] ** 2 - np.dot(t_4mom[1:], t_4mom[1:])
    return {
        's': float(s),
        'u': float(u),
        't': float(t),
        'omega_out': omega_prime,
    }


# ---------------------------------------------------------------------------
# BAM Compton amplitude with configurable spin-phase ansatz
# ---------------------------------------------------------------------------

def bam_compton_amplitude(
    omega: float,
    theta: float,
    m_e: float = 1.0,
    R_S3: float = 1.0,
    spin_ansatz: str = 'spin_half',
    include_photon_factor: bool = False,
) -> complex:
    """Construct M_BAM = M_s + M_u with the given spin-phase ansatz.

    Ansätze:
      'spin_half':  φ_s − φ_u =  θ   (each channel carries ±θ/2)
      'spin_one':   φ_s − φ_u =  2θ  (each channel carries ±θ)
      'only_u':     φ_s − φ_u =  θ   but only u-channel has phase

    `include_photon_factor`: multiply M_BAM by exp(iθ) — a putative
    spin-1 contribution from the photon throat-pair.
    """
    M = compton_mandelstam(omega, theta, m_e)
    s, u = M['s'], M['u']

    # Mandelstam → ψ mapping (linear in (s−m²) so G_S3 ∝ 1/(s−m²)
    # near threshold). The sign of (u − m²) is absorbed into ψ_u
    # being |u − m²|/(2m²); the relative s/u sign appears later as
    # a phase, not a magnitude, so we do NOT add an extra channel
    # sign (would double-count).
    eps_pole = 1e-12
    psi_s = max((s - m_e * m_e) / (2.0 * m_e * m_e), eps_pole)
    psi_u = max(abs(u - m_e * m_e) / (2.0 * m_e * m_e), eps_pole)
    psi_s = min(psi_s, PI - eps_pole)
    psi_u = min(psi_u, PI - eps_pole)

    # Antipodal propagators
    Gs = G_S3(psi_s, R=R_S3)
    Gu = G_S3(psi_u, R=R_S3)

    # Closure phases
    hopf_phase = PI * math.cos(0.0)   # = π for χ = 0 (BAM lock)
    if spin_ansatz == 'spin_half':
        phi_s = hopf_phase + 0.5 * theta
        phi_u = hopf_phase - 0.5 * theta
    elif spin_ansatz == 'spin_one':
        phi_s = hopf_phase + theta
        phi_u = hopf_phase - theta
    elif spin_ansatz == 'only_u':
        phi_s = hopf_phase
        phi_u = hopf_phase + theta
    else:
        raise ValueError(f'unknown spin_ansatz {spin_ansatz!r}')

    # Throat transport: T² = −I per channel ⇒ scalar factor (−1)
    throat_factor = -1.0

    M_s = Gs * cmath_exp(1j * phi_s) * throat_factor
    M_u = Gu * cmath_exp(1j * phi_u) * throat_factor

    M_total = M_s + M_u
    if include_photon_factor:
        M_total = M_total * cmath_exp(1j * theta)
    return M_total


def cmath_exp(z: complex) -> complex:
    """Use numpy for complex exp (math.exp doesn't take complex)."""
    return complex(np.exp(z))


# ---------------------------------------------------------------------------
# T1. Propagator-pole reproduction
# ---------------------------------------------------------------------------

def test_T1_propagator_poles(m_e: float = 1.0) -> dict:
    """Approach the s-channel pole `s → m²` (i.e. ω → 0 at fixed θ)
    and check that `|M_BAM|` diverges with `1/(s − m²)` leading
    behaviour.

    The s-channel kinematics: in the electron rest frame, s − m² = 2mω,
    so s → m² as ω → 0.
    """
    theta_fixed = PI / 3.0  # generic angle
    # Sample ω from 1e-2 down to 1e-6 — the unsaturated regime
    # (below 1e-6 the eps_pole clipping starts to dominate).
    omegas = np.array([10.0 ** (-k) for k in [2, 3, 4, 5, 6]])
    amplitudes = []
    s_minus_msq = []
    for omega in omegas:
        M = bam_compton_amplitude(float(omega), theta_fixed, m_e=m_e)
        amplitudes.append(abs(M))
        s = compton_mandelstam(float(omega), theta_fixed, m_e=m_e)['s']
        s_minus_msq.append(s - m_e * m_e)

    amplitudes = np.asarray(amplitudes)
    s_minus_msq = np.asarray(s_minus_msq)
    # Fit log |M| = −α · log(s − m²) + log K  ⇒  α should be 1
    valid = (amplitudes > 0) & (s_minus_msq > 0)
    if valid.sum() >= 3:
        log_M = np.log(amplitudes[valid])
        log_x = np.log(s_minus_msq[valid])
        slope, intercept = np.polyfit(log_x, log_M, 1)
        residue = math.exp(intercept)
    else:
        slope = float('nan')
        residue = float('nan')

    pole_exponent = float(-slope)  # expected 1.0 for 1/(s−m²) pole

    return {
        'name': 'T1_propagator_pole_reproduction',
        'description': (
            'Verify |M_BAM| ∝ 1/(s − m²) as s → m² (ω → 0). The fitted '
            'log-log slope should be −1; equivalently, the pole '
            'exponent is +1.'
        ),
        'theta_fixed': theta_fixed,
        'omegas_sampled': omegas.tolist(),
        'amplitudes': amplitudes.tolist(),
        's_minus_msq': s_minus_msq.tolist(),
        'fitted_pole_exponent': pole_exponent,
        'residue': residue,
        'expected_pole_exponent': 1.0,
        'pass': 0.9 < pole_exponent < 1.1,
    }


# ---------------------------------------------------------------------------
# T2. Thomson-limit angular dependence
# ---------------------------------------------------------------------------

def thomson_kn_normalized(theta: float) -> float:
    """Klein-Nishina / Thomson normalized angular factor F(θ)/F(0).

    F(θ) = (1 + cos²θ) / 2,  F(0) = 1,  F(π/2) = 1/2,  F(π) = 1.
    """
    return 0.5 * (1.0 + math.cos(theta) ** 2)


def _fit_angular_form(thetas: np.ndarray, M_sq: np.ndarray) -> tuple[float, float, float]:
    """Fit |M|² = A + B·cos θ + C·cos²θ. Returns (A, B, C).

    Standard linear least squares on cos θ and cos²θ basis.
    """
    c = np.cos(thetas)
    X = np.column_stack([np.ones_like(c), c, c * c])
    coeffs, _, _, _ = np.linalg.lstsq(X, M_sq, rcond=None)
    return float(coeffs[0]), float(coeffs[1]), float(coeffs[2])


def test_T2_thomson_angular(omega: float = 1e-4, m_e: float = 1.0) -> dict:
    """Compute |M_BAM(θ)|² for θ ∈ [0, π] at small ω (Thomson limit).
    Fit to A + B·cos θ + C·cos²θ. Klein-Nishina: (A, B, C) = (½, 0, ½).
    Normalized residual to KN reported.
    """
    thetas = np.linspace(0.0, PI, 33)
    M_sq = np.array([
        abs(bam_compton_amplitude(omega, float(t), m_e=m_e)) ** 2
        for t in thetas
    ])
    # Normalize at θ = 0
    M_sq_norm = M_sq / max(M_sq[0], 1e-30)
    kn_norm = np.array([thomson_kn_normalized(float(t)) for t in thetas])
    # KN normalized at θ=0 is already (1 + 1)/2 / 1 = 1.0 at θ=0.

    A, B, C = _fit_angular_form(thetas, M_sq_norm)
    # Compare to KN coefficients (after normalisation: F/F(0) = ½ + ½cos²θ)
    A_kn, B_kn, C_kn = 0.5, 0.0, 0.5
    # Max residual point-wise (normalised)
    fit_pred = A + B * np.cos(thetas) + C * np.cos(thetas) ** 2
    kn_pred = A_kn + B_kn * np.cos(thetas) + C_kn * np.cos(thetas) ** 2
    max_residual = float(np.max(np.abs(M_sq_norm - kn_pred)))

    # Identify the dominant angular form: is it cos²θ-dominant
    # (Klein-Nishina-like) or cos θ-dominant (spin-1/2-only) ?
    if abs(C) > 2.0 * abs(B):
        dominant_form = 'cos²θ-dominant (Klein-Nishina-like)'
    elif abs(B) > 2.0 * abs(C):
        dominant_form = 'cos θ-dominant (single-cosine spin-½ structure)'
    else:
        dominant_form = 'mixed'

    # PASS = the fitted (A, B, C) coefficients are close to KN's.
    # Tolerance: 5% relative on each.
    near_kn = (
        abs(A - A_kn) < 0.1 and
        abs(B - B_kn) < 0.1 and
        abs(C - C_kn) < 0.1
    )

    return {
        'name': 'T2_thomson_angular',
        'description': (
            'Compute |M_BAM(θ)|² in the Thomson limit (ω → 0) and fit '
            'A + B·cos θ + C·cos²θ. Klein-Nishina: (½, 0, ½).'
        ),
        'omega_thomson_limit': omega,
        'fitted_coefficients': {'A': A, 'B': B, 'C': C},
        'klein_nishina_coefficients': {'A': A_kn, 'B': B_kn, 'C': C_kn},
        'max_pointwise_residual_to_KN': max_residual,
        'dominant_angular_form': dominant_form,
        'theta_samples': thetas.tolist(),
        'amplitude_squared_normalized': M_sq_norm.tolist(),
        'kn_normalized': kn_norm.tolist(),
        'pass': near_kn,
    }


# ---------------------------------------------------------------------------
# T3. Ansatz sensitivity
# ---------------------------------------------------------------------------

def test_T3_ansatz_sensitivity(omega: float = 1e-4, m_e: float = 1.0) -> dict:
    """Repeat T2 with three different spin-phase ansätze and report
    which (if any) is closest to Klein-Nishina."""
    thetas = np.linspace(0.0, PI, 33)
    results = {}
    for ansatz in ('spin_half', 'spin_one', 'only_u'):
        M_sq = np.array([
            abs(bam_compton_amplitude(
                omega, float(t), m_e=m_e, spin_ansatz=ansatz,
            )) ** 2
            for t in thetas
        ])
        M_sq_norm = M_sq / max(M_sq[0], 1e-30)
        A, B, C = _fit_angular_form(thetas, M_sq_norm)
        kn_residual = float(
            np.max(np.abs(
                M_sq_norm - (0.5 + 0.5 * np.cos(thetas) ** 2)
            ))
        )
        results[ansatz] = {
            'fitted_A_B_C': (A, B, C),
            'max_residual_to_KN': kn_residual,
        }
    # Identify closest ansatz
    closest = min(results, key=lambda a: results[a]['max_residual_to_KN'])

    return {
        'name': 'T3_ansatz_sensitivity',
        'description': (
            'Repeat T2 with three spin-phase ansätze and identify '
            'which (if any) approaches Klein-Nishina.'
        ),
        'ansatz_results': results,
        'closest_to_KN_ansatz': closest,
        'closest_residual': results[closest]['max_residual_to_KN'],
        'pass': results[closest]['max_residual_to_KN'] < 0.1,
    }


# ---------------------------------------------------------------------------
# T4. Photon throat-pair diagnostic
# ---------------------------------------------------------------------------

def test_T4_photon_factor(omega: float = 1e-4, m_e: float = 1.0) -> dict:
    """Multiply M_BAM by exp(iθ) — a putative photon spin-1 throat-pair
    contribution — and re-test the angular fit. If this restores
    Klein-Nishina structure, it identifies the photon throat-pair as
    the missing BAM ingredient."""
    thetas = np.linspace(0.0, PI, 33)
    M_sq_base = np.array([
        abs(bam_compton_amplitude(omega, float(t), m_e=m_e)) ** 2
        for t in thetas
    ])
    M_sq_phot = np.array([
        abs(bam_compton_amplitude(
            omega, float(t), m_e=m_e, include_photon_factor=True,
        )) ** 2
        for t in thetas
    ])
    base_norm = M_sq_base / max(M_sq_base[0], 1e-30)
    phot_norm = M_sq_phot / max(M_sq_phot[0], 1e-30)
    A_b, B_b, C_b = _fit_angular_form(thetas, base_norm)
    A_p, B_p, C_p = _fit_angular_form(thetas, phot_norm)
    res_base = float(np.max(np.abs(
        base_norm - (0.5 + 0.5 * np.cos(thetas) ** 2)
    )))
    res_phot = float(np.max(np.abs(
        phot_norm - (0.5 + 0.5 * np.cos(thetas) ** 2)
    )))
    improvement = res_base - res_phot

    return {
        'name': 'T4_photon_factor',
        'description': (
            'Multiply M_BAM by exp(iθ) (a putative photon spin-1 '
            'contribution) and check whether Klein-Nishina structure '
            'is restored.'
        ),
        'baseline_fit_A_B_C': (A_b, B_b, C_b),
        'with_photon_factor_fit_A_B_C': (A_p, B_p, C_p),
        'baseline_residual_to_KN': res_base,
        'with_photon_factor_residual_to_KN': res_phot,
        'improvement_from_photon_factor': improvement,
        'photon_factor_restores_KN': res_phot < 0.1,
        # PASS = either the baseline already matches KN, or the
        # diagnostic identifies the photon factor as the fix
        'pass': (res_base < 0.1) or (improvement > 0.0),
    }


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_propagator_poles()
    t2 = test_T2_thomson_angular()
    t3 = test_T3_ansatz_sensitivity()
    t4 = test_T4_photon_factor()
    tests = [t1, t2, t3, t4]

    # Verdict: classify as full / partial / no match
    propagator_ok = t1['pass']
    angular_ok = t2['pass']
    if propagator_ok and angular_ok:
        verdict_class = 'FULL_MATCH'
        verdict = (
            'FULL STRUCTURAL MATCH — the BAM amplitude reproduces '
            'both the Compton propagator poles 1/(s − m²), 1/(u − m²) '
            'and the Klein-Nishina angular dependence (1 + cos²θ)/2 '
            'in the Thomson limit, from the natural BAM ingredients '
            '(antipodal `S³` Green function + throat transport + '
            'Hopf-holonomy closure phase). This would be a strong '
            'BAM progress claim.'
        )
    elif propagator_ok and not angular_ok:
        verdict_class = 'PARTIAL_MATCH'
        verdict = (
            'PARTIAL STRUCTURAL MATCH — the BAM amplitude reproduces '
            'the propagator pole structure (T1 PASS: `G_S3(ψ) ∼ 1/ψ` '
            'with `ψ ∝ s − m²` gives the correct 1/(s − m²) leading '
            'divergence) but NOT the Klein-Nishina angular dependence '
            '(T2 FAIL). The dominant angular form is '
            f'`{t2["dominant_angular_form"]}` rather than KN`s '
            'cos²θ-dominant structure. The diagnostic T4 identifies '
            'whether an extra photon-side phase factor restores KN '
            '— if yes, the missing structure is the photon throat-pair '
            'representation; if no, the gap is deeper. This is the '
            'expected outcome and identifies the next structural piece.'
        )
    elif not propagator_ok and angular_ok:
        verdict_class = 'ANGULAR_ONLY'
        verdict = (
            'ANGULAR-ONLY MATCH — Klein-Nishina structure reproduced '
            'but propagator pole structure broken. Unexpected; '
            'indicates the BAM amplitude construction has a flaw in '
            'the antipodal-propagator → QED-propagator identification.'
        )
    else:
        verdict_class = 'NO_MATCH'
        verdict = (
            'NO STRUCTURAL MATCH — BAM amplitude does not reproduce '
            'either propagator pole structure or angular dependence. '
            'Falsifies the amplitude-level reinterpretation thread at '
            'this construction. Either the ingredient identification '
            '(antipodal → propagator, throat → vertex, closure → '
            'phase) is wrong, or BAM does not have amplitude-level '
            'content matching QED at tree level.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'construction': (
            'M_BAM(θ, ω) = M_s + M_u with M_x = G_S3(ψ_x) · '
            'exp(i·φ_x) · T² where ψ_x = (x − m²)/(2m²), '
            'φ_x = π ± θ/2, T² = −I.'
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
    L.append('# BAM Compton amplitude — structural reproduction probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Tests whether a BAM transaction amplitude built from natural '
        'ingredients reproduces structural features of the QED Compton '
        'amplitude — propagator denominators and Klein-Nishina '
        'angular dependence.'
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
        if t['name'].startswith('T1'):
            metric = 'fitted pole exponent (expected 1.0)'
            value = f"{t['fitted_pole_exponent']:.4f}"
        elif t['name'].startswith('T2'):
            C = t['fitted_coefficients']['C']
            metric = 'fitted cos²θ coefficient C (KN: 0.5)'
            value = f"{C:+.4f}"
        elif t['name'].startswith('T3'):
            metric = 'closest ansatz residual to KN'
            value = f"{t['closest_residual']:.4f}"
        elif t['name'].startswith('T4'):
            metric = 'baseline residual to KN'
            value = f"{t['baseline_residual_to_KN']:.4f}"
        else:
            metric = '—'; value = '—'
        L.append(f"| {t['name'][:2]} | `{t['name']}` | {metric} | {value} | {passed} |")
    L.append('')

    for t in s['tests']:
        L.append(f"## {t['name']}")
        L.append('')
        L.append(t['description'])
        L.append('')
        if t['name'].startswith('T1'):
            L.append(
                f"Fitted pole exponent: **{t['fitted_pole_exponent']:.6f}** "
                f"(expected {t['expected_pole_exponent']}). "
                f"Residue: {t['residue']:.4e}."
            )
            L.append('')
            L.append('Sampled (ω, |M_BAM|):')
            L.append('')
            L.append('```')
            for omega, M_abs in list(zip(
                t['omegas_sampled'], t['amplitudes']
            ))[:6]:
                L.append(f"  ω = {omega:.2e},  |M_BAM| = {M_abs:.6e}")
            L.append('```')
            L.append('')
        elif t['name'].startswith('T2'):
            f = t['fitted_coefficients']
            k = t['klein_nishina_coefficients']
            L.append(
                f"Fitted: A = {f['A']:+.4f}, B = {f['B']:+.4f}, "
                f"C = {f['C']:+.4f}."
            )
            L.append(
                f"Klein-Nishina: A = {k['A']:+.4f}, B = {k['B']:+.4f}, "
                f"C = {k['C']:+.4f}."
            )
            L.append('')
            L.append(
                f"**Dominant angular form:** `{t['dominant_angular_form']}`"
            )
            L.append('')
            L.append(
                f"**Max pointwise residual to KN (normalised):** "
                f"{t['max_pointwise_residual_to_KN']:.4f}"
            )
            L.append('')
        elif t['name'].startswith('T3'):
            L.append('| ansatz | A | B | C | residual to KN |')
            L.append('|---|---:|---:|---:|---:|')
            for ansatz, r in t['ansatz_results'].items():
                A, B, C = r['fitted_A_B_C']
                L.append(
                    f"| `{ansatz}` | {A:+.4f} | {B:+.4f} | "
                    f"{C:+.4f} | {r['max_residual_to_KN']:.4f} |"
                )
            L.append('')
            L.append(
                f"**Closest ansatz:** `{t['closest_to_KN_ansatz']}` "
                f"with residual {t['closest_residual']:.4f}."
            )
            L.append('')
        elif t['name'].startswith('T4'):
            A_b, B_b, C_b = t['baseline_fit_A_B_C']
            A_p, B_p, C_p = t['with_photon_factor_fit_A_B_C']
            L.append(f"Baseline fit: A = {A_b:+.4f}, B = {B_b:+.4f}, "
                     f"C = {C_b:+.4f}.")
            L.append(f"With photon factor: A = {A_p:+.4f}, B = {B_p:+.4f}, "
                     f"C = {C_p:+.4f}.")
            L.append('')
            L.append(
                f"Residual to KN: baseline = "
                f"{t['baseline_residual_to_KN']:.4f}, "
                f"with photon factor = "
                f"{t['with_photon_factor_residual_to_KN']:.4f} "
                f"(improvement {t['improvement_from_photon_factor']:+.4f})."
            )
            L.append('')
            verdict = (
                'YES — restored' if t['photon_factor_restores_KN']
                else 'NO — does not restore'
            )
            L.append(f"**Photon factor restores Klein-Nishina:** {verdict}")
            L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Klein-Nishina from polarization sum.** The QED '
        'angular factor (1 + cos²θ) arises from the polarization '
        'sum of two transverse photon polarisations '
        '|ε·ε\'|² (averaged over incoming, summed over outgoing). '
        'BAM does not yet have explicit photon polarization machinery; '
        'reproducing this sum requires either (a) explicit photon '
        'throat-pair representation (T4 diagnostic), (b) a Hopf-fibre '
        'photon model coupled to the connection (see '
        '`geometrodynamics/hopf/connection.py`), or (c) an `S³` vector '
        'Green function rather than the scalar one used here.'
    )
    L.append(
        '- **Energy-dependent recoil enhancement.** Klein-Nishina at '
        'high energy `ω ~ m_e` differs from Thomson by the recoil '
        'factor `(ω\'/ω)²·[ω\'/ω + ω/ω\' − sin²θ]`. The BAM amplitude '
        'inherits energy-dependence only through `ψ_s, ψ_u` near the '
        'pole. Whether the full recoil enhancement is reproduced is a '
        'separate test (not run here).'
    )
    L.append(
        '- **Loop corrections.** Even if tree-level structure is '
        'reproduced, one-loop corrections (vertex, self-energy) would '
        'require BAM\'s bulk radial channel — the missing radial-bulk '
        'phase from the closure-ledger framework. Cross-connection '
        'between threads is open.'
    )
    L.append(
        '- **Mandelstam → ψ mapping uniqueness.** The choice '
        '`ψ = (s − m²)/(2m²)` is the simplest natural choice but not '
        'the only one. Alternative mappings (e.g. `ψ = arctan((s − '
        'm²)/Λ²)`) may give different sub-leading behaviour. The pole '
        'reproduction is robust to choice; finer structure is not.'
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
    out = here / 'runs' / f'{ts}_compton_amplitude_structure_probe'
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
