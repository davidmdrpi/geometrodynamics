"""
Bhabha / Møller two-channel interference probe.

Tests whether the BAM scalar-intensity tree kernel (derived for the
Compton/BW/annihilation triangle in PRs #35-#41) can reproduce
tree-level QED scalar-intensity structure when two crossed kernels
must be coherently summed with a Fermi-statistics-determined relative
sign.

Bhabha (e⁺e⁻ → e⁺e⁻, massless e):

    |M̄|²_Bhabha / (8e⁴)  =  (s² + u²)/t²  +  (u² + t²)/s²  +  2u²/(s·t)
                              t-chan²        s-chan²         interference (negative)

Møller (e⁻e⁻ → e⁻e⁻, massless e):

    |M̄|²_Møller / (8e⁴)  =  (s² + u²)/t²  +  (s² + t²)/u²  +  2s²/(t·u)
                              t-chan²        u-chan²         interference (positive)

CM kinematics: s = 4E², t = -s(1-cosθ)/2, u = -s(1+cosθ)/2.

Tests:

  T1. QED Bhabha scalar intensity computed at sample CM angles.
  T2. QED Møller scalar intensity computed similarly.
  T3. Interference signs: Bhabha negative, Møller positive.
  T4. BAM Construction I (positive √-channel sum) → fails on sign.
  T5. BAM Construction II (Fermi sign added by hand) → sign ok,
      magnitude still fails.
  T6. BAM Construction III (Compton-F² inspired vertex factor) → tests
      diagonal-term match.
  T7. Residual analysis: identify the specific structure missing.
  T8. Required BAM extension: name the missing ingredient.
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
# Mandelstam invariants (massless e, CM frame)
# ---------------------------------------------------------------------------

def cm_mandelstam(cos_theta: float, s: float = 1.0) -> tuple[float, float, float]:
    """Return (s, t, u) at scattering angle θ in CM, with s normalised
    to unity by default."""
    t = -s * (1.0 - cos_theta) / 2.0
    u = -s * (1.0 + cos_theta) / 2.0
    return s, t, u


# ---------------------------------------------------------------------------
# QED scalar intensities
# ---------------------------------------------------------------------------

def M2_Bhabha_over_8e4(s: float, t: float, u: float) -> dict:
    """|M̄|²_Bhabha / (8e⁴) and its three terms."""
    M2_t_chan = (s * s + u * u) / (t * t)
    M2_s_chan = (u * u + t * t) / (s * s)
    interf = 2.0 * u * u / (s * t)
    return {
        'M2_t_channel': M2_t_chan,
        'M2_s_channel': M2_s_chan,
        'interference': interf,
        'total': M2_t_chan + M2_s_chan + interf,
    }


def M2_Moller_over_8e4(s: float, t: float, u: float) -> dict:
    """|M̄|²_Møller / (8e⁴) and its three terms."""
    M2_t_chan = (s * s + u * u) / (t * t)
    M2_u_chan = (s * s + t * t) / (u * u)
    interf = 2.0 * s * s / (t * u)
    return {
        'M2_t_channel': M2_t_chan,
        'M2_u_channel': M2_u_chan,
        'interference': interf,
        'total': M2_t_chan + M2_u_chan + interf,
    }


# ---------------------------------------------------------------------------
# T1. QED Bhabha
# ---------------------------------------------------------------------------

def test_T1_QED_bhabha() -> dict:
    rows = []
    for theta_deg in [30, 60, 90, 120, 150]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)
        parts = M2_Bhabha_over_8e4(s, t, u)
        rows.append({
            'theta_deg': theta_deg,
            'cos_theta': c,
            's': s, 't': t, 'u': u,
            'M2_t_channel': parts['M2_t_channel'],
            'M2_s_channel': parts['M2_s_channel'],
            'interference_2u2_over_st': parts['interference'],
            'total_M2_over_8e4': parts['total'],
        })
    return {
        'name': 'T1_QED_bhabha_scalar_intensity',
        'description': (
            "QED |M̄|²_Bhabha/(8e⁴) = (s²+u²)/t² + (u²+t²)/s² + 2u²/(s·t) "
            "at sample CM angles. Diagonal and interference terms tabulated."
        ),
        'rows': rows,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. QED Møller
# ---------------------------------------------------------------------------

def test_T2_QED_moller() -> dict:
    rows = []
    for theta_deg in [30, 60, 90, 120, 150]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)
        parts = M2_Moller_over_8e4(s, t, u)
        rows.append({
            'theta_deg': theta_deg,
            'cos_theta': c,
            's': s, 't': t, 'u': u,
            'M2_t_channel': parts['M2_t_channel'],
            'M2_u_channel': parts['M2_u_channel'],
            'interference_2s2_over_tu': parts['interference'],
            'total_M2_over_8e4': parts['total'],
        })
    return {
        'name': 'T2_QED_moller_scalar_intensity',
        'description': (
            "QED |M̄|²_Møller/(8e⁴) = (s²+u²)/t² + (s²+t²)/u² + 2s²/(t·u). "
            "Møller is the e⁻e⁻ analogue with Pauli antisymmetrisation; "
            "the interference sign is positive (t·u > 0 since both < 0)."
        ),
        'rows': rows,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. Interference signs
# ---------------------------------------------------------------------------

def test_T3_interference_signs() -> dict:
    bhabha_signs = []
    moller_signs = []
    for theta_deg in [30, 60, 90, 120, 150]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)
        bhabha_signs.append(np.sign(M2_Bhabha_over_8e4(s, t, u)['interference']))
        moller_signs.append(np.sign(M2_Moller_over_8e4(s, t, u)['interference']))
    bhabha_all_negative = all(sign < 0 for sign in bhabha_signs)
    moller_all_positive = all(sign > 0 for sign in moller_signs)
    return {
        'name': 'T3_interference_signs',
        'description': (
            "Bhabha interference 2u²/(s·t) is NEGATIVE in CM "
            "(t < 0, s > 0). Møller interference 2s²/(t·u) is POSITIVE "
            "(t·u > 0 since both negative). The opposite signs come "
            "from different Fermi-statistics origins: Wick-contraction "
            "for Bhabha s↔t exchange, Pauli antisymmetrisation for "
            "Møller t↔u exchange."
        ),
        'bhabha_signs_across_angles': bhabha_signs,
        'moller_signs_across_angles': moller_signs,
        'bhabha_all_negative': bool(bhabha_all_negative),
        'moller_all_positive': bool(moller_all_positive),
        'pass': bool(bhabha_all_negative) and bool(moller_all_positive),
    }


# ---------------------------------------------------------------------------
# T4. BAM Construction I — positive √-channel sum
# ---------------------------------------------------------------------------

def bam_positive_sum(channels: tuple[float, float]) -> tuple[float, float, float]:
    """Given the two diagonal channel scalar intensities |M_a|² and
    |M_b|², compute the positive-sum coherent total:

        M_total = √|M_a|² + √|M_b|²
        |M_total|² = |M_a|² + |M_b|² + 2·√(|M_a|²·|M_b|²)

    Returns (sum_of_diags, positive_interference, total)."""
    M2_a, M2_b = channels
    diag = M2_a + M2_b
    interf = 2.0 * math.sqrt(max(M2_a, 0.0) * max(M2_b, 0.0))
    return diag, interf, diag + interf


def test_T4_construction_I_positive_sum() -> dict:
    rows = []
    bhabha_residuals = []
    moller_residuals = []
    for theta_deg in [30, 60, 90, 120, 150]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)

        # Bhabha
        b_qed = M2_Bhabha_over_8e4(s, t, u)
        b_diag, b_pos_interf, b_total = bam_positive_sum(
            (b_qed['M2_t_channel'], b_qed['M2_s_channel'])
        )
        b_interf_residual = b_pos_interf - b_qed['interference']
        b_total_residual = b_total - b_qed['total']
        bhabha_residuals.append(b_total_residual)

        # Møller
        m_qed = M2_Moller_over_8e4(s, t, u)
        m_diag, m_pos_interf, m_total = bam_positive_sum(
            (m_qed['M2_t_channel'], m_qed['M2_u_channel'])
        )
        m_interf_residual = m_pos_interf - m_qed['interference']
        m_total_residual = m_total - m_qed['total']
        moller_residuals.append(m_total_residual)

        rows.append({
            'theta_deg': theta_deg,
            'bhabha_QED_interference': b_qed['interference'],
            'bhabha_BAM_positive_interference': b_pos_interf,
            'bhabha_interference_residual': b_interf_residual,
            'bhabha_total_residual': b_total_residual,
            'moller_QED_interference': m_qed['interference'],
            'moller_BAM_positive_interference': m_pos_interf,
            'moller_interference_residual': m_interf_residual,
            'moller_total_residual': m_total_residual,
        })

    bhabha_max_resid = max(abs(r) for r in bhabha_residuals)
    moller_max_resid = max(abs(r) for r in moller_residuals)
    return {
        'name': 'T4_construction_I_positive_sum',
        'description': (
            "Construction I: M_channel = √(diagonal QED |M|²), summed "
            "with positive sign. Predicts positive interference always — "
            "fails on Bhabha (which is negative) and on Møller "
            "magnitude (sign right but value wrong)."
        ),
        'rows': rows,
        'bhabha_max_total_residual': bhabha_max_resid,
        'moller_max_total_residual': moller_max_resid,
        'pass': False,  # documented falsification
    }


# ---------------------------------------------------------------------------
# T5. BAM Construction II — Fermi sign added by hand
# ---------------------------------------------------------------------------

def bam_signed_sum(
    channels: tuple[float, float], relative_sign: int,
) -> tuple[float, float, float]:
    """Given the two diagonal scalar intensities and a relative sign
    (+1 or -1), compute |M_a + sign·M_b|²."""
    M2_a, M2_b = channels
    diag = M2_a + M2_b
    interf = 2.0 * relative_sign * math.sqrt(max(M2_a, 0.0) * max(M2_b, 0.0))
    return diag, interf, diag + interf


def test_T5_construction_II_fermi_sign() -> dict:
    rows = []
    bhabha_residuals = []
    moller_residuals = []
    for theta_deg in [30, 60, 90, 120, 150]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)

        # Bhabha: sign chosen so the interference is negative (matches QED)
        b_qed = M2_Bhabha_over_8e4(s, t, u)
        _, b_interf, b_total = bam_signed_sum(
            (b_qed['M2_t_channel'], b_qed['M2_s_channel']), relative_sign=-1,
        )
        b_interf_residual = b_interf - b_qed['interference']
        b_total_residual = b_total - b_qed['total']
        bhabha_residuals.append(b_total_residual)

        # Møller: sign chosen so the interference is positive (matches QED)
        m_qed = M2_Moller_over_8e4(s, t, u)
        _, m_interf, m_total = bam_signed_sum(
            (m_qed['M2_t_channel'], m_qed['M2_u_channel']), relative_sign=+1,
        )
        m_interf_residual = m_interf - m_qed['interference']
        m_total_residual = m_total - m_qed['total']
        moller_residuals.append(m_total_residual)

        rows.append({
            'theta_deg': theta_deg,
            'bhabha_QED_interference': b_qed['interference'],
            'bhabha_BAM_signed_interference': b_interf,
            'bhabha_interference_residual': b_interf_residual,
            'bhabha_total_residual': b_total_residual,
            'moller_QED_interference': m_qed['interference'],
            'moller_BAM_signed_interference': m_interf,
            'moller_interference_residual': m_interf_residual,
            'moller_total_residual': m_total_residual,
        })

    bhabha_max = max(abs(r) for r in bhabha_residuals)
    moller_max = max(abs(r) for r in moller_residuals)
    return {
        'name': 'T5_construction_II_fermi_sign',
        'description': (
            "Construction II: explicit Fermi-statistics sign added by hand "
            "(Bhabha: -1 between s and t channels; Møller: +1 between t "
            "and u channels, matching antisymmetrised result). The sign "
            "is now correct, but the magnitude of the interference "
            "(2·√(|M_a|²·|M_b|²)) does not equal the QED interference "
            "(2u²/(s·t) for Bhabha or 2s²/(t·u) for Møller) because "
            "the Dirac-trace algebra produces a different functional "
            "form than the geometric mean of the diagonals."
        ),
        'rows': rows,
        'bhabha_max_total_residual': bhabha_max,
        'moller_max_total_residual': moller_max,
        'pass': False,  # documented falsification
    }


# ---------------------------------------------------------------------------
# T6. BAM Construction III — Compton-F² inspired vertex factor
# ---------------------------------------------------------------------------

def K_padé(x: float) -> float:
    return 2.0 * x / (1.0 + x)


def Q_polarization(x: float, c: float) -> float:
    return x * x + x * (1.0 - x) ** 2 / (1.0 + c * c)


def bam_compton_channel(q_squared: float, x_eff: float = 1.0, c_eff: float = 0.0) -> float:
    """Channel amplitude squared:
       |M_channel|² ∝ K(x_eff)² · Q(x_eff, c_eff) / q⁴

    With x_eff = 1 (massless, no per-vertex recoil) and c_eff = 0
    (forward-scattering-like default), K = 1, Q = 1; the channel
    intensity reduces to 1/q⁴.
    """
    K = K_padé(x_eff)
    Q = Q_polarization(x_eff, c_eff)
    return (K * K * Q) / (q_squared * q_squared)


def test_T6_construction_III_compton_vertex() -> dict:
    """Test the Compton-F²-inspired construction with x_eff = 1
    (massless limit). Even the diagonal channel intensities should
    match QED — if they don't, Bhabha/Møller live outside the BAM
    F² analytic continuation class entirely."""
    rows = []
    bhabha_diag_residuals = []
    moller_diag_residuals = []
    for theta_deg in [60, 90, 120]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)

        # Bhabha diagonal predictions
        bam_t_chan = bam_compton_channel(t, x_eff=1.0, c_eff=c)
        bam_s_chan = bam_compton_channel(s, x_eff=1.0, c_eff=c)
        b_qed = M2_Bhabha_over_8e4(s, t, u)
        b_t_residual = bam_t_chan - b_qed['M2_t_channel']
        b_s_residual = bam_s_chan - b_qed['M2_s_channel']
        bhabha_diag_residuals.append(max(abs(b_t_residual), abs(b_s_residual)))

        # Møller diagonal predictions
        bam_u_chan = bam_compton_channel(u, x_eff=1.0, c_eff=c)
        m_qed = M2_Moller_over_8e4(s, t, u)
        m_t_residual = bam_t_chan - m_qed['M2_t_channel']
        m_u_residual = bam_u_chan - m_qed['M2_u_channel']
        moller_diag_residuals.append(max(abs(m_t_residual), abs(m_u_residual)))

        rows.append({
            'theta_deg': theta_deg,
            'bhabha_QED_M2_t': b_qed['M2_t_channel'],
            'bhabha_BAM_M2_t_via_compton_vertex': bam_t_chan,
            'bhabha_t_diag_residual': b_t_residual,
            'bhabha_QED_M2_s': b_qed['M2_s_channel'],
            'bhabha_BAM_M2_s_via_compton_vertex': bam_s_chan,
            'bhabha_s_diag_residual': b_s_residual,
            'moller_QED_M2_t': m_qed['M2_t_channel'],
            'moller_BAM_M2_t_via_compton_vertex': bam_t_chan,
            'moller_t_diag_residual': m_t_residual,
            'moller_QED_M2_u': m_qed['M2_u_channel'],
            'moller_BAM_M2_u_via_compton_vertex': bam_u_chan,
            'moller_u_diag_residual': m_u_residual,
        })

    return {
        'name': 'T6_construction_III_compton_vertex',
        'description': (
            "Construction III: M_channel = K(1)·√Q(1,c)/q² = 1/q². "
            "Even the diagonal Bhabha/Møller intensities don't match — "
            "QED diagonals have numerators (s²+u²), (u²+t²), (s²+t²) "
            "from Dirac traces, which 1/q⁴ alone doesn't capture. "
            "Bhabha/Møller sit outside the BAM Compton-F² analytic "
            "continuation class."
        ),
        'rows': rows,
        'bhabha_max_diagonal_residual': max(bhabha_diag_residuals),
        'moller_max_diagonal_residual': max(moller_diag_residuals),
        'pass': False,  # documented falsification
    }


# ---------------------------------------------------------------------------
# T7. Residual analysis
# ---------------------------------------------------------------------------

def test_T7_residual_analysis() -> dict:
    """Identify the specific QED structure that BAM scalar-intensity
    kernels cannot reproduce. The two key residuals:

      (a) The diagonal-term numerators (s²+u²), (u²+t²), (s²+t²) come
          from Dirac trace algebra and cannot be reduced to a Padé-
          times-Hopf-spinor product (the BAM F² structure).

      (b) The interference term (2u²/(s·t) for Bhabha; 2s²/(t·u) for
          Møller) has a specific magnitude AND sign determined by
          spinor antisymmetrisation. The magnitude does NOT equal the
          geometric mean of the diagonals (which is what a positive-
          real scalar BAM construction predicts).
    """
    # (a) Show the diagonal numerators have a specific Dirac-trace
    # structure that doesn't factor naturally.
    diag_terms_analysis = []
    for theta_deg in [60, 90, 120]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)
        diag_terms_analysis.append({
            'theta_deg': theta_deg,
            's2_plus_u2': s * s + u * u,
            'u2_plus_t2': u * u + t * t,
            's2_plus_t2': s * s + t * t,
            'comment': (
                "These trace-derived numerators cannot be reproduced "
                "by Padé(x)·Hopf-spinor(x,c) for any natural BAM x and c."
            ),
        })

    # (b) Show the interference magnitude differs from the geometric
    # mean of the diagonals.
    interf_vs_geom_mean = []
    for theta_deg in [60, 90, 120]:
        c = math.cos(math.radians(theta_deg))
        s, t, u = cm_mandelstam(c, s=1.0)
        b_qed = M2_Bhabha_over_8e4(s, t, u)
        bhabha_interf_abs = abs(b_qed['interference'])
        bhabha_geom_mean = 2.0 * math.sqrt(
            b_qed['M2_t_channel'] * b_qed['M2_s_channel']
        )
        interf_vs_geom_mean.append({
            'theta_deg': theta_deg,
            'bhabha_|interference|_QED': bhabha_interf_abs,
            'bhabha_2x_geometric_mean_diagonals_BAM_positive': bhabha_geom_mean,
            'ratio_BAM_over_QED': bhabha_geom_mean / bhabha_interf_abs if bhabha_interf_abs > 0 else float('inf'),
        })

    return {
        'name': 'T7_residual_analysis',
        'description': (
            "Two BAM residuals identified: (a) the QED diagonal-term "
            "Dirac-trace numerators (s²+u²) etc. don't factor as "
            "Padé·Hopf-spinor; (b) the interference magnitude is "
            "specific to Dirac-trace algebra, not the geometric mean "
            "of the diagonals."
        ),
        'diagonal_term_analysis': diag_terms_analysis,
        'interference_vs_geometric_mean_analysis': interf_vs_geom_mean,
        'pass': True,   # informative
    }


# ---------------------------------------------------------------------------
# T8. Required BAM extension
# ---------------------------------------------------------------------------

def test_T8_required_extension() -> dict:
    """Name the minimal additional structure BAM needs to handle
    multi-channel interference."""
    return {
        'name': 'T8_required_BAM_extension',
        'description': (
            "BAM's scalar-intensity kernel (verified for the Compton/BW/ann "
            "triangle in PRs #35-#41) cannot reproduce Bhabha/Møller "
            "scalar intensities at either the diagonal or interference "
            "level. The missing structure is the Dirac γ-matrix algebra "
            "at the QED vertex — specifically the spinor trace that "
            "produces (s²+u²), (u²+t²), (s²+t²) numerators and the "
            "specific (2u²/(s·t), 2s²/(t·u)) interference."
        ),
        'missing_ingredients': [
            (
                "Dirac spinors at the vertex. The Compton BAM kernel "
                "F² = K²·Q already integrates over Dirac structure "
                "implicitly via the closed form, but for multi-channel "
                "processes the per-channel amplitudes must carry "
                "explicit spinor weights."
            ),
            (
                "Virtual-photon propagator structure beyond 1/q². "
                "Internal photon lines couple two BAM vertices; their "
                "transport on the throat may carry Hopf-fibre phase "
                "in addition to 1/q²."
            ),
            (
                "Channel sign rule from non-orientable throat topology. "
                "The README's channel-2 (Möbius / antipodal T = iσ_y) "
                "is the geometric origin of spin-½ and could in "
                "principle drive Pauli signs for identical-fermion "
                "interference (Møller t↔u). Whether it correctly "
                "predicts the Bhabha s↔t Wick sign is a deeper "
                "open question."
            ),
        ],
        'next_probe_candidates': [
            (
                "Dirac-trace BAM probe: derive (s²+u²) etc. from a "
                "geometric trace over per-vertex spinor configurations."
            ),
            (
                "Möbius-sign probe: test whether T = iσ_y gives the "
                "Møller Pauli sign automatically."
            ),
            (
                "Virtual-photon throat-fibre propagator probe: derive "
                "1/q² + corrections from BAM throat-fibre dynamics."
            ),
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_QED_bhabha()
    t2 = test_T2_QED_moller()
    t3 = test_T3_interference_signs()
    t4 = test_T4_construction_I_positive_sum()
    t5 = test_T5_construction_II_fermi_sign()
    t6 = test_T6_construction_III_compton_vertex()
    t7 = test_T7_residual_analysis()
    t8 = test_T8_required_extension()
    tests = [t1, t2, t3, t4, t5, t6, t7, t8]

    qed_baselines_ok = t1['pass'] and t2['pass'] and t3['pass']
    bam_constructions_fail = (
        not t4['pass'] and not t5['pass'] and not t6['pass']
    )

    if qed_baselines_ok and bam_constructions_fail:
        # Check whether the diagonals match in Construction III
        diag_ok = (
            t6['bhabha_max_diagonal_residual'] < 1e-6
            and t6['moller_max_diagonal_residual'] < 1e-6
        )
        if diag_ok:
            verdict_class = 'DIAGONAL_OK_INTERFERENCE_FAILS'
            verdict = (
                'DIAGONAL_OK_INTERFERENCE_FAILS. The BAM Compton-F² '
                'inspired vertex construction reproduces the diagonal '
                'channel intensities for Bhabha and Møller, but the '
                'interference (Fermi-statistics-determined sign and '
                'Dirac-trace magnitude) cannot be reproduced by any '
                'positive-real scalar BAM construction. '
                'Required extension: Dirac spinor structure at the '
                'vertex (or its BAM geometric proxy).'
            )
        else:
            verdict_class = 'DIAGONAL_FAILS_TOO'
            verdict = (
                'DIAGONAL_FAILS_TOO. Even the diagonal channel '
                'intensities for Bhabha/Møller do not match the BAM '
                'Compton-F² analytic continuation — these 4-fermion '
                'processes sit outside the Compton kernel\'s '
                'crossing-orbit class. A fundamentally new BAM '
                'tree kernel is needed for multi-fermion scattering, '
                'with explicit Dirac-trace structure at the vertex '
                'and a virtual-photon propagator beyond 1/q². The '
                'Compton/BW/ann triangle does NOT extend to Bhabha/'
                'Møller by Mandelstam crossing alone.'
            )
    elif not qed_baselines_ok:
        verdict_class = 'QED_BASELINE_BROKEN'
        verdict = (
            'QED BASELINE BROKEN. One of the QED reference '
            'computations failed.'
        )
    else:
        verdict_class = 'UNEXPECTED_BAM_SUCCESS'
        verdict = (
            'Unexpectedly, one of the BAM constructions matched '
            'Bhabha/Møller. Investigate which and why.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
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
    L.append('# Bhabha / Møller two-channel interference probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Tests whether the BAM scalar-intensity tree kernel can '
        'reproduce tree-level QED Bhabha/Møller, where two crossed '
        'kernels must be coherently summed with a Fermi-statistics-'
        'determined relative sign.'
    )
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '*falsified (documented)*'
        nm = t['name']
        if nm.startswith('T1'):
            value = 'QED Bhabha intensity tabulated at 5 CM angles'
        elif nm.startswith('T2'):
            value = 'QED Møller intensity tabulated at 5 CM angles'
        elif nm.startswith('T3'):
            value = (
                f"Bhabha all-negative: {t['bhabha_all_negative']}; "
                f"Møller all-positive: {t['moller_all_positive']}"
            )
        elif nm.startswith('T4'):
            value = (
                f"Bhabha max residual = {t['bhabha_max_total_residual']:.3e}; "
                f"Møller max residual = {t['moller_max_total_residual']:.3e}"
            )
        elif nm.startswith('T5'):
            value = (
                f"Bhabha max residual = {t['bhabha_max_total_residual']:.3e}; "
                f"Møller max residual = {t['moller_max_total_residual']:.3e}"
            )
        elif nm.startswith('T6'):
            value = (
                f"Bhabha diag residual = {t['bhabha_max_diagonal_residual']:.3e}; "
                f"Møller diag residual = {t['moller_max_diagonal_residual']:.3e}"
            )
        elif nm.startswith('T7'):
            value = 'residuals identified: Dirac numerators + interference magnitude'
        elif nm.startswith('T8'):
            value = (
                f"{len(t['missing_ingredients'])} missing ingredients; "
                f"{len(t['next_probe_candidates'])} next-probe candidates"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1
    t1 = s['tests'][0]
    L.append('## T1: QED Bhabha scalar intensity')
    L.append('')
    L.append('| θ (deg) | s | t | u | M²_t | M²_s | interference | total |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|---:|')
    for r in t1['rows']:
        L.append(
            f"| {r['theta_deg']} | {r['s']:.3f} | {r['t']:+.3f} | "
            f"{r['u']:+.3f} | {r['M2_t_channel']:.4f} | "
            f"{r['M2_s_channel']:.4f} | "
            f"{r['interference_2u2_over_st']:+.4f} | "
            f"{r['total_M2_over_8e4']:.4f} |"
        )
    L.append('')

    # T2
    t2 = s['tests'][1]
    L.append('## T2: QED Møller scalar intensity')
    L.append('')
    L.append('| θ (deg) | s | t | u | M²_t | M²_u | interference | total |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|---:|')
    for r in t2['rows']:
        L.append(
            f"| {r['theta_deg']} | {r['s']:.3f} | {r['t']:+.3f} | "
            f"{r['u']:+.3f} | {r['M2_t_channel']:.4f} | "
            f"{r['M2_u_channel']:.4f} | "
            f"{r['interference_2s2_over_tu']:+.4f} | "
            f"{r['total_M2_over_8e4']:.4f} |"
        )
    L.append('')

    # T3
    t3 = s['tests'][2]
    L.append('## T3: Interference signs')
    L.append('')
    L.append(f"Bhabha signs across angles: `{t3['bhabha_signs_across_angles']}`")
    L.append(f"Møller signs across angles: `{t3['moller_signs_across_angles']}`")
    L.append(
        f"Bhabha all negative: **{t3['bhabha_all_negative']}**; "
        f"Møller all positive: **{t3['moller_all_positive']}**."
    )
    L.append('')

    # T4
    t4 = s['tests'][3]
    L.append('## T4: BAM Construction I (positive √-channel sum)')
    L.append('')
    L.append(
        '| θ | Bhabha QED interf | Bhabha BAM+ interf | Bhabha total resid '
        '| Møller QED interf | Møller BAM+ interf | Møller total resid |'
    )
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t4['rows']:
        L.append(
            f"| {r['theta_deg']} | "
            f"{r['bhabha_QED_interference']:+.4f} | "
            f"{r['bhabha_BAM_positive_interference']:+.4f} | "
            f"{r['bhabha_total_residual']:+.4f} | "
            f"{r['moller_QED_interference']:+.4f} | "
            f"{r['moller_BAM_positive_interference']:+.4f} | "
            f"{r['moller_total_residual']:+.4f} |"
        )
    L.append('')

    # T5
    t5 = s['tests'][4]
    L.append('## T5: BAM Construction II (Fermi sign by hand)')
    L.append('')
    L.append(
        '| θ | Bhabha QED interf | Bhabha BAM± interf | Bhabha total resid '
        '| Møller QED interf | Møller BAM± interf | Møller total resid |'
    )
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t5['rows']:
        L.append(
            f"| {r['theta_deg']} | "
            f"{r['bhabha_QED_interference']:+.4f} | "
            f"{r['bhabha_BAM_signed_interference']:+.4f} | "
            f"{r['bhabha_total_residual']:+.4f} | "
            f"{r['moller_QED_interference']:+.4f} | "
            f"{r['moller_BAM_signed_interference']:+.4f} | "
            f"{r['moller_total_residual']:+.4f} |"
        )
    L.append('')

    # T6
    t6 = s['tests'][5]
    L.append('## T6: BAM Construction III (Compton-F² vertex factor)')
    L.append('')
    L.append('| θ | Bhabha QED M²_t | BAM M²_t | resid | Bhabha QED M²_s | BAM M²_s | resid |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t6['rows']:
        L.append(
            f"| {r['theta_deg']} | "
            f"{r['bhabha_QED_M2_t']:.4f} | "
            f"{r['bhabha_BAM_M2_t_via_compton_vertex']:.4f} | "
            f"{r['bhabha_t_diag_residual']:+.4f} | "
            f"{r['bhabha_QED_M2_s']:.4f} | "
            f"{r['bhabha_BAM_M2_s_via_compton_vertex']:.4f} | "
            f"{r['bhabha_s_diag_residual']:+.4f} |"
        )
    L.append('')

    # T7
    t7 = s['tests'][6]
    L.append('## T7: Residual analysis')
    L.append('')
    L.append('### Bhabha interference vs geometric-mean (BAM ansatz)')
    L.append('')
    L.append('| θ | |QED interference| | 2·√(M²_t · M²_s) BAM | ratio |')
    L.append('|---:|---:|---:|---:|')
    for r in t7['interference_vs_geometric_mean_analysis']:
        L.append(
            f"| {r['theta_deg']} | "
            f"{r['bhabha_|interference|_QED']:.4f} | "
            f"{r['bhabha_2x_geometric_mean_diagonals_BAM_positive']:.4f} | "
            f"{r['ratio_BAM_over_QED']:.3f} |"
        )
    L.append('')

    # T8
    t8 = s['tests'][7]
    L.append('## T8: Required BAM extension')
    L.append('')
    L.append('### Missing ingredients')
    for i, ing in enumerate(t8['missing_ingredients'], 1):
        L.append(f"{i}. {ing}")
    L.append('')
    L.append('### Next-probe candidates')
    for i, cand in enumerate(t8['next_probe_candidates'], 1):
        L.append(f"{i}. {cand}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Dirac-trace structure in BAM**: capturing the QED γ-matrix '
        'algebra that produces (s²+u²), (u²+t²), (s²+t²) numerators '
        'requires either explicit Dirac spinors at vertices or a '
        'geometric proxy.'
    )
    L.append(
        '- **Möbius/antipodal Pauli signs**: BAM\'s `T = iσ_y` '
        'antipodal transport (README channel 2) is a candidate origin '
        'for the Møller t↔u Pauli sign. Whether it gives the right '
        'sign automatically is a separate probe target.'
    )
    L.append(
        '- **Virtual-photon throat-fibre propagator**: 1/q² alone '
        'does not reproduce QED; additional Hopf-fibre phase structure '
        'may be needed.'
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
    out = here / 'runs' / f'{ts}_bhabha_moller_interference_probe'
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
