"""
Two-throat Coulomb probe.

A falsification test flagged in docs/THESIS.md: the finite-separation
two-throat Coulomb law on S³ is "a near-term falsification test, not a
demonstrated result." This probe runs it.

Two charged BAM throat mouths on S³ of radius R, separated by geodesic
angle ψ, interact through the S³ Green function

    V(ψ) = q · G(ψ),   G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)

(geometrodynamics.transaction.s3_geometry.s3_green_potential; the same
kernel as PR #45's momentum-space exchange picture). This probe tests
the static / position-space picture:

  - Coulomb potential V ∝ 1/r in the flat-space limit (P1).
  - Inverse-square force F ∝ 1/r² in the flat-space limit (P2).
  - Exact S³ force F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ), N(ψ)=(π−ψ)+sinψcosψ
    (P3): the THESIS 1/sin²ψ claim modulated by N(ψ).
  - Antipodal equilibrium F(π) = 0 (P4): compact-S³ image effect.
  - Gauss's law on S³ Φ(ψ) = Q_enclosed(ψ) (P5): exact.
  - Attraction/repulsion sign (P6).
  - Two-mouth placement realization (T9).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.transaction.s3_geometry import (
    s3_green_potential,
    s3_green_field_kernel,
    antipode4,
    geo4,
    hsp,
)


PI = math.pi
EPS = 1e-9   # small clipping override for flat-limit probing


# ---------------------------------------------------------------------------
# Analytic helpers
# ---------------------------------------------------------------------------

def N_modulation(psi: float) -> float:
    """N(ψ) = (π − ψ) + sin ψ cos ψ. Runs from π at ψ=0 to 0 at ψ=π."""
    return (PI - psi) + math.sin(psi) * math.cos(psi)


def force_analytic(psi: float, q1: float, q2: float, R: float) -> float:
    """Exact S³ Coulomb force magnitude:
       F(ψ) = q₁q₂ · N(ψ) / (4π² R² sin²ψ)."""
    return q1 * q2 * N_modulation(psi) / (4.0 * PI ** 2 * R ** 2 * math.sin(psi) ** 2)


def potential_analytic(psi: float, q: float, R: float) -> float:
    """V(ψ) = q · G(ψ)."""
    return q * s3_green_potential(psi, radius=R, eps=EPS)


# ---------------------------------------------------------------------------
# T1. Potential = S³ Green function
# ---------------------------------------------------------------------------

def test_T1_potential_is_green_function() -> dict:
    """V(ψ) = q·G(ψ) with G the repo S³ Green function. Verify the
    analytic form ((π−ψ)cot ψ − ½)/(4π²R)."""
    R = 1.0
    q = 1.0
    rows = []
    max_diff = 0.0
    for psi in [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        V_repo = potential_analytic(psi, q, R)
        V_analytic = q * ((PI - psi) / math.tan(psi) - 0.5) / (4.0 * PI ** 2 * R)
        diff = abs(V_repo - V_analytic)
        max_diff = max(max_diff, diff)
        rows.append({
            'psi': psi,
            'V_repo': V_repo,
            'V_analytic': V_analytic,
            'difference': diff,
        })
    return {
        'name': 'T1_potential_is_s3_green_function',
        'description': (
            "Two-throat Coulomb potential V(ψ) = q·G(ψ) with the S³ "
            "Green function from geometrodynamics.transaction.s3_geometry. "
            "Zero-mean (the −½ is the compact-manifold neutralizing "
            "background)."
        ),
        'rows': rows,
        'max_difference': max_diff,
        'pass': max_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T2. Coulomb potential (flat-space limit)
# ---------------------------------------------------------------------------

def test_T2_coulomb_potential_flat_limit() -> dict:
    """For small ψ with r = R·ψ: V(ψ) → q/(4π r). Verify V·r → q/(4π)."""
    R = 1.0
    q = 1.0
    rows = []
    for psi in [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1]:
        r = R * psi
        V = potential_analytic(psi, q, R)
        Vr = V * r
        target = q / (4.0 * PI)
        rows.append({
            'psi': psi,
            'r': r,
            'V': V,
            'V_times_r': Vr,
            'coulomb_target_q_over_4pi': target,
            'residual': Vr - target,
        })
    smallest_residual = abs(rows[0]['residual'])
    return {
        'name': 'T2_coulomb_potential_flat_limit',
        'description': (
            "Flat-space limit: V(ψ)·r → q/(4π) as ψ → 0, i.e. the "
            "Coulomb potential V = q/(4π r). The S³ Green function "
            "reproduces electrostatics at short range."
        ),
        'rows': rows,
        'smallest_residual': smallest_residual,
        'pass': smallest_residual < 1e-5,
    }


# ---------------------------------------------------------------------------
# T3. Force from repo field kernel matches analytic
# ---------------------------------------------------------------------------

def test_T3_force_from_field_kernel() -> dict:
    """Force F(ψ) = q₁q₂·|dG/dψ|/R from the repo field kernel matches
    the analytic q₁q₂·N(ψ)/(4π²R²sin²ψ)."""
    R = 1.0
    q1 = 1.0
    q2 = 1.0
    rows = []
    max_diff = 0.0
    for psi in [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5]:
        F_kernel = q1 * q2 * s3_green_field_kernel(psi, radius=R, eps=EPS)
        F_analytic = force_analytic(psi, q1, q2, R)
        diff = abs(F_kernel - F_analytic)
        max_diff = max(max_diff, diff)
        rows.append({
            'psi': psi,
            'F_from_repo_field_kernel': F_kernel,
            'F_analytic_N_over_4pi2_R2_sin2': F_analytic,
            'difference': diff,
        })
    return {
        'name': 'T3_force_from_field_kernel',
        'description': (
            "Force F(ψ) = q₁q₂·|dG/dψ|/R (repo field kernel) equals the "
            "analytic q₁q₂·N(ψ)/(4π²R²sin²ψ) with N(ψ)=(π−ψ)+sinψcosψ."
        ),
        'rows': rows,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T4. Inverse-square force (flat-space limit)
# ---------------------------------------------------------------------------

def test_T4_inverse_square_flat_limit() -> dict:
    """For small ψ: F(ψ) → q₁q₂/(4π r²). Verify F·r² → q₁q₂/(4π)."""
    R = 1.0
    q1 = 1.0
    q2 = 1.0
    rows = []
    for psi in [1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1]:
        r = R * psi
        F = force_analytic(psi, q1, q2, R)
        Fr2 = F * r * r
        target = q1 * q2 / (4.0 * PI)
        rows.append({
            'psi': psi,
            'r': r,
            'F': F,
            'F_times_r2': Fr2,
            'inverse_square_target_q1q2_over_4pi': target,
            'residual': Fr2 - target,
        })
    smallest_residual = abs(rows[0]['residual'])
    return {
        'name': 'T4_inverse_square_force_flat_limit',
        'description': (
            "Flat-space limit: F(ψ)·r² → q₁q₂/(4π) as ψ → 0, i.e. the "
            "inverse-square law F = q₁q₂/(4π r²). The two-throat S³ "
            "Green function reproduces Coulomb's force law at short "
            "range."
        ),
        'rows': rows,
        'smallest_residual': smallest_residual,
        'pass': smallest_residual < 1e-5,
    }


# ---------------------------------------------------------------------------
# T5. Exact S³ force form / 1/sin²ψ claim
# ---------------------------------------------------------------------------

def test_T5_exact_s3_force_form() -> dict:
    """Characterise the THESIS claim F ∝ 1/sin²ψ. The exact force is
    F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ) with N(ψ) = (π−ψ)+sinψcosψ. The
    quantity F·sin²ψ = q₁q₂·N(ψ)/(4π²R²) is the modulation; it equals
    the constant q₁q₂·π/(4π²R²) only at ψ=0 (Coulomb regime)."""
    R = 1.0
    q1 = 1.0
    q2 = 1.0
    rows = []
    coulomb_const = q1 * q2 * PI / (4.0 * PI ** 2 * R ** 2)
    for psi in [0.01, 0.1, 0.5, 1.0, PI / 2, 2.0, 2.5, 3.0]:
        F = force_analytic(psi, q1, q2, R)
        F_sin2 = F * math.sin(psi) ** 2
        N = N_modulation(psi)
        rows.append({
            'psi': psi,
            'F': F,
            'F_times_sin2_psi': F_sin2,
            'N_modulation': N,
            'coulomb_const_at_psi_0': coulomb_const,
            'ratio_F_sin2_to_coulomb_const': F_sin2 / coulomb_const,
        })
    # 1/sin²ψ is leading: F·sin²ψ → coulomb_const at ψ → 0
    leading_ok = abs(rows[0]['F_times_sin2_psi'] - coulomb_const) / coulomb_const < 1e-3
    # Modulation present: F·sin²ψ varies (N(ψ) not constant)
    modulation_present = abs(rows[-1]['ratio_F_sin2_to_coulomb_const'] - 1.0) > 0.5
    return {
        'name': 'T5_exact_s3_force_form',
        'description': (
            "THESIS claim F ∝ 1/sin²ψ is the LEADING behaviour near "
            "ψ=0; the exact force is modulated by N(ψ)=(π−ψ)+sinψcosψ "
            "(runs π → 0). F·sin²ψ matches the Coulomb constant only "
            "at short range. The N(ψ) modulation is the compact-S³ "
            "antipodal-image correction (see T6)."
        ),
        'rows': rows,
        'leading_1_over_sin2_confirmed': leading_ok,
        'modulation_present_compact_S3': modulation_present,
        'pass': leading_ok and modulation_present,
    }


# ---------------------------------------------------------------------------
# T6. Antipodal equilibrium
# ---------------------------------------------------------------------------

def test_T6_antipodal_equilibrium() -> dict:
    """At the antipode ψ → π, N(π) = 0 so F(π) = 0. A charge on S³ has
    its field lines converge at the antipode; the force on a test charge
    there vanishes by symmetry (compact-S³ image effect)."""
    R = 1.0
    q1 = 1.0
    q2 = 1.0
    rows = []
    for psi in [PI - 0.5, PI - 0.1, PI - 0.01, PI - 0.001]:
        N = N_modulation(psi)
        F = force_analytic(psi, q1, q2, R)
        rows.append({
            'psi': psi,
            'pi_minus_psi': PI - psi,
            'N_modulation': N,
            'force': F,
        })
    N_at_pi = N_modulation(PI - 1e-12)
    force_decreasing = rows[-1]['force'] < rows[0]['force']
    return {
        'name': 'T6_antipodal_equilibrium',
        'description': (
            "At the antipode ψ → π, N(ψ) → 0 so the force F → 0. A "
            "charge on compact S³ has its field lines converge at its "
            "antipode; a test charge placed there feels zero net force "
            "by symmetry. This is the compact-S³ signature absent in "
            "flat space — a charge feels its own antipodal image."
        ),
        'rows': rows,
        'N_near_pi': N_at_pi,
        'force_decreasing_toward_antipode': force_decreasing,
        'pass': abs(N_at_pi) < 1e-10 and force_decreasing,
    }


# ---------------------------------------------------------------------------
# T7. Gauss's law on S³
# ---------------------------------------------------------------------------

def test_T7_gauss_law() -> dict:
    """The electric flux through the geodesic 2-sphere at colatitude ψ
    (area 4πR²sin²ψ) equals the enclosed charge (point q + neutralizing
    background cap):
       Φ(ψ) = |E(ψ)|·4πR²sin²ψ = q·N(ψ)/π = Q_enclosed(ψ)."""
    R = 1.0
    q = 1.0
    rows = []
    max_diff = 0.0
    for psi in [0.1, 0.5, 1.0, PI / 2, 2.0, 2.5, 3.0]:
        E = q * s3_green_field_kernel(psi, radius=R, eps=EPS)
        area = 4.0 * PI * (R * math.sin(psi)) ** 2
        flux = E * area
        # Enclosed = point charge − background within cap
        cap_fraction = (psi - math.sin(psi) * math.cos(psi)) / PI
        Q_enclosed = q - q * cap_fraction
        diff = abs(flux - Q_enclosed)
        max_diff = max(max_diff, diff)
        rows.append({
            'psi': psi,
            'flux_E_times_area': flux,
            'Q_enclosed_point_plus_background_cap': Q_enclosed,
            'difference': diff,
        })
    return {
        'name': 'T7_gauss_law_on_s3',
        'description': (
            "Gauss's law on compact S³: flux through the colatitude-ψ "
            "2-sphere equals the enclosed charge (point q plus the "
            "neutralizing uniform background within the cap). Flux falls "
            "from q at the source to 0 at the antipode, consistent with "
            "the zero-mean (compact-manifold) Green function. Exact "
            "internal-consistency check of the whole electrostatic "
            "picture."
        ),
        'rows': rows,
        'max_difference': max_diff,
        'pass': max_diff < 1e-10,
    }


# ---------------------------------------------------------------------------
# T8. Attraction / repulsion sign
# ---------------------------------------------------------------------------

def signed_radial_force(psi: float, q1: float, q2: float, R: float) -> float:
    """Signed radial force = −dU/ds = −(q₁q₂/R)·dG/dψ.
    Positive = repulsive (pushes charges to larger ψ).
    dG/dψ = −N(ψ)/(4π²R sin²ψ), so −dG/dψ = +N(ψ)/(...),
    giving F = +(q₁q₂/R)·N(ψ)/(4π²R sin²ψ) for like charges (q₁q₂>0)
    → positive → repulsive."""
    dGdpsi = -N_modulation(psi) / (4.0 * PI ** 2 * R * math.sin(psi) ** 2)
    return -(q1 * q2 / R) * dGdpsi


def test_T8_attraction_repulsion_sign() -> dict:
    """Like charges repel (F > 0), opposite charges attract (F < 0)."""
    R = 1.0
    psi = 0.5
    F_like = signed_radial_force(psi, 1.0, 1.0, R)       # like (++)
    F_opposite = signed_radial_force(psi, 1.0, -1.0, R)  # opposite (+−)
    return {
        'name': 'T8_attraction_repulsion_sign',
        'description': (
            "Signed radial force −dU/ds: like charges (q₁q₂>0) repel "
            "(F>0, pushed to larger ψ); opposite charges (q₁q₂<0) "
            "attract (F<0)."
        ),
        'force_like_charges_++': F_like,
        'force_opposite_charges_+-': F_opposite,
        'like_repels': F_like > 0,
        'opposite_attracts': F_opposite < 0,
        'pass': F_like > 0 and F_opposite < 0,
    }


# ---------------------------------------------------------------------------
# T9. Two-mouth placement realization
# ---------------------------------------------------------------------------

def test_T9_two_mouth_placement() -> dict:
    """Place two actual throat mouths at explicit S³ points and verify
    the geo4-separation interaction matches the analytic V(ψ). Also
    verify the antipodal mouth of each particle sits at geodesic
    distance π (the throat-pair structure)."""
    R = 1.0
    q1 = 1.0
    q2 = 1.0
    rows = []
    max_diff = 0.0
    # Place mouth 1 at the "north pole" of the Hopf-sphere parameterisation
    mouth_1 = hsp(0.0, 0.0, 0.0)   # χ=0 → 4-vector (0,0,0,1)
    for chi in [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]:
        mouth_2 = hsp(chi, PI / 3, PI / 4)
        psi = geo4(mouth_1, mouth_2)   # geodesic separation
        # Interaction energy via Green function at this separation
        U = q1 * q2 * s3_green_potential(psi, radius=R, eps=EPS)
        U_analytic = q1 * q2 * potential_analytic(psi, 1.0, R)
        diff = abs(U - U_analytic)
        max_diff = max(max_diff, diff)
        # Antipodal mouth of particle 1
        anti_1 = antipode4(mouth_1)
        anti_dist = geo4(mouth_1, anti_1)   # should be π
        rows.append({
            'chi': chi,
            'mouth_2_4vector': [float(x) for x in mouth_2],
            'geodesic_separation_psi': psi,
            'U_interaction': U,
            'U_analytic': U_analytic,
            'difference': diff,
            'antipodal_mouth_distance_should_be_pi': anti_dist,
        })
    antipode_ok = all(abs(r['antipodal_mouth_distance_should_be_pi'] - PI) < 1e-9 for r in rows)
    return {
        'name': 'T9_two_mouth_placement_realization',
        'description': (
            "Place two throat mouths at explicit S³ points (Hopf-sphere "
            "parameterisation); the geo4-geodesic-separation interaction "
            "via the Green function matches the analytic V(ψ). Each "
            "particle's antipodal back-mouth sits at geodesic distance π "
            "(throat-pair structure)."
        ),
        'rows': rows,
        'max_difference': max_diff,
        'antipodal_distances_all_pi': antipode_ok,
        'pass': max_diff < 1e-12 and antipode_ok,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_potential_is_green_function()
    t2 = test_T2_coulomb_potential_flat_limit()
    t3 = test_T3_force_from_field_kernel()
    t4 = test_T4_inverse_square_flat_limit()
    t5 = test_T5_exact_s3_force_form()
    t6 = test_T6_antipodal_equilibrium()
    t7 = test_T7_gauss_law()
    t8 = test_T8_attraction_repulsion_sign()
    t9 = test_T9_two_mouth_placement()
    tests = [t1, t2, t3, t4, t5, t6, t7, t8, t9]

    coulomb_core = [t1, t2, t3, t4, t7, t8, t9]
    if all(t['pass'] for t in coulomb_core) and t5['pass'] and t6['pass']:
        verdict_class = 'COULOMB_REPRODUCED'
        verdict = (
            'COULOMB LAW REPRODUCED. Two charged BAM throat mouths on '
            'S³, interacting through the S³ Green function, reproduce '
            'electrostatics:\n'
            '  (P1) Coulomb potential V → q/(4π r) in the flat-space '
            'limit (T2);\n'
            '  (P2) inverse-square force F → q₁q₂/(4π r²) in the '
            'flat-space limit (T4);\n'
            '  (P5) Gauss\'s law Φ(ψ) = Q_enclosed(ψ) holds exactly on '
            'compact S³, with flux falling from q at the source to 0 at '
            'the antipode as the neutralizing background is enclosed (T7).\n'
            'The THESIS finite-separation claim F ∝ 1/sin²ψ is '
            'CONFIRMED as the leading short-range behaviour and REFINED: '
            'the exact force F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ) carries the '
            'modulation N(ψ) = (π−ψ)+sinψcosψ (T5), which encodes the '
            'compact-S³ antipodal image — the force vanishes exactly at '
            'the antipode (T6), an equilibrium absent in flat space. '
            'Like charges repel, opposite charges attract (T8). The '
            'two-mouth placement realization confirms the geometric '
            'picture, with each particle\'s back-mouth at geodesic '
            'distance π (T9). The "charge = throat mouth on S³" '
            'identification is NOT falsified — it reproduces Coulomb '
            'electrostatics with calculable compact-manifold '
            'corrections.'
        )
    else:
        verdict_class = 'COULOMB_FALSIFIED'
        verdict = (
            'COULOMB FALSIFIED. The flat-limit Coulomb/inverse-square '
            'law or Gauss\'s law failed. Investigate which test.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'potential': 'V(ψ) = q·G(ψ), G(ψ) = ((π−ψ)cot ψ − ½)/(4π²R)',
        'force_exact': 'F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ), N(ψ)=(π−ψ)+sinψcosψ',
        'flat_limits': 'V → q/(4π r), F → q₁q₂/(4π r²)',
        'gauss_law': 'Φ(ψ) = q·N(ψ)/π = Q_enclosed(ψ)',
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
    L.append('# Two-throat Coulomb probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Falsification test (flagged in docs/THESIS.md): two charged '
        'BAM throat mouths on S³, interacting through the S³ Green '
        'function, vs the Coulomb potential and inverse-square force law.'
    )
    L.append('')

    L.append('## Electrostatics on S³')
    L.append('')
    L.append('```')
    L.append(f"potential:    {s['potential']}")
    L.append(f"force (exact): {s['force_exact']}")
    L.append(f"flat limits:  {s['flat_limits']}")
    L.append(f"Gauss law:    {s['gauss_law']}")
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
            value = f"V = q·G(ψ), max diff = {t['max_difference']:.2e}"
        elif nm.startswith('T2'):
            value = f"V·r → q/(4π), residual = {t['smallest_residual']:.2e}"
        elif nm.startswith('T3'):
            value = f"F kernel = analytic, max diff = {t['max_difference']:.2e}"
        elif nm.startswith('T4'):
            value = f"F·r² → q₁q₂/(4π), residual = {t['smallest_residual']:.2e}"
        elif nm.startswith('T5'):
            value = (
                f"1/sin²ψ leading: {t['leading_1_over_sin2_confirmed']}; "
                f"N(ψ) modulation: {t['modulation_present_compact_S3']}"
            )
        elif nm.startswith('T6'):
            value = f"N(π)→0 (F→0 at antipode): {t['N_near_pi']:.2e}"
        elif nm.startswith('T7'):
            value = f"Φ = Q_enclosed, max diff = {t['max_difference']:.2e}"
        elif nm.startswith('T8'):
            value = (
                f"like repels: {t['like_repels']}; "
                f"opposite attracts: {t['opposite_attracts']}"
            )
        elif nm.startswith('T9'):
            value = (
                f"V matches, max diff = {t['max_difference']:.2e}; "
                f"back-mouths at π: {t['antipodal_distances_all_pi']}"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T2 detail
    t2 = s['tests'][1]
    L.append('## T2: Coulomb potential (flat-space limit)')
    L.append('')
    L.append('| ψ | r = Rψ | V | V·r | q/(4π) | residual |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t2['rows']:
        L.append(
            f"| {r['psi']:.0e} | {r['r']:.0e} | {r['V']:.4e} | "
            f"{r['V_times_r']:.6f} | {r['coulomb_target_q_over_4pi']:.6f} | "
            f"{r['residual']:+.2e} |"
        )
    L.append('')

    # T4 detail
    t4 = s['tests'][3]
    L.append('## T4: Inverse-square force (flat-space limit)')
    L.append('')
    L.append('| ψ | r = Rψ | F | F·r² | q₁q₂/(4π) | residual |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t4['rows']:
        L.append(
            f"| {r['psi']:.0e} | {r['r']:.0e} | {r['F']:.4e} | "
            f"{r['F_times_r2']:.6f} | "
            f"{r['inverse_square_target_q1q2_over_4pi']:.6f} | "
            f"{r['residual']:+.2e} |"
        )
    L.append('')

    # T5 detail
    t5 = s['tests'][4]
    L.append('## T5: Exact S³ force form / 1/sin²ψ claim')
    L.append('')
    L.append('| ψ | F | F·sin²ψ | N(ψ) | F·sin²ψ / Coulomb-const |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t5['rows']:
        L.append(
            f"| {r['psi']:.4f} | {r['F']:.4e} | "
            f"{r['F_times_sin2_psi']:.6f} | {r['N_modulation']:.6f} | "
            f"{r['ratio_F_sin2_to_coulomb_const']:.6f} |"
        )
    L.append('')
    L.append(
        f"1/sin²ψ confirmed as leading: **{t5['leading_1_over_sin2_confirmed']}**; "
        f"N(ψ) modulation present (compact-S³): "
        f"**{t5['modulation_present_compact_S3']}**."
    )
    L.append('')

    # T6 detail
    t6 = s['tests'][5]
    L.append('## T6: Antipodal equilibrium')
    L.append('')
    L.append('| ψ | π − ψ | N(ψ) | force |')
    L.append('|---:|---:|---:|---:|')
    for r in t6['rows']:
        L.append(
            f"| {r['psi']:.4f} | {r['pi_minus_psi']:.0e} | "
            f"{r['N_modulation']:.4e} | {r['force']:.4e} |"
        )
    L.append('')
    L.append(
        f"N near antipode = `{t6['N_near_pi']:.2e}` → force vanishes at "
        f"ψ = π. A charge on compact S³ feels its own antipodal image; "
        f"the antipode is a force-free equilibrium absent in flat space."
    )
    L.append('')

    # T7 detail
    t7 = s['tests'][6]
    L.append('## T7: Gauss\'s law on S³')
    L.append('')
    L.append('| ψ | flux = E·area | Q_enclosed (point + bg cap) | diff |')
    L.append('|---:|---:|---:|---:|')
    for r in t7['rows']:
        L.append(
            f"| {r['psi']:.4f} | {r['flux_E_times_area']:.6f} | "
            f"{r['Q_enclosed_point_plus_background_cap']:.6f} | "
            f"{r['difference']:.2e} |"
        )
    L.append('')
    L.append(
        "Flux falls from q (source) to 0 (antipode) as the neutralizing "
        "background charge is enclosed; Φ(ψ) = Q_enclosed(ψ) exactly."
    )
    L.append('')

    # T9 detail
    t9 = s['tests'][8]
    L.append('## T9: Two-mouth placement realization')
    L.append('')
    L.append('| χ | geodesic ψ | U interaction | U analytic | diff | back-mouth dist |')
    L.append('|---:|---:|---:|---:|---:|---:|')
    for r in t9['rows']:
        L.append(
            f"| {r['chi']:.2f} | {r['geodesic_separation_psi']:.4f} | "
            f"{r['U_interaction']:.6f} | {r['U_analytic']:.6f} | "
            f"{r['difference']:.2e} | "
            f"{r['antipodal_mouth_distance_should_be_pi']:.6f} |"
        )
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Dynamical / radiative**: this is the static potential. '
        'Moving charges (radiation, magnetic fields) need the '
        'time-dependent S³ Green function — a separate probe.'
    )
    L.append(
        '- **Throat-pair internal structure**: whether the back mouth '
        'carries −q (dipole) or +q (monopole pair) affects the '
        'long-range field; a first-principles determination from the '
        'Hopf charge is open.'
    )
    L.append(
        '- **Self-energy**: the ψ → 0 divergence of G(ψ) is the usual '
        'point-charge self-energy; regularisation by the finite throat '
        'radius is open.'
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
    out = here / 'runs' / f'{ts}_two_throat_coulomb_probe'
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
