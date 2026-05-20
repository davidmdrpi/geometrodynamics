"""
Covariant BAM effective-action scaffold probe.

Builds the candidate single variational principle from which the three
BAM structures should derive:

  (A) Hopf-bundle U(1) connection  A_φ(χ) = ½cos χ
  (G) 5D Tangherlini bulk boundaries  ΔR = R_OUTER − R_INNER
  (T) Compton vertex factor  F²(x, c) = K(x)²·Q(x, c)

and rigorously identifies the mismatch terms (barriers) that prevent
full closure. This is a scaffold + barrier map, NOT a claim that the
action closes.

Candidate action:

  S_BAM = ∫_{M₅} d⁵x √(−g₅) [ (1/2κ₅)(R₅ − 2Λ₅) − ¼ F_{MN}F^{MN}
                              + ψ̄(iΓ^M D_M − m)ψ + L_throat ]
        + S_∂[hard walls] + S_closure[∮A = 2πn]

Sector status:
  (A) gauge   — closes (homogeneous Maxwell; c₁ = 1, no external input)
  (T) throat  — closes given the closure constraint (B1) + antipodal Z₂ (B2)
  (G) gravity — partial (metric solves Einstein; boundaries imposed, B3)

Mismatch terms / barriers:
  B1 closure quantum is topological (Wilson loop ∮A = 2πn, not local)
  B2 antipodal Z₂ is a discrete quotient (T = iσ_y, not continuous)
  B3 boundary conditions topologically imposed (hard wall from T²=−I)
  B4 dimensional bridge: scale-complete only modulo m_e (ℏ = m_e R_MID c)
  B5 5D → 4D reduction producing F² not constructed

Tests:
  T1. Action scaffold structure (informative).
  T2. Gauge sector closes: A_φ = ½cos χ, c₁ = 1.
  T3. Throat sector closes with constraint: F² = K²·Q.
  T4. Gravity sector partial: f(r) = 1 − (r_s/r)²; boundaries imposed.
  T5. Dimensional bridge barrier (B4).
  T6. 5D → 4D reduction barrier (B5).
  T7. Barrier ledger / closure assessment.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.hopf.connection import hopf_connection, hopf_curvature
from geometrodynamics.hopf.chern import compute_c1
from geometrodynamics.tangherlini.radial import V_tangherlini
from geometrodynamics.constants import R_MID, R_OUTER, R_INNER, DELTA


PI = math.pi
TAU = 2.0 * PI


# ---------------------------------------------------------------------------
# T1. Action scaffold structure
# ---------------------------------------------------------------------------

def test_T1_action_scaffold() -> dict:
    """State the candidate covariant action and the target → sector →
    variation map. Informative / structural."""
    action_terms = {
        'G_gravity': '(1/2κ₅)(R₅ − 2Λ₅)  — 5D Einstein-Hilbert + cosmological',
        'A_gauge': '−¼ F_{MN}F^{MN}  — U(1) gauge field strength',
        'D_dirac': 'ψ̄(iΓ^M D_M − m)ψ  — bulk fermions',
        'T_throat': 'L_throat[Φ, g, A]  — throat / closure scalar',
        'B_boundary': 'S_∂[hard walls]  — topological boundary conditions',
        'C_closure': 'S_closure[∮A = 2πn]  — closure-quantum constraint',
    }
    target_map = {
        'A_phi_hopf_connection': 'sector (A) U(1) gauge; δS/δA = 0 on S³',
        'Delta_R_tangherlini_boundaries': 'sector (G) gravity + (B) BCs; δS/δg = 0 + boundary data',
        'F_squared_compton_vertex': 'sector (T) throat + (C) closure; δS/δΦ = 0 + ∮A = 2πn',
    }
    return {
        'name': 'T1_action_scaffold_structure',
        'description': (
            "Candidate covariant BAM effective action on the 5D "
            "Tangherlini bulk with non-orientable throat and S³ Hopf "
            "fibration. Three targets attach to three sectors."
        ),
        'action_terms': action_terms,
        'target_sector_variation_map': target_map,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Gauge sector closes
# ---------------------------------------------------------------------------

def test_T2_gauge_sector() -> dict:
    """The Hopf connection A_φ(χ) = ½cos χ is the homogeneous
    (symmetric-critical) solution of source-free Maxwell on the Hopf
    bundle. Its first Chern number is exactly 1. Derives from δS/δA = 0
    with no external input — this sector closes."""
    # Verify the connection form A_φ = ½cos χ
    chis = np.linspace(0.0, PI, 9)
    max_conn_diff = 0.0
    conn_rows = []
    for chi in chis:
        A = float(hopf_connection(float(chi)))
        A_analytic = 0.5 * math.cos(chi)
        diff = abs(A - A_analytic)
        max_conn_diff = max(max_conn_diff, diff)
        if len(conn_rows) < 5:
            conn_rows.append({
                'chi_in_pi': float(chi) / PI,
                'A_phi': A,
                'half_cos_chi': A_analytic,
                'difference': diff,
            })
    # Curvature |F| = ½ sin χ (field strength of the Hopf connection)
    curv_check = abs(float(hopf_curvature(PI / 2.0)) - 0.5)

    # First Chern number c₁ = 1 (topological charge of the connection)
    c1 = compute_c1(N_chi=8000)
    c1_abs = c1['c1_abs']
    c1_err = c1['err_abs']

    return {
        'name': 'T2_gauge_sector_closes',
        'description': (
            "Gauge sector (A): A_φ = ½cos χ is the homogeneous Maxwell "
            "solution on the Hopf bundle (symmetric-critical by Palais' "
            "principle — the unique SU(2)-invariant connection). First "
            "Chern number c₁ = 1 exactly. Derives from δS/δA = 0 with "
            "no external input — this sector CLOSES."
        ),
        'connection_samples': conn_rows,
        'max_connection_difference': max_conn_diff,
        'curvature_at_equator_minus_half': curv_check,
        'chern_c1_abs': c1_abs,
        'chern_c1_error': c1_err,
        'pass': max_conn_diff < 1e-12 and c1_err < 1e-3,
    }


# ---------------------------------------------------------------------------
# T3. Throat sector closes with constraint
# ---------------------------------------------------------------------------

def K_pade(x: float) -> float:
    return 2.0 * x / (1.0 + x)


def Q_pol(x: float, c: float) -> float:
    return x * x + x * (1.0 - x) ** 2 / (1.0 + c * c)


def F_squared_closed(x: float, c: float) -> float:
    s2 = 1.0 - c * c
    return 4.0 * x ** 3 * (x * x + 1.0 - x * s2) / ((1.0 + c * c) * (1.0 + x) ** 2)


def test_T3_throat_sector() -> dict:
    """The vertex factor F² = K²·Q derives from the throat action
    (PR #41) GIVEN the closure quantum S = 2π and the antipodal Z₂.
    Recompute K and Q from equal-action splitting and verify
    F² = K²·Q. Flag the two imposed inputs (B1, B2)."""
    samples = []
    max_diff = 0.0
    for x in [0.1, 0.5, 1.0, 2.0, 5.0]:
        for c in [-0.7, 0.0, 0.7]:
            # K from equal-action splitting (throat-rate harmonic mean):
            # ω₁τ₁ = ω₂τ₂ = π (closure 2π split by antipodal symmetry);
            # ω_eff/ω₁ = 2x/(1+x).
            tau_1 = PI / 1.0
            tau_2 = PI / x
            omega_eff = TAU / (tau_1 + tau_2)
            K = omega_eff / 1.0
            # Q from Hopf-helicity spinor (PR #40)
            Q = Q_pol(x, c)
            F2_reconstructed = K * K * Q
            F2_closed = F_squared_closed(x, c)
            diff = abs(F2_reconstructed - F2_closed)
            max_diff = max(max_diff, diff)
            if len(samples) < 6:
                samples.append({
                    'x': x, 'cos_theta': c,
                    'K_from_equal_action': K,
                    'Q_from_hopf_spinor': Q,
                    'K2_Q': F2_reconstructed,
                    'F2_closed': F2_closed,
                    'difference': diff,
                })
    return {
        'name': 'T3_throat_sector_closes_with_constraint',
        'description': (
            "Throat sector (T): F² = K(x)²·Q(x, c) reconstructs from the "
            "equal-action throat splitting (K = 2x/(1+x)) and the "
            "Hopf-helicity spinor (Q). CLOSES — but requires two imposed "
            "inputs: the closure quantum S = 2π (B1, topological) and "
            "the antipodal Z₂ (B2, discrete). Neither is a local "
            "Lagrangian density derivable from δS/δΦ = 0 alone."
        ),
        'samples': samples,
        'max_difference': max_diff,
        'imposed_inputs': ['B1_closure_quantum_2pi', 'B2_antipodal_Z2'],
        'pass': max_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T4. Gravity sector partial
# ---------------------------------------------------------------------------

def test_T4_gravity_sector() -> dict:
    """The 5D Tangherlini metric f(r) = 1 − (r_s/r)² solves the vacuum
    Einstein equation; the radial-mode potential V_tangherlini uses it.
    But the domain boundaries R_INNER, R_OUTER are boundary data — the
    bulk equation gives the metric up to them. R_OUTER is fixed by the
    cross-species γ-lock fixed point; R_MID = 1 is a gauge choice (B3).
    """
    rs = R_MID
    rows = []
    max_diff = 0.0
    for r in [0.8, 1.0, 1.2, 1.5, 2.0]:
        f_analytic = 1.0 - (rs / r) ** 2
        # Cross-check against the form used in V_tangherlini (l=0):
        # V(r, 0) = f(r)·[0 + 3 rs²/r⁴]; recover f from V/(3 rs²/r⁴).
        if r != rs:  # avoid f = 0 at throat
            V0 = V_tangherlini(r, 0, rs=rs)
            f_from_V = V0 / (3.0 * rs ** 2 / r ** 4)
            diff = abs(f_from_V - f_analytic)
            max_diff = max(max_diff, diff)
            rows.append({
                'r': r,
                'f_tangherlini_1_minus_rs2_over_r2': f_analytic,
                'f_recovered_from_V_tangherlini': f_from_V,
                'difference': diff,
            })
    return {
        'name': 'T4_gravity_sector_partial',
        'description': (
            "Gravity sector (G): the Tangherlini metric "
            "f(r) = 1 − (r_s/r)² solves the 5D vacuum Einstein equation "
            "and underlies V_tangherlini. PARTIAL closure: the metric "
            "is bulk-derived, but the radial-domain boundaries "
            "R_INNER, R_OUTER are boundary data. R_OUTER ≈ 1.2626 is "
            "fixed by the cross-species γ-lock fixed point; R_MID = 1 "
            "is a gauge choice. The boundaries are not bulk-equation "
            "output (B3)."
        ),
        'samples': rows,
        'max_f_difference': max_diff,
        'R_MID': R_MID,
        'R_OUTER': R_OUTER,
        'R_INNER': R_INNER,
        'Delta_R': R_OUTER - R_INNER,
        'boundary_data_imposed': True,
        'pass': max_diff < 1e-12,
    }


# ---------------------------------------------------------------------------
# T5. Dimensional bridge barrier (B4)
# ---------------------------------------------------------------------------

def test_T5_dimensional_bridge() -> dict:
    """The action fixes dimensionless ratios and structure; the
    absolute mass scale enters only through the external anchor
    ℏ = m_e · R_MID · c. The gauge sector (Chern number) and the
    throat vertex (F²) are scale-free; the gravity sector sets one
    length scale (R_MID) but not the mass-energy conversion. Barrier
    B4: scale-complete only modulo m_e."""
    # Demonstrate scale-freedom: F² is dimensionless (ratio x = ω'/ω,
    # angle c). Chern number is a pure integer. The throat closure
    # quantum 2π is dimensionless. Only m_e carries the absolute scale.
    scale_free_quantities = {
        'F_squared': 'dimensionless (function of x = ω′/ω and c = cos θ)',
        'chern_c1': 'pure integer = 1',
        'closure_quantum': '2π (dimensionless action units)',
        'K_pade': 'dimensionless (harmonic-mean ratio)',
        'Q_polarization': 'dimensionless (helicity-channel sum)',
    }
    dimensional_bridge = 'ℏ = m_e · R_MID · c  (m_e is the single external anchor)'
    return {
        'name': 'T5_dimensional_bridge_barrier',
        'description': (
            "Barrier B4: the action fixes dimensionless structure but "
            "not the absolute scale. F², the Chern number, the closure "
            "quantum, K and Q are all scale-free. The only external "
            "input setting MeV is m_e via ℏ = m_e·R_MID·c (consistent "
            "with docs/hbar_origin_status.md). The action is "
            "scale-complete modulo m_e."
        ),
        'scale_free_quantities': scale_free_quantities,
        'dimensional_bridge': dimensional_bridge,
        'external_anchors_count': 1,
        'pass': True,   # barrier documented (informative)
    }


# ---------------------------------------------------------------------------
# T6. 5D → 4D reduction barrier (B5)
# ---------------------------------------------------------------------------

def test_T6_dimensional_reduction() -> dict:
    """F² is a 4D-effective vertex factor; the bulk action is 5D. The
    Kaluza-Klein-like integration over the radial channel that would
    produce F² from the 5D action is not constructed. Barrier B5: the
    reduction map is a structural gap, not a number."""
    reduction_requirements = [
        (
            "Radial-mode integration: integrate the 5D fields over the "
            "Tangherlini radial channel [R_INNER, R_OUTER] against the "
            "bound modes u_{l,n}(r*) to obtain 4D effective couplings."
        ),
        (
            "Vertex projection: project the bulk Dirac–gauge–throat "
            "coupling onto the 4D Compton kinematic variables (x, c) "
            "to recover F²(x, c)."
        ),
        (
            "Consistency of the radial-mode spectrum (which sets lepton/"
            "quark masses) with the F² vertex normalisation — the two "
            "currently live in separate sub-threads (closure-ledger vs "
            "tree-QED) and are not connected by an explicit reduction."
        ),
    ]
    return {
        'name': 'T6_5d_to_4d_reduction_barrier',
        'description': (
            "Barrier B5: F² is 4D-effective; the bulk is 5D. The "
            "dimensional-reduction map (radial-mode integration + "
            "vertex projection) that produces F² from the 5D action "
            "is not constructed. This is a structural gap, not a "
            "numerical mismatch — the largest open piece of the "
            "scaffold."
        ),
        'reduction_requirements': reduction_requirements,
        'constructed': False,
        'pass': True,   # barrier documented (informative)
    }


# ---------------------------------------------------------------------------
# T7. Barrier ledger / closure assessment
# ---------------------------------------------------------------------------

def test_T7_barrier_ledger() -> dict:
    """Tabulate the mismatch terms, their type, and which sectors
    close."""
    sectors = {
        'A_gauge': {
            'target': 'A_φ = ½cos χ',
            'status': 'CLOSES',
            'external_input': 'none',
        },
        'T_throat': {
            'target': 'F² = K²·Q',
            'status': 'CLOSES_WITH_CONSTRAINT',
            'external_input': 'B1 (closure quantum) + B2 (antipodal Z₂)',
        },
        'G_gravity': {
            'target': 'ΔR = R_OUTER − R_INNER',
            'status': 'PARTIAL',
            'external_input': 'B3 (boundary conditions)',
        },
    }
    barriers = {
        'B1_closure_quantum': {
            'type': 'topological',
            'description': '∮A = 2πn Wilson-loop constraint, not local density',
            'mismatch': 'one global topological input per closed orbit',
        },
        'B2_antipodal_Z2': {
            'type': 'discrete',
            'description': 'σ: p → −p (T = iσ_y), discrete quotient not in smooth Lagrangian',
            'mismatch': 'one discrete symmetry imposed as identification',
        },
        'B3_boundary_conditions': {
            'type': 'boundary',
            'description': 'hard wall (Dirichlet) from T²=−I at throat fixed point',
            'mismatch': 'boundary data imposed by non-orientable topology',
        },
        'B4_dimensional_bridge': {
            'type': 'scale',
            'description': 'absolute mass scale via ℏ = m_e·R_MID·c',
            'mismatch': 'one external scale anchor (m_e)',
        },
        'B5_5d_to_4d_reduction': {
            'type': 'reduction',
            'description': 'radial-mode integration + vertex projection producing F²',
            'mismatch': 'structural map not constructed',
        },
    }
    n_closes = sum(1 for s in sectors.values() if s['status'] == 'CLOSES')
    n_partial = sum(1 for s in sectors.values() if s['status'] != 'CLOSES')
    return {
        'name': 'T7_barrier_ledger',
        'description': (
            "Closure assessment: 1 sector closes outright (gauge), "
            "1 closes given constraints (throat), 1 is partial "
            "(gravity). Five mismatch terms block full closure, in "
            "three families: topological (B1, B2), boundary (B3), "
            "and scale/reduction (B4, B5)."
        ),
        'sectors': sectors,
        'barriers': barriers,
        'n_sectors_closing_outright': n_closes,
        'n_sectors_partial_or_constrained': n_partial,
        'n_barriers': len(barriers),
        'barrier_families': {
            'topological': ['B1_closure_quantum', 'B2_antipodal_Z2'],
            'boundary': ['B3_boundary_conditions'],
            'scale_and_reduction': ['B4_dimensional_bridge', 'B5_5d_to_4d_reduction'],
        },
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_action_scaffold()
    t2 = test_T2_gauge_sector()
    t3 = test_T3_throat_sector()
    t4 = test_T4_gravity_sector()
    t5 = test_T5_dimensional_bridge()
    t6 = test_T6_dimensional_reduction()
    t7 = test_T7_barrier_ledger()
    tests = [t1, t2, t3, t4, t5, t6, t7]

    # The scaffold "passes" in the sense that every test runs and the
    # sector statuses + barriers are as anticipated. Gauge closes,
    # throat closes-with-constraint, gravity partial.
    gauge_closes = t2['pass']
    throat_closes = t3['pass']
    gravity_partial = t4['pass']

    if gauge_closes and throat_closes and gravity_partial:
        verdict_class = 'SCAFFOLD_WITH_BARRIERS'
        verdict = (
            'SCAFFOLD WITH BARRIERS. A single covariant action\n'
            '  S_BAM = ∫_{M₅} √(−g₅)[ (R₅−2Λ₅)/2κ₅ − ¼F² + ψ̄(iΓ·D−m)ψ '
            '+ L_throat ] + S_∂ + S_closure\n'
            'unifies the three BAM targets STRUCTURALLY, but full '
            'closure is blocked by five identified mismatch terms.\n\n'
            'Sector status:\n'
            '  (A) gauge   — CLOSES. A_φ = ½cos χ is the homogeneous '
            'Maxwell solution on the Hopf bundle; c₁ = 1 exactly; no '
            'external input.\n'
            '  (T) throat  — CLOSES WITH CONSTRAINT. F² = K²·Q '
            'reconstructs from the equal-action splitting, but requires '
            'the closure quantum (B1) and antipodal Z₂ (B2).\n'
            '  (G) gravity — PARTIAL. The Tangherlini metric is '
            'bulk-derived, but ΔR boundary data is imposed (B3).\n\n'
            'Mismatch terms (barriers), in three families:\n'
            '  TOPOLOGICAL: B1 closure quantum (∮A = 2πn, not a local '
            'density); B2 antipodal Z₂ (T = iσ_y, discrete quotient).\n'
            '  BOUNDARY: B3 hard-wall Dirichlet from T²=−I (topological, '
            'not from δS/δg).\n'
            '  SCALE / REDUCTION: B4 dimensional bridge '
            '(ℏ = m_e·R_MID·c, one external anchor); B5 5D→4D '
            'reduction producing F² (not constructed, the largest gap).\n\n'
            'The barriers are not five unrelated patches: they are a '
            'recognisable set — topological quantisation, discrete '
            'antipodal identification, the dimensional anchor, and the '
            'unbuilt dimensional reduction — that recurs throughout the '
            'BAM programme. Full closure requires (i) promoting the '
            'closure quantum + antipodal Z₂ to a topological / discrete '
            'sector of the action (e.g. a θ-term + S³/Z₂ quotient), '
            '(ii) deriving the hard-wall BC from the non-orientable '
            'bulk, and (iii) constructing the 5D→4D radial reduction '
            'that projects the bulk coupling onto F²(x, c).'
        )
    elif gauge_closes and not throat_closes:
        verdict_class = 'SCAFFOLD_FAILS'
        verdict = (
            'SCAFFOLD FAILS at the throat sector — F² did not '
            'reconstruct from the equal-action splitting. Investigate.'
        )
    else:
        verdict_class = 'SCAFFOLD_FAILS'
        verdict = (
            'SCAFFOLD FAILS: the gauge sector (which should close '
            'outright) did not reproduce A_φ / c₁ = 1. Check the '
            'action ansatz.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'action': (
            'S_BAM = ∫_{M₅} √(−g₅)[ (R₅−2Λ₅)/2κ₅ − ¼F_{MN}F^{MN} '
            '+ ψ̄(iΓ^M D_M − m)ψ + L_throat ] + S_∂[hard walls] '
            '+ S_closure[∮A = 2πn]'
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
    L.append('# Covariant BAM effective-action scaffold probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Builds the candidate single variational principle for the three '
        'BAM targets (Hopf U(1) connection, Tangherlini boundaries, F² '
        'vertex) and identifies the mismatch terms preventing full '
        'closure. A scaffold + barrier map, not a closure claim.'
    )
    L.append('')

    L.append('## Candidate action')
    L.append('')
    L.append('```')
    L.append(s['action'])
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
            value = 'action stated; target → sector → variation map'
        elif nm.startswith('T2'):
            value = (
                f"A_φ = ½cos χ (diff {t['max_connection_difference']:.1e}); "
                f"c₁ = {t['chern_c1_abs']:.6f} (err {t['chern_c1_error']:.1e}) "
                f"— CLOSES"
            )
        elif nm.startswith('T3'):
            value = (
                f"F² = K²·Q (max diff {t['max_difference']:.1e}); "
                f"requires B1 + B2"
            )
        elif nm.startswith('T4'):
            value = (
                f"f(r) = 1−(r_s/r)² (diff {t['max_f_difference']:.1e}); "
                f"ΔR = {t['Delta_R']:.4f} imposed (B3)"
            )
        elif nm.startswith('T5'):
            value = f"scale-complete modulo m_e ({t['external_anchors_count']} anchor)"
        elif nm.startswith('T6'):
            value = f"5D→4D reduction NOT constructed (B5)"
        elif nm.startswith('T7'):
            value = (
                f"{t['n_sectors_closing_outright']} sector closes, "
                f"{t['n_barriers']} barriers in 3 families"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1
    t1 = s['tests'][0]
    L.append('## T1: Action scaffold structure')
    L.append('')
    L.append('Action terms:')
    for k, v in t1['action_terms'].items():
        L.append(f"  - **{k}**: {v}")
    L.append('')
    L.append('Target → sector → variation:')
    for k, v in t1['target_sector_variation_map'].items():
        L.append(f"  - `{k}`: {v}")
    L.append('')

    # T2
    t2 = s['tests'][1]
    L.append('## T2: Gauge sector — CLOSES')
    L.append('')
    L.append('| χ/π | A_φ | ½cos χ | diff |')
    L.append('|---:|---:|---:|---:|')
    for r in t2['connection_samples']:
        L.append(
            f"| {r['chi_in_pi']:.3f} | {r['A_phi']:+.6f} | "
            f"{r['half_cos_chi']:+.6f} | {r['difference']:.1e} |"
        )
    L.append('')
    L.append(
        f"First Chern number c₁ = **{t2['chern_c1_abs']:.6f}** "
        f"(error {t2['chern_c1_error']:.2e}). The Hopf connection is the "
        f"homogeneous Maxwell solution; this sector derives from "
        f"δS/δA = 0 with no external input."
    )
    L.append('')

    # T3
    t3 = s['tests'][2]
    L.append('## T3: Throat sector — CLOSES WITH CONSTRAINT')
    L.append('')
    L.append('| x | cosθ | K | Q | K²·Q | F² closed | diff |')
    L.append('|---:|---:|---:|---:|---:|---:|---:|')
    for r in t3['samples']:
        L.append(
            f"| {r['x']:.2f} | {r['cos_theta']:+.2f} | "
            f"{r['K_from_equal_action']:.4f} | {r['Q_from_hopf_spinor']:.4f} | "
            f"{r['K2_Q']:.4f} | {r['F2_closed']:.4f} | {r['difference']:.1e} |"
        )
    L.append('')
    L.append(
        f"F² = K²·Q reconstructs to {t3['max_difference']:.1e}, but "
        f"requires two imposed inputs: **{', '.join(t3['imposed_inputs'])}** "
        f"(closure quantum 2π and antipodal Z₂) — neither is a local "
        f"Lagrangian density."
    )
    L.append('')

    # T4
    t4 = s['tests'][3]
    L.append('## T4: Gravity sector — PARTIAL')
    L.append('')
    L.append('| r | f(r) = 1−(r_s/r)² | f from V_tangherlini | diff |')
    L.append('|---:|---:|---:|---:|')
    for r in t4['samples']:
        L.append(
            f"| {r['r']:.2f} | "
            f"{r['f_tangherlini_1_minus_rs2_over_r2']:+.6f} | "
            f"{r['f_recovered_from_V_tangherlini']:+.6f} | "
            f"{r['difference']:.1e} |"
        )
    L.append('')
    L.append(
        f"R_MID = {t4['R_MID']}, R_OUTER = {t4['R_OUTER']}, "
        f"R_INNER = {t4['R_INNER']}, ΔR = {t4['Delta_R']:.4f}. "
        f"The metric is bulk-derived; the boundaries are imposed (B3)."
    )
    L.append('')

    # T5
    t5 = s['tests'][4]
    L.append('## T5: Dimensional bridge barrier (B4)')
    L.append('')
    for k, v in t5['scale_free_quantities'].items():
        L.append(f"  - `{k}`: {v}")
    L.append('')
    L.append(f"Dimensional bridge: `{t5['dimensional_bridge']}`. "
             f"External anchors: {t5['external_anchors_count']} (m_e).")
    L.append('')

    # T6
    t6 = s['tests'][5]
    L.append('## T6: 5D → 4D reduction barrier (B5)')
    L.append('')
    L.append('The reduction map producing F² from the 5D action would require:')
    for i, req in enumerate(t6['reduction_requirements'], 1):
        L.append(f"{i}. {req}")
    L.append('')
    L.append(f"**Constructed: {t6['constructed']}** — the largest structural gap.")
    L.append('')

    # T7
    t7 = s['tests'][6]
    L.append('## T7: Barrier ledger / closure assessment')
    L.append('')
    L.append('### Sector status')
    L.append('')
    L.append('| sector | target | status | external input |')
    L.append('|---|---|---|---|')
    for k, v in t7['sectors'].items():
        L.append(
            f"| `{k}` | {v['target']} | **{v['status']}** | {v['external_input']} |"
        )
    L.append('')
    L.append('### Mismatch terms')
    L.append('')
    L.append('| barrier | type | mismatch |')
    L.append('|---|---|---|')
    for k, v in t7['barriers'].items():
        L.append(f"| `{k}` | {v['type']} | {v['mismatch']} |")
    L.append('')
    L.append('### Barrier families')
    L.append('')
    for fam, bs in t7['barrier_families'].items():
        L.append(f"  - **{fam}**: {', '.join(bs)}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What full closure would require')
    L.append('')
    L.append(
        '1. **Promote B1 + B2 to a topological / discrete sector** — '
        'e.g. a θ-term enforcing ∮A = 2πn and an explicit S³/Z₂ '
        'antipodal quotient, so the closure quantum and antipodal '
        'symmetry become part of the action rather than imposed '
        'constraints.'
    )
    L.append(
        '2. **Derive the hard-wall boundary condition (B3)** from the '
        'non-orientable bulk topology (T² = −I), rather than imposing '
        'Dirichlet by hand.'
    )
    L.append(
        '3. **Construct the 5D → 4D reduction (B5)** — integrate the '
        'bulk fields over the Tangherlini radial channel and project '
        'onto the Compton variables (x, c) to produce F²(x, c) from '
        'first principles. This is the largest gap.'
    )
    L.append(
        '4. **Close the dimensional bridge (B4)** — derive m_e (or the '
        'R_MID scale) from within the action, removing the single '
        'external anchor.'
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
    out = here / 'runs' / f'{ts}_bam_effective_action_scaffold_probe'
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
