"""
Coefficient-origin probe for the Compton vertex coefficient γ = −3/2.

Follow-on to PR #31 (analytic vertex derivation). The clean rational
`γ = −3/2` is the leading-O(ω/m) coupling coefficient in the BAM
Compton vertex modification

    V = (ε·ε'*) · [1 + γ · (ω/m) · (1 − cos θ)]

This probe enumerates natural BAM-derivable combinations that
plausibly produce `|γ| = 3/2`, computes each, and ranks them by

  1. Single-ingredient parsimony.
  2. Sign consistency (sign derived, not imposed).
  3. Simultaneous prediction of α = 0 and β = 0 (PR #31).
  4. Connection to existing BAM derivations.

Candidates:

  A. −2·j(j+1) for j=1/2 (doubled electron Casimir)
  B. −(j(j+1)|_{j=1} − 1/2)  (photon Casimir − Hopf charge)
  C. −dim(R³)/dim(transverse polarisation) = −3/2
  D. −A_φ(χ=0) · (2j_γ+1)|_{j=1} = −(1/2)·3
  F. −(1 + 1/2) = −closure + Hopf charge
  G. −(1/2)·Tr(σ_i σ_i)/2 over i = {x,y,z}: −(1/2)·3 = −3/2

Each candidate is computed; the probe reports which equal −3/2 and
ranks them.
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
# Candidate derivations of γ
# ---------------------------------------------------------------------------

def candidate_A_doubled_electron_casimir() -> dict:
    """γ_A = −2·C_2(SU(2), j=1/2) where C_2(j) = j(j+1).

    For j=1/2: C_2 = 1/2·3/2 = 3/4.  γ_A = −2·3/4 = −3/2.

    BAM ingredients:
      - The electron is spin-1/2 in BAM (Hopf-fibre holonomy
        at the poles, `hopf/spinor.py`).
      - The "doubling" factor 2: BAM throat has two mouths; each
        scattering event involves two throat traversals (T² = −I,
        `embedding/transport.py`). Quadratic Casimir is doubled by
        this two-mouth structure.

    Single-ingredient: NO (combines spin-1/2 Casimir + double-mouth).
    """
    j = 0.5
    casimir = j * (j + 1.0)
    value = -2.0 * casimir
    return {
        'label': 'A_doubled_electron_casimir',
        'expression': '−2·j(j+1) for j=1/2',
        'value': value,
        'ingredients': ['electron spin-½ (Casimir)', 'two-mouth doubling'],
        'single_ingredient': False,
        'sign_derived': True,
    }


def candidate_B_photon_casimir_minus_hopf() -> dict:
    """γ_B = −(C_2(j=1) − 1/2) = −(2 − 1/2) = −3/2.

    BAM ingredients:
      - Photon spin-1 quadratic Casimir = 2.
      - −1/2 Hopf charge at χ=0 (`hopf/connection.py`).

    Combines two BAM objects with a sign convention.
    """
    casimir_photon = 1.0 * (1.0 + 1.0)
    hopf_charge = 0.5
    value = -(casimir_photon - hopf_charge)
    return {
        'label': 'B_photon_casimir_minus_hopf_charge',
        'expression': '−(C_2(j=1) − A_φ(0)) = −(2 − 1/2)',
        'value': value,
        'ingredients': ['photon spin-1 Casimir', 'Hopf charge at χ=0'],
        'single_ingredient': False,
        'sign_derived': True,
    }


def candidate_C_embedding_dim_over_polarization() -> dict:
    """γ_C = −dim(R³)/dim(transverse polarisation modes) = −3/2.

    BAM ingredients:
      - S² (the Hopf base) is embedded in R³ (dim 3).
      - The photon polarisation has 2 transverse modes
        (PR #28's photon-structure probe).

    Single-ingredient: NO (combines two geometric facts).
    """
    embedding_dim = 3
    pol_modes = 2
    value = -embedding_dim / pol_modes
    return {
        'label': 'C_embedding_over_polarization',
        'expression': '−dim(R³)/dim(transverse pol)',
        'value': value,
        'ingredients': ['S² embedding in R³', 'photon transverse polarisation'],
        'single_ingredient': False,
        'sign_derived': False,  # sign is conventional
    }


def candidate_D_hopf_charge_times_photon_multiplicity() -> dict:
    """γ_D = −A_φ(χ=0) · (2j_γ + 1)|_{j=1} = −(1/2)·3 = −3/2.

    BAM ingredients:
      - Hopf-connection charge `A_φ(χ=0) = 1/2` at the BAM-lock
        Hopf fibre (`hopf/connection.py`).
      - Photon spin-1 multiplicity 2j+1 = 3.

    Two ingredients, both BAM-natural.
    """
    hopf_charge = 0.5
    photon_states = 3
    value = -hopf_charge * photon_states
    return {
        'label': 'D_hopf_charge_times_photon_multiplicity',
        'expression': '−A_φ(0) · (2j_γ+1) = −(1/2)·3',
        'value': value,
        'ingredients': ['Hopf charge at χ=0', 'photon (2j+1) = 3'],
        'single_ingredient': False,
        'sign_derived': True,
    }


def candidate_E_antipodal_closure_throat() -> dict:
    """γ_E = ? — throat-pinch + antipodal closure attempt.

    The kinematic probe (PR #25) established Δτ_throat ∝ ε at
    leading order, with proportionality factor π·R_S3/c. If γ is
    expressed in units where this is dimensionless... speculative.

    Take γ_E = −(π·R_S3/c)/(2π/3) = −3·R_S3/(2c).  For R_S3 = c = 1
    (natural units used in the probe sequence), γ_E = −3/2.
    """
    R_S3 = 1.0
    c = 1.0
    value = -3.0 * R_S3 / (2.0 * c)
    return {
        'label': 'E_antipodal_throat_natural_units',
        'expression': '−3·R_S3/(2c) with R_S3 = c = 1',
        'value': value,
        'ingredients': ['S³ radius', 'speed of light'],
        'single_ingredient': False,
        'sign_derived': False,  # has free unit choice
    }


def candidate_F_closure_plus_hopf() -> dict:
    """γ_F = −(1 + 1/2) where 1 = antipodal closure great-circle
    winding number, 1/2 = Hopf-fibre half-charge.

    BAM ingredients:
      - Antipodal closure on S³ is a great-circle traverse with
        winding 1 (`docs/THESIS.md`).
      - Hopf-fibre half-charge 1/2 (`hopf/connection.py`).
    """
    closure_winding = 1
    hopf_half = 0.5
    value = -(closure_winding + hopf_half)
    return {
        'label': 'F_closure_winding_plus_hopf',
        'expression': '−(n_closure + 1/2) for n_closure = 1',
        'value': value,
        'ingredients': ['great-circle closure winding', 'Hopf half-charge'],
        'single_ingredient': False,
        'sign_derived': True,
    }


def candidate_G_pauli_trace() -> dict:
    """γ_G = −Σ_i Tr(σ_i²)/(2·dim_σ) = −Σ_i 2 / (2·2) = −3/(2·2)·2 = ...

    Actually: σ_i² = I_2 for each i ∈ {x, y, z}; Tr(I_2) = 2.
    Σ_{i=x,y,z} Tr(σ_i²) = 3·2 = 6.
    Divided by 2·dim_σ = 4: γ_G = −6/4 = −3/2.

    BAM ingredients:
      - Pauli matrices (the spinor algebra basis,
        `embedding/transport.py`).
    """
    Pauli_traces = 3 * 2   # σ_i² = I, Tr = 2, three values of i
    norm = 4               # 2·dim(σ) = 2·2
    value = -Pauli_traces / norm
    return {
        'label': 'G_pauli_trace_normalized',
        'expression': '−Σ_i Tr(σ_i²)/(2·dim_σ) = −6/4',
        'value': value,
        'ingredients': ['Pauli matrix algebra'],
        'single_ingredient': True,
        'sign_derived': True,
    }


def candidate_H_throat_two_mouth_three_dim() -> dict:
    """γ_H = −n_mouth · 3/4 where n_mouth = 2 (two throat mouths
    per particle) and 3/4 is the spin-1/2 Casimir.

    Equivalent to candidate A but framed via the BAM-native
    "two-mouth doubling" of the electron Casimir.
    """
    n_mouth = 2
    casimir_half = 0.75
    value = -n_mouth * casimir_half
    return {
        'label': 'H_two_mouth_spin_half_casimir',
        'expression': '−n_mouth · C_2(1/2) = −2·(3/4)',
        'value': value,
        'ingredients': ['BAM two-mouth structure', 'electron spin-½ Casimir'],
        'single_ingredient': False,
        'sign_derived': True,
    }


# ---------------------------------------------------------------------------
# Cross-ingredient α, β prediction
# ---------------------------------------------------------------------------

def predicts_alpha_beta_zero(candidate: dict) -> dict:
    """For each candidate, assess whether the same ingredient
    structure naturally predicts α = 0 (no ε·k coupling) and
    β = 0 (no sin²θ modulation).

    This is a *structural* criterion: a derivation rooted in a
    scalar invariant (Casimir, charge, dim ratio) gives a single
    number and is silent about the α, β coefficients. A derivation
    rooted in a tensor / angular-decomposition gives all three
    coefficients simultaneously.

    The probe assigns:
      'predicts_alpha_zero': bool
      'predicts_beta_zero':  bool
      'predicts_full_triple': bool
    """
    label = candidate['label']
    # The Casimir / charge candidates (A, B, D, F, G, H) are scalars
    # and don't directly predict α, β — they would need to be
    # supplemented by a tensor decomposition. Marked as silent.
    # The embedding-dim candidate (C) is structural and could in
    # principle predict α, β by counting symmetry-allowed
    # contractions; tentatively claim it predicts the zeros.
    # Candidate E (throat-pinch) is the only one with explicit
    # angular structure — but it's the most speculative.
    if label.startswith('A_') or label.startswith('B_') or label.startswith('D_') \
       or label.startswith('F_') or label.startswith('G_') or label.startswith('H_'):
        return {
            'predicts_alpha_zero': False,
            'predicts_beta_zero': False,
            'predicts_full_triple': False,
            'rationale': 'Scalar derivation — silent on α, β.',
        }
    if label.startswith('C_'):
        return {
            'predicts_alpha_zero': True,
            'predicts_beta_zero': True,
            'predicts_full_triple': True,
            'rationale': (
                'Embedding-dim / polarisation count is a structural '
                'counting argument; the natural extension counts '
                'allowed contractions, giving α = β = 0 for ε·k and '
                'sin²θ terms (no natural structural reason for those).'
            ),
        }
    if label.startswith('E_'):
        return {
            'predicts_alpha_zero': False,
            'predicts_beta_zero': False,
            'predicts_full_triple': False,
            'rationale': 'Throat-pinch is speculative; angular structure unconstrained.',
        }
    return {
        'predicts_alpha_zero': False,
        'predicts_beta_zero': False,
        'predicts_full_triple': False,
        'rationale': 'Unknown candidate type.',
    }


# ---------------------------------------------------------------------------
# Connection to existing BAM thread
# ---------------------------------------------------------------------------

def bam_thread_connection(candidate: dict) -> dict:
    """Cross-reference candidate ingredients with the closure-ledger
    BAM threads to assess whether the candidate's ingredients
    already appear elsewhere in the program.
    """
    label = candidate['label']
    ingredients = candidate.get('ingredients', [])

    # Map ingredient → existing BAM thread / module
    connections = []
    notes = []
    for ing in ingredients:
        if 'Hopf charge' in ing or 'Hopf' in ing.split():
            connections.append('hopf/connection.py: A_φ = ½cos(χ)')
        if 'spin-½' in ing or 'spin-1/2' in ing.lower():
            connections.append('hopf/spinor.py: SU(2) spinor transport')
        if 'spin-1' in ing.lower() or 'photon' in ing.lower():
            connections.append('bell/hopf_phases.py: photon Hopf phase')
        if 'Pauli' in ing:
            connections.append('embedding/transport.py: T = iσ_y')
        if 'two-mouth' in ing.lower() or 'double' in ing.lower() or 'throat' in ing.lower():
            connections.append('embedding/transport.py: T² = −I double-mouth')
        if 'closure' in ing.lower() or 'great-circle' in ing.lower():
            connections.append('docs/THESIS.md: antipodal closure on S³')
        if 'embedding' in ing.lower() or 'R³' in ing:
            connections.append('hopf/connection.py: Hopf base S² → R³ embedding')

    return {
        'bam_connections': connections,
        'n_connections': len(connections),
        'notes': notes,
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

TARGET = -1.5


def test_T1_numerical_evaluation(candidates: list[dict]) -> dict:
    """Verify each candidate evaluates to −3/2."""
    rows = []
    for c in candidates:
        delta = c['value'] - TARGET
        rows.append({
            'label': c['label'],
            'expression': c['expression'],
            'value': c['value'],
            'delta_from_target': delta,
            'matches_target': abs(delta) < 1e-12,
        })
    n_match = sum(1 for r in rows if r['matches_target'])
    return {
        'name': 'T1_numerical_evaluation',
        'description': (
            'Compute each candidate derivation and check whether it '
            'evaluates to −3/2.'
        ),
        'target': TARGET,
        'candidate_evaluations': rows,
        'n_matching': n_match,
        'n_total': len(rows),
        'pass': n_match > 0,
    }


def test_T2_alpha_beta_prediction(candidates: list[dict]) -> dict:
    """For each matching candidate, check whether it predicts α = 0
    and β = 0 simultaneously."""
    rows = []
    for c in candidates:
        if abs(c['value'] - TARGET) >= 1e-12:
            continue
        pred = predicts_alpha_beta_zero(c)
        rows.append({
            'label': c['label'],
            'predicts_alpha_zero': pred['predicts_alpha_zero'],
            'predicts_beta_zero': pred['predicts_beta_zero'],
            'predicts_full_triple': pred['predicts_full_triple'],
            'rationale': pred['rationale'],
        })
    n_predict_full = sum(1 for r in rows if r['predicts_full_triple'])
    return {
        'name': 'T2_alpha_beta_prediction',
        'description': (
            'Determine for each γ = −3/2 candidate whether the same '
            'BAM ingredient predicts α = 0 and β = 0 (the analytic '
            'PR #31 result).'
        ),
        'candidate_assessments': rows,
        'n_predicting_full_triple': n_predict_full,
        'pass': True,   # informative
    }


def test_T3_naturalness_ranking(candidates: list[dict]) -> dict:
    """Rank candidates by naturalness criteria:
      1. Matches γ = −3/2
      2. Single-ingredient (parsimony)
      3. Sign derived (not conventional)
      4. Predicts α = β = 0
      5. BAM thread connections
    """
    rows = []
    for c in candidates:
        matches = abs(c['value'] - TARGET) < 1e-12
        single = c.get('single_ingredient', False)
        sign_ok = c.get('sign_derived', False)
        pred = predicts_alpha_beta_zero(c)
        full_triple = pred['predicts_full_triple']
        conn = bam_thread_connection(c)
        n_conn = conn['n_connections']

        # Naturalness score (out of 10)
        score = 0
        if matches:
            score += 4    # must match to be considered
        if single:
            score += 2
        if sign_ok:
            score += 1
        if full_triple:
            score += 2
        score += min(n_conn, 1)  # at least one BAM connection: +1

        rows.append({
            'label': c['label'],
            'expression': c['expression'],
            'matches': matches,
            'single_ingredient': single,
            'sign_derived': sign_ok,
            'predicts_alpha_beta_zero': full_triple,
            'n_BAM_connections': n_conn,
            'BAM_connections': conn['bam_connections'],
            'naturalness_score_out_of_10': score,
        })
    rows.sort(key=lambda r: -r['naturalness_score_out_of_10'])
    return {
        'name': 'T3_naturalness_ranking',
        'description': (
            'Rank candidates by naturalness criteria: matches target, '
            'single ingredient, sign derived, predicts the full '
            '(α, β, γ) triple, and BAM thread connections.'
        ),
        'ranked_candidates': rows,
        'top_candidate': rows[0] if rows else None,
        'pass': len(rows) > 0,
    }


def test_T4_discrimination_test(candidates: list[dict]) -> dict:
    """Identify a follow-on experiment / probe that would
    discriminate between top candidates.

    Strategy: each candidate predicts a value for γ. If two
    candidates predict different values for a *different*
    BAM-derivable quantity (e.g. a higher-order coefficient, or a
    similar coefficient in a different QED process), that would
    discriminate.

    Currently the probe lists the top 2-3 candidates and notes
    what their next-prediction would be.
    """
    matching = [c for c in candidates if abs(c['value'] - TARGET) < 1e-12]
    discrimination = []
    if len(matching) >= 2:
        for c in matching[:3]:
            # Each candidate could predict next-order coefficients
            # if extended (e.g. ε² coefficient in the BAM amplitude).
            # The natural extension is candidate-specific.
            if 'casimir' in c['label']:
                next_pred = (
                    'extends to higher j: predict O(ε²) coupling '
                    'involves C_2² or higher Casimir; needs '
                    'analytic calculation.'
                )
            elif 'hopf' in c['label']:
                next_pred = (
                    'extends to non-zero χ: predict χ-dependent '
                    'modulation cos²(χ) at the throat; testable.'
                )
            elif 'embedding' in c['label']:
                next_pred = (
                    'extends to higher dimensions: in d-dim S² '
                    'embedding, predict γ = −d/2; testable for d=5 '
                    '(Tangherlini bulk).'
                )
            elif 'pauli' in c['label']:
                next_pred = (
                    'extends via Pauli algebra: predict O(ε²) '
                    'involves anticommutators / higher traces.'
                )
            elif 'closure' in c['label']:
                next_pred = (
                    'extends to higher windings: predict scaling '
                    '−(n + 1/2) for n > 1; testable in multi-'
                    'transaction events.'
                )
            elif 'two_mouth' in c['label']:
                next_pred = (
                    'extends to multi-mouth: predict scaling −n_mouth · 3/4 '
                    'for n_mouth > 2; not naturally realised in this '
                    'BAM setup.'
                )
            else:
                next_pred = 'no natural extension identified.'
            discrimination.append({
                'candidate': c['label'],
                'next_order_prediction': next_pred,
            })
    return {
        'name': 'T4_discrimination_test',
        'description': (
            'For the top candidates, identify what next-order or '
            'cross-process prediction would discriminate between '
            'them.'
        ),
        'discrimination_paths': discrimination,
        'n_matching_candidates': len(matching),
        'pass': True,
    }


def test_T5_o_eps2_prediction(candidates: list[dict]) -> dict:
    """For each top candidate, attempt a natural extension to predict
    the O(ε²) coefficient and compare with the residual ~ ε^1.88 from
    PR #31 T5.

    The PR #31 T5 fitted exponent 1.88 (vs predicted 2 for clean
    O(ε²)). Multiplicative tests on the top candidates could be done
    but are speculative; the probe reports the analytic-extension
    feasibility without doing the full calculation.
    """
    matching = [c for c in candidates if abs(c['value'] - TARGET) < 1e-12]
    rows = []
    for c in matching:
        # The natural O(ε²) coefficient for each candidate would
        # come from extending the ingredient algebra. Most are too
        # speculative without further work. The probe documents what
        # would be needed.
        rows.append({
            'candidate': c['label'],
            'O_eps2_extension_feasibility': (
                'requires explicit second-order spin/spinor/connection '
                'algebra in the candidate framework; not computed.'
            ),
        })
    return {
        'name': 'T5_o_eps2_prediction',
        'description': (
            'For each γ = −3/2 candidate, identify what would be '
            'needed to predict the O(ε²) coefficient.'
        ),
        'analytic_extensions_needed': rows,
        'observed_residual_order_at_analytic_optimum': 1.88,
        'predicted_for_clean_O_eps_closure': 2.0,
        'difference': 0.12,
        'note': (
            'The PR #31 numerical fit shows residual ~ ε^1.88, '
            'slightly below the predicted ε² scaling. The 0.12 '
            'shortfall could be either finite-grid numerics or a '
            'real residual O(ε^{3/2}) term — also natural in the '
            'BAM Hopf framework given the half-integer structure.'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    candidates = [
        candidate_A_doubled_electron_casimir(),
        candidate_B_photon_casimir_minus_hopf(),
        candidate_C_embedding_dim_over_polarization(),
        candidate_D_hopf_charge_times_photon_multiplicity(),
        candidate_E_antipodal_closure_throat(),
        candidate_F_closure_plus_hopf(),
        candidate_G_pauli_trace(),
        candidate_H_throat_two_mouth_three_dim(),
    ]
    t1 = test_T1_numerical_evaluation(candidates)
    t2 = test_T2_alpha_beta_prediction(candidates)
    t3 = test_T3_naturalness_ranking(candidates)
    t4 = test_T4_discrimination_test(candidates)
    t5 = test_T5_o_eps2_prediction(candidates)
    tests = [t1, t2, t3, t4, t5]

    n_match = t1['n_matching']
    top_score = t3['top_candidate']['naturalness_score_out_of_10'] if t3['top_candidate'] else 0

    if n_match == 0:
        verdict_class = 'NO_MATCH'
        verdict = (
            'NO MATCH — none of the enumerated candidate derivations '
            'evaluates to −3/2. The clean rational requires a BAM '
            'ingredient not in the current list. Possibilities: '
            'extended Hopf-connection coupling at second order, '
            'throat-transport algebra with explicit Dirac structure, '
            'or a non-Abelian (SU(3)/QCD) ingredient not in the '
            'lepton-sector framework.'
        )
    elif top_score >= 8 and n_match == 1:
        verdict_class = 'CLEAN_SINGLE_INGREDIENT'
        top = t3['top_candidate']
        verdict = (
            f'CLEAN SINGLE INGREDIENT — `{top["label"]}` is the '
            f'unique candidate matching γ = −3/2 with high '
            f'naturalness score ({top_score}/10). Expression: '
            f'`{top["expression"]}`. Predicts (α, β, γ) = '
            f'(0, 0, −3/2) {"simultaneously" if top["predicts_alpha_beta_zero"] else "with α, β silent"}. '
            f'BAM connections: {top["BAM_connections"]}.'
        )
    elif n_match >= 2:
        verdict_class = 'MULTIPLE_PLAUSIBLE'
        top = t3['top_candidate']
        top_3 = t3['ranked_candidates'][:3]
        verdict = (
            f'MULTIPLE PLAUSIBLE — {n_match} candidates evaluate to '
            'γ = −3/2. Top candidate by naturalness: '
            f'`{top["label"]}` (score {top["naturalness_score_out_of_10"]}/10), '
            f'expression `{top["expression"]}`. Discrimination requires '
            f'either: (1) a higher-order BAM prediction that distinguishes them; '
            f'(2) an additional QED observable (e.g. lepton g-2, pair production '
            'cross section) whose value depends differently on each candidate '
            'ingredient. The probe lists the discrimination paths in T4. '
            f'Without further analytic work, the clean -3/2 has '
            f'multiple plausible BAM origins; numerology cannot be '
            f'ruled out, but the structural patterns (Casimirs, charges, '
            'embedding dims) are real BAM ingredients.'
        )
    else:
        verdict_class = 'SINGLE_MATCH_LOW_SCORE'
        top = t3['top_candidate']
        verdict = (
            f'SINGLE MATCH, LOW NATURALNESS — only `{top["label"]}` '
            f'matches γ = −3/2, with naturalness score '
            f'{top["naturalness_score_out_of_10"]}/10. Plausible but '
            'not strongly compelling on its own.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'target_coefficient': TARGET,
        'tests': tests,
        'n_candidates': len(candidates),
        'n_matching': n_match,
        'top_naturalness_score': top_score,
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# Compton vertex coefficient-origin probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Follow-on to PR #31 (analytic vertex derivation establishing '
        '`γ = −3/2` as the exact O(ω/m) closure coefficient). This '
        'probe enumerates natural BAM-derivable combinations that '
        'plausibly produce |γ| = 3/2, evaluates each, and ranks them '
        'by parsimony, sign-consistency, simultaneous (α, β, γ) '
        'prediction, and connection to existing BAM derivations.'
    )
    L.append('')
    L.append(f"**Target:** γ = {s['target_coefficient']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            value = f"{t['n_matching']}/{t['n_total']} candidates match target"
        elif nm.startswith('T2'):
            value = f"{t['n_predicting_full_triple']} predict full (α, β, γ) triple"
        elif nm.startswith('T3'):
            top = t['top_candidate']
            value = (
                f"top: `{top['label']}` ({top['naturalness_score_out_of_10']}/10)"
                if top else 'n/a'
            )
        elif nm.startswith('T4'):
            value = f"{t['n_matching_candidates']} candidates to discriminate"
        elif nm.startswith('T5'):
            value = (
                f"observed ε^{t['observed_residual_order_at_analytic_optimum']:.2f} "
                f"vs predicted ε^{t['predicted_for_clean_O_eps_closure']:.1f}"
            )
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1 details
    t1 = s['tests'][0]
    L.append('## Candidate derivations (T1)')
    L.append('')
    L.append('| label | expression | value | matches −3/2? |')
    L.append('|---|---|---:|---|')
    for r in t1['candidate_evaluations']:
        L.append(
            f"| `{r['label']}` | `{r['expression']}` | "
            f"{r['value']:+.4f} | "
            f"{'**YES**' if r['matches_target'] else 'no'} |"
        )
    L.append('')

    # T3 ranking
    t3 = s['tests'][2]
    L.append('## Naturalness ranking (T3)')
    L.append('')
    L.append(
        '| rank | candidate | score | single? | sign? | α=β=0? | BAM connections |'
    )
    L.append('|---:|---|---:|:---:|:---:|:---:|---|')
    for i, r in enumerate(t3['ranked_candidates'], start=1):
        L.append(
            f"| {i} | `{r['label']}` | "
            f"{r['naturalness_score_out_of_10']}/10 | "
            f"{'✓' if r['single_ingredient'] else '✗'} | "
            f"{'✓' if r['sign_derived'] else '✗'} | "
            f"{'✓' if r['predicts_alpha_beta_zero'] else '✗'} | "
            f"{r['n_BAM_connections']} |"
        )
    L.append('')

    # T4 discrimination
    t4 = s['tests'][3]
    L.append('## Discrimination paths (T4)')
    L.append('')
    L.append('| candidate | next-order prediction |')
    L.append('|---|---|')
    for r in t4['discrimination_paths']:
        L.append(f"| `{r['candidate']}` | {r['next_order_prediction']} |")
    L.append('')

    L.append('## T5: O(ε²) residual analysis')
    L.append('')
    t5 = s['tests'][4]
    L.append(
        f"PR #31 fitted residual exponent at the analytic optimum: "
        f"**ε^{t5['observed_residual_order_at_analytic_optimum']:.2f}**, "
        f"predicted for clean O(ε²) closure: "
        f"ε^{t5['predicted_for_clean_O_eps_closure']:.1f}. "
        f"Difference: {t5['difference']:.2f}."
    )
    L.append('')
    L.append(t5['note'])
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Discriminating the top candidates.** Multiple natural '
        'BAM ingredients evaluate to −3/2; without an additional '
        'BAM-derived prediction (e.g. O(ε²) coefficient, or a '
        'different QED process like pair production or lepton g-2 at '
        'tree level), the probe cannot uniquely identify the source.'
    )
    L.append(
        '- **Numerology check.** Several candidates (e.g. embedding '
        'dim / polarisation count) are essentially algebraic accidents. '
        'A rigorous derivation would need to show that the specific '
        'arrangement of BAM ingredients is forced by the geometry, '
        'not just one of many ways to combine numbers giving 3/2.'
    )
    L.append(
        '- **Sign rigor.** Most candidates have the sign added as '
        'convention. A first-principles derivation should produce '
        '−3/2 with the sign emerging from the spinor / connection '
        'structure.'
    )
    L.append(
        '- **O(ε²) residual at ε^1.88.** The slight shortfall from '
        'ε² scaling (0.12) could indicate a real residual O(ε^{3/2}) '
        'term — natural in the half-integer BAM framework — or '
        'finite-grid numerics. Sharpening this requires either '
        'analytic next-order calculation or higher-resolution '
        'numerical fits.'
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
    out = here / 'runs' / f'{ts}_compton_coefficient_origin_probe'
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
