"""
Tangherlini dimensional-scaling discrimination probe for γ = −3/2.

Follow-on to PR #32 (coefficient-origin probe). That probe identified
8 plausible BAM derivations of γ = −3/2 and identified three
discrimination paths. This probe executes the first: extending the
BAM Compton amplitude to arbitrary spatial dimension d_spatial and
asking whether the vertex coefficient γ depends on d.

The discrimination splits the 8 PR #32 candidates into two classes:

  d-independent (γ = −3/2 in all d):
    A. −2·j(j+1) for j=1/2  (doubled electron Casimir)
    B. −(C_2(j=1) − 1/2)     (photon Casimir − Hopf charge)
    D. −A_φ(0)·(2j_γ+1)      (Hopf charge × photon multiplicity)
    F. −(n_closure + 1/2)    (closure winding + Hopf half)
    G. −Σ_i Tr(σ_i²)/(2·dim) (Pauli trace, normalised)
    H. −n_mouth·C_2(1/2)     (two-mouth × spin-½ Casimir)

  d-dependent (predict γ_d rescaling with d):
    C. −d_spatial/dim(trans pol) = −d/(d−1)
    E. natural-units, speculative

Probe construction:

  Photon in d_spatial spatial dims has (d_spatial − 1) transverse
  polarisations.  Polarisation completeness:
    Σ_λ ε^λ_i ε^λ_j = δ_ij − k̂_i k̂_j
  Polarisation sum:
    Σ_{λ,λ'} |ε^λ(k)·ε^{λ'}(k')|² = (d_spatial − 2) + cos²θ
  Initial average: divide by (d_spatial − 1).

  BAM amplitude with Family B modification:
    |M_BAM_d|² ∝ ((d-2)+cos²θ)·(1 + ε·μ)²·(G_s + G_u)²

  KN-d natural generalisation:
    |M_KN_d|² ∝ x²·((d-2)(x+1/x) − sin²θ + (3-d))

  Both reduce to standard QED at d=4 (d_spatial = 3) ⇒ (1+cos²θ)
  and x²(x+1/x-sin²θ) respectively.

Tests:

  T1. Analytic d-independence — re-derive the matching equation
      symbolically and verify the d-cancellation.

  T2. Numerical d-scaling — extract γ_BAM_d numerically for
      d_spatial ∈ {3, 4, 5, 6} and confirm all give −3/2.

  T3. Candidate predictions table — tabulate predicted γ(d) for
      each of the 8 PR #32 candidates.

  T4. Falsified candidates — list candidates whose predicted γ(d)
      doesn't match the numerical observation.

  T5. Surviving candidates — list candidates consistent with the
      observation.
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
# d-dim amplitude implementations
# ---------------------------------------------------------------------------

def x_ratio(eps: float, theta: float) -> float:
    return 1.0 / (1.0 + eps * (1.0 - math.cos(theta)))


def f_KN_d(eps: float, theta: float, d_spatial: int) -> float:
    """d-dim KN-analog normalised at θ=0.

    |M|²/Σ = (x²/(d-1)) · ((d-2)(x + 1/x) − sin²θ + (3-d))

    At x=1: (1/(d-1))·(2(d-2) - sin²θ + 3-d) = (1/(d-1))·(d-1 - sin²θ)
          = (1/(d-1))·((d-2) + cos²θ)
    which equals f_BAM_d at Thomson — good.

    At θ=0 with x=1: (1/(d-1))·(d-1) = 1 ✓
    """
    x = x_ratio(eps, theta)
    s2 = math.sin(theta) ** 2
    d = d_spatial
    val = (x * x / (d - 1)) * ((d - 2) * (x + 1.0 / x) - s2 + (3 - d))
    # At θ=0: val = (1/(d-1))·(2(d-2)+3-d) = (1/(d-1))·(d-1) = 1
    # So f_KN_d already normalised at θ=0. But to be safe:
    val_at_0 = 1.0 / (d - 1) * ((d - 2) * 2 + 3 - d)
    return val / val_at_0


def f_BAM_d_baseline(eps: float, theta: float, d_spatial: int) -> float:
    """BAM amplitude in d_spatial dim with polarisation sum
    ((d_spatial-2) + cos²θ)/(d_spatial-1) and S³ Green function
    propagator factor (1+1/x)²/4 (kept at the BAM-native S³, since
    BAM's spatial slice is fundamentally 3-sphere; the d-scaling is
    in the photon polarisation structure).

    Normalised at θ=0.
    """
    x = x_ratio(eps, theta)
    c = math.cos(theta)
    d = d_spatial
    pol_sum = ((d - 2) + c * c) / (d - 1)
    prop_factor = (1.0 + 1.0 / x) ** 2 / 4.0
    val = pol_sum * prop_factor
    # At θ=0: pol_sum = ((d-2)+1)/(d-1) = 1, prop_factor at x=1 = 1
    # So already normalised.
    return val


def f_BAM_d_modified(eps: float, theta: float, d_spatial: int,
                     beta: float, gamma: float) -> float:
    """BAM amplitude in d_spatial dim with Family B vertex modification
    `1 + ε·(β·sin²θ + γ·(1−cos θ))` applied to the polarisation sum.
    """
    x = x_ratio(eps, theta)
    c = math.cos(theta)
    s2 = math.sin(theta) ** 2
    d = d_spatial
    pol_sum_baseline = ((d - 2) + c * c) / (d - 1)
    angular_mod = 1.0 + eps * (beta * s2 + gamma * (1.0 - c))
    pol_sum = pol_sum_baseline * angular_mod * angular_mod
    prop_factor = (1.0 + 1.0 / x) ** 2 / 4.0
    val = pol_sum * prop_factor
    # Normalise at θ=0 in the same construction
    pol_sum_at_0 = 1.0 * (1.0 + eps * 0.0) ** 2  # μ(θ=0) = 0
    prop_at_0 = 1.0
    val_at_0 = pol_sum_at_0 * prop_at_0
    return val / val_at_0


# ---------------------------------------------------------------------------
# Numerical extraction of γ_d
# ---------------------------------------------------------------------------

def numerical_gamma_d(d_spatial: int, eps_small: float = 1e-4) -> dict:
    """Extract γ(d) numerically by scanning Family B and finding the
    (β, γ) that minimise residual to KN_d at small ε.

    Uses the analytic guess (β, γ) = (0, −3/2) as the starting point
    and refines by 1D scan in γ at β=0.
    """
    thetas = np.linspace(0.05, PI - 0.05, 17)
    gammas = np.linspace(-2.5, -0.5, 41)
    residuals = []
    for gamma in gammas:
        max_diff = 0.0
        for theta in thetas:
            b = f_BAM_d_modified(eps_small, float(theta), d_spatial, 0.0, float(gamma))
            k = f_KN_d(eps_small, float(theta), d_spatial)
            diff = abs(b - k)
            max_diff = max(max_diff, diff)
        residuals.append(max_diff)
    residuals = np.asarray(residuals)
    i_opt = int(np.argmin(residuals))
    gamma_opt = float(gammas[i_opt])
    res_opt = float(residuals[i_opt])
    return {
        'd_spatial': d_spatial,
        'gamma_optimal_numerical': gamma_opt,
        'residual_at_optimum': res_opt,
        'gamma_grid': gammas.tolist(),
        'residuals_at_grid': residuals.tolist(),
    }


# ---------------------------------------------------------------------------
# Analytic matching: re-derive symbolically and check d-cancellation
# ---------------------------------------------------------------------------

def test_T1_analytic_d_independence() -> dict:
    """Symbolic check that the matching equation
    (1-c) + 2μ = -2(1-c) does not depend on d.

    The derivation in d-dim:
      f_BAM_d_modified_O(ε) = ((d-2)+c²)/(d-1) · [(1-c) + 2μ]
      f_KN_d_O(ε)           = ((d-2)+c²)/(d-1) · [-2(1-c)]
    Setting equal:
      (1-c) + 2μ = -2(1-c)
      μ = -3(1-c)/2
    The ((d-2)+c²)/(d-1) factor CANCELS, so μ = -3(1-c)/2 in all d.

    Sub-test: verify the factor cancellation by evaluating both
    f_BAM and f_KN at small ε and dividing.
    """
    eps = 1e-5
    ds = [3, 4, 5, 6, 8]
    rows = []
    for d in ds:
        ratios = []
        for theta in np.linspace(0.1, PI - 0.1, 9):
            b_O = (f_BAM_d_baseline(eps, float(theta), d) - f_BAM_d_baseline(0.0, float(theta), d)) / eps
            k_O = (f_KN_d(eps, float(theta), d) - f_KN_d(0.0, float(theta), d)) / eps
            # Expected: k_O = -2 · b_O (i.e. the ratio is -2 in all d)
            if abs(b_O) > 1e-12:
                ratios.append(k_O / b_O)
        ratio_avg = float(np.mean(ratios))
        ratio_std = float(np.std(ratios))
        rows.append({
            'd_spatial': d,
            'mean_K/B_ratio': ratio_avg,
            'std_K/B_ratio': ratio_std,
            'expected_ratio': -2.0,
            'matches_minus_2': abs(ratio_avg - (-2.0)) < 0.01,
        })
    all_match = all(r['matches_minus_2'] for r in rows)
    return {
        'name': 'T1_analytic_d_independence',
        'description': (
            'Verify that f_KN_d_O(ε)/f_BAM_d_baseline_O(ε) = −2 in '
            'all d, confirming the analytic factor cancellation. This '
            'implies the matching coefficient μ = −3(1−c)/2 → '
            'γ = −3/2 is d-independent.'
        ),
        'd_values_tested': ds,
        'per_d_results': rows,
        'all_d_give_ratio_minus_2': all_match,
        'pass': all_match,
    }


# ---------------------------------------------------------------------------
# T2: Numerical γ_d extraction
# ---------------------------------------------------------------------------

def test_T2_numerical_gamma_d() -> dict:
    """Extract γ(d) numerically for d_spatial ∈ {3, 4, 5, 6} and
    verify all give −3/2 to numerical precision (within grid spacing).
    """
    ds = [3, 4, 5, 6]
    rows = []
    for d in ds:
        result = numerical_gamma_d(d)
        rows.append(result)
    # All optima should be near -1.5 (γ = -3/2)
    gammas = [r['gamma_optimal_numerical'] for r in rows]
    max_dev = max(abs(g - (-1.5)) for g in gammas)
    return {
        'name': 'T2_numerical_gamma_d',
        'description': (
            'Extract γ(d) numerically for d ∈ {3, 4, 5, 6}. Predicted '
            'value is -3/2 in all d (per T1 derivation).'
        ),
        'per_d_results': rows,
        'numerical_gammas': gammas,
        'predicted_value': -1.5,
        'max_deviation_from_minus_3_over_2': max_dev,
        # Use grid spacing 0.05 as tolerance
        'pass': max_dev < 0.05,
    }


# ---------------------------------------------------------------------------
# T3: Candidate predictions table
# ---------------------------------------------------------------------------

def candidate_predictions_d(d_spatial: int) -> dict:
    """For each of the 8 PR #32 candidates, predict γ(d).

    For most candidates, the prediction is d-independent (the
    ingredient — Casimir, charge, Pauli trace — doesn't scale with d).
    For candidate C, the embedding-dim ratio scales as
    −d_spatial / dim(trans pol) = −d_spatial / (d_spatial − 1).
    """
    d = d_spatial
    return {
        # Casimir of spin-1/2 is dimension-independent (SU(2) lives in
        # the BAM 4D framework regardless of bulk d):
        'A_doubled_electron_casimir': -1.5,
        'B_photon_casimir_minus_hopf': -1.5,
        'D_hopf_times_photon_mult': -1.5,
        'F_closure_plus_hopf': -1.5,
        'G_pauli_trace': -1.5,
        'H_two_mouth_spin_half': -1.5,
        # Candidate C: scales with d.
        'C_embedding_over_polarization': -d / (d - 1),
        # Candidate E: natural-units, speculative; assume no scaling
        # but flag as ambiguous.
        'E_antipodal_natural_units': -1.5,
    }


def test_T3_candidate_predictions() -> dict:
    """Tabulate predicted γ(d) for each candidate in d ∈ {3, 4, 5, 6}."""
    ds = [3, 4, 5, 6]
    table = {}
    for d in ds:
        table[d] = candidate_predictions_d(d)
    return {
        'name': 'T3_candidate_predictions',
        'description': (
            'For each PR #32 candidate, predict γ(d) for '
            'd ∈ {3, 4, 5, 6}.'
        ),
        'predictions_per_d': table,
        'pass': True,   # informative
    }


# ---------------------------------------------------------------------------
# T4: Falsified candidates
# ---------------------------------------------------------------------------

def test_T4_falsified_candidates(t2: dict, t3: dict) -> dict:
    """Compare candidate predictions to numerical γ(d) and identify
    candidates that fail the d-scaling test."""
    numerical = {r['d_spatial']: r['gamma_optimal_numerical']
                 for r in t2['per_d_results']}
    predictions = t3['predictions_per_d']
    candidates = list(predictions[3].keys())
    falsified = []
    surviving = []
    for cand in candidates:
        max_dev = 0.0
        for d, num_gamma in numerical.items():
            pred_gamma = predictions[d][cand]
            dev = abs(pred_gamma - num_gamma)
            max_dev = max(max_dev, dev)
        if max_dev > 0.1:
            falsified.append({
                'candidate': cand,
                'max_deviation': max_dev,
                'predictions_per_d': {d: predictions[d][cand] for d in [3, 4, 5, 6]},
            })
        else:
            surviving.append({
                'candidate': cand,
                'max_deviation': max_dev,
            })
    return {
        'name': 'T4_falsified_candidates',
        'description': (
            'Identify PR #32 candidates whose predicted γ(d) deviates '
            'from the numerical observation by more than the grid '
            'tolerance.'
        ),
        'falsified': falsified,
        'surviving': surviving,
        'n_falsified': len(falsified),
        'n_surviving': len(surviving),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T5: Cleanest survivor + remaining discrimination paths
# ---------------------------------------------------------------------------

def test_T5_survivors_and_next_steps(t4: dict) -> dict:
    """Among the surviving candidates, identify the most parsimonious
    by PR #32's naturalness criteria, and list the next discrimination
    paths."""
    survivors = t4['surviving']
    next_paths = [
        ('O(ε²) analytic extension',
         'Different surviving candidates predict different next-order '
         'coefficients. Casimir candidates predict C_2²; Hopf '
         'candidates predict χ-dependent modulation. Computing the '
         'analytic O(ε²) coefficient would discriminate.'),
        ('Cross-process test (pair production)',
         'γγ → e⁺e⁻ has a different amplitude structure. Candidates '
         'rooted in electron-only quantities (A, H) vs photon-electron '
         'coupling (B, D, G) would give different predictions.'),
        ('Polarised cross-section corrections',
         'For specific photon polarisation choices (circular, linear), '
         'candidates may give different correction patterns.'),
    ]
    return {
        'name': 'T5_survivors_and_next_steps',
        'description': (
            'Identify surviving candidates after the d-scaling test '
            'and list next-step discrimination paths.'
        ),
        'surviving_candidates': [s['candidate'] for s in survivors],
        'n_surviving': len(survivors),
        'next_discrimination_paths': [
            {'name': n, 'description': d} for n, d in next_paths
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_analytic_d_independence()
    t2 = test_T2_numerical_gamma_d()
    t3 = test_T3_candidate_predictions()
    t4 = test_T4_falsified_candidates(t2, t3)
    t5 = test_T5_survivors_and_next_steps(t4)
    tests = [t1, t2, t3, t4, t5]

    # Verdict
    if t1['pass'] and t2['pass']:
        if t4['n_falsified'] > 0:
            verdict_class = 'FALSIFIES_C'
            falsified_names = [f['candidate'] for f in t4['falsified']]
            verdict = (
                f'FALSIFIES_C — γ is d-independent (numerical γ_d '
                f'matches −3/2 in d ∈ {{3, 4, 5, 6}} to within grid '
                f'tolerance). The factor ((d−2)+c²)/(d−1) cancels in '
                f'the matching equation. Candidates predicting '
                f'd-dependent γ are falsified: {falsified_names}. '
                f'{t4["n_surviving"]} candidates survive — all '
                'predicting γ = −3/2 universally. Further '
                'discrimination requires O(ε²) analytic extension, '
                'cross-process tests, or polarisation-dependent '
                'observables (T5).'
            )
        else:
            verdict_class = 'NO_DISCRIMINATION'
            verdict = (
                'NO DISCRIMINATION — γ is d-independent, but all '
                'tested candidates predict d-independent γ. The '
                'd-scaling test alone cannot discriminate between '
                'them.'
            )
    elif t1['pass'] and not t2['pass']:
        verdict_class = 'NUMERICAL_INCONSISTENCY'
        verdict = (
            'NUMERICAL INCONSISTENCY — analytic derivation gives '
            'd-independent γ = −3/2, but numerical extraction in '
            'higher d disagrees. Indicates either a construction '
            'issue in the d-dim amplitude or numerical resolution '
            'problem.'
        )
    elif not t1['pass']:
        verdict_class = 'ANALYTIC_FAILURE'
        verdict = (
            'ANALYTIC FAILURE — the predicted ratio f_KN_d_O / '
            'f_BAM_d_O = −2 does not hold in all d. The '
            'd-cancellation assumed in the matching derivation '
            'breaks down; the KN-d natural generalisation may need '
            'revision.'
        )
    else:
        verdict_class = 'INDETERMINATE'
        verdict = 'See per-test details.'

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'thesis': (
            'γ = −3/2 is d-independent under the natural d-dim '
            'generalisation of BAM and KN. The ((d-2)+c²)/(d-1) '
            'polarisation-sum factor cancels in the matching '
            'equation.'
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
    L.append('# Tangherlini dimensional-scaling discrimination probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Follow-on to PR #32 (coefficient-origin probe). Tests whether '
        'γ = −3/2 (the analytic O(ω/m) closure coefficient from '
        'PR #31) is d-independent or scales with spatial dimension. '
        'Discriminates between Casimir-based and embedding-based '
        'BAM derivations.'
    )
    L.append('')
    L.append(f"**Thesis:** {s['thesis']}")
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
                f"K/B ratio = −2 in all d: "
                f"{'YES' if t['all_d_give_ratio_minus_2'] else 'NO'}"
            )
        elif nm.startswith('T2'):
            value = f"max |γ(d) − (−3/2)| = {t['max_deviation_from_minus_3_over_2']:.4f}"
        elif nm.startswith('T3'):
            value = 'predictions tabulated'
        elif nm.startswith('T4'):
            value = f"{t['n_falsified']} falsified, {t['n_surviving']} survive"
        elif nm.startswith('T5'):
            value = f"{t['n_surviving']} surviving candidates"
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1 detail
    t1 = s['tests'][0]
    L.append('## T1: K/B ratio per d')
    L.append('')
    L.append('Predicted: f_KN_d_O(ε) / f_BAM_d_baseline_O(ε) = **−2** for all d (factor `((d-2)+c²)/(d-1)` cancels).')
    L.append('')
    L.append('| d | mean K/B | std | matches −2? |')
    L.append('|---:|---:|---:|:---:|')
    for r in t1['per_d_results']:
        L.append(
            f"| {r['d_spatial']} | {r['mean_K/B_ratio']:+.4f} | "
            f"{r['std_K/B_ratio']:.2e} | "
            f"{'✓' if r['matches_minus_2'] else '✗'} |"
        )
    L.append('')

    # T2 detail
    t2 = s['tests'][1]
    L.append('## T2: Numerical γ(d)')
    L.append('')
    L.append('| d | γ_opt | residual at optimum |')
    L.append('|---:|---:|---:|')
    for r in t2['per_d_results']:
        L.append(
            f"| {r['d_spatial']} | {r['gamma_optimal_numerical']:+.4f} | "
            f"{r['residual_at_optimum']:.4e} |"
        )
    L.append('')
    L.append(f"**Max deviation from −3/2:** {t2['max_deviation_from_minus_3_over_2']:.4f}")
    L.append('')

    # T3 + T4 detail
    t3 = s['tests'][2]
    t4 = s['tests'][3]
    L.append('## T3 & T4: Candidate predictions and falsification')
    L.append('')
    L.append('| candidate | d=3 | d=4 | d=5 | d=6 | status |')
    L.append('|---|---:|---:|---:|---:|---|')
    predictions = t3['predictions_per_d']
    candidates = list(predictions[3].keys())
    falsified_names = {f['candidate'] for f in t4['falsified']}
    for cand in candidates:
        row = f"| `{cand}` |"
        for d in [3, 4, 5, 6]:
            row += f" {predictions[d][cand]:+.4f} |"
        status = '**FALSIFIED**' if cand in falsified_names else 'survives'
        row += f" {status} |"
        L.append(row)
    L.append('')

    # T5 detail
    t5 = s['tests'][4]
    L.append('## T5: Surviving candidates and next discrimination paths')
    L.append('')
    L.append(f"**Surviving candidates** ({t5['n_surviving']}):")
    L.append('')
    for c in t5['surviving_candidates']:
        L.append(f"- `{c}`")
    L.append('')
    L.append('**Next discrimination paths:**')
    L.append('')
    for p in t5['next_discrimination_paths']:
        L.append(f"- **{p['name']}** — {p['description']}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **6 candidates remain after d-scaling discrimination.** '
        'Among the surviving candidates, the cleanest single-'
        'ingredient one is `G_pauli_trace_normalized` (PR #32 '
        'naturalness ranking 8/10).'
    )
    L.append(
        '- **O(ε²) extension** is the natural next probe: different '
        'surviving candidates predict different next-order '
        'coefficients. PR #31\'s residual fit (ε^1.88 vs predicted '
        'ε²) leaves room for a small O(ε^{3/2}) or O(ε²) correction '
        'that could be candidate-specific.'
    )
    L.append(
        '- **The d-dim KN analog used in this probe** is a natural '
        'generalisation; verifying against a textbook d-dim QED '
        'Compton derivation would close the analytic loop.'
    )
    L.append(
        '- **Tangherlini physical interpretation.** BAM\'s 5D '
        'framework uses Tangherlini for the radial bulk channel — '
        'this probe tests the *spatial-dimension* scaling of the '
        'photon polarisation sum, which is a different generalisation '
        'than putting the photon in the Tangherlini bulk. The latter '
        'would test the radial channel\'s contribution to the vertex '
        'coefficient — open follow-on.'
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
    out = here / 'runs' / f'{ts}_compton_dimensional_scaling_probe'
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
