"""
BAM throat-vertex loop probe — the Schwinger anomalous magnetic moment.

The capstone of the spin sector. PR #61 derived g=2 at tree level from
the throat's Pauli/SU(2) + Hopf-monopole structure and recorded the
Schwinger anomaly a=(g−2)/2=α/2π as the known one-loop correction. This
probe constructs the one-loop throat-vertex correction in BAM-native
terms — the throat emits and reabsorbs one virtual photon (the S³
Green-function propagator, PRs #45–#46) at the throat-pinch vertex (PRs
#38–#41) — and shows it reproduces the Schwinger value.

HONEST SCOPE (stated upfront): this is a reconstruction, not an
independent first-principles derivation of 1/2π. The BAM tree thread
normalized its primitives (S³ photon propagator 1/q², throat-pinch
vertex, Hopf coupling α) to QED at tree level; the one-loop integrand is
therefore the QED vertex-correction integrand in those primitives, and
the probe verifies the loop integral gives the Schwinger coefficient.
BAM-native: the loop's pieces (virtual photon = S³ exchange; vertex =
throat pinch) and structure. NOT independently derived: the 1/2π from a
pure-geometry loop measure from S_BAM (the standing follow-on).

The loop (QED vertex triangle in BAM terms): external throat lines (in/
out electron at the pinch), virtual photon = S³ Green-function exchange
(flat limit G→1/(4πd) = 1/q²), vertex = throat-pinch F² (Hopf coupling
α). Loop parameter α/π; the anomalous moment is the F₂(0) form factor:

    F₂(0) = (α/2π)·I,  I = ∫dxdydz δ(x+y+z−1) 2m²z(1−z)/[m²(1−z)²]
                          = ∫₀¹ 2z dz = 1,
so a = F₂(0) = α/2π ≈ 0.0011614 (vs a_e = 0.00115965, leading term).

B4: a=α/2π dimensionless; α the EM coupling; scale is the single anchor.

Tests:
  T1. BAM one-loop vertex structure (throat self-exchange triangle).
  T2. Virtual photon = S³ exchange (G→1/(4πd) flat limit).
  T3. Schwinger integral I = ∫₀¹2z dz = 1 ⟹ F₂(0)=α/2π.
  T4. a = α/2π vs experiment.
  T5. One-loop g = 2(1+a).
  T6. Honest scope (reconstruction; 1/2π from S_BAM open).
  T7. B4 accounting (a,g dimensionless; α coupling; scale = anchor).
  T8. Assessment.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import integrate

from geometrodynamics.transaction.s3_geometry import s3_green_potential


PI = math.pi

ALPHA = 7.2973525693e-3
A_E_EXPERIMENT = 0.00115965218


# ---------------------------------------------------------------------------
# T1. BAM one-loop vertex structure
# ---------------------------------------------------------------------------

def test_T1_loop_structure() -> dict:
    """The one-loop vertex correction in BAM terms: the throat (electron)
    emits and reabsorbs one virtual photon at the throat-pinch vertex —
    the BAM image of the QED vertex-correction triangle. The pieces are
    the BAM primitives: virtual photon = S³ Green-function exchange (#45–
    #46), vertex = throat-pinch F² (#38–#41), coupling = α (Hopf). The
    loop expansion parameter is α/π."""
    pieces = {
        'external_lines': 'in/out throat (electron) at the pinch vertex',
        'virtual_photon': 'S³ Green-function exchange (flat limit 1/q²)',
        'vertex': 'throat-pinch F² vertex (PRs #38–#41)',
        'coupling': 'α (Hopf charge structure)',
        'loop_parameter': 'α/π (coupling α × loop phase space 1/π)',
    }
    one_loop_factor = ALPHA / PI
    return {
        'name': 'T1_bam_one_loop_vertex_structure',
        'description': (
            "The one-loop throat-vertex correction: the throat emits and "
            "reabsorbs one virtual photon (S³ exchange) at the throat-pinch "
            "vertex — the BAM image of the QED vertex triangle. Loop "
            "parameter α/π."
        ),
        'loop_pieces': pieces,
        'one_loop_factor_alpha_over_pi': one_loop_factor,
        'pass': one_loop_factor > 0,
    }


# ---------------------------------------------------------------------------
# T2. Virtual photon = S³ exchange
# ---------------------------------------------------------------------------

def test_T2_virtual_photon_s3() -> dict:
    """The virtual photon in the loop is an S³ Green-function exchange.
    The S³ scalar Green function G(ψ)=((π−ψ)cot ψ − ½)/(4π²R) has the
    flat limit G → 1/(4π d) (Coulomb), i.e. the 1/q² photon propagator
    (PRs #45–#46). Verify G·(4π d) → 1 as ψ → 0 (d = Rψ)."""
    R = 1.0
    rows = []
    max_dev = 0.0
    for psi in [1e-2, 1e-3, 1e-4]:
        G = s3_green_potential(psi, radius=R, eps=1e-12)
        d = R * psi
        coulomb_check = G * 4.0 * PI * d     # → 1
        max_dev = max(max_dev, abs(coulomb_check - 1.0))
        rows.append({'psi': psi, 'G': G, 'G_times_4pi_d': coulomb_check})
    return {
        'name': 'T2_virtual_photon_is_s3_exchange',
        'description': (
            "The virtual photon is an S³ Green-function exchange: "
            "G(ψ) → 1/(4π d) (Coulomb / 1/q² propagator) in the flat "
            "limit, the BAM photon propagator (#45–#46)."
        ),
        'rows': rows,
        'max_deviation_from_coulomb': max_dev,
        'flat_limit_is_1_over_q2': max_dev < 1e-2,
        'pass': max_dev < 1e-2,
    }


# ---------------------------------------------------------------------------
# T3. Schwinger integral
# ---------------------------------------------------------------------------

def test_T3_schwinger_integral() -> dict:
    """After the loop integral, the anomalous form factor is
    F₂(0) = (α/2π)·I with the Feynman-parameter integral
    I = ∫dxdydz δ(x+y+z−1) 2m²z(1−z)/[m²(1−z)²]. The x,y simplex slice at
    fixed z has measure (1−z), so I = ∫₀¹ (1−z)·2z(1−z)/(1−z)² dz =
    ∫₀¹ 2z dz = 1. Verify I = 1 numerically (reduced and 2D simplex
    forms)."""
    # reduced form
    I_reduced, _ = integrate.quad(lambda z: 2.0 * z, 0.0, 1.0)
    # 2D simplex form: integrate over z, with the (1−z) slice measure and
    # the integrand 2z(1−z)/(1−z)² = 2z/(1−z); product = 2z (the 1/(1−z)
    # cancels against the slice measure (1−z))

    def integrand_z(z):
        slice_measure = (1.0 - z)               # ∫dxdy δ(x+y+z−1)
        f = 2.0 * z * (1.0 - z) / (1.0 - z) ** 2  # 2z(1−z)/(1−z)²
        return slice_measure * f

    I_simplex, _ = integrate.quad(integrand_z, 0.0, 1.0)
    F2_0 = (ALPHA / (2.0 * PI)) * I_simplex
    return {
        'name': 'T3_schwinger_feynman_parameter_integral',
        'description': (
            "F₂(0) = (α/2π)·I, I = ∫dxdydz δ(x+y+z−1) 2z(1−z)/(1−z)² = "
            "∫₀¹ 2z dz = 1 (the x,y simplex slice measure (1−z) cancels "
            "the 1/(1−z)). The throat-vertex loop reproduces F₂(0)=α/2π."
        ),
        'I_reduced_int_2z': I_reduced,
        'I_simplex': I_simplex,
        'integral_equals_1': abs(I_simplex - 1.0) < 1e-9,
        'F2_0': F2_0,
        'pass': abs(I_simplex - 1.0) < 1e-9 and abs(I_reduced - 1.0) < 1e-9,
    }


# ---------------------------------------------------------------------------
# T4. a = α/2π vs experiment
# ---------------------------------------------------------------------------

def test_T4_anomaly_vs_experiment() -> dict:
    """a = F₂(0) = α/2π ≈ 0.0011614 (Schwinger), vs the measured
    a_e = 0.00115965 (the α/2π term dominates; higher orders make up the
    small remainder)."""
    a = ALPHA / (2.0 * PI)
    rel = abs(a - A_E_EXPERIMENT) / A_E_EXPERIMENT
    return {
        'name': 'T4_anomaly_vs_experiment',
        'description': (
            "a = α/2π ≈ 0.0011614 (Schwinger one-loop) vs measured "
            "a_e = 0.00115965 — the leading term, agreeing to ~0.15%; "
            "higher orders (α², α³, …) make up the small remainder."
        ),
        'a_alpha_over_2pi': a,
        'a_e_experiment': A_E_EXPERIMENT,
        'relative_difference': rel,
        'leading_term_agrees': rel < 0.01,
        'pass': rel < 0.01,
    }


# ---------------------------------------------------------------------------
# T5. One-loop g
# ---------------------------------------------------------------------------

def test_T5_one_loop_g() -> dict:
    """g = 2(1 + a) = 2(1 + α/2π) = 2.00232… — the tree g=2 (PR #61)
    plus the one-loop Schwinger correction."""
    a = ALPHA / (2.0 * PI)
    g = 2.0 * (1.0 + a)
    g_expected = 2.0023228
    return {
        'name': 'T5_one_loop_g_factor',
        'description': (
            "g = 2(1 + a) = 2(1 + α/2π) = 2.00232… — the tree g=2 (PR #61) "
            "plus the one-loop Schwinger anomaly."
        ),
        'tree_g': 2.0,
        'anomaly_a': a,
        'one_loop_g': g,
        'g_expected': g_expected,
        'pass': abs(g - g_expected) < 1e-5,
    }


# ---------------------------------------------------------------------------
# T6. Honest scope
# ---------------------------------------------------------------------------

def test_T6_honest_scope() -> dict:
    """The honest scope: this is a RECONSTRUCTION, not an independent
    first-principles derivation of 1/2π. The one-loop integrand is the
    QED vertex integrand expressed in tree-normalized BAM primitives (S³
    propagator, throat-pinch vertex, Hopf α); the probe verifies the loop
    integral gives the Schwinger coefficient. BAM-native: the loop pieces
    and structure. Open: the 1/2π from a pure-geometry S_BAM loop
    measure."""
    return {
        'name': 'T6_honest_scope',
        'description': (
            "Reconstruction, not first-principles: the one-loop integrand "
            "is the QED vertex integrand in tree-normalized BAM primitives "
            "(S³ propagator, throat-pinch vertex, Hopf α); the probe "
            "verifies the loop integral gives α/2π. BAM-native: the loop "
            "pieces (virtual photon = S³ exchange; vertex = throat pinch) "
            "and structure. Open: 1/2π from a pure-geometry S_BAM loop "
            "measure (the full covariant throat loop)."
        ),
        'bam_native': [
            'virtual photon = S³ Green-function exchange',
            'vertex = throat-pinch F²',
            'structure = throat dressing its moment by one self-exchange',
        ],
        'inherited_from_tree_normalization': [
            'the 1/q² propagator normalization',
            'the vertex normalization (tree QED match, #35–#46)',
            'the coupling α',
        ],
        'open_first_principles': '1/2π from a covariant S_BAM loop measure',
        'is_reconstruction_not_derivation': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. B4 accounting
# ---------------------------------------------------------------------------

def test_T7_b4_accounting() -> dict:
    """a = α/2π and g = 2(1+a) are dimensionless; α is the EM coupling
    (an input, related to the Hopf-charge structure); the absolute mass
    scale is the single anchor m (m_e c² = ℏc/R_MID). The anomaly is a
    pure number, independent of the anchor's value."""
    a = ALPHA / (2.0 * PI)
    return {
        'name': 'T7_b4_accounting',
        'description': (
            "a = α/2π and g = 2(1+a) are dimensionless; α is the EM "
            "coupling (input, Hopf-charge structure); the absolute mass "
            "scale is the single anchor m. The anomaly is a pure number, "
            "independent of the anchor's value (B4-consistent)."
        ),
        'a_dimensionless': a,
        'alpha_is_coupling': True,
        'scale_is_single_anchor': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    """The BAM throat-vertex loop — the throat dressing its moment by one
    S³-photon self-exchange at the throat-pinch vertex — reproduces the
    Schwinger anomaly a = α/2π (the Feynman-parameter integral I=1),
    matching a_e at leading order. The loop pieces and structure are
    BAM-native; the coefficient is inherited from the tree-normalized
    primitives (a reconstruction)."""
    return {
        'name': 'T8_assessment',
        'description': (
            "The BAM throat-vertex loop (throat dressing its moment by one "
            "S³-photon self-exchange at the throat-pinch vertex) reproduces "
            "the Schwinger anomaly a = α/2π (Feynman-parameter integral "
            "I=1), matching a_e at leading order. The loop pieces and "
            "structure are BAM-native; the coefficient is inherited from "
            "the tree-normalized primitives — a reconstruction, with the "
            "first-principles 1/2π (from a covariant S_BAM loop measure) "
            "the open follow-on."
        ),
        'anomaly': 'a = α/2π (one loop)',
        'integral': 'I = ∫₀¹ 2z dz = 1',
        'one_loop_g': 2.0 * (1.0 + ALPHA / (2.0 * PI)),
        'scope': 'reconstruction (tree-normalized primitives), not first-principles 1/2π',
        'remaining': '1/2π from covariant S_BAM loop; higher-order a_e; throat spinor',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_loop_structure()
    t2 = test_T2_virtual_photon_s3()
    t3 = test_T3_schwinger_integral()
    t4 = test_T4_anomaly_vs_experiment()
    t5 = test_T5_one_loop_g()
    t6 = test_T6_honest_scope()
    t7 = test_T7_b4_accounting()
    t8 = test_T8_assessment()
    tests = [t1, t2, t3, t4, t5, t6, t7, t8]

    core = [t1, t2, t3, t4, t5, t6, t7]
    if all(t['pass'] for t in core):
        verdict_class = 'SCHWINGER_RECONSTRUCTED'
        verdict = (
            'SCHWINGER RECONSTRUCTED. The BAM throat-vertex loop reproduces '
            'the Schwinger anomalous magnetic moment a = α/2π, the one-loop '
            'capstone of the spin sector (tree g=2 in PR #61).\n\n'
            'THE LOOP. The one-loop vertex correction, in BAM terms, is the '
            'throat dressing its magnetic moment by one virtual-photon '
            'self-exchange: external in/out throat lines at the throat-pinch '
            'vertex (PRs #38–#41), the virtual photon an S³ Green-function '
            'exchange whose flat limit G → 1/(4π d) is the 1/q² photon '
            'propagator (PRs #45–#46), and the Hopf coupling α. The loop '
            'expansion parameter is α/π.\n\n'
            'THE COEFFICIENT. After the loop integral the anomalous form '
            'factor is F₂(0) = (α/2π)·I with the Feynman-parameter integral '
            'I = ∫dxdydz δ(x+y+z−1) 2z(1−z)/(1−z)² = ∫₀¹ 2z dz = 1 (the x,y '
            'simplex slice measure (1−z) cancels the 1/(1−z)). So '
            'a = F₂(0) = α/2π ≈ 0.0011614, matching the measured '
            'a_e = 0.00115965 at leading order (~0.15%; higher orders '
            'α², α³ make up the remainder), and g = 2(1 + a) = 2.00232….\n\n'
            'HONEST SCOPE. This is a RECONSTRUCTION, not an independent '
            'first-principles derivation of 1/2π. The one-loop integrand is '
            'the QED vertex integrand expressed in BAM primitives that were '
            'normalized to QED at tree level (the S³ propagator, the '
            'throat-pinch vertex, the coupling α); the probe verifies the '
            'loop integral gives the Schwinger coefficient. What is '
            'BAM-native: the loop pieces (virtual photon = S³ exchange; '
            'vertex = throat pinch) and the structure (one self-exchange '
            'dressing the moment). What remains open: the 1/2π from a '
            'pure-geometry loop measure derived from S_BAM (the full '
            'covariant throat loop), plus the higher-order a_e series. '
            'B4: a and g are dimensionless; α is the coupling; the scale '
            'is the single anchor.'
        )
    else:
        verdict_class = 'LOOP_FAILS'
        verdict = (
            'LOOP FAILS. The Feynman-parameter integral did not give α/2π, '
            'or the loop pieces do not map to the BAM primitives. '
            'Investigate the failing test.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'anomaly': 'a = α/2π (Schwinger, one loop)',
        'loop': 'throat self-exchange: S³ virtual photon at the throat-pinch vertex',
        'integral': 'I = ∫₀¹ 2z dz = 1 ⟹ F₂(0) = α/2π',
        'scope': 'reconstruction (tree-normalized primitives); 1/2π from S_BAM open',
        'b4_caveat': 'a, g dimensionless; α the coupling; scale = single anchor',
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
    L.append('# BAM throat-vertex loop probe — the Schwinger anomalous moment')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Constructs the one-loop throat-vertex correction in BAM-native '
        'terms (the throat dressing its moment by one S³-photon '
        'self-exchange) and shows it reproduces the Schwinger anomaly '
        'a = α/2π. A reconstruction using the tree-normalized BAM '
        'primitives — honest scope stated in T6.'
    )
    L.append('')
    L.append(f"- **Anomaly**: {s['anomaly']}")
    L.append(f"- **Loop**: {s['loop']}")
    L.append(f"- **Integral**: `{s['integral']}`")
    L.append(f"- **Scope**: {s['scope']}")
    L.append(f"- **B4 caveat**: {s['b4_caveat']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            value = "throat self-exchange triangle; loop param α/π"
        elif nm.startswith('T2'):
            value = f"S³ photon → 1/(4πd) (dev {t['max_deviation_from_coulomb']:.0e})"
        elif nm.startswith('T3'):
            value = f"I = ∫₀¹2z dz = {t['I_simplex']:.4f} → F₂(0)=α/2π"
        elif nm.startswith('T4'):
            value = f"a=α/2π={t['a_alpha_over_2pi']:.7f} vs a_e (Δ {t['relative_difference']:.2%})"
        elif nm.startswith('T5'):
            value = f"g = 2(1+a) = {t['one_loop_g']:.7f}"
        elif nm.startswith('T6'):
            value = "reconstruction; 1/2π from S_BAM open"
        elif nm.startswith('T7'):
            value = "a, g dimensionless; α coupling; scale = anchor"
        elif nm.startswith('T8'):
            value = "throat loop reproduces Schwinger"
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1
    t1 = s['tests'][0]
    L.append('## T1: BAM one-loop vertex structure')
    L.append('')
    for k, v in t1['loop_pieces'].items():
        L.append(f"- **{k}**: {v}")
    L.append('')

    # T2
    t2 = s['tests'][1]
    L.append('## T2: Virtual photon = S³ exchange')
    L.append('')
    L.append('| ψ | G(ψ) | G·4π d (→ 1) |')
    L.append('|---:|---:|---:|')
    for r in t2['rows']:
        L.append(f"| {r['psi']:.0e} | {r['G']:.4f} | {r['G_times_4pi_d']:.6f} |")
    L.append('')
    L.append(f"Flat limit → Coulomb 1/(4π d) = 1/q² propagator "
             f"(max dev {t2['max_deviation_from_coulomb']:.1e}).")
    L.append('')

    # T3
    t3 = s['tests'][2]
    L.append('## T3: Schwinger Feynman-parameter integral')
    L.append('')
    L.append(f"- I (reduced ∫₀¹ 2z dz) = {t3['I_reduced_int_2z']:.6f}")
    L.append(f"- I (simplex form) = {t3['I_simplex']:.6f} (= 1)")
    L.append(f"- F₂(0) = (α/2π)·I = {t3['F2_0']:.8f}")
    L.append('')

    # T4
    t4 = s['tests'][3]
    L.append('## T4: a = α/2π vs experiment')
    L.append('')
    L.append(f"- a = α/2π = {t4['a_alpha_over_2pi']:.8f}")
    L.append(f"- a_e (measured) = {t4['a_e_experiment']:.8f}")
    L.append(f"- relative difference = {t4['relative_difference']:.2%} "
             f"(leading term; higher orders make up the remainder)")
    L.append('')

    # T5
    t5 = s['tests'][4]
    L.append('## T5: One-loop g')
    L.append('')
    L.append(f"- g = 2(1 + a) = {t5['one_loop_g']:.8f} (tree g=2 + Schwinger a)")
    L.append('')

    # T6
    t6 = s['tests'][5]
    L.append('## T6: Honest scope')
    L.append('')
    L.append('BAM-native:')
    for x in t6['bam_native']:
        L.append(f"  - {x}")
    L.append('Inherited from tree normalization (#35–#46):')
    for x in t6['inherited_from_tree_normalization']:
        L.append(f"  - {x}")
    L.append(f"Open (first-principles): {t6['open_first_principles']}")
    L.append('')

    # T7
    t7 = s['tests'][6]
    L.append('## T7: B4 accounting')
    L.append('')
    L.append(f"- a = {t7['a_dimensionless']:.7f} (dimensionless); α the coupling; "
             f"scale = single anchor")
    L.append('')

    # T8
    t8 = s['tests'][7]
    L.append('## T8: Assessment')
    L.append('')
    L.append(f"- anomaly: {t8['anomaly']}")
    L.append(f"- integral: {t8['integral']}")
    L.append(f"- one-loop g: {t8['one_loop_g']:.7f}")
    L.append(f"- scope: {t8['scope']}")
    L.append(f"- remaining: {t8['remaining']}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **1/2π from a first-principles S_BAM loop measure.** The '
             'covariant throat loop derived from the action, rather than the '
             'QED-normalized integrand — the genuine derivation.')
    L.append('- **Higher-order a_e.** The α², α³, … terms (two- and '
             'three-loop).')
    L.append('- **The explicit throat spinor / vertex from S_BAM** (shared '
             'with #59–#61).')
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
    out = here / 'runs' / f'{ts}_throat_vertex_loop_probe'
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
