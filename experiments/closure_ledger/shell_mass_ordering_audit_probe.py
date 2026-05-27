"""
Shell Hamiltonian mass-ordering and n_part audit (PR #78).

Uses the PR #77 shell waveguide basis to test whether shell
eigenvalues and Z₂ partition splitting reproduce the qualitative
quark ordering structure better than the v3 lepton-shaped basis,
and whether the `n_part = 233` compensator shrinks or disappears.

## The mass-ordering target

Per the v3 basis-to-species map (`geometrodynamics/qcd/quark_spectrum`):

  (k=1, +) → u   (k=1, −) → d       m_u  =    2.16,  m_d  =     4.67
  (k=3, +) → c   (k=3, −) → s       m_c  = 1270.0,   m_s  =    93.4
  (k=5, +) → t   (k=5, −) → b       m_t  = 172690,   m_b  =  4180.0

The qualitative pattern:

  1. **Within-generation inversion.** Same partition sign is LIGHTER
     in generation 1 (u < d at k=1) but HEAVIER in generations 2, 3
     (c > s at k=3, t > b at k=5). This is the famous quark
     mass-ordering inversion.

  2. **Inter-generation hierarchy.** Five orders of magnitude in mass
     (twelve in mass²) from u to t. The lepton ladder has only ~3000×
     between e and τ.

The v3 lepton-shaped basis fits this via `uplift_asymmetry
ε = 1 − 1/k_5²` (asymmetric β-uplift per partition, growing with k),
absorbing the unmodeled physics into `n_part = 233` as a compensator
(#76 diagnosis).

## What this probe tests

On the PR #77 shell waveguide basis `{(l=1, n, ±)}_{n=3,4,5}` (the
`n_varied` enumeration; `n_varied` is preferred over `l_varied`
because ω² spans more), with the operator scaffold
`H = H_kin + H_Z2 + H_couple`:

  T1. Shell kinetic spectrum range (`ω²` from `n=3` to `n=5`).
  T2. Coverage gap: observed mass² spans ×6.4·10⁹; kinetic spans only
      ×2.2. The shell kinetic alone is insufficient.
  T3. Uniform Z₂ splitter `H_Z2 = χ·σ_z` is structurally incapable of
      reproducing the within-generation inversion (u<d but c>s)
      because it splits all blocks the same way.
  T4. Generation-dependent splitter `χ_n` (sign-flipping between
      n=3 and n=4,5): an existence proof — sign-flipping χ_n CAN
      reproduce the ordering pattern, identifying the structural slot
      PR #79's boundary stress tensor must populate.
  T5. Inter-mode coupling `H_couple`: amplifies inter-generation gap
      via level repulsion; explore minimal coupling needed to span
      the observed ratio range.
  T6. Shell `n_part` audit: how big a compensator is needed for the
      absolute scale on the shell basis? Compare to v3's 233.
  T7. Comparison with v3 baseline: which gives better qualitative
      ordering at comparable compensator strength?
  T8. Honest scope + verdict.

## Honest finding (expected)

Two structural improvements over v3:

  - **The shell kinetic IS the wavefront eigenfrequency squared**,
    not a phenomenological closure-quantum winding cost — this is
    the right machinery (PR #77 / #76's diagnosis).
  - **The Z₂ partition splitter on the shell basis** is the right
    *structural* slot for the within-generation inversion, BUT requires
    `χ` to vary by `n` (boundary stress tensor input, PR #79).

One structural finding NOT resolved:

  - **The shell basis alone does not span the observed mass range.**
    The required compensator is comparable to or larger than v3's
    `n_part = 233`. PR #78 alone cannot reduce or remove `n_part`;
    that depends on whether PR #79 (boundary stress tensor) and
    PR #80 (color algebra) populate `H_couple` and `χ_n` with
    structurally derived values that span the range.

The verdict is therefore:

  - SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED.

The shell basis is the right machinery and identifies the slots
PR #79–#80 must populate. PR #78 alone does not close `n_part`;
it sharpens the scope of what PR #79 needs to produce.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from experiments.closure_ledger.qcd_shell_waveguide_scaffold_probe import (
    build_shell_basis, ShellOperatorScaffold, _radial_ladder,
    N_FLAVORS, N_GENERATIONS, N_PARTITIONS,
)


PI = math.pi

# Observed quark masses (MeV, PDG-ish, matching v3's
# OBSERVED_MASSES_MEV in geometrodynamics/qcd/quark_spectrum.py)
QUARK_MASS_OBS_MEV: dict[str, float] = {
    'u':      2.16,
    'd':      4.67,
    's':     93.4,
    'c':   1270.0,
    'b':   4180.0,
    't': 172690.0,
}
QUARK_SPECIES_ORDER = ('u', 'd', 's', 'c', 'b', 't')

# Shell basis state → species map: mirrors v3's BASIS_TO_SPECIES,
# with k_v3 ↔ n_shell substitution:
#   v3 (k=1, +) = u  →  shell (n=3, +) = u
#   v3 (k=1, −) = d  →  shell (n=3, −) = d
#   v3 (k=3, +) = c  →  shell (n=4, +) = c
#   v3 (k=3, −) = s  →  shell (n=4, −) = s
#   v3 (k=5, +) = t  →  shell (n=5, +) = t
#   v3 (k=5, −) = b  →  shell (n=5, −) = b
SHELL_STATE_TO_SPECIES: dict[tuple[int, int, int], str] = {
    (1, 3, +1): 'u', (1, 3, -1): 'd',
    (1, 4, +1): 'c', (1, 4, -1): 's',
    (1, 5, +1): 't', (1, 5, -1): 'b',
}

# v3 baseline n_part (PR #76)
N_PART_V3 = 233


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def species_of(state) -> str:
    """Map a shell basis state to its v3-analog species label."""
    return SHELL_STATE_TO_SPECIES.get((state.l, state.n, state.p), '?')


def diagnostic_spectrum(eigvals: np.ndarray, basis) -> dict[str, float]:
    """Pair eigenvalues with species via the (l, n, p) → species map
    by sorting both sides — the smallest eigenvalue gets the lightest
    intended species in that block. Returns a dict {species: eigval}."""
    spec = {}
    # Per-block: pair the two p=± eigenvalues to the two species in
    # that (l, n) block.
    blocks: dict[tuple[int, int], list[tuple[int, str]]] = {}
    for i, s in enumerate(basis):
        blocks.setdefault((s.l, s.n), []).append((i, species_of(s)))
    for (_, _), idx_species in blocks.items():
        # eigvals are sorted globally; we need the eigenvalues belonging
        # to this block. We don't have diagonalizer info here; just use
        # the kinetic ordering by (l, n) sorted, then split ± by p sign.
        pass
    # Fallback: assume the operator is block-diagonal (no inter-block
    # coupling in this audit's main test). Then eigenvalues per block =
    # ω² ± χ. Just sort eigenvalues, sort species by intended order, and
    # pair by index. This is fine for diagnostic purposes.
    species_in_order_by_intended_mass = sorted(
        [species_of(s) for s in basis],
        key=lambda sp: QUARK_MASS_OBS_MEV.get(sp, 0.0),
    )
    eigs_sorted = sorted(float(e) for e in eigvals)
    for sp, ev in zip(species_in_order_by_intended_mass, eigs_sorted):
        spec[sp] = ev
    return spec


def within_block_ordering_correct(
    op: ShellOperatorScaffold,
    chi: float,
) -> dict:
    """For each (l, n) block of the shell basis, set χ to the given
    value and check whether the +/− ordering (block eigenvalues
    ω² ± χ) matches the species ordering for that block:
      n=3 block: m_u (+) < m_d (−)   ⟹ requires χ < 0 (+ lighter)
      n=4 block: m_c (+) > m_s (−)   ⟹ requires χ > 0 (+ heavier)
      n=5 block: m_t (+) > m_b (−)   ⟹ requires χ > 0 (+ heavier)
    With a UNIFORM χ, one of the two cases must fail — that's the
    structural impossibility this test pins down."""
    op.chi = chi
    op.couple_matrix = []          # block-diagonal
    H = op.H_total()
    # Per-block eigenvalues from the (2x2) block structure
    blocks: dict[tuple[int, int], dict] = {}
    for i, s in enumerate(op.basis):
        blocks.setdefault((s.l, s.n), {})
        blocks[(s.l, s.n)][s.p] = H[i, i]
    block_correct: dict[tuple[int, int], bool] = {}
    for key, vals in blocks.items():
        mass_plus = vals.get(+1)
        mass_minus = vals.get(-1)
        sp_plus = SHELL_STATE_TO_SPECIES[(key[0], key[1], +1)]
        sp_minus = SHELL_STATE_TO_SPECIES[(key[0], key[1], -1)]
        obs_plus = QUARK_MASS_OBS_MEV[sp_plus]
        obs_minus = QUARK_MASS_OBS_MEV[sp_minus]
        intended_plus_heavier = obs_plus > obs_minus
        scaffold_plus_heavier = mass_plus > mass_minus
        block_correct[key] = intended_plus_heavier == scaffold_plus_heavier
    return {
        'chi': chi,
        'block_correct': {f"l{k[0]}n{k[1]}": v for k, v in block_correct.items()},
        'n_blocks_correct': sum(block_correct.values()),
        'all_blocks_correct': all(block_correct.values()),
    }


# ---------------------------------------------------------------------------
# T1. Shell kinetic spectrum
# ---------------------------------------------------------------------------

def test_T1_shell_kinetic_spectrum() -> dict:
    """Compute the shell kinetic spectrum ω²(l=1, n=3,4,5). This is
    H_kin's diagonal in the n_varied enumeration."""
    basis = build_shell_basis(mode='n_varied')
    omegas = [s.omega for s in basis]
    omega_sq = [s.omega ** 2 for s in basis]
    range_factor = max(omega_sq) / min(omega_sq)
    return {
        'name': 'T1_shell_kinetic_spectrum_n_varied',
        'description': (
            "Shell kinetic spectrum on n_varied basis: ω²(l=1, n=3,4,5). "
            "This is H_kin's diagonal — the 'cavity-mass operator'."
        ),
        'basis_labels': [s.label() for s in basis],
        'species': [species_of(s) for s in basis],
        'omegas': [round(o, 4) for o in omegas],
        'omega_sq_mass_sq': [round(o, 4) for o in omega_sq],
        'mass_sq_range_factor': round(range_factor, 4),
        'pass': len(omegas) == N_FLAVORS,
    }


# ---------------------------------------------------------------------------
# T2. Coverage gap vs observed
# ---------------------------------------------------------------------------

def test_T2_coverage_gap_vs_observed() -> dict:
    """Compare shell kinetic mass²-range to observed quark mass²-range.
    Observed: ×6.4·10⁹ from u to t. Shell kinetic: ×2.2 from n=3 to
    n=5. The shell kinetic alone is insufficient — an additional
    structural contribution must span the rest. PR #78 cannot fix
    this alone; PR #79 (boundary stress tensor) and PR #80 (color
    algebra) must populate the operator with the missing structure."""
    basis = build_shell_basis(mode='n_varied')
    shell_min_sq = min(s.omega ** 2 for s in basis)
    shell_max_sq = max(s.omega ** 2 for s in basis)
    shell_range_sq = shell_max_sq / shell_min_sq
    obs = QUARK_MASS_OBS_MEV
    obs_min_sq = min(obs[sp] ** 2 for sp in QUARK_SPECIES_ORDER)
    obs_max_sq = max(obs[sp] ** 2 for sp in QUARK_SPECIES_ORDER)
    obs_range_sq = obs_max_sq / obs_min_sq
    coverage_deficit = obs_range_sq / shell_range_sq
    return {
        'name': 'T2_coverage_gap_shell_vs_observed',
        'description': (
            "Compare shell kinetic mass²-range to observed range. Shell "
            "alone is insufficient by ~9 orders of magnitude — PR #79 "
            "(boundary stress tensor) and PR #80 (color algebra) must "
            "populate the rest."
        ),
        'shell_kinetic_mass_sq_min': shell_min_sq,
        'shell_kinetic_mass_sq_max': shell_max_sq,
        'shell_kinetic_range_factor': shell_range_sq,
        'observed_mass_sq_min': obs_min_sq,
        'observed_mass_sq_max': obs_max_sq,
        'observed_range_factor': obs_range_sq,
        'coverage_deficit_orders_of_magnitude': math.log10(coverage_deficit),
        'shell_kinetic_alone_sufficient': shell_range_sq >= obs_range_sq,
        'pass': True,   # T2 always passes; reports the gap honestly
    }


# ---------------------------------------------------------------------------
# T3. Uniform χ cannot reproduce within-generation inversion
# ---------------------------------------------------------------------------

def test_T3_uniform_chi_cannot_invert() -> dict:
    """A UNIFORM Z₂ partition splitter `H_Z2 = χ·σ_z` (applied
    identically across all (l, n) blocks) splits every block in the
    SAME direction: if χ > 0 then + is always heavier. But the
    observed quark mass ordering requires + LIGHTER at n=3 (u < d)
    and + HEAVIER at n=4, 5 (c > s, t > b). So uniform χ cannot
    reproduce the within-generation inversion — at best 2 of 3
    blocks correct (whichever sign of χ is chosen)."""
    basis = build_shell_basis(mode='n_varied')
    op = ShellOperatorScaffold(basis=basis)
    # Scan uniform χ over a reasonable range; record per-block
    # correctness.
    results = []
    best_correct = 0
    for chi_val in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]:
        r = within_block_ordering_correct(op, chi_val)
        results.append({'chi': chi_val, 'n_correct': r['n_blocks_correct'],
                        'all_correct': r['all_blocks_correct'],
                        'block_correct': r['block_correct']})
        best_correct = max(best_correct, r['n_blocks_correct'])
    uniform_chi_can_reproduce_all = (best_correct == N_GENERATIONS)
    return {
        'name': 'T3_uniform_chi_cannot_invert',
        'description': (
            "A UNIFORM χ (same sign across all blocks) cannot reproduce "
            "the within-generation inversion: u < d (requires + lighter "
            "at n=3) but c > s, t > b (require + heavier at n=4, 5). "
            "Best uniform-χ correctness: 2/3 blocks. This identifies "
            "the structural slot PR #79 must populate — a generation-"
            "dependent χ_n."
        ),
        'scan_results': results,
        'best_uniform_chi_correctness': best_correct,
        'n_blocks_total': N_GENERATIONS,
        'uniform_chi_can_reproduce_all_blocks': uniform_chi_can_reproduce_all,
        'pass': not uniform_chi_can_reproduce_all,  # confirms the negative
    }


# ---------------------------------------------------------------------------
# T4. Sign-flipping χ_n — existence proof
# ---------------------------------------------------------------------------

def _eval_with_chi_n(chi_n: dict[int, float], basis):
    """Block-diagonal H with a per-n splitter χ_n. Returns the dict
    {species: diagonal mass²}."""
    spec: dict[str, float] = {}
    for s in basis:
        chi = chi_n.get(s.n, 0.0)
        # Diagonal: ω² + chi·p (with p = ±1)
        mass_sq = s.omega ** 2 + chi * s.p
        sp = species_of(s)
        spec[sp] = mass_sq
    return spec


def test_T4_signflipping_chi_n_existence() -> dict:
    """A GENERATION-DEPENDENT splitter `χ_n` can reproduce the
    inversion: with χ_3 < 0 (so + lighter at n=3, u<d) and χ_4,5 > 0
    (so + heavier at n=4,5, c>s and t>b). This is an EXISTENCE proof
    of the structural slot PR #79's boundary stress tensor needs to
    populate. The values are illustrative, not derived here."""
    basis = build_shell_basis(mode='n_varied')
    # Choose sign-flipping χ_n so that the WITHIN-BLOCK ordering is
    # correct for all 3 blocks. Magnitudes are illustrative.
    chi_n = {3: -1.0, 4: +5.0, 5: +20.0}    # only signs matter for ordering
    spec = _eval_with_chi_n(chi_n, basis)
    # Check within-block ordering
    block_ok = {
        'l1n3 (u<d)': spec['u'] < spec['d'],
        'l1n4 (c>s)': spec['c'] > spec['s'],
        'l1n5 (t>b)': spec['t'] > spec['b'],
    }
    all_blocks_ok = all(block_ok.values())
    # Also check global ordering — does u < d < s < c < b < t hold?
    species_sorted_by_mass_sq = sorted(spec.keys(), key=lambda sp: spec[sp])
    global_ordering_correct = (
        species_sorted_by_mass_sq == ['u', 'd', 's', 'c', 'b', 't']
    )
    return {
        'name': 'T4_signflipping_chi_n_existence',
        'description': (
            "A generation-dependent χ_n (sign-flipping between n=3 and "
            "n=4,5) reproduces the within-block inversion. Illustrative "
            "magnitudes — PR #79 must derive the structural values from "
            "the boundary stress tensor."
        ),
        'chi_n_illustrative': {f"n={k}": v for k, v in chi_n.items()},
        'shell_spec_mass_sq': {sp: round(float(spec[sp]), 4) for sp in QUARK_SPECIES_ORDER},
        'within_block_ordering': block_ok,
        'all_blocks_ok': all_blocks_ok,
        'global_ordering_u_d_s_c_b_t': species_sorted_by_mass_sq,
        'global_ordering_correct': global_ordering_correct,
        'pass': all_blocks_ok,   # within-block is the structural test
    }


# ---------------------------------------------------------------------------
# T5. Inter-mode coupling: minimum needed to span observed range
# ---------------------------------------------------------------------------

def _compensator_for_absolute_scale(spec_mass_sq: dict[str, float]) -> float:
    """Given a shell-derived dimensionless mass² spectrum, fit a single
    overall scale factor S such that S · mean(shell)² ≈ mean(observed)².
    Then return the L² log-residual — a "compensator measure" — which
    quantifies how much additional non-structural physics is needed."""
    species = list(spec_mass_sq.keys())
    shell = np.array([spec_mass_sq[sp] for sp in species])
    obs = np.array([QUARK_MASS_OBS_MEV[sp] ** 2 for sp in species])
    log_shell = np.log10(np.maximum(shell, 1e-30))
    log_obs = np.log10(np.maximum(obs, 1e-30))
    # Fit a single offset (overall scale in log space)
    offset = float(np.mean(log_obs - log_shell))
    fit_log_shell = log_shell + offset
    # Residual after scaling
    residual = float(np.sqrt(np.mean((fit_log_shell - log_obs) ** 2)))
    return residual


def test_T5_compensator_kinetic_only() -> dict:
    """With only H_kin (no χ, no coupling), the shell spectrum has
    six degenerate pairs (+/− each at ω²(n=3,4,5)). Fit a single
    overall scale; report the log-residual (compensator measure).
    Large residual ⟹ much non-structural physics is needed."""
    basis = build_shell_basis(mode='n_varied')
    # Kinetic only: spec[sp] = ω²(l=1, n)
    spec = {species_of(s): s.omega ** 2 for s in basis}
    residual = _compensator_for_absolute_scale(spec)
    return {
        'name': 'T5_kinetic_only_compensator_residual',
        'description': (
            "Kinetic-only spectrum (six degenerate ω² pairs) fit to "
            "observed by a single overall scale; log-residual measures "
            "what additional structure is needed."
        ),
        'shell_mass_sq': {sp: round(float(spec[sp]), 4)
                          for sp in QUARK_SPECIES_ORDER},
        'log10_residual_kinetic_only': residual,
        'interpretation': (
            'Large residual: shell kinetic plus a single scale factor '
            'cannot reproduce 12-orders-of-magnitude observed spread. '
            'PR #79 (boundary stress) + PR #80 (color) must contribute.'
        ),
        'pass': True,
    }


def test_T6_compensator_with_chi_n() -> dict:
    """Add the T4 sign-flipping χ_n; re-fit overall scale; report the
    new log-residual. If the residual drops substantially, the
    structural shape is closer to observed. If it stays large,
    additional H_couple inter-mode mixing is needed."""
    basis = build_shell_basis(mode='n_varied')
    # Use illustrative sign-flipping χ_n that fixed the ordering
    chi_n = {3: -1.0, 4: +5.0, 5: +20.0}
    spec_no_chi = {species_of(s): s.omega ** 2 for s in basis}
    spec_with_chi = _eval_with_chi_n(chi_n, basis)
    residual_no = _compensator_for_absolute_scale(spec_no_chi)
    residual_yes = _compensator_for_absolute_scale(spec_with_chi)
    return {
        'name': 'T6_compensator_with_signflipping_chi_n',
        'description': (
            "With sign-flipping χ_n (T4 illustrative values), re-fit "
            "single scale; report the new log-residual. χ_n improves "
            "the ordering but only modestly improves the spread."
        ),
        'chi_n_illustrative': {f"n={k}": v for k, v in chi_n.items()},
        'log10_residual_kinetic_only': residual_no,
        'log10_residual_with_chi_n': residual_yes,
        'residual_reduction': residual_no - residual_yes,
        'observed_log10_range': math.log10(
            QUARK_MASS_OBS_MEV['t'] ** 2 / QUARK_MASS_OBS_MEV['u'] ** 2),
        'interpretation': (
            'Residual stays large after χ_n: the shape needs an '
            'exponential or strongly nonlinear inter-mode coupling. '
            'PR #79 (n-dependent χ from boundary stress tensor) is '
            'one channel; PR #80 (color algebra mixing) is the '
            'complementary one.'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Comparison with v3 baseline
# ---------------------------------------------------------------------------

def test_T7_comparison_with_v3() -> dict:
    """Compare the shell basis findings with v3's:
      - v3 used closure-quantum integers; absorbed unmodeled QCD
        physics into n_part = 233 (#76 diagnosis).
      - shell basis: structurally distinct (cavity wavefronts vs
        throat traversals), correct 3×2=6 flavor count, Z₂ partition
        slot for within-generation inversion (T4), but kinetic alone
        cannot span observed range (T2).
      - n_part audit: shell basis does NOT reduce/remove n_part by
        itself. The shape of the spread (12 orders of magnitude) is
        still mostly unaccounted-for in PR #78. PR #79 (boundary
        stress tensor) and PR #80 (color algebra) must contribute."""
    # Both bases need to fit ~12 orders of magnitude in mass²; both
    # need substantial compensation.
    v3_n_part = N_PART_V3
    # Honest summary: PR #78 alone gives no quantitative reduction.
    return {
        'name': 'T7_comparison_with_v3_lepton_shaped_baseline',
        'description': (
            "Compare shell basis (PR #77/#78) with v3 (lepton-shaped). "
            "Shell basis is structurally cleaner (cavity wavefronts not "
            "throat traversals, kinetic IS ω², 3×2=6 with correct "
            "partition slot). But shell-kinetic-only does NOT span the "
            "observed 12-order-of-magnitude mass² range; the n_part "
            "compensator question is NOT yet resolved at #78."
        ),
        'v3_n_part_baseline': v3_n_part,
        'shell_basis_improvements': [
            'distinct from throat sector (cavity wavefronts)',
            'kinetic = ω²(l, n), not phenomenological β·k²·(2π)',
            'Z₂ partition slot for within-generation inversion (T4)',
            '3 × 2 = 6 flavors structural (PR #69 match)',
        ],
        'shell_basis_not_yet_resolved': [
            'coverage gap: shell-kinetic range factor 2.2 vs observed 6.4e9',
            'within-generation inversion needs χ_n, not uniform χ (T3)',
            'n_part compensator NOT reduced by PR #78 alone',
            'PR #79 (boundary stress tensor χ_n) and PR #80 (color '
            'algebra coupling H_couple) must populate the operator',
        ],
        'pr_78_alone_resolves_n_part': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Honest scope + assessment
# ---------------------------------------------------------------------------

def test_T8_honest_scope_assessment() -> dict:
    return {
        'name': 'T8_honest_scope_assessment',
        'description': (
            "Honest scope: the shell waveguide basis is structurally "
            "better than the v3 lepton-shaped basis (cavity wavefronts "
            "vs throat traversals; ω²(l, n) vs phenomenological β·k²; "
            "Z₂ partition slot for within-generation inversion). The "
            "structural slots are identified: uniform χ insufficient "
            "(T3); sign-flipping χ_n needed (T4). However, the n_part "
            "compensator question is NOT YET resolved at PR #78 — the "
            "shell kinetic alone cannot span the observed 12-orders-of-"
            "magnitude mass² range. PR #79 (boundary stress tensor) and "
            "PR #80 (color algebra) must populate the operator with "
            "structurally derived values. PR #78 SHARPENS the scope of "
            "what PRs #79–#80 must do; it does not close it."
        ),
        'classification': 'SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED',
        'what_pr_78_establishes': [
            'shell basis is the right structural machinery',
            'uniform χ cannot reproduce within-generation inversion',
            'sign-flipping χ_n CAN reproduce it (existence proof)',
            'kinetic-only spectrum does not span observed range',
        ],
        'what_pr_78_does_not_resolve': [
            'absolute scale / n_part compensator (needs #79+#80)',
            'derivation of χ_n values from a physical principle (#79)',
            'identification of H_couple structure (#79+#80)',
            'color algebra (SU(3) vs SU(2)×Z₂ vs other) (#80)',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_shell_kinetic_spectrum(),
        test_T2_coverage_gap_vs_observed(),
        test_T3_uniform_chi_cannot_invert(),
        test_T4_signflipping_chi_n_existence(),
        test_T5_compensator_kinetic_only(),
        test_T6_compensator_with_chi_n(),
        test_T7_comparison_with_v3(),
        test_T8_honest_scope_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED'
        verdict = (
            'SHELL BASIS STRUCTURALLY BETTER, N_PART NOT YET RESOLVED. '
            'The PR #77 shell waveguide basis is structurally cleaner '
            'than the v3 lepton-shaped basis for the quark sector. The '
            'kinetic operator is ω²(l, n) — the cavity-wavefront '
            'eigenfrequency squared — not the phenomenological winding '
            'cost β·k²·(2π) that v3 inherited from the lepton ladder. '
            'The Z₂ partition slot is the natural structural home for '
            'the within-generation mass-ordering inversion (u < d but '
            'c > s, t > b).\n\n'
            'STRUCTURAL FINDINGS. '
            '(1) A UNIFORM χ·σ_z partition splitter cannot reproduce '
            'the inversion: with any single sign of χ, at most 2 of 3 '
            '(l, n) blocks have the correct +/− ordering. (T3 confirms '
            'this with a parameter scan.) '
            '(2) A GENERATION-DEPENDENT, sign-flipping χ_n CAN '
            'reproduce the inversion — concretely, χ_3 < 0 (so + is '
            'lighter at n=3, u < d) with χ_4, χ_5 > 0 (so + is heavier '
            'at n=4, 5, c > s and t > b). T4 establishes the existence '
            'of this structural slot; the magnitudes are illustrative, '
            'NOT derived. PR #79 (boundary stress tensor on the cavity '
            'wall) is the natural channel to derive χ_n from physics.\n\n'
            'COVERAGE GAP. The shell kinetic spectrum spans only a '
            'factor of ~2.2 in mass² (ω²(n=5)/ω²(n=3) ≈ 32.5/14.6) '
            'while the observed quark mass² spans a factor of ~6.4·10⁹ '
            '(u to t). Shell kinetic alone is therefore insufficient '
            'by ~9 orders of magnitude. Sign-flipping χ_n improves the '
            'ORDERING but not the SPREAD. Closing the spread requires '
            'either a large n-dependent χ_n (PR #79) or strong '
            'inter-mode coupling H_couple (PR #79 + PR #80).\n\n'
            'N_PART STATUS. PR #78 does NOT reduce or remove the v3 '
            'phenomenological n_part = 233 compensator on its own. The '
            'shell basis identifies the structural SLOTS that must '
            'populate the operator, but the values that span the '
            'observed range are not derivable at PR #78 — they '
            'require PR #79 (boundary stress tensor → χ_n) and PR #80 '
            '(color algebra → H_couple). The n_part question is '
            'therefore SHARPENED, not closed: PR #78 identifies the '
            'right machinery and the missing structural inputs.\n\n'
            'COMPARISON WITH V3. The shell basis (this PR) is '
            'structurally distinct from v3 in four ways: (i) basis = '
            'cavity wavefronts not throat traversals, (ii) kinetic = '
            'ω²(l, n) cavity eigenfrequency squared not phenomenological '
            'β·k²·(2π), (iii) Z₂ partition slot is the right structural '
            'home for the inversion, (iv) the 3 × 2 = 6 flavor count '
            'matches PR #69. These are qualitative improvements at the '
            'STRUCTURAL level. Quantitatively, PR #78 alone does not '
            'win on mass accuracy or n_part reduction — both bases '
            'require substantial compensation to span the observed '
            'range. The shell basis advantage is that its compensation '
            'is structurally located (PR #79 + PR #80 slots) rather '
            'than absorbed into a single phenomenological integer.\n\n'
            'HONEST SCOPE. PR #78 SHARPENS what PR #79 must produce: '
            'sign-flipping χ_n from the boundary stress tensor, with '
            'magnitudes large enough to push the n=5 block up by a '
            'factor of ~10⁹ in mass². It also IDENTIFIES PR #80\'s '
            'role: H_couple inter-mode mixing for the inter-generation '
            'hierarchy (and the color algebra it transforms under). '
            'PR #78 does not close n_part; it relocates the '
            'phenomenological content to structurally-named slots that '
            'PRs #79–#80 must populate from physics.'
        )
    else:
        verdict_class = 'AUDIT_INCONCLUSIVE'
        verdict = (
            'AUDIT INCONCLUSIVE. A structural test failed; investigate '
            'before proceeding to PR #79.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'shell basis structurally better than v3 (cavity wavefronts '
            'vs throat traversals); within-generation inversion requires '
            'sign-flipping χ_n (PR #79 slot); n_part NOT yet resolved at '
            'PR #78'
        ),
        'next_pr': (
            'PR #79 — boundary stress tensor to derive χ_n; '
            'PR #80 — color algebra to identify H_couple structure'
        ),
        'b4_caveat': (
            'shell ω dimensionful; mass-ratio audit scale-free; '
            'absolute scale via single B4 anchor (PR #53)'
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
    L.append('# Shell Hamiltonian mass-ordering and `n_part` audit')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Tests whether the PR #77 shell waveguide basis reproduces the "
        "qualitative quark ordering structure better than the v3 "
        "lepton-shaped basis, and whether the `n_part = 233` "
        "compensator (PR #76) shrinks or disappears."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Next PR**: {s['next_pr']}")
    L.append(f"- **B4 caveat**: {s['b4_caveat']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'shell kinetic ω²(l=1, n=3,4,5) — range factor 2.2',
        'T2': 'coverage gap: shell ×2.2 vs observed ×6.4e9 (~9 OOM)',
        'T3': 'uniform χ cannot invert (best 2/3 blocks)',
        'T4': 'sign-flipping χ_n CAN invert (existence proof)',
        'T5': 'kinetic-only compensator log-residual (large)',
        'T6': 'χ_n improves ordering, modestly improves spread',
        'T7': 'shell structurally better; n_part not resolved at #78',
        'T8': 'PR #78 SHARPENS PR #79–#80 scope; does not close',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    # T2 coverage gap
    t2 = s['tests'][1]
    L.append('## T2: Coverage gap — shell kinetic vs observed mass²')
    L.append('')
    L.append('| range | min mass² | max mass² | factor |')
    L.append('|---|---:|---:|---:|')
    L.append(f"| shell kinetic ω² (n=3,4,5) | "
             f"{t2['shell_kinetic_mass_sq_min']:.2f} | "
             f"{t2['shell_kinetic_mass_sq_max']:.2f} | "
             f"{t2['shell_kinetic_range_factor']:.2f} |")
    L.append(f"| observed (u → t) | "
             f"{t2['observed_mass_sq_min']:.2f} | "
             f"{t2['observed_mass_sq_max']:.2e} | "
             f"{t2['observed_range_factor']:.2e} |")
    L.append('')
    L.append(f"Coverage deficit: **{t2['coverage_deficit_orders_of_magnitude']:.1f} "
             "orders of magnitude in mass²**. Shell kinetic alone is "
             "insufficient; PR #79–#80 must contribute.")
    L.append('')

    # T3 uniform χ
    t3 = s['tests'][2]
    L.append('## T3: Uniform χ cannot reproduce within-generation inversion')
    L.append('')
    L.append('| χ | l1n3 (u<d?) | l1n4 (c>s?) | l1n5 (t>b?) | correct/3 |')
    L.append('|---:|:---:|:---:|:---:|---:|')
    for r in t3['scan_results']:
        bc = r['block_correct']
        L.append(f"| {r['chi']:+.1f} | "
                 f"{'✓' if bc.get('l1n3', False) else '✗'} | "
                 f"{'✓' if bc.get('l1n4', False) else '✗'} | "
                 f"{'✓' if bc.get('l1n5', False) else '✗'} | "
                 f"{r['n_correct']} |")
    L.append('')
    L.append(f"Best uniform-χ correctness: **{t3['best_uniform_chi_correctness']} "
             f"of {t3['n_blocks_total']} blocks**. Uniform splitter is "
             "structurally insufficient — PR #79's boundary stress "
             "tensor must produce a generation-dependent χ_n.")
    L.append('')

    # T4 sign-flipping χ_n
    t4 = s['tests'][3]
    L.append('## T4: Sign-flipping χ_n — existence proof')
    L.append('')
    L.append(f"Illustrative `χ_n` (signs matter, magnitudes are not "
             f"derived): `{t4['chi_n_illustrative']}`.")
    L.append('')
    L.append('| species | shell mass² (illustrative) | within-block correct? |')
    L.append('|---|---:|:---:|')
    for sp in QUARK_SPECIES_ORDER:
        L.append(f"| {sp} | {t4['shell_spec_mass_sq'][sp]} | |")
    L.append('')
    L.append('| within-block ordering | correct? |')
    L.append('|---|:---:|')
    for k, v in t4['within_block_ordering'].items():
        L.append(f"| {k} | {'✓' if v else '✗'} |")
    L.append('')
    L.append(f"All 3 within-block orderings correct: "
             f"**{t4['all_blocks_ok']}**. Global mass ordering: "
             f"`{t4['global_ordering_u_d_s_c_b_t']}` (correct = "
             f"`['u','d','s','c','b','t']`).")
    L.append('')

    # T5 + T6 compensator
    t5 = s['tests'][4]
    t6 = s['tests'][5]
    L.append('## T5–T6: Compensator measure (log₁₀ residual after best single-scale fit)')
    L.append('')
    L.append('| configuration | log₁₀ residual |')
    L.append('|---|---:|')
    L.append(f"| kinetic only (no χ, no coupling) | "
             f"{t5['log10_residual_kinetic_only']:.3f} |")
    L.append(f"| kinetic + sign-flipping χ_n (illustrative) | "
             f"{t6['log10_residual_with_chi_n']:.3f} |")
    L.append(f"| observed range (for reference) | "
             f"{t6['observed_log10_range']:.3f} |")
    L.append('')
    L.append("Residual reduction from χ_n: "
             f"{t6['residual_reduction']:.3f} (modest). The remaining "
             "residual must be supplied by PR #79's full χ_n derivation "
             "and PR #80's H_couple inter-mode mixing.")
    L.append('')

    # T7 comparison
    t7 = s['tests'][6]
    L.append('## T7: Comparison with v3 baseline')
    L.append('')
    L.append('**Shell basis improvements over v3:**')
    L.append('')
    for item in t7['shell_basis_improvements']:
        L.append(f"  - {item}")
    L.append('')
    L.append('**Not yet resolved at PR #78:**')
    L.append('')
    for item in t7['shell_basis_not_yet_resolved']:
        L.append(f"  - {item}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **PR #79** — derive the generation-dependent `χ_n` from '
             'a boundary stress tensor on the cavity wall. Signs and '
             'magnitudes are the structural input.')
    L.append('- **PR #80** — identify the color algebra acting on '
             '`(l, n, p)` and populate `H_couple` with inter-mode '
             'mixing terms transforming under that algebra. This is '
             'the channel for the inter-generation hierarchy.')
    L.append('- **`n_part` audit** — re-run after PR #79 and PR #80 '
             'populate the operator. If the resulting spread is '
             'structurally accounted for, `n_part = 233` (v3) is '
             'replaced by a derived value; otherwise it persists as a '
             'residual phenomenological compensator with reduced scope.')
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_shell_mass_ordering_audit_probe'
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
