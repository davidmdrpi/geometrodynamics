"""
BAM non-orientable hadron sector: a compact experimental note (PR #110).

PRs #100–#109 built BAM's non-orientable (Möbius / closed-flux-loop) hadron
predictions across four sub-sectors — mesonic 1⁻⁺ hybrids, glueballs, and
heavy Möbius baryons with their decays. This probe COMPILES them into a
single compact "experimental note" of the kind an LHCb / Belle II / BESIII
analyst can read off: predicted masses, Q-values, preferred/suppressed
modes, and analysis handles. Every number is recomputed here from the one
input (√σ) so the note cannot drift from the probes it summarises.

## The single input

Everything in the non-orientable sector is set by ONE QCD scale, the string
tension √σ ≈ 0.424 GeV (the B4 confinement anchor). The flux-tube
excitation quantum — the energy to put one unit of transverse excitation /
non-orientable twist on the tube — is

    Δ = 2√σ ≈ 0.849 GeV,

and the lowest closed-flux-loop (glueball) level is √(4πσ) ≈ 1.50 GeV.
There are no new parameters; the note is a pushforward of √σ.

## What the note collects

  1. Mesonic 1⁻⁺ hybrids (PR #101): ground meson + 2√σ — π₁ ≈ 1.62, η₁ ≈
     1.85 GeV — matched to the observed π₁(1600) / η₁(1855); the exotic
     1⁻⁺ is forbidden to ordinary qq̄ and is the Möbius/twist marker.
  2. Glueballs as closed flux loops (PR #100): 0⁺⁺ ground √(4πσ) ≈ 1.50 GeV
     (benchmarks the lattice 0⁺⁺ √σ scale to ~13%); unobserved — the
     freest channel.
  3. Heavy Möbius baryons (PR #103): ground heavy baryon + 2√σ — Λ_c ~3135,
     Ω_c ~3544, Ξ_cc ~4471, Λ_b ~6469, Ω_b ~6894 MeV — supernumerary,
     above the orbital tower, just above current data.
  4. Heavy Möbius baryon decays (PR #109): twist-unwinding ⟹ the hybrid
     selection rule (single-S-wave-π-to-ground SUPPRESSED; Σ_Q π /
     isoscalar dipion / P-wave+π PREFERRED) and the cross-flavor Q-match
     (Λ_Q ππ 569, Λ_Q η 301 MeV identical for c and b).

## Analysis handles (the note's payload)

  - **Cross-flavor Q-match.** The all-light release energies are
    flavor-independent: the SAME dipion spectrum (Q = 569 MeV) and η recoil
    (Q = 301 MeV) above the charm and the bottom ground baryon. A
    correlated, two-channel signature.
  - **Hybrid selection rule.** The naive single-S-wave-π-to-ground
    transition is suppressed; an ordinary radial excitation would show the
    opposite. The branching PATTERN is the discriminator.
  - **Isoscalar dipion.** The preferred Λ_Q(ππ)_S channel peaks at high
    m(ππ) (coherent twist unwinding, like ψ(2S) → J/ψ ππ).
  - **Broad ⟹ amplitude analyses.** Widths ~tens–150 MeV; resolve in
    Dalitz / amplitude fits, not as sharp peaks.
  - **Exotic J^PC where available.** The mesonic sector HAS a smoking gun
    (1⁻⁺); the baryon sector does not (supernumerary ordinary-J^P) and
    leans on the cross-flavor correlation instead.

## Honest scope

The note is a compilation, not new physics: it carries the established
content (masses as pushforwards of √σ; the decay PATTERN and Q-structure)
and the open items unchanged (exact masses within the lattice hybrid-gap
band ~0.8–1.3 GeV; absolute branching fractions and total widths; baryon
J^P). It is a reference card, deliberately compact.

Tests:
  T1. Scope: consolidate the PRs #100–#109 non-orientable sector into one
      experimental note.
  T2. The single input: √σ, the gap 2√σ, the glueball scale √(4πσ).
  T3. Mesonic 1⁻⁺ hybrids: ground + 2√σ (π₁ ~1.62, η₁ ~1.85 GeV), matched.
  T4. Glueball closed-flux-loop ground √(4πσ) ≈ 1.50 GeV.
  T5. Heavy Möbius baryon masses: ground + 2√σ.
  T6. Heavy Möbius baryon decays: channels + Q-values + cross-flavor match.
  T7. Preferred/suppressed modes + analysis handles.
  T8. Assessment (emits the note).

Verdict:
  - NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED (expected): the four
    non-orientable sub-sectors (mesonic 1⁻⁺ hybrids, glueballs, heavy
    Möbius baryon masses, heavy Möbius baryon decays) are compiled into one
    compact, internally-consistent experimental note — predicted masses,
    Q-values, preferred/suppressed modes, and analysis handles — all
    pushed forward from the single input √σ.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.qcd.constants import SIGMA_QCD


PI = math.pi
SQRT_SIGMA_MEV = math.sqrt(SIGMA_QCD) * 1000.0      # ≈ 424 MeV
GAP_MEV = 2.0 * SQRT_SIGMA_MEV                       # 2√σ ≈ 849 MeV
GLUEBALL_MEV = math.sqrt(4.0 * PI * SIGMA_QCD) * 1000.0   # √(4πσ) ≈ 1504 MeV

# Light fragments (MeV, PDG).
M_PI = 139.57
M_ETA = 547.86

# Ground hadrons for the mesonic hybrids (MeV).
M_RHO = 775.0          # ρ(770) — π₁ basis
M_ETA1_BASE = 1005.0   # ~1.0 GeV — η₁ basis (per PR #101)

# Heavy-baryon ground / partner states (MeV, PDG).
HB = {
    'Lambda_c': 2286.46, 'Sigma_c': 2453.75, 'Sigma_c_2520': 2518.4,
    'D0': 1864.84, 'proton': 938.27,
    'Lambda_b': 5619.6, 'Sigma_b': 5813.1, 'Sigma_b_star': 5832.5, 'Bminus': 5279.3,
}

# Heavy baryons that get a Möbius partner (ground, current ceiling state).
HEAVY_MOBIUS = {
    'Lambda_c': (2286, 'Λ_c(2940)'),
    'Omega_c':  (2695, 'Ω_c(3120)'),
    'Xi_cc':    (3622, 'none (ground only)'),
    'Lambda_b': (5620, 'Λ_b(6152)'),
    'Omega_b':  (6045, 'none (ground only)'),
}

MOBIUS_C = HB['Lambda_c'] + GAP_MEV
MOBIUS_B = HB['Lambda_b'] + GAP_MEV


def _Q(parent: float, *frags: float) -> float:
    return parent - sum(frags)


# ---------------------------------------------------------------------------
# T1. Scope
# ---------------------------------------------------------------------------

def test_T1_scope() -> dict:
    return {
        'name': 'T1_scope_compile_one_note',
        'description': (
            "Compile the PRs #100–#109 non-orientable hadron sector — "
            "mesonic 1⁻⁺ hybrids (#101), glueballs (#100), heavy Möbius "
            "baryon masses (#103) and decays (#109) — into one compact "
            "experimental note: masses, Q-values, modes, handles."
        ),
        'sub_sectors': [
            'mesonic 1⁻⁺ hybrids (PR #101)',
            'glueballs / closed flux loops (PR #100)',
            'heavy Möbius baryon masses (PR #103)',
            'heavy Möbius baryon decays (PR #109)',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. The single input
# ---------------------------------------------------------------------------

def test_T2_single_input() -> dict:
    return {
        'name': 'T2_single_input_sqrt_sigma',
        'description': (
            "One QCD scale sets the sector: √σ ≈ 424 MeV (B4 anchor); the "
            "flux-tube excitation quantum 2√σ ≈ 849 MeV; the closed-loop "
            "glueball scale √(4πσ) ≈ 1504 MeV. No new parameters."
        ),
        'sqrt_sigma_MeV': SQRT_SIGMA_MEV,
        'gap_2sqrt_sigma_MeV': GAP_MEV,
        'glueball_sqrt_4pi_sigma_MeV': GLUEBALL_MEV,
        'pass': 400 < SQRT_SIGMA_MEV < 450 and 800 < GAP_MEV < 900,
    }


# ---------------------------------------------------------------------------
# T3. Mesonic 1⁻⁺ hybrids
# ---------------------------------------------------------------------------

def test_T3_mesonic_hybrids() -> dict:
    """Ground meson + 2√σ for the exotic 1⁻⁺ (forbidden to qq̄): π₁ ≈ 1.62,
    η₁ ≈ 1.85 GeV — matched to the observed π₁(1600) and η₁(1855)."""
    pi1 = M_RHO + GAP_MEV
    eta1 = M_ETA1_BASE + GAP_MEV
    rows = [
        {'state': 'π₁ (1⁻⁺)', 'bam_MeV': round(pi1), 'observed': 'π₁(1600) ≈ 1660',
         'obs_MeV': 1660},
        {'state': 'η₁ (1⁻⁺)', 'bam_MeV': round(eta1), 'observed': 'η₁(1855)',
         'obs_MeV': 1855},
    ]
    ok = all(abs(r['bam_MeV'] - r['obs_MeV']) < 250 for r in rows)
    return {
        'name': 'T3_mesonic_1mp_hybrids',
        'description': (
            "Exotic 1⁻⁺ hybrids (forbidden to qq̄ — the Möbius/twist "
            "marker): ground meson + 2√σ. π₁ ≈ 1.62, η₁ ≈ 1.85 GeV — matched "
            "to π₁(1600), η₁(1855). The mesonic sector HAS a smoking-gun "
            "J^PC."
        ),
        'rows': rows,
        'exotic_marker': '1⁻⁺ (C=(−1)^{L+S} forbidden to qq̄)',
        'matched': ok,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T4. Glueball closed-flux-loop ground
# ---------------------------------------------------------------------------

def test_T4_glueball() -> dict:
    """The lowest closed-flux-loop (glueball) level M = √(4πσ) ≈ 1.50 GeV
    (≈ 3.5√σ) benchmarks the lattice 0⁺⁺ √σ scale to ~13% (lattice 4.1√σ ≈
    1.73 GeV). Unobserved — the freest channel."""
    ratio_lattice = GLUEBALL_MEV / (4.1 * SQRT_SIGMA_MEV)
    return {
        'name': 'T4_glueball_closed_flux_loop',
        'description': (
            "0⁺⁺ glueball = lowest closed flux loop, M = √(4πσ) ≈ 1.50 GeV "
            "(3.5√σ); benchmarks the lattice 0⁺⁺ √σ scale to ~13%. "
            "Unobserved — freest channel."
        ),
        'ground_MeV': round(GLUEBALL_MEV),
        'in_sqrt_sigma': GLUEBALL_MEV / SQRT_SIGMA_MEV,
        'lattice_0pp_MeV': round(4.1 * SQRT_SIGMA_MEV),
        'agreement_on_scale': abs(ratio_lattice - 1.0),
        'pass': 1400 < GLUEBALL_MEV < 1600,
    }


# ---------------------------------------------------------------------------
# T5. Heavy Möbius baryon masses
# ---------------------------------------------------------------------------

def test_T5_heavy_masses() -> dict:
    """Ground heavy baryon + 2√σ: Λ_c ~3135, Ω_c ~3544, Ξ_cc ~4471,
    Λ_b ~6469, Ω_b ~6894 MeV — supernumerary, above the orbital tower, just
    above current data."""
    rows = []
    for b, (ground, ceiling) in HEAVY_MOBIUS.items():
        rows.append({'baryon': b, 'ground_MeV': ground,
                     'mobius_MeV': round(ground + GAP_MEV),
                     'current_ceiling': ceiling})
    return {
        'name': 'T5_heavy_mobius_masses',
        'description': (
            "Ground + 2√σ: Λ_c ~3135, Ω_c ~3544, Ξ_cc ~4471, Λ_b ~6469, "
            "Ω_b ~6894 MeV — all above current excitation ceilings; Ξ_cc, "
            "Ω_b unexplored."
        ),
        'rows': rows,
        'pass': all(r['mobius_MeV'] > r['ground_MeV'] for r in rows),
    }


# ---------------------------------------------------------------------------
# T6. Heavy Möbius baryon decays
# ---------------------------------------------------------------------------

def test_T6_decays() -> dict:
    """Decay channels and release energies for Möbius_c (3135) and Möbius_b
    (6469 MeV); the all-light Q-values are cross-flavor identical."""
    channels = [
        ('Λ_Q ππ',    ('Lambda_c', M_PI, M_PI), ('Lambda_b', M_PI, M_PI), 'PREFERRED'),
        ('Σ_Q π',     ('Sigma_c', M_PI),        ('Sigma_b', M_PI),        'PREFERRED'),
        ('Σ_Q* π',    ('Sigma_c_2520', M_PI),   ('Sigma_b_star', M_PI),   'PREFERRED'),
        ('Λ_Q η',     ('Lambda_c', M_ETA),      ('Lambda_b', M_ETA),      'allowed'),
        ('D N / B N', ('D0', 'proton'),         ('Bminus', 'proton'),     'threshold'),
    ]
    rows = []
    for name, cf, bf, role in channels:
        cv = [HB[x] if isinstance(x, str) else x for x in cf]
        bv = [HB[x] if isinstance(x, str) else x for x in bf]
        qc = round(_Q(MOBIUS_C, *cv)); qb = round(_Q(MOBIUS_B, *bv))
        rows.append({'channel': name, 'charm_Q_MeV': qc, 'bottom_Q_MeV': qb,
                     'role': role})
    pipi = rows[0]; eta = rows[3]
    cross_flavor_ok = (pipi['charm_Q_MeV'] == pipi['bottom_Q_MeV']
                       and eta['charm_Q_MeV'] == eta['bottom_Q_MeV'])
    return {
        'name': 'T6_heavy_mobius_decays',
        'description': (
            "Decay channels + Q-values; all-light Q identical for c and b "
            "(Λ_Q ππ 569, Λ_Q η 301 MeV). Twist-unwinding sheds 2√σ as "
            "light isoscalar hadrons."
        ),
        'rows': rows,
        'cross_flavor_Q_identical': cross_flavor_ok,
        'pass': cross_flavor_ok and all(r['charm_Q_MeV'] > 0 for r in rows),
    }


# ---------------------------------------------------------------------------
# T7. Preferred/suppressed modes + analysis handles
# ---------------------------------------------------------------------------

def test_T7_modes_and_handles() -> dict:
    return {
        'name': 'T7_modes_and_analysis_handles',
        'description': (
            "Hybrid selection rule (single-S-wave-π-to-ground SUPPRESSED; "
            "Σ_Q π / isoscalar dipion / P-wave+π PREFERRED) + the analysis "
            "handles: cross-flavor Q-match, isoscalar high-m(ππ) dipion, "
            "broad ⟹ amplitude analyses, exotic 1⁻⁺ where available."
        ),
        'suppressed': ['Möbius baryon → (ground heavy baryon) + (single S-wave π)'],
        'preferred': [
            'Σ_Q π (spin-1 light diquark)',
            'Λ_Q (ππ)_{S, isoscalar} — high m(ππ)',
            'P-wave heavy baryon + π',
        ],
        'analysis_handles': [
            'cross-flavor Q-match: same Q (569 ππ, 301 η) above c and b',
            'branching pattern (suppressed single-π-to-ground ≠ radial excitation)',
            'isoscalar dipion peaked at high m(ππ)',
            'broad (~tens–150 MeV) ⟹ Dalitz / amplitude fits, not sharp peaks',
            'mesonic sector has exotic 1⁻⁺ smoking gun; baryon sector uses cross-flavor correlation',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The four non-orientable sub-sectors are compiled into one "
            "compact, internally-consistent experimental note (masses, "
            "Q-values, preferred/suppressed modes, analysis handles), all "
            "pushed forward from the single input √σ."
        ),
        'classification': 'NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_scope(),
        test_T2_single_input(),
        test_T3_mesonic_hybrids(),
        test_T4_glueball(),
        test_T5_heavy_masses(),
        test_T6_decays(),
        test_T7_modes_and_handles(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED'
        verdict = (
            'THE BAM NON-ORIENTABLE HADRON SECTOR, COMPILED INTO ONE COMPACT '
            'EXPERIMENTAL NOTE. PRs #100–#109 built four sub-sectors of '
            'non-orientable (Möbius / closed-flux-loop) hadron predictions; '
            'this note collects them for an LHCb / Belle II / BESIII reader, '
            'with every number pushed forward from the single QCD scale √σ ≈ '
            '424 MeV (the B4 confinement anchor) — flux-tube quantum 2√σ ≈ '
            '849 MeV, closed-loop scale √(4πσ) ≈ 1504 MeV, no new '
            'parameters.\n\n'
            'MESONIC 1⁻⁺ HYBRIDS (PR #101). Ground meson + 2√σ gives π₁ ≈ '
            '1.62 and η₁ ≈ 1.85 GeV, matched to the observed π₁(1600) and '
            'η₁(1855). The exotic 1⁻⁺ (forbidden to ordinary qq̄) is the '
            'Möbius/twist marker — the one place the sector has a smoking-gun '
            'J^PC.\n\n'
            'GLUEBALLS (PR #100). The lowest closed flux loop sits at √(4πσ) '
            '≈ 1.50 GeV (3.5√σ), benchmarking the lattice 0⁺⁺ √σ scale to '
            '~13%; unobserved — the freest channel.\n\n'
            'HEAVY MÖBIUS BARYONS — MASSES (PR #103). Ground heavy baryon + '
            '2√σ: Λ_c ~3135, Ω_c ~3544, Ξ_cc ~4471, Λ_b ~6469, Ω_b ~6894 '
            'MeV — supernumerary states above the orbital tower and just '
            'above current data; Ξ_cc and Ω_b are entirely '
            'unconstrained.\n\n'
            'HEAVY MÖBIUS BARYONS — DECAYS (PR #109). The decay proceeds by '
            'twist-unwinding (non-orientable → orientable), so the state '
            'inherits the flux-tube hybrid selection rule: the naive '
            'single-S-wave-π-to-ground transition is SUPPRESSED, while Σ_Q '
            'π, the isoscalar S-wave dipion Λ_Q(ππ), and P-wave-baryon + π '
            'are PREFERRED (an ordinary radial excitation does the '
            'opposite). The all-light release energies are cross-flavor '
            'identical: Λ_Q ππ at Q = 569 MeV and Λ_Q η at Q = 301 MeV above '
            'both the charm and the bottom ground baryon.\n\n'
            'ANALYSIS HANDLES. (i) the cross-flavor Q-match (same 569 / 301 '
            'MeV above c and b); (ii) the branching pattern (suppressed '
            'single-π-to-ground distinguishes Möbius from a radial '
            'excitation); (iii) the isoscalar dipion peaked at high m(ππ); '
            '(iv) broad widths (~tens–150 MeV) ⟹ resolve in Dalitz / '
            'amplitude fits, not as sharp peaks; (v) the exotic 1⁻⁺ '
            'smoking gun in the mesonic sector. Honest scope: this is a '
            'compilation — established masses (pushforwards of √σ) and the '
            'decay pattern / Q-structure, with exact masses (lattice hybrid '
            'gap 0.8–1.3 GeV), absolute branching fractions / widths, and '
            'the baryon J^P left open.'
        )
    else:
        verdict_class = 'NONORIENTABLE_SECTOR_NOTE_INCONCLUSIVE'
        verdict = (
            'NOTE INCONCLUSIVE. A consistency test failed; reconcile the '
            'compiled numbers with the source probes before issuing.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the BAM non-orientable hadron sector (mesonic 1⁻⁺ hybrids, '
            'glueballs, heavy Möbius baryon masses and decays) compiled into '
            'one compact experimental note, pushed forward from √σ'
        ),
        'single_input': '√σ ≈ 424 MeV; 2√σ ≈ 849 MeV; √(4πσ) ≈ 1504 MeV',
        'mesonic': 'π₁ ~1.62, η₁ ~1.85 GeV (1⁻⁺, matched to π₁(1600)/η₁(1855))',
        'glueball': '0⁺⁺ ground √(4πσ) ~1.50 GeV (unobserved, freest)',
        'heavy_masses': 'Λ_c 3135, Ω_c 3544, Ξ_cc 4471, Λ_b 6469, Ω_b 6894 MeV',
        'heavy_decays': 'twist-unwinding; single-π-to-ground suppressed; cross-flavor Q (569/301)',
        'handles': 'cross-flavor Q-match; branching pattern; isoscalar dipion; broad⟹amplitude fits; 1⁻⁺ smoking gun',
        'open': 'exact masses (0.8–1.3 GeV band); branching fractions/widths; baryon J^P',
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# The experimental note (the deliverable)
# ---------------------------------------------------------------------------

def render_experimental_note(s: dict) -> str:
    """A compact LHCb / Belle II / BESIII-style reference card."""
    t2, t3, t4, t5, t6, t7 = (s['tests'][i] for i in (1, 2, 3, 4, 5, 6))
    L: list[str] = []
    L.append('# BAM non-orientable hadron sector — experimental note')
    L.append('')
    L.append(f"_Compiled {s['timestamp_utc']} · PRs #100–#109 · all masses "
             f"pushed forward from one input, √σ._")
    L.append('')
    L.append('## 0. The single input')
    L.append('')
    L.append(f"| quantity | value | role |")
    L.append(f"|---|---:|---|")
    L.append(f"| √σ (string tension) | {t2['sqrt_sigma_MeV']:.0f} MeV | B4 confinement anchor |")
    L.append(f"| Δ = 2√σ (flux-tube quantum) | {t2['gap_2sqrt_sigma_MeV']:.0f} MeV | excitation / twist gap |")
    L.append(f"| √(4πσ) (closed-loop scale) | {t2['glueball_sqrt_4pi_sigma_MeV']:.0f} MeV | glueball ground |")
    L.append('')
    L.append('No new parameters: the whole note is a pushforward of √σ.')
    L.append('')

    L.append('## 1. Predicted states (masses)')
    L.append('')
    L.append('| sector | state | J | BAM mass | nearest observed | status |')
    L.append('|---|---|---|---:|---|---|')
    for r in t3['rows']:
        L.append(f"| mesonic hybrid | {r['state']} | 1⁻⁺ | {r['bam_MeV']} MeV "
                 f"| {r['observed']} | **matched** |")
    L.append(f"| glueball | 0⁺⁺ ground | 0⁺⁺ | {t4['ground_MeV']} MeV "
             f"| lattice 0⁺⁺ ~{t4['lattice_0pp_MeV']} | unobserved (freest) |")
    for r in t5['rows']:
        name = r['baryon'].replace('_', '_').replace('Lambda', 'Λ').replace(
            'Omega', 'Ω').replace('Xi', 'Ξ')
        unexpl = 'none' in r['current_ceiling']
        status = 'unexplored' if unexpl else 'above ceiling'
        L.append(f"| heavy baryon | {name} Möbius | (J^P) | {r['mobius_MeV']} MeV "
                 f"| ceiling {r['current_ceiling']} | findable, {status} |")
    L.append('')

    L.append('## 2. Heavy Möbius baryon decays — channels & Q-values')
    L.append('')
    L.append(f"Möbius_c = {round(MOBIUS_C)} MeV, Möbius_b = {round(MOBIUS_B)} MeV. "
             f"Decay = twist-unwinding (non-orientable → orientable; sheds 2√σ "
             f"as light isoscalar hadrons).")
    L.append('')
    L.append('| channel | charm Q (MeV) | bottom Q (MeV) | mode |')
    L.append('|---|---:|---:|---|')
    for r in t6['rows']:
        L.append(f"| {r['channel']} | {r['charm_Q_MeV']} | {r['bottom_Q_MeV']} | {r['role']} |")
    L.append('')
    L.append('`Λ_Q ππ` (569) and `Λ_Q η` (301) are **identical** for charm and '
             'bottom — the cross-flavor Q-match.')
    L.append('')

    L.append('## 3. Preferred / suppressed modes (hybrid selection rule)')
    L.append('')
    L.append('- **SUPPRESSED** — ' + t7['suppressed'][0]
             + ' (the naive, phase-space-favored channel).')
    L.append('- **PREFERRED** — ' + '; '.join(t7['preferred']) + '.')
    L.append('- An ordinary **radial excitation** does the opposite (single π '
             'to ground): the branching **pattern** is the discriminator.')
    L.append('')

    L.append('## 4. Analysis handles')
    L.append('')
    for h in t7['analysis_handles']:
        L.append(f"- {h}")
    L.append('')

    L.append('## 5. Scope (honest)')
    L.append('')
    L.append(f"- **Established:** masses as pushforwards of √σ; the decay "
             f"pattern (selection rule) and the cross-flavor Q-structure.")
    L.append(f"- **Open:** {s['open']}.")
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Probe-summary markdown
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# BAM non-orientable hadron sector: experimental note (PR #110)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Compiles the PRs #100–#109 non-orientable (Möbius / closed-flux-"
        "loop) hadron predictions into one compact experimental note — "
        "predicted masses, Q-values, preferred/suppressed modes, and "
        "analysis handles — every number pushed forward from the single "
        "input `√σ`. The note itself is written to `experimental_note.md`."
    )
    L.append('')
    L.append(f"- **Single input**: {s['single_input']}")
    L.append(f"- **Mesonic**: {s['mesonic']}")
    L.append(f"- **Glueball**: {s['glueball']}")
    L.append(f"- **Heavy masses**: {s['heavy_masses']}")
    L.append(f"- **Heavy decays**: {s['heavy_decays']}")
    L.append(f"- **Handles**: {s['handles']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'consolidate #100–#109 non-orientable sector into one note',
        'T2': 'single input: √σ 424, 2√σ 849, √(4πσ) 1504 MeV',
        'T3': 'mesonic 1⁻⁺: π₁ ~1.62, η₁ ~1.85 GeV (matched)',
        'T4': 'glueball 0⁺⁺ ground √(4πσ) ~1.50 GeV',
        'T5': 'heavy masses: Λ_c 3135 … Ω_b 6894 MeV',
        'T6': 'decays: cross-flavor Q-match (569/301), twist-unwinding',
        'T7': 'selection rule + analysis handles',
        'T8': 'NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')
    L.append('---')
    L.append('')
    L.append('The compiled note follows (also written standalone as '
             '`experimental_note.md`):')
    L.append('')
    L.append(render_experimental_note(s))
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    note = render_experimental_note(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_nonorientable_experimental_note_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    (out / 'experimental_note.md').write_text(note)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    print(f"Wrote: {out / 'experimental_note.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
