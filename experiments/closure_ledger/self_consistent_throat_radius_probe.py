"""
Self-consistent throat radius / finite self-energy probe.

Targets the remaining BAM scale anchor: the throat radius R_MID (≈ the
invariant bulk separation ΔR), currently imposed in constants.py. The
THESIS open problem asks for R_MID determined dynamically as the
equilibrium throat radius, with finite self-energy.

This probe builds that equilibrium and reports — honestly, against the
B4 scale-modulus theorem (PR #52) — what it does and does not fix.

The B4 constraint: the closure-ledger machinery is scale-free (fixes
ratios, not an absolute length). A self-energy built from scale-free
terms (all ∝ 1/R) has NO stationary point in R. To pin an absolute R*,
the energy must contain two terms of different R-scaling, one carrying a
dimensionful coupling. The probe makes this explicit.

Finite self-energy (BAM regularization): a point charge has divergent
U = e²/(8πε₀ r) → ∞. The throat removes the divergence geometrically:
  - short distance: the wormhole throat is the inner boundary (no
    r < R_MID), capping the field → U_EM = α ℏc/(2 R_MID) finite;
  - long distance: on compact S³ the Green function G(ψ) is regular at
    the antipode with a zero-mean background → no far-field divergence.
So U_EM/(m c²) = α/2 ≈ 0.0036 — finite, no UV divergence.

Equilibrium: E(R) = A/R + B R² (EM repulsion ∝ 1/R + cohesion ∝ R²) →
unique stable minimum R* = (A/2B)^(1/3). The throat radius is a
stationary point, not imposed — but R* rides on A/B (rescaling B → B/λ³
sends R* → λ R*), the one scale modulus, consistent with B4.

Pure-EM relocation: demanding m c² = U_EM is R-independent (both ∝1/R)
→ fixes g = 2/α, relating geometry to the fine-structure constant
rather than fixing R — the B4 obstruction in self-energy language.

Tests:
  T1. Finite self-energy (throat caps the field) vs divergent point charge.
  T2. S³ far-field regular (antipode finite; zero-mean background).
  T3. Scale-free obstruction (EM-only 1/R has no equilibrium).
  T4. Self-consistent equilibrium E=A/R+BR² → stable R*=(A/2B)^(1/3).
  T5. Scale modulus: B→B/λ³ ⟹ R*→λR* (absolute value = one anchor).
  T6. Finite, renormalization-free: U_EM/(m c²) = α/2.
  T7. Pure-EM relocation: m c²=U_EM is R-independent ⟹ g=2/α.
  T8. Assessment: equilibrium recasts the anchor; value still one input.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.transaction.s3_geometry import s3_green_potential
from geometrodynamics.constants import R_MID


PI = math.pi

# Physical constants (SI)
HBAR = 1.054571817e-34       # J·s
M_E = 9.1093837015e-31       # kg
C_LIGHT = 2.99792458e8       # m/s
ALPHA = 7.2973525693e-3      # fine-structure constant (dimensionless)

LAMBDA_C = HBAR / (M_E * C_LIGHT)     # reduced Compton wavelength, m
M_E_C2 = M_E * C_LIGHT ** 2           # rest energy, J


# ---------------------------------------------------------------------------
# Self-energy helpers
# ---------------------------------------------------------------------------

def coulomb_self_energy(cutoff_radius: float) -> float:
    """Coulomb self-energy of a charge with field capped below
    `cutoff_radius`: U = e²/(8πε₀ R) = α ℏc/(2R). Diverges as R → 0."""
    return ALPHA * HBAR * C_LIGHT / (2.0 * cutoff_radius)


def throat_energy(R: float, A: float, B: float) -> float:
    """Throat energy functional: EM repulsion (∝1/R) + cohesion (∝R²)."""
    return A / R + B * R ** 2


# ---------------------------------------------------------------------------
# T1. Finite self-energy (short distance) — throat caps the field
# ---------------------------------------------------------------------------

def test_T1_finite_self_energy() -> dict:
    """A point charge has divergent self-energy as the cutoff → 0. The
    BAM throat is the inner boundary (no r < R_MID), so the field is
    capped at the throat and U_EM = α ℏc/(2 R_MID) is finite — the
    cutoff is the throat radius, not an arbitrary ε."""
    # point charge: self-energy blows up as cutoff shrinks
    cutoffs = [1e-13, 1e-15, 1e-17, 1e-19]
    point_energies = [coulomb_self_energy(c) for c in cutoffs]
    diverges = point_energies[-1] > 100.0 * point_energies[0]
    # throat: capped at R_MID (= λ_C in proper units)
    R_throat = LAMBDA_C
    U_throat = coulomb_self_energy(R_throat)
    finite = math.isfinite(U_throat) and U_throat > 0
    return {
        'name': 'T1_finite_self_energy_throat_cutoff',
        'description': (
            "A point charge's Coulomb self-energy U = α ℏc/(2R) diverges "
            "as the cutoff R → 0. The BAM throat is the inner boundary "
            "(f(r)=1−(rs/r)²=0 at r=R_MID; no r<R_MID), capping the field "
            "at the throat → U_EM = α ℏc/(2 R_MID) finite. The cutoff is "
            "the throat radius itself, not an arbitrary ε."
        ),
        'point_charge_cutoffs_m': cutoffs,
        'point_charge_self_energies_J': point_energies,
        'point_charge_diverges': diverges,
        'throat_radius_m': R_throat,
        'throat_self_energy_J': U_throat,
        'throat_self_energy_finite': finite,
        'pass': diverges and finite,
    }


# ---------------------------------------------------------------------------
# T2. S³ far-field regular (long distance)
# ---------------------------------------------------------------------------

def test_T2_s3_farfield_regular() -> dict:
    """On compact S³ the Green function G(ψ)=((π−ψ)cot ψ − ½)/(4π²R) is
    regular at the antipode (G(π) finite) with a zero-mean neutralizing
    background (the −½) — no long-distance divergence (unlike a flat-space
    image-charge construction)."""
    R = 1.0
    # near the antipode ψ → π: G stays finite
    psis = [PI - 1e-2, PI - 1e-3, PI - 1e-4]
    G_near_antipode = [s3_green_potential(p, radius=R, eps=1e-9) for p in psis]
    antipode_finite = all(math.isfinite(g) for g in G_near_antipode)
    # the values converge (regular), not blow up
    bounded = max(abs(g) for g in G_near_antipode) < 1.0
    # analytic antipode limit: (π−ψ)cot ψ → −1, so G(π) = (−1 − ½)/(4π²R)
    G_antipode_analytic = (-1.0 - 0.5) / (4.0 * PI ** 2 * R)
    matches = abs(G_near_antipode[-1] - G_antipode_analytic) < 1e-2
    return {
        'name': 'T2_s3_farfield_regular',
        'description': (
            "On compact S³ the Green function G(ψ) is regular at the "
            "antipode (G(π) = (−1−½)/(4π²R), finite) with a zero-mean "
            "background (the −½ neutralizing charge). No far-field "
            "divergence — the long-distance end of the self-energy is "
            "regularized by S³ compactness."
        ),
        'psis_near_antipode': psis,
        'G_near_antipode': G_near_antipode,
        'G_antipode_analytic': G_antipode_analytic,
        'antipode_finite': antipode_finite,
        'bounded': bounded,
        'matches_analytic': matches,
        'pass': antipode_finite and bounded and matches,
    }


# ---------------------------------------------------------------------------
# T3. Scale-free obstruction: EM-only 1/R has no equilibrium
# ---------------------------------------------------------------------------

def test_T3_scale_free_obstruction() -> dict:
    """With only the EM self-energy E(R)=A/R (scale-free, ∝1/R), the
    energy is monotone decreasing — no stationary point, no equilibrium.
    A scale-free 1/R energy cannot fix R (the B4 obstruction in
    self-energy form)."""
    A = ALPHA * HBAR * C_LIGHT / 2.0
    Rs = np.array([0.5, 1.0, 2.0, 4.0, 8.0]) * LAMBDA_C
    E = [A / R for R in Rs]
    dE = [(-A / R ** 2) for R in Rs]
    monotone_decreasing = all(d < 0 for d in dE)
    has_stationary = any(abs(d) < 1e-60 for d in dE)
    return {
        'name': 'T3_scale_free_obstruction',
        'description': (
            "EM-only energy E(R)=A/R is monotone decreasing (dE/dR=−A/R²<0 "
            "everywhere): no stationary point, no equilibrium. A scale-free "
            "1/R energy cannot fix R — the B4 scale-modulus obstruction in "
            "self-energy language. A second term of different R-scaling is "
            "required."
        ),
        'A_coupling_J_m': A,
        'radii_m': [float(r) for r in Rs],
        'energies_J': [float(e) for e in E],
        'dE_dR': [float(d) for d in dE],
        'monotone_no_equilibrium': monotone_decreasing and not has_stationary,
        'pass': monotone_decreasing and not has_stationary,
    }


# ---------------------------------------------------------------------------
# T4. Self-consistent equilibrium
# ---------------------------------------------------------------------------

def test_T4_self_consistent_equilibrium() -> dict:
    """Adding a cohesive term B·R² (∝R², wants to shrink) to the EM
    repulsion gives E(R)=A/R+B R² with a unique stable minimum
    R*=(A/2B)^(1/3). Verify the stationary point numerically (analytic
    match + d²E/dR²>0) and that a scan of E(R) has its minimum at R*."""
    A = ALPHA * HBAR * C_LIGHT / 2.0
    # choose B so that R* lands at λ_C (illustrative; B is the anchor)
    R_target = LAMBDA_C
    B = A / (2.0 * R_target ** 3)            # from R*³ = A/2B
    R_star_analytic = (A / (2.0 * B)) ** (1.0 / 3.0)
    # numeric minimization via dense scan
    Rs = np.linspace(0.2 * R_target, 5.0 * R_target, 200001)
    E = A / Rs + B * Rs ** 2
    R_star_numeric = float(Rs[int(np.argmin(E))])
    # second derivative at R*: d²E/dR² = 2A/R³ + 2B > 0
    d2E = 2.0 * A / R_star_analytic ** 3 + 2.0 * B
    stable = d2E > 0
    rel_err = abs(R_star_numeric - R_star_analytic) / R_star_analytic
    return {
        'name': 'T4_self_consistent_equilibrium',
        'description': (
            "E(R)=A/R+B R² (EM repulsion ∝1/R + cohesion ∝R²) has a "
            "unique stable minimum R*=(A/2B)^(1/3) (d²E/dR²=2A/R³+2B>0). "
            "The throat radius is a self-consistent stationary point, not "
            "an imposed constant."
        ),
        'A_coupling': A,
        'B_coupling': B,
        'R_star_analytic_m': R_star_analytic,
        'R_star_numeric_m': R_star_numeric,
        'relative_error': rel_err,
        'second_derivative_positive': stable,
        'equilibrium_exists_and_stable': stable and rel_err < 1e-3,
        'pass': stable and rel_err < 1e-3,
    }


# ---------------------------------------------------------------------------
# T5. Scale modulus — absolute R* rides on one coupling
# ---------------------------------------------------------------------------

def test_T5_scale_modulus() -> dict:
    """R*=(A/2B)^(1/3): rescaling the cohesive coupling B → B/λ³ sends
    R* → λ R*. The absolute equilibrium radius rides on the single
    dimensionful coupling B (or A/B) — consistent with the B4 theorem
    (the self-consistency fixes R* given the coupling, not in absolute
    terms)."""
    A = ALPHA * HBAR * C_LIGHT / 2.0
    B0 = A / (2.0 * LAMBDA_C ** 3)
    rows = []
    ok = True
    R_star_0 = (A / (2.0 * B0)) ** (1.0 / 3.0)
    for lam in [1.0, 2.0, 0.5, 10.0]:
        B = B0 / lam ** 3
        R_star = (A / (2.0 * B)) ** (1.0 / 3.0)
        predicted = lam * R_star_0
        dev = abs(R_star - predicted) / predicted
        ok = ok and dev < 1e-12
        rows.append({
            'lambda': lam,
            'B_coupling': B,
            'R_star_m': R_star,
            'predicted_lambda_R0': predicted,
            'deviation': dev,
        })
    return {
        'name': 'T5_scale_modulus',
        'description': (
            "R*=(A/2B)^(1/3): rescaling B → B/λ³ gives R* → λ R*. The "
            "absolute equilibrium radius is set by the single dimensionful "
            "coupling B (the cohesive stiffness) — one scale modulus, as "
            "the B4 theorem requires. The self-consistency fixes R* given "
            "the coupling; it does not produce an absolute length from "
            "scale-free ingredients."
        ),
        'rows': rows,
        'R_star_scales_linearly_with_lambda': ok,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T6. Finite, renormalization-free: U_EM/(m c²) = α/2
# ---------------------------------------------------------------------------

def test_T6_renormalization_free() -> dict:
    """At the throat (R_MID = λ_C in proper units), the EM self-energy is
    a finite fraction of the rest energy: U_EM/(m c²) = α/2 ≈ 0.0036 — a
    small finite correction with no UV divergence, in contrast to QED's
    divergent electron self-energy (which needs renormalization)."""
    R_throat = LAMBDA_C
    U_EM = coulomb_self_energy(R_throat)
    fraction = U_EM / M_E_C2
    expected = ALPHA / 2.0
    rel_err = abs(fraction - expected) / expected
    return {
        'name': 'T6_renormalization_free',
        'description': (
            "U_EM/(m c²) = [α ℏc/(2R)]/[ℏc/R] = α/2 ≈ 0.0036 at the "
            "throat (R_MID = λ_C): a finite, small EM self-energy "
            "correction — no UV divergence, no renormalization, unlike "
            "QED's divergent point-electron self-energy."
        ),
        'throat_self_energy_J': U_EM,
        'rest_energy_J': M_E_C2,
        'fraction_U_over_mc2': fraction,
        'expected_alpha_over_2': expected,
        'relative_error': rel_err,
        'finite_and_small': fraction < 0.01,
        'pass': rel_err < 1e-6 and fraction < 0.01,
    }


# ---------------------------------------------------------------------------
# T7. Pure-EM relocation: m c² = U_EM is R-independent → g = 2/α
# ---------------------------------------------------------------------------

def test_T7_pure_em_relocation() -> dict:
    """The BAM-native balance "rest energy = EM self-energy" with a
    geometric factor g, m c² = g·α ℏc/(2R), and the bridge m c² = ℏc/R,
    gives 1 = g α/2 ⟹ g = 2/α — both sides scale as 1/R, so the balance
    is R-INDEPENDENT. It fixes a geometry–α relation, not the length. The
    self-consistency relocates the scale question to α (still one input);
    it does not fix R."""
    # mc² = ℏc/R  (bridge);  U_EM = g·α ℏc/(2R)
    # equate at two radii to show R-independence
    g_required = 2.0 / ALPHA
    rows = []
    r_independent = True
    base_g = None
    for R in [0.5 * LAMBDA_C, LAMBDA_C, 2.0 * LAMBDA_C]:
        mc2 = HBAR * C_LIGHT / R
        # solve mc2 = g·α ℏc/(2R) for g
        g = mc2 / (ALPHA * HBAR * C_LIGHT / (2.0 * R))
        if base_g is None:
            base_g = g
        r_independent = r_independent and abs(g - base_g) < 1e-9
        rows.append({'R_m': R, 'mc2_J': mc2, 'g_solved': g})
    g_matches = abs(base_g - g_required) < 1e-9
    return {
        'name': 'T7_pure_em_relocation',
        'description': (
            "Demanding rest energy = EM self-energy, m c² = g·α ℏc/(2R), "
            "with the bridge m c² = ℏc/R, gives g = 2/α — and since both "
            "sides ∝ 1/R the balance is R-INDEPENDENT. The self-consistency "
            "fixes a geometry–α relation (g = 2/α ≈ 274), not the length. "
            "The scale question is relocated to α (still one external "
            "input), not eliminated — the B4 obstruction in self-energy "
            "form."
        ),
        'g_required_2_over_alpha': g_required,
        'rows': rows,
        'g_is_R_independent': r_independent,
        'g_equals_2_over_alpha': g_matches,
        'pass': r_independent and g_matches,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    """The throat radius is realized as a finite-self-energy stable
    equilibrium (not imposed); the self-energy is finite (throat caps the
    field; no UV divergence). The absolute value still rides on one
    dimensionful coupling (B4); the self-consistency recasts the anchor as
    an equilibrium condition and relates it to α, without deriving the
    value."""
    return {
        'name': 'T8_assessment',
        'description': (
            "The remaining BAM scale anchor (R_MID ≈ ΔR) is recast: from "
            "imposed → invariant geometric length ΔR (#53) → "
            "finite-self-energy stable equilibrium R*=(A/2B)^(1/3) (this "
            "probe). The self-energy is finite (throat caps the field; no "
            "UV divergence; U_EM/mc²=α/2). The absolute value still rides "
            "on one dimensionful coupling (the cohesive stiffness B, "
            "equivalently a relation to α), as the B4 theorem requires. "
            "Progress: an explicit stable equilibrium + finite self-energy "
            "+ the α-relocation; not a derivation of the absolute scale."
        ),
        'relocation_chain': [
            'imposed R_MID (constants.py)',
            'ΔR invariant geometric length (#53)',
            'finite-self-energy stable equilibrium R*=(A/2B)^(1/3) (this probe)',
        ],
        'self_energy_finite': True,
        'equilibrium_stable': True,
        'absolute_value_derived': False,
        'remaining_prize': 'pin B (or α-relation) to a second fixed scale (e.g. Planck via closure quantum)',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_finite_self_energy()
    t2 = test_T2_s3_farfield_regular()
    t3 = test_T3_scale_free_obstruction()
    t4 = test_T4_self_consistent_equilibrium()
    t5 = test_T5_scale_modulus()
    t6 = test_T6_renormalization_free()
    t7 = test_T7_pure_em_relocation()
    t8 = test_T8_assessment()
    tests = [t1, t2, t3, t4, t5, t6, t7, t8]

    core = [t1, t2, t3, t4, t5, t6, t7]
    if all(t['pass'] for t in core):
        verdict_class = 'SELF_CONSISTENT_THROAT_EQUILIBRIUM'
        verdict = (
            'SELF-CONSISTENT THROAT EQUILIBRIUM. The remaining BAM scale '
            'anchor — the throat radius R_MID (≈ the invariant bulk '
            'separation ΔR) — is recast from an imposed constant into a '
            'finite-self-energy stable equilibrium.\n\n'
            'FINITE SELF-ENERGY. A point charge has divergent Coulomb '
            'self-energy U = α ℏc/(2R) → ∞ as R → 0. The BAM throat '
            'removes the divergence geometrically: (short distance) the '
            'wormhole throat is the inner boundary (no r < R_MID), capping '
            'the field at the throat → U_EM = α ℏc/(2 R_MID) finite; (long '
            'distance) on compact S³ the Green function G(ψ) is regular at '
            'the antipode with a zero-mean background. The result is '
            'U_EM/(m c²) = α/2 ≈ 0.0036 — a finite, small correction with '
            'NO UV divergence (no renormalization), unlike QED.\n\n'
            'SELF-CONSISTENT EQUILIBRIUM. With only the scale-free EM '
            'energy (∝1/R) there is no stationary point — a scale-free 1/R '
            'energy cannot fix R (the B4 obstruction in self-energy form). '
            'Adding a cohesive term (∝R²) gives E(R)=A/R+B R² with a '
            'unique stable minimum R*=(A/2B)^(1/3) (d²E/dR²>0): the throat '
            'radius is a stationary point, not imposed.\n\n'
            'HONEST CAVEAT (B4). The absolute R* rides on the single '
            'dimensionful coupling B: rescaling B → B/λ³ sends R* → λ R* '
            '(one scale modulus). Equivalently, the BAM-native balance '
            '"rest energy = EM self-energy" is R-independent (both ∝1/R) '
            'and fixes a geometry–α relation g = 2/α ≈ 274, not the length '
            '— relocating the scale question to α. So the self-consistency '
            'recasts the anchor as a finite-self-energy equilibrium '
            'condition (genuine progress) and relates it to α, but it does '
            'NOT derive the absolute value — by the scale-modulus theorem, '
            'that needs one external dimensionful coupling.\n\n'
            'Relocation chain: imposed R_MID → ΔR invariant geometric '
            'length (#53) → finite-self-energy stable equilibrium (this '
            'probe). The remaining prize: pin B (or the α-relation) to a '
            'second fixed scale, e.g. a closure-quantum relation to the '
            'Planck length.'
        )
    else:
        verdict_class = 'EQUILIBRIUM_FAILS'
        verdict = (
            'EQUILIBRIUM FAILS. No stable stationary point exists, or the '
            'self-energy is not finite. Investigate the failing test.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'target': 'remaining BAM scale anchor: throat radius R_MID (≈ ΔR)',
        'finite_self_energy': 'U_EM = α ℏc/(2 R_MID); throat caps the field; U/mc² = α/2',
        'equilibrium': 'E(R)=A/R+B R² → stable R*=(A/2B)^(1/3)',
        'b4_caveat': 'absolute R* rides on one dimensionful coupling (or a relation to α)',
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
    L.append('# Self-consistent throat radius / finite self-energy probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Targets the remaining BAM scale anchor — the throat radius '
        'R_MID (≈ the invariant bulk separation ΔR) — recasting it from '
        'an imposed constant into a finite-self-energy stable equilibrium, '
        'honestly against the B4 scale-modulus theorem (PR #52).'
    )
    L.append('')
    L.append(f"- **Target**: {s['target']}")
    L.append(f"- **Finite self-energy**: {s['finite_self_energy']}")
    L.append(f"- **Equilibrium**: `{s['equilibrium']}`")
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
            value = "throat caps field → U_EM finite; point charge diverges"
        elif nm.startswith('T2'):
            value = "S³ Green regular at antipode (zero-mean background)"
        elif nm.startswith('T3'):
            value = "EM-only 1/R monotone → no equilibrium (B4)"
        elif nm.startswith('T4'):
            value = f"E=A/R+BR² → stable R*=(A/2B)^⅓ (err {t['relative_error']:.0e})"
        elif nm.startswith('T5'):
            value = "B→B/λ³ ⟹ R*→λR* (one modulus)"
        elif nm.startswith('T6'):
            value = f"U_EM/mc² = α/2 = {t['fraction_U_over_mc2']:.4f} (finite, no UV div)"
        elif nm.startswith('T7'):
            value = f"mc²=U_EM R-independent ⟹ g=2/α={t['g_required_2_over_alpha']:.1f}"
        elif nm.startswith('T8'):
            value = "equilibrium recasts anchor; value still one input"
        else:
            value = '—'
        L.append(f"| {nm[:2]} | `{nm}` | {value} | {passed} |")
    L.append('')

    # T1
    t1 = s['tests'][0]
    L.append('## T1: Finite self-energy — throat caps the field')
    L.append('')
    L.append('| point-charge cutoff (m) | self-energy (J) |')
    L.append('|---:|---:|')
    for c, e in zip(t1['point_charge_cutoffs_m'], t1['point_charge_self_energies_J']):
        L.append(f"| {c:.0e} | {e:.3e} |")
    L.append('')
    L.append(f"Throat (capped at R_MID = {t1['throat_radius_m']:.3e} m): "
             f"U_EM = {t1['throat_self_energy_J']:.3e} J (finite).")
    L.append('')

    # T2
    t2 = s['tests'][1]
    L.append('## T2: S³ far-field regular')
    L.append('')
    L.append(f"- G near antipode: {[f'{g:.5f}' for g in t2['G_near_antipode']]}")
    L.append(f"- analytic G(π) = (−1−½)/(4π²R) = {t2['G_antipode_analytic']:.5f} "
             f"(finite; zero-mean −½ background)")
    L.append('')

    # T3
    t3 = s['tests'][2]
    L.append('## T3: Scale-free obstruction (EM-only 1/R)')
    L.append('')
    L.append('| R/λ_C | E=A/R (J) | dE/dR |')
    L.append('|---:|---:|---:|')
    for R, e, d in zip(t3['radii_m'], t3['energies_J'], t3['dE_dR']):
        L.append(f"| {R/LAMBDA_C:.1f} | {e:.3e} | {d:.3e} |")
    L.append('')
    L.append('dE/dR < 0 everywhere → no stationary point → no equilibrium '
             'from scale-free 1/R alone (the B4 obstruction).')
    L.append('')

    # T4
    t4 = s['tests'][3]
    L.append('## T4: Self-consistent equilibrium')
    L.append('')
    L.append(f"- E(R) = A/R + B·R², A = {t4['A_coupling']:.3e}, B = {t4['B_coupling']:.3e}")
    L.append(f"- R* analytic = (A/2B)^⅓ = {t4['R_star_analytic_m']:.4e} m")
    L.append(f"- R* numeric (scan minimum) = {t4['R_star_numeric_m']:.4e} m "
             f"(rel err {t4['relative_error']:.1e})")
    L.append(f"- d²E/dR² > 0 (stable minimum): {t4['second_derivative_positive']}")
    L.append('')

    # T5
    t5 = s['tests'][4]
    L.append('## T5: Scale modulus — absolute R* rides on one coupling')
    L.append('')
    L.append('| λ | B coupling | R* (m) | predicted λ·R₀ | deviation |')
    L.append('|---:|---:|---:|---:|---:|')
    for r in t5['rows']:
        L.append(
            f"| {r['lambda']:.1f} | {r['B_coupling']:.3e} | {r['R_star_m']:.4e} | "
            f"{r['predicted_lambda_R0']:.4e} | {r['deviation']:.1e} |"
        )
    L.append('')
    L.append('Rescaling the cohesive coupling rescales R* linearly — the one '
             'scale modulus (B4-consistent).')
    L.append('')

    # T6
    t6 = s['tests'][5]
    L.append('## T6: Finite, renormalization-free')
    L.append('')
    L.append(f"- U_EM (throat) = {t6['throat_self_energy_J']:.3e} J")
    L.append(f"- rest energy m c² = {t6['rest_energy_J']:.3e} J")
    L.append(f"- U_EM/(m c²) = {t6['fraction_U_over_mc2']:.6f} = α/2 = "
             f"{t6['expected_alpha_over_2']:.6f} (finite, small; no UV divergence)")
    L.append('')

    # T7
    t7 = s['tests'][6]
    L.append('## T7: Pure-EM relocation → g = 2/α')
    L.append('')
    L.append('| R (m) | m c² (J) | g solved |')
    L.append('|---:|---:|---:|')
    for r in t7['rows']:
        L.append(f"| {r['R_m']:.3e} | {r['mc2_J']:.3e} | {r['g_solved']:.4f} |")
    L.append('')
    L.append(f"g is R-independent: {t7['g_is_R_independent']}; g = 2/α = "
             f"{t7['g_required_2_over_alpha']:.4f}. The balance fixes a "
             f"geometry–α relation, not the length — the scale relocates to α.")
    L.append('')

    # T8
    t8 = s['tests'][7]
    L.append('## T8: Assessment')
    L.append('')
    L.append('Relocation chain:')
    for step in t8['relocation_chain']:
        L.append(f"  - {step}")
    L.append('')
    L.append(f"- self-energy finite: {t8['self_energy_finite']}; "
             f"equilibrium stable: {t8['equilibrium_stable']}")
    L.append(f"- absolute value derived: {t8['absolute_value_derived']}")
    L.append(f"- remaining prize: {t8['remaining_prize']}")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **The absolute value.** Still one dimensionful coupling (B, or '
             'a relation to α). A genuine derivation would pin it to a second '
             'fixed scale (e.g. a closure-quantum relation to the Planck length).')
    L.append('- **The cohesive term from first principles.** B·R² is a '
             'Poincaré-stress-like cohesion; deriving it from the BAM throat '
             'action is the follow-on.')
    L.append('- **Pair-production threshold.** Identified structurally '
             '(2 m_e c² = nucleating a throat–antithroat pair at the lowest '
             'stable R*); a dynamical nucleation calculation is future work.')
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
    out = here / 'runs' / f'{ts}_self_consistent_throat_radius_probe'
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
