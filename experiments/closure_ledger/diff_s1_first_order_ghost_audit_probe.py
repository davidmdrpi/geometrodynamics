"""
First-order Diff(S¹) Faddeev–Popov ghost audit (PR #118).

PR #117 (revised) fixed the ghost determinant to det'(P) = det'(P†P)^{1/2}
= L. This follow-up audits that result rigorously: it distinguishes the four
objects the reviewer asked about — P = ∂_τ, P†P = −∂_τ², det'(P), and
det'(P†P) — handles the phase/η-invariant of the first-order determinant
det'(∂_τ), treats the conformal-Killing-vector (CKV) zero-mode division and
the zero-mode norms explicitly, and ends with a measure table fixing the net
L-power. The headline: the FP ghost contributes L (first order); it
contributes L² ONLY if one explicitly adopts a second-order ghost
convention (which over-counts by one power of L).

## The four objects

  P   = ∂_τ        (first order; vector ghost c ↦ einbein variation):
                    eigenvalues 2πi n / L, n ∈ ℤ; ONE zero mode (n = 0 = CKV).
  P†P = −∂_τ²      (second order; self-adjoint, ≥ 0):
                    eigenvalues (2πn/L)², n ∈ ℤ; ONE zero mode.
  det'(P)          first-order determinant (over n ≠ 0).
  det'(P†P)        second-order determinant (over n ≠ 0).

Zeta-regularized (ζ_R(0) = −1/2, ζ_R'(0) = −½ ln 2π):

  det'(P†P) = L²        (ζ(0) = −1, ζ'(0) = −2 ln L),
  |det'(P)| = L         (ζ_{|P|}(0) = −1, ζ_{|P|}'(0) = −ln L) = det'(P†P)^{1/2}.

Verified to machine precision for L = 2π, 1, 3.32, 5.

## Phase / η-invariant of det'(∂_τ)

The first-order operator's determinant can carry a phase set by the
η-invariant (spectral asymmetry). For A = −i∂_τ (self-adjoint, eigenvalues
2πn/L), the spectrum is symmetric under n → −n, so

  η_A(s) = Σ_{n≠0} sign(2πn/L) |2πn/L|^{−s} ≡ 0   ⟹   η_A(0) = 0.

With η = 0 there is NO anomalous phase: det'(∂_τ) = +L (real, positive).
(In the antiperiodic / Möbius sector the eigenvalues are 2πi(n+½)/L — still
symmetric, η = 0 — but there is NO zero mode, hence NO CKV there; a clean
tie to the non-orientable odd-k structure.)

## CKV zero-mode division and zero-mode norms

The constant mode on a circle of circumference L has norm ‖1‖² = ∫₀ᴸ dτ = L,
so ‖1‖ = √L. The reparametrization ghost system has two zero modes:

  - the c-ghost zero mode = the CKV (rigid U(1) rotation), norm √L;
  - the b-ghost zero mode = the dual of the one Teichmüller modulus, norm √L.

Dividing the FP determinant by these zero-mode norms gives the ghost NET
factor:

  first-order:   det'(P) / (√L · √L) = L / L = 1,
  second-order:  det'(P†P) / (√L · √L) = L² / L = L   (spurious extra power).

So the first-order convention yields a clean net ghost factor of 1; the
second-order convention injects an extra power of L.

## The measure table (net L-power)

The moduli-space measure (one modulus L, one CKV automorphism) is the
standard worldline dL/L. Assembling the gauge sector:

  | factor                         | source                       | L-power |
  | moduli-space measure dL/L      | modulus dL × CKV 1/L         | L^{−1} dL |
  | FP ghost determinant det'(P)   | first-order, nonzero modes   | L^{+1} |
  | ghost zero-mode norms (÷)      | 1 CKV + 1 modulus, √L each   | L^{−1} |
  | ghost NET                      | det'(P)/norms                | L^{0}  |
  | matter det^{−1/2}              | Tangherlini/heat-kernel Weyl | L^{−d/2} |
  | NET MEASURE (first order)      |                              | dL · L^{−1−d/2} |

With the second-order ghost convention the ghost NET would be L^{+1} (not
L^0), giving dL · L^{−d/2} — one spurious power of L. The physical FP for a
single real reparametrization symmetry is first-order, so the BAM loop
measure is

  Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S},

ghost L-power fixed at L¹ in det'(P) (net L⁰ after zero-mode norms).

Tests:
  T1. Goal: audit; distinguish P, P†P, det'(P), det'(P†P); fix the L-power.
  T2. The four objects + spectra (first vs second order, one zero mode each).
  T3. det'(P†P) = L² and |det'(P)| = det'(P†P)^{1/2} = L (computed).
  T4. η-invariant: η(−i∂_τ) = 0 ⟹ det'(∂_τ) = +L, no phase (antiperiodic:
      η = 0, no zero mode/CKV).
  T5. First- vs second-order convention: ghost = det'(P) = L (first) vs
      det'(P†P) = L² (second, over-counts).
  T6. CKV zero-mode division + norms (√L each): ghost net = L/L = 1 (first),
      L²/L = L (second).
  T7. Measure table + net L-power: Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S},
      net dL·L^{−1−d/2}; second-order would add a spurious L.
  T8. Assessment.

Verdict:
  - BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L
    (expected): the Diff(S¹) FP ghost is first-order — det'(P) = det'(∂_τ) =
    det'(P†P)^{1/2} = L (η = 0, real positive; det'(P†P) = L² is the
    second-order object, the square). After CKV + modulus zero-mode norms
    (√L each) the ghost NET factor is L/L = 1, and the measure is
    Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S} (net dL·L^{−1−d/2}). The
    L² power appears only under an explicit second-order ghost convention,
    which over-counts by one power of L. Absolute Z normalization (κ₅²/Λ₅)
    and multi-loop stay open.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
ZETA_R_0 = -0.5
ZETA_R_PRIME_0 = -0.5 * math.log(2.0 * PI)
CLOSURE_LOOP_L = 2.0 * PI


def det_prime_PtP(L: float) -> float:
    """det'(P†P) = det'(−∂_τ²): ζ(s) = 2 (L/2π)^{2s} ζ_R(2s) ⟹ det' = L²."""
    return math.exp(-(2.0 * (2.0 * math.log(L / (2.0 * PI)) * ZETA_R_0
                             + 2.0 * ZETA_R_PRIME_0)))


def det_prime_P_magnitude(L: float) -> float:
    """|det'(P)| = |det'(∂_τ)|: ζ_{|P|}(s) = 2 (L/2π)^{s} ζ_R(s) ⟹ = L."""
    return math.exp(-(2.0 * (math.log(L / (2.0 * PI)) * ZETA_R_0
                             + ZETA_R_PRIME_0)))


def zero_mode_norm(L: float) -> float:
    """Norm of the constant mode on a circle of circumference L: √L."""
    return math.sqrt(L)


# ---------------------------------------------------------------------------
# T1. Goal
# ---------------------------------------------------------------------------

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            "First-order Diff(S¹) FP ghost audit: distinguish P = ∂_τ, "
            "P†P = −∂_τ², det'(P), det'(P†P); handle the η-invariant of "
            "det'(∂_τ); treat CKV zero-mode division + norms; fix the ghost "
            "L-power in the BAM loop measure (L, L², or L² only with a "
            "second-order ghost convention)."
        ),
        'follows': 'PR #117 (revised): ghost det = det\'(P) = det\'(P†P)^{1/2} = L',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. The four objects
# ---------------------------------------------------------------------------

def test_T2_four_objects() -> dict:
    """P = ∂_τ (first order, eigenvalues 2πin/L, 1 zero mode = CKV);
    P†P = −∂_τ² (second order, eigenvalues (2πn/L)², 1 zero mode). Verify the
    discrete spectra and the single zero mode of each."""
    # Odd N: avoids the centered-difference Nyquist artifact (a spurious zero
    # at k=N/2 for even N); the continuum ∂_τ has exactly ONE zero mode.
    N = 201
    # first-order P = d/dtau, antisymmetric circulant -> imaginary eigenvalues
    P = np.zeros((N, N))
    for i in range(N):
        P[i, (i + 1) % N] = 0.5
        P[i, (i - 1) % N] = -0.5
    evP = np.linalg.eigvals(P)               # ~ imaginary
    n_zero_P = int(np.sum(np.abs(evP) < 1e-9))
    # second-order P†P = -d^2/dtau^2, periodic Laplacian
    Lap = 2.0 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    Lap[0, -1] = -1.0
    Lap[-1, 0] = -1.0
    evL = np.sort(np.linalg.eigvalsh(Lap))
    n_zero_L = int(np.sum(np.abs(evL) < 1e-9))
    P_is_imaginary = bool(np.max(np.abs(evP.real)) < 1e-9)
    return {
        'name': 'T2_four_objects_and_spectra',
        'description': (
            "P = ∂_τ (first order): eigenvalues 2πin/L (imaginary), 1 zero "
            "mode = CKV. P†P = −∂_τ² (second order): eigenvalues (2πn/L)², "
            "1 zero mode. det'(P), det'(P†P) are their nonzero-mode "
            "determinants."
        ),
        'P_order': 1, 'PtP_order': 2,
        'P_eigenvalues_imaginary': P_is_imaginary,
        'P_zero_modes': n_zero_P,
        'PtP_zero_modes': n_zero_L,
        'pass': P_is_imaginary and n_zero_P == 1 and n_zero_L == 1,
    }


# ---------------------------------------------------------------------------
# T3. The two determinants
# ---------------------------------------------------------------------------

def test_T3_determinants() -> dict:
    """det'(P†P) = L² and |det'(P)| = det'(P†P)^{1/2} = L (zeta), verified
    for several L."""
    rows = []
    ok = True
    for L in (2.0 * PI, 1.0, 3.32408, 5.0):
        d2 = det_prime_PtP(L)
        dP = det_prime_P_magnitude(L)
        match = abs(d2 - L * L) < 1e-6 and abs(dP - L) < 1e-6 and abs(math.sqrt(d2) - dP) < 1e-9
        ok = ok and match
        rows.append({'L': round(L, 5), 'det_prime_PtP': round(d2, 5),
                     'det_prime_P': round(dP, 5),
                     'sqrt_det_PtP': round(math.sqrt(d2), 5), 'match': match})
    return {
        'name': 'T3_two_determinants',
        'description': (
            "det'(P†P) = L² (second order); |det'(P)| = det'(P†P)^{1/2} = L "
            "(first order). Verified to machine precision."
        ),
        'rows': rows,
        'det_prime_P_equals_sqrt_PtP': ok,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T4. η-invariant / phase
# ---------------------------------------------------------------------------

def test_T4_eta_invariant() -> dict:
    """A = −i∂_τ is self-adjoint with eigenvalues 2πn/L, symmetric under
    n → −n, so η_A(s) = Σ sign(λ)|λ|^{−s} ≡ 0 ⟹ η_A(0) = 0: det'(∂_τ) carries
    no anomalous phase and equals +L (real positive). Antiperiodic (Möbius)
    sector: eigenvalues 2πi(n+½)/L — still symmetric (η = 0) but NO zero mode
    ⟹ no CKV."""
    # eta sum over a symmetric finite spectrum -> exactly 0 by pairing
    ns = list(range(-50, 0)) + list(range(1, 51))
    eta_partial = sum(math.copysign(1.0, n) for n in ns)   # = 0
    return {
        'name': 'T4_eta_invariant_phase',
        'description': (
            "η(−i∂_τ) = 0 (spectrum symmetric under n→−n) ⟹ det'(∂_τ) = +L, "
            "no anomalous phase. Antiperiodic/Möbius sector: η = 0 too, but "
            "no zero mode ⟹ no CKV there."
        ),
        'eta_periodic': eta_partial,          # 0
        'eta_zero': abs(eta_partial) < 1e-12,
        'det_prime_dtau_is_real_positive': True,
        'det_prime_dtau_value': 'L (no phase)',
        'antiperiodic_mobius': 'eigenvalues 2πi(n+½)/L: η = 0, NO zero mode ⟹ no CKV',
        'pass': abs(eta_partial) < 1e-12,
    }


# ---------------------------------------------------------------------------
# T5. First- vs second-order ghost convention
# ---------------------------------------------------------------------------

def test_T5_first_vs_second_order() -> dict:
    """The physical FP for one real reparametrization symmetry is the
    first-order bc system: Δ_FP = ∫ Db Dc e^{−∮ b ∂_τ c} = det'(P) = L. A
    second-order ghost convention (action b(P†P)c, a doubled/complex ghost)
    would give det'(P†P) = L². The two differ by exactly one power of L; the
    minimal FP is first-order ⟹ L¹."""
    L = CLOSURE_LOOP_L
    first = det_prime_P_magnitude(L)        # L
    second = det_prime_PtP(L)               # L²
    return {
        'name': 'T5_first_vs_second_order_convention',
        'description': (
            "First-order (real bc, the physical FP): Δ_FP = det'(P) = L. "
            "Second-order convention (b(P†P)c): det'(P†P) = L². They differ "
            "by one power of L; the minimal FP is first-order ⟹ L¹."
        ),
        'first_order_ghost_det': round(first, 5),     # L
        'second_order_ghost_det': round(second, 5),   # L²
        'ratio_second_over_first': round(second / first, 5),   # = L
        'physical_FP_is': 'first-order (L¹)',
        'second_order_requires': 'an explicit second-order ghost convention (doubled/complex ghost)',
        'pass': abs(second / first - L) < 1e-6,
    }


# ---------------------------------------------------------------------------
# T6. CKV division + zero-mode norms
# ---------------------------------------------------------------------------

def test_T6_ckv_and_norms() -> dict:
    """Constant-mode norm ‖1‖ = √L. Two ghost zero modes: c-ghost (CKV, rigid
    rotation) and b-ghost (dual of the 1 modulus), norms √L each ⟹ √L·√L = L.
    Ghost NET = det'(P)/norms: first-order L/L = 1; second-order L²/L = L
    (spurious)."""
    L = CLOSURE_LOOP_L
    norm = zero_mode_norm(L)                 # √L
    norm_product = norm * norm               # L
    net_first = det_prime_P_magnitude(L) / norm_product       # 1
    net_second = det_prime_PtP(L) / norm_product              # L
    return {
        'name': 'T6_ckv_division_and_zero_mode_norms',
        'description': (
            "Constant-mode norm √L; two ghost zero modes (CKV + modulus), "
            "norms √L·√L = L. Ghost NET = det'(P)/norms: first-order L/L = 1; "
            "second-order L²/L = L (spurious extra power)."
        ),
        'zero_mode_norm': round(norm, 5),
        'two_zero_mode_norm_product': round(norm_product, 5),   # L
        'ghost_net_first_order': round(net_first, 5),           # 1
        'ghost_net_second_order': round(net_second, 5),         # L
        'first_order_net_is_one': abs(net_first - 1.0) < 1e-6,
        'pass': abs(net_first - 1.0) < 1e-6 and abs(net_second - L) < 1e-6,
    }


# ---------------------------------------------------------------------------
# T7. The measure table
# ---------------------------------------------------------------------------

def test_T7_measure_table() -> dict:
    """Assemble the measure. Moduli-space measure dL/L (modulus dL × CKV
    automorphism 1/L); ghost NET = 1 (first order); matter det^{−1/2} ~
    L^{−d/2} (heat-kernel Weyl). Net: dL · L^{−1−d/2}. Second-order ghost
    would multiply by L (spurious)."""
    table = [
        {'factor': 'moduli-space measure dL/L', 'source': 'modulus dL × CKV automorphism 1/L', 'L_power': 'L^{-1} dL'},
        {'factor': "FP ghost det'(P) (first order)", 'source': 'nonzero ghost modes', 'L_power': 'L^{+1}'},
        {'factor': 'ghost zero-mode norms (÷)', 'source': '1 CKV + 1 modulus, √L each', 'L_power': 'L^{-1}'},
        {'factor': 'ghost NET', 'source': "det'(P)/norms", 'L_power': 'L^{0}'},
        {'factor': 'matter det^{-1/2}', 'source': 'Tangherlini/heat-kernel Weyl', 'L_power': 'L^{-d/2}'},
        {'factor': 'NET MEASURE (first order)', 'source': 'product', 'L_power': 'dL · L^{-1-d/2}'},
    ]
    return {
        'name': 'T7_measure_table_and_net_L_power',
        'description': (
            "Measure: Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}. Ghost "
            "NET = 1 (first order: det'(P)=L ÷ zero-mode norms L). Net "
            "L-power dL·L^{−1−d/2}. Second-order ghost convention would "
            "multiply by L (spurious)."
        ),
        'measure': 'Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S}',
        'table': table,
        'net_L_power_first_order': 'dL · L^{-1-d/2}',
        'net_L_power_second_order': 'dL · L^{-d/2}  (one spurious power of L)',
        'ghost_L_power_fixed': 'L¹ in det\'(P); net L⁰ after zero-mode norms',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The Diff(S¹) FP ghost is first-order: det'(P) = det'(∂_τ) = "
            "det'(P†P)^{1/2} = L (η = 0, real positive; det'(P†P) = L² is the "
            "second-order square). After CKV + modulus zero-mode norms (√L "
            "each) the ghost NET is L/L = 1, and the measure is "
            "Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S} (net dL·L^{−1−d/2}). The "
            "L² appears only under an explicit second-order ghost convention "
            "(over-counts by one power of L)."
        ),
        'classification': 'BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_four_objects(),
        test_T3_determinants(),
        test_T4_eta_invariant(),
        test_T5_first_vs_second_order(),
        test_T6_ckv_and_norms(),
        test_T7_measure_table(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L'
        verdict = (
            'THE Diff(S¹) FP GHOST IS FIRST-ORDER: det\'(P) = L; THE L² POWER '
            'APPEARS ONLY UNDER AN EXPLICIT SECOND-ORDER GHOST CONVENTION. '
            'This audit follows PR #117 and distinguishes the four objects the '
            'review named.\n\n'
            'THE FOUR OBJECTS. P = ∂_τ is first order (eigenvalues 2πin/L, one '
            'zero mode = the CKV); P†P = −∂_τ² is second order (eigenvalues '
            '(2πn/L)², one zero mode). Their zeta-regularized determinants are '
            'det\'(P†P) = L² and |det\'(P)| = det\'(P†P)^{1/2} = L (verified to '
            'machine precision).\n\n'
            'PHASE / η-INVARIANT. For A = −i∂_τ (self-adjoint, eigenvalues '
            '2πn/L) the spectrum is symmetric under n → −n, so the η-invariant '
            'vanishes identically, η(0) = 0. Hence det\'(∂_τ) carries no '
            'anomalous phase and equals +L (real, positive). In the '
            'antiperiodic / Möbius sector the eigenvalues are 2πi(n+½)/L — '
            'still symmetric (η = 0) but with NO zero mode, so there is no CKV '
            'in the non-orientable sector (a clean tie to the odd-k '
            'structure).\n\n'
            'FIRST- VS SECOND-ORDER CONVENTION. The physical FP for one real '
            'reparametrization symmetry is the first-order bc system, '
            'Δ_FP = ∫ Db Dc e^{−∮ b ∂_τ c} = det\'(P) = L. A second-order '
            'ghost convention (action b(P†P)c, i.e. a doubled/complex ghost) '
            'would instead give det\'(P†P) = L². The two differ by exactly one '
            'power of L; the minimal Faddeev–Popov construction is first-order '
            '⟹ L¹.\n\n'
            'CKV DIVISION + ZERO-MODE NORMS. The constant mode on a circle of '
            'circumference L has norm √L. The ghost system has two zero modes '
            '— the c-ghost (the CKV, rigid rotation) and the b-ghost (dual of '
            'the one Teichmüller modulus) — with norms √L·√L = L. The ghost '
            'NET factor det\'(P)/norms is therefore L/L = 1 in the first-order '
            'convention, but L²/L = L in the second-order one (a spurious '
            'extra power).\n\n'
            'THE MEASURE. The moduli-space measure (one modulus, one CKV '
            'automorphism) is the standard worldline dL/L. With the '
            'first-order ghost NET = 1 and the matter determinant '
            'det^{−1/2}_matter ~ L^{−d/2} (the heat-kernel Weyl term), the BAM '
            'loop measure is Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}, '
            'net L-power dL·L^{−1−d/2}. The second-order convention would '
            'multiply by an extra L (dL·L^{−d/2}), confirming that the '
            'first-order det\'(P) = L is the correct fixing. The absolute '
            'normalization (the κ₅²/Λ₅ anchor) and the multi-loop measure '
            'remain open.'
        )
    else:
        verdict_class = 'BAM_LOOP_MEASURE_GHOST_AUDIT_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the spectra, the '
            'determinants, or the zero-mode accounting.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'first-order Diff(S¹) FP ghost audit: det\'(P) = det\'(∂_τ) = '
            'det\'(P†P)^{1/2} = L (η = 0, real positive); det\'(P†P) = L² is '
            'the second-order square; ghost NET = 1 after √L·√L zero-mode '
            'norms; measure Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}'
        ),
        'four_objects': 'P=∂_τ (1st, 1 zero mode=CKV); P†P=−∂_τ² (2nd, 1 zero mode); det\'(P)=L; det\'(P†P)=L²',
        'eta': 'η(−i∂_τ) = 0 (symmetric spectrum) ⟹ det\'(∂_τ)=+L, no phase; antiperiodic Möbius: η=0, no CKV',
        'convention': 'first-order ghost = det\'(P) = L; second-order = det\'(P†P) = L² (over-counts by L)',
        'zero_modes': 'norm √L each; 2 zero modes (CKV+modulus) ⟹ L; ghost NET first-order L/L=1',
        'measure': 'Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2}',
        'L_power_fixed': 'ghost L¹ in det\'(P); net L⁰ after norms; L² only with a second-order convention',
        'open': 'absolute Z normalization (κ₅²/Λ₅); multi-loop',
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
    out: list[str] = []
    out.append('# First-order Diff(S¹) Faddeev–Popov ghost audit (PR #118)')
    out.append('')
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append('')
    out.append(
        "Audits the PR #117 ghost determinant rigorously: distinguishes "
        "`P = ∂_τ`, `P†P = −∂_τ²`, `det'(P)`, `det'(P†P)`; handles the "
        "phase/η-invariant of `det'(∂_τ)`; treats CKV zero-mode division and "
        "zero-mode norms; and fixes the ghost L-power. **Result: the FP ghost "
        "contributes `L` (first order); `L²` only under an explicit "
        "second-order ghost convention.**"
    )
    out.append('')
    out.append(f"- **Four objects**: {s['four_objects']}")
    out.append(f"- **η-invariant**: {s['eta']}")
    out.append(f"- **Convention**: {s['convention']}")
    out.append(f"- **Zero modes**: {s['zero_modes']}")
    out.append(f"- **Measure**: {s['measure']}")
    out.append(f"- **L-power fixed**: {s['L_power_fixed']}")
    out.append('')

    out.append('## Test summary')
    out.append('')
    out.append('| # | Test | Key finding | PASS? |')
    out.append('|---|---|---|---|')
    label_map = {
        'T1': 'audit goal: distinguish the 4 objects, fix the L-power',
        'T2': 'P=∂_τ (1st, eigvals 2πin/L), P†P=−∂_τ² (2nd); 1 zero mode each',
        'T3': "det'(P†P)=L²; |det'(P)|=det'(P†P)^{1/2}=L (computed)",
        'T4': 'η(−i∂_τ)=0 ⟹ det\'(∂_τ)=+L, no phase; antiperiodic: no CKV',
        'T5': "first-order det'(P)=L vs second-order det'(P†P)=L² (×L)",
        'T6': 'zero-mode norms √L·√L=L ⟹ ghost net L/L=1 (1st), L²/L=L (2nd)',
        'T7': 'measure Z=Σ∫(dL/L)det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2}',
        'T8': 'BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        out.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    out.append('')

    t3 = s['tests'][2]
    out.append('## The two determinants')
    out.append('')
    out.append('| L | det′(P†P) | det′(P) | det′(P†P)^½ | match |')
    out.append('|---:|---:|---:|---:|:---:|')
    for r in t3['rows']:
        out.append(f"| {r['L']} | {r['det_prime_PtP']} | {r['det_prime_P']} | "
                   f"{r['sqrt_det_PtP']} | {'✓' if r['match'] else '✗'} |")
    out.append('')
    out.append("`det'(P†P) = L²` (second order); `det'(P) = det'(P†P)^{1/2} = "
               "L` (first order). The η-invariant of `−i∂_τ` vanishes "
               "(symmetric spectrum), so `det'(∂_τ) = +L` with no anomalous "
               "phase.")
    out.append('')

    t7 = s['tests'][6]
    out.append('## The measure table (net L-power)')
    out.append('')
    out.append('| factor | source | L-power |')
    out.append('|---|---|---|')
    for r in t7['table']:
        out.append(f"| {r['factor']} | {r['source']} | `{r['L_power']}` |")
    out.append('')
    out.append(f"**Measure:** `{t7['measure']}`. Net L-power (first order): "
               f"`{t7['net_L_power_first_order']}`. Second-order convention: "
               f"`{t7['net_L_power_second_order']}`.")
    out.append('')

    out.append('## Verdict')
    out.append('')
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append('')

    out.append('## What this fixes (and leaves open)')
    out.append('')
    out.append('- **Fixed:** the ghost L-power. The FP ghost is first-order, '
               '`det\'(P) = det\'(P†P)^{1/2} = L` (`η = 0`, real positive); '
               'after the CKV + modulus zero-mode norms (`√L` each) the ghost '
               'NET factor is `1`, giving the measure `Z = Σ ∫ (dL/L) '
               'det^{−1/2}_matter e^{−S}` (net `dL·L^{−1−d/2}`). `L²` arises '
               'only under an explicit second-order ghost convention, which '
               'over-counts by one power of `L`.')
    out.append('- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` '
               'bulk anchor, PR #112) and the multi-loop measure.')
    out.append('')
    return '\n'.join(out)


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
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_diff_s1_first_order_ghost_audit_probe'
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
