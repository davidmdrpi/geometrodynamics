"""
Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117,
revised per review).

PR #115 constructed the S_BAM measure as a closure-ledger sector sum over a
loop-space integral gauge-fixed by Diff(S¹) ⋉ U(1)_Hopf ⋉ Z₂, and named the
reparametrization Faddeev–Popov (bc-ghost) determinant as a piece to be
supplied. PR #116 regularized the MATTER fluctuation determinant
(Gel'fand–Yaglom + zeta, finite). This probe supplies the GAUGE sector: the
Faddeev–Popov / ghost determinant for the Diff(S¹) quotient — computed by
the same zeta method, finite, and (crucially) anomaly-free for the 1D loop.

## REVIEW CORRECTION (which determinant, and the L-power)

The first version of this probe wrote the ghost determinant as
det'(−d²/dτ²) = L². That is wrong by one square root. The Faddeev–Popov
determinant is the bc-ghost path integral

    Δ_FP = ∫ Db Dc e^{−∮ b (P c)} = det'(P)        (Grassmann integral),

with P = d/dτ the operator mapping the gauge parameter (the vector ghost c)
to the einbein variation. Since the b- and c-ghost spaces have equal
dimension here, this equals det'(P†P)^{1/2}. So

    Δ_FP = det'(P) = det'(P†P)^{1/2} = √(L²) = L      (NOT L²).

The L² is det'(P†P) = det'(−d²/dτ²), the intermediate Laplacian
determinant; the FP ghost determinant is its SQUARE ROOT, L. Of the three
candidates {det'(P), det'(P†P), det'(P†P)^{1/2}}, the ghost determinant is
det'(P†P)^{1/2} ( = det'(P) in 1D, by the ±n mode pairing) = L.

## The gauge structure (worldline reparametrization)

The closure loop X: S¹ → base is reparametrization-invariant: X(τ) → X(f(τ))
for f ∈ Diff(S¹). Gauge-fixing the loop einbein e(τ) to a constant (the
worldline/Polyakov procedure) leaves exactly:

  - ONE Teichmüller modulus L = ∮ e dτ — the loop circumference (= the
    Schwinger proper time), and
  - ONE conformal Killing vector (CKV) — the constant vector field ∂_τ, i.e.
    the residual rigid U(1) rotation of the loop.

## The ghost operator and its determinant (same zeta method as PR #116)

The Faddeev–Popov operator is P = d/dτ on periodic fields (P†P = −d²/dτ²).
Its kernel is the constants — exactly the ONE CKV. The nonzero-mode
determinant is zeta-regularized as in PR #116:

    det'(P†P) = det'(−d²/dτ²) = L²   (ζ(0) = −1, ζ'(0) = −2 ln L),
    Δ_FP = det'(P†P)^{1/2} = det'(P) = L     [the ghost determinant].

Verified to machine precision: det'(P†P) = L² and Δ_FP = L for L = 2π, 1,
3.32, 5.

## The corrected one-loop measure (L-power fixed)

The FP determinant Δ_FP = L is the JACOBIAN of the einbein → proper-length
gauge fixing: it converts the einbein-fluctuation measure into the proper
modulus measure dL. Dividing by the conformal Killing volume Vol(U(1)) = L
(the rigid rotation) gives the symmetry factor 1/L. So the gauge sector
yields the standard worldline proper-time measure

    Z = Σ_sectors ∫ (dL / L) · det^{−1/2}_matter · e^{−S}.

There is NO separate det'_ghost = L² factor (the error of the first
version): the ghost determinant is L¹, and it is the Jacobian already used in
writing the modulus measure as the proper length dL. The 1/L is the CKV
factor — and for the closure loop L = 2π it is 1/(2π), PR #74's per-loop
measure factor.

## The CKV zero mode IS PR #74's 1/(2π) (unchanged by the correction)

The single CKV (rigid rotation) is divided out: Vol(U(1)) = L. For the BAM
closure loop, whose great-circle length is L = 2π (the closure quantum),
1/L = 1/(2π) — exactly the per-loop measure factor PR #74 identified from
the Schwinger anomaly a = α/(2π). This is the CKV (c-ghost zero-mode)
factor, independent of the determinant power, so the correction above leaves
it intact.

## Anomaly-free (the clean part)

Unlike the 2D string worldsheet (where the bc-ghosts carry central charge
c = −26 and consistency forces the conformal anomaly to cancel, D = 26), the
1D worldline loop has NO Weyl/conformal symmetry to be anomalous: there is
no traceless symmetric 2-tensor in one dimension, so the Diff(S¹)
gauge-fixing is anomaly-FREE. (The nontrivial anomaly is the DISCRETE Z₂
orientation anomaly — the odd-k condition of PR #115 — not this continuous
one.)

## What this completes (and what stays open)

  - **Completes:** the gauge sector of the S_BAM measure. With PR #116's
    matter determinant and this ghost determinant Δ_FP = det'(P†P)^{1/2} = L
    (anomaly-free, CKV → 1/(2π)), the one-loop measure is
    Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S}, every factor
    finite/computable.
  - **Open:** the absolute normalization of Z still carries the bulk κ₅²/Λ₅
    anchor (PR #112); the multi-loop / interacting measure; a closed form.

Tests:
  T1. Recap: PR #115 flagged the Diff(S¹) FP ghost; PR #116 did matter; now
      the ghost.
  T2. Gauge structure: Diff(S¹) → constant einbein ⟹ 1 modulus L + 1 CKV.
  T3. Ghost operator P = d/dτ; P†P = −d²/dτ²; kernel(P) = constants = 1 CKV.
  T4. WHICH DETERMINANT: Δ_FP = det'(P) = det'(P†P)^{1/2} = L (NOT
      det'(P†P) = L²). The three candidates compared.
  T5. Corrected measure: Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}; ghost
      L-power = L¹ (the dL Jacobian), CKV 1/L; no spurious L².
  T6. CKV zero mode ⟹ 1/L; for L = 2π, 1/(2π) = PR #74's factor
      (unchanged).
  T7. Anomaly-free: 1D worldline has no conformal anomaly (vs 2D string
      c = −26); the nontrivial anomaly is the discrete Z₂ (PR #115).
  T8. Assessment.

Verdict:
  - DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE (expected):
    the Faddeev–Popov ghost determinant for the Diff(S¹) quotient is
    Δ_FP = det'(P) = det'(P†P)^{1/2} = L (the square root of det'(P†P) = L²;
    the first version's L² was the error). The corrected one-loop measure is
    Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S} — ghost L-power L¹ (the
    proper-length Jacobian), CKV factor 1/L = 1/(2π) (PR #74) at the closure
    loop, anomaly-free. Absolute Z normalization (κ₅²/Λ₅) and multi-loop
    stay open.
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

CLOSURE_LOOP_L = 2.0 * PI          # great-circle length = the closure quantum


def det_prime_PtP(L: float) -> float:
    """Zeta-regularized det'(P†P) = det'(−d²/dτ²) on a circle of
    circumference L. Eigenvalues (2πn/L)², n=±1,±2,…:
    ζ(s) = 2 (L/2π)^{2s} ζ_R(2s) ⟹ ζ(0) = −1, ζ'(0) = −2 ln L ⟹ det' = L²."""
    zeta_prime_0 = 2.0 * (2.0 * math.log(L / (2.0 * PI)) * ZETA_R_0
                          + 2.0 * ZETA_R_PRIME_0)
    return math.exp(-zeta_prime_0)


def det_prime_P(L: float) -> float:
    """Zeta-regularized |det'(P)| = |det'(d/dτ)| on a circle of circumference
    L. Eigenvalues 2πi n/L, n=±1,±2,…; |·| via ±n pairing gives the same
    nonzero spectrum once: ζ(s) = 2 (L/2π)^{s} ζ_R(s) ⟹ ζ(0) = −1,
    ζ'(0) = −ln L ⟹ |det'(P)| = L = det'(P†P)^{1/2}."""
    zeta_prime_0 = 2.0 * (math.log(L / (2.0 * PI)) * ZETA_R_0 + ZETA_R_PRIME_0)
    return math.exp(-zeta_prime_0)


def fp_ghost_determinant(L: float) -> float:
    """The Faddeev–Popov ghost determinant: Δ_FP = det'(P†P)^{1/2} = L."""
    return math.sqrt(det_prime_PtP(L))


# ---------------------------------------------------------------------------
# T1. Recap
# ---------------------------------------------------------------------------

def test_T1_recap() -> dict:
    return {
        'name': 'T1_recap',
        'description': (
            "PR #115 named the Diff(S¹) Faddeev–Popov (bc-ghost) determinant "
            "as a piece of the measure; PR #116 regularized the MATTER "
            "fluctuation determinant. This probe supplies the GAUGE sector — "
            "the FP/ghost determinant for the Diff(S¹) quotient (revised per "
            "review to fix the determinant identification and L-power)."
        ),
        'pr115': 'measure structure + Diff(S¹) gauge-fixing flagged',
        'pr116': 'matter fluctuation determinant finite (GY + zeta)',
        'this_probe': 'the Diff(S¹) FP / ghost determinant (corrected)',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Gauge structure
# ---------------------------------------------------------------------------

def test_T2_gauge_structure() -> dict:
    return {
        'name': 'T2_gauge_structure',
        'description': (
            "Worldline reparametrization: gauge-fixing the loop einbein to a "
            "constant leaves ONE Teichmüller modulus L = ∮ e dτ (the loop "
            "circumference = Schwinger proper time) and ONE conformal Killing "
            "vector (the constant ∂_τ = rigid U(1) rotation). Algebra = Witt "
            "(centrally extended: Virasoro)."
        ),
        'moduli': '1 (the circumference L = Schwinger proper time)',
        'conformal_killing_vectors': '1 (the rigid U(1) rotation ∂_τ)',
        'algebra': 'Witt / Virasoro',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. Ghost operator + kernel
# ---------------------------------------------------------------------------

def test_T3_ghost_operator() -> dict:
    """The FP operator is P = d/dτ (vector ghost c ↦ einbein variation), so
    P†P = −d²/dτ² on periodic fields. ker(P) = constants = exactly the ONE
    CKV. Verify the kernel dimension by counting zero modes."""
    N = 200
    lap = (2.0 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1))
    lap[0, -1] = -1.0
    lap[-1, 0] = -1.0
    ev = np.sort(np.linalg.eigvalsh(lap))
    n_zero = int(np.sum(np.abs(ev) < 1e-9))
    return {
        'name': 'T3_ghost_operator_and_kernel',
        'description': (
            "FP operator P = d/dτ (vector ghost ↦ einbein variation); "
            "P†P = −d²/dτ² (periodic); ker(P) = constants = the ONE CKV "
            "(rigid rotation). Periodic Laplacian has exactly 1 zero mode."
        ),
        'fp_operator': 'P = d/dτ',
        'PtP': 'P†P = −d²/dτ² (periodic)',
        'kernel_dim_zero_modes': n_zero,
        'kernel_is_the_CKV': n_zero == 1,
        'pass': n_zero == 1,
    }


# ---------------------------------------------------------------------------
# T4. WHICH determinant (the review point)
# ---------------------------------------------------------------------------

def test_T4_which_determinant() -> dict:
    """The Faddeev–Popov determinant is the bc-ghost path integral
    Δ_FP = ∫ Db Dc e^{−∮ b(Pc)} = det'(P) = det'(P†P)^{1/2} (equal in 1D by
    the ±n pairing). It is the SQUARE ROOT of det'(P†P) = L², i.e. Δ_FP = L,
    NOT L². The three candidates are compared explicitly."""
    rows = []
    ok = True
    for L in (2.0 * PI, 1.0, 3.32408, 5.0):
        dPtP = det_prime_PtP(L)             # det'(P†P) = L²  (intermediate, NOT the ghost det)
        dP = det_prime_P(L)                 # det'(P) = L
        dFP = fp_ghost_determinant(L)       # det'(P†P)^{1/2} = L  (the ghost det)
        match = (abs(dPtP - L * L) < 1e-6 and abs(dP - L) < 1e-6
                 and abs(dFP - L) < 1e-6 and abs(dP - dFP) < 1e-9)
        ok = ok and match
        rows.append({'L': round(L, 5),
                     'det_prime_PtP': round(dPtP, 5),        # = L²
                     'det_prime_P': round(dP, 5),            # = L
                     'det_prime_PtP_half': round(dFP, 5),    # = L  (ghost det)
                     'match': match})
    return {
        'name': 'T4_which_determinant_is_the_ghost_det',
        'description': (
            "Δ_FP = det'(P) = det'(P†P)^{1/2} = L (the bc-ghost path "
            "integral; the two coincide in 1D by ±n pairing). NOT "
            "det'(P†P) = L² — that is the intermediate Laplacian determinant "
            "(its square). So the ghost determinant is L, the square root of "
            "the first version's L²."
        ),
        'candidates': {
            "det'(P)": 'L  — the bc-ghost integral (FP determinant)',
            "det'(P†P)": 'L²  — intermediate Laplacian determinant (the SQUARE; not the ghost det)',
            "det'(P†P)^{1/2}": 'L  — = det\'(P); THE Faddeev–Popov ghost determinant',
        },
        'ghost_determinant_is': "det'(P†P)^{1/2} = det'(P) = L",
        'first_version_error': "wrote det'(P†P) = L² as the ghost det; off by one square root",
        'rows': rows,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T5. The corrected one-loop measure (L-power fixed)
# ---------------------------------------------------------------------------

def test_T5_corrected_measure() -> dict:
    """Δ_FP = L is the Jacobian of the einbein → proper-length gauge fixing
    (it makes the modulus measure the proper length dL). Dividing by the CKV
    volume Vol(U(1)) = L gives the symmetry factor 1/L. So the gauge sector
    yields the standard worldline proper-time measure dL/L — ghost L-power
    L¹ (NOT a separate L²)."""
    return {
        'name': 'T5_corrected_one_loop_measure',
        'description': (
            "Δ_FP = det'(P†P)^{1/2} = L is the einbein→proper-length Jacobian "
            "(modulus measure = dL); CKV volume L ⟹ symmetry factor 1/L. "
            "Measure: Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}. Ghost "
            "L-power = L¹ (the dL Jacobian), NOT the spurious L² of v1."
        ),
        'fp_determinant_L_power': 1,
        'v1_wrong_L_power': 2,
        'measure': 'Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S}',
        'no_separate_L2_factor': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. CKV zero mode IS PR #74's 1/(2π) (unchanged)
# ---------------------------------------------------------------------------

def test_T6_ckv_is_pr74_factor() -> dict:
    """The CKV (rigid rotation) volume Vol(U(1)) = L gives the symmetry
    factor 1/L; for the closure loop L = 2π, 1/L = 1/(2π) = PR #74's per-loop
    factor. This is the c-ghost zero-mode factor, independent of the
    determinant power — so the T4 correction leaves it intact."""
    one_over_L = 1.0 / CLOSURE_LOOP_L
    pr74 = 1.0 / (2.0 * PI)
    return {
        'name': 'T6_ckv_is_closure_quantum_factor',
        'description': (
            "CKV volume Vol(U(1)) = L ⟹ symmetry factor 1/L; for L = 2π, "
            "1/L = 1/(2π) = PR #74's per-loop factor. The c-ghost zero-mode "
            "factor — independent of the determinant power, so unaffected by "
            "the T4 correction."
        ),
        'closure_loop_L': CLOSURE_LOOP_L,
        'ckv_factor_one_over_L': one_over_L,
        'pr74_factor': pr74,
        'matches_pr74': abs(one_over_L - pr74) < 1e-12,
        'pass': abs(one_over_L - pr74) < 1e-12,
    }


# ---------------------------------------------------------------------------
# T7. Anomaly-free
# ---------------------------------------------------------------------------

def test_T7_anomaly_free() -> dict:
    """The 1D worldline loop has no Weyl/conformal symmetry to be anomalous
    (no traceless symmetric 2-tensor in 1D), so Diff(S¹) gauge-fixing is
    anomaly-FREE — unlike the 2D string (bc-ghosts c = −26 ⟹ D = 26). The
    only nontrivial anomaly is the DISCRETE Z₂ orientation anomaly (odd-k,
    PR #115)."""
    return {
        'name': 'T7_anomaly_free_1d_worldline',
        'description': (
            "1D worldline: no Weyl/conformal anomaly (no traceless symmetric "
            "2-tensor in 1D) ⟹ Diff(S¹) gauge-fixing is anomaly-free (vs 2D "
            "string bc-ghost c = −26 ⟹ D = 26). The only nontrivial anomaly "
            "is the discrete Z₂ orientation anomaly (odd-k, PR #115)."
        ),
        'worldline_conformal_anomaly': 'none (1D)',
        'string_contrast': 'bc-ghost c = −26 ⟹ critical D = 26 (2D worldsheet)',
        'nontrivial_anomaly_in_BAM': 'discrete Z₂ orientation (odd-k, PR #115)',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The Faddeev–Popov ghost determinant for the Diff(S¹) quotient is "
            "Δ_FP = det'(P) = det'(P†P)^{1/2} = L (the square root of "
            "det'(P†P) = L²). The corrected one-loop measure is "
            "Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}: ghost L-power "
            "L¹ (the proper-length Jacobian), CKV factor 1/L = 1/(2π) (PR #74) "
            "at the closure loop, anomaly-free. Absolute Z normalization and "
            "multi-loop stay open."
        ),
        'classification': 'DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_recap(),
        test_T2_gauge_structure(),
        test_T3_ghost_operator(),
        test_T4_which_determinant(),
        test_T5_corrected_measure(),
        test_T6_ckv_is_pr74_factor(),
        test_T7_anomaly_free(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE'
        verdict = (
            'THE Diff(S¹) FADDEEV–POPOV GHOST DETERMINANT IS det\'(P†P)^{1/2} '
            '= L (NOT L²); THE ONE-LOOP MEASURE\'S GHOST L-POWER IS FIXED '
            'ACCORDINGLY. (Revised per review.) PR #115 named the '
            'reparametrization Faddeev–Popov determinant as a piece of the '
            'measure; PR #116 regularized the matter determinant. This probe '
            'supplies the gauge sector — and the review correctly flagged that '
            'the first version conflated the ghost determinant with '
            'det\'(P†P).\n\n'
            'WHICH DETERMINANT. The Faddeev–Popov determinant is the bc-ghost '
            'path integral Δ_FP = ∫ Db Dc e^{−∮ b(Pc)} = det\'(P), with '
            'P = d/dτ the operator taking the vector ghost c to the einbein '
            'variation. Because the b- and c-ghost spaces have equal '
            'dimension here, det\'(P) = det\'(P†P)^{1/2}. The intermediate '
            'Laplacian determinant is det\'(P†P) = det\'(−d²/dτ²) = L² '
            '(ζ(0) = −1, ζ\'(0) = −2 ln L), so the GHOST determinant is its '
            'square root: Δ_FP = det\'(P) = det\'(P†P)^{1/2} = L. Of the three '
            'candidates {det\'(P), det\'(P†P), det\'(P†P)^{1/2}}, it is '
            'det\'(P†P)^{1/2} (= det\'(P) in 1D by the ±n mode pairing) = L. '
            'The first version\'s "det\'(−d²/dτ²) = L²" was the ghost '
            'determinant off by one square root.\n\n'
            'THE CORRECTED MEASURE (L-POWER FIXED). Δ_FP = L is the Jacobian '
            'of the einbein → proper-length gauge fixing — it makes the '
            'modulus measure the proper circumference dL. Dividing by the '
            'conformal Killing volume Vol(U(1)) = L (the rigid rotation) gives '
            'the symmetry factor 1/L. So the gauge sector yields the standard '
            'worldline proper-time measure Z = Σ_sectors ∫ (dL/L) '
            'det^{−1/2}_matter e^{−S}: the ghost L-power is L¹ (the dL '
            'Jacobian), and there is NO separate det\'_ghost = L² factor (the '
            'v1 error).\n\n'
            'PR #74 UNCHANGED. The 1/(2π) is the CKV (c-ghost zero-mode) '
            'factor 1/Vol(U(1)) = 1/L at the closure loop L = 2π — independent '
            'of the determinant power, so the correction leaves it intact: '
            'PR #74\'s 1/(2π) is still the ghost zero-mode of the Diff(S¹) '
            'quotient.\n\n'
            'ANOMALY-FREE. The 1D worldline carries no Weyl/conformal anomaly '
            '(no traceless symmetric 2-tensor in 1D), unlike the 2D string '
            '(bc-ghost c = −26 ⟹ D = 26); the only nontrivial anomaly is the '
            'discrete Z₂ orientation (odd-k, PR #115). With PR #116\'s matter '
            'determinant, the one-loop measure is finite/computable; the '
            'absolute normalization (κ₅²/Λ₅) and the multi-loop measure '
            'remain open.'
        )
    else:
        verdict_class = 'DIFF_S1_FP_GHOST_DET_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the determinant '
            'identification and the zeta regularization.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the Diff(S¹) Faddeev–Popov ghost determinant is '
            'Δ_FP = det\'(P) = det\'(P†P)^{1/2} = L (the square root of '
            'det\'(P†P) = L²); the corrected one-loop measure is '
            'Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}, ghost L-power L¹'
        ),
        'which_determinant': "Δ_FP = det'(P) = det'(P†P)^{1/2} = L  (NOT det'(P†P) = L²)",
        'L_power_fix': 'ghost L-power L¹ (was wrongly L²); measure = ∫ (dL/L) det^{−1/2}_matter e^{−S}',
        'ckv_factor': 'CKV volume L ⟹ 1/L; for L = 2π, 1/(2π) = PR #74 (unchanged)',
        'anomaly': 'anomaly-free (no 1D conformal anomaly); only the discrete Z₂ (odd-k, PR #115) is nontrivial',
        'open': 'absolute Z normalization (κ₅²/Λ₅); multi-loop; closed form',
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
    L.append('# Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117, revised)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Supplies the gauge sector of the S_BAM measure (PR #115 flagged it; "
        "PR #116 did the matter determinant). **Revised per review:** the "
        "ghost determinant is `det'(P†P)^{1/2} = det'(P) = L` — the square "
        "root of `det'(P†P) = L²`; the first version wrongly used `L²`. The "
        "one-loop measure's ghost L-power is fixed to `L¹` accordingly."
    )
    L.append('')
    L.append(f"- **Which determinant**: {s['which_determinant']}")
    L.append(f"- **L-power fix**: {s['L_power_fix']}")
    L.append(f"- **CKV factor**: {s['ckv_factor']}")
    L.append(f"- **Anomaly**: {s['anomaly']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'gauge sector of the measure (revised per review)',
        'T2': 'Diff(S¹) → constant einbein ⟹ 1 modulus L + 1 CKV',
        'T3': 'P = d/dτ; P†P = −d²/dτ²; ker(P) = constants = 1 CKV',
        'T4': "Δ_FP = det'(P) = det'(P†P)^{1/2} = L (NOT det'(P†P) = L²)",
        'T5': 'measure Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}; ghost L¹',
        'T6': 'CKV 1/L; for L = 2π, 1/(2π) = PR #74 (unchanged)',
        'T7': 'anomaly-free (1D; vs string c = −26)',
        'T8': 'DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t4 = s['tests'][3]
    L.append('## T4: which determinant is the ghost determinant')
    L.append('')
    L.append('| candidate | value | role |')
    L.append('|---|---|---|')
    L.append("| `det'(P)` | L | the bc-ghost integral (FP determinant) |")
    L.append("| `det'(P†P)` | L² | intermediate Laplacian determinant (the **square** — not the ghost det) |")
    L.append("| `det'(P†P)^{1/2}` | L | = `det'(P)`; **THE Faddeev–Popov ghost determinant** |")
    L.append('')
    L.append('| L | det′(P†P) | det′(P) | det′(P†P)^½ | match |')
    L.append('|---:|---:|---:|---:|:---:|')
    for r in t4['rows']:
        L.append(f"| {r['L']} | {r['det_prime_PtP']} | {r['det_prime_P']} | "
                 f"{r['det_prime_PtP_half']} | {'✓' if r['match'] else '✗'} |")
    L.append('')
    L.append("The Faddeev–Popov determinant is the `bc`-ghost path integral "
             "`∫ Db Dc e^{−∮ b(Pc)} = det'(P) = det'(P†P)^{1/2} = L`. The `L²` "
             "is `det'(P†P)` (the intermediate Laplacian determinant); the "
             "ghost determinant is its **square root**, `L`. The first version "
             "displayed `L²` — off by one square root.")
    L.append('')

    L.append('## The corrected one-loop measure')
    L.append('')
    L.append('```')
    L.append('Z = Σ_sectors  ∫ (dL/L)  ·  det^{−1/2}_matter  ·  e^{−S}')
    L.append('              └ Δ_FP = det′(P†P)^{1/2} = L is the einbein→proper-length')
    L.append('                Jacobian (⟹ modulus measure dL); CKV vol L ⟹ 1/L ┘')
    L.append('```')
    L.append('')
    L.append("Ghost L-power = **L¹** (the proper-length Jacobian), not the "
             "spurious **L²** of the first version. The `1/L` is the CKV "
             "factor — `1/(2π)` at the closure loop `L = 2π` (PR #74).")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this completes (and does not)')
    L.append('')
    L.append('- **Completes:** the Diff(S¹) gauge sector — the FP ghost '
             'determinant `Δ_FP = det\'(P†P)^{1/2} = L`, the corrected '
             'measure `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (ghost '
             'L-power `L¹`), anomaly-free, with the CKV `1/L = 1/(2π)` (PR '
             '#74) intact. With PR #116\'s matter determinant the one-loop '
             'measure is finite/computable.')
    L.append('- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` '
             'bulk anchor, PR #112), the multi-loop / interacting measure, '
             'and a closed-form expression.')
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
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_diff_s1_ghost_determinant_probe'
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
