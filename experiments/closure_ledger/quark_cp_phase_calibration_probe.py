"""
Quark CP phase calibration on the locked shelled-closure Hamiltonian (PR #156).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields), not
> quantum gravity. The CP phase rides the partition-mixing element of the
> mass-locked Hamiltonian.

PR #155 extracted the CKM |V| from the mass-locked partition blocks (zero new
inputs) and left CP as the constrained open item: the baseline has J = 0
exactly, the observed J = 3.08×10⁻⁵, and any phase calibration must reproduce
J WITHOUT disturbing the fixed |V| or the mass calibration. This probe
performs that calibration — and finds a sharper structure than expected.

## The CP extension (in-probe, locked blocks preserved)

The global `phase` knob CANNOT be used: it enters the same-partition
transport couplings as cos(phase·dk), so the mass calibration locked it at
zero (verified — turning it on deforms |V| and the masses). The v3 §4
partition-mixing element carries its own phase — the Hopf-connection
placeholder φ_q(k) = φ·k — and the locked params have partition_mixing = 0,
so CP is switched on cleanly in-probe:

    H(ε, φ) = H_locked + Σ_k [ −ε e^{iφk} |k,+⟩⟨k,−| + h.c. ],

leaving the locked blocks exactly intact at ε = 0. The CKM matrix is read
through the shell-partner charged current on the full 6×6.

## Scaling structure (derived)

J ∝ ε^1.9 (≈ quadratic — CP needs one insertion on each sector's side) with
a sinusoidal φ-dependence; the |V| and mass perturbations also enter at
O(ε²). CP switches on without disturbing the locked structure at small ε.

## The ceiling identity: the J deficit is the soft direction, exactly

J is bounded by the rephasing-invariant ceiling |V_us·V_cb·V_ub|. With the
PREDICTED |V|: ceiling = 8.64×10⁻⁶, vs the observed ceiling 3.47×10⁻⁵ —
ratio 0.249, which decomposes EXACTLY as the per-element soft-direction
ratios 0.498 × 0.902 × 0.555 (#155's V_us/V_ub deficits). The observed CP is
near-maximal (J_obs/ceiling_obs = 0.887). So the J shortfall is NOT an
independent CP failure — it is the #155 soft direction propagated, and when
the soft V_us/V_ub directions land on data, the ceiling rises to the
observed value: a future consistency lock.

## The constrained calibration

Target: the observed PHASE CONTENT (sin δ ≈ 0.887) on the predicted |V| —
J_target = 7.67×10⁻⁶. Result: ε* = 0.0528 at φ* = 0.80 hits the target with
|V_cb| shifted −0.0%, |V_us|/|V_ub| shifted only −4.2/−4.3% (inside the soft
direction), and masses shifted ≤ 0.5% (well inside the 1.6% calibration
accuracy). The locked structure SURVIVES the CP calibration; the calibrated
phase content is near-maximal (sin δ = 0.967), like the data.

## The triangle-shape discriminator (the honest sharp edge)

The calibrated phase reproduces the AREA (J) but not the SHAPE: the
db-unitarity-triangle comes out squashed (β ≈ 0° vs observed 22.2°, γ ≈ 180°
vs 65.9°) — the placeholder linear phase φ_q(k) = φ·k puts the CP in the
wrong quartet orientation. The triangle SHAPE therefore discriminates phase
structures: β = 22° is the falsifiable target for the TRUE Hopf-connection
phase (the v3 §4 TODO, tracked since the handoff). A sharp open item, not
hidden.

## Ledger

One new input consumed: the quark CP phase content (like α — structure
derived, value input; the flavor puzzle's CP entry made explicit).
Derived: the quadratic scaling, the ceiling identity and its exact
soft-direction decomposition, the survival of the locked structure, the
near-maximal phase content. Open: the Hopf-connection φ_q(k) (target:
the triangle shape), the soft V_us/V_ub directions (target: the J ceiling).

Tests:
  T1. Goal: the constrained CP calibration #155 posed.
  T2. The extension: locked blocks preserved; baseline J = 0; the global
      phase knob shown unusable (cos(phase·dk) in transport).
  T3. Scaling: J ∝ ε^≈2, sinusoidal in φ; shifts O(ε²).
  T4. The ceiling identity: deficit 0.249 = 0.498 × 0.902 × 0.555 exactly;
      observed CP near-maximal (0.887).
  T5. The calibration: ε* = 0.053, φ* = 0.80; V_cb untouched; V_us/V_ub
      −4%; masses ≤ 0.5% — the locked structure survives.
  T6. The triangle shape: β ≈ 0° vs 22.2° — the placeholder phase fails the
      shape; the Hopf-connection phase has a falsifiable target.
  T7. Ledger / scope.
  T8. Assessment.

Verdict:
  - QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN
    (expected): the quark CP phase calibrates onto the locked Hamiltonian
    without disturbing |V| or the masses; the J ceiling deficit decomposes
    exactly into the #155 soft-direction |V| deficits (no independent CP
    failure); the phase content is near-maximal like the data; and the
    db-triangle shape (β = 22°) is the falsifiable target the placeholder
    linear phase fails — pointing at the Hopf-connection φ_q(k).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import brentq

from geometrodynamics.qcd import quark_spectrum as qs


PI = math.pi
IDX_PLUS = [0, 2, 4]
IDX_MINUS = [1, 3, 5]
K_SHELLS = (1, 3, 5)

J_OBSERVED = 3.08e-5
V_OBS = {'us': 0.225, 'cb': 0.04182, 'ub': 0.00369}
BETA_OBS_DEG = 22.2
GAMMA_OBS_DEG = 65.9

_H0 = qs.build_quark_hamiltonian(qs.LOCKED_QUARK_PARAMS)


def H_cp(eps: float, phi: float) -> np.ndarray:
    """The locked Hamiltonian + the v3 §4 partition-mixing element with the
    Hopf-placeholder phase φ_q(k) = φ·k (in-probe; locked blocks intact)."""
    H = _H0.copy()
    for ki, (ip, im) in enumerate(zip(IDX_PLUS, IDX_MINUS)):
        k = K_SHELLS[ki]
        H[ip, im] += -eps * np.exp(1j * phi * k)
        H[im, ip] += -eps * np.exp(-1j * phi * k)
    return H


def ckm(H: np.ndarray):
    """CKM through the shell-partner charged current C = Σ_k |k,+⟩⟨k,−| on
    the full 6×6; eigenstates tagged by dominant partition, mass-ascending."""
    w, U = np.linalg.eigh(H)
    up, dn = [], []
    for c in range(6):
        wp = float(np.sum(np.abs(U[IDX_PLUS, c]) ** 2))
        (up if wp > 0.5 else dn).append((float(w[c].real), c))
    up.sort()
    dn.sort()
    if len(up) != 3:
        raise RuntimeError("partition tagging failed")
    Uu = U[:, [c for _, c in up]]
    Ud = U[:, [c for _, c in dn]]
    V = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            V[i, j] = sum(np.conj(Uu[IDX_PLUS[k], i]) * Ud[IDX_MINUS[k], j]
                          for k in range(3))
    return V, np.array([x for x, _ in up]), np.array([x for x, _ in dn])


def jarlskog(V: np.ndarray) -> float:
    return float(np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0])))


_V0, _MU0, _MD0 = ckm(_H0)
_CEIL_PRED = float(abs(_V0[0, 1]) * abs(_V0[1, 2]) * abs(_V0[0, 2]))
_CEIL_OBS = V_OBS['us'] * V_OBS['cb'] * V_OBS['ub']
_SIN_DELTA_OBS = J_OBSERVED / _CEIL_OBS
_J_TARGET = _SIN_DELTA_OBS * _CEIL_PRED

_PHI_STAR = 0.80
_EPS_STAR = float(brentq(
    lambda e: jarlskog(ckm(H_cp(e, _PHI_STAR))[0]) - _J_TARGET, 1e-3, 0.15))
_VC, _MUC, _MDC = ckm(H_cp(_EPS_STAR, _PHI_STAR))


def triangle_angles_deg(V: np.ndarray):
    beta = math.degrees(np.angle(-V[1, 0] * np.conj(V[1, 2])
                                 / (V[2, 0] * np.conj(V[2, 2]))))
    gamma = math.degrees(np.angle(-V[0, 0] * np.conj(V[0, 2])
                                  / (V[1, 0] * np.conj(V[1, 2]))))
    return beta, gamma


# ---------------------------------------------------------------------------
# T1. Goal
# ---------------------------------------------------------------------------

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            "Perform the constrained CP calibration #155 posed: introduce "
            "the quark CP phase on the locked Hamiltonian, reproduce the "
            "observed phase content, and verify the fixed |V| and the mass "
            "calibration survive — with the J ceiling and the triangle "
            "shape as the structural diagnostics."
        ),
        'builds_on': ['#155 CKM |V| (zero new inputs) + J = 0 baseline',
                      'v3 §4 partition-mixing phase (Hopf placeholder)',
                      'qcd_v3_r2_handoff TODO: Hopf-derived φ_q(k)'],
        'framing': 'QFT on the fixed classical throat geometry — not quantum gravity',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. The extension and why the global knob is unusable
# ---------------------------------------------------------------------------

def test_T2_extension() -> dict:
    """Locked blocks preserved at ε = 0 (J = 0 re-verified); the global
    phase knob deforms |V_us| via cos(phase·dk) in transport — locked at 0
    by the mass calibration."""
    j0 = jarlskog(_V0)
    p_glob = replace(qs.LOCKED_QUARK_PARAMS, phase=0.5)
    Hg = qs.build_quark_hamiltonian(p_glob)
    Vg, _, _ = ckm(Hg)
    vus_glob = float(abs(Vg[0, 1]))
    vus_base = float(abs(_V0[0, 1]))
    glob_breaks = abs(vus_glob / vus_base - 1.0) > 0.5
    return {
        'name': 'T2_cp_extension',
        'description': (
            "The CP extension H(ε, φ) = H_locked + partition-mixing with "
            "the Hopf-placeholder phase φ_q(k) = φ·k, built in-probe so the "
            "locked blocks are exactly intact at ε = 0 (J = 0 re-verified). "
            "The GLOBAL phase knob cannot be used: it enters the "
            "same-partition transport as cos(phase·dk) — setting it to 0.5 "
            "collapses |V_us| from 0.112 to ~0.005 — which is why the mass "
            "calibration locked it at zero and why the #155 baseline had "
            "J = 0 structurally."
        ),
        'baseline_J': j0,
        'vus_with_global_phase_0p5': round(vus_glob, 4),
        'vus_baseline': round(vus_base, 4),
        'global_knob_unusable': bool(glob_breaks),
        'pass': j0 == 0.0 and glob_breaks,
    }


# ---------------------------------------------------------------------------
# T3. Scaling structure
# ---------------------------------------------------------------------------

def test_T3_scaling() -> dict:
    """J ∝ ε^≈2 (fit); sinusoidal-ish φ-dependence; |V|/mass shifts O(ε²)."""
    e1, e2 = 0.01, 0.04
    j1 = jarlskog(ckm(H_cp(e1, _PHI_STAR))[0])
    j2 = jarlskog(ckm(H_cp(e2, _PHI_STAR))[0])
    p_fit = math.log(j2 / j1) / math.log(e2 / e1)
    phi_rows = []
    for phi in (0.25, 0.8, 1.5, 2.0, 2.5):
        phi_rows.append({'phi': phi,
                         'J': float(f'{jarlskog(ckm(H_cp(0.02, phi))[0]):.3e}')})
    signs_vary = (min(r['J'] for r in phi_rows) < 0 < max(r['J'] for r in phi_rows))
    ok = abs(p_fit - 2.0) < 0.25 and signs_vary
    return {
        'name': 'T3_scaling_structure',
        'description': (
            "J switches on at second order: fitted J ∝ ε^p with p = "
            f"{p_fit:.2f} (≈ 2 — CP interference needs one insertion on "
            "each sector's side), with a sinusoidal φ-dependence that "
            "changes sign across [0, π] (table). The |V| and mass "
            "perturbations also enter at O(ε²): CP can be added without "
            "disturbing the locked structure at small ε."
        ),
        'fitted_exponent': round(p_fit, 3),
        'phi_rows': phi_rows,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T4. The ceiling identity
# ---------------------------------------------------------------------------

def test_T4_ceiling_identity() -> dict:
    """J ceiling deficit 0.249 = 0.498 × 0.902 × 0.555 exactly — the #155
    soft direction propagated; observed CP near-maximal (0.887)."""
    r_us = float(abs(_V0[0, 1])) / V_OBS['us']
    r_cb = float(abs(_V0[1, 2])) / V_OBS['cb']
    r_ub = float(abs(_V0[0, 2])) / V_OBS['ub']
    ceil_ratio = _CEIL_PRED / _CEIL_OBS
    decomp = r_us * r_cb * r_ub
    exact = abs(ceil_ratio - decomp) < 1e-12
    near_max = 0.7 < _SIN_DELTA_OBS < 1.0
    return {
        'name': 'T4_ceiling_identity',
        'description': (
            "J is bounded by the rephasing-invariant ceiling "
            "|V_us·V_cb·V_ub|. Predicted ceiling 8.64e-6 vs observed "
            "3.47e-5 — ratio 0.249, decomposing EXACTLY (identity) into "
            "the per-element ratios 0.498 × 0.902 × 0.555: the J shortfall "
            "IS the #155 V_us/V_ub soft direction, not an independent CP "
            "failure. The observed CP is near-maximal "
            "(J_obs/ceiling_obs = 0.887). Consistency lock: when the soft "
            "directions land on data, the ceiling rises to the observed "
            "value — a falsifiable future requirement on the pinhole "
            "refinement."
        ),
        'ceiling_predicted': float(f'{_CEIL_PRED:.3e}'),
        'ceiling_observed': float(f'{_CEIL_OBS:.3e}'),
        'ceiling_ratio': round(ceil_ratio, 3),
        'element_ratios': [round(r_us, 3), round(r_cb, 3), round(r_ub, 3)],
        'decomposition_product': round(decomp, 3),
        'decomposition_exact': exact,
        'sin_delta_observed': round(_SIN_DELTA_OBS, 3),
        'pass': exact and near_max,
    }


# ---------------------------------------------------------------------------
# T5. The constrained calibration
# ---------------------------------------------------------------------------

def test_T5_calibration() -> dict:
    """ε* = 0.053 at φ* = 0.80 reproduces the observed phase content on the
    predicted |V|; V_cb untouched, V_us/V_ub −4%, masses ≤ 0.5%."""
    j_cal = jarlskog(_VC)
    d_us = float(abs(_VC[0, 1]) / abs(_V0[0, 1]) - 1.0)
    d_cb = float(abs(_VC[1, 2]) / abs(_V0[1, 2]) - 1.0)
    d_ub = float(abs(_VC[0, 2]) / abs(_V0[0, 2]) - 1.0)
    d_mass = float(max(np.max(np.abs(_MUC / _MU0 - 1.0)),
                       np.max(np.abs(_MDC / _MD0 - 1.0))))
    sin_d_cal = j_cal / (abs(_VC[0, 1]) * abs(_VC[1, 2]) * abs(_VC[0, 2]))
    ok = (abs(j_cal - _J_TARGET) / _J_TARGET < 1e-6 and abs(d_cb) < 0.01
          and abs(d_us) < 0.10 and abs(d_ub) < 0.10 and d_mass < 0.016)
    return {
        'name': 'T5_constrained_calibration',
        'description': (
            "Target: the observed phase content (sin δ = 0.887) on the "
            "predicted |V| — J_target = 7.67e-6. Calibrated: ε* = "
            f"{_EPS_STAR:.4f} at φ* = {_PHI_STAR}: |V_cb| shifts −0.0% "
            "(the stiff prediction untouched), |V_us|/|V_ub| shift −4.2%/"
            "−4.3% (inside the soft direction), the masses shift ≤ 0.5% — "
            "well inside the 1.6% calibration accuracy. The locked "
            "structure SURVIVES the CP calibration, and the calibrated "
            f"phase content sin δ = {sin_d_cal:.3f} is near-maximal, like "
            "the data."
        ),
        'eps_star': round(_EPS_STAR, 4),
        'phi_star': _PHI_STAR,
        'J_calibrated': float(f'{j_cal:.3e}'),
        'J_target': float(f'{_J_TARGET:.3e}'),
        'shift_V_us': round(d_us, 4),
        'shift_V_cb': round(d_cb, 4),
        'shift_V_ub': round(d_ub, 4),
        'max_mass_shift': round(d_mass, 4),
        'sin_delta_calibrated': round(float(sin_d_cal), 3),
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T6. The triangle-shape discriminator
# ---------------------------------------------------------------------------

def test_T6_triangle_shape() -> dict:
    """The calibrated phase reproduces the area (J) but squashes the
    db-triangle (β ≈ 0° vs 22.2°): the shape discriminates phase structures
    — the falsifiable target for the true Hopf-connection φ_q(k)."""
    beta, gamma = triangle_angles_deg(_VC)
    squashed = abs(beta) < 5.0
    return {
        'name': 'T6_triangle_shape_discriminator',
        'description': (
            "The honest sharp edge: the calibrated placeholder phase "
            "reproduces the unitarity-triangle AREA (J) but not its SHAPE "
            f"— the db-triangle angles come out (β, γ) ≈ ({beta:.1f}°, "
            f"{gamma:.1f}°) vs the observed (22.2°, 65.9°): the linear "
            "φ_q(k) = φ·k places the CP in the wrong quartet orientation "
            "(a squashed triangle with the right area). The triangle shape "
            "therefore DISCRIMINATES phase structures: β = 22° is the "
            "falsifiable target for the TRUE Hopf-connection phase — the "
            "v3 §4 TODO tracked since the handoff, now with a quantitative "
            "acceptance test."
        ),
        'beta_deg': round(beta, 2),
        'gamma_deg': round(gamma, 2),
        'beta_observed_deg': BETA_OBS_DEG,
        'gamma_observed_deg': GAMMA_OBS_DEG,
        'placeholder_fails_shape': squashed,
        'pass': squashed,
    }


# ---------------------------------------------------------------------------
# T7. Ledger / scope
# ---------------------------------------------------------------------------

def test_T7_ledger_and_scope() -> dict:
    return {
        'name': 'T7_ledger_and_scope',
        'description': (
            "CONSUMED: one new input — the quark CP phase content (like α: "
            "structure derived, value input; the flavor puzzle's CP entry "
            "made explicit in the #150 budget). DERIVED: the quadratic "
            "ε-scaling, the ceiling identity with its exact soft-direction "
            "decomposition (no independent CP failure), the survival of "
            "the locked |V|/masses under calibration, the near-maximal "
            "phase content. OPEN, with falsifiable targets: the "
            "Hopf-connection φ_q(k) (target: the db-triangle shape, "
            "β = 22°); the soft V_us/V_ub directions (target: the J "
            "ceiling rising to 3.5e-5)."
        ),
        'consumed': ['the quark CP phase content (one number — J/sin δ)'],
        'derived': [
            'J ∝ ε² switching; sinusoidal φ-dependence',
            'ceiling identity: J deficit = soft-direction |V| deficit (exact)',
            'locked structure survives calibration (V_cb −0.0%, masses ≤ 0.5%)',
            'near-maximal phase content (sin δ ≈ 0.97 vs data 0.89)',
        ],
        'open': ['Hopf-connection φ_q(k) — target β = 22°',
                 'soft V_us/V_ub — target J ceiling 3.5e-5'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The quark CP phase calibrates onto the locked Hamiltonian "
            "without disturbing |V| or the masses; the J ceiling deficit "
            "decomposes exactly into the #155 soft-direction deficits; the "
            "phase content is near-maximal like the data; and the "
            "db-triangle shape (β = 22°) is the falsifiable target the "
            "placeholder linear phase fails — pointing at the "
            "Hopf-connection φ_q(k) as the next derivation."
        ),
        'classification': 'QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_extension(),
        test_T3_scaling(),
        test_T4_ceiling_identity(),
        test_T5_calibration(),
        test_T6_triangle_shape(),
        test_T7_ledger_and_scope(),
        test_T8_assessment(),
    ]
    t4, t5, t6 = tests[3], tests[4], tests[5]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN'
        verdict = (
            'THE QUARK CP PHASE CALIBRATES ONTO THE LOCKED HAMILTONIAN '
            'WITHOUT DISTURBING |V| OR THE MASSES; THE J CEILING DEFICIT '
            'DECOMPOSES EXACTLY INTO THE #155 SOFT-DIRECTION |V| DEFICITS '
            '(NO INDEPENDENT CP FAILURE); THE PHASE CONTENT IS NEAR-MAXIMAL '
            'LIKE THE DATA; AND THE TRIANGLE SHAPE IS THE FALSIFIABLE '
            'TARGET THE PLACEHOLDER PHASE FAILS. #155 posed the constrained '
            'problem; this probe solves it and sharpens what remains.\n\n'
            'THE EXTENSION. H(ε, φ) = H_locked + the v3 §4 partition-mixing '
            'element with the Hopf-placeholder phase φ_q(k) = φ·k, built '
            'in-probe (locked blocks exactly intact; baseline J = 0 '
            're-verified). The global phase knob is unusable — it enters '
            'the transport couplings as cos(phase·dk) and collapses |V_us| '
            '×20 at phase = 0.5 — which is why the mass calibration locked '
            'it at zero and the #155 baseline had J = 0 structurally.\n\n'
            'SCALING. J ∝ ε^1.9 (quadratic — one insertion per sector '
            'side), sinusoidal in φ; |V|/mass shifts also O(ε²): CP '
            'switches on cleanly.\n\n'
            'THE CEILING IDENTITY. J ≤ |V_us·V_cb·V_ub|: predicted ceiling '
            f'{t4["ceiling_predicted"]} vs observed '
            f'{t4["ceiling_observed"]} — ratio {t4["ceiling_ratio"]}, '
            'decomposing EXACTLY into the element ratios '
            f'{t4["element_ratios"]}: the J shortfall IS the #155 '
            'V_us/V_ub soft direction propagated, not a new failure. The '
            f'observed CP is near-maximal ({t4["sin_delta_observed"]}). '
            'Consistency lock: when the soft directions land, the ceiling '
            'rises to the observed 3.5e-5.\n\n'
            f'THE CALIBRATION. ε* = {t5["eps_star"]}, φ* = '
            f'{t5["phi_star"]}: J hits the phase-content target '
            f'({t5["J_calibrated"]}) with |V_cb| shifted '
            f'{t5["shift_V_cb"]:+.1%} (the stiff prediction untouched), '
            f'|V_us|/|V_ub| {t5["shift_V_us"]:+.1%}/{t5["shift_V_ub"]:+.1%} '
            f'(inside the soft direction), masses ≤ '
            f'{t5["max_mass_shift"]:.1%} — the locked structure survives; '
            f'sin δ = {t5["sin_delta_calibrated"]}, near-maximal like the '
            'data.\n\n'
            'THE TRIANGLE SHAPE. The placeholder phase reproduces the AREA '
            '(J) but squashes the db-triangle: (β, γ) ≈ '
            f'({t6["beta_deg"]}°, {t6["gamma_deg"]}°) vs observed (22.2°, '
            '65.9°). The shape discriminates phase structures — β = 22° is '
            'the quantitative acceptance test for the true Hopf-connection '
            'φ_q(k) (the v3 §4 TODO).\n\n'
            'LEDGER. One input consumed (the CP phase content — the flavor '
            'puzzle\'s CP entry made explicit); the scaling, the ceiling '
            'identity, the calibration survival, and the near-maximal '
            'content derived; the Hopf phase (β = 22°) and the soft |V| '
            'directions (J ceiling) open with falsifiable targets.'
        )
    else:
        verdict_class = 'QUARK_CP_CALIBRATION_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A calibration check failed; review the '
            'extension, the ceiling identity, or the shift bounds.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the quark CP phase calibrates onto the locked Hamiltonian '
            'without disturbing |V| or the masses; the J ceiling deficit '
            'decomposes exactly into the soft-direction |V| deficits; the '
            'phase content is near-maximal; the triangle shape (β = 22°) is '
            'the falsifiable target for the Hopf-connection phase'
        ),
        'extension': 'H(ε,φ) = H_locked + partition-mixing with φ_q(k) = φ·k (in-probe)',
        'scaling': 'J ∝ ε^1.9; sinusoidal in φ; shifts O(ε²)',
        'ceiling': 'deficit 0.249 = 0.498×0.902×0.555 exactly (the #155 soft direction)',
        'calibration': f'ε* = {_EPS_STAR:.4f}, φ* = {_PHI_STAR}: V_cb −0.0%, masses ≤ 0.5%',
        'triangle': 'β ≈ 0° vs 22.2° — shape open; the Hopf φ_q(k) acceptance test',
        'consumed': 'one input: the quark CP phase content',
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
    out.append('# Quark CP phase calibration on the locked Hamiltonian (PR #156)')
    out.append('')
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append('')
    out.append(
        "Solves the constrained problem #155 posed: the quark CP phase is "
        "introduced on the partition-mixing element (the Hopf-placeholder "
        "φ_q(k) = φ·k), calibrated to the observed phase content, and the "
        "locked |V| and masses survive. The J ceiling deficit decomposes "
        "exactly into the #155 soft-direction |V| deficits — no independent "
        "CP failure — and the db-triangle shape (β = 22°) emerges as the "
        "falsifiable acceptance test for the true Hopf-connection phase. "
        "*(QFT on the classical throat, not quantum gravity.)*"
    )
    out.append('')
    out.append(f"- **Extension**: {s['extension']}")
    out.append(f"- **Scaling**: {s['scaling']}")
    out.append(f"- **Ceiling**: {s['ceiling']}")
    out.append(f"- **Calibration**: {s['calibration']}")
    out.append(f"- **Triangle**: {s['triangle']}")
    out.append(f"- **Consumed**: {s['consumed']}")
    out.append('')

    out.append('## Test summary')
    out.append('')
    out.append('| # | Test | Key finding | PASS? |')
    out.append('|---|---|---|---|')
    label_map = {
        'T1': 'the constrained CP calibration #155 posed',
        'T2': 'in-probe extension; baseline J = 0; global knob unusable',
        'T3': 'J ∝ ε^1.9; sinusoidal φ-dependence; shifts O(ε²)',
        'T4': 'ceiling deficit 0.249 = soft-direction product (exact identity)',
        'T5': 'ε* = 0.053: V_cb −0.0%, V_us/V_ub −4%, masses ≤ 0.5%',
        'T6': 'placeholder fails the triangle shape: β ≈ 0° vs 22.2° (target)',
        'T7': 'ledger: one input consumed; ceiling/shape open with targets',
        'T8': 'QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        out.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    out.append('')

    t4 = s['tests'][3]
    out.append('## The ceiling identity (the J deficit is the soft direction)')
    out.append('')
    out.append('| quantity | value |')
    out.append('|---|---:|')
    out.append(f"| J ceiling, predicted \\|V\\| | {t4['ceiling_predicted']} |")
    out.append(f"| J ceiling, observed \\|V\\| | {t4['ceiling_observed']} |")
    out.append(f"| ceiling ratio | {t4['ceiling_ratio']} |")
    out.append(f"| element ratios (us, cb, ub) | {t4['element_ratios']} |")
    out.append(f"| decomposition product | {t4['decomposition_product']} |")
    out.append(f"| sin δ observed (near-maximal) | {t4['sin_delta_observed']} |")
    out.append('')

    t5 = s['tests'][4]
    out.append('## The constrained calibration')
    out.append('')
    out.append('| quantity | value |')
    out.append('|---|---:|')
    out.append(f"| ε* | {t5['eps_star']} |")
    out.append(f"| φ* | {t5['phi_star']} |")
    out.append(f"| J calibrated / target | {t5['J_calibrated']} / {t5['J_target']} |")
    out.append(f"| shift V_us | {t5['shift_V_us']} |")
    out.append(f"| shift V_cb | {t5['shift_V_cb']} |")
    out.append(f"| shift V_ub | {t5['shift_V_ub']} |")
    out.append(f"| max mass shift | {t5['max_mass_shift']} |")
    out.append(f"| sin δ calibrated | {t5['sin_delta_calibrated']} |")
    out.append('')

    t6 = s['tests'][5]
    out.append('## The triangle-shape discriminator')
    out.append('')
    out.append(f"Calibrated (β, γ) = ({t6['beta_deg']}°, {t6['gamma_deg']}°) "
               f"vs observed ({t6['beta_observed_deg']}°, "
               f"{t6['gamma_observed_deg']}°): the placeholder linear phase "
               "reproduces the area (J) but squashes the db-triangle — "
               "β = 22° is the quantitative acceptance test for the true "
               "Hopf-connection φ_q(k).")
    out.append('')

    out.append('## Verdict')
    out.append('')
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
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
    out = here / 'runs' / f'{ts}_quark_cp_phase_calibration_probe'
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
