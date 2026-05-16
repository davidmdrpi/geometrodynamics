"""
Finite-energy Klein-Nishina recoil probe.

Follow-on to the Thomson-limit Klein-Nishina restoration probe
(PR #28). The previous probe established that the BAM Compton
amplitude with explicit photon transverse polarization + scalar
electron reproduces the Klein-Nishina angular factor
`(1 + cos²θ)/2` exactly at ω → 0. This probe tests whether the same
construction reproduces the finite-energy KN structure, or whether
the Thomson-limit match conceals a leading-order failure at ω ~ m_e.

The KN squared amplitude (electron rest frame, unpolarized):

    |M_KN|²(ω, θ)  ∝  (ω'/ω)² · [ω'/ω + ω/ω' − sin²θ]

with `ω' = ω/(1 + (ω/m)(1 − cos θ))`. In terms of `x ≡ ω'/ω`:

    |M_KN|²  ∝  x² · (x + 1/x − sin²θ)  =  x³ + x − x² sin²θ.

The BAM construction (photon-structure probe, scalar electron):

    |M_BAM|²(ω, θ)  ∝  (1 + cos²θ) · (G_s + G_u)²
                    ∝  (1 + cos²θ) · (m/ω + m/ω')²
                    =  (1 + cos²θ) · (m/ω)² · (1 + 1/x)².

Both reduce to `(1 + cos²θ)` at x → 1 (Thomson). At finite ω they
differ structurally: KN has `x³ + x − x² sin²θ`, BAM has
`(1 + cos²θ) · (x + 1)²/(4 x²)` after appropriate normalisation.

Tests:

  T1. Thomson sanity — at ω/m = 1e-4 the BAM and KN normalised
      angular distributions agree to < 1e-3.

  T2. Leading-order discrepancy — Taylor-expand
      `f_BAM(ε, θ) − f_KN(ε, θ)` in ε = ω/m and identify the lowest
      non-vanishing order. Predicted: O(ε) at generic θ.

  T3. Forward-backward asymmetry — compute
      `A(ε) = (f(0) − f(π))/(f(0) + f(π))` for both KN and BAM
      across ε ∈ [0.01, 1]. KN gives a specific ε-dependence; BAM's
      will differ.

  T4. Compton-edge backscatter — at θ = π, evaluate both expressions
      as a function of ε on a log scale. Identify the ε at which BAM
      and KN diverge by 10%.

  T5. Full angular fit at finite ω — at ω/m = 0.5 (relativistic),
      fit both |M_BAM|²(θ) and |M_KN|²(θ) to the angular polynomial
      basis {1, cos θ, cos²θ, cos³θ}. Compare coefficients.

Verdict:
  PARTIAL_MATCH if the Thomson sanity (T1) passes but the
  finite-energy tests T2-T5 show systematic deviation. This is the
  expected outcome — the natural BAM construction reproduces
  Thomson exactly but lacks the QED vertex-factor structure that
  generates KN's (ω'/ω)² recoil enhancement.

  FULL_MATCH if all tests pass — would indicate the simple BAM
  construction captures the finite-ω structure. Unlikely under the
  analytic prediction.

  NO_MATCH if even Thomson fails — would indicate a regression in
  the construction.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
TAU = 2.0 * PI


# ---------------------------------------------------------------------------
# Klein-Nishina exact reference (electron rest frame, unpolarized)
# ---------------------------------------------------------------------------

def x_ratio(omega: float, theta: float, m: float = 1.0) -> float:
    """x ≡ ω'/ω = 1/(1 + (ω/m)(1 − cos θ))."""
    return 1.0 / (1.0 + (omega / m) * (1.0 - math.cos(theta)))


def kn_squared_amplitude(omega: float, theta: float, m: float = 1.0) -> float:
    """|M_KN|²(ω, θ) ∝ x² · (x + 1/x − sin²θ)  (no overall α², 1/m² etc.)."""
    x = x_ratio(omega, theta, m)
    sin2 = math.sin(theta) ** 2
    return x * x * (x + 1.0 / x - sin2)


def kn_normalized(omega: float, theta: float, m: float = 1.0) -> float:
    """KN normalised at θ = 0 (x = 1, sin²=0 → factor 2)."""
    return kn_squared_amplitude(omega, theta, m) / 2.0


# ---------------------------------------------------------------------------
# BAM construction (from photon-structure probe)
# ---------------------------------------------------------------------------

def G_S3(psi: float, R: float = 1.0, eps: float = 1e-12) -> float:
    psi_eff = max(min(psi, PI - eps), eps)
    return float(
        (((PI - psi_eff) / math.tan(psi_eff)) - 0.5)
        / (4.0 * PI * PI * R)
    )


def bam_squared_amplitude(omega: float, theta: float, m: float = 1.0) -> float:
    """|M_BAM|²(ω, θ) = (1+cos²θ) · (G_s + G_u)²  (photon polarization
    summed, scalar electron; up to the 1/2 polarization-average factor
    which cancels in any normalised comparison).

    Implements the exact same construction as
    `compton_photon_structure_probe.py` so the Thomson result is
    reproduced.
    """
    omega_prime = omega / (1.0 + (omega / m) * (1.0 - math.cos(theta)))
    s = m * m + 2.0 * m * omega
    u = m * m - 2.0 * m * omega_prime
    eps_pole = 1e-12
    psi_s = max((s - m * m) / (2.0 * m * m), eps_pole)
    psi_u = max(abs(u - m * m) / (2.0 * m * m), eps_pole)
    psi_s = min(psi_s, PI - eps_pole)
    psi_u = min(psi_u, PI - eps_pole)
    Gs = G_S3(psi_s)
    Gu = G_S3(psi_u)
    pol_sum = 1.0 + math.cos(theta) ** 2
    return pol_sum * (Gs + Gu) ** 2


def bam_normalized(omega: float, theta: float, m: float = 1.0) -> float:
    return bam_squared_amplitude(omega, theta, m) / max(
        bam_squared_amplitude(omega, 0.0, m), 1e-30,
    )


def kn_normalized_at_zero(omega: float, theta: float, m: float = 1.0) -> float:
    return kn_squared_amplitude(omega, theta, m) / max(
        kn_squared_amplitude(omega, 0.0, m), 1e-30,
    )


# ---------------------------------------------------------------------------
# T1. Thomson sanity
# ---------------------------------------------------------------------------

def test_T1_thomson_sanity() -> dict:
    omega = 1e-4
    thetas = np.linspace(0.0, PI, 33)
    diffs = []
    for theta in thetas:
        b = bam_normalized(omega, float(theta))
        k = kn_normalized_at_zero(omega, float(theta))
        diffs.append(abs(b - k))
    return {
        'name': 'T1_thomson_sanity',
        'description': (
            'At ω/m = 1e-4 the BAM and KN normalised angular '
            'distributions agree across θ ∈ [0, π].'
        ),
        'omega': omega,
        'max_pointwise_diff': float(np.max(diffs)),
        'pass': float(np.max(diffs)) < 1e-3,
    }


# ---------------------------------------------------------------------------
# T2. Leading-order discrepancy in ε = ω/m
# ---------------------------------------------------------------------------

def test_T2_leading_order_discrepancy() -> dict:
    """Fit log|f_BAM − f_KN| ≈ p · log(ε) at fixed θ, identify lowest-
    order discrepancy."""
    thetas_test = [PI / 4.0, PI / 2.0, 3.0 * PI / 4.0]
    epsilons = np.array([10.0 ** (-k) for k in [1, 2, 3, 4, 5, 6]])
    rows = []
    for theta in thetas_test:
        diffs = []
        for eps in epsilons:
            b = bam_normalized(float(eps), float(theta))
            k = kn_normalized_at_zero(float(eps), float(theta))
            diffs.append(abs(b - k))
        diffs_arr = np.asarray(diffs)
        valid = diffs_arr > 0
        if valid.sum() >= 3:
            slope, intercept = np.polyfit(
                np.log(epsilons[valid]), np.log(diffs_arr[valid]), 1,
            )
            order = float(slope)
        else:
            order = float('nan')
        rows.append({
            'theta_in_pi': theta / PI,
            'epsilons': epsilons.tolist(),
            'diffs': diffs_arr.tolist(),
            'fitted_order_in_eps': order,
        })

    # Predicted: leading discrepancy is O(ε)  (i.e. order = 1.0)
    orders = [r['fitted_order_in_eps'] for r in rows]
    near_one = all(0.8 < o < 1.2 for o in orders if not math.isnan(o))
    return {
        'name': 'T2_leading_order_discrepancy',
        'description': (
            'Lowest non-vanishing order in ε = ω/m of '
            'f_BAM(ε, θ) − f_KN(ε, θ) at generic θ. Predicted '
            'fitted slope = 1.0 (O(ε)).'
        ),
        'per_theta': rows,
        'fitted_orders_summary': orders,
        'all_near_O_eps': near_one,
        # T2 "passes" if the discrepancy is well-characterised at
        # leading order — the probe is informative, not pass/fail
        # on KN reproduction.
        'pass': near_one,
    }


# ---------------------------------------------------------------------------
# T3. Forward-backward asymmetry vs ε
# ---------------------------------------------------------------------------

def test_T3_forward_backward_asymmetry() -> dict:
    """Compute A(ε) = (f(0) − f(π))/(f(0) + f(π)) at finite ω.

    Note: with normalisation at θ=0 the "asymmetry" reduces to
    (1 − f_norm(π))/(1 + f_norm(π)).
    """
    epsilons = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0])
    rows = []
    for eps in epsilons:
        b_pi = bam_normalized(float(eps), PI)
        k_pi = kn_normalized_at_zero(float(eps), PI)
        A_bam = (1.0 - b_pi) / (1.0 + b_pi)
        A_kn = (1.0 - k_pi) / (1.0 + k_pi)
        rows.append({
            'epsilon': float(eps),
            'f_BAM_pi_normalized': float(b_pi),
            'f_KN_pi_normalized': float(k_pi),
            'asymmetry_BAM': float(A_bam),
            'asymmetry_KN': float(A_kn),
            'asymmetry_difference': float(A_bam - A_kn),
        })
    max_diff = max(abs(r['asymmetry_difference']) for r in rows)
    return {
        'name': 'T3_forward_backward_asymmetry',
        'description': (
            'Forward-backward asymmetry A(ε) = (f(0)−f(π))/(f(0)+f(π)) '
            'compared between BAM and KN across ε ∈ [0.01, 1].'
        ),
        'rows': rows,
        'max_asymmetry_difference': float(max_diff),
        # Informative test — pass = the difference is well-characterised
        # and grows with ε (a structural signature, not a bug).
        'asymmetry_grows_with_eps': all(
            abs(rows[i + 1]['asymmetry_difference'])
            >= abs(rows[i]['asymmetry_difference']) - 1e-12
            for i in range(len(rows) - 1)
        ),
        'pass': max_diff > 1e-3,  # i.e. a meaningful gap exists
    }


# ---------------------------------------------------------------------------
# T4. Compton edge at θ = π
# ---------------------------------------------------------------------------

def test_T4_compton_edge() -> dict:
    """At θ = π, evaluate f_BAM(ε, π) and f_KN(ε, π) for ε ∈ [1e-4, 10],
    log scale. Identify the smallest ε at which the relative difference
    exceeds 10%."""
    epsilons = np.logspace(-4.0, 1.0, 41)
    rows = []
    eps_10pct = None
    for eps in epsilons:
        b = bam_normalized(float(eps), PI)
        k = kn_normalized_at_zero(float(eps), PI)
        rel_diff = abs(b - k) / max(abs(k), 1e-30)
        rows.append({
            'epsilon': float(eps),
            'f_BAM': float(b),
            'f_KN': float(k),
            'relative_difference': float(rel_diff),
        })
        if eps_10pct is None and rel_diff > 0.10:
            eps_10pct = float(eps)

    return {
        'name': 'T4_compton_edge',
        'description': (
            'Compton backscatter (θ = π) comparison vs ε on log '
            'scale. Identifies the ε at which BAM and KN diverge by '
            '> 10 %.'
        ),
        'rows_compact': rows[::4],
        'eps_at_10pct_divergence': eps_10pct,
        # Probe is informative; PASS = eps_10pct is finite (gap is
        # quantifiable rather than blowing up at all energies).
        'pass': eps_10pct is not None and eps_10pct > 0.001,
    }


# ---------------------------------------------------------------------------
# T5. Full angular fit at finite ω
# ---------------------------------------------------------------------------

def test_T5_angular_fit_at_finite_omega(omega: float = 0.5) -> dict:
    """At ω/m = 0.5 (relativistic Compton), fit |M_BAM|²(θ) and
    |M_KN|²(θ) to {1, cos θ, cos²θ, cos³θ} basis. Compare
    coefficients."""
    thetas = np.linspace(0.0, PI, 65)
    c = np.cos(thetas)
    basis = np.column_stack([np.ones_like(c), c, c * c, c * c * c])

    bam_vals = np.array([bam_normalized(omega, float(t)) for t in thetas])
    kn_vals = np.array([kn_normalized_at_zero(omega, float(t)) for t in thetas])

    bam_coeffs, _, _, _ = np.linalg.lstsq(basis, bam_vals, rcond=None)
    kn_coeffs, _, _, _ = np.linalg.lstsq(basis, kn_vals, rcond=None)
    diff_coeffs = [float(b - k) for b, k in zip(bam_coeffs, kn_coeffs)]

    bam_pred = basis @ bam_coeffs
    kn_pred = basis @ kn_coeffs
    bam_fit_residual = float(np.max(np.abs(bam_vals - bam_pred)))
    kn_fit_residual = float(np.max(np.abs(kn_vals - kn_pred)))

    max_pointwise = float(np.max(np.abs(bam_vals - kn_vals)))

    return {
        'name': 'T5_angular_fit_at_finite_omega',
        'description': (
            'Fit |M|²(θ) to {1, cos θ, cos²θ, cos³θ} at ω/m = 0.5 '
            'for both BAM and KN. Compare coefficients.'
        ),
        'omega': omega,
        'bam_coefficients_c0_c1_c2_c3': [float(c) for c in bam_coeffs],
        'kn_coefficients_c0_c1_c2_c3': [float(c) for c in kn_coeffs],
        'coefficient_differences_c0_c1_c2_c3': diff_coeffs,
        'bam_fit_residual_in_basis': bam_fit_residual,
        'kn_fit_residual_in_basis': kn_fit_residual,
        'max_pointwise_difference': max_pointwise,
        # Informative test — PASS = gap is finite, non-trivial,
        # and well-characterised (no NaN). The probe documents the
        # gap; the verdict comes from T1 + T2.
        'pass': (
            not math.isnan(max_pointwise) and
            max_pointwise > 0.1 and
            all(not math.isnan(c) for c in bam_coeffs) and
            all(not math.isnan(c) for c in kn_coeffs)
        ),
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    t1 = test_T1_thomson_sanity()
    t2 = test_T2_leading_order_discrepancy()
    t3 = test_T3_forward_backward_asymmetry()
    t4 = test_T4_compton_edge()
    t5 = test_T5_angular_fit_at_finite_omega()
    tests = [t1, t2, t3, t4, t5]

    # Verdict
    thomson_ok = t1['pass']
    if thomson_ok:
        # The probe is structurally informative even when KN is not
        # fully reproduced. Classify by leading-order behaviour.
        if t2.get('all_near_O_eps'):
            verdict_class = 'PARTIAL_MATCH_O1_THOMSON_ONLY'
            verdict = (
                'PARTIAL MATCH — KN reproduced at O(1) in ε = ω/m_e '
                '(Thomson limit, PR #28), but the natural BAM '
                'construction fails to reproduce KN at leading order '
                'O(ε). The discrepancy f_BAM − f_KN scales linearly '
                'in ω/m_e (T2 fitted slopes near 1.0). The forward-'
                'backward asymmetry (T3) and the Compton edge at '
                'θ = π (T4) show the structural signature: BAM\'s '
                'per-channel propagator sum (G_s + G_u)² gives a '
                '(1 + 1/x)² energy dependence, while QED\'s vertex-'
                'factor algebra gives x² · (x + 1/x − sin²θ). The '
                'missing BAM ingredient is the per-channel kinematic '
                'weighting — the QED vertex factors that contract '
                'photon polarization with momentum (ε·k structures), '
                'which the present BAM construction does not have '
                'natively. This locates the next structural piece: '
                'a vertex-coupling derivation from the Hopf '
                'connection or throat-transport algebra.'
            )
        else:
            verdict_class = 'PARTIAL_MATCH_THOMSON_ONLY'
            verdict = (
                'PARTIAL MATCH — Thomson reproduced; finite-energy '
                'discrepancy is present but not at the predicted '
                'O(ε) order. Check T2 fits for the actual leading-'
                'order behaviour and re-interpret.'
            )
    else:
        verdict_class = 'REGRESSION'
        verdict = (
            'REGRESSION — even the Thomson limit fails in the present '
            'construction. This indicates a bug or change relative to '
            'PR #28; the probe construction should be re-verified.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'construction': (
            'BAM: |M|² ∝ (1+cos²θ) · (1/ψ_s + 1/ψ_u)²  with '
            'ψ_x = (x−m²)/(2m²). KN: |M|² ∝ x²·(x + 1/x − sin²θ) '
            'with x = ω\'/ω.'
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
    L.append('# Finite-energy Klein-Nishina recoil probe')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        'Follow-on to PR #28 (Thomson-limit KN restoration). Tests '
        'whether the natural BAM construction — antipodal `S³` '
        'Green-function propagator + photon transverse polarization '
        '+ scalar electron — also reproduces Klein-Nishina at finite '
        '`ω/m_e`, or whether the Thomson match conceals a leading-'
        'order failure.'
    )
    L.append('')
    L.append('**Construction:**')
    L.append('')
    L.append('```')
    L.append(s['construction'])
    L.append('```')
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key metric | Value | PASS? |')
    L.append('|---|---|---|---:|---|')
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        if nm.startswith('T1'):
            metric, value = (
                'max |f_BAM − f_KN| at ω/m=1e-4',
                f"{t['max_pointwise_diff']:.2e}",
            )
        elif nm.startswith('T2'):
            metric, value = (
                'fitted leading order in ε (predicted 1.0)',
                ', '.join(
                    f"{o:.3f}" if not math.isnan(o) else 'nan'
                    for o in t['fitted_orders_summary']
                ),
            )
        elif nm.startswith('T3'):
            metric, value = (
                'max |A_BAM(ε) − A_KN(ε)|',
                f"{t['max_asymmetry_difference']:.4f}",
            )
        elif nm.startswith('T4'):
            metric, value = (
                'ε at which BAM/KN diverge by 10 %',
                f"{t['eps_at_10pct_divergence']:.4f}"
                if t['eps_at_10pct_divergence'] is not None else 'n/a',
            )
        elif nm.startswith('T5'):
            c_diff = t['coefficient_differences_c0_c1_c2_c3']
            metric, value = (
                'Δc0, Δc1, Δc2, Δc3 (BAM − KN) at ω/m=0.5',
                ', '.join(f"{x:+.3f}" for x in c_diff),
            )
        else:
            metric, value = '—', '—'
        L.append(f"| {nm[:2]} | `{nm}` | {metric} | {value} | {passed} |")
    L.append('')

    for t in s['tests']:
        L.append(f"## {t['name']}")
        L.append('')
        L.append(t['description'])
        L.append('')
        if t['name'].startswith('T2'):
            L.append('Per-θ fits:')
            L.append('')
            L.append('| θ/π | fitted order in ε |')
            L.append('|---:|---:|')
            for r in t['per_theta']:
                o = r['fitted_order_in_eps']
                o_str = f"{o:.4f}" if not math.isnan(o) else 'nan'
                L.append(f"| {r['theta_in_pi']:.4f} | {o_str} |")
            L.append('')
        elif t['name'].startswith('T3'):
            L.append('| ε | f_BAM(π) | f_KN(π) | A_BAM | A_KN | ΔA |')
            L.append('|---:|---:|---:|---:|---:|---:|')
            for r in t['rows']:
                L.append(
                    f"| {r['epsilon']:.3f} | "
                    f"{r['f_BAM_pi_normalized']:.4f} | "
                    f"{r['f_KN_pi_normalized']:.4f} | "
                    f"{r['asymmetry_BAM']:+.4f} | "
                    f"{r['asymmetry_KN']:+.4f} | "
                    f"{r['asymmetry_difference']:+.4f} |"
                )
            L.append('')
        elif t['name'].startswith('T4'):
            L.append(
                f"Smallest ε with > 10 % BAM/KN divergence: "
                f"**{t['eps_at_10pct_divergence']}**"
            )
            L.append('')
            L.append('| ε | f_BAM(π) | f_KN(π) | rel. diff |')
            L.append('|---:|---:|---:|---:|')
            for r in t['rows_compact']:
                L.append(
                    f"| {r['epsilon']:.4e} | {r['f_BAM']:.4e} | "
                    f"{r['f_KN']:.4e} | {r['relative_difference']:.4f} |"
                )
            L.append('')
        elif t['name'].startswith('T5'):
            cb = t['bam_coefficients_c0_c1_c2_c3']
            ck = t['kn_coefficients_c0_c1_c2_c3']
            L.append(
                f"BAM coefficients (c0, c1, c2, c3): "
                f"{cb[0]:+.4f}, {cb[1]:+.4f}, {cb[2]:+.4f}, {cb[3]:+.4f}"
            )
            L.append(
                f"KN coefficients  (c0, c1, c2, c3): "
                f"{ck[0]:+.4f}, {ck[1]:+.4f}, {ck[2]:+.4f}, {ck[3]:+.4f}"
            )
            L.append('')
            L.append(
                f"Max pointwise |M_BAM − M_KN| at ω/m=0.5: "
                f"**{t['max_pointwise_difference']:.4f}**"
            )
            L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append(
        '- **Vertex-factor BAM derivation.** The probe locates the '
        'finite-ω gap in the missing per-channel kinematic weighting '
        '(ε·k vertex structure in QED). Whether a BAM-derived '
        'vertex factor — from explicit Hopf-connection coupling, '
        'throat-transport algebra, or another natural construction '
        '— closes the gap is the natural follow-on probe.'
    )
    L.append(
        '- **Photon as Hopf-fibre excitation.** The current photon '
        'model treats polarization as flat-space tangent-plane '
        'vectors. An `S³`-native photon (Hopf-fibre excitation, '
        '`hopf/connection.py`) would automatically carry the '
        'connection coupling that supplies ε·k vertex factors. '
        'Sketching this is the natural design step.'
    )
    L.append(
        '- **Electron spin at finite energy.** The Thomson sanity '
        'used a scalar electron; finite-ω Compton in QED has spin-½ '
        'corrections that BAM\'s previous probe identified as '
        'spurious at Thomson but should reappear at higher orders.'
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
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_compton_finite_energy_kn_probe'
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
