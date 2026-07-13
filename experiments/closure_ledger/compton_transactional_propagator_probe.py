"""
The transactional Compton propagator: the Feynman propagator, its two
completion orderings, and their coherent relative phase, derived from a
frozen classical bulk geometry (PR #213).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

WHAT WAS IMPORTED, UNTIL NOW
----------------------------
The Compton arc (#35, #45, #46, #211) derived the propagator's SPATIAL
structure from geometry: the pole from the antipodal S^3 Green function
(#35), the magnitude 1/q^2 from the flat limit (#45), the Lorentz tensor
from the Hopf connection (#46).  But the propagator's TIME structure -
the Feynman i*epsilon prescription, i.e. the split of G_F into the two
time-ordering completions and the assertion that they add coherently
with relative phase +1 - was imported from QED, not derived.  In the
standard theory that structure encodes "virtual particles running
forward and backward in time"; here it must come from the one place the
program allows: complete classical histories on the frozen bulk.

THE DERIVATION (Wheeler-Feynman / transactional, on the frozen bulk)
--------------------------------------------------------------------
Work per mode (frequency omega); conventions:
    (d^2/dt^2 + omega^2) G(t) = delta(t)
    G_ret(t) =  theta(t)  sin(omega t)/omega
    G_adv(t) = -theta(-t) sin(omega t)/omega
    FT: Gt(W) = int e^{+iWt} G(t) dt.

1. TIME SYMMETRY IS NOT A CHOICE.  The bulk is static (frozen): t -> -t
   is an isometry, so the elementary field of a source is the
   time-symmetric Gbar = (G_ret + G_adv)/2 - the Wheeler-Feynman
   propagator.  Neither retardation nor advancement can be preferred by
   the geometry.

2. COMPLETE ABSORPTION IS A GEOMETRIC FACT, NOT AN ASSUMPTION.  On a
   CLOSED S^3 spatial section every retarded wavefront refocuses at the
   antipode at t = pi R (#166's 1/sin(psi) focusing) and returns to the
   source at t = 2 pi R: nothing escapes, every offer meets an absorber.
   Wheeler-Feynman's "perfect absorber" hypothesis - contingent on
   distant matter in flat space - is a THEOREM of the closed geometry.

3. THE ABSORBER RESPONSE SELECTS G_F UNIQUELY.  Complete histories
   leave no free remnants: the total field must be purely
   positive-frequency (e^{-i omega t}) in the far future (every emitted
   quantum is a confirmed transaction) and purely negative-frequency in
   the far past.  Adding the general homogeneous response
   a*cos(omega t) + b*sin(omega t) to Gbar, that boundary condition is
   a 2x2 linear system with the UNIQUE solution a = i/(2 omega), b = 0:
       G_F(t) = Gbar(t) + (i/2 omega) cos(omega t)
              = (i/2 omega) e^{-i omega |t|}.
   The Feynman propagator IS the time-symmetric field plus the
   absorber response of the closed universe.

4. THE TWO COMPLETION ORDERINGS.  G_F(t) splits into theta(t) and
   theta(-t) pieces - the offer segment (excitation after the first
   vertex) and the confirmation segment (excitation before it).  Their
   half-line Fourier transforms are the two old-fashioned
   perturbation-theory energy denominators
       I+ = -(1/2w) / (Delta - w + i eps)
       I- = +(1/2w) / (Delta + w - i eps)
   each individually regulator-dependent (finite-T truncations never
   converge) and individually non-covariant (not functions of Delta^2);
   their coherent sum is exactly the covariant pole
       I+ + I- = -(w - i eps)/(w (Delta^2 - (w - i eps)^2))
              -> -1/(Delta^2 - w^2 + i eps).

5. THE RELATIVE PHASE IS FORCED.  The t -> -t isometry maps the two
   segments into each other (G_ret(t) = G_adv(-t) on the S^3 tower),
   so they enter every complete history with equal weight: relative
   phase e^{i phi} with phi = 0.  Deforming phi away from 0 destroys
   the covariant pole (Delta^2-dependence broken, pole form broken)
   and fails the Compton denominators.  The repo's transactional
   engine (retro_phase_match) already scores phase closure with its
   maximum exactly at the coherent point.

Tests:
  T1. Goal.
  T2. The transactional construction: Gbar + unique complete-absorption
      response = G_F (jumps, homogeneity, the 2x2 uniqueness, the
      closed form; the frequency-splitting identity with its epsilon
      scaling).
  T3. The two completion orderings: half-line transforms = OFPT energy
      denominators; finite-T truncations individually unstable, damped
      orderings converge; the coherent sum = the exact covariant pole;
      each ordering alone breaks Delta -> -Delta evenness at O(1).
  T4. The frozen bulk supplies both conditions: the t -> -t identity
      G_ret(t) = G_adv(-t) on the conformal tower omega_l = (l+1)/R;
      complete absorption as geometry - antipodal refocus at t = pi R
      and full return at t = 2 pi R.
  T5. The coherent relative phase: phi = 0 exact; the deform test
      e^{i phi} breaks covariance and the pole form at O(1); the
      engine's phase closure peaks at the coherent point (with the
      pi branch = the #48 exchange sign, recorded honestly).
  T6. The Compton tie-in: the OFPT two-ordering sums equal 1/(s - m^2)
      and 1/(u - m^2) exactly across lab kinematics; the deformed
      phase fails them at O(1).  The denominators #35/#211 assumed
      are now derived.
  T7. Honest scope.
  T8. Assessment.

Verdict:
  FEYNMAN_PROPAGATOR_FROM_COMPLETE_HISTORIES_THE_TWO_ORDERINGS_ARE_
  OFFER_AND_CONFIRMATION_THE_COHERENT_PHASE_FORCED_BY_THE_FROZEN_BULK
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

_CACHE: dict = {}

# ========================================================================
# SECTION A - single-mode Green functions (conventions in the docstring)
# ========================================================================

_W0 = 1.7          # reference mode frequency for the mode-level checks
_EPS = 1e-3        # complete-absorption regulator for exact identities


def g_ret(t: np.ndarray, w: float = _W0) -> np.ndarray:
    return np.where(t > 0, np.sin(w * t) / w, 0.0)


def g_adv(t: np.ndarray, w: float = _W0) -> np.ndarray:
    return np.where(t < 0, -np.sin(w * t) / w, 0.0)


def g_feyn(t: np.ndarray, w: float = _W0) -> np.ndarray:
    return (1j / (2 * w)) * np.exp(-1j * w * np.abs(t))


def ordering_plus(delta, w: float = _W0, eps: float = _EPS):
    """theta(t) piece of G_F, half-line transformed: the 'offer' ordering."""
    return -(1.0 / (2 * w)) / (delta - w + 1j * eps)


def ordering_minus(delta, w: float = _W0, eps: float = _EPS):
    """theta(-t) piece of G_F: the 'confirmation' ordering."""
    return +(1.0 / (2 * w)) / (delta + w - 1j * eps)


def covariant_pole(delta, w: float = _W0, eps: float = _EPS):
    """Exact finite-eps covariant pole: the coherent sum, closed form."""
    return -(w - 1j * eps) / (w * (delta ** 2 - (w - 1j * eps) ** 2))


# ========================================================================
# SECTION B - the S^3 conformal tower (R = 1): omega_l = l + 1
# ========================================================================


def tower_kernel(t: float, psi: np.ndarray, sign: int = +1,
                 nmax: int = 400, sig: float = 0.02) -> np.ndarray:
    """Mode-sum retarded (sign=+1, support t>0) or advanced (sign=-1,
    support t<0) conformal-scalar Green kernel on S^3 x R at angular
    separation psi, Gaussian-smeared in mode number for summability.

    G(t, psi) = theta(sign*t) * sign *
                sum_n e^{-(n sig)^2/2} sin(n t) sin(n psi)
                / (2 pi^2 sin psi),          n = l + 1 = omega_l R.
    """
    if sign * t <= 0:
        return np.zeros_like(psi)
    n = np.arange(1, nmax + 1)
    damp = np.exp(-0.5 * (n * sig) ** 2)
    s = np.sum(damp * np.sin(n * t)[None, :] * np.sin(np.outer(psi, n)),
               axis=1)
    return sign * s / (2 * math.pi ** 2 * np.sin(psi))


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Derive the Feynman propagator, its two completion '
            'orderings, and their coherent relative phase from the '
            'frozen classical bulk: static geometry forces the '
            'time-symmetric Wheeler-Feynman field, the closed S^3 '
            'makes complete absorption a geometric fact, and the '
            'complete-history (transactional) boundary condition '
            'selects G_F uniquely - closing the one element of the '
            'Compton amplitude (#35/#45/#46/#211) that was imported '
            'rather than derived: the i*epsilon time structure.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The transactional construction of G_F
# ========================================================================


def test_T2_transactional_construction() -> dict:
    w = _W0
    h = 1e-7

    # (i) jump conditions: d/dt G jumps by +1 across t = 0 for all three
    jump_ret = (math.cos(w * h)) - 0.0                     # ret: 0 for t<0
    jump_adv = 0.0 - (-math.cos(-w * h))                   # adv: 0 for t>0
    dgf = lambda t: (1j / (2 * w)) * (-1j * w * np.sign(t)) \
        * np.exp(-1j * w * abs(t))
    jump_feyn = complex(dgf(h) - dgf(-h))
    jumps_ok = (abs(jump_ret - 1) < 1e-6 and abs(jump_adv - 1) < 1e-6
                and abs(jump_feyn - 1) < 1e-6)

    # (ii) the absorber response: G_F - Gbar = (i/2w) cos(wt), homogeneous
    t = np.linspace(-30, 30, 200001)
    gbar = 0.5 * (g_ret(t, w) + g_adv(t, w))
    resp = g_feyn(t, w) - gbar
    target = (1j / (2 * w)) * np.cos(w * t)
    decomp_err = float(np.max(np.abs(resp - target)))
    dt = t[1] - t[0]
    d2 = (resp[2:] - 2 * resp[1:-1] + resp[:-2]) / dt ** 2
    homog_err = float(np.max(np.abs(d2 + w ** 2 * resp[1:-1])))

    # (iii) uniqueness: complete absorption (pure e^{-iwt} future, pure
    # e^{+iwt} past) on Gbar + a cos + b sin is a 2x2 system
    A = np.array([[0.5, 1 / 2j], [0.5, -1 / 2j]], dtype=complex)
    rhs = np.array([-1 / (4j * w), -1 / (4j * w)], dtype=complex)
    a, b = np.linalg.solve(A, rhs)
    uniq_err = float(abs(a - 1j / (2 * w)) + abs(b))
    cond = float(np.linalg.cond(A))          # well-posed, unique

    # (iv) frequency-splitting: Gt_F(W) = Gt_ret(W) for W>0, Gt_adv(W)
    # for W<0, exact as eps -> 0 (linear rate)
    mism = {}
    W = np.linspace(-6, 6, 4001)
    mask = (np.abs(np.abs(W) - w) > 0.2) & (np.abs(W) > 0.05)
    for eps in (1e-1, 1e-2, 1e-3):
        gr = 1.0 / (w ** 2 - (W + 1j * eps) ** 2)
        ga = 1.0 / (w ** 2 - (W - 1j * eps) ** 2)
        gf = ordering_plus(W, w, eps) + ordering_minus(W, w, eps)
        split = np.where(W > 0, gr, ga)
        mism[eps] = float(np.max(np.abs(gf - split)[mask]))
    scaling = mism[1e-1] / mism[1e-2], mism[1e-2] / mism[1e-3]
    split_ok = (mism[1e-3] < 5e-4
                and 8 < scaling[0] < 12 and 8 < scaling[1] < 12)

    ok = (jumps_ok and decomp_err < 1e-12 and homog_err < 1e-6
          and uniq_err < 1e-14 and cond < 10 and split_ok)
    return {
        'name': 'T2_transactional_construction',
        'description': (
            'G_F = Gbar + absorber response, with the response the '
            'unique homogeneous field satisfying complete absorption'
        ),
        'jump_conditions_ok': bool(jumps_ok),
        'decomposition_error': decomp_err,
        'response_homogeneity': homog_err,
        'uniqueness_error': uniq_err,
        'uniqueness_condition_number': cond,
        'response_coefficients': f"a = i/(2w) = {a.imag:.6f}i, b = {abs(b):.1e}",
        'frequency_split_mismatch': {f"{k:.0e}": v for k, v in mism.items()},
        'split_scaling_per_decade': [float(s) for s in scaling],
        'pass': bool(ok),
    }


# ========================================================================
# T3. The two completion orderings and the covariant pole
# ========================================================================


def test_T3_two_orderings() -> dict:
    w, eps = _W0, _EPS

    # (i) finite-T truncation of the offer ordering: undamped never
    # converges (late-time spread stays at the oscillation amplitude);
    # the eps-damped version converges to the closed form
    delta = 2.3
    Ts = np.linspace(50, 500, 400)
    trunc = (1j / (2 * w)) * (np.exp(1j * (delta - w) * Ts) - 1) \
        / (1j * (delta - w))
    late = trunc[-100:]
    spread = float(np.max(np.abs(late - np.mean(late))))
    osc_amp = 1 / (2 * w * abs(delta - w))
    eps_d = 0.05
    trunc_d = (1j / (2 * w)) \
        * (np.exp(1j * (delta - w + 1j * eps_d) * Ts) - 1) \
        / (1j * (delta - w + 1j * eps_d))
    damped_err = float(abs(trunc_d[-1] - ordering_plus(delta, w, eps_d)))

    # (ii) the exact covariant-pole identity on a random (Delta, w) grid
    rng = np.random.default_rng(7)
    D = rng.uniform(-5, 5, 400)
    ww = rng.uniform(0.3, 4.0, 400)
    S = ordering_plus(D, ww, eps) + ordering_minus(D, ww, eps)
    pole = covariant_pole(D, ww, eps)
    pole_err = float(np.max(np.abs(S - pole) / np.abs(pole)))

    # (iii) covariance: the sum is exactly even in Delta -> -Delta;
    # each ordering alone violates evenness at O(1)
    Dg = np.linspace(0.2, 4.0, 50)
    S_p = ordering_plus(Dg, w, eps) + ordering_minus(Dg, w, eps)
    S_m = ordering_plus(-Dg, w, eps) + ordering_minus(-Dg, w, eps)
    sum_even = float(np.max(np.abs(S_p - S_m) / np.abs(S_p)))
    one = ordering_plus(Dg, w, eps)
    one_m = ordering_plus(-Dg, w, eps)
    single_viol = float(np.min(np.abs(one - one_m) / np.abs(one)))

    ok = (spread > 0.5 * osc_amp and damped_err < 1e-9
          and pole_err < 1e-12 and sum_even < 1e-12
          and single_viol > 0.2)
    return {
        'name': 'T3_two_orderings',
        'description': (
            'the theta(t)/theta(-t) split gives the two OFPT energy '
            'denominators; individually regulator-dependent and '
            'non-covariant, coherently summing to the exact pole'
        ),
        'undamped_truncation_spread': spread,
        'oscillation_amplitude': float(osc_amp),
        'damped_convergence_error': damped_err,
        'covariant_pole_identity_error': pole_err,
        'coherent_sum_evenness': sum_even,
        'single_ordering_evenness_violation': single_viol,
        'pass': bool(ok),
    }


# ========================================================================
# T4. The frozen bulk supplies both conditions
# ========================================================================


def test_T4_frozen_bulk_conditions() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    psi = np.linspace(0.05, math.pi - 0.05, 2000)

    # (i) the t -> -t isometry on the conformal tower omega_l = (l+1)/R:
    # G_ret(t, psi) = G_adv(-t, psi), checked mode-sum vs mode-sum
    tr_err = 0.0
    for t in (0.4, 1.1, 2.0, math.pi - 0.2):
        gr = tower_kernel(t, psi, sign=+1)
        ga = tower_kernel(-t, psi, sign=-1)
        tr_err = max(tr_err, float(np.max(np.abs(gr - ga))))

    # (ii) complete absorption as geometry: the retarded front peaks at
    # psi = t (launch), refocuses at the antipode at t = pi R with O(1)
    # of the |G|^2 mass concentrated there, and returns at t = 2 pi R
    peak_launch = float(psi[np.argmax(np.abs(tower_kernel(0.6, psi)))])
    prof_refocus = np.abs(tower_kernel(math.pi, psi)) ** 2
    peak_refocus = float(psi[np.argmax(prof_refocus)])
    conc_refocus = float(np.sum(prof_refocus[psi > math.pi - 0.3])
                         / np.sum(prof_refocus))
    prof_mid = np.abs(tower_kernel(math.pi / 2, psi)) ** 2
    conc_mid = float(np.sum(prof_mid[psi > math.pi - 0.3])
                     / np.sum(prof_mid))
    peak_return = float(psi[np.argmax(np.abs(
        tower_kernel(2 * math.pi - 0.6, psi)))])

    ok = (tr_err < 1e-12
          and abs(peak_launch - 0.6) < 0.05
          and peak_refocus > math.pi - 0.11
          and conc_refocus > 0.4 and conc_mid < 1e-3
          and abs(peak_return - 0.6) < 0.05)
    out = {
        'name': 'T4_frozen_bulk_conditions',
        'description': (
            'static bulk: G_ret(t) = G_adv(-t) on the tower (time '
            'symmetry); closed bulk: antipodal refocus + full return '
            '(complete absorption is geometric, not assumed)'
        ),
        'time_reversal_identity_error': tr_err,
        'launch_peak_at_psi': peak_launch,
        'refocus_peak_at_psi': peak_refocus,
        'refocus_concentration_near_antipode': conc_refocus,
        'midtime_concentration_near_antipode': conc_mid,
        'return_peak_at_psi': peak_return,
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The coherent relative phase and the deform test
# ========================================================================


def test_T5_coherent_phase() -> dict:
    w, eps = _W0, _EPS
    D = np.linspace(-5, 5, 801)
    mask = np.abs(np.abs(D) - w) > 0.25
    Ip, Im_ = ordering_plus(D, w, eps), ordering_minus(D, w, eps)

    def diagnostics(phi):
        S = Ip + np.exp(1j * phi) * Im_
        even = float(np.max((np.abs(S - S[::-1]) / np.abs(S))[mask]))
        pf = float(np.max(np.abs(
            S * (D ** 2 - (w - 1j * eps) ** 2)
            * (-w / (w - 1j * eps)) - 1)[mask]))
        return even, pf

    even0, pf0 = diagnostics(0.0)
    deform = {phi: diagnostics(phi) for phi in (0.3, math.pi / 2, math.pi)}
    deform_broken = all(e > 0.3 and p > 0.3 for e, p in deform.values())

    # the engine already enforces the coherent point: retro_phase_match
    # weight is maximal exactly at zero total phase (and at the pi
    # branch - the #48 exchange sign, a separate discrete sector)
    from geometrodynamics.transaction.handshake import retro_phase_match
    phis = np.linspace(-math.pi, math.pi, 721)
    weights = np.array([retro_phase_match(np.exp(1j * 0.4),
                                          np.exp(1j * (p - 0.4)))[1]
                        for p in phis])
    peaks = phis[weights > 1 - 1e-12]
    peak_at_zero = bool(np.any(np.abs(peaks) < 1e-6))
    peak_branches = sorted(set(np.round(np.abs(peaks) / math.pi).astype(int)))
    monotone = bool(np.all(np.diff(
        weights[(phis >= 0) & (phis <= math.pi / 2 + 1e-9)]) <= 1e-12))

    ok = (even0 < 1e-12 and pf0 < 1e-12 and deform_broken
          and peak_at_zero and monotone)
    return {
        'name': 'T5_coherent_phase',
        'description': (
            'phi = 0 forced by the t -> -t isometry; any deformed '
            'e^{i phi} breaks covariance and the pole form at O(1); '
            'the engine phase closure peaks at the coherent point'
        ),
        'phi0_evenness': even0,
        'phi0_pole_form_deviation': pf0,
        'deform_diagnostics': {
            f"{phi:.3f}": {'evenness_violation': e, 'pole_form_dev': p}
            for phi, (e, p) in deform.items()},
        'engine_peak_at_zero': peak_at_zero,
        'engine_peak_branches_x_pi': [int(b) for b in peak_branches],
        'engine_weight_monotone_from_peak': monotone,
        'pass': bool(ok),
    }


# ========================================================================
# T6. The Compton tie-in: the derived denominators
# ========================================================================


def test_T6_compton_denominators() -> dict:
    m = 1.0
    err_s, err_u, deform_rel = [], [], []
    for wg in (0.1, 0.5, 1.0, 3.0):
        for th in np.linspace(0.1, math.pi - 0.1, 7):
            wp = wg / (1 + (wg / m) * (1 - math.cos(th)))
            s = m ** 2 + 2 * m * wg
            u = m ** 2 - 2 * m * wp
            # s-channel: Delta = m + wg through the internal line with
            # 3-momentum k, E_n = sqrt(wg^2 + m^2)
            En = math.hypot(wg, m)
            ofpt_s = (1 / (2 * En)) * (1 / (m + wg - En)
                                       - 1 / (m + wg + En))
            err_s.append(abs(ofpt_s - 1 / (s - m ** 2)))
            # u-channel: Delta = m - wp, internal 3-momentum -k',
            # E_n' = sqrt(wp^2 + m^2); u < m^2 always (no pole), the
            # identity holds off-shell
            Enp = math.hypot(wp, m)
            ofpt_u = (1 / (2 * Enp)) * (1 / (m - wp - Enp)
                                        - 1 / (m - wp + Enp))
            err_u.append(abs(ofpt_u - 1 / (u - m ** 2)))
            # the deformed phase fails the s-channel denominator
            phi = math.pi / 2
            ofpt_def = (1 / (2 * En)) * (1 / (m + wg - En)
                                         - np.exp(1j * phi)
                                         / (m + wg + En))
            deform_rel.append(abs(ofpt_def - 1 / (s - m ** 2))
                              * (s - m ** 2))
    ok = (max(err_s) < 1e-12 and max(err_u) < 1e-12
          and min(deform_rel) > 0.02)
    return {
        'name': 'T6_compton_denominators',
        'description': (
            'the two-ordering coherent sums equal the 1/(s - m^2) and '
            '1/(u - m^2) denominators the Compton arc assumed - exact '
            'across lab kinematics; the deformed phase fails them'
        ),
        's_channel_identity_error': float(max(err_s)),
        'u_channel_identity_error': float(max(err_u)),
        'deformed_phase_min_relative_error': float(min(deform_rel)),
        'kinematics_grid': '4 photon energies x 7 angles, electron at rest',
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'DERIVED here: the time structure of the propagator - the '
        'Feynman i*epsilon, the two ordering completions, and their '
        'coherent relative phase - from the static + closed frozen '
        'bulk with the complete-history boundary condition.',
        'NOT rederived here: the spinor numerators and vertex factors '
        '(#37-#44, #46); this PR supplies the denominators those '
        'results multiply.',
        'The epsilon -> 0 limit is the complete-absorption '
        'idealization: on the closed bulk, epsilon is a physical '
        'inverse absorption time, not a formal regulator - but its '
        'limit is still taken.',
        'The analysis is per-mode on the free (tree-level) propagator; '
        'no loop or interacting-history statement is made.',
        'The geometry is frozen: no backreaction of the exchanged '
        'field on the bulk (consistent with the program framing).',
        'The engine phase closure also admits the pi branch - the #48 '
        'Moebius exchange sign, a discrete topological sector distinct '
        'from the continuous phase derived here; recorded, not hidden.',
    ]
    return {
        'name': 'T7_honest_scope',
        'description': 'what this PR does and does not establish',
        'scope': scope,
        'pass': True,
    }


# ========================================================================
# T8. Assessment
# ========================================================================


def test_T8_assessment() -> dict:
    t2 = test_T2_transactional_construction()
    t3 = test_T3_two_orderings()
    t4 = test_T4_frozen_bulk_conditions()
    t5 = test_T5_coherent_phase()
    t6 = test_T6_compton_denominators()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The transactional derivation closes: on a frozen bulk that is '
        'static (time symmetry forced) and closed (complete absorption '
        'geometric), the complete-history boundary condition selects '
        'the Feynman propagator uniquely; its two ordering completions '
        'are the offer and confirmation segments of one classical '
        'history, and their relative phase is forced to zero by the '
        't -> -t isometry - the deform test shows any other phase '
        'destroys covariance, the pole, and the Compton denominators. '
        'What QED postulates as the i*epsilon prescription is, on this '
        'geometry, a theorem about complete histories.'
        if core else
        'INCONCLUSIVE - a core identity failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the transactional derivation',
        'core_green': bool(core),
        'assessment': assessment,
        'pass': bool(core),
    }


# ========================================================================
# Probe runner
# ========================================================================


def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_transactional_construction(),
        test_T3_two_orderings(),
        test_T4_frozen_bulk_conditions(),
        test_T5_coherent_phase(),
        test_T6_compton_denominators(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "FEYNMAN_PROPAGATOR_FROM_COMPLETE_HISTORIES_THE_TWO_"
            "ORDERINGS_ARE_OFFER_AND_CONFIRMATION_THE_COHERENT_PHASE_"
            "FORCED_BY_THE_FROZEN_BULK"
        )
        verdict = (
            "DERIVED (the argument is in "
            "docs/compton_transactional_propagator.md).\n\n"
            "THE CONSTRUCTION. The static bulk forces the "
            "time-symmetric Wheeler-Feynman field; the closed S^3 "
            "makes complete absorption a geometric fact (antipodal "
            f"refocus concentration {t4['refocus_concentration_near_antipode']:.2f} "
            f"at t = pi R vs {t4['midtime_concentration_near_antipode']:.0e} "
            "at mid-flight, full return at 2 pi R); the "
            "complete-history boundary condition is a well-posed 2x2 "
            f"system (uniqueness error {t2['uniqueness_error']:.0e}) "
            "whose unique solution is the Feynman propagator, "
            "G_F = Gbar + (i/2w) cos(wt).\n\n"
            "THE TWO ORDERINGS. The theta(+t)/theta(-t) split of G_F "
            "is offer and confirmation; each half-line transform is an "
            "OFPT energy denominator - individually regulator-"
            f"dependent (truncation spread {t3['undamped_truncation_spread']:.2f}) "
            f"and non-covariant (evenness violation {t3['single_ordering_evenness_violation']:.2f}) "
            "- and the coherent sum is the exact covariant pole "
            f"(identity to {t3['covariant_pole_identity_error']:.0e}).\n\n"
            "THE PHASE. The t -> -t isometry (tower identity to "
            f"{t4['time_reversal_identity_error']:.0e}) forces relative "
            f"phase 0; deforming it breaks evenness and the pole form "
            "at O(1) and fails the Compton denominators (min error "
            f"{t6['deformed_phase_min_relative_error']:.2f} at phi = "
            "pi/2). The s- and u-channel denominators the Compton arc "
            "assumed are recovered exactly "
            f"({t6['s_channel_identity_error']:.0e} / "
            f"{t6['u_channel_identity_error']:.0e}): the i*epsilon "
            "prescription is, on this geometry, a theorem about "
            "complete histories."
        )
    else:
        verdict_class = "TRANSACTIONAL_PROPAGATOR_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core identity failed; re-examine before "
            "quoting the derivation."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The transactional Compton propagator: the Feynman "
            "propagator, its two completion orderings (offer and "
            "confirmation), and their coherent relative phase, derived "
            "from the frozen classical bulk - static implies time "
            "symmetry, closed implies complete absorption, and the "
            "complete-history boundary condition selects G_F uniquely"
        ),
        "executes": (
            "the next Compton frontier: replace the imported "
            "i*epsilon time structure with a derivation"
        ),
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return "<array>"
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The transactional Compton propagator (PR #213)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/compton_transactional_propagator.md` "
        "- the Feynman propagator, its two completion orderings, and "
        "their coherent relative phase from the frozen bulk. *(QFT on "
        "the fixed classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: derive the imported time structure",
        "T2": "G_F = Gbar + unique absorber response",
        "T3": "the two orderings -> the exact covariant pole",
        "T4": "static + closed: both conditions are geometric",
        "T5": "the relative phase forced; deform test refutes others",
        "T6": "the Compton denominators derived, not assumed",
        "T7": "honest scope",
        "T8": "assessment",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**Class:** `{s['verdict_class']}`")
    out.append("")
    out.append(s["verdict"])
    out.append("")
    return "\n".join(out)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_compton_transactional_propagator_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
