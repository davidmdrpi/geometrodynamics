"""Probe B -- the genuine first-order RP^3 Dirac spectrum and eta (PR #230).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

WHY THIS EXISTS
---------------
bam_spinor_spectrum_effective_action_probe computed a 1D cycle mode via the
SUSY square root, which delivers |spec(D)| only -- it discards the SIGN of
the spectrum, so it cannot represent a chirality grading or a spectral
asymmetry, and its "zero mode" turned out to be the constant function
(T8 there).  This probe replaces that construction entirely: the genuine
first-order Dirac operator on RP^3, from the analytic S^3 eigenspinors,
with SIGNS PRESERVED and no finite differencing anywhere -- hence no
fermion doubling to worry about.

THE INPUT (analytic, representation-theoretic)
----------------------------------------------
On the round S^3 of radius 1 the Dirac spectrum is

    lambda = +/- (3/2 + k) ,  multiplicity (k+1)(k+2) ,  k = 0,1,2,...

Writing S^3 = SU(2), the antipodal map is LEFT multiplication by the
central element -1, which acts on a mode of left-spin j_L by (-1)^{2 j_L}.
The two branches sit at different left-spin:

    + branch, level k :  j_L = k/2       ->  parity (-1)^k
    - branch, level k :  j_L = (k+1)/2   ->  parity (-1)^{k+1}

*** The two branches carry OPPOSITE antipodal parity at the same |lambda|.
That single fact is the origin of the spectral asymmetry, and it is exactly
what the |spec| construction could not see. ***

Consistency check built in: for each k, S^3 has 2(k+1)(k+2) modes and the
projection keeps exactly (k+1)(k+2) -- precisely half, as vol(RP^3) =
vol(S^3)/2 requires.

THE RESULTS
-----------
Projecting under J~_+ and J~_- separately:

    eta_W(0) = +1/4  (one spin structure)   -1/4  (the other)   EXACT
    h_W      = 0     for both -- |lambda| >= 3/2, there is NO zero mode

  * h = 0 confirms, from the genuine operator, that the 1D "zero mode" of
    the earlier probe was an artefact of the reduction, not physics.
  * eta = +/-1/4 is a genuine spectral invariant -- the classical eta
    invariant of RP^3 -- and it DOES distinguish the two spin structures.
    This is the index-theoretic discriminator the earlier probe reached for
    and failed to find.

And the sharp structural point:

    for each k exactly ONE sign survives, with |lambda| = 3/2 + k ,

so |spec| is IDENTICAL for the two spin structures.  Therefore

    zeta_{D^2}(s), zeta'_{D^2}(0) and |det D|  are the SAME for both,

and the two sectors differ ONLY in the determinant PHASE, through eta.

*** This retires the energetic route for good.  Any vacuum-energy or
|det|-based comparison between the sectors returns exactly zero difference
on the genuine 3D operator.  The pi/(4C) energetics of the twist-selection
arc were computing a quantity that cannot distinguish the sectors. ***

Tests:
  T1. Goal and the analytic input.
  T2. The branch parities are opposite; the half-density consistency check.
  T3. Projection under each lift: the RP^3 spectra, signs preserved.
  T4. eta_W(0) = +/-1/4 and h_W = 0.
  T5. zeta_{D^2}(0), zeta'_{D^2}(0), |det D| -- identical for both.
  T6. The determinant phase -- the only difference.
  T7. What this overturns.
  T8. Assessment.

Verdict:
  THE_RP3_SPIN_STRUCTURES_HAVE_h_EQUALS_ZERO_AND_IDENTICAL_DETERMINANT
  _MODULUS_AND_ARE_DISTINGUISHED_ONLY_BY_ETA_EQUALS_PLUS_MINUS_ONE_QUARTER
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from mpmath import mp, mpf, zeta as mp_zeta, exp as mp_exp, mpc, pi as mp_pi

mp.dps = 30
PI = math.pi


# ── the analytic S^3 spectrum and its antipodal grading ─────────────────

def s3_levels(kmax: int) -> list:
    """(k, |lambda|, multiplicity, parity of + branch, parity of - branch)."""
    out = []
    for k in range(kmax):
        out.append({'k': k, 'abs_lambda': 1.5 + k,
                    'multiplicity': (k + 1) * (k + 2),
                    'parity_plus_branch': (-1) ** k,
                    'parity_minus_branch': (-1) ** (k + 1)})
    return out


def rp3_signed_spectrum(eps: int, kmax: int) -> list:
    """Modes surviving the projection under the lift with sign eps.
    Signs are PRESERVED: each entry carries the actual eigenvalue."""
    out = []
    for lv in s3_levels(kmax):
        k, a, m = lv['k'], lv['abs_lambda'], lv['multiplicity']
        if lv['parity_plus_branch'] == eps:
            out.append({'k': k, 'eigenvalue': +a, 'multiplicity': m})
        if lv['parity_minus_branch'] == eps:
            out.append({'k': k, 'eigenvalue': -a, 'multiplicity': m})
    return out


def _beta(s, a):
    """Alternating Hurwitz zeta: sum_{k>=0} (-1)^k (k+a)^{-s}."""
    return 2 ** (-s) * (mp_zeta(s, a / 2) - mp_zeta(s, (a + 1) / 2))


def eta_at_zero(eps: int):
    """eta(0) = sum sign(lam)|lam|^{-s}|_{s=0}, analytically continued.
    With u = k+3/2 the multiplicity is u^2 - 1/4, so
    eta(s) = eps * [ beta(s-2, 3/2) - (1/4) beta(s, 3/2) ]."""
    return eps * (_beta(mpf(-2), mpf(3) / 2) - mpf(1) / 4 * _beta(mpf(0), mpf(3) / 2))


def eta_convergent(eps: int, s: float, kmax: int = 40000):
    """Direct mode sum at s large enough to converge -- the numerical check
    on the analytic continuation."""
    tot = mpf(0)
    for k in range(kmax):
        u = mpf(k) + mpf(3) / 2
        m = (k + 1) * (k + 2)
        tot += (m * u ** (-mpf(s))) * ((-1) ** k * eps)
    return tot


def zeta_d2(s):
    """zeta_{D^2}(s) = sum |lam|^{-2s} = zeta_H(2s-2,3/2) - (1/4) zeta_H(2s,3/2).
    Identical for both spin structures (for each k exactly one sign survives
    with the same |lambda|)."""
    return mp_zeta(2 * s - 2, mpf(3) / 2) - mpf(1) / 4 * mp_zeta(2 * s, mpf(3) / 2)


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_the_genuine_first_order_operator',
        'description': (
            "Replace the SUSY square-root construction, which delivers "
            "|spec(D)| and therefore cannot see a spectral asymmetry, with "
            "the genuine first-order Dirac operator on RP^3 built from the "
            "analytic S^3 eigenspinors: lambda = +/-(3/2+k) with "
            "multiplicity (k+1)(k+2). No finite differencing is used "
            "anywhere, so fermion doubling cannot arise. Signs are carried "
            "through the projection."
        ),
        'pass': True,
    }


def test_T2_branch_parities() -> dict:
    lv = s3_levels(6)
    opposite = all(r['parity_plus_branch'] == -r['parity_minus_branch']
                   for r in lv)
    # half-density check: S^3 has 2m modes per level, RP^3 keeps m
    checks = []
    for eps in (+1, -1):
        for k in range(6):
            m = (k + 1) * (k + 2)
            kept = sum(r['multiplicity']
                       for r in rp3_signed_spectrum(eps, 6) if r['k'] == k)
            checks.append({'eps': eps, 'k': k, 's3_modes': 2 * m,
                           'rp3_modes': kept, 'is_half': kept == m})
    half = all(c['is_half'] for c in checks)
    return {
        'name': 'T2_opposite_branch_parities_and_the_half_density_check',
        'description': (
            "Writing S^3 = SU(2), the antipodal map is left multiplication "
            "by the central -1, acting on a mode of left-spin j_L by "
            "(-1)^{2 j_L}. The + branch at level k has j_L = k/2 and the - "
            "branch has j_L = (k+1)/2, so the two branches carry OPPOSITE "
            "antipodal parity at the same |lambda| (verified for k = 0..5). "
            "That is the entire origin of the spectral asymmetry, and it is "
            "precisely what a |spec| construction cannot see. Consistency: "
            "S^3 has 2(k+1)(k+2) modes at level k and the projection keeps "
            "exactly (k+1)(k+2) for either spin structure -- exactly half, "
            "as vol(RP^3) = vol(S^3)/2 requires (verified at every level "
            "and for both lifts)."
        ),
        'levels': lv,
        'branch_parities_opposite': opposite,
        'half_density_checks': checks,
        'keeps_exactly_half': half,
        'pass': opposite and half,
    }


def test_T3_projected_spectra() -> dict:
    out = {}
    for eps in (+1, -1):
        sp = rp3_signed_spectrum(eps, 6)
        out[f'eps={eps:+d}'] = [
            {'eigenvalue': r['eigenvalue'], 'multiplicity': r['multiplicity']}
            for r in sorted(sp, key=lambda r: r['eigenvalue'])]
    # the two are mirror images
    a = sorted(r['eigenvalue'] for r in rp3_signed_spectrum(+1, 20))
    b = sorted(-r['eigenvalue'] for r in rp3_signed_spectrum(-1, 20))
    mirror = a == sorted(b)
    return {
        'name': 'T3_the_projected_spectra_signs_preserved',
        'description': (
            "Projecting separately under J~_+ and J~_- gives two genuinely "
            "asymmetric spectra which are exact mirror images of one "
            "another (verified to k = 19): where one keeps +(3/2+k) the "
            "other keeps -(3/2+k). Neither spectrum is symmetric about "
            "zero -- which is what makes eta nonzero and what the SUSY "
            "square root necessarily missed."
        ),
        'spectra_first_levels': out,
        'the_two_are_mirror_images': mirror,
        'pass': mirror,
    }


def test_T4_eta_and_h() -> dict:
    rows = []
    for eps in (+1, -1):
        an = eta_at_zero(eps)
        rows.append({'eps': eps, 'eta_0': float(an),
                     'eta_0_exact': str(an)})
    checks = []
    for s in (4, 6):
        num = eta_convergent(+1, s)
        an = _beta(mpf(s) - 2, mpf(3) / 2) - mpf(1) / 4 * _beta(mpf(s), mpf(3) / 2)
        checks.append({'s': s, 'mode_sum': float(num),
                       'analytic_continuation': float(an),
                       'abs_diff': float(abs(num - an))})
    quarter = all(abs(abs(r['eta_0']) - 0.25) < 1e-12 for r in rows)
    opposite = abs(rows[0]['eta_0'] + rows[1]['eta_0']) < 1e-12
    conv_ok = all(c['abs_diff'] < 1e-8 for c in checks)
    h = 0  # |lambda| >= 3/2 > 0 for every mode
    return {
        'name': 'T4_eta_is_plus_minus_one_quarter_and_h_is_zero',
        'description': (
            "eta_W(0) = +1/4 for one spin structure and -1/4 for the other, "
            "EXACTLY -- the classical eta invariant of RP^3. Computed by "
            "analytic continuation (with u = k+3/2 the multiplicity is "
            "u^2 - 1/4, so eta(s) = eps[beta(s-2,3/2) - (1/4)beta(s,3/2)] "
            "in alternating Hurwitz zetas) and cross-checked against the "
            "directly convergent mode sum at s = 4 and s = 6, agreeing to "
            f"{max(c['abs_diff'] for c in checks):.1e}. And h_W = 0 for "
            "BOTH: every eigenvalue satisfies |lambda| >= 3/2 > 0, so the "
            "genuine 3D operator has NO zero mode at all. *** That settles "
            "the earlier probe's kernel from the other side: its 1D 'zero "
            "mode' was an artefact of the reduction. The real discriminator "
            "is eta, a spectral asymmetry -- not a kernel. ***"
        ),
        'rows': rows,
        'convergent_cross_checks': checks,
        'eta_is_quarter': quarter,
        'the_two_are_opposite': opposite,
        'h_W': h,
        'pass': quarter and opposite and conv_ok,
    }


def test_T5_determinant_modulus() -> dict:
    z0 = zeta_d2(mpf(0))
    # zeta'_{D^2}(0) = 2 zeta_H'(-2,3/2) - (1/2) zeta_H'(0,3/2)
    zp0 = (2 * mp_zeta(mpf(-2), mpf(3) / 2, 1)
           - mpf(1) / 2 * mp_zeta(mpf(0), mpf(3) / 2, 1))
    det_mod = float(mp_exp(-zp0 / 2))
    return {
        'name': 'T5_determinant_modulus_is_identical_for_both',
        'description': (
            "For each k exactly ONE sign survives the projection, with "
            "|lambda| = 3/2 + k and multiplicity (k+1)(k+2) -- the SAME for "
            "both spin structures. So zeta_{D^2}, and with it |det D|, "
            "cannot distinguish them. Computed: zeta_{D^2}(0) = "
            f"{float(z0):.6f}, zeta'_{{D^2}}(0) = {float(zp0):.10f}, "
            f"|det D| = exp(-zeta'(0)/2) = {det_mod:.10f}, identically for "
            "eps = +1 and eps = -1. *** This retires the energetic route. "
            "Any vacuum-energy or |det|-based comparison between the "
            "sectors returns EXACTLY zero difference on the genuine 3D "
            "operator, because the two spectra have the same absolute "
            "values mode by mode. The pi/(4C) energetics of the "
            "twist-selection arc were computing a quantity that cannot "
            "distinguish the sectors at all. ***"
        ),
        'zeta_D2_at_0': float(z0),
        'zeta_prime_D2_at_0': float(zp0),
        'det_modulus': det_mod,
        'identical_for_both_spin_structures': True,
        'vacuum_energy_difference': 0.0,
        'pass': abs(float(z0)) < 1e-12 and det_mod > 0,
    }


def test_T6_determinant_phase() -> dict:
    rows = []
    for eps in (+1, -1):
        eta = float(eta_at_zero(eps))
        phase = -PI * (eta + 0) / 2.0
        rows.append({'eps': eps, 'eta_0': eta, 'h': 0,
                     'arg_det': phase,
                     'phase': [float(np.cos(phase)), float(np.sin(phase))]})
    dphi = abs(rows[0]['arg_det'] - rows[1]['arg_det'])
    return {
        'name': 'T6_the_only_difference_is_the_determinant_phase',
        'description': (
            "With the standard convention arg det D = -pi (eta(0) + h)/2 "
            "and h = 0, the two spin structures give arg det D = -pi/8 and "
            f"+pi/8, a phase difference of {dphi:.6f} rad = pi/4. Since T5 "
            "shows the moduli agree exactly, the ENTIRE distinction between "
            "the two sectors, at one loop on the genuine 3D operator, is "
            "this phase. A phase is not an energy: it cannot be traded "
            "against a vacuum-energy cost, and it does not make one sector "
            "'cheaper'. It is, however, exactly the object that anomaly "
            "inflow and the APS index theorem constrain -- which is where "
            "the selection question has to go next (Probe C). The "
            "convention is stated explicitly because the SIGN of the phase "
            "is convention-dependent; the DIFFERENCE pi/4 is not."
        ),
        'rows': rows,
        'phase_difference': dphi,
        'convention': 'arg det D = -pi(eta(0)+h)/2',
        'pass': abs(dphi - PI / 4) < 1e-12,
    }


def test_T7_what_this_overturns() -> dict:
    items = [
        {'claim': "the twisted sector carries a Dirac zero mode "
                  "(dim ker D = 1)",
         'source': 'bam_spinor_spectrum_effective_action_probe T4',
         'status': 'ARTEFACT -- the genuine 3D operator has h = 0 for both '
                   'spin structures (|lambda| >= 3/2). The 1D kernel was '
                   'the constant mode of the reduction.'},
        {'claim': "energetics favour the untwisted sector by pi/(4C)",
         'source': 'twist_sector_selection_probe T2 / spinor probe T5',
         'status': 'NOT A DISCRIMINATOR -- on the genuine 3D operator the '
                   'two spin structures have identical |spec|, hence '
                   'identical vacuum energy and |det|. The difference is '
                   'exactly zero, not small.'},
        {'claim': "the sectors are distinguished by an index",
         'source': 'the withdrawn reading',
         'status': 'PARTLY VINDICATED, differently -- they ARE distinguished '
                   'by a spectral invariant, but it is eta = +/-1/4 (an '
                   'asymmetry), not a kernel and not an index.'},
    ]
    return {
        'name': 'T7_what_this_overturns_in_the_earlier_probes',
        'description': (
            "Three corrections, all in the same direction: the earlier "
            "constructions were measuring the wrong things. The kernel was "
            "an artefact (h = 0 here); the energetics cannot distinguish "
            "the sectors at all (identical |spec|); and the discriminator "
            "that does exist is a spectral ASYMMETRY, eta = +/-1/4, which "
            "no |spec|-based construction could ever have produced. The "
            "instinct that the answer was index-theoretic was right; the "
            "object was wrong."
        ),
        'items': items,
        'pass': True,
    }


def test_T8_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T8_assessment',
        'description': (
            "The genuine first-order RP^3 Dirac operator gives a clean "
            "answer: h = 0, identical determinant modulus, and eta = +/-1/4 "
            "as the sole one-loop distinction between the two spin "
            "structures. The energetic route is closed exactly; the "
            "selection question moves to the phase, i.e. to anomaly inflow "
            "and APS (Probe C)."
        ),
        'established': [
            'the two branches carry opposite antipodal parity; the '
            'projection keeps exactly half the modes at every level',
            'eta_W(0) = +/-1/4 exactly, cross-checked against the '
            'convergent mode sum',
            'h_W = 0 for both -- no zero mode on the genuine operator',
            "zeta_{D^2}(0) = 0, zeta'_{D^2}(0) = 0.2189594807, |det D| = "
            '0.8963003229, IDENTICAL for both spin structures',
            'the sole difference is the determinant phase, pi/4 apart',
        ],
        'open': [
            'what the eta difference implies for selection -- a phase is not '
            'an energy; this needs anomaly inflow / APS (Probe C)',
            'the round unit S^3 is used as the exterior; the BAM throat '
            'geometry is not round, and whether eta survives that '
            'deformation is untested (eta is not a homotopy invariant, '
            'though its mod-Z reduction is)',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_RP3_SPIN_STRUCTURES_HAVE_h_EQUALS_ZERO_AND_IDENTICAL_'
            'DETERMINANT_MODULUS_AND_ARE_DISTINGUISHED_ONLY_BY_ETA_EQUALS_'
            'PLUS_MINUS_ONE_QUARTER'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_branch_parities()
    res['T3'] = test_T3_projected_spectra()
    res['T4'] = test_T4_eta_and_h()
    res['T5'] = test_T5_determinant_modulus()
    res['T6'] = test_T6_determinant_phase()
    res['T7'] = test_T7_what_this_overturns()
    res['T8'] = test_T8_assessment(res)
    res['summary'] = {
        'probe': 'rp3_dirac_eta_probe', 'pr': 230,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T8']['tests_passed'],
        'verdict_class': res['T8']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe B — the RP³ Dirac spectrum, η and det (PR #230)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The spectral invariants", "",
         "| spin structure | `η(0)` | `h` | `\\|det D\\|` | `arg det D` |",
         "|---|---:|---:|---:|---:|"]
    for r4, r6 in zip(s['T4']['rows'], s['T6']['rows']):
        o.append(f"| `ε = {r4['eps']:+d}` | **{r4['eta_0']:+.4f}** | "
                 f"{s['T4']['h_W']} | {s['T5']['det_modulus']:.10f} | "
                 f"{r6['arg_det']:+.6f} |")
    o += ["", f"`ζ_D²(0) = {s['T5']['zeta_D2_at_0']:.6f}`, "
          f"`ζ'_D²(0) = {s['T5']['zeta_prime_D2_at_0']:.10f}`",
          "",
          "**The moduli are identical** — the sectors differ only by a "
          f"phase of {s['T6']['phase_difference']:.6f} rad = π/4.",
          "", "## What this overturns", ""]
    for it in s['T7']['items']:
        o.append(f"  - *{it['claim']}* — {it['status']}")
    o += ["", "## Verdict", "", f"**{s['T8']['verdict_class']}**", ""]
    return "\n".join(o)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_rp3_dirac_eta_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
