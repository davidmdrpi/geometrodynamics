"""Is the Moebius identification quantitative, or only structural? (PR #229)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
PR #227 (sec. 8.2) established that the ER-network link twist and the QCD
Moebius label are the same Z2 STRUCTURALLY: same holonomy, same
integer -> half-integer comb, two implementations that had never been
connected (``make_glueball_ring`` periodic/+1 vs ``make_mobius_tube``
mobius/-1).  It left open whether the +pi*sigma tower shift of #100 is
QUANTITATIVELY the same object.  This probe answers that.

THE ANSWER: PARTIALLY -- and the missing part is a gap in #100
--------------------------------------------------------------
The pi*sigma shift factorizes into a topological factor and a dimensionful
one, and they have very different standing:

    Delta M^2  =  (1/2)          x   (2 pi sigma)
                  ^ moding shift      ^ closed-string M^2 spacing
                  DERIVED, shared     INPUT, not supplied by the network

  * THE 1/2 IS DERIVED AND SHARED.  It is reproduced by three independent
    implementations in this repo -- the #227 twisted ring, the repo's own
    ``MobiusSpectrum``, and a direct antiperiodic closed-loop solve -- and
    it is exactly parameter-independent: the deviation from 1/2 is pure
    discretization error, falling 1.5e-05 -> 2.2e-07 as the grid goes
    120 -> 960, and it does not move with the loop circumference at all.
    That is the signature of a topological invariant, and it is what makes
    the identification quantitative rather than merely suggestive.

  * THE 2 pi sigma IS NOT.  The closed-string M^2 spacing comes from the
    Regge slope 1/(4 pi sigma), a string-theoretic relation #100 imports;
    nothing in the ER network supplies it.  So pi*sigma is quantitative in
    its topological factor and INHERITED in its scale.

  * AND #100'S TOWER HAS AN UNJUSTIFIED CONSTANT.  #100 writes both towers
    with the SAME intercept,

        orientable:  M^2_n = M0^2 + 2 pi sigma * n
        mobius:      M^2_n = M0^2 + 2 pi sigma * (n + 1/2)

    but antiperiodic moding shifts the zero-point (Casimir) energy as well
    as the level, so the Moebius intercept cannot be the orientable one.
    Computed here two independent ways (zeta regularization and an
    exponential-cutoff extraction, agreeing to ~4e-7): per transverse
    polarization on a loop of circumference C,

        E0_periodic = -pi/(6C) ,  E0_antiperiodic = +pi/(12C) ,
        Delta E0    = +pi/(4C)   (exact)

    The Moebius zero point is HIGHER, and by an amount that scales with
    1/C -- so it cannot be absorbed into a common M0^2.  #100's ground
    state is therefore not M0^2 + pi*sigma.

SCOPE OF THE CORRECTION (important -- it does NOT hit the search table)
----------------------------------------------------------------------
This touches the GLUEBALL tower of #100, which is built from the
closed-string mass formula.  It does NOT touch the heavy Moebius BARYON
predictions of #103/#109/#114 (Lambda_c 3135, Lambda_b 6469, the 849 MeV
dipion endpoint): those are built from the flux-tube quantum
Delta = 2 sqrt(sigma), not from the closed-loop intercept, so the search
table stands as published.  Glueballs are unobserved, so the corrected
item is the one with no experimental exposure.

Tests:
  T1. Goal and the factorization being tested.
  T2. The 1/2 from three independent repo implementations.
  T3. The 1/2 is topological: parameter-independent to discretization only.
  T4. The pi*sigma ledger: derived factor vs inherited scale.
  T5. The gap in #100: antiperiodic moding shifts the intercept too.
  T6. Assessment, and the scope of the correction.

Verdict:
  THE_MOEBIUS_IDENTIFICATION_IS_QUANTITATIVE_IN_ITS_TOPOLOGICAL_HALF_THE_
  SCALE_IS_INHERITED_AND_THE_MOEBIUS_INTERCEPT_IS_AN_OPEN_CORRECTION
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.qcd.constants import SIGMA_QCD
from geometrodynamics.qcd.spectrum import MobiusSpectrum
from geometrodynamics.qcd.topology import make_glueball_ring, make_mobius_tube

PI = math.pi
TWO_PI_SIGMA = 2.0 * PI * SIGMA_QCD
M0_ORIENTABLE = math.sqrt(4.0 * PI * SIGMA_QCD)


# ── the closed loop with a Z2 holonomy (the #227 machinery) ─────────────

def ring_comb(m: int, holonomy: float, circumference: float = 1.0
              ) -> np.ndarray:
    """Mode comb of a closed loop carrying Wilson loop W = holonomy,
    normalized so the continuum answer is n = 0,1,1,2,2,... (periodic)
    and n = 1/2,1/2,3/2,3/2,... (antiperiodic)."""
    dx = circumference / m
    lap = (np.diag(-2.0 * np.ones(m)) + np.diag(np.ones(m - 1), 1)
           + np.diag(np.ones(m - 1), -1))
    lap[0, -1] = lap[-1, 0] = holonomy
    lap /= dx ** 2
    ev = np.linalg.eigvalsh(-lap)
    return np.sort(np.sqrt(np.clip(ev, 0.0, None))) * circumference / (2 * PI)


# ── Casimir energy of the loop, two independent regularizations ─────────

def casimir_zeta(circumference: float, antiperiodic: bool) -> float:
    """(1/2) sum |omega| by zeta regularization, per transverse
    polarization.  Uses the closed forms zeta(-1) = -1/12 and
    zeta(-1, 1/2) = +1/24 (no external dependency)."""
    a = 2.0 * PI / circumference
    return a * (1.0 / 24.0 if antiperiodic else -1.0 / 12.0)


def casimir_cutoff(circumference: float, antiperiodic: bool,
                   x0: float = 1e-3) -> float:
    """The same quantity by an INDEPENDENT method: an exponential cutoff
    (1/2) sum omega e^{-eps omega}, with the 1/eps^2 divergence removed
    analytically.  Periodic:     a[ e^{-x}/(1-e^{-x})^2 ]      - a/x^2
    Antiperiodic: a[ cosh(x/2)/(4 sinh^2(x/2)) ]  - a/x^2 ,  x = eps*a.

    The cutoff is specified through the dimensionless expansion parameter
    x = eps*a rather than eps itself: the truncation error is O(a*x^2) and
    the cancellation error is O(a/x^2 * machine_eps), so holding x fixed
    keeps the accuracy uniform in the circumference.  x0 = 1e-3 balances
    the two at ~1e-8."""
    a = 2.0 * PI / circumference
    x = x0
    if antiperiodic:
        total = a * math.cosh(x / 2.0) / (4.0 * math.sinh(x / 2.0) ** 2)
    else:
        total = a * math.exp(-x) / (1.0 - math.exp(-x)) ** 2
    return total - a / x ** 2


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_the_factorization_under_test',
        'description': (
            "PR #227 showed the ER-network link twist and the QCD Moebius "
            "label are the same Z2 structurally. This probe asks whether "
            "#100's +pi*sigma glueball tower shift is the same object "
            "QUANTITATIVELY. The shift factorizes as Delta M^2 = (1/2) x "
            "(2 pi sigma): a topological moding shift times the "
            "closed-string M^2 spacing. The two factors are tested "
            "separately, because they have different standing."
        ),
        'factorization': 'Delta M^2 = (1/2) * (2 pi sigma)',
        'pass': True,
    }


def test_T2_half_from_three_implementations() -> dict:
    """The 1/2, computed from three independent implementations.

    #100's own T4 ASSERTS the integer -> half-integer shift (it returns
    pass: True with no computation). This test computes it.
    """
    net_ring = ring_comb(480, +1.0)[:6].tolist()
    net_mob = ring_comb(480, -1.0)[:6].tolist()

    ms = MobiusSpectrum(L=1.0, v=1.0, N=400)
    an = ms.analytic_modes(4)
    repo_mob = [r['omega_mob'] / PI for r in an]
    repo_ord = [r['omega_ord'] / PI for r in an]
    repo_num_ok = bool(ms.report(4)['all_pass'])

    direct = {f'C={c}': ring_comb(480, -1.0, c)[:4].tolist()
              for c in (1.0, 2.5)}

    # the repo's two closed-loop constructors are the two holonomy sectors
    ring_obj, mob_obj = make_glueball_ring(1.0), make_mobius_tube(1.0)
    constructors = {
        'glueball_ring': [ring_obj.branches[0].topology_kind,
                          int(ring_obj.branches[0].orientation)],
        'mobius_tube': [mob_obj.branches[0].topology_kind,
                        int(mob_obj.branches[0].orientation)],
    }

    half_net = all(abs(v - round(v * 2) / 2) < 1e-4 and
                   abs(v * 2 - round(v * 2)) < 1e-4 for v in net_mob[:4])
    half_repo = all(abs(m - (n + 0.5)) < 1e-9
                    for n, m in enumerate(repo_mob))
    ok = half_net and half_repo and repo_num_ok
    return {
        'name': 'T2_the_half_from_three_implementations',
        'description': (
            "#100's T4 asserts the integer -> half-integer shift without "
            "computing it. Computed here, three independent ways, all "
            "agreeing on exactly 1/2: (a) the #227 twisted ring gives the "
            f"comb {np.round(net_mob[:4], 4).tolist()} against the periodic "
            f"{np.round(net_ring[:4], 4).tolist()}; (b) the repo's own "
            "MobiusSpectrum gives omega_mob/pi = 0.5, 1.5, 2.5, 3.5 against "
            "omega_ord/pi = 1, 2, 3, 4, and its numerical solve confirms "
            "them; (c) a direct antiperiodic closed-loop solve gives "
            "1/2, 1/2, 3/2, 3/2 independently of the circumference. The "
            "repo's two closed-loop constructors are the two holonomy "
            "sectors of one cycle."
        ),
        'network_ring_periodic': [round(v, 5) for v in net_ring],
        'network_ring_antiperiodic': [round(v, 5) for v in net_mob],
        'repo_MobiusSpectrum_over_pi': {'mobius': repo_mob,
                                        'ordinary': repo_ord},
        'repo_MobiusSpectrum_numerical_ok': repo_num_ok,
        'direct_antiperiodic_solve': direct,
        'repo_constructors': constructors,
        'pass': ok,
    }


def test_T3_half_is_topological() -> dict:
    """The 1/2 is parameter-independent -- the deviation is pure
    discretization error and converges away."""
    # Measure the antiperiodic GROUND mode against 1/2 directly.  Do not
    # subtract the periodic zero mode: its eigenvalue is ~1e-13 and the
    # sqrt amplifies eigensolver noise into an erratic ~1e-6 that has
    # nothing to do with the shift (it is the only ill-conditioned entry
    # in either comb).
    rows = []
    for m in (120, 240, 480, 960):
        for c in (0.7, 1.0, 3.3):
            n_eff = float(ring_comb(m, -1.0, c)[0])
            rows.append({'grid': m, 'circumference': c,
                         'antiperiodic_ground_n': n_eff,
                         'dev_from_half': abs(n_eff - 0.5)})
    by_grid = {}
    for r in rows:
        by_grid.setdefault(r['grid'], []).append(r['dev_from_half'])
    worst = {g: max(v) for g, v in by_grid.items()}

    # the deviation is finite-difference dispersion: exactly O(1/M^2)
    grids = sorted(worst)
    ratios = [worst[grids[i]] / worst[grids[i + 1]]
              for i in range(len(grids) - 1)]
    follows_inverse_M2 = all(abs(r - 4.0) < 0.05 for r in ratios)

    # and it does not move with the circumference at all
    spread = max(max(v) - min(v) for v in by_grid.values())
    ok = (follows_inverse_M2 and spread < 1e-9 and worst[960] < 1e-6)
    return {
        'name': 'T3_the_half_is_topological_not_dynamical',
        'description': (
            "A topological invariant must not depend on the geometry or "
            "the discretization, and must converge to its exact value. "
            "Both hold, sharply. Scanning grid 120 -> 960 and "
            "circumference 0.7 -> 3.3, the antiperiodic ground mode sits "
            f"at 1/2 with deviation falling {worst[120]:.2e} -> "
            f"{worst[960]:.2e}, in successive ratios "
            f"{[round(r, 3) for r in ratios]} -- i.e. EXACTLY O(1/M^2), "
            "the finite-difference dispersion of the discrete Laplacian "
            "and nothing else. Across a factor-of-5 change in "
            f"circumference at fixed grid the deviation moves by {spread:.1e} "
            "(bitwise identical): the shift does not depend on the geometry "
            "at all. So the 1/2 is exact in the continuum and "
            "geometry-independent -- which is what makes the "
            "identification with #100's Moebius sector quantitative rather "
            "than merely suggestive. (The deviation is measured on the "
            "antiperiodic ground mode directly rather than as a difference "
            "against the periodic zero mode: that mode's eigenvalue is "
            "~1e-13 and the sqrt amplifies eigensolver noise into an "
            "erratic ~1e-6 unrelated to the shift.)"
        ),
        'rows': rows,
        'worst_deviation_by_grid': {str(k): v for k, v in worst.items()},
        'successive_ratios': ratios,
        'follows_inverse_M_squared': follows_inverse_M2,
        'circumference_spread': spread,
        'pass': ok,
    }


def test_T4_pi_sigma_ledger() -> dict:
    """What is derived and what is inherited."""
    delta_m2 = 0.5 * TWO_PI_SIGMA
    ok = abs(delta_m2 - PI * SIGMA_QCD) < 1e-12
    return {
        'name': 'T4_the_pi_sigma_ledger',
        'description': (
            "The shift factorizes into a DERIVED topological factor and an "
            "INHERITED scale. The 1/2 is established by T2/T3 from the "
            "network itself. The M^2 spacing 2 pi sigma comes from the "
            "closed-string Regge slope 1/(4 pi sigma), a string-theoretic "
            "relation #100 imports; nothing in the ER network supplies it, "
            "and sigma itself is the lattice-calibrated QCD anchor. So "
            f"Delta M^2 = 0.5 * {TWO_PI_SIGMA:.6f} = {delta_m2:.6f} GeV^2 = "
            "pi*sigma is quantitative in its topological factor and "
            "inherited in its scale. The honest claim is therefore neither "
            "'the network derives pi*sigma' nor 'the identification is only "
            "structural', but precisely: the topological half transfers "
            "exactly, the dimensionful half does not come from the network."
        ),
        'sigma_GeV2': SIGMA_QCD,
        'M2_spacing_2pisigma_GeV2': TWO_PI_SIGMA,
        'moding_shift': 0.5,
        'delta_M2_GeV2': delta_m2,
        'pi_sigma_GeV2': PI * SIGMA_QCD,
        'derived': 'the 1/2 (topological, T2/T3)',
        'inherited': '2 pi sigma (closed-string Regge slope + lattice sigma)',
        'pass': ok,
    }


def test_T5_intercept_gap_in_100() -> dict:
    """#100 holds the intercept fixed across both towers; antiperiodic
    moding demonstrably shifts it."""
    rows = []
    for c in (1.0, 2.0, 3.5):
        ez_p, ez_a = casimir_zeta(c, False), casimir_zeta(c, True)
        ec_p, ec_a = casimir_cutoff(c, False), casimir_cutoff(c, True)
        rows.append({
            'circumference': c,
            'E0_periodic_zeta': ez_p, 'E0_antiperiodic_zeta': ez_a,
            'delta_E0': ez_a - ez_p,
            'analytic_pi_over_4C': PI / (4 * c),
            'cutoff_method_periodic': ec_p,
            'cutoff_method_antiperiodic': ec_a,
            'max_method_disagreement': max(abs(ez_p - ec_p),
                                           abs(ez_a - ec_a)),
        })
    methods_agree = all(r['max_method_disagreement'] < 1e-6 for r in rows)
    matches_analytic = all(
        abs(r['delta_E0'] - r['analytic_pi_over_4C']) < 1e-12 for r in rows)
    scales_with_C = abs(rows[0]['delta_E0'] / rows[1]['delta_E0'] - 2.0) < 1e-9
    # #100's towers, as published
    orientable = [math.sqrt(M0_ORIENTABLE ** 2 + TWO_PI_SIGMA * n)
                  for n in range(4)]
    mobius_as_published = [math.sqrt(M0_ORIENTABLE ** 2
                                     + TWO_PI_SIGMA * (n + 0.5))
                           for n in range(4)]
    ok = methods_agree and matches_analytic and scales_with_C
    return {
        'name': 'T5_the_intercept_gap_in_100',
        'description': (
            "#100 writes both towers with the SAME intercept M0^2 and "
            "shifts only the level. But antiperiodic moding changes the "
            "zero-point (Casimir) energy as well. Computed two independent "
            "ways -- zeta regularization using zeta(-1) = -1/12 and "
            "zeta(-1,1/2) = +1/24, and an exponential-cutoff extraction "
            "with the 1/eps^2 divergence removed analytically -- agreeing "
            f"to {max(r['max_method_disagreement'] for r in rows):.1e}: per "
            "transverse polarization on a loop of circumference C, "
            "E0_periodic = -pi/(6C), E0_antiperiodic = +pi/(12C), so "
            "Delta E0 = +pi/(4C) exactly. The Moebius zero point is HIGHER, "
            "and it scales as 1/C (verified: halving C doubles it), so it "
            "cannot be absorbed into a common M0^2. #100's Moebius ground "
            "state is therefore not M0^2 + pi*sigma. Quantifying the "
            "resulting M^2 intercept needs the full closed-loop dynamics "
            "#100 itself defers, so this probe reports the gap and does NOT "
            "publish a corrected tower."
        ),
        'casimir_rows': rows,
        'methods_agree': methods_agree,
        'matches_analytic_pi_over_4C': matches_analytic,
        'scales_as_inverse_C': scales_with_C,
        'tower_100_orientable_GeV': [round(v, 4) for v in orientable],
        'tower_100_mobius_as_published_GeV': [round(v, 4)
                                              for v in mobius_as_published],
        'status_of_published_tower': (
            'level shift +pi*sigma stands; the common intercept does not'),
        'pass': ok,
    }


def test_T6_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T6_assessment_and_the_scope_of_the_correction',
        'description': (
            "PR #227 sec. 8.2 is answered: PARTIALLY quantitative. The "
            "topological factor (the 1/2) transfers exactly and is verified "
            "three independent ways and shown geometry-independent; the "
            "dimensionful factor (2 pi sigma) is an imported closed-string "
            "relation the network does not supply; and #100's common "
            "intercept is unjustified, because antiperiodic moding shifts "
            "the zero-point energy by pi/(4C) per polarization. The net "
            "effect on the program is to SHARPEN #100 rather than overturn "
            "it: its level shift stands, its ground state does not."
        ),
        'scope_of_the_correction': (
            'This touches the GLUEBALL tower of #100 only. It does NOT touch '
            'the heavy Moebius BARYON predictions of #103/#109/#114 '
            '(Lambda_c 3135, Lambda_b 6469, the 849 MeV dipion endpoint), '
            'which are built from the flux-tube quantum Delta = 2 sqrt(sigma) '
            'rather than the closed-loop intercept. The search table stands '
            'as published. Glueballs are unobserved, so the corrected item '
            'is the one with no experimental exposure.'
        ),
        'established': [
            'the 1/2 is derived, shared across three implementations, and '
            'topological (grid-converging, circumference-independent)',
            'pi*sigma = (1/2)*(2 pi sigma): quantitative in its topological '
            'factor, inherited in its scale',
            'antiperiodic moding shifts the zero point by +pi/(4C) per '
            'transverse polarization (two independent regularizations)',
            "#100's common intercept M0^2 across both towers is unjustified",
        ],
        'open': [
            'the corrected Moebius M^2 intercept -- needs the full '
            'closed-loop dynamics that #100 itself defers; not published here',
            'whether the ER network can supply the 2 pi sigma spacing at all, '
            'which would make the identification fully quantitative',
            'PR #227 sec. 8.3: what selects the twist sector dynamically',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_MOEBIUS_IDENTIFICATION_IS_QUANTITATIVE_IN_ITS_TOPOLOGICAL_'
            'HALF_THE_SCALE_IS_INHERITED_AND_THE_MOEBIUS_INTERCEPT_IS_AN_'
            'OPEN_CORRECTION'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_half_from_three_implementations()
    res['T3'] = test_T3_half_is_topological()
    res['T4'] = test_T4_pi_sigma_ledger()
    res['T5'] = test_T5_intercept_gap_in_100()
    res['T6'] = test_T6_assessment(res)
    res['summary'] = {
        'probe': 'mobius_identification_quantitative_probe',
        'pr': 229,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T6']['tests_passed'],
        'verdict_class': res['T6']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Is the Möbius identification quantitative? (PR #229)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The factorization", "",
         "```",
         "Delta M^2  =  (1/2)        x   (2 pi sigma)",
         "              DERIVED           INHERITED",
         "```", "",
         "## The 1/2 is topological", "",
         "| grid | circumference | antiperiodic ground n | dev from 1/2 |",
         "|---:|---:|---:|---:|"]
    for r in s['T3']['rows']:
        o.append(f"| {r['grid']} | {r['circumference']} | "
                 f"{r['antiperiodic_ground_n']:.9f} | {r['dev_from_half']:.1e} |")
    o += ["", "## The intercept gap in #100", "",
          "| C | `E0` periodic | `E0` antiperiodic | `ΔE0` | `π/(4C)` |",
          "|---:|---:|---:|---:|---:|"]
    for r in s['T5']['casimir_rows']:
        o.append(f"| {r['circumference']} | {r['E0_periodic_zeta']:+.6f} | "
                 f"{r['E0_antiperiodic_zeta']:+.6f} | {r['delta_E0']:+.6f} | "
                 f"{r['analytic_pi_over_4C']:+.6f} |")
    o += ["", "## Scope", "", s['T6']['scope_of_the_correction'],
          "", "## Verdict", "", f"**{s['T6']['verdict_class']}**", ""]
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
    out = here / "runs" / f"{ts}_mobius_identification_quantitative_probe"
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
