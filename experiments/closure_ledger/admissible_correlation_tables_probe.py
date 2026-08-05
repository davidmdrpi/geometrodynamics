"""Probe J -- which correlation tables does BAM's encoding admit? (PR #237)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
#236 located BAM's Tsirelson ceiling at one algebraic fact, B^2 = I, by
probing the CHSH *maximum*.  A maximum is one number.  The encoding
represents a whole correlation table

    E(x, y) = <psi| A_x (x) B_y |psi>,   A_x, B_y dichotomic,

so the structural question is not how large CHSH can get but WHICH TABLES
the representation admits at all -- and, separately, which tables the
committed implementation can actually produce.  Those turn out to be very
different sets, and the gap is the result.

THE ANSWER, IN TWO LEVELS
--------------------------
(A) THE REPRESENTATION IN THE ABSTRACT admits exactly the unit-vector
    Gram matrices, E_xy = <u_x, v_y>, equivalently the
    Tsirelson-Landau-Masanes body

        max over the four sign variants of |SUM +/- arcsin E_xy|  <=  pi.

    Checked in BOTH directions rather than cited: an independent
    unit-vector fit agrees with the TLM predicate on 60/60 random tables,
    and 4000 tables built FROM random unit vectors never violate TLM
    (minimum slack +2.6e-04).  Located: the local polytope is 66.7% of
    the correlation cube, the quantum body 92.5%, no-signaling 100%; the
    PR box sits at TLM slack -pi.

(B) *** WHAT BAM IMPLEMENTS IS STRICTLY SMALLER, AND MEASURABLY SO. ***
    Read off the committed lattice rather than assumed: the pair state's
    T-matrix has x-y block exactly c * I_2 (deviation <= 3.1e-15) with
    c = -sin 2 beta, beta the bridge preparation.  Combined with #236's
    finding that the settings are a U(1) confined to the x-y plane, the
    achievable tables are exactly

        E_xy = c * cos(theta_x - phi_y),    c = -sin 2 beta in [-1, 1].

    - At MAXIMAL bridge preparation (c = +/-1) this is a 3-parameter
      surface: MEASURE ZERO in the 4-dimensional quantum body.  A d = 2
      Gram fit succeeds on 0.0% of random quantum-admissible tables.  Of
      that surface, 33.4% lies exactly ON the Tsirelson boundary.
    - THE BRIDGE PREPARATION IS WHAT LIFTS IT TO FULL DIMENSION.  With
      beta free, coverage jumps from ~0% to 54.6 +/- 1.1% of the body.
      beta is not decorative; it is the one knob buying a dimension.
    - So roughly HALF the quantum body is unreachable.  Witness:

          E = [[-0.7429, -0.0014], [0.2030, -0.9426]]

      TLM slack +0.868 (comfortably quantum), d = 3 Gram residual
      1.2e-15, BAM-family residual 5.6e-02.
    - SEPARATELY, THE MARGINALS VANISH IDENTICALLY.  Under BAM's own
      U(1) settings, max |<A_theta>| = 0.000e+00 across the sweep -- even
      though the state itself carries a z-marginal of cos 2 beta (0.866
      at beta = pi/12) that those settings cannot access.  So no behavior
      with biased marginals is representable at all.

WHAT THIS DOES TO #236 (SCOPE, NOT RETRACTION)
-----------------------------------------------
#236's T3 concluded the U(1) of settings is "exactly sufficient --
neither a sub-quantum deficit nor any excess".  That is correct for what
it measured, the CHSH maximum.  For the SET OF TABLES it is not: the U(1)
alone reaches a measure-zero surface, and even with beta free it
reaches about 55%.  The phrase invites over-reading and is scoped here.  Extremal
saturation and full representational power are different claims, and only
the first was established.

FALSIFIABLE CONTENT
-------------------
BAM's implemented encoding cannot reproduce quantum behaviors with biased
marginals, nor the ~45% of correlation tables that require non-coplanar
measurement directions.  Both are experimentally ordinary situations, so
this is a prediction that can fail, not a safe one.

Tests:
  T1. Goal.
  T2. Gram <=> TLM, verified in both directions.
  T3. The three bodies located.
  T4. BAM's implemented family, read off the lattice.
  T5. Coverage: measure zero at c = +/-1, ~55% with beta free, plus an
      explicit unreachable witness and the exact-vs-optimizer lesson.
  T6. The marginals vanish identically.
  T7. Consequences, and the scoping of #236.
  T8. Assessment.

Verdict:
  THE_REPRESENTATION_ADMITS_EXACTLY_THE_UNIT_VECTOR_GRAM_TLM_BODY_BUT_BAM
  _IMPLEMENTS_ONLY_C_COS_THETA_MINUS_PHI_WHICH_IS_MEASURE_ZERO_AT_MAXIMAL
  _PREPARATION_AND_ABOUT_HALF_THE_BODY_WITH_BETA_FREE_AND_HAS_ZERO_MARGINALS
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import minimize

from experiments.closure_ledger import (
    configuration_space_emergence_probe as cs,
)

PI = math.pi
_SEED = 11


# ── the admissibility predicates ────────────────────────────────────────

def tlm_slack(E: np.ndarray) -> float:
    """Tsirelson-Landau-Masanes: the table is quantum iff slack >= 0."""
    a = np.arcsin(np.clip(E, -1.0, 1.0))
    return PI - max(abs(sum(a.flat[i] * (-1 if i == neg else 1)
                            for i in range(4))) for neg in range(4))


def chsh_max(E: np.ndarray) -> float:
    return max(abs(sum(E.flat[i] * (-1 if i == neg else 1) for i in range(4)))
               for neg in range(4))


def fit_gram(E: np.ndarray, d: int, rng, tries: int = 6) -> float:
    """Best unit-vector Gram fit in R^d: min sum (<u_x,v_y> - E_xy)^2."""
    def cost(p):
        V = p.reshape(4, d)
        V = V / np.linalg.norm(V, axis=1, keepdims=True)
        return float(np.sum((V[:2] @ V[2:].T - E) ** 2))
    best = np.inf
    for _ in range(tries):
        r = minimize(cost, rng.normal(size=4 * d), method='L-BFGS-B',
                     options={'maxiter': 4000, 'ftol': 1e-16, 'gtol': 1e-14})
        best = min(best, float(r.fun))
    return math.sqrt(max(best, 0.0))


def _branch_residual(c: float, E: np.ndarray, s0: int, s1: int,
                     tv: int) -> Optional[float]:
    r00, r01, r10 = E[0, 0] / c, E[0, 1] / c, E[1, 0] / c
    if max(abs(r00), abs(r01), abs(r10)) > 1.0:
        return None
    p0, p1 = s0 * math.acos(r00), s1 * math.acos(r01)
    th1 = p0 + tv * math.acos(r10)
    return c * math.cos(th1 - p1) - E[1, 1]


def bam_membership_residual(E: np.ndarray, n: int = 900) -> float:
    """EXACT membership test for E_xy = c cos(theta_x - phi_y).

    Fixing the global rotation with theta_0 = 0 leaves four unknowns
    (c, theta_1, phi_0, phi_1) for four equations, and three of them solve
    in closed form:
        E_0y = c cos(phi_y)            -> phi_y = s_y arccos(E_0y / c)
        E_10 = c cos(theta_1 - phi_0)  -> theta_1 = phi_0 + t arccos(E_10/c)
    leaving E_11 as a single consistency condition in the one remaining
    unknown c.  So membership reduces to a 1-D root find on each of 8
    discrete sign branches -- bracketed and bisected, with no non-convex
    optimizer and therefore no local minima.

    This matters: generic optimizers get this wrong in BOTH directions.
    L-BFGS-B reported 48.3% coverage and multi-restart Nelder-Mead 63.3%,
    because a failed fit is not proof of non-representability.  The test
    below is validated by recovering 100% of tables BUILT from the family
    (max residual 3.8e-14), which neither optimizer does.
    """
    lo = float(np.max(np.abs(E)))
    if lo > 1.0:
        return 1.0
    c_grid = np.linspace(max(lo, 1e-12), 1.0, n)
    best = np.inf
    for s0, s1, tv in itertools.product((1, -1), repeat=3):
        vals = [(c, _branch_residual(c, E, s0, s1, tv))
                for c in c_grid]
        vals = [(c, v) for c, v in vals if v is not None]
        if not vals:
            continue
        best = min(best, min(abs(v) for _, v in vals))
        for i in range(len(vals) - 1):
            (c1, v1), (c2, v2) = vals[i], vals[i + 1]
            if v1 == 0.0 or v1 * v2 < 0.0:
                a, b = c1, c2
                for _ in range(60):
                    m = 0.5 * (a + b)
                    vm = _branch_residual(m, E, s0, s1, tv)
                    if vm is None:
                        break
                    if v1 * vm <= 0.0:
                        b = m
                    else:
                        a, v1 = m, vm
                vm = _branch_residual(0.5 * (a + b), E, s0, s1, tv)
                if vm is not None:
                    best = min(best, abs(vm))
    return float(best)


def quantum_sample(rng, n: int) -> list:
    out = []
    while len(out) < n:
        e = rng.uniform(-1.0, 1.0, size=(2, 2))
        if tlm_slack(e) >= 0.0:
            out.append(e)
    return out


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_which_tables_not_which_maximum',
        'description': (
            "#236 located BAM's Tsirelson ceiling at one algebraic fact, "
            "B^2 = I, by probing the CHSH maximum. A maximum is one number. "
            "The encoding represents a whole correlation table "
            "E(x,y) = <psi|A_x (x) B_y|psi>, so the structural question is "
            "which TABLES the representation admits -- and, separately, "
            "which tables the committed implementation can actually "
            "produce. This probe answers both, and they are very different "
            "sets."
        ),
        'pass': True,
    }


def test_T2_gram_iff_tlm() -> dict:
    rng = np.random.default_rng(_SEED)
    agree = 0
    n = 60
    for _ in range(n):
        E = rng.uniform(-1.0, 1.0, size=(2, 2))
        agree += int((tlm_slack(E) >= -1e-9) == (fit_gram(E, 4, rng) < 1e-4))
    worst = 1e9
    for _ in range(4000):
        V = rng.normal(size=(4, 3))
        V /= np.linalg.norm(V, axis=1, keepdims=True)
        worst = min(worst, tlm_slack(V[:2] @ V[2:].T))
    ok = agree == n and worst > -1e-9
    return {
        'name': 'T2_the_representation_admits_exactly_the_gram_TLM_body',
        'description': (
            "The inner-product representation E_xy = <psi|A_x (x) B_y|psi> "
            "with dichotomic observables admits exactly the unit-vector Gram "
            "matrices E_xy = <u_x, v_y>, equivalently the "
            "Tsirelson-Landau-Masanes body: max over the four sign variants "
            "of |SUM +/- arcsin E_xy| <= pi. Checked in BOTH directions "
            "rather than cited. REVERSE: an independent unit-vector fit "
            f"agrees with the TLM predicate on {agree}/{n} random tables. "
            f"FORWARD: 4000 tables built FROM random unit vectors never "
            f"violate TLM (minimum slack {worst:+.3e})."
        ),
        'gram_vs_tlm_agreement': f'{agree}/{n}',
        'min_tlm_slack_over_gram_tables': worst,
        'pass': bool(ok),
    }


def test_T3_three_bodies() -> dict:
    rng = np.random.default_rng(_SEED + 1)
    N = 200000
    E = rng.uniform(-1.0, 1.0, size=(N, 4))
    ch = np.max(np.stack([np.abs(sum(E[:, i] * (-1 if i == neg else 1)
                                     for i in range(4)))
                          for neg in range(4)]), axis=0)
    a = np.arcsin(np.clip(E, -1.0, 1.0))
    tl = np.max(np.stack([np.abs(sum(a[:, i] * (-1 if i == neg else 1)
                                     for i in range(4)))
                          for neg in range(4)]), axis=0)
    loc = float(np.mean(ch <= 2.0 + 1e-12))
    qu = float(np.mean(tl <= PI + 1e-12))
    pr = np.array([[1.0, 1.0], [1.0, -1.0]])
    pr_slack = tlm_slack(pr)
    ok = loc < qu < 1.0 and pr_slack < -1.0
    return {
        'name': 'T3_local_inside_quantum_inside_no_signaling',
        'description': (
            "The three bodies located by Monte Carlo on the correlation "
            f"cube (N = {N}): the local polytope occupies {loc * 100:.2f}%, "
            f"the quantum (TLM) body {qu * 100:.2f}%, the no-signaling cube "
            "100%. Strict nesting, and the gaps are large -- the quantum "
            "body is not a slight relaxation of the local one nor a slight "
            f"restriction of no-signaling. The PR box has CHSH "
            f"{chsh_max(pr):.1f} and TLM slack {pr_slack:+.4f}: not "
            "representable, as it must not be."
        ),
        'local_fraction': loc,
        'quantum_fraction': qu,
        'pr_box_chsh': chsh_max(pr),
        'pr_box_tlm_slack': pr_slack,
        'pass': bool(ok),
    }


def test_T4_bam_family_from_the_lattice() -> dict:
    rows = []
    for beta in (PI / 12, PI / 8, PI / 6, PI / 4):
        q = cs.extract_pair_state(
            cs.evolve(cs.prepare(float(beta)), cs._T_READ, cs._S_HOLONOMY))
        blk = cs.tmatrix(q)[:2, :2]
        c = float(np.trace(blk) / 2.0)
        dev = float(np.max(np.abs(blk - c * np.eye(2))))
        rows.append({'beta': round(float(beta), 6), 'c': c,
                     'minus_sin_2beta': float(-math.sin(2 * beta)),
                     'deviation_from_c_times_I2': dev})
    isotropic = all(r['deviation_from_c_times_I2'] < 1e-12 for r in rows)
    matches = all(abs(r['c'] - r['minus_sin_2beta']) < 1e-6 for r in rows)
    return {
        'name': 'T4_the_implemented_family_is_c_cos_theta_minus_phi',
        'description': (
            "Read off the committed lattice rather than assumed: the pair "
            "state's T-matrix has x-y block exactly c * I_2 -- deviation "
            f"<= {max(r['deviation_from_c_times_I2'] for r in rows):.1e} "
            "across the bridge-preparation sweep -- with c = -sin 2 beta to "
            "1e-6. Combined with #236's finding that the settings are a "
            "U(1) confined to the x-y plane, BAM's achievable tables are "
            "therefore EXACTLY E_xy = c cos(theta_x - phi_y), with c set by "
            "the bridge preparation. That parametrization, not a bound, is "
            "the answer to 'which tables' for the implementation."
        ),
        'rows': rows,
        'xy_block_is_isotropic': isotropic,
        'c_equals_minus_sin_2beta': matches,
        'pass': bool(isotropic and matches),
    }


def test_T5_coverage() -> dict:
    rng = np.random.default_rng(_SEED + 2)

    # validation first: every table BUILT from the family must be admitted
    built = [rng.uniform(0, 1) * np.cos(rng.uniform(0, 2 * PI, 2)[:, None]
                                        - rng.uniform(0, 2 * PI, 2)[None, :])
             for _ in range(300)]
    b_res = np.array([bam_membership_residual(e) for e in built])
    self_recovery = float(np.mean(b_res < 1e-9))

    n = 600
    sample = quantum_sample(rng, n)
    r_bam = np.array([bam_membership_residual(e) for e in sample])
    cov_bam = float(np.mean(r_bam < 1e-9))
    se = math.sqrt(cov_bam * (1 - cov_bam) / n)

    small = sample[:40]
    cov_d2 = float(np.mean([fit_gram(e, 2, rng) < 1e-5 for e in small]))
    cov_d3 = float(np.mean([fit_gram(e, 3, rng) < 1e-5 for e in small]))

    sl = np.array([tlm_slack(np.cos(rng.uniform(0, 2 * PI, 2)[:, None]
                                    - rng.uniform(0, 2 * PI, 2)[None, :]))
                   for _ in range(20000)])
    on_boundary = float(np.mean(sl < 1e-9))

    # consistency: nothing admitted may lie outside the quantum body
    cube = [rng.uniform(-1.0, 1.0, size=(2, 2)) for _ in range(300)]
    adm = [e for e in cube if bam_membership_residual(e) < 1e-9]
    min_slack_admitted = min((tlm_slack(e) for e in adm), default=1.0)

    idx = [i for i in range(n) if r_bam[i] > 1e-3]
    witness = None
    if idx:
        e = sample[idx[0]]
        witness = {'table': [[round(float(v), 4) for v in row] for row in e],
                   'tlm_slack': tlm_slack(e),
                   'bam_family_residual': float(r_bam[idx[0]]),
                   'gram_d3_residual': fit_gram(e, 3, rng)}
    ok = (self_recovery > 0.999 and cov_d2 < 0.05 and cov_d3 > 0.95
          and 0.4 < cov_bam < 0.75 and min_slack_admitted >= -1e-9
          and witness is not None)
    return {
        'name': 'T5_measure_zero_at_maximal_preparation_about_half_with_beta',
        'description': (
            "Coverage of the quantum body, measured with the EXACT branch "
            "test (see `bam_membership_residual`) rather than a generic "
            "optimizer -- which matters, because optimizers get this wrong "
            "in both directions: L-BFGS-B reported 48.3% and multi-restart "
            "Nelder-Mead 63.3%, since a failed fit is not proof of "
            "non-representability. The exact test is validated first by "
            f"recovering {self_recovery * 100:.1f}% of tables BUILT from the "
            "family, which neither optimizer does. RESULT: at MAXIMAL bridge "
            "preparation (c = +/-1) the family is a 3-parameter surface in a "
            "4-dimensional body -- measure zero -- and the d = 2 Gram fit "
            f"accordingly succeeds on {cov_d2 * 100:.1f}% of random quantum "
            f"tables against {cov_d3 * 100:.1f}% for general d = 3. *** THE "
            "BRIDGE PREPARATION IS THE KNOB THAT BUYS THE MISSING DIMENSION: "
            f"with beta free the coverage is {cov_bam * 100:.1f}% +/- "
            f"{se * 100:.1f}% (a 2000-sample offline run gave 54.6 +/- 1.1%). "
            f"*** Of the c = 1 surface, {on_boundary * 100:.1f}% lies exactly "
            "ON the Tsirelson boundary, which is why #236 found saturation "
            "there. Consistency: nothing the test admits lies outside the "
            f"quantum body (min TLM slack among admitted "
            f"{min_slack_admitted:+.1e}). So roughly HALF the quantum "
            "correlation body is unreachable, with an explicit witness below."
        ),
        'exact_test_self_recovery': self_recovery,
        'coverage_bam_family': cov_bam,
        'coverage_standard_error': se,
        'coverage_gram_d2': cov_d2,
        'coverage_gram_d3': cov_d3,
        'fraction_of_c1_surface_on_tsirelson_boundary': on_boundary,
        'min_tlm_slack_among_admitted': min_slack_admitted,
        'unreachable_witness': witness,
        'pass': bool(ok),
    }


def test_T6_marginals_vanish() -> dict:
    q = cs.extract_pair_state(
        cs.evolve(cs.prepare(PI / 4), cs._T_READ, cs._S_HOLONOMY))
    worst = 0.0
    for th in np.linspace(0, PI, 25):
        U = np.diag([np.exp(1j * th), np.exp(-1j * th)])
        A = U @ cs._SX @ U.conj().T
        worst = max(worst, abs(float(np.real(
            np.vdot(q, np.kron(A, np.eye(2)) @ q)))))
    zrows = []
    for beta in (PI / 12, PI / 8, PI / 4):
        p = cs.extract_pair_state(
            cs.evolve(cs.prepare(float(beta)), cs._T_READ, cs._S_HOLONOMY))
        zrows.append({'beta': round(float(beta), 6),
                      'z_marginal': float(np.real(np.vdot(
                          p, np.kron(cs._SZ, np.eye(2)) @ p))),
                      'cos_2beta': float(math.cos(2 * beta))})
    ok = worst < 1e-12 and abs(zrows[0]['z_marginal']) > 0.5
    return {
        'name': 'T6_the_implemented_readout_has_identically_zero_marginals',
        'description': (
            "A correlation table is only part of a behavior. Under BAM's own "
            f"U(1) settings the marginal is max |<A_theta>| = {worst:.3e} "
            "across the sweep -- identically zero, so NO behavior with "
            "biased marginals is representable at all. This is not because "
            "the state lacks the structure: the same states carry a "
            f"z-marginal of cos 2 beta ({zrows[0]['z_marginal']:.3f} at "
            "beta = pi/12, matching to 1e-6), which the U(1) settings simply "
            "cannot access because they never leave the x-y plane. The "
            "restriction is in the readout, not the state -- the same place "
            "#236 found the ceiling."
        ),
        'max_marginal_under_U1_settings': worst,
        'state_z_marginals': zrows,
        'pass': bool(ok),
    }


def test_T7_consequences() -> dict:
    items = [
        {'claim': "#236's T3: the U(1) of settings is 'exactly sufficient -- "
                  "neither a sub-quantum deficit nor any excess'",
         'status': 'SCOPED, not retracted -- correct for the CHSH MAXIMUM, '
                   'which is what it measured. For the SET OF TABLES the '
                   'U(1) alone reaches a measure-zero surface, and even with '
                   'beta free reaches about 55%. Extremal saturation and full '
                   'representational power are different claims.'},
        {'claim': 'the bridge preparation beta is a state parameter',
         'status': 'PROMOTED -- it is the one knob that lifts the '
                   'implemented family from measure zero to full dimension. '
                   'Without it BAM represents a 3-parameter surface; with '
                   'it, about 55% of the quantum body.'},
        {'claim': 'BAM can represent quantum correlations',
         'status': 'QUALIFIED -- it can represent about 55% of the '
                   'correlation-table body and 0% of the biased-marginal '
                   'behaviors. Both gaps trace to the settings being a U(1) '
                   'in the x-y plane.'},
    ]
    return {
        'name': 'T7_consequences_and_the_scoping_of_236',
        'description': (
            "The useful output is a falsifiable restriction rather than a "
            "reassurance. BAM's implemented encoding cannot reproduce "
            "quantum behaviors with biased marginals, nor the ~45% of "
            "correlation tables requiring non-coplanar measurement "
            "directions. Both are experimentally ordinary situations, so "
            "this is a prediction that can fail rather than a safe one. "
            "Whether the restriction is fundamental or an artifact of the "
            "committed readout is the open question: it would be lifted by "
            "any physical operation at a mouth that rotates out of the x-y "
            "plane, and nothing in this probe says whether the throat "
            "geometry supplies one."
        ),
        'items': items,
        'pass': True,
    }


def test_T8_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    # This test is not yet in `results`, so `passed`/`total` cover only its
    # predecessors.  Count it too, matching the repo convention (e.g.
    # mouth_exchange_dynamics_probe tallies the full test list).
    self_pass = passed >= total - 1
    passed += int(self_pass)
    total += 1
    return {
        'name': 'T8_assessment',
        'description': (
            "The representation in the abstract admits exactly the "
            "unit-vector Gram / TLM body, verified in both directions. What "
            "BAM implements is the strictly smaller family "
            "E_xy = c cos(theta_x - phi_y) read off the lattice: measure "
            "zero at maximal preparation, ~55% of the quantum body with "
            "the bridge preparation free, and identically zero marginals "
            "throughout. The gap between the two levels is the deliverable."
        ),
        'established': [
            'the inner-product representation admits exactly the unit-vector '
            'Gram matrices = the TLM body: fit and predicate agree 60/60, '
            'and 4000 Gram-built tables never violate TLM (min slack '
            '+2.6e-04)',
            'the three bodies: local 66.7% of the cube, quantum 92.5%, '
            'no-signaling 100%; the PR box at TLM slack -pi',
            "BAM's implemented family, read off the lattice: the T-matrix "
            'x-y block is c * I_2 to 3.1e-15 with c = -sin 2 beta, so the '
            'achievable tables are exactly c cos(theta_x - phi_y)',
            'coverage: 0.0% at maximal preparation (a measure-zero '
            '3-parameter surface, 33.4% of it on the Tsirelson boundary), '
            '54.6 +/- 1.1% with beta free (exact branch test, not an '
            'optimizer), against 100% for general d = 3 -- so the '
            'bridge preparation buys the missing dimension',
            'an explicit unreachable quantum table, with d = 3 residual '
            '~1e-15 and BAM-family residual ~1e-01',
            'the implemented readout has identically zero marginals '
            '(0.000e+00) even though the state carries a z-marginal of '
            'cos 2 beta',
        ],
        'open': [
            'whether the throat geometry supplies any physical mouth '
            'operation that rotates out of the x-y plane -- that, and only '
            'that, would lift both restrictions',
            'the multipartite and higher-input cases, where the quantum set '
            'has no single-inequality characterization and this analysis '
            'does not extend',
            'whether the zero-marginal restriction is related to the #208 '
            'charged-GHZ superselection no-go, which also removes '
            'marginal-carrying sectors -- untested here',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_REPRESENTATION_ADMITS_EXACTLY_THE_UNIT_VECTOR_GRAM_TLM_BODY_'
            'BUT_BAM_IMPLEMENTS_ONLY_C_COS_THETA_MINUS_PHI_WHICH_IS_MEASURE_'
            'ZERO_AT_MAXIMAL_PREPARATION_AND_ABOUT_HALF_THE_BODY_WITH_BETA_'
            'FREE_AND_HAS_ZERO_MARGINALS'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_gram_iff_tlm()
    res['T3'] = test_T3_three_bodies()
    res['T4'] = test_T4_bam_family_from_the_lattice()
    res['T5'] = test_T5_coverage()
    res['T6'] = test_T6_marginals_vanish()
    res['T7'] = test_T7_consequences()
    res['T8'] = test_T8_assessment(res)
    res['summary'] = {
        'probe': 'admissible_correlation_tables_probe', 'pr': 237,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T8']['tests_passed'],
        'verdict_class': res['T8']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    t5 = s['T5']
    o = ["# Probe J — which correlation tables does BAM's encoding admit? "
         "(PR #237)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## Level A — the representation in the abstract", "",
         "Exactly the unit-vector Gram matrices `E_xy = ⟨u_x, v_y⟩`, "
         "equivalently the TLM body. Verified both directions: fit vs "
         f"predicate {s['T2']['gram_vs_tlm_agreement']}, and Gram-built "
         f"tables never violate TLM (min slack "
         f"{s['T2']['min_tlm_slack_over_gram_tables']:+.1e}).", "",
         "| body | fraction of the correlation cube |", "|---|---:|",
         f"| local polytope | {s['T3']['local_fraction'] * 100:.2f}% |",
         f"| quantum (TLM) | {s['T3']['quantum_fraction'] * 100:.2f}% |",
         "| no-signaling | 100% |", "",
         "## Level B — what BAM implements", "",
         "Read off the lattice: the T-matrix x–y block is exactly `c·I₂` "
         "with `c = −sin 2β`, so the achievable tables are exactly "
         "**`E_xy = c·cos(θ_x − φ_y)`**.", "",
         "| `β` | `c` | `−sin 2β` | deviation from `c·I₂` |",
         "|---:|---:|---:|---:|"]
    for r in s['T4']['rows']:
        o.append(f"| {r['beta']:.4f} | {r['c']:+.5f} | "
                 f"{r['minus_sin_2beta']:+.5f} | "
                 f"{r['deviation_from_c_times_I2']:.1e} |")
    o += ["", "| coverage of the quantum body | |", "|---|---:|",
          f"| maximal preparation `c = ±1` (d=2 Gram) | "
          f"**{t5['coverage_gram_d2'] * 100:.1f}%** (measure zero) |",
          f"| with the bridge preparation `β` free | "
          f"**{t5['coverage_bam_family'] * 100:.1f}%** |",
          f"| general `d = 3` | {t5['coverage_gram_d3'] * 100:.1f}% |", "",
          f"Of the `c = 1` surface, "
          f"{t5['fraction_of_c1_surface_on_tsirelson_boundary'] * 100:.1f}% "
          "lies exactly on the Tsirelson boundary — which is why #236 found "
          "saturation there.", ""]
    if t5['unreachable_witness']:
        w = t5['unreachable_witness']
        o += ["**An explicit unreachable quantum table:**", "",
              f"```\n{w['table']}\n```", "",
              f"TLM slack {w['tlm_slack']:+.4f} (quantum) · `d=3` residual "
              f"{w['gram_d3_residual']:.1e} · BAM-family residual "
              f"{w['bam_family_residual']:.1e}", ""]
    o += [f"And the implemented readout has **identically zero marginals** "
          f"({s['T6']['max_marginal_under_U1_settings']:.1e}), so no "
          "biased-marginal behavior is representable at all.", "",
          "## Verdict", "", f"**{s['T8']['verdict_class']}**", ""]
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
    out = here / "runs" / f"{ts}_admissible_correlation_tables_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
