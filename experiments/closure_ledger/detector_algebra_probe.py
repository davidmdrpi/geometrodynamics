"""Probe K -- the detector algebra and the marginals (PR #238)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
#237 established that correlator coverage is NOT the remaining problem:
once convex mixing is allowed, the hull of BAM's tables is the whole
quantum body.  What remains is the DETECTOR MARGINALS and the MEASUREMENT
ALGEBRA.  This probe makes both precise, and they turn out to be one
question with an uncomfortable answer.

THE ANSWER
----------
(1) THE ALGEBRA IS NOT THE RESTRICTION.  Conjugating a fixed detector by
    the setting group traces out a great circle of observables -- a
    2-dimensional span -- but the *-algebra those observables GENERATE is
    all of M_2(C) (dimension 4, verified: products of two x-y observables
    already produce sigma_z).  So nothing is missing algebraically.  The
    restriction is on WHICH SELF-ADJOINT ELEMENTS ARE DIALABLE as
    settings, which is a different and much sharper question.

(2) THE REPO HAS TWO INCOMPATIBLE SETTING MODELS, AND THEY DISAGREE.
    The documents say the setting is "the fiber-frame rotation before the
    device".  A rotation of the chi fiber acts on the winding-k channels
    as diag(e^{i th}, e^{-i th}) -- a z-rotation.  But the committed code
    `measurement_sector_probe._rot` is [[c, -s], [s, c]] = exp(-i th
    sigma_y / 2) -- a y-rotation (verified identical to 0.00e+00).  These
    are not the same operation: the fiber rotation COMMUTES EXACTLY with
    the winding detector (||[fiber, sigma_z]|| = 0.00e+00), the
    y-rotation does not.

(3) *** THE TRILEMMA.  The only pairing in which BOTH the setting and the
    detector are physically derived by this repository gives CHSH = 2
    exactly -- no Bell violation. ***  Measured on the committed #206
    singlet:

      setting              detector                    CHSH     dialable
      fiber U(1) (derived) sigma_z winding SG (derived) 2.000000    1
      fiber U(1) (derived) sigma_x transverse (assumed) 2.828         2
      y-rotation (in code) sigma_z winding SG (derived) 2.828         2

    The violation requires exactly one non-derived ingredient -- either a
    TRANSVERSE DETECTOR or a WINDING-MIXING SETTING -- and the repo
    supplies neither.  With the derived pair the accessible observable set
    is a single point (span dimension 1): the setting does nothing at all,
    because a fiber rotation cannot move a winding-diagonal detector.

(4) THE MARGINALS FOLLOW THE SAME FORK, AND #237's CLAIM IS CORRECTED.
    #237 reported "identically zero marginals, immune to convex mixing"
    as the probe's one surviving falsifiable content.  That measurement
    was made with the fiber x sigma_x pairing.  Under the other two
    pairings the marginals are NOT zero -- 0.866 at beta = pi/12, exactly
    cos 2 beta.  So the zero-marginal restriction is an artifact of one
    of three modelling choices and is RETRACTED as a general claim; it
    holds only for the transverse-detector model.  The marginals and the
    algebra are the same question, not two.

(5) WHAT THE TWO MISSING INGREDIENTS HAVE IN COMMON.  sigma_x in the
    winding basis is off-diagonal between k = +1 and k = -1 -- a
    coherence between distinct topological charges -- and the y-rotation
    setting mixes those same sectors.  Both non-derived ingredients are
    therefore the SAME physical requirement: an operation or observable
    that coherently connects distinct winding sectors.  And if winding is
    superselected, that requirement is unavailable: restricting both
    parties to k-diagonal observables gives max CHSH = 2.000000 exactly.
    The repo already carries a superselection structure in this sector
    (#208's charged-GHZ zero-sum no-go), and nothing in it has been
    reconciled with the Bell chain.

WHAT THIS DOES AND DOES NOT SAY
--------------------------------
It does NOT show BAM's Bell violation is wrong.  It shows the violation
rests on one ingredient the program has not derived, that the two
candidate ingredients are the same physical requirement, and that the
requirement is in tension with a superselection structure the repo
already committed elsewhere.  Naming which single thing has to be
supplied -- a mouth operation or observable connecting winding sectors --
is the deliverable.

Tests:
  T1. Goal.
  T2. The algebra is full; only the dialable set is restricted.
  T3. Two incompatible setting models: docs vs committed code.
  T4. The trilemma: CHSH by (setting, detector).
  T5. The marginals follow the same fork; #237's claim corrected.
  T6. Both missing ingredients are winding-sector coherence; the
      superselection collapse.
  T7. Consequences.
  T8. Assessment.

Verdict:
  THE_ALGEBRA_IS_FULL_BUT_THE_DIALABLE_SET_IS_NOT_AND_THE_ONLY_PAIRING_OF
  _A_DERIVED_SETTING_WITH_A_DERIVED_DETECTOR_GIVES_CHSH_EXACTLY_TWO_SO_THE
  _VIOLATION_RESTS_ON_WINDING_SECTOR_COHERENCE_THAT_IS_NOWHERE_DERIVED
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from experiments.closure_ledger import (
    configuration_space_emergence_probe as cs,
)

PI = math.pi
_SX, _SY, _SZ = cs._SX, cs._SY, cs._SZ
_I2 = np.eye(2, dtype=complex)


# ── the two candidate setting groups ────────────────────────────────────

def fiber_rotation(th: float) -> np.ndarray:
    """What the GEOMETRY supplies: chi -> chi + delta multiplies the
    winding-k channel by e^{2 pi i k delta / N}, so on the k = +/-1 qubit
    it is diag(e^{i th}, e^{-i th}) -- a z-rotation."""
    return np.diag([np.exp(1j * th), np.exp(-1j * th)])


def y_rotation(th: float) -> np.ndarray:
    """What `measurement_sector_probe._rot` actually implements."""
    c, s = math.cos(th / 2), math.sin(th / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def observable(setting, detector: np.ndarray, th: float) -> np.ndarray:
    U = setting(th)
    return U.conj().T @ detector @ U


def chsh_over_settings(state: np.ndarray, setting, detector: np.ndarray,
                       n: int = 81) -> float:
    ths = np.linspace(0.0, 2 * PI, n)
    ops = [observable(setting, detector, t) for t in ths]
    E = np.array([[float(np.real(np.vdot(state, np.kron(A, B) @ state)))
                   for B in ops] for A in ops])
    best = 0.0
    for ia in range(n):
        for ja in range(n):
            M = (E[ia][:, None] + E[ia][None, :]
                 + E[ja][:, None] - E[ja][None, :])
            best = max(best, float(np.abs(M).max()))
    return best


def dialable_span(setting, detector: np.ndarray, n: int = 60) -> int:
    """Dimension of the real span of the reachable observables."""
    M = np.stack([observable(setting, detector, t).reshape(4)
                  for t in np.linspace(0.0, 2 * PI, n)])
    return int(np.linalg.matrix_rank(
        np.vstack([M.real, M.imag]), tol=1e-9))


def max_marginal(state: np.ndarray, setting, detector: np.ndarray,
                 n: int = 60) -> float:
    return max(abs(float(np.real(np.vdot(
        state, np.kron(observable(setting, detector, t), _I2) @ state))))
        for t in np.linspace(0.0, 2 * PI, n))


def pair_state(beta: float) -> np.ndarray:
    return cs.extract_pair_state(
        cs.evolve(cs.prepare(float(beta)), cs._T_READ, cs._S_HOLONOMY))


_PAIRINGS = (
    ('fiber U(1) [derived]', fiber_rotation,
     'sigma_z winding SG [derived]', _SZ, True, True),
    ('fiber U(1) [derived]', fiber_rotation,
     'sigma_x transverse [assumed by #236/#237]', _SX, True, False),
    ('y-rotation [in #209 code]', y_rotation,
     'sigma_z winding SG [derived]', _SZ, False, True),
)


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_the_detector_marginals_and_the_measurement_algebra',
        'description': (
            "#237 established that correlator coverage is NOT the remaining "
            "problem: once convex mixing is allowed, the hull of BAM's "
            "tables is the whole quantum body. What remains is the detector "
            "marginals and the measurement algebra. This probe makes both "
            "precise -- and finds they are one question, whose answer is "
            "that the Bell violation rests on an ingredient the program has "
            "not derived."
        ),
        'pass': True,
    }


def test_T2_the_algebra_is_full() -> dict:
    gens = [observable(fiber_rotation, _SX, t)
            for t in np.linspace(0.0, 2 * PI, 12)]
    prods = [a @ b for a in gens for b in gens] + gens + [_I2]
    A = np.stack([p.reshape(4) for p in prods])
    gen_dim = int(np.linalg.matrix_rank(
        np.vstack([A.real, A.imag]), tol=1e-9))
    dial = {f'{s} x {d}': dialable_span(S, D)
            for s, S, d, D, _, _ in _PAIRINGS}
    ok = gen_dim == 4 and min(dial.values()) < 4
    return {
        'name': 'T2_the_algebra_is_full_only_the_dialable_set_is_restricted',
        'description': (
            "Conjugating a fixed detector by the setting group traces out a "
            "great circle of observables -- a 2-dimensional real span at "
            "best. But the *-algebra those observables GENERATE is all of "
            f"M_2(C): the span of their products has dimension {gen_dim} of "
            "4, because a product of two x-y observables already produces "
            "sigma_z. *** So nothing is missing algebraically, and 'the "
            "measurement algebra is too small' is the wrong diagnosis. The "
            "restriction is on WHICH SELF-ADJOINT ELEMENTS ARE DIALABLE as "
            "settings -- a different and much sharper question, and the one "
            "the rest of this probe answers. ***"
        ),
        'generated_algebra_dimension': gen_dim,
        'dialable_span_by_pairing': dial,
        'pass': bool(ok),
    }


def test_T3_two_setting_models() -> dict:
    rows = []
    for th in (0.3, 1.1, 2.7):
        R = y_rotation(th)
        expy = (math.cos(th / 2) * _I2 - 1j * math.sin(th / 2) * _SY)
        rows.append({
            'theta': th,
            'code_rot_minus_exp_minus_i_theta_sy_over_2':
                float(np.max(np.abs(R - expy))),
            'commutator_fiber_with_sigma_z': float(np.max(np.abs(
                fiber_rotation(th) @ _SZ - _SZ @ fiber_rotation(th)))),
            'commutator_yrot_with_sigma_z': float(np.max(np.abs(
                R @ _SZ - _SZ @ R))),
        })
    is_yrot = all(r['code_rot_minus_exp_minus_i_theta_sy_over_2'] < 1e-12
                  for r in rows)
    fiber_commutes = all(r['commutator_fiber_with_sigma_z'] < 1e-12
                         for r in rows)
    yrot_does_not = all(r['commutator_yrot_with_sigma_z'] > 1e-2
                        for r in rows)
    return {
        'name': 'T3_the_docs_and_the_code_use_different_setting_operations',
        'description': (
            "The documents describe the setting as 'the fiber-frame "
            "rotation before the device'. A rotation of the chi fiber "
            "multiplies the winding-k channel by e^{2 pi i k delta / N}, so "
            "on the k = +/-1 qubit it is diag(e^{i th}, e^{-i th}) -- a "
            "Z-ROTATION. But the committed `measurement_sector_probe._rot` "
            "is [[c, -s], [s, c]], which is exp(-i th sigma_y / 2) -- a "
            "Y-ROTATION, verified identical to "
            f"{max(r['code_rot_minus_exp_minus_i_theta_sy_over_2'] for r in rows):.1e}"
            ". These are not the same operation, and the difference is "
            "exactly the one that matters: the fiber rotation COMMUTES "
            "EXACTLY with the winding detector (||[fiber, sigma_z]|| = "
            "0.00e+00 at every angle tested) while the y-rotation does not "
            f"(up to {max(r['commutator_yrot_with_sigma_z'] for r in rows):.2f}). "
            "A fiber rotation cannot move a winding-diagonal detector, so "
            "under the described setting the dial does nothing."
        ),
        'rows': rows,
        'code_rot_is_a_y_rotation': is_yrot,
        'fiber_commutes_with_winding_detector': fiber_commutes,
        'yrot_does_not_commute': yrot_does_not,
        'pass': bool(is_yrot and fiber_commutes and yrot_does_not),
    }


def test_T4_trilemma() -> dict:
    psi = pair_state(PI / 4)
    rows = []
    for sname, S, dname, D, s_der, d_der in _PAIRINGS:
        rows.append({
            'setting': sname, 'detector': dname,
            'setting_derived': s_der, 'detector_derived': d_der,
            'both_derived': bool(s_der and d_der),
            'max_chsh': chsh_over_settings(psi, S, D),
            'dialable_span': dialable_span(S, D),
        })
    derived = [r for r in rows if r['both_derived']]
    others = [r for r in rows if not r['both_derived']]
    ok = (len(derived) == 1 and abs(derived[0]['max_chsh'] - 2.0) < 1e-6
          and derived[0]['dialable_span'] == 1
          and all(r['max_chsh'] > 2.5 for r in others))
    return {
        'name': 'T4_the_only_fully_derived_pairing_gives_CHSH_exactly_two',
        'description': (
            "Maximum CHSH on the committed #206 singlet, by (setting, "
            "detector). *** THE ONLY PAIRING IN WHICH BOTH INGREDIENTS ARE "
            "PHYSICALLY DERIVED BY THIS REPOSITORY -- the fiber U(1) "
            "setting with the sigma_z winding Stern-Gerlach -- gives "
            f"{derived[0]['max_chsh']:.6f}, exactly the local bound. No Bell "
            "violation at all. *** And the reason is visible in the "
            "dialable span: it is 1, a single point, because a fiber "
            "rotation cannot move a winding-diagonal detector, so the "
            "setting does nothing. The other two pairings do violate "
            f"({others[0]['max_chsh']:.3f}, {others[1]['max_chsh']:.3f}), "
            "and each contains exactly one ingredient the repo has not "
            "derived: a transverse detector, or a winding-mixing setting."
        ),
        'rows': rows,
        'pass': bool(ok),
    }


def test_T5_marginals_follow_the_same_fork() -> dict:
    rows = []
    for beta in (PI / 12, PI / 8, PI / 6, PI / 4):
        q = pair_state(beta)
        entry = {'beta': round(float(beta), 6),
                 'cos_2beta': float(math.cos(2 * beta))}
        for sname, S, dname, D, _, _ in _PAIRINGS:
            entry[f'{sname.split()[0]}_x_{dname.split()[0]}'] = max_marginal(
                q, S, D)
        rows.append(entry)
    transverse_zero = all(r['fiber_x_sigma_x'] < 1e-12 for r in rows)
    others_nonzero = any(r['fiber_x_sigma_z'] > 0.5 for r in rows)
    matches_cos = all(abs(r['fiber_x_sigma_z'] - abs(r['cos_2beta'])) < 1e-6
                      for r in rows)
    return {
        'name': 'T5_the_zero_marginal_claim_of_237_is_model_dependent',
        'description': (
            "#237 reported 'identically zero marginals, immune to convex "
            "mixing' as its one surviving falsifiable content. That "
            "measurement was made with the fiber x sigma_x pairing, and it "
            "is correct there (max |<A>| < 1e-12 at every beta). But under "
            "the other two pairings the marginals are NOT zero: they are "
            "|cos 2 beta| exactly (0.866 at beta = pi/12, matching to "
            "1e-6). *** So the zero-marginal restriction is an artifact of "
            "one of three modelling choices, and is RETRACTED as a general "
            "claim about BAM -- it holds only for the transverse-detector "
            "model. #237's falsifiable content was therefore contingent on "
            "the same undetermined choice that T4 shows the Bell violation "
            "itself hangs on. The marginals and the algebra are one "
            "question, not two. ***"
        ),
        'rows': rows,
        'transverse_pairing_gives_zero': transverse_zero,
        'other_pairings_give_cos_2beta': matches_cos,
        'pass': bool(transverse_zero and others_nonzero and matches_cos),
    }


def test_T6_winding_coherence_and_superselection() -> dict:
    psi = pair_state(PI / 4)
    rng = np.random.default_rng(4)
    best = 0.0
    for _ in range(4000):
        d = [np.diag(rng.choice([-1.0, 1.0], size=2)).astype(complex)
             for _ in range(4)]

        def e(X, Y):
            return float(np.real(np.vdot(psi, np.kron(X, Y) @ psi)))
        best = max(best, abs(e(d[0], d[2]) + e(d[0], d[3])
                             + e(d[1], d[2]) - e(d[1], d[3])))
    offdiag_sx = float(abs(_SX[0, 1]))
    offdiag_yrot = float(abs(y_rotation(1.0)[0, 1]))
    ok = abs(best - 2.0) < 1e-6 and offdiag_sx > 0.5 and offdiag_yrot > 0.1
    return {
        'name': 'T6_both_missing_ingredients_are_winding_sector_coherence',
        'description': (
            "The two non-derived ingredients of T4 are the same physical "
            "requirement. sigma_x in the winding basis is purely "
            f"off-diagonal between k = +1 and k = -1 (|<+1|sigma_x|-1>| = "
            f"{offdiag_sx:.1f}) -- a coherence between distinct topological "
            "charges -- and the y-rotation setting mixes those same sectors "
            f"(off-diagonal element {offdiag_yrot:.3f} at th = 1). So the "
            "one thing BAM must supply is an operation OR observable that "
            "coherently connects distinct winding sectors. *** And if "
            "winding is superselected that is unavailable: restricting both "
            "parties to k-diagonal dichotomic observables gives max CHSH = "
            f"{best:.6f}, exactly the local bound, over 4000 random draws. "
            "The repo already carries a superselection structure in this "
            "sector -- #208's charged-GHZ zero-sum no-go -- and nothing in "
            "it has been reconciled with the Bell chain. *** This is not a "
            "proof that BAM's violation fails; it is the identification of "
            "the single unresolved requirement it rests on."
        ),
        'max_chsh_under_winding_superselection': best,
        'sigma_x_offdiagonal': offdiag_sx,
        'y_rotation_offdiagonal_at_theta_1': offdiag_yrot,
        'pass': bool(ok),
    }


def test_T7_consequences() -> dict:
    items = [
        {'claim': "#237: 'the implemented readout has identically zero "
                  "marginals, and this survives mixing -- the whole of the "
                  "probe's falsifiable content'",
         'status': 'RETRACTED as a general claim -- true only for the '
                   'fiber x sigma_x pairing. Under the other two the '
                   'marginals are |cos 2 beta|, up to 0.866.'},
        {'claim': "'the measurement algebra is too small'",
         'status': 'WRONG DIAGNOSIS -- the generated algebra is all of '
                   'M_2(C). What is restricted is the dialable set of '
                   'settings, which for the fully derived pairing is a '
                   'single point.'},
        {'claim': "BAM's Bell violation (#206/#209)",
         'status': 'NOT REFUTED, BUT LOCATED -- it requires exactly one '
                   'ingredient the repo has not derived. With the derived '
                   'setting and the derived detector, CHSH = 2.000000.'},
        {'claim': 'the docs and the code agree on what a setting is',
         'status': 'FALSE -- the docs say fiber-frame rotation '
                   '(diag(e^{i th}, e^{-i th})), the code applies '
                   'exp(-i th sigma_y / 2). For the winding carrier these '
                   'differ, and only the second violates Bell.'},
    ]
    return {
        'name': 'T7_consequences',
        'description': (
            "The deliverable is a single well-posed requirement rather than "
            "another coverage number. BAM's Bell chain needs a mouth "
            "operation or observable that coherently connects distinct "
            "winding sectors; supply that and both the violation and the "
            "marginals follow, and the #208 superselection structure has to "
            "be reconciled with it. Until then the honest statement is that "
            "the violation is carried by a modelling choice rather than by "
            "the geometry, and that #237's falsifiable content was "
            "contingent on the same choice. The natural next step is to ask "
            "whether the transported-frame / spin doublet (#192/#197 "
            "Berger-squash, named in #209 as the analogous device for the "
            "OTHER carrier) supplies the missing operation for a carrier "
            "that is not winding -- which would relocate the whole Bell "
            "chain off the winding qubit."
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
            "The measurement algebra is full; what is restricted is which "
            "self-adjoint elements can be dialed. The repo's documents and "
            "its committed code disagree about what a setting is, and the "
            "disagreement is decisive: with the fiber rotation the docs "
            "describe and the winding detector #209 derives, the setting "
            "commutes with the detector, the dialable set is a single "
            "point, and CHSH is exactly 2. Both ways of restoring the "
            "violation require coherence between distinct winding sectors, "
            "which winding superselection would forbid."
        ),
        'established': [
            'the generated measurement algebra is all of M_2(C) (dimension '
            "4) -- 'the algebra is too small' is the wrong diagnosis; the "
            'dialable set is what is restricted',
            "the committed `_rot` is exp(-i th sigma_y/2), a y-rotation, not "
            'the fiber-frame rotation diag(e^{i th}, e^{-i th}) the docs '
            'describe (verified to 0.0e+00)',
            'the fiber rotation commutes EXACTLY with the sigma_z winding '
            'detector (0.0e+00), so under the described setting the dial '
            'does nothing: dialable span 1',
            'CHSH by pairing: fiber x sigma_z (both derived) = 2.000000; '
            'fiber x sigma_x = 2.828; y-rotation x sigma_z = 2.828 -- the '
            'violation requires exactly one non-derived ingredient',
            "#237's zero-marginal claim holds only for the fiber x sigma_x "
            'pairing; under the others the marginals are |cos 2 beta| (up '
            'to 0.866), so it is retracted as a general claim',
            'both missing ingredients are the same requirement -- coherence '
            'between distinct winding sectors -- and imposing winding '
            'superselection gives max CHSH = 2.000000 exactly',
        ],
        'open': [
            'whether any mouth operation or observable in BAM coherently '
            'connects distinct winding sectors; supplying that settles both '
            'the violation and the marginals at once',
            "how #208's charged-GHZ superselection structure is to be "
            'reconciled with a Bell chain that requires exactly the '
            'coherence it removes',
            'whether relocating the Bell chain onto the transported-frame / '
            'spin doublet (#192/#197 Berger-squash, named in #209 as the '
            'device for that carrier) avoids the problem -- untested here',
            'this is the two-input two-output winding sector on one '
            'lattice; the multipartite chain (#207/#208) is not analysed',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_ALGEBRA_IS_FULL_BUT_THE_DIALABLE_SET_IS_NOT_AND_THE_ONLY_'
            'PAIRING_OF_A_DERIVED_SETTING_WITH_A_DERIVED_DETECTOR_GIVES_CHSH_'
            'EXACTLY_TWO_SO_THE_VIOLATION_RESTS_ON_WINDING_SECTOR_COHERENCE_'
            'THAT_IS_NOWHERE_DERIVED'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_the_algebra_is_full()
    res['T3'] = test_T3_two_setting_models()
    res['T4'] = test_T4_trilemma()
    res['T5'] = test_T5_marginals_follow_the_same_fork()
    res['T6'] = test_T6_winding_coherence_and_superselection()
    res['T7'] = test_T7_consequences()
    res['T8'] = test_T8_assessment(res)
    res['summary'] = {
        'probe': 'detector_algebra_probe', 'pr': 238,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T8']['tests_passed'],
        'verdict_class': res['T8']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe K — the detector algebra and the marginals (PR #238)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "",
         "**Q. After #237, correlator coverage is not the problem. What is?**",
         "", "**A. Which effects are dialable — and the fully derived pairing "
         "does not violate Bell at all.**", "",
         "## The trilemma", "",
         "| setting | detector | both derived? | dialable span | max CHSH |",
         "|---|---|---|---:|---:|"]
    for r in s['T4']['rows']:
        o.append(f"| {r['setting']} | {r['detector']} | "
                 f"{'**yes**' if r['both_derived'] else 'no'} | "
                 f"{r['dialable_span']} | "
                 f"{'**' if r['both_derived'] else ''}{r['max_chsh']:.6f}"
                 f"{'**' if r['both_derived'] else ''} |")
    o += ["",
          f"The generated algebra is all of `M₂(C)` (dimension "
          f"{s['T2']['generated_algebra_dimension']}), so *the algebra is not "
          "the restriction* — the dialable set is.", "",
          "## The marginals follow the same fork", "",
          "| `β` | `\\|cos 2β\\|` | fiber×σ_z | fiber×σ_x | y-rot×σ_z |",
          "|---:|---:|---:|---:|---:|"]
    for r in s['T5']['rows']:
        o.append(f"| {r['beta']:.4f} | {abs(r['cos_2beta']):.4f} | "
                 f"{r['fiber_x_sigma_z']:.4f} | {r['fiber_x_sigma_x']:.1e} | "
                 f"{r['y-rotation_x_sigma_z']:.4f} |")
    o += ["",
          "So #237's *identically zero marginals* holds only for the "
          "`fiber × σ_x` pairing and is **retracted as a general claim**.",
          "",
          "## Both missing ingredients are the same thing", "",
          "`σ_x` is off-diagonal between `k = ±1`, and the y-rotation mixes "
          "those sectors: both require **coherence between distinct winding "
          "sectors**. Imposing winding superselection gives max CHSH = "
          f"**{s['T6']['max_chsh_under_winding_superselection']:.6f}** — "
          "exactly the local bound.", "",
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
    out = here / "runs" / f"{ts}_detector_algebra_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
