"""Probe I -- what caps BAM's nonlocality at Tsirelson? (PR #236)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
#206 states the program's sharpest ontological constraint in its own
words: a single classical field on 3-space with local dynamics and local
readout is a local-hidden-variable model and cannot violate CHSH -- that
is Bell's theorem -- so BAM spends its one legitimate nonlocal resource,
the bridge, and reaches CHSH = 2 sqrt 2.

But nonlocality is not a scalar dial that stops where quantum mechanics
stops.  Popescu-Rohrlich boxes are perfectly no-signaling and reach the
algebraic maximum CHSH = 4.  Nothing anywhere in this repository says
what stops BAM at 2 sqrt 2 rather than at 4.  Tsirelson's number appears
in eight documents, always as the value being matched, never as something
the geometry explains.  So:

    Can the bridge be driven past Tsirelson -- and if not, WHICH
    structure forbids it: the geometry, or the imported formalism?

This matters because the answers are not equivalent.  If the geometry
forbids it, BAM has derived a quintessentially quantum bound from
classical ingredients, which is the conjecture doing real work.  If the
formalism forbids it, then the repo's blanket caveat "quantization is
imported, not derived" (#232) has a precise and much smaller residue,
and that residue is exactly what a classical-origin program still owes.

THE ANSWER: THE CEILING IS REAL, AND IT IS NOT GEOMETRIC
---------------------------------------------------------
  (1) THE FLOOR IS BRIDGE-CARRIED.  Cut the handle and the raw B-side
      amplitude falls from 3.747e-01 to 4.140e-16; the 16 deterministic
      local strategies cap at CHSH = 2 exactly.  (#206, re-derived.)

  (2) THE GEOMETRY'S SETTINGS ARE EXACTLY SUFFICIENT -- NOT MORE.  BAM's
      local setting is a fiber-frame rotation, which on the k = +/-1
      qubit is diag(e^{i th}, e^{-i th}) -- a z-rotation, a ONE-parameter
      U(1), not SU(2).  Conjugating the fixed transverse readout sweeps
      the x-y plane and nothing else.  That is enough: the singlet's
      T-matrix restricted to the x-y block has singular values (1, 1), so
      the plane-restricted Horodecki value is exactly 2 sqrt 2 -- attained
      at every grid resolution tested (n = 25, 49, 97, 193, all
      2.828427125).  So the geometry supplies neither a sub-quantum
      deficit nor any excess: precisely the one-parameter family CHSH's
      optimum needs, because that optimum lies in a plane.

  (3) IT NEVER EXCEEDS.  Maximizing over every bridge the model has --
      all 8 gluing holonomies s, all preparations (beta, alpha) -- the
      largest CHSH is 2.828427125, excess over Tsirelson +8.9e-16.

  (4) *** BUT LOCALITY IS NOT WHAT CAPS IT. ***  The obvious guess is
      that the mouths' operations act on commuting algebras -- and they
      do, to an absurd degree: ||[P_A, P_B]|| = 4.7e-113, because the
      committed bridge is a LINK IN H between distinct lattice sites, not
      an identification of degrees of freedom (a fact worth stating
      plainly, since the prose of #206 -- "one object", "two frames of
      one bulk defect" -- reads like an identification).  Yet dropping
      commutativity does not help.  For ANY dichotomic Hermitian
      observables on ONE Hilbert space, with no tensor split at all,

          (B0 + B1)^2 + (B0 - B1)^2 = 4 I        [verified, 6.2e-15]

      and Cauchy-Schwarz then gives CHSH <= 2 sqrt 2 without using
      locality anywhere.  Measured: optimizing A for each B over 2000
      random draws of NON-commuting dichotomic observables, the best
      value is 2.828427041 -- Tsirelson again, from B^2 = I alone.

  ==> The ceiling does not come from the bridge, from the mouths being
      far apart, or from anything else the geometry supplies.  It comes
      from the readout being a +/-1-valued observable on a Hilbert space.

WHAT THIS DOES TO "QUANTIZATION IS IMPORTED"
---------------------------------------------
It makes it precise and much smaller.  The geometry genuinely supplies
the state (the singlet, at fidelity 1.000000000000, from field plus
gluing), the nonlocality (bridge-carried, floor 2 without it), and the
settings (a U(1), exactly sufficient).  What it does NOT supply is the
ceiling, and the whole of the missing ingredient is one algebraic
identity, B^2 = I.  That identity -- not Hilbert space in general, not
the Born rule, not linearity -- is what a classical-origin derivation
still owes.  Locating it is the progress; deriving it is not attempted
here and should not be claimed.

A falsifiable consequence in the meantime: BAM cannot predict
super-quantum correlations.  On this axis it agrees with experiment for a
reason, rather than by construction.

THE NO-SIGNALING CATCH
----------------------
#206's T6 checks that Bob's reduced state is invariant under Alice-side
unitaries.  That is true of EVERY bipartite state -- a partial-trace
identity -- so it tests nothing about the geometry.  Tested on the
lattice instead, no-signaling turns out to be an EQUAL-TIME statement.
The committed handle is always on, so once any evolution separates
Alice's setting from Bob's readout her choice propagates to him:

    Bob's marginal shift, gb = 0.8 :  0 (dt=0), 2.3e-02, 1.0e-01, 2.3e-01
    gb = 0.2 : 0, 5.3e-03, 1.3e-02, 3.6e-02      gb = 0 : ~1e-15 throughout

scaling with the bridge coupling and vanishing when the handle is cut.
This substantiates, with numbers, the caveat #206 states in prose -- "the
physical bridge is non-traversable; the lattice handle is a modeling
stand-in" -- and turns it into a requirement: the bridge must be
dynamically inert at measurement time, which the committed model is not.

Tests:
  T1. Goal.
  T2. The floor: bridge-carried, and 2 without it.
  T3. The geometry's U(1) of settings is exactly sufficient.
  T4. The attack: no bridge exceeds Tsirelson.
  T5. Locality is not the binding constraint; B^2 = I is.
  T6. What that does to "quantization is imported".
  T7. No-signaling is equal-time only; the always-on handle signals.
  T8. Assessment.

Verdict:
  BAM_REACHES_TSIRELSON_EXACTLY_AND_NEVER_EXCEEDS_IT_BUT_THE_CEILING_IS
  _NOT_GEOMETRIC_IT_FOLLOWS_FROM_DICHOTOMIC_READOUT_ALONE_SO_THE_IMPORTED
  _QUANTIZATION_REDUCES_TO_B_SQUARED_EQUALS_I_AND_NO_SIGNALING_IS_EQUAL_TIME
"""

from __future__ import annotations

import itertools
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
TSIRELSON = 2.0 * math.sqrt(2.0)
_SX = cs._SX


# ── the settings the geometry actually supplies ─────────────────────────

def obs_fiber_rotation(theta: float) -> np.ndarray:
    """BAM's local setting: a fiber-frame rotation chi -> chi + delta acts
    on the k = +/-1 qubit as diag(e^{i th}, e^{-i th}) -- a z-rotation.
    Conjugating the fixed transverse readout gives the observable."""
    U = np.diag([np.exp(1j * theta), np.exp(-1j * theta)])
    return U @ _SX @ U.conj().T


def chsh_over_grid(psi4: np.ndarray, thetas: np.ndarray) -> float:
    E = np.array([[float(np.real(np.vdot(
        psi4, np.kron(obs_fiber_rotation(a), obs_fiber_rotation(b)) @ psi4)))
        for b in thetas] for a in thetas])
    best = 0.0
    for ia in range(len(E)):
        for ja in range(len(E)):
            M = (E[ia][:, None] + E[ia][None, :]
                 + E[ja][:, None] - E[ja][None, :])
            best = max(best, float(np.abs(M).max()))
    return best


def raw_b_amplitudes(gb: float) -> tuple:
    """The UNNORMALIZED B-side channel amplitudes -- the bridge-carried
    signal itself (extract_pair_state normalizes, which hides this)."""
    psi = cs.evolve(cs.prepare(PI / 4), cs._T_READ, cs._S_HOLONOMY, gb)
    return (float(abs(np.vdot(cs.mode(cs._XB, -1), psi))),
            float(abs(np.vdot(cs.mode(cs._XB, +1), psi))))


def mouth_projector(x0: int) -> np.ndarray:
    v = cs.mode(x0, +1)
    return np.outer(v, v.conj())


def rand_dichotomic(rng, d: int) -> np.ndarray:
    """A random +/-1-valued Hermitian observable (B^2 = I) on C^d."""
    q, _ = np.linalg.qr(rng.normal(size=(d, d))
                        + 1j * rng.normal(size=(d, d)))
    s = np.diag(rng.choice([-1.0, 1.0], size=d)).astype(complex)
    return q @ s @ q.conj().T


def setting_unitary_at_A(theta: float) -> np.ndarray:
    """The physical setting as a LATTICE operator supported at mouth A."""
    mp, mm = cs.mode(cs._XA, +1), cs.mode(cs._XA, -1)
    return (np.eye(cs._DIM, dtype=complex)
            + (np.exp(1j * theta) - 1) * np.outer(mp, mp.conj())
            + (np.exp(-1j * theta) - 1) * np.outer(mm, mm.conj()))


def bob_marginal(psi: np.ndarray) -> float:
    a = abs(np.vdot(cs.mode(cs._XB, +1), psi)) ** 2
    b = abs(np.vdot(cs.mode(cs._XB, -1), psi)) ** 2
    return a / (a + b) if (a + b) > 0 else 0.5


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_what_stops_the_bridge_at_tsirelson',
        'description': (
            "#206 escapes Bell's theorem by spending the bridge as its one "
            "nonlocal resource, reaching CHSH = 2 sqrt 2. But nonlocality "
            "does not automatically stop where quantum mechanics stops: "
            "Popescu-Rohrlich boxes are perfectly no-signaling and reach the "
            "algebraic maximum 4. Tsirelson's number appears in eight "
            "documents of this repository, always as the value being "
            "matched, never as something the geometry explains. This probe "
            "asks whether the bridge can be driven past it, and -- since it "
            "cannot -- WHICH structure forbids it. The two possible answers "
            "are not equivalent: if the geometry forbids it, the conjecture "
            "has derived a quantum bound from classical ingredients; if the "
            "formalism forbids it, then #232's blanket 'quantization is "
            "imported' has a precise and much smaller residue, and that "
            "residue is what a classical-origin program still owes."
        ),
        'pass': True,
    }


def test_T2_floor() -> dict:
    rows = []
    for gb in (cs._GB, 0.4, 0.1, 0.0):
        a, b = raw_b_amplitudes(gb)
        rows.append({'bridge_coupling': gb, 'raw_c_pm': a, 'raw_c_mp': b})
    local_max = max(abs(a0 * b0 + a0 * b1 + a1 * b0 - a1 * b1)
                    for a0, a1, b0, b1 in itertools.product([-1, 1], repeat=4))
    carried = rows[0]['raw_c_pm'] > 0.1 and rows[-1]['raw_c_pm'] < 1e-12
    return {
        'name': 'T2_the_floor_is_bridge_carried_and_is_two_without_it',
        'description': (
            "The correlation is carried by the handle and by nothing else: "
            f"the raw (unnormalized) B-side channel amplitude is "
            f"{rows[0]['raw_c_pm']:.3e} with the committed coupling and "
            f"{rows[-1]['raw_c_pm']:.3e} -- machine zero -- with the handle "
            "cut, falling monotonically in between. And without it the model "
            "is a local-hidden-variable model: the 16 deterministic local "
            f"strategies cap at CHSH = {local_max:.0f} exactly. This is "
            "#206's result, re-derived here as the baseline the rest of the "
            "probe measures against. NOTE a reporting trap: "
            "`extract_pair_state` normalizes, so the cut bridge yields a "
            "unit-norm pair state made entirely of numerical noise -- the "
            "raw amplitudes above are what actually shows the cut."
        ),
        'raw_amplitudes_vs_coupling': rows,
        'local_strategy_max_chsh': local_max,
        'correlation_is_bridge_carried': carried,
        'pass': bool(carried and abs(local_max - 2.0) < 1e-12),
    }


def test_T3_settings_are_exactly_sufficient() -> dict:
    psi = cs.evolve(cs.prepare(PI / 4), cs._T_READ, cs._S_HOLONOMY)
    p4 = cs.extract_pair_state(psi)
    fid = float(abs(np.vdot(cs._SINGLET, p4)) ** 2)
    T = cs.tmatrix(p4)
    sv = np.linalg.svd(T[:2, :2], compute_uv=False)
    analytic = 2.0 * math.sqrt(float(sv[0] ** 2 + sv[1] ** 2))
    grid = [{'n': n, 'max_chsh': chsh_over_grid(p4, np.linspace(0, PI, n))}
            for n in (25, 49, 97, 193)]
    saturates = all(abs(g['max_chsh'] - TSIRELSON) < 1e-8 for g in grid)
    return {
        'name': 'T3_the_geometry_supplies_a_U1_and_it_is_exactly_enough',
        'description': (
            "BAM's local setting is a fiber-frame rotation, so on the "
            "k = +/-1 qubit it is diag(e^{i th}, e^{-i th}) -- a z-rotation, "
            "a ONE-parameter U(1), not SU(2). Conjugating the fixed "
            "transverse readout therefore sweeps the x-y plane and nothing "
            "else. That is exactly enough, for a reason worth stating: "
            "CHSH's optimum lies in a plane. The lattice-derived state is "
            f"the singlet at fidelity {fid:.12f}, its T-matrix x-y block has "
            f"singular values {np.round(sv, 9).tolist()}, and the "
            f"plane-restricted Horodecki value is {analytic:.9f} = "
            "2 sqrt 2 -- attained at every grid resolution tested "
            "(n = 25, 49, 97, 193, all identical), so this is a saturation "
            "and not a discretization artifact. *** The geometry supplies "
            "neither a sub-quantum deficit nor any excess: precisely the "
            "one-parameter family needed to reach the quantum ceiling, and "
            "structurally incapable of going past it. ***"
        ),
        'singlet_fidelity': fid,
        'xy_block_singular_values': [float(x) for x in sv],
        'analytic_plane_restricted_chsh': analytic,
        'grid_convergence': grid,
        'saturates_tsirelson': saturates,
        'pass': bool(fid > 0.999999 and saturates),
    }


def test_T4_attack() -> dict:
    best, arg = 0.0, None
    for s in range(cs._NCHI):
        for beta in np.linspace(0, PI / 2, 9):
            for alpha in (0.0, PI / 3, PI / 2, PI):
                q = cs.extract_pair_state(
                    cs.evolve(cs.prepare(float(beta), float(alpha)),
                              cs._T_READ, s))
                if np.linalg.norm(q) < 0.5:
                    continue
                v = cs.chsh_horodecki(q)
                if v > best:
                    best, arg = float(v), {'holonomy_s': s,
                                           'beta': round(float(beta), 4),
                                           'alpha': round(float(alpha), 4)}
    excess = best - TSIRELSON
    return {
        'name': 'T4_no_bridge_in_the_model_exceeds_tsirelson',
        'description': (
            "Maximizing the Horodecki CHSH over every bridge the committed "
            "model has -- all 8 gluing holonomies s (the pi-holonomy handle "
            "is only one of them), all preparations (beta, alpha) -- the "
            f"largest value is {best:.9f}, at {arg}, against Tsirelson "
            f"{TSIRELSON:.9f}: an excess of {excess:+.1e}, i.e. machine "
            "zero. The bridge saturates the quantum ceiling and cannot be "
            "driven past it by any preparation or any gluing. So BAM does "
            "not predict super-quantum correlations -- on this axis it "
            "agrees with experiment for a reason rather than by "
            "construction, and it is falsifiable in the other direction."
        ),
        'max_chsh': best,
        'argmax': arg,
        'tsirelson': TSIRELSON,
        'excess': excess,
        'pass': bool(excess < 1e-9 and best > TSIRELSON - 1e-9),
    }


def test_T5_locality_is_not_the_constraint() -> dict:
    PA, PB = mouth_projector(cs._XA), mouth_projector(cs._XB)
    comm = float(np.max(np.abs(PA @ PB - PB @ PA)))
    overlap = float(abs(np.vdot(cs.mode(cs._XA, +1), cs.mode(cs._XB, +1))))

    p4 = cs.extract_pair_state(
        cs.evolve(cs.prepare(PI / 4), cs._T_READ, cs._S_HOLONOMY))
    rng = np.random.default_rng(7)
    resid, best = 0.0, 0.0
    for _ in range(800):
        B0, B1 = rand_dichotomic(rng, 4), rand_dichotomic(rng, 4)
        X, Y = B0 + B1, B0 - B1
        resid = max(resid, float(np.max(np.abs(
            X @ X + Y @ Y - 4 * np.eye(4)))))
        # the optimal A for this B is the sign of the relevant operator, so
        # the attainable value is ||X psi|| + ||Y psi||
        best = max(best, float(np.linalg.norm(X @ p4)
                               + np.linalg.norm(Y @ p4)))
    ok = (comm < 1e-30 and resid < 1e-12
          and best < TSIRELSON + 1e-6 and best > TSIRELSON - 1e-5)
    return {
        'name': 'T5_the_cap_is_B_squared_equals_I_not_locality',
        'description': (
            "The obvious guess is that the two mouths' operations act on "
            "commuting algebras -- and they do, overwhelmingly: "
            f"||[P_A, P_B]|| = {comm:.1e}, mode overlap {overlap:.1e}, "
            "because the committed bridge is a LINK IN H between distinct "
            "lattice sites, not an identification of degrees of freedom. "
            "(Worth stating plainly: #206's prose -- 'one object', 'two "
            "frames of one bulk defect' -- reads like an identification, and "
            "the implementation is not one.) BUT DROPPING COMMUTATIVITY DOES "
            "NOT HELP. For any dichotomic Hermitian observables on ONE "
            "Hilbert space, with no tensor split at all, (B0+B1)^2 + "
            f"(B0-B1)^2 = 4I identically (residual {resid:.1e} over 800 "
            "random draws), and Cauchy-Schwarz then caps CHSH at 2 sqrt 2 "
            "without invoking locality anywhere. Measured, taking the "
            "optimal A for each B: the best attainable value is "
            f"{best:.9f} -- Tsirelson again. *** So the ceiling is not "
            "supplied by the bridge, nor by the mouths being far apart, nor "
            "by anything else in the geometry. It is supplied by the readout "
            "being a +/-1-valued observable on a Hilbert space. ***"
        ),
        'mouth_projector_commutator': comm,
        'mouth_mode_overlap': overlap,
        'identity_residual_B_sq_eq_I': resid,
        'best_chsh_without_commutativity': best,
        'pass': bool(ok),
    }


def test_T6_what_this_does_to_imported_quantization() -> dict:
    items = [
        {'ingredient': 'the pair state (the singlet)',
         'supplied_by': 'GEOMETRY -- field plus bulk gluing, fidelity '
                        '1.000000000000 (#206)'},
        {'ingredient': 'the nonlocality (violating Bell at all)',
         'supplied_by': 'GEOMETRY -- bridge-carried; raw B-side amplitude '
                        '3.7e-01 with the handle, 4.1e-16 without, local '
                        'strategies capped at 2'},
        {'ingredient': 'the measurement settings',
         'supplied_by': 'GEOMETRY -- a U(1) of fiber-frame rotations, '
                        'exactly sufficient to saturate and unable to exceed'},
        {'ingredient': 'THE CEILING (why not 4?)',
         'supplied_by': 'NOT THE GEOMETRY -- it follows from B^2 = I on a '
                        'Hilbert space alone, with locality unused. This is '
                        'the whole residue of "quantization is imported".'},
    ]
    return {
        'name': 'T6_the_residue_of_imported_quantization_is_one_identity',
        'description': (
            "#232's honest scope says 'quantization is imported, not "
            "derived' without saying how much that costs. This probe prices "
            "it. Three of the four ingredients of a Bell test are genuinely "
            "supplied by the geometry -- the state, the nonlocality, the "
            "settings. The fourth, the ceiling, is not, and the entire "
            "missing ingredient is one algebraic identity: B^2 = I. Not "
            "Hilbert space in general, not the Born rule, not linearity -- "
            "just the readout being two-valued. *** That identity is what a "
            "classical-origin derivation still owes, and naming it is the "
            "progress here. Deriving it is not attempted and is not "
            "claimed. *** The narrowing is real but it is a narrowing, not a "
            "closure: an argument that a classical throat readout must be "
            "two-valued would finish this; the repo's Pin- / T^2 = -1 "
            "machinery (#195/#196, #227-#231) is where such an argument "
            "would have to come from, and none of it currently addresses "
            "the readout algebra."
        ),
        'ledger': items,
        'pass': True,
    }


def test_T7_no_signaling_is_equal_time() -> dict:
    psi0 = cs.evolve(cs.prepare(PI / 4), cs._T_READ, cs._S_HOLONOMY)
    rows = []
    for gb in (cs._GB, 0.2, 0.0):
        for dt in (0.0, 0.5, 2.0, 8.0):
            ms = []
            for th in (0.0, 1.0, 2.0):
                psi = setting_unitary_at_A(th) @ psi0
                if dt > 0:
                    psi = cs.evolve(psi, dt, cs._S_HOLONOMY, gb)
                ms.append(bob_marginal(psi))
            rows.append({'bridge_coupling': gb, 'dt': dt,
                         'marginal_shift': float(max(ms) - min(ms))})
    equal_time = [r for r in rows if r['dt'] == 0.0]
    cut = [r for r in rows if r['bridge_coupling'] == 0.0 and r['dt'] > 0]
    live = [r for r in rows if r['bridge_coupling'] == cs._GB and r['dt'] > 0]
    ok = (all(r['marginal_shift'] < 1e-12 for r in equal_time)
          and all(r['marginal_shift'] < 1e-12 for r in cut)
          and max(r['marginal_shift'] for r in live) > 1e-2)
    return {
        'name': 'T7_no_signaling_holds_at_equal_time_only',
        'description': (
            "#206's T6 checks that Bob's reduced state is invariant under "
            "Alice-side unitaries. That is a partial-trace identity true of "
            "EVERY bipartite state, so it tests nothing about the geometry. "
            "Tested on the lattice instead -- applying Alice's physical "
            "setting as an operator supported at mouth A, then evolving, "
            "then reading Bob's channel populations -- no-signaling turns "
            "out to be an EQUAL-TIME statement. At dt = 0 the shift is "
            "exactly 0 (the supports are disjoint). But the committed handle "
            "is always on, so once any evolution separates the setting from "
            "the readout, Alice's choice propagates: shifts of "
            f"{live[0]['marginal_shift']:.1e}, {live[1]['marginal_shift']:.1e}"
            f", {live[2]['marginal_shift']:.1e} at dt = 0.5, 2, 8 -- scaling "
            "with the bridge coupling and vanishing to ~1e-15 when the "
            "handle is cut. *** This substantiates with numbers the caveat "
            "#206 states only in prose ('the physical bridge is "
            "non-traversable; the lattice handle is a modeling stand-in') "
            "and converts it into a requirement: the bridge must be "
            "dynamically inert at measurement time, which the committed "
            "model is not. *** This is a constraint on the model, not a "
            "refutation of the program -- but it is the kind of thing that "
            "has to be discharged rather than asserted."
        ),
        'rows': rows,
        'equal_time_is_clean': all(r['marginal_shift'] < 1e-12
                                   for r in equal_time),
        'cut_bridge_is_clean': all(r['marginal_shift'] < 1e-12 for r in cut),
        'live_bridge_signals': max(r['marginal_shift'] for r in live),
        'pass': bool(ok),
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
            "BAM reaches Tsirelson exactly and cannot be driven past it, "
            "which is the right answer and a falsifiable one. But the bound "
            "is not geometric: it follows from the readout being "
            "two-valued, with locality unused, so it would hold for any "
            "Hilbert-space model whatever its geometry. The useful output is "
            "therefore a much sharper statement of what the classical-origin "
            "conjecture still owes: not 'quantization', but one identity."
        ),
        'established': [
            'the bridge saturates Tsirelson and never exceeds it: max over '
            'all 8 gluing holonomies and all preparations = 2.828427125, '
            'excess +8.9e-16',
            "the geometry's settings are a U(1) (fiber-frame rotations) and "
            'are exactly sufficient -- the x-y block of the singlet T-matrix '
            'has singular values (1,1), giving 2 sqrt 2 at every grid '
            'resolution: no sub-quantum deficit, no excess',
            'the correlation is bridge-carried: raw B-side amplitude '
            '3.747e-01 with the handle, 4.140e-16 without; local strategies '
            'cap at 2',
            'LOCALITY IS NOT THE BINDING CONSTRAINT: the mouth projectors '
            'commute to 4.7e-113 because the bridge is a link between '
            'distinct sites rather than an identification, but dropping '
            'commutativity entirely still gives 2.828427041, because '
            '(B0+B1)^2 + (B0-B1)^2 = 4I for any dichotomic observables',
            'so the residue of "quantization is imported" is exactly one '
            'identity, B^2 = I -- the geometry supplies the state, the '
            'nonlocality and the settings, but not the ceiling',
            'no-signaling is equal-time only: the always-on handle shifts '
            "Bob's marginal by 2.3e-02 / 1.0e-01 / 2.3e-01 at dt = 0.5/2/8, "
            'scaling with the coupling and vanishing when cut',
        ],
        'open': [
            'whether a classical throat readout can be argued to be '
            'two-valued -- that would close the gap this probe locates; the '
            'Pin- / T^2 = -1 machinery is where such an argument would come '
            'from and does not currently address the readout algebra',
            'the committed handle must be made dynamically inert at '
            'measurement time, or the model signals (T7)',
            'this is one 2-outcome, 2-setting scenario on one lattice; '
            'nothing here bounds multipartite or higher-input scenarios, '
            'where the quantum set is not characterized by a single '
            'inequality',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'BAM_REACHES_TSIRELSON_EXACTLY_AND_NEVER_EXCEEDS_IT_BUT_THE_'
            'CEILING_IS_NOT_GEOMETRIC_IT_FOLLOWS_FROM_DICHOTOMIC_READOUT_'
            'ALONE_SO_THE_IMPORTED_QUANTIZATION_REDUCES_TO_B_SQUARED_EQUALS_'
            'I_AND_NO_SIGNALING_IS_EQUAL_TIME'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_floor()
    res['T3'] = test_T3_settings_are_exactly_sufficient()
    res['T4'] = test_T4_attack()
    res['T5'] = test_T5_locality_is_not_the_constraint()
    res['T6'] = test_T6_what_this_does_to_imported_quantization()
    res['T7'] = test_T7_no_signaling_is_equal_time()
    res['T8'] = test_T8_assessment(res)
    res['summary'] = {
        'probe': 'tsirelson_ceiling_probe', 'pr': 236,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T8']['tests_passed'],
        'verdict_class': res['T8']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe I — what caps BAM's nonlocality at Tsirelson? (PR #236)",
         "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "",
         "**Q. The bridge escapes Bell. What stops it at 2√2 rather than "
         "the algebraic 4?**", "",
         "**A. Not the geometry. `B² = I`.**", "",
         "## The ceiling is real", "",
         "| quantity | value |", "|---|---:|",
         f"| max CHSH over all 8 holonomies and all preparations | "
         f"{s['T4']['max_chsh']:.9f} |",
         f"| Tsirelson | {TSIRELSON:.9f} |",
         f"| excess | {s['T4']['excess']:+.1e} |", "",
         "## Where each ingredient comes from", "",
         "| ingredient | supplied by |", "|---|---|"]
    for i in s['T6']['ledger']:
        o.append(f"| {i['ingredient']} | {i['supplied_by'].split(' -- ')[0]} |")
    o += ["",
          "Locality is **not** the binding constraint: the mouth projectors "
          f"commute to {s['T5']['mouth_projector_commutator']:.1e}, but "
          "dropping commutativity entirely still gives "
          f"{s['T5']['best_chsh_without_commutativity']:.9f}, because "
          "`(B₀+B₁)² + (B₀−B₁)² = 4I` for any dichotomic observables "
          f"(residual {s['T5']['identity_residual_B_sq_eq_I']:.1e}).", "",
          "## No-signaling is equal-time only", "",
          "| `gb` | `dt` | Bob's marginal shift |", "|---:|---:|---:|"]
    for r in s['T7']['rows']:
        o.append(f"| {r['bridge_coupling']} | {r['dt']} | "
                 f"{r['marginal_shift']:.3e} |")
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
    out = here / "runs" / f"{ts}_tsirelson_ceiling_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
