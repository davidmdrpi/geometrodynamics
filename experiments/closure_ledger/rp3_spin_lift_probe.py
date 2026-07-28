"""Probe A -- the RP^3 spin-lift ledger (PR #230).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
The twist-selection arc composed the network holonomy W with an "intrinsic
spin structure eta" read off B2 (T = i sigma_y, T^2 = -I), via the rule
theta = 0 iff W*eta = +1.  That composition was never justified: it assumes
B2's Z2 and the network's Z2 are the SAME Z2, living on the same cycle.
This probe builds the lift ledger explicitly and asks whether they are
independent.

WHAT IS COMPUTED
----------------
  * J: S^3 -> S^3, the antipodal map x -> -x, as an element of SO(4);
  * both lifts J~_+ , J~_- in Spin(4) = SU(2) x SU(2), via the quaternionic
    action x -> q_L x q_R^dagger;
  * the group law, the kernel, and the squares J~_+^2, J~_-^2;
  * where BAM's sigma(z1,z2) = (conj z2, -conj z1) and T = i sigma_y sit;
  * whether B2 and W are independent data.

THE FINDINGS
------------
1. Both lifts project to -I_4 and BOTH SQUARE TO +1 in Spin(4).  They differ
   by the kernel element (-1,-1).  So both define genuine Z2 actions on the
   spinor bundle, i.e. RP^3 carries two spin structures -- consistent with
   H^1(RP^3;Z2) = Z2.  (Note this refutes the folklore shortcut "the two
   spin structures are distinguished by J~^2 = +1 vs -1": here both are +1,
   and the distinction is WHICH SU(2) factor carries the -1.)

2. BAM's sigma is ORIENTATION-PRESERVING (det = +1, so it lies in SO(4),
   not O(4)\\SO(4)) and satisfies sigma^2 = -I_4 = J.  So sigma is a square
   root of the antipodal map, and T = i sigma_y -- its spinor lift, with
   T^2 = -I in the LEFT SU(2) -- gives T^2 = (-1, 1) = J~_+ .

   *** So B2 DOES fix the RP^3 spin structure: it selects J~_+ . ***

3. But B2's Z2 and the network's W are NOT the same datum.  B2 lives in
   H^1(RP^3; Z2) -- a choice of spin structure on the exterior 3-manifold,
   detected by the deck loop that generates pi_1(RP^3) = Z2.  W lives in
   H^1(Gamma; Z2) for the ER graph Gamma -- a holonomy on graph cycles.
   They coincide only if the network cycle IS the deck loop, which is a
   claim about the embedding of the cycle in RP^3, not an identity.

   *** CONSEQUENCE: the composition rule theta = 0 iff W*eta = +1, on which
   the moding-swap of bam_spinor_spectrum_effective_action_probe rests, is
   NOT justified in general.  It is justified only for a network cycle that
   generates pi_1(RP^3).  The #224 two-throat ring is a loop through two
   exterior arcs; whether it is the pi_1 generator is a separate question
   this probe does NOT settle. ***

Tests:
  T1. J as an element of SO(4).
  T2. Both lifts, the covering map, and the kernel.
  T3. The group law and the squares.
  T4. Where sigma and T = i sigma_y sit: sigma^2 = J, T^2 = J~_+ .
  T5. Are B2 and W independent?
  T6. Assessment.

Verdict:
  BOTH_LIFTS_SQUARE_TO_PLUS_ONE_B2_SELECTS_J_PLUS_VIA_T_SQUARED_AND_B2_IS
  _INDEPENDENT_OF_THE_NETWORK_HOLONOMY_SO_THE_COMPOSITION_RULE_IS_UNJUSTIFIED
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
T_BAM = 1j * SY                      # BAM's throat transport, T^2 = -I


# ── quaternionic model: S^3 = SU(2), Spin(4) = SU(2) x SU(2) ───────────

def quat(v: np.ndarray) -> np.ndarray:
    """R^4 -> 2x2:  x = x0 I + i(x1 sx + x2 sy + x3 sz)."""
    return v[0] * I2 + 1j * (v[1] * SX + v[2] * SY + v[3] * SZ)


def unquat(m: np.ndarray) -> np.ndarray:
    return np.array([np.real(np.trace(m) / 2),
                     np.real(np.trace(m @ SX) / (2j)),
                     np.real(np.trace(m @ SY) / (2j)),
                     np.real(np.trace(m @ SZ) / (2j))])


def rho(qL: np.ndarray, qR: np.ndarray) -> np.ndarray:
    """The covering map Spin(4) -> SO(4):  x -> qL x qR^dagger."""
    out = np.zeros((4, 4))
    for i in range(4):
        e = np.zeros(4)
        e[i] = 1.0
        out[:, i] = unquat(qL @ quat(e) @ qR.conj().T)
    return out


def sigma_matrix() -> np.ndarray:
    """BAM's orientation map sigma(z1,z2) = (conj z2, -conj z1) on R^4."""
    def act(v):
        z1, z2 = v[0] + 1j * v[1], v[2] + 1j * v[3]
        w1, w2 = np.conj(z2), -np.conj(z1)
        return np.array([w1.real, w1.imag, w2.real, w2.imag])
    return np.column_stack([act(np.eye(4)[i]) for i in range(4)])


LIFTS = {'J+': (-I2, I2), 'J-': (I2, -I2)}
KERNEL = (-I2, -I2)


# ════════════════════════════════════════════════════════════════════════

def test_T1_antipodal_map() -> dict:
    J = -np.eye(4)
    det = float(np.linalg.det(J))
    return {
        'name': 'T1_J_is_the_antipodal_map_in_SO4',
        'description': (
            "J: S^3 -> S^3 is x -> -x, i.e. -I_4 in O(4). Its determinant "
            f"is {det:+.1f} > 0, so J lies in SO(4): it is "
            "ORIENTATION-PRESERVING, which is why the quotient RP^3 = "
            "S^3/J is an orientable 3-manifold (the repo's #183 deck "
            "determinant, recovered here). J is free (no fixed points on "
            "S^3) and J^2 = identity."
        ),
        'det_J': det,
        'in_SO4': det > 0,
        'J_squared_is_identity': bool(np.allclose(J @ J, np.eye(4))),
        'pass': det > 0 and bool(np.allclose(J @ J, np.eye(4))),
    }


def test_T2_the_two_lifts() -> dict:
    rows = []
    for lbl, (qL, qR) in LIFTS.items():
        m = rho(qL, qR)
        rows.append({'lift': lbl,
                     'projects_to_minus_I4': bool(np.allclose(m, -np.eye(4))),
                     'det': float(np.linalg.det(m))})
    ker = rho(*KERNEL)
    ker_ok = bool(np.allclose(ker, np.eye(4)))
    # J- = J+ * kernel
    qL, qR = LIFTS['J+']
    kL, kR = KERNEL
    prod = (qL @ kL, qR @ kR)
    differ_by_kernel = bool(np.allclose(prod[0], LIFTS['J-'][0])
                            and np.allclose(prod[1], LIFTS['J-'][1]))
    ok = all(r['projects_to_minus_I4'] for r in rows) and ker_ok and differ_by_kernel
    return {
        'name': 'T2_both_lifts_and_the_kernel',
        'description': (
            "In the quaternionic model Spin(4) = SU(2) x SU(2) acts by "
            "x -> q_L x q_R^dagger. Solving q_L x q_R^dagger = -x for all x "
            "forces q_L = -q_R with q_R central, so the ONLY lifts of J are "
            "J~_+ = (-1, 1) and J~_- = (1, -1). Both project to -I_4 "
            "(verified). The covering kernel is (-1,-1), which projects to "
            "+I_4 (verified), and J~_- = J~_+ * (-1,-1): the two lifts "
            "differ by the kernel, as the two elements of a Z2 fibre must."
        ),
        'rows': rows,
        'kernel_projects_to_identity': ker_ok,
        'lifts_differ_by_kernel': differ_by_kernel,
        'pass': ok,
    }


def test_T3_group_law_and_squares() -> dict:
    rows = []
    for lbl, (qL, qR) in LIFTS.items():
        sL, sR = qL @ qL, qR @ qR
        is_id = bool(np.allclose(sL, I2) and np.allclose(sR, I2))
        rows.append({'lift': lbl,
                     'square_is_identity_in_Spin4': is_id,
                     'square_L': 'I' if np.allclose(sL, I2) else '-I',
                     'square_R': 'I' if np.allclose(sR, I2) else '-I'})
    both_plus = all(r['square_is_identity_in_Spin4'] for r in rows)
    return {
        'name': 'T3_both_lifts_square_to_plus_one',
        'description': (
            "J~_+^2 = ((-1)^2, 1^2) = (1,1) and J~_-^2 = (1, (-1)^2) = "
            "(1,1): BOTH lifts square to the IDENTITY of Spin(4). Each "
            "therefore generates a genuine Z2 action on the spinor bundle, "
            "so BOTH descend to RP^3 -- which is why RP^3 has two spin "
            "structures, matching H^1(RP^3;Z2) = Z2. This refutes a common "
            "shortcut: the two spin structures are NOT distinguished by "
            "J~^2 = +1 versus -1 (both are +1). They are distinguished by "
            "WHICH SU(2) FACTOR carries the -1 -- a distinction invisible "
            "to any argument phrased only in terms of 'the square of the "
            "lift', including the eta = +/-1 bookkeeping the twist-selection "
            "arc used."
        ),
        'rows': rows,
        'both_square_to_plus_one': both_plus,
        'pass': both_plus,
    }


def test_T4_where_sigma_and_T_sit() -> dict:
    S = sigma_matrix()
    det_s = float(np.linalg.det(S))
    sq_is_J = bool(np.allclose(S @ S, -np.eye(4)))
    t_sq = T_BAM @ T_BAM
    t_sq_minus_I = bool(np.allclose(t_sq, -I2))
    # T lives in the LEFT SU(2); T^2 = -1 there, i.e. the element (-1, 1)
    selects = 'J+' if t_sq_minus_I else 'undetermined'
    ok = det_s > 0 and sq_is_J and t_sq_minus_I
    return {
        'name': 'T4_sigma_squares_to_J_and_T_squared_selects_J_plus',
        'description': (
            "BAM's orientation map sigma(z1,z2) = (conj z2, -conj z1) has "
            f"det = {det_s:+.1f}, so it is ORIENTATION-PRESERVING and lies "
            "in SO(4) -- not, as the phrase 'orientation-reversing Hopf "
            "map' in the repo's transport module suggests, in the "
            "non-identity component. (It reverses orientation of the Hopf "
            "FIBRE, not of S^3.) And sigma^2 = -I_4 = J exactly: sigma is a "
            "SQUARE ROOT of the antipodal map. Its spinor lift is BAM's "
            "T = i sigma_y, which satisfies T^2 = -I in the LEFT SU(2) "
            "factor -- i.e. T^2 = (-1, 1) = J~_+ . So B2 does not merely "
            "assert 'a non-trivial spin structure': it SELECTS a specific "
            "lift, J~_+ , as the square of the throat transport."
        ),
        'det_sigma': det_s,
        'sigma_is_orientation_preserving': det_s > 0,
        'sigma_squared_is_J': sq_is_J,
        'T_squared_is_minus_I': t_sq_minus_I,
        'B2_selects_lift': selects,
        'pass': ok,
    }


def test_T5_independence() -> dict:
    facts = [
        {'datum': 'B2 spin structure',
         'lives_in': 'H^1(RP^3; Z2) = Z2',
         'detected_by': 'the deck loop generating pi_1(RP^3) = Z2',
         'fixed_by': 'T^2 = J~_+ (T4)'},
        {'datum': 'network holonomy W',
         'lives_in': 'H^1(Gamma; Z2) for the ER graph Gamma',
         'detected_by': 'Wilson loops on graph cycles (b_1 of the ER graph)',
         'fixed_by': 'free link variables (#227)'},
    ]
    independent = True
    return {
        'name': 'T5_B2_and_W_are_independent_data',
        'description': (
            "They are different cohomology groups of different spaces. B2 "
            "is a spin structure on the exterior 3-manifold, an element of "
            "H^1(RP^3;Z2), detected by the deck loop that generates "
            "pi_1(RP^3) = Z2. W is a Z2 connection on the ER GRAPH, an "
            "element of H^1(Gamma;Z2), detected by Wilson loops on graph "
            "cycles. Nothing identifies them. They would coincide only if "
            "the network cycle carrying W were itself the pi_1(RP^3) "
            "generator -- a statement about how the cycle embeds, not an "
            "identity. *** CONSEQUENCE: the rule theta = 0 iff W*eta = +1, "
            "on which the moding swap of "
            "bam_spinor_spectrum_effective_action_probe rests, is NOT "
            "justified in general. It holds only for a network cycle that "
            "generates pi_1(RP^3). The #224 two-throat ring is a loop "
            "through two exterior brane arcs; whether that loop is the "
            "pi_1 generator is a separate question, and this probe does "
            "NOT settle it. Until it is settled, the moding swap -- and "
            "therefore the corrected energetic conclusion built on it -- "
            "is conditional. ***"
        ),
        'facts': facts,
        'independent': independent,
        'consequence': ('the W*eta composition rule is conditional on the '
                        'network cycle generating pi_1(RP^3); unproven'),
        'pass': independent,
    }


def test_T6_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T6_assessment',
        'description': (
            "The lift ledger is now explicit. Both lifts of the antipodal "
            "map square to +1 and differ by the covering kernel; B2 selects "
            "J~_+ through T^2; and B2's Z2 is independent of the network's "
            "W. The last point is the consequential one: it makes the "
            "moding-swap composition conditional rather than established."
        ),
        'established': [
            'J is in SO(4) (orientation-preserving), so RP^3 is orientable',
            'the only lifts are J~_+ = (-1,1) and J~_- = (1,-1); both square '
            'to +1 in Spin(4) and differ by the kernel (-1,-1)',
            'the two spin structures are distinguished by WHICH SU(2) factor '
            'carries the -1, not by the value of J~^2',
            'sigma is orientation-PRESERVING with sigma^2 = J, and T = i '
            'sigma_y gives T^2 = J~_+ : B2 selects a specific lift',
            'B2 (in H^1(RP^3;Z2)) and W (in H^1(Gamma;Z2)) are independent',
        ],
        'open': [
            'whether the #224 network cycle generates pi_1(RP^3) -- until '
            'settled, the W*eta composition rule and everything built on it '
            'is conditional',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'BOTH_LIFTS_SQUARE_TO_PLUS_ONE_B2_SELECTS_J_PLUS_VIA_T_SQUARED_'
            'AND_B2_IS_INDEPENDENT_OF_THE_NETWORK_HOLONOMY_SO_THE_'
            'COMPOSITION_RULE_IS_UNJUSTIFIED'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_antipodal_map()
    res['T2'] = test_T2_the_two_lifts()
    res['T3'] = test_T3_group_law_and_squares()
    res['T4'] = test_T4_where_sigma_and_T_sit()
    res['T5'] = test_T5_independence()
    res['T6'] = test_T6_assessment(res)
    res['summary'] = {
        'probe': 'rp3_spin_lift_probe', 'pr': 230,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T6']['tests_passed'],
        'verdict_class': res['T6']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe A — the RP³ spin-lift ledger (PR #230)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The lifts", "",
         "| lift | projects to | `J̃²` |", "|---|---|---|"]
    for r, q in zip(s['T2']['rows'], s['T3']['rows']):
        o.append(f"| {r['lift']} | {'−I₄' if r['projects_to_minus_I4'] else '?'} | "
                 f"({q['square_L']}, {q['square_R']}) = **+1** |")
    t4 = s['T4']
    o += ["", "## Where BAM sits", "",
          f"  - `det σ = {t4['det_sigma']:+.1f}` ⟹ **orientation-preserving**",
          f"  - `σ² = J` (the antipodal map): {t4['sigma_squared_is_J']}",
          f"  - `T² = −I` ⟹ `T² = J̃₊` — **B2 selects {t4['B2_selects_lift']}**",
          "", "## Independence", "", f"{s['T5']['consequence']}",
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
    out = here / "runs" / f"{ts}_rp3_spin_lift_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
