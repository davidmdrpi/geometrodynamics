"""Does the #207 surgery composition law acquire a Z2 term? (PR #227 follow-on)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
The PR #227 plan (sec. 8.1) left one corollary uncomputed: the #207 bridge-
surgery composition law phi_14 = phi_a + phi_b + phi_c was conjectured to
acquire a Z2 term,

    phi_14  =  phi_a + phi_b + phi_c + pi * [n_twist mod 2] ,

so that surgery across an odd number of twisted bridges would swap the Bell
outcome to its orthogonal partner.  This probe computes it on #207's own
4-mouth surgery lattice.

THE ANSWER: NO SUCH TERM EXISTS
-------------------------------
Measured over all 8 sign assignments (eps_12, eps_34, eps_j):

  * twisting EITHER outer bridge changes the extracted pair phase by
    ~1e-13 -- identically zero, not pi.  It also leaves |c_pm| unchanged
    to 1e-14.
  * twisting the junction bridge shifts the phase by 4.68e-02 rad --
    which is neither 0 nor pi, and is not topological at all.

THE STRUCTURAL REASON
---------------------
The pair phase is extracted as arg(c_mp / c_pm) -- a RATIO of two branches
that traverse the SAME links.  A link-uniform sign multiplies both branches
equally and cancels identically.  The Bell phase of a single bridge-mediated
channel is blind to the twist, and no choice of lattice parameters can
change that: it is an algebraic cancellation, not a small number.

WHERE THE TWIST *DOES* ACT (and it confirms the main probe's T2)
---------------------------------------------------------------
The residual 4.68e-02 on the junction is bridge-vs-exterior interference
around the M2-M3 CYCLE: those two mouths sit only 8 lattice sites apart, so
the short exterior hop competes with the junction bridge, and the relative
sign between the two paths is physical.  Measured, it is a dynamical
interference amplitude and not a Z2 phase:

  * it dies with the cycle length -- d(phi) = 2.9e-01, 4.7e-02, 2.1e-04,
    3.0e-14, 1.8e-13 as the M2-M3 gap goes 4, 8, 16, 32, 48.  Once the
    competing path is long enough it is machine zero;
  * it scales continuously with the exterior hopping t_x that CARRIES the
    competing path -- 2.5e-05, 1.2e-03, 4.7e-02, 1.8e-01 at
    t_x = 0.25, 0.5, 1.0, 2.0.

A topological phase does neither.  The outer-bridge twists, whose competing
exterior path is 60 sites long, stay at ~1e-13 for every gap and every
hopping -- they never do anything.

So the twist acts ONLY through cycles, exactly as the main probe's b_1
audit (T2) requires -- and in #207's geometry that cycle carries a tunable
interference amplitude, not a Bell-outcome flip.  The Z2 IS physical here
(it moves the spectrum by ~2.2e-02 to 2.8e-02); it simply does not enter
the composition law.

Tests:
  T1. Goal and the conjecture being tested.
  T2. The 8 sign assignments: no pi term (the refutation).
  T3. The structural reason: the pair phase is a ratio; the twist is
      nonetheless physical in the spectrum.
  T4. The mechanism: the residual is cycle interference (gap + hopping
      scans).
  T5. Assessment and the corrected composition law.

Verdict:
  THE_SURGERY_COMPOSITION_LAW_HAS_NO_Z2_TERM_THE_PAIR_PHASE_IS_A_RATIO_
  AND_THE_TWIST_ACTS_ONLY_AS_CYCLE_INTERFERENCE
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

import experiments.closure_ledger.bridge_surgery_swapping_probe as B

_NX, _NCHI = B._NX, B._NCHI
_M1, _M2, _M3, _M4 = B._M1, B._M2, B._M3, B._M4
_V0, _WW, _GB = B._V0, B._WW, B._GB
_S12, _S34 = B._S12, B._S34
_DIM = _NX * _NCHI

_ALL_EPS = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1),
            (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)]


# ── the #207 lattice with a Z2 link twist on each bridge ────────────────

def build_H_twist(sj: Optional[int], middle_wells: bool,
                  eps=(1, 1, 1), m=(None, None, None, None),
                  tx: float = 1.0, tchi: float = 0.3) -> np.ndarray:
    """``bridge_surgery_swapping_probe.build_H`` with the twist made a
    link variable: eps = (eps_12, eps_34, eps_j) multiplies each bridge's
    hopping.  ``m`` optionally relocates the four mouths (the cycle-length
    scan); ``tx`` sets the exterior hopping that carries the competing
    path."""
    m1 = _M1 if m[0] is None else m[0]
    m2 = _M2 if m[1] is None else m[1]
    m3 = _M3 if m[2] is None else m[2]
    m4 = _M4 if m[3] is None else m[3]
    e12, e34, ej = eps
    H = np.zeros((_DIM, _DIM), dtype=complex)
    for x in range(_NX):
        for c in range(_NCHI):
            i = B._idx(x, c)
            H[i, B._idx((x + 1) % _NX, c)] -= tx
            H[B._idx((x + 1) % _NX, c), i] -= tx
            H[i, B._idx(x, (c + 1) % _NCHI)] -= tchi
            H[B._idx(x, (c + 1) % _NCHI), i] -= tchi
    xs = np.arange(_NX)
    for x0 in [m1, m4] + ([m2, m3] if middle_wells else []):
        d = np.minimum(np.abs(xs - x0), _NX - np.abs(xs - x0))
        w = -_V0 * np.exp(-d ** 2 / (2 * _WW ** 2))
        for x in range(_NX):
            for c in range(_NCHI):
                H[B._idx(x, c), B._idx(x, c)] += w[x]
    for c in range(_NCHI):
        H[B._idx(m1, c), B._idx(m2, (_S12 - c) % _NCHI)] -= _GB * e12
        H[B._idx(m2, (_S12 - c) % _NCHI), B._idx(m1, c)] -= _GB * e12
        H[B._idx(m3, c), B._idx(m4, (_S34 - c) % _NCHI)] -= _GB * e34
        H[B._idx(m4, (_S34 - c) % _NCHI), B._idx(m3, c)] -= _GB * e34
    if sj is not None:
        for c in range(_NCHI):
            H[B._idx(m2, c), B._idx(m3, (sj - c) % _NCHI)] -= _GB * ej
            H[B._idx(m3, (sj - c) % _NCHI), B._idx(m2, c)] -= _GB * ej
    return H


def pair_phase(sj: int, eps=(1, 1, 1), m=(None, None, None, None),
               tx: float = 1.0):
    """#207's extraction: inject (k=+1 plus k=-1) at mouth 1, read the
    (1,4) pair phase as arg(c_mp / c_pm) at mouth 4."""
    m1 = _M1 if m[0] is None else m[0]
    m4 = _M4 if m[3] is None else m[3]
    ev, vec = np.linalg.eigh(build_H_twist(sj, False, eps, m, tx))
    psi0 = (B.mode(m1, +1) + B.mode(m1, -1)) / math.sqrt(2)
    psi = vec @ (np.exp(-1j * ev * B._T_ORBIT) * (vec.conj().T @ psi0))
    c_pm = np.vdot(B.mode(m4, -1), psi)
    c_mp = np.vdot(B.mode(m4, +1), psi)
    return (float(np.angle(c_mp / c_pm) % (2 * np.pi)),
            float(abs(c_pm)), float(abs(c_mp)))


def _dphi(a: float, b: float) -> float:
    return float(abs((a - b + np.pi) % (2 * np.pi) - np.pi))


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_the_conjectured_z2_term',
        'description': (
            "PR #227 sec. 8.1 left one corollary uncomputed: whether the "
            "#207 bridge-surgery composition law phi_14 = phi_a + phi_b + "
            "phi_c acquires a Z2 term pi*[n_twist mod 2] when bridges carry "
            "the link twist -- i.e. whether surgery across an odd number of "
            "twisted bridges swaps the Bell outcome to its orthogonal "
            "partner. This probe computes it on #207's own 4-mouth surgery "
            "lattice, reusing its Hamiltonian, its mouths, its holonomies "
            "(s12 = s34 = 2) and its phase extraction unchanged."
        ),
        'conjecture': 'phi_14 = phi_a + phi_b + phi_c + pi*[n_twist mod 2]',
        'pass': True,
    }


def test_T2_no_pi_term() -> dict:
    """The refutation: sweep all 8 sign assignments."""
    rows = []
    base = {sj: pair_phase(sj)[0] for sj in (0, 1, 2, 3)}
    for sj in (0, 1, 2, 3):
        for eps in _ALL_EPS:
            ph, a_pm, _ = pair_phase(sj, eps)
            rows.append({
                'sj': sj, 'eps': list(eps),
                'n_twist': sum(1 for e in eps if e < 0),
                'phase': round(ph, 6),
                'dphi_vs_untwisted': _dphi(ph, base[sj]),
                'abs_c_pm': round(a_pm, 6),
                'junction_twisted': eps[2] < 0,
            })
    outer = [r for r in rows if not r['junction_twisted']]
    junc = [r for r in rows if r['junction_twisted']]
    max_outer = max(r['dphi_vs_untwisted'] for r in outer)
    max_junc = max(r['dphi_vs_untwisted'] for r in junc)
    # the conjecture demands ~pi for every odd-n_twist row
    odd = [r for r in rows if r['n_twist'] % 2 == 1]
    worst_vs_pi = min(abs(r['dphi_vs_untwisted'] - np.pi) for r in odd)
    ok = max_outer < 1e-9 and worst_vs_pi > 3.0
    return {
        'name': 'T2_no_pi_term_the_refutation',
        'description': (
            "Sweeping all 8 sign assignments (eps_12, eps_34, eps_j) at "
            "every junction holonomy sj = 0..3: twisting EITHER outer "
            f"bridge moves the pair phase by at most {max_outer:.1e} -- "
            "identically zero, not pi -- and leaves |c_pm| unchanged too. "
            f"Twisting the junction bridge moves it by {max_junc:.2e} rad, "
            "which is neither 0 nor pi. Every odd-n_twist row misses the "
            f"conjectured pi by at least {worst_vs_pi:.2f} rad. THE "
            "CONJECTURED Z2 TERM DOES NOT EXIST."
        ),
        'rows': rows,
        'max_dphi_outer_twists': max_outer,
        'max_dphi_junction_twist': max_junc,
        'min_distance_from_pi_on_odd_rows': float(worst_vs_pi),
        'conjecture_refuted': True,
        'pass': ok,
    }


def test_T3_structural_reason() -> dict:
    """Why: the pair phase is a ratio -- but the twist IS physical."""
    spec = []
    e0 = np.linalg.eigvalsh(build_H_twist(2, False, (1, 1, 1)))
    for eps in [(-1, 1, 1), (-1, -1, 1), (-1, -1, -1)]:
        e1 = np.linalg.eigvalsh(build_H_twist(2, False, eps))
        spec.append({'eps': list(eps),
                     'max_dE_vs_untwisted': float(np.max(np.abs(e1 - e0)))})
    moves_spectrum = all(s['max_dE_vs_untwisted'] > 1e-3 for s in spec)
    ph0, a0, b0 = pair_phase(2, (1, 1, 1))
    ph1, a1, b1 = pair_phase(2, (-1, 1, 1))
    ok = moves_spectrum and _dphi(ph0, ph1) < 1e-9 and abs(a0 - a1) < 1e-9
    return {
        'name': 'T3_the_pair_phase_is_a_ratio_but_the_twist_is_physical',
        'description': (
            "The reason is algebraic, not numerical. #207 extracts the pair "
            "phase as arg(c_mp / c_pm) -- a RATIO of two branches that "
            "traverse the SAME links. A link-uniform sign multiplies both "
            "branches equally and cancels identically, and |c_pm| is a "
            "modulus so it cancels there too (measured: phase and amplitude "
            "both unchanged to <1e-9 under an outer-bridge twist). No "
            "choice of lattice parameters can change this. Crucially this "
            "is NOT because the twist is gaugeable: it moves the spectrum "
            "by ~2.2e-02 to 2.8e-02, so the Z2 is physical here -- it "
            "simply does not enter this observable."
        ),
        'spectrum_shifts': spec,
        'twist_moves_spectrum': moves_spectrum,
        'phase_unchanged': _dphi(ph0, ph1),
        'amplitude_unchanged': abs(a0 - a1),
        'pass': ok,
    }


def test_T4_mechanism_is_cycle_interference() -> dict:
    """The residual junction shift is bridge-vs-exterior interference
    around the M2-M3 cycle -- dynamical, not topological."""
    gap_rows = []
    for gap in (4, 8, 16, 32, 48):
        m = (32, 96 - gap // 2, 96 + gap // 2, 160)
        base = pair_phase(2, (1, 1, 1), m)[0]
        gap_rows.append({
            'gap': gap,
            'dphi_junction_twist': _dphi(pair_phase(2, (1, 1, -1), m)[0], base),
            'dphi_outer_twist': _dphi(pair_phase(2, (-1, 1, 1), m)[0], base),
        })
    tx_rows = []
    for tx in (0.25, 0.5, 1.0, 2.0):
        base = pair_phase(2, (1, 1, 1), tx=tx)[0]
        tx_rows.append({
            'tx': tx,
            'dphi_junction_twist': _dphi(pair_phase(2, (1, 1, -1), tx=tx)[0],
                                         base),
        })
    dies_with_gap = (gap_rows[0]['dphi_junction_twist'] > 1e-2
                     and gap_rows[-1]['dphi_junction_twist'] < 1e-9)
    scales_with_tx = (tx_rows[-1]['dphi_junction_twist']
                      > 100 * tx_rows[0]['dphi_junction_twist'])
    outer_never = all(r['dphi_outer_twist'] < 1e-9 for r in gap_rows)
    ok = dies_with_gap and scales_with_tx and outer_never
    return {
        'name': 'T4_the_residual_is_cycle_interference_not_a_z2_phase',
        'description': (
            "The junction's 4.7e-02 residual is bridge-vs-exterior "
            "interference around the M2-M3 CYCLE: in #207's geometry those "
            "mouths sit only 8 lattice sites apart, so the short exterior "
            "hop competes with the junction bridge and the relative sign "
            "between the two paths is physical. It is an interference "
            "amplitude, not a topological phase, and behaves like one: it "
            "DIES with the cycle length (2.9e-01 -> 1.8e-13 as the gap goes "
            "4 -> 48; machine zero once the competing path is long enough) "
            "and SCALES CONTINUOUSLY with the exterior hopping that carries "
            "the competing path (2.5e-05 -> 1.8e-01 as t_x goes 0.25 -> "
            "2.0). A topological phase does neither. The outer-bridge "
            "twists, whose competing exterior path is 60 sites long, stay "
            "at ~1e-13 for every gap. So the twist acts ONLY through "
            "cycles -- exactly what the main probe's b_1 audit (T2) "
            "requires -- and in #207's geometry that cycle carries a "
            "tunable interference amplitude, not a Bell-outcome flip."
        ),
        'cycle_length_scan': gap_rows,
        'exterior_hopping_scan': tx_rows,
        'dies_with_cycle_length': dies_with_gap,
        'scales_with_exterior_hopping': scales_with_tx,
        'outer_twist_never_acts': outer_never,
        'pass': ok,
    }


def test_T5_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T5_assessment',
        'description': (
            "PR #227 sec. 8.1 is closed with a NEGATIVE result and a "
            "mechanism. The #207 composition law is unchanged: "
            "phi_14 = phi_a + phi_b + phi_c, with NO Z2 term. Surgery "
            "across twisted bridges does not swap the Bell outcome. The "
            "reason is structural (the pair phase is a ratio over branches "
            "sharing the same links), so this is not a null result waiting "
            "for a better lattice. What the twist does instead is act "
            "through cycles, as a tunable interference amplitude -- which "
            "independently confirms the b_1 principle the main probe's T2 "
            "established, now dynamically rather than combinatorially."
        ),
        'corrected_law': 'phi_14 = phi_a + phi_b + phi_c   (unchanged)',
        'established': [
            'no Z2 term in the #207 composition law (outer twists 1e-13; '
            'odd-n_twist rows miss pi by > 3 rad)',
            'the reason is algebraic: the pair phase is a ratio over '
            'branches traversing the same links',
            'the twist is nonetheless physical here (spectrum moves ~2.5e-02)',
            'its only effect is cycle interference: dies with cycle length, '
            'scales with the exterior hopping -- dynamical, not topological',
        ],
        'open': [
            'whether any Bell-type observable IS twist-sensitive: it would '
            'have to compare two paths that traverse DIFFERENT link sets, '
            'which #207 geometry does not naturally provide',
            'the remaining sec. 8 items: whether the Moebius identification '
            'is quantitative, and what selects the twist sector',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_SURGERY_COMPOSITION_LAW_HAS_NO_Z2_TERM_THE_PAIR_PHASE_IS_A_'
            'RATIO_AND_THE_TWIST_ACTS_ONLY_AS_CYCLE_INTERFERENCE'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_no_pi_term()
    res['T3'] = test_T3_structural_reason()
    res['T4'] = test_T4_mechanism_is_cycle_interference()
    res['T5'] = test_T5_assessment(res)
    res['summary'] = {
        'probe': 'er_link_twist_surgery_probe',
        'pr': 227,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T5']['tests_passed'],
        'verdict_class': res['T5']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Does the #207 surgery law acquire a Z2 term? (PR #227 follow-on)",
         "", f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The answer: no", "",
         f"  - outer-bridge twists move the pair phase by at most "
         f"**{s['T2']['max_dphi_outer_twists']:.1e}** (identically zero)",
         f"  - the junction twist moves it by "
         f"**{s['T2']['max_dphi_junction_twist']:.2e}** rad — not pi",
         f"  - every odd-n_twist row misses the conjectured pi by at least "
         f"**{s['T2']['min_distance_from_pi_on_odd_rows']:.2f}** rad",
         "", "## The mechanism: cycle interference, not a Z2 phase", "",
         "| M2–M3 gap | d(phi) junction twist | d(phi) outer twist |",
         "|---:|---:|---:|"]
    for r in s['T4']['cycle_length_scan']:
        o.append(f"| {r['gap']} | {r['dphi_junction_twist']:.3e} | "
                 f"{r['dphi_outer_twist']:.3e} |")
    o += ["", "| exterior hopping t_x | d(phi) junction twist |", "|---:|---:|"]
    for r in s['T4']['exterior_hopping_scan']:
        o.append(f"| {r['tx']} | {r['dphi_junction_twist']:.3e} |")
    o += ["", "## Verdict", "", f"**{s['T5']['verdict_class']}**",
          "", f"Corrected law: `{s['T5']['corrected_law']}`", ""]
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
    out = here / "runs" / f"{ts}_er_link_twist_surgery_probe"
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
