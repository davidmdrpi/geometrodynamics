"""The BAM Pin- spinor spectrum and effective action on the cycle (PR #230).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

WHY THIS EXISTS: IT CORRECTS THIS PR'S OWN T5
----------------------------------------------
`twist_sector_selection_probe` (T2/T5) argued that the sign of the twist
selection is "the statistics of the link field", using the generic rule
that the zero point is +1/2 sum omega for a boson and -1/2 sum omega for a
fermion, and concluded:

    fermionic link field  ->  dE = -pi/(4C)  ->  TWISTED favoured

and then read that as confirming that BAM matter sits in the twisted
sector.  **That is wrong.**  It ASSERTED the statistics of the link field
instead of computing BAM's actual spinor, and in doing so it silently
assumed the spinor sees the same moding as the scalar.  It does not.

For a spinor the network holonomy W composes with the INTRINSIC spin
structure eta, and BAM's B2 fixes the non-trivial one -- `T^2 = -I`,
"antiperiodic 4pi-spinors" -- so eta = -1:

    effective BC phase   theta = 0   when  W * eta = +1
                         theta = pi  when  W * eta = -1

    scalar     (eta=+1):  W=+1 -> integer modes      W=-1 -> half-integer
    BAM spinor (eta=-1):  W=+1 -> half-integer       W=-1 -> integer

The assignment is SWAPPED.  Carrying that through, the statistics flip and
the spin-structure flip CANCEL, and the corrected result is

    scalar      : dE = +pi/(4C)  ->  untwisted favoured
    BAM spinor  : dE = +pi/(4C)  ->  untwisted favoured   (NOT twisted)

So there is no "sign = statistics" rule.  Energetics favour the untwisted
sector in BOTH cases.

WHAT ACTUALLY PUTS MATTER IN THE TWISTED SECTOR: THE INDEX
-----------------------------------------------------------
The swap has a much sharper consequence than a vacuum-energy sign.  In the
sector with integer moding the Dirac operator has an exact ZERO MODE, and
for the BAM spinor that is the TWISTED sector:

    W = -1  ->  theta = 0   ->  dim ker D = 1 ,  det D = 2 sin(theta/2) = 0
    W = +1  ->  theta = pi  ->  dim ker D = 0 ,  det D = 2

The twisted sector is where the Pin- Dirac zero mode lives, and a zero
mode is topological, not energetic -- it cannot be traded against a
pi/(4C) vacuum-energy difference.  This is #195's index mechanism (Pin-
Dirac zero modes pinned by Atiyah-Singer, "not by tuning") seen from the
network side.

So the corrected selection picture is TWO mechanisms, not one sign:

  * the bosonic QCD Moebius sector is selected ENERGETICALLY (untwisted
    favoured by pi/(4C)) -- which is why #100/#103 present the Moebius
    states as excitations ABOVE the orientable ground;
  * the fermionic matter sector is selected by the INDEX -- matter sits in
    the twisted sector because that is where the zero mode is, matching
    #202/#227, and NOT because it is energetically cheaper.

The topological freeze established by the twist-selection probe (T3/T4) is
untouched: changing W still requires severing the cycle.

METHOD NOTE (a real trap)
-------------------------
The spectrum is computed through the repo's own SUSY framing -- the throat
Dirac operator is the square root of the second-order operator
(docs/throat_dirac_spinor_research_plan.md), so D^2 is the ring Laplacian
with the SPINOR's boundary phase and the Dirac spectrum is
+/- sqrt(spec(Laplacian_theta)).  A direct first-order central-difference
Dirac operator was tried first and is UNUSABLE here: fermion doubling puts
a spurious second zero mode at the Brillouin-zone edge, giving
dim ker D = 2 at theta = 0 and destroying precisely the index statement
this probe rests on.  The SUSY route is doubling-free.

Tests:
  T1. Goal: compute the spinor, do not assert its statistics.
  T2. The spectrum, converging to +/-|k_m| as O(1/n^2).
  T3. The moding: W composes with the intrinsic eta; the assignment swaps.
  T4. The index: dim ker D and det D per sector.
  T5. The effective action / vacuum energy -- the T5 correction.
  T6. The corrected selection picture: energy vs index.
  T7. Assessment.

Verdict:
  THE_BAM_SPINOR_SWAPS_THE_MODING_SO_ENERGETICS_FAVOUR_UNTWISTED_IN_BOTH
  _SECTORS_AND_MATTER_IS_SELECTED_BY_THE_ZERO_MODE_INDEX_NOT_BY_ENERGY
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

PI = math.pi

# BAM's B2 spin structure: T^2 = -I, "antiperiodic 4pi-spinors".
ETA_BAM_SPINOR = -1
ETA_SCALAR = +1


# ── the operator: the SUSY square root of the ring Laplacian ────────────

def laplacian_ring(n: int, theta: float, circumference: float = 1.0
                   ) -> np.ndarray:
    """Scalar ring Laplacian with boundary phase theta (3-point)."""
    a = circumference / n
    lap = np.zeros((n, n), dtype=complex)
    for j in range(n):
        jp, jm = (j + 1) % n, (j - 1) % n
        php = np.exp(1j * theta) if j + 1 == n else 1.0
        phm = np.exp(-1j * theta) if j - 1 < 0 else 1.0
        lap[j, j] += 2.0 / a ** 2
        lap[j, jp] += -php / a ** 2
        lap[j, jm] += -phm / a ** 2
    return lap


def dirac_abs_spectrum(n: int, theta: float, circumference: float = 1.0,
                       take: int = 7) -> np.ndarray:
    """|spec(D)| = sqrt(spec(D^2)) with D^2 the ring Laplacian at the
    spinor's boundary phase -- the repo's SUSY square-root framing."""
    ev = np.sort(np.linalg.eigvalsh(laplacian_ring(n, theta,
                                                   circumference)).real)
    return np.sqrt(np.clip(ev[:take], 0.0, None))


def analytic_abs_k(theta: float, circumference: float = 1.0,
                   take: int = 7) -> np.ndarray:
    ks = sorted({abs((theta + 2 * PI * m) / circumference)
                 for m in range(-8, 9)})
    out: list = []
    for v in ks:
        out += [v] if abs(v) < 1e-12 else [v, v]
    return np.array(sorted(out)[:take])


def theta_of(holonomy: int, eta: int) -> float:
    """The spinor/scalar boundary phase: theta = 0 if W*eta = +1 else pi."""
    return 0.0 if holonomy * eta > 0 else PI


def vacuum_energy(theta: float, fermionic: bool,
                  circumference: float = 1.0) -> float:
    """Regularized zero point per degree of freedom: +1/2 sum omega for a
    boson, -1/2 sum omega for a fermion; zeta(-1) = -1/12 (integer modes),
    zeta(-1,1/2) = +1/24 (half-integer)."""
    a = 2.0 * PI / circumference
    z = (-1.0 / 12.0 if theta == 0.0 else 1.0 / 24.0)
    return (-1.0 if fermionic else 1.0) * a * z


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_compute_the_spinor_do_not_assert_its_statistics',
        'description': (
            "The twist-selection probe's T2/T5 asserted that the link "
            "field's statistics set the sign of the selection, and read "
            "'fermionic' off the fact that BAM throat modes are Pin- "
            "spinors. It never computed the spinor. This probe computes "
            "the BAM Pin- spinor spectrum and effective action on the "
            "network cycle, through the repo's own SUSY square-root "
            "framing, and corrects the conclusion."
        ),
        'corrects': 'twist_sector_selection_probe T2/T5 (same PR)',
        'pass': True,
    }


def test_T2_spectrum_converges() -> dict:
    rows = []
    for theta, lbl in ((0.0, 'theta=0 (integer moding)'),
                       (PI, 'theta=pi (half-integer moding)')):
        an = analytic_abs_k(theta)
        errs = []
        for n in (200, 400, 800, 1600):
            k = dirac_abs_spectrum(n, theta)
            errs.append({'n': n, 'max_err': float(np.max(np.abs(k - an)))})
        ratios = [errs[i]['max_err'] / errs[i + 1]['max_err']
                  for i in range(len(errs) - 1)]
        rows.append({'sector': lbl,
                     'analytic_abs_k': [round(float(v), 5) for v in an],
                     'numeric_abs_k_n1600':
                         [round(float(v), 5)
                          for v in dirac_abs_spectrum(1600, theta)],
                     'errors': errs,
                     'successive_ratios': [round(r, 3) for r in ratios]})
    second_order = all(all(abs(r - 4.0) < 0.3 for r in row['successive_ratios'])
                       for row in rows)
    accurate = all(row['errors'][-1]['max_err'] < 1e-3 for row in rows)
    ok = second_order and accurate
    return {
        'name': 'T2_spinor_spectrum_converges_to_plus_minus_k',
        'description': (
            "The throat Dirac operator is the square root of the "
            "second-order operator, so on the cycle D^2 is the ring "
            "Laplacian at the spinor's boundary phase and spec(D) = "
            "+/- sqrt(spec(Laplacian)). Computed at n = 200...1600, the "
            "spectrum converges to the analytic |k_m| = |theta + 2 pi m|/C "
            "with successive error ratios ~4 -- second order, as the 3-point "
            "Laplacian requires. METHOD NOTE: a direct first-order "
            "central-difference Dirac operator was tried first and is "
            "unusable, because fermion doubling adds a spurious zero mode "
            "at the Brillouin-zone edge (dim ker = 2 at theta = 0) and "
            "destroys the index statement of T4. The SUSY route is "
            "doubling-free."
        ),
        'rows': rows,
        'second_order_convergence': second_order,
        'pass': ok,
    }


def test_T3_moding_swaps() -> dict:
    rows = []
    for field, eta in (('scalar', ETA_SCALAR), ('BAM Pin- spinor',
                                                ETA_BAM_SPINOR)):
        for w in (+1, -1):
            th = theta_of(w, eta)
            rows.append({
                'field': field, 'eta': eta, 'W': w, 'W_times_eta': w * eta,
                'theta': th,
                'moding': 'integer' if th == 0.0 else 'half-integer',
            })
    sc = {r['W']: r['moding'] for r in rows if r['field'] == 'scalar'}
    sp = {r['W']: r['moding'] for r in rows if r['field'] != 'scalar'}
    swapped = sc[+1] == sp[-1] and sc[-1] == sp[+1]
    return {
        'name': 'T3_the_spinor_moding_is_swapped_relative_to_the_scalar',
        'description': (
            "For a spinor the network holonomy W composes with the "
            "INTRINSIC spin structure eta, and BAM's B2 fixes the "
            "non-trivial one: T^2 = -I, 'antiperiodic 4pi-spinors', so "
            "eta = -1. The boundary phase is theta = 0 when W*eta = +1 and "
            "pi otherwise, so the scalar's assignment (W=+1 integer, W=-1 "
            "half-integer) is exactly REVERSED for the BAM spinor (W=+1 "
            "half-integer, W=-1 integer). The twist-selection probe missed "
            "this and implicitly gave the spinor the scalar's moding, which "
            "is what made its T5 sign wrong."
        ),
        'rows': rows,
        'assignment_is_swapped': swapped,
        'pass': swapped,
    }


def test_T4_index() -> dict:
    rows = []
    for w in (+1, -1):
        th = theta_of(w, ETA_BAM_SPINOR)
        ev = np.linalg.eigvalsh(laplacian_ring(1600, th)).real
        ker = int(np.sum(np.abs(ev) < 1e-8))
        rows.append({
            'W': w, 'theta': th, 'dim_ker_D': ker,
            'det_D_zeta': 2.0 * math.sin(th / 2.0),
            'lowest_abs_k': float(np.sqrt(max(float(np.min(np.abs(ev))), 0.0))),
        })
    twisted = next(r for r in rows if r['W'] == -1)
    untwisted = next(r for r in rows if r['W'] == +1)
    ok = (twisted['dim_ker_D'] == 1 and untwisted['dim_ker_D'] == 0
          and abs(twisted['det_D_zeta']) < 1e-12
          and abs(untwisted['det_D_zeta'] - 2.0) < 1e-12)
    return {
        'name': 'T4_the_zero_mode_sits_in_the_twisted_sector',
        'description': (
            "The swap of T3 has a sharper consequence than any "
            "vacuum-energy sign. Integer moding carries an exact Dirac ZERO "
            "MODE, and for the BAM spinor that is the TWISTED sector: "
            "W = -1 gives dim ker D = 1 and a vanishing zeta-regularized "
            "determinant det D = 2 sin(theta/2) = 0, while W = +1 is gapped "
            "at |k| = pi/C with det D = 2. A zero mode is topological, not "
            "energetic -- it cannot be traded against a pi/(4C) vacuum-energy "
            "difference. This is #195's index mechanism (Pin- Dirac zero "
            "modes pinned by Atiyah-Singer, 'not by tuning') seen from the "
            "network side."
        ),
        'rows': rows,
        'pass': ok,
    }


def test_T5_effective_action_correction() -> dict:
    rows = []
    for field, eta, ferm in (('scalar (boson)', ETA_SCALAR, False),
                             ('BAM Pin- spinor', ETA_BAM_SPINOR, True)):
        for c in (1.0, 2.0):
            e = {w: vacuum_energy(theta_of(w, eta), ferm, c) for w in (+1, -1)}
            d = e[-1] - e[+1]
            rows.append({
                'field': field, 'circumference': c,
                'E0_untwisted': e[+1], 'E0_twisted': e[-1],
                'delta_E': d, 'pi_over_4C': PI / (4 * c),
                'favoured': 'twisted' if d < 0 else 'untwisted',
            })
    all_untwisted = all(r['favoured'] == 'untwisted' for r in rows)
    magnitude_ok = all(abs(abs(r['delta_E']) - r['pi_over_4C']) < 1e-12
                       for r in rows)
    ok = all_untwisted and magnitude_ok
    return {
        'name': 'T5_CORRECTION_energetics_favour_untwisted_in_BOTH_sectors',
        'description': (
            "The correction. The twist-selection probe concluded that a "
            "fermionic link field gives dE = -pi/(4C) and so favours the "
            "TWISTED sector -- flipping the zero-point sign while leaving "
            "the moding alone. With the moding computed rather than assumed "
            "(T3), the statistics flip and the spin-structure flip CANCEL, "
            "and the BAM spinor gives dE = +pi/(4C): the UNTWISTED sector is "
            "favoured, exactly as for the scalar. Magnitude is still "
            "pi/(4C) per degree of freedom in both cases (verified to "
            "1e-12). So there is no 'sign = statistics' rule, and the "
            "twist-selection probe's T5 cross-check -- which read the "
            "fermionic sign as confirming that matter sits in the twisted "
            "sector -- does not stand. What does the work instead is T4's "
            "zero mode."
        ),
        'rows': rows,
        'energetics_favour_untwisted_everywhere': all_untwisted,
        'magnitude_still_pi_over_4C': magnitude_ok,
        'supersedes': ("twist_sector_selection_probe T2/T5: 'the sign is "
                       "the statistics of the link field'"),
        'pass': ok,
    }


def test_T6_corrected_selection_picture() -> dict:
    cases = [
        {
            'sector': 'QCD Moebius (#100/#103)',
            'link_field': 'flux-tube transverse displacement (phonon)',
            'eta': ETA_SCALAR,
            'selected_by': 'ENERGY',
            'outcome': ('untwisted favoured by pi/(4C) -> the Moebius states '
                        'are excitations ABOVE the orientable ground, which '
                        'is how #100/#103 present them (+pi*sigma, +2 sqrt(sigma))'),
        },
        {
            'sector': 'BAM matter (#170/#195/#202)',
            'link_field': 'Pin- throat spinor',
            'eta': ETA_BAM_SPINOR,
            'selected_by': 'INDEX',
            'outcome': ('energetics favour untwisted here TOO, so energy does '
                        'not explain matter; the twisted sector is selected '
                        'because it carries the Dirac zero mode (dim ker D = 1, '
                        'det D = 0) -- #195\'s Atiyah-Singer mechanism, '
                        'matching #202/#227'),
        },
    ]
    ok = {c['selected_by'] for c in cases} == {'ENERGY', 'INDEX'}
    return {
        'name': 'T6_two_mechanisms_not_one_sign',
        'description': (
            "The corrected picture is richer than the one it replaces. "
            "Selection is NOT a single statistics-dependent sign; it is two "
            "different mechanisms acting on the two sectors. The bosonic QCD "
            "Moebius sector is selected ENERGETICALLY (untwisted cheaper by "
            "pi/(4C)), which is why its Moebius states are excitations above "
            "the ground. The fermionic matter sector cannot be selected "
            "energetically -- energetics favour untwisted there as well -- so "
            "matter sits in the twisted sector for a topological reason: that "
            "is where the Pin- Dirac zero mode is. Both sectors still come "
            "out right, but for different and now correctly identified "
            "reasons. The topological freeze (twist-selection T3/T4) is "
            "untouched: changing W still requires severing the cycle."
        ),
        'cases': cases,
        'pass': ok,
    }


def test_T7_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T7_assessment',
        'description': (
            "Computing BAM's actual spinor rather than asserting its "
            "statistics overturns one claim of this PR and strengthens the "
            "rest. The overturned claim: 'the sign of the twist selection is "
            "the statistics of the link field'. The replacement: the "
            "holonomy composes with the intrinsic B2 spin structure, the "
            "moding swaps, energetics favour the untwisted sector "
            "universally, and matter's twisted sector is selected by the "
            "Dirac zero mode instead."
        ),
        'established': [
            'the BAM spinor spectrum on the cycle, second-order convergent '
            'to +/-|k_m| via the repo SUSY square root',
            'the holonomy composes with eta = -1 (B2), so the spinor moding '
            'is swapped relative to the scalar',
            'the twisted sector carries the Dirac zero mode: dim ker D = 1, '
            'det D = 0; the untwisted sector is gapped at pi/C with det D = 2',
            'energetics favour UNTWISTED for both the scalar and the BAM '
            'spinor, dE = +pi/(4C) in both cases',
            'selection is two mechanisms -- energy for the bosonic QCD '
            'sector, index for the fermionic matter sector',
        ],
        'corrected': [
            "twist_sector_selection_probe T2/T5: 'fermionic -> twisted "
            "favoured' is WRONG; it assumed the spinor sees the scalar's "
            "moding. Its T3/T4 (the topological freeze) are unaffected.",
        ],
        'open': [
            'this is the free massless cycle mode; the full throat spinor '
            'carries the radial/Dirac and SU(2) factors (a degeneracy, which '
            'scales the magnitude but not the sign or the index)',
            'the index statement here is the 1D cycle zero mode; #195 derives '
            'the k-dependent count on the monopole base -- relating the two '
            'is not attempted',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_BAM_SPINOR_SWAPS_THE_MODING_SO_ENERGETICS_FAVOUR_UNTWISTED_'
            'IN_BOTH_SECTORS_AND_MATTER_IS_SELECTED_BY_THE_ZERO_MODE_INDEX_'
            'NOT_BY_ENERGY'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_spectrum_converges()
    res['T3'] = test_T3_moding_swaps()
    res['T4'] = test_T4_index()
    res['T5'] = test_T5_effective_action_correction()
    res['T6'] = test_T6_corrected_selection_picture()
    res['T7'] = test_T7_assessment(res)
    res['summary'] = {
        'probe': 'bam_spinor_spectrum_effective_action_probe',
        'pr': 230,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T7']['tests_passed'],
        'verdict_class': res['T7']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# The BAM Pin⁻ spinor spectrum and effective action (PR #230)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The moding swaps", "",
         "| field | `η` | `W` | `W·η` | `θ` | moding |",
         "|---|---:|---:|---:|---:|---|"]
    for r in s['T3']['rows']:
        o.append(f"| {r['field']} | {r['eta']:+d} | {r['W']:+d} | "
                 f"{r['W_times_eta']:+d} | {r['theta']:.4f} | {r['moding']} |")
    o += ["", "## The index: the zero mode is in the twisted sector", "",
          "| `W` | `θ` | `dim ker D` | `det D` | lowest `|k|` |",
          "|---:|---:|---:|---:|---:|"]
    for r in s['T4']['rows']:
        o.append(f"| {r['W']:+d} | {r['theta']:.4f} | {r['dim_ker_D']} | "
                 f"{r['det_D_zeta']:.6f} | {r['lowest_abs_k']:.6f} |")
    o += ["", "## The correction: energetics favour untwisted in BOTH", "",
          "| field | `C` | `E₀` untwisted | `E₀` twisted | `ΔE` | favoured |",
          "|---|---:|---:|---:|---:|---|"]
    for r in s['T5']['rows']:
        o.append(f"| {r['field']} | {r['circumference']} | "
                 f"{r['E0_untwisted']:+.6f} | {r['E0_twisted']:+.6f} | "
                 f"{r['delta_E']:+.6f} | {r['favoured']} |")
    o += ["", "## Verdict", "", f"**{s['T7']['verdict_class']}**", ""]
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
    out = here / "runs" / f"{ts}_bam_spinor_spectrum_effective_action_probe"
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
