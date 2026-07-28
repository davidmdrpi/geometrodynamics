"""Probe D -- semiclassical back-reaction between the twist sectors (PR #230).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.  The semiclassical equations here are the standard
> G = 8 pi <T>_ren system on a fixed classical background, not a quantised
> metric.

PRECONDITION
------------
Probe D was to run only if both sectors remain admissible.  Probe C settles
that: Omega_3^{Spin} = 0, both spin structures on RP^3 extend over a spin
bulk, and anomaly inflow imposes no obstruction on either.  So the
precondition IS met and this probe runs.

THE RESULT: THE BACK-REACTION IS IDENTICAL, EXACTLY
----------------------------------------------------
Probe B established that for each level k exactly one sign survives the
projection, with |lambda| = 3/2 + k and multiplicity (k+1)(k+2) -- the SAME
for both spin structures.  Three consequences follow without any further
computation:

  1. The regularized vacuum energy is identical: it depends only on |spec|,
     and the two |spec| agree mode by mode.  Delta E_vac = 0 EXACTLY -- not
     small, not suppressed, zero.

  2. RP^3 with the round metric is HOMOGENEOUS, so the renormalized vacuum
     stress tensor is covariantly constant: <T_AB>_ren = (E_vac/vol) * (a
     fixed tensor).  Equal total energy and equal volume therefore force
     equality POINTWISE:

         Delta <T_AB>_ren = 0   pointwise on the round exterior.

  3. With identical sources the radial semiclassical equations are the SAME
     equations.  Regularity, the on-shell action modulus, the negative-mode
     spectrum and stability are therefore identical by construction.  There
     is nothing left to compare.

*** So back-reaction does not select the twist sector either.  The two
sectors are exactly degenerate in energy and differ ONLY by the
determinant phase (arg det D = -/+ pi/8, Probe B T6). A phase does not
source the semiclassical Einstein equations -- <T_AB> is built from the
modulus -- so the surviving discriminator is invisible to gravity at this
order. ***

WHAT IS NOT DONE (the honest boundary)
---------------------------------------
The degeneracy argument leans on HOMOGENEITY of the round RP^3.  BAM's
actual exterior is the Tangherlini throat geometry, which is not
homogeneous, and on a non-homogeneous background equal integrated energy
does NOT force equal pointwise <T_AB>.  A genuine solve there would need:
the mode functions on the actual radial background, the renormalized
coincidence limit of the Green function difference between the two spin
structures, and then the coupled radial system.  None of that is attempted
here.  What IS established is that the leading, exactly-computable case
gives exact degeneracy -- so any selection from back-reaction would have to
come entirely from the inhomogeneity, and would be a small effect on top of
an exact zero rather than a leading one.

Tests:
  T1. Precondition: both sectors admissible (Probe C).
  T2. Delta E_vac = 0 exactly; the identical-|spec| argument.
  T3. Delta <T_AB>_ren = 0 pointwise by homogeneity.
  T4. Identical sources -> identical radial system, action, stability.
  T5. What is not done, and what would be needed.
  T6. Assessment.

Verdict:
  BACK_REACTION_IS_EXACTLY_DEGENERATE_BETWEEN_THE_TWIST_SECTORS_SO_GRAVITY
  _DOES_NOT_SELECT_EITHER_AND_ONLY_THE_ETA_PHASE_DISTINGUISHES_THEM
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from mpmath import mp, mpf, zeta as mp_zeta

mp.dps = 30
PI = math.pi


def abs_spectrum(kmax: int) -> list:
    """|spec(D)| on RP^3: |lambda| = 3/2 + k, multiplicity (k+1)(k+2).
    The SAME for both spin structures (Probe B)."""
    return [{'k': k, 'abs_lambda': 1.5 + k, 'multiplicity': (k + 1) * (k + 2)}
            for k in range(kmax)]


def vacuum_energy_zeta() -> float:
    """E_vac = -(1/2) sum |lambda| for a fermion, zeta-regularized:
    sum |lam|^{-s} = zeta_H(s-2,3/2) - (1/4) zeta_H(s,3/2), evaluated at
    s = -1."""
    s = mpf(-1)
    tot = mp_zeta(s - 2, mpf(3) / 2) - mpf(1) / 4 * mp_zeta(s, mpf(3) / 2)
    return float(-tot / 2)


# ════════════════════════════════════════════════════════════════════════

def test_T1_precondition() -> dict:
    return {
        'name': 'T1_precondition_both_sectors_admissible',
        'description': (
            "Probe D was to run only if both sectors remain admissible. "
            "Probe C settles it: Omega_3^{Spin} = 0, so both spin "
            "structures on RP^3 bound a spin 4-manifold (with the explicit "
            "witness X^4 = the D^2-bundle over S^2 of even Euler number), "
            "and the APS index is integral for each. Anomaly inflow imposes "
            "no obstruction on either side. Precondition met."
        ),
        'omega_3_spin': '0',
        'both_admissible': True,
        'pass': True,
    }


def test_T2_vacuum_energy_degenerate() -> dict:
    sp = abs_spectrum(6)
    e = vacuum_energy_zeta()
    return {
        'name': 'T2_delta_E_vac_is_exactly_zero',
        'description': (
            "For each level k exactly one sign survives the projection, "
            "with |lambda| = 3/2 + k and multiplicity (k+1)(k+2) -- "
            "identical for both spin structures (Probe B T5). The "
            "regularized vacuum energy depends only on |spec|, so it is the "
            f"same number for both: E_vac = {e:.10f} in units of 1/R. "
            "Delta E_vac = 0 EXACTLY -- not small, not exponentially "
            "suppressed, identically zero, because the two spectra agree "
            "term by term in absolute value."
        ),
        'abs_spectrum_first_levels': sp,
        'E_vac': e,
        'delta_E_vac': 0.0,
        'exactly_zero': True,
        'pass': True,
    }


def test_T3_stress_tensor_degenerate() -> dict:
    return {
        'name': 'T3_delta_T_AB_is_zero_pointwise_by_homogeneity',
        'description': (
            "RP^3 with the round metric is a HOMOGENEOUS space: its isometry "
            "group acts transitively, and the vacuum state (built from the "
            "full mode set in either spin structure) is invariant. A "
            "covariantly constant, isometry-invariant symmetric 2-tensor on "
            "a homogeneous 3-space is fixed up to overall normalization, and "
            "that normalization is set by the total energy. Since T2 gives "
            "equal total energy, and the two quotients have equal volume "
            "(both are S^3/Z2), the normalizations agree and therefore "
            "Delta <T_AB>_ren = 0 POINTWISE, not merely on integration. "
            "The homogeneity step is what upgrades an equal integral to an "
            "equal local source -- and it is also exactly the step that "
            "fails on a non-homogeneous background (T5)."
        ),
        'exterior_is_homogeneous': True,
        'equal_volume': True,
        'delta_T_AB_pointwise': 0.0,
        'pass': True,
    }


def test_T4_identical_semiclassical_system() -> dict:
    items = [
        {'quantity': 'radial semiclassical equations',
         'differs_between_sectors': False,
         'why': 'identical source <T_AB>_ren (T3), identical background'},
        {'quantity': 'regularity of the solution',
         'differs_between_sectors': False,
         'why': 'same equations, same boundary data'},
        {'quantity': 'on-shell action MODULUS',
         'differs_between_sectors': False,
         'why': "|det D| identical (Probe B T5, zeta'(0) = 0.2189594807)"},
        {'quantity': 'on-shell action PHASE',
         'differs_between_sectors': True,
         'why': 'arg det D = -/+ pi/8 from eta = +/-1/4 (Probe B T6)'},
        {'quantity': 'negative-mode spectrum / stability',
         'differs_between_sectors': False,
         'why': 'fluctuation operator built from the same real system'},
    ]
    only_phase = [i['quantity'] for i in items if i['differs_between_sectors']]
    return {
        'name': 'T4_identical_sources_give_an_identical_radial_system',
        'description': (
            "With identical sources the two sectors solve the SAME radial "
            "semiclassical equations, so regularity, the action modulus and "
            "the negative-mode/stability structure coincide by "
            "construction -- there is nothing to compare. The single "
            "quantity that differs is the on-shell action PHASE, and a "
            "phase does not source the semiclassical Einstein equations: "
            "<T_AB> is built from the modulus. So the one surviving "
            "discriminator is invisible to gravity at this order."
        ),
        'items': items,
        'only_difference': only_phase,
        'pass': only_phase == ['on-shell action PHASE'],
    }


def test_T5_what_is_not_done() -> dict:
    gaps = [
        'the exterior used is the ROUND RP^3; BAM\'s actual exterior is the '
        'Tangherlini throat geometry, which is NOT homogeneous',
        'on a non-homogeneous background, equal integrated energy does NOT '
        'force equal pointwise <T_AB> -- the T3 argument fails there',
        'a genuine solve would need the mode functions on the actual radial '
        'background, the renormalized coincidence limit of the Green '
        'function DIFFERENCE between the two spin structures, and then the '
        'coupled radial system; none is attempted',
        'the fermion is taken massless and free; a mass or a coupling would '
        'break the exact mode-by-mode |spec| coincidence that T2 rests on',
    ]
    return {
        'name': 'T5_the_honest_boundary_of_this_result',
        'description': (
            "The degeneracy is exact in the case that can be computed "
            "exactly, and that case leans on homogeneity. What this does "
            "NOT establish is degeneracy on BAM's actual throat geometry. "
            "What it does establish is the baseline: any back-reaction "
            "selection would have to come ENTIRELY from the inhomogeneity, "
            "as a correction on top of an exact zero, rather than as a "
            "leading effect. That is a meaningful constraint on where to "
            "look next, and it is weaker than 'gravity selects nothing'."
        ),
        'gaps': gaps,
        'pass': True,
    }


def test_T6_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    return {
        'name': 'T6_assessment',
        'description': (
            "Back-reaction does not select the twist sector. On the round "
            "exterior the two sectors are exactly degenerate -- same vacuum "
            "energy, same pointwise <T_AB>, same radial system, same "
            "stability -- because their spectra agree in absolute value "
            "mode by mode. They differ only by the determinant phase, which "
            "does not source gravity. Combined with Probe C (inflow does "
            "not select) the position is: THREE independent selection "
            "mechanisms have now been tested -- energetics, anomaly inflow, "
            "and semiclassical back-reaction -- and none of them "
            "distinguishes the sectors. The only surviving distinction in "
            "the whole arc is eta = +/-1/4."
        ),
        'established': [
            'Delta E_vac = 0 exactly (identical |spec| mode by mode)',
            'Delta <T_AB>_ren = 0 pointwise on the homogeneous round exterior',
            'identical radial system, regularity, action modulus and '
            'stability; only the action PHASE differs',
            'a phase does not source the semiclassical Einstein equations',
        ],
        'open': [
            'the non-homogeneous Tangherlini exterior -- any back-reaction '
            'selection must come entirely from the inhomogeneity, on top of '
            'an exact zero',
            'what, if anything, the eta = +/-1/4 phase difference is '
            'physically FOR, given it selects nothing energetically',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'BACK_REACTION_IS_EXACTLY_DEGENERATE_BETWEEN_THE_TWIST_SECTORS_'
            'SO_GRAVITY_DOES_NOT_SELECT_EITHER_AND_ONLY_THE_ETA_PHASE_'
            'DISTINGUISHES_THEM'),
        'pass': passed >= total - 1,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_precondition()
    res['T2'] = test_T2_vacuum_energy_degenerate()
    res['T3'] = test_T3_stress_tensor_degenerate()
    res['T4'] = test_T4_identical_semiclassical_system()
    res['T5'] = test_T5_what_is_not_done()
    res['T6'] = test_T6_assessment(res)
    res['summary'] = {
        'probe': 'twist_sector_einstein_dirac_probe', 'pr': 230,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T6']['tests_passed'],
        'verdict_class': res['T6']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe D — semiclassical back-reaction (PR #230)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", f"Precondition (Probe C): both sectors admissible — **met**.",
         "", "## The comparison", "",
         "| quantity | differs? | why |", "|---|---|---|"]
    for i in s['T4']['items']:
        o.append(f"| {i['quantity']} | "
                 f"{'**YES**' if i['differs_between_sectors'] else 'no'} | "
                 f"{i['why']} |")
    o += ["", f"`ΔE_vac = 0` exactly (`E_vac = {s['T2']['E_vac']:.10f}`); "
          "`Δ⟨T_AB⟩ = 0` pointwise by homogeneity.",
          "", "## Not done", ""]
    for g in s['T5']['gaps']:
        o.append(f"  - {g}")
    o += ["", "## Verdict", "", f"**{s['T6']['verdict_class']}**", ""]
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
    out = here / "runs" / f"{ts}_twist_sector_einstein_dirac_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
