"""
Independent identifiability audit for the proposed O(1) lepton coefficient.

This companion probe deliberately keeps the definitions used by PR #202 and
PR #203 separate from the cavity definition introduced by PR #221:

* ``sigma_mode`` in #202/#203 is the infrared localization scale supplied by
  the soliton sector.  The vacuum throat problem has no bound state and only
  determines the near-neck suppression law.
* ``sigma_antinode`` in #221 is the distance from the cross-cap to the first
  antinode of an odd cavity mode.

The audit asks whether the existing equations derive a unique weld between
those two lengths.  A PASS means the logical status has been identified
correctly; it does not mean the mass coefficient has been derived.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from experiments.closure_ledger import (
    coupled_5d_soliton_solve_probe as p203,
    lepton_o1_coefficient_probe as p221,
    psi_phi_q_soliton_hardening_probe as p180,
)

_CACHE: dict = {}


def _weighted_radius(r: np.ndarray, density: np.ndarray, fraction: float) -> float:
    """Radius enclosing ``fraction`` of a positive radial measure."""
    dr = float(r[1] - r[0])
    weights = 4.0 * math.pi * r**2 * density * dr
    cumulative = np.cumsum(weights)
    target = fraction * float(cumulative[-1])
    return float(np.interp(target, cumulative, r))


def soliton_localization_scales() -> dict:
    """Independent IR scales from the one locked #180 soliton."""
    if "soliton" in _CACHE:
        return _CACHE["soliton"]

    sol = p180.relax(3.5, 0.05)
    r = np.asarray(sol["r"], float)
    psi = np.asarray(sol["psi"], float)
    density = psi**2
    dr = float(r[1] - r[0])
    norm = float(np.sum(4.0 * math.pi * r**2 * density) * dr)
    rms = math.sqrt(
        float(np.sum(4.0 * math.pi * r**4 * density) * dr) / norm
    )

    locked = p203.locked_scales()
    out = {
        "r50": _weighted_radius(r, density, 0.50),
        "r90": _weighted_radius(r, density, 0.90),
        "r95": _weighted_radius(r, density, 0.95),
        "rms": float(rms),
        "Rstar_A": float(locked["r_star_A"]),
        "Rstar_B": float(locked["r_star_B"]),
        "core_half": float(locked["r_q_half"]),
        "core_threshold": float(locked["r_rhoc"]),
    }
    _CACHE["soliton"] = out
    return out


def cavity_data() -> dict:
    """The robust cavity invariant established by #221."""
    if "cavity" in _CACHE:
        return _CACHE["cavity"]
    t3 = p221.test_T3_measurement()
    t4 = p221.test_T4_robustness()
    t5 = p221.test_T5_eigenhistory_orbit()
    out = {
        "omega": float(t3["odd"]["w"]),
        "sigma_antinode": float(t3["odd"]["sigma_match"]),
        "X_antinode": float(t3["odd"]["X_match"]),
        "X_amplitude_rms": float(t3["odd"]["X_u"]),
        "X_energy_rms": float(t3["odd"]["X_E"]),
        "depth_band": [float(x) for x in t4["X_match_depth_band"]],
        "q_star": float(t5["q_star"]),
        "p_star": float(t5["p_star"]),
        "half_amplitude_residual": float(
            t5["half_amplitude_periodicity_residual"]
        ),
    }
    _CACHE["cavity"] = out
    return out


def candidate_cross_sector_coefficients() -> dict:
    """Apply the cavity frequency to pre-existing soliton length definitions.

    This is not asserted to be physically correct.  It is the identifiability
    test: if the current theory supplied a unique cross-sector weld, the
    admissible definitions would collapse to one value before comparison with
    the measured lepton masses.
    """
    s = soliton_localization_scales()
    omega = cavity_data()["omega"]
    names = ("r50", "r90", "r95", "rms", "Rstar_A", "Rstar_B")
    return {name: float(omega * s[name]) for name in names}


# ========================================================================
# T1. Preserve the historical definitions
# ========================================================================


def test_T1_definition_ledger() -> dict:
    return {
        "name": "T1_definition_ledger",
        "description": (
            "The audit preserves the upstream ledger: #202 defines "
            "sigma_mode as the IR localization scale supplied by the soliton "
            "sector, because the vacuum winding problem has no bound state; "
            "#203 measures that scale from the locked #180 soliton.  The "
            "first-antinode distance introduced in #221 is recorded as a "
            "separate cavity observable, sigma_antinode."
        ),
        "sigma_mode_source": "#180/#203 localized soliton",
        "sigma_antinode_source": "#221 odd cavity mode",
        "definitions_identical_by_assumption": False,
        "pass": True,
    }


# ========================================================================
# T2. Verify the quarter-wave result, but classify it correctly
# ========================================================================


def test_T2_quarter_wave_benchmark() -> dict:
    c = cavity_data()
    hard_wall_error = abs(c["X_antinode"] - math.pi / 2.0)
    definition_spread = max(
        c["X_antinode"], c["X_amplitude_rms"], c["X_energy_rms"]
    ) - min(c["X_antinode"], c["X_amplitude_rms"], c["X_energy_rms"])
    ok = (
        hard_wall_error < 0.08
        and definition_spread > 0.25
        and c["q_star"] < 1e-12
        and c["p_star"] < 1e-12
    )
    return {
        "name": "T2_quarter_wave_benchmark",
        "description": (
            "The odd cavity mode genuinely carries a robust quarter-wave "
            "invariant: omega*sigma_antinode is close to pi/2.  However, "
            "amplitude-RMS and energy-RMS lengths on the same mode give "
            "different O(1) values, and the orbit is exactly source-decoupled. "
            "The result is therefore classified as a cavity benchmark, not "
            "yet as the soliton mass-map coefficient."
        ),
        "X_antinode": c["X_antinode"],
        "pi_over_2": math.pi / 2.0,
        "hard_wall_error": float(hard_wall_error),
        "X_amplitude_rms": c["X_amplitude_rms"],
        "X_energy_rms": c["X_energy_rms"],
        "definition_spread": float(definition_spread),
        "source_decoupled": bool(c["q_star"] < 1e-12 and c["p_star"] < 1e-12),
        "pass": bool(ok),
    }


# ========================================================================
# T3. Measure the original soliton localization family
# ========================================================================


def test_T3_soliton_localization_family() -> dict:
    s = soliton_localization_scales()
    wave = [s[k] for k in ("r50", "r90", "r95", "rms", "Rstar_A", "Rstar_B")]
    spread = max(wave) / min(wave)
    ok = (
        0 < s["r50"] < s["r90"] < s["r95"]
        and s["rms"] > 0
        and spread > 2.0
    )
    return {
        "name": "T3_soliton_localization_family",
        "description": (
            "The locked #180 soliton supplies several standard, independently "
            "computable IR localization lengths.  They are not numerically "
            "identical: quantile, RMS, and exterior-overlap radii span an "
            "O(1) family.  A physical matching functional is therefore needed "
            "to select which scale enters the #202 suppression law."
        ),
        "scales": s,
        "wave_scale_max_over_min": float(spread),
        "pass": bool(ok),
    }


# ========================================================================
# T4. Direct identifiability test
# ========================================================================


def test_T4_cross_sector_nonuniqueness() -> dict:
    candidates = candidate_cross_sector_coefficients()
    values = list(candidates.values())
    ratio = max(values) / min(values)
    distance_to_quarter_wave = {
        k: float(abs(v - math.pi / 2.0)) for k, v in candidates.items()
    }
    ok = ratio > 2.0 and max(distance_to_quarter_wave.values()) > 0.5
    return {
        "name": "T4_cross_sector_nonuniqueness",
        "description": (
            "Applying the #221 cavity frequency to the pre-existing #203 "
            "soliton length definitions does not produce one coefficient. "
            "The candidate X values retain order-one definition dependence. "
            "Therefore the present equations do not yet identify "
            "sigma_mode with the first cavity antinode."
        ),
        "candidate_X": candidates,
        "candidate_max_over_min": float(ratio),
        "distance_from_pi_over_2": distance_to_quarter_wave,
        "unique_coefficient_selected": False,
        "pass": bool(ok),
    }


# ========================================================================
# T5. The missing dimensionful weld
# ========================================================================


def test_T5_scale_weld_falsification() -> dict:
    """A coordinate-scale test for the absent cross-sector calibration.

    Normalized soliton shape observables are unchanged if its radial coordinate
    is expressed in units rescaled by ``a``.  Until an equation ties those units
    to the Tangherlini radius used by #221, omega*sigma_mode scales with ``a``.
    """
    base = soliton_localization_scales()["rms"]
    omega = cavity_data()["omega"]
    factors = (0.5, 1.0, 2.0)
    mapped = {str(a): float(omega * a * base) for a in factors}
    exact_scaling = abs(mapped["2.0"] / mapped["1.0"] - 2.0)
    unresolved = max(mapped.values()) / min(mapped.values())
    ok = exact_scaling < 1e-12 and unresolved > 3.9
    return {
        "name": "T5_scale_weld_falsification",
        "description": (
            "The soliton and cavity sectors currently use separately normalized "
            "length units.  In the absence of an explicit field equation or "
            "matching condition fixing their ratio, a rescaling of the soliton "
            "radial unit rescales X linearly while leaving normalized soliton "
            "shape data unchanged.  This is the missing cross-sector weld."
        ),
        "coordinate_scale_factors": list(factors),
        "mapped_X_from_rms": mapped,
        "linearity_error": float(exact_scaling),
        "range_ratio": float(unresolved),
        "cross_sector_scale_fixed_by_existing_equation": False,
        "pass": bool(ok),
    }


# ========================================================================
# T6. What the #220 stability information can and cannot do
# ========================================================================


def test_T6_stability_role() -> dict:
    c = cavity_data()
    linear = c["half_amplitude_residual"] < 1e-9
    source_decoupled = c["q_star"] < 1e-12 and c["p_star"] < 1e-12
    ok = linear and source_decoupled
    return {
        "name": "T6_stability_role",
        "description": (
            "The #220 monodromy is valuable as a rejection test for unstable "
            "branches.  In the branch used by #221, however, parity makes the "
            "source vanish and the mode is exactly linear under amplitude "
            "rescaling.  Stability therefore does not normalize the soliton "
            "length or supply the missing 3D-to-5D matching coefficient."
        ),
        "exactly_linear_sector": bool(linear),
        "source_decoupled": bool(source_decoupled),
        "stability_selects_normalization": False,
        "pass": bool(ok),
    }


# ========================================================================
# T7. Executable requirements for a genuine derivation
# ========================================================================


def test_T7_successor_contract() -> dict:
    requirements = [
        "Solve one coupled localized Pin-Dirac/soliton state on the Tangherlini bridge.",
        "Use the rho^3 five-dimensional radial measure in the mode normalization.",
        "Derive the 3D soliton radius to 5D bridge-coordinate map from the action.",
        "Compute the neck-to-asymptotic overlap functional without defining it by an antinode.",
        "Show convergence under radial-grid, outer-boundary, and matching-surface refinement.",
        "Use Floquet stability only to reject unstable branches.",
        "Lock the coefficient before comparing with the observed lepton masses.",
    ]
    return {
        "name": "T7_successor_contract",
        "description": (
            "These are the minimum executable conditions under which the "
            "alpha-to-lepton coefficient would become identifiable rather "
            "than definition-selected."
        ),
        "requirements": requirements,
        "pass": True,
    }


# ========================================================================
# T8. Assessment
# ========================================================================


def test_T8_assessment() -> dict:
    core = all(
        test()["pass"]
        for test in (
            test_T2_quarter_wave_benchmark,
            test_T3_soliton_localization_family,
            test_T4_cross_sector_nonuniqueness,
            test_T5_scale_weld_falsification,
            test_T6_stability_role,
        )
    )
    assessment = (
        "The quarter-wave number is validated as a robust cavity invariant, "
        "but the original #202 soliton localization scale remains a distinct "
        "object.  Existing equations do not supply a unique dimensional weld "
        "or matching functional between the soliton and cavity sectors.  The "
        "O(1) lepton coefficient is therefore not yet identifiable from the "
        "current model; PR #221 should retain pi/2 as a benchmark and narrow "
        "its mass-unification claim until the successor contract is executed."
        if core
        else "INCONCLUSIVE: an audit invariant failed; inspect before quoting."
    )
    return {
        "name": "T8_assessment",
        "description": "logical status of the proposed coefficient",
        "core_green": bool(core),
        "assessment": assessment,
        "pass": bool(core),
    }


def run_probe() -> dict:
    tests = [
        test_T1_definition_ledger(),
        test_T2_quarter_wave_benchmark(),
        test_T3_soliton_localization_family(),
        test_T4_cross_sector_nonuniqueness(),
        test_T5_scale_weld_falsification(),
        test_T6_stability_role(),
        test_T7_successor_contract(),
        test_T8_assessment(),
    ]
    all_ok = all(t["pass"] for t in tests)
    verdict_class = (
        "QUARTER_WAVE_INVARIANT_VALIDATED_BUT_LEPTON_COEFFICIENT_NOT_IDENTIFIABLE_UNTIL_THE_SOLITON_TO_5D_SCALE_WELD_IS_DERIVED"
        if all_ok
        else "LEPTON_O1_IDENTIFIABILITY_AUDIT_INCONCLUSIVE"
    )
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "Independent audit separating the #221 quarter-wave cavity invariant "
            "from the #202/#203 soliton localization scale"
        ),
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": tests[-1]["assessment"],
    }


def main(output_dir: Optional[Path] = None) -> dict:
    result = run_probe()
    print(json.dumps(result, indent=2))
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "probe.json").write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()
