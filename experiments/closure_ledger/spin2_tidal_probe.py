"""
Spin 0 against spin 2 on the same S²: what a tensor wave does differently

> Framing: a spin-2 **field on a fixed** S² — the tensor analogue of the
> scalar wave, not a solution of the 4-D linearised Einstein equations.
> General relativity in 2+1 dimensions has no propagating tensor modes at all.

THE CORRECTION THIS PROBE ENCODES
─────────────────────────────────
`warped_sphere_geometry_probe` displays a SCALAR u(t, θ, φ) extrinsically, as a
radial height.  A metric perturbation is not that kind of object: h_ab is
symmetric and trace-free — spin 2 — and it does not push a surface outward.  It
shears it.  Everything below is the list of consequences, each measured.

THE FIELD AND ITS EQUATION
──────────────────────────
An axisymmetric spin-s field on the unit sphere obeys

    ∂²_t h = ∂²_d h + cot d ∂_d h − (s²/sin²d) h,

with eigenvalue −ℓ(ℓ+1) on ₛY_ℓ0.  So the tensor shares the scalar's dispersion
ω_ℓ = √(ℓ(ℓ+1)) and its t = π refocus; what the centrifugal term does instead is
force h → 0 at both poles.  Writing h = sin²d · q removes that term exactly,

    ∂²_t q = (1/sin⁵d) ∂_d( sin⁵d ∂_d q ) − 6 q,

which is integrated in conservative form on a cell-centred grid: the poles
carry no flux and never appear in a denominator.

WHAT IS CHECKED
───────────────
T2  EXACT MODES.  q = 1 is ℓ = 2, q = cos d is ℓ = 3, q = 7cos²d − 1 is ℓ = 4,
    at ω = √6, √12, √20.  Shape and invariant after many periods.

T3  NO MONOPOLE, NO DIPOLE, NO AMPLITUDE AT A POLE.  h = sin²d·q vanishes at
    both poles for every q.  So the smallest source a spin-2 field admits is a
    RING, and its focus is a ring around the antipode rather than a peak on it
    — the exact opposite of the scalar, which piles up there.  Both measured
    on one clock.

T4  IT IS A TENSOR.  Symmetric, trace-free, and the deformation it produces is
    area-preserving: the ellipse has semi-axes 1 ± εh, so the first-order area
    change vanishes identically and what is left is the ε²h² term, reported
    against its closed form.  A scalar height perturbation changes area at
    first order; that is the visible difference.

T5  SPIN WEIGHT 2.  The strain is h₊cos2β + h_ˣsin2β, so rotating the frame by
    β rotates the pattern by 2β: identical after 180°, inverted after 90°.
    Measured on the field itself, not asserted from the formula.

T6  THE CAUSTIC IS A QUARTER TURN, NOT A FLIP.  One passage through the
    antipodal focus does NOT swap the stretch and compression axes.  It shifts
    the phase by 90° — the outbound waveform is the Hilbert transform of the
    inbound one.  This is the Gouy shift, Maslov index 1, and the probe
    distinguishes it from a flip by correlating against all four candidates.

T7  BUT THE ROUND TRIP DOES.  Two focal passages — the antipode at t = π, home
    at t = 2π — compose to a half turn, and the field returns as minus itself.
    That is where the polarisation axes really do swap.  It is not exact,
    because ω_ℓ = √(ℓ(ℓ+1)) only approaches ℓ + ½; the residual is reported.

SCOPE
─────
The Gouy shift belongs to the wave equation, not to the spin: the scalar has
it too.  What belongs to the spin is the pole behaviour, the ℓ ≥ 2 content, the
two polarisations, and the area-preserving shear.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
    TidalField,
    measure_area_preservation,
    measure_caustic_phase,
    measure_exact_modes,
    measure_node_at_the_focus,
    measure_round_trip_inversion,
    measure_scalar_contrast,
    measure_spin_weight,
)

N_RADIAL = 1200
PULSE_WIDTH = 0.18


def _field() -> TidalField:
    return TidalField(sim=Spin2WaveSim(n=N_RADIAL, pulse_width=PULSE_WIDTH))


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("show what changes when the field on S² is a tensor rather "
                  "than a scalar — spin 2, displayed as tidal shear"),
        "equation": "∂²_t h = ∂²_d h + cot d ∂_d h − 4h/sin²d",
        "substitution": "h = sin²d · q  ⇒  ∂²_t q = (sin⁵d q')'/sin⁵d − 6q",
        "dispersion": "ω_ℓ = √(ℓ(ℓ+1)), the same as the scalar",
        "pass": True,
    }


def t2_exact_modes() -> dict:
    r = measure_exact_modes(n=N_RADIAL, periods=10.0)
    return {
        "name": "T2_exact_modes",
        "periods": r["periods"],
        "modes": r["modes"],
        "worst_shape_error": r["worst_shape_error"],
        "pass": bool(r["worst_shape_error"] < 2e-3
                     and all(m["energy_drift"] < 1e-9 for m in r["modes"])),
    }


def t3_no_amplitude_at_a_pole() -> dict:
    f = _field()
    node = measure_node_at_the_focus(f, frames=220)
    contrast = measure_scalar_contrast(f, frames=220)
    return {
        "name": "T3_no_amplitude_at_a_pole",
        **{k: node[k] for k in ("peak_distance", "peak_time",
                                "ring_radius_at_focus",
                                "amplitude_at_the_antipode",
                                "is_a_ring_not_a_peak")},
        "scalar_peak_distance": contrast["scalar_peak_distance"],
        "scalar_sits_on_the_antipode": contrast["scalar_sits_on_the_antipode"],
        "lowest_multipole": 2,
        "pass": bool(node["is_a_ring_not_a_peak"]
                     and contrast["scalar_sits_on_the_antipode"]
                     and node["amplitude_at_the_antipode"]
                     < 0.02 * node["peak_amplitude"]),
    }


def t4_area_preserving_shear() -> dict:
    r = measure_area_preservation(_field(), eps=0.05)
    return {
        "name": "T4_area_preserving_shear",
        **{k: r[k] for k in ("trace", "symmetric", "amplitude", "area_ratio",
                             "first_order_area_change",
                             "second_order_prediction",
                             "area_preserved_to_first_order")},
        "pass": bool(r["area_preserved_to_first_order"]
                     and abs(r["trace"]) < 1e-14),
    }


def t5_spin_weight() -> dict:
    r = measure_spin_weight(_field())
    return {
        "name": "T5_spin_weight",
        **{k: r[k] for k in ("distance", "repeats_after_180_deg",
                             "differs_after_90_deg", "consistent")},
        "pass": bool(r["consistent"]),
    }


def t6_caustic_is_a_quarter_turn() -> dict:
    r = measure_caustic_phase(_field(), frames=2400)
    return {
        "name": "T6_caustic_is_a_quarter_turn",
        **{k: r[k] for k in ("ring_radius", "correlations", "best_match",
                             "best_correlation",
                             "is_a_quarter_turn_not_a_flip", "maslov_index")},
        "pass": bool(r["is_a_quarter_turn_not_a_flip"]),
    }


def t7_round_trip_inverts() -> dict:
    r = measure_round_trip_inversion(_field())
    return {
        "name": "T7_round_trip_inverts",
        **{k: r[k] for k in ("corr_after_one_round_trip_with_minus_start",
                             "corr_after_one_round_trip_with_start",
                             "corr_after_two_round_trips_with_start",
                             "inversion_residual", "inverts")},
        "pass": bool(r["inverts"]),
    }


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_exact_modes(), t3_no_amplitude_at_a_pole(),
             t4_area_preserving_shear(), t5_spin_weight(),
             t6_caustic_is_a_quarter_turn(), t7_round_trip_inverts()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t5, t6, t7 = tests[1], tests[2], tests[3], tests[4], tests[5], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "A_TENSOR_WAVE_CANNOT_SIT_ON_ITS_OWN_FOCUS"
        verdict = (
            "A TENSOR WAVE CANNOT SIT ON ITS OWN FOCUS. Replacing the scalar "
            "with a spin-2 field on the same sphere changes the picture "
            "structurally, not cosmetically, and the changes are measured "
            "rather than drawn.\n\n"
            "THE SOLVER IS EXACT WHERE IT CAN BE CHECKED. Three closed-form "
            "modes — q = 1 at ℓ = 2, q = cos d at ℓ = 3, q = 7cos²d − 1 at "
            f"ℓ = 4 — keep their shape to {t2['worst_shape_error']:.1e} after "
            f"{t2['periods']:.0f} periods at ω = √6, √12, √20, with the "
            "invariant conserved to round-off.\n\n"
            "NO AMPLITUDE AT A POLE. h = sin²d·q vanishes at both poles for "
            "every q, so a spin-2 field cannot be nonzero where the frame "
            "degenerates. At the focus it is therefore a RING of radius "
            f"{t3['ring_radius_at_focus']:.3f} about the antipode, with "
            f"{t3['amplitude_at_the_antipode']:.1e} on the antipode itself — "
            "while the scalar, on the same clock, piles up exactly there at "
            f"d = {t3['scalar_peak_distance']:.4f}. The same fact at the other "
            "end says the smallest source a spin-2 field admits is a ring: "
            "there is no such thing as a point source of tidal shear.\n\n"
            "PURE SHEAR, NOT BREATHING. The tensor is symmetric and trace-free "
            f"({t4['trace']:.1e}), so the ellipse it makes has semi-axes "
            "1 ± εh and its first-order area change vanishes identically: "
            f"measured {t4['first_order_area_change']:.2e} against the "
            f"second-order prediction {t4['second_order_prediction']:.2e}. A "
            "radial height perturbation changes area at first order — that is "
            "the difference you can see.\n\n"
            "SPIN WEIGHT 2, DIRECTLY. The strain pattern is identical after a "
            f"180° rotation of the frame ({t5['repeats_after_180_deg']:.1e}) "
            f"and inverted after 90° ({t5['differs_after_90_deg']:.2f}× the "
            "amplitude). That factor of two in the angle is the spin weight, "
            "visible without any decomposition.\n\n"
            "THE CAUSTIC IS A QUARTER TURN, NOT A FLIP. The obvious guess — "
            "that passing the antipodal focus swaps the stretch and "
            "compression axes — is not what happens. Correlating the outbound "
            "waveform against all four candidates gives "
            f"{t6['correlations']['phase_-90']:+.3f} for a −90° phase shift "
            f"against {t6['correlations']['inverted']:+.3f} for an inversion: "
            "the outbound front is the HILBERT TRANSFORM of the inbound one. "
            "That is the Gouy shift, Maslov index 1, and it is a property of "
            "passing through a focus rather than of the spin — the scalar has "
            "it too.\n\n"
            "THE ROUND TRIP IS WHERE THE AXES SWAP. Two focal passages, the "
            "antipode and then home, compose to a half turn: at t = 2π the "
            "field is minus what it started as, correlation "
            f"{t7['corr_after_one_round_trip_with_minus_start']:.4f}, and "
            "after two round trips it is itself again "
            f"({t7['corr_after_two_round_trips_with_start']:.4f}). The "
            f"residual {t7['inversion_residual']:.4f} is the dispersion left "
            "over because ω_ℓ = √(ℓ(ℓ+1)) only approaches ℓ + ½.\n\n"
            "SCOPE. A spin-2 field on a FIXED S², not a solution of the 4-D "
            "linearised Einstein equations — 2+1 gravity has no propagating "
            "tensor modes at all. What is faithful is the polarisation "
            "structure: two components, spin weight 2, ℓ ≥ 2 only, "
            "area-preserving shear, and the behaviour at a focus."
        )
    else:
        verdict_class = "SPIN2_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the exact modes, the "
                   "pole behaviour, the trace, the spin weight, or the "
                   "caustic correlations.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a spin-2 wave on S² displayed as tidal shear, against the scalar "
            "displayed as radial height, on one clock"
        ),
        "the_poles": "h = sin²d·q vanishes there — no point source, no peak at the focus",
        "the_shear": "symmetric, trace-free, area-preserving to first order",
        "the_spin": "strain ∝ cos 2β — identical after 180°, inverted after 90°",
        "the_caustic": "a 90° Gouy shift per passage, not a polarisation flip",
        "the_round_trip": "two passages compose to −1: there the axes do swap",
        "geometry": {
            "n_radial": N_RADIAL, "pulse_width": PULSE_WIDTH,
            "antipodal_time": ANTIPODAL_TIME, "return_time": RETURN_TIME,
            "lowest_multipole": 2,
        },
        "tests": tests,
        "n_passed": tests[-1]["n_passed"],
        "n_total": tests[-1]["n_total"],
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    t = {x["name"]: x for x in s["tests"]}
    out: List[str] = []
    out.append("# Spin 0 against spin 2 on the same S²\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> A spin-2 **field on a fixed** `S²` — the tensor analogue of "
               "the scalar wave, not linearised GR on a spacetime.\n")

    c = t["T2_exact_modes"]
    out.append("## The solver, against exact modes\n")
    out.append("| `q` | `ℓ` | `ω` | shape error | invariant drift |")
    out.append("|---|---:|---:|---:|---:|")
    names = {2: "`1`", 3: "`cos d`", 4: "`7cos²d − 1`"}
    for m in c["modes"]:
        out.append(f"| {names[m['ell']]} | {m['ell']} | {m['omega']:.4f} | "
                   f"`{m['shape_error']:.1e}` | `{m['energy_drift']:.1e}` |")
    out.append(f"\nAfter {c['periods']:.0f} periods.\n")

    c = t["T3_no_amplitude_at_a_pole"]
    out.append("## A spin-2 field cannot sit on a pole\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| tensor peak at distance | {c['peak_distance']:.4f} |")
    out.append(f"| — i.e. a ring of radius | {c['ring_radius_at_focus']:.4f} |")
    out.append(f"| amplitude *on* the antipode | `{c['amplitude_at_the_antipode']:.1e}` |")
    out.append(f"| **scalar** peak at distance | {c['scalar_peak_distance']:.4f} |")
    out.append(f"| lowest multipole | ℓ = {c['lowest_multipole']} |\n")

    c = t["T4_area_preserving_shear"]
    out.append("## Pure shear, not breathing\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| trace | `{c['trace']:.1e}` |")
    out.append(f"| area ratio | {c['area_ratio']:.10f} |")
    out.append(f"| first-order area change | `{c['first_order_area_change']:.2e}` |")
    out.append(f"| `ε²h²/2` prediction | `{c['second_order_prediction']:.2e}` |\n")

    c = t["T5_spin_weight"]
    out.append("## Spin weight 2\n")
    out.append("| rotation of the frame | strain |")
    out.append("|---|---:|")
    out.append(f"| 180° | identical, `{c['repeats_after_180_deg']:.1e}` |")
    out.append(f"| 90° | inverted, {c['differs_after_90_deg']:.2f}× amplitude |\n")

    c = t["T6_caustic_is_a_quarter_turn"]
    out.append("## The caustic: a quarter turn, not a flip\n")
    out.append("| the outbound front is the inbound one… | correlation |")
    out.append("|---|---:|")
    for k, v in c["correlations"].items():
        out.append(f"| {k.replace('_', ' ')} | {v:+.4f} |")
    out.append(f"\nBest match: **{c['best_match']}** — the Gouy shift, Maslov "
               f"index {c['maslov_index']}.\n")

    c = t["T7_round_trip_inverts"]
    out.append("## The round trip, where the axes do swap\n")
    out.append("| | correlation |")
    out.append("|---|---:|")
    out.append(f"| `h(2π)` with `−h(0)` | {c['corr_after_one_round_trip_with_minus_start']:+.4f} |")
    out.append(f"| `h(2π)` with `h(0)` | {c['corr_after_one_round_trip_with_start']:+.4f} |")
    out.append(f"| `h(4π)` with `h(0)` | {c['corr_after_two_round_trips_with_start']:+.4f} |")
    out.append(f"| inversion residual | {c['inversion_residual']:.4f} |\n")

    out.append("## Verdict\n")
    out.append(f"**{s['n_passed']}/{s['n_total']} checks passed.**\n")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_spin2_tidal_probe"
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
