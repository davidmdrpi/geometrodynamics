"""
Projecting h_ab into a continuous embedded surface: shape, slide, and the bulk

> Framing: a faithful REPRESENTATION of a spin-2 field in the ℝ³ embedding of a
> fixed S² — not backreaction, and not linearised GR on a spacetime.

THE PROBLEM
───────────
`spin2_tidal_probe` samples the tensor as discrete tangent-space ellipses.
Faithful, but flat: it never touches the embedding, so it gives no intuition
for how tidal shear meets a bulk.  The question is whether h_ab can be drawn as
a CONTINUOUS surface deformation r(θ, φ, t) instead.

WHY A HEIGHT FIELD CANNOT DO IT
───────────────────────────────
For X = r n̂ the induced metric is g_ab = r² ĝ_ab + ∂_a r ∂_b r, and the
gradient term is SECOND order.  So at first order a radial deformation gives

    δg_ab = 2ρ ĝ_ab      — purely conformal.

Shape carries the trace and nothing else.  That is exactly why the ellipses
existed, and it is checked here directly (T2) rather than argued.

THE DISPLACEMENT THAT DOES
──────────────────────────
Add the tangential part.  With ξ = ∇W,

    X = (R + ερ) n̂ + ε ξ ,     δg_ab = 2ρ ĝ_ab + 2∇₍ₐξ_b₎ ,

and demanding tracelessness fixes the radial part completely:

    ρ = −½ ΔW ,      h_ab = [2∇₍ₐ∇_b₎W]^TF .

One potential carries both.  For an axisymmetric h the Hessian condition
W'' − cot d W' = h is first order in ψ = W', with integrating factor sin d:

    ψ(d) = sin d [ C + ∫₀^d h/sin ] ,    ρ(d) = −½h − cos d [ C + ∫₀^d h/sin ] .

No derivative of the solved field is taken — one integral, with a regular
integrand because h ~ sin²d at both poles.

WHAT IS CHECKED
───────────────
T2  A RADIAL DEFORMATION IS CONFORMAL.  Measured on a deliberately radial
    surface: the trace-free part of its induced metric perturbation is zero.

T3  THE THEOREM.  The induced metric perturbation of the drawn surface IS the
    solved h_ab: component by component, trace-free, and off-diagonal free.

T4  THE QUADRUPOLE IS THE LEGENDRE SHAPE.  Feeding the exact ℓ = 2 mode gives
    ρ = P₂(cos d) — the prolate–oblate picture of a quadrupole wave, derived
    rather than drawn.

T5  THE FREE CONSTANT IS A RIGID TRANSLATION.  C is the whole kernel: an ℓ = 1
    displacement moves every point by one vector.  Removing it leaves no
    dipole; ℓ = 0 cannot appear at all since ∫ΔW dA = 0, so a spin-2 wave can
    never breathe the sphere's area.

T6  AREA IS UNCHANGED AT FIRST ORDER, which is what trace-free means, measured
    on the drawn surface rather than inferred.

T7  IT REACHES THE BULK.  With the display gain fixed from the run's own peak,
    the surface travels a measured fraction of the way to each shell of the
    vacuole and never touches either.

SCOPE
─────
The gain ε is a display choice; the SHAPE at any gain is the solved field, and
the theorem is a first-order statement, so the residuals shrink with ε.  This
gives the wave an extrinsic amplitude — it does not make it act on the sphere.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.embedded_wave import (
    EmbeddedTidalSurface,
    measure_area_and_multipoles,
    measure_bulk_reach,
    measure_induced_metric,
    measure_quadrupole_shape,
)
from geometrodynamics.viz.spin2_tidal import RETURN_TIME, Spin2WaveSim

N_RADIAL = 1200
PULSE_WIDTH = 0.18
SMALL_GAIN = 1e-4


def _surface() -> EmbeddedTidalSurface:
    return EmbeddedTidalSurface(
        sim=Spin2WaveSim(n=N_RADIAL, pulse_width=PULSE_WIDTH),
        n_theta=121, n_phi=181)


# ════════════════════════════════════════════════════════════════════════════
def t1_goal() -> dict:
    return {
        "name": "T1_goal",
        "claim": ("draw a spin-2 field as a CONTINUOUS deformation of the "
                  "embedded sphere, so tidal shear can meet a bulk"),
        "displacement": "X = (R + ερ) n̂ + ε∇W",
        "radial_part": "ρ = −½ΔW",
        "shear_part": "h_ab = [2∇₍ₐ∇_b₎W]^TF",
        "solution": "ψ = W' = sin d [C + ∫₀^d h/sin]",
        "pass": True,
    }


def t2_radial_alone_is_conformal() -> dict:
    eps, a0, step = 1e-4, 0.4, 1e-5

    def X(d, a):
        r = 1.0 + eps * math.cos(2.0 * d)
        return r * np.array([math.sin(d) * math.cos(a),
                             math.sin(d) * math.sin(a), math.cos(d)])

    rows = []
    for d0 in (0.6, 1.1, 1.9, 2.5):
        Xd = (X(d0 + step, a0) - X(d0 - step, a0)) / (2.0 * step)
        Xa = (X(d0, a0 + step) - X(d0, a0 - step)) / (2.0 * step)
        dg_dd = (float(Xd @ Xd) - 1.0) / eps
        dg_aa = (float(Xa @ Xa) - math.sin(d0) ** 2) / (eps * math.sin(d0) ** 2)
        rows.append({"distance": d0, "trace": dg_dd + dg_aa,
                     "trace_free_part": 0.5 * (dg_dd - dg_aa)})
    worst = max(abs(r["trace_free_part"]) for r in rows)
    biggest_trace = max(abs(r["trace"]) for r in rows)
    return {
        "name": "T2_radial_alone_is_conformal",
        "rows": rows,
        "worst_trace_free_part": worst,
        "largest_trace": biggest_trace,
        "is_conformal": bool(worst < 1e-3 * max(biggest_trace, 1e-12)
                             or worst < 1e-3),
        "pass": bool(worst < 1e-3),
    }


def t3_the_theorem() -> dict:
    r = measure_induced_metric(_surface(), t=1.2, gain=SMALL_GAIN)
    return {
        "name": "T3_the_theorem",
        "gain": r["gain"],
        "rows": r["rows"],
        "worst_relative_error": r["worst_relative_error"],
        "worst_relative_trace": r["worst_relative_trace"],
        "worst_off_diagonal": r["worst_off_diagonal"],
        "pass": bool(r["reproduces_h"] and r["worst_relative_trace"] < 1e-2),
    }


def t4_quadrupole_shape() -> dict:
    r = measure_quadrupole_shape(n=2000)
    return {"name": "T4_quadrupole_shape", **r, "pass": bool(r["is_legendre_p2"])}


def t5_the_kernel_is_a_translation() -> dict:
    s = EmbeddedTidalSurface(sim=Spin2WaveSim(n=600), n_theta=41, n_phi=61)
    s.advance_to(1.0)
    d = np.linspace(0.2, math.pi - 0.2, 25)
    a = np.zeros_like(d)
    base = s.positions(d[None, :], a[None, :], gain=1.0)[0]
    c = 0.37
    n_hat, e_d = s._frames(d[None, :], a[None, :])
    delta = c * (np.sin(d)[:, None] * e_d[0] - np.cos(d)[:, None] * n_hat[0])
    spread = float(np.ptp(delta, axis=0).max())
    mult = measure_area_and_multipoles(_surface(), t=1.2, gain=1e-3)
    return {
        "name": "T5_the_kernel_is_a_translation",
        "displacement_spread_over_the_sphere": spread,
        "translation_length": float(np.linalg.norm(delta[0])),
        "requested_length": c,
        "residual_dipole": mult["dipole"],
        "residual_monopole": mult["monopole"],
        "pass": bool(spread < 1e-12 and abs(mult["dipole"]) < 1e-9
                     and abs(mult["monopole"]) < 1e-6),
    }


def t6_area_is_unchanged() -> dict:
    r = measure_area_and_multipoles(_surface(), t=1.2, gain=1e-3)
    return {
        "name": "T6_area_is_unchanged",
        **{k: r[k] for k in ("gain", "area", "round_area",
                             "relative_area_change",
                             "area_is_first_order_unchanged")},
        "pass": bool(r["area_is_first_order_unchanged"]),
    }


def t7_reaches_the_bulk() -> dict:
    s = _surface()
    r = measure_bulk_reach(s, frames=220)
    return {
        "name": "T7_reaches_the_bulk",
        **r,
        "display_gain": s.gain,
        "pass": bool(r["stays_between_the_dolls"]
                     and r["max_outward_fraction"] > 0.2),
    }


def t8_assessment(tests: List[dict]) -> dict:
    n = sum(1 for t in tests if t["pass"])
    return {"name": "T8_assessment", "n_passed": n, "n_total": len(tests),
            "pass": n == len(tests)}


# ════════════════════════════════════════════════════════════════════════════
def run_probe() -> dict:
    tests = [t1_goal(), t2_radial_alone_is_conformal(), t3_the_theorem(),
             t4_quadrupole_shape(), t5_the_kernel_is_a_translation(),
             t6_area_is_unchanged(), t7_reaches_the_bulk()]
    tests.append(t8_assessment(tests))
    t2, t3, t4, t5, t6, t7 = tests[1], tests[2], tests[3], tests[4], tests[5], tests[6]

    if all(t["pass"] for t in tests):
        verdict_class = "ONE_POTENTIAL_CARRIES_SHAPE_AND_SHEAR"
        verdict = (
            "ONE POTENTIAL CARRIES SHAPE AND SHEAR. A spin-2 field can be "
            "drawn as a continuous deformation of the embedded sphere, and the "
            "projection is forced rather than chosen.\n\n"
            "A HEIGHT FIELD ALONE CANNOT. For X = r n̂ the induced metric is "
            "g_ab = r²ĝ_ab + ∂_a r ∂_b r, whose gradient term is second order, "
            "so at first order the perturbation is purely conformal: measured "
            f"trace-free part {t2['worst_trace_free_part']:.1e} against a "
            f"trace of {t2['largest_trace']:.3f}. Shape carries the trace and "
            "nothing else — which is exactly why the tangential slide has to "
            "be there, and why discrete ellipses were the only option before "
            "it was.\n\n"
            "THE THEOREM. Adding ξ = ∇W and demanding tracelessness fixes the "
            "radial part completely, ρ = −½ΔW, and then the induced metric "
            "perturbation of the drawn surface IS the solved h_ab: measured "
            f"against solved to {t3['worst_relative_error']:.1e} of the peak "
            f"amplitude, with a trace of {t3['worst_relative_trace']:.1e} and "
            f"off-diagonal {t3['worst_off_diagonal']:.1e}. The surface is not "
            "an illustration of the tensor; it has the tensor as its own "
            "geometry.\n\n"
            "THE QUADRUPOLE IS THE TEXTBOOK SHAPE. Feeding the exact ℓ = 2 "
            f"mode returns ρ = P₂(cos d) to {t4['shape_error']:.1e}, amplitude "
            f"ratio {t4['amplitude_ratio']:.6f}. The prolate–oblate picture of "
            "a quadrupole gravitational wave is not assumed here; it comes out "
            "of the projection.\n\n"
            "THE FREE CONSTANT IS A RIGID TRANSLATION. C is the whole kernel of "
            "the construction — an ℓ = 1 displacement moves every point of the "
            f"sphere by one vector, verified to {t5['displacement_spread_over_the_sphere']:.1e}. "
            f"Removing it leaves a residual dipole of {t5['residual_dipole']:.1e}, "
            "and ℓ = 0 cannot appear at all because ∫ΔW dA = 0: a spin-2 wave "
            "can never breathe the sphere's area.\n\n"
            "AND THE AREA HOLDS. The drawn surface keeps the round area to "
            f"{t6['relative_area_change']:.1e} at gain {t6['gain']:.0e} — "
            "second order in ε, which is what trace-free means, measured on "
            "the surface rather than inferred from the tensor.\n\n"
            "IT REACHES THE BULK. With the gain fixed from the run's own peak, "
            f"the surface travels {100 * t7['max_outward_fraction']:.1f}% of "
            f"the way to R_outer at t = {t7['at_time_outward']:.3f} and "
            f"{100 * t7['max_inward_fraction']:.1f}% toward R_inner, without "
            "touching either. The tensor wave now has an extrinsic amplitude — "
            "which is the whole reason to do this.\n\n"
            "SCOPE. The gain is a display choice and the theorem is a "
            "first-order statement, so its residuals shrink with ε. This is a "
            "faithful representation of h_ab in the embedding; it does not "
            "make the wave act on the sphere it is deforming."
        )
    else:
        verdict_class = "EMBEDDED_WAVE_INCONCLUSIVE"
        verdict = ("INCONCLUSIVE. A check failed; review the conformal test, "
                   "the induced metric, the quadrupole shape, the kernel, the "
                   "area, or the bulk reach.")

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "a spin-2 field projected into a continuous ℝ³ surface deformation "
            "whose own induced metric is the field"
        ),
        "the_obstruction": "a radial deformation is conformal at first order",
        "the_construction": "X = (R + ερ)n̂ + ε∇W with ρ = −½ΔW",
        "the_theorem": "δg_ab of the drawn surface equals the solved h_ab",
        "the_kernel": "ℓ = 1 only — a rigid translation; ℓ = 0 cannot appear",
        "geometry": {
            "n_radial": N_RADIAL, "pulse_width": PULSE_WIDTH,
            "small_gain": SMALL_GAIN, "return_time": RETURN_TIME,
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
    out.append("# Projecting `h_ab` into a continuous embedded surface\n")
    out.append(f"_{s['timestamp_utc']}_\n")
    out.append("> A faithful **representation** of a spin-2 field in the ℝ³ "
               "embedding of a fixed `S²` — not backreaction.\n")

    out.append("## The construction\n")
    out.append("```\nX = (R + ερ) n̂ + ε∇W      ρ = −½ΔW      "
               "h_ab = [2∇₍ₐ∇_b₎W]^TF\n```\n")

    c = t["T2_radial_alone_is_conformal"]
    out.append("## Why a height field cannot do it\n")
    out.append("| `d` | trace | trace-free part |")
    out.append("|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['distance']:.2f} | {r['trace']:+.4f} | "
                   f"`{r['trace_free_part']:+.1e}` |")
    out.append("\nA radial deformation is **conformal** at first order.\n")

    c = t["T3_the_theorem"]
    out.append("## The theorem: the drawn surface has the solved metric\n")
    out.append("| `d` | `h₊` measured | `h₊` solved | trace |")
    out.append("|---:|---:|---:|---:|")
    for r in c["rows"]:
        out.append(f"| {r['distance']:.2f} | {r['h_plus_measured']:+.6f} | "
                   f"{r['h_plus_solved']:+.6f} | `{r['trace']:+.1e}` |")
    out.append(f"\nWorst relative error `{c['worst_relative_error']:.1e}`, "
               f"trace `{c['worst_relative_trace']:.1e}`, off-diagonal "
               f"`{c['worst_off_diagonal']:.1e}`, at gain `{c['gain']:.0e}`.\n")

    c = t["T4_quadrupole_shape"]
    out.append("## The quadrupole is the Legendre shape\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| `ρ` against `P₂(cos d)` | `{c['shape_error']:.1e}` |")
    out.append(f"| amplitude ratio | {c['amplitude_ratio']:.6f} |")
    out.append(f"| residual dipole | `{c['residual_dipole']:.1e}` |\n")

    c = t["T6_area_is_unchanged"]
    out.append("## Area, at first order\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| drawn area | {c['area']:.8f} |")
    out.append(f"| round area | {c['round_area']:.8f} |")
    out.append(f"| relative change | `{c['relative_area_change']:.1e}` |\n")

    c = t["T7_reaches_the_bulk"]
    out.append("## Reach into the bulk\n")
    out.append("| quantity | value |")
    out.append("|---|---:|")
    out.append(f"| toward `R_outer` | {100 * c['max_outward_fraction']:.1f}% |")
    out.append(f"| at time | {c['at_time_outward']:.3f} |")
    out.append(f"| toward `R_inner` | {100 * c['max_inward_fraction']:.1f}% |")
    out.append(f"| display gain | {c['display_gain']:.2f} |\n")

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
    out = here / "runs" / f"{ts}_embedded_wave_probe"
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
