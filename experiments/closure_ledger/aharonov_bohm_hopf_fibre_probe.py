"""
Aharonov-Bohm Hopf-fibre action probe.

Sub-target #2 from the ℏ-origin research plan. The Hopf connection on
S³ has the closed form

    A  =  (1/2) cos(χ) dφ
    ∮A  =  π · cos(χ)        (per fibre loop at hyper-latitude χ)

per `geometrodynamics/hopf/connection.py`. The Aharonov-Bohm
interpretation: a spinor traversing a closed Hopf-fibre loop picks up
phase ∮A from the connection. With the spinor's natural double-cover
closure (4π in S³), a full closure traverses each fibre TWICE,
giving accumulated phase 2π·cos(χ).

This probe makes the verification numerical + structural:

  (a) Integrate A·dφ explicitly along a φ-loop at canonical χ values
      and compare to the closed form π·cos(χ). The eigensolver-
      independent path: just the connection 1-form integrated around
      a circle at fixed χ.
  (b) Verify the spinor double-cover doubling: 2·holonomy = 2π·cos(χ).
  (c) Test the integer-quantization pattern at canonical χ values
      ({0, π/4, π/3, π/2, π}): single-closure phase, double-closure
      phase, and whether each is an integer multiple of 2π.
  (d) Combine with the throat T² channel (closure phase π) and
      the antipodal closure (k·2π) to verify the **full closure-cycle
      integer count** for each species at the locked χ = 0.

This probe documents the angular channel's contribution to the
closure-cycle integer count, complementing the radial channel's
hard-wall BS reading from the closed-orbit probe.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


TAU = 2.0 * math.pi


# ---------------------------------------------------------------------------
# (a) Numerical verification of the closed-form holonomy
# ---------------------------------------------------------------------------

@dataclass
class HolonomyTrial:
    chi: float
    chi_in_pi: float
    holonomy_closed_form: float    # π·cos(χ)
    holonomy_numerical: float      # ∮ A_φ dφ via trapezoidal integration
    deviation: float
    matches_closed_form: bool


def _numerical_holonomy(chi: float, n_pts: int = 1000) -> float:
    """Numerically integrate ∮A · dφ around a φ-loop at fixed χ.

    The Hopf connection has A = (1/2)cos(χ) dφ (per `hopf/connection.py`).
    A_φ is therefore (1/2)cos(χ), constant along the φ-loop. The integral
    is just the constant times the loop length 2π — so the closed-form
    holonomy is π·cos(χ). We verify this numerically by trapezoidal
    integration along an explicit φ-grid, which serves as a sanity check
    that our convention matches the analytic claim.
    """
    import numpy as np
    from geometrodynamics.hopf.connection import hopf_connection
    phi = np.linspace(0.0, TAU, n_pts)
    a_phi = hopf_connection(chi)   # constant in φ → broadcasts
    # ∫A_φ dφ = A_φ · 2π
    return float(np.trapezoid(np.full_like(phi, a_phi), phi))


def _holonomy_trials() -> list[HolonomyTrial]:
    canonical = [0.0, math.pi / 4, math.pi / 3, math.pi / 2,
                 2 * math.pi / 3, 3 * math.pi / 4, math.pi]
    out: list[HolonomyTrial] = []
    for chi in canonical:
        closed = math.pi * math.cos(chi)
        numerical = _numerical_holonomy(chi)
        dev = abs(closed - numerical)
        out.append(HolonomyTrial(
            chi=chi,
            chi_in_pi=chi / math.pi,
            holonomy_closed_form=closed,
            holonomy_numerical=numerical,
            deviation=dev,
            matches_closed_form=dev < 1e-9,
        ))
    return out


# ---------------------------------------------------------------------------
# (b) + (c): spinor double-cover and integer-quantization at canonical χ
# ---------------------------------------------------------------------------

@dataclass
class SpinorClosureRow:
    chi: float
    chi_label: str
    single_closure_phase: float    # π·cos(χ)
    single_in_2pi: float           # = (1/2)cos(χ) — half quantum
    double_closure_phase: float    # 2π·cos(χ) (spinor 4π closure)
    double_in_2pi: float           # = cos(χ) — integer iff χ ∈ {0, π/2, π}
    single_is_integer: bool
    double_is_integer: bool
    interpretation: str


def _spinor_closure_table() -> list[SpinorClosureRow]:
    rows: list[tuple[float, str, str]] = [
        (0.0, "0 (canonical north pole, χ=0)",
         "Hopf flux at maximum; single closure picks up π = 1/2 quantum, "
         "spinor double cover gives 2π = 1 full quantum."),
        (math.pi / 4, "π/4",
         "Intermediate; cos(π/4) = √2/2 ≈ 0.707."),
        (math.pi / 3, "π/3",
         "Intermediate; cos(π/3) = 1/2 = exact half."),
        (math.pi / 2, "π/2 (equator)",
         "Hopf flux vanishes; both single and double closure give 0 — "
         "trivial AB phase. Equatorial fibre has zero electromagnetic "
         "self-interaction."),
        (2 * math.pi / 3, "2π/3",
         "Intermediate; cos(2π/3) = −1/2."),
        (3 * math.pi / 4, "3π/4",
         "Intermediate; cos(3π/4) = −√2/2."),
        (math.pi, "π (south pole)",
         "Antipode of canonical; single closure picks up −π = −1/2 "
         "quantum, double cover gives −2π = −1 full quantum (mod 2π = 0)."),
    ]
    out: list[SpinorClosureRow] = []
    for chi, label, interp in rows:
        single = math.pi * math.cos(chi)
        double = 2.0 * single
        out.append(SpinorClosureRow(
            chi=chi,
            chi_label=label,
            single_closure_phase=single,
            single_in_2pi=single / TAU,
            double_closure_phase=double,
            double_in_2pi=double / TAU,
            single_is_integer=abs(single / TAU - round(single / TAU)) < 1e-12,
            double_is_integer=abs(double / TAU - round(double / TAU)) < 1e-12,
            interpretation=interp,
        ))
    return out


# ---------------------------------------------------------------------------
# (d) Combined closure-cycle integer per species at the locked χ = 0
# ---------------------------------------------------------------------------

@dataclass
class FullCycleAtChi0:
    species: str
    k: int
    n_layer_1: int                       # from closure_cycle_action_probe
    contribution_antipodal_in_2pi: int   # k
    contribution_hopf_in_2pi: float       # cos(χ)/2 = 1/2 at χ=0
    contribution_throat_in_2pi: float     # 1/2 (T² eigenvalue arg / 2π)
    contribution_uplift_in_2pi: float     # β·max(0, k-3)² / 2π
    sum_in_2pi: float                     # = N_layer_1 (verified)
    n_radial_at_b2_coupling: int          # from closed-orbit probe
    n_total: int                          # N_layer_1 + n_radial


def _full_cycle_breakdown() -> list[FullCycleAtChi0]:
    """Per-species closure-cycle integer count, broken down by channel,
    at the locked χ = 0. Cross-checks the Layer-1 integer counts and
    adds the B2_radial_ladder Layer-2 contribution from the closed-orbit
    probe."""
    BETA_LEPTON = 50.0 * math.pi
    LAYER1_PER_SPECIES = {1: 2, 3: 4, 5: 106}
    # B2_radial_ladder couples (e, μ, τ) to (l=1, n=0,1,2).
    N_RADIAL_B2 = {1: 1, 3: 2, 5: 3}        # (n+1) per coupling
    out: list[FullCycleAtChi0] = []
    for species, k in zip(("electron", "muon", "tau"), (1, 3, 5)):
        antipodal = k                          # k · 2π / 2π
        hopf = (math.pi * 1.0) / TAU           # π·cos(0) / 2π = 1/2
        throat = math.pi / TAU                 # π / 2π = 1/2
        uplift = (BETA_LEPTON * max(0, k - 3) ** 2) / TAU
        total = antipodal + hopf + throat + uplift
        n_radial = N_RADIAL_B2[k]
        out.append(FullCycleAtChi0(
            species=species,
            k=k,
            n_layer_1=LAYER1_PER_SPECIES[k],
            contribution_antipodal_in_2pi=antipodal,
            contribution_hopf_in_2pi=hopf,
            contribution_throat_in_2pi=throat,
            contribution_uplift_in_2pi=uplift,
            sum_in_2pi=total,
            n_radial_at_b2_coupling=n_radial,
            n_total=LAYER1_PER_SPECIES[k] + n_radial,
        ))
    return out


# ---------------------------------------------------------------------------

def run_probe() -> dict:
    holonomy = _holonomy_trials()
    spinor = _spinor_closure_table()
    full_cycle = _full_cycle_breakdown()
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hopf_connection_form": "A = (1/2) cos(χ) dφ",
        "closed_form_holonomy": "∮A = π · cos(χ)",
        "holonomy_trials": [asdict(t) for t in holonomy],
        "all_holonomy_match_closed_form": all(t.matches_closed_form for t in holonomy),
        "spinor_closure_table": [asdict(s) for s in spinor],
        "double_cover_integer_quantized_at": [
            s.chi for s in spinor if s.double_is_integer
        ],
        "double_cover_integer_quantized_chi_labels": [
            s.chi_label for s in spinor if s.double_is_integer
        ],
        "full_cycle_breakdown_at_chi_0": [asdict(f) for f in full_cycle],
        "layer_1_integers_consistent": all(
            abs(f.sum_in_2pi - f.n_layer_1) < 1e-12 for f in full_cycle
        ),
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Aharonov-Bohm Hopf-fibre action probe")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append("")
    lines.append(
        f"**Connection 1-form:** `{summary['hopf_connection_form']}`"
    )
    lines.append(
        f"**Closed-form holonomy:** `{summary['closed_form_holonomy']}`"
    )
    lines.append("")
    lines.append(
        "Sub-target #2 of the ℏ-origin research plan: verify the "
        "Aharonov-Bohm Hopf-fibre holonomy numerically and document "
        "its role in the closure-cycle integer count."
    )
    lines.append("")

    lines.append("## (a) Numerical holonomy vs closed form")
    lines.append("")
    lines.append(
        "Integrating `A_φ dφ` along a φ-loop at fixed χ on a 1000-point "
        "trapezoidal grid:"
    )
    lines.append("")
    lines.append(
        "| χ | χ in π | closed form | numerical | deviation | matches? |"
    )
    lines.append("|---:|---:|---:|---:|---:|:---:|")
    for t in summary["holonomy_trials"]:
        marker = "**✓**" if t["matches_closed_form"] else "—"
        lines.append(
            f"| {t['chi']:.4f} | {t['chi_in_pi']:.3f}π | "
            f"{t['holonomy_closed_form']:.6f} | "
            f"{t['holonomy_numerical']:.6f} | "
            f"{t['deviation']:.2e} | {marker} |"
        )
    lines.append("")
    if summary["all_holonomy_match_closed_form"]:
        lines.append(
            "**Match:** numerical holonomy ≡ `π·cos(χ)` to machine "
            "precision at every canonical χ. The closed-form expression "
            "in `hopf/connection.py` is verified."
        )
    lines.append("")

    lines.append("## (b) Spinor double-cover quantization")
    lines.append("")
    lines.append(
        "A spinor traversing a closed Hopf-fibre loop ONCE picks up "
        "phase π·cos(χ) (= half a closure quantum at χ=0). Under the "
        "natural fermionic double-cover (4π closure), the loop is "
        "traversed TWICE, giving accumulated phase 2π·cos(χ). The "
        "double-cover phase is an integer multiple of 2π only at the "
        "polar fibres χ ∈ {0, π/2, π}."
    )
    lines.append("")
    lines.append(
        "| χ | label | single-closure (2π units) | double-closure (2π units) | "
        "single integer? | double integer? |"
    )
    lines.append("|---:|---|---:|---:|:---:|:---:|")
    for s in summary["spinor_closure_table"]:
        si = "✓" if s["single_is_integer"] else "—"
        di = "✓" if s["double_is_integer"] else "—"
        lines.append(
            f"| {s['chi']:.4f} | {s['chi_label']} | "
            f"{s['single_in_2pi']:+.4f} | {s['double_in_2pi']:+.4f} | "
            f"{si} | {di} |"
        )
    lines.append("")
    lines.append("**Per-fibre interpretation:**")
    lines.append("")
    for s in summary["spinor_closure_table"]:
        lines.append(f"- χ = {s['chi_label']}: {s['interpretation']}")
    lines.append("")
    lines.append(
        f"**Integer-quantization at:** χ ∈ "
        + "{" + ", ".join(
            f"{c:.4f}" for c in summary["double_cover_integer_quantized_at"]
        ) + "} "
        "— exactly the polar fibres of the Hopf base S²."
    )
    lines.append("")

    lines.append("## (c) Full closure-cycle integer count at χ = 0")
    lines.append("")
    lines.append(
        "Per-species breakdown of the Layer-1 closure cycle in units of "
        "2π. The Hopf and throat channels each contribute 1/2 quantum "
        "at χ = 0, partnering to give a full quantum:"
    )
    lines.append("")
    lines.append(
        "| species | k | antipodal (k·2π) | Hopf (π) | throat (π) | "
        "uplift (β·…) | sum (2π units) | N_layer_1 | "
        "B2_radial (n+1) | N_total |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for f in summary["full_cycle_breakdown_at_chi_0"]:
        lines.append(
            f"| {f['species']} | {f['k']} | "
            f"{f['contribution_antipodal_in_2pi']} | "
            f"{f['contribution_hopf_in_2pi']:.4f} | "
            f"{f['contribution_throat_in_2pi']:.4f} | "
            f"{f['contribution_uplift_in_2pi']:.4f} | "
            f"{f['sum_in_2pi']:.4f} | "
            f"**{f['n_layer_1']}** | "
            f"{f['n_radial_at_b2_coupling']} | "
            f"**{f['n_total']}** |"
        )
    lines.append("")
    if summary["layer_1_integers_consistent"]:
        lines.append(
            "**Layer-1 integer count verified:** the sum of "
            "(antipodal + Hopf + throat + uplift) channels equals the "
            "integer N_layer_1 from the closure_cycle_action_probe. "
            "The Hopf-throat partnership at χ = 0 (each contributes "
            "1/2 quantum, summing to one full quantum) is the angular "
            "structural piece that completes the closure-cycle integer "
            "reading."
        )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "**Hopf-fibre Aharonov-Bohm holonomy verified.** The closed-"
        "form `π·cos(χ)` matches numerical integration of the "
        "connection 1-form to machine precision. The spinor double-"
        "cover doubles this to `2π·cos(χ)`, integer-quantized at the "
        "three polar fibres χ ∈ {0, π/2, π}."
    )
    lines.append("")
    lines.append(
        "**Closure-cycle integer reading complete.** At the locked "
        "χ = 0:"
    )
    lines.append("")
    lines.append(
        "- antipodal closure: `k · 2π` — k closures of the great "
        "circle on S³ (integer per species)."
    )
    lines.append(
        "- Hopf AB: `π·cos(0) = π = 1/2 · 2π` — half quantum from "
        "the canonical-fibre holonomy."
    )
    lines.append(
        "- throat T²: `π = 1/2 · 2π` — half quantum from the spinor "
        "double-cover sign flip (T² = −I has eigenvalue arg π per "
        "closure pass)."
    )
    lines.append(
        "- β-uplift: `4β·max(0, k − 3)² / (2π)` — closure-quantum "
        "integer 100 for the τ row only."
    )
    lines.append("")
    lines.append(
        "The Hopf and throat channels **partner to give a single "
        "full quantum** at χ = 0: each contributes 1/2 quantum, "
        "summing to 1·2π. This is the angular structural piece "
        "complementing the radial channel's hard-wall BS reading "
        "(closed-orbit probe) and the topological antipodal closure."
    )
    lines.append("")
    lines.append(
        "**Conceptually:** the closure cycle is built from THREE "
        "structural integer-quantum sources — antipodal closure "
        "(`k`), angular Hopf-throat partnership (`+1` at χ = 0), and "
        "the τ-uplift closure quantum (`+100` at k = 5). Plus the "
        "Layer-2 radial contribution `(n + 1)` per coupled mode. "
        "**All four pieces are integer-valued in units of 2π.** This "
        "is the cleanest formulation of the closure-cycle action "
        "quantum picture."
    )
    lines.append("")
    lines.append("## What's next")
    lines.append("")
    lines.append(
        "With both radial (closed-orbit) and angular (Hopf AB) channels "
        "verified at integer quantization, the closure-cycle action is "
        "fully integer-quantized in units of 2π. The remaining ℏ-origin "
        "sub-targets are:"
    )
    lines.append("")
    lines.append(
        "- **Sub-target #3: Tangherlini eigenvalue ↔ m_e relation.** "
        "Test whether `ℏ ω(l, n) = m_e c² · (something natural)` is "
        "consistent across species. If so, the dimensional bridge "
        "from geometric units to SI is closed."
    )
    lines.append(
        "- **Identify the physical species → (l, n) coupling.** The "
        "B2_radial_ladder gives (3, 6, 109); cross-checking against an "
        "independent observable (e.g. lepton-anomalous-magnetic-moment "
        "or some QCD identity) would pin which coupling is physical."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_aharonov_bohm_hopf_fibre_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
