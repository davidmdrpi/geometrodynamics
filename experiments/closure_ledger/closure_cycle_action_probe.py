"""
Closure-cycle action quantum probe.

The first concrete probe of the ℏ-origin research plan
(`docs/hbar_origin_research_plan.md`). Tests prediction P1:

  Φ_cycle  =  ∮(L / ℏ) dt  =  2π · N(species)

with N a small integer derivable from species-level inputs (k, color,
Z₂ partition). The Layer-1 ledger is already known to close mod 2π
universally for odd k (odd-k closure lemma). This probe asks the
sharper question: is each species' Layer-1 ledger sum already an
INTEGER multiple of 2π (not just zero mod 2π), and what is that
integer?

If yes: the closure cycle has a per-species action-quantum count
N(species) which IS the geometric content of the closure structure.
The conversion to physical ℏ requires the dimensional bridge
`ℏ_SI = (m_e R_MID c) / (2π)` per single closure quantum, which is
prediction P2 (the reduced-Compton-wavelength identification).

This probe also computes the leading-order numerical conversion for
the canonical R_MID = ℏ/(m_e c) identification and reports whether
it gives self-consistent ℏ in CGS units.

Layer 1 ledger sum at depth k under the locked baseline:

  Φ_avail(k)  =  k · 2π  +  π·cos(χ=0)  +  2·(π/2)  +  β · max(0, k−3)²
              =  2π · (k + 1)  +  50π · max(0, k − 3)²

Per species:

  electron (k=1):  Φ_avail = 2π · 2  +  0           = 4π     = 2 · 2π
  muon (k=3):      Φ_avail = 2π · 4  +  0           = 8π     = 4 · 2π
  tau (k=5):       Φ_avail = 2π · 6  +  50π · 4     = 212π   = 106 · 2π

So each species' Layer-1 ledger sum IS an integer multiple of 2π —
the values being 2, 4, 106 respectively. This is the leading-order
test of P1.

The probe also explores what Layer 2 contributions would preserve or
break the integer-quantum reading.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


TAU = 2.0 * math.pi
ACTION_BASE = TAU                       # 2π — S³ great-circle action
HOPF_HOLONOMY_AT_CHI_0 = math.pi        # π · cos(0)
THROAT_T2_PHASE = math.pi               # T² = -I → eigenvalue arg per closure
BETA_LEPTON = 50.0 * math.pi            # closure quantum: 4β = 100·(2π)
LEPTON_DEPTHS = (1, 3, 5)
LEPTON_LABELS = ("electron", "muon", "tau")

# CGS reference values for the ℏ identification check.
M_E_GRAMS = 9.1093837015e-28            # electron rest mass (g)
HBAR_CGS = 1.054571817e-27              # ℏ (erg·s)
C_CGS = 2.99792458e10                   # speed of light (cm/s)
LAMBDA_E_REDUCED_CM = HBAR_CGS / (M_E_GRAMS * C_CGS)  # reduced Compton wavelength


@dataclass
class SpeciesCycle:
    """Per-species Layer-1 closure-cycle action contribution."""
    label: str
    k: int
    antipodal: float        # k · action_base
    hopf: float             # Hopf holonomy at χ=0
    throat: float           # T² closure phase
    uplift: float           # β · max(0, k-3)²
    phi_total: float
    n_quanta: float         # phi_total / (2π); integer for clean quantization
    is_integer_quantized: bool
    n_quanta_rounded: int


def _layer_1_per_species() -> list[SpeciesCycle]:
    rows: list[SpeciesCycle] = []
    for k, label in zip(LEPTON_DEPTHS, LEPTON_LABELS):
        antipodal = k * ACTION_BASE
        hopf = HOPF_HOLONOMY_AT_CHI_0
        throat = THROAT_T2_PHASE
        uplift = BETA_LEPTON * max(0, k - 3) ** 2
        phi_total = antipodal + hopf + throat + uplift
        n_quanta = phi_total / TAU
        n_rounded = round(n_quanta)
        is_int = abs(n_quanta - n_rounded) < 1e-12
        rows.append(SpeciesCycle(
            label=label,
            k=k,
            antipodal=antipodal,
            hopf=hopf,
            throat=throat,
            uplift=uplift,
            phi_total=phi_total,
            n_quanta=n_quanta,
            is_integer_quantized=is_int,
            n_quanta_rounded=n_rounded,
        ))
    return rows


@dataclass
class HBarConversionCheck:
    """Test of the R_MID = ℏ/(m_e c) identification (prediction P2)."""
    species_label: str
    n_quanta: int
    R_MID_geom_units: float              # always 1.0 by construction
    R_MID_physical_cm: float             # = N · ℏ/(m_e c) under P2
    closure_cycle_geometric_length: float  # = N · 2π in geometric units
    closure_cycle_physical_cm: float
    closure_cycle_physical_compton_wavelengths: float
    self_consistent: bool


def _hbar_conversion_check(rows: list[SpeciesCycle]) -> list[HBarConversionCheck]:
    """For each species, compute the physical length of the closure cycle
    under the canonical R_MID = ℏ/(m_e c) identification."""
    out: list[HBarConversionCheck] = []
    for r in rows:
        N = r.n_quanta_rounded
        # Under P2: each "1 unit" of geometric length is the reduced
        # Compton wavelength of the electron.
        R_MID_cm = LAMBDA_E_REDUCED_CM
        cycle_geom = N * TAU   # phi_cycle / 1 (taking geometric units = 1)
        cycle_cm = cycle_geom * R_MID_cm
        cycle_in_compton = cycle_cm / LAMBDA_E_REDUCED_CM
        out.append(HBarConversionCheck(
            species_label=r.label,
            n_quanta=N,
            R_MID_geom_units=1.0,
            R_MID_physical_cm=R_MID_cm,
            closure_cycle_geometric_length=cycle_geom,
            closure_cycle_physical_cm=cycle_cm,
            closure_cycle_physical_compton_wavelengths=cycle_in_compton,
            self_consistent=abs(cycle_in_compton - N * TAU) < 1e-12,
        ))
    return out


@dataclass
class Layer2EffectAnalysis:
    """How does adding a Layer-2 contribution affect integer quantization?"""
    layer_2_candidate: str
    description: str
    delta_phi_in_pi: list[float]        # per species, in units of π
    new_phi_in_2pi: list[float]         # per species, in units of 2π
    preserves_integer_quantization: bool
    new_residue_circular_spread_rad: float


def _circular_spread_rad(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    s = sorted(values)
    gaps = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    wrap_gap = TAU - (s[-1] - s[0])
    return TAU - max(*gaps, wrap_gap)


def _layer_2_effect_analysis(rows: list[SpeciesCycle]) -> list[Layer2EffectAnalysis]:
    """
    For each candidate Layer-2 contribution from prior probes, compute
    whether adding it preserves or breaks the integer-quantum reading.
    """
    # Per-species Layer-2 candidate values, in units of π.
    # These are the residues observed in the closure-ledger experiment.
    candidates = [
        (
            "Layer 2 absent (ledger universal at 0 mod 2π)",
            "Layer 1 only — the odd-k lemma's domain of validity.",
            [0.0, 0.0, 0.0],   # no contribution
        ),
        (
            "C1 (eigenvector-weighted B1 modes)",
            "Best scalar S(k) candidate from the closure-ledger sweep; "
            "circular spread 0.326 rad.",
            [0.864195, 0.788396, 0.760506],
        ),
        (
            "D1 (operator-valued V_j-V_i Hermitian matrix element)",
            "Best operator-valued candidate from the D-family probe; "
            "circular spread 0.577 rad.",
            [0.097827, 1.914264, 1.987910],
        ),
    ]
    out: list[Layer2EffectAnalysis] = []
    for name, desc, dphi_pi in candidates:
        new_phis_2pi = []
        new_residues = []
        for r, d_pi in zip(rows, dphi_pi):
            new_phi_2pi = (r.phi_total + d_pi * math.pi) / TAU
            new_phis_2pi.append(new_phi_2pi)
            # Reduce to mod 1 (i.e. mod 2π in phi units) for spread.
            new_residues.append((new_phi_2pi % 1.0) * TAU)
        # Preserves integer quantization iff every new_phi_2pi is integer.
        is_int = all(abs(p - round(p)) < 1e-12 for p in new_phis_2pi)
        spread = _circular_spread_rad(new_residues)
        out.append(Layer2EffectAnalysis(
            layer_2_candidate=name,
            description=desc,
            delta_phi_in_pi=list(dphi_pi),
            new_phi_in_2pi=new_phis_2pi,
            preserves_integer_quantization=is_int,
            new_residue_circular_spread_rad=spread,
        ))
    return out


# ---------------------------------------------------------------------------

def run_probe() -> dict:
    rows = _layer_1_per_species()
    conversions = _hbar_conversion_check(rows)
    layer2_effects = _layer_2_effect_analysis(rows)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "constants": {
            "action_base_2pi": ACTION_BASE,
            "hopf_at_chi_0": HOPF_HOLONOMY_AT_CHI_0,
            "throat_T2": THROAT_T2_PHASE,
            "beta_lepton_50pi": BETA_LEPTON,
            "M_e_grams": M_E_GRAMS,
            "hbar_cgs": HBAR_CGS,
            "c_cgs": C_CGS,
            "lambda_e_reduced_cm": LAMBDA_E_REDUCED_CM,
        },
        "layer_1_per_species": [asdict(r) for r in rows],
        "all_species_integer_quantized": all(r.is_integer_quantized for r in rows),
        "p1_status": (
            "PASS" if all(r.is_integer_quantized for r in rows)
            else "FAIL"
        ),
        "n_quanta_per_species": {
            r.label: r.n_quanta_rounded for r in rows
        },
        "hbar_conversion_check": [asdict(c) for c in conversions],
        "layer_2_effect_analysis": [asdict(e) for e in layer2_effects],
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Closure-cycle action quantum probe")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append("")
    lines.append(
        "First concrete probe of the ℏ-origin research plan "
        "(`docs/hbar_origin_research_plan.md`). Tests prediction P1: "
        "the per-species Layer-1 closure-cycle phase is an integer "
        "multiple of 2π, and that integer N is the geometric content "
        "of the closure structure."
    )
    lines.append("")

    lines.append("## Layer-1 ledger contributions per species")
    lines.append("")
    lines.append(
        "All values in units of π. `Φ_total / 2π` is the integer-"
        "quantum count if Layer 1 is exactly quantized."
    )
    lines.append("")
    lines.append(
        "| species | k | antipodal | Hopf | throat | uplift | "
        "Φ_total | Φ_total / 2π | integer? |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|:---:|")
    for r in summary["layer_1_per_species"]:
        f = lambda v: f"{v / math.pi:.1f}π"
        marker = "**✓**" if r["is_integer_quantized"] else "—"
        lines.append(
            f"| {r['label']} | {r['k']} | {f(r['antipodal'])} | "
            f"{f(r['hopf'])} | {f(r['throat'])} | {f(r['uplift'])} | "
            f"{f(r['phi_total'])} | {r['n_quanta']:.6f} | {marker} |"
        )
    lines.append("")
    p1 = summary["p1_status"]
    n_quanta = summary["n_quanta_per_species"]
    lines.append(
        f"**Prediction P1:** {p1}. Per-species integer quantum counts: "
        + ", ".join(f"`{lbl}: N = {n}`" for lbl, n in n_quanta.items())
        + "."
    )
    lines.append("")
    if p1 == "PASS":
        lines.append(
            "Each species' Layer-1 ledger sum is an integer multiple of "
            "`2π` (the action quantum candidate). The N counts grow as "
            "`{k+1, k+1, k+1+100·(k=5)} = {2, 4, 106}` — the first two "
            "are simply (k + 1) closure passes (electron and muon have "
            "no β·max(0, k-3)² contribution), and the τ row picks up "
            "the closure-quantum integer 100 from the 4β = 100·(2π) "
            "lock."
        )
    lines.append("")

    lines.append("## R_MID = ℏ/(m_e c) identification (prediction P2)")
    lines.append("")
    lines.append(
        f"Reduced Compton wavelength of the electron: "
        f"λ_e_reduced = ℏ/(m_e c) = "
        f"{summary['constants']['lambda_e_reduced_cm']:.6e} cm."
    )
    lines.append("")
    lines.append(
        "Under P2, BAM's geometric `R_MID = 1` corresponds to "
        "λ_e_reduced. The closure cycle in physical units:"
    )
    lines.append("")
    lines.append(
        "| species | N quanta | cycle (geom units) | cycle (cm) | "
        "cycle (Compton wavelengths) |"
    )
    lines.append("|---|---:|---:|---:|---:|")
    for c in summary["hbar_conversion_check"]:
        lines.append(
            f"| {c['species_label']} | {c['n_quanta']} | "
            f"{c['closure_cycle_geometric_length']:.4f} | "
            f"{c['closure_cycle_physical_cm']:.4e} | "
            f"{c['closure_cycle_physical_compton_wavelengths']:.4f} |"
        )
    lines.append("")
    lines.append(
        "By construction, the cycle length in Compton wavelengths "
        "equals N · 2π (since R_MID = λ_e_reduced is the conversion). "
        "This is **self-consistent** — but it does not yet PREDICT "
        "ℏ; it just identifies that under the canonical R_MID "
        "convention, the closure cycle traverses N · 2π reduced "
        "Compton wavelengths. Predicting ℏ requires R_MID to be "
        "geometrically determined (THESIS.md \"Self-consistent throat "
        "radius\" target), which is open."
    )
    lines.append("")

    lines.append("## Layer-2 effects on integer quantization")
    lines.append("")
    lines.append(
        "Each row tests whether adding a candidate Layer-2 contribution "
        "preserves integer-quantum closure. (None of the prior "
        "probes' candidates did — but the table records by how much "
        "they break it.)"
    )
    lines.append("")
    lines.append(
        "| Layer 2 candidate | Δφ per species (π) | new φ/2π per species | "
        "preserves integer? | new spread (rad) |"
    )
    lines.append("|---|---|---|:---:|---:|")
    for e in summary["layer_2_effect_analysis"]:
        dphi = "[" + ", ".join(f"{x:.3f}" for x in e["delta_phi_in_pi"]) + "]"
        nphi = "[" + ", ".join(f"{x:.3f}" for x in e["new_phi_in_2pi"]) + "]"
        marker = "✓" if e["preserves_integer_quantization"] else "—"
        lines.append(
            f"| {e['layer_2_candidate']} | {dphi} | {nphi} | {marker} | "
            f"{e['new_residue_circular_spread_rad']:.4f} |"
        )
    lines.append("")
    lines.append(
        "Layer 1 alone preserves integer quantization (P1 PASS). C1 "
        "and D1 — the best Layer-2 scalar/operator candidates from "
        "the closure-ledger sweep — both BREAK integer quantization "
        "(they shift the per-species N by non-integer amounts). This "
        "is consistent with their FAIL on Layer-2 universality: a "
        "candidate that broke the per-species integer count but still "
        "made the residue universal would imply a NEW action quantum "
        "different from `2π`, which the lemma's structure does not "
        "allow."
    )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if p1 == "PASS":
        lines.append(
            "**P1 PASS.** Each species' Layer-1 ledger sum is exactly "
            "`N · 2π` for integer N ∈ {2, 4, 106}. The geometric "
            "action quantum candidate is `2π` per closure unit; "
            "summing the four wired channels in the locked baseline "
            "lands every species on an integer multiple, with no "
            "fractional residue."
        )
        lines.append("")
        lines.append(
            "**P2 self-consistent, not yet predictive.** Identifying "
            "BAM's geometric R_MID = 1 with the reduced Compton "
            "wavelength of the electron (λ_e_reduced ≈ 3.86 × 10⁻¹¹ "
            "cm) makes the conversion `ℏ_SI = m_e R_MID c` a tautology. "
            "Predicting ℏ requires R_MID to be geometrically determined "
            "from a self-consistency condition (e.g. equilibrium "
            "throat radius for the locked mass spectrum), which is "
            "open."
        )
        lines.append("")
        lines.append(
            "**Layer 2 status.** No Layer-2 candidate from the closure-"
            "ledger sweep preserves integer quantization at the per-"
            "species level. This refines the earlier Layer-2 verdict: "
            "any Layer-2 closure form must contribute exactly an "
            "integer multiple of `2π` per species (NOT just be "
            "universal mod 2π) to preserve the action-quantum reading. "
            "The current candidates do neither."
        )
    else:
        lines.append(
            "**P1 FAIL.** Layer-1 ledger does not lock to integer "
            "multiples of `2π` per species under the locked baseline. "
            "This contradicts the odd-k closure lemma's universality "
            "claim and would require re-examining the locked constants."
        )
    lines.append("")
    lines.append("## What's next")
    lines.append("")
    lines.append(
        "The ℏ-origin research plan identifies four sub-targets in "
        "decreasing tractability. P1 PASS narrows the next sub-probe "
        "scope:"
    )
    lines.append("")
    lines.append(
        "- **(1) Layer 2 closure as integer-quantum contribution.** "
        "The next probe should test whether a worldline-integral "
        "Layer-2 form (rather than the per-mode WKB sum used by C1/D1) "
        "yields integer-multiple-of-2π per species. The natural "
        "candidate is the full `∮p·dq` over the closed orbit on the "
        "Tangherlini grid, with the closure cycle treated as a single "
        "integral rather than a sum of channels."
    )
    lines.append(
        "- **(2) Aharonov-Bohm form.** Compute the AB phase explicitly "
        "for a Hopf-fibre loop and verify it equals `2π` per spinor "
        "double-cover closure. This is the conceptually-clean "
        "geometric realization of the action quantum and should be "
        "verifiable in 50 lines of code given `geometrodynamics/hopf/"
        "connection.py`."
    )
    lines.append(
        "- **(3) Tangherlini eigenvalue ω_0 = m_e c² / ℏ?** Test "
        "whether `ℏ ω_0 = m_e c²` is consistent across species. If "
        "the lowest Tangherlini eigenvalue equals the electron rest "
        "energy in geometric units, that's an ℏ relationship."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_closure_cycle_action_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
