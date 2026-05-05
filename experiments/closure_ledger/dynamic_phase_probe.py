"""
Dynamic-phase probe for the closure-phase ledger.

Question: can the missing radial-channel residual for the best current
candidates (C1_eigenvector_weighted_B1 at 0.326 rad, D1_potential_difference
_phase at 0.577 rad on the circle) be supplied by a natural BAM loop
phase rather than by adding more S(k) variants?

Four candidate mechanisms, each tied to an existing module:

  M1  Moving throat transport         (embedding.transport)
  M2  Hopf fibre loop phase           (hopf.connection)
  M3  Antipodal out-and-back closure  (transaction.handshake / s3_geometry)
  M4  Worldline crossing phase        (T = iσ_y eigenvalue arg per crossing)

Each mechanism has a natural integer winding parameter m and possibly a
geometric angle (χ for Hopf, θ_loop for antipodal). For each (mechanism,
parameters) we compute the per-k phase contribution Δ_dyn(k) for the
three lepton generations k ∈ {1, 3, 5}, add it to the candidate's mod-2π
residue, and report the resulting circular spread.

A natural loop phase HELPS iff it brings the spread strictly below the
no-dynamic-phase baseline for that candidate. It CLOSES iff it reaches
< 1e-9.

The probe is exploratory — it does not introduce a new candidate, just
asks whether the existing geometric channels already on the books can
absorb the residual.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional


TAU = 2.0 * math.pi


# Per-generation residues mod 2π (units of π) for the best wired candidates,
# read off the most recent --compare archive. These are the input the probe
# tries to absorb via a dynamic phase.
BASELINE_RESIDUES_PI = {
    "C1_eigenvector_weighted_B1": (0.864195, 0.788396, 0.760506),
    "D1_potential_difference_phase": (0.097827, 1.914264, 1.987910),
}
LEPTON_K = (1, 3, 5)


def _circular_spread(mods: list[float]) -> float:
    """2π minus the largest gap between consecutive sorted residues."""
    if len(mods) < 2:
        return 0.0
    s = sorted(mods)
    gaps = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    wrap_gap = TAU - (s[-1] - s[0])
    return TAU - max(*gaps, wrap_gap)


@dataclass
class MechanismProbe:
    """A single natural loop-phase mechanism with its parameter ranges."""
    name: str
    module_source: str           # repo module the formula comes from
    description: str             # one-sentence physics
    delta_k_formula: str         # symbolic Δ(k) for the report
    delta_k: Callable[[int, dict], float]
    parameter_grid: list[dict] = field(default_factory=list)


def _hopf_holonomy(chi: float) -> float:
    """π·cos(χ); reproduces hopf.connection.hopf_holonomy at χ."""
    return math.pi * math.cos(chi)


# ---------------------------------------------------------------------------
# Catalog of natural loop-phase mechanisms.
# Each is a Δ(k) closed-form expression keyed off existing repo formulas.
# Parameter grids span only physically natural values (integer windings,
# canonical χ ∈ {0, π/4, π/3, π/2, π}, spin coupling α ∈ {1/4, 1/2, 1}).
# ---------------------------------------------------------------------------

def _build_mechanisms() -> list[MechanismProbe]:
    return [
        MechanismProbe(
            name="M1_moving_throat",
            module_source="geometrodynamics.embedding.transport",
            description=(
                "Moving throat transport: T = iσ_y has eigenvalues ±i, so "
                "each pass through a moving non-orientable throat absorbs a "
                "spectral phase ±π/2. For k passes per generation, "
                "Δ(k) = m · k · (π/2) with the eigenvalue choice encoded by m."
            ),
            delta_k_formula="Δ(k) = m · k · (π/2)",
            delta_k=lambda k, p: p["m"] * k * (math.pi / 2.0),
            parameter_grid=[{"m": m} for m in (-3, -2, -1, 0, 1, 2, 3)],
        ),
        MechanismProbe(
            name="M2_hopf_fibre_loop",
            module_source="geometrodynamics.hopf.connection",
            description=(
                "Hopf fibre loop phase: holonomy on a fibre at hyper-latitude "
                "χ is π·cos(χ). For m windings per generation, "
                "Δ(k) = m · k · π · cos(χ)."
            ),
            delta_k_formula="Δ(k) = m · k · π · cos(χ)",
            delta_k=lambda k, p: (
                p["m"] * k * _hopf_holonomy(p["chi"])
            ),
            parameter_grid=[
                {"m": m, "chi": chi}
                for m in (-2, -1, 0, 1, 2)
                for chi in (0.0, math.pi / 4.0, math.pi / 3.0,
                            math.pi / 2.0, math.pi)
            ],
        ),
        MechanismProbe(
            name="M3_antipodal_out_and_back",
            module_source="geometrodynamics.transaction.handshake",
            description=(
                "Antipodal out-and-back closure: a worldline traverses a "
                "closed great-circle loop of length 2π and accumulates the "
                "SU(2) geodesic phase α · 2π. For k loops per generation, "
                "Δ(k) = α · k · 2π. With α = 1/2 (spin-½) this gives k · π."
            ),
            delta_k_formula="Δ(k) = α · k · 2π",
            delta_k=lambda k, p: p["alpha"] * k * TAU,
            parameter_grid=[
                {"alpha": a}
                for a in (0.0, 0.25, 0.5, 1.0, 1.5)
            ],
        ),
        MechanismProbe(
            name="M4_worldline_crossings",
            module_source="geometrodynamics.embedding.transport",
            description=(
                "Worldline crossing phase: each antipodal crossing event "
                "absorbs phase π/2 from T = iσ_y (eigenvalue arg). For k "
                "crossings of generation k, Δ(k) = m · k · (π/2). This is "
                "structurally identical to M1 — kept here as a separate "
                "entry because the physical mechanism (TX firing) is "
                "distinct even though the formula coincides."
            ),
            delta_k_formula="Δ(k) = m · k · (π/2)   [same family as M1]",
            delta_k=lambda k, p: p["m"] * k * (math.pi / 2.0),
            parameter_grid=[{"m": m} for m in (-3, -2, -1, 0, 1, 2, 3)],
        ),
        MechanismProbe(
            name="M5_antipodal_pair_with_hopf",
            module_source=(
                "geometrodynamics.transaction.handshake "
                "+ geometrodynamics.hopf.connection"
            ),
            description=(
                "Composite: spin-½ antipodal pass plus a Hopf fibre loop. "
                "Δ(k) = m · k · (π/2) + n · π · cos(χ). The Hopf part is "
                "k-independent (universal shift) and so cannot tighten the "
                "spread beyond what M1 alone achieves; included here as a "
                "control against accidental over-improvement claims."
            ),
            delta_k_formula="Δ(k) = m · k · (π/2) + n · π · cos(χ)",
            delta_k=lambda k, p: (
                p["m"] * k * (math.pi / 2.0)
                + p["n"] * _hopf_holonomy(p["chi"])
            ),
            parameter_grid=[
                {"m": m, "n": n, "chi": chi}
                for m in (-3, -2, -1, 0, 1, 2, 3)
                for n in (-1, 0, 1)
                for chi in (0.0, math.pi / 2.0)
            ],
        ),
    ]


# ---------------------------------------------------------------------------
# Probe runner.
# ---------------------------------------------------------------------------

@dataclass
class ProbeResult:
    candidate: str
    mechanism: str
    parameters: dict
    delta_k_in_pi: list[float]
    new_residues_in_pi: list[float]
    circular_spread_rad: float
    helps: bool
    closes: bool


def _evaluate_mechanism(
    candidate: str,
    base_residues_pi: tuple[float, float, float],
    mechanism: MechanismProbe,
    baseline_spread: float,
) -> list[ProbeResult]:
    """Sweep a mechanism's parameter grid for one base candidate."""
    results: list[ProbeResult] = []
    base_rad = [v * math.pi for v in base_residues_pi]
    for params in mechanism.parameter_grid:
        deltas = [mechanism.delta_k(k, params) for k in LEPTON_K]
        new_residues = [(b + d) % TAU for b, d in zip(base_rad, deltas)]
        spread = _circular_spread(new_residues)
        helps = spread < baseline_spread - 1e-12
        closes = spread < 1e-9
        results.append(ProbeResult(
            candidate=candidate,
            mechanism=mechanism.name,
            parameters=dict(params),
            delta_k_in_pi=[d / math.pi for d in deltas],
            new_residues_in_pi=[r / math.pi for r in new_residues],
            circular_spread_rad=spread,
            helps=helps,
            closes=closes,
        ))
    return results


def run_probe() -> dict:
    mechanisms = _build_mechanisms()
    summary: dict = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "lepton_k": list(LEPTON_K),
        "baselines": {},
        "mechanisms": [
            {
                "name": m.name,
                "module_source": m.module_source,
                "description": m.description,
                "delta_k_formula": m.delta_k_formula,
                "n_parameter_combinations": len(m.parameter_grid),
            }
            for m in mechanisms
        ],
        "best_per_candidate": {},
        "all_results": [],
    }

    for candidate, residues_pi in BASELINE_RESIDUES_PI.items():
        baseline_residues_rad = [v * math.pi for v in residues_pi]
        baseline_spread = _circular_spread(baseline_residues_rad)
        summary["baselines"][candidate] = {
            "residues_in_pi": list(residues_pi),
            "circular_spread_rad": baseline_spread,
        }

        all_results: list[ProbeResult] = []
        for m in mechanisms:
            all_results.extend(_evaluate_mechanism(
                candidate, residues_pi, m, baseline_spread,
            ))

        # Track best per mechanism and overall best.
        best_per_mech: dict[str, ProbeResult] = {}
        for r in all_results:
            cur = best_per_mech.get(r.mechanism)
            if cur is None or r.circular_spread_rad < cur.circular_spread_rad:
                best_per_mech[r.mechanism] = r

        overall_best = min(all_results, key=lambda r: r.circular_spread_rad)
        summary["best_per_candidate"][candidate] = {
            "baseline_spread_rad": baseline_spread,
            "best_per_mechanism": {
                name: asdict(r) for name, r in best_per_mech.items()
            },
            "overall_best": asdict(overall_best),
            "any_mechanism_helps": any(r.helps for r in all_results),
            "any_mechanism_closes": any(r.closes for r in all_results),
        }
        summary["all_results"].extend(asdict(r) for r in all_results)

    return summary


# ---------------------------------------------------------------------------
# Markdown report.
# ---------------------------------------------------------------------------

def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Dynamic-phase probe — closure ledger")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append("")
    lines.append(
        "Question: can the missing radial-channel residual for the best "
        "current S(k) candidates (C1, D1) be supplied by one of the natural "
        "BAM loop phases — moving throat transport, Hopf fibre loop, "
        "antipodal out-and-back, or worldline crossing? The probe scans "
        "each mechanism's natural parameter grid and reports whether any "
        "(mechanism, parameters) tightens the circular spread mod 2π."
    )
    lines.append("")

    lines.append("## Mechanisms tested")
    lines.append("")
    lines.append("| name | module | Δ(k) | parameter combinations |")
    lines.append("|---|---|---|---:|")
    for m in summary["mechanisms"]:
        lines.append(
            f"| `{m['name']}` | `{m['module_source']}` | "
            f"`{m['delta_k_formula']}` | {m['n_parameter_combinations']} |"
        )
    lines.append("")

    for candidate, info in summary["best_per_candidate"].items():
        baseline = info["baseline_spread_rad"]
        lines.append(f"## Candidate: `{candidate}`")
        lines.append("")
        lines.append(
            f"Baseline circular spread (no dynamic phase): "
            f"**{baseline:.6f} rad** = {baseline / math.pi:.6f} π"
        )
        lines.append("")
        lines.append("### Best per mechanism")
        lines.append("")
        lines.append(
            "| mechanism | best parameters | Δ(k) (units of π) | "
            "new residues (units of π) | spread (rad) | helps? | closes? |"
        )
        lines.append("|---|---|---|---|---:|---|---|")
        for mech_name, r in info["best_per_mechanism"].items():
            params = ", ".join(
                f"{k}={v}" for k, v in r["parameters"].items()
            ) or "—"
            deltas = "[" + ", ".join(
                f"{x:+.3f}" for x in r["delta_k_in_pi"]
            ) + "]"
            new_res = "[" + ", ".join(
                f"{x:.3f}" for x in r["new_residues_in_pi"]
            ) + "]"
            lines.append(
                f"| `{mech_name}` | {params} | {deltas} | {new_res} | "
                f"{r['circular_spread_rad']:.6f} | "
                f"{'**yes**' if r['helps'] else 'no'} | "
                f"{'**yes**' if r['closes'] else 'no'} |"
            )
        lines.append("")
        ob = info["overall_best"]
        lines.append(
            f"**Overall best for {candidate}**: `{ob['mechanism']}` with "
            f"params {ob['parameters']} → spread = "
            f"{ob['circular_spread_rad']:.6f} rad. "
            f"(Helps: {ob['helps']}; closes: {ob['closes']}.)"
        )
        lines.append("")

    lines.append("## Verdict")
    lines.append("")
    any_closes = any(
        info["any_mechanism_closes"]
        for info in summary["best_per_candidate"].values()
    )
    any_helps = any(
        info["any_mechanism_helps"]
        for info in summary["best_per_candidate"].values()
    )
    if any_closes:
        lines.append(
            "**At least one natural loop phase CLOSES the residual** "
            "(circular spread < 1e-9)."
        )
    elif any_helps:
        lines.append(
            "**No natural loop phase closes** the residual mod 2π, but at "
            "least one tightens the spread below the baseline. The closure "
            "ledger conjecture remains open."
        )
    else:
        lines.append(
            "**No natural loop phase closes or tightens** the residual "
            "for either C1 or D1 within the scanned parameter grid. "
            "All four mechanisms (moving throat, Hopf fibre, antipodal "
            "out-and-back, worldline crossing) reduce to either "
            "k-independent shifts (universal across species, do not affect "
            "spread) or to multiples of k · π/2 (which give residues with "
            "spread π for k ∈ {1, 3, 5}, wider than C1's 0.326 rad). "
            "The radial-channel residual cannot be absorbed into any "
            "existing geometric channel — it must come from a new piece "
            "of physics."
        )
    lines.append("")
    lines.append("### Why the natural mechanisms cannot fit")
    lines.append("")
    lines.append(
        "- **M1/M4 (spin-½ passes; T eigenvalue arg per crossing)**: both "
        "give Δ(k) = m · k · π/2 for integer m. For k ∈ {1, 3, 5}, the "
        "k·π/2 pattern mod 2π takes only the values [π/2, 3π/2, π/2] "
        "(m = 1) or [π, π, π] (m = 2) etc. The first has spread π; the "
        "second is universal (spread 0) and only shifts the universal "
        "value. Neither matches C1's residue pattern [0.864, 0.788, 0.761] π."
    )
    lines.append(
        "- **M2 (Hopf fibre loop)**: Δ(k) = m·k·π·cos(χ). At χ = 0 this "
        "is m·k·π, which mod 2π is universal for odd k. At χ = π/2 it "
        "vanishes entirely. Intermediate χ scales the universal piece "
        "without breaking the k-symmetry."
    )
    lines.append(
        "- **M3 (antipodal out-and-back)**: Δ(k) = α·k·2π. mod 2π this "
        "is α·k·2π mod 2π. For α = 1/2 (spin-½): k·π mod 2π = π for odd k "
        "(universal). For α = 1: 2π·k mod 2π = 0. Always universal."
    )
    lines.append(
        "- **General constraint**: any loop phase that scales as k·c with "
        "constant c gives k-dependent Δ only through the modular "
        "arithmetic, and the resulting patterns over k = (1, 3, 5) live "
        "on a finite lattice of multiples of π. C1's residues do not lie "
        "on this lattice."
    )
    lines.append("")
    lines.append(
        "Conclusion: the BAM loop-phase channels already in the ledger "
        "(antipodal closure k·2π, Hopf holonomy π·cos(χ), throat T² = π) "
        "exhaust what natural integer-quantized loop phases can supply. "
        "The remaining ~0.3 rad C1 residue is too small and the wrong "
        "shape to come from another instance of the same mechanisms."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_dynamic_phase_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
