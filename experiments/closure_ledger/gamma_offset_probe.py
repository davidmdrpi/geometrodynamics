"""
γ-offset probe.

Closes in on the residual ~2.2% gap between the bare radial geometric
quantity Σ_{l=1..5} V_max(l) ≈ 22.008 and the locked γ_lepton = 22.5.
Probes principled augmentations of the QCD-style pinhole formula:

  (1) Extending the angular sum to include l = 0 (the 5D-specific
      centrifugal-free 3·rs²/r⁴ barrier).
  (2) Alternative angular weightings (l → 2l+1, l → l(l+1), etc.)
  (3) Higher-cutoff sums (l up to 6, 7).
  (4) V evaluated at the classical turning point of each ω(l, 0).
  (5) Symmetric combinations of inner and outer barrier maxima.

For each candidate the probe reports the value, %Δ vs γ_lepton and
γ_quark, and the muon-mass error obtained when the locked block is
re-evaluated with that γ. A candidate "explains" the offset if it
brings γ within ≤ 1% of the lock AND reduces the muon-mass error to
< 5%.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


GAMMA_LEPTON = 22.5
GAMMA_QUARK = 22.25


@dataclass
class OffsetCandidate:
    name: str
    formula: str
    value: float
    pct_diff_lepton: float
    pct_diff_quark: float
    muon_mass_err_pct: float
    tau_mass_err_pct: float
    closes_gamma_within_1pct: bool
    closes_muon_within_5pct: bool


def _build_grid():
    """Canonical Chebyshev tortoise grid matching `solve_radial_modes` defaults."""
    import numpy as np
    from geometrodynamics.tangherlini.radial import r_to_rstar, rstar_to_r
    from geometrodynamics.constants import R_MID, R_OUTER
    rs = float(R_MID)
    N = 80
    rs_min = r_to_rstar(rs + 5e-4, rs)
    rs_max = r_to_rstar(R_OUTER - 5e-4, rs)
    x = np.cos(math.pi * np.arange(N + 1) / N)
    L = (rs_max - rs_min) / 2.0
    rsg = rs_min + L * (1.0 - x)
    rg = np.array([rstar_to_r(s, rs) for s in rsg])
    return rg, rs


def _v_max(l: int, rg, rs) -> float:
    import numpy as np
    from geometrodynamics.tangherlini.radial import V_tangherlini
    return float(np.max(V_tangherlini(rg, l, rs)))


def _v_at_turning_point(l: int, omega: float, rs: float) -> float:
    """V at the outer classical turning point of the (l, n=0) ground state."""
    import numpy as np
    from geometrodynamics.tangherlini.radial import V_tangherlini
    from geometrodynamics.constants import R_OUTER
    rgrid = np.linspace(rs * 1.001, R_OUTER * 0.999, 100000)
    Vg = V_tangherlini(rgrid, l, rs)
    diff = Vg - omega ** 2
    sc = np.where(np.diff(np.sign(diff)) != 0)[0]
    if len(sc) == 0:
        return float("nan")
    return float(Vg[sc[-1]])


def _build_candidates() -> list[tuple[str, str, float]]:
    """Return list of (name, formula_str, value) tuples."""
    from geometrodynamics.tangherlini.radial import solve_radial_modes
    rg, rs = _build_grid()

    v_max_per_l = {l: _v_max(l, rg, rs) for l in range(0, 8)}
    omegas = {
        l: float(solve_radial_modes(l, N=80, n_modes=2)[0][0])
        for l in range(1, 6)
    }
    v_tp = {l: _v_at_turning_point(l, omegas[l], rs) for l in range(1, 6)}

    cands: list[tuple[str, str, float]] = []

    # Bare reference (the QCD pinhole formula).
    cands.append((
        "ref_QCD_Sum_l_1to5_V_max",
        "Σ_{l=1..5} V_max(l)",
        sum(v_max_per_l[l] for l in range(1, 6)),
    ))

    # Extension to l = 0 (the headline candidate).
    cands.append((
        "extend_Sum_l_0to5_V_max",
        "Σ_{l=0..5} V_max(l)   [adds 5D centrifugal-free l=0 barrier]",
        sum(v_max_per_l[l] for l in range(0, 6)),
    ))

    # Higher cutoffs.
    cands.append((
        "extend_Sum_l_1to6_V_max",
        "Σ_{l=1..6} V_max(l)",
        sum(v_max_per_l[l] for l in range(1, 7)),
    ))
    cands.append((
        "extend_Sum_l_0to6_V_max",
        "Σ_{l=0..6} V_max(l)",
        sum(v_max_per_l[l] for l in range(0, 7)),
    ))

    # Multiplicity-weighted sums (S^3 spherical-harmonic degeneracies).
    cands.append((
        "weight_Sum_2lp1_V_max_l_1to5",
        "Σ_{l=1..5} (2l+1) V_max(l)",
        sum((2 * l + 1) * v_max_per_l[l] for l in range(1, 6)),
    ))
    cands.append((
        "weight_Sum_lp1_V_max_l_0to5",
        "Σ_{l=0..5} (l+1) V_max(l)",
        sum((l + 1) * v_max_per_l[l] for l in range(0, 6)),
    ))
    cands.append((
        "weight_Sum_lpls2_V_max_l_1to5",
        "Σ_{l=1..5} l(l+2) V_max(l)   [S^3 angular eigenvalue weighting]",
        sum(l * (l + 2) * v_max_per_l[l] for l in range(1, 6)),
    ))

    # Turning-point evaluations.
    cands.append((
        "tp_Sum_l_1to5_V_at_outer_tp",
        "Σ_{l=1..5} V(r_tp_outer, l)",
        sum(v_tp[l] for l in range(1, 6)),
    ))

    # Composite: bare + half ω² average (a sub-leading correction guess).
    bare = sum(v_max_per_l[l] for l in range(1, 6))
    cands.append((
        "composite_bare_plus_min_omega_sq",
        "Σ_{l=1..5} V_max(l) + min_l ω(l, 0)²",
        bare + min(omegas[l] ** 2 for l in range(1, 6)),
    ))
    cands.append((
        "composite_bare_plus_V_max_l_eq_0",
        "Σ_{l=1..5} V_max(l) + V_max(l=0)",
        bare + v_max_per_l[0],
    ))

    # Mid-grid V_max (extended R_outer to capture analytic peaks).
    import numpy as np
    from geometrodynamics.tangherlini.radial import V_tangherlini
    rg_dense = np.linspace(rs * 1.001, 2.0 * rs, 200000)
    v_max_dense = {l: float(np.max(V_tangherlini(rg_dense, l, rs))) for l in range(0, 6)}
    cands.append((
        "ext_R_to_2rs_Sum_l_1to5",
        "Σ_{l=1..5} V_max(l) on r ∈ [rs, 2·rs]",
        sum(v_max_dense[l] for l in range(1, 6)),
    ))
    cands.append((
        "ext_R_to_2rs_Sum_l_0to5",
        "Σ_{l=0..5} V_max(l) on r ∈ [rs, 2·rs]",
        sum(v_max_dense[l] for l in range(0, 6)),
    ))

    return cands


def _mass_sensitivity_for_gamma(gamma: float) -> tuple[float, float]:
    """Plug γ into the locked block and return (m_μ_err_pct, m_τ_err_pct)."""
    import numpy as np
    from scipy.linalg import eigh
    from geometrodynamics.tangherlini.lepton_spectrum import (
        _build_generation_block,
        LEPTON_BASELINE_DEPTHS,
        LEPTON_BASELINE_PHASE,
        LEPTON_BASELINE_TRANSPORT,
        LEPTON_BASELINE_RESISTANCE,
        S3_ACTION_BASE,
        TAU_BETA_50PI,
    )
    M_E = 0.5109989461
    OBS = {3: 105.6583745, 5: 1776.86}
    h = _build_generation_block(
        depths=LEPTON_BASELINE_DEPTHS,
        phase_per_pass=LEPTON_BASELINE_PHASE,
        transport_strength=LEPTON_BASELINE_TRANSPORT,
        resistance_model="exponential",
        resistance_scale=LEPTON_BASELINE_RESISTANCE,
        hard_pinhole_gamma=gamma,
        action_base=S3_ACTION_BASE,
        action_slope=0.5,
        depth_cost_mode="tunnel_only",
        winding_mode="max",
        k_uplift_beta=TAU_BETA_50PI,
    )
    w = sorted(float(x) for x in eigh(h, eigvals_only=True) if x > 0)
    if len(w) < 3:
        return float("nan"), float("nan")
    scale = M_E / w[0]
    m_mu = w[1] * scale
    m_tau = w[2] * scale
    return (
        100.0 * abs(m_mu - OBS[3]) / OBS[3],
        100.0 * abs(m_tau - OBS[5]) / OBS[5],
    )


def run_probe() -> dict:
    raw = _build_candidates()
    out: list[OffsetCandidate] = []
    for name, formula, value in raw:
        pl = 100.0 * (value - GAMMA_LEPTON) / GAMMA_LEPTON
        pq = 100.0 * (value - GAMMA_QUARK) / GAMMA_QUARK
        m_mu_err, m_tau_err = _mass_sensitivity_for_gamma(value)
        out.append(OffsetCandidate(
            name=name,
            formula=formula,
            value=value,
            pct_diff_lepton=pl,
            pct_diff_quark=pq,
            muon_mass_err_pct=m_mu_err,
            tau_mass_err_pct=m_tau_err,
            closes_gamma_within_1pct=abs(pl) <= 1.0,
            closes_muon_within_5pct=m_mu_err < 5.0,
        ))

    # Best by γ-closeness; best by muon-match.
    best_gamma = min(out, key=lambda c: abs(c.pct_diff_lepton))
    best_muon = min(out, key=lambda c: c.muon_mass_err_pct)
    joint_winners = [
        c for c in out
        if c.closes_gamma_within_1pct and c.closes_muon_within_5pct
    ]

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "gamma_lepton_target": GAMMA_LEPTON,
        "gamma_quark_target": GAMMA_QUARK,
        "candidates": [asdict(c) for c in out],
        "best_by_gamma_closeness": asdict(best_gamma),
        "best_by_muon_match": asdict(best_muon),
        "joint_winners_within_1pct_gamma_and_5pct_muon": [
            asdict(c) for c in joint_winners
        ],
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# γ-offset probe — closing the residual ~2.2% pinhole gap")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append(
        f"**Targets:** γ_lepton = {summary['gamma_lepton_target']}, "
        f"γ_quark = {summary['gamma_quark_target']}"
    )
    lines.append("")
    lines.append(
        "Bare radial geometric pinhole Σ_{l=1..5} V_max(l) = 22.0082 "
        "sits 2.2% below the locked γ_lepton = 22.5. This probe asks "
        "whether a principled augmentation of the QCD-style formula "
        "closes the offset, and whether the resulting γ also reproduces "
        "the muon mass within a few percent (the strict empirical test)."
    )
    lines.append("")

    lines.append("## All candidates")
    lines.append("")
    lines.append(
        "| candidate | formula | value | %Δ vs γ_lep | %Δ vs γ_qrk | "
        "m_μ err | m_τ err |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for c in summary["candidates"]:
        lines.append(
            f"| `{c['name']}` | `{c['formula']}` | "
            f"{c['value']:.4f} | "
            f"{c['pct_diff_lepton']:+.3f}% | {c['pct_diff_quark']:+.3f}% | "
            f"{c['muon_mass_err_pct']:.3f}% | {c['tau_mass_err_pct']:.3f}% |"
        )
    lines.append("")

    bg = summary["best_by_gamma_closeness"]
    bm = summary["best_by_muon_match"]
    lines.append("## Best by γ-closeness")
    lines.append("")
    lines.append(
        f"`{bg['name']}` = {bg['value']:.4f} ({bg['pct_diff_lepton']:+.3f}% "
        f"vs γ_lepton). Formula: `{bg['formula']}`. "
        f"Muon mass error at this γ: **{bg['muon_mass_err_pct']:.3f}%**."
    )
    lines.append("")
    lines.append("## Best by muon-mass match")
    lines.append("")
    lines.append(
        f"`{bm['name']}` = {bm['value']:.4f}. Muon error: "
        f"**{bm['muon_mass_err_pct']:.3f}%**. "
        f"%Δ vs γ_lepton: {bm['pct_diff_lepton']:+.3f}%."
    )
    lines.append("")

    winners = summary["joint_winners_within_1pct_gamma_and_5pct_muon"]
    lines.append("## Joint winners (within 1% of γ AND within 5% on m_μ)")
    lines.append("")
    if not winners:
        lines.append("(none)")
    else:
        lines.append("| candidate | formula | value | %Δ γ_lep | m_μ err |")
        lines.append("|---|---|---:|---:|---:|")
        for c in winners:
            lines.append(
                f"| `{c['name']}` | `{c['formula']}` | {c['value']:.4f} | "
                f"{c['pct_diff_lepton']:+.3f}% | {c['muon_mass_err_pct']:.3f}% |"
            )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if winners:
        top = winners[0]
        lines.append(
            f"**The γ offset is explained by `{top['name']}`.** "
            f"Formula: `{top['formula']}` evaluates to {top['value']:.4f} "
            f"({top['pct_diff_lepton']:+.3f}% vs γ_lepton, "
            f"{top['muon_mass_err_pct']:.3f}% on the muon mass)."
        )
        lines.append("")
        lines.append(
            "**Implication.** The lepton pinhole is a Tangherlini "
            "barrier-spectrum sum on the canonical tortoise grid, with "
            "a definite l-range that includes one more channel than the "
            "QCD pinhole's `Σ_{l=1..5}` formula. The full lepton "
            "diagonal is now constructed from geometric ingredients "
            "alone (closure quantum 4β = 100·(2π) at τ; barrier sum "
            f"`{top['formula']}` at the muon row), with no remaining "
            "residual that requires an extra fit parameter."
        )
    else:
        lines.append(
            "**The γ offset is NOT closed within the scanned candidate "
            "list.** No single formula in this probe brings γ within "
            "1% of the lock AND reproduces the muon mass within 5%. "
            "The locked γ = 22.5 retains a small calibration on top of "
            "the radial barrier-sum scale."
        )
        # Best near-miss
        lines.append("")
        lines.append(
            f"Best near-miss (γ-closeness): `{bg['name']}` at "
            f"{bg['pct_diff_lepton']:+.3f}% γ-error, "
            f"{bg['muon_mass_err_pct']:.3f}% muon-error. "
            "Closing the remaining gap is the next concrete probe target."
        )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_gamma_offset_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
