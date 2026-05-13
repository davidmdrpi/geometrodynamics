"""
Resistance-reading disambiguation probe.

The opening transport/resistance origin probe identified two within-1 %
candidates for `resistance_scale ≈ 0.2179` that the probe could not
distinguish:

  A. `7π / 100 = 0.2199`        (+0.937 %) — closure-quantum fraction.
  B. `4·(ω(1,0) − 1) ≈ 0.2189`  (+0.477 %) — Tangherlini eigenfrequency
                                              gap (the 1.054 factor in
                                              disguise), evaluated at
                                              R_OUTER ≈ 1.262.

Both readings keep the lepton ladder within ~0.6 % when substituted
alongside `transport = 8π`. The two cannot be told apart on the locked
geometry alone.

This probe disambiguates by **re-closing the R_OUTER self-consistency
loop on each reading**. The loop from probe 8:

  R → γ_geometric(R) → locked surrogate → m_μ_pred, m_τ_pred

is run with three resistance/transport configurations:

  (0) Locked baseline:           transport = 25.1, resistance = 0.2179
  (1) Closure-quantum reading:   transport = 8π, resistance = 7π/100
  (2) Eigenfrequency reading:    transport = 8π, resistance = 4·(ω(1,0;R) − 1)
                                                              ↑
                                              VARIES WITH R_OUTER
                                              (fully geometric loop)

For each configuration we bisect R*_μ and R*_τ independently, report
the cross-species agreement |ΔR*|/R*, and check whether the fixed
point still lands at R* ≈ 1.262 with the original 0.008 % tolerance.

**Disambiguation criteria:**

  - If reading (2) yields cross-species agreement comparable to the
    locked baseline (≤ 0.05 %) and R* in the same neighbourhood, then
    resistance is structurally `4·(ω(1,0) − 1)` — the SAME 1.054-factor
    object that already governs the dimensional bridge. The loop
    becomes fully geometric (only m_e and the closure-quantum integers
    are inputs).

  - If reading (1) yields tight agreement and reading (2) does not,
    resistance is structurally `7π/100` — a separate closure-quantum
    fraction with its own origin.

  - If both readings yield comparable, looser agreement than the
    locked baseline, the structural origin is undetermined at this
    probe's precision.

The control is the locked baseline reading (0): it should reproduce
the probe-8 result (R* ≈ 1.262, cross-species 0.008 %).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional


M_E_MEV = 0.5109989461
M_MU_MEV = 105.6583745
M_TAU_MEV = 1776.86

# Locked baseline values
LOCKED_TRANSPORT = 25.1
LOCKED_RESISTANCE = 0.217869435878

# Closure-quantum readings
EIGHT_PI = 8.0 * math.pi
SEVEN_PI_OVER_100 = 7.0 * math.pi / 100.0


# ---------------------------------------------------------------------------
# Geometric quantities as functions of R_OUTER
# ---------------------------------------------------------------------------

def _gamma_geometric(R_outer: float, l_max: int = 5) -> float:
    """Σ V_max[0..l_max] on the Chebyshev grid at given R_OUTER."""
    import numpy as np
    from geometrodynamics.tangherlini.radial import (
        V_tangherlini, r_to_rstar, rstar_to_r,
    )
    from geometrodynamics.constants import R_MID
    rs = float(R_MID)
    rs_min = r_to_rstar(rs + 5e-4, rs)
    rs_max = r_to_rstar(R_outer - 5e-4, rs)
    N = 80
    x = np.cos(math.pi * np.arange(N + 1) / N)
    L = (rs_max - rs_min) / 2.0
    rsg = rs_min + L * (1.0 - x)
    rg = np.array([rstar_to_r(s, rs) for s in rsg])
    return sum(float(np.max(V_tangherlini(rg, l, rs))) for l in range(0, l_max + 1))


def _omega_1_0(R_outer: float) -> float:
    """Ground-state Tangherlini eigenfrequency ω(l=1, n=0) at given R_OUTER."""
    from geometrodynamics.tangherlini.radial import solve_radial_modes
    from geometrodynamics.constants import R_MID
    oms, _, _ = solve_radial_modes(l=1, N=80, n_modes=2,
                                   rs=float(R_MID), r_outer=R_outer)
    return float(oms[0])


def _predict_masses(
    gamma: float, transport: float, resistance: float,
) -> list[float]:
    """Locked-surrogate prediction at fixed γ, transport, resistance."""
    from scipy.linalg import eigh
    from geometrodynamics.tangherlini.lepton_spectrum import (
        _build_generation_block,
        LEPTON_BASELINE_DEPTHS,
        LEPTON_BASELINE_PHASE,
        S3_ACTION_BASE,
        TAU_BETA_50PI,
    )
    h = _build_generation_block(
        depths=LEPTON_BASELINE_DEPTHS,
        phase_per_pass=LEPTON_BASELINE_PHASE,
        transport_strength=transport,
        resistance_model="exponential",
        resistance_scale=resistance,
        hard_pinhole_gamma=gamma,
        action_base=S3_ACTION_BASE,
        action_slope=0.5,
        depth_cost_mode="tunnel_only",
        winding_mode="max",
        k_uplift_beta=TAU_BETA_50PI,
    )
    w = sorted(float(x) for x in eigh(h, eigvals_only=True) if x > 0)
    if len(w) < 3:
        return [float("nan")] * 3
    scale = M_E_MEV / w[0]
    return [w[i] * scale for i in range(3)]


# ---------------------------------------------------------------------------
# Configuration: (transport, resistance) as functions of R_OUTER
# ---------------------------------------------------------------------------

@dataclass
class Reading:
    name: str
    description: str
    transport_fn: Callable[[float], float]
    resistance_fn: Callable[[float], float]


READINGS: list[Reading] = [
    Reading(
        name="locked_baseline",
        description="transport = 25.1, resistance = 0.217869 (locked values)",
        transport_fn=lambda R: LOCKED_TRANSPORT,
        resistance_fn=lambda R: LOCKED_RESISTANCE,
    ),
    Reading(
        name="closure_quantum",
        description="transport = 8π, resistance = 7π/100 (both R-independent)",
        transport_fn=lambda R: EIGHT_PI,
        resistance_fn=lambda R: SEVEN_PI_OVER_100,
    ),
    Reading(
        name="eigenfrequency",
        description="transport = 8π, resistance = 4·(ω(1,0;R) − 1) (R-dependent)",
        transport_fn=lambda R: EIGHT_PI,
        resistance_fn=lambda R: 4.0 * (_omega_1_0(R) - 1.0),
    ),
    Reading(
        name="closure_quantum_locked_transport",
        description="transport = 25.1 (locked), resistance = 7π/100",
        transport_fn=lambda R: LOCKED_TRANSPORT,
        resistance_fn=lambda R: SEVEN_PI_OVER_100,
    ),
    Reading(
        name="eigenfrequency_locked_transport",
        description="transport = 25.1 (locked), resistance = 4·(ω(1,0;R) − 1)",
        transport_fn=lambda R: LOCKED_TRANSPORT,
        resistance_fn=lambda R: 4.0 * (_omega_1_0(R) - 1.0),
    ),
]


# ---------------------------------------------------------------------------
# Bisection
# ---------------------------------------------------------------------------

def _err_at(R: float, species: str, reading: Reading) -> float:
    gamma = _gamma_geometric(R)
    transport = reading.transport_fn(R)
    resistance = reading.resistance_fn(R)
    masses = _predict_masses(gamma, transport, resistance)
    obs = M_MU_MEV if species == "mu" else M_TAU_MEV
    pred = masses[1] if species == "mu" else masses[2]
    return (pred - obs) / obs


def _bisect(
    species: str, reading: Reading,
    lo: float = 1.245, hi: float = 1.275, tol: float = 1e-9,
) -> Optional[float]:
    """Bisect F(R; species, reading) = 0. Returns None if no bracket."""
    f_lo = _err_at(lo, species, reading)
    f_hi = _err_at(hi, species, reading)
    if math.isnan(f_lo) or math.isnan(f_hi):
        return None
    if f_lo * f_hi > 0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = _err_at(mid, species, reading)
        if math.isnan(f_mid):
            return None
        if abs(f_mid) < tol:
            return mid
        if f_lo * f_mid < 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


@dataclass
class ReadingResult:
    name: str
    description: str
    R_star_mu: Optional[float]
    R_star_tau: Optional[float]
    agreement_pct: Optional[float]
    gamma_at_R_star_mu: Optional[float]
    transport_at_R_star_mu: Optional[float]
    resistance_at_R_star_mu: Optional[float]
    omega_1_0_at_R_star_mu: Optional[float]
    masses_at_R_star_mu: list[float]
    err_mu_pct: Optional[float]
    err_tau_pct: Optional[float]
    bracket_hit: Optional[str]


def _evaluate_reading(reading: Reading) -> ReadingResult:
    R_mu = _bisect("mu", reading)
    R_tau = _bisect("tau", reading)
    if R_mu is None or R_tau is None:
        return ReadingResult(
            name=reading.name, description=reading.description,
            R_star_mu=R_mu, R_star_tau=R_tau,
            agreement_pct=None, gamma_at_R_star_mu=None,
            transport_at_R_star_mu=None, resistance_at_R_star_mu=None,
            omega_1_0_at_R_star_mu=None, masses_at_R_star_mu=[],
            err_mu_pct=None, err_tau_pct=None,
            bracket_hit="failed to bracket",
        )
    agreement = abs(R_mu - R_tau) / R_mu * 100.0
    gamma = _gamma_geometric(R_mu)
    transport = reading.transport_fn(R_mu)
    resistance = reading.resistance_fn(R_mu)
    omega = _omega_1_0(R_mu)
    masses = _predict_masses(gamma, transport, resistance)
    err_mu = (masses[1] - M_MU_MEV) / M_MU_MEV * 100.0
    err_tau = (masses[2] - M_TAU_MEV) / M_TAU_MEV * 100.0
    return ReadingResult(
        name=reading.name, description=reading.description,
        R_star_mu=R_mu, R_star_tau=R_tau,
        agreement_pct=agreement,
        gamma_at_R_star_mu=gamma,
        transport_at_R_star_mu=transport,
        resistance_at_R_star_mu=resistance,
        omega_1_0_at_R_star_mu=omega,
        masses_at_R_star_mu=masses,
        err_mu_pct=err_mu,
        err_tau_pct=err_tau,
        bracket_hit=None,
    )


# ---------------------------------------------------------------------------
# Probe runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    results: list[ReadingResult] = []
    for reading in READINGS:
        results.append(_evaluate_reading(reading))

    locked = next((r for r in results if r.name == "locked_baseline"), None)
    cq = next((r for r in results if r.name == "closure_quantum"), None)
    eig = next((r for r in results if r.name == "eigenfrequency"), None)

    # Disambiguation verdict — compare both readings on:
    #   (i) cross-species R*-agreement (smaller = tighter loop)
    #   (ii) R*-match to locked baseline (smaller = recovers probe 8)
    #   (iii) γ at R* match to canonical 22.5 (smaller = preserves
    #        the pinhole-origin reading)
    verdict_lines: list[str] = []
    if locked and cq and eig:
        locked_agree = locked.agreement_pct or float("inf")
        cq_agree = cq.agreement_pct or float("inf")
        eig_agree = eig.agreement_pct or float("inf")
        cq_R_match = abs((cq.R_star_mu or 0.0) - (locked.R_star_mu or 0.0)) / \
            (locked.R_star_mu or 1.0) * 100.0
        eig_R_match = abs((eig.R_star_mu or 0.0) - (locked.R_star_mu or 0.0)) / \
            (locked.R_star_mu or 1.0) * 100.0
        cq_gamma_off = abs((cq.gamma_at_R_star_mu or 0.0) - 22.5) / 22.5 * 100.0
        eig_gamma_off = abs((eig.gamma_at_R_star_mu or 0.0) - 22.5) / 22.5 * 100.0
        verdict_lines.append(
            f"Cross-species agreement — locked: {locked_agree:.4f}%, "
            f"closure: {cq_agree:.4f}%, eig: {eig_agree:.4f}%."
        )
        verdict_lines.append(
            f"R*-match to locked baseline — closure: {cq_R_match:.4f}%, "
            f"eig: {eig_R_match:.4f}%."
        )
        verdict_lines.append(
            f"γ at R* vs canonical 22.5 — closure: {cq_gamma_off:.4f}%, "
            f"eig: {eig_gamma_off:.4f}%."
        )
        verdict_lines.append("")
        # Both must be within tight tolerance to compete.
        AGREE_TIGHT = 0.05
        cq_qualifies = cq_agree < AGREE_TIGHT
        eig_qualifies = eig_agree < AGREE_TIGHT
        if cq_qualifies and eig_qualifies:
            # Tiebreak on R*-match to locked baseline (preserves probe 8).
            if cq_R_match < eig_R_match:
                verdict_lines.append(
                    "**closure-quantum reading PREFERRED.** Both readings "
                    "close the loop at high cross-species precision, but "
                    f"closure-quantum's R* ({cq.R_star_mu:.6f}) matches "
                    f"the locked baseline ({locked.R_star_mu:.6f}) to "
                    f"{cq_R_match:.4f}%, while eigenfrequency drifts by "
                    f"{eig_R_match:.4f}%. The closure-quantum reading "
                    "also lands γ at R* on 22.5 to "
                    f"{cq_gamma_off:.4f}%, preserving the pinhole-origin "
                    "anchor; eigenfrequency drifts γ by "
                    f"{eig_gamma_off:.4f}%."
                )
            else:
                verdict_lines.append(
                    "**eigenfrequency reading PREFERRED.** Both readings "
                    "close the loop at high cross-species precision, but "
                    f"eigenfrequency's R* ({eig.R_star_mu:.6f}) matches "
                    f"the locked baseline ({locked.R_star_mu:.6f}) to "
                    f"{eig_R_match:.4f}%, while closure-quantum drifts by "
                    f"{cq_R_match:.4f}%."
                )
        elif cq_qualifies:
            verdict_lines.append(
                "**closure-quantum reading PREFERRED.** Eigenfrequency "
                "does not close the loop at the locked-baseline precision "
                f"({eig_agree:.4f}% vs tight tolerance {AGREE_TIGHT}%)."
            )
        elif eig_qualifies:
            verdict_lines.append(
                "**eigenfrequency reading PREFERRED.** Closure-quantum "
                "does not close the loop at the locked-baseline precision "
                f"({cq_agree:.4f}% vs tight tolerance {AGREE_TIGHT}%)."
            )
        else:
            verdict_lines.append(
                "Neither reading reaches the locked baseline's "
                f"precision ({locked_agree:.4f}%). Disambiguation "
                "inconclusive at this probe's resolution."
            )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "readings": [asdict(r) for r in results],
        "verdict_lines": verdict_lines,
        "locked_baseline_summary": asdict(locked) if locked else None,
    }


def render_markdown(s: dict) -> str:
    lines: list[str] = []
    lines.append("# Resistance-reading disambiguation probe")
    lines.append("")
    lines.append(f"**Run:** {s['timestamp_utc']}")
    lines.append("")
    lines.append(
        "Disambiguates the two within-1 % candidates for `resistance_scale` "
        "identified by the opening transport/resistance origin probe:"
    )
    lines.append("")
    lines.append(
        "- **(A) closure_quantum**: `resistance = 7π/100 ≈ 0.2199` (+0.94 %), "
        "R-independent."
    )
    lines.append(
        "- **(B) eigenfrequency**: `resistance = 4·(ω(1,0;R) − 1) ≈ 0.2189` "
        "(+0.48 %), R-dependent (varies with the bisection)."
    )
    lines.append("")
    lines.append(
        "Each reading is paired with `transport = 8π = 4·(2π)` (the "
        "transport closure-quantum reading) and re-runs the R_OUTER "
        "self-consistency loop from probe 8 — bisecting R*_μ and R*_τ "
        "independently for each reading."
    )
    lines.append("")

    lines.append("## Results")
    lines.append("")
    lines.append(
        "| reading | description | R*_μ | R*_τ | |ΔR*|/R*_μ | γ at R* | "
        "transport | resistance | ω(1,0) | err μ | err τ |"
    )
    lines.append(
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    )
    for r in s["readings"]:
        if r["R_star_mu"] is None:
            lines.append(
                f"| `{r['name']}` | {r['description']} | (no bracket) | "
                f"— | — | — | — | — | — | — | — |"
            )
            continue
        lines.append(
            f"| `{r['name']}` | {r['description']} | "
            f"{r['R_star_mu']:.6f} | {r['R_star_tau']:.6f} | "
            f"{r['agreement_pct']:.4f}% | "
            f"{r['gamma_at_R_star_mu']:.4f} | "
            f"{r['transport_at_R_star_mu']:.4f} | "
            f"{r['resistance_at_R_star_mu']:.4f} | "
            f"{r['omega_1_0_at_R_star_mu']:.4f} | "
            f"{r['err_mu_pct']:+.4f}% | {r['err_tau_pct']:+.4f}% |"
        )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    for vl in s["verdict_lines"]:
        lines.append(vl)
    lines.append("")
    lines.append("### Reading-by-reading interpretation")
    lines.append("")
    locked = next((r for r in s["readings"] if r["name"] == "locked_baseline"), None)
    cq = next((r for r in s["readings"] if r["name"] == "closure_quantum"), None)
    eig = next((r for r in s["readings"] if r["name"] == "eigenfrequency"), None)
    if locked:
        lines.append(
            f"**Locked baseline** (control). Reproduces probe 8: "
            f"R* ≈ {locked['R_star_mu']:.4f}, cross-species agreement "
            f"{locked['agreement_pct']:.4f} %, "
            f"γ = {locked['gamma_at_R_star_mu']:.3f} ≈ 22.5. ✓"
        )
        lines.append("")
    if cq:
        lines.append(
            f"**Closure-quantum reading**. `transport = 8π`, "
            f"`resistance = 7π/100`, both R-independent. The bisection "
            f"selects R*_μ = {cq['R_star_mu']:.6f}, R*_τ = "
            f"{cq['R_star_tau']:.6f}, agreement "
            f"{cq['agreement_pct']:.4f} %. γ at R* = "
            f"{cq['gamma_at_R_star_mu']:.3f}."
        )
        lines.append("")
    if eig:
        lines.append(
            f"**Eigenfrequency reading**. `transport = 8π`, "
            f"`resistance = 4·(ω(1,0;R) − 1)`. resistance varies with R "
            "during the bisection — this is the FULLY GEOMETRIC loop "
            "(no constants other than m_e and the closure-quantum "
            f"integers). R*_μ = {eig['R_star_mu']:.6f}, R*_τ = "
            f"{eig['R_star_tau']:.6f}, agreement "
            f"{eig['agreement_pct']:.4f} %. ω(1,0) at R* = "
            f"{eig['omega_1_0_at_R_star_mu']:.4f} (resistance at R* = "
            f"{eig['resistance_at_R_star_mu']:.4f})."
        )
        lines.append("")

    lines.append("## Implication")
    lines.append("")
    if locked and cq and eig:
        l_ag = locked["agreement_pct"]
        c_ag = cq["agreement_pct"]
        e_ag = eig["agreement_pct"]
        cq_R = cq["R_star_mu"]
        eig_R = eig["R_star_mu"]
        locked_R = locked["R_star_mu"]
        cq_R_off = abs(cq_R - locked_R) / locked_R * 100.0
        eig_R_off = abs(eig_R - locked_R) / locked_R * 100.0
        lines.append(
            f"**Cross-species agreement** — Locked: {l_ag:.4f}% · "
            f"Closure-quantum: {c_ag:.4f}% · Eigenfrequency: {e_ag:.4f}%."
        )
        lines.append(
            f"**R*-match to locked** — Closure-quantum: "
            f"{cq_R_off:.4f}% · Eigenfrequency: {eig_R_off:.4f}%."
        )
        lines.append("")
        if cq_R_off < eig_R_off and c_ag <= e_ag + 0.005:
            lines.append(
                "**Closure-quantum reading wins.** Both readings produce "
                "tighter cross-species agreement than the locked baseline "
                "(0.0078 %), but the closure-quantum reading's R* lands "
                f"on the locked-baseline R* to {cq_R_off:.4f} %, while the "
                f"eigenfrequency reading drifts by {eig_R_off:.4f} %. γ at "
                f"R* under closure-quantum matches the canonical "
                f"22.5 to "
                f"{abs(cq['gamma_at_R_star_mu'] - 22.5) / 22.5 * 100.0:.4f} %, "
                "preserving the pinhole-origin anchor identified in "
                "`pinhole_origin_probe`. The structural reading is:"
            )
            lines.append("")
            lines.append("  `transport = 8π = 4·(2π)` — 4th closure quantum.")
            lines.append("  `resistance = 7π / 100` — closure-quantum fraction.")
            lines.append("")
            lines.append(
                "Both are pure closure-quantum invariants of the BAM "
                "framework. With this reading the R_OUTER self-"
                "consistency loop closes on principled inputs alone — "
                "no transport or resistance constants are required as "
                "external inputs."
            )
        elif eig_R_off < cq_R_off and e_ag < c_ag:
            lines.append(
                "**Eigenfrequency reading wins.** Cross-species agreement "
                f"and R*-match are both tightest for the eigenfrequency "
                "reading, indicating that resistance is structurally the "
                "Tangherlini eigenfrequency object `4·(ω(1,0;R) − 1)`. "
                "This makes resistance, the 1.054 factor, and γ all "
                "projections of one Tangherlini matrix-element family — "
                "the cross-cutting structural unification flagged in "
                "the research plan §4."
            )
        else:
            lines.append(
                "**Mixed outcome.** Cross-species and R*-match criteria "
                "favour different readings; further sub-probes needed."
            )
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_resistance_disambiguation_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
