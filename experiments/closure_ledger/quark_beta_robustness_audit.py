"""
Quark β: §8 robustness audit sub-probe.

The decomposition probe reported a candidate structural reading where
the closure-quantum count splits as

    N_q  =  m_q · k_5  +  δ_q       with  m_q = 93,  δ_q = +1.

That probe (and the boundary-correction probe before it) used a
SYNTHETIC ablation table with extrapolated N values for the c+b
combined drift. This audit uses the EXACT N values logged in
`docs/quark_axioms.md` §8 (commit `f7c90b5`, the decisive N-ablation
table) and re-evaluates the (m, δ) decomposition.

Targets:

  Are the previous probes' structural claims preserved when the audit
  uses the actual logged N values rather than extrapolated ones?

Structural claims to audit:

  (1) δ = +1 holds for the BASELINE and the uniform-mass-scale runs
      (the "no-physics-change" perturbations).
  (2) δ = +1 is preserved across per-species mass perturbations.
  (3) δ = +1 is preserved across anchor-species changes.
  (4) The δ-invariance rate is ≥ 50% across all logged ablations.

The audit reports the actual (m, δ) for each ablation and grades the
four claims. A claim is FALSE if the data contradicts it; UNVERIFIED
if it would have required data not logged in §8.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


K_5 = 5
N_QUARK_BASELINE = 466


@dataclass
class Ablation:
    """One row of the §8 N-stability table."""
    label: str
    N: int
    fit_residual_pct: float
    perturbation_class: str   # "baseline" | "uniform_scale" | "anchor"
                              # | "single_species" | "all_species" | "convention"


# §8 N-stability ablation table from `docs/quark_axioms.md`. These are
# the exact logged values from `experiment_n_ablation.py`.
ABLATIONS = [
    Ablation("baseline (anchor=d, PDG, min_eig)", 466, 1.6, "baseline"),
    Ablation("PDG × 1.10 (uniform scale)",        466, 1.6, "uniform_scale"),
    Ablation("PDG × 0.90 (uniform scale)",        466, 1.6, "uniform_scale"),
    Ablation("anchor = s",                        476, 1.1, "anchor"),
    Ablation("anchor = c",                        474, 1.0, "anchor"),
    Ablation("anchor = b",                        474, 1.9, "anchor"),
    Ablation("anchor = t",                        482, 1.6, "anchor"),
    Ablation("c × 1.10",                          432, 1.8, "single_species"),
    Ablation("b × 1.10",                          494, 4.2, "single_species"),
    Ablation("t × 1.10",                          494, 5.5, "single_species"),
    Ablation("t × 0.90",                          440, 3.7, "single_species"),
    Ablation("all ±5% (deterministic)",           510, 4.6, "all_species"),
    # spectrum_zero=second_min row excluded — fit fails.
]


@dataclass
class AblationDecomposition:
    label: str
    perturbation_class: str
    N: int
    m: int                       # N // K_5
    delta_unsigned: int          # N % K_5  (in [0, k_5))
    delta_signed: int            # in (−k_5/2, k_5/2]
    delta_drift_from_baseline: int
    fit_residual_pct: float


def _decompose(a: Ablation, baseline_signed: int) -> AblationDecomposition:
    m = a.N // K_5
    d_unsigned = a.N % K_5
    d_signed = d_unsigned if d_unsigned <= K_5 // 2 else d_unsigned - K_5
    return AblationDecomposition(
        label=a.label,
        perturbation_class=a.perturbation_class,
        N=a.N,
        m=m,
        delta_unsigned=d_unsigned,
        delta_signed=d_signed,
        delta_drift_from_baseline=d_signed - baseline_signed,
        fit_residual_pct=a.fit_residual_pct,
    )


@dataclass
class ClaimVerdict:
    name: str
    description: str
    expected: str
    observed: str
    verdict: str              # "TRUE" | "FALSE" | "UNVERIFIED"


def _audit_claims(
    decomps: list[AblationDecomposition],
) -> list[ClaimVerdict]:
    """Grade the four structural claims against the actual decompositions."""
    by_class: dict[str, list[AblationDecomposition]] = {}
    for d in decomps:
        by_class.setdefault(d.perturbation_class, []).append(d)

    results: list[ClaimVerdict] = []

    # Claim 1: δ = +1 holds for baseline + uniform_scale.
    no_physics = (
        by_class.get("baseline", []) + by_class.get("uniform_scale", [])
    )
    delta_1_count = sum(1 for d in no_physics if d.delta_signed == 1)
    rate = delta_1_count / max(len(no_physics), 1)
    results.append(ClaimVerdict(
        name="Claim 1: δ = +1 in baseline / uniform-scale runs",
        description=(
            "Under perturbations that do NOT change the physics (just "
            "the MeV anchor or a uniform mass rescaling), δ should "
            "equal the baseline +1."
        ),
        expected=f"all {len(no_physics)} runs at δ = +1",
        observed=(
            f"{delta_1_count}/{len(no_physics)} at δ = +1 "
            f"(rate {100*rate:.0f}%)"
        ),
        verdict="TRUE" if rate == 1.0 else (
            "FALSE" if rate < 0.5 else "PARTIAL"
        ),
    ))

    # Claim 2: δ = +1 preserved under single-species perturbations.
    single = by_class.get("single_species", [])
    delta_1_single = sum(1 for d in single if d.delta_signed == 1)
    rate_s = delta_1_single / max(len(single), 1)
    deltas_single = [d.delta_signed for d in single]
    results.append(ClaimVerdict(
        name="Claim 2: δ = +1 across single-species perturbations",
        description=(
            "Under single-species ±10% perturbations (c, b, t each "
            "×1.10 or ×0.90), δ should remain at +1 if the structural "
            "reading is robust."
        ),
        expected=f"all {len(single)} runs at δ = +1",
        observed=(
            f"{delta_1_single}/{len(single)} at δ = +1; "
            f"observed δ values: {sorted(deltas_single)}"
        ),
        verdict="TRUE" if rate_s == 1.0 else (
            "FALSE" if rate_s < 0.5 else "PARTIAL"
        ),
    ))

    # Claim 3: δ = +1 preserved under anchor-species changes.
    anchors = by_class.get("anchor", [])
    delta_1_anchor = sum(1 for d in anchors if d.delta_signed == 1)
    rate_a = delta_1_anchor / max(len(anchors), 1)
    deltas_anchor = [d.delta_signed for d in anchors]
    results.append(ClaimVerdict(
        name="Claim 3: δ = +1 across anchor-species changes",
        description=(
            "Anchor-species choice (s, c, b, t instead of d) is a "
            "convention; structural reading expects δ = +1 throughout."
        ),
        expected=f"all {len(anchors)} runs at δ = +1",
        observed=(
            f"{delta_1_anchor}/{len(anchors)} at δ = +1; "
            f"observed δ values: {sorted(deltas_anchor)}"
        ),
        verdict="TRUE" if rate_a == 1.0 else (
            "FALSE" if rate_a < 0.5 else "PARTIAL"
        ),
    ))

    # Claim 4: overall δ-invariance rate ≥ 50%.
    overall_rate = sum(1 for d in decomps if d.delta_signed == 1) / len(decomps)
    results.append(ClaimVerdict(
        name="Claim 4: overall δ-invariance rate ≥ 50%",
        description=(
            "Aggregated over every logged ablation, at least half "
            "should sit at δ = +1 if the boundary correction is a "
            "structural invariant."
        ),
        expected="≥ 50% of runs at δ = +1",
        observed=f"{100*overall_rate:.0f}% of runs at δ = +1",
        verdict="TRUE" if overall_rate >= 0.5 else "FALSE",
    ))

    return results


# ---------------------------------------------------------------------------
# Probe runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    baseline = next(a for a in ABLATIONS if a.label.startswith("baseline"))
    base_d = _decompose(baseline, baseline_signed=0)
    decomps = [_decompose(a, baseline_signed=base_d.delta_signed) for a in ABLATIONS]
    claims = _audit_claims(decomps)

    counter = Counter(d.delta_signed for d in decomps)
    counter_unsigned = Counter(d.delta_unsigned for d in decomps)
    n_at_plus_1 = counter.get(1, 0)
    n_total = len(decomps)
    plus_one_rate = n_at_plus_1 / n_total

    # Drift envelope of m: how much does m wander?
    m_values = [d.m for d in decomps]
    m_min, m_max = min(m_values), max(m_values)
    m_baseline = base_d.m

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "k_5": K_5,
        "baseline_N": N_QUARK_BASELINE,
        "ablation_decompositions": [asdict(d) for d in decomps],
        "delta_signed_distribution": dict(counter),
        "delta_unsigned_distribution": dict(counter_unsigned),
        "plus_one_rate": plus_one_rate,
        "m_envelope": {
            "baseline_m": m_baseline,
            "m_min": m_min,
            "m_max": m_max,
            "m_drift_pm": max(abs(m_max - m_baseline), abs(m_baseline - m_min)),
        },
        "claim_verdicts": [asdict(c) for c in claims],
        "n_claims_true": sum(1 for c in claims if c.verdict == "TRUE"),
        "n_claims_partial": sum(1 for c in claims if c.verdict == "PARTIAL"),
        "n_claims_false": sum(1 for c in claims if c.verdict == "FALSE"),
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Quark β: §8 robustness audit sub-probe")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append(
        f"**k_5 = {summary['k_5']}**, "
        f"**baseline N = {summary['baseline_N']}** "
        f"(m = {summary['m_envelope']['baseline_m']}, "
        f"δ_signed = +1)."
    )
    lines.append("")
    lines.append(
        "Audits the structural reading `N_q = m·k_5 + δ` against the "
        "EXACT §8 N-ablation values logged in `docs/quark_axioms.md`. "
        "Earlier probes used extrapolated N values for some perturbation "
        "classes; this audit drops them and reports only what the §8 "
        "log actually documents."
    )
    lines.append("")

    lines.append("## Per-ablation (m, δ) decompositions")
    lines.append("")
    lines.append(
        "| ablation | class | N | m | δ_signed | drift | fit residual |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for d in summary["ablation_decompositions"]:
        lines.append(
            f"| {d['label']} | {d['perturbation_class']} | {d['N']} | "
            f"{d['m']} | {d['delta_signed']:+d} | "
            f"{d['delta_drift_from_baseline']:+d} | "
            f"{d['fit_residual_pct']:.1f}% |"
        )
    lines.append("")

    lines.append("## δ distribution across all ablations")
    lines.append("")
    lines.append("| signed δ | count | fraction |")
    lines.append("|---:|---:|---:|")
    total = len(summary["ablation_decompositions"])
    for d_val, n in sorted(summary["delta_signed_distribution"].items()):
        lines.append(f"| {d_val:+d} | {n} | {100*n/total:.0f}% |")
    lines.append("")

    env = summary["m_envelope"]
    lines.append("## m envelope")
    lines.append("")
    lines.append(
        f"baseline m = {env['baseline_m']}, m range = "
        f"[{env['m_min']}, {env['m_max']}], "
        f"drift envelope ±{env['m_drift_pm']} units."
    )
    lines.append("")

    lines.append("## Claim verdicts")
    lines.append("")
    for c in summary["claim_verdicts"]:
        marker = {
            "TRUE": "✓",
            "PARTIAL": "~",
            "FALSE": "✗",
            "UNVERIFIED": "?",
        }.get(c["verdict"], "?")
        lines.append(f"### {c['name']}")
        lines.append(f"- {c['description']}")
        lines.append(f"- **Expected:** {c['expected']}")
        lines.append(f"- **Observed:** {c['observed']}")
        lines.append(f"- **Verdict:** {marker} **{c['verdict']}**")
        lines.append("")

    n_true = summary["n_claims_true"]
    n_partial = summary["n_claims_partial"]
    n_false = summary["n_claims_false"]
    n_total_claims = n_true + n_partial + n_false

    lines.append("## Verdict summary")
    lines.append("")
    lines.append(
        f"- **TRUE:** {n_true}/{n_total_claims}\n"
        f"- **PARTIAL:** {n_partial}/{n_total_claims}\n"
        f"- **FALSE:** {n_false}/{n_total_claims}\n"
    )
    lines.append("")

    lines.append("## Headline correction to prior probes")
    lines.append("")
    plus_rate = summary["plus_one_rate"]
    if plus_rate < 0.5:
        lines.append(
            f"**The earlier boundary-correction probe overstated the "
            f"+1 invariance.** Using the actual §8 ablation N values, "
            f"only **{100*plus_rate:.0f}%** of runs sit at δ = +1; the "
            f"remainder are distributed across δ ∈ "
            f"{{{', '.join(str(d) for d in sorted(summary['delta_signed_distribution']))}}}. "
            f"The δ-distribution is approximately uniform on the integers "
            f"in (−k_5/2, k_5/2] = (−2, 2], consistent with the §8 "
            f"conclusion that 'N is a compensator, not a topological "
            f"invariant'."
        )
        lines.append("")
        lines.append(
            "**Implications for the structural decomposition:**"
        )
        lines.append("")
        lines.append(
            "- The structural reading `N_q = ((k_5−1)·k_5 + 2·k_5(k_5+2) "
            "+ N_c)·k_5 + (N_c − 2)` correctly fits the *baseline* "
            "value 466 and the *uniform-scale* invariance, but does NOT "
            "predict the perturbation behaviour that §8 documents. "
            "Under per-species or anchor changes, both m and δ drift."
        )
        lines.append(
            "- The +1 = N_c − 2 boundary-correction interpretation is "
            "consistent at the baseline but is NOT a robust topological "
            "invariant. It survives uniform-mass-scale perturbations "
            "(those don't change the underlying physics, just the MeV "
            "anchor) but breaks under anchor-species and per-species "
            "perturbations."
        )
        lines.append(
            "- This is consistent with the §8 conclusion: **N (and "
            "therefore m and δ both) is a fit compensator. The "
            "structural decomposition is a useful descriptor of the "
            "BASELINE locked value, not a derivation of N.**"
        )
    else:
        lines.append(
            f"The +1 invariance holds for {100*plus_rate:.0f}% of "
            "ablations — sufficient to support the structural reading "
            "as a robust feature."
        )
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append(
        "The structural decomposition `m_q = (k_5−1)·k_5 + 2·k_5(k_5+2) "
        "+ N_c` and boundary correction `δ = N_c − 2` remain "
        "**descriptively useful for the baseline** and the "
        "uniform-scale runs — these read off the locked N = 466 in a "
        "geometrically meaningful way. They are NOT a derivation of N "
        "in the strong sense: under per-species and anchor "
        "perturbations, both m and δ wander, exactly as §8 reports. "
        "The repo's own conclusion stands: **β_quark is "
        "phenomenological**, and the cleanest structural reading is a "
        "post-hoc identification of the locked value, not a "
        "first-principles prediction."
    )
    lines.append("")
    lines.append(
        "The next concrete sub-target is to ask whether **a different "
        "perturbation-robust quantity** (e.g. the NN_l + N_c piece "
        "alone, or the lepton-block sub-piece) is what's actually "
        "topologically locked, with the rest being the compensator."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_quark_beta_robustness_audit"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
