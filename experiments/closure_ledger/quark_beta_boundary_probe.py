"""
Quark β: focused boundary-correction probe.

The previous quark_beta_origin_probe surfaced two principled near-miss
patterns:

  18 · (k_5² + 1) = 468  vs  N_q = 466   (off by  +2)
  14 · (k_5² + 1) = 364  vs  ΔN  = 366   (off by  −2)

This probe asks the focused question: what is the *cleanest* unit u
and boundary-correction δ such that the locked closure-quantum
integers fit u-multiplet form

    N(target)  =  m · u  +  δ

simultaneously across all three targets — N_l = 100, N_q = 466, and
ΔN = 366 = N_q − N_l? A clean reading requires δ to come from a
small natural set ({0, ±1, ±2, ±k_5, ±(k_5+1)}, etc.) AND to be
consistent (or at least reduced) across targets.

Structural units scanned:

    k_5             = 5
    k_5²            = 25
    k_5² + 1        = 26
    k_5² − 1        = 24
    (k_5 + 1)²      = 36
    (k_5 − 1)²      = 16
    k_5 (k_5 + 2)   = 35   (S³ angular eigenvalue at ℓ = k_5)
    2 k_5           = 10
    (2 k_5)²        = 100
    4 (k_5² + 1)    = 104

For each unit, we fit `m = round(N / u)` and report the residual
δ = N − m · u for each of the three targets. Rank by joint cleanness
score `|δ_l| + |δ_q| + |δ_ΔN|` and flag candidates whose δ values are
also drawn from a small natural set.

Robustness side-test (per `docs/quark_axioms.md` §8): N_q drifts by
~90 units under per-species mass perturbations. If the structural
reading is `m · k_5 + 1` with k_5 = 5, then a 90-unit drift in N_q
corresponds to 18 units of drift in m (with δ fixed at +1), which is
a quantitative prediction the quark code's existing PDG-perturbation
runs can be re-read against.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


N_LEPTON = 100
N_QUARK = 466
DELTA_N = N_QUARK - N_LEPTON   # 366
K_5 = 5

# Natural boundary-correction set: small integers tied to k_5 or to
# elementary structural counts (Z₂ partition residue ±1, double-cover ±2,
# half-shell ±k_5, etc.).
NATURAL_DELTAS = {
    0,
    1, -1,
    2, -2,
    K_5, -K_5,                      # ±5
    K_5 + 1, -(K_5 + 1),            # ±6
    K_5 - 1, -(K_5 - 1),            # ±4
    2 * K_5, -2 * K_5,              # ±10
}


@dataclass
class FitResult:
    unit_name: str
    unit_value: int
    m_lepton: int
    delta_lepton: int
    m_quark: int
    delta_quark: int
    m_delta: int
    delta_gap: int
    sum_abs_delta: int
    deltas_in_natural_set: bool
    deltas_consistent_across_quark_and_gap: bool
    deltas_zero_for_lepton: bool


def _fit_unit(name: str, u: int) -> FitResult:
    if u <= 0:
        raise ValueError(f"unit must be positive: {u}")
    def _fit(N):
        m = round(N / u)
        d = N - m * u
        return m, d
    m_l, d_l = _fit(N_LEPTON)
    m_q, d_q = _fit(N_QUARK)
    m_d, d_d = _fit(DELTA_N)
    in_natural = (
        d_l in NATURAL_DELTAS
        and d_q in NATURAL_DELTAS
        and d_d in NATURAL_DELTAS
    )
    consistent = (d_q == d_d)
    zero_lepton = (d_l == 0)
    return FitResult(
        unit_name=name,
        unit_value=u,
        m_lepton=m_l, delta_lepton=d_l,
        m_quark=m_q, delta_quark=d_q,
        m_delta=m_d, delta_gap=d_d,
        sum_abs_delta=abs(d_l) + abs(d_q) + abs(d_d),
        deltas_in_natural_set=in_natural,
        deltas_consistent_across_quark_and_gap=consistent,
        deltas_zero_for_lepton=zero_lepton,
    )


def _structural_units() -> list[tuple[str, int]]:
    return [
        ("k_5",                K_5),
        ("k_5²",               K_5 ** 2),
        ("k_5² + 1",           K_5 ** 2 + 1),
        ("k_5² − 1",           K_5 ** 2 - 1),
        ("(k_5 + 1)²",         (K_5 + 1) ** 2),
        ("(k_5 − 1)²",         (K_5 - 1) ** 2),
        ("k_5(k_5 + 2)",       K_5 * (K_5 + 2)),
        ("2 k_5",              2 * K_5),
        ("(2 k_5)²",           (2 * K_5) ** 2),
        ("4 (k_5² + 1)",       4 * (K_5 ** 2 + 1)),
        ("k_5 + 1",            K_5 + 1),
        ("k_5(k_5+1)/2",       K_5 * (K_5 + 1) // 2),
    ]


# ---------------------------------------------------------------------------
# Robustness test: under N_q drift, does m or δ move?
# ---------------------------------------------------------------------------

@dataclass
class RobustnessSnapshot:
    """Re-read of the locked quark spectrum's N stability per quark_axioms §8."""
    perturbation: str
    N_q_at_perturbation: int
    fit_unit: str
    m_at_perturbation: int
    delta_at_perturbation: int
    delta_drift_from_baseline: int


def _robustness_snapshots(unit_name: str, u: int) -> list[RobustnessSnapshot]:
    """
    Document what the existing N-stability ablation says, restated under
    the (m, δ) decomposition.

    Source: `docs/quark_axioms.md` §8 N-stability table (12 logged
    ablation points from `experiment_n_ablation.py`). All N values are
    EXACT logged values; no extrapolation.

    The decisive N-ablation conclusion of §8 is that N is a compensator,
    not a topological invariant. This snapshot routine just restates
    the table under the (m, δ) decomposition; the audit-grade analysis
    is in `quark_beta_robustness_audit.py`.
    """
    snapshots: list[tuple[str, int]] = [
        ("baseline (anchor=d, PDG, min_eig)", 466),
        ("PDG × 1.10 (uniform scale)", 466),
        ("PDG × 0.90 (uniform scale)", 466),
        ("anchor = s", 476),
        ("anchor = c", 474),
        ("anchor = b", 474),
        ("anchor = t", 482),
        ("c × 1.10", 432),
        ("b × 1.10", 494),
        ("t × 1.10", 494),
        ("t × 0.90", 440),
        ("all ±5% (deterministic)", 510),
    ]
    out: list[RobustnessSnapshot] = []
    baseline_delta = N_QUARK - round(N_QUARK / u) * u
    for label, N in snapshots:
        m = round(N / u)
        d = N - m * u
        out.append(RobustnessSnapshot(
            perturbation=label,
            N_q_at_perturbation=N,
            fit_unit=unit_name,
            m_at_perturbation=m,
            delta_at_perturbation=d,
            delta_drift_from_baseline=d - baseline_delta,
        ))
    return out


# ---------------------------------------------------------------------------
# Probe runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    units = _structural_units()
    fits = [_fit_unit(name, u) for name, u in units]

    # Rank by sum |δ|, with tie-break by the cleanness flags.
    ranked = sorted(
        fits,
        key=lambda f: (
            f.sum_abs_delta,
            not f.deltas_in_natural_set,
            not f.deltas_consistent_across_quark_and_gap,
            not f.deltas_zero_for_lepton,
        ),
    )
    best = ranked[0]
    robustness = _robustness_snapshots(best.unit_name, best.unit_value)

    # Check the structural prediction: under N_q drifts, does δ stay
    # invariant (m absorbs the drift) or does δ drift (m stays)?
    deltas_at_drifts = [s.delta_at_perturbation for s in robustness]
    delta_invariance_count = sum(
        1 for d in deltas_at_drifts if d == best.delta_quark
    )
    delta_invariance_rate = delta_invariance_count / max(len(deltas_at_drifts), 1)

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "targets": {
            "N_lepton": N_LEPTON,
            "N_quark": N_QUARK,
            "delta_N": DELTA_N,
        },
        "k_5": K_5,
        "natural_delta_set": sorted(NATURAL_DELTAS),
        "all_fits": [asdict(f) for f in fits],
        "ranked_fits_by_joint_cleanness": [asdict(f) for f in ranked],
        "best_fit": asdict(best),
        "robustness_snapshots_under_best_fit": [
            asdict(s) for s in robustness
        ],
        "delta_invariance_rate_under_documented_drifts": delta_invariance_rate,
        "structural_prediction": (
            "If N_q = m · u + δ is the correct decomposition with δ "
            "fixed by structure, then per-species mass perturbations "
            "should leave δ invariant and drift m proportionally. The "
            "quark_axioms §8 ~90-unit N drift translates to "
            f"{round(90 / max(best.unit_value, 1))} units of m drift "
            f"under {best.unit_name} = {best.unit_value} as the structural "
            "increment."
        ),
    }


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Quark β: focused boundary-correction probe")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append(
        f"**Targets:** N_l = {summary['targets']['N_lepton']}, "
        f"N_q = {summary['targets']['N_quark']}, "
        f"ΔN = {summary['targets']['delta_N']}; "
        f"k_5 = {summary['k_5']}."
    )
    lines.append(
        "**Natural δ set:** "
        f"{summary['natural_delta_set']}."
    )
    lines.append("")

    lines.append("## All structural-unit fits, ranked by joint cleanness")
    lines.append("")
    lines.append(
        "| unit | u | N_l = m·u+δ | N_q = m·u+δ | ΔN = m·u+δ | Σ|δ| | "
        "δ ∈ natural? | δ_q == δ_ΔN? | δ_l = 0? |"
    )
    lines.append("|---|---:|---|---|---|---:|---|---|---|")
    for f in summary["ranked_fits_by_joint_cleanness"]:
        l_str = f"{f['m_lepton']}·{f['unit_value']} {f['delta_lepton']:+d}"
        q_str = f"{f['m_quark']}·{f['unit_value']} {f['delta_quark']:+d}"
        d_str = f"{f['m_delta']}·{f['unit_value']} {f['delta_gap']:+d}"
        nat = "✓" if f["deltas_in_natural_set"] else "—"
        cons = "✓" if f["deltas_consistent_across_quark_and_gap"] else "—"
        zlep = "✓" if f["deltas_zero_for_lepton"] else "—"
        lines.append(
            f"| `{f['unit_name']}` | {f['unit_value']} | {l_str} | {q_str} | "
            f"{d_str} | {f['sum_abs_delta']} | {nat} | {cons} | {zlep} |"
        )
    lines.append("")

    best = summary["best_fit"]
    lines.append("## Best fit")
    lines.append("")
    lines.append(
        f"**Unit:** `{best['unit_name']}` = {best['unit_value']}."
    )
    lines.append("")
    lines.append("**Decomposition:**")
    lines.append("")
    lines.append(
        f"- `N_lepton  = {best['m_lepton']} · {best['unit_value']} "
        f"{best['delta_lepton']:+d}` "
        f"⇒ residue δ = **{best['delta_lepton']}**"
    )
    lines.append(
        f"- `N_quark   = {best['m_quark']} · {best['unit_value']} "
        f"{best['delta_quark']:+d}` "
        f"⇒ residue δ = **{best['delta_quark']}**"
    )
    lines.append(
        f"- `ΔN_q-l   = {best['m_delta']} · {best['unit_value']} "
        f"{best['delta_gap']:+d}` "
        f"⇒ residue δ = **{best['delta_gap']}**"
    )
    lines.append("")
    if best["deltas_zero_for_lepton"] and best["deltas_consistent_across_quark_and_gap"]:
        lines.append(
            "**Structural cleanness criteria met:**"
            " δ_lepton = 0 (lepton sector closes exactly on a multiple"
            f" of {best['unit_name']}) AND"
            f" δ_quark = δ_gap = {best['delta_quark']} (the quark sector"
            f" carries a single fixed boundary correction +{best['delta_quark']}"
            " across both targets)."
        )
        lines.append("")
        lines.append(
            f"This is consistent with reading the quark closure-quantum"
            f" count as `m · k_5 + 1`, where m is the integer winding"
            f" number that absorbs all per-species mass-perturbation drift,"
            f" and the +1 is a single sector-level boundary correction"
            f" tied to the orientation-reversing throat closure (parallel"
            f" to the lepton-sector exact reading `N_l = 20 · k_5`)."
        )
        lines.append("")

    lines.append("## Robustness re-read under documented N-drifts")
    lines.append("")
    lines.append(
        f"From `docs/quark_axioms.md` §8 N-stability ablation, restated"
        f" under the (m, δ) decomposition with u = `{best['unit_name']}`:"
    )
    lines.append("")
    lines.append("| perturbation | N_q | m | δ | δ-drift |")
    lines.append("|---|---:|---:|---:|---:|")
    for s in summary["robustness_snapshots_under_best_fit"]:
        lines.append(
            f"| {s['perturbation']} | {s['N_q_at_perturbation']} | "
            f"{s['m_at_perturbation']} | {s['delta_at_perturbation']:+d} | "
            f"{s['delta_drift_from_baseline']:+d} |"
        )
    lines.append("")
    rate = summary["delta_invariance_rate_under_documented_drifts"]
    lines.append(
        f"**δ-invariance rate under documented drifts:** "
        f"{100 * rate:.0f}% of perturbations leave δ at its baseline "
        f"value of {best['delta_quark']:+d}."
    )
    lines.append("")
    if rate >= 0.5:
        lines.append(
            "**Most documented drifts are absorbed by m**, with δ staying"
            " near its baseline. Consistent with the structural reading."
        )
    else:
        lines.append(
            f"**Only {100 * rate:.0f}% of documented drifts leave δ at"
            " baseline**; under per-species and anchor perturbations,"
            " both m and δ wander. The structural reading is descriptively"
            " useful for the BASELINE locked value but is NOT a robust"
            " topological invariant — consistent with §8's headline"
            " conclusion that N (and therefore both m and δ) is a fit"
            " compensator. See `quark_beta_robustness_audit.py` for the"
            " full audit."
        )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if (
        best["deltas_zero_for_lepton"]
        and best["deltas_consistent_across_quark_and_gap"]
        and abs(best["delta_quark"]) <= 1
    ):
        lines.append(
            f"**The cleanest structural reading is `m · {best['unit_name']}"
            f" + {best['delta_quark']}`** with the lepton sector exact"
            f" (δ = 0) and the quark sector carrying a single fixed"
            f" boundary correction δ = {best['delta_quark']:+d}. The"
            f" pattern is suggestive but not yet derived: `+1` could be"
            " a Z₂ partition residue (orientation-flip closure), an l = 0"
            " s-wave closure quantum, or a single-mouth boundary"
            " contribution. Identifying which is the next probe target."
        )
    elif best["sum_abs_delta"] <= 4:
        lines.append(
            f"**Cleanest unit `{best['unit_name']}` gives Σ|δ|"
            f" = {best['sum_abs_delta']}**, but the δ values are not"
            " uniformly small or consistent across targets. The"
            " structural pattern is suggestive but not yet pinned."
        )
    else:
        lines.append(
            f"**No structural unit in the catalog gives a clean joint"
            f" decomposition.** Best is `{best['unit_name']}` with"
            f" Σ|δ| = {best['sum_abs_delta']}; further candidates beyond"
            " this catalog (e.g. structural units constructed from"
            " multiple shells) would be needed to make progress."
        )

    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_quark_beta_boundary_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
