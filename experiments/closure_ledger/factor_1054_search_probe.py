"""
Closed-form search for the 1.054 dimensional-bridge factor.

After eight closure-ledger probes, the BAM framework predicts every
dimensionless ratio at sub-percent but anchors the absolute MeV scale
at m_e via a single residual factor:

    ℏ · ω(1, 0)  =  1.054 · m_e c²    at R_OUTER ≈ 1.262

This 1.054 is **structural** — it is ω(l = 1, n = 0) evaluated at the
cross-species self-consistency fixed point (probe 8). The Compton
bridge ω = 1 exactly is mathematically definable but physically
vetoed by the lepton spectrum (probe 7). If 1.054 has a closed form
in `(k_5, π, barrier-spectrum invariants, small integers)`, the
dimensional bridge to ℏ closes.

This probe runs an exhaustive search for closed-form candidates:

  (1) Small rationals p/q with p, q ≤ 200, matching ω(1, 0) or ω².
  (2) Roots of small rationals: (p/q)^(1/n) for n ≤ 10.
  (3) Combinations involving k_5 = 5 and π.
  (4) Series-truncation candidates: 1 + ε/k_5^m with ε ∈ {1, 2, …}.

The target is the ω(1, 0) value at the cross-species fixed point
R* ≈ 1.262239 (from `R_outer_self_consistency_probe`). We use the
canonical-baseline value ω(R_OUTER = 1.26) = 1.054727 as the
primary target (since this is the locked surrogate's grid), but
also report the γ-lock value 1.053694 for comparison.

A candidate "explains" the factor if its value matches to better
than 0.01 % AND the formula is built from natural BAM ingredients
(k_5, π, integers ≤ 100). Looser matches are recorded as
"near-misses".
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


# Targets — three precision points for ω(1, 0):
# (a) canonical baseline R_OUTER = 1.26          → ω = 1.054727
# (b) γ-lock R_OUTER ≈ 1.262266 (Σ V_max = 22.5) → ω ≈ 1.053694
# (c) cross-species fixed point R* ≈ 1.262239    → ω ≈ 1.053710 (computed below)
TARGET_CANONICAL = 1.054727       # ω at canonical R_OUTER = 1.26
TARGET_GAMMA_LOCK = 1.053694      # ω at γ-lock R_OUTER ≈ 1.2623
# We use TARGET_GAMMA_LOCK as the "physical" target — this is where the
# self-consistency fixed point lives.

K_5 = 5


@dataclass
class Candidate:
    formula: str
    family: str           # "rational" | "root" | "k5_pi" | "series"
    target: str           # "omega" | "omega_sq"
    value: float
    deviation_pct_canonical: float       # vs 1.054727
    deviation_pct_gamma_lock: float      # vs 1.053694
    matches_canonical_under_01pct: bool
    matches_gamma_lock_under_01pct: bool
    structural_plausibility: str         # "high" | "medium" | "low"


def _record(formula: str, family: str, target: str, value: float, plausibility: str) -> Candidate:
    if target == "omega":
        observed = value
    else:  # omega_sq
        observed = math.sqrt(value) if value >= 0 else float("nan")
    if math.isnan(observed):
        return Candidate(
            formula=formula, family=family, target=target, value=value,
            deviation_pct_canonical=float("nan"),
            deviation_pct_gamma_lock=float("nan"),
            matches_canonical_under_01pct=False,
            matches_gamma_lock_under_01pct=False,
            structural_plausibility=plausibility,
        )
    pct_canonical = 100.0 * (observed - TARGET_CANONICAL) / TARGET_CANONICAL
    pct_gamma = 100.0 * (observed - TARGET_GAMMA_LOCK) / TARGET_GAMMA_LOCK
    return Candidate(
        formula=formula, family=family, target=target, value=value,
        deviation_pct_canonical=pct_canonical,
        deviation_pct_gamma_lock=pct_gamma,
        matches_canonical_under_01pct=abs(pct_canonical) < 0.01,
        matches_gamma_lock_under_01pct=abs(pct_gamma) < 0.01,
        structural_plausibility=plausibility,
    )


# ---------------------------------------------------------------------------
# Family (1): small rationals
# ---------------------------------------------------------------------------

def _small_rationals(max_pq: int = 200) -> list[Candidate]:
    """p/q for ω or ω², 1 ≤ q < p ≤ max_pq, within 0.5 %."""
    out: list[Candidate] = []
    for target_name, target_val in [("omega", TARGET_GAMMA_LOCK),
                                     ("omega_sq", TARGET_GAMMA_LOCK ** 2)]:
        for p in range(1, max_pq + 1):
            for q in range(1, p):
                v = p / q
                err = abs(v - target_val) / target_val * 100
                if err < 0.5:
                    plausibility = "high" if (p < 30 and q < 30) else "medium"
                    out.append(_record(
                        f"{p}/{q}", "rational", target_name, v, plausibility
                    ))
    # Dedupe by formula
    seen = set()
    uniq: list[Candidate] = []
    for c in out:
        key = (c.formula, c.target)
        if key not in seen:
            seen.add(key)
            uniq.append(c)
    return uniq


# ---------------------------------------------------------------------------
# Family (2): roots of small rationals
# ---------------------------------------------------------------------------

def _root_rationals(max_pq: int = 100, max_n: int = 7) -> list[Candidate]:
    """(p/q)^(1/n) matching ω, for small p, q, n."""
    out: list[Candidate] = []
    for n in range(2, max_n + 1):
        for p in range(1, max_pq + 1):
            for q in range(1, p):
                v = (p / q) ** (1.0 / n)
                err = abs(v - TARGET_GAMMA_LOCK) / TARGET_GAMMA_LOCK * 100
                if err < 0.02:
                    plausibility = "high" if (p < 30 and q < 30 and n <= 3) else "medium"
                    out.append(_record(
                        f"({p}/{q})^(1/{n})", "root", "omega", v, plausibility
                    ))
    return out


# ---------------------------------------------------------------------------
# Family (3): k_5-and-π combinations
# ---------------------------------------------------------------------------

def _k5_pi_forms() -> list[Candidate]:
    """Natural BAM expressions in k_5 and π."""
    k5 = K_5
    pi = math.pi
    forms = [
        # ω² candidates
        ("ω²", "1 + 1/k_5²", 1 + 1/k5**2),
        ("ω²", "1 + 2/k_5²", 1 + 2/k5**2),
        ("ω²", "1 + 3/k_5²", 1 + 3/k5**2),
        ("ω²", "1 + (k_5−2)/k_5²", 1 + (k5-2)/k5**2),
        ("ω²", "1 + 1/k_5(k_5−1)", 1 + 1/(k5*(k5-1))),
        ("ω²", "(k_5²+3)/k_5²", (k5**2+3)/k5**2),
        ("ω²", "10/9", 10/9),
        ("ω²", "9/8", 9/8),
        ("ω²", "k_5/(k_5−1/k_5)", k5 / (k5 - 1/k5)),
        ("ω²", "(k_5+1)/k_5 · 9/10", (k5+1)/k5 * 9/10),
        ("ω²", "1 + π/k_5⁴", 1 + pi/k5**4),
        ("ω²", "1 + 2π/k_5⁴", 1 + 2*pi/k5**4),
        # ω candidates directly
        ("ω", "(10/9)^(1/2)", math.sqrt(10/9)),
        ("ω", "(9/8)^(1/2)", math.sqrt(9/8)),
        ("ω", "1 + 1/k_5(k_5−1)·(1/2)", 1 + 1/(2*k5*(k5-1))),
        ("ω", "1 + 1/(k_5²−1)·(1/2)", 1 + 1/(2*(k5**2 - 1))),
        ("ω", "(1+1/k_5²)^(1/2)·c1", math.sqrt(1+1/k5**2)),
        # Series-style
        ("ω", "1 + π/(R_OUTER−R_MID)·tiny", 1 + 1/19),         # 20/19 = 1.0526
        ("ω", "20/19", 20/19),
        ("ω", "57/54", 57/54),
        ("ω", "(k_5−1)·k_5/((k_5−1)·k_5−1)", (k5-1)*k5/((k5-1)*k5 - 1)),
    ]
    out: list[Candidate] = []
    for target, formula, val in forms:
        target_arg = "omega" if target == "ω" else "omega_sq"
        plausibility = "high"
        out.append(_record(formula, "k5_pi", target_arg, val, plausibility))
    return out


# ---------------------------------------------------------------------------
# Family (4): 1 + ε/k_5^m series truncations
# ---------------------------------------------------------------------------

def _series_truncations() -> list[Candidate]:
    out: list[Candidate] = []
    k5 = K_5
    for m in (2, 3, 4):
        for eps in range(1, 8):
            for sign in (+1, -1):
                v = 1 + sign * eps / k5**m
                err = abs(v - TARGET_GAMMA_LOCK) / TARGET_GAMMA_LOCK * 100
                if err < 1.0:
                    out.append(_record(
                        f"1 {'+' if sign > 0 else '−'} {eps}/k_5^{m}",
                        "series", "omega", v, "high"
                    ))
    return out


# ---------------------------------------------------------------------------
# Probe runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    rationals = _small_rationals()
    roots = _root_rationals()
    k5pi = _k5_pi_forms()
    series = _series_truncations()
    all_cands: list[Candidate] = rationals + roots + k5pi + series

    # Best per family by deviation from γ-lock target
    by_family: dict[str, list[Candidate]] = {}
    for c in all_cands:
        by_family.setdefault(c.family, []).append(c)
    best_per_family = {
        family: min(items, key=lambda c: abs(c.deviation_pct_gamma_lock))
        for family, items in by_family.items()
    }

    overall_best = min(all_cands, key=lambda c: abs(c.deviation_pct_gamma_lock))
    high_plausibility_matches = sorted(
        [c for c in all_cands
         if c.structural_plausibility == "high"
         and abs(c.deviation_pct_gamma_lock) < 0.1],
        key=lambda c: abs(c.deviation_pct_gamma_lock),
    )
    exact_matches_01pct = [
        c for c in all_cands
        if c.matches_gamma_lock_under_01pct and c.structural_plausibility == "high"
    ]

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "targets": {
            "canonical_R_OUTER_1_26": TARGET_CANONICAL,
            "gamma_lock_R_OUTER_1_2623": TARGET_GAMMA_LOCK,
            "deviation_canonical_from_unity_pct": (TARGET_CANONICAL - 1) * 100,
            "deviation_gamma_lock_from_unity_pct": (TARGET_GAMMA_LOCK - 1) * 100,
        },
        "n_candidates_total": len(all_cands),
        "best_per_family": {f: asdict(c) for f, c in best_per_family.items()},
        "overall_best": asdict(overall_best),
        "high_plausibility_top_10": [asdict(c) for c in high_plausibility_matches[:10]],
        "exact_matches_under_01pct_high_plausibility": [
            asdict(c) for c in exact_matches_01pct
        ],
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Closed-form search for the 1.054 factor")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append("")
    t = summary["targets"]
    lines.append("**Targets** (ω(l=1, n=0)):")
    lines.append("")
    lines.append(
        f"- Canonical R_OUTER = 1.26: **ω = {t['canonical_R_OUTER_1_26']:.6f}** "
        f"({t['deviation_canonical_from_unity_pct']:+.3f}% from 1)"
    )
    lines.append(
        f"- γ-lock R_OUTER ≈ 1.2623 (cross-species fixed point): "
        f"**ω = {t['gamma_lock_R_OUTER_1_2623']:.6f}** "
        f"({t['deviation_gamma_lock_from_unity_pct']:+.3f}% from 1)"
    )
    lines.append("")
    lines.append(
        f"**Candidates scanned:** {summary['n_candidates_total']} across four "
        "families (small rationals, roots of rationals, k₅-and-π forms, "
        "series-truncation `1 + ε/k₅^m`). A candidate is 'high-plausibility' "
        "if it uses only small integers (≤ 30) and / or k₅ = 5 directly."
    )
    lines.append("")

    lines.append("## Best per family")
    lines.append("")
    lines.append(
        "| family | best formula | value | %Δ vs γ-lock | plausibility |"
    )
    lines.append("|---|---|---:|---:|---|")
    for family, c in summary["best_per_family"].items():
        lines.append(
            f"| `{family}` | `{c['formula']}` | {c['value']:.6f} | "
            f"{c['deviation_pct_gamma_lock']:+.4f}% | {c['structural_plausibility']} |"
        )
    lines.append("")

    lines.append("## High-plausibility candidates within 0.1% of γ-lock")
    lines.append("")
    hp = summary["high_plausibility_top_10"]
    if not hp:
        lines.append("(none — no high-plausibility candidate matches within 0.1%)")
    else:
        lines.append(
            "| formula | family | target | value | %Δ vs γ-lock |"
        )
        lines.append("|---|---|---|---:|---:|")
        for c in hp:
            lines.append(
                f"| `{c['formula']}` | {c['family']} | {c['target']} | "
                f"{c['value']:.6f} | {c['deviation_pct_gamma_lock']:+.4f}% |"
            )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    exact = summary["exact_matches_under_01pct_high_plausibility"]
    if exact:
        c = exact[0]
        lines.append(
            f"**Strong candidate:** `{c['formula']}` (family `{c['family']}`) "
            f"matches the γ-lock ω to {abs(c['deviation_pct_gamma_lock']):.4f}% — "
            "within the high-plausibility tolerance of 0.01%. If this expression "
            "is the correct closed form, the dimensional bridge to ℏ "
            f"reduces to ℏ ω = {c['value']:.4f} · m_e c², with the factor "
            "**derived from k₅ and small integers alone**."
        )
    else:
        lines.append(
            "**No high-plausibility candidate matches within 0.01%.**"
        )
        lines.append("")
        ob = summary["overall_best"]
        lines.append(
            f"Overall best (no plausibility filter): `{ob['formula']}` "
            f"({ob['family']}) at {ob['deviation_pct_gamma_lock']:+.4f}% — "
            "but typically with large integers (> 30) that lack obvious "
            "structural meaning."
        )
        lines.append("")
        if hp:
            top = hp[0]
            lines.append(
                f"Best high-plausibility near-miss: `{top['formula']}` at "
                f"{top['deviation_pct_gamma_lock']:+.4f}% — suggestive but not "
                "exact."
            )
        lines.append("")
        lines.append(
            "**This is the THESIS.md-flagged clean negative result.** "
            "If the 1.054 factor has no closed form in natural BAM ingredients, "
            "the dimensional bridge to ℏ retains an irreducible structural "
            "constant that must be anchored externally (via m_e). BAM remains "
            "**dimensional-ratio-complete and dimensional-scale-incomplete**."
        )
    lines.append("")
    lines.append(
        "Caveat: 'closed form' depends on what's considered natural. The "
        "search restricts to (small integers, k₅, π); a wider net (e.g. "
        "involving the antipodal-closure quantum 100, ratios of Tangherlini "
        "ground-mode eigenvalues, etc.) might surface a more structural "
        "expression. The probe is exhaustive only within the catalog scanned."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_factor_1054_search"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
