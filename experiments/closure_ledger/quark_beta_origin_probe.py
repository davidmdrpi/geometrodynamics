"""
Quark β-origin probe.

The locked quark β = 466 · π/2 corresponds to a closure-quantum integer
4β_q / (2π) = 466. Per `geometrodynamics/qcd/quark_spectrum.py` and
`docs/quark_axioms.md` §8, this integer is empirically robust under
*uniform* mass scaling but drifts by ~90 units under per-species
perturbations — meaning its precise numerical value is not a clean
topological invariant. Per THESIS.md, the research target is to either
derive 466 from a principled enumeration OR record a clean negative
result (no enumeration in the framework's natural list lands at 466).

This probe scans candidate integer enumerations and reports nearest
matches to:

  - N_q = 466                        (the locked quark winding number)
  - ΔN  = N_q − N_l = 466 − 100 = 366   (the lepton↔quark gap)

The candidate enumerations are drawn from the structures the BAM
framework actually exposes:

  S³ angular harmonics      (l+1)² for l = 0, 1, 2, … and partial sums
  S² angular harmonics      (2l+1) for l = 0, 1, 2, … and partial sums
  Tangherlini barrier sums  Σ V_max(l) (already validated for the
                            lepton pinhole, scaled to integer counts)
  S³ angular eigenvalues    l(l+2) for l = 0, 1, 2, … and partial sums
  SU(3) representation dim  fundamental (3), sextet (6), adjoint (8),
                            decuplet (10), 15, 21, 27, 35, … and sums
  Torus-knot crossings      (p−1)(q−1) for small (p, q)
  Closure shell composites  k(k+2) at k_5 = 5 (= 35), N · 35,
                            (k_5² + 1) · m, etc.
  Triangular / tetrahedral  T_n = n(n+1)/2, and Σ T_n = n(n+1)(n+2)/6
  Multi-shell sums          Σ_{k=1,3,5} k² = 35, weighted variants

A candidate "explains" the target if its value is exact (≤ 0.5 unit
tolerance, i.e. rounds to N_q or ΔN). A candidate "approximates" if
within ±5%. We report exact matches (if any), best near-misses by
category, and the joint distribution.

The probe is exploratory — its purpose is to either find the natural
home of 466 OR document the absence of one in this catalog.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


N_QUARK = 466
N_LEPTON = 100
DELTA_N = N_QUARK - N_LEPTON   # 366

EXACT_TOL = 0.5          # ≤ 0.5 unit ⇒ rounds to target ⇒ "exact"
APPROX_PCT = 5.0         # ≤ 5% relative ⇒ "approximate"


@dataclass
class CandidateInteger:
    name: str
    category: str
    formula: str
    value: int
    pct_diff_466: float
    pct_diff_366: float
    diff_abs_466: int
    diff_abs_366: int
    explains_quark: bool
    explains_gap: bool


def _record(name: str, category: str, formula: str, v: int) -> CandidateInteger:
    p466 = 100.0 * (v - N_QUARK) / N_QUARK
    p366 = 100.0 * (v - DELTA_N) / DELTA_N
    return CandidateInteger(
        name=name,
        category=category,
        formula=formula,
        value=int(v),
        pct_diff_466=p466,
        pct_diff_366=p366,
        diff_abs_466=int(abs(v - N_QUARK)),
        diff_abs_366=int(abs(v - DELTA_N)),
        explains_quark=abs(v - N_QUARK) <= EXACT_TOL,
        explains_gap=abs(v - DELTA_N) <= EXACT_TOL,
    )


# ---------------------------------------------------------------------------
# Enumeration generators (each returns a list of (formula_str, value) tuples).
# ---------------------------------------------------------------------------

def _s3_angular_harmonics() -> list[tuple[str, int]]:
    """(l+1)² for l = 0, 1, 2, … plus partial sums."""
    out: list[tuple[str, int]] = []
    deg = [(l + 1) ** 2 for l in range(0, 30)]
    for l, d in enumerate(deg):
        if 50 <= d <= 700:
            out.append((f"S³ harmonic dim (l+1)² at l={l}", d))
    cum = 0
    cum_pairs: list[tuple[int, int]] = []
    for l, d in enumerate(deg):
        cum += d
        cum_pairs.append((l, cum))
        if 50 <= cum <= 700:
            out.append((f"Σ_{{l'=0..{l}}} (l'+1)²", cum))
    # Bracketed sums (l = a..b for a > 0) as a separate flavor.
    for a in (1, 2):
        for b in range(a + 1, 25):
            s = sum((l + 1) ** 2 for l in range(a, b + 1))
            if 50 <= s <= 700:
                out.append((f"Σ_{{l={a}..{b}}} (l+1)²", s))
    return out


def _s2_angular_harmonics() -> list[tuple[str, int]]:
    """(2l+1) for l = 0, 1, 2, … plus partial sums."""
    out: list[tuple[str, int]] = []
    cum = 0
    for l in range(0, 30):
        d = 2 * l + 1
        cum += d
        if 50 <= d <= 700:
            out.append((f"(2l+1) at l={l}", d))
        if 50 <= cum <= 700:
            out.append((f"Σ_{{l'=0..{l}}} (2l'+1) = (l+1)²", cum))
    # (2l+1)² for l = 0..N
    cum2 = 0
    for l in range(0, 25):
        d = (2 * l + 1) ** 2
        cum2 += d
        if 50 <= cum2 <= 700:
            out.append((f"Σ_{{l'=0..{l}}} (2l'+1)²", cum2))
    return out


def _s3_angular_eigenvalues() -> list[tuple[str, int]]:
    """l(l+2) for l = 0, 1, … plus weighted variants."""
    out: list[tuple[str, int]] = []
    for l in range(0, 25):
        v = l * (l + 2)
        if 50 <= v <= 700:
            out.append((f"l(l+2) at l={l}", v))
    cum = 0
    for l in range(0, 30):
        cum += l * (l + 2)
        if 50 <= cum <= 700:
            out.append((f"Σ_{{l'=0..{l}}} l'(l'+2)", cum))
    for k_5 in (5,):
        for n in range(2, 25):
            v = n * k_5 * (k_5 + 2)   # n · 35
            if 50 <= v <= 700:
                out.append((f"{n} · k_5(k_5+2) (k_5={k_5})", v))
            v2 = n * (k_5 ** 2)
            if 50 <= v2 <= 700:
                out.append((f"{n} · k_5² (k_5={k_5})", v2))
            v3 = n * (k_5 ** 2 + 1)   # n · 26
            if 50 <= v3 <= 700:
                out.append((f"{n} · (k_5² + 1)", v3))
    return out


def _su3_representation_dimensions() -> list[tuple[str, int]]:
    """Dimensions of low-lying SU(3) irreps (Dynkin/Young-tableaux)."""
    # SU(3) irrep dimensions in lexicographic order of (p, q) labels.
    irreps: list[tuple[tuple[int, int], int]] = []
    for p in range(0, 12):
        for q in range(0, 12):
            d = (p + 1) * (q + 1) * (p + q + 2) // 2
            irreps.append(((p, q), d))
    out: list[tuple[str, int]] = []
    for (p, q), d in irreps:
        if 50 <= d <= 700:
            out.append((f"SU(3) irrep ({p},{q}) dim", d))
    # Sums over the first N irreps in dimension-ascending order (small reps).
    by_dim = sorted(irreps, key=lambda x: x[1])
    cum = 0
    for i, (lab, d) in enumerate(by_dim[:30]):
        cum += d
        if 50 <= cum <= 700:
            out.append((f"Σ first {i+1} SU(3) irreps (dim-ascending)", cum))
    return out


def _torus_knot_crossings() -> list[tuple[str, int]]:
    """(p−1)(q−1) for small coprime (p, q) and small multiples."""
    out: list[tuple[str, int]] = []
    for p in range(2, 30):
        for q in range(p + 1, 30):
            if math.gcd(p, q) != 1:
                continue
            c = (p - 1) * (q - 1)
            if 50 <= c <= 700:
                out.append((f"T({p},{q}) crossing # = (p−1)(q−1)", c))
            # Genus also: (p−1)(q−1)/2
            g = c // 2 if c % 2 == 0 else None
            if g and 50 <= g <= 700:
                out.append((f"T({p},{q}) genus = (p−1)(q−1)/2", g))
    return out


def _composite_arithmetic() -> list[tuple[str, int]]:
    """Suggestive small-integer arithmetic combinations."""
    out: list[tuple[str, int]] = []
    # k_5² + 1 = 26
    out.append(("(k_5² + 1) · 18 = 26·18", 26 * 18))
    out.append(("(k_5² + 1) · 18 - 2", 26 * 18 - 2))
    # 4 · 91 = 364 (4 generations × Σ_{l=0..5} (l+1)²)
    out.append(("4 · Σ_{l=0..5} (l+1)²  [= 4·91]", 4 * 91))
    out.append(("5 · 91 + 11", 5 * 91 + 11))
    # 100 + (gap)
    for n in (366, 364, 365, 367, 350, 360, 372):
        out.append((f"100 + {n}", 100 + n))
    # k_5 · 93 + 1 = 466
    out.append(("5 · 93 + 1", 5 * 93 + 1))
    out.append(("2 · 233 (233 prime)", 2 * 233))
    # Sum of first N primes
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
              59, 61, 67, 71, 73, 79, 83, 89]
    cum = 0
    for i, p in enumerate(primes):
        cum += p
        if 50 <= cum <= 700:
            out.append((f"Σ first {i+1} primes", cum))
    # k(k+2) at lepton k = 1,3,5: 3, 15, 35; sum 53. weighted variants:
    out.append(("Σ_{k=1,3,5} k(k+2)·(k+2)", 3*3 + 15*5 + 35*7))
    out.append(("Σ_{k=1,3,5} k²·(k+2)", 1*3 + 9*5 + 25*7))
    return out


def _tangherlini_barrier_integers() -> list[tuple[str, int]]:
    """Tangherlini-derived integer counts (Σ V_max scaled, etc.)."""
    out: list[tuple[str, int]] = []
    # The lepton pinhole γ = 22.5 ≈ Σ_{l=0..5} V_max(l). The corresponding
    # closure-quantum count would be 4γ/(2π) ≈ 14.32 (NOT integer). Test
    # whether multiples of small integer ratios land near 466.
    target_lep = 22.5
    target_qrk_pinhole = 22.25
    for n in range(1, 50):
        # n · γ_lepton / π ≈ ?
        v = round(n * target_lep / math.pi)
        if 50 <= v <= 700:
            out.append((f"round({n} · γ_lep / π)", v))
        v2 = round(n * target_qrk_pinhole / math.pi)
        if 50 <= v2 <= 700:
            out.append((f"round({n} · γ_qrk / π)", v2))
    return out


# ---------------------------------------------------------------------------
# Probe runner
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    generators = {
        "S³_harmonics": _s3_angular_harmonics,
        "S²_harmonics": _s2_angular_harmonics,
        "S³_eigenvalues_l(l+2)": _s3_angular_eigenvalues,
        "SU(3)_rep_dimensions": _su3_representation_dimensions,
        "torus_knot_invariants": _torus_knot_crossings,
        "tangherlini_barrier_integers": _tangherlini_barrier_integers,
        "composite_arithmetic": _composite_arithmetic,
    }

    by_cat: dict[str, list[CandidateInteger]] = {}
    all_cands: list[CandidateInteger] = []
    for cat, gen in generators.items():
        for name, val in gen():
            c = _record(name, cat, name, val)
            by_cat.setdefault(cat, []).append(c)
            all_cands.append(c)

    # Best per category by 466.
    best_466_per_cat = {
        cat: min(items, key=lambda c: c.diff_abs_466)
        for cat, items in by_cat.items()
    }
    best_366_per_cat = {
        cat: min(items, key=lambda c: c.diff_abs_366)
        for cat, items in by_cat.items()
    }
    exact_466 = [c for c in all_cands if c.explains_quark]
    exact_366 = [c for c in all_cands if c.explains_gap]
    # `composite_arithmetic` is a sandbox for structural-numerology
    # checks — NOT a principled enumeration family. Filter it out for
    # the verdict-grade match list so we don't credit `100 + 366` etc.
    PRINCIPLED_CATS = {
        "S³_harmonics",
        "S²_harmonics",
        "S³_eigenvalues_l(l+2)",
        "SU(3)_rep_dimensions",
        "torus_knot_invariants",
        "tangherlini_barrier_integers",
    }
    principled_exact_466 = [
        c for c in exact_466 if c.category in PRINCIPLED_CATS
    ]
    principled_exact_366 = [
        c for c in exact_366 if c.category in PRINCIPLED_CATS
    ]
    principled_near_466 = sorted(
        [
            c for c in all_cands
            if c.category in PRINCIPLED_CATS and abs(c.pct_diff_466) <= APPROX_PCT
        ],
        key=lambda c: abs(c.pct_diff_466),
    )
    principled_near_366 = sorted(
        [
            c for c in all_cands
            if c.category in PRINCIPLED_CATS and abs(c.pct_diff_366) <= APPROX_PCT
        ],
        key=lambda c: abs(c.pct_diff_366),
    )
    near_466 = sorted(
        [c for c in all_cands if abs(c.pct_diff_466) <= APPROX_PCT],
        key=lambda c: abs(c.pct_diff_466),
    )
    near_366 = sorted(
        [c for c in all_cands if abs(c.pct_diff_366) <= APPROX_PCT],
        key=lambda c: abs(c.pct_diff_366),
    )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "targets": {
            "N_quark": N_QUARK,
            "N_lepton": N_LEPTON,
            "delta_N_quark_lepton_gap": DELTA_N,
        },
        "n_candidates_total": len(all_cands),
        "best_per_category_for_466": {
            cat: asdict(c) for cat, c in best_466_per_cat.items()
        },
        "best_per_category_for_366": {
            cat: asdict(c) for cat, c in best_366_per_cat.items()
        },
        "exact_matches_466": [asdict(c) for c in exact_466],
        "exact_matches_366": [asdict(c) for c in exact_366],
        "near_matches_466_within_5pct": [asdict(c) for c in near_466],
        "near_matches_366_within_5pct": [asdict(c) for c in near_366],
        "principled_categories": sorted(PRINCIPLED_CATS),
        "principled_exact_matches_466": [asdict(c) for c in principled_exact_466],
        "principled_exact_matches_366": [asdict(c) for c in principled_exact_366],
        "principled_near_matches_466_within_5pct": [
            asdict(c) for c in principled_near_466
        ],
        "principled_near_matches_366_within_5pct": [
            asdict(c) for c in principled_near_366
        ],
        "all_candidates": [asdict(c) for c in all_cands],
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Quark β-origin probe — searching principled enumerations for N = 466")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append(
        f"**Targets:** N_quark = {summary['targets']['N_quark']} "
        f"(closure-quantum count for β_quark = N · π/2), "
        f"ΔN = {summary['targets']['delta_N_quark_lepton_gap']} "
        "(quark/lepton gap)."
    )
    lines.append(
        f"**Candidates scanned:** {summary['n_candidates_total']} integers "
        "across seven principled-enumeration categories (S³/S² harmonics, "
        "S³ angular eigenvalues, SU(3) rep dimensions, torus-knot "
        "invariants, Tangherlini barrier counts, suggestive arithmetic "
        "composites)."
    )
    lines.append("")

    lines.append("## Best per category against N = 466")
    lines.append("")
    lines.append("| category | candidate | value | %Δ | |value−466| |")
    lines.append("|---|---|---:|---:|---:|")
    for cat, c in summary["best_per_category_for_466"].items():
        lines.append(
            f"| `{cat}` | `{c['name']}` | {c['value']} | "
            f"{c['pct_diff_466']:+.3f}% | {c['diff_abs_466']} |"
        )
    lines.append("")

    lines.append("## Best per category against ΔN = 366")
    lines.append("")
    lines.append("| category | candidate | value | %Δ | |value−366| |")
    lines.append("|---|---|---:|---:|---:|")
    for cat, c in summary["best_per_category_for_366"].items():
        lines.append(
            f"| `{cat}` | `{c['name']}` | {c['value']} | "
            f"{c['pct_diff_366']:+.3f}% | {c['diff_abs_366']} |"
        )
    lines.append("")

    lines.append("## Exact integer matches")
    lines.append("")
    if summary["exact_matches_466"]:
        lines.append("**N = 466 hits:**")
        lines.append("")
        for c in summary["exact_matches_466"]:
            lines.append(f"- `{c['name']}` (category: `{c['category']}`)")
        lines.append("")
    else:
        lines.append("- No exact integer enumeration in the catalog "
                     "produces 466.")
        lines.append("")
    if summary["exact_matches_366"]:
        lines.append("**ΔN = 366 hits:**")
        lines.append("")
        for c in summary["exact_matches_366"]:
            lines.append(f"- `{c['name']}` (category: `{c['category']}`)")
        lines.append("")
    else:
        lines.append("- No exact integer enumeration in the catalog "
                     "produces 366.")
        lines.append("")

    lines.append("## Near matches within ±5%")
    lines.append("")
    lines.append("**N = 466 (within ±5%):**")
    lines.append("")
    if not summary["near_matches_466_within_5pct"]:
        lines.append("(none)")
    else:
        lines.append("| candidate | category | value | %Δ vs 466 |")
        lines.append("|---|---|---:|---:|")
        for c in summary["near_matches_466_within_5pct"][:20]:
            lines.append(
                f"| `{c['name']}` | {c['category']} | {c['value']} | "
                f"{c['pct_diff_466']:+.3f}% |"
            )
    lines.append("")
    lines.append("**ΔN = 366 (within ±5%):**")
    lines.append("")
    if not summary["near_matches_366_within_5pct"]:
        lines.append("(none)")
    else:
        lines.append("| candidate | category | value | %Δ vs 366 |")
        lines.append("|---|---|---:|---:|")
        for c in summary["near_matches_366_within_5pct"][:20]:
            lines.append(
                f"| `{c['name']}` | {c['category']} | {c['value']} | "
                f"{c['pct_diff_366']:+.3f}% |"
            )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    lines.append(
        "Principled-enumeration categories: "
        f"{', '.join(f'`{c}`' for c in summary['principled_categories'])}. "
        "The `composite_arithmetic` category is a sandbox for "
        "structural-numerology checks (`100 + 366`, `5·93 + 1`, etc.) "
        "and is excluded from the verdict — its exact hits are "
        "hand-constructed identities, not derived enumerations."
    )
    lines.append("")
    p_exact_466 = summary["principled_exact_matches_466"]
    p_exact_366 = summary["principled_exact_matches_366"]
    p_near_466 = summary["principled_near_matches_466_within_5pct"]
    p_near_366 = summary["principled_near_matches_366_within_5pct"]
    if p_exact_466:
        winner = p_exact_466[0]
        lines.append(
            f"**N = 466 is reproduced exactly by `{winner['name']}` "
            f"(category `{winner['category']}`).** Candidate topological "
            "reading for the locked quark β."
        )
    elif p_exact_366:
        winner = p_exact_366[0]
        lines.append(
            f"**ΔN = 366 (the quark/lepton gap) is reproduced exactly "
            f"by `{winner['name']}` (category `{winner['category']}`).** "
            "Suggests N_q = 100 + 366 = lepton closure quantum + "
            f"sector-specific contribution from `{winner['name']}`."
        )
    else:
        lines.append(
            "**No principled enumeration in the catalog produces N = 466 "
            "or ΔN = 366 exactly.** This is the clean negative result "
            "THESIS.md flagged as a valid research outcome."
        )
        lines.append("")
        if p_near_466 or p_near_366:
            lines.append(
                "Closest principled near-misses (within ±5%):"
            )
            lines.append("")
            if p_near_466:
                top = p_near_466[:5]
                lines.append("**For N = 466:**")
                lines.append("")
                for c in top:
                    lines.append(
                        f"- `{c['name']}` = {c['value']} "
                        f"({c['pct_diff_466']:+.3f}%, off by "
                        f"{c['diff_abs_466']})"
                    )
                lines.append("")
            if p_near_366:
                top = p_near_366[:5]
                lines.append("**For ΔN = 366:**")
                lines.append("")
                for c in top:
                    lines.append(
                        f"- `{c['name']}` = {c['value']} "
                        f"({c['pct_diff_366']:+.3f}%, off by "
                        f"{c['diff_abs_366']})"
                    )
                lines.append("")

    # Highlight the (k_5² + 1) = 26 multiplet pattern if present.
    k5sq1_466 = [
        c for c in p_near_466 if "k_5² + 1" in c["name"] or "k_5²+1" in c["name"]
    ]
    k5sq1_366 = [
        c for c in p_near_366 if "k_5² + 1" in c["name"] or "k_5²+1" in c["name"]
    ]
    if k5sq1_466 and k5sq1_366:
        lines.append("### Suggestive structural pattern: multiples of (k_5² + 1) = 26")
        lines.append("")
        lines.append(
            "Both targets are bracketed by clean multiples of "
            "`(k_5² + 1) = 26`:"
        )
        lines.append("")
        for c in k5sq1_466[:3]:
            lines.append(
                f"- `{c['name']}` = {c['value']}, off by "
                f"{c['diff_abs_466']} from N = 466"
            )
        for c in k5sq1_366[:3]:
            lines.append(
                f"- `{c['name']}` = {c['value']}, off by "
                f"{c['diff_abs_366']} from ΔN = 366"
            )
        lines.append("")
        lines.append(
            "If `(k_5² + 1)` is the natural quark-sector closure-quantum "
            "increment, both 466 and 366 sit within 2 of integer "
            "multiples (18·26 = 468, 14·26 = 364). The unexplained "
            "±2 is small enough to plausibly come from a parity- or "
            "boundary-correction term, but is not derived in this probe. "
            "Whether `(k_5² + 1) · m ± 2` is the correct structural form "
            "is the next concrete sub-target."
        )
        lines.append("")

    lines.append("### Caveat from `docs/quark_axioms.md` §8")
    lines.append("")
    lines.append(
        "Per the existing N-stability analysis, N_q = 466 itself drifts "
        "by ~90 units under per-species mass perturbations while staying "
        "invariant under uniform mass scaling. A topological reading "
        "must not only land on 466 numerically but also reproduce that "
        "robustness signature. None of the candidates in this catalog "
        "has been tested for that — only for raw integer agreement."
    )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_quark_beta_origin_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
