"""
Quark β: sub-block stability sub-probe.

The §8 robustness audit established that no logged ablation preserves
the +1 boundary correction across the full 12-point table — only 33%
of runs sit at δ = +1, with the rest distributed across {-1, 0, +2}.
This sub-probe asks the next concrete question: **does any sub-block
of the baseline decomposition

    N_q  =  (k_5−1)·k_5²  +  2·k_5²(k_5+2)  +  N_c·k_5  +  (N_c − 2)
         =  100           +  350            +  15        +  1
         =  466

remain perturbation-stable?**

Strategy. For each candidate sub-block S drawn from the baseline
decomposition (and several composite sub-blocks), compute the residual
R = N_q − S across all 12 logged ablations. Report:

  - the residual range and width
  - the residual mean and standard deviation
  - whether R is structurally cleaner than the bare N_q (smaller spread
    or tighter modular structure)

A sub-block is "differentially stable" if subtracting it tightens the
residual envelope — i.e. if the variance of R across ablations is
strictly smaller than the variance of N_q. (Constant-shift subtraction
preserves variance, so the test is non-trivial only if S itself
carries some of the drift.)

Independently, the probe scans modular-arithmetic invariants (does
N_q satisfy N_q ≡ 0 mod m for some m across all logged ablations?).
A modular invariant that holds across all 12 ablations is a genuine
structural lock — even if it's coarser than a specific value.
"""

from __future__ import annotations

import json
import math
import statistics
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


K_5 = 5
N_C_QUARK = 3
N_QUARK_BASELINE = 466

# §8 logged N-stability ablation values, copied verbatim from
# `quark_beta_robustness_audit.py`.
ABLATION_N_VALUES = [
    ("baseline (anchor=d, PDG, min_eig)", 466),
    ("PDG × 1.10 (uniform scale)",        466),
    ("PDG × 0.90 (uniform scale)",        466),
    ("anchor = s",                        476),
    ("anchor = c",                        474),
    ("anchor = b",                        474),
    ("anchor = t",                        482),
    ("c × 1.10",                          432),
    ("b × 1.10",                          494),
    ("t × 1.10",                          494),
    ("t × 0.90",                          440),
    ("all ±5% (deterministic)",           510),
]


@dataclass
class SubBlockTest:
    name: str
    formula: str
    baseline_value: int
    residuals: list[int]
    residual_range: tuple[int, int]
    residual_width: int
    residual_mean: float
    residual_stdev: float
    differentially_stable: bool   # variance lower than bare-N_q variance?


@dataclass
class ModularInvariantTest:
    name: str
    modulus: int
    formula: str
    baseline_residue: int
    residues_observed: list[int]
    is_invariant_across_all_ablations: bool
    structural_interpretation: str


# Candidate sub-blocks built from the baseline decomposition.
LEPTON_BLOCK = (K_5 - 1) * K_5 ** 2          # 100
SHELL_BLOCK = 2 * K_5 ** 2 * (K_5 + 2)        # 350
COLOR_K5 = N_C_QUARK * K_5                     # 15
BOUNDARY = N_C_QUARK - 2                       # 1


def _build_subblock_candidates() -> list[tuple[str, str, int]]:
    return [
        ("S_zero",                        "0 (no subtraction)",                            0),
        ("S_lepton",                      "(k_5−1)·k_5²",                                  LEPTON_BLOCK),
        ("S_lepton_plus_boundary",        "(k_5−1)·k_5² + (N_c−2)",                        LEPTON_BLOCK + BOUNDARY),
        ("S_shell",                       "2·k_5²(k_5+2)",                                 SHELL_BLOCK),
        ("S_color_x_k5",                  "N_c·k_5",                                       COLOR_K5),
        ("S_color_plus_boundary",         "N_c·k_5 + (N_c−2)",                             COLOR_K5 + BOUNDARY),
        ("S_lepton_plus_shell",           "(k_5−1)·k_5² + 2·k_5²(k_5+2)",                  LEPTON_BLOCK + SHELL_BLOCK),
        ("S_lepton_plus_color",           "(k_5−1)·k_5² + N_c·k_5",                        LEPTON_BLOCK + COLOR_K5),
        ("S_shell_plus_color",            "2·k_5²(k_5+2) + N_c·k_5",                       SHELL_BLOCK + COLOR_K5),
        ("S_all_but_boundary",            "(k_5−1)·k_5² + 2·k_5²(k_5+2) + N_c·k_5",        LEPTON_BLOCK + SHELL_BLOCK + COLOR_K5),
        ("S_full_baseline",               "(k_5−1)·k_5² + 2·k_5²(k_5+2) + N_c·k_5 + (N_c−2)", LEPTON_BLOCK + SHELL_BLOCK + COLOR_K5 + BOUNDARY),
    ]


def _evaluate_subblock(name: str, formula: str, S: int, n_values: list[int]) -> SubBlockTest:
    residuals = [n - S for n in n_values]
    width = max(residuals) - min(residuals)
    return SubBlockTest(
        name=name,
        formula=formula,
        baseline_value=S,
        residuals=residuals,
        residual_range=(min(residuals), max(residuals)),
        residual_width=width,
        residual_mean=statistics.mean(residuals),
        residual_stdev=statistics.stdev(residuals) if len(residuals) > 1 else 0.0,
        differentially_stable=False,   # filled in after computing baseline width
    )


def _scan_modular_invariants(n_values: list[int]) -> list[ModularInvariantTest]:
    """Test whether N_q ≡ const (mod m) across all ablations for various m."""
    interpretations = {
        2:  "Z₂ partition class — two partition sectors (p = ±) contribute equally to β",
        3:  "SU(3) color triplet (N_c = 3)",
        4:  "spinor double cover or 4-fold winding",
        5:  "k_5 modulus (heaviest shell label)",
        6:  "(k_5 + 1) modulus — closure-quantum on shell+1",
        10: "2 · k_5 (Z₂ × k_5 combined)",
        15: "N_c · k_5",
        25: "k_5²",
        26: "k_5² + 1",
        100:"lepton closure quantum",
    }
    out: list[ModularInvariantTest] = []
    for m, interp in interpretations.items():
        residues = [n % m for n in n_values]
        baseline_res = n_values[0] % m
        invariant = all(r == baseline_res for r in residues)
        out.append(ModularInvariantTest(
            name=f"N_q mod {m}",
            modulus=m,
            formula=f"N_q mod {m}",
            baseline_residue=baseline_res,
            residues_observed=residues,
            is_invariant_across_all_ablations=invariant,
            structural_interpretation=interp,
        ))
    return out


def run_probe() -> dict:
    n_values = [n for _, n in ABLATION_N_VALUES]

    # Sub-block subtractions (sanity check: width should be invariant
    # under constant subtraction; this serves as the null hypothesis).
    sub_blocks = _build_subblock_candidates()
    sub_block_tests = [
        _evaluate_subblock(name, formula, S, n_values)
        for name, formula, S in sub_blocks
    ]
    bare_width = next(
        t.residual_width for t in sub_block_tests if t.name == "S_zero"
    )
    for t in sub_block_tests:
        # Constant-shift subtraction can never reduce variance below bare.
        # So differentially_stable is always False for this family. We
        # report the test for completeness — the result should be
        # identical width across all rows, confirming the math.
        t.differentially_stable = t.residual_width < bare_width

    # Modular-arithmetic invariants.
    modular_tests = _scan_modular_invariants(n_values)
    invariants = [t for t in modular_tests if t.is_invariant_across_all_ablations]

    # Half-N_q diagnostic: if N is always even, what's the per-partition count?
    n_half_values = [n // 2 for n in n_values]
    half_baseline = N_QUARK_BASELINE // 2

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "k_5": K_5,
        "N_c_quark": N_C_QUARK,
        "baseline_decomposition": {
            "lepton_block_(k_5-1)·k_5²": LEPTON_BLOCK,
            "shell_block_2·k_5²(k_5+2)": SHELL_BLOCK,
            "color_x_k_5": COLOR_K5,
            "boundary_N_c-2": BOUNDARY,
            "total": LEPTON_BLOCK + SHELL_BLOCK + COLOR_K5 + BOUNDARY,
        },
        "ablation_N_values": dict(ABLATION_N_VALUES),
        "subblock_tests": [asdict(t) for t in sub_block_tests],
        "modular_invariant_tests": [asdict(t) for t in modular_tests],
        "preserved_modular_invariants": [
            asdict(t) for t in invariants
        ],
        "half_partition_diagnostic": {
            "all_N_even": all(n % 2 == 0 for n in n_values),
            "half_baseline": half_baseline,
            "half_values": n_half_values,
            "half_range": (min(n_half_values), max(n_half_values)),
            "half_drift_pm_baseline": max(
                abs(max(n_half_values) - half_baseline),
                abs(half_baseline - min(n_half_values)),
            ),
        },
    }


def render_markdown(summary: dict) -> str:
    lines: list[str] = []
    lines.append("# Quark β: sub-block stability sub-probe")
    lines.append("")
    lines.append(f"**Run:** {summary['timestamp_utc']}")
    lines.append(
        f"**k_5 = {summary['k_5']}**, **N_c = {summary['N_c_quark']}**, "
        f"**baseline N_q = {summary['baseline_decomposition']['total']}** "
        f"= {summary['baseline_decomposition']['lepton_block_(k_5-1)·k_5²']} "
        f"+ {summary['baseline_decomposition']['shell_block_2·k_5²(k_5+2)']} "
        f"+ {summary['baseline_decomposition']['color_x_k_5']} "
        f"+ {summary['baseline_decomposition']['boundary_N_c-2']}."
    )
    lines.append("")
    lines.append(
        "Asks the focused question: **does any sub-block of the baseline "
        "decomposition remain perturbation-stable across the §8 ablation "
        "table?** Two test families are run: (a) constant-shift "
        "subtraction (sanity-check that subtracting a baseline-locked "
        "constant cannot reduce variance — confirmed below), and (b) "
        "modular-arithmetic invariants (which IS a non-trivial test, "
        "since residue-class membership can survive ablations even when "
        "the absolute value drifts)."
    )
    lines.append("")

    # Section 1: sub-block subtraction.
    lines.append("## Constant-shift subtraction (null-hypothesis check)")
    lines.append("")
    lines.append(
        "Each row subtracts a baseline sub-block S from N_q across "
        "every ablation and reports the residual envelope. Constant-"
        "shift subtraction cannot reduce variance (Var(N − c) = Var(N)), "
        "so this is a sanity check that no sub-block is differentially "
        "stable in the trivial sense."
    )
    lines.append("")
    lines.append(
        "| sub-block S | formula | S value | residual range | width | std dev |"
    )
    lines.append("|---|---|---:|---|---:|---:|")
    for t in summary["subblock_tests"]:
        rng = t["residual_range"]
        lines.append(
            f"| `{t['name']}` | `{t['formula']}` | {t['baseline_value']} | "
            f"[{rng[0]}, {rng[1]}] | {t['residual_width']} | "
            f"{t['residual_stdev']:.2f} |"
        )
    lines.append("")
    bare = next(t for t in summary["subblock_tests"] if t["name"] == "S_zero")
    same_width_count = sum(
        1 for t in summary["subblock_tests"]
        if t["residual_width"] == bare["residual_width"]
    )
    lines.append(
        f"As expected: **all {same_width_count} sub-block subtractions "
        f"give the same drift width ({bare['residual_width']} units)**. "
        "No sub-block is differentially stable in the constant-shift "
        "sense — they all carry the same drift envelope. This rules out "
        "the simple reading 'subtract the lepton sub-block to expose a "
        "tighter compensator.'"
    )
    lines.append("")

    # Section 2: modular invariants.
    lines.append("## Modular-arithmetic invariants")
    lines.append("")
    lines.append(
        "Tests whether N_q ≡ baseline_residue (mod m) across all 12 "
        "logged ablations. Each row is one modulus m, with the structural "
        "interpretation of m. A modular invariance that holds across "
        "every ablation is a genuine sub-block stability — even if the "
        "absolute value of N_q wanders."
    )
    lines.append("")
    lines.append(
        "| modulus m | residues across 12 ablations | invariant? | "
        "interpretation |"
    )
    lines.append("|---:|---|:---:|---|")
    for t in summary["modular_invariant_tests"]:
        marker = "**✓**" if t["is_invariant_across_all_ablations"] else "—"
        residues_str = ", ".join(str(r) for r in t["residues_observed"])
        lines.append(
            f"| {t['modulus']} | [{residues_str}] | {marker} | "
            f"{t['structural_interpretation']} |"
        )
    lines.append("")

    invs = summary["preserved_modular_invariants"]
    if invs:
        lines.append("### Preserved modular invariants")
        lines.append("")
        for t in invs:
            lines.append(
                f"- **N_q ≡ {t['baseline_residue']} (mod {t['modulus']})** "
                f"holds across all 12 logged ablations. Interpretation: "
                f"{t['structural_interpretation']}."
            )
        lines.append("")
    else:
        lines.append("No modular invariant holds across all ablations.")
        lines.append("")

    # Section 3: half-partition diagnostic.
    half = summary["half_partition_diagnostic"]
    lines.append("## Half-partition diagnostic")
    lines.append("")
    if half["all_N_even"]:
        lines.append(
            f"**N_q is even across all 12 logged §8 ablations.** "
            f"Per-partition closure-quantum count n_part = N_q / 2:"
        )
        lines.append("")
        lines.append(
            f"- baseline n_part = {half['half_baseline']}"
        )
        lines.append(
            f"- §8 envelope n_part ∈ {list(half['half_range'])}, "
            f"drift ±{half['half_drift_pm_baseline']} units"
        )
        lines.append("")
        lines.append(
            "The factor-of-2 reflects the v3 ansatz's two-partition "
            "Hamiltonian basis {(k, +), (k, −)}: each Z₂ partition "
            "sector contributes equally to the total β, forcing N_q to "
            "be twice an integer count. The two-ness IS the topological "
            "invariant; the per-partition count itself drifts under "
            "perturbations as the fit compensator."
        )
    else:
        lines.append(
            "N_q is NOT always even — the partition-symmetry reading "
            "needs more nuance."
        )
    lines.append("")

    # Verdict.
    lines.append("## Verdict")
    lines.append("")
    if invs:
        invariant_names = ", ".join(
            f"mod {t['modulus']}" for t in invs
        )
        lines.append(
            f"**The only perturbation-stable structural property of the "
            f"baseline decomposition is parity** "
            f"(`N_q ≡ 0 (mod 2)`, i.e. `{invariant_names}` invariance). "
            "All other modular structures and all constant-shift "
            "sub-block subtractions vary across §8 ablations."
        )
        lines.append("")
        lines.append(
            "Structurally: the two-Z₂-partition multiplicity is the "
            "topological invariant; the specific N_q value within the "
            "even integers is the compensator. This refines the "
            "decomposition probe's reading from"
        )
        lines.append("")
        lines.append(
            "  `N_q = ((k_5−1)·k_5 + 2·k_5(k_5+2) + N_c)·k_5 + (N_c−2)`"
        )
        lines.append("")
        lines.append("to the more honest")
        lines.append("")
        lines.append(
            "  `N_q = 2 · n_part`,   where `n_part` is a phenomenological "
            "per-partition closure-quantum count whose specific value "
            "(233 at baseline) is NOT topologically locked."
        )
        lines.append("")
        lines.append(
            "The 2× partition factor is the sole structural piece; the "
            "rest is fit. This is the cleanest summary of the §8 "
            "compensator behavior under a sub-block-stability lens."
        )
    else:
        lines.append(
            "**No perturbation-stable sub-block found.** The §8 ablation "
            "table moves N_q across {-1, 0, +1, +2} mod 5 and across "
            "an envelope of [432, 510] in absolute value, with no "
            "modular structure preserved across all ablations. The "
            "structural decomposition of the baseline value 466 is a "
            "post-hoc reading; β_quark is fully phenomenological."
        )
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    print(render_markdown(summary))
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_quark_beta_subblock_stability"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "probe.md").write_text(render_markdown(summary))
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
