"""
Experiment runner: orchestrates Layer 1 (ledger) and Layer 2 (blocker
report), produces a single ExperimentResult that can be serialized to
JSON and rendered to a human-readable markdown summary.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from experiments.closure_ledger.blockers import (
    SkBridgeBlocker,
    layer2_blocker_report,
)
from experiments.closure_ledger.ledger import (
    LedgerRow,
    RepoConstants,
    compute_lepton_ledger,
    compute_quark_sector_summary,
)
from experiments.closure_ledger.sk_bridge import DEFAULT_SK_CANDIDATE


TAU = 2.0 * math.pi


@dataclass
class ExperimentResult:
    timestamp_utc: str
    experiment_name: str
    transport_convention: str          # "T2_sign_flip" | "T1_diagnostic"
    chi: float                          # Hopf fibre angle used
    sk_candidate: str                   # S(k) bridge candidate ("none" = Layer 1)
    constants: dict                     # RepoConstants as dict
    import_errors: list[str]            # any failed repo imports
    rows: list[dict]                    # each lepton ledger row
    universality_check: dict            # lepton-sector closure mod 2π summary
    quark_sector: dict                  # quark β-lock structural summary
    sk_bridge_blocker: dict             # Layer 2 blocker report
    overall_status: str                 # one-line conclusion

    def to_dict(self) -> dict:
        return asdict(self)

    def write_json(self, path: Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str))

    def write_summary_markdown(self, path: Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self._render_markdown())

    def _render_markdown(self) -> str:
        lines: list[str] = []
        lines.append("# Closure-phase ledger — run summary")
        lines.append("")
        lines.append(f"**Run:** {self.timestamp_utc}")
        lines.append(f"**Experiment:** {self.experiment_name}")
        lines.append(f"**Transport convention:** `{self.transport_convention}`")
        lines.append(f"**χ (Hopf fibre):** {self.chi}")
        lines.append(f"**S(k) candidate:** `{self.sk_candidate}`")
        lines.append("")
        lines.append(f"**Overall status:** {self.overall_status}")
        lines.append("")

        # Constants section
        lines.append("## Constants used")
        lines.append("")
        lines.append("| name | value | source |")
        lines.append("|---|---|---|")
        c = self.constants
        sources = c.get("source", {})
        lines.append(f"| `action_base` | {c['action_base']:.6f} "
                     f"({c['action_base'] / math.pi:.3f}π) | "
                     f"{sources.get('action_base', 'unknown')} |")
        lines.append(f"| `beta_lepton` | {c['beta_lepton']:.6f} "
                     f"({c['beta_lepton'] / math.pi:.3f}π) | "
                     f"{sources.get('beta_lepton', 'unknown')} |")
        lines.append(f"| `beta_quark` | {c['beta_quark']:.6f} "
                     f"({c['beta_quark'] / math.pi:.3f}π) | "
                     f"{sources.get('beta_quark', 'unknown')} |")
        lines.append(f"| `lepton_quanta` (4β/2π) | {c['lepton_quanta']} | derived |")
        lines.append(f"| `quark_quanta` (4β/2π) | {c['quark_quanta']} | derived |")
        lines.append("")

        # Import errors, if any
        if self.import_errors:
            lines.append("### Import errors (fell back to README values)")
            lines.append("")
            for err in self.import_errors:
                lines.append(f"- `{err}`")
            lines.append("")

        # Per-lepton ledger
        lines.append("## Per-lepton ledger")
        lines.append("")
        lines.append("Phase contributions in units of π. `mod 2π` column "
                     "shows the available-term sum reduced to [0, 2π).")
        lines.append("")
        lines.append("| lepton | k | antipodal | Hopf | throat | uplift | "
                     "radial | moving | total | mod 2π | status |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for row in self.rows:
            term_by_name = {t["name"]: t for t in row["terms"]}

            def fmt(name: str) -> str:
                t = term_by_name.get(name)
                if t is None or t["value"] is None:
                    return "—"
                return f"{t['value'] / math.pi:.3f}π"

            total_pi = (
                f"{row['available_total'] / math.pi:.3f}π"
                if row["available_total"] is not None else "—"
            )
            mod_pi = (
                f"{row['available_total_mod_2pi'] / math.pi:.6f}π"
                if row["available_total_mod_2pi"] is not None else "—"
            )
            lines.append(
                f"| {row['label']} | {row['k']} | "
                f"{fmt('antipodal_closure')} | "
                f"{fmt('hopf_holonomy')} | "
                f"{fmt(f'throat_transport_T{2 if self.transport_convention.startswith(chr(84) + chr(50)) else 1}')} | "
                f"{fmt('lepton_uplift')} | "
                f"{fmt('radial_bulk_phase')} | "
                f"{fmt('moving_mouth_phase')} | "
                f"{total_pi} | {mod_pi} | {row['closure_status']} |"
            )
        lines.append("")

        # Per-mode radial breakdown when the S(k) bridge is active.
        if self.sk_candidate != "none":
            lines.append("## Radial bulk channel — per-mode breakdown")
            lines.append("")
            lines.append(
                "| lepton | k | (l, n) | ω(l, n) | Φ(l, n) | "
                "weight | weight·Φ | N_tp | Δ_maslov | status |"
            )
            lines.append("|---|---|---|---|---|---|---|---|---|---|")
            for row in self.rows:
                detail = row.get("radial_detail")
                if not detail or not detail.get("modes"):
                    continue
                for m in detail["modes"]:
                    omega = (
                        f"{m['omega']:.6f}" if m["omega"] is not None else "—"
                    )
                    phi = (
                        f"{m['phi'] / math.pi:.6f}π"
                        if m["phi"] is not None else "—"
                    )
                    weight = m.get("weight", 1.0)
                    weighted_phi = (
                        f"{(weight * m['phi']) / math.pi:.6f}π"
                        if m["phi"] is not None else "—"
                    )
                    n_tp = m.get("n_turning_points", 0)
                    delta = m.get("maslov_correction", 0.0)
                    delta_str = (
                        f"{delta / math.pi:+.6f}π"
                        if delta else "0"
                    )
                    lines.append(
                        f"| {row['label']} | {row['k']} | "
                        f"(l={m['l']}, n={m['n']}) | {omega} | {phi} | "
                        f"{weight:.6f} | {weighted_phi} | "
                        f"{n_tp} | {delta_str} | "
                        f"{m['status']} |"
                    )
            lines.append("")

        # Universality
        u = self.universality_check
        layer_label = (
            "Layer 1" if self.sk_candidate == "none" else "Layer 2"
        )
        lines.append(f"## Universality check ({layer_label})")
        lines.append("")
        lines.append(
            "- Per-lepton totals mod 2π (in units of π): "
            f"{[f'{v / math.pi:.6f}' for v in u['per_lepton_mod_2pi']]}"
        )
        lines.append(f"- Spread across leptons: {u['spread']:.2e}")
        lines.append(f"- Universal mod 2π (within 1e-9): "
                     f"**{u['universal']}**")
        lines.append(f"- Universal value: "
                     f"{u['universal_value'] / math.pi:.6f}π" if u["universal"]
                     else f"- Universal value: N/A (not universal)")
        lines.append("")

        # Quark sector
        q = self.quark_sector
        lines.append("## Quark sector (structural)")
        lines.append("")
        lines.append(f"- Lepton lock quanta: {q['lepton_lock_quanta']}")
        lines.append(f"- Quark lock quanta: {q['quark_lock_quanta']}")
        lines.append(f"- Lock quanta gap: {q['lock_quanta_gap']}")
        lines.append("")
        lines.append(f"_{q['interpretation']}_")
        lines.append("")

        # Layer 2 blocker
        b = self.sk_bridge_blocker
        lines.append("## Layer 2 blocker — S(k) bridge")
        lines.append("")
        lines.append(f"**Verdict:** {b['verdict']}")
        lines.append("")
        lines.append("### Evidence")
        for e in b["evidence"]:
            lines.append(f"- {e}")
        lines.append("")
        lines.append("### Candidate S(k) maps")
        lines.append("")
        for cand in b["candidates"]:
            lines.append(f"**{cand['name']}**")
            lines.append("")
            lines.append(f"- Formula: `{cand['formula']}`")
            lines.append(f"- Physical picture: {cand['physical_picture']}")
            lines.append(f"- Advantages: {cand['advantages']}")
            lines.append(f"- Open questions: {cand['open_questions']}")
            lines.append("")
        lines.append("### Next steps")
        for ns in b["next_steps"]:
            lines.append(f"- {ns}")
        lines.append("")
        lines.append("### Downgraded predictions")
        for dp in b["downgraded_predictions"]:
            lines.append(f"- {dp}")
        lines.append("")

        return "\n".join(lines)


def _build_universality_check(
    rows: list[LedgerRow],
) -> dict:
    """Summarize whether the available-term sums agree mod 2π across leptons."""
    mods = [
        r.available_total_mod_2pi
        for r in rows
        if r.available_total_mod_2pi is not None
    ]
    if len(mods) != len(rows) or not mods:
        return {
            "per_lepton_mod_2pi": [],
            "spread": float("nan"),
            "universal": False,
            "universal_value": None,
        }
    spread = max(mods) - min(mods)
    # Spread might also be ~2π if values straddle the wrap point.
    universal = spread < 1e-9 or abs(spread - TAU) < 1e-9
    universal_value = mods[0] if universal else None
    return {
        "per_lepton_mod_2pi": mods,
        "spread": spread,
        "universal": universal,
        "universal_value": universal_value,
    }


def _overall_status(
    universality: dict,
    blocker: SkBridgeBlocker,
    import_errors: list[str],
    sk_candidate: str,
    any_radial_missing: bool,
) -> str:
    parts: list[str] = []
    if sk_candidate == "none":
        if universality["universal"]:
            parts.append(
                f"Layer 1 PASS: lepton ledger universal mod 2π at "
                f"{universality['universal_value'] / math.pi:.3f}π."
            )
        else:
            parts.append(
                "Layer 1 FAIL: lepton ledger does NOT close universally."
            )
        parts.append("Layer 2 DISABLED (sk_candidate='none').")
    else:
        if any_radial_missing:
            parts.append(
                f"Layer 2 PARTIAL: candidate '{sk_candidate}' wired but at "
                f"least one row's radial channel did not resolve."
            )
        elif universality["universal"]:
            parts.append(
                f"Layer 2 PASS: candidate '{sk_candidate}' yields universal "
                f"closure mod 2π at {universality['universal_value'] / math.pi:.3f}π."
            )
        else:
            parts.append(
                f"Layer 2 FALSIFIES candidate '{sk_candidate}': lepton "
                f"ledger does NOT close universally mod 2π "
                f"(spread {universality['spread']:.3e})."
            )
    if import_errors:
        parts.append(
            f"({len(import_errors)} repo import(s) fell back to README values.)"
        )
    return " ".join(parts)


def run_experiment(
    chi: float = 0.0,
    transport_power: int = 2,
    sk_candidate: str = DEFAULT_SK_CANDIDATE,
) -> ExperimentResult:
    """
    Run the closure-phase ledger experiment and return a complete
    ExperimentResult. No I/O — call .write_json(...) and
    .write_summary_markdown(...) on the result to persist.

    sk_candidate selects the S(k) → {(l, n)} bridge. Defaults to the
    wired candidate (`A_lowest_radial_per_l`); pass "none" to reproduce
    the Layer-1 ledger.
    """
    rows, constants, import_errors = compute_lepton_ledger(
        chi=chi,
        transport_power=transport_power,
        sk_candidate=sk_candidate,
    )

    universality = _build_universality_check(rows)
    quark_sector = compute_quark_sector_summary(constants)
    blocker = layer2_blocker_report(implemented_candidate=sk_candidate)

    convention = (
        "T2_sign_flip" if transport_power == 2
        else "T1_diagnostic"
    )

    any_radial_missing = any(
        any(
            t.name == "radial_bulk_phase" and t.status != "available"
            for t in r.terms
        )
        for r in rows
    )

    experiment_name = (
        "closure_ledger.layer1" if sk_candidate == "none"
        else "closure_ledger.layer2"
    )

    return ExperimentResult(
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        experiment_name=experiment_name,
        transport_convention=convention,
        chi=chi,
        sk_candidate=sk_candidate,
        constants=asdict(constants),
        import_errors=import_errors,
        rows=[r.to_dict() for r in rows],
        universality_check=universality,
        quark_sector=quark_sector,
        sk_bridge_blocker=blocker.to_dict(),
        overall_status=_overall_status(
            universality, blocker, import_errors,
            sk_candidate, any_radial_missing,
        ),
    )


def _classify_universality(uc: dict, sk_candidate: str) -> str:
    """Map a universality_check dict to the PASS/FAIL/OPEN status table."""
    if sk_candidate == "none":
        return "PASS" if uc["universal"] else "FAIL"
    return "PASS" if uc["universal"] else "FAIL"


def run_comparison(
    chi: float = 0.0,
    transport_power: int = 2,
    candidates: Optional[tuple[str, ...]] = None,
) -> dict:
    """
    Run the closure ledger across multiple S(k) candidates plus the
    Layer-1-only baseline (`sk_candidate='none'`), and return a comparison
    dict suitable for serialization.

    Default candidate set: ('none',) + sk_bridge.WIRED_CANDIDATES.
    """
    from experiments.closure_ledger.sk_bridge import WIRED_CANDIDATES, s_k_membership

    if candidates is None:
        candidates = ("none",) + WIRED_CANDIDATES

    runs: dict[str, ExperimentResult] = {}
    for cand in candidates:
        runs[cand] = run_experiment(
            chi=chi, transport_power=transport_power, sk_candidate=cand,
        )

    # Status table: per-candidate PASS/FAIL with mod-2π residues.
    status_rows: list[dict] = []
    for cand, res in runs.items():
        uc = res.universality_check
        residues_pi = [
            v / math.pi for v in uc["per_lepton_mod_2pi"]
        ] if uc["per_lepton_mod_2pi"] else []
        status_rows.append({
            "candidate": cand,
            "experiment_layer": res.experiment_name.split(".")[-1],
            "result": _classify_universality(uc, cand),
            "per_lepton_mod_2pi_in_pi": residues_pi,
            "spread": uc["spread"],
            "universal_value_in_pi": (
                uc["universal_value"] / math.pi
                if uc.get("universal_value") is not None else None
            ),
            "overall_status": res.overall_status,
        })

    # Mode-map comparison: what does each candidate select per generation?
    mode_map_rows: list[dict] = []
    for cand in candidates:
        if cand == "none":
            continue
        for k in (1, 3, 5):
            members = s_k_membership(k, cand)
            mode_map_rows.append({
                "candidate": cand,
                "k": k,
                "modes": [(m.l, m.n) for m in members],
            })

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "transport_convention": (
            "T2_sign_flip" if transport_power == 2 else "T1_diagnostic"
        ),
        "chi": chi,
        "candidates_run": list(candidates),
        "status_table": status_rows,
        "mode_map_comparison": mode_map_rows,
        "runs": {cand: res.to_dict() for cand, res in runs.items()},
    }


def render_comparison_markdown(comparison: dict) -> str:
    """Render the multi-candidate comparison as a single markdown report."""
    lines: list[str] = []
    lines.append("# Closure-phase ledger — candidate comparison")
    lines.append("")
    lines.append(f"**Run:** {comparison['timestamp_utc']}")
    lines.append(
        f"**Transport convention:** `{comparison['transport_convention']}`"
    )
    lines.append(f"**χ (Hopf fibre):** {comparison['chi']}")
    lines.append("")

    lines.append("## Status table")
    lines.append("")
    lines.append(
        "| candidate | layer | result | per-lepton mod 2π (units of π) | "
        "spread | universal value |"
    )
    lines.append("|---|---|---|---|---|---|")
    for r in comparison["status_table"]:
        residues = (
            "[" + ", ".join(f"{v:.6f}" for v in r["per_lepton_mod_2pi_in_pi"]) + "]"
            if r["per_lepton_mod_2pi_in_pi"] else "—"
        )
        univ = (
            f"{r['universal_value_in_pi']:.6f}π"
            if r["universal_value_in_pi"] is not None else "—"
        )
        lines.append(
            f"| `{r['candidate']}` | {r['experiment_layer']} | **{r['result']}** | "
            f"{residues} | {r['spread']:.3e} | {univ} |"
        )
    lines.append("")

    lines.append("## Mode-map comparison")
    lines.append("")
    lines.append("| candidate | k=1 | k=3 | k=5 |")
    lines.append("|---|---|---|---|")
    by_cand: dict[str, dict[int, list]] = {}
    for row in comparison["mode_map_comparison"]:
        by_cand.setdefault(row["candidate"], {})[row["k"]] = row["modes"]
    for cand, ks in by_cand.items():
        cells = []
        for k in (1, 3, 5):
            modes = ks.get(k, [])
            cells.append(
                "{" + ", ".join(f"({l},{n})" for l, n in modes) + "}"
                if modes else "—"
            )
        lines.append(f"| `{cand}` | {cells[0]} | {cells[1]} | {cells[2]} |")
    lines.append("")

    lines.append("## Per-candidate one-liners")
    lines.append("")
    for r in comparison["status_table"]:
        lines.append(f"- `{r['candidate']}`: {r['overall_status']}")
    lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "- Layer 1 (`sk_candidate='none'`): closed-form topological "
        "ledger is universal mod 2π across leptons. PASS."
    )
    lines.append(
        "- Layer 2A (`A_lowest_radial_per_l`, WKB convention): "
        "cumulative odd-l ground modes do not preserve closure mod 2π. "
        "FAIL — rejects this radial-mode interpretation of lepton depth."
    )
    lines.append(
        "- Layer 2B1 (`B1_single_angular_mode`, WKB convention): "
        "single l=k ground mode per generation. The per-row Φ values are "
        "exactly the candidate-A summands; if A's cumulative residues are "
        "non-universal, the same per-mode numbers cannot be universal "
        "alone unless they happen to coincide mod 2π. FAIL."
    )
    lines.append(
        "- Layer 2B2 (`B2_single_radial_excitation`, WKB convention): "
        "single l=1 ladder. Asymptotic WKB gives Φ(1, n) → (n + 1) π, "
        "which produces a parity pattern across generations rather than "
        "a single universal value. FAIL."
    )
    lines.append(
        "- Layer 2C1 (`C1_eigenvector_weighted_B1`, WKB convention): "
        "B1 modes weighted by |v_species,i|² from the locked lepton "
        "generation block. The eigenvector mixing tightens the spread "
        "relative to B1 (the {1,3} mixing pulls e and μ residues toward "
        "each other) but does NOT close them to a universal value. "
        "FAIL — under this WKB convention, the surrogate Hamiltonian's "
        "own eigenvectors do not supply the missing bridge."
    )
    lines.append(
        "- Layer 2C2 (`C2_eigenvector_weighted_B2`, WKB convention): "
        "B2 modes weighted by the same eigenvectors. The (n+1)π "
        "spacing between B2 modes dominates the eigenvector mixing, "
        "leaving residues distributed across [0, 2π). FAIL."
    )
    lines.append(
        "- Layer 2 C1+Maslov (`C1_maslov_standard`, Bohr-Sommerfeld "
        "convention): same C1 mode set and weights, but each integrated "
        "phase is shifted by −π/2 per detected classical turning point "
        "(sign change of ω² − V_eff inside the tortoise grid). When the "
        "turning-point count is uniform across the B1 ground modes the "
        "Maslov shift is a uniform offset and the spread inherits from C1; "
        "when the count varies across modes the correction can in "
        "principle redistribute residues. The radial-detail table in each "
        "sub-run reports the per-mode `n_turning_points` so the regime can "
        "be read directly."
    )
    lines.append(
        "- Layer 2 B2+Maslov (`B2_maslov_standard`, Bohr-Sommerfeld "
        "convention): single-mode B2 ladder (l = 1, n = (k−1)/2) with the "
        "standard −π/2 turning-point shift. The B2 ladder is the natural "
        "test for differential Maslov: the (l=1, n=0) ground mode has a "
        "centrifugal-barrier soft turning point while the n ≥ 1 "
        "excitations sit above the barrier and have no interior turning "
        "points. Whether the differential shift collapses the (n+1)π "
        "parity pattern is read directly off the per-row Φ + Δ_maslov "
        "in the breakdown table."
    )
    lines.append(
        "- Layer 2 C2+Maslov (`C2_maslov_standard`, Bohr-Sommerfeld "
        "convention): C2's eigenvector-weighted B2 ladder with the same "
        "standard Maslov shift. Combines C-family weighting (lifts the "
        "B2 single-mode degeneracy) with the per-mode turning-point "
        "correction. Differential N_turning across the B2 ladder lets "
        "the eigenvector weights redistribute the Maslov shift across "
        "species — the most aggressive composition of the two known "
        "structural levers."
    )
    lines.append("")

    return "\n".join(lines)


def write_comparison_outputs(comparison: dict, out_dir: Path) -> None:
    """Write the comparison JSON, markdown, and per-candidate sub-runs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "comparison.json").write_text(
        json.dumps(comparison, indent=2, default=str)
    )
    (out_dir / "status_table.md").write_text(render_comparison_markdown(comparison))
    for cand, run in comparison["runs"].items():
        sub = out_dir / cand
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "result.json").write_text(json.dumps(run, indent=2, default=str))
        # Reconstruct an ExperimentResult-like object for markdown rendering.
        # We use the same fields verbatim; render via a minimal helper.
        result = ExperimentResult(**run)
        (sub / "summary.md").write_text(result._render_markdown())
