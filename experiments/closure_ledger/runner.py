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


TAU = 2.0 * math.pi


@dataclass
class ExperimentResult:
    timestamp_utc: str
    experiment_name: str
    transport_convention: str          # "T2_sign_flip" | "T1_diagnostic"
    chi: float                          # Hopf fibre angle used
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

        # Universality
        u = self.universality_check
        lines.append("## Universality check (Layer 1)")
        lines.append("")
        lines.append(f"- All available-term sums mod 2π: "
                     f"{[f'{v:.6f}' for v in u['per_lepton_mod_2pi']]}")
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
) -> str:
    parts: list[str] = []
    if universality["universal"]:
        parts.append(
            f"Layer 1 PASS: lepton ledger universal mod 2π at "
            f"{universality['universal_value'] / math.pi:.3f}π."
        )
    else:
        parts.append("Layer 1 FAIL: lepton ledger does NOT close universally.")
    parts.append("Layer 2 BLOCKED: S(k) bridge undefined.")
    if import_errors:
        parts.append(
            f"({len(import_errors)} repo import(s) fell back to README values.)"
        )
    return " ".join(parts)


def run_experiment(
    chi: float = 0.0,
    transport_power: int = 2,
) -> ExperimentResult:
    """
    Run the closure-phase ledger experiment and return a complete
    ExperimentResult. No I/O — call .write_json(...) and
    .write_summary_markdown(...) on the result to persist.
    """
    rows, constants, import_errors = compute_lepton_ledger(
        chi=chi,
        transport_power=transport_power,
    )

    universality = _build_universality_check(rows)
    quark_sector = compute_quark_sector_summary(constants)
    blocker = layer2_blocker_report()

    convention = (
        "T2_sign_flip" if transport_power == 2
        else "T1_diagnostic"
    )

    return ExperimentResult(
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        experiment_name="closure_ledger.layer1",
        transport_convention=convention,
        chi=chi,
        constants=asdict(constants),
        import_errors=import_errors,
        rows=[r.to_dict() for r in rows],
        universality_check=universality,
        quark_sector=quark_sector,
        sk_bridge_blocker=blocker.to_dict(),
        overall_status=_overall_status(universality, blocker, import_errors),
    )
