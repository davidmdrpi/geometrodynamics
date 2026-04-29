"""
CLI entry point for the closure-phase ledger experiment.

Usage (from repo root):

    python -m experiments.closure_ledger
    python -m experiments.closure_ledger --chi 0.0
    python -m experiments.closure_ledger --transport-power 1     # T¹ diagnostic
    python -m experiments.closure_ledger --no-write               # stdout only

Default behaviour: run the experiment with T² convention and χ = 0,
write `result.json` and `summary.md` under
`experiments/closure_ledger/runs/<UTC timestamp>/`, and print the
summary to stdout.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

from experiments.closure_ledger.runner import run_experiment


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m experiments.closure_ledger",
        description=(
            "Run the BAM closure-phase ledger experiment. Layer 1 (the "
            "wired channels) runs to completion; Layer 2 (radial bulk "
            "phase) emits a structured blocker citing the missing S(k) "
            "bridge in lepton_spectrum.py."
        ),
    )
    parser.add_argument(
        "--chi",
        type=float,
        default=0.0,
        help="Hopf fibre angle χ (default: 0, canonical fibre, holonomy = π).",
    )
    parser.add_argument(
        "--transport-power",
        type=int,
        default=2,
        choices=[1, 2],
        help=(
            "Throat transport power: 2 for T² closure convention "
            "(default; validated in embedding/transport.py), "
            "1 for T¹ single-pass diagnostic."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Where to write run artifacts. Default: "
            "experiments/closure_ledger/runs/<UTC timestamp>/"
        ),
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Skip writing files; print summary to stdout only.",
    )
    args = parser.parse_args(argv)

    result = run_experiment(
        chi=args.chi,
        transport_power=args.transport_power,
    )

    print(result._render_markdown())
    print()

    if not args.no_write:
        if args.output_dir is None:
            ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            here = Path(__file__).resolve().parent
            out = here / "runs" / ts
        else:
            out = Path(args.output_dir)
        result.write_json(out / "result.json")
        result.write_summary_markdown(out / "summary.md")
        print(f"Wrote: {out / 'result.json'}")
        print(f"Wrote: {out / 'summary.md'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
