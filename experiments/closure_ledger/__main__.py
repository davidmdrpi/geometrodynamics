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

from experiments.closure_ledger.runner import (
    render_comparison_markdown,
    run_comparison,
    run_experiment,
    write_comparison_outputs,
)
from experiments.closure_ledger.sk_bridge import (
    DEFAULT_SK_CANDIDATE,
    SK_CANDIDATES,
)


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
        "--sk-candidate",
        choices=SK_CANDIDATES,
        default=DEFAULT_SK_CANDIDATE,
        help=(
            "S(k) → {(l, n)} bridge candidate. Default wires the radial "
            "channel via candidate A (lowest radial mode per odd l up to k). "
            "Pass 'none' to reproduce the Layer-1 ledger."
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
    parser.add_argument(
        "--compare",
        action="store_true",
        help=(
            "Run all wired S(k) candidates plus the Layer-1 baseline and "
            "write a comparison directory with a status table, mode-map "
            "comparison, and per-candidate sub-runs."
        ),
    )
    args = parser.parse_args(argv)

    if args.compare:
        comparison = run_comparison(
            chi=args.chi, transport_power=args.transport_power,
        )
        print(render_comparison_markdown(comparison))
        print()
        if not args.no_write:
            if args.output_dir is None:
                ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                here = Path(__file__).resolve().parent
                out = here / "runs" / f"{ts}_comparison"
            else:
                out = Path(args.output_dir)
            write_comparison_outputs(comparison, out)
            print(f"Wrote: {out}/")
        return 0

    result = run_experiment(
        chi=args.chi,
        transport_power=args.transport_power,
        sk_candidate=args.sk_candidate,
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
