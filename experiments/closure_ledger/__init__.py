"""
Closure-phase ledger experiment for Bulk Antipodal Mechanics (BAM).

Tests whether stable particle states share a universal closure-phase
invariant modulo 2π. See README.md for the experiment's motivation,
the layered design (Layer 1 runnable now, Layer 2 blocked on the S(k)
bridge), and the falsification predictions.

Public API:

    from experiments.closure_ledger import run_experiment, ExperimentResult

    result = run_experiment()
    result.write_json(path)
    result.write_summary_markdown(path)
"""

from experiments.closure_ledger.ledger import (
    PhaseTerm,
    LedgerRow,
    compute_lepton_ledger,
    compute_quark_sector_summary,
)
from experiments.closure_ledger.blockers import (
    SkBridgeBlocker,
    layer2_blocker_report,
)
from experiments.closure_ledger.runner import (
    ExperimentResult,
    run_experiment,
)

__all__ = [
    "PhaseTerm",
    "LedgerRow",
    "compute_lepton_ledger",
    "compute_quark_sector_summary",
    "SkBridgeBlocker",
    "layer2_blocker_report",
    "ExperimentResult",
    "run_experiment",
]
