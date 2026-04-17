"""
History subpackage — closed-history framework.

The unifying principle: only globally closed, conservation-respecting
histories become physical events.  Connects Bell, QED transactions,
and black-hole thermodynamics through a shared closure condition.
"""

from geometrodynamics.history.closure import (
    EventType,
    Event,
    Worldline,
    History,
    ClosureResult,
    BranchResult,
    make_bell_history,
    make_transaction_history,
    enumerate_bell_branches,
)

__all__ = [
    "EventType",
    "Event",
    "Worldline",
    "History",
    "ClosureResult",
    "BranchResult",
    "make_bell_history",
    "make_transaction_history",
    "enumerate_bell_branches",
]
