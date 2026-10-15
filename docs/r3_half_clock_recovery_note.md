# Recovery of interrupted execution

2026-10-10, after 20:30 UTC. The original producer execution disappeared
after printing successful convergence for 63 of 84 phases. It had saved
`known.json` but no primary scan, noise checks or result. The individual
primary obstruction values were not inspected or persisted and cannot be
recovered. No finding is assigned to that incomplete attempt.

On the user's instruction to continue, preserve `known.json` and execute
all 84 phase calculations again with the same equations, seeds, precision,
tolerances, Newton limits and registered scoring code. This is an execution
recovery, not a retry after an unfavorable scientific result. It is a second
attempt at executing the calculation, not an independent replication.

The recovery harness schedules four independent phases at a time and saves
each result immediately. Its numerical per-phase body was extracted from
the frozen producer. An AST comparison, excluding only scheduling and print
statements, verifies its identity before any calculation. On interruption,
completed source-bound checkpoints are reused; failed numerical points are
not retried. The completed scan is ordered by the original phase index.

The original source and protocol remain unchanged. The recovery implementation
and this note are published before restarting production. An additional
`execution.json` binds this harness, the frozen inputs, and hashes of the
per-phase checkpoints. The checkpoints are redundant with `scan.json`;
the aggregate archive contains all their numerical point records. The
previously published supplemental validation is still required for confirmation.
