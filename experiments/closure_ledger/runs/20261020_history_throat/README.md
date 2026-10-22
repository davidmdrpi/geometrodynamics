# History-coupled throat replay data

CSV and JSON files are directly readable. The complete NPZ archives are split
into lossless parts of at most 480,000 bytes so interrupted uploads can resume.
`replay_parts.json` records every part hash and the original archive hash.

From the repository root:

```bash
python -m experiments.closure_ledger.history_throat_report --assemble
```

This verifies and reassembles the exact NPZ files; it does not rerun a trajectory.
The solver and test data are classical. Geometric event tables are ray predictions,
not measured tagged pulses. Read the protocol and results before interpreting them.
