# `experiments/closure_ledger`

**Closure-phase ledger experiment for Bulk Antipodal Mechanics (BAM).**

This experiment tests whether the locked stable particle states share a
universal closure-phase invariant modulo `2π`. It is a falsification
test for the BAM thesis that quantum particles are self-consistent
topological boundary conditions — under that thesis, every stable
state's worldline closure should accumulate the same total geometric
phase (modulo `2π`) from the four channels BAM identifies (antipodal
closure, Hopf holonomy, throat transport, bulk eigenmodes).

## Run

From the repo root:

```bash
python -m experiments.closure_ledger
```

Writes `result.json` and `summary.md` to a timestamped directory under
`experiments/closure_ledger/runs/`.

Options:

```bash
python -m experiments.closure_ledger --chi 0.0          # Hopf fibre angle (default 0)
python -m experiments.closure_ledger --transport-power 1 # T¹ diagnostic instead of T² closure
python -m experiments.closure_ledger --no-write          # stdout only, no run dir
```

## What the experiment computes

For each lepton generation `k ∈ {1, 3, 5}` (electron, muon, tau), the
ledger sums six phase contributions. Four are **wired into the repo**
and produce values; two are **structurally missing** and produce
`null`:

| term                  | source                                        | status        |
|-----------------------|-----------------------------------------------|---------------|
| antipodal closure     | `k · action_base` (= `k · 2π`)                | available     |
| Hopf holonomy         | `∮A = π·cos(χ)`                                | available     |
| throat transport `T²` | `T = iσ_y`, `T² = −I` (phase π per closure)   | available     |
| β-uplift              | `β · max(0, k−3)²`, `4β = 100·(2π)`            | available     |
| **radial bulk phase** | `Σ ∫ k_local(r*) dr*` over `S(k)` mode set     | **missing**   |
| moving-mouth phase    | time-dependent `T(t)` along closed trajectory  | not modeled   |

The four available terms are computed and reported. The two missing
terms produce a structured Layer 2 blocker report (see below).

## Layered design

- **Layer 1 (runs to completion)**: ledger from the four wired
  channels. Tests whether their sum closes mod `2π` universally
  across `e, μ, τ`. Under the `T²` convention with `χ = 0`, the
  answer is **yes, all three lepton ledgers close to 0 mod 2π**.

- **Layer 2 (structured blocker)**: explicit "no S(k) bridge"
  report citing repo evidence and listing three principled
  candidate maps. Defining `S(k)` is the next theoretical target,
  per the repo audit.

## Why Layer 2 is blocked

Per the repo audit:

- `geometrodynamics/tangherlini/lepton_spectrum.py` is an
  **instanton-transition surrogate** operating on depth labels
  `k ∈ {1, 3, 5}`. The radial-eigenmode parameters (`l`, `n_points`,
  `rs`, `r_outer`) are accepted for API compatibility but explicitly
  unused in `compute_knotted_lepton_spectrum`.

- The locked diagonal Hamiltonian
  `H_kk = action_base + resistance_scale · k² + res_diag(k) +
  pinhole(k ∈ {3, 5}) + β · max(0, k−3)²` has no eigenmode index
  and no Tangherlini potential evaluation.

- The quark residual sector reads scalars off the tortoise grid
  (transport, pinhole, resistance) but those are integrated quantities
  for the whole sector, not a per-generation `S(k) → (l, n)` map.

The bridge from generation depth `k` to Tangherlini eigenmodes
`(l, n)` — call it `S(k)` — does not yet exist. Until it does,
`Φ_radial(k)` is not computable and the full closure-phase invariant
remains a partial result.

## Predictions

Three falsification-level predictions, plus their current status:

| Prediction                                                                       | Status                       |
|----------------------------------------------------------------------------------|------------------------------|
| **P1**. Lepton universal closure: `Φ_total(e) ≡ Φ_total(μ) ≡ Φ_total(τ) (mod 2π)` | partial pass (4 of 6 terms) |
| **P2**. Quark sector universal closure                                            | blocked on `S(k)`            |
| **P3**. Lepton-quark closure-quanta gap = 466 − 100 = 366                         | hypothesis (S(k)-conditional)|

The Layer 1 result (P1 partial pass) is consistent with the
universal-closure conjecture: the four wired channels do agree
mod `2π` across all three leptons. Whether the radial channel
also agrees mod `2π` — making the conjecture's verdict definitive —
requires `S(k)` to be defined.

## Convention defaults

The experiment defaults to:

- **Throat transport convention `T² = −I`** (phase π per closure).
  This is the validated closure signature in `embedding/transport.py`.
  The `T¹` single-throat-pass interpretation (phase π/2 per closure)
  is available as a diagnostic via `--transport-power 1`. Both
  conventions yield universal closure; they differ in the universal
  value (0 vs 3π/2 mod 2π).

- **Hopf fibre `χ = 0`**, the canonical fibre at the equator of the
  Hopf base S². Holonomy at this fibre is `∮A = π · cos(0) = π`.

## Repo dependencies

The experiment **prefers to import constants from the repo** to stay
in sync with canonical definitions, but **falls back to README-published
values** on any `ImportError`. Fallbacks are recorded in the result
JSON's `import_errors` field and in the markdown summary.

Imports the runner attempts (Claude Code: verify these against the
actual modules on commit and adjust the import paths in `ledger.py`
if any name differs):

- `from geometrodynamics.tangherlini.lepton_spectrum import S3_ACTION_BASE, TAU_BETA_50PI`
- `from geometrodynamics.qcd.quark_spectrum import BETA_QUARK_LOCK` (or `LOCKED_BETA`)
- `from geometrodynamics.hopf.connection import hopf_holonomy`
- `from geometrodynamics.embedding.transport import T_matrix` (or `T`)

Where the name doesn't match, the ledger falls back to the closed-form
value. The fallback is correct for Layer 1 because all four phase
contributions are analytically expressible from the README.

## Testing

```bash
pytest experiments/closure_ledger/test_closure_ledger.py
```

The tests do not require the geometrodynamics package to be
importable — they verify that the experiment runs end-to-end with
fallback values, that Layer 1 closes mod 2π universally, that
Layer 2 produces a structured blocker, and that JSON/markdown
serialization round-trips.

When Claude Code verifies the real import paths, additional
integration tests against actual `geometrodynamics.*` symbols can
be added.

## Output layout

```
experiments/closure_ledger/runs/<UTC timestamp>/
├── result.json     # full structured output (LedgerRow dicts, blocker, etc.)
└── summary.md      # human-readable per-lepton table + universality + blocker
```

Both are committable artifacts; the JSON is machine-readable for
future cross-run analysis.

## What success looks like

A passing run should report:

- All four available terms valued for each lepton.
- Available-term sums universal mod 2π across `e, μ, τ` (universal
  value 0 under `T²`).
- Quark lock quanta gap = 366.
- Layer 2 blocker non-empty with three candidate `S(k)` maps and a
  next-steps list.
- `import_errors` empty (after Claude Code verifies the imports) or
  populated with the symbols that need adapter shims.

## What failure modes mean

- **Available-term sums NOT universal mod 2π**: Layer 1 fails. Either
  `T²` is the wrong convention, or one of the four wired channels
  has been miscomputed. Falsifies `P1` at the partial level.

- **Constants don't match README values when imported**: the lock
  values themselves drift in the repo. The fallback values are correct
  per README; if the repo has changed, update both the README and the
  `_FALLBACK_*` constants in `ledger.py`.

- **Quark lock quanta ≠ 466**: the quark β-lock has been changed in
  the repo. Update `_FALLBACK_BETA_QUARK` and the README accordingly.

## See also

- `docs/THESIS.md` — BAM conjecture, three mechanisms, falsification roadmap.
- `docs/lepton_axioms.md` — locked lepton baseline, β-lock derivation log.
- `docs/quark_axioms.md` — locked quark spectrum, residual-sector
  geometrization, β = 466·π/2 phenomenological status.
