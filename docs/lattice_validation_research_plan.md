# High-resolution lattice validation of the BAM measure operators (PR #120)

The measure arc (PRs #116–#119) computed determinants and the η-invariant
from **continuous analytic** derivations — the Tangherlini Gel'fand–Yaglom
determinant, the ghost determinant `det'(−∂_τ²) = L²`, the first-order
`det'(∂_τ) = L`, the antiperiodic `det = 2`, and `η = 0`. This PR is a
**validation**: it checks that the **discrete software** (the
finite-difference / lattice operators actually used in the code) reproduces
those continuum results — **exactly** for the structural/symmetry quantities
and with the expected **`O(1/N²)`** convergence for the finite-difference
quantities, at high lattice resolution.

## What is validated, and how

  - **Eigenvalues** of the discrete `−∂_τ²` → continuum `(2πk/L)²`, relative
    error `O(1/N²)` (ratio ≈ 16 per `N×4`).
  - **Ghost determinant** (periodic): the lattice log-determinant
    `Σ_k log[2 − 2cos(2πk/N) + (mh)²]` → log of the continuum
    `(2 sinh(mL/2))²`, `O(1/N²)`; cross-checked by the transfer-matrix closed
    form `2(cosh Nα − 1)` `[2cosh α = 2 + m²h²]` at `N = 10⁶`.
  - **Antiperiodic determinant** → continuum `(2 cosh(mL/2))²`, `O(1/N²)`.
  - **`det'(−∂_τ²) = L²`** and **`det(∂_τ)_antiperiodic = 2`**: the continuum
    `m → 0` limits (`(2 sinh)²/m² → L²`, `(2 cosh)² → 4`), with the discrete
    matching the continuum at fixed `m` (PR #117/#119).
  - **η-invariant = 0** and the **zero-mode count**: the discrete centered
    `∂_τ` (odd `N`) has a spectrum symmetric under `k ↔ N−k`, so `Σ sign = 0`
    **exactly** at finite `N`, with exactly one zero mode and a purely
    imaginary spectrum (PR #118/#119).
  - **Tangherlini Gel'fand–Yaglom** `det(H)/det(H_free) → 1.574370` (PR #116)
    at high resolution, stable to ~1e-7 by `N = 2000`.

Structural/symmetry quantities (`η`, zero-mode count, spectrum symmetry)
match **exactly** at finite `N`; finite-difference quantities converge as
`O(1/N²)`. So the discrete implementation behaves as the continuum analytic
derivation in the high-resolution limit.

## Numerical-robustness notes

  - The bare eigenvalue product `Π(4/h²)sin²(...)` overflows from the `1/h²`
    prefactor; the physical object is the **dimensionless** lattice
    determinant `Π[2 − 2cos(2πk/N) + (mh)²]`, computed as a **log-determinant**
    (sum of logs) to avoid overflow/underflow at high `N`.
  - The **odd-`N`** centered difference avoids the Nyquist zero-mode artifact
    (even `N` has a spurious zero at `k = N/2`), so the first-order operator
    shows exactly one zero mode (cf. PR #118).

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | validate discrete ≡ continuum (PRs #116–#119) |
| T2 | eigenvalues | `−∂_τ²` → `(2πk/L)²`, `O(1/N²)` (ratio 16) |
| T3 | ghost det | lattice log-det → `(2sinh(mL/2))²`, `O(1/N²)`; TM at `N=10⁶` |
| T4 | antiperiodic | lattice → `(2cosh(mL/2))²`, `O(1/N²)` |
| T5 | `m→0` limits | `det'(−∂_τ²)=L²`, `det_AP=2` |
| T6 | η exact | `η = 0` EXACT at finite `N` + 1 zero mode (centered `∂_τ`, odd `N`) |
| T7 | Tangherlini GY | `det(H)/det(H_free) → 1.574370` (PR #116), high-`N` |
| T8 | assessment | `LATTICE_VALIDATION_DISCRETE_MATCHES_CONTINUUM` |

## Established and open

  - **Established (BAM-native):** the discrete lattice operators reproduce
    every continuum analytic result of PRs #116–#119 — eigenvalues and
    determinants converge `O(1/N²)` (ghost → `(2sinh)²`, antiperiodic →
    `(2cosh)²`, Tangherlini GY → 1.574370), and the structural quantities
    (`η = 0`, one zero mode, spectrum symmetry) match **exactly** at finite
    `N`. The software behaves exactly as the continuous derivation.

  - **Open (unchanged):** the analytic open pieces of the measure arc remain
    as before (absolute `Z` normalization `κ₅²/Λ₅`, multi-loop, the
    intermediate-holonomy η-phase). This probe validates the numerics, not
    those.

## Cross-references

  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116.
  - `docs/diff_s1_ghost_determinant_research_plan.md` — PR #117.
  - `docs/diff_s1_first_order_ghost_audit_research_plan.md` — PR #118.
  - `docs/detprime_dtau_eta_invariant_phase_research_plan.md` — PR #119.

## Run

```
python -m experiments.closure_ledger.lattice_validation_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_lattice_validation_probe/`.
Expected verdict: `LATTICE_VALIDATION_DISCRETE_MATCHES_CONTINUUM`, 8/8 PASS.
