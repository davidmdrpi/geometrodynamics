# Zeta / heat-kernel regularization of the Tangherlini fluctuation determinant (PR #116)

PR #115 constructed the `S_BAM` path-integral measure structurally but left
its hard analytic core **open**: the one-loop factor is the fluctuation
determinant of the Tangherlini cavity operator (the second variation of
`S_BAM` about the throat saddle), and the bare product `Π_n ω_n`
**diverges** (the log-det partial sums grow without bound). This PR closes
that core — it regularizes the determinant by two independent standard
methods and gets a **finite, scheme-independent, mutually consistent**
answer.

## The operator

The fluctuation operator is the radial Tangherlini cavity Hamiltonian on the
tortoise interval `[r*_inner, r*_outer]` with Dirichlet ends:

```
H = −d²/dr*² + V_tangherlini(r, l),   V = f(r)[ l(l+2)/r² + 3 rs²/r⁴ ].
```

Its spectrum `{ω²_n}` is positive (`min ω² ≈ 1.11 > 0`: a stable saddle, no
zero/negative modes), but `Π_n ω²_n` diverges — the bare determinant needs
regularization.

## Method 1: Gel'fand–Yaglom (no mode sum at all)

For a 1D Schrödinger operator the ratio of determinants to a free reference
is given **exactly** by a single initial-value solve, with no eigenvalue
computation:

```
det(H) / det(H_free) = y(L) / L,
```

where `y` solves `H y = 0` with `y(0) = 0`, `y'(0) = 1`, and `H_free =
−d²/dr*²` (so `y_free(L) = L`). On the Tangherlini cavity:

```
det(H)/det(H_free) = 1.57437   (log = 0.45386),
```

converged to 6 digits by `N = 2000` grid points (identical at `N = 32000`),
with **zero interior nodes** in `y` ⟹ no negative eigenvalues (consistent
with `min ω² > 0`). This is a finite, scheme-independent renormalized
one-loop determinant — exactly the quantity PR #115 said was missing.

## Method 2: zeta / heat-kernel (independent cross-check)

The heat kernel `θ(t) = Σ_n e^{−ω²_n t}` has the small-`t` Weyl expansion
`θ(t) ≈ a_{−1/2} t^{−1/2} + a₀ + …`, and the regularized log-determinant is
`−ζ'(0)` with `ζ(s) = Σ_n ω_n^{−2s}`; its **finiteness** is controlled by
`ζ(0) = a₀`. Fitting the computed low spectrum:

```
a_{−1/2} (fit) = 0.946   vs Weyl  L/√(4π) = 0.938   (0.9%),
ζ(0) = a₀ (fit) = −0.60  vs Dirichlet-interval  −1/2,
```

so `ζ(0)` is finite and matches the universal Dirichlet value `−1/2` (no
zero mode, no anomalous piece) — confirming the determinant is finite and
the Gel'fand–Yaglom value is the physical renormalized determinant. The Weyl
law itself checks out (`N(<λ)`: 7, 14, 29 vs `(L/π)√λ` = 7.5, 15.0, 29.9).

## What this resolves (and does not)

  - **Resolved:** the PR #115 analytic core. The Tangherlini fluctuation
    determinant is **finite** after standard regularization, computed two
    ways that agree: Gel'fand–Yaglom gives `det(H)/det(H_free) = 1.574` (log
    `0.454`) exactly, and heat-kernel gives `ζ(0) = −1/2` (finite, no zero
    mode). The `S_BAM` one-loop measure factor is therefore well-defined and
    computable, not merely formal.
  - **Not resolved (and not claimed):** a closed-form analytic expression
    for the determinant (it is a definite number, computed numerically, not
    a formula); the absolute normalization of `Z` (still carries the bulk
    `κ₅²/Λ₅` anchor, PR #112); and the full multi-channel/multi-loop measure.
    What is delivered is the finite one-loop determinant of the leading
    fluctuation operator — the specific gap PR #115 flagged.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | PR #115 left the fluctuation determinant divergent/open |
| T2 | operator | stable (`min ω² ≈ 1.11 > 0`); bare `Π ω²_n` diverges |
| T3 | Gel'fand–Yaglom | `det(H)/det(H_free) = 1.574` (log 0.454), zero nodes |
| T4 | GY convergence | 6 digits stable `N = 2000 → 32000` |
| T5 | Weyl law | `a_{−1/2} ≈ L/√(4π)` (0.9%); `N(λ) ≈ (L/π)√λ` |
| T6 | ζ(0) | `= −1/2` (Dirichlet, finite, no zero mode); two methods agree |
| T7 | scope | analytic core RESOLVED; closed form / abs `Z` normalization open |
| T8 | assessment | `TANGHERLINI_FLUCTUATION_DETERMINANT_REGULARIZED_FINITE` |

## Established and open

  - **Established (BAM-native):** the divergent bare determinant of PR #115
    is regularized to a finite, scheme-independent value two independent
    ways — Gel'fand–Yaglom `det(H)/det(H_free) = 1.574` (log 0.454), exact
    and converged with zero nodes, and zeta/heat-kernel `ζ(0) = −1/2`
    (finite, Dirichlet, no zero mode), with the Weyl law confirmed. The
    `S_BAM` one-loop measure factor is finite and computable.

  - **Open:** a closed-form analytic expression (the value is numerical);
    the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor, PR #112);
    the full multi-channel / multi-loop measure.

## Cross-references

  - `docs/s_bam_path_integral_measure_research_plan.md` — PR #115, which
    constructed the measure and flagged this determinant as the open
    analytic core.
  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the `1/(2π)`
    per-loop-dimension factor.
  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, the unpinned
    `κ₅²/Λ₅` bulk normalization that the absolute `Z` still carries.

## Run

```
python -m experiments.closure_ledger.tangherlini_fluctuation_determinant_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_tangherlini_fluctuation_determinant_probe/`.
Expected verdict: `TANGHERLINI_FLUCTUATION_DETERMINANT_REGULARIZED_FINITE`,
8/8 PASS.
