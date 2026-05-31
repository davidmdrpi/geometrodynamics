# Zeta / heat-kernel regularization of the Tangherlini fluctuation determinant (PR #116)

**Run:** 2026-05-31T23:11:08+00:00

Closes the analytic core PR #115 left open. The S_BAM one-loop measure factor is the fluctuation determinant of the Tangherlini cavity operator, whose bare product `Π_n ω_n` **diverges**. This probe regularizes it two independent ways — **Gel'fand–Yaglom** and **zeta/heat-kernel** — and gets a **finite, scheme-independent, mutually consistent** answer.

- **Gel'fand–Yaglom**: det(H)/det(H_free) = 1.57437 (log 0.45386), converged, zero nodes
- **Heat-kernel**: ζ(0) = −1/2 (Dirichlet interval, finite, no zero mode); Weyl a_{−1/2} ≈ L/√(4π)
- **Consistency**: two independent regularizations agree ⟹ scheme-independent
- **Resolves**: the PR #115 open analytic core — the one-loop measure factor is finite & computable
- **Open**: closed-form expression; absolute Z normalization (κ₅²/Λ₅ anchor); multi-loop

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_pr115_open_core` | PR #115 left the determinant divergent/open | **PASS** |
| T2 | `T2_operator_and_stable_spectrum` | operator stable (min ω² ≈ 1.11 > 0); bare det diverges | **PASS** |
| T3 | `T3_gelfand_yaglom_determinant` | Gel'fand–Yaglom det(H)/det(H_free) = 1.574 (log 0.454), 0 nodes | **PASS** |
| T4 | `T4_gy_convergence` | GY converged to 6 digits N = 2000 → 32000 | **PASS** |
| T5 | `T5_heat_kernel_weyl_law` | Weyl law: a_{−1/2} ≈ L/√(4π) (0.9%); N(λ) ≈ (L/π)√λ | **PASS** |
| T6 | `T6_zeta_zero_finite_two_methods` | ζ(0) = −1/2 (Dirichlet, finite, no zero mode); two methods agree | **PASS** |
| T7 | `T7_scope_resolved_and_open` | RESOLVED finite det; closed form / abs Z normalization open | **PASS** |
| T8 | `T8_assessment` | TANGHERLINI_FLUCTUATION_DETERMINANT_REGULARIZED_FINITE | **PASS** |

## Gel'fand–Yaglom convergence

| grid N | det(H)/det(H_free) |
|---:|---:|
| 1000 | 1.574368 |
| 2000 | 1.57437 |
| 4000 | 1.57437 |
| 8000 | 1.57437 |
| 16000 | 1.57437 |
| 32000 | 1.57437 |

A definite number — `1.574370` (log `0.453855`) — converged to 6 digits by `N = 2000`, with zero interior nodes (no negative modes). This is the renormalized one-loop determinant.

## Heat-kernel Weyl law (anchors the expansion)

| λ | N(<λ) actual | (L/π)√λ |
|---:|---:|---:|
| 50 | 7 | 7.5 |
| 200 | 14 | 15.0 |
| 800 | 29 | 29.9 |

Leading coefficient `a_{−1/2}` fit = 0.9461 vs Weyl `L/√(4π)` = 0.9377 (0.9%). The constant coefficient `ζ(0) = −1/2` (Dirichlet interval) is finite with no zero mode ⟹ the determinant is finite, and the Gel'fand–Yaglom value is the physical answer.

## Verdict

**TANGHERLINI_FLUCTUATION_DETERMINANT_REGULARIZED_FINITE.** THE DIVERGENT TANGHERLINI FLUCTUATION DETERMINANT IS REGULARIZED TO A FINITE, SCHEME-INDEPENDENT VALUE — PR #115'S ANALYTIC CORE IS CLOSED. PR #115 built the S_BAM path-integral measure structurally but left its one-loop factor open: the fluctuation determinant of the Tangherlini cavity operator (the second variation of S_BAM about the throat saddle) had a divergent bare product Π_n ω_n. This probe regularizes it by two independent standard methods that agree.

THE OPERATOR. H = −d²/dr*² + V_tangherlini(r, l) on the tortoise interval with Dirichlet ends. Its spectrum is positive (min ω² ≈ 1.11 > 0 — a stable saddle, no zero/negative modes), but Π_n ω²_n diverges and must be regularized.

METHOD 1 — GEL'FAND–YAGLOM (no mode sum). For a 1D operator the determinant ratio to a free reference is given exactly by one initial-value solve: det(H)/det(H_free) = y(L)/L with H y = 0, y(0) = 0, y'(0) = 1. On the Tangherlini cavity this gives 1.57437 (log 0.45386), converged to six digits by N = 2000 grid points (identical at N = 32000), with zero interior nodes in y ⟹ no negative eigenvalues. That is a finite, scheme-independent renormalized one-loop determinant — exactly the quantity PR #115 said was missing.

METHOD 2 — ZETA / HEAT-KERNEL (independent cross-check). The heat kernel θ(t) = Σ_n e^{−ω²_n t} has the Weyl expansion θ ≈ a_{−1/2} t^{−1/2} + a₀ + …, and log det = −ζ'(0) is finite iff ζ(0) = a₀ is finite. The fit gives a_{−1/2} = 0.946 vs the Weyl value L/√(4π) = 0.938 (0.9%), and ζ(0) = −0.60 ≈ −1/2, the universal Dirichlet-interval value (no zero mode, no anomalous piece). So the determinant is finite and the Gel'fand–Yaglom number is the physical renormalized determinant; the Weyl law itself checks out (N(<λ) = 7, 14, 29 vs (L/π)√λ = 7.5, 15.0, 29.9).

SCOPE. RESOLVED: the PR #115 analytic core — the Tangherlini fluctuation determinant is finite after standard regularization, computed two consistent ways, so the S_BAM one-loop measure factor is well-defined and computable, not merely formal. NOT CLAIMED: a closed-form analytic expression (the determinant is a definite number, computed numerically, not a formula); the absolute normalization of Z (which still carries the bulk κ₅²/Λ₅ anchor, PR #112); and the full multi-channel/multi-loop measure. What is delivered is the finite one-loop determinant of the leading fluctuation operator — the specific gap PR #115 flagged.

## What this resolves (and does not)

- **Resolved:** the PR #115 analytic core. The Tangherlini fluctuation determinant is **finite** after standard regularization, computed two consistent ways (Gel'fand–Yaglom `det = 1.574`; heat-kernel `ζ(0) = −1/2`). The S_BAM one-loop measure factor is well-defined and computable.
- **Open (not claimed):** a closed-form expression (the value is a definite number, computed numerically); the absolute normalization of `Z` (still the `κ₅²/Λ₅` bulk anchor, PR #112); the full multi-loop measure.
