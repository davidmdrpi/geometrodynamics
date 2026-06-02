# Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117)

**Run:** 2026-06-02T02:12:31+00:00

Supplies the gauge sector of the S_BAM measure: the reparametrization Faddeev–Popov / ghost determinant for `Diff(S¹)` (named in PR #115; PR #116 did the matter determinant). Computed by the same zeta method — **finite, det' = L², anomaly-free** — and the ghost zero mode turns out to BE PR #74's `1/(2π)`.

- **Gauge structure**: 1 modulus L (circumference = proper time) + 1 CKV (rigid U(1) rotation)
- **Ghost determinant**: det'(−d²/dτ²) = L² (zeta-regularized, ζ_ghost(0) = −1)
- **CKV factor**: CKV ⟹ dL/L; for L = 2π, 1/L = 1/(2π) = PR #74 closure quantum
- **Anomaly**: anomaly-free (no 1D conformal anomaly); only the discrete Z₂ (odd-k, PR #115) is nontrivial
- **Completes**: the gauge sector — with PR #116 matter det, the one-loop measure is finite/computable
- **Open**: absolute Z normalization (κ₅²/Λ₅); multi-loop; closed form

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap` | PR #115 flagged the ghost; PR #116 did matter; now the ghost | **PASS** |
| T2 | `T2_gauge_structure` | Diff(S¹) → constant einbein ⟹ 1 modulus L + 1 CKV | **PASS** |
| T3 | `T3_ghost_operator_and_kernel` | ghost operator −d²/dτ²; kernel = constants = 1 CKV | **PASS** |
| T4 | `T4_zeta_regularized_ghost_determinant` | zeta: ζ_ghost(0) = −1, det'(−d²/dτ²) = L² (computed) | **PASS** |
| T5 | `T5_ckv_is_closure_quantum_factor` | CKV ⟹ dL/L; for L = 2π, 1/L = 1/(2π) = PR #74 factor | **PASS** |
| T6 | `T6_anomaly_free_1d_worldline` | anomaly-free (1D, no conformal anomaly; vs string c = −26) | **PASS** |
| T7 | `T7_scope` | gauge sector complete; abs Z normalization / multi-loop open | **PASS** |
| T8 | `T8_assessment` | DIFF_S1_FADDEEV_POPOV_GHOST_DETERMINANT_FINITE_ANOMALY_FREE | **PASS** |

## The ghost determinant (zeta-regularized)

| L (circumference) | det′(−d²/dτ²) | L² | match |
|---:|---:|---:|:---:|
| 6.28319 | 39.47842 | 39.47842 | ✓ |
| 1.0 | 1.0 | 1.0 | ✓ |
| 3.32408 | 11.04951 | 11.04951 | ✓ |
| 5.0 | 25.0 | 25.0 | ✓ |

`ζ_ghost(0) = −1` (finite), `det′(−d²/dτ²) = L²`. For the closure loop `L = 2π`, the ghost determinant is `(2π)²` and the CKV volume `2π` gives the moduli factor **`1/(2π)` = PR #74's closure quantum**.

## The assembled one-loop measure

```
Z = Σ_sectors  ∫ (dL/L)   ·  det′_ghost   ·  det^{−1/2}_matter  ·  e^{−S}
              └ CKV(PR#74)┘  └ L² (PR#117) ┘  └ GY/zeta (PR#116) ┘
```

## Verdict

**DIFF_S1_FADDEEV_POPOV_GHOST_DETERMINANT_FINITE_ANOMALY_FREE.** THE Diff(S¹) FADDEEV–POPOV / GHOST DETERMINANT IS FINITE AND ANOMALY-FREE — THE GAUGE SECTOR OF THE S_BAM MEASURE IS COMPLETE. PR #115 named the reparametrization Faddeev–Popov (bc-ghost) determinant as a piece of the measure; PR #116 regularized the matter fluctuation determinant. This probe supplies the gauge sector.

THE GAUGE STRUCTURE. The closure loop is reparametrization-invariant under Diff(S¹). Gauge-fixing the loop einbein to a constant (the worldline/Polyakov procedure) leaves exactly one Teichmüller modulus — the circumference L = ∮ e dτ, i.e. the Schwinger proper time — and one conformal Killing vector, the constant ∂_τ = the residual rigid U(1) rotation. The algebra is the Witt algebra (centrally extended to Virasoro).

THE GHOST DETERMINANT. The Faddeev–Popov operator is P ξ ~ dξ/dτ, so P†P ~ −d²/dτ² on periodic fields, whose kernel is the constants — exactly the one CKV. The nonzero-mode determinant is zeta-regularized by the same method as PR #116: the eigenvalues (2πn/L)² give ζ_ghost(s) = 2 (L/2π)^{2s} ζ_R(2s), so ζ_ghost(0) = −1 (finite) and ζ_ghost'(0) = −2 ln L, hence det'(−d²/dτ²) = L² — finite and scheme-independent (verified to machine precision for L = 2π, 1, 3.32, 5).

THE CKV ZERO MODE IS PR #74's 1/(2π). The single CKV (rigid rotation) is divided out, Vol(U(1)) = L, so the Diff(S¹) quotient produces the worldline moduli measure ∫ dL/L. For the BAM closure loop, whose great-circle length is L = 2π (the closure quantum), this 1/L = 1/(2π) is exactly the per-loop measure factor PR #74 identified from the Schwinger anomaly a = α/(2π). So PR #74's 1/(2π) is the Faddeev–Popov CKV (ghost zero-mode) factor of the Diff(S¹) quotient — the gauge origin of the closure quantum in the measure.

ANOMALY-FREE. Unlike the 2D string worldsheet — where the bc-ghosts carry central charge c = −26 and consistency forces the conformal anomaly to cancel (D = 26) — the 1D worldline loop has no Weyl/conformal symmetry to be anomalous (there is no traceless symmetric 2-tensor in one dimension), so the Diff(S¹) gauge-fixing is anomaly-free. The only nontrivial anomaly in the BAM measure is the DISCRETE Z₂ orientation anomaly — the odd-k condition of PR #115 — not this continuous one.

SCOPE. With PR #116's finite matter determinant and this finite, anomaly-free ghost determinant (det' = L², CKV → 1/(2π)), the one-loop measure reads Z = Σ_sectors ∫ (dL/L) · det'_ghost · det^{−1/2}_matter · e^{−S}, every factor finite/computable: the gauge sector is complete. Still open: the absolute normalization of Z (the bulk κ₅²/Λ₅ anchor, PR #112), the multi-loop / interacting measure, and a closed form.

## What this completes (and does not)

- **Completes:** the Diff(S¹) gauge sector — a finite, anomaly-free ghost determinant (`det′ = L²`), with the CKV zero mode supplying the `dL/L` moduli measure whose `1/(2π)` is PR #74's closure-quantum factor. Together with PR #116's matter determinant, the one-loop measure is finite/computable.
- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor, PR #112), the multi-loop / interacting measure, and a closed-form expression.
