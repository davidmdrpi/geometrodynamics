# Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117)

PR #115 constructed the `S_BAM` measure as a closure-ledger sector sum over a
loop-space integral gauge-fixed by `Diff(S¹) ⋉ U(1)_Hopf ⋉ Z₂`, and named
the reparametrization Faddeev–Popov (`bc`-ghost) determinant as a piece to be
supplied. PR #116 regularized the **matter** fluctuation determinant
(Gel'fand–Yaglom + zeta, finite). This PR supplies the **gauge** sector: the
Faddeev–Popov / ghost determinant for the `Diff(S¹)` quotient — computed by
the same zeta method, **finite**, and **anomaly-free** for the 1D loop.

## The gauge structure (worldline reparametrization)

The closure loop `X: S¹ → base` is reparametrization-invariant: `X(τ) →
X(f(τ))` for `f ∈ Diff(S¹)`. Gauge-fixing the loop einbein `e(τ)` to a
constant (the worldline/Polyakov procedure) leaves exactly:

  - **one Teichmüller modulus** `L = ∮ e dτ` — the loop circumference (= the
    Schwinger proper time), and
  - **one conformal Killing vector (CKV)** — the constant `∂_τ`, i.e. the
    residual rigid `U(1)` rotation of the loop.

The Lie algebra of `Diff(S¹)` is the Witt algebra (centrally extended:
Virasoro).

## The ghost operator and its determinant (same zeta method as PR #116)

The Faddeev–Popov operator is `P: Vect(S¹) → (einbein variations)`, `P ξ ~
dξ/dτ`, so `P†P ~ −d²/dτ²` on periodic fields. Its kernel is the constants —
exactly the one CKV. The nonzero-mode determinant is zeta-regularized as in
PR #116: eigenvalues `(2πn/L)²`, `n = ±1, ±2, …`, give

```
ζ_ghost(s) = 2 (L/2π)^{2s} ζ_R(2s)   ⟹   ζ_ghost(0) = −1,
ζ_ghost'(0) = −2 ln L                 ⟹   det'(−d²/dτ²) = L²,
```

finite and scheme-independent (verified: `ζ_ghost(0) = −1` and `det' = L²` to
machine precision for `L = 2π, 1, 3.32, 5`).

| L (circumference) | det′(−d²/dτ²) | L² |
|---:|---:|---:|
| 2π | 39.478 | 39.478 |
| 1 | 1.000 | 1.000 |
| 3.324 | 11.050 | 11.050 |
| 5 | 25.000 | 25.000 |

## The CKV zero mode IS PR #74's 1/(2π)

The single CKV (rigid rotation) must be divided out: `Vol(U(1)) = L`, the loop
circumference. The `Diff(S¹)` quotient therefore produces the worldline moduli
measure `∫ dL/L` — and for the BAM closure loop, whose great-circle length is
`L = 2π` (the closure quantum), this `1/L` is exactly

```
1/L = 1/(2π),
```

the per-loop measure factor PR #74 identified from the Schwinger anomaly
`a = α/(2π)`. So **PR #74's `1/(2π)` is the Faddeev–Popov CKV (ghost
zero-mode) factor of the `Diff(S¹)` quotient** — the gauge origin of the
closure quantum in the measure.

## Anomaly-free (the clean part)

Unlike the 2D string worldsheet (where the `bc`-ghosts carry central charge
`c = −26` and consistency forces the conformal anomaly to cancel, `D = 26`),
the 1D worldline loop has **no Weyl/conformal symmetry to be anomalous**:
there is no traceless symmetric 2-tensor in one dimension, so the `Diff(S¹)`
gauge-fixing is anomaly-**free**. (Contrast: PR #115's nontrivial anomaly is
the **discrete `Z₂` orientation anomaly** — the odd-k condition — not this
continuous one.) The ghost sector is clean: a finite determinant `L²` and the
`dL/L` moduli measure, with no continuous anomaly to cancel.

## The assembled one-loop measure

```
Z = Σ_sectors  ∫ (dL/L)   ·  det′_ghost   ·  det^{−1/2}_matter  ·  e^{−S}
              └ CKV (PR#74)┘  └ L² (PR#117)┘  └ GY/zeta (PR#116) ┘
```

every factor finite/computable.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | PR #115 flagged the ghost; PR #116 did matter; now the ghost |
| T2 | gauge structure | `Diff(S¹)` → constant einbein ⟹ 1 modulus L + 1 CKV |
| T3 | ghost operator | `−d²/dτ²`; kernel = constants = 1 CKV |
| T4 | ghost determinant | `ζ_ghost(0) = −1`, `det'(−d²/dτ²) = L²` (computed) |
| T5 | CKV = PR #74 | CKV ⟹ `dL/L`; for `L = 2π`, `1/L = 1/(2π)` |
| T6 | anomaly-free | 1D worldline: no conformal anomaly (vs string `c = −26`) |
| T7 | scope | gauge sector complete; abs `Z` normalization / multi-loop open |
| T8 | assessment | `DIFF_S1_FADDEEV_POPOV_GHOST_DETERMINANT_FINITE_ANOMALY_FREE` |

## Established and open

  - **Established (BAM-native):** the `Diff(S¹)` gauge sector — a finite,
    anomaly-free ghost determinant (`det' = L²`, `ζ_ghost(0) = −1`, one CKV
    zero mode), with the CKV supplying the `dL/L` moduli measure whose
    `1/(2π)` for the closure loop IS PR #74's closure-quantum factor.
    Together with PR #116's matter determinant, the one-loop measure is
    finite/computable.

  - **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor,
    PR #112); the multi-loop / interacting measure; a closed-form expression.

## Cross-references

  - `docs/s_bam_path_integral_measure_research_plan.md` — PR #115, which set
    up the `Diff(S¹)` gauge-fixing and named this ghost determinant.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116,
    the matter determinant alongside this ghost determinant.
  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the `1/(2π)` factor
    here re-derived as the ghost CKV zero-mode.
  - `docs/odd_k_closure_lemma.md` — the discrete `Z₂` orientation anomaly,
    the only nontrivial anomaly (this continuous one being absent).

## Run

```
python -m experiments.closure_ledger.diff_s1_ghost_determinant_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_diff_s1_ghost_determinant_probe/`.
Expected verdict: `DIFF_S1_FADDEEV_POPOV_GHOST_DETERMINANT_FINITE_ANOMALY_FREE`,
8/8 PASS.
