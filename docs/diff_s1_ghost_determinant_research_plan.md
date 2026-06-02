# Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117, revised)

PR #115 constructed the `S_BAM` measure as a closure-ledger sector sum over a
loop-space integral gauge-fixed by `Diff(S¹) ⋉ U(1)_Hopf ⋉ Z₂`, and named
the reparametrization Faddeev–Popov (`bc`-ghost) determinant as a piece to be
supplied. PR #116 regularized the **matter** fluctuation determinant
(Gel'fand–Yaglom + zeta, finite). This PR supplies the **gauge** sector: the
Faddeev–Popov / ghost determinant for the `Diff(S¹)` quotient — computed by
the same zeta method, **finite**, and **anomaly-free** for the 1D loop.

## Review correction (which determinant, and the L-power)

The first version wrote the ghost determinant as `det'(−d²/dτ²) = L²`. That is
**wrong by one square root**. The Faddeev–Popov determinant is the `bc`-ghost
path integral

```
Δ_FP = ∫ Db Dc e^{−∮ b (P c)} = det'(P),
```

with `P = d/dτ` the operator mapping the gauge parameter (the vector ghost
`c`) to the einbein variation. Because the `b`- and `c`-ghost spaces have
equal dimension here, `det'(P) = det'(P†P)^{1/2}`. Hence

```
Δ_FP = det'(P) = det'(P†P)^{1/2} = √(L²) = L      (NOT L²).
```

The `L²` is `det'(P†P) = det'(−d²/dτ²)`, the intermediate Laplacian
determinant (the **square**); the FP ghost determinant is its **square root**,
`L`. Of the three candidates `{det'(P), det'(P†P), det'(P†P)^{1/2}}`, the
ghost determinant is `det'(P†P)^{1/2}` ( `= det'(P)` in 1D, by the `±n` mode
pairing) `= L`.

| L | det′(P†P) | det′(P) | det′(P†P)^½ |
|---:|---:|---:|---:|
| 2π | 39.478 | 6.283 | 6.283 |
| 1 | 1.000 | 1.000 | 1.000 |
| 3.324 | 11.050 | 3.324 | 3.324 |
| 5 | 25.000 | 5.000 | 5.000 |

## The gauge structure (worldline reparametrization)

Gauge-fixing the loop einbein `e(τ)` to a constant leaves exactly:

  - **one Teichmüller modulus** `L = ∮ e dτ` — the loop circumference (= the
    Schwinger proper time), and
  - **one conformal Killing vector (CKV)** — the constant `∂_τ` (rigid `U(1)`
    rotation).

The FP operator is `P = d/dτ`; `P†P = −d²/dτ²`; `ker(P) =` constants `=` the
one CKV.

## The corrected one-loop measure (L-power fixed)

`Δ_FP = L` is the **Jacobian** of the einbein → proper-length gauge fixing:
it converts the einbein-fluctuation measure into the proper modulus measure
`dL`. Dividing by the conformal Killing volume `Vol(U(1)) = L` gives the
symmetry factor `1/L`. So the gauge sector yields the standard worldline
proper-time measure

```
Z = Σ_sectors ∫ (dL / L) · det^{−1/2}_matter · e^{−S}.
```

There is **no separate `det'_ghost = L²` factor** (the first version's error):
the ghost determinant is `L¹`, the Jacobian already used in writing the
modulus measure as the proper length `dL`. The `1/L` is the CKV factor.

## The CKV zero mode IS PR #74's 1/(2π) (unchanged by the correction)

The CKV (rigid rotation) volume is `Vol(U(1)) = L`. For the BAM closure loop,
whose great-circle length is `L = 2π` (the closure quantum), `1/L = 1/(2π)` —
exactly the per-loop measure factor PR #74 identified from `a = α/(2π)`. This
is the `c`-ghost zero-mode factor, independent of the determinant power, so
the correction above leaves it intact.

## Anomaly-free

The 1D worldline has no Weyl/conformal symmetry to be anomalous (no traceless
symmetric 2-tensor in one dimension), so the `Diff(S¹)` gauge-fixing is
anomaly-**free** — unlike the 2D string (`bc`-ghost `c = −26 ⟹ D = 26`). The
only nontrivial anomaly is the **discrete `Z₂` orientation anomaly** (odd-k,
PR #115).

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | gauge sector of the measure (revised per review) |
| T2 | gauge structure | `Diff(S¹)` → constant einbein ⟹ 1 modulus L + 1 CKV |
| T3 | ghost operator | `P = d/dτ`; `P†P = −d²/dτ²`; `ker(P)` = constants = 1 CKV |
| T4 | which determinant | `Δ_FP = det'(P) = det'(P†P)^{1/2} = L` (NOT `det'(P†P) = L²`) |
| T5 | corrected measure | `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}`; ghost L-power `L¹` |
| T6 | CKV = PR #74 | CKV `1/L`; for `L = 2π`, `1/(2π)` (unchanged) |
| T7 | anomaly-free | 1D worldline: no conformal anomaly (vs string `c = −26`) |
| T8 | assessment | `DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE` |

## Established and open

  - **Established (BAM-native):** the `Diff(S¹)` gauge sector — the FP ghost
    determinant `Δ_FP = det'(P) = det'(P†P)^{1/2} = L` (the square root of
    `det'(P†P) = L²`), the corrected measure `Z = Σ ∫ (dL/L)
    det^{−1/2}_matter e^{−S}` (ghost L-power `L¹`), anomaly-free, with the CKV
    `1/L = 1/(2π)` (PR #74) intact. With PR #116's matter determinant the
    one-loop measure is finite/computable.

  - **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor,
    PR #112); the multi-loop / interacting measure; a closed-form expression.

## Cross-references

  - `docs/s_bam_path_integral_measure_research_plan.md` — PR #115, which set
    up the `Diff(S¹)` gauge-fixing and named this ghost determinant.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116,
    the matter determinant alongside this ghost determinant.
  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the `1/(2π)` factor
    here re-derived as the ghost CKV zero-mode.
  - `docs/odd_k_closure_lemma.md` — the discrete `Z₂` orientation anomaly.

## Run

```
python -m experiments.closure_ledger.diff_s1_ghost_determinant_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_diff_s1_ghost_determinant_probe/`.
Expected verdict: `DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE`,
8/8 PASS.
