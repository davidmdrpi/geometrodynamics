# The factorized BAM sector sum Z (PR #122)

PRs #74 and #115–#121 built and validated every ingredient of the BAM loop
measure. This PR **assembles** them into the full factorized sector sum

```
Z = Σ_{k odd, c₁∈ℤ, n_part}  (−1)^k  ∫₀^∞ (dL/L)  det^{−1/2}_matter(L)
                                      · e^{i(π/2)(1−2a)} · e^{−S_BAM},
```

and shows it **factorizes** into a discrete (topological, Z₂-signed) sum
times a continuous (analytic, η-phased) moduli integral — with a concrete
payoff: the Z₂ grading cancels the leading UV divergence between the
orientable and Möbius sectors.

## The assembled measure (each factor validated)

| factor | meaning | source |
|---|---|---|
| `Σ_{k,c₁,n_part}` | closure-ledger sector sum | #115 |
| `(−1)^k` | discrete Z₂ orientation sign | #115/#118/#121 |
| `∫ (dL/L)` | gauge-fixed moduli measure (CKV = closure quantum `1/L`) | #74/#117/#118 |
| `det^{−1/2}_matter(L)` | matter fluctuation det (finite, Gel'fand–Yaglom) | #116 |
| `det'(P) = L` (ghost) | first-order ghost det | #117/#118 |
| `e^{i(π/2)(1−2a)}` | continuous η-phase (holonomy `a`) | #119/#121 |
| `e^{−S_BAM}` | leading bounce saddle | #87–#90 |

## Why it factorizes

The discrete Z₂ orientation sign `(−1)^k` is a **sector-constant** — it
depends only on the winding parity, not on the continuous moduli `L` or
holonomy `a` — so it pulls **out** of the continuous integral:

```
Z = Σ_{discrete sectors} (−1)^k × [ ∫(dL/L) · det · η-phase · e^{−S} ].
```

So `Z` is a discrete Z₂-signed sum of continuous moduli integrals. The
continuous part carries the η-phase (confined to the right half-circle,
PR #121) and the finite determinants; the discrete part carries the
orientation signs. No double-counting (PR #121: the U(1)-valued η-phase and
the Z₂ sign are independent).

## The Z₂-graded UV cancellation (the payoff)

Grouped by orientation, `Z` is the Z₂-graded combination of the orientable
(periodic, +) and Möbius (antiperiodic, −) contributions. The leading
heat-kernel (Weyl) coefficient `a_{−1/2} = L/√(4π)` is a **bulk** quantity,
**independent of the boundary condition**, so it is identical in both
sectors and **cancels** in the graded difference. Concretely the heat traces
`θ_per(t) = Σ_n e^{−(2πn/L)²t}` and `θ_anti(t) = Σ_n e^{−(2π(n+½)/L)²t}` each
diverge as `L/√(4πt)` as `t → 0`, but

```
θ_per(t) − θ_anti(t)  ~  e^{−π²/t}  →  0
```

is **UV-finite** (the polynomial UV divergence cancels to all orders,
leaving only exponentially small instanton terms). So the orientation Z₂
grading renders the bulk UV of the sector sum finite.

| t | θ_per | θ_anti | Weyl L/√(4πt) | θ_per − θ_anti |
|---:|---:|---:|---:|---:|
| 0.2 | 3.96333 | 3.96333 | 3.96333 | 0 |
| 0.1 | 5.60499 | 5.60499 | 5.60499 | ~1e-15 |
| 0.05 | 7.92666 | 7.92666 | 7.92666 | 0 |
| 0.02 | 12.53314 | 12.53314 | 12.53314 | ~2e-15 |

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | assemble factorized Z from PRs #74, #115–#121 |
| T2 | formula | the assembled measure (each factor + source PR) |
| T3 | discrete sum | closure ledger; `(−1)^k` sector-constant pulls out |
| T4 | continuous integral | `dL/L` + finite dets + η-phase; IR-convergent |
| T5 | factorization | `Z = (discrete Z₂ sum) ⊗ (continuous integral)`; no double-count |
| T6 | graded UV | `θ_per − θ_anti ~ e^{−π²/t} → 0` (Weyl term cancels) |
| T7 | scope | assembled; absolute normalization / non-pert sum open |
| T8 | assessment | `BAM_FACTORIZED_SECTOR_SUM_Z_DISCRETE_Z2_TIMES_CONTINUOUS_ETA_GRADED_UV_CANCELS` |

## Established and open

  - **Established (BAM-native):** the full factorized one-loop sector sum
    `Z` — a discrete Z₂-signed (topological) sum of continuous η-phased
    (analytic) moduli integrals, every factor finite/validated (PRs #74,
    #116–#121), with the Z₂ grading cancelling the leading UV
    (`θ_per − θ_anti ~ e^{−π²/t} → 0`).

  - **Open:** the absolute normalization (the bulk `κ₅²/Λ₅` anchor,
    PR #112); the full non-perturbative convergence of the sector sum; the
    multi-loop measure. The assembly organizes the structure; it does not fix
    the overall scale.

## Cross-references

  - `docs/s_bam_path_integral_measure_research_plan.md` — PR #115, the
    measure structure assembled here.
  - `docs/bam_sector_phase_ledger_research_plan.md` — PR #121, the
    continuous-η / discrete-Z₂ separation (no double-counting).
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` /
    `docs/diff_s1_first_order_ghost_audit_research_plan.md` — PRs #116/#118,
    the finite matter and ghost determinants.
  - `docs/lattice_validation_research_plan.md` — PR #120, the lattice
    validation of all the pieces.

## Run

```
python -m experiments.closure_ledger.bam_factorized_sector_sum_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_bam_factorized_sector_sum_probe/`.
Expected verdict:
`BAM_FACTORIZED_SECTOR_SUM_Z_DISCRETE_Z2_TIMES_CONTINUOUS_ETA_GRADED_UV_CANCELS`,
8/8 PASS.
