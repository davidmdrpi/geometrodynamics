# Alpha normalization ledger for the gauge–matter coupling (PR #143)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The gauge coupling is
> the EM coupling on the classical antipodal throat; this ledger audits its
> normalisation `α`.

PRs #141 and #142 derived the gauge–matter coupling **structure** at the
antipodal throat (minimal coupling, the Σl-even vertex, current conservation, the
Ward identity, photon masslessness) but left, each time, the same flagged input:
the coupling **strength** `α` (the "137 problem", #105). This PR is the
consolidating **ledger** for `α` — parallel to the bulk-scale ledger (#133) for
`κ₅²/Λ₅`. It separates what the geometry **derives** about the EM coupling (the
charge quantum, the `1/2π` loop measure, the coupling structure, the running)
from the one irreducible **input**: the **value** `α ≈ 1/137`.

## How α enters

Every EM observable is a function of `α` and the geometry:

  - the EM amplitude `A_EM = α·ℏc/2` (#105);
  - the gauge–matter vertex strength `∝ c₁² α` (the #141 coupling squared);
  - the Schwinger anomaly `a = α/2π` (#74) — the one-loop `g−2`.

So `α` is a **single** dimensionless number feeding the whole EM sector.

## The charge QUANTUM is derived (geometric, integer)

The Hopf charge is the integer Hopf number, `|c₁| = 1` (#58/#74) — charge
quantisation is topological, not input. The charge **unit** is fixed by the
geometry; only the coupling **strength** (how strongly that unit charge couples)
is `α`.

## The 1/2π loop-measure factor is derived (the closure quantum)

In the Schwinger anomaly `a = α/2π`, the `2π` is the BAM closure-quantum loop
measure (#74) — derived. So of the famous `g−2` result, the geometry fixes the
`1/2π` and leaves only `α` as the input prefactor: **BAM derives the measure, not
the coupling.**

## The running is derived; the value is not

The RG flow of `α` — the vacuum polarisation, transverse by the #142 Ward
identity — is derived structurally (BAM derives **how** `α` runs). The boundary
value `α(μ_0) ≈ 1/137` is the input (BAM does not derive **where** it starts).
This is the #105 classification, sharpened: the running derived, the value
residual.

## The value α ≈ 1/137 is the one EM input residual (the 137 problem)

A fit-independent scan of `α⁻¹ = 137.036` against the BAM closure numbers (`2π`,
`k₅`, `β_lepton = 50π`, …) finds **no clean match**:

| candidate | value | % off | needs ad-hoc term? |
|---|---:|---:|:---:|
| 2π | 6.28 | −95% | — |
| β_lepton = 50π | 157.08 | +15% | — |
| k₅³ + 2π | 131.28 | −4.2% | — |
| 50π − 20 | 137.08 | +0.0% | ✗ (fit) |
| 4·k₅² + 37 | 137.0 | −0.0% | ✗ (fit) |
| 8π·k₅ | 125.66 | −8.3% | — |

The sub-% near-misses (`50π − 20`, `4·k₅² + 37`) each require an ad-hoc additive
`O(20–37)` integer — fits, not derivations, the same reverse-engineering failure
mode #107/#108 documented for `√σ/m_e`. So `α` is plausibly **irreducible**, like
`√σ/m_e` (#108): the EM sector contributes exactly **one** dimensionless
residual, the value `α`.

## Tie to the input budget

`α` joins the program's handful of dimensionless residuals — `{n_part, √σ/m_e,
ε, α}` (#104/#108) — as the EM one. The charge quantum, the `1/2π` measure, the
coupling structure, and the running are derived; the value `α` is the single EM
input, the 137 problem.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | α normalization ledger for the gauge–matter coupling |
| T2 | how α enters | `A_EM = α·ℏc/2`, vertex `∝ c₁²α`, `a = α/2π` — one number |
| T3 | charge quantum derived | `|c₁| = 1` (integer Hopf number, #58/#74) |
| T4 | 1/2π measure derived | `a = α/2π`, the `2π` = closure quantum (#74) |
| T5 | running derived, value not | RG flow (#142) derived; `α(μ_0)` input (#105) |
| T6 | value = one residual | no clean closure match (137 problem, #108) |
| T7 | ledger / input budget | derived (quantum/measure/structure/running) vs input (value α) |
| T8 | assessment | `ALPHA_NORMALIZATION_LEDGER_..._VALUE_ONE_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** the geometry derives the EM coupling's charge
    quantum (`|c₁| = 1`), its `1/2π` loop measure (the closure quantum), the
    coupling structure (#141/#142), and the running of `α`; the **value** `α ≈
    1/137` is the one EM input residual (the 137 problem), with no clean closure
    match — plausibly irreducible like `√σ/m_e` (#108).

  - **Does not / open:** the ledger does **not** derive `α` (the 137 problem
    stays open) or fix the EM normalisation absolutely. The `α` (#105/#108),
    bulk-scale (#133), and flavor (#134) residuals stand.

## Cross-references

  - `docs/gauge_matter_coupling_research_plan.md` /
    `docs/gauge_ward_identity_research_plan.md` — #141/#142, the coupling
    structure and the running this ledger normalises.
  - `docs/alpha_G_ledger_classification_research_plan.md` — #105, `α` as the
    universal residual (value input, running derived).
  - `docs/lepton_qcd_ratio_legitimate_search_research_plan.md` /
    `docs/ratio_832_npart_recycling_research_plan.md` — #108/#107, the
    fit-independent search methodology (ad-hoc-term failure mode).
  - `docs/bulk_scale_ledger_research_plan.md` — #133, the parallel bulk-scale
    ledger.
  - `docs/program_synthesis_research_plan.md` — #104, the input budget.

## Run

```
python -m experiments.closure_ledger.alpha_normalization_ledger_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_alpha_normalization_ledger_probe/`.
Expected verdict:
`ALPHA_NORMALIZATION_LEDGER_CHARGE_QUANTUM_AND_2PI_MEASURE_DERIVED_VALUE_ONE_RESIDUAL`, 8/8 PASS.
