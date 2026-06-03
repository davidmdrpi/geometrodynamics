# The APS quark partition index from the factorized sector sum (PR #123)

PR #122 assembled the BAM loop measure into the factorized sector sum
`Z = Σ_{k odd, c₁, n_part} (−1)^k ∫(dL/L) det^{−1/2}_matter e^{i(π/2)(1−2a)}
e^{−S_BAM}` — a discrete Z₂-signed (topological) sum × a continuous η-phased
(analytic) integral. A Z₂-graded partition sum has a **Witten / APS index**:
the graded trace `Tr(−1)^k` is a topological invariant, computed via the
Atiyah–Patodi–Singer η-invariant boundary term. This PR extracts that index
for the **quark** sector and shows what it does — and does not — fix about
the quark partition `n_part`.

## The index from the Z₂-grading

The orientation sign `(−1)^k` of the factorized sum is a Z₂ grading, so

```
I = Tr (−1)^k e^{−βH}      (the Witten / APS index)
```

is a topological invariant: the orientable (+) and Möbius (−) nonzero modes
pair, leaving only the net graded content. `Z(β)` is β-dependent, but `I` is
not.

## The APS η-invariant boundary term

For the first-order closure operator `∂_τ` with holonomy `a` (eigenvalues
`2πi(n+a)/L`), the Atiyah–Patodi–Singer ξ-invariant — the η-boundary
correction to the index — is

```
ξ(a) = (η_A(0) + h)/2 = 1/2 − a = ζ_H(0,a)      (twisted, h = 0),
```

continuous in `a` (`ξ(0)=1/2`, `ξ(1/2)=0`, `ξ(1⁻)=−1/2`). This is the
η-invariant machinery of PRs #119–#121.

| a | η(a) = 1−2a | ξ(a) = 1/2 − a |
|---:|---:|---:|
| →0 | +1 | +0.5 |
| 1/4 | +0.5 | +0.25 |
| 1/2 | 0 | 0 |
| 3/4 | −0.5 | −0.25 |
| →1 | −1 | −0.5 |

## The index is an integer (spectral flow)

As the holonomy winds once, `a : 0 → 1`, exactly one eigenvalue `2π(n+a)/L`
crosses zero (the `n = 0` mode), so the APS index — the spectral flow — is

```
spectral flow = ξ(0⁺) − ξ(1⁻) = 1/2 − (−1/2) = 1      (an INTEGER).
```

The fractional `ξ` is the boundary η-correction; the spectral flow is the
integer topological index — the discrete Z₂ content of the sector-phase
ledger (PR #121), made into an index.

## Apply to the quark partition

The quark sector's closure count (PR #76/#97) is `N_q = 2·n_part = 466`
(`n_part = 233`, `N_lepton = 4 k₅² = 100`). The factor of 2 — the **even
doubling** `N_q = 2·n_part` — is precisely the Z₂-graded structure: the
orientation index pairs the modes, doubling the count. So the APS index
extracts the **topological (mod-2 / doubling) content** of the quark
partition, which is §8-stable.

## Topological vs residual (the honest split)

The APS index formalises exactly the split PRs #97/#107 found empirically:

  - **§8-stable (topological):** the index / the doubling `N_q = 2·n_part`
    (even across all twelve `quark_axioms` §8 ablations) and the spectral-flow
    integer — the Z₂ mod-2 content, the Atiyah–Patodi–Singer topological
    invariant, protected against the §8 ablations.
  - **Residual (non-topological):** the bare value `n_part` itself — the
    continuous, ξ-type (η-boundary) part — which drifts `216–255` across the
    §8 ablations (the phenomenological compensator, PR #97/#107).

So the APS quark partition index **derives the topological structure** (the
graded doubling, the integer spectral flow) but **not the value of `n_part`**,
exactly as the compensator status requires.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | derive the APS quark partition index from the factorized sum |
| T2 | index | `I = Tr(−1)^k` from the Z₂ grading (topological, β-indep) |
| T3 | APS ξ | `ξ(a) = (η+h)/2 = 1/2 − a` (η-boundary correction) |
| T4 | spectral flow | integer index `= ξ(0⁺) − ξ(1⁻) = 1` |
| T5 | quark partition | `N_q = 2·n_part` (the Z₂-graded doubling) |
| T6 | topological vs residual | doubling/index §8-stable; `n_part` value drifts |
| T7 | scope | index fixes the structure, not the `n_part` value |
| T8 | assessment | `APS_QUARK_PARTITION_INDEX_TOPOLOGICAL_DOUBLING_STABLE_NPART_VALUE_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** the factorized sector sum's Z₂ grading
    defines a Witten/APS index — a topological spectral-flow integer (`= 1`
    per holonomy cycle), with the APS ξ-invariant `ξ(a) = 1/2 − a` the
    η-boundary correction. For the quark sector it extracts the §8-stable
    **topological content** — the graded doubling `N_q = 2·n_part` (even) —
    formalising the PR #97/#107 split.

  - **Open / residual:** the bare value `n_part = 233` — the non-topological
    residual, the §8-drifting compensator (PR #97/#107). The index derives
    the structure, not the value.

## Cross-references

  - `docs/bam_factorized_sector_sum_research_plan.md` — PR #122, the
    factorized sector sum this index is read off.
  - `docs/bam_sector_phase_ledger_research_plan.md` /
    `docs/detprime_dtau_eta_invariant_phase_research_plan.md` — PR #121/#119,
    the Z₂ / η machinery.
  - `docs/ratio_832_npart_recycling_research_plan.md` /
    `docs/npart_dynamical_hierarchy_research_plan.md` — PR #107/#97, the
    `n_part` compensator and its §8 drift.

## Run

```
python -m experiments.closure_ledger.aps_quark_partition_index_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_aps_quark_partition_index_probe/`.
Expected verdict:
`APS_QUARK_PARTITION_INDEX_TOPOLOGICAL_DOUBLING_STABLE_NPART_VALUE_RESIDUAL`,
8/8 PASS.
