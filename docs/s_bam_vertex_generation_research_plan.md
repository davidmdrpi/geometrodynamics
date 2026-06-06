# S_BAM vertex generation (PR #140)

> **Framing (to avoid a category error).** `S_BAM` is the action of the
> **matter** field theory on the *fixed classical* throat geometry — this is
> QFT-on-a-curved-GR-background, **not** a quantum theory of gravity. BAM derives
> QFT *from* continuous GR (geometry → fields); it does not quantise the metric.
> "Generating the vertices from S_BAM" means reading the matter self-interaction
> off the classical geometry, the opposite direction from quantum gravity.

PRs #137–#139 built the antipodal matter interaction with cubic and quartic
vertices but flagged, each time, the same open item: the vertices were
**modelled**, not derived — "whether the S_BAM measure (#115–#122) actually
generates the cubic/quartic term" was left open. This PR closes that item
**structurally**: it shows the vertices are the Taylor coefficients of the S_BAM
action expanded about the throat background, that their `Σl`-even selection rule
(#137/#138) is the **antipodal symmetry** of S_BAM (a Ward identity), and that
the positive quartic sign (#138) is the **measure-consistency condition** (#122).
What stays inherited is only the coupling **magnitudes** (the action form and the
#133 scale).

## Vertices = Taylor coefficients of the action

Expanding S_BAM about the classical throat background `φ_cl` in the fluctuation
`φ`,

```
S_BAM[φ_cl + φ] = S_cl + S_2[φ] + S_3[φ] + S_4[φ] + …,   S_n = (1/n!) ∫ (δⁿS/δφⁿ) φⁿ,
```

  - **S_2** is the quadratic fluctuation action — the #116 Tangherlini
    determinant / the #135 free propagator;
  - **S_3 = (g/3!) ∫ φ³ √g** generates the cubic vertex `g_{knm} = ∫ ψ_k ψ_n ψ_m`
    (#137);
  - **S_4 = (λ/4!) ∫ φ⁴ √g** generates the quartic vertex `g_4 = ∫ ψ⁴` (#138).

So the vertices are not added by hand: they are the higher functional
derivatives of the S_BAM action. A geometric (non-quadratic) S_BAM generates the
whole tower; a free (purely quadratic) action would have none.

## The Σl-even selection rule is the antipodal symmetry of S_BAM (a Ward identity)

The S_BAM measure carries the loop quotient `Diff(S¹) ⋉ U(1)_Hopf ⋉ Z₂` (#74),
whose `Z₂` is the antipodal map `A : x → −x` (the throat ↔ antithroat C-swap #63,
#128). Under `A`, a fluctuation mode of angular momentum `l` carries the harmonic
parity, so its amplitude transforms as

```
a_{lm} → (−1)^l a_{lm} .
```

Therefore a vertex of `n` modes picks up `∏_i (−1)^{l_i} = (−1)^{Σl}`, and is
`A`-invariant only if `Σl` is even. Because S_BAM is `A`-invariant (`A` is a
gauge symmetry of the measure), every generated vertex must have `Σl` even —
exactly the selection rule #137/#138 found.

| l | `a_l →` | `(−1)^l` |
|---:|---:|---:|
| 0 | `+a_l` | `+1` |
| 1 | `−a_l` | `−1` |
| 2 | `+a_l` | `+1` |
| 3 | `−a_l` | `−1` |

| vertex modes | Σl | A-invariant? | S³ integral ≠ 0? |
|---|---:|:---:|:---:|
| (0,0,0) | 0 | ✓ | ✓ |
| (1,1,0) | 2 | ✓ | ✓ |
| (1,1,1) | 3 | ✗ | ✗ |
| (1,1,2) | 4 | ✓ | ✓ |
| (1,1,1,1) | 4 | ✓ | ✓ |
| (1,1,1,0) | 3 | ✗ | ✗ |

The selection rule is a **Ward identity** of the antipodal symmetry, not a
modelling choice (the `A`-invariance condition coincides with the explicit S³
vanishing of odd-`Σl` integrals).

## The positive quartic sign is the measure-consistency condition

The S_BAM measure `∫ Dμ e^{−S}` exists — it is reflection-positive and gives the
unitary kernel (#135), and it converges non-perturbatively (#122) — only if the
action is bounded below, which (for a potential-type quartic) requires the
quartic coupling positive. So the positive sign of the quartic (#138) is not a
free choice: it is fixed by the measure's existence. The geometric overlap
`∫ψ⁴ > 0` (#138) realises it.

## What stays inherited

The coupling **magnitudes** (`g, λ`) are the numerical values of the higher
functional derivatives of S_BAM — set by the specific action form and the overall
S_BAM normalisation, which carries the `κ₅²/Λ₅` bulk scale (#133). So the
structure (the vertices' existence, their `Σl`-even selection, the positive
quartic sign) is derived from the action's symmetry and the measure's
consistency; the magnitudes inherit the standing #133 scale residual.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | derive the S_BAM vertices (close the #137–#139 open item) |
| T2 | Taylor coefficients | vertices = `S_n` of S_BAM (S_2/#116, S_3/#137, S_4/#138) |
| T3 | non-quadratic required | a free action has no vertices; S_BAM is geometric |
| T4 | Ward identity | `Σl`-even = the antipodal Z₂ Ward identity (`a_l → (−1)^l a_l`) |
| T5 | quartic sign | positive = measure consistency (#122) |
| T6 | ledger | derived (existence/selection/sign) vs inherited (magnitudes/#133) |
| T7 | scope | structure derived; action form / magnitudes / #133 open |
| T8 | assessment | `S_BAM_VERTICES_FROM_ACTION_EXPANSION_Z2_WARD_SELECTION_MEASURE_POSITIVITY` |

## Established and open

  - **Established (BAM-native):** the antipodal matter vertices are the Taylor
    coefficients of the S_BAM action; their `Σl`-even selection rule (#137/#138)
    is the antipodal `Z₂` Ward identity of S_BAM (the `Z₂` of the #74 loop
    quotient), and the positive quartic sign (#138) is the measure-consistency
    condition (#122). The vertex structure is *generated*, not modelled.

  - **Does not / open:** the exact S_BAM functional form (Polyakov/Nambu-Goto-type
    choice) and the coupling **magnitudes** (which carry the `κ₅²/Λ₅` scale, #133)
    are not fixed; higher vertices follow the same `Σl`-even Ward identity. The
    bulk-scale (#133) and flavor (#134) residuals stand.

## Cross-references

  - `docs/cubic_vertex_ledger_research_plan.md` /
    `docs/quartic_vertex_bounded_interaction_research_plan.md` — #137/#138, the
    vertices generated here.
  - `docs/antipodal_matter_interaction_synthesis_research_plan.md` — #139, the
    arc synthesis that flagged this open item.
  - `docs/s_bam_path_integral_measure_research_plan.md` /
    `docs/bam_factorized_sector_sum_research_plan.md` — #115/#122, the S_BAM
    measure and its convergence.
  - `docs/charge_conjugation_swap_research_plan.md` — #63, the C-swap (the Z₂).
  - `docs/bulk_scale_ledger_research_plan.md` — #133, the κ₅²/Λ₅ scale the
    magnitudes inherit.

## Run

```
python -m experiments.closure_ledger.s_bam_vertex_generation_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_s_bam_vertex_generation_probe/`.
Expected verdict:
`S_BAM_VERTICES_FROM_ACTION_EXPANSION_Z2_WARD_SELECTION_MEASURE_POSITIVITY`, 8/8 PASS.
