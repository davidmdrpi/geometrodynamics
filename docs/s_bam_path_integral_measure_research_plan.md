# The hard S_BAM path-integral measure: the full loop-measure construction (PR #115)

PR #74 (`s_bam_loop_measure_probe`) identified the per-loop-dimension
measure **factor** — the `1/(2π)` of the Schwinger anomaly `a = α/(2π)` — as
the BAM closure quantum, but explicitly flagged the "full covariant
`(2π)^d` path-integral derivation from `S_BAM`" as open. This sprint takes
up that hard open work: constructing the full path-integral **measure**
around that factor,

```
Z = Σ_sectors ∫ Dμ[X] e^{−S_BAM[X]},
```

in a loop-measure formalism — and reporting honestly which parts are
structurally fixed and which remain analytically open. (BAM's quantization
to date is otherwise saddle-point/instanton only: the bounce actions of
PRs #87–#90.)

## The arena: loop space

A BAM throat is characterised by its **closure loop** — a map `X: S¹ →`
base of the Hopf fibration carrying winding `k` and Hopf charge `c₁`. The
natural configuration space is therefore **loop space** `LS³` (with the
antipodal identification for the non-orientable sector), and the measure is
a Wiener-like measure on loops, gauge-fixed for the loop's symmetries:

```
Z = Σ_{k odd, c₁∈ℤ, n_part}  ∫_{LS³ / (Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂)}  Dμ[X]  e^{−S_BAM[X]},
       └──── closure ledger ────┘  └─────── gauge-fixed loop integral ───────┘
   with PR #74:  Dμ ~ Π dk/(2π)   (one closure quantum per loop dimension).
```

## What is structurally fixed (and computable)

  1. **Closure quantum = loop holonomy.** The Hopf/Wilson holonomy
     `e^{i∮A} = e^{ikχ}` is single-valued ⟹ the holonomy period is `2π` —
     the closure quantum (the PR #74 factor, here the measure period).
  2. **Sectors = closure ledger.** The measure's superselection sectors are
     the homotopy classes: winding `k` (`π₁` of the Hopf fibre), Hopf charge
     `c₁ ∈ π₃(S²) = ℤ`, and the partition `n_part`. `Σ_sectors` IS the
     ledger.
  3. **Odd-k = orientation anomaly condition.** The Möbius loop space is
     non-orientable; the integrand is a section of a twisted line bundle,
     anti-periodic under the antipodal `χ → χ+π`, so consistency requires
     `e^{ikπ} = −1 ⟹ k odd` (even `k` closes on the torus cover only). The
     kinematic odd-k lemma (`docs/odd_k_closure_lemma.md`) is thereby
     **upgraded** to the measure's `Z₂` orientation-anomaly-freedom
     condition.
  4. **Leading saddle = the bounce.** Stationary phase reproduces
     `e^{−S_BAM[saddle]}` — the bounce actions used in PRs #87–#90. The
     measure adds the prefactor (Faddeev–Popov × fluctuation determinant);
     the exponential, on which all prior results rest, is the leading term.

## The hard part: the one-loop measure factor

Reparametrization invariance `Diff(S¹)` must be gauge-fixed ⟹ a
Faddeev–Popov (`bc`-ghost) determinant; the measure factor is then
`(FP-det) × (fluctuation-det)`. The fluctuation operator is the second
variation of `S_BAM` about the throat saddle — the Tangherlini cavity
operator whose spectrum the cavity probes already compute. We check it: its
low spectrum is **positive** (`min ω² ≈ 1.11 > 0`), so the saddle is stable
and the Gaussian measure is real and well-defined pointwise (in the
tunneling sector the single Coleman negative mode supplies `Im Z` = the
decay rate; the zero modes are collective coordinates, traded for the moduli
measure — both standard).

## What remains analytically OPEN

The bare fluctuation determinant `Π_n ω_n` **diverges** — the log-det
partial sums grow without bound:

| modes N | `Σ ln ω_n` |
|---:|---:|
| 10 | 14.75 |
| 50 | 145.85 |
| 100 | 358.10 |
| 200 | 850.43 |
| 400 | 1964.17 |

— so it needs zeta/heat-kernel regularization, and the normalisation `Z` is
regularization-dependent and not yet rigorously constructed. So the
loop-measure formalism gives a **consistent structural definition** of the
`S_BAM` measure (sector sum × gauge-fixed loop integral, with PR #74's
`1/(2π)` the per-dimension factor, `2π` the loop holonomy, odd-k the
orientation anomaly, and the bounce as leading saddle), but the **hard
analytic core** — a finite, regularization-independent fluctuation
determinant / a constructive measure — stays open. Crucially the prior
saddle-point results are the leading `e^{−S}` term and do not depend on the
unresolved normalisation, so they stand.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal/status | PR #74 found `1/(2π)`; construct the full measure around it |
| T2 | arena | loop space `LS³ / (Diff S¹ ⋉ U(1) ⋉ Z₂)`; `Dμ ~ Π dk/(2π)` |
| T3 | holonomy | closure quantum = loop holonomy period `2π` |
| T4 | sectors | `Σ_sectors` = closure ledger (homotopy `k`, `c₁`, `n_part`) |
| T5 | odd-k anomaly | `e^{ikπ} = −1 ⟹ k odd` (Z₂ orientation anomaly) |
| T6 | one-loop | FP(`bc`-ghost) × fluctuation-det; spectrum positive (stable) |
| T7 | hard core | bare det diverges ⟹ needs regularization; `Z` not constructed |
| T8 | assessment | `S_BAM_PATH_INTEGRAL_MEASURE_STRUCTURALLY_DEFINED_ANALYTIC_CORE_OPEN` |

## Established and open

  - **Established (structural):** the `S_BAM` measure as a closure-ledger
    sector sum over a `Diff(S¹)`-gauge-fixed loop-space integral; PR #74's
    `1/(2π)` as the per-dimension factor; the closure quantum `2π` as the
    loop holonomy; the odd-k lemma upgraded to the `Z₂` orientation-anomaly
    condition; the PRs #87–#90 bounces as the leading saddle; the
    fluctuation operator stable.

  - **Open (analytic):** a finite, regularization-independent fluctuation
    determinant / a constructive measure — the bare determinant diverges, so
    the normalisation `Z` is not yet rigorously defined. Prior saddle-point
    results stand (leading `e^{−S}`, normalisation-independent).

## Cross-references

  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the `1/(2π)`
    per-loop-dimension measure factor this construction is built around.
  - `docs/odd_k_closure_lemma.md` — the odd-k lemma, here re-derived as the
    measure's `Z₂` orientation-anomaly condition.
  - `docs/throat_action_derivation_research_plan.md` /
    `docs/brane_tension_tuning_research_plan.md` — the `S_BAM` terms and the
    closure quantum `2π` this measure integrates.
  - `docs/majorana_bounce_action_research_plan.md` — PR #88, the bounce that
    is the leading saddle of this measure.

## Run

```
python -m experiments.closure_ledger.s_bam_path_integral_measure_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_s_bam_path_integral_measure_probe/`.
Expected verdict:
`S_BAM_PATH_INTEGRAL_MEASURE_STRUCTURALLY_DEFINED_ANALYTIC_CORE_OPEN`, 8/8 PASS.
