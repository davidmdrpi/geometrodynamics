# The hard S_BAM path-integral measure: the full loop-measure construction (PR #115)

**Run:** 2026-05-31T19:53:22+00:00

PR #74 found the per-loop-dimension measure FACTOR `1/(2π)` (the Schwinger `a = α/(2π)`) and flagged the full covariant path-integral measure as open. This sprint takes up that hard work: constructing `Z = Σ_sectors ∫ Dμ[X] e^{−S_BAM[X]}` around it, in a loop-measure formalism. **Result: the measure is defined STRUCTURALLY (sector sum × gauge-fixed loop integral), but the hard analytic core (a finite fluctuation determinant) stays OPEN.**

- **Arena**: loop space LS³ / (Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂); Σ_sectors = the closure ledger
- **Closure quantum**: 2π = the Hopf/Wilson loop holonomy period (PR #74 factor)
- **Odd-k**: orientation anomaly condition e^{ikπ} = −1 ⟹ k odd
- **One-loop**: FP(bc-ghost) × fluctuation-det; fluctuation spectrum positive (stable saddle)
- **Open**: bare determinant diverges ⟹ needs zeta/heat-kernel reg; Z not rigorously constructed
- **Saddle results**: leading e^{−S} term (PRs #87–#90) — unaffected by the normalisation

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_and_status` | PR #74 found 1/(2π); construct the full measure around it | **PASS** |
| T2 | `T2_loop_space_arena_and_formal_measure` | arena = loop space LS³ / (Diff S¹ ⋉ U(1) ⋉ Z₂) | **PASS** |
| T3 | `T3_closure_quantum_is_loop_holonomy` | closure quantum = loop holonomy period 2π | **PASS** |
| T4 | `T4_sectors_are_the_closure_ledger` | sectors = closure ledger (homotopy k, c₁, n_part) | **PASS** |
| T5 | `T5_odd_k_is_orientation_anomaly_condition` | odd-k = Z₂ orientation anomaly: e^{ikπ} = −1 ⟹ k odd | **PASS** |
| T6 | `T6_gauge_fixing_and_one_loop_factor` | FP(bc-ghost) × fluctuation-det; spectrum positive (stable) | **PASS** |
| T7 | `T7_hard_analytic_core_open` | bare det diverges ⟹ needs reg; Z not rigorously constructed | **PASS** |
| T8 | `T8_assessment` | S_BAM_PATH_INTEGRAL_MEASURE_STRUCTURALLY_DEFINED_ANALYTIC_CORE_OPEN | **PASS** |

## The formal measure

```
Z = Σ_{k odd, c₁∈ℤ, n_part}  ∫_{LS³ / (Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂)}  Dμ[X]  e^{−S_BAM[X]}
     └── closure ledger ──┘   └────── gauge-fixed loop integral ──────┘
   with PR #74:  Dμ ~ Π dk/(2π)   (one closure quantum per loop dimension)
```

## T5: odd-k as the orientation-anomaly condition

| k | e^{ikπ} | closes? |
|---|---:|---|
| 1 | -1 | Möbius (non-orientable) |
| 2 | +1 | torus cover only |
| 3 | -1 | Möbius (non-orientable) |
| 4 | +1 | torus cover only |
| 5 | -1 | Möbius (non-orientable) |
| 6 | +1 | torus cover only |

The integrand is a twisted-bundle section, anti-periodic under the antipodal `χ → χ+π`, so `e^{ikπ} = −1` ⟹ **k odd**. The kinematic odd-k lemma is the measure's Z₂ orientation-anomaly condition.

## T7: the hard analytic core (open)

| modes N | log-det partial sum Σ ln ω_n |
|---:|---:|
| 10 | 14.75 |
| 50 | 145.85 |
| 100 | 358.1 |
| 200 | 850.43 |
| 400 | 1964.17 |

The fluctuation spectrum is positive (stable saddle), but the **bare determinant diverges** — the log-det grows without bound ⟹ zeta/heat-kernel regularization is needed and `Z` is regularization-dependent, not yet rigorously constructed. The prior saddle results are the leading `e^{−S}` term and are unaffected.

## Verdict

**S_BAM_PATH_INTEGRAL_MEASURE_STRUCTURALLY_DEFINED_ANALYTIC_CORE_OPEN.** THE LOOP-MEASURE FORMALISM CONSTRUCTS THE S_BAM PATH-INTEGRAL MEASURE STRUCTURALLY; THE HARD ANALYTIC CORE REMAINS OPEN. PR #74 identified the per-loop-dimension measure factor 1/(2π) (the Schwinger a = α/(2π)) as the BAM closure quantum but flagged the full covariant path-integral measure as open; BAM's only other quantization was saddle-point/instanton (the PRs #87–#90 bounces). This sprint constructs the full measure Z = Σ_sectors ∫ Dμ[X] e^{−S_BAM[X]} around that factor.

THE ARENA. A throat is its closure loop X: S¹ → base of the Hopf fibration (winding k, Hopf charge c₁), so the configuration space is loop space LS³ (antipodal for the non-orientable sector), and Z = Σ_{k odd, c₁∈ℤ, n_part} ∫_{LS³/(Diff S¹ ⋉ U(1)_Hopf ⋉ Z₂)} Dμ[X] e^{−S_BAM[X]}: the discrete sum is the closure ledger, the continuous integral the gauge-fixed loop integral, and PR #74's 1/(2π) is the per-loop-dimension normalization (dk/(2π)).

WHAT IS STRUCTURALLY FIXED (and computable). (1) The closure quantum is the loop holonomy: e^{i∮A} = e^{ikχ} single-valued ⟹ period 2π. (2) The measure's superselection sectors are the homotopy classes — winding k (π₁ of the Hopf fibre), Hopf charge c₁ ∈ π₃(S²) = ℤ, partition n_part — so Σ_sectors IS the closure ledger. (3) The odd-k lemma is UPGRADED to the measure's Z₂ orientation-anomaly condition: the Möbius loop space is non-orientable, the integrand a twisted-bundle section anti-periodic under the antipodal χ → χ+π, so consistency forces e^{ikπ} = −1 ⟹ k odd (even k closes on the torus cover only). (4) Stationary phase reproduces e^{−S_BAM[saddle]} — the bounce actions of PRs #87–#90 — as the leading term.

THE HARD PART. Reparametrization invariance Diff(S¹) must be gauge-fixed ⟹ a Faddeev–Popov (bc-ghost) determinant, so the measure factor is (FP-det) × (fluctuation-det). The fluctuation operator is the second variation of S_BAM — the Tangherlini cavity operator — and its low spectrum is positive (min ω² ≈ 1.11 > 0), so the saddle is stable and the Gaussian measure is real and well-defined pointwise (the tunneling sector's single Coleman negative mode supplies Im Z = the decay rate; the zero modes are collective coordinates, traded for the moduli measure — both standard).

WHAT REMAINS OPEN. The bare fluctuation determinant Π_n ω_n diverges — the log-det partial sums grow without bound (≈ 15, 146, 358, 850, 1964 at N = 10, 50, 100, 200, 400 modes) — so it needs zeta/heat-kernel regularization, and the normalisation Z is regularization-dependent and not yet rigorously constructed. So the loop-measure formalism gives a consistent STRUCTURAL definition of the S_BAM measure, but the HARD analytic core — a finite, regularization-independent fluctuation determinant / a constructive measure — stays open. Crucially the prior saddle-point results are the leading e^{−S} term and do not depend on the unresolved normalisation, so they stand.

## What this establishes (and does not)

- **Established (structural):** the S_BAM measure as a closure-ledger sector sum over a Diff(S¹)-gauge-fixed loop-space integral; PR #74's 1/(2π) as the per-dimension factor; the closure quantum 2π as the loop holonomy; the odd-k lemma upgraded to the Z₂ orientation-anomaly condition; the PRs #87–#90 bounces as the leading saddle; the fluctuation operator stable.
- **Open (analytic):** a finite, regularization-independent fluctuation determinant / a constructive measure — the bare determinant diverges, so the normalisation `Z` is not yet rigorously defined. Prior saddle-point results stand (leading `e^{−S}`, normalisation-independent).
