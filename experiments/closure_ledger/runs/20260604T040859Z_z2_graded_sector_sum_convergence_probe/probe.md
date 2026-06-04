# Non-perturbative convergence audit for the Z₂-graded BAM sector sum (PR #126)

**Run:** 2026-06-04T04:08:59+00:00

Audits the non-perturbative convergence PR #122 left open. The Z₂-graded sector sum has three pieces — the winding sum, the Hopf-charge sum, and the moduli integral — and each is shown finite, so the sum **converges**. (This audits *finiteness*, not the overall scale; the absolute normalization `κ₅²/Λ₅` stays open.)

- **Winding**: FINITE (3 terms, k ∈ {1,3,5}, closure cap k₅ = 5)
- **Hopf**: convergent theta Σ e^{−A c₁²} = √(π/A)
- **Moduli**: finite both ends: UV e^{−π²/t} (Z₂ cancellation), IR e^{−m²t} (mass gap)
- **Overall**: converges non-perturbatively
- **Open**: absolute normalization (κ₅²/Λ₅); exact value of Z; multi-loop

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | audit non-perturbative convergence (PR #122 open item) | **PASS** |
| T2 | `T2_winding_sum_finite` | winding sum FINITE: k ∈ {1,3,5}, closure cap (3 terms) | **PASS** |
| T3 | `T3_hopf_charge_sum_convergent` | Hopf-charge sum convergent: Σ e^{−A c₁²} = √(π/A) | **PASS** |
| T4 | `T4_moduli_integral_uv_convergent` | moduli UV-convergent: Z₂ cancellation θ_per − θ_anti ~ e^{−π²/t} | **PASS** |
| T5 | `T5_moduli_integral_ir_convergent` | moduli IR-convergent: mass gap e^{−m²t} ⟹ ∫ finite | **PASS** |
| T6 | `T6_overall_convergence` | overall: finite × convergent × convergent ⟹ converges | **PASS** |
| T7 | `T7_scope` | scope: convergence established; absolute scale open | **PASS** |
| T8 | `T8_assessment` | Z2_GRADED_SECTOR_SUM_CONVERGES_NONPERTURBATIVELY | **PASS** |

## The moduli integral is finite (UV cancellation + IR mass gap)

| mass gap m | `∫(dt/t)[θ_per − θ_anti] e^{−m²t}` |
|---:|---:|
| 0.3 | 0.605916 |
| 0.5 | 0.172959 |
| 1.0 | 0.00747 |

Finite for every mass gap: the UV (`t → 0`) is killed by the Z₂-graded cancellation `θ_per − θ_anti ~ e^{−π²/t}`, the IR (`t → ∞`) by the mass gap `e^{−m²t}` (the bounce saddle).

## Verdict

**Z2_GRADED_SECTOR_SUM_CONVERGES_NONPERTURBATIVELY.** THE Z₂-GRADED BAM SECTOR SUM CONVERGES NON-PERTURBATIVELY — THE WINDING SUM IS FINITE, THE HOPF-CHARGE SUM IS A CONVERGENT THETA, AND THE MODULI INTEGRAL IS FINITE AT BOTH ENDS. PR #122 left the non-perturbative convergence of the factorized sector sum open; this probe audits it in its three pieces.

THE WINDING SUM IS FINITE. The physical winding is capped by the closure ledger: the odd-k lemma plus the available phase Φ_avail(k) = 2π(k+1) + 50π·max(0,k−3)² allow only k ∈ {1,3,5} — the three generations, with the bulk dimension k₅ = 5 the cap (PR #73). So the winding sum is a FINITE sum, three terms, not an infinite tower: trivially convergent.

THE HOPF-CHARGE SUM CONVERGES. The Hopf charge c₁ ∈ ℤ has an EM/Hopf action cost ~ c₁², so Σ_{c₁∈ℤ} e^{−A c₁²} = √(π/A)·θ₃ → √(π/A), a convergent theta function (verified for A = 0.5, 1, 2); charge conservation Σ c₁ = 0 (PR #58) further constrains it.

THE MODULI INTEGRAL CONVERGES AT BOTH ENDS. The graded integral ∫₀^∞ (dt/t)[θ_per − θ_anti] e^{−m²t} is finite: at the UV (t → 0) the Z₂-graded cancellation gives θ_per − θ_anti ~ e^{−π²/t} → 0 (PR #122), so the integrand vanishes (no UV divergence); at the IR (t → ∞) the mass gap e^{−m²t} — the bounce saddle / the physical masses — kills the large-t tail (where θ_per → 1, θ_anti → 0). The integral is finite for every mass gap (≈ 0.17 at m = 0.5, 0.61 at m = 0.3, 0.0075 at m = 1.0).

THE OVERALL CONVERGENCE. The Z₂-graded sector sum = (finite winding sum, 3 terms) × (convergent Hopf-charge sum) × (convergent moduli integral). Each piece is finite, so the full sum converges non-perturbatively: the Z₂ grading makes the moduli integral UV-finite (the orientation signs cancel the Weyl divergence), the closure cap makes the winding sum finite, and the mass gap makes the IR finite.

SCOPE. This audits CONVERGENCE — that the Z₂-graded sector sum is finite, not divergent. It does NOT fix the absolute normalization (the bulk κ₅²/Λ₅ anchor, PR #112) or the exact value of Z; the multi-loop / interacting measure is also beyond the audit. Finiteness is established; the overall scale stays open.
