# Quartic vertex ledger and bounded interaction audit (PR #138)

**Run:** 2026-06-05T02:21:34+00:00

Extends the #137 cubic-vertex ledger to the quartic vertex and audits whether the BAM matter interaction is bounded below (a stable vacuum). The quartic carries the same antipodal Z₂ selection rule; its geometric self-overlap ∫ψ⁴ > 0 is positive, so the effective potential is bounded below — the S_BAM measure-convergence condition (#122).

- **Factorisation**: V_4 = λ_4 · (angular ∫YYYY) · (radial ∫ψψψψ)
- **Angular rule**: DERIVED: Σl even (same antipodal Z₂ as #137) + SO(4) common channel
- **Positive overlap**: g_4 = ∫ψ_k⁴ dr* > 0 (manifest)
- **Bounded**: a⁴ coeff λ_4 g_4/24 > 0 ⟹ V → +∞ ⟹ bounded below ⟹ stable vacuum
- **Measure tie**: boundedness = convergence of ∫ Dμ e^{−S} (#122)
- **Input**: coupling magnitudes λ_3, λ_4 (sign λ_4 > 0 from #122)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | quartic vertex ledger + bounded interaction audit | **PASS** |
| T2 | `T2_quartic_factorisation` | factorisation V_4 = λ_4 · (angular ∫YYYY) · (radial ∫ψψψψ) | **PASS** |
| T3 | `T3_quartic_angular_selection_rule` | quartic selection rule: Σl even (same Z₂ as #137) + SO(4) channel | **PASS** |
| T4 | `T4_positive_quartic_self_overlap` | positive quartic self-overlap g_4 = ∫ψ_k⁴ > 0 (manifest) | **PASS** |
| T5 | `T5_bounded_interaction_stable_vacuum` | bounded interaction: a⁴ coeff > 0 ⟹ V bounded ⟹ stable vacuum | **PASS** |
| T6 | `T6_boundedness_is_measure_convergence` | boundedness = S_BAM measure convergence (#122); stability thread | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger/scope: rule + overlap + boundedness derived; couplings input | **PASS** |
| T8 | `T8_assessment` | QUARTIC_VERTEX_LEDGER_Z2_SELECTION_POSITIVE_OVERLAP_BOUNDED_STABLE_VACUUM | **PASS** |

## The quartic angular selection rule (Σl even, the same antipodal Z₂)

| (l₁,l₂,l₃,l₄) | Σl even? | SO(4) channel? | ∫YYYY | nonzero? |
|---|:---:|:---:|---:|:---:|
| (0,0,0,0) | ✓ | ✓ | 1.0 | ✓ |
| (1,1,0,0) | ✓ | ✓ | 0.25 | ✓ |
| (1,1,1,1) | ✓ | ✓ | 0.041667 | ✓ |
| (2,1,1,0) | ✓ | ✓ | 0.041667 | ✓ |
| (1,1,1,0) | ✗ | ✓ | 0.0 | ✗ |
| (1,0,0,0) | ✗ | ✗ | 0.0 | ✗ |

## Positive self-overlap ⟹ bounded interaction (stable vacuum)

| mode k | g_4 = ∫ψ_k⁴ dr* |
|---:|---:|
| 0 | 1.0335 |
| 1 | 1.0205 |
| 2 | 1.0204 |
| 3 | 1.0204 |

| λ₃ (cubic) | a⁴ coeff | V(+10⁴) | V(−10⁴) | bounded? |
|---:|---:|---:|---:|:---:|
| 0.0 | 0.0431 | 431000000000000.0 | 431000000000000.0 | ✓ |
| 5.0 | 0.0431 | 430000000000000.0 | 431000000000000.0 | ✓ |
| 50.0 | 0.0431 | 422000000000000.0 | 439000000000000.0 | ✓ |
| 200.0 | 0.0431 | 397000000000000.0 | 464000000000000.0 | ✓ |

`g_4 = ∫ψ⁴ > 0` ⟹ the a⁴ coefficient is positive ⟹ `V → +∞` for any cubic ⟹ bounded below, a **stable vacuum**. Boundedness is the condition for the S_BAM measure `∫ Dμ e^{−S}` to converge (#122).

## Verdict

**QUARTIC_VERTEX_LEDGER_Z2_SELECTION_POSITIVE_OVERLAP_BOUNDED_STABLE_VACUUM.** THE QUARTIC VERTEX CARRIES THE SAME ANTIPODAL Z₂ SELECTION RULE, AND ITS POSITIVE GEOMETRIC OVERLAP BOUNDS THE INTERACTION BELOW — A STABLE VACUUM. PR #137 drew the cubic-vertex ledger; this probe extends it to the quartic and audits the stability the cubic alone left open.

THE QUARTIC VERTEX FACTORISES. V_4 = λ_4 · [∫_{S³} Y_{l1}Y_{l2}Y_{l3}Y_{l4} dΩ] · [∫ ψ_k ψ_l ψ_m ψ_n dr*] — angular integral × radial overlap × coupling, the same factorisation as the cubic.

THE QUARTIC SELECTION RULE IS DERIVED. The angular integral is nonzero only if (a) l1+l2+l3+l4 is EVEN — the SAME antipodal parity Z₂ as the cubic (#137): under x → −x (the throat ↔ antithroat C-swap #63), Y_l → (−1)^l Y_l, so (−1)^{Σl} = +1 over the inversion-symmetric S³ — and (b) a common SO(4) intermediate channel exists, L ∈ [|l1−l2|,l1+l2] ∩ [|l3−l4|,l3+l4]. The antipodal Z₂ persists from the cubic to the quartic (verified exactly via the S³ monomial integral: odd-Σl → 0).

THE QUARTIC SELF-OVERLAP IS POSITIVE. The diagonal quartic radial overlap g_4 = ∫ ψ_k⁴ dr* > 0 is manifestly positive — an integral of a fourth power — for every mode.

THE INTERACTION IS BOUNDED BELOW — A STABLE VACUUM. The single-mode effective potential V(a) = ½ ω_k² a² + (λ_3 g_3/6) a³ + (λ_4 g_4/24) a⁴ is bounded below iff its a⁴ coefficient λ_4 g_4/24 > 0. With g_4 > 0 and λ_4 > 0, V → +∞ as |a| → ∞ for ANY cubic — the vacuum is stable, the cubic only tilting the minimum, never unbounding it. The positive geometric quartic bounds the cubic.

BOUNDEDNESS IS THE MEASURE-CONVERGENCE CONDITION (#122). A bounded-below action is exactly the condition for the path-integral measure ∫ Dμ e^{−S} to converge — established non-perturbatively for the Z₂-graded sector sum in #122. So the positive quartic is not an extra assumption: a bounded interaction is required by, and consistent with, the measure's existence. This extends the program's stability thread — free modes stable (#130), one-loop self-energy unitarity-preserving (#136), and now the full interacting vacuum bounded below.

SCOPE. Audits the quartic vertex STRUCTURE (the Z₂ + SO(4) selection rule, the positive geometric overlap, the boundedness it implies) and ties boundedness to #122. It does NOT derive the coupling magnitudes λ_3, λ_4 from S_BAM (the sign λ_4 > 0 is required by #122 convergence; the magnitude is input), nor address quintic / higher vertices. The bulk-scale (#133) and flavor (#134) residuals stand.
