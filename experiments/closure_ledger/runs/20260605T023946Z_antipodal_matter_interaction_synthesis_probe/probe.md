# Antipodal matter interaction ledger synthesis (PR #139)

**Run:** 2026-06-05T02:39:46+00:00

Capstone of the antipodal matter-interaction arc — PRs #129–#138. Re-verifies a keystone from each arc member together, shows the whole arc is two threads from one postulate (the antipodal Z₂ and the unitary mirror), and lays out the epistemic ledger.

- **Thread A**: antipodal Z₂ (−1)^l: BC (#129) + kernel grading (#135) + vertex selection (#137/#138)
- **Thread B**: unitary mirror: spectrum (#130) + propagator (#135) + self-energy (#136) + vacuum (#138/#122)
- **One postulate**: the real l-parity BC (#129) is both threads; from the antipodal identification (#128)
- **Derived**: Z₂ selection structure; unitary stable propagator/self-energy/vacuum
- **Postulated**: the antipodal identification (#128) — self-consistent, not forced
- **Input**: coupling magnitudes λ_3, λ_4 (sign λ_4 > 0 from #122)
- **Open**: S_BAM vertex generation; higher loops/vertices; scale (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | synthesise the antipodal matter-interaction arc (#129–#138) | **PASS** |
| T2 | `T2_arc_layer_by_layer` | the arc layer by layer (#129→#130→#135→#136→#137/#138) | **PASS** |
| T3 | `T3_keystones_reverified_together` | keystones re-verified together (cross-arc consistency) | **PASS** |
| T4 | `T4_thread_A_antipodal_z2` | Thread A — antipodal Z₂ (−1)^l: BC, kernel, vertices | **PASS** |
| T5 | `T5_thread_B_unitarity_stability` | Thread B — unitarity/stability: spectrum→propagator→Σ→vacuum | **PASS** |
| T6 | `T6_two_threads_one_postulate` | two threads, one postulate (the real l-parity BC is both) | **PASS** |
| T7 | `T7_epistemic_ledger` | epistemic ledger: derived / postulated / input / open | **PASS** |
| T8 | `T8_assessment` | ANTIPODAL_MATTER_INTERACTION_SYNTHESIS_Z2_SELECTION_AND_UNITARY_STABILITY | **PASS** |

## The arc keystones, re-verified together

| PR | keystone | value |
|---|---|---|
| #129 | Y_l antipodal parity | [1, -1, 1, -1] = (−1)^l |
| #130 | antipodal fundamental (real) | 1.166 |
| #130 | absorbing fundamental (complex) | 1.893-1.141i |
| #135 | kernel reciprocity \|G−Gᵀ\|/\|G\| | 1e-13 |
| #135 | lowest pole ω² | 1.36 |
| #136 | lightest-mode Im Σ (stable) | -0.003 |
| #138 | ∫ψ⁴ overlap (bounded vacuum) | 1.033 |

All mutually consistent in one run.

## Two threads, one postulate

| | Thread A (Z₂ selection) | Thread B (unitarity/stability) |
|---|---|---|
| #129 BC | even-l N / odd-l D (parity) | unitary mirror (zero flux) |
| #130 | — | real stable spectrum |
| #135 | kernel parity grading | reciprocal unitary propagator |
| #136 | — | stable lightest mode |
| #137/#138 | Σl-even vertex selection | bounded vacuum (#122) |

Both threads flow from the single antipodal identification (#128): the real l-parity boundary condition (#129) is at once the Z₂ grading and the unitary mirror.

## Verdict

**ANTIPODAL_MATTER_INTERACTION_SYNTHESIS_Z2_SELECTION_AND_UNITARY_STABILITY.** THE ANTIPODAL MATTER-INTERACTION ARC IS TWO THREADS FROM ONE POSTULATE: THE ANTIPODAL Z₂ SELECTION AND THE UNITARY MIRROR. PRs #129–#138 built the BAM matter interaction layer by layer on the antipodal horizon boundary data; this capstone re-verifies the arc's keystones together and organises them.

THE ARC, LAYER BY LAYER. The boundary condition (#129, antipodal l-parity, a unitary mirror) → the spectrum (#130, real and undamped vs absorbing ringdown) → the free propagator (#135, the reciprocal unitary resolvent, a mode-sum over the stable modes) → the one-loop self-energy (#136, a finite real mass shift with an exactly-stable lightest mode) → the cubic and quartic vertices (#137/#138, the Z₂ selection rule and the bounded-below vacuum).

THE KEYSTONES RE-VERIFIED TOGETHER. In one run: Y_l(−x) = (−1)^l Y_l (#129); the exchange kernel reciprocal with real poles (#135); the lightest-mode Im Σ = 0 (#136); ∫ψ⁴ > 0 (#138); the antipodal fundamental real (≈ 1.17) vs the absorbing one complex (≈ 1.89 − 1.16i, #130). All mutually consistent.

THREAD A — THE ANTIPODAL Z₂ (−1)^l. The C-swap inversion x → −x (#63) carries Y_l → (−1)^l Y_l. This one parity fixes the boundary condition (#129), grades the exchange kernel (#135), and selects the cubic (#137) and quartic (#138) vertices (Σl even). One Z₂ threading the entire interaction structure.

THREAD B — UNITARITY / STABILITY. The antipodal BC is a unitary mirror (#129), giving a real, stable spectrum (#130), a unitary reciprocal propagator (#135), a unitarity-preserving self-energy with an exactly-stable lightest mode (#136), and a bounded-below interacting vacuum (#138, the #122 measure-convergence condition). BAM matter is stable at every order because the throat reflects antipodally.

TWO THREADS, ONE POSTULATE. The threads are not independent: the real l-parity boundary condition IS both the Z₂ grading (Thread A) and the unitary mirror (Thread B). One antipodal identification (#128), two faces — the selection structure and the stability.

THE EPISTEMIC LEDGER. DERIVED (given the antipodal BC): the Z₂ selection structure of the propagator and vertices; the unitary, reciprocal, real-pole propagator; the unitarity-preserving self-energy and stable lightest mode; the bounded-below stable vacuum. POSTULATED: the antipodal identification itself (#128) — shown self-consistent (a unitary, stable, bounded interacting theory), not forced. INPUT: the coupling magnitudes λ_3, λ_4 (the sign λ_4 > 0 required by #122; the magnitudes not derived). OPEN: the S_BAM generation of the vertices; higher loops / vertices; the bulk-scale (#133) and flavor (#134) residuals.

SCOPE. A synthesis/consistency capstone: it re-verifies the arc's keystones together and organises them into the two threads and the ledger. It does not add new derivations, remove any open item, or strengthen the antipodal postulate from "self-consistent" to "forced".
