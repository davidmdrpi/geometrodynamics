# The index mechanism: a Pin/Dirac zero mode for the k=1 sector (PR #195)

**Run:** 2026-07-02T04:04:00+00:00

Answers #194's mechanism question: the Pin⁻ spinor (which BAM already requires) on the #193 monopole reduction carries an Atiyah–Singer index that pins exactly k chiral zero modes in winding sector k. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: answer #194 — mechanism or hierarchy problem? | **PASS** |
| T2 | `T2_candidate` | the candidate: Pin⁻ ⟹ spinor ⟹ Dirac + index | **PASS** |
| T3 | `T3_count_verified` | {1,3,5} zero modes in sectors {1,3,5}; exact towers | **PASS** |
| T4 | `T4_protection_certificates` | certificates: zero PINNED (1e-10/1e-15) vs scalar moves | **PASS** |
| T5 | `T5_flux_change_control` | count moves only with a flux quantum (1 → 3 → 5) | **PASS** |
| T6 | `T6_natural_mass` | natural mass: one-mouth forbidden; two-mouth linear ε·o | **PASS** |
| T7 | `T7_honest_scope` | scope: mechanism established, ladder not re-derived | **PASS** |
| T8 | `T8_assessment` | MECHANISM FOUND — the spinor supplies the protection | **PASS** |

## The index, verified (D²₊ = L_{q−½} − (q−½), D²₋ = L_{q+½} + (q+½))

| k | zero modes | predicted | worst residual | 1st excited (+) | lowest (−) | exact gap |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1.5e-10 | 2.0 | 2.0 | 2.0 |
| 3 | 3 | 3 | 2.88e-08 | 4.0 | 4.0 | 4.0 |
| 5 | 5 | 5 | 9.34e-08 | 5.999999 | 6.0 | 6.0 |

## The protection certificates

- gauge wobble (flux fixed, ε = 0.05): zero energy pinned in **[0, 6e-10]** (wavefunction deforms, energy does not)
- metric deformation: ker D conformally rigid — certificate **9e-16**, Ω-independent
- the scalar on the same deformed metric: 1.5 → 1.401226 (moves **-0.0988**) — **8 orders** of contrast

## The natural mass (two-mouth pairing)

- within one mouth: first-order lift **forbidden** (no opposite-chirality j = q−½ partner; lowest D²₋ = 2.0)
- two mouths (±k winding): |E_e| = ε·o, o = **1.0** — linear, multiplicative, sign-stable ('t Hooft-natural)

## Verdict

**K1_ZERO_MODE_IS_INDEX_PROTECTED_THE_PIN_DIRAC_STRUCTURE_SUPPLIES_THE_MECHANISM_THE_SCALAR_SURROGATE_LACKS.** MECHANISM FOUND — #194's question is answered YES, with no new ingredients.

THE MECHANISM. BAM throats are Pin⁻ (#183/#188), so the throat mode is a SPINOR; on the #193 sector reduction it is a spinor on the base S² with monopole charge q = k/2, and the Atiyah–Singer index pins exactly 2q = k chiral zero modes in winding sector k — verified: sectors {1,3,5} carry {1,3,5} zero modes (residuals ≤ 1e-7), opposite chirality gapped at 2q+1, towers matching the exact Dirac spectrum. The 'geometric identity tying the diagonal to the repulsion' is the SUSY factorization D² = A†A — a perfect square that cancels on the kernel, available only to the spinor.

THE PROTECTION (the discriminator vs #194). Under a flux-preserving gauge wobble the zero energy is certified pinned in [0, 6e-10]; under a metric deformation of the base, ker D is conformally rigid (certificate 9e-16, Ω-independent) while the SCALAR ground on the same metric moves by -0.099 — 8 orders of contrast: energy-PINNED vs energy-TUNED. The count changes only with a flux quantum (1 → 3 → 5).

THE NATURAL MASS. Within one mouth a first-order lift is forbidden by angular momentum (no opposite-chirality j = q−½ partner exists); the mass comes from pairing the throat's TWO MOUTHS (±k winding, opposite chirality — the BAM wormhole supplies the Dirac pair): |E_e| = ε·o with o = 1.000 — linear, multiplicative, sign-stable, 't Hooft-natural. The surrogate's dialing (#194: a sign-flipping difference of two O(7) numbers, Δ = 74.7) is an artifact of treating a spinor problem with scalar dynamics. FOLLOW-UP: rebuild the mass ladder on the Dirac tower; compute the mouth coupling ε from the throat overlap machinery (#185/#190).
