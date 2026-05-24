# Charge conjugation from the inner/outer swap

**Run:** 2026-05-24T04:00:23+00:00

Promotes C-symmetry from a postulate to a geometric statement: the inner/outer reflection across the throat (S: r ↦ 2R_MID − r) is an involution under which the eigenmodes are odd (B3) and the integrated Hopf curvature flips sign (c₁ → −c₁), taking a throat to its antithroat.

- **C operation**: `S : r ↦ 2 R_MID − r (inner/outer reflection)`
- **Effect**: c₁ → −c₁ (throat → antithroat)
- **Charge**: `c₁ = (1/2π)∮F = ±1 (integrated Hopf curvature)`
- **Consistency**: antipodal Z₂ / T=iσ_y (B2); pair production (#58)
- **B4 caveat**: c₁=±1 dimensionless topological integer; C scale-independent

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_charge_is_integrated_hopf_curvature` | c₁ = ±1 (|c₁|=1.0000) | **PASS** |
| T2 | `T2_swap_involution_modes_odd` | S²=id; fixes throat; modes odd (B3) | **PASS** |
| T3 | `T3_swap_reverses_mouth_orientation` | n̂=±r̂ opposite (inner/outer); swap reverses dr | **PASS** |
| T4 | `T4_integrated_curvature_flips` | c₁→−c₁ (-1.00 / +1.00) | **PASS** |
| T5 | `T5_C_is_inner_outer_swap` | C=swap: throat→antithroat (#58, B2) | **PASS** |
| T6 | `T6_C_squared_identity` | C²=id (involution) | **PASS** |
| T7 | `T7_falsification_b4` | swap flips c₁ (not invariant); topological integer | **PASS** |
| T8 | `T8_assessment` | C geometric, not a postulate | **PASS** |

## T1: Charge = integrated Hopf curvature

- c₁ = (1/2π)∮F, |c₁| = 1.000000 (analytic -1); unit charge: True

## T2: The inner/outer swap is an involution; modes odd (B3)

- fixes throat R_MID: True
- exchanges R_INNER ↔ R_OUTER: True
- involution S²=id: True
- modes odd under S (B3 antisymmetric extension): True
- Dirichlet node at throat: True

## T3: The swap reverses the mouth orientation

- outward normal n̂·r̂: outer = +1, inner = -1 (opposite: True)
- swap reverses dr (d(2R_MID−r)=−dr): True

## T4: The integrated curvature flips

- c₁ outer mouth (dχ∧dφ) = -1.000000
- c₁ inner mouth (reversed) = +1.000000
- sum = +0.00e+00 → c₁ → −c₁ under the swap: True

## T5: C = swap (the geometric statement)

- throat charge = +1; C image (antithroat) = -1; C flips charge: True
- matches pair-production antithroat (#58): True
- matches antipodal Z₂ (B2): True

## T6: C² = id; discrete-symmetry consistency

- c₁: +1 →(C) -1 →(C) +1; C²=id: True
- coordinate involution S²=id: True

## T7: Falsification / B4

- charge flips under swap: True (a swap-invariant charge would falsify)
- c₁ is a topological integer (scale-independent): True
- **BAM passes: True**

## T8: Assessment

- C operation: S : r ↦ 2 R_MID − r (inner/outer reflection)
- effect: c₁ → −c₁ (throat → antithroat)
- involution: C² = id
- consistency: antipodal Z₂ / T=iσ_y (B2); pair production (#58)
- remaining: full CPT from S_BAM; C on the Dirac spinor (ψ → C ψ̄ᵀ)

## Verdict

**C_IS_INNER_OUTER_SWAP.** C IS THE INNER/OUTER SWAP. Charge conjugation is promoted from a postulate to a geometric statement: the inner/outer reflection across the throat is C.

THE SWAP. The two wormhole regions r < R_MID (inner) and r > R_MID (outer), with R_INNER, R_OUTER symmetric about the throat, are exchanged by the reflection S: r ↦ 2R_MID − r — fixing R_MID, exchanging R_INNER ↔ R_OUTER, an involution (S²=id). It is the reflection of the B3 hard-wall odd extension: the throat modes are odd under S (u(2R_MID−r)=−u(r)).

THE CHARGE FLIPS. The charge is the integrated Hopf curvature c₁ = (1/2π)∮F = ±1 (first Chern number). The mouth's induced orientation is set by its outward normal n̂=±r̂ (outer +r̂, inner −r̂); the normals point oppositely, so the inner and outer mouths carry opposite orientation and c₁(inner) = −c₁(outer) — exactly compute_c1's c1_chiphi=−1, c1_phichi=+1. The swap reverses the mouth orientation, flipping the integrated curvature: c₁ → −c₁. With the modes odd under S, the swap takes a throat (c₁=+1) to its antithroat (c₁=−1).

C IS GEOMETRIC. Charge conjugation — particle → antiparticle — is realized as C = S, with C: c₁ → −c₁ and C²=id. It is no longer a postulate but the throat-reflection involution, consistent with the antipodal Z₂ / T=iσ_y (B2) and the C-conjugate antithroat of pair production (#58). B4: c₁=±1 is a dimensionless topological integer, C a discrete geometric involution — scale-independent. Remaining: the full CPT statement from S_BAM, and the explicit action of S on the throat Dirac spinor (ψ → C ψ̄ᵀ) beyond the charge sign.

## What this leaves open

- **The full CPT theorem from S_BAM.** C (this probe), P (parity), and T (T=iσ_y, B2) as one geometric CPT statement.
- **C on the Dirac spinor.** The explicit action of S on the throat spinor (ψ → C ψ̄ᵀ), beyond the charge sign.
