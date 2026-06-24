# The odd-k ladder: forced, rigid, unique to the non-orientable 5D geometry (PR #174)

**Run:** 2026-06-24T05:18:12+00:00

Takes the cleanest discrete feature — the odd-k charged-lepton ladder {1,3,5} — and asks whether the geometry forces it and whether anything but this geometry could produce it, using the #173 direction analysis (active / null / mixed). *(QFT on the classical throat, not quantum gravity.)*

- **Origin**: T=iσ_y, T²=−I: odd k off-diagonal (fermion/non-orientable), even k diagonal (boson)
- **Active**: rank-10 subspace moves observables linearly (exponent ≈ 1.0)
- **Null**: 10-dim kernel flat to first order (exponent ≈ 2.0; ~10⁴× smaller)
- **Mixed**: active-dominated (≈ 1.0): nonlinearity does not break the rank story
- **Forced**: labels are integer winding + ℤ₂ grading; no continuous knob moves them
- **Unique**: orientable geometry → even/bosonic; {1,3,5} needs k≤5=D_bulk

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the two questions for the odd-k ladder | **PASS** |
| T2 | `T2_discrete_feature_and_grading` | the discrete feature + the T^k orientability grading | **PASS** |
| T3 | `T3_active_singular_directions` | active directions (10): move observables linearly | **PASS** |
| T4 | `T4_null_compensator_directions` | null directions (10): flat to first order (quadratic) | **PASS** |
| T5 | `T5_mixed_directions_nonlinearity` | mixed directions: nonlinearity does not break the rank story | **PASS** |
| T6 | `T6_forced_rigid_discreteness` | FORCED: the ladder is rigid against all continuous deformation | **PASS** |
| T7 | `T7_unique_to_non_orientable_5d` | UNIQUE: only the non-orientable 5D geometry gives odd-{1,3,5} | **PASS** |
| T8 | `T8_assessment` | ODD_K_LADDER_FORCED_RIGID_UNIQUE | **PASS** |

## The monodromy grading  (T = iσ_y, T² = −I)

| k | T^k class | sector |
|---|---|---|
| 1 | off_diagonal_fermion | odd → fermion (non-orientable) |
| 2 | diagonal_boson | even → boson (orientable) |
| 3 | off_diagonal_fermion | odd → fermion (non-orientable) |
| 4 | diagonal_boson | even → boson (orientable) |
| 5 | off_diagonal_fermion | odd → fermion (non-orientable) |
| 6 | diagonal_boson | even → boson (orientable) |

(charged leptons = spin-½ fermions ⟹ odd k; k ≤ k_5 = 5 ⟹ {1,3,5})

## The three direction sets

| set | dimension | scaling exponent | reading |
|---|---:|---:|---|
| active | 10 | 1.03 | linear — moves observables |
| null | 10 | 2.0 | quadratic — flat to 1st order |
| mixed | — | 1.02 | active-dominated — rank story holds |

## Verdict

**ODD_K_LADDER_FORCED_RIGID_UNIQUE_TO_NON_ORIENTABLE_5D_GEOMETRY.** FORCED, RIGID, UNIQUE. The cleanest discrete feature — the odd-k charged-lepton ladder {1,3,5} — is geometry, not bookkeeping.

THE ORIGIN. The throat monodromy T = iσ_y (T² = −I) makes T^k off-diagonal for odd k (the orientation-reversing, non-orientable closure of a spin-½ fermion) and diagonal for even k (orientable, bosonic). Charged leptons are fermions ⟹ odd k; the bulk boundary k ≤ k_5 = 5 = D_bulk ⟹ {1,3,5} = three generations.

THE THREE DIRECTION SETS. ACTIVE (the rank-10 subspace): perturbing along it moves the observables LINEARLY (exponent 1.03) — the 10 directions that actually move the masses and the CKM. NULL (the 10-dim kernel): flat to first order (exponent 2.00, quadratic; ~10⁴× smaller response than active). MIXED (active+null): active-dominated (exponent 1.02) — nonlinear effects do NOT break the local rank story; the null leakage stays quadratic.

FORCED & RIGID. The odd-k labels and the generation count are integer winding plus the ℤ₂ orientability grading (T² = −I) — discrete topological data outside the entire continuous deformation manifold (there is no generation-number knob). The continuous geometry deforms only the masses and the CKM; the ladder is rigid against every active, null, and mixed deformation, linear and nonlinear.

UNIQUE. An orientable geometry (T² = +I) gives the even/bosonic sector, not an odd-only fermion ladder; the specific {1,3,5} needs k ≤ 5 = the bulk dimension. So odd-{1,3,5} is the joint signature of the non-orientable antipodal spin structure and the 5D bulk — an exclusion/signature argument within BAM, not a no-go against every conceivable alternative.
