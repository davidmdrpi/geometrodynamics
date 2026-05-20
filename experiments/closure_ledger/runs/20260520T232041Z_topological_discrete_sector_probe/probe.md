# Topological/discrete action-sector probe

**Run:** 2026-05-20T23:20:41+00:00

Tests whether the scaffold barriers B1 (closure quantum) and B2 (antipodal Z₂) can be promoted from imposed constraints to terms in a topological action.

## The topological/discrete sector

```
RP³ = S³/Z₂  +  non-trivial spin structure (T² = −I)  +  2π winding θ-term
```

- **B1**: closure quantum → winding θ-term S_top = 2π·n
- **B2**: antipodal Z₂ → RP³ quotient + non-trivial spin structure
- **Unification**: great circle on S³ (2π = closure quantum, B1) = double traversal of non-contractible RP³ loop (π, π₁ generator, B2)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_antipodal_involution_Z2` | σ² = id (residual 0.0e+00); free action → smooth RP³ | **PASS** |
| T2 | `T2_double_cover_closure_quantum` | great circle 2π = 2× RP³ loop π (diff 0.0e+00) | **PASS** |
| T3 | `T3_pi1_rp3_is_Z2` | π₁(RP³) = Z₂: True | **PASS** |
| T4 | `T4_spin_structure_T2_minus_I` | T² = −I (non-trivial spin structure), residual 0.0e+00 | **PASS** |
| T5 | `T5_winding_theta_term` | ∮dφ = 2πn integer: True; metric-independent: True | **PASS** |
| T6 | `T6_variational_consistency` | δS_top under local variation = 1.2e-13 | **PASS** |
| T7 | `T7_K_Q_from_topological_sector` | F² = K²·Q from topological sector (diff 2.8e-14) | **PASS** |
| T8 | `T8_promotion_assessment` | both promoted: True; barriers 5 → 3 | **PASS** |

## T2: Double cover S³ → RP³

| RP³ loop length (π) | great circle (2π) | closure quantum |
|---:|---:|---:|
| 3.141593 | 6.283185 | 6.283185 |
| 3.141593 | 6.283185 | 6.283185 |
| 3.141593 | 6.283185 | 6.283185 |
| 3.141593 | 6.283185 | 6.283185 |

## T4: Spin structure T² = −I

```
  +0 +1
  -1 +0
```
T² = −I residual: `0.00e+00`; det = 1 residual: `0.00e+00`. RP³ has 2 spin structures (H¹(RP³,Z₂)=Z₂); T²=−I selects the non-trivial (antiperiodic) one.

## T5: Winding θ-term ∮dφ = 2πn

| intended n | measured winding | S_top = 2πn | matches 2πn |
|---:|---:|---:|:---:|
| 1 | 1.000000 | 6.283185 | True |
| 2 | 2.000000 | 12.566371 | True |
| 3 | 3.000000 | 18.849556 | True |
| -1 | -1.000000 | -6.283185 | True |
| 5 | 5.000000 | 31.415927 | True |

All windings integer: **True**; winding metric-independent (reparam-invariant): **True**.

## T7: K and Q from the topological sector

| x | cosθ | K (winding + Z₂) | Q | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | -0.70 | 0.1818 | 0.0644 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.10 | +0.00 | 0.1818 | 0.0910 | 0.0030 | 0.0030 | 4.3e-19 |
| 0.10 | +0.70 | 0.1818 | 0.0644 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.50 | -0.70 | 0.6667 | 0.3339 | 0.1484 | 0.1484 | 2.8e-17 |
| 0.50 | +0.00 | 0.6667 | 0.3750 | 0.1667 | 0.1667 | 0.0e+00 |
| 0.50 | +0.70 | 0.6667 | 0.3339 | 0.1484 | 0.1484 | 2.8e-17 |

## T8: Promotion assessment / barrier reduction

### Promotions

| barrier | promoted to | topological datum | status |
|---|---|---|---|
| `B1_closure_quantum` | winding θ-term S_top = 2π·n | winding number of phase map S¹ → U(1) | **PROMOTED** |
| `B2_antipodal_Z2` | RP³ = S³/Z₂ quotient + non-trivial spin structure | π₁(RP³) = Z₂ deck transformation + spin structure (T²=−I) | **PROMOTED** |

### Residual barriers (3 of the original 5)

  - **B3_boundary_conditions**: hard-wall Dirichlet from T²=−I (consistent with non-trivial spin structure, but BC derivation from bulk still separate)
  - **B4_dimensional_bridge**: 2π winding is metric-free, but physical length 2πR and mass scale still need m_e anchor
  - **B5_5d_to_4d_reduction**: unaffected; radial reduction producing F² remains unconstructed

Scaffold barrier count: **5 → 3**.

## Verdict

**PROMOTION_SUCCEEDS.** PROMOTION SUCCEEDS. Both B1 (closure quantum) and B2 (antipodal Z₂) are promoted from imposed constraints to terms in a single topological/discrete action sector:

    RP³ = S³/Z₂  +  non-trivial spin structure (T² = −I)  +  2π winding θ-term

B2 → quotient + spin structure: σ is a free involution (σ² = id, no fixed points), so S³/Z₂ = RP³ is a smooth manifold with π₁(RP³) = Z₂ (the non-contractible loop squared lifts to a contractible S³ path). The spinor transport T = iσ_y with T² = −I selects the non-trivial of the two RP³ spin structures (4π-periodic / antiperiodic spinors). The antipodal Z₂ is now the deck transformation plus a spin-structure choice — topological data, not an imposed symmetry.

B1 → winding θ-term: the closure quantum is the winding number of the phase map S¹ → U(1), ∮dφ = 2π·n, a topological total-derivative term S_top = 2π·n. It is integer-quantized by single-valuedness and metric-independent (reparam-invariant), and being a total derivative it does not modify the local EOM — variationally consistent with the smooth Einstein–Maxwell–Dirac–throat action.

The two are unified by the S³ → RP³ double cover: a great circle on S³ (length 2π = closure quantum, B1) is exactly a double traversal of the non-contractible RP³ loop (length π, the π₁(RP³) = Z₂ generator, B2). With both as action data, K(x) = 2x/(1+x) and Q(x, c) — hence F² = K²·Q — follow from the topological sector + stationary action to machine precision, with the closure quantum and antipodal symmetry no longer imposed.

Scaffold barrier count reduced 5 → 3. Residual: B3 (boundary conditions), B4 (dimensional bridge — the 2π winding is metric-free but the physical scale still needs m_e), B5 (5D→4D reduction, unaffected).

## What this leaves open

- **B3 — boundary conditions**: the hard-wall Dirichlet at the throat is consistent with the non-trivial spin structure identified here, but deriving the BC from the bulk action is a separate step.
- **B4 — dimensional bridge**: the 2π winding is metric-independent (purely topological), but the physical great-circle length 2πR and the mass scale still need the m_e anchor.
- **B5 — 5D → 4D reduction**: unaffected by this promotion; the radial reduction producing F² remains unconstructed.
