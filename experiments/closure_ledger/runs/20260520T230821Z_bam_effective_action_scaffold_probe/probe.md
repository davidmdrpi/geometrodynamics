# Covariant BAM effective-action scaffold probe

**Run:** 2026-05-20T23:08:21+00:00

Builds the candidate single variational principle for the three BAM targets (Hopf U(1) connection, Tangherlini boundaries, F² vertex) and identifies the mismatch terms preventing full closure. A scaffold + barrier map, not a closure claim.

## Candidate action

```
S_BAM = ∫_{M₅} √(−g₅)[ (R₅−2Λ₅)/2κ₅ − ¼F_{MN}F^{MN} + ψ̄(iΓ^M D_M − m)ψ + L_throat ] + S_∂[hard walls] + S_closure[∮A = 2πn]
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_action_scaffold_structure` | action stated; target → sector → variation map | **PASS** |
| T2 | `T2_gauge_sector_closes` | A_φ = ½cos χ (diff 0.0e+00); c₁ = 1.000000 (err 1.3e-08) — CLOSES | **PASS** |
| T3 | `T3_throat_sector_closes_with_constraint` | F² = K²·Q (max diff 2.8e-14); requires B1 + B2 | **PASS** |
| T4 | `T4_gravity_sector_partial` | f(r) = 1−(r_s/r)² (diff 0.0e+00); ΔR = 0.5200 imposed (B3) | **PASS** |
| T5 | `T5_dimensional_bridge_barrier` | scale-complete modulo m_e (1 anchor) | **PASS** |
| T6 | `T6_5d_to_4d_reduction_barrier` | 5D→4D reduction NOT constructed (B5) | **PASS** |
| T7 | `T7_barrier_ledger` | 1 sector closes, 5 barriers in 3 families | **PASS** |

## T1: Action scaffold structure

Action terms:
  - **G_gravity**: (1/2κ₅)(R₅ − 2Λ₅)  — 5D Einstein-Hilbert + cosmological
  - **A_gauge**: −¼ F_{MN}F^{MN}  — U(1) gauge field strength
  - **D_dirac**: ψ̄(iΓ^M D_M − m)ψ  — bulk fermions
  - **T_throat**: L_throat[Φ, g, A]  — throat / closure scalar
  - **B_boundary**: S_∂[hard walls]  — topological boundary conditions
  - **C_closure**: S_closure[∮A = 2πn]  — closure-quantum constraint

Target → sector → variation:
  - `A_phi_hopf_connection`: sector (A) U(1) gauge; δS/δA = 0 on S³
  - `Delta_R_tangherlini_boundaries`: sector (G) gravity + (B) BCs; δS/δg = 0 + boundary data
  - `F_squared_compton_vertex`: sector (T) throat + (C) closure; δS/δΦ = 0 + ∮A = 2πn

## T2: Gauge sector — CLOSES

| χ/π | A_φ | ½cos χ | diff |
|---:|---:|---:|---:|
| 0.000 | +0.500000 | +0.500000 | 0.0e+00 |
| 0.125 | +0.461940 | +0.461940 | 0.0e+00 |
| 0.250 | +0.353553 | +0.353553 | 0.0e+00 |
| 0.375 | +0.191342 | +0.191342 | 0.0e+00 |
| 0.500 | +0.000000 | +0.000000 | 0.0e+00 |

First Chern number c₁ = **1.000000** (error 1.29e-08). The Hopf connection is the homogeneous Maxwell solution; this sector derives from δS/δA = 0 with no external input.

## T3: Throat sector — CLOSES WITH CONSTRAINT

| x | cosθ | K | Q | K²·Q | F² closed | diff |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | -0.70 | 0.1818 | 0.0644 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.10 | +0.00 | 0.1818 | 0.0910 | 0.0030 | 0.0030 | 4.3e-19 |
| 0.10 | +0.70 | 0.1818 | 0.0644 | 0.0021 | 0.0021 | 0.0e+00 |
| 0.50 | -0.70 | 0.6667 | 0.3339 | 0.1484 | 0.1484 | 2.8e-17 |
| 0.50 | +0.00 | 0.6667 | 0.3750 | 0.1667 | 0.1667 | 0.0e+00 |
| 0.50 | +0.70 | 0.6667 | 0.3339 | 0.1484 | 0.1484 | 2.8e-17 |

F² = K²·Q reconstructs to 2.8e-14, but requires two imposed inputs: **B1_closure_quantum_2pi, B2_antipodal_Z2** (closure quantum 2π and antipodal Z₂) — neither is a local Lagrangian density.

## T4: Gravity sector — PARTIAL

| r | f(r) = 1−(r_s/r)² | f from V_tangherlini | diff |
|---:|---:|---:|---:|
| 0.80 | -0.562500 | -0.562500 | 0.0e+00 |
| 1.20 | +0.305556 | +0.305556 | 0.0e+00 |
| 1.50 | +0.555556 | +0.555556 | 0.0e+00 |
| 2.00 | +0.750000 | +0.750000 | 0.0e+00 |

R_MID = 1.0, R_OUTER = 1.26, R_INNER = 0.74, ΔR = 0.5200. The metric is bulk-derived; the boundaries are imposed (B3).

## T5: Dimensional bridge barrier (B4)

  - `F_squared`: dimensionless (function of x = ω′/ω and c = cos θ)
  - `chern_c1`: pure integer = 1
  - `closure_quantum`: 2π (dimensionless action units)
  - `K_pade`: dimensionless (harmonic-mean ratio)
  - `Q_polarization`: dimensionless (helicity-channel sum)

Dimensional bridge: `ℏ = m_e · R_MID · c  (m_e is the single external anchor)`. External anchors: 1 (m_e).

## T6: 5D → 4D reduction barrier (B5)

The reduction map producing F² from the 5D action would require:
1. Radial-mode integration: integrate the 5D fields over the Tangherlini radial channel [R_INNER, R_OUTER] against the bound modes u_{l,n}(r*) to obtain 4D effective couplings.
2. Vertex projection: project the bulk Dirac–gauge–throat coupling onto the 4D Compton kinematic variables (x, c) to recover F²(x, c).
3. Consistency of the radial-mode spectrum (which sets lepton/quark masses) with the F² vertex normalisation — the two currently live in separate sub-threads (closure-ledger vs tree-QED) and are not connected by an explicit reduction.

**Constructed: False** — the largest structural gap.

## T7: Barrier ledger / closure assessment

### Sector status

| sector | target | status | external input |
|---|---|---|---|
| `A_gauge` | A_φ = ½cos χ | **CLOSES** | none |
| `T_throat` | F² = K²·Q | **CLOSES_WITH_CONSTRAINT** | B1 (closure quantum) + B2 (antipodal Z₂) |
| `G_gravity` | ΔR = R_OUTER − R_INNER | **PARTIAL** | B3 (boundary conditions) |

### Mismatch terms

| barrier | type | mismatch |
|---|---|---|
| `B1_closure_quantum` | topological | one global topological input per closed orbit |
| `B2_antipodal_Z2` | discrete | one discrete symmetry imposed as identification |
| `B3_boundary_conditions` | boundary | boundary data imposed by non-orientable topology |
| `B4_dimensional_bridge` | scale | one external scale anchor (m_e) |
| `B5_5d_to_4d_reduction` | reduction | structural map not constructed |

### Barrier families

  - **topological**: B1_closure_quantum, B2_antipodal_Z2
  - **boundary**: B3_boundary_conditions
  - **scale_and_reduction**: B4_dimensional_bridge, B5_5d_to_4d_reduction

## Verdict

**SCAFFOLD_WITH_BARRIERS.** SCAFFOLD WITH BARRIERS. A single covariant action
  S_BAM = ∫_{M₅} √(−g₅)[ (R₅−2Λ₅)/2κ₅ − ¼F² + ψ̄(iΓ·D−m)ψ + L_throat ] + S_∂ + S_closure
unifies the three BAM targets STRUCTURALLY, but full closure is blocked by five identified mismatch terms.

Sector status:
  (A) gauge   — CLOSES. A_φ = ½cos χ is the homogeneous Maxwell solution on the Hopf bundle; c₁ = 1 exactly; no external input.
  (T) throat  — CLOSES WITH CONSTRAINT. F² = K²·Q reconstructs from the equal-action splitting, but requires the closure quantum (B1) and antipodal Z₂ (B2).
  (G) gravity — PARTIAL. The Tangherlini metric is bulk-derived, but ΔR boundary data is imposed (B3).

Mismatch terms (barriers), in three families:
  TOPOLOGICAL: B1 closure quantum (∮A = 2πn, not a local density); B2 antipodal Z₂ (T = iσ_y, discrete quotient).
  BOUNDARY: B3 hard-wall Dirichlet from T²=−I (topological, not from δS/δg).
  SCALE / REDUCTION: B4 dimensional bridge (ℏ = m_e·R_MID·c, one external anchor); B5 5D→4D reduction producing F² (not constructed, the largest gap).

The barriers are not five unrelated patches: they are a recognisable set — topological quantisation, discrete antipodal identification, the dimensional anchor, and the unbuilt dimensional reduction — that recurs throughout the BAM programme. Full closure requires (i) promoting the closure quantum + antipodal Z₂ to a topological / discrete sector of the action (e.g. a θ-term + S³/Z₂ quotient), (ii) deriving the hard-wall BC from the non-orientable bulk, and (iii) constructing the 5D→4D radial reduction that projects the bulk coupling onto F²(x, c).

## What full closure would require

1. **Promote B1 + B2 to a topological / discrete sector** — e.g. a θ-term enforcing ∮A = 2πn and an explicit S³/Z₂ antipodal quotient, so the closure quantum and antipodal symmetry become part of the action rather than imposed constraints.
2. **Derive the hard-wall boundary condition (B3)** from the non-orientable bulk topology (T² = −I), rather than imposing Dirichlet by hand.
3. **Construct the 5D → 4D reduction (B5)** — integrate the bulk fields over the Tangherlini radial channel and project onto the Compton variables (x, c) to produce F²(x, c) from first principles. This is the largest gap.
4. **Close the dimensional bridge (B4)** — derive m_e (or the R_MID scale) from within the action, removing the single external anchor.
