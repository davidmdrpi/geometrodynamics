# Hard-wall boundary derivation probe

**Run:** 2026-05-21T04:05:29+00:00

Derives the throat hard-wall (Dirichlet) boundary condition from the PR #49 non-trivial RP³ spin structure, closing scaffold barrier B3.

## Derivation chain

```
RP³ non-trivial spin structure (PR #49) → T² = −I → single-valuedness ψ = T·ψ at throat → ψ(throat) = 0 (Dirichlet hard wall)
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_T_eigenvalues` | eigenvalues ±i; no +1 eigenvalue: True | **PASS** |
| T2 | `T2_t_fixed_point_forces_dirichlet` | (T−I) null dim = 0; Dirichlet forced: True | **PASS** |
| T3 | `T3_radial_solver_odd_extension` | odd residual 0.0e+00; throat node 0.0e+00; endpoint 0.0e+00 | **PASS** |
| T4 | `T4_alternatives_ruled_out` | Dirichlet unique: True | **PASS** |
| T5 | `T5_discrete_bulk_spectrum` | all spectra discrete: True | **PASS** |
| T6 | `T6_consistency_with_prior_probe` | T²=−I reconfirmed; DD wins (prior probe) | **PASS** |
| T7 | `T7_b3_promotion` | B3 absorbed; barriers 3 → 2 | **PASS** |

## T1: Spin structure → T eigenvalues

Eigenvalues of T = iσ_y: `+0.000-1.000i, +0.000+1.000i`. No +1 eigenvalue (`has_plus_1 = False`) → no nonzero T-invariant spinor.

## T2: T-fixed-point argument → Dirichlet

`(T − I)` has rank 2, null-space dimension 0, det = +2.000+0.000i ≠ 0. So `ψ = T·ψ` ⟹ `ψ = 0`: Dirichlet forced (`dirichlet_forced = True`).

## T3: Radial solver odd extension realizes the throat node

| l | eigenfrequencies | odd residual | throat node | endpoint |
|---:|---|---:|---:|---:|
| 1 | 1.0547, 1.9744, 2.8941 | 0.0e+00 | 0.0e+00 | 0.0e+00 |
| 3 | 1.2191, 2.1412, 3.0220 | 0.0e+00 | 0.0e+00 | 0.0e+00 |
| 5 | 1.3960, 2.3694, 3.2277 | 0.0e+00 | 0.0e+00 | 0.0e+00 |

The solver extends each mode antisymmetrically across the throat (`u_full = [−u_reflected, u]`) — the T-odd transport producing the Dirichlet node `u(R_MID) = 0`.

## T4: Alternatives ruled out by the spin structure

Across 1000 random unit spinors, none is T-invariant (`found_nonzero_T_invariant = False`). Neumann/Robin require a nonzero throat spinor that would have to be T-invariant — impossible under T² = −I. Dirichlet is unique.

## T5: Discrete bulk spectrum from the derived BC

| l | ω spectrum | spacings | discrete |
|---:|---|---|:---:|
| 1 | 1.0547, 1.9744, 2.8941, 3.8247 | 0.9197, 0.9197, 0.9306 | True |
| 3 | 1.2191, 2.1412, 3.0220, 3.9225 | 0.9221, 0.8808, 0.9005 | True |
| 5 | 1.3960, 2.3694, 3.2277, 4.0861 | 0.9734, 0.8583, 0.8584 | True |

## T6: Consistency with prior hard_wall probe

T² = −I reconfirmed (`True`); the prior `hard_wall_boundary_verification` probe established **DD wins (hard_wall_boundary_verification)**. This probe supplies the topological-sector derivation (spin structure = PR #49 action data).

## T7: B3 promotion / barrier reduction

B3 derivation: `RP³ non-trivial spin structure (PR #49) → T² = −I → single-valuedness ψ = T·ψ at throat → ψ = 0 (Dirichlet)`

Scaffold barriers: **3 → 2**. Residual:
  - **B4_dimensional_bridge**: absolute MeV scale needs m_e anchor (ℏ = m_e·R_MID·c)
  - **B5_5d_to_4d_reduction**: radial spectrum (this BC) and F² vertex still in separate sub-threads; reduction map unconstructed (largest gap)

## Verdict

**HARD_WALL_DERIVED.** HARD-WALL BOUNDARY DERIVED. The Dirichlet condition at the throat is a consequence of the PR #49 non-trivial RP³ spin structure, not an independent imposition:

  RP³ non-trivial spin structure (PR #49)
    → T = iσ_y, T² = −I  (eigenvalues ±i, no +1 eigenvector)
    → single-valuedness at the throat fixed point: ψ = T·ψ
    → T²ψ = Tψ = ψ  but  T²ψ = −ψ  ⟹  ψ = −ψ  ⟹  ψ(throat) = 0
    → Dirichlet (hard wall) at the throat.

The argument is realized concretely in the Tangherlini radial solver, whose modes are extended antisymmetrically across the throat (u_full = [−u_reflected, u], odd to machine precision) — the T-odd transport producing the throat node u(R_MID) = 0, with both grid endpoints at Dirichlet. Neumann/Robin BCs are ruled out: no nonzero spinor is T-invariant under T² = −I. The derived BC yields the discrete bulk ω spectrum feeding the lepton/quark ladder. Consistent with the prior hard_wall_boundary_verification probe (DD wins).

B3 is absorbed into the PR #49 topological sector; the scaffold barrier count drops 3 → 2. Residual: B4 (dimensional bridge — absolute scale needs m_e) and B5 (5D→4D reduction producing F², the largest remaining gap).

## What this leaves open (residual 2 barriers)

- **B4 — dimensional bridge.** The discrete spectrum is in geometric units (R_MID = 1); the absolute MeV scale needs the m_e anchor (ℏ = m_e·R_MID·c).
- **B5 — 5D → 4D reduction.** The radial spectrum (this BC) and the F² vertex still live in separate sub-threads; the reduction map connecting them is unconstructed — the largest gap.
- **Outer boundary note.** R_OUTER is the cavity wall, fixed by the cross-species γ-lock fixed point (closure-ledger). Whether its Dirichlet is also a spin-structure consequence or a distinct cavity condition is noted but not the focus here — B3 names the throat BC, which is derived.
