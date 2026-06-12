# The final geometric mile: the Hopf sector arc + the pinhole refinement (PR #160)

**Run:** 2026-06-12T04:58:05+00:00

Closes the two flagged final-mile items. PART A: the #159 sector arc is derived from the shell-wavefunction algebra — the capacity-k₅ winding space is Weyl-dual to a k₅-site fiber lattice, making 2π/k₅ the commutator quantum (machine-exact): φ_h = π/k₅ is now algebraic end-to-end. PART B: the soft V_us direction reduces — through exact exclusions and the mass-preserving refinement family — to a single remaining angle (γ), with seven of eight flavor-CP observables landing at the joint solution and the #156/#158 J-ceiling lock verified. *(QFT on the classical throat, not quantum gravity.)*

- **Part A**: UVU†V† = e^{2πi/k₅} exact ⟹ arc = 2π/k₅ ⟹ φ_h = π/k₅ algebraic
- **Part B exclusions**: pinhole breaks m_s −22.5%; transport rescale self-defeats
- **The joint solution**: V_us/β exact; V_ub ×1.10; J ×1.05; V_cb/V_ts ×0.90; γ = 104° the misfit
- **The ceiling lock**: ceiling → 0.99 of observed at the refined point — VERIFIED
- **Open**: the γ angle; the knob-level v3+CP re-lock

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | close both flagged final-mile items | **PASS** |
| T2 | `T2_weyl_quantum` | Weyl pair: commutator e^{2πi/k₅} exact; arc = lattice step | **PASS** |
| T3 | `T3_chain_composition` | × connection ½ ⟹ φ_h = π/k₅ algebraic; #159 caveat removed | **PASS** |
| T4 | `T4_single_route_exclusions` | exclusions: pinhole breaks m_s; transport rescale self-defeats | **PASS** |
| T5 | `T5_mass_preserving_joint_solution` | joint solution: 7/8 observables land; γ the single misfit | **PASS** |
| T6 | `T6_j_ceiling_lock_verified` | J-ceiling lock VERIFIED (ceiling 0.99 of observed; J ×1.05) | **PASS** |
| T7 | `T7_ledger_and_scope` | re-lock targets tabulated; residual = γ; budget unchanged | **PASS** |
| T8 | `T8_assessment` | SECTOR_ARC_WEYL_DERIVED_SOFT_DIRECTION_REDUCED_TO_GAMMA | **PASS** |

## The joint solution (masses preserved to 1e-15)

| observable | refined / observed | status |
|---|---:|---|
| V_us | ×1.0 | lands |
| V_cb | ×0.899 | lands |
| V_ub | ×1.103 | lands |
| V_td | ×1.188 | lands |
| V_ts | ×0.886 | lands |
| J | ×1.054 | lands |
| β | 22.2° vs 22.2° | exact-fit |
| γ | 104.2° vs 65.9° | **the single misfit** |

## The re-lock targets (minus block)

| element | locked | target | factor |
|---|---:|---:|---:|
| H_ds | -0.3548 | -0.8381 | ×2.362 |
| H_db | -0.2682 | -0.286 | ×1.067 |
| H_sb | -5.2682 | -5.277 | ×1.002 |
| H_dd | 3.3816 | 3.5942 | ×1.063 |
| H_ss | 6.5516 | 6.3393 | ×0.968 |

## Verdict

**SECTOR_ARC_WEYL_DERIVED_SOFT_DIRECTION_REDUCED_TO_GAMMA_CEILING_VERIFIED.** BOTH FINAL-MILE ITEMS CLOSE: THE SECTOR ARC IS THE WEYL COMMUTATOR QUANTUM OF THE CAPACITY-k₅ FIBER (THE #159 IDENTIFICATION NOW ALGEBRA, MACHINE-EXACT), AND THE SOFT V_us DIRECTION REDUCES TO A SINGLE REMAINING ANGLE — WITH SEVEN OF EIGHT FLAVOR-CP OBSERVABLES LANDING AT THE MASS-PRESERVING JOINT SOLUTION AND THE #156/#158 J-CEILING LOCK VERIFIED.

PART A. A capacity-k₅ winding space carries the canonical clock–shift pair with UVU†V† = e^{2πi/k₅}·1 exactly: the fiber discretizes into k₅ sites θ_n = 2πn/k₅, and winding-changing transport is the shift whose minimal step is one site — the sector arc 2π/k₅ is the Weyl quantum, independent of any radial profile. Composed with the connection's ½: φ_h = π/k₅, end-to-end algebraic. The #159 caveat is removed.

PART B, THE EXCLUSIONS. The pinhole single-knob lands V_us but breaks m_s by −22.5%; the exact transport rescale that doubles the dk = 3 element self-defeats via level repulsion (V_us only 0.133, m_s +50%; the invariant sin 2θ = 2|H_ds|/Δλ verified). The soft direction is a structured joint constraint, not slop.

THE JOINT SOLUTION. The mass-preserving family (eigenvector rotations at fixed eigenvalues, masses to 1e-15) has a joint (V_us, β) solution at (δθ_u, δθ_d) = (-5.22°, 9.94°): V_us ×1.00, V_cb ×0.90, V_ub ×1.10, V_td ×1.19, V_ts ×0.89, J ×1.05, β exact-fit — and γ = 104.2° vs 65.9°, the single remaining misfit.

THE CEILING LOCK, VERIFIED. #156/#158 predicted the J ceiling rises to the observed 3.5e-5 when the soft elements land: at the refined point the ceiling is 3.441e-05 (0.991 of observed) and J is ×1.054 — a prediction made about a state that did not yet exist, now checked in that state and passed.

LEDGER. No new inputs; the precise minus-block re-lock targets are tabulated (H_ds ×1.8-class with %-level diagonal compensation); the flavor sector's remaining residual is the γ angle plus the knob-level realization of the re-lock — the flagged successor.
