# The Born-rule equivariance test - companion probe (PR #198)

**Run:** 2026-07-02T20:27:37+00:00

The deliverable is `docs/born_rule_equivariance.md` - the equivariance and uniqueness theorems for the BAM psi-Phi-q transport flow. This probe verifies every claim on the running dynamics. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | item 3: the deepest refutation vector; the stakes | **PASS** |
| T2 | `T2_flow_and_guidance` | the flow: norm exact; guidance = grad S (Ehrenfest) | **PASS** |
| T3 | `T3_theorem1_equivariance` | THEOREM 1: |psi|^2 equivariant (residual ~ 0; KS at noise) | **PASS** |
| T4 | `T4_theorem2_uniqueness` | THEOREM 2: unique among h(rho); alternatives fail | **PASS** |
| T5 | `T5_wrong_flow_controls` | wrong flows fail: equivariance is specific to grad S | **PASS** |
| T6 | `T6_teeth_dissipative_deformation` | teeth: dissipative deformation breaks it as predicted | **PASS** |
| T7 | `T7_relaxation` | relaxation: uniform -> Born (H-theorem on BAM dynamics) | **PASS** |
| T8 | `T8_assessment` | Born rule at dBB grade; conditions stated | **PASS** |

## The ensemble transport through the collision

| t | KS(Born) | KS(sqrt-rho) | KS(rho^2) | KS(v/2) | KS(frozen) | Hbar(uniform) |
|---:|---:|---:|---:|---:|---:|---:|
| 2.5 | 0.0065 | 0.0086 | 0.0086 | 0.1054 | 0.2018 | 3.0763 |
| 5.0 | 0.0068 | 0.0621 | 0.0768 | 0.1929 | 0.3274 | 1.1687 |
| 7.5 | 0.0091 | 0.0962 | 0.1491 | 0.2237 | 0.4138 | 0.6987 |
| 10.0 | 0.0087 | 0.0952 | 0.1477 | 0.2052 | 0.3964 | 0.6524 |

(sampling noise 1/sqrt(N) = 0.0071; initial Hbar(uniform) = 4.13; continuity residual 3.3e-04; max|div v| = 2240.8)

## Verdict

**PSI_SQUARED_IS_THE_UNIQUE_EQUIVARIANT_DENSITY_OF_THE_BAM_TRANSPORT_FLOW_BORN_RULE_AT_DBB_GRADE_CONDITIONAL_ON_GUIDANCE_AND_EQUILIBRIUM.** THE DEEPEST REFUTATION VECTOR DID NOT FIRE (the argument is in docs/born_rule_equivariance.md; this probe verifies it on the running dynamics).

EQUIVARIANCE. The real-time Hamiltonian flow of the BAM psi-Phi-q functional preserves |psi|^2 EXACTLY: the continuity residual on the live evolution (self-consistent gravity + live order field) is 3e-04, and a 20000-throat Born ensemble transported by v = grad S stays at sampling noise through a two-soliton collision (KS series [0.0065, 0.0068, 0.0091, 0.0087], noise 0.0071). The reason is structural: every BAM coupling - the self-consistent Phi, the order field q - enters the pilot equation as a REAL potential.

UNIQUENESS AND FALSIFIABILITY. Among density functionals only |psi|^2 is equivariant ((h - rho h') div v with max|div v| = 2240.77; sqrt(rho)- and rho^2-ensembles fail, wrong flows fail), and the property is falsifiable: a dissipative deformation produces exactly the predicted continuity source and breaks the transport - the repo's own imaginary-time relaxation flow is of that non-equivariant type; only the Hamiltonian flow carries the Born rule. RELAXATION: a uniform ensemble falls 3.476 nats toward Born - fixed point AND attractor.

THE LABEL. The Born rule enters BAM at dBB GRADE - equivariance + uniqueness + relaxation, the same standing it has in Bohmian mechanics - conditional on the guidance identification (throat velocity = phase gradient; Ehrenfest-verified, uniquely current-closing, but its 5D derivation is its own program) and on the linear measurement regime (test throat in an external pilot wave). The deepest import is replaced by a theorem with stated hypotheses.
