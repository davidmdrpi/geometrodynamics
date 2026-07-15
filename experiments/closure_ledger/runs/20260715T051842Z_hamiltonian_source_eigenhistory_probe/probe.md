# The Hamiltonian source eigenhistory (PR #219)

**Run:** 2026-07-15T05:18:42+00:00

The deliverable is `docs/hamiltonian_source_eigenhistory.md` - the imposed source law replaced by an explicit Hamiltonian, and U(X)X = X solved with total-energy closure. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: derive the source, solve jointly | **PASS** |
| T2 | `T2_source_derived` | the source from H: unitary/reactive by construction | **PASS** |
| T3 | `T3_ring_modes` | the ring: the gap resolved; a mode BRANCH w*(A) | **PASS** |
| T4 | `T4_joint_solve` | U(X)X = X + energy closure: residuals at 1e-14 | **PASS** |
| T5 | `T5_energy_closure` | flux constant; zero source power; 1e4-pass persistence | **PASS** |
| T6 | `T6_stability` | stability eigenvalues: all on the unit circle | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_EIGENHISTORY_SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_CLOSURE_AND_ITS_STABILITY_SPECTRUM_IS_MARGINAL`

ESTABLISHED (the argument is in docs/hamiltonian_source_eigenhistory.md).

THE SOURCE, DERIVED. H = p^2/2 + w0^2 q^2/2 + mu q^4/4 + g q u(0): unitary scattering as a CONSEQUENCE (2e-16), reactive (0e+00), and #218's imposed law recovered as the weak-coupling limit (ratios [0.9999999884665239, 0.9999999997047431]).

THE RING. The reflecting source makes the loop a true two-direction ring: real trace, the tr > 2 gap segment resolved (2.7324-2.7390, max 2.0006 - cleaning up #217's coarse-scan tangency), and the homogeneous condition defining a nonlinear mode BRANCH w*(A) (span 9.8e-05): homogeneity alone cannot fix the amplitude.

THE JOINT SOLVE. U(X)X = X TOGETHER WITH total-energy closure: (w*, U0*) = (2.732375, 0.912168), a* = -0.1963; FIXED-POINT RESIDUALS: Newton (-9e-16, 5e-14), U(X)X = X 2e-14, source slaving 1e-16; energy partition E_field = 11.0064 + E_source = 0.1706 = E0 to 5e-14. Flux constant around the ring (2e-16); 1e4-pass persistence (drift 4e-11).

THE STABILITY EIGENVALUES: [1.000000453707063, 0.0] , [0.9999995465499416, 0.0] , [0.9999881564514221, 0.004866903856260777] , [0.9999881564514221, -0.004866903856260777] - moduli within 5e-07 of the unit circle (the marginal pair a Jordan block at 1, reciprocal product to 3e-10) - the conservative eigenhistory is marginal: Novikov-passive, no runaway; perturbations shear at most linearly (log-log slope 1.00, late ratio 1.82) - polynomial, never exponential.
