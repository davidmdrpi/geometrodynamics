# The Hamiltonian source eigenhistory (PR #219)

**Run:** 2026-07-15T06:20:01+00:00

The deliverable is `docs/hamiltonian_source_eigenhistory.md` - the imposed source law replaced by an explicit Hamiltonian, and U(X)X = X solved with total-energy closure. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: derive the source, solve jointly | **PASS** |
| T2 | `T2_source_derived` | the source from H: unitary/reactive by construction | **PASS** |
| T3 | `T3_ring_modes` | the ring: the gap resolved; a mode BRANCH w*(A) | **PASS** |
| T4 | `T4_joint_solve` | U(X)X = X + energy closure: residuals at 1e-14 | **PASS** |
| T5 | `T5_energy_closure` | flux constant; zero source power; 1e4-pass persistence | **PASS** |
| T6 | `T6_stability` | the FULL Hamiltonian Floquet spectrum ((q,p) evolved) | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_EIGENHISTORY_SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_CLOSURE_AND_ITS_STABILITY_SPECTRUM_IS_MARGINAL`

ESTABLISHED (the argument is in docs/hamiltonian_source_eigenhistory.md).

THE SOURCE, DERIVED. H = p^2/2 + w0^2 q^2/2 + mu q^4/4 + g q u(0): unitary scattering as a CONSEQUENCE (2e-16), reactive (0e+00), and #218's imposed law recovered as the weak-coupling limit (ratios [0.9999999884665239, 0.9999999997047431]).

THE RING. The reflecting source makes the loop a true two-direction ring: real trace, the tr > 2 gap segment resolved (2.7324-2.7390, max 2.0006 - cleaning up #217's coarse-scan tangency), and the homogeneous condition defining a nonlinear mode BRANCH w*(A) (span 9.8e-05): homogeneity alone cannot fix the amplitude.

THE JOINT SOLVE, CORRECTED LEDGER. U(X)X = X TOGETHER WITH E_field + E_source + <g q u(0)> = E0: (w*, U0*) = (2.732375, 0.914367), a* = -0.1967; FIXED-POINT RESIDUALS: Newton (-2e-16, 1e-13), U(X)X = X 2e-14, source slaving 0e+00; energy partition E_field = 11.0596 + E_source = 0.1714 + E_int = -0.05397 = E0 to 1e-13. Flux constant around the ring (2e-16); 1e4-pass persistence (drift 3e-11).

THE FULL HAMILTONIAN STABILITY SPECTRUM. The Duffing (q, p) evolved as INDEPENDENT variables in the 4x4 variational monodromy about the shooting-refined periodic orbit (residual 2e-14; NNM frequency consistent with the full ring to 8.8e-05): eigenvalues [0.9999999999999949, 7.603780428470208e-09] (double - the Floquet-trivial pair, to 8e-09) and [0.4540579773262343, 0.8909721394220808] - THE SOURCE PAIR, rotating at its dressed frequency 3.2102 (bare w0 = 3.2), invisible to the slaved harmonic-balance map; ALL moduli within 2e-14 of the unit circle; det M = 1 to 6e-14, symplectic to 5e-14, energy drift 6e-11: the conservative eigenhistory is marginal in the FULL Hamiltonian sense - Novikov-passive, no runaway.
