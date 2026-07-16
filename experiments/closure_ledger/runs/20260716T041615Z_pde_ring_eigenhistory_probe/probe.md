# The PDE-ring eigenhistory (PR #220)

**Run:** 2026-07-16T04:16:15+00:00

The deliverable is `docs/pde_ring_eigenhistory.md` - the full Hamiltonian successor to #219: the explicit time-domain source-field system, its periodic orbit, and the complete monodromy. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: the full Hamiltonian successor | **PASS** |
| T2 | `T2_hamiltonian_system` | H with g q u(0) conserved; the term load-bearing | **PASS** |
| T3 | `T3_periodic_orbit` | the literal periodic-orbit solve, machine residual | **PASS** |
| T4 | `T4_monodromy` | the complete monodromy: all unit circle, symplectic | **PASS** |
| T5 | `T5_pairs` | the source pair confirms #219 at 0.07% | **PASS** |
| T6 | `T6_energy_ledger` | the time-resolved ledger closes | **PASS** |
| T7 | `T7_spatial_convergence` | spatial convergence: O(dx^2), Richardson continuum | **PASS** |
| T8 | `T8_honest_scope` | honest scope | **PASS** |
| T9 | `T9_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_FULL_PDE_EIGENHISTORY_EXISTS_ITS_COMPLETE_MONODROMY_IS_UNIT_CIRCLE_SYMPLECTIC_AND_THE_219_REDUCTION_IS_CONFIRMED_BY_THE_FIELD_ITSELF`

ESTABLISHED (the argument is in docs/pde_ring_eigenhistory.md).

THE SYSTEM. H = H_field + p^2/2 + w0^2 q^2/2 + mu q^4/4 + g q u(0) on the periodic two-barrier ring, evolved by symplectic leapfrog: H INCLUDING THE INTERACTION TERM is conserved (bound 2.1e-05 over 100 periods, non-secular, O(dt^2) scaling ratio 4.0), while the ledger without it fails 99x worse in a single period - the term is load-bearing.

THE ORBIT. Gauss-Newton on the literal [X(T) - X(0); H - E0; p(0)] over the complete state: quadratic convergence to 2e-13; E0 = 20.834 hit to 0e+00; w_orbit = 2.669296 (nonlinear pull -3.00e-03 from the linear mode), dt-converged to 3e-07.

THE COMPLETE MONODROMY. 386 x 386: ALL eigenvalues on the unit circle to 4e-14 - no parametric instability of any field mode; det = 1 to 2e-13; the dx-weighted symplectic form preserved to 4e-15; the trivial pair to 9e-11.

THE PAIRS. The source pair (qp-weight 0.36 vs 0.008 for the next) rotates at 3.2124 - within 0.1% of #219's reduced-model 3.2102 (bare w0 = 3.2): the reduction is CONFIRMED by the field itself; the low field modes sit at their folded angles (worst gap 1.8e-02).

THE LEDGER. Time-resolved on the orbit: E_int oscillates about -0.0467 (negative), equal to the #219 harmonic-balance value to 0.0%; the full H is constant along the orbit to 5e-05.

THE CONVERGENCE. The complete solve repeated at N = 128/192/256/384 (ring, barriers, coupling, E0 fixed): the orbit frequency, q and u(0), the partition, the interpolated profile, the source pair, and the low Floquet angles ALL converge O(dx^2) - measured diff-to-finest ratios (exact prediction 6.40 / 2.40): w_orbit 6.39 / 2.40, profile 6.39 / 2.40; every monodromy stays on the unit circle (worst 1e-10); Richardson continuum w_orbit = 2.673464, source frequency = 3.2132 - 0.09% from #219's 3.2102.
