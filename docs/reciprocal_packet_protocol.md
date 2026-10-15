# Reciprocal characteristic packet: production protocol

Written after the analytic preflight #334 and before implementing or running
this evolution. The null-constraint compatibility issue is already anticipated;
this is not a blind discovery protocol. Preserve this file after production.

## Problem and orders

Use two independent vacuum Tangherlini exterior patches, b=1, c=1,
G_ab=kappa5 T_ab with kappa5=1. Propagate u=r^(3/2) Phi to first order in
amplitude A. Compute the spherical Misner-Sharp mass parameter response at
O(A²), with m=r²(1+4 r_U r_V/F) and M=3 Omega3 m/(2 kappa5).
The spherical scalar is constant on S3, not normalized to a unit harmonic;
physical energies include Omega3=2 pi². l=1 is a separate test-field control:
no exact nonspherical Misner-Sharp mass is assigned to it.

Incoming rectangle: U in [-8,0], V in [.5,4], with zero data on V=.5 and
u(-8,V)=A exp(-z²/2) exp(1-1/(1-y²)) for |y|<1, zero otherwise,
z=(V-2)/.25, y=V-2. This is a smooth compact packet on V in [1,3].
A=.01. A second solution with boundary value (V-2)u supplies an independent
symplectic variation. Do not claim that the evolved horizon trace is compact.

Read a(V)=u(0,V), then prescribe u_out(-1/V,0)=(-1)^l a(V).
The second outgoing characteristic is U_out=-2 with zero field. Its corner
agrees with the trace. Outgoing V runs [0,80]. Extract the field at r_obs=2
only where that worldline lies inside the computed domain. This is a
prescribed surface-transfer law, not a derived material interaction.

## Measurements and failure rules

1. Solve u_UV + (F V_l/(4f))u=0 using the corrected operator. Check the
   characteristic solver against a nonzero-potential exact solution and
   zero-potential characteristic addition. Compare resolutions N=64,128,256,512;
   use N incoming U intervals and 2N V intervals, and 2N intervals in each
   outgoing direction. Report measured orders; expect second-order convergence.
2. For l=0 compute both components of the Einstein mass one-form and integrate
   around two rectangular paths. Require relative path discrepancy <1e-3 on
   the finest grid and decreasing with refinement. If necessary one additional
   doubling may resolve a failing numerical tolerance; do not change the gate.
3. Check characteristic symplectic conservation using two independent real
   solutions and all four edges. Require relative residual <1e-4 at the finest
   grid, with refinement. Report reciprocal seam energy and oriented pairing
   separately from full-domain balances.
4. Apply passive boosts lambda=.5,1,2,4 to grids, data, map and clock origin.
   Compare the same physical observer field; relative L2 difference <1e-10.
   Holding physical grid points fixed in this comparison tests covariance,
   not a second numerical discretization.
5. Check A -> A/2 and A -> 2A. The linear field and leading mass must scale
   as A and A². These are perturbative identities, not nonlinear validation.
6. Integrate the two horizon null Einstein constraints independently, with
   equal initial radius and zero expansion before the packet. Compare induced
   radii at matched generators. A nonzero O(A²) difference surviving refinement
   is a **physical failure of a perturbatively matched fixed reciprocal null
   continuation** with the canonical scalar, not a numerical solver failure.

The packet and mass response on each side may be computed even if step 6 fails.
Do not call their pair a self-consistent backreacted spacetime. Do not interpret
the quasi-local mass-parameter difference as an independently derived shell
Hamiltonian or a solved surface constitutive law. No finite-thickness resolution,
nonlinear completion, or quantum event is claimed. Preserve #334's original
protocol and results. Archive numerical summaries, waveforms, and reproducible
source/version hashes, distinguishing numerical gates from physical admission.
