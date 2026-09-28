# Two gravitating regions: prospective initial-data and tracking test

Baseline: merged PR #311, bf74cba62dd39e17b53ba7ca629f8b9a39769b62.
This specification is published in a draft PR before the new solves. No
measurement has been run. The following analytic formulas are design work.

## Scope

Add two separately seeded, separately measured regions to one nonlinear
Einstein/quartet constraint solution. Compute their instantaneous momentum
rates from the coupled field equations and a separate stress balance. This
is a time-symmetric initial slice and its first time jet, not a trajectory,
finite-time momentum-exchange experiment, two black holes, or two mouths.
The previous spherical evolution cannot evolve these data. Traversability
is not required. No result here establishes emergent quantum mechanics.

Use the same Einstein-frame four-scalar action as #310/#311:
G_AB=delta_AB/f+phi_A phi_B/(6 f^2), f=1-|phi|^2/6,
U=3/(2 f^2), Einstein equation G_munu=T_munu. On the unit S3 cover set
phi=q(x)x, q=(sqrt(3)/2)[1+e1 b(c1.x)+e2 b(c2.x)],
b(z)=cosh(8z)/cosh(8), c1=(1,0,0,0),
c2=(cos(1.2),sin(1.2),0,0). Set scalar normal velocity and K_ij to zero.
The two locations are distinct on RP3; each has two antipodal cover images.
Metric and stress descend. The odd quartet descends only with the internal
phi -> -phi twist (a target-space isometry), not as four ordinary even
scalars on RP3. All raw solves are on S3; paired invariant integrals are
halved when reporting quotient inventories. This is not the old handle.

With g=psi^4 gamma, solve the full nonlinear Hamiltonian equation

    -8 Delta_gamma psi + (6-Sbar) psi - 2 U psi^5 = 0,
    Sbar=|grad q|^2/f^2+3 q^2/f.

The momentum constraint is identically satisfied by these time-symmetric
data. It is not a nonzero-momentum constraint test. The round solution is
psi=f0^(1/4), f0=7/8. No delta functions or external probe masses are added.
Newton solves retain the homogeneous mode; no stabilizing projection.

## Numerical schedule and controls

Exploit SO(2) invariance rotating x2,x3, not spherical symmetry. Use real
disk Zernike eigenfunctions of Delta_S3, eigenvalues -n(n+2), even n through
12,20,28. Gauss-Legendre in r^2 and periodic angle quadrature use
(20,64),(28,88),(36,112); independent validation uses (64,192).
Newton starts from the round solution, positivity-preserving line search,
at most 30 iterations, projected residual <1e-11. Stop on failure or
f<=0.1 or psi<=0.1; do not silently drop a case.

Cases: round (0,0), A-only (.02,0), B-only (0,.03), pair (.02,.03),
A-perturbed (.021,.03), B-perturbed (.02,.031). Solve all at each resolution.
Compare fields on the independent grid and record unprojected physical
Hamiltonian residual. Finest normalized max residual <1e-7 and successive
field differences <1e-7 are the constraint gates. Also evaluate physical
curvature through finite differences of the coordinate metric at off-grid
points, independent of the spectral Laplacian, using successively halved
steps. These are numerical checks, not rigorous existence certificates.

## Independent tracking and momentum convention

Assign permanent IDs A and B to disjoint paired caps of radius .45 around
+/-c1 and +/-c2. Their smooth fixed windows are w=max(h,0)^4,
h=((c.x)^2-cos(.45)^2)/sin(.45)^2. They never share an oriented cut.
Record separately proper window volume, energy inventory, signed excess
relative to the round background with the same window, and a positive
energy centroid using the lift sign(c.x)x. Centroids are coordinate/frame
observables, not gauge-invariant worldlines or material boundaries.
Record all six momentum inventories P_ab=int w j_i X_ab^i dV/2,
X_ab=x_a partial_b-x_b partial_a, using the same specified S3 Killing fields
for both regions. They are rotation charges, not ADM linear momenta. Do
not impose P_B=-P_A. All initial P values vanish because normal velocities
vanish; their rates need not vanish or be equal and opposite.

Compute jdot from the sigma-model acceleration in unit lapse, zero shift;
independently compute each rate using

    Pdot_X = int [S^j_i X^i partial_j w
                  + 2 w X(log psi) tr_g S] dV/2.

The first term is the smooth-window stress flux; the second is metric work.
Also measure the complementary window 1-wA-wB so the matter ledger includes
the exterior. Require the two rate calculations to agree within 1e-5
absolute on the validation grid, and publish quadrature refinement. A matter
ledger with metric work is not a local gravitational momentum density or a
complete gravitational recoil balance. There are no inner boundaries here.

Separate seed perturbations must produce separately measured responses;
report the 2x2 finite-difference inventory response and metric response in
the other cap. An elliptic cross-response compares different initial data;
it is NOT a causal propagation signal. Do not select gates on response sign.

At K=0, every sphere has future expansions theta_+=H, theta_-=-H, so no
strictly future-trapped sphere exists on this slice. This analytic exclusion
does not change the #309 result for different, non-time-symmetric neck data.

## Deliverables and next gate

Publish solver, snapshot tracker, reproducible runner, raw coefficients,
residuals, controls, source hashes, tests, and an honest verdict. This adds
two gravitating initial regions with independent observables. A nonspherical
constraint-monitored evolution, transported moving windows, and an exterior
causal response experiment remain required before claiming two-body recoil.
