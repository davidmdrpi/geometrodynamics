# Vacuum momentum completion on a compact twisted handle: prospective test

Date: 2026-09-15. Baseline: merged #300,
`4cd86541d3b836b35561b0c4a3a54629d28851cd`.
This specification is to be published before implementing or running its
new symbolic checks, elliptic solves or numerical controls. Analytic
predictions made while designing the experiment are disclosed below.

## Question and scope

Can a compact slice with a spatial antipodal handle carry a nonspherical,
nonzero extrinsic curvature satisfying the vacuum momentum constraint and
the nonlinear Hamiltonian constraint simultaneously?

This is a new vacuum Einstein initial-data benchmark, not an extension of
#300's homogeneous four-scalar ansatz. No new matter or apparatus law is
introduced. It does not retrofit a mouth into #300, reproduce the earlier
round-S3 exterior, establish traversability, or evolve a worldline crossing.
A nontrivial momentum solve on this topology is the intended limited advance.
The tensor seed is wave-like free gravitational initial data; propagation,
directional radiation flux and particle momentum transfer require evolution.

## Topology, action and conventions

Use the Einstein-Hilbert action with positive cosmological constant and no
matter. Units fix the conformal sphere radius to one. The compact spatial
manifold is the mapping torus

    (s+2pi, n) ~ (s, -n), n in S2,
    gbar = ds^2 + dtheta^2 + sin(theta)^2 dphi^2, Rbar=2.

It is the nonorientable S2 bundle over S1; its orientable double cover has
s period 4pi. This is a spatial twist, with an unchanged timelike normal.
It is not a derivation of a time-reversing wormhole identification. A cut
through a fiber exposes two sphere faces that are glued by the antipodal
map. There is no physical boundary, thin shell, imposed force or surface
stress. Smooth tensor descent replaces interface source conditions.
This realizes a compact twisted handle topology, not a geometrically round
S3 with small antipodal holes and a separately resolved exterior.

Let physical data be

    gamma_ij = psi^4 gbar_ij,
    K_ij = psi^-2 Abar_ij, tr_gamma K=0, psi>0.

Use K_ij=-1/2 L_n gamma_ij. The vacuum constraints are

    R(gamma)-|K|_gamma^2=2 Lambda,
    D_j K^j_i=0.

Their conformal equations are

    divbar Abar=0,
    Delta_bar psi - psi/4 + |Abar|^2 psi^-7/8
                            + Lambda psi^5/4=0.

The s=0 fiber will be tested as a reflection-symmetric minimal section;
its area and nearby section areas are diagnostics, not a proof of stable
minimality against every deformation. At H=0 its future null expansions
are both -tr_surface K. Record these: do not call a minimal section a
traversable mouth, an apparent horizon or a measurement event.

## A nontrivial momentum equation

On the unit sphere put u=cos(theta),

    V_phi=3u(1-u^2), V_theta=0,
    S_AB=D_A V_B+D_B V_A.

V is the axial l=2 vector harmonic associated with P2(u). Predicted identities:
div V=0, tr S=0, div S=-4V, |V|^2=9u^2(1-u^2),
S_theta_phi=-3 sin(theta)^3, |S|^2=18(1-u^2)^2.
Both V and S change sign under antipodal pullback.

Take the background TT tensor

    A0_ss=2C, A0_AB=-C Omega_AB, C=1/2, Lambda=1/4.

It has |A0|^2=6C^2=3/2 and the exact solution psi=1.
The scalar linearization at this solution is Delta_bar-5/4, hence invertible
on the compact domain. This predicts a local continuous family by the
implicit-function theorem; it does not prove a numerical finite-amplitude
existence interval or stability under time evolution.

Start with an explicitly NOT transverse seed

    T_sA=epsilon p(s)V_A, T_AB=T_ss=0, p(s)=sin(k s).

Solve for W_A=epsilon w(s)V_A, W_s=0 using the conformal Killing operator
(LW)_ij=D_i W_j+D_j W_i-(2/3)gbar_ij div W:

    w''-4w=-p',
    Abar=A0+T+LW,
    Abar_sA=epsilon [p+w'] V_A, Abar_AB=A0_AB+epsilon w S_AB.

The primary twisted case has k=1/2, so p and w are antiperiodic over 2pi;
their product with the odd angular tensors descends smoothly. The primary
analytic reference is w=k cos(ks)/(k^2+4), p+w'=4 sin(ks)/(k^2+4).
This reference is checked against an independently assembled numerical
one-dimensional solve. The nonzero seed divergence and its cancellation
must be measured. Setting K=0 or merely inserting a known TT tensor cannot
pass the momentum-completion gate.

Orthogonality of the axial perturbation to A0 predicts

    |Abar|^2=6C^2+2 epsilon^2(p+w')^2 |V|^2
                       +epsilon^2 w^2 |S|^2.

It is even in u and 2pi-periodic in s. The positive branch connected to
psi=1 is solved on that scalar fundamental domain. Tensor seam checks must
use the pullback including the angular Jacobian, not scalar periodicity.

## Prospective numerical schedule

No random search. Primary epsilon values: 0, .02, .05, .1, .2.
Fourier(s) x Gauss-Legendre(u) grids: (24,12), (40,20), (64,32).
Newton iteration with positive line search, at most 30 iterations, residual
target 1e-11; solve the linearized system by a spectral preconditioned
iterative or direct method to 1e-12 relative accuracy. Start at psi=1 or
continue from the preceding amplitude; record which. No clipping a failed
solution or discarding a failed amplitude. Record actual residuals and
iterations, including failures.

Solve the momentum ODE independently on a 4pi periodic double cover with
32,64,128 nodes using second-order centered differences; compare to the
closed form and demand convergence. The exact spectral solution may then
be used in the Hamiltonian solve; disclose this distinction.

Validation points: s=(.23,.77,1.39,2.17,2.81,4.13),
theta=(.41,.83,1.21,1.87,2.39), phi=.37 (30 Cartesian combinations).
Evaluate independently reconstructed coordinate metric and extrinsic
curvature on these off-grid points. Obtain Christoffels, Ricci scalar and
physical covariant divergence using centered coordinate differences of
the metric/tensor with h=1e-3,5e-4,2.5e-4. This is independent of the
elliptic residual formulas; record refinement, do not report algebraic
reuse of the PDE as an independent curvature check.

Evaluate the finest scalar interpolant on an independent (96,48) grid for
the PDE residual. Compare successive numerical solutions on the common
(64,32) grid. For all positive amplitudes record area(s), minimum psi,
max |K|, and the null expansions at the s=0 section. Check H=0 there by
reflection and numerical derivative. Test section areas at s=0,.05,.1,pi;
do not substitute coordinate radius for invariant area.

Controls at the medium grid, epsilon=.1 unless specified:

* Zero-wave control: epsilon=0, psi=1, full nonzero A0 retained.
* Amplitude reversal: epsilon=-.1; psi unchanged and axial K reversed.
* Full time reversal: C=-.5 and epsilon=-.1, same Lambda; psi unchanged,
  full K and both future null expansions reversed.
* Untwisted S2 x S1 control: k=1, ordinary periodic tensor seam. This
  tests whether solving the constraints distinguishes the chosen twist;
  success in both sectors is not evidence selecting antipodal topology.
* Bad tensor seam: use k=1 with the antipodal identification. Require a
  nonzero seam mismatch even if its local PDE residual is small.
* Uncorrected seed: w=0. Require nonzero momentum residual.
* Damaged completion: multiply the solved w by 1.1. Require nonzero
  momentum residual.
* Wrong conformal tensor weight: K=psi^+2 Abar instead of psi^-2 Abar.
  Independent physical constraint checks must distinguish it at epsilon=.2.

## Frozen gates and verdict dependencies

1. Symbolic angular/divergence/norm identities: exact zero residuals.
2. Momentum finite-difference solution: last max error <1e-3 and successive
   error ratios in [3.8,4.2]; analytic residual <1e-12; uncorrected and
   damaged residuals >1e-4.
3. Smooth correct tensor descent: max mismatch <1e-10 at the off-grid
   angular points; wrong-seam mismatch >1e-3. Check tensor values and first
   s derivatives. Scalar antipodal parity mismatch <1e-10.
4. Nonlinear solutions: all frozen grids/amplitudes converge, psi>0,
   on-grid residual <1e-10; finest off-grid residual <1e-8; successive
   solutions' finest max difference <1e-7. Failure remains failure.
5. Independent physical constraints at epsilon=.2: normalized Hamiltonian
   and momentum residuals <1e-5 at the finest coordinate step, with last
   error ratio in [2.5,5.5] whenever preceding error exceeds 1e-8; wrong
   conformal weight must exceed the correct error by at least 10x.
   Normalize each residual by max(1,sum of norms of its constituent terms).
6. Minimal-section diagnostic: max |H|<1e-8 at s=0; nonzero wave produces
   positive area differences Area(.05)-Area(0), Area(.1)-Area(0),
   Area(pi)-Area(0). All quantities are measured; a failure withdraws the
   neck interpretation without erasing a valid compact constraint solve.
7. Sign and topology controls: scalar agreement under amplitude/full time
   reversal <1e-9; untwisted local constraints meet the same on-grid target.
8. Evidence integrity: missing cases, nonfinite values, corrupted solutions,
   removed kinetic data, missing/unknown gate names and stale output cannot
   yield an affirmative combined verdict. Recompute gates from raw evidence.

`COMPACT_TWISTED_VACUUM_CONSTRAINT_DATA` needs gates 1-5,7,8.
`MINIMAL_SECTION_WITH_NONTRIVIAL_MOMENTUM_COMPLETION` additionally needs 6.
Keep unsupported claims explicitly false: round-S3 mouth embedding,
four-scalar interface, traversability, moving-mouth response, radiative
momentum transfer, crossing events, discrete action and quantum statistics.
An unexpected failure requires a new prospective addendum before additional
scientific measurements; do not relax gates or relabel failed cases.

## Context and deliverables

Inherited source inventory: `waves/physical_throat.py`, `waves/areal.py`,
`waves/nonlinear_supported_tt.py` and their corresponding documents. The
earlier neck is time-symmetric; #300's nontrivial completion is homogeneous
and has no mouth. Open #302 separately localizes the absent scalar interface.

Standard-method references (not claims of novelty):

* Chrusciel and Gicquaud, *Bifurcating solutions of the Lichnerowicz equation*,
  https://arxiv.org/abs/1506.00101 (compact S1 x S2, positive Lambda).
* Isenberg, Mazzeo and Pollack, *Gluing and wormholes for the Einstein
  constraint equations*, https://arxiv.org/abs/gr-qc/0109045 (constraint
  gluing; this experiment does not claim to implement their theorem).
* Corvino, *Constructing initial data for the Einstein equations*,
  https://www.univie.ac.at/AGESI_2017/school/corvino_ESI3.pdf (constraints).

Deliver a module, reproducible probe/raw arrays, tests and a result document.
Preserve this public freeze unchanged and report its commit. The scientific
contribution sought is a repository implementation and checked initial-data
bridge, not a new general existence theorem or quantization derivation.
