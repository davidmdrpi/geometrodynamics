# Specification: R3 beyond the local LRS coefficient

Date: 2026-09-29. Parent: `claude/geometrodynamics-qft-audit-vpktax` at `0100c74` (PR #317).
Publish this specification before any cubic-order computation for the five
polarisations and before any new ESU invariant circle. The measurement code
is committed with this document. Corrections are dated notes only.

## 1. What can and cannot be ruled out here

#316 established SHIFT_AWAY_FROM_TARGET prospectively, and #317 cross-checks
it. The claim is local: the linearly polarised (LRS) n=2 tensor at small
action. Four caveats remain. This extension settles what the exact
homogeneous system can settle, and labels the rest NOT_TESTED.

| caveat | treatment here |
|---|---|
| a turn at larger amplitude | **Part A**: continue the LRS invariant-circle family until it ends |
| other polarisations | **Part B**: leading twist for every polarisation of the homogeneous n=2 tensor (all 5 components) |
| other closure conditions | **Part D**: which rational rotations the family actually reaches |
| other sectors (inhomogeneous n >= 3 tensor, vector, scalar) | **Part C: NOT_TESTED**. See section 5 |

A negative result below rules out the R3 refocusing closure only inside the
stated scope. It is not a proof about the full theory.

## 2. Part B: all five homogeneous polarisations

**System.** Use the full homogeneous matrix system: the `conformal_rhs`
equations, with q = (q,0,0,0) and a general unimodular shape
M = exp(2 beta), where beta = sum_k x_k e_k over an orthonormal STF basis.
A pre-freeze check confirms that the equations match `conformal_rhs` to
2e-15. The section is q = 0 with q' < 0; the Hamiltonian constraint is
solved for q'. The section coordinates are
z = (A, p_A, x_1, p_1, ..., x_5, p_5), with p_k = H tr(beta' e_k).

**Canonical coordinates.** These coordinates are canonical at the fixed
point, and the leading twist needs nothing more. A non-canonical
symplectic-form correction linear in z integrates to zero over the
centrally symmetric leading-order circle. Higher corrections enter the
action at O(|w|^4). The LRS block agrees with #316's canonical coordinates
exactly.

**Momentum constraint.** Frame-rotating polarisations carry gravitational
angular momentum of order eps^2. The momentum constraint then requires a
matter current, supplied by quartet components dq of order eps^2. Those
components enter the (A, q, M, L) equations only through |q|^2, q·q' and
|q'|^2, which are O(eps^4). So they change the section map only at degree
>= 4, and the cubic normal form computed with uncompensated tensor data is
the physical one.

**Normal form.** SO(3) symmetry makes the five elliptic pairs 1:1-resonant.
The order-3 normal form (`r3_extension.multimode_normal_form`):
- diagonalise into u, s, w_1..w_5 and their conjugates;
- remove all quadratic terms, which are nonresonant;
- read the resonant cubic tensor T_d(w_a w_b wbar_c)/mu.

For polarisation w, the mean phase-advance shift per unit action is

    nu(w) = Im<w, G(w)> / (2π c |w|^4),   G_d(w) = sum T[d,a,b,c] w_a w_b wbar_c.

For relative equilibria this is exact. For any other motion, the time
average of nu lies between min nu and max nu over the sphere. Two
quantities are reported:
- max and min of nu over the unit sphere of C^5, found by 64 BFGS restarts;
- the Sym^2 eigenvalue bound. If it is negative, it certifies that nu < 0
  for every w.

The method reproduces a toy map with a known quartic to 1e-11. The toy
includes a relative-phase coupling term.

**Checks.** Failing any gives UNRESOLVED.

| id | requirement |
|---|---|
| B1 | fixed-point error <= 1e-10 |
| B2 | \|theta - theta0(#310)\| <= 1e-9; off-block linear entries <= 1e-9; the five tensor blocks agree to 1e-9 |
| B3 | dissipative (Hermitian) part of the quartic form <= 1e-6 × twist scale |
| B4 | SO(3) covariance: \|nu(D(R)w) - nu(w)\| <= 1e-7 × scale, for 5 random rotations × 5 random w |
| B5 | nu(LRS) equals #316's value, -0.9501180968942993, within 1e-8 |
| B6 | all real (linear) polarisations give equal nu within 1e-7 × scale (an SO(3) prediction, since tr b^4 = (tr b^2)^2/2) |
| B7 | Richardson (1024,2048) vs (2048,4096): max, min and LRS nu agree within 1e-7 × scale |
| B8 | every quadratic divisor >= .1 *relative* to the larger multiplier involved; \|mu^k - 1\| >= .1 for k = 1..4 |
| B9 | image-constraint and section-q jets <= 1e-8 × max coefficient |

B8 is written in relative form. The #317 absolute C6 threshold was
mis-specified, and this avoids repeating that error. The change is
prospective here; #317's archived label is not rescored.

**Labels.**
- **ALL_POLARISATIONS_SHIFT_AWAY:** max nu < 0, resolved as |max nu| >= 100 × the B7 difference.
- **SOME_POLARISATION_SHIFTS_TOWARD:** max nu > 0, resolved on the same criterion. The maximising polarisation is reported.
- **UNRESOLVED:** otherwise.

## 3. Part A: the LRS family at large amplitude

**Circles.** Use invariant circles of the #317 map (`esu_map`, full
system, DOP853 at 1e-12/1e-14) on a 63-point grid. Solve by Gauss–Newton
with a finite-difference Jacobian computed in 4 processes, to a tolerance
of 1e-11, with at most 15 iterations.

**Ladder.** a_k = .004 × 2^(k/4), up to a = 2. Each circle starts from the
previous one, rescaled.

**Acceptance.** A circle is accepted if its residual is <= 1e-10, its
Fourier tail (|k| >= 29) is <= 1e-10, and the constraint at every image is
<= 1e-10.

**Stopping.** On the first failure, retry once at the half step
(× 2^(1/8)). If the retry also fails, the family ends there, and the
reason is recorded.

**Labels.** Rotation is ordered by action.
- **CROSSING_IN_FAMILY:** some omega >= π.
- **TURN_IN_FAMILY:** omega increases by more than 1e-9 between consecutive circles.
- **NO_TURN_IN_FAMILY:** neither, with at least 8 accepted circles.
- **UNRESOLVED:** otherwise.

**Scope.** A family that ends does not exclude a turn beyond its end,
where invariant circles break up. That case is reported, not claimed.

## 4. Part D: other closure conditions

Along a twist family, every rational rotation p/q in the family's range
has closed histories (Poincaré–Birkhoff), at isolated actions for each p/q.
The union over all rationals is dense. Closure alone therefore selects no
discrete set; only a physically fixed rational would.

Part D lists:
- every p/q with q <= 12 inside the measured range of frac(rho);
- which of the natural low-order closures (0, 1/4, 1/3, 1/2, 2/3, 3/4, 1)
  fall inside it.

It is descriptive and assigns no verdict.

## 5. Part C: sectors not tested, and why

- **Inhomogeneous n >= 3 tensor modes.** A cubic frequency shift needs
  second-order GR perturbations on the breathing ESU with n-dependent
  harmonic couplings. The code here evolves only homogeneous
  (spatially uniform) geometries.
- **Scalar n = 2.** It is linearly hyperbolic, with #310 multiplier 1.2475.
  There is no elliptic centre, so the R3 closure is not a small-amplitude
  question there.
- **Vector modes.** Unresolved already at the linear level in #310.
- **Homogeneous quartet rotations.** Neutral: multiplier 1, the SO(4)
  symmetry directions. They have no twist to measure.

These stay open, and no conclusion here extends to them.

## 6. Pre-freeze disclosure

Before this freeze, only the following were run:
- toy-map tests of the multi-mode normal form;
- the equation, parametrisation-round-trip and SO(3)-representation tests;
- the order-1 five-polarisation map: all five blocks match #310 to 5e-12,
  off-block entries are exactly 0, and the multipliers are 85.02 and .01176;
- a reproduction of the already archived a = .004 circle (omega = 3.0451116429)
  by the parallel solver.

No cubic five-polarisation jet and no new-amplitude circle has been
computed. Crashes may be fixed under a dated note, without changing
parameters or rules.
