# Prospective constrained R3 return-map and leading-twist test

Date: 2026-09-29. Baseline: `c0c5f59bddbba7f5072265df05329085aa577e6f`
(PR #315). Publish this specification before implementing or measuring the
new nonlinear coefficient. Prior R3 runs and their UNRESOLVED verdict are
unchanged. The analytic reduction below is disclosed design work.

## Question and scope

Compute the leading action-frequency coefficient nu of the local LRS n=2
centre dynamics around the breathing ESU:

    rho(I) = rho0 + nu I + higher-order normal-form terms.

I is a local canonical action per unit S3 coordinate volume, with its
normalization fixed below. It is not absorbed action, a quantum, or a
spatial receiver output. The only possible physical conclusions here are
SHIFT_TOWARD_TARGET, SHIFT_AWAY_FROM_TARGET, or UNRESOLVED for the local
shift toward rho=3/2. No sign establishes a crossing, full-state closure,
global centre manifold, or physical selection of histories.

## Exact constraint and clock reduction

Use the unchanged `nonlinear_supported_tt` equations with a=kappa=1,
q=(q,0,0,0), beta=x b0, b0=diag(1,1,-2)/sqrt(6), and v=x'.
Write u=A', k=1/sqrt(6),

    T=tr(M^-1)=2 exp(-2kx)+exp(4kx),
    r=8 exp(-2kx)-2 exp(-8kx), ell=v^2,
    q=R cos(phi), q'=-2R sin(phi), R>0,
    R^2=[3u^2+A^2(r-v^2)/2-3A^4/2]
        /[2 sin(phi)^2+cos(phi)^2(T/2+(r-v^2)/12)],
    Q=R^2 cos(phi)^2, H=A^2-Q/6,
    H'=2Au+(2/3)R^2 sin(phi)cos(phi),
    Omega^2=T+(r+v^2)/6,
    phi'=2 sin(phi)^2+(Omega^2/2)cos(phi)^2.

The scalar amplitude is solved from the Hamiltonian constraint. Integrate
the exact inherited A,u,x,v equations divided by phi', plus eta'=1/phi',
from phi=pi/2 to 5pi/2. This takes the negative-velocity q=0 section to
its next negative-velocity crossing. It includes the state-dependent
return time, without event-time finite differencing or trajectory resets.
R^2,H,A,phi' must stay positive. A direct 29-state conformal-time solve
with actual section events is the independent formulation check.

At q=0 the canonical one-form per unit volume is

    Theta = -6u dA + A^2 v dx,
    z=(A-1,p_A,x,p_x), p_A=-6u, p_x=A^2 v.

The resulting four-dimensional return map P is symplectic on its domain.
No kicking, constraint projection during integration, damping, or long-time
Birkhoff estimator is used.

## Derivatives, centre graph and normal form

Integrate multivariate Taylor/variational jets in the four initial section
coordinates through total degree three, using ordinary monomial coefficients
(derivative divided by multi-index factorial). Arithmetic is truncated only
in polynomial degree. Record every coefficient and the return-time jet.
Use primary DOP853 rtol=2e-12, atol=2e-14, max_step=.025 in scalar phase;
repeat with RK45 rtol=2e-12, atol=2e-14, max_step=.0125.

The linear map splits into the hyperbolic (A,p_A) and elliptic (x,p_x)
blocks. Canonically normalize the elliptic block using its positive
invariant quadratic form, with transformation determinant +1. Fix rho0's
branch from the already checked linear phase orientation/winding, not from
new nonlinear data; expected rho0 approximately 1.484666408416.

Solve polynomial centre-graph invariance through degree three. Hyperbolic
coordinates are functions h2+h3 of the two elliptic coordinates; there is
no assumption that the hyperbolic variables can simply be set to zero.
Include the induced symplectic area form on the graph and its cubic
Darboux correction. Remove quadratic terms by a near-identity homological
solve. For zeta=(Q-iP)/sqrt(2), I=|zeta|^2, extract the resonant coefficient:

    zeta_new = lambda zeta + a21 zeta^2 conjugate(zeta) + ...,
    nu = Im(a21/lambda)/(2pi), lambda=exp(2pi i rho0).

The real part of a21/lambda must vanish numerically. This is a finite-order
local normal form, not proof of a convergent all-order normal form or of
an invariant torus at resonance. Report conditioning of all homological
solves. No fit to the earlier R3 data or fitted action unit is permitted.

## Independent finite-data validation

Use eight equally spaced angles and radii .008,.004,.002 in the canonically
normalized linear centre coordinates, embedded with h2+h3. These are
validation points, not an amplitude search or new frequency estimator.
For each, compare the cubic map with an unmodified 29-state DOP853 solve
(rtol=2e-12, atol=2e-14, max_step=.01) to the first negative q crossing
after the intervening positive crossing. Also compare against the exact
four-state phase formulation at the same tolerances. Save initial/returned
states and times, sampled constraints/chart minima, and all errors.

Use Euclidean canonical-coordinate norm. For each angle, the cubic-map
remainder and the centre-graph invariance remainder must fall with observed
order at least 3.5 on each radius halving unless the finer error is already
below the 1e-10 reporting floor. Require their finest-radius errors
below 1e-6. This is a numerical finite-order check, not a rigorous radius
of existence. Record all points, including failures.

## Gates and interpretation

All gates are required:

- G1: background return and linear off-block errors <1e-9; linear tensor
  trace agrees with #310 within 1e-9. The canonical linear symplectic
  residual is <1e-8. Linear rotation branch matches the disclosed value
  within 1e-9; no refitting that branch.
- G2: polynomial symplectic residual through degree two has max coefficient
  <1e-4 and normalized residual <1e-9 (normalizer 1 plus maximum coefficient
  of the sum of absolute products forming DP^T J DP).
- G3: centre-graph and normal-form homological residual max coefficient
  <1e-7; condition numbers <1e10; abs(Re(a21/lambda)) <1e-7.
- G4: primary/secondary nu difference <1e-5 max(1,abs(nu)); the sign is
  resolved only when abs(nu)>100 max(the difference,1e-9). Report full jet
  differences as well, without using a sign-only comparison.
- G5: full/phase return-state and return-time discrepancies <1e-9 at every
  validation point; sampled full-system constraint residual <1e-9; chart
  remains positive. Both remainder requirements above pass.

Any failure gives UNRESOLVED. Otherwise compare sign(nu) with
sign(3/2-rho0) to report SHIFT_TOWARD_TARGET or SHIFT_AWAY_FROM_TARGET.
Near-zero nu is UNRESOLVED, not no nonlinearity. If toward, a separate
prospective continuation/closure experiment is needed; if away, stop this
local search without excluding later turns, other sectors, or R2.

## Outputs and controls

Publish immutable freeze provenance, source hashes, environment, raw map
jets, intermediate reductions, full validation states and a replayable
report. Include synthetic canonical maps with known positive and negative
twist, a zero-twist control, finite-difference checks of jet arithmetic,
and tests rejecting missing/nonfinite/altered evidence and changed verdicts.
Archive numerical failures without loosening gates. Corrections require
dated notes; additional measurements require a prospective specification.
