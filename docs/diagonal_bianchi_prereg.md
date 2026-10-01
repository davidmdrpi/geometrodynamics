# Prospective diagonal Bianchi IX continuation and two-return closure test

Date: 2026-09-29 UTC. Parent: `590cf031a6fdd40df1d0179fb99a2821b7b015b6`
(merged #316 on `codex/r3-prereg-review`). Publish before implementing or
measuring new-amplitude diagonal circles or periodic orbits.

## Prior information and scope

#317 reports positive leading twist for W proportional to diag(1,omega,omega^2).
A review-time diagonal scalar-phase cubic calculation reproduced
nu=+0.847971923049527 using #317's normal-form reduction. That calculation
and #316's negative LRS coefficient are disclosed design inputs. No new
finite-amplitude diagonal circle or closure measurement has been made.
The leading estimate I~0.018 motivates the amplitude window below; it is
not evidence of a crossing. Earlier freezes, data and labels stay unchanged.

Test the exact diagonal homogeneous Einstein–quartet subsystem, with
a=kappa=1, q=(q,0,0,0), beta=x E0+y E1,
E0=diag(1,1,-2)/sqrt(6), E1=diag(1,-1,0)/sqrt(2), and M=exp(2 beta).
Diagonal M,L have zero gravitational angular momentum; the matter current
also vanishes, so all momentum constraints hold without compensators.
This is neither an inhomogeneous test nor a test of selected action.

## Exact map

Use the scalar-phase clock q=R cos(phi), q'=-2R sin(phi), R>0.
Write u=A', v=x', w=y', T=sum exp(-2 beta_i),
r=2 sum(2 exp(-2 beta_i)-exp(4 beta_i)), ell=v^2+w^2. Then

    R^2=[3u^2+A^2(r-ell)/2-3A^4/2]
        /[2 sin(phi)^2+cos(phi)^2(T/2+(r-ell)/12)],
    H=A^2-R^2 cos(phi)^2/6,
    H'=2Au+(2/3)R^2 sin(phi)cos(phi),
    phi'=2 sin(phi)^2+[T+(r+ell)/6]cos(phi)^2/2.

Integrate the inherited A,u,x,v,y,w equations and eta'=1, divided by phi',
from pi/2 to 5pi/2. The shape acceleration in direction Ek is

    -H'/H velocity_k + sum_i Ek_ii [(-4+q^2/H) exp(-2 beta_i)-4 exp(4 beta_i)].

Section coordinates z=(A,p_A,x,p_x,y,p_y) are canonical,
p_A=-6u, (p_x,p_y)=A^2(v,w) at q=0. Initial q' is the negative Hamiltonian
root. No resets, kicks or evolving constraint projections.
Primary map: DOP853, rtol=2e-12, atol=2e-14, max_step=.025 in phase.
Batched integration and complex-step map Jacobians are allowed. Validate
the Jacobian by centred differences before any amplitude continuation.
The full-system check uses unchanged 29-state conformal_rhs with actual
section events: stop at the positive q crossing, then at the negative one.
Use DOP853 rtol=2e-12, atol=2e-14, max_step=.01; maximum total time 8.
Check sampled chart positivity and all absolute constraints at 257 states.

## Circular-family continuation

Solve P(K(theta))=K(theta+omega) on 63 Fourier nodes, fixing the first
x harmonic to a/2 real. Seed the first circle with the linear elliptic
eigenvector and quadrature y harmonic (w_y=i w_x). Use the already fixed
positive rotation branch omega0=2pi*.484666408416, followed continuously.
Start a=.004, multiply by 2^(1/4), and stop before exceeding .2.
Use Newton/Gauss–Newton, at most 20 iterations and 10 step halvings per
iteration, with target maximum residual 5e-11. Scale the preceding circle
to initialize the next. On failure retry once at the geometric half step;
if this fails, stop as CONTINUATION_FAILED, never a physical family endpoint.

Accept only finite circles with positive canonical action magnitude,
grid residual <=1e-9, off-grid residual <=1e-9 on 63 half-offset nodes,
Fourier tail at |k|>=29 <=1e-9, and amplitude/phase errors <=1e-10.
The signed configuration-space area of (x,y) must retain the seed's sign.
At four quarter-period points compare the phase map with the full-system
event map: canonical-state/time discrepancies <=1e-8, sampled absolute
constraints <=1e-9, and positive A,H,M. Keep every attempted circle and
failure, including trial states and diagnostic exceptions when available.

Action is |(1/2pi) integral(p_A dA+p_x dx+p_y dy)|, evaluated spectrally.
Do not use signed displacement as action. Report rho=1+omega/(2pi).
Every accepted point is a measurement; none proves a continuous family.

## Crossing and refinement

On the first pair of accepted circles bracketing pi, stop the outer ladder.
Re-solve both endpoints on 95 nodes (same solver settings), require
grid and half-offset residuals <=1e-9 and tail at |k|>=45 <=1e-9.
For each endpoint define numerical uncertainty
delta=max(|omega95-omega63|,1e-9). Require delta<=1e-7 and endpoint
distance from pi >10 delta, with opposite signs. Only then label
CROSSING_BRACKETED_NUMERICALLY. This is an empirical resolution check,
not an interval proof. Otherwise label CROSSING_UNRESOLVED.
If no bracket appears, report NO_CROSSING_ON_ACCEPTED_LADDER and the
termination reason; do not exclude a crossing between or beyond samples.

## Direct closure by multiple shooting

Attempt closure if a bracket is found OR an accepted circle comes within
.03 radians of pi. Use the accepted circle closest to pi. For eight
equally spaced seed angles, take two nodes at theta and theta+pi and solve

    P(z0)-z1=0, P(z1)-z0=0.

Use least-squares with the map Jacobian, xtol=ftol=gtol=1e-12, max_nfev=80.
No fitted action, fixed frequency, amplitude constraint or trajectory reset
is used in this solve. Require max canonical residual <=1e-9,
|P(z0)-z0|_2>1e-5, and each node's distance from the ESU point greater
than one quarter of the minimum seed-circle distance. This rejects collapse
to the background or a one-return fixed point.

Verify candidates independently with two full 29-state event returns, both
DOP853 and Radau (same tolerances/max_step). Require full-state closure
max norm <=1e-8, excluding elapsed eta, plus phase/full intermediate-node
agreement <=1e-8, sampled absolute constraints <=1e-9 and positive charts.
Only then label NONTRIVIAL_TWO_RETURN_ORBIT_VERIFIED_NUMERICALLY. Otherwise
CLOSURE_UNRESOLVED, or NOT_ATTEMPTED if the proximity condition was not met.
Record all seeds, failures, nodes, residuals, return times and full histories.
No uniqueness, stability, physical selection or quantization follows.

## Evidence and controls

Archive immutable source hashes, environment, freeze commit, all raw K
arrays and map images, off-grid images, full validation histories, solver
traces, refinement records and shooting attempts. Refuse overwrites.
Replay must authenticate archives and reconstruct actions, acceptance,
crossing and closure decisions from the saved evidence, without trusting
stored ok flags. Missing/nonfinite/altered records must be rejected.
Add toy circular twist/period-two controls, exact/full vector-field and
map/Jacobian comparisons, and mutation tests. Preserve negative outcomes;
any implementation correction requires a dated disclosure. No changes to
the registered ladder, thresholds or physical decision rules after results.
