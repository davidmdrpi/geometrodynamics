# Prospective nonlinear test with an explicit intervention budget

Date: 2026-10-05 UTC. Parent: #320 at `80d894a9304b9c1cbe4db51b955865a0083d6f4b`.
Publish before any new nonlinear trajectory. Freeze criteria, cases and
controller before observing outcomes; bind the specification and producer
sources in the evidence. This PR is stacked on #320 while it remains open.

## Hypothesis, scope and fine tuning

Test whether four specified diagonal action-detuned initial conditions stay
within a prescribed tube around #319's family until a twist-based horizon,
using only the corrective-reset budget below. This is a controlled nonlinear
robustness test, not autonomous stability, a physical selection mechanism,
or persistence of exact resonant closure. No general off-diagonal or
inhomogeneous nonlinear claim will be made.

The Einstein-static multiplier is approximately 7242 per two returns.
Initial data are deliberately fine-tuned: the scale coordinates A,p_A are
adjusted to remove the *linear* stable/unstable components relative to the
inherited family. This does not put them on a proven nonlinear centre
manifold. For N steps, linear autonomous shadowing would require an unstable
initial error of order tube_radius * 7242^(-N). Report this precision estimate
in decimal logarithms rather than implying the tuning is physically generic.

Corrective resets are external interventions, not terms derived from GR.
They alter A,p_A and the Hamiltonian-completed clock momentum. Every reset,
including the proposed reset that first exceeds the budget, must be recorded.
An over-budget reset is never applied. Re-centering or resetting cannot be
hidden in the phase fit, integrator or output processing.

## Coordinate norm and reference family

Restrict the unchanged full 29-state conformal_rhs to diagonal tensor data
and quartet q=(q0,0,0,0). At descending q0=0 sections use

    w=(A,p_A,x,p_x,y,p_y,q0'), p_A=-6A', (p_x,p_y)=A^2(x',y').

Use the unweighted Euclidean norm of these seven dimensionless coordinates.
It mixes coordinate and momentum units fixed by kappa=a=1 and is not an
invariant physical distance. The reset cost is ||w_after-w_before||_2.
The clock momentum is included, so constraint completion is not free.

Use all ordered z0 points of #319 stage F, periodic cubic interpolation in
normalized chord arclength theta in [0,2pi). Complete each reference section
point with the negative Hamiltonian root for q0'. Distance to the family is
min_theta ||w-K7(theta)||_2. Locate the minimum from all local minima of a
fixed 256-point phase grid followed by bounded scalar minimization on each
neighbouring interval (xatol=1e-12). Save the fitted phase and residual.
Phase fitting changes no state. Record actual wrapped and unwrapped phase;
only the distance calculation may quotient phase motion away.

## Fixed initial conditions and tuning algorithm

Use delta_I = -0.001, -0.0005, +0.0005, +0.001, starting at #319 node 0.
Scale tensor coordinates and momenta on the entire reference loop by s.
For every loop point adjust A,p_A so the displacement annihilates the two
hyperbolic left coordinates defined below. Choose s by a scalar root solve
on [0.9,1.1] so the canonical loop integral differs from the unperturbed
loop by delta_I (absolute residual <=1e-10). Use #319's loop_action routine.
The scaled, tuned loop is NOT asserted to be invariant. Its action change
is a precisely defined input proxy; it is not an independently measured
conserved action of each perturbed trajectory.

For the 12 inherited stage-S matrices, form the real spectral projector
onto eigenvalues of modulus <.1 or >10. Interpolate its entries periodically
in the same phase coordinate. At any fitted phase use the two leading right
singular vectors as rows L; solve L[:,0:2] delta_(A,p_A) = -L[:,2:6] delta_tensor.
Require that 2-by-2 matrix's condition number <=100. This supplies the same
fixed linear tuning/controller rule at every phase; no new local monodromy
fits or adaptive controller gains are allowed. Record the initial tuning
vector and its norm separately from subsequent reset costs.

Let d0 be the tuned initial state's distance to the reference family. If
no root, valid clock completion or nonzero d0>1e-6 exists, report INVALID_SETUP
and do not silently drop the case.

## Twist-based horizon

Prior information: quadratic interpolation of I(omega) through #319's three
archived circles gives tau = d omega/dI at omega=pi approximately 4.90858.
Recompute it from those authenticated numbers. The predicted drift per
two-return step is 2*tau*delta_I. Set

    N(delta_I) = ceil(1 / (2*abs(tau*delta_I))).

This is 102 steps for |delta_I|=.001 and 204 for |delta_I|=.0005, targeting
one radian of detuning drift. The empirical twist is approximate, and the
trial loop is not proved invariant. Therefore a positive result additionally
requires measured absolute unwrapped phase drift >=0.5 radians by N, with
each fitted-phase increment <pi/2. If it does not, report
DRIFT_HORIZON_NOT_RESOLVED, not success. A safety cap of 256 steps must never
be reported as completion if the formula demands more. No horizon shortening
or extension after observing the nonlinear response.

## Evolution and corrective-reset budget

For each case run both an unforced control and a reset-enabled arm. Evolve
full states consecutively through two actual descending q0 sections per step.
The unforced arm never resets any component. At each completed step:

1. Save the pre-reset state, time, constraints, fitted phase and distance.
2. Fail immediately if pre-reset distance >2*d0 (the tube radius).
3. For the controlled arm, compute the A,p_A target with the fixed L rule
   relative to the nearest family point, retaining the current tensor data.
   Reconstruct the full state and solve the Hamiltonian clock root. Preserve
   elapsed time. Save the proposed full state and seven-coordinate cost.
4. Apply only if BOTH per-step cost <=0.02*d0 and cumulative cost <=0.25*d0.
   Exceeding either is CONTROL_BUDGET_EXCEEDED; record but do not apply it.
5. Check the post-reset distance <=2*d0. Post-reset tube escape is also failure.

The cumulative budget is a sum of magnitudes, never a signed cancellation.
The initial fine-tuning is disclosed separately; it is not evidence that the
unforced system selected a special initial condition. The budgets and tube
are operational choices for this experiment, not constants derived from GR.

## Pre-stated failure and numerical validity

Run every case/arm with DOP853 at rtol=2e-12, atol=2e-14, max_step=.025 and
independently with Radau at the same settings. Two returns use terminal
section events after a short initial departure from the starting section.
Check positive chart and all four absolute constraint residuals <=1e-8 on
33 uniform samples per return, as well as all event and reset states.
Save sampled full histories and endpoints.

A controlled case fails the hypothesis on TUBE_ESCAPE or
CONTROL_BUDGET_EXCEEDED, even at the first step. A chart exit also prevents
a positive result. Solver failure, violated constraint accuracy, phase-fit
failure, or integrator disagreement gives NUMERICALLY_UNRESOLVED. Require
matching terminal reason and step between integrators, with matching
recorded section states to 1e-4*d0 and distance/cost discrepancies to 1e-4*d0.
Both integrators must independently reach the same physical/budget failure
for it to count as a resolved negative. Diagnostics must be finite.

A controlled case passes only if it completes N steps, stays in the tube
both before and after resets, stays within both budgets, has valid constraints
and chart throughout sampled histories, resolves >=0.5 radian phase drift,
and agrees between integrators. Aggregate labels:

- CONTROLLED_NONLINEAR_BOUND_FAILED if at least one controlled case has an
  integrator-confirmed tube/budget failure. Report any other unresolved cases.
- CONTROLLED_NONLINEAR_BOUND_SUPPORTED_ON_TESTED_HORIZON only if ALL four
  controlled cases pass.
- NONLINEAR_STUDY_UNRESOLVED otherwise, including insufficient phase drift.

Always report unforced outcomes and the total tuning/reset intervention.
An early failure may save computation; it cannot turn into a success label.
No outcome establishes autonomous nonlinear stability or action selection.

## Evidence and replay

Authenticate #319/#320 inputs, bind source hashes and this freeze, and retain
all successful and failed cases in readable archives. Replay recomputes the
trial-loop action, horizons, nearest-phase distances, reset proposals/costs,
budgets, constraints, per-case and aggregate decisions from raw states; it
must not trust saved pass flags. Pin archive hashes. Tests must include a
budget-exhausting kick, a tube escape, insufficient drift, a cumulative
budget overflow, and corrupted/missing/reordered evidence. Preserve previous
freezes and measurements. Document implementation fixes without changing
criteria or erasing completed attempts.
