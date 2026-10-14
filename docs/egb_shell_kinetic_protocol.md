# Action-level admission protocol for the EGB shell

Written after the analytic completion/constraint calculation, before numerical
checks. This is not blind preregistration. The original #330 producer and
archives remain frozen. Use R=1, a=.3, alpha=-.06, eta=.1, C=.76. No kicks,
stabilization, new radius fitting or post-result change of these parameters.

The target is the physical canonical norm, not the sign of the coefficient D
in a pressure equation. A physical norm may be assigned only to perturbations
that satisfy the linearized bulk and junction constraints of a specified action.

Specify a minimally coupled, shift-symmetric scalar chi throughout the retained
bulk, continuous at the shell, with surface action P(Y). Define X=-grad(chi)^2/2,
Y=-D(chi)^2/2 and chi=q t in the static background. Choose the canonical bulk
P_b=X-V, q^2=3C/R^2 and V=(9/R^2+12 alpha/R^4)/2. Choose the power-law-plus-
constant surface action to reproduce the existing affine EOS exactly, for
0<eta<=1. It is an added material postulate, not derived BAM matter.

Also test P_b,lambda=X-V+lambda*(X-X0)^2/(2X0), with lambda=-1,0,1. This family
has identical background stress and frozen-bulk flux but different quadratic
bulk kinetic coefficients. Lambda=-1 is an explicit rejection control, not an
accepted theory. The lambda=0 completion is the nominated candidate.

Numerical verification gates, fixed before running:

1. Match static bulk stress, surface stress and affine EOS to relative 1e-10.
   Verify the normal scalar current junction at equilibrium.
2. Independently differentiate the actual matter Lagrangians in temporal and
   spatial field velocities. Recover the bare quadratic kinetic/gradient
   coefficients to normalized 1e-5 at two step sizes, and recover c_surface^2=eta.
   The lambda=-1 control must have negative bulk time-kinetic coefficient.
3. Check the geometric density slope against finite differences to normalized
   1e-7, and the linear scalar-current obstruction/closed total shift-charge
   derivative to normalized 1e-7. Include flux signs and both mouths once each.
4. Test turning-point displacements epsilon=1e-3,1e-4,1e-5,1e-6 with bulk metric
   and chi=q t frozen. Divide the density junction residual by epsilon, so a
   mere small-amplitude limit cannot make a first-order constraint failure pass.
   A constraint-satisfying tangent requires the dimensionless residual divided
   by epsilon to approach zero (threshold 1e-6 at the last step). Otherwise the
   old radial oscillator is not an admissible mode of this completion.
5. Evaluate the same mismatch and scalar-current residual on the original
   archived trajectory, using its prescribed acceleration. This is a residual
   audit, not an evolution under the new action; do not relabel those arrays.

If numerical identities fail, label NUMERICALLY_UNRESOLVED. If they pass but
the constrained tangent fails, label STATIC_ACTION_MATCH_FROZEN_MODE_REJECTED
and leave physical_radial_kinetic_norm null. That outcome does not reject
coupled solutions of the new action, or #330's stated restricted radial window.

Positive bare matter Hessians do not pass the physical canonical-norm gate.
That requires the bulk metric/scalar perturbations, lapse and shift constraints,
linearized junction and scalar current conditions, and the physical quadratic
form with its gravitational boundary contributions. Do not assign a sign by
fitting an oscillator mass, eliminating support fields while keeping their
stress prescribed, or treating normal-to-boundary canonical momentum as the
time-evolution kinetic norm. Smooth regularity and nonspherical stability remain
separate gates. If the constrained tangent passes, the physical norm remains
UNCOMPUTED until this reduction is actually performed.
