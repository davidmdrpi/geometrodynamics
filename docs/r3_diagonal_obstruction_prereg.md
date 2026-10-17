# High-precision diagonal obstruction: specification before production

11 October 2026. Follow-up to #326, frozen after off-family and already-known
LRS/synthetic controls, before any new high-precision diagonal-loop obstruction.
The loop, its low-precision null, and the order-twelve selection rule are known.
The new signal, its amplitude, zero pattern and isolated roots are held out.
Historical archives and producers remain unchanged.

## Question and distinct outcomes

Does the #319 labelled circular component resolve into isolated periodic points
when the two-node obstruction is measured below its old 1e-11 resolution?
The selection rule permits k=12, but does not imply a nonzero coefficient.
A leading twelfth term suggests 24 simple zeros on the labelled cover. This
pattern prediction can fail independently of detecting nonzero breaking.

Use the historical coordinate ordering (A,p_A,x1,p1,x2,p2), dimensionless maximum
norm, the archived periodic splines, Euclidean unit tangent phase slice, and
unit symplectic-dual unfolding direction from r3_breaking. Preserve those
normalizations; do not interpret lambda as energy. The two node seed splines
have separate chord parameters, as in the old producer. They are only Newton
seeds, not an assertion of a common invariant-circle angle. The solved phase
slice is anchored to node zero only.

At 192 fixed phases phi=2*pi*j/192 solve the 13-dimensional augmented system
P(z0)-z1=0, P(z1)-z0-lambda*g=0, t dot (z0-c)=0. No parameters, background
geometry or amplitude are retuned. There are no corrective kicks in the
physical calculation. This is periodic shooting, not long-time stability or
autonomous action selection. The physical ESU hyperbolic direction remains.

## Independent numerical formulations

Primary: eight-state canonical-momentum Taylor integration, 224 bits (over
64 decimal digits), order 48, local coefficient scale 1e-56. Clock crossing
uses Newton on its Taylor polynomial. Check: eight-state shear-velocity Taylor
integration, 336 bits (over 96 decimal digits), order 64, scale 1e-80; crossing
uses bracketed bisection. The implementations share scalar convolution and
exponential recurrence machinery, but integrate different state equations and
locate crossings differently. Both use a direct full return, not H squared.
Off-family agreement is tested against the original double-precision matrix
map. The method is numerical, not certified interval integration.

The primary augmented solve requires returned-node residual <=1e-45; the check
requires <=1e-65. A double-precision Jacobian supplies chord corrections, while
states, residuals and updates retain arbitrary precision. Recompute residuals
after the last update even on iteration exhaustion. Maximum 24 iterations;
failed convergence fails the gate. Record constraint errors throughout solves.

For each phase let E_j = |lambda_primary-lambda_check| +
2*max(||row_lambda(J^-1)||_1)*(residual_primary+residual_check).
This is a conservative *numerical error estimate*, not a rigorous nonlinear
bound. Set E=max(E_j) and F=max(E,1e-30). The accuracy gate requires E<=1e-30,
so accepted production uses a fixed resolution floor F=1e-30 and a detection
threshold 100F=1e-28. This avoids labeling uncalibrated extra arithmetic digits
as physical evidence. Also require all solves to pass, condition(J)<=1e7,
seed displacement <=.01, all high-precision map constraints <=1e-48, and all
positive/null controls to pass. Failing any gate yields NUMERICALLY_UNRESOLVED.

With the gates passing:
- max|lambda|>100F gives RESOLVED_OBSTRUCTION, pending root verification.
- max|lambda|<=100F gives ADDITIONAL_SUPPRESSION_CANDIDATE. It is a sampled
  upper limit, not exact integrability or a calibrated anomaly.
- Successful distinct, numerically isolated roots at every reliable sign-change
  bracket upgrade the first outcome to RESOLVED_BREAKING.
- Failed root verification retains RESOLVED_OBSTRUCTION_ROOTS_UNRESOLVED.

Record the raw Fourier spectrum descriptively. Nonmultiples of twelve are not
a symmetry falsifier in this nonuniform phase and normalization. Separately
report whether 24 simple roots with no unresolved grid phases support the
leading-twelfth pattern. A different resolved count rejects that pattern;
it does not invalidate the exact clock/spatial symmetries.

The prospective new amplitude window is 1e-28 to 1e-11, or effective exponent
-ln(max|lambda|)/12 approximately 2.11 to 5.37. The earlier limit is sampled.
The LRS family's measured exponents do not calibrate this circular family's
coefficient. A new null would warrant investigation, not prove a hidden integral.

## Controls and root verification

Before production, compare three phase-distinct LRS 3/7 slices (j=0,17,41)
with the pinned 60-phase 48-digit archive at e7afdeb46e274799bd1c0ca749e85e52f70dc50c.
Use the *same* saved phase slices and augmented equations. Require lambda
agreement <=1e-32 with the archive, and <=1e-40 between the new formulations.
Archive hashes are checked before use. Its high-precision solutions are known
controls, not held-out evidence for the diagonal result.

A synthetic symplectic kick-drift map has an exact lambda=-epsilon*sin(12*phi)
under its phase slice. At epsilon=1e-24 recover it to <=1e-55 and count 24
sign changes on 192 offset phases. Also pass the epsilon=0 null control.
Check the numerical isolation diagnostic on an exact synthetic root with an
ill-conditioned inverse (~1e24). The synthetic map is a solver test only.

For each cyclic adjacent phase pair with opposite signs and both magnitudes
>10F, interpolate the saved slice anchors/tangents, normalizing the tangent at
high precision. Use bracketed secants on that slice parameter, with at most
48 trials. At each trial redo the double seed solve and high-precision shooting.
Require slice residual <=1e-56 and |lambda|<1e-51. Do not regard merely small
closure residual as proof of isolation.

At the resulting two-node solution, form the unaugmented 12x12 closure
Jacobian by centered high-precision differences: primary step 1e-22, independent
check step 1e-24. Invert the latter at 336 bits with partial pivoting. Require:
- inverse-norm times formulation/step Jacobian difference <.01;
- independent closure residual <1e-45;
- inferred Newton correction maximum norm <1e-12;
- inverse multiplication residual <1e-40.

These checks establish numerical nonsingularity at the tested resolution;
they are not an interval existence theorem. Distinct node-zero roots must be
separated by >1e-6. Report root equivalence under P and C, and separately H
and C, with a matching tolerance 1e-8. All labelled representatives and all
matching errors remain visible. Do not claim 24 physically distinct objects,
classify elliptic/hyperbolic stability, or infer an action spectrum.

## Execution and provenance

Freeze sources, protocol, tests, and control results publicly before scanning.
Eight numerical workers; no delegated agents. Atomically save each primary
solve before its check, each completed phase, and each completed root. On
resumption verify source/input fingerprints and reuse completed jobs. Never
silently replace a failed phase, refine a grid, alter a threshold or change
formulation after reading results: any repair must be disclosed separately.

The producer prints only progress indices until the complete scan is scored.
Run root verification after that score, without changing its rules. Publish
both numerical stages, source-bound archives, authenticated replay and the
result, including failure labels. If execution is interrupted, preserve the
per-job checkpoints and resume the frozen calculation.
