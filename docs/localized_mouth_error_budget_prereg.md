# Prospective independent physical error budget

Date: 2026-09-22. Parent evidence: PR #306, commit
`491fce5d47a2b4f634c0eb76476c598ef1be0b36`.

Publish this specification before running the new measurements. This is a
separate extension, not a rescore that turns the historical 6/8 or 7/8 into
a pass. The previous convergence test remains failed. No new BVP solve,
interpolant, profile, branch search, or altered source is permitted.

## Question

Does the finite-difference curvature approach a small nonzero residual of
the saved interpolant, rather than zero? A ratio of total residual norms
cannot separate those errors. Nor can a fitted constant from maxima at
different points. Measure the signed residual of the saved fields directly
and then test the differentiation error against that independent reference.

The reference may differentiate stored polynomial coefficients and elementary
coordinate functions. It must NOT evaluate the Hamiltonian ODE, substitute
an equation for psi'', or use the conformal scalar-curvature identity.
The reconstruction itself remains ODE-informed; this is independent
verification of its physical constraints, not an independent BVP solution.

## Inputs and reference

1. Verify the original decompressed archive SHA256
   `e510ba44aeb7d8d8d25cfc28cfb54a6a402d45ba43d373d5f6b48b6bfdffd227`
   and stable reconstruction decompressed SHA256
   `b9fd7d2e652ccdae2012c636c7fc249fc8b5cca05cb872df5b32576f591cb3ae`.
   Replay both historical scorers and require respectively 6/8 and 7/8,
   with precisely their recorded failed gates. Inherit the seven passing
   refinement gates only if this evidence check passes.
2. Use the archived L=5.5, eta=.3 finest reconstructed solution and the
   same 20 physical points. Preserve the binary64 coordinates, polynomial
   coefficients and knots exactly when importing into arbitrary precision.
   Use exact rational model constants (f=7/8, Lambda=3/2, q=sqrt(3)/2).
3. Evaluate the diagonal Jordan metric, extrinsic curvature, quartet and
   normal momenta from those polynomials and elementary functions. Compute
   coordinate first/second derivatives by high-precision local differentiation
   of those functions. Assemble Christoffels, their derivatives, Ricci and
   covariant momentum divergence by coordinate tensor contractions. Evaluate
   scalar gradients and all matter contractions independently as well.
4. Repeat the reference at 60 and 80 decimal digits. Every returned geometry,
   source, signed constraint and normalization component must agree to
   1e-40 times max(1, absolute 80-digit component). Reject nonfinite values.
   Archive both references with at least 50 significant decimal digits.

## Fixed differentiation experiment

The saved representation is piecewise polynomial, not globally analytic.
For each point find the distance d to the nearest knot among the psi,
theta and tensor profiles. Reject a point on a knot. Set h0=min(.008,d/8),
then use exactly h0,h0/2,h0/4, at 80 digits. This geometric rule is fixed
before any residual is measured; no adaptive tuning is allowed.

Construct fourth-order centered first/second coordinate differences and
tensor-product mixed differences from function VALUES only. Recompute
the physical constraints using the same coordinate tensor contractions.
All stencil points must remain in their respective polynomial intervals.
Use the 80-digit reference H and M scales for all three errors at that
point. Record signed H and all three signed M components before norms.

At EVERY point require:

- Reference and finest finite-difference |H|/Hscale and norm(M)/Mscale <1e-5.
  These are the existing absolute tolerances, applied to both estimates.
- For H, E(h)=abs(H_FD(h)-H_reference)/Hscale. For M use the norm of the
  vector difference divided by Mscale. BOTH consecutive E ratios must be
  in [8,24], unless the preceding E is <=1e-40. The tiny-error exception
  does not waive absolute accuracy. All reference and difference values
  must be finite. Require at least one active ratio for H and for M.
- Finest differentiation error E <1e-8 for H and M. Subtracting a large
  residual floor cannot generate a pass because absolute residuals are
  checked separately.
- Wrong-f Hamiltonian residual (replace f only on the geometric side by 1),
  normalized with the reference scale, exceeds ten times the larger of
  the correctly weighted reference H and M residuals.

Also evaluate the original .008,.004,.002 fourth-order stencils at 80 digits
as DIAGNOSTICS, including signed pointwise differences from the reference.
These may cross knots; no local asymptotic convergence claim is made for
them. Retain the float64 historical results and their failure unchanged.

## Controls, evidence and verdict

Before data measurement, validate the coordinate engine against: the unit
round S3 metric (R=6); the product cylinder ds^2+dOmega^2 (R=2); and that
cylinder with K=diag(0,s,s sin(t)^2), whose trace is 2s, norm squared is
2s^2 and momentum constraint geometric divergence is (-2,0,0). Require
reference relative/scaled error <1e-40 and demonstrate fourth-order
finite-difference curvature convergence for the sphere at steps
.008,.004,.002, ratios [8,24]. These calibrate geometry, not the GR data.

Archive input hashes, public freeze commit, source hashes, dependency
versions, all signed measurements, errors and ratios. Rescoring must
recompute the measurement rows from the hash-verified inputs and compare
numeric evidence with scaled tolerance 1e-40. Missing/changed points,
steps, source records, K, polynomial coefficients or constraint samples
must invalidate evidence. Tests must reject a fabricated constant floor
that exceeds the absolute tolerance even if its differentiation converges.

The new physical gate is the conjunction of calibration, precision,
absolute accuracy, differentiation convergence, differentiation accuracy
and wrong-f controls. Combine it with the seven inherited gates only
after evidence replay succeeds. Publish separate extension verdicts for
FOUR_SCALAR_HANDLE_CONSTRAINT_DATA (gates 1-5,7,8) and
LOCALIZED_BULK_MOUTH_INITIAL_DATA (also gate 6). Report failure and stop if
any criterion fails; do not change resolution or thresholds to obtain a pass.

Passing certifies the registered numerical initial-data checks only. The
chosen sign bundle, matter quartet and preparation remain assumptions.
Traversability, evolved crossings, momentum-transfer events, discrete
action and quantum statistics remain unestablished.
