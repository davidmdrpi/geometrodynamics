# Prospective derivative-reconstruction refinement

Date: 2026-09-16. Original public freeze: #305, `66bce68`.

The original experiment remains **6/8**. All 27 registered BVP solves
converged and six gates passed. The largest off-grid second-order equation
residual was 4.85645e-7 against 1e-7. At the finest registered coordinate
step the normalized physical Hamiltonian residual was 1.71258e-4 against
1e-5, with last error ratio 1.6115. Neither failed gate is relabeled.
Momentum and localization checks passed, including independent physical
momentum residual 9.97e-12 and 1.94% bulk radius deviation at L=5.5.

Original decompressed raw SHA256:
`e510ba44aeb7d8d8d25cfc28cfb54a6a402d45ba43d373d5f6b48b6bfdffd227`.
Original module SHA256:
`424026b1caf05e15ddb423a49d66938654388bf953e6c2f4a773988cf3229f17`.
The exact original module/probe snapshots accompany the raw archive; the
hashes in this document must be checked against them before publication.

## Hypothesis, not yet a result

The cubic collocation state interpolant is C1. Checking a second-order ODE
by differentiating its scalar component twice exposes a lower-accuracy,
piecewise second derivative. Conformal curvature multiplies this error by
psi^-5 at a small neck. The independent finite-difference curvature check
also suffers second-order coordinate truncation. The original fine scalar
values agreed with medium values to about 1e-10, but that alone does not
establish correct curvature or prove this explanation.

## Frozen refinement

No new BVP branch search, new scalar/momentum profile, amplitude, cap size,
action, or relaxed physical/error threshold is authorized by this extension.
Use only the nine archived finest solutions from the original run.

1. On each original BVP mesh interval, construct the unique quintic Hermite
   interpolant matching psi and psi' from the stored state at both endpoints
   and psi'' from the registered ODE at both endpoints. Use these as interpolation
   data, not as a replacement for evaluating an off-grid equation. This is
   C2 reconstruction of the same numerical trajectory. Store its polynomial
   coefficients and independently differentiate them in all validation.
2. Check all nine reconstructed data sets at the original 1001 off-grid
   points. Demand max |H|<1e-7 as before, max |M|<1e-9 and psi>0. Check all
   interpolated scalar values against original fine data to relative 1e-8;
   knot values/slopes must agree to absolute 1e-12. Retain original
   medium-to-fine comparisons and boundary tolerances.
3. Use exactly the same 20 physical off-grid coordinates at L=5.5,eta=.3.
   Build metric jets with fourth-order centered first/second differences;
   mixed derivatives use tensor products of the fourth-order first-derivative
   stencil. Coordinate steps: .008,.004,.002. Form Christoffels, their
   derivatives, Ricci scalar, momentum divergence and four-field source
   contractions from these jets. Do not substitute the Hamiltonian equation
   into curvature or use analytic conformal curvature as the verification.
4. Require normalized physical H and M residuals <1e-5 at the finest step.
   Require last error ratios in [8,24] when preceding error exceeds 1e-8.
   Otherwise mark the ratio as roundoff/small-error limited, without using
   that designation to waive the absolute residual criterion. The wrong-f
   physical control must remain >10x worse.
5. Repeat seam and localization diagnostics on each reconstruction with the
   original bounds. Reuse the original time-reversal control and confirm
   that changing interpolation does not change its sign relations. Repeat
   null contractions with the same action; no traversability gate is added.
6. Tamper with a saved reconstructed polynomial and with an archived source
   or physical K/source sample: dependent verdicts must fail. The extension
   must verify the original raw SHA256 and original 6/8 gate pattern first.
   A changed original run cannot be treated as the frozen experiment.

The refinement verdict separately combines original passing action,
momentum, transport, localization and control evidence with the new checked
derivative reconstruction and physical residuals. The original raw archive,
6/8 verdict and nonzero original command exit remain unchanged. Publish
this extension before constructing or measuring any refined interpolant.
If the refined checks fail, report failure; do not keep changing resolution
until a desired headline appears.
