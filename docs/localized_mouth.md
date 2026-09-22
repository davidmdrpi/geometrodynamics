# Localized four-scalar handle: results and remaining failure

The localized initial-data experiment is implemented, but its full registered
certification is **not passed**. The original experiment passes **6/8** gates.
The prospectively registered derivative reconstruction passes **7/8**. Its
independent physical Hamiltonian convergence-rate gate remains failed. Both
`FOUR_SCALAR_HANDLE_CONSTRAINT_DATA` and `LOCALIZED_BULK_MOUTH_INITIAL_DATA`
remain false. Small absolute errors do not override the frozen rate criterion.

This extends the compact vacuum handle benchmark in #304 by retaining the
repository's existing four conformal scalars and their round-S3 bulk profile,
with nonzero scalar momentum, a solved momentum constraint and a nonlinear
Hamiltonian BVP. It does not derive those four fields from vacuum GR.

## Registration and evidence history

- Original specification: [localized_mouth_prereg.md](localized_mouth_prereg.md),
  published in #305 at `66bce68ce5c9e07f4369283b18cb1ac4d3e94f73` before runs.
- All 27 registered BVP solves converged. The original off-grid Hamiltonian
  and independent physical gates failed; the original archive and 6/8 verdict
  are retained without replacement.
- Derivative reconstruction was then frozen publicly at
  `dbec68f2b6de8f745c269167a2ac7a41a38f4654` in
  [localized_mouth_refinement_prereg.md](localized_mouth_refinement_prereg.md).
  It uses the nine original finest trajectories, with no further BVP solves.
- The first implementation formed Bernstein polynomials and converted them
  to power coefficients. Endpoint slope cancellation gave a 3.28e-12 knot
  mismatch against 1e-12. That attempt, its 6/8 report and module snapshot
  are retained as `bernstein_*` artifacts. Direct scaled Hermite coefficient
  arithmetic fixes this implementation defect while preserving the unique
  registered interpolant, original knots, schedule and acceptance bounds.
- The stable implementation gives 7/8. The frozen extension says to stop if
  it fails; no extra resolution or revised threshold was tried afterward.

Artifacts are in
[`experiments/closure_ledger/runs/20260916_localized_mouth/`](../experiments/closure_ledger/runs/20260916_localized_mouth/).
`source_manifest.json` records code, snapshots and compressed/decompressed
archive hashes. The original raw SHA256 is
`e510ba44aeb7d8d8d25cfc28cfb54a6a402d45ba43d373d5f6b48b6bfdffd227`.
The stable refinement raw SHA256 is
`b9fd7d2e652ccdae2012c636c7fc249fc8b5cca05cb872df5b32576f591cb3ae`.

Publication note, 2026-09-20: the initial publication of #305 carried the
hashes but not their corresponding evidence. They were therefore not
externally verifiable then. The implementation now supplies the original
snapshots and all three raw archives as lossless `.partNNN` files, with
`.parts.json` manifests. Concatenation in manifest order recovers the exact
original compressed bytes. The probe readers verify both individual parts
and the assembled archive automatically; their command-line paths remain
`probe.json.gz` and `refinement.json.gz`. No numerical samples were changed.

## What was solved

For the original Jordan-frame action, with kappa=1 and Lambda=3/2,

\[
I=\int\sqrt{-g}\left[\frac{R-2\Lambda}{2}
-\frac12\sum_A(\nabla\phi_A)^2-\frac{R}{12}\sum_A\phi_A^2\right],
\]

the initial quartet has constant norm 3/4 and normal momenta tangent to its
field-space sphere. Thus f=1-sum(phi^2)/6=7/8 and n(f)=0 on this slice.
This restriction is not imposed on subsequent evolution. In the Einstein
frame, gE=f gJ, G_AB=delta_AB/f+phi_A phi_B/(6f^2), U=Lambda/f^2.

On the compact cylinder with `(L,n)~(-L,-n)`, let
phi=q(sin(theta),cos(theta)n), q=sqrt(3)/2, and
PiE=psi^-6 p dphi/dtheta. The smooth cutoff leaves theta=asin(tanh(s))
in the bulk and makes theta constant near the seam. The scalar fields and
their momenta use a specified sign line bundle across the identification:
phi(L,n)=-phi(-L,-n). This is an additional global sector choice, not a
consequence of the constraints or an ordinary scalar periodicity rule.

With gammaE=psi^4(ds^2+dOmega^2), trace KE=0 and
KE=psi^-2 diag(2a,-a Omega), the coupled equations are

\[
2a'=-\alpha p\theta',\qquad
\psi''-\frac{2-\alpha(\theta'^2+2\cos^2\theta)}8\psi
+\frac{6a^2+\alpha p^2}{8}\psi^{-7}+\frac U4\psi^5=0,
\qquad \alpha=6/7.
\]

The metric obeys Neumann conditions at the center and seam. The odd momentum
source satisfies the compact compatibility condition; deliberately making
it even fails that condition. Removing the momentum correction also fails.
The collar, initial normal momenta and throat integration constant are free
initial data. Their preparation has not been dynamically explained.

For these constant-norm, tangent-momentum data the original improved stress
has the ADM projections

\[
f[R_J+(\operatorname{tr}K_J)^2-|K_J|^2]
=|\Pi_J|^2+|D\phi|_J^2+2\Lambda,
\qquad
f[D_jK_J{}^j{}_i-D_i\operatorname{tr}K_J]
=-\sum_A\Pi_{J,A}\partial_i\phi_A.
\]

The coordinate check forms Christoffels, Ricci curvature, extrinsic-curvature
divergence and all four field gradients directly from metric/K/field samples.
It does not replace curvature by the Hamiltonian equation. Reconstruction
uses the ODE only to specify endpoint second derivatives; off-grid residuals
differentiate the saved polynomial independently.

This reconstruction is an equation-informed consistency check of a new
interpolant, not an independent replication of the original cubic residual
test. The ODE fixes its endpoint accelerations, so residuals at those knots
vanish by construction, up to arithmetic error. Away from knots the check
is substantive, but the unchanged numeric bound does not give it the same
evidential independence. The original 6/8 failure is retained, and the
independent coordinate-curvature gate remains mandatory and failed. No
threshold has been tightened retrospectively to answer the review.

## Registered outcomes

| Quantity | Original | Stable reconstruction | Requirement |
|---|---:|---:|---|
| Gates passed | 6/8 | 7/8 | All required gates |
| Maximum off-grid Hamiltonian residual | 4.86e-7 | 4.40e-10 | <1e-7 |
| Maximum off-grid momentum residual | 3.04e-14 | 3.04e-14 | <1e-9 |
| Finest normalized physical Hamiltonian residual | 1.71e-4 | 6.02e-8 | <1e-5 |
| Finest normalized physical momentum residual | 9.97e-12 | 1.13e-15 | <1e-5 |
| Last physical Hamiltonian error ratio | 1.61 | 4.71 | Original [2.5,5.5]; extension [8,24] |
| Seam mismatch | 2.78e-17 | 2.78e-17 | <1e-7 |

The stable physical Hamiltonian errors at h=.008,.004,.002 are
4.10573e-6, 2.83434e-7 and 6.02144e-8. The preceding error exceeds 1e-8,
so the final ratio 4.70709 must meet [8,24]; it does not. Although the momentum
errors have a ratio about 16, they are below the small-error trigger and do
not require that rate for a pass. The wrong-f Hamiltonian control is 0.0714.
Reconstruction changes scalar values by at most 2.17e-12 relative and
matches knot values/slopes to 5.55e-17.

A diagnostic using differentiated saved polynomials finds a small nonzero
physical Hamiltonian residual near the narrow neck. The conformal residual
is multiplied by -8 f^2 psi^-5 when converted to this Jordan constraint.
At s=.93L for L=5.5, eta=.3, this gives about -2.69e-7 before normalization.
This provides an error-floor explanation for the stalled coordinate rate;
it is a diagnosis, not an independent validation or a waiver of the gate.

For eta=.3, the geometric localization diagnostics are:

| L | Maximum bulk radius-factor error on abs(s)<=1 | Neck/central radius |
|---|---:|---:|
| 3.5 | 15.04% | 0.25226 |
| 4.5 | 5.36% | 0.09031 |
| 5.5 | 1.94% | 0.03284 |

At L=5.5 the central Jordan areal radius is 0.980613 and the neck radius
0.0322038. Nearby sections have larger radii. This supports a local section
minimum, not a general minimal-surface stability theorem. The bulk scalar
profile is retained; the solved bulk metric is close to, not exactly, round S3.
For eta=.3 both future null expansions at the neck are about +0.00140133:
the section is anti-trapped. For eta=0 both vanish to numerical precision.
Neither case proves an outermost horizon or a traversable mouth. The regular
Einstein-frame kinetic matrix is positive and sampled null contractions
are nonnegative. No negative-energy support channel has appeared.

## Reproduction and software validation

From the repository root, with project dependencies installed:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.localized_mouth_probe --output-dir /tmp/localized-original
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.localized_mouth_refinement_probe --original experiments/closure_ledger/runs/20260916_localized_mouth/probe.json.gz --output-dir /tmp/localized-refinement
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_localized_mouth.py tests/test_mouth_momentum.py
```

Both scientific commands intentionally exit **1** for their failed gates.
The refinement also accepts `--rescore <archive>` to evaluate saved evidence.
An independent replay of the stable refinement reproduced its raw archive
byte for byte. The targeted suite passes **73 tests**, including source/frame
corruption, missing evidence, polynomial tampering, preservation of both
failed scientific verdicts, stale-output withdrawal and Python 3.10 grammar.
Seven additional packaging checks verify exact compressed/decompressed hashes and reject missing, reordered or corrupted evidence parts. Software test success does not change the scientific gate results.

## Reproducibility repair, 2026-09-22

The first GitHub CI run failed two tests because it rescored the refinement
as 6/8 instead of 7/8. This was a software evidence-validation failure, not
a newly failed constraint. With NumPy 2.5.3 and SciPy 1.18.1, disabling the
AVX-512 dispatch groups locally reproduced it: the scientific residuals
and controls retained their results, but exact
dictionary equality against freshly reconstructed polynomial coefficients
failed the evidence gate. The underlying physical gate remained false.

All nine reconstructed coefficient records differed in their last bits.
Across the 1001 comparison points their scalar values were identical; the
largest second-derivative difference was 1.67e-16. This explains why matching
library versions alone did not reproduce the discrepancy: CPU-dispatched
floating-point operations can take different arithmetic paths.

The repaired validator keeps metadata, mesh and polynomial layout exact.
It bounds differences over every complete interval, for psi through its
second derivative and the stored velocity through its first derivative.
For coefficients c_k and interval width h, the derivative-order-r error is
bounded by the sum of abs(delta c_k) k!/(k-r)! h^(k-r). That sum bounds the
error at every point of the interval, including between validation samples.
The arithmetic budget is 32 float64 eps times the larger of one and the
corresponding expected coefficient sum. This roundoff budget is distinct
from, and far below, the unchanged knot and PDE acceptance bounds.

The reproducing AVX-512-disabled case uses 0.789 of that arithmetic budget
and now scores 7/8. Native dispatch also scores 7/8. Regression tests run
the scorer in a fresh process with the current NumPy build's optional CPU
dispatch groups disabled. Further tests retain rejection of changed mesh,
metadata, array shape, velocity and curvature; the curvature control keeps
endpoint values/slopes within 1e-12 and is still rejected by the interval
derivative bound. CI checks these tests on Python 3.10 and 3.12 before the
full suite.

The solver, all raw evidence, both original reports and both public freeze
files are unchanged. Original 6/8 and refined 7/8 remain their recorded
outcomes, and both milestone verdicts remain false. Re-scoring now reports
the fraction of the arithmetic reconstruction budget used; this diagnostic
is not a new physical acceptance criterion. The source manifest preserves
the pre-repair hashes alongside current source hashes.

## Consequence for the quantum-foundations claim

There is now a concrete numerical momentum solve with nonzero quartet
current on a slice carrying a localized antipodal handle, accompanied by
independent physical momentum checks. Full localized initial-data
certification remains unresolved because of the Hamiltonian rate failure.
An evolved throat/worldline crossing, reciprocal momentum-transfer event,
traversability, selection of the sign sector, discrete action and quantum
statistics remain unestablished. There is no event detector or threshold
that could manufacture measurement discreteness in this experiment.

The immediate numerical follow-up would require a new prospective
specification that separates solver/interpolation error from coordinate
truncation in the physical curvature test. It must not retroactively change
these verdicts. Crossing dynamics remains a subsequent, distinct milestone.
