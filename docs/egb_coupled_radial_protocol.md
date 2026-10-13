# Coupled radial validation protocol

This is a post-derivation validation protocol, not blind preregistration.
Exploratory shooting and a finite-element check already showed negative kinetic
signs and growing modes when this document was written. No selection threshold
is being presented as fixed before that exploration.

Keep #331's action and benchmark R=1, a=.3, alpha=-.06, eta=.1, bulk lambda=0.
Do not change the material law, add control kicks, or refit the geometry. Solve
linear, spherically symmetric perturbations of both parity sectors on the
retained compact bulk, including moving-shell matching. Nonspherical and
nonlinear evolution are outside this calculation.

The physical kinetic form must come from the action's symplectic form, including
both gravitational boundary faces and one surface scalar. An auxiliary spectral
pencil mass or the old pressure coefficient alone is not an acceptable norm.
Remove the constant scalar shift and the exceptional odd coordinate mode.

Validation gates:

- Check the full orthonormal Einstein/GB tensor and full Davis tensor by
  independent contractions on the low modes (relative 1e-6).
- Check density, pressure, scalar-current matching and total delta Q on each
  computed even mode (absolute residual 1e-7 with z=1 normalization).
- Match the physical quadratic form, differentiated Robin form and direct
  fixed-shell symplectic form in two arbitrary radial gauges (relative 1e-7).
- Recover the gravitational kinetic term from the full warped-curvature action
  plus Myers terms using decreasing amplitudes; final relative error <2e-7.
- Refine independent linear finite elements through 80,160,320 intervals. For
  the first four even eigenvalues the error must decrease by a factor <.3
  per doubling and the final absolute error must be <4e-4. Do not claim the
  entire infinite spectrum numerically resolved from this finite matrix.
- Repeat shooting with tighter tolerances; low eigenvalues agree within 1e-7.
- Exclude a spurious root at the rational Robin pole by using a pole-free
  determinant. Treat h=0 separately, not by division or discarding every
  negative eigenvalue as gauge.
- Track kinetic inertia separately from squared frequency: a real-frequency
  negative-norm mode fails the positive-energy gate; a negative squared
  frequency fails linear stability even with positive norm. Complex roots
  would also fail stability.

Failure of the numerical identities gives NUMERICALLY_UNRESOLVED. If checks
pass and a physical negative-norm mode or growing mode survives, report
COUPLED_RADIAL_ADMISSION_FAILED with the particular reasons. A positive result
would require positive physical norm and no unstable radial modes, not merely
a frequency close to #330's oscillator. No conclusion about every alternative
focus material, finite-thickness regularity or emergence of quantum events is
licensed by this test.
