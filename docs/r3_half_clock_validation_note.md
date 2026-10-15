# Additional phase coverage for independent validation

2026-10-10. Added after freeze `82c832f6`, while the production scan was
running (23 of 84 phases had printed convergence flags). No new lambda value,
spectrum, noise estimate or registered score had been inspected. Existing
2/5 retrospective results had been inspected. The frozen producer and all
its gates remain unchanged.

Reviewing the validation grid reveals that j=1+12m on 84 samples always has
the same leading seventh-harmonic phase. Its seven checks can therefore
miss a phase-dependent map error. This is a limitation of the registered
noise estimate, even though they are different section points.

Additional checks, specified here before inspecting the new scan results:

- Chord re-solve **every one of the 84 phases** using the full matrix
  equations, direct full-return event detection, DOP853 1e-13/1e-15.
- Also re-solve j=**0,14,28,42,56,70,83** with full-matrix Radau
  1e-12/1e-14. These cover seven distinct seventh-harmonic phases.
- Use six chord iterations and tolerance 2e-13, with each original Jacobian.
- Record the maximum difference-plus-residual, the complete alternative
  spectrum, and the alternative correlation. No missing-point FFT.
- Report both the original frozen label and whether replacing its resolution
  by max(original resolution, 10 times the supplemental noise) would retain
  the prediction. Also require the full alternative scan to have dominant
  harmonic 7, 14 cyclic sign changes, and harmonic-7 amplitude at least ten
  times that combined resolution. If either check fails, the registered label
  is retained for provenance but the scientific finding is **not confirmed**.
  Confirmation also requires the registered support label and supplemental
  noise no larger than the original 1e-10 numerical cap.

This is a validation amendment after production began, not an original
preregistered gate. It adds no parameter search, rerun of the primary scan,
or change in the dynamical equations. `r3_half_clock_validation.py` writes
the checks to a separate exclusive archive with source and input hashes.
