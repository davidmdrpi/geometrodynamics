# Disclosed follow-up after the registered seed-locality failure

11 October 2026. Written **after the complete 192-phase scan was inspected**,
and before any supplemental root search. The original protocol, producer,
archives and score remain unchanged. Registered result:
**NUMERICALLY_UNRESOLVED**, because the seed-locality gate failed.

All numerical convergence, constraint, error, conditioning and control gates
pass. The independent formulations give max|lambda|=2.1600485904e-18,
a dominant raw harmonic of 12 and 24 reliable sign-change brackets. The
estimated error is 1.1071e-47; the registered detection threshold remains
1e-28. These numerical observations do not override the failed locality gate.

## Post-inspection diagnosis

The historical `control_seeds` splines parameterize each six-dimensional node
curve by its **own** Euclidean chord length. Equal spline parameter therefore
does not generally select an actual paired return. The original protocol
recognized that the splines were only seeds, but its blanket .01 displacement
cap was not calibrated to this mismatch. This is a design error in the new
locality gate, not a defect newly discovered in the historical GR equations.

- 131/192 phases fail the registered .01 seed-displacement gate.
- Maximum original displacement: 0.0278719303644, from node one.
- Maximum node-zero displacement: 3.4888656861e-6.
- Reusing node zero's chord grid for **both** node coordinates gives maximum
  displacement 3.4888656861e-6. The node-zero anchor and tangent are exactly
  unchanged in this audit, so the augmented equations are identical.
- All solved nodes lie within 9.7506e-10 of the unsymmetrized degree-20 loop fit,
  inside #326's previously registered 1.9919e-8 membership gate.

These are post-hoc diagnostics, not a replacement registered locality test.
No high-precision phase is rerun with the more favorable seeds, and no score
is silently recomputed with a relaxed gate. This check does not turn the
original result into a preregistered success.

## Additional root calculation, fixed before its execution

Run the **unchanged frozen root algorithm** on all 24 already observed sign
brackets, using `r3_diagonal_obstruction_supplement.py`. This is an additional
analysis selected after seeing the scan, rather than the registered conditional
root stage (which correctly skipped the roots after the failed gate).

All frozen root thresholds are retained: augmented residual <=1e-56,
|lambda|<1e-51; two formulations/finite-difference steps for the 12x12
unaugmented Jacobian; inverse-norm times Jacobian discrepancy <.01;
independent closure <1e-45; inferred correction <1e-12; inverse residual <1e-40;
and distinct node-zero solutions separated by >1e-6. Retain failures and record
symmetry equivalence classes and matching errors. No elliptic/hyperbolic
classification is claimed.

Only all isolated-root checks passing allows the separate label
SUPPLEMENTAL_ROOTS_VERIFIED; otherwise SUPPLEMENTAL_ROOTS_UNRESOLVED. This
label never changes `result.json` or its NUMERICALLY_UNRESOLVED verdict.
The supplemental archive binds this note, its driver and the original input
archives. The complete registered data, its failure, the seed audit and this
follow-up are published before starting the supplemental roots.
