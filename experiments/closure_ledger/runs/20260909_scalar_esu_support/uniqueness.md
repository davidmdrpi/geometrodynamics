# Post-freeze scalar-support uniqueness

This extension was not a prediction in freeze de55f3f.

**ONLY_HOMOGENEOUS_SUPPORT_IN_STATED_CLASS**

Alignment requires B != 0. At B = 0 only Bd = 0 follows from the pointwise matrix equation.
The independent zero-momentum route also excludes every nonconstant smooth support.

| Supplementary gate | Pass |
|---|---|
| coefficient_expansion | True |
| rank_and_alignment | True |
| regularity_contradiction | True |
| independent_momentum_route | True |
| degenerate_controls | True |
| failure_paths | True |

The full component-boundary and zero-slice continuation proofs are in
docs/scalar_esu_uniqueness.md. Numerical controls are not a proof by search.

The even homogeneous family remains outside the imposed odd sector.
Support selection and Phi selection are NOT_DERIVED; the causality gate remains OPEN.
