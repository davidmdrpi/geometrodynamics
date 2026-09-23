# Independent physical error budget: 8/8

The prospective extension passes all eight gates. Its separate verdicts are
`FOUR_SCALAR_HANDLE_CONSTRAINT_DATA = true` and
`LOCALIZED_BULK_MOUTH_INITIAL_DATA = true`. This completes the registered
numerical initial-data checks for the chosen four-scalar handle family.

The original experiment remains 6/8 and its derivative reconstruction remains
7/8. Neither historical verdict has been replaced. The change is a new
physical verification experiment, not a new solution or a relaxed absolute
accuracy threshold.

## Public freeze and method

[PR #307](https://github.com/davidmdrpi/geometrodynamics/pull/307) published
[the specification](localized_mouth_error_budget_prereg.md) at
`f9025668b2d6ad4bd94a2f704aa699e8849831a7` before measurements.
There are no new BVP solves, changed profiles or interpolants.

The checker imports the original binary64 polynomial coefficients exactly,
evaluates the Jordan metric, K and four fields with arbitrary precision,
and differentiates these coordinate functions directly. Christoffels,
their derivatives, Ricci curvature and momentum divergence are assembled
by coordinate contractions. Curvature never substitutes the Hamiltonian
equation or the conformal scalar-curvature identity. The underlying
reconstruction remains ODE-informed; this is an independent check of that
saved representation, not an independently solved BVP.

Both 60- and 80-digit references agree within the registered scaled 1e-40
bound. A second calculation uses fourth-order differences of function
values. Each point's three steps are fixed by its distance to the nearest
saved polynomial knot: h0=min(.008,d/8), followed by h0/2 and h0/4.
Thus all local stencils stay within smooth polynomial pieces. The five
h0 values range from 6.71e-6 to 4.70e-5; no residual-dependent step selection
or resolution search was used.

For each of the same 20 physical points, the test compares signed H and M
against the independent reference BEFORE taking norms. Both successive
ratios of differentiation errors must lie in [8,24], apart from the fixed
1e-40 small-error exception. Absolute residuals of BOTH the reference and
the finest difference estimate must remain below 1e-5.

## Results

| Quantity, maximum across registered points | Result | Bound |
|---|---:|---:|
| Normalized reference Hamiltonian residual | 3.95608074265e-8 | <1e-5 |
| Normalized finest-difference Hamiltonian residual | 3.95608074456e-8 | <1e-5 |
| Normalized reference momentum residual | 5.38387849960e-16 | <1e-5 |
| Normalized finest-difference momentum residual | 5.38387849811e-16 | <1e-5 |
| Finest Hamiltonian differentiation error | 1.91046665744e-17 | <1e-8 |
| Finest momentum differentiation error | 1.68551259348e-25 | <1e-8 |
| All Hamiltonian differentiation ratios | 15.9999999851 to 16.0000000017 | [8,24] |
| All momentum differentiation ratios | 15.9999999988 to 16.0000001981 | [8,24] |

Every point passes the wrong-f control and precision, absolute accuracy,
differentiation accuracy and convergence checks. Calibration reproduces
R=6 for the round S3, R=2 for the cylinder, and nonzero geometric momentum
divergence (-2,0,0) for the prescribed cylinder tensor. Additional tests
verify the exact round quartet constraints and a nonzero scalar-current sign.

## What caused the plateau

At s=.93L, the directly measured signed Hamiltonian residual is
-2.69411881737e-7, or -3.95608074265e-8 with the registered normalization.
This independently confirms a nonzero residual of the saved representation.
It is small enough for absolute accuracy but prevents the total residual
from continuing to decrease toward zero at fourth order.

The original larger stencils were also evaluated at 80 digits as diagnostics.
Their maximum normalized Hamiltonian residuals are 4.10706062508e-6,
2.84703997405e-7 and 5.48086669946e-8. The original binary64 values remain
4.10573192702e-6, 2.83434492388e-7 and 6.02143793503e-8. Thus the old sequence
also contains finite-precision effects; it was not a clean measurement of
one constant floor. Larger stencils can cross polynomial knots. No unique
floor was inferred from ratios of maxima at changing points.

The independent review of #307 confirmed the result and corrected its
earlier knot table: the step rule uses all three meshes, not just psi.
The distances below follow directly from the archived h0 values (d=8h0):

| s | Nearest knot across psi, theta and tensor | h0 | Original reach 2h=.004 crosses a knot? |
|---|---:|---:|---|
| .385 | 7.16146e-5 | 8.95182e-6 | yes |
| 1.265 | 1.07422e-4 | 1.34277e-5 | yes |
| 2.805 | 5.37109e-5 | 6.71387e-6 | yes |
| 4.345 | 2.14844e-4 | 2.68555e-5 | yes |
| 5.115 | 3.75977e-4 | 4.69971e-5 | yes |

All five radial locations therefore cross a polynomial boundary in the
original finest stencil. That invalidates an unqualified smooth-piece
fourth-order truncation argument; it does not by itself measure the size
of that contribution. The measured nonzero solution residual and finite-
precision differences are separate effects. Every new H and M ratio is
active at every point: none passes through the 1e-40 exception.

For subsequent evolution the review correctly requires independent scalar
amplitudes: fixed initial norm is not preserved by the field equations.
The matched round background breathes, and its homogeneous gravitational
mode is unstable. These must remain explicit controls in an evolution test.

The registered solution residual and differentiation error are now measured
separately. Subtracting the reference cannot hide an inaccurate solution:
absolute constraints must pass independently, and a regression test rejects
a large fabricated constant residual despite fourth-order differentiation.

## Reproduce and audit

Install the package with development dependencies, then run from its root:

```sh
python -m experiments.closure_ledger.localized_mouth_error_budget_probe \
  --input-dir experiments/closure_ledger/runs/20260916_localized_mouth \
  --rescore experiments/closure_ledger/runs/20260922_localized_mouth_error_budget/error_budget.json \
  --output-dir /tmp/localized-mouth-error-budget-replay
pytest -q tests/test_localized_mouth_error_budget.py tests/test_localized_mouth.py tests/test_mouth_momentum.py
```

Omit `--rescore` to measure the same registered data afresh. Rescoring
recomputes every measurement from hash-verified input archives, verifies
the source hashes, and requires scaled 1e-40 agreement with saved evidence.
Changing a point, step, K, source, polynomial, or constraint sample invalidates
evidence. Failed CLI runs clear any previous affirmative output.

The combined targeted suite passes **89 tests**, including a full replay of
all saved measurements and preservation of the two historical failed verdicts.
CI runs this suite on Python 3.10 and 3.12 before the full repository suite.

Raw measurements, source and input hashes, dependency versions and the
pointwise verdict are in
[`runs/20260922_localized_mouth_error_budget/`](../experiments/closure_ledger/runs/20260922_localized_mouth_error_budget/).
The original raw archives and both older specifications are unchanged.

## Scope

The chosen quartet, sign bundle, collar and initial momenta remain assumptions.
Localization is already present at zero momentum and is not evidence for
momentum-driven throat formation. These results do not establish derivation
from vacuum GR, traversability, evolved worldline crossings, operational
momentum-transfer events, action discreteness or quantum statistics.
