# Diagonal-loop group-action test (11 October 2026)

This additive study follows #325 and the independent review of its results.
The #319 two-return loop and its unresolved breaking at approximately 1e-13 are
already known. No new loop-action evaluations have been inspected before this
protocol. This is a new test on historical data, not a new independent discovery
of the loop. Frozen historical producers and archives remain unchanged.

## Hypothesis and outcome

On the particular archived diagonal-circular component, the half-clock H
preserves the component and has order four; cyclic axis permutation C preserves
it and has order three. They commute. A reflection exchanges the two chiral
components. If supported, a common circle angle makes the orientation-preserving
group cyclic of order twelve. The first *permitted* scalar resonant angular
harmonic is then 12, not 4 or 6. This predicts no absolute amplitude, no resolved
twelfth harmonic at double precision, and no extra integral.

We do not assume which oriented lift, 1/4 or 3/4, occurs: orientation is the
geometric azimuth atan2(x2,x1), unlike the LRS Fourier convention. Determine it
from four successive H steps. Variation of individual azimuth steps is allowed;
geometric azimuth is not asserted to be a uniform normal-form angle.

## Fixed computation

Use all 110 continuation points in stage_F.json to fit the first six coordinates;
its concatenated twelve entries are **two six-dimensional nodes**, not one full
homogeneous state. Do not append a duplicate endpoint or enforce any symmetry
on the fit. Compare Fourier degrees 20 and 26; also fit the even and odd halves
separately at degree 20 and predict the withheld half. Degree selection and
thresholds precede the new action evaluation.

Sample indices floor(110*j/24), j=0,...,23. For each, compute four successive
half returns, one direct full return, H(Cz), C(Hz), the reflected point, and an
independent full-matrix half return. At j divisible by four also use full-matrix
Radau (six checks). Primary: explicit eight-state diagonal equations with
DOP853; all methods rtol=1e-13, atol=1e-15, max_step=.05. Checkpoint each job.
Use the unmodified full Einstein-quartet equations for matrix checks.

The coordinate norm is the maximum absolute component in dimensionless
(A,p_A,x1,p1,x2,p2). Interpolation uncertainty e is the maximum of training,
even/odd withheld errors, and degree-20/26 disagreement on 1024 azimuths.
The membership gate is max(20e,1e-9), subject to e<=5e-9 and gate<=1e-7.
H-images (including H^4) and C-images must lie within that gate of the fitted
component at their own azimuth. This is a sampled, finite-resolution membership
test, not a proof of an exact invariant continuum.

Additional fixed gates:
- H^4 identity <=1e-7; H and H^2 and C each move the point by >1e-2.
- Four positive modulo-one azimuth increments sum to the same integer 1 or 3
  at all samples, within 1e-6; each increment is within .10 of that integer/4.
- H^2=P, P matches the paired archived node, and HC=CH, each <=1e-10.
- Reduced/matrix DOP853 and reduced/matrix Radau differences <=1e-10;
  absolute Hamiltonian constraint residual <=1e-10.
- x1*p2-x2*p1 has one sign over the archived component, with magnitude >1e-3;
  the reflected samples have the opposite sign with the same margin. This is
  a coordinate chirality diagnostic, **not** physical SO(3) angular momentum
  or a conserved quantity.

If interpolation, independent-map, or constraint gates fail, report
NUMERICALLY_UNRESOLVED and identify the failed gates. If those checks pass but
any proposed group-action gate fails, report PROPOSED_GROUP_ACTION_FAILED.
Only all gates passing permits ORDER_TWELVE_SELECTION_SUPPORTED. These labels
refer to the group action; no label about measured breaking is issued.

## Interpretation and next falsifier

In a common conjugate angle, invariant scalar Fourier terms must satisfy both
4|k and 3|k. Reflection relates coefficients on opposite components; it does
not by itself force the k=12 coefficient to vanish on either one. Raw Euclidean
phase-slice lambda and geometric-azimuth spectra can have sidebands, so their
nonmultiples of twelve do not falsify this normal-form selection rule.

After measuring the actions, specify a high-precision test that can distinguish
a nonzero order-12 coefficient from additional suppression, including numerical
floor and control. Such a test requires a covariant obstruction or a common
symmetry angle and must not infer an absolute coefficient from selection alone.
No arbitrary corrective kicks, new fine tuning, or long-time evolution are
introduced: the horizon is exactly four half-clock transits of historical
periodic initial data. The existing physical hyperbolic instability remains.
