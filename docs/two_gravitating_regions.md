# Two separately measured gravitating regions

PR #312 follows the merged finite-packet experiment (#311). It adds **two
separate source regions in a single nonlinear Einstein/quartet initial-data
solve**, with independent region IDs, centroids, energy inventories, and six
rotation-momentum observables. It computes the instantaneous momentum rates
and checks them against a stress balance including the exterior.

The measurement windows are disjoint, but the scalar seed profiles overlap.
The initial momentum-rate signs do **not** isolate gravitational attraction.
The post-review controls below show that field-profile overlap dominates their
toward-the-other-region signs; the metric change opposes those signs.

It does **not** yet evolve these regions. These weak density perturbations
are not demonstrated self-bound objects, wormhole mouths, or black holes.
The four-scalar support remains assumed matter. No quantum milestone follows.

## Public freeze and numerical verdict

The [prospective specification](two_gravitating_regions_prereg.md) was
published as commit `27bf9c6d5840506bd1c786b4cddbf1d0ee2f79a6` in draft
[PR #312](https://github.com/davidmdrpi/geometrodynamics/pull/312), created
2026-09-27 19:48:25 UTC, before the new implementation and measurements.
The raw record is
[`initial_data.json`](../experiments/closure_ledger/runs/20260927_two_regions/initial_data.json).
It includes coefficients for all 18 solves, source hashes, controls, and
all registered gates. No case was dropped.

The registered verdict is **REGISTERED_GATE_FAILURE**, because the strict
all-successive-grid field-difference gate fails on the coarsest refinement.
The largest degree-12 to degree-20 change is 1.32470e-7, above 1e-7. The
largest degree-20 to degree-28 change is 1.51041e-10. The finest constraints,
instantaneous stress balance, independent metric responses, and independent
coordinate-curvature refinement pass. The failed gate was not loosened.
This is useful two-region infrastructure and numerically well-resolved finest
initial data, not an all-gates-passed experiment milestone.

## Geometry and physical meaning

On the unit S3 cover, take phi=q(x)x, with two independent smooth seeds in q.
The centers have quotient separation 1.2, with window radius .45 each.
Each region has two antipodal images. A and B are distinct locations even
after quotienting; they are not opposite orientations of one cut. The solve
retains SO(2) symmetry in the transverse plane, which permits independent
structure along the plane containing both centers. It is not spherically
symmetric, but does not support arbitrary three-dimensional dynamics.

The metric/stress are antipodally even. The quartet is odd and requires the
internal phi -> -phi twist to descend to RP3; it is not four ordinary even
scalar functions there. Solves use the full S3 cover. Reported invariant
inventories use half-cover integrals. The background and source formulas
are stated in the freeze; this new geometry has no handle or inner boundary.

With zero extrinsic curvature and scalar normal velocity, the momentum
constraint vanishes exactly. The Hamiltonian equation is solved nonlinearly
for g=psi^4 gamma. Newton iteration from the round solution selects the nearby
branch. The linearization there is -8 Delta - 96/7 and is invertible on the
retained even degrees; no global uniqueness of the focusing nonlinear
Lichnerowicz equation is asserted. The paired solution is not the sum of the two one-seed
metrics: max|psi_AB-psi_A-psi_B+psi_round|=7.51790e-6. The initial slice is
not static: the sigma-model acceleration and the Einstein evolution have
nonzero time derivatives beyond first order in the metric.

The time-symmetric choice also fixes the trapped-surface distinction:
for every sphere, theta_plus=H and theta_minus=-H. Thus this slice contains
no strictly future-trapped sphere. This does not reclassify the different,
non-time-symmetric neck data audited in #309.

## Separate measurements

The fixed smooth windows partition the cover into A, B, and the remaining
exterior. These are measurement windows, not material boundaries. The table
uses the finest field solution and registered (64,192) validation quadrature.
Energy includes the supporting field and potential; the excess subtracts
the exact round-background inventory with the same coordinate window.
These quantities are not ADM masses.

| Measurement | Region A | Region B |
|---|---:|---:|
| Seed amplitude | .02 | .03 |
| Proper energy inventory | .1117520072 | .1133888487 |
| Excess over round background | .0020903849 | .0037272297 |
| Initial rotation momenta (all six) | 0 | 0 |
| Initial Pdot_01 | +8.2697506e-5 | -6.4482017e-5 |

All momenta use the same six specified ambient rotation Killing fields
X_ab=x_a partial_b-x_b partial_a. This is a shared frame convention, not
an identification of tangent vectors at distinct points or ADM translation
momentum. Zero initial momenta follow from time symmetry; no P_B=-P_A
constraint is imposed. A synthetic localized-current test also verifies that
A can have nonzero measured momentum while B has zero momentum.

Independent seed increments .001 give the following finite-difference
energy response. Rows are measured regions, columns are perturbed seeds:

| | Seed A | Seed B |
|---|---:|---:|
| Region A | +.13977476 | -.02253229 |
| Region B | -.02241031 | +.14299561 |

Each seed also changes the solved metric in the other cap. Those responses
compare different elliptic initial-data solutions; they are not instantaneous
signals in a spacetime and are not evidence of causal transmission. The registered nonzero metric-response
gate is a sensitivity/bookkeeping check: it does not distinguish direct scalar
interaction from metric-mediated response. A nonzero elliptic response alone
is weak evidence for the intended mechanism. A finite threshold can still fail
for sufficiently small or cancelling responses; it is not literally guaranteed
to pass for every elliptic problem.

## Momentum-rate ledger and validation

The rate is computed directly from the sigma-model normal acceleration,
jdot_i=-G_AB Pdot^A partial_i phi^B. An independent integrated stress
calculation supplies a window-gradient flux plus metric work. The other five
rotation rates vanish by transverse SO(2) symmetry, not by an imposed
pairwise momentum identity.

The exterior rate is -1.8469931e-5. A+B alone gives +1.8215489e-5 and does
not close. Including the exterior gives -2.5444231e-7, equal to the integrated
metric-work term to floating-point accuracy. Matter momentum need not be
conserved against a spatial vector field that is not Killing for the
physical metric. This is **not a completed gravitational recoil ledger**:
no local gravitational momentum density or quasilocal gravitational boundary
charge is supplied, and no finite-time motion is measured.

| Numerical check | Result |
|---|---:|
| Largest finest normalized Hamiltonian residual, all six cases | 3.61e-11 |
| Largest pair rate-balance discrepancy, (32,96) quadrature | 2.31e-5 |
| Same, registered (64,192) quadrature | 2.97e-7 |
| Same, (96,288) quadrature | 1.28e-9 |
| Independent coordinate-metric Ricci residual at finest difference step | <=8.94e-6 |

The independent Ricci calculation finite-differences the coordinate metric,
its Christoffels, and their derivatives; it does not reuse the spectral
Laplacian or replace curvature by the constraint's matter source. Errors
shrink approximately fourfold on each halving of the difference step at all
three registered off-grid points. These checks are numerical, not rigorous
continuum existence/error certificates.

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.two_regions_probe
python -m pytest -q tests/test_two_regions.py
```

The experiment command intentionally exits 1 after saving the complete
record because the registered coarse-grid gate fails. Unit/regression tests
check the identities and implementation and can pass while that scientific
gate fails. To retain the original record, supply `--output /tmp/two-regions.json`.

## Post-review controls: profile overlap and metric response

The [review](https://github.com/davidmdrpi/geometrodynamics/pull/312#issuecomment-5861469221)
correctly identifies a missing mechanism control. These additions are explicitly
**post-review diagnostics**, not new preregistered gates. They reuse the saved
finest coefficients; the original solver, runner, freeze, archive and failed
verdict are unchanged. Reproduce with:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.two_regions_controls
python -m pytest -q tests/test_two_regions.py tests/test_two_regions_controls.py
```

The separate record is
[`review_controls.json`](../experiments/closure_ledger/runs/20260927_two_regions/review_controls.json).
It contains both (64,192) and (96,288) quadratures, source hashes, and each
control's Hamiltonian residual. The following table uses (96,288):

| Field / metric used | A: Pdot_01 | B: Pdot_01 |
|---|---:|---:|
| Paired field / solved pair metric | +8.269751e-5 | -6.448187e-5 |
| Paired field / fixed round metric | +1.529712e-4 | -1.117383e-4 |
| A-only field / pair metric | -6.982385e-5 | — |
| B-only field / pair metric | — | +4.725753e-5 |

The fixed-round field dynamics already produce larger rates with the original
signs. Subtracting this control from the paired solution gives -7.027371e-5
for A and +4.725643e-5 for B: the metric change opposes those signs. With only
each window's own seed retained in the pair metric, its rate also has the
opposite sign. These controls support the review's central interpretation:
the reported signs are dominated by scalar-profile effects and cannot be
reported as gravitational attraction. Their direction refers to rotation
momentum relative to the centers, not measured center acceleration.

Except for the paired field/paired metric row, these source/metric combinations
are **off shell**: their normalized Hamiltonian residuals are about .013–.020.
They are controlled evaluations of the matter equation on specified metrics,
not alternative constrained Einstein solutions. Their subtraction defines a
useful comparison convention, not a unique, gauge-invariant split into scalar
and gravitational forces or proof of repulsive two-body trajectories. Nonlinear
field/metric effects need not sum as independent forces.

The seed tails are substantial. The exact supremum of either normalized seed
inside the other radius-.45 cap is

    cosh(8 cos(1.2-.45))/cosh(8) = .1168947963.

The supremum is the boundary limit of the open window support. Both directions
have the same value by geometry. Different maxima such as 10.2% and 11.5% can
result from grid sampling; they are not different continuum separations.
Separate IDs and seed parameters mean independent bookkeeping, not isolated
field profiles or dynamically independent bodies.

The fixed-metric energy control makes the metric-mediated response clearer:

| Energy response | A inventory / seed B | B inventory / seed A |
|---|---:|---:|
| Constrained metric changes with seed | -.0225322880 | -.0224103105 |
| Metric fixed at round solution | +.0012285534 | +.0012317504 |
| Difference under this convention | -.0237608414 | -.0236420609 |

Thus the sign and magnitude of the off-diagonal energy response depend strongly
on solving for the metric. This is an initial-data metric effect on a proper
matter inventory, not a gravitational binding energy or causal exchange.

For any sigma-model normal velocity Pi, put T=G_AB Pi^A Pi^B and let Z denote
the physical squared spatial gradient. Directly from its stress tensor,

    rho = (T+Z)/2+U,
    tr S = 3T/2-Z/2-3U,
    rho+tr S = 2T-2U.

At this time-symmetric slice T=0, so rho+tr S=-2U<0; the spatial gradient
terms cancel even when the stress is anisotropic. Writing rho+3p is valid if
p means tr S/3, not an assumption of isotropic local stress. On the paired
data the contraction ranges from -3.9876924 to -3.9184426; the numerical
identity error is below 9e-16. This is negative active density for Ricci
focusing. It is consistent with the opposing metric contribution in these
controls, but does not by itself determine a regional force: spatial stress,
background subtraction, geometry and the measurement convention also matter.

For the exact breathing background, R=(sqrt(3)/2) cos(2 eta),
f=1-R^2/6 and Pi=R'/sqrt(f), the contraction is

    rho+tr S = 2 R'^2/f^3 - 3/f^2.

It is -192/49 at eta=0 and +3 at eta=pi/4, so the background Ricci-focusing
source changes sign. This is an analytic background prediction, not an evolved
prediction for the two regions' recoil signs.

## Next required experiment

Before evolution, preregister and solve new data with compactly supported seed
perturbations, or bound both their field and gradient tails in the other window
(the review suggests a field-tail target below 1e-6 of peak). Include fixed-metric
and field/metric controls, phase-resolved active-density diagnostics and the
breathing-background sign prediction above. Separate perturbation support does
not remove the common supporting background field. Do not silently replace the
seeds in the original frozen experiment.

Evolve the new constrained data with a nonspherical Einstein/quartet solver,
monitor both constraints, and retain independently identified moving regions
with the corresponding window-motion terms. Then compare controlled source
changes after their light-travel time through the exterior, include field
and gravitational boundary terms in a specified momentum convention, and
resolve the unstable homogeneous background response. Initial elliptic
cross-response and initial time jets cannot replace that experiment. Adding
strong-field mouths, if desired, also needs a new trapped-surface assessment;
traversability is not required.
