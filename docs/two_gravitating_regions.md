# Two independently measured gravitating regions

PR #312 follows the merged finite-packet experiment (#311). It adds **two
separate source regions in a single nonlinear Einstein/quartet initial-data
solve**, with independent region IDs, centroids, energy inventories, and six
rotation-momentum observables. It computes the instantaneous momentum rates
and checks them against a stress balance including the exterior.

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
for g=psi^4 gamma. The paired solution is not the sum of the two one-seed
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
signals in a spacetime and are not evidence of causal transmission.

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

## Next required experiment

Evolve these constrained data with a nonspherical Einstein/quartet solver,
monitor both constraints, and retain independently identified moving regions
with the corresponding window-motion terms. Then compare controlled source
changes after their light-travel time through the exterior, include field
and gravitational boundary terms in a specified momentum convention, and
resolve the unstable homogeneous background response. Initial elliptic
cross-response and initial time jets cannot replace that experiment. Adding
strong-field mouths, if desired, also needs a new trapped-surface assessment;
traversability is not required.
