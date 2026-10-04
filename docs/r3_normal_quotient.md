# Normal response of the two-return family

Measured 2026-10-03 UTC, following merged #319. The [specification](r3_normal_quotient_prereg.md)
was published at `c5c8c6be168408a82b86c9208c5182983c22bffc` before the new
quotient and perturbation calculations. All six registered gates pass:
**BOUNDED_CENTER_QUOTIENT_NUMERICALLY**.

The secular drift seen in #319 is removed from the tested linear normal
response when distance is measured to the entire family, allowing three
principal-axis rotations and displacement along the diagonal family loop.
Four neutral directions remain in the quotient centre; angular momentum
has not been discarded. The physical Einstein-static stable/unstable pair
is excluded only to define the conditional centre sector, and is reported
as **UNREDUCED_HYPERBOLIC_PAIR_PRESENT**.

This establishes numerical evidence for conditional linear normal
boundedness in the stated coordinates. It does not establish nonlinear
orbital stability, stability of an individual periodic history, an exact
symmetry generating the diagonal family, or action selection. The full
unreduced homogeneous system is still unstable.

## Interpretation after review (2026-10-04)

The [review of `36a8277`](https://github.com/davidmdrpi/geometrodynamics/pull/320#issuecomment-5976783259)
independently reproduced the label and the full-system perturbation response.
It also identified why boundedness of this quotient is largely an algebraic
consequence of #319's measured Jordan structure, rather than additional
evidence for dynamical stability of individual histories.

For a 2-by-2 unit Jordan block, removing its eigenvector leaves the identity
on the remaining generalized direction. A post-review check at sample 0,
using an ordered real Schur decomposition to exclude the physical hyperbolic
pair, gives numerical nullities 4, 8 and 8 for `(C10-I)`, `(C10-I)^2` and
`(C10-I)^3` at an absolute singular-value cutoff of 1e-7. This agrees with
four 2-by-2 unit Jordan blocks and the elliptic pair, with no resolved larger
unit block. These are numerical rank observations, not an exact Jordan-form
proof or a new registered measurement.

The substantive checks are that no larger secular chain is resolved and
that the removed subspace matches the rotation and family eigendirections.
Given that structure and the elliptic pair, the bounded quotient follows;
Q3/Q4 provide numerical checks of that consequence. Invariance in Q2 alone
would not identify a complete unit eigenspace or rule out larger blocks.

Bounded distance to the family does **not** imply persistence of exact
two-return closure. The retained neutral directions include angular-momentum
changes and the generalized direction conjugate to family phase (the local
action/detuning direction). They can persist while their accumulated rotation
or phase displacement is absorbed by the quotient. Along a detuned nearby
invariant circle, rotation need not remain one-half per return, so the exact
3/2 closure is lost. Thus the quotient result does not make resonant closure
robust to action perturbations. Local codimension-one language assumes a
nonzero frequency/action slope; this PR does not prove a global foliation,
an exact action-coordinate identification, or uniqueness of the resonant
family.

## The quotient and its norm

The section has 12 dimensionless coordinates `(A,p_A,x_1,p_1,...,x_5,p_5)`.
Let R be an orthonormal basis for the three analytic SO(3) tangents plus
the numerically identified diagonal-family tangent. Its complement Q gives
the norm `min_a ||delta z - R a||_2 = ||Q^T delta z||_2`. Invariance of R is
checked before forming `B = Q^T M2 Q`, an eight-dimensional quotient map.

This Euclidean norm mixes scale-factor, momentum and tensor coordinates.
The full section coordinates are not asserted to be globally canonical.
The reported amplification 2.852 and condition number 2.929 depend on their
scaling and on the chosen similarity basis; they are not invariant physical
amplification factors. Boundedness of powers of a finite-dimensional exact
linear map is preserved under a fixed invertible coordinate change, whereas
these numerical constants are not.

An ordered real Schur decomposition selects a six-dimensional invariant
centre subspace U of B. The two excluded multipliers at sample 0 are about
7242.22684 and 0.00013807908. They are physical growing/decaying directions,
not coordinate symmetries. The centre map `C = U^T B U` contains four
neutral directions and an elliptic pair with trace 1.3418618576.

At every sample a real similarity W brings C numerically to four identity
directions and a planar rotation. The resulting positive metric
`G = W^(-T) W^(-1)` is approximately preserved. This is a numerical
tangent-space construction; no nonlinear symplectic reduction is claimed.

## Registered measurements

All 12 points use #319's authenticated arclength sample schedule. Worst
values over those points, unless otherwise stated:

| Quantity | Result |
|---|---:|
| Removed-subspace leakage, absolute operator norm | 6.78e-10 |
| Centre-subspace leakage | 1.01e-11 |
| Family tangent / neighbouring chord sine | 0.008008 |
| Neutral centre dimensions | 4 at every point |
| Similarity-basis condition number | 2.929 |
| Relative similarity defect | 6.98e-10 |
| Relative positive-metric invariance defect | 3.77e-10 |
| Largest Euclidean power norm, n=1,2,4,...,1024 | 2.852 |
| Sample-0 derivative refinement difference from archive | 6.49e-12 |
| Elliptic-trace refinement difference | 5.45e-12 |
| Power-envelope relative refinement difference | 6.50e-10 |
| Full-system quotient response error, epsilon=2e-6 | 1.59e-5 |
| Full-system quotient response error, epsilon=1e-6 | 3.18e-5 |
| Absolute initial/final constraint residual | 2.39e-13 |

The full-system checks cover all six independent centre quotient
directions at both epsilon values, with both signs. All four initial
constraints are solved together, including matter-current compensation.
Both clock returns are integrated consecutively, then the response is
projected into the normal quotient. The smaller epsilon gives a larger
finite-difference error, consistent with numerical subtraction sensitivity;
both pass the preregistered 1e-4 threshold. No extrapolated error bound is
claimed. Q5 contains one resolution comparison and one reproducibility
check: `(1024,2048)` differs from the archived derivative by 6.49e-12 relative;
`(2048,4096)` is #319's original default and reproduces it exactly in the
measurement environment, with zero matrix, trace and power differences.
The latter supplies no additional resolution beyond the archived baseline.

The finite power ladder and approximate invariant metric do not constitute
an infinite-time proof: residual floating-point errors can accumulate.
They support the registered numerical label in its conditional scope.

## Reproduction and evidence

The producer is `experiments/closure_ledger/r3_normal_quotient_probe.py`;
the portable replay is `experiments/closure_ledger/r3_normal_quotient_replay.py`.
The readable [archive](../experiments/closure_ledger/runs/20261003_r3_normal_quotient/)
contains the inherited sample matrices, new derivative refinements, all
signed initial/final 29-state perturbation endpoints, quotient bases and
metrics, power norms and decisions. The manifest SHA-256 is pinned in replay:

    132a99d9f051d7e96f96c5b8d8673d7e9deeffe045342ea1563c66ebe1a2f163

The raw evidence binds the published specification, measured source files
and #319's original archive hashes. Replay authenticates the inventory,
recomputes all six gates, reconstructs perturbation responses from the saved
states, and checks the results against the saved decisions. It does not
reintegrate the trajectories during routine replay.

The bound source inventory intentionally includes
`experiments/closure_ledger/r3_family_replay.py`. Consequently, even a later
change to that #319 replay makes this historical #320 provenance check fail.
Reproduce from the published experiment revision, or retain the original
source in a versioned compatibility layer for future development. Do not
silently refresh historical hashes. This documents the coupling without
changing the measured source inventory or archive.

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_normal_quotient_replay
    OPENBLAS_NUM_THREADS=1 pytest -q tests/test_r3_normal_quotient.py

Controls cover unremoved Jordan growth, hyperbolicity, noninvariant removed
directions, analytic rotation tangents and invariance under orthogonal
coordinate changes. Mutations reject changed matrices, omitted samples or
directions, nonfinite values, altered frames and endpoints, changed labels
and archive bytes.

## Post-measurement replay corrections

The measurement producer originally reconstructed its Schur frame during
replay. SciPy 1.17 and 1.18 can return different orthonormal bases for the
same repeated neutral eigenspace, so matching a recorded perturbation to a
newly chosen basis can fail despite identical subspaces. A separate portable
replay layer now uses the authenticated recorded frame after checking its
orthonormality and projector agreement with the freshly computed quotient
and centre. Every diagnostic and gate is recomputed; thresholds are unchanged.
An orthogonal-frame control exercises this equivalence explicitly.

The replay also enforces a strictly relative Q6 error by dividing by the
predicted response norm, without the producer's unit floor. The smallest
predicted column norm is 0.806. This raises the worst reported response error
from 2.84e-5 to 3.18e-5; every case still passes the frozen 1e-4 threshold.
The table above uses the stricter replay values. Original raw states and
historical producer diagnostics are preserved, and no acceptance threshold
or categorical decision changed.

The producer, measured dynamics, specification and archived results remain
unchanged. The portable replay reproduces the label with both NumPy 2.5.3 /
SciPy 1.18.1 and NumPy 2.3.5 / SciPy 1.17.0.
The combined normal-quotient and #319 suite passes 31 tests in the primary
environment; the 17 new tests also pass in the alternate environment.

## Remaining question

This closes the conditional linear-normal-response test proposed after
#319. A next experiment would need to establish nonlinear evolution in a
defined symmetry-reduced neighbourhood, including how the physical unstable
direction is treated. The family tangent has not been promoted to a proven
exact symmetry. Neither the present quotient nor its bounded response
supplies a mechanism that selects the measured family action.
