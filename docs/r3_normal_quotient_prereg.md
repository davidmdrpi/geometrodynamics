# Prospective test: normal response of the two-return family

Date: 2026-10-03 UTC. Parent: merged #319, `3304cc5c10ccebf549f7dee734724aef9e3c25e8`.
Publish this specification before evaluating any new quotient or perturbation
measurements. The implementation and evidence will bind this file's SHA-256.

## Question and prior information

#319 found a period-two family and unit Jordan blocks, with unreduced secular
drift. Its archived 12 matrices, loop nodes, elliptic transverse pair, large
Einstein-static pair and review validation are prior information. No blind
discovery claim is made. The new question is whether the linear normal
response is numerically bounded after quotienting the three rotation tangents
and one family tangent, conditional on excluding the hyperbolic pair.

This is a tangent-space quotient with an explicitly stated norm, not a
nonlinear symplectic reduction or proof of an exact continuous symmetry
generating the diagonal family. It retains angular-momentum perturbations.
It tests distance to the whole family, not stability of one periodic orbit.

## Construction

Use the unchanged #319 12-coordinate section and all 12 registered samples
in their original order. Authenticate original archives with #319's pinned
hashes, verify their source hashes and reproduce its replay before proceeding.
Use the Euclidean norm in the existing dimensionless section coordinates.

For each sample with two-return derivative M:
1. Form three analytic rotation generators from commutators [Omega,beta]
   and [Omega,p], projected on the five STF basis tensors. At the section
   p=A^2 beta'; SO(3) acts on both tensors by conjugation.
2. Extract the diagonal family tangent from the smallest right singular
   vector of M_diag-I. Check its null defect and alignment (up to sign) with
   the neighbouring archived loop-node chord. It is a numerical tangent,
   not an independently proved exact symmetry.
3. Orthonormalize those four columns into R; let Q span their orthogonal
   complement. Check invariance before forming B=Q^T M Q (dimension 8).
   The quotient norm is min_a ||delta z-R a||_2 = ||Q^T delta z||_2.
4. Use ordered real Schur decomposition of B, selecting eigenvalues with
   0.1 < |lambda| < 10. Require six selected directions and one excluded
   stable/unstable pair. Let U be the orthonormal invariant centre basis
   and C=U^T B U. Report the excluded pair, never label it gauge.
5. Require four near-unit semisimple directions and one elliptic pair in C.
   Build a real similarity W from four null vectors of C-I and the real and
   imaginary parts of an elliptic eigenvector. Compare C W with W D, where
   D is identity on the four neutral directions and a unit planar rotation
   on the remaining pair. Record cond(W), the positive metric
   G=W^(-T)W^(-1), and its relative invariance defect ||C^T G C-G||/||G||.
   Evaluate ||C^n||_2 at n=1,2,4,...,1024. These finite, floating-point checks
   do not prove boundedness for infinite time.

## Independent checks

At sample 0, repeat the derivative with #319 DP at step pairs (1024,2048)
and (2048,4096); compare with the archived matrix and repeat reduction.
Compare elliptic trace, power envelope and positive-metric defect.

At sample 0 lift each of the six centre quotient basis vectors by Q U.
For epsilon=2e-6 and 1e-6, integrate plus/minus perturbations for two actual
clock returns with unchanged full conformal_rhs (DOP853, rtol=1e-12,
atol=1e-14, max_step=.025). Solve all four initial constraints jointly for
q_1,q_2,q_3,q'_0, with q'_1..3=0, starting from first-order momentum
compensation. This implements a section of the constraint surface, not
an assumption that matter currents vanish. Use event states consecutively.
Compare Q^T times the central-difference response with B U, including any
leakage outside the centre space. Save both initial and final full states.

## Gates and labels (fixed before evaluation)

All numbers must be finite; require exact sample count/indices, expected
matrix shapes and matching source/archive hashes. Failures yield UNRESOLVED.

| Gate | Requirement |
|---|---|
| Q1 | rotation rank 3 (smallest singular value >1e-6); family tangent residual <=1e-7 and chord sine <=1e-2; four independent removed directions |
| Q2 | absolute removed-subspace leakage ||Q^T M R||_2 <=1e-7; centre leakage <=1e-7 |
| Q3 | 8 quotient and 6 centre dimensions; excluded moduli <.1 and >10; four singular values of C-I <=1e-7 and other two >1e-3 |
| Q4 | elliptic modulus error <=1e-7; cond(W)<=1e4; relative similarity and metric defects <=1e-7; all finite power norms <=1.01 cond(W) |
| Q5 | sample-0 derivative relative differences <=1e-7; elliptic trace differences <=1e-6; power-envelope relative differences <=1e-4; Q1-Q4 also hold at both new resolutions |
| Q6 | all 12 signed finite differences (6 directions x 2 epsilon values) relative quotient errors <=1e-4; all initial/final absolute constraint residuals <=1e-8 |

If all gates pass: BOUNDED_CENTER_QUOTIENT_NUMERICALLY. Otherwise:
NORMAL_RESPONSE_UNRESOLVED, with each failure recorded. No parameter tuning
or new branch search after observing results. Numerical runtime/serialization
fixes must be disclosed without changing these criteria.

Always report UNREDUCED_HYPERBOLIC_PAIR_PRESENT when observed; nonlinear and
inhomogeneous stability, exact symmetry, and physical action selection remain
NOT_ESTABLISHED. Four neutral quotient directions must be retained rather
than discarded as gauge. A positive result is conditional linear normal
stability evidence in the specified norm and centre sector only.

## Evidence and controls

Archive nodes, matrices, bases, singular values, gates, power norms, refined
derivatives and signed full-system perturbation endpoints in readable JSON.
Bind all producer/dynamics/specification sources and original archives with
SHA-256. Replay authenticates, recomputes quotient diagnostics and gates,
and reconstructs perturbation checks from saved states; optional live checks
may rerun the integrations. Preserve historical sources and labels.

Controls: a toy unit Jordan block loses secular growth only when its genuine
invariant tangent is removed; an additional unremoved Jordan block and an
extra hyperbolic centre mode must fail the bounded-response test. Rotating
orthonormal quotient coordinates must preserve eigenvalues and power norms.
Reject changed matrices, omitted samples/directions, NaNs, modified endpoint
states, and altered labels. Never use saved acceptance booleans as evidence.
