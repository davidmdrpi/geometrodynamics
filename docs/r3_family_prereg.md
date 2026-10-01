# Specification: action and transverse stability of the two-return family

Date: 2026-10-01. Parent: `claude/geometrodynamics-qft-audit-vpktax` at `b1b265b` (PR #317).
Publish before tracing the family loop or computing any multiplier at a
family point. The code is committed with this document. Corrections are
dated notes only.

## 1. Questions

#318 verified two-return closure for the diagonal-circular mode. Reviewing
it, I found that the closures are not isolated: they lie on a
one-parameter family (#318 review comment). Two questions follow.

1. **Family action.** Does the family close into a loop, and if so, at
   what action? A single closed loop has one action. If it is the resonant
   invariant circle, closure at 3/2 fixes the action and leaves the phase
   free. The test also checks whether the loop's action agrees with the
   action at which the circle family crosses ω = π.
2. **Transverse stability.** Within the homogeneous sector, are the family
   orbits linearly stable in the directions transverse to the family? The
   directions are the diagonal transverse pair and the three off-diagonal
   tensor pairs. The Einstein-static pair is hyperbolic for every
   homogeneous orbit, including the ESU, and is reported separately.

**Scope.** Homogeneous perturbations only, with the quartet along
(q,0,0,0). Inhomogeneous perturbations are NOT_TESTED. No result here
establishes action selection or quantisation.

## 2. Maps and linearisation (`geometrodynamics/waves/r3_family.py`)

**Section map.** The 12-coordinate section map
z = (A, p_A, x_1, p_1, ..., x_5, p_5) is evaluated with the unchanged
`conformal_rhs`: DOP853, 1e-12/1e-14, event q = 0 with q' < 0.
- The shape and its rate map exactly through the eigenbasis: matrix
  exp/log and their divided-difference Fréchet derivatives.
- The diagonal subsystem is z[:6] (E0, E1), identical to #318's.
- Pre-freeze checks: the map agrees with my earlier diagonal event map to
  7e-18, and the round trip is exact to 4e-16.

**Exact linearisation.** Order-1 jets of the 22-component reduced state
(`r3_extension.rhs`, which equals `conformal_rhs` for q = (q,0,0,0) to
1e-15).
- Fixed-step RK4 with N = 2048 and 4096, Richardson-extrapolated.
- The return time is solved as a jet, so the section-time correction is exact.
- Section-coordinate maps are differentiated by centred differences, h = 1e-6.

Pre-freeze validation of the linearisation:
- **At the ESU:** multipliers 85.019695 and .011762, and all five tensor
  blocks match #310's trace to 1.7e-10.
- **At an arbitrary point off the family:** it agrees with centred
  finite differences of the event map to 2.1e-9 (12D) and 2.9e-9 (6D).

**Matter compensation.** Off-diagonal perturbations carry angular momentum
at first order. The momentum constraint is met by quartet components dq_i
of the same order. In the (A, q, M, L) equations these enter only through
|q|^2, q·q' and |q'|^2, which are second order. So the tensor
linearisation computed with q = (q,0,0,0) is the physical one. The dq_i
themselves obey the clock's Hill equation. They are the neutral SO(4)
symmetry directions (multiplier 1) and are not stability modes.

## 3. Stage F: the family loop

- **Start.** Newton from #318's published seed-0 nodes, used as an initial
  guess only.
- **Continuation.** Two-node closure F = [P(z0) - z1, P(z1) - z0] in the
  diagonal 6D subsystem. Pseudo-arclength continuation along the null
  tangent, with step .02 in the 12-vector (z0, z1). Finite-difference
  Jacobians are used for Newton. Corrector tolerance: max residual 1e-12.
  At most 400 steps.
- **Stop.** After travelling more than 10 steps, stop when within 1.5
  steps of the start. Then land on the start along its own tangent row.
- **Loop action.** I_loop = (1/2π)|∮ p_A dA + p_x dx + p_y dy| over the
  ordered z0 points, using a periodic cubic spline in chord length. The
  uncertainty is the difference between all points and every second point.
- **Circle interpolation.** Circles at a = .0640, .0761 and .0905 on 63
  nodes, continued from a = .004 with the same solver as #318's. I*_circ is
  the action at ω = π by quadratic interpolation of I(ω), with the linear
  fit over the two nearest circles for comparison.

## 4. Stage S: stability at 12 loop points

At 12 points equally spaced in arclength:
- the exact-jet DP for both nodes, in all 12 coordinates, giving M2 = DP(z1) DP(z0);
- the diagonal 6×6 block, classified into three reciprocal pairs:
  - the Einstein-static pair (largest |trace|);
  - the family (unit) pair (trace closest to 2);
  - the diagonal transverse pair;
- the off-diagonal 6×6 block, classified into three pairs, each ELLIPTIC
  (|tr| < 2 - 1e-6), HYPERBOLIC (|tr| > 2 + 1e-6) or MARGINAL;
- the family-tangent defect |(M2_diag - I) t|, with t the null vector of the
  exact two-node Jacobian, plus that Jacobian's smallest two singular values;
- Radau closure over two returns.

At the first sample, each off-diagonal eigen-direction is also checked
against a direct nonlinear perturbation: P∘P at ±1e-6.

## 5. Checks and labels

| id | requirement |
|---|---|
| F1 | every loop point has closure residual <= 1e-10 |
| F2 | the loop returned to the start, landing distance <= 1e-8 |
| F3 | family-tangent defect <= 1e-6 at all 12 samples |
| F4 | Radau two-return closure <= 1e-9 at all 12 samples |
| F5 | \|I(all) - I(half)\| <= 1e-6 I |
| S1 | diagonal/off-diagonal coupling blocks of M2 <= 1e-8 × max entry (exact by reflection symmetry) |
| S2 | \|det M2 - 1\| <= 1e-6 |
| S3 | direct perturbation agrees with M2 to 1e-3 (relative) |
| S4 | exact-jet two-return closure <= 1e-9 at all samples |

**Family labels.**
- **CLOSED_FAMILY_LOOP_NUMERICALLY:** F1–F5 all pass.
- **FAMILY_UNRESOLVED:** otherwise.
- **ACTION_CONSISTENCY** (descriptive): CONSISTENT if
  |I_loop - I*_circ| <= 3(|quadratic - linear| + action uncertainty);
  otherwise INCONSISTENT.

**Stability labels** (require S1–S4):
- DIAGONAL_TRANSVERSE and OFF_DIAGONAL: **ELLIPTIC / ALL_ELLIPTIC**,
  **HYPERBOLIC_PRESENT** or **MARGINAL**;
- **UNRESOLVED** if any S-check fails.

**Interpretation.** CLOSED_FAMILY_LOOP with ALL_ELLIPTIC and ELLIPTIC would
mean the family is a single closed loop at one action, linearly stable in
every homogeneous direction except the Einstein-static one shared with the
ESU. It would still not establish a selection mechanism, nonlinear
stability, or inhomogeneous stability.

## 6. Prior information (disclosed)

From my #318 review, all exploratory:
- the family exists over ±76° of phase;
- eight seeds give four distinct orbits;
- finite-difference P^2 multipliers ≈ e^{±0.835i} plus a pair at 1 ± 8e-6;
- circles at a = .0761 and .0905 with I = .01692 and .02408 and
  omega - π = -.0095 and +.0254;
- interpolated I* ≈ .019.

Pre-freeze method checks are in `experiments/closure_ledger/r3_family_prefreeze.py`
and `tests/test_r3_family.py`. No loop tracing and no multiplier at a
family point beyond those exploratory ones has been computed.
