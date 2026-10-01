# Two-return family: action and transverse stability (results)

Date: 2026-10-01.
- Specification and code: [`1fcccf9`](r3_family_prereg.md), pushed 01:37 UTC, before any loop tracing or any multiplier at a family point.
- Archives: `experiments/closure_ledger/runs/20261001_r3_family/` (`stage_F.json`, `stage_S.json`, `result.json`), each binding the five source hashes.
- Scoring: `experiments/closure_ledger/r3_family_score.py` (see §5).

## 1. Registered labels

| label | result |
|---|---|
| FAMILY | **CLOSED_FAMILY_LOOP_NUMERICALLY** |
| loop action | **I = 0.0188474 ± 1.6e-9** |
| ACTION_CONSISTENCY | **CONSISTENT**: circle interpolation gives 0.0188475 (quadratic), 0.0188681 (linear) |
| DIAGONAL_TRANSVERSE | **ELLIPTIC** |
| OFF_DIAGONAL | **MARGINAL**, explained by symmetry in §3 |

Every check F1–F5 and S1–S4 passes:

| id | check | value |
|---|---|---|
| F1 | loop closure residual | <= 1.0e-12 over 111 points |
| F2 | landing distance on return to start | 6.0e-13, after 110 steps |
| F3 | family-tangent defect | <= 5.2e-11 |
| F4 | Radau two-return closure | <= 1.0e-10 |
| F5 | action, all points vs half the points | 1.6e-9 |
| S1, S2 | block decoupling, determinant | pass |
| S3 | direct nonlinear perturbation vs M2 | <= 6.1e-10 |
| S4 | jet closure | <= 2.0e-12 |

The exact two-node Jacobian has one null singular value (3e-13 to 7e-12) at every sample.

## 2. The family is one closed loop at the resonant action

The two-return solutions trace a single closed loop. Its action, 0.0188474, agrees to 1.6e-7 with the action at which the invariant-circle family crosses omega = π, from quadratic interpolation through a = .0640, .0761 and .0905. The loop is therefore the resonant invariant circle itself. It survives as a closed curve of period-two points, unbroken at the 1e-12 level.

In the diagonal subsystem, the 3/2 closure condition fixes **one action**, I* = 0.0188474. It leaves the phase along the loop free, and also the orientation of the principal axes (§3). This is a numerical statement in units of action per unit S3 coordinate volume (κ = a = 1). It is not an action quantum, and it is not evidence of a physical selection mechanism.

**Constant multipliers (descriptive, not registered).** The multipliers are identical at all 12 loop points to about 1e-10:
- diagonal transverse trace 1.3418618576;
- Einstein-static multiplier 7242.2269805.

Distinct periodic orbits normally have different multipliers. Equal multipliers along a whole family suggest a continuous symmetry or an additional integral of the reduced dynamics, at least on this torus. The diagonal Bianchi IX potential has only the discrete D3 symmetry, and I have not identified the mechanism. It is an open question, and it bears on whether the family is exact or only extremely weakly broken.

## 3. Transverse stability in the homogeneous sector

At each of the 12 samples the multipliers of P² split into three groups.

| directions | multipliers | meaning |
|---|---|---|
| Einstein-static (A, p_A) | 7242.2 and 1/7242.2 | the homogeneous instability, shared with the ESU itself (about 85 per return) |
| diagonal family pair | 1 (Jordan) | motion along the loop |
| diagonal transverse | e^{±0.837i} (trace 1.34186) | **elliptic: linearly stable** |
| off-diagonal (3 pairs) | all 1 | SO(3) rotations and conserved angular momentum (below) |

**Why the off-diagonal multipliers are all 1.** At a diagonal, anisotropic configuration, off-diagonal shape perturbations are rotations of the principal axes together with their conjugate angular momenta. The check:
- (M2_off - I) has rank 3;
- its null space is exactly the three rotation generators, fixed to 2e-14 to 1e-10;
- the angular-momentum directions map into rotations.

So an off-diagonal perturbation, with its angular momentum compensated by a quartet current, makes the axes precess at a constant rate. That is secular, linear drift, not exponential growth. The registered label MARGINAL is retained. The specification described these as "off-diagonal tensor pairs" and did not anticipate that, at an anisotropic base point, they are symmetry directions. They have no oscillation of their own to classify.

**Answer to "does it survive general homogeneous perturbations?"**
- No exponential instability appears beyond the Einstein-static direction, which every homogeneous orbit has, including the ESU.
- The diagonal transverse pair is elliptic.
- The off-diagonal directions are neutral by symmetry.
- Within the homogeneous centre directions the family is therefore linearly stable, with neutral drift.
- Nonlinear stability is NOT_TESTED.
- Inhomogeneous perturbations are NOT_TESTED.
- The Einstein-static instability remains a physical caveat for any claim built on this orbit, as it does for the static universe.

## 4. What this establishes and what it does not

Within the exact diagonal homogeneous system:
- **Established numerically:** the refocusing closure at 3/2 occurs on one closed loop, at one action, I* = 0.0188474; it is linearly stable in the diagonal transverse direction and neutral in the rotation directions.
- **Not established:**
  - a physical reason why histories must close;
  - that closure at other rationals does not select other actions;
  - stability against inhomogeneous perturbations;
  - any connection to an action quantum.

## 5. Disclosures

- **Scoring wrapper.** The frozen probe's `score` stage crashed after computing, while printing, because numpy bools are not JSON-serialisable. Editing the probe would change its bound source hash. Instead, `r3_family_score.py` calls the unchanged `score()` and converts numpy scalars. It changes no number or rule.
- **Interpretation added after the result.** The rotation-generator analysis in §3 and the constant-multiplier observation in §2 were computed after the labels were fixed. They are interpretation, not registered tests.
- **Prior information.** It is listed in §6 of the specification.
