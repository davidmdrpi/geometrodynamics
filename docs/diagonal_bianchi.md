# Diagonal Bianchi IX: numerical crossing and two-return closure

Run date: 2026-09-29 UTC. Publication follow-up: 2026-09-30 UTC.
Prospective freeze: `0c51f9b5b528432c012fc15a4765af8829bb3867`, published in
PR #318 before implementation and new-amplitude measurements. The
[specification](diagonal_bianchi_prereg.md) and all historical labels are unchanged.

The exact diagonal Einstein–quartet subsystem passes both registered tests:

- **CROSSING_BRACKETED_NUMERICALLY**.
- **NONTRIVIAL_TWO_RETURN_ORBIT_VERIFIED_NUMERICALLY**.

Action selection remains **NOT_ESTABLISHED**. These are numerical statements,
not interval existence proofs, uniqueness claims, or quantum action levels.

## Crossing evidence

Nineteen circles were accepted on the frozen amplitude ladder, with no
failed attempts. The run stopped at its first bracket. The two endpoints
were re-solved on 95 nodes, in addition to the original 63-node solves:

| Endpoint amplitude | Canonical circle action | rho | omega - pi |
| --- | ---: | ---: | ---: |
| 0.07610925536 | 0.0169197308683 | 1.49848794405117 | -0.00950052772 |
| 0.09050966799 | 0.0240755208901 | 1.50404130601999 | +0.02539227461 |

The two grids agree in omega at the recorded precision. The frozen
uncertainty floor is 1e-9 radians; both endpoints exceed its tenfold
separation requirement. This is a resolved numerical bracket across
accepted circles, not a proof of a continuous invariant-circle foliation
through resonance. The direct shooting experiment supplies the separate
closure evidence.

Across all 19 circles and both refinements, the maximum grid residual is
3.94e-11 and the maximum half-offset residual is 3.96e-11, below the frozen
1e-9 limits. Four points on each circle were also integrated using the
unchanged full 29-state conformal-time equations with actual section events.
The maximum phase/full state-or-time discrepancy is 1.24e-12.

## Direct two-return shooting

The accepted circle closest to pi supplied eight equally spaced seed
angles. Each solve used two independent six-component section nodes and
solved P(z0)=z1, P(z1)=z0 without imposing an action or amplitude. All eight
candidates passed the registered nontriviality and closure checks.

| Diagnostic | Worst value across eight seeds |
| --- | ---: |
| Maximum component of shooting residual | 1.23e-13 |
| Full-state closure error, both DOP853 and Radau | 1.29e-11 |
| Absolute constraint residual, every saved validation state | 8.04e-13 |

The full-system checks evolve both returns consecutively, preserving the
first returned state as the second segment's initial state. Only the
elapsed-time coordinate is excluded from the closure norm. They stop at
the required section crossings rather than integrating beyond the return.
Both methods use rtol=2e-12, atol=2e-14 and max_step=.01. Positivity and
constraints are checked on 257 samples per return; these are sampled
checks, not rigorous bounds between samples.

For seed 0, the two section nodes `(A,p_A,x,p_x,y,p_y)` are approximately:

    z0 = (1.000503231205612, -0.000025712699647,
          0.077190977191945,  0.000802940697183,
          0.000364201695803, -0.254452104496236)
    z1 = (1.000436620468819,  0.000025726421217,
         -0.083326512003473, -0.001163768798227,
         -0.000312623403504,  0.211950511588411)

Use the archived full-precision nodes for reproduction. Its total conformal
period is 6.271615250392913 (DOP853) or 6.271615250392921 (Radau). The
one-return node separation is about 0.4933, so the solution is not the
background or a one-return fixed point.

Eight successful seeds do not establish eight distinct isolated orbits.
The experiment does not classify equivalences, degeneracy, stability or
isolation, and does not attach the bracketing circles' actions to the
period-two solutions as an exact orbit action.

## Post-publication review: numerical family

The [independent review of `03afa4e`](https://github.com/davidmdrpi/geometrodynamics/pull/318#issuecomment-5922620742)
on 2026-10-01 reproduced both registered labels, replay and all 12 diagonal
tests. Its separate event-located map using the unchanged 29-state equations
reproduced the bracket and closed all eight archived candidates to 8.5e-13
or better.

The review also reports exploratory evidence that the two-return solutions
are **numerically non-isolated**, forming a continuous one-parameter family:
opposite seeds exchange the two nodes; the shooting Jacobian has one
near-null singular value (4.3e-10, versus 0.139 for the next); and continuation
along that direction remains closed to 6e-13 or better over the traced arc.
Its finite-difference multipliers include a near-unit pair and the large
Einstein-static unstable pair. These are attributed review findings, not
new registered measurements or additions to the immutable archive. The
review's continuation data are not archived here, and exact non-isolation,
a global closed family and an exact unit Jordan block are not proved.

Carry this numerical family forward as prior information for the next
freeze. Neither the family nor the successful closures establish action
selection. In particular, a common action on one family would not by itself
establish a discrete spectrum or a mechanism selecting that family.

## Implementation, provenance and replay

`geometrodynamics/waves/diagonal_bianchi.py` evolves the exact scalar-phase
reduction with diagonal beta=x E0+y E1. The momentum constraints vanish
identically because both M and L are diagonal and the quartet has zero
current. Canonical section momenta are p_A=-6A' and (p_x,p_y)=A^2(x',y').

The batched DOP853 map uses rtol=2e-12, atol=2e-14 and phase max_step=.025.
Complex-step derivatives propagate through the same numerical flow. They
were checked against centred differences before the amplitude run. Their
tiny imaginary parts are not separately controlled by the absolute
solver tolerance; the fixed maximum step and independent map checks
therefore remain important. Toy controls cover a known circular twist
and a known nonzero period-two orbit. A toy-circle assertion was corrected
before measurement to use its solved action: that symmetric toy permits
the second mode's amplitude to change while the first harmonic is fixed.
No measured source changed during or after this run.

The immutable evidence is in
[`20260929_diagonal_bianchi`](../experiments/closure_ledger/runs/20260929_diagonal_bianchi/).
Each circle/refinement archive includes K, map images, off-grid images,
full validation histories and Newton traces. Each shooting archive includes
all trial-node traces, final nodes and both full-system integrator histories.
`provenance.json` binds the freeze and sources; `report.json` records the
recomputed decisions. `manifest.json` binds every evidence file with SHA-256.
Its pinned digest is:

    eebc287e12240f2d293f99a7206c10465a089825dd93b115a05c4bd1b6536dda

Replay authenticates the manifest, files and measured sources, reconstructs
actions and residuals, checks the continuation and seed schedules, and
requires exact categorical agreement. It does not trust saved acceptance
flags. Mutation controls reject missing curves, changed images, nonfinite
data, missing ladder points, changed histories and altered reports; a
background two-cycle is rejected as trivial.

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.diagonal_bianchi_replay
    OPENBLAS_NUM_THREADS=1 pytest -q tests/test_diagonal_bianchi.py

The recorded environment is Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1.
All 12 new tests pass there and with NumPy 2.3.5 / SciPy 1.17.0. The
combined diagonal and #316 phase-map suite passes 28 tests.
The producer refuses an existing output directory. A new run can use:

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.diagonal_bianchi_probe --output /tmp/diagonal-fresh

## Next candidate: family action and transverse stability

The next prospective specification should disclose the review's family
finding and test its extent and action rather than treating isolation as
an unexplored question. Continue with a fixed phase condition around the
putative full family, report closure residuals and their convergence, and
measure the canonical action along the family with a stated uncertainty.
Distinguish that section-family loop integral from the time integral along
one physical two-return history; the latter requires the appropriate full
canonical one-form. Do not substitute interpolation of the two bracket
actions for either measurement.

Use variational equations for the constraint-compatible full-period
monodromy and closure Jacobian, with convergence and direct-perturbation
checks. Distinguish family, gauge and symmetry directions from physical
growth and the known homogeneous background instability. Report the full
unstable spectrum separately from stability within a specified centre
sector. General tensor perturbations require the full 12-dimensional
section map, including off-diagonal directions; angular-momentum-carrying
perturbations require the second-order matter compensation described in
the #317 extension specification. A diagonal spectrum alone cannot settle
that broader stability question.

That test can establish or reject robustness within a stated sector. Even
an isolated stable orbit would still require a physical selection mechanism
before an action-quantization claim. Inhomogeneous perturbations remain a
separate open test. No stability result is asserted here.
