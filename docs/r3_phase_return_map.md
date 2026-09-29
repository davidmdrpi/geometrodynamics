# Constrained R3 phase return map: measured local shift

Date: 2026-09-29. PR #316. The prospective specification is
[r3_return_map_prereg.md](r3_return_map_prereg.md), published in
`259ab9c08b774af5e76819789a40a379172a318d` before this implementation and
measurement. Its text, thresholds and prior R3 evidence remain unchanged.

The registered outcome is **SHIFT_AWAY_FROM_TARGET**, with all G1–G5 passing:

    rho0 = 1.484666408416047
    nu   = -0.9501180968942993  turns per unit canonical action

This is the leading coefficient in the local finite-order normal form
`rho(I) = rho0 + nu I + ...`, with `I = |(Q-iP)/sqrt(2)|^2` per unit S3
coordinate volume. Since rho0 is below 3/2, the negative coefficient points
away from that target. Following the frozen decision rule, this stops the
local R3 search. It does not exclude a later turn at larger amplitude,
other sectors, or R2. Resonance crossing and full-state closure are
NOT_TESTED; dynamical action selection is NOT_ESTABLISHED.

## Calculation and evidence

`geometrodynamics/waves/r3_phase_return_map.py` implements the exact
constraint reduction and scalar-phase clock from the freeze. The section
coordinates are `(A-1, -6 A', x, A^2 x')`. The phase interval is pi/2 to
5pi/2; the conformal return time is integrated as a fifth variable. No
trajectory resets, kicks or constraint projections are used.

`taylor_jets.py` carries all ordinary monomial coefficients through total
degree three in four initial coordinates. The return map includes the
hyperbolic variables. The reduction solves the quadratic and cubic centre
graph, includes its induced area form and cubic Darboux correction, and
removes quadratic centre terms before extracting the resonant coefficient.
The elliptic canonical basis is approximately diag(0.5899638217, 1.6950191915),
with determinant one; the action scale is fixed by the symplectic form.

The immutable run is in
[`20260929_r3_phase_return_map`](../experiments/closure_ledger/runs/20260929_r3_phase_return_map/).
`raw.json.gz.b64` contains both complete cubic maps, both return-time jets,
and the 24 full-system validation histories (257 states each), initial and
returned states, section times and phase-formulation endpoints.
`report.json` contains the derived maps, graph, Darboux and normal-form
coefficients, diagnostics, every remainder and the decision.

| Check | Measured value |
| --- | ---: |
| DOP853 nu | -0.9501180968942993 |
| RK45 nu | -0.9501180968937700 |
| Absolute nu difference | 5.29e-13 |
| Maximum difference over all map coefficients | 1.80e-8 |
| Maximum linear tensor-trace discrepancy | 5.08e-13 |
| Maximum linear symplectic residual | 1.83e-12 |
| Polynomial symplectic residual, absolute / scaled | 3.73e-9 / 1.30e-17 |
| Maximum homological condition number | 1720.38 |
| Maximum centre-graph algebraic residual | 7.22e-15 |
| Maximum absolute Re(a21/lambda) | 1.61e-11 |
| Maximum full/phase state-or-time discrepancy | 6.59e-13 |
| Maximum sampled normalized constraint residual | 9.29e-15 |
| Maximum sampled absolute constraint residual | 5.58e-14 |
| Minimum sampled chart diagnostic | 0.8750005261 |
| Largest finest-radius map / graph error | 1.11e-9 / 7.11e-11 |

The archive's `constraint_max` is the inherited normalized residual. The
absolute residual above was additionally recomputed from all archived
states during review; both are below the frozen 1e-9 requirement. Chart
and constraint checks on these histories are sampled, not interval bounds.
The finite-radius experiment checks the cubic map and centre graph; it is
not an independent measurement of nonlinear frequency on invariant circles.

## Review clarification: adaptive error control

The solver state has 5 times 35 = **175 scalar components**. Every jet
coefficient, including degree three and the return-time coefficients,
participates in SciPy's adaptive error calculation. For each component j,
the scale is `atol + rtol max(abs(y_j), abs(y_new_j))`. RK45 takes the RMS
of its scaled embedded error. DOP853 combines the squared norms of its
scaled fifth- and third-order estimators with the same dimension
normalization (SciPy `_ivp/rk.py`, `_estimate_error_norm`). This is aggregate
error control, not a strict maximum-component or global error bound.

The frozen tolerances are rtol=2e-12, atol=2e-14. Phase max_step is 0.025
for DOP853 and 0.0125 for RK45. The respective runs used 3038 and 5786 RHS
evaluations. The full-map coefficient discrepancy is reported separately
from the much smaller difference in nu; they are not interchangeable
accuracy estimates. The independent 29-state event-return comparisons use
the original conformal-time equations with max_step=0.01.

## Review clarification: reporting-floor use

All 16 map-error halving orders lie between 3.99581 and 4.00361. Two of
their finer errors are below 1e-10, but neither needs the floor exemption
to pass the 3.5 order requirement.

Of the 16 graph-error halvings, ten have finer error below 1e-10. Exactly
one requires the exemption: at angle pi, radius 0.008 to 0.004, the error
falls from 2.43850e-10 to 2.24012e-11, giving order **3.44435**. Its next
halving gives order 3.96551 and error 1.43394e-12. The other 15 graph
orders exceed 3.5. All eight finest graph errors are below the reporting
floor. Thus the graph check includes one explicitly floor-qualified pass;
it does not demonstrate fourth-order scaling above the floor at every
angle. No radius, floor or gate was changed after measurement.

## Independent PR #317 and dated namespace note

The [review of the specification](https://github.com/davidmdrpi/geometrodynamics/pull/316#issuecomment-5883261108)
reports compatible coefficients from #317's event-return jets, invariant
circles, and a separately implemented phase-clock check. These are external
cross-checks, not inputs to this run or its G1–G5 decision. #317's registered
UNRESOLVED result under its own failed C6 gate remains unchanged. Earlier
R3 UNRESOLVED evidence likewise remains unchanged. This implementation and
measurement completed locally before that review was retrieved; only the
specification was public when the review was posted.

On 2026-09-29, after measurement, the unpublished module, probe, replay,
test and run-directory names were changed from `r3_return_map` to
`r3_phase_return_map` to coexist with #317. The engine bytes and both
archived artifacts are unchanged. In the producer, only its module import,
source-path entries and run-directory path changed. The freeze's original
path and bytes are unchanged.

The report deliberately retains its original producer paths and hashes.
Replay maps the two historical engine/producer paths to their new paths
and reverses the producer's namespace substitution before hashing. This
recovers the exact measured source bytes; it does not regenerate historical
hashes or silently rebind the evidence to modified physics. New measurements
use the current namespace and a fresh output directory.

Pinned SHA-256 fingerprints:

    raw.json.gz.b64
    016877441c066b4c491c89c88b4ac54507439284049ffbe4958668232fa20587
    report.json
    0154641c9e6476f5d9e4954bf2c2c860e2c6fd394b7456101af7849c5ab087da

## Reproduction and verification

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_phase_return_map_replay
OPENBLAS_NUM_THREADS=1 pytest -q tests/test_r3_phase_return_map.py
# New measurement; existing evidence cannot be overwritten:
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_phase_return_map_probe --output-dir /tmp/r3-phase-fresh
```

The recorded environment is Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1.
All 16 new tests pass there and with NumPy 2.3.5 / SciPy 1.17.0. The
combined R3 historical, review-control and phase-map suite passes 41 tests.
Controls cover known positive, negative and zero twists, nontrivial
canonical centre embeddings, Taylor remainder scaling, the inherited
full vector field, missing/nonfinite/altered evidence and changed verdicts.

Replay authenticates both artifacts and measured sources, rebuilds the
reductions and every gate from raw evidence, and requires exact categorical
agreement. It compares remainder errors at 1e-10 and other derived numbers
at 1e-8 using the inherited scaled comparison. Ratios of near-floor errors
are numerically ill-conditioned, so the original order values are
authenticated by the report fingerprint while their pass/fail decisions
are recomputed exactly; replay does not require roundoff-identical ratios.
