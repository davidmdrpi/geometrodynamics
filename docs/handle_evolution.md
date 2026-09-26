# Evolved neck crossings: continuous momentum transfer

The registered short-time evolution passes its validation and gives two
affirmative results: `EVOLVED_TEST_WORLDLINE_CROSSING` and
`COVARIANT_MATTER_MOMENTUM_BALANCE`. `FINITE_CROSSING_IMPULSE` and
`DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE` retain their historical Boolean
value false. The impulse is a measured negative; reciprocal exchange is
**NOT_TESTED**, because this experiment has no independent gravitating
receiver. Crossing a smooth evolved neck did not produce a finite momentum
jump in these preparations. The [subsequent saved-field audit](handle_surface_audit.md)
finds future-trapped neck spheres at t=.01 in both metric frames and
distinguishes cut-face flux cancellation from two-object recoil.

The full specification was public in
[PR #308](https://github.com/davidmdrpi/geometrodynamics/pull/308) at
`cdb045eddfa0f01bf560559c5ea71b90887e3d0f` before implementation and measurement.
The prior initial-data extension remains 8/8; the original 6/8 and
reconstruction 7/8 remain unchanged.

## Coupled evolution and independent checks

The Einstein-frame metric is -dt^2+A^2 ds^2+B^2 dOmega^2, with evolved
extrinsic-curvature eigenvalues k,l. The four scalars are represented as
(u,v n); both independent amplitudes and their momenta evolve. The norm
is not fixed, and the field stress backreacts on both metric functions.
All equations and normalizations are in
[the preregistration](localized_mouth_evolution_prereg.md).

Metric compatibility of the full four-dimensional target-space Christoffel
connection and equivariant angular closure were checked symbolically.
The evolution agrees with exact round Jordan-universe breathing data to
9.55e-14, and spatial Ricci/divergence agree with an independent coordinate
engine to 8.89e-16. The analytic initial Einstein-time norm-acceleration
prediction is
-3.958973274443148; it is not a separate numerical measurement. The
round-breathing right-hand-side check above validates the active norm
equation that fixed-norm evolution would discard.

Both archived eta=0 and eta=.3 data were evolved with N=512,1024,2048 and
200,400,800 RK4 steps, respectively, over Einstein proper time [0,.01].
There was no filtering, constraint projection, added damping or resolution
search. The smallest metric factors stayed above .029 and f stayed >=.875.

| Quantity | eta=0 | eta=.3 |
|---|---:|---:|
| Finest maximum normalized Hamiltonian constraint | 2.45e-9 | 2.08e-9 |
| Finest maximum normalized momentum constraint | 1.62e-8 | 1.70e-8 |
| Finest maximum absolute Hamiltonian constraint | 2.14e-6 | 2.27e-6 |
| Finest maximum absolute momentum constraint | 2.39e-8 | 2.22e-8 |
| Medium/fine field difference | 6.08e-9 | 6.08e-9 |
| Coarse/medium versus medium/fine field-difference ratio | 14.87 | 14.87 |
| Matter momentum balance defect | 5.80e-16 | 5.80e-16 |
| Defect if geometric source is omitted | 2.11e-7 | 2.11e-7 |

Constraints decrease with refinement. Evolution uses the explicitly
registered Ricci-component normalization, which differs from #307's
physical initial-data normalization; the smaller normalized numbers are
not a claim that the archived initial solution became more accurate.
The radius-ratio differences are already below the 1e-8 convergence trigger.
The registered field comparison covers A, B, u, v; it does not impose the
same bound on k, l, P, Q. In particular it must not be quoted as convergence
of all eight state components.

## Background and neck are distinct on this interval

The exact homogeneous Jordan scale stays one. The two constraint-satisfying
perturbed controls start at 1+/-1e-4 and end at 1.00010152323 and
.999898476573. Their growing mode is retained; their Friedmann defects
stay below 3.34e-16. Einstein-frame breathing is explicitly distinguished
from Jordan-frame expansion.

For eta=.3 the Jordan neck radius decreases **2.77564%**, while the central
bulk radius decreases **0.000162917%**. The zero-momentum case is nearly
identical: **2.77640%** and **0.000162986%**. The neck remains a local
minimum of the Jordan areal radius. The bulk field norm changes by
-2.02285e-4. No homogeneous mode was subtracted from the evolved equations.

This separates the measured local neck response from the small central
bulk drift over the fixed short interval. It does not establish long-time
stability, a stationary throat, or momentum-driven throat formation.

## Crossings and the shrinking-window test

Three Jordan-frame timelike geodesics start at s=L-.05 with local speeds
.25,.5,.75; their reflected partners start at -L+.05. All six cross the
neck in both preparations. For eta=.3 the positive-side crossing times
are .00610591821, .00302184367 and .00201081605. Reflected trajectories
agree to 7.55e-15, and the medium/fine crossing comparisons pass.
Coordinates are unwrapped: the antipodal scalar bundle transition never
inserts a radial momentum sign flip.

These probes start close to the neck; they are not a demonstration of
travel from the round bulk through a globally traversable wormhole. Their
stress is neglected, so they cannot serve as independent recoil receivers.

For eta=.3, the signed local-frame momentum changes in symmetric windows are:

| Initial local speed | Half-width .000625 | Half-width .0003125 | Half-width .00015625 |
|---|---:|---:|---:|
| .25 | -.00211871159 | -.00105955078 | -.000529799771 |
| .5 | -.00238632353 | -.00119340744 | -.000596734440 |
| .75 | -.00312797953 | -.00156435314 | -.000782222009 |

Successive absolute ratios are approximately **2**, not the registered
finite-impulse window [.8,1.2]. The changes tend to zero with window width.
At the symmetric seam,
A_s=f_s=0 and dp_s/dt=0. The leading local-frame change is
`d p_hat/dt = p_hat (k + f_t/(2f))`: gravitational redshift due to the
changing radial frame, not an independently measured recoil receiver.
The signed covariant radial momentum changes shrink approximately **eightfold**
per halving, consistent with a cubic window dependence at a reflection-symmetric
neck where the instantaneous radial force vanishes. Neither momentum
definition yields a finite impulse candidate. Reflected probes have opposite
signed changes, and eta=0 gives the same continuous behavior.

This is also the local regularity prediction. The Hamiltonian geodesic ODE
has a bounded smooth force while A and f remain positive, so a momentum
change on a window of width 2delta is bounded by 2delta times that force's
supremum and vanishes as delta tends to zero. Smooth topology transport
does not supply an extra delta-function force. This tests the instantaneous
impulse interpretation; it does not exclude all possible finite-duration,
discretely selected exchanges in a different coupled preparation.

## Matter balance and the missing receiver

The one-sided tube [L-.2,L] has a nontrivial covariant matter momentum
balance. Its final momentum is -4.18972e-9 for eta=.3 and -3.96034e-9 for
eta=0. Flux and the geometric source were independently integrated through
the RK4 stages; comparison with momentum reconstructed from evolved fields
gives the defects above. Omitting the geometric source fails the negative
control by many orders of magnitude.

The radial vector X=partial_s is not a Killing vector. This is a specified
matter-current balance with gravitational exchange, not a total conserved
ADM momentum, an invariant particle four-momentum at infinity, or a
mouth/absorber recoil ledger. The code reports discrete reciprocal exchange
as unestablished because no independent gravitating receiver is present,
in addition to the failed finite-impulse test. It does not conceal that
missing ingredient behind the two affirmative verdicts.
The inherited preparation also contains no separately injected localized
incident wave packet; this is the evolution of the certified quartet data.

The next substantive requirement is a constraint-completed, dynamically
localized receiver and a reciprocal transfer observable, with finite-duration
event/action selection tested independently. Simply detecting this crossing
or integrating a force over a chosen window would not derive a quantum.

## Evidence and reproduction

Raw final fields, all 201 diagnostic times, worldlines, signed windows,
convergence comparisons, input/source hashes and separate verdicts are in
[`runs/20260923_handle_evolution/`](../experiments/closure_ledger/runs/20260923_handle_evolution/).
Raw JSON SHA256:
`8a858b6d7ecce0de3047f92a8baa77fd9c262ec9475d6886d539eac26cadf815`.

```sh
python -m experiments.closure_ledger.restore_handle_evidence
python -m experiments.closure_ledger.handle_evolution_probe \
  --rescore experiments/closure_ledger/runs/20260923_handle_evolution/evolution.json \
  --output-dir /tmp/handle-evolution-replay
pytest -q tests/test_handle_evolution.py
```

Replay recomputes the prior 8/8 certificate, equation checks, background
controls and all six evolutions. It compares recorded evidence with scaled
1e-11 tolerance, far below the scientific acceptance bounds. Changed sources,
fields, momenta, missing records and nonfinite evidence clear the verdicts.
CLI exit zero means a validated evolution experiment completed; the separate
impulse and reciprocal-exchange verdicts remain false. Failed execution
withdraws stale affirmative output.

The new tests cover full numerical replay, topology, nonconstant norm,
exact controls, corruption and failed-output paths. A positive control
confirms that the impulse classifier recognizes a finite, resolution-
independent jump; the actual evolved trajectories fail that criterion.
All 103 targeted tests pass locally: the combined 102-test run plus the
additional positive impulse-classifier control. CI includes all 103 on
Python 3.10 and 3.12 before the full suite.

### Publication recovery and lossless transport

The workspace rollback was repaired by reconstructing the recorded patches
and repeating the unchanged experiment. Both the original raw SHA256 above
and the complete original tree `a6fec369cbfde77fe73b3117cf5044544073457b`
were reproduced exactly. All 103 tests passed again in one run.

Large publication requests repeatedly failed while small GitHub writes
succeeded. The raw JSON is therefore transported as 16 small gzip/base64
parts plus a hash manifest. The restoration command above verifies every
part, the compressed stream and the registered raw hash before creating
`evolution.json`. CI and the evolution tests restore it automatically.
The physics modules, frozen specification, raw bytes and verdicts are
unchanged. The final publication tree differs only through this packaging,
its restoration tests and documentation; it is not the original tree.
