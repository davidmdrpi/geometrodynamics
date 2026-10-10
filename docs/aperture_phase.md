# Controlled aperture phase family: measured result

The 57 trajectories were run once after publication of freeze
`70736ccaa7fb3693bbabb404ba87af09a5d5b3be` in [PR #324](https://github.com/davidmdrpi/geometrodynamics/pull/324).
The [preregistered](aperture_phase_prereg.md) numerical and phase-coverage gates
are **True** and **True**, respectively.

| Registered claim | Frozen result |
|---|---|
| Every designated high-phase case retains >=10% of round capture | **FAILED_IN_DECLARED_FAMILY** |
| Both specified fits obey the approximate inverse-phase law | **SUPPORTED_IN_DECLARED_FAMILY** |

9 of 10 designated high-phase cases fall below 10% retention.
The smallest is 0.02426. This rejects uniform 10%
retention in the declared family. It does not show zero transport, nor does it
exclude other port models or geometries. The inverse-phase result is a
three-carrier, finite-range test at b=.8, not an asymptotic scaling theorem.
Neither result establishes quantum mechanics or self-consistent feedback.

![Capture and retention](figures/aperture_phase.png)

## What was controlled

Primary a*w=4.8; the second coherent-port family has a*w=7.2. The carrier
values are 12, 24 and 48; gamma/w=2/3 and the compact source contains three
cycles at each carrier. This fixes coordinate aperture/wavelength, relative
pulse bandwidth and dimensionless damping strength. All captures use the
complete incident energy and the same stop time 1.75 in round-transit units.
The horizon is not enlarged after observing a poor capture. Changes in bulk
spectral density and the geometry remain part of the question, rather than
being normalized away. The numerical implementation removes only exact dark
and +/-m degeneracies of the #323 port representation.

These distributed L2-normalized coherent transducers remain assumptions.
They are not physical holes with a derived boundary-matching law. Gamma is
scaled as a controlled parameter, not determined by GR. The corrective-kick
budget is zero. Static field-plus-lead energy is accounted for throughout.
There is no metric pump, dynamical R3 background, moving mouth or reinjection.

## Primary family, a*w=4.8

Capture and remaining energy are percentages of complete incident energy.
Retention is relative to the matching round case. Exact phase spread uses the
registered round free-source modal weights, not measured receiver amplitudes.

| b | w | Proxy Phi | Exact weighted phase std | Capture (%) | Retention | Remaining bulk (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.8 | 12 | 10.603 | 2.700 | 4.0845 | 0.0964 | 59.967 |
| 0.8 | 24 | 21.206 | 5.464 | 2.1398 | 0.0506 | 62.409 |
| 0.8 | 48 | 42.412 | 10.961 | 1.0266 | 0.0243 | 63.681 |
| 0.9 | 12 | 4.422 | 1.189 | 14.6524 | 0.3459 | 48.514 |
| 0.9 | 24 | 8.843 | 2.408 | 4.4245 | 0.1045 | 59.135 |
| 0.9 | 48 | 17.686 | 4.832 | 2.3762 | 0.0562 | 61.187 |
| 0.98 | 12 | 0.777 | 0.217 | 40.9074 | 0.9657 | 21.466 |
| 0.98 | 24 | 1.554 | 0.439 | 35.9546 | 0.8495 | 26.394 |
| 0.98 | 48 | 3.109 | 0.882 | 21.6863 | 0.5125 | 40.650 |
| 1 | 12 | 0.000 | 0.000 | 42.3607 | 1.0000 | 19.695 |
| 1 | 24 | 0.000 | 0.000 | 42.3241 | 1.0000 | 19.700 |
| 1 | 48 | 0.000 | 0.000 | 42.3149 | 1.0000 | 19.697 |
| 1.02 | 12 | 0.732 | 0.207 | 40.4585 | 0.9551 | 21.276 |
| 1.02 | 24 | 1.464 | 0.421 | 35.7335 | 0.8443 | 25.962 |
| 1.02 | 48 | 2.928 | 0.844 | 22.0190 | 0.5204 | 39.664 |
| 1.1 | 12 | 3.271 | 0.955 | 17.0178 | 0.4017 | 43.408 |
| 1.1 | 24 | 6.543 | 1.938 | 5.5647 | 0.1315 | 54.808 |
| 1.1 | 48 | 13.086 | 3.891 | 3.0513 | 0.0721 | 57.307 |
| 1.2 | 12 | 5.760 | 1.736 | 6.2604 | 0.1478 | 52.571 |
| 1.2 | 24 | 11.519 | 3.526 | 3.4774 | 0.0822 | 55.300 |
| 1.2 | 48 | 23.038 | 7.079 | 1.7206 | 0.0407 | 57.043 |

## Second aperture family, a*w=7.2

| b | w | Proxy Phi | Exact weighted phase std | Capture (%) | Retention | Remaining bulk (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.8 | 12 | 10.603 | 2.128 | 1.9118 | 0.1062 | 34.264 |
| 0.8 | 24 | 21.206 | 4.330 | 1.2034 | 0.0674 | 35.696 |
| 0.8 | 48 | 42.412 | 8.697 | 0.5605 | 0.0315 | 36.477 |
| 1 | 12 | 0.000 | 0.000 | 18.0078 | 1.0000 | 11.552 |
| 1 | 24 | 0.000 | 0.000 | 17.8441 | 1.0000 | 11.594 |
| 1 | 48 | 0.000 | 0.000 | 17.8034 | 1.0000 | 11.570 |

## Preregistered inverse-phase fits

Fit log(retention)=intercept-beta*log(Phi), only b=.8, at both footprints.
The proxy spans 10.60–42.41 rad. The registered acceptance region is
beta in [.5,1.5] and maximum absolute log residual <=.25 for **both** fits.
No alternative b slice, peak-flux exponent or time window replaces this test.

| a*w | Fine beta | Coarse beta | Fine max log residual | Within declared law bounds |
|---:|---:|---:|---:|---|
| 4.8 | 0.995334 | 0.995412 | 0.029531 | True |
| 7.2 | 0.876802 | 0.876967 | 0.102671 | True |

These deterministic fits have only three points each. A passing law is only
supported in this finite interval. The exact source-weighted phase standard
deviation passes the separately stated coverage gate, but a broad asymptotic
regime and extrapolation to throat-scale resolution remain unestablished.
The near-round deformation cases report behavior outside the fit range;
they are not added as convenient points to improve the fitted law.

## Verification and energy accounting

| Diagnostic | Largest observed | Frozen limit |
|---|---:|---:|
| energy_error | 4.56498e-13 | 1e-09 |
| final_energy_error | 0 | 1e-09 |
| port_error | 4.85723e-15 | 1e-08 |
| reconstructed_energy_error | 3.45121e-15 | 1e-08 |
| state_error | 3.7115e-15 | 1e-08 |
| early_flux | 4.5192e-19 | 1e-05 |
| cap_tail | 3.79184e-06 | 0.0001 |

Largest paired/refinement comparison is 0.00786581
in source-normalized two-port output L2, against .03. Separate time, mode and
extension comparisons are 7.29039e-05,
1.59296e-09 and 0.
The extension adds 4.20197e-11 of incident energy to B capture after 1.75 in
the one designated case. Its later flux is reported, not folded into the
primary measurement or either verdict. Reflected A energy plus captured B
energy plus remaining bulk energy closes the ledger in every trajectory.

Independent reconstruction uses archived forces incoming-outgoing, rebuilding
all field states and energies without re-solving the coupled port evolution.
Source formulas allow absolute roundoff <=1e-14 while grid, compact support
and inactive port remain exact. Stored incoming samples are used in physics
ledgers after validation. The frozen thresholds and source files are intact.

All 58 focused tests pass. A complete authenticated replay with AVX2, FMA3
and AVX512F NumPy features disabled reproduces both labels, both exponents
and every archived capture measurement. All 57 NPZ archives were inspected
with pickle disabled: 342 finite float64 arrays, with the declared configuration
metadata, and every manifest hash verified. CI runs the focused tests and full
archive reconstruction before the repository-wide suite.

## Evidence and chronology

- Run start: `2026-10-10T03:49:20.810600+00:00`.
- Run finish: `2026-10-10T03:50:44.676926+00:00`.
- Runtime: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1.
- [All 57 raw trajectories, result, manifest and timestamped provenance](../experiments/closure_ledger/runs/20261010_aperture_phase/).
- Manifest SHA-256: `02b8bd392282292cbf0f82a42dfdf17ba2263a84b709bebf38217fa394037a6d`.

Per-case start and finish UTC are recorded directly; no NPZ ZIP timestamp is
used as a run clock. The manifest authenticates the original arrays before
reconstruction. Replay checks all source hashes, measurements and both
verdicts. It never launches new production trajectories.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.aperture_phase_probe --manifest-sha 02b8bd392282292cbf0f82a42dfdf17ba2263a84b709bebf38217fa394037a6d
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_aperture_phase.py tests/test_aperture_transfer.py tests/test_mty_packet.py
python -m experiments.closure_ledger.aperture_phase_report
```

## Consequence for the next physical test

This family measures how much the static result depends on phase accumulation
under declared source and port scaling. It supplies no reason to extrapolate
static suppression unchanged to an evolving triaxial R3 geometry: anisotropy
can change during a transit, and mode coupling and metric work must then be
included. A future R3 study needs a new freeze, independently sampled geometry
phases, an explicit metric/port-work ledger and a pre-stated capture failure
criterion. The present study does not claim to have carried out that test.
