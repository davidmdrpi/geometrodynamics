# Constrained R3 return map: leading frequency-shift coefficient

Date: 2026-09-29.
- Specification and code: [`cff0ec0`](r3_return_map_prereg.md), pushed 02:04:23 UTC, before any cubic-order ESU computation.
- Archive: `experiments/closure_ledger/runs/20260929_r3_return_map/return_map.json`. It binds the SHA-256 of all five sources and stores the Method 1 jet coefficients, both normal forms, and every circle (K, omega, I, residuals).
- Replay: `tests/test_r3_return_map.py::test_return_map_archive_rescores_and_binds_sources`.

## 1. Registered label: UNRESOLVED

Check **C6 (nonresonance) failed**. Every other check passed. By the frozen
rule the label is UNRESOLVED, and it is not relabelled.

| id | result | value |
|---|---|---|
| C1 fixed point | pass | Method 1: 0; Method 2: 1.0e-12 |
| C2 linear | pass | theta_M1 = theta0 = 3.045248856 (#310 trace); the finite-difference eigen-angle agrees |
| C3 symplectic | pass | jet defect 2.4e-8 on coefficients up to 5.0e5; finite-difference defect 3.6e-6 |
| C4 dissipative | pass | Re(g/mu) = 2.1e-14, against Im(g/mu) = -3.7066 |
| C5 constraint | pass | jet 2.5e-8 (scale 5e5); section q 2.7e-15; circles <= 9.7e-13 |
| **C6 nonresonance** | **fail** | three quadratic divisors below .1 (section 3) |
| C7 Method 1 convergence | pass | 1.1e-13 |
| D1 invariance | pass | <= 1.8e-12 |
| D2 Fourier tail | pass | <= 1.1e-14 |
| D3 Radau | pass | 1.2e-12 at a = .008 |
| A1 agreement | pass | see section 2 |
| A2 remainder | pass | slopes 2.0004 and 2.0011 |

## 2. Measurements

| Quantity | Value |
|---|---|
| nu_M1 (jets + normal form, Richardson 2048/4096) | **-0.950118097** turns per unit action |
| nu_M1, Richardson 1024/2048 | same to 1.1e-13 |
| nu_M2 (invariant circles, I -> 0) | **-0.950118084**, uncertainty 4.5e-9 |
| \|nu_M1 - nu_M2\| | 1.3e-8 (1.4e-8 relative) |
| Ablation: fixed time pi, no return-time correction | -0.8086 (not a section map) |

The invariant circles, with rho0 frac = theta0/2π = .48466641:

| a | I | omega | (omega - theta0)/(2π I) | remainder after nu_M1 I |
|---|---|---|---|---|
| .004 | 2.2985e-5 | 3.0451116429 | -0.9501230 | -7.1e-10 |
| .008 | 9.1937e-5 | 3.0447000038 | -0.9501379 | -1.14e-8 |
| .016 | 3.6772e-4 | 3.0430534630 | -0.9501972 | -1.83e-7 |
| .032 | 1.4705e-3 | 3.0364675477 | -0.9504350 | -2.93e-6 |

Once nu_M1 I is subtracted, the remainder is quadratic in I to three digits.
Its coefficient is about -1.35 in radians. The return-time jet has
quadratic coefficients -0.450 in x^2 and -0.0549 in p_x^2, and none in A.
The clock frequency Omega^2 depends only on the shape. The return-time
correction changes nu by 15%, which the ablation row shows.

## 3. Why C6 failed: a specification error

This section is a post-hoc diagnostic, computed from the archived jet.

C6 requires every quadratic divisor to be at least .1 in absolute value.
The three that fail all belong to the contracting hyperbolic component s,
whose multiplier is lambda_s = .01176:

| component | monomial | \|divisor\| | relative to the larger multiplier |
|---|---|---|---|
| s | s^2 | .0116 | .988 |
| s | s w | .0235 | 1.998 |
| s | s wbar | .0235 | 1.998 |

These divisors are small only because lambda_s is small. Measured relative
to the multipliers involved, they are order one, far from any resonance.
The specification's own §7 describes C6 as the nonresonance condition for
the elliptic pair, with smallest divisor |mu^2 - 1| ≈ .19. The divisors
that actually govern that condition, those of the w and wbar components,
are all at least .988, and |mu^k - 1| >= .192 for k = 1..4. An absolute
threshold applied to every component was the wrong way to write that
condition. The error is mine, and it was avoidable: lambda_s was known
before the freeze.

## 4. What the numbers indicate, without a verdict

- **The two methods agree.** They are independent: separate equations,
  integrators and algorithms, and a jet expansion against invariant circles
  on the full system. Both give nu = -0.950118 to 1.4e-8.
- **The expansion is validated.** The remainder scales as I^2 exactly as
  predicted, and every symplectic, constraint and convergence check passes.
- **The resolution criterion is met.** |nu| exceeds every error estimate
  by more than 1e7.
- **If C6 had been written as §7 describes it, the frozen rules would give
  SHIFT_AWAY_FROM_TARGET.** Near the breathing ESU, the LRS n=2 tensor's
  rotation moves away from the refocusing resonance, from frac .4847
  toward smaller values.

**This is not adopted.** Any change to C6 now would be made after the
result was seen. It could be recorded only as a retrospective correction,
not as a prospective confirmation. Whether to record one is the user's
decision.

Even under that reading, the conclusion stays narrow. It says nothing about
large amplitudes beyond this ladder, where the largest circle has
I = 1.5e-3 and frac = .4833. It says nothing about other polarisations, the
n >= 3 sectors, or other closure conditions. It also does not bear on the
separate localised-receiver route.

The two archived numbers are not directly comparable. The kicked run's
descriptive coefficient was taken per signed displacement squared, and this
one is per action. At the smallest circle, I/a^2 = 1.437.

## 5. Disclosures

- **Early look at intermediate values.** The run took about 45 minutes,
  longer than estimated, because each circle's Newton solve stalls at
  about 1.6e-12. That is the integration floor, just above the 1e-12
  stopping tolerance, so every circle used all 15 iterations. While it
  ran, I inspected the live process with `py-spy dump --locals`. That
  exposed omega for the .016 circle and an intermediate omega for the
  .032 circle before the run finished. No code, parameter or rule changed
  afterwards.
- **Pre-freeze work** is listed in §8 of the specification and in
  `experiments/closure_ledger/r3_return_map_prefreeze.py`.

## 6. Reproduction

    python -m pytest -q tests/test_r3_return_map.py      # toy exactness, equations, linear map, archive re-score
    python -m experiments.closure_ledger.r3_return_map_probe --output /tmp/rerun.json   # about 45 min
