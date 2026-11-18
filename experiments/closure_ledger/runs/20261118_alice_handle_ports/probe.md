# Alice handle ports: probe against freeze `36ba9c2`

A real scalar on the unit sphere with one handle joining antipodal mouths on the z-axis, glued torus-type (twist c) or Klein-type. Fixed geometry, separable finite-element solver per azimuthal number, evolution by mode sum; #338's mesh as a cross-check. Classical; not GR, not S3.

**Verdict: COUNTEREXAMPLE** (P4 split)

* also missed: P5 identity (not shrinking)

| statement | N = 2000 | N = 4000 |
|---|---|---|
| P1 | pass (0) | pass (0) |
| P2 degeneracy | pass (1.31e-14) | pass (1.92e-14) |
| P2 conservation | pass (1.67e-10) | pass (2.07e-10) |
| P3 | pass (1.19e-14) | pass (2.99e-14) |
| P4 split | FAIL (1.01e-08) | FAIL (1.01e-08) |
| P4 timing | pass (0.318) | pass (0.318) |
| P5 identity | FAIL (7.03e-07) | FAIL (8.59e-06) |
| P5 non-vacuous | pass (0.989) | pass (0.989) |
| P6 sign | pass (7.82e-08) | pass (1.7e-08) |
| P6 tolerance | pass (1.56e-07) | pass (3.4e-08) |

P7 geometry: True; spectrum: True
* torus: largest m among the lowest 20 = 3; max relative difference 6000/48: 0.00333, 12000/64: 0.0018
* klein: largest m among the lowest 20 = 3; max relative difference 6000/48: 0.00333, 12000/64: 0.0018

Lowest Klein splittings Δ/ω̄ (first resolution):

* a = 0.3, m = 1: 0.00316, 0.00536, 0.00771, 0.0103, 0.0132, 0.0164  (expectation 0.023)
* a = 0.3, m = 2: 5.51e-06, 1.96e-05, 4.61e-05, 8.98e-05, 0.000156, 0.000252  (expectation 0.00052)
* a = 0.3, m = 3: 1.01e-08, 5.52e-08, 1.83e-07, 4.71e-07, 1.04e-06, 2.08e-06  (expectation 1.2e-05)
* a = 0.5, m = 1: 0.0262, 0.0462, 0.0694, 0.0934, 0.113, 0.122  (expectation 0.065)
* a = 0.5, m = 2: 0.000379, 0.0014, 0.00351, 0.00743, 0.0142, 0.025  (expectation 0.0043)
* a = 0.5, m = 3: 5.72e-06, 3.24e-05, 0.000113, 0.000311, 0.000751, 0.00166  (expectation 0.00028)

Assembled in 0 s from stages cached by earlier runs; their run times are not recorded here.
