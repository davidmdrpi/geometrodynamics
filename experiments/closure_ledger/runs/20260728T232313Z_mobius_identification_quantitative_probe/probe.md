# Is the Möbius identification quantitative? (PR #229)

_Run 2026-07-28T23:23:13.833436+00:00 · 6/6 PASS_

## The factorization

```
Delta M^2  =  (1/2)        x   (2 pi sigma)
              DERIVED           INHERITED
```

## The 1/2 is topological

| grid | circumference | antiperiodic ground n | dev from 1/2 |
|---:|---:|---:|---:|
| 120 | 0.7 | 0.499985721 | 1.4e-05 |
| 120 | 1.0 | 0.499985721 | 1.4e-05 |
| 120 | 3.3 | 0.499985721 | 1.4e-05 |
| 240 | 0.7 | 0.499996430 | 3.6e-06 |
| 240 | 1.0 | 0.499996430 | 3.6e-06 |
| 240 | 3.3 | 0.499996430 | 3.6e-06 |
| 480 | 0.7 | 0.499999108 | 8.9e-07 |
| 480 | 1.0 | 0.499999108 | 8.9e-07 |
| 480 | 3.3 | 0.499999108 | 8.9e-07 |
| 960 | 0.7 | 0.499999777 | 2.2e-07 |
| 960 | 1.0 | 0.499999777 | 2.2e-07 |
| 960 | 3.3 | 0.499999777 | 2.2e-07 |

## The intercept gap in #100

| C | `E0` periodic | `E0` antiperiodic | `ΔE0` | `π/(4C)` |
|---:|---:|---:|---:|---:|
| 1.0 | -0.523599 | +0.261799 | +0.785398 | +0.785398 |
| 2.0 | -0.261799 | +0.130900 | +0.392699 | +0.392699 |
| 3.5 | -0.149600 | +0.074800 | +0.224399 | +0.224399 |

## Scope

This touches the GLUEBALL tower of #100 only. It does NOT touch the heavy Moebius BARYON predictions of #103/#109/#114 (Lambda_c 3135, Lambda_b 6469, the 849 MeV dipion endpoint), which are built from the flux-tube quantum Delta = 2 sqrt(sigma) rather than the closed-loop intercept. The search table stands as published. Glueballs are unobserved, so the corrected item is the one with no experimental exposure.

## Verdict

**THE_MOEBIUS_IDENTIFICATION_IS_QUANTITATIVE_IN_ITS_TOPOLOGICAL_HALF_THE_SCALE_IS_INHERITED_AND_THE_MOEBIUS_INTERCEPT_IS_AN_OPEN_CORRECTION**
