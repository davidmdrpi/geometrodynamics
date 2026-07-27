# The Z2 link twist as a gauge field over the ER network (PR #227)

_Run 2026-07-27T02:51:58.824787+00:00 · 8/8 PASS_

## The headline: a twisted link freezes mouth exchange

| r_s | Dw periodic | Dw twisted | ratio | T_exchange periodic | twisted |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 9.464e-04 | 1.9e-13 | 2.1e-10 | 3319 | 1.62e+13 |
| 0.3 | 1.714e-03 | 8.3e-14 | 4.9e-11 | 1833 | 3.76e+13 |
| 0.4 | 4.747e-03 | 2.4e-13 | 5.2e-11 | 662 | 1.28e+13 |
| 0.5 | 1.079e-02 | 2.4e-13 | 2.2e-11 | 291 | 1.31e+13 |

The S^2 = W degeneracy theorem (max intra-pair gap, low spectrum):

  - periodic (W=+1): 2.08e-03
  - twisted  (W=-1): 7.10e-12

## The b_1 audit (the plan's primary falsifier)

| topology | V | E | b_1 | reading |
|---|---:|---:|---:|---|
| #224 two-throat ring (A,B + 2 exterior arcs) | 2 | 2 | 1 | as built |
| #223 single bridge + exterior arc | 2 | 2 | 1 | as built |
| #207 perfect matching, 6 mouths (bridges only) | 6 | 3 | 0 | bridges only -- the naive falsifier |
| #207 perfect matching + shared S^3 exterior | 6 | 9 | 4 | bridges + exterior -- the physical graph |
| #208 Y-junction (3 mouths + core), bare | 4 | 3 | 0 | bare junction -- a tree |
| #208 Y-junction + shared S^3 exterior | 4 | 6 | 3 | junction + exterior -- the physical graph |

## The quotient revival

| sector | t | \|f-f0\| | \|f+f0\| | \|f+Pf0\| |
|---|---|---:|---:|---:|
| full | pi R | 1.41e+00 | 1.41e+00 | 1.68e-15 |
| full | 2 pi R | 0.00e+00 | 2.00e+00 | 1.41e+00 |
| untwisted_even_l | pi R | 2.00e+00 | 0.00e+00 | 1.69e-15 |
| untwisted_even_l | 2 pi R | 0.00e+00 | 2.00e+00 | 2.00e+00 |
| twisted_odd_l | pi R | 0.00e+00 | 2.00e+00 | 1.65e-15 |
| twisted_odd_l | 2 pi R | 0.00e+00 | 2.00e+00 | 1.65e-15 |

## Corrections to the research plan

  - NO factor-of-4 energy advantage: energy conserved exactly, per-sector caustic gains agree to <1%, the 2 m_e c^2 threshold of #58/#166 is untouched
  - the vertex rule Sum l == n_twist (mod 2) is a TAUTOLOGY (eps = (-1)^l is forced); the real content is the reinterpretation of #137 plus the network junction rule

## Verdict

**ER_LINK_TWIST_IS_A_Z2_GAUGE_FIELD_THAT_FREEZES_MOUTH_EXCHANGE_AND_HALVES_THE_QUOTIENT_REVIVAL_BUT_LEAVES_THE_NUCLEATION_THRESHOLD_UNTOUCHED**
