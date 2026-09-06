# Round 9 — can positive counting reach the quantum correlation?

**9/9 checks pass.**

## Verdict

* **Q1_positivity_bound** — `NO_BOUND_ALGEBRAIC_MAXIMUM_REACHED`
* **Q1_sup_CHSH** — `4.0`
* **Q2_quantum_law** — `QUANTUM_LAW_ATTAINABLE_BY_POSITIVE_COUNTING`
* **Q2_witness** — `Phi(D) = D^2 (1 - D/5)`
* **Q2_worst_error_vs_minus_cos** — `1.1102230246251565e-16`
* **no_go_requires_a_dual_certificate** — `True`
* **withdrawn** — `round 6's 'the distance to quantum mechanics is exactly the distance from |D| to D' — counting suffices`
* **remaining** — `nothing in the geometry selects Phi; and by sector_sum_identity, attaining the singlet law does not close the source-readout hazard`

## Checks

* PASS  O1  reduction D = t + sqrt(2t) cos psi holds as a law
* PASS  O5  marginals are 1/2 for every member of the class
* PASS  Q1  a nonnegative threshold weight reaches CHSH = 4
* PASS  Q2  an explicit nonnegative Phi reproduces -cos gamma
* PASS  Q2  and gives exactly Tsirelson
* PASS  S1  degree 3 is the minimal nonnegative solution
* PASS  S2  no polynomial solution is globally nonnegative
* PASS  S3  G_n has degree n-1 with leading coefficient 1
* PASS  R5  nothing in rounds 5-8 or #283 moved

## The class spans everything

| weight `Φ(D)` | nonnegative? | standard-angle CHSH |
|---|---|---:|
| `|D| (round 5)` | yes | 2.1422831632 |
| `D (round 6, signed)` | **no** (comparison) | 2.8284271247 |
| `D^2(1-D/5) (round 9)` | yes | 2.8284271247 |
| threshold bump at `1+√2` | yes | 4.0000000000 |

Every **nonnegative** row has marginals exactly `1/2`; the signed row is listed only to show the quantum value was previously reached only that way.

## The explicit quantum witness `Φ(D) = D²(1 − D/5)`

| `γ` | `E` | `−cos γ` | error |
|---|---:|---:|---:|
| 0.3 | -0.955336489126 | -0.955336489126 | 1.1e-16 |
| 0.5 | -0.877582561890 | -0.877582561890 | 1.1e-16 |
| 1.0 | -0.540302305868 | -0.540302305868 | 1.1e-16 |
| 1.5 | -0.070737201668 | -0.070737201668 | 1.1e-16 |
| 2.0 | +0.416146836547 | +0.416146836547 | 5.6e-17 |
| 2.5 | +0.801143615547 | +0.801143615547 | 1.1e-16 |
| 3.0 | +0.989992496600 | +0.989992496600 | 0.0e+00 |

* minimal degree is 3: `G_2` odd part `1.0`, `G_3` odd part `5.0`, condition `a_2 + 5 a_3 = 0`
* combined `G(u) = 1.2 + (-0.2) u²` — even about `t = 1`
* `min Φ` on `[-1/2, 4]` = `2.50e-11`
* top degree is always odd, so no globally nonnegative polynomial solution exists

## Dependency ledger

* `class restriction: W_s = int Phi(D_s) dsigma, Phi >= 0` — **chosen** (excludes (int D)^2, x-dependence beyond D, sector-dependent Phi)
* `geodesic itinerary and realignment` — **chosen** (inherited from rounds 5-6)
* `conditioning variable = phase` — **chosen** (inherited from round 8)
* `equal sector prior` — **chosen** (inherited; the class gives 1/2 marginals for any paired prior)
* `nonnegativity of Phi on [-1/2, 4]` — **derived** (the range the model evaluates; global nonnegativity is impossible for polynomials (no_global_polynomial))
* `which Phi the physics selects` — **open** (THE remaining question; #283's equilibrium selects |.|)
