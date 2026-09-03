# Round 7 — does a classical BAM action select the oriented branch?

**16/16 checks pass.**

## Five independent verdicts

* **A_action** — `HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION`
* **B_sectors** — `LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE`
* **C_readout** — `CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY`
* **D_compatibility** — `HISTORY_ACTION_INDEPENDENTLY_POSTULATED`
* **E_causality** — `SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS`
* **headline** — `NOT MET — no field of the pre-registered conjunction holds`

## Checks

* PASS  O1-O4  the frozen oracle reproduces
* PASS  O4     saddle magnitude IS the positive coarea density
* PASS  A-ctl  every F(cos theta) shares the critical set
* PASS  A1     theta additive, S_H not
* PASS  A2     the additive functional is nowhere stationary
* PASS  F1     saddle magnitudes never separate the branches
* PASS  F2     ratio real iff 4kappa/pi odd
* PASS  B1     two sector orbits at every angle but pi/2
* PASS  B2     r not forced at any CHSH angle
* PASS  B3     fibre symmetries cannot change sector weights
* PASS  C1     every BAM coupling is degree 2
* PASS  C2     a quadratic readout is superquantum
* PASS  D1     the radial action is a genuine one-form integral
* PASS  E1     source density is exactly antipodally even
* PASS  E2     odd observables blind, even observables signal
* PASS  E3     non-coplanar supports are mutually singular

## The candidate, and why it is not an action

* `-1/2 Tr G = -cos(Omega/2) = -D/sqrt(D^2+N^2)` to `3.3e-16`
* closure set is the critical manifold: `|grad S_H|` <= `1.2e-06`
* transverse curvature `sgn(D)|uxv|^2/D^2`, rel `2.6e-06`; index is `sgn D`: True
* **saddle magnitude = positive coarea density**, rel `1.3e-06`

* `theta` additive to `3.3e-16`; `S_H` additivity defect `1.805`
* `min |grad theta|` on closure = `0.1620` > 0

| kappa | 4k/pi | real? | selects |
|---|---|---|---|
| pi/4 | 1.0000 | True | positive count (+1) |
| 3pi/4 | 3.0000 | True | oriented (-1) |
| 1 (hbar=1) | 1.2732 | False | neither: complex relative weight |
| pi/2 | 2.0000 | False | neither: complex relative weight |

## B — sector orbits

| gamma | group order | orbits | r forced |
|---|---|---|---|
| 0.3000 | 4 | 2 | False |
| 0.7854 | 4 | 2 | False |
| 1.0000 | 4 | 2 | False |
| 1.5708 | 8 | 1 | True |
| 2.3562 | 4 | 2 | False |

## C — measured homogeneity of every BAM coupling

| coupling | where | degree |
|---|---|---|
| `null_stress_minimal_scalar` | `source_audit.py:115` | 2.000000 |
| `null_stress_complex_order_field` | `source_audit.py:126` | 2.000000 |
| `null_stress_maxwell` | `source_audit.py:137` | 2.000000 |
| `null_stress_perfect_fluid` | `source_audit.py:150` | 2.000000 |
| `null_stress_nonminimal_at_node` | `source_audit.py:157` | 2.000000 |
| `mouth_flux Im(q* A q)` | `throat_operator.py:676` | 2.000000 |

* linear readout: `S_max = 2.828427` (Tsirelson)
* quadratic readout: `S_max = 3.771236` = `8 sqrt2/3` at `gamma = 0.7854` — **superquantum**
* quadratic marginals stay `1/2` to `5.6e-17`

## E — causality gate

* conditioned density antipodally even to `3.4e-20`
* odd observables blind (spread `2.8e-17`)
* even observables signal (spread `0.0103`)
* non-coplanar total variation `1.0` — mutually singular, a one-shot signal

## Dependency ledger

* `Hopf bundle = P_Spin(S^2)` — **derived** (round 3 (PR #279))
* `reduction of the chosen itinerary to x -> u -> -v -> x` — **derived** (round 6 (PR #280))
* `geodesic-realignment ansatz at the detectors` — **chosen** (round 5)
* `the source -> A -> J -> B -> source itinerary` — **chosen** (round 5)
* `physicality of the pair direction x` — **chosen** (round 5; collides with gate E)
* `equal sector prior r = 1` — **open** (question B: free at every angle but pi/2)
* `S_H = -1/2 Tr G as the functional` — **chosen** (question A: stationary but not additive)
* `kappa, the normalisation in e^{i kappa S_H}` — **open** (question A/F3: no repository source)
* `orientation convention in the Maslov factor` — **chosen** (shifts kappa by pi/2; F2 is invariant)
* `linear current-to-frequency readout` — **open** (question C: every BAM coupling is degree 2)
* `antipodal scalar BC, eta, quotient-vs-cover` — **not used** (C1: not in this strand)
