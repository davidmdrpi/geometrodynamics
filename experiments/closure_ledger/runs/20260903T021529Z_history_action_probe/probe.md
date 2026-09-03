# Round 7 — does a classical BAM action select the oriented branch?

**21/21 checks pass.**

## Five independent verdicts

* **A_action** — `HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION`
* **B_sectors** — `LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE`
* **C_readout** — `NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_READOUT`
* **D_compatibility** — `HISTORY_ACTION_INDEPENDENTLY_POSTULATED`
* **E_causality** — `SETTING_INFORMATION_IS_PRESENT_AT_SOURCE_READOUT_DYNAMICS_OPEN`
* **headline** — `NOT MET — no field of the pre-registered conjunction holds`

## Checks

* PASS  O1-O4  the frozen oracle reproduces
* PASS  O4     saddle magnitude IS the positive coarea density
* PASS  A-ctl  every F(cos theta) shares the critical set
* PASS  A1     theta additive, S_H not
* PASS  A2     the additive functional is nowhere stationary
* PASS  A3     Crit(S_H) = Gamma: no critical points OFF closure either
* PASS  F1     component masses are unequal; only the PHASE is undetermined
* PASS  F1b    the masses reproduce BOTH candidate aggregations exactly
* PASS  F2     phase factor real iff 4kappa/pi odd
* PASS  B1     two sector orbits at every angle but pi/2
* PASS  B2     r not forced at any CHSH angle
* PASS  B3     fibre symmetries cannot change sector weights
* PASS  C1     every existing BAM observable is degree 2
* PASS  C2     but two ordinary quadratics disagree: no readout is named
* PASS  C3     <D_s^2> = (1+c)(2+c) closed form
* PASS  D1     the radial action is a genuine one-form integral
* PASS  E1     source density is exactly antipodally even
* PASS  E2     odd observables blind; SOME even ones separate
* PASS  E2b    but NOT every even observable separates (constants, x.x)
* PASS  E4     no operational readout is claimed
* PASS  E3     non-coplanar supports are mutually singular

## The candidate, and why it is not an action

* `-1/2 Tr G = -cos(Omega/2) = -D/sqrt(D^2+N^2)` to `3.3e-16`
* closure set is the critical manifold: `|grad S_H|` <= `1.2e-06`
* transverse curvature `sgn(D)|uxv|^2/D^2`, rel `2.6e-06`; index is `sgn D`: True
* **saddle magnitude = positive coarea density**, rel `1.3e-06`

* `theta` additive to `3.3e-16`; `S_H` additivity defect `1.805`
* `min |grad theta|` on closure = `0.1620` > 0


* component masses `M_0 = 11.670806`, `M_pi = 0.169512`, ratio `0.014524` — **not** 1
* `(M_0 - M_pi)|uxv| = int D` (residual `1.8e-15`) and `(M_0 + M_pi)|uxv| = int |D|` (residual `0.0e+00`): stationary phase supplies both candidate magnitudes; only the relative phase is open
* no critical points off closure: `min |grad theta| = 0.2850`, and the only candidate needs `x_p = -sec(gamma/2)`, outside the sphere

| kappa | 4k/pi | phase real? | selects |
|---|---|---|---|
| pi/4 | 1.0000 | True | positive count  M_0 + M_pi |
| 3pi/4 | 3.0000 | True | oriented sum  M_0 - M_pi |
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
* quadratic, square-of-integral: `S_max = 3.771236` = `8 sqrt2/3`
* quadratic, integral-of-square: `S_max = 3.394113`
* both keep marginals at exactly `1/2`; **they disagree**, so degree-2 homogeneity does not name a readout

## E — causality gate

* conditioned density antipodally even to `3.4e-20`
* odd observables blind (spread `2.8e-17`)
* *some* even functions separate (spread `0.0103`); constants and `x.x` do not (spread `0.0e+00`)
* no map from `x` to field configurations, and no two-boundary-compatible source readout, is constructed — a hazard, not a channel
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
* `current-to-frequency readout` — **open** (question C: every existing observable is degree 2, but two ordinary quadratics disagree and none is a derived detector coupling)
* `a map from the source variable x to field configurations` — **open** (gate E: needed before degree in phi says anything about parity in x)
* `a source-local readout that respects the two-boundary problem` — **open** (gate E: not constructed; measurement is itself a BC)
* `antipodal scalar BC, eta, quotient-vs-cover` — **not used** (C1: not in this strand)
