# Scale-bridge regularization probe

**Run:** 2026-05-13T01:23:05+00:00

Asks the regularization-convergence question for the 1.054 factor (ω(l=1, n=0) at the closure-quantum R*). The 5D Tangherlini radial problem has a singular inner boundary at the throat r = r_s; the closure-ledger code regularizes by truncating the grid at r = r_s + ε with ε = 5×10⁻⁴. If the lowest eigenvalue depends on ε, then the 1.054 is not a converged structural constant of the Tangherlini geometry but a property of this specific regularization scheme.

## (1) Grid-resolution convergence (sanity check)

ω at the locked (R* = 1.262636, ε = 5×10⁻⁴) for increasing Chebyshev-grid resolution N. The closure-ledger probes use N = 80 throughout.

| N | ω(1, 0) |
|---:|---:|
| 60 | 1.0535269266 |
| 80 | 1.0535269266 |
| 120 | 1.0535269266 |
| 160 | 1.0535269266 |
| 200 | 1.0535269266 |
| 240 | 1.0535269266 |

## (2) Regularization-ε convergence (the test)

ω(1, 0) at fixed R = 1.262636 for a sweep of inner-boundary regularization ε. The locked baseline uses ε = 5×10⁻⁴.

| ε | ω(1, 0) | (ω − 1)·100 |
|---:|---:|---:|
| 1e-02 | 1.895296 | +89.5296 |
| 5e-03 | 1.597802 | +59.7802 |
| 2e-03 | 1.326942 | +32.6942 |
| 1e-03 | 1.175363 | +17.5363 |
| 5e-04 | 1.053527 | +5.3527 |
| 2e-04 | 0.924745 | -7.5255 |
| 1e-04 | 0.845361 | -15.4639 |
| 5e-05 | 0.777802 | -22.2198 |
| 1e-05 | 0.654427 | -34.5573 |

**ω spread across the ε sweep:** 1.2409. ω at the smallest ε (1e-05) is 0.6544.

**ω is NOT converged in ε at fixed R.** The lowest eigenvalue depends strongly on where the inner regularization is placed. The 1.054 value is the eigenvalue at ε = 5×10⁻⁴ specifically — not a structural constant of the Tangherlini boundary problem.

**Why this happens.** In tortoise coordinates `r* = r + (r_s/2)·ln|(r − r_s)/(r + r_s)|`, the throat r → r_s maps to r* → −∞. Putting a hard wall at r = r_s + ε corresponds to a hard wall at r* finite but log-divergent in ε: as ε → 0 the box becomes infinite in extent. The lowest eigenvalue tracks the box width and is not a converged Sturm-Liouville eigenvalue of the bare radial operator.

## (3) Self-consistent R*(ε) under the closure-quantum reading

For each ε, re-bisect R*_μ and R*_τ under the closure-quantum reading (transport = 8π, resistance = 7π/100). Report R*(ε), cross-species agreement, γ at R*, and ω at R*.

| ε | R*_μ | R*_τ | agreement | γ at R* | ω at R* |
|---:|---:|---:|---:|---:|---:|
| 5e-03 | 1.267136 | 1.267163 | 0.0021% | 22.5076 | 1.591606 |
| 2e-03 | 1.264136 | 1.264163 | 0.0021% | 22.5076 | 1.325684 |
| 1e-03 | 1.263137 | 1.263163 | 0.0021% | 22.5076 | 1.175060 |
| 5e-04 | 1.262636 | 1.262662 | 0.0021% | 22.5076 | 1.053527 |
| 2e-04 | 1.262338 | 1.262365 | 0.0021% | 22.5076 | 0.924840 |
| 1e-04 | 1.262236 | 1.262263 | 0.0021% | 22.5076 | 0.845462 |

**Self-consistent ω spread across the ε sweep:** 0.7461. R* spread: 0.3867%. γ spread: 0.0000%.

**The closure-quantum self-consistency loop does NOT pin** **ω at R* across ε.** Both R* and ω drift with the regularization. The 1.054 factor inherits its value from the chosen ε = 5×10⁻⁴; smaller ε produces a different ω at a different R*.

## (4) Clean-target ε candidates

ε values at which ω lands within 0.1 % of a simple closed-form target. A clean hit identifies a specific regularization at which the dimensional bridge would close.

### Fixed R = 1.262636

| ε | ω | target | target value | %Δ |
|---:|---:|---|---:|---:|
| 5e-04 | 1.053527 | `(10/9)^(1/2)` | 1.054093 | -0.0537% |

### Self-consistent R*(ε)

| ε | ω | target | target value | %Δ |
|---:|---:|---|---:|---:|
| 5e-04 | 1.053527 | `(10/9)^(1/2)` | 1.054093 | -0.0536% |

## (5) Compton-bridge ε (ω = 1 exactly)

ω at fixed R = R*_locked is decreasing in ε. The bisection finds the ε at which ω(1, 0; R*, ε) = 1 — the Compton-bridge condition `ℏ = m_e R_MID c` would close exactly at that regularization.

**ε at which ω = 1:** `3.5087e-04` (ω at this ε: 1.00000000).

Natural ε candidates derived from BAM ingredients, evaluated at R*_locked:

| candidate | ε | ω(1, 0) | (ω − 1)·100 |
|---|---:|---:|---:|
| `1 / (100·π)` | 3.1831e-03 | 1.451682 | +45.1682 |
| `1 / (1000·π)` | 3.1831e-04 | 0.986158 | -1.3842 |
| `1 / (8π·100)` | 3.9789e-04 | 1.018412 | +1.8412 |
| `(R* − 1) / 1000` | 2.6264e-04 | 0.959856 | -4.0144 |
| `(R* − 1) / 500` | 5.2527e-04 | 1.061404 | +6.1404 |
| `(R* − 1)² · 1e-2` | 6.8978e-04 | 1.106989 | +10.6989 |
| `π / 10⁴` | 3.1416e-04 | 0.984320 | -1.5680 |
| `1 / (2π)³` | 4.0314e-03 | 1.524542 | +52.4542 |

**Closest natural ε candidate to ω = 1:** `1 / (1000·π)` (ε = 3.1831e-04, ω = 0.9862, -1.384% from 1).

## Implications for the scale bridge

Three structural facts emerge from this probe:

1. **(R*, γ) are ε-invariant.** Across 6 ε values spanning 1e-04 to 5e-03, the closure-quantum cross-species fixed point R* shifts by 0.3867 %, and γ at R* by 0.0000 %. The mass ratios m_μ/m_e and m_τ/m_e are predicted by the closure-quantum machinery at the SAME precision regardless of the inner regularization.

2. **ω(1, 0) at R* is NOT ε-invariant.** ω drifts from 1.5916 to 0.8455 across the same ε sweep. The 1.054 reading at the closure-ledger default ε = 5×10⁻⁴ is therefore NOT a converged Sturm-Liouville eigenvalue of the Tangherlini geometry — it is the eigenvalue at that specific inner regularization, and the closure-quantum machinery does NOT pin it.

3. **The Compton bridge is restorable.** At ε* = 3.5087e-04, ω(1, 0; R*, ε*) = 1 exactly. The dimensional bridge `ℏ = m_e · R_MID · c` would close with NO 1.054 factor at this regularization. The closure-quantum machinery still predicts the mass ratios (self-consistency holds at ε = 2×10⁻⁴ and below). So the Compton-bridge geometry, vetoed by the lepton spectrum in probe 7 under the canonical ε, is **recovered** at the Compton-bridge ε.

**Reframing of the m_e / 1.054 scale problem.** The original factor-1054 search asked: 'is 1.054 expressible in closed form in (k_5, π, integers)?' That probe returned a clean negative result. This probe shows the question was the wrong one. The 1.054 factor is not a structural Tangherlini eigenvalue — it is the lowest eigenvalue at a specific (R, ε) point that the closure-ledger code chose by numerical convenience. The proper structural object is the regularization ε itself.

Two readings are now on the table:

- **Compton-bridge reading.** Take ε = ε* ≈ 3.509e-04 (the regularization at which ω = 1). The dimensional bridge is `ℏ = m_e R_MID c` with no remaining factor. BAM is dimensional-scale-incomplete with the m_e anchor as the unique external input. The 1.054 factor of the prior framing is absorbed into the regularization choice.
- **Default-regularization reading.** Take ε = 5×10⁻⁴ (closure-ledger default). The dimensional bridge has the form `ℏ · ω(1, 0) = 1.054 · m_e c²`. Both ε and the 1.054 are external choices, but they are not independent — fixing ε determines ω.

The nearest natural-BAM candidate ε is `1 / (1000·π)` = 3.183e-04, giving ω = 0.9862 (-1.384 % from 1). This is suggestive — `1/(1000·π)` involves only π and the τ-uplift quantum 100 scaled by a factor of 10 — but the 1.4 % gap from ω = 1 is large compared to the closure-quantum precisions established in PR #16 (transport 0.13 %, resistance 0.94 %, γ 0.03 %). The Compton-bridge ε does not have a closed form in BAM ingredients at this probe's precision.

**What this leaves open.** The dimensional bridge to ℏ requires choosing the inner-boundary regularization. Two follow-up directions:

1. **Derive ε structurally.** Identify a natural BAM ingredient (closure-quantum integers, throat geometry, Hopf invariants) that selects ε* without external input. If achievable, the Compton bridge closes and BAM is ratio-and-scale-complete modulo m_e.
2. **Replace hard-wall with quasi-regular boundary.** The tortoise-coordinate hard wall at finite ε is a numerical convenience. A boundary condition derived from the throat dynamics (THESIS.md 'Self-consistent throat radius') would remove the regularization dependence entirely.

Either direction makes the question SHARPER than the factor-1054 framing: the residual external input has been reduced from 'find a closed form for 1.054' to 'derive the inner regularization ε from BAM ingredients'.