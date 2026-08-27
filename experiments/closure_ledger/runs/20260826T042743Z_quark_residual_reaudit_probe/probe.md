# The quark residual sector, re-derived from the corrected operator

**8/8 checks pass.**

PR #271 corrected the 5D radial scalar operator and reopened the **lepton** sector's geometric closure. The frozen v3 **quark** lock derives three residual knobs from the same eigensolver, so the expectation was that they would break the same way. They do not.

## Q1 — all three residuals move toward their locked values

| residual | locked | legacy | corrected | legacy off | corrected off | drift |
|--|--|--|--|--|--|--|
| `pinhole` | `22.25` | `22.008240` | `22.331195` | `-1.09%` | **`+0.36%`** | `+1.467%` |
| `transport` | `0.54` | `0.544746` | `0.543800` | `+0.88%` | **`+0.70%`** | `-0.174%` |
| `resistance` | `0.14` | `0.140683` | `0.139965` | `+0.49%` | **`-0.02%`** | `-0.510%` |

*The readme's resistance residual (-0.43%) used the locked transport 0.54, not the derived one; both conventions are carried in `derived[*]['resistance*']`.*

## Q7 — why the two sectors part company

One barrier feeds both, and the correction moves it the same way for both. What differs is how hard each Hamiltonian leans on it.

| sector | scalar | corrected residual | elasticity | linearised | **measured** |
|--|--|--|--|--|--|
| lepton (`m_mu`) | gamma = Sum V_max[1..5] vs the locked 22.5 | `-0.75%` | `-17.471` | `+14.01%` | **`+15.16%`** |
| quark (`m_s`) | pinhole = Sum V_max[1..5] vs the fitted 22.25 | `+0.36%` | `+4.798` | `+1.77%` | **`+1.79%`** |

The elasticity ratio is `3.64`. The same barrier feeds both sectors and the correction moves it the same way; the lepton hamiltonian amplifies its residual by 17.5 and the quark ladder by 4.8, so a comparable geometric agreement is a 15% miss in one sector and a 2% miss in the other.

> **A percentage agreement between a geometric quantity and a fitted knob means nothing until multiplied by the elasticity of what it feeds.**

**A number this round corrects.** The README row written in #271 quoted `9%` for the muon error at the corrected `γ`. The elasticity and the residual in that row are both right; the `9%` is not — it came from a generic illustration in the #271 probe (*a half-percent error in γ is a nine-percent error in the muon*) applied to a residual that is not a half percent. Linearising gives `+14.01%`; the locked block returns `+15.16%`, and #271's own A/B/C table carried that value the whole time.

## Q2 — the composite is not the sum of its parts

Every substitution below goes into the **locked** quark Hamiltonian with **nothing retuned**, `d`-anchored, max relative error over `{s, c, b, t}`.

| substitution | legacy-derived | corrected-derived |
|--|--|--|
| `pinhole` alone | `3.61%` | `3.40%` |
| `transport` alone | `2.10%` | `2.00%` |
| `resistance` alone | `1.70%` | `1.61%` |
| **all three** | `3.44%` | **`3.78%`** |
| *the fitted v3 lock* | `1.61%` | `1.61%` |

The ordering **reverses**: each corrected residual is individually closer to its lock, yet the corrected triple scores worse. The legacy triple enjoys a partial cancellation its pinhole error alone would not earn; the corrected triple has none.

> **"Each knob is derived to within 1 %" and "the derived knobs reproduce the ladder" are different claims, and only the first was ever established** — under either operator. The correction did not create that gap, it exposed it.

## Q3 — the elasticity, and the pinhole the ladder actually wants

| species | local `d ln m/d ln pinhole` at the lock | secant across the two derived pinholes |
|--|--|--|
| `s` | `+4.7976` | `+4.8188` |
| `c` | `-1.5064` | `-1.6322` |
| `b` | `-1.8610` | `-1.9865` |
| `t` | `-2.0130` | `-2.1380` |

Scanning the locked ladder over the pinhole alone, it is minimised at `22.228022` (`1.348%`), not at the fitted `22.25` (`1.610%`). Measured against the ladder's **own** optimum rather than a fitted knob, the correction roughly halves the miss:

| pinhole | distance from what the ladder wants |
|--|--|
| v3 lock `22.25` | `+0.099%` |
| legacy `22.008` | `-0.989%` |
| corrected `22.331` | **`+0.464%`** |

## Q4 — the one thing that genuinely improves: cross-sector `R_OUTER`

Each sector, read as a demand on `R_OUTER` through its own locked scalar:

| operator | lepton (`γ=22.5`) | quark (`pinhole=22.25`) | split | brackets `1.26`? |
|--|--|--|--|--|
| legacy | `1.28737` (`+2.17%`) | `1.27229` (`+0.98%`) | `1.179%` | no |
| corrected | `1.26788` (`+0.63%`) | `1.25645` (`-0.28%`) | `0.906%` | **yes** |

Under the legacy operator both sectors demand **more** than `1.26`, putting the canonical value outside the bracket they define. Under the corrected operator they **straddle** it — the quark sector wants less, the lepton sector more, and `1.26` sits at `0.31` of the way across.

This is evidence *for* `R_OUTER = 1.26`, and it is **not** the evidence #271 reopened: that was a single-sector fixed point, this is a two-sector bracket. It is also weak — a 0.9% bracket admits any value inside it; this supports R_OUTER = 1.26 but does not derive it. Recorded as suggestive, not as a restored derivation.

## Q5 — the last lepton claim that reads a radial eigenvalue

README: *"Resistance = 7π/100 — selected over `4·(ω−1)` by `R_OUTER` bisection."*

| candidate | value | off the locked `0.217869` |
|--|--|--|
| `7π/100` (operator-independent) | `0.219911` | `+0.937%` |
| `4(ω−1)`, legacy | `0.218908` | `+0.477%` |
| `4(ω−1)`, corrected | `0.223306` | `+2.495%` |

Under the legacy operator the **rejected** candidate fitted *better* than the selected one, so the choice rested entirely on the `R_OUTER` bisection — which #271 reopened. Under the correction the competitor degrades to `+2.50%` and `7π/100` wins on proximity too.

> **Conclusion survives, stated reason does not: the r_outer-bisection selector is reopened by pr #271, and under the legacy operator it had selected the worse-fitting of the two candidates.**

## Q6 — `N` still drifts

| residuals pinned to | `N` range | width |
|--|--|--|
| baseline | `[430, 494]` | `64` |
| legacy | `[434, 480]` | `46` |
| corrected | `[438, 496]` | `58` |

Unchanged by the correction: beta remains the model's last fit knob. Pinning the residuals to *corrected* geometry does not stabilise `N` any more than pinning them to legacy geometry did. The v3 lock's own conclusion stands unchanged; only the digits move, `N` at PDG masses included.

## Q8 — the ledger

| claim | verdict | evidence |
|--|--|--|
| quark shell-index axioms (eps, eta, chi, phase) in k_5 = 5 | **EXACTLY INVARIANT** | expressible in k_5 alone; reads no radial operator |
| transport operator V_l2 - V_l1 as an algebraic object | **EXACTLY INVARIANT** | the correction 3A^2/4r^2 carries no l |
| pinhole = Sum V_max[1..5] vs the fitted 22.25 | **NUMERICALLY SHIFTED, AND IMPROVED** | -1.09% -> +0.36% |
| transport = mean off-diagonal vs the fitted 0.54 | **NUMERICALLY SHIFTED, AND IMPROVED** | +0.88% -> +0.70% |
| resistance = transport * ln(alpha_5/alpha_1) vs the fitted 0.14 | **NUMERICALLY SHIFTED, AND IMPROVED** | +0.49% -> -0.02% |
| the derived residuals reproduce the quark ladder | **INTERPRETATION CHANGED** | never established under either operator: composite 3.44% (legacy) and 3.78% (corrected) against the fitted lock's 1.61%; per-knob agreement is not ladder agreement |
| N = 466 as a stable integer of the residual-pinned fit | **NUMERICALLY SHIFTED** | still a compensator; N at PDG masses moves under the correction and the drift width stays O(50) |
| lepton resistance = 7pi/100 selected over 4(omega-1) | **INTERPRETATION CHANGED** | conclusion survives on proximity (+0.94% vs +2.50%) but the R_OUTER-bisection selector that chose it is reopened by PR #271 |
| R_OUTER = 1.26 from a single-sector fixed point | **INTERPRETATION CHANGED** | reopened by PR #271 and not restored here |
| R_OUTER = 1.26 bracketed by two independent sectors | **NEW, AND WEAK** | legacy puts 1.26 outside the lepton/quark bracket; corrected straddles it at 0.31 of the way across a 0.91% window — suggestive, not a derivation |

**Headline.** The quark residual sector survives the operator correction and improves; the lepton sector's closure did not, and the difference is elasticity: d ln m_s/d ln pinhole = +4.798 against the lepton's -17.471.

**Not re-run, and why.** The four quark shell-index axioms are expressible in `k₅ = 5` alone, the Hopf transport phase is topological, and the flavor-CP layer inherits the v3 eigenvalues by construction. None of them reads the radial scalar operator. Proximity is not dependence.

**Still open.** Nothing here derives `22.25`, `0.54` or `0.14` — each remains a fitted knob that a geometric quantity lands near. The composite gap (Q2) is the sharp successor target: a residual set that is individually right and jointly wrong points at a missing correlation between the three, not at three independent near-misses.
