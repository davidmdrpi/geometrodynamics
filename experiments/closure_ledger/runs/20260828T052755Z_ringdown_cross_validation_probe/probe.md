# Settling PR #270's ringdown cross-validation

**7/7 checks pass.**

PR #270 built two time-domain codes for a test scalar on a fixed `D = 5` Tangherlini background. Both stable, both converged, and they disagreed — real parts within `0.3%`, **damping rates differing by 37%**. #270 refused to quote a frequency, correctly, and named its own prime suspect: the Kerr–Schild inner cut.

> ## The Kerr–Schild code was right. The tortoise code's damping was wrong.

| source | `ω` at `ℓ = 1` | gap to this round |
|--|--|--|
| **this round** (characteristic, `h → 0`) | `1.01612 -0.36244i` | — |
| #270 Kerr–Schild | `1.01622 -0.36231i` | **`0.035%`** damping |
| #270 tortoise | `1.01876 -0.26404i` | **`37.3%`** damping |

An independent Gundlach–Price–Pullin evolution, written from scratch and sharing no code with either, agrees with Kerr–Schild to `~1e-4` and misses the tortoise damping by 37%. **#270's suspicion pointed at the wrong code.**

## R0 — why `D = 5` makes the answer checkable

Two exact facts do most of the work.

**No logarithmic tail.** `r* − r` at `r = 50, 200, 1000` is `2.0e-02, 5.0e-03, 1.0e-03` — unlike 4D, `r* → r`, so the far field is a pure Hankel wave with no Coulomb phase.

**The tail is exactly Bessel.** `V·r²` against `(ℓ+1)² − ¼`:

| `ℓ` | limit | `V·r²` at `r = 1000` | rel. error |
|--|--|--|--|
| 0 | `0.75` | `0.750001` | `2.0e-06` |
| 1 | `3.75` | `3.749998` | `4.0e-07` |
| 2 | `8.75` | `8.749993` | `7.4e-07` |
| 3 | `15.75` | `15.749986` | `8.6e-07` |

*V -> [(l+1)^2 - 1/4]/r^2 is the flat-limit identity PR #271 used to settle which radial operator was correct; it is reused here as a boundary condition.*

## R1 — an exact asymptote to judge every solver against

| quantity | value | exact? |
|--|--|--|
| `r_ph²` | `2.000000000000000` | **2 exactly** |
| `Ω_c` | `0.500000000000000` | **1/2 exactly** |
| `λ` | `0.707107` | **1/√2** |

> `omega -> 0.5 (l+1) - 0.353553i  as l -> infinity`

## R2 — the independent solver converges

Gundlach–Price–Pullin on the null grid `u = t − r*`, `v = t + r*`. **No spatial boundary conditions are applied or needed** — the domain of dependence is the null diamond, so the horizon and infinity are limits rather than boundaries. That is exactly why this construction is immune to the excision question that both the frequency-domain shoot and the Kerr–Schild inner cut raise.

| `ℓ` | `h = 0.1` | `h = 0.05` | `h = 0.025` | successive `\|Δ\|` |
|--|--|--|--|--|
| 1 | `1.0162114 -0.3621571i` | `1.0161856 -0.3623949i` | `1.0161189 -0.3624352i` | `2.4e-04`, `7.8e-05` |
| 2 | `1.5110469 -0.3571935i` | `1.5107645 -0.3574379i` | `1.5105884 -0.3575318i` | `3.7e-04`, `2.0e-04` |

## R3 — the frequencies, against the exact asymptote

| `ℓ` | characteristic | eikonal | 1st-order WKB |
|--|--|--|--|
| 0 | `0.535380 -0.384175i` | `0.500000 -0.353553i` | `0.805398 -0.378527i` |
| 1 | `1.016186 -0.362395i` | `1.000000 -0.353553i` | `1.168269 -0.360949i` |
| 2 | `1.510765 -0.357438i` | `1.500000 -0.353553i` | `1.613545 -0.356117i` |
| 3 | `2.008183 -0.355527i` | `2.000000 -0.353553i` | `2.085488 -0.354728i` |

**Real parts approach 0.5(l+1) from above and damping approaches -0.353553; both happen.** *The l = 0 barrier is weakest and the power-law tail contaminates its fit earliest, so its uncertainty is wider than the others.*

## R4 — four lines of evidence, all against the tortoise damping

| source | `Im ω` at `ℓ = 1` |
|--|--|
| characteristic (this round) | `-0.36244` |
| Kerr-Schild (PR #270) | `-0.36231` |
| first-order WKB | `-0.36095` |
| exact eikonal asymptote | `-0.35355` |
| tortoise (PR #270) | `-0.26404`  ← **excluded** |

> **The Kerr-Schild code was right and the tortoise code's damping was wrong. An independent characteristic evolution sharing no code with either agrees with Kerr-Schild to ~1e-4 and misses the tortoise damping by 37%. PR #270 named the Kerr-Schild inner cut as the prime suspect; that suspicion pointed at the wrong code.**

**What this round cannot do.** Neither #270 code was landed in the tree, only their reported numbers, so there is no autopsy of WHICH line of the tortoise code produced the wrong damping -- only the demonstration that it did.

## R5 — what did not work

**Frequency-domain shooting.** REPRODUCED PR #270's NON-CONVERGENCE. Roots across `ε`: `1.229 - 0.152i`, `1.173 - 0.102i`, `1.155 - 0.214i`. The QNM boundary-value problem is exponentially ill-conditioned in real r: for Im w < 0 the outgoing piece grows like e^{|Im w| R} and swamps the incoming coefficient one is trying to zero.

**Sixth-order WKB.** UNUSABLE BY FINITE DIFFERENCES. Under refinement: `9.01 + 8.97i` → `18.63 + 18.61i` → `623.09 + 623.09i`. The Iyer-Will formula needs V^(6) at the barrier peak, and central differences amplify roundoff as h^-6, so refining the grid makes the answer DIVERGE.

**First-order WKB accuracy**, stated rather than assumed: damping 0.4% -- good, real part 13% -- poor at `ℓ = 1`, 6.8% -- improving at `ℓ = 2`. Low-l WKB is simply inaccurate on the real part; this is a known limitation of the method and NOT a discrepancy between solvers.

## R6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| PR #270's two time-domain codes disagreed by 37% in damping | **CONFIRMED, AND NOW RESOLVED** | an independent characteristic evolution agrees with Kerr-Schild to ~1e-4 and excludes the tortoise damping |
| the Kerr-Schild inner cut was the prime suspect | **WRONG SUSPECT** | that code's frequency is confirmed; the fault was in the tortoise evolution |
| a quasinormal frequency can now be quoted | **YES, FOR l = 1, 2, 3** | converged in step size, window-stable, and consistent with the exact eikonal asymptote |
| the l = 0 frequency is equally well determined | **NO** | its barrier is weakest and the power-law tail contaminates the fit earliest; quoted with a wider uncertainty |
| frequency-domain shooting can settle this | **NO -- REPRODUCED THE FAILURE** | the root moves with every numerical knob; the problem is exponentially ill-conditioned in real r |
| higher-order WKB can settle this | **NOT BY FINITE DIFFERENCES** | V^(6) by central differences diverges under refinement |
| the tortoise code's error can be diagnosed | **NO** | neither #270 code was landed, only their numbers; this round shows WHICH is wrong, not WHY |

**The standing lesson held.** PR #270 refused to quote a frequency from two converged codes that disagreed. That was right, and the way out was a third implementation sharing no code with either -- plus an exact asymptote to judge all three against.

**Still open.** Overtones, backreaction, and the retarded outer-to-inner transfer function that #270 also deferred; the transfer function is a ratio of the same two signals and is now unblocked, since the signals can be trusted.
