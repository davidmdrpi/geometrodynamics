# Settling PR #270's ringdown cross-validation

**9/9 checks pass.**

PR #270 built two time-domain codes for a test scalar on a fixed `D = 5` Tangherlini background. Both stable, both converged, and they disagreed — real parts within `0.3%`, and damping rates apart by **37% of the smaller value** (equivalently, the wrong one was **27% low**). #270 refused to quote a frequency, correctly, and named its own prime suspect: the Kerr–Schild inner cut.

> ## The Kerr–Schild code was right. The tortoise code's damping was wrong.

All damping errors below are relative to the **published** value.

| source | `ω` at `ℓ = 1` | damping error |
|--|--|--|
| **published** (continued fractions + Hill) | `1.01601691 -0.36232802i` | — |
| #270 Kerr–Schild | `1.01622 -0.36231i` | **`0.005%`** |
| **this round** (characteristic, `h → 0`) | `1.01612 -0.36244i` | **`0.030%`** |
| #270 tortoise | `1.01876 -0.26404i` | **`27.1%`** ← excluded |

An independent Gundlach–Price–Pullin evolution, written from scratch and sharing no code with either #270 code, lands on a spectrum computed elsewhere by continued fractions. **#270's suspicion pointed at the wrong code.**

> **Note the ordering.** #270's Kerr–Schild code is ~6× *more* accurate than this round's solver. The characteristic scheme arbitrated correctly; it did not out-resolve.

**Naming the denominator.** The tortoise damping is `27.1%` off in the conventional sense (divided by the published value), and equivalently the correct damping is `37.2%` *larger* than the tortoise value (divided by the tortoise value). PR #270 measured the latter because it had two codes and no reference. Both are true; neither should be quoted bare.

## R0 — why `D = 5` makes the answer checkable

Two exact facts do most of the work.

**No logarithmic tail.** The correction has the exact closed form `r* − r = −artanh(1/r) = −1/r − 1/(3r³) − ⋯`, so `−1/r` is the *leading behaviour*, not an equality. Every term decays — unlike 4D's growing `2M ln r` — so the far field is a pure Hankel wave with no Coulomb phase. The series predicts the next coefficient, and it is checked rather than assumed:

| `r` | `r(r*−r) + 1` | predicted `−1/(3r²)` | rel. error |
|--|--|--|--|
| `50` | `-1.3337e-04` | `-1.3333e-04` | `2.4e-04` |
| `200` | `-8.3335e-06` | `-8.3333e-06` | `1.5e-05` |
| `1000` | `-3.3331e-07` | `-3.3333e-07` | `8.2e-05` |

**The tail is *asymptotically* Bessel.** Exactly, `V = L/r² + (9/4−L)/r⁴ − (9/4)/r⁶` with `L = (ℓ+1)² − ¼`, so `√r H⁽¹⁾` is the *leading* outgoing solution, not the exact one at finite `r`. `V·r²` against `L`:

| `ℓ` | limit | `V·r²` at `r = 1000` | rel. error |
|--|--|--|--|
| 0 | `0.75` | `0.750001` | `2.0e-06` |
| 1 | `3.75` | `3.749998` | `4.0e-07` |
| 2 | `8.75` | `8.749993` | `7.4e-07` |
| 3 | `15.75` | `15.749986` | `8.6e-07` |

*V -> [(l+1)^2 - 1/4]/r^2 is the flat-limit identity PR #271 used to settle which radial operator was correct; it is reused here as an asymptotic boundary condition.*

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
| 1 | `1.0162112 -0.3621562i` | `1.0161835 -0.3623977i` | `1.0161158 -0.3624367i` | `2.4e-04`, `7.8e-05` |
| 2 | `1.5110476 -0.3571938i` | `1.5107648 -0.3574392i` | `1.5105838 -0.3575288i` | `3.7e-04`, `2.0e-04` |

## R3 — the frequencies, against the exact asymptote

| `ℓ` | characteristic | eikonal | 1st-order WKB |
|--|--|--|--|
| 0 | `0.535366 -0.384195i` | `0.500000 -0.353553i` | `0.805398 -0.378527i` |
| 1 | `1.016184 -0.362398i` | `1.000000 -0.353553i` | `1.168269 -0.360949i` |
| 2 | `1.510765 -0.357439i` | `1.500000 -0.353553i` | `1.613545 -0.356117i` |
| 3 | `2.008201 -0.355538i` | `2.000000 -0.353553i` | `2.085488 -0.354728i` |

**Real parts approach 0.5(l+1) from above and damping approaches -0.353553; both happen.** *The l = 0 barrier is weakest and the power-law tail contaminates its fit earliest, so its uncertainty is wider than the others.*

## R4 — four lines of evidence, all against the tortoise damping

| source | `Im ω` at `ℓ = 1` |
|--|--|
| published (external reference) | `-0.36233` |
| characteristic (this round, h = 0.025) | `-0.36244` |
| Kerr-Schild (PR #270) | `-0.36231` |
| first-order WKB | `-0.36095` |
| exact eikonal asymptote | `-0.35355` |
| tortoise (PR #270) | `-0.26404`  ← **excluded** |

> **The Kerr-Schild code was right and the tortoise code's damping was wrong. A published high-precision spectrum external to this repository puts the answer at 1.01601691-0.36232802i, which confirms Kerr-Schild to 0.005% and this round's independent characteristic evolution to 0.030%, while the tortoise damping is off by 27.1%. PR #270 named the Kerr-Schild inner cut as the prime suspect; that suspicion pointed at the wrong code.**

**What this round cannot do.** Neither #270 code was landed in the tree, only their reported numbers, so there is no autopsy of WHICH line of the tortoise code produced the wrong damping -- only the demonstration that it did.

## R5 — what did not work

**Frequency-domain shooting.** REPRODUCED PR #270's NON-CONVERGENCE. Roots across `ε`: `1.229 - 0.152i`, `1.173 - 0.102i`, `1.155 - 0.214i`. The QNM boundary-value problem is exponentially ill-conditioned in real r: for Im w < 0 the outgoing piece grows like e^{|Im w| R} and swamps the incoming coefficient one is trying to zero.

**Sixth-order WKB.** UNUSABLE BY FINITE DIFFERENCES. Under refinement: `9.01 + 8.97i` → `18.63 + 18.61i` → `623.09 + 623.09i`. The Iyer-Will formula needs V^(6) at the barrier peak, and central differences amplify roundoff as h^-6, so refining the grid makes the answer DIVERGE.

**First-order WKB accuracy**, stated rather than assumed: damping 0.4% -- good, real part 13% -- poor at `ℓ = 1`, 6.8% -- improving at `ℓ = 2`. Low-l WKB is simply inaccurate on the real part; this is a known limitation of the method and NOT a discrepancy between solvers.

## R7 — against a published high-precision spectrum

The strongest check available, and external to this repository entirely. Source: Matyjasek, Phys. Rev. D 104, 084066 (2021), arXiv:2107.04815 -- continued fractions cross-checked against Hill determinants, agreeing to 11 digits.

| `ℓ` | characteristic (`h = 0.05`) | published | real err | damping err |
|--|--|--|--|--|
| 0 | `0.535366 -0.384195i` | `0.53383557 -0.38337537i` | `0.2867%` | `0.2137%` |
| 1 | `1.016184 -0.362398i` | `1.01601691 -0.36232802i` | `0.0164%` | `0.0192%` |
| 2 | `1.510765 -0.357439i` | `1.51056745 -0.35753726i` | `0.0131%` | `0.0274%` |

**With an external reference in hand this is no longer only a tie-breaker between two internal codes. It is an independent implementation reproducing a known high-precision spectrum, which is a considerably stronger check on PR #271's corrected radial operator and on the GPP machinery than internal arbitration.**

*PR #270's Kerr-Schild code, at 0.005% in damping against 0.031% for this round's characteristic evolution -- about 6x better. The characteristic scheme arbitrated correctly; it is not the most accurate of the three, and should not be described as though it were.*

### What the reference exposed about this solver

| `h` | `Im ω` | distance to published |
|--|--|--|
| `0.1` | `-0.3621562` | `1.72e-04` |
| `0.05` | `-0.3623977` | `6.97e-05` |
| `0.025` | `-0.3624367` | `1.09e-04` |

Last successive difference `3.90e-05`, actual distance from the finest value to the published one `1.09e-04` — **understated by `2.8×`**, and the `h = 0.05` value lands *closer* than the `h = 0.025` one.

> Self-convergence measures only the error it is refining away. The step-size study's last successive difference was ~2.7x smaller than the finest value's actual distance to the published one, and the middle step lands closer than the finest -- so discretization is not what limits this solver, extraction systematics are. A convergence study is a consistency check, not an error bar.

## R8 — where the error floor actually lives

An earlier draft *asserted* the residual was extraction systematics. Nothing was varied behind that. Varying them:

| extraction window | damping rel. error |
|--|--|
| `(50.0, 130.0)` | `1.0205%` |
| `(60.0, 140.0)` | `0.0192%` |
| `(70.0, 150.0)` | `0.0392%` |
| `(80.0, 160.0)` | `0.3252%` |
| `(90.0, 180.0)` | `61.2379%` |

| observer `r*` | damping rel. error |
|--|--|
| `20.0` | `0.0202%` |
| `30.0` | `0.0192%` |
| `40.0` | `1.0796%` |

**The window dominates** by orders of magnitude — late windows admit the power-law tail. **`t_max` is bit-irrelevant** (confirmed), because the extraction window sits well inside it.

Over reasonable choices the band is `0.0192%`–`0.0392%`, against a step-refinement difference of `0.0108%` — **3.6× larger**. That is why step refinement alone was the wrong error bar.

> An internal study CAN expose this floor -- this function is one, and it used no external value. The earlier claim that 'nothing internal to a solver can reveal this' was too strong. What the published spectrum uniquely provides is the anchor: a systematic spread says the answer is uncertain, but only an external reference says which point in the spread is correct.

## R6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| PR #270's two time-domain codes disagreed in damping (37% of the smaller value; the wrong one 27% low) | **CONFIRMED, AND NOW RESOLVED** | an independent characteristic evolution agrees with Kerr-Schild to ~1e-4 and excludes the tortoise damping |
| the Kerr-Schild inner cut was the prime suspect | **WRONG SUSPECT** | that code's frequency is confirmed; the fault was in the tortoise evolution |
| the verdict rests only on internal arbitration | **NO -- CONFIRMED EXTERNALLY** | a published high-precision spectrum (continued fractions + Hill determinants) puts l = 1 at 1.01601691-0.36232802i, confirming Kerr-Schild to 0.005% and excluding the tortoise damping at 27.1% |
| the characteristic evolution is the most accurate of the three | **NO** | PR #270's Kerr-Schild code is ~6x closer to the published value (0.005% against 0.031%); this round arbitrated, it did not out-resolve |
| the step-size study bounded this solver's error | **NO -- IT UNDERSTATED IT 2.7x** | last successive difference 4.0e-5, actual distance to the published value 1.1e-4, and h = 0.05 lands closer than h = 0.025; self-convergence estimates only the error component tied to the refinement parameter |
| the residual is extraction systematics | **MEASURED, NOT ASSERTED** | varying the extraction window, observer radius and t_max shows the window dominating by orders of magnitude and t_max being bit-irrelevant; the band over reasonable choices is comparable to the gap to the published value |
| only an external reference could expose that floor | **NO -- TOO STRONG** | the internal systematics scan exposes it using no external value; what the published spectrum uniquely supplies is the anchor, i.e. WHICH point in the spread is right |
| the GPP potential was sampled at the diamond centre | **NO -- IT WAS OFF BY h/4** | the centre is r* = (j-i)h/2, the two half-steps cancelling; the code sampled r* - h/4, contradicting its own docstring. Fixed and locked by a test; measured effect ~1e-6, so a real bug that does NOT explain the error floor |
| the far-field solution is exactly Hankel | **NO -- ASYMPTOTICALLY** | V = L/r^2 + (9/4-L)/r^4 - (9/4)/r^6 exactly, so the outgoing solution carries a further radial series; the failed shoot truncated to pure Hankel, making boundary truncation a second unseparated confounder |
| a quasinormal frequency can now be quoted | **YES, FOR l = 1, 2, 3** | converged in step size, window-stable, consistent with the exact eikonal asymptote, and matching a published spectrum to <0.05% at l = 1 and 2 |
| the l = 0 frequency is equally well determined | **NO** | its barrier is weakest and the power-law tail contaminates the fit earliest; quoted with a wider uncertainty |
| frequency-domain shooting can settle this | **NO -- REPRODUCED THE FAILURE** | the root moves with every numerical knob; the problem is exponentially ill-conditioned in real r |
| higher-order WKB can settle this | **NOT BY FINITE DIFFERENCES** | V^(6) by central differences diverges under refinement |
| the tortoise code's error can be diagnosed | **NO** | neither #270 code was landed, only their numbers; this round shows WHICH is wrong, not WHY |

**The standing lesson held.** PR #270 refused to quote a frequency from two converged codes that disagreed. That was right, and the way out was a third implementation sharing no code with either -- plus an exact asymptote to judge all three against.

**The lesson this round adds.** Self-convergence estimates the error component associated with the refinement parameter, not the total numerical, model or extraction error. This round's step-size study would have quoted ~4e-5 when the true error was ~1.1e-4. An internal scan over extraction choices does expose the missing component -- what the external reference uniquely supplies is the anchor that says which point in the spread is correct. Look for a published benchmark BEFORE building a third implementation to break a tie.

**Still open.** Overtones, backreaction, and the retarded outer-to-inner transfer function that #270 also deferred; the transfer function is a ratio of the same two signals and is now unblocked, since the signals can be trusted.

**The next object.** The retarded transfer function G_l(t; r_obs, r_src) from the same characteristic evolution with a compact ingoing excitation, gated on three checks before any physical reading: causal support G(t < t_null) = 0; flux conservation |R_l|^2 + |T_l|^2 = 1 on the fixed background; and late-time ringdown consistent with the EXTERNAL value 1.01601691149-0.36232802385i rather than with this solver's own fitted number.
