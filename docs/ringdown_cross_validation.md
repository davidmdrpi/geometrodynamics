# Which of PR #270's two ringdown codes was right?

*Module `geometrodynamics/tangherlini/ringdown.py`, probe
`experiments/closure_ledger/ringdown_cross_validation_probe.py` (7/7), tests
`tests/test_ringdown.py` (24/24).*

---

## The answer

**The Kerr–Schild code was right. The tortoise code's damping was wrong.**

PR #270 built two horizon-penetrating time-domain codes for a test scalar on a
fixed `D = 5` Tangherlini background. Both were stable, both converged in their
own refinement parameter, and they disagreed: real parts within `0.3 %` at
`ℓ = 1`, damping rates differing by `37 %`. #270 declined to quote a frequency
— correctly, because *a converged number is not a correct number* — and named
its own prime suspect: the Kerr–Schild operator's inner cut.

That suspicion pointed at the wrong code.

| source | `ω` at `ℓ = 1` | gap to this round (damping) |
|--|--|--|
| **this round** — characteristic, `h → 0` | `1.01612 − 0.36244i` | — |
| #270 Kerr–Schild | `1.01622 − 0.36231i` | **`0.035 %`** |
| #270 tortoise | `1.01876 − 0.26404i` | **`37.3 %`** |

The `37.3 %` this round measures against the tortoise number is the same `37 %`
#270 measured between its own two codes — which is the consistency check that
the third implementation is landing on one of them and not somewhere new.

---

## Why a third code, and why this one

Two converged codes that disagree cannot be adjudicated by refining either of
them further. The only way out is an implementation that shares no code with
either, plus an exact asymptote that all three can be judged against
independently of any of them.

The third implementation is a **Gundlach–Price–Pullin characteristic
evolution** on the null grid `u = t − r*`, `v = t + r*`:

```
ψ_N = ψ_W + ψ_E − ψ_S − (h²/8) V_c (ψ_W + ψ_E)
```

The decisive structural property: **no spatial boundary conditions are applied
or needed.** The domain of dependence of the initial null diamond is the
diamond itself, so the horizon and infinity are reached only as limits, never
as boundaries. That is precisely why this construction is immune to the
excision question that both the frequency-domain shoot and the Kerr–Schild
inner cut raise — the question that #270 suspected. It cannot be wrong for
that reason, so it can arbitrate it.

`V` at a diamond centre depends only on `r* = h(j − i)/2`, so the potential is
a one-dimensional lookup indexed by `j − i`, and the inner `j`-recursion is
linear and first order,

```
a_j = A_j a_{j−1} + B_j        A_j = 1 − (h²/8)V_j
```

which is evaluated in closed form by `cumprod`/`cumsum` rather than a Python
loop. Without that, the finest run (`t_max = 400`, `h = 0.025`, so
`16000 × 16000`) is ~`2.6e8` scalar iterations; with it the whole sweep is
seconds.

---

## What makes `D = 5` checkable at all

Two exact facts, and they are why this round can be settled where a generic
background would leave it open.

**The tortoise correction is a decaying power, not a log.**
`r* = r + ½ ln((r−1)/(r+1))`, whose correction is exactly `−1/r`. In 4D it is
`2M ln r`, which grows, and the far field carries a Coulomb-like phase. Here
`r* → r`, so the far field is a pure Hankel wave with nothing to unwind.

| `r` | `r* − r` | `−1/r` |
|--|--|--|
| 50 | `−2.000267e-02` | `−2.0e-02` |
| 200 | `−5.000042e-03` | `−5.0e-03` |
| 1000 | `−1.000000e-03` | `−1.0e-03` |

The test does not check the value at one radius — it checks the *power law*,
that `r(r* − r) + 1` falls like `1/r²` (`1.33e-4 → 8.33e-6 → 3.33e-7`). A
threshold at a single large radius would have passed for a wrong tail too.

**The potential's tail is exactly Bessel.** `V → [(ℓ+1)² − ¼]/r²`, so the
outgoing solution is exactly `√r H⁽¹⁾_{ℓ+1}(ωr)`.

| `ℓ` | `(ℓ+1)² − ¼` | `V·r²` at `r = 1000` | rel. error |
|--|--|--|--|
| 0 | `0.75` | `0.750001` | `2.0e-06` |
| 1 | `3.75` | `3.749998` | `4.0e-07` |
| 2 | `8.75` | `8.749993` | `7.4e-07` |
| 3 | `15.75` | `15.749986` | `8.6e-07` |

This is the *same identity* PR #271 used to settle which radial operator was
correct. It is reused here as a boundary condition, so the operator this round
evolves is #271's corrected one, taken from `dynamics` rather than restated.

---

## The exact asymptote every solver is judged against

The photon sphere solves `d/dr (f/r²) = 0`, which at `n = 3` gives `r_ph² = 2`
**exactly**, and with it

| quantity | value | exact |
|--|--|--|
| `r_ph²` | `2.000000000000000` | `2` |
| `Ω_c = √f(r_ph)/r_ph` | `0.500000000000000` | `1/2` |
| `λ` | `0.707107` | `1/√2` |

```
ω → Ω_c(ℓ+1) − iλ/2 = 0.5(ℓ+1) − 0.353553i
```

Nothing in this asymptote comes from any of the three codes. It is the fixed
point all of them must approach, and it is what turns "two numbers disagree"
into "one of them is on the wrong side of an exact limit."

---

## The solver earns its number

Converged in step size, at two `ℓ`:

| `ℓ` | `h = 0.1` | `h = 0.05` | `h = 0.025` | successive `\|Δ\|` |
|--|--|--|--|--|
| 1 | `1.0162114 − 0.3621571i` | `1.0161856 − 0.3623949i` | `1.0161189 − 0.3624352i` | `2.4e-04`, `7.8e-05` |
| 2 | `1.5110469 − 0.3571935i` | `1.5107645 − 0.3574379i` | `1.5105884 − 0.3575318i` | `3.7e-04`, `2.0e-04` |

and the frequencies track the exact asymptote in the way a correct solver must
— real parts approaching `0.5(ℓ+1)` **from above**, damping approaching
`−0.353553`:

| `ℓ` | characteristic | eikonal | 1st-order WKB |
|--|--|--|--|
| 0 | `0.535380 − 0.384175i` | `0.500000 − 0.353553i` | `0.805398 − 0.378527i` |
| 1 | `1.016186 − 0.362395i` | `1.000000 − 0.353553i` | `1.168269 − 0.360949i` |
| 2 | `1.510765 − 0.357438i` | `1.500000 − 0.353553i` | `1.613545 − 0.356117i` |
| 3 | `2.008183 − 0.355527i` | `2.000000 − 0.353553i` | `2.085488 − 0.354728i` |

---

## Four lines of evidence, all on one side

```
characteristic evolution   −0.36244    (this round, converged in h)
#270 Kerr–Schild           −0.36231
first-order WKB            −0.36095
exact eikonal asymptote    −0.35355
──────────────────────────────────
#270 tortoise              −0.26404    ← excluded
```

The four agreeing lines span `0.0089`. The tortoise damping sits `0.098` away,
on the other side of the eikonal asymptote, and shows no sign of approaching
it. Three of the four are independent of this round's solver entirely.

---

## What this round cannot do

**There is no autopsy.** Neither #270 code was landed in the tree — only their
reported numbers. This round establishes *which* number is right and by how
much the other is excluded. It cannot say which line of the unlanded tortoise
code produced the wrong damping.

**Scope.** A test scalar on a **fixed** exact Tangherlini background, no
backreaction. Fundamental (`n = 0`) mode only; no overtones. `ℓ = 0` is quoted
with a wider uncertainty than the rest because its barrier is weakest and the
power-law tail contaminates its fit earliest.

---

## What did not work

Both are reported as negatives rather than quietly dropped.

**Frequency-domain shooting reproduced #270's non-convergence.** Matching an
ingoing Frobenius series at the horizon to `√r H⁽¹⁾` at large `r` gives a root
that moves with every numerical knob:

| knob | roots |
|--|--|
| horizon offset `ε` | `1.229 − 0.152i`, `1.173 − 0.102i`, `1.155 − 0.214i` |
| matching radius | `1.204 − 0.209i`, `1.173 − 0.102i`, `1.166 − 0.105i` |

This is not a bug to be found. The QNM boundary-value problem is exponentially
ill-conditioned in real `r`: for `Im ω < 0` the outgoing piece grows like
`e^{|Im ω|R}` and swamps the incoming coefficient one is trying to zero. #270
reported the same failure; that diagnosis stands, and it is the second reason
the characteristic scheme is the right tool.

**Sixth-order WKB by finite differences is unusable.** The Iyer–Will formula
needs `V⁽⁶⁾` at the barrier peak, and central differences amplify roundoff as
`h⁻⁶`. Refining the grid makes the answer *diverge*: `9.01 → 18.63 → 623.09`.
First-order WKB is well conditioned and is what the tables above use, with its
accuracy **stated rather than assumed** — `0.4 %` on the damping at `ℓ = 1`
(good), `13 %` on the real part (poor), improving to `6.8 %` at `ℓ = 2`.
Low-`ℓ` WKB is simply not accurate on the real part; that is a known limitation
of the method, not a discrepancy between solvers, and the test asserts both
halves of that so the claim cannot rot.

---

## Ledger

| claim | verdict | evidence |
|--|--|--|
| #270's two time-domain codes disagreed by `37 %` in damping | **CONFIRMED, AND NOW RESOLVED** | independent evolution agrees with Kerr–Schild to `~1e-4`, excludes tortoise |
| the Kerr–Schild inner cut was the prime suspect | **WRONG SUSPECT** | that code's frequency is confirmed; the fault was in the tortoise evolution |
| a quasinormal frequency can now be quoted | **YES, for `ℓ = 1, 2, 3`** | converged in step size, window-stable, consistent with the exact asymptote |
| the `ℓ = 0` frequency is equally well determined | **NO** | weakest barrier, earliest tail contamination; wider uncertainty |
| frequency-domain shooting can settle this | **NO — reproduced the failure** | the root moves with every knob; exponentially ill-conditioned in real `r` |
| higher-order WKB can settle this | **NOT BY FINITE DIFFERENCES** | `V⁽⁶⁾` by central differences diverges under refinement |
| the tortoise code's error can be diagnosed | **NO** | neither #270 code was landed; this shows WHICH is wrong, not WHY |

---

## Two things learned, added to the standing list

**A converged number is not a correct number** — already on the list, and #270
applied it correctly by refusing to quote. What this round adds is the way
*out*: when two converged codes disagree, refining either one further is wasted
work. What breaks the tie is a third implementation that shares no code, plus
an exact limit none of the three can influence.

**Name the property that makes the arbitrator immune to the suspected fault.**
The characteristic scheme was not chosen because it is more accurate. It was
chosen because it applies no spatial boundary condition at all, and the
suspected fault was an inner cut. An arbitrator that could fail the same way
would have settled nothing.

---

## Still open

Overtones, backreaction, and the retarded outer-to-inner transfer function
#270 also deferred. The transfer function is a ratio of the same two signals,
so it is now unblocked — the signals can be trusted.
