# Which of PR #270's two ringdown codes was right?

*Module `geometrodynamics/tangherlini/ringdown.py`, probe
`experiments/closure_ledger/ringdown_cross_validation_probe.py` (8/8), tests
`tests/test_ringdown.py` (32/32).*

---

## The answer

**The Kerr–Schild code was right. The tortoise code's damping was wrong.**

PR #270 built two horizon-penetrating time-domain codes for a test scalar on a
fixed `D = 5` Tangherlini background. Both were stable, both converged in their
own refinement parameter, and they disagreed at `ℓ = 1`: real parts within
`0.3 %`, damping rates apart by `37 %` of the smaller value. #270 declined to
quote a frequency — correctly, because *a converged number is not a correct
number* — and named its own prime suspect: the Kerr–Schild operator's inner cut.

That suspicion pointed at the wrong code.

Two independent things establish it. A characteristic evolution written here
from scratch, sharing no code with either; and — decisively — a **published
high-precision spectrum external to this repository entirely**: Matyjasek,
*Phys. Rev. D* **104**, 084066 (2021), [arXiv:2107.04815](https://arxiv.org/abs/2107.04815),
computed by continued fractions and cross-checked against Hill determinants,
the two agreeing to 11 digits.

The paper tabulates the scaled frequency `ω̃ = ω/T_H` with `T_H = 1/(2π)`; at
`ℓ = 1` it gives `6.38382253011 − 2.27657411582i`, so for `r_h = 1`

```
ω = 1.01601691149 − 0.36232802385i
```

All errors below are relative to that value.

| source | `ω` at `ℓ = 1` | damping error |
|--|--|--|
| **published** (CF / Hill) | `1.01601691 − 0.36232802i` | — |
| #270 Kerr–Schild | `1.01622 − 0.36231i` | **`0.005 %`** |
| **this round** (characteristic, `h → 0`) | `1.01612 − 0.36244i` | **`0.031 %`** |
| #270 tortoise | `1.01876 − 0.26404i` | **`27.1 %`** ← excluded |

**Note the ordering.** #270's Kerr–Schild code is about **six times more
accurate** than this round's characteristic evolution. The characteristic
scheme arbitrated correctly; it did not out-resolve, and describing it as the
best of the three would misread what it was for.

### What this round actually is

With an external reference in hand, this stops being merely a tie-breaker
between two internal codes. It is **an independent implementation reproducing a
known high-precision spectrum** — a considerably stronger check on PR #271's
corrected radial operator and on the GPP machinery than any amount of internal
arbitration could be.

### Naming the denominator

"`X %` off" is ambiguous until the denominator is stated, and earlier drafts of
this round quoted the two conventions in different places. Both, once:

| quantity | value | meaning |
|--|--|--|
| `\|Δ\| / \|published\|` | **`27.1 %`** | the tortoise result's relative error — *the number to quote* |
| `\|Δ\| / \|tortoise\|` | **`37.3 %`** | how much *larger* the correct damping is |

Both are true. `37.3 %` is what #270 measured, because it was comparing two
codes to each other with no reference available — precisely the situation this
round removes. `27.1 %` is the conventional relative error against truth.

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

**The tortoise correction is a decaying power, not a log.** The exact closed
form is

```
r* − r = ½ ln((r−1)/(r+1)) = −artanh(1/r) = −1/r − 1/(3r³) − 1/(5r⁵) − ⋯
```

so `−1/r` is the **leading asymptotic behaviour, not an exact equality** — an
earlier draft of this document said "exactly `−1/r`", which the series above
contradicts. What matters physically survives that correction untouched: *every*
term decays. In 4D the correction is `2M ln r`, which grows, and the far field
carries a Coulomb-like phase; here there is nothing to unwind.

Better still, the series **predicts the next coefficient**, so the test checks a
number rather than a trend — `r(r*−r) + 1` against `−1/(3r²)`:

| `r` | `r* − r` | `r(r*−r) + 1` | predicted `−1/(3r²)` | rel. error |
|--|--|--|--|--|
| 50 | `−2.000267e-02` | `−1.3337e-04` | `−1.3333e-04` | `2.4e-04` |
| 200 | `−5.000042e-03` | `−8.3335e-06` | `−8.3333e-06` | `1.5e-05` |
| 1000 | `−1.000000e-03` | `−3.3331e-07` | `−3.3333e-07` | `8.2e-05` |

A threshold at a single large radius would have passed for a wrong tail; so
would "it shrinks". Matching the predicted coefficient would not.

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

Those `|Δ|` are complex moduli; the damping-only differences used in the
comparison further down are smaller (`2.4e-4`, `4.0e-5` at `ℓ = 1`). **And this
table is exactly the thing the next section shows to be over-optimistic** — it
demonstrates that the scheme converges, which is not the same as bounding how
far the converged value sits from the truth.

The frequencies track the exact asymptote in the way a correct solver must
— real parts approaching `0.5(ℓ+1)` **from above**, damping approaching
`−0.353553`:

| `ℓ` | characteristic | eikonal | 1st-order WKB |
|--|--|--|--|
| 0 | `0.535380 − 0.384175i` | `0.500000 − 0.353553i` | `0.805398 − 0.378527i` |
| 1 | `1.016186 − 0.362395i` | `1.000000 − 0.353553i` | `1.168269 − 0.360949i` |
| 2 | `1.510765 − 0.357438i` | `1.500000 − 0.353553i` | `1.613545 − 0.356117i` |
| 3 | `2.008183 − 0.355527i` | `2.000000 − 0.353553i` | `2.085488 − 0.354728i` |

---

## Against the published spectrum

Three modes, computed independently here, against values obtained by continued
fractions and Hill determinants:

| `ℓ` | characteristic (this round) | published | real err | damping err |
|--|--|--|--|--|
| 0 | `0.535380 − 0.384175i` | `0.53383557 − 0.38337537i` | `0.2893 %` | `0.2086 %` |
| 1 | `1.016186 − 0.362395i` | `1.01601691 − 0.36232802i` | `0.0166 %` | `0.0184 %` |
| 2 | `1.510765 − 0.357438i` | `1.51056745 − 0.35753726i` | `0.0130 %` | `0.0278 %` |

`ℓ = 0` is the loosest by an order of magnitude, which independently vindicates
having given it a visibly wider extraction uncertainty rather than quoting it
alongside the others.

### What the reference exposed about this solver

This is the part that could not have been learned from inside.

| `h` | `Im ω` | distance to published |
|--|--|--|
| `0.1` | `−0.3621571` | `1.71e-04` |
| `0.05` | `−0.3623949` | `6.68e-05` |
| `0.025` | `−0.3624352` | `1.07e-04` |

The step-size study's last successive difference is `4.0e-5`. The finest value's
actual distance to the published one is `1.07e-4` — **2.7× larger** — and the
`h = 0.05` value lands *closer* to the truth than the `h = 0.025` one.

So discretization is **not** what limits this solver. Extraction systematics are:
window placement, Prony order, observer radius, finite `t_max`, power-law tail
contamination. Had this round quoted an error bar from its own refinement
sequence, that bar would have been almost three times too small.

> **Self-convergence measures only the error it is refining away. It is a
> consistency check, not an error bar.**

That is a sharper form of #270's own lesson, and nothing internal to a solver
can reveal it.

---

## Five lines of evidence, all on one side

```
published (external)       −0.36233    ← the reference
characteristic evolution   −0.36244    (this round, converged in h)
#270 Kerr–Schild           −0.36231
first-order WKB            −0.36095
exact eikonal asymptote    −0.35355
──────────────────────────────────
#270 tortoise              −0.26404    ← excluded
```

The five agreeing lines span `0.0089`. The tortoise damping sits `0.098` away,
on the other side of the eikonal asymptote, and shows no sign of approaching
it. Four of the five are independent of this round's solver entirely.

---

## Ledger

| claim | verdict | evidence |
|--|--|--|
| #270's two time-domain codes disagreed in damping | **CONFIRMED, AND NOW RESOLVED** | independent evolution agrees with Kerr–Schild, excludes tortoise |
| the Kerr–Schild inner cut was the prime suspect | **WRONG SUSPECT** | that code's frequency is confirmed; the fault was in the tortoise evolution |
| the verdict rests only on internal arbitration | **NO — CONFIRMED EXTERNALLY** | published CF/Hill spectrum confirms Kerr–Schild to `0.005 %`, excludes tortoise at `27.1 %` |
| the characteristic evolution is the most accurate of the three | **NO** | #270's Kerr–Schild is ~6× closer to the published value |
| the step-size study bounded this solver's error | **NO — UNDERSTATED 2.7×** | last `Δ` `4.0e-5`, true error `1.1e-4`; `h = 0.05` beats `h = 0.025` |
| a quasinormal frequency can now be quoted | **YES, for `ℓ = 1, 2, 3`** | converged, window-stable, consistent with the exact asymptote, matching published to `<0.05 %` |
| the `ℓ = 0` frequency is equally well determined | **NO** | `0.21 %` against published — an order of magnitude looser |
| frequency-domain shooting can settle this | **NO — reproduced the failure** | the root moves with every knob; exponentially ill-conditioned in real `r` |
| higher-order WKB can settle this | **NOT BY FINITE DIFFERENCES** | `V⁽⁶⁾` by central differences diverges under refinement |
| the tortoise code's error can be diagnosed | **NO** | neither #270 code was landed; this shows WHICH is wrong, not WHY |

---

## What this round adds to the standing list

**A converged number is not a correct number** — already on the list, and #270
applied it correctly by refusing to quote. What this round adds is the way
*out*: when two converged codes disagree, refining either one further is wasted
work. What breaks the tie is a third implementation that shares no code, plus a
reference none of them can influence.

**Name the property that makes the arbitrator immune to the suspected fault.**
The characteristic scheme was not chosen because it is more accurate — it
demonstrably is not, being ~6× worse than #270's Kerr–Schild code. It was chosen
because it applies no spatial boundary condition at all, and the suspected fault
was an inner cut. An arbitrator that could fail the same way would have settled
nothing.

**Self-convergence is a consistency check, not an error bar.** This round's own
step-size study would have quoted `~4e-5` when the true error was `~1.1e-4`.
Only an external reference could show that.

**Look for a published benchmark before building a third implementation.** Had
the external spectrum been found first, the arbitration would have been
unnecessary — the two #270 codes could have been scored against it directly. The
characteristic evolution earned its place anyway, because it is now validated
machinery for the transfer function rather than a one-off tie-breaker. But the
search order was wrong, and it is the cheap step.

---

## Scope, and what this round cannot do

**There is no autopsy.** Neither #270 code was landed in the tree — only their
reported numbers. This round establishes *which* number is right and by how much
the other is excluded. It cannot say which line of the unlanded tortoise code
produced the wrong damping.

**Fixed background.** A test scalar on an exact Tangherlini metric, no
backreaction. Fundamental (`n = 0`) mode only; no overtones.

---

## What did not work

Both are reported as negatives rather than quietly dropped, and the published
reference now supplies what both were reaching for.

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
reported the same failure; that diagnosis stands. The published continued-fraction
calculation is the frequency-domain reference this attempt was trying to be, and
there is no value in forcing the shoot to work merely to accumulate methods.

**Sixth-order WKB by finite differences is unusable.** The Iyer–Will formula
needs `V⁽⁶⁾` at the barrier peak, and central differences amplify roundoff as
`h⁻⁶`. Refining the grid makes the answer *diverge*: `9.01 → 18.63 → 623.09`.
First-order WKB is well conditioned and is what the tables above use, with its
accuracy **stated rather than assumed** — `0.4 %` on the damping at `ℓ = 1`
(good), `13 %` on the real part (poor), improving to `6.8 %` at `ℓ = 2`.
Low-`ℓ` WKB is simply not accurate on the real part; that is a known limitation
of the method, not a discrepancy between solvers, and the test asserts both
halves so the claim cannot rot.

---

## The next object

The retarded transfer function is now unblocked, and it is the right next
calculation — but it should be defined more sharply than #270 did:

```
G_ℓ^ret(t; r_obs, r_src)
```

from the same characteristic evolution with a compact **purely ingoing**
excitation, gated on three independent checks before any physical reading:

1. **Causal support** — `G(t < t_null) = 0`.
2. **Flux conservation** — `|R_ℓ(ω)|² + |T_ℓ(ω)|² = 1` on the fixed background.
3. **Late-time ringdown** consistent with `1.01601691149 − 0.36232802385i`.

Check 3 must use the **external** value, not this solver's own fitted number.
Scoring a solver against a frequency extracted from the same solver is not a
check, and this round has just shown that its own refinement sequence would have
been the wrong yardstick.

Once those pass, the kernel is trustworthy enough to compare against the rigid
exchange kernels used elsewhere in BAM — which asks whether the actual causal 5D
geometry produces the propagation structure the later model assumes. That is a
more interesting question than another QNM calculation.

**Separately blocked:** the nonlinear backreaction oracle, by the `C(v)` /
inner-flux issue, which is unrelated to anything here.
