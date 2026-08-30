# The retarded outer→inner transfer kernel on the `D = 5` background

*Module `geometrodynamics/tangherlini/transfer_kernel.py`, probe
`experiments/closure_ledger/transfer_kernel_probe.py` (14/14), tests
`tests/test_transfer_kernel.py`.*

---

## The answer

**The transfer is not rigid, and not marginally.**

PR #270 deferred this object because it is a ratio of two signals it could not
trust. PR #274 settled which signal was right, against a published spectrum.
This round builds the kernel itself.

A wave sent in from the far region reaches the horizon filtered. That filter is
the transmission amplitude `T_ℓ(ω)`; in time it is a convolution kernel,

```
ψ_transmitted(v) = (K_ℓ ⋆ ψ_incident)(v) ,     v = t + r*
```

Because `T_ℓ(ω) → 1` at high frequency — fast waves fly over the barrier — the
kernel splits into an instantaneous part and a memory:

```
K_ℓ(t) = δ(t) + K_ℓ^reg(t)
```

| quantity | measured | meaning |
|--|--|--|
| `∫K_reg dt` | `−0.999996` | exact value **`−1`**, a sum rule, not a fit |
| `∫\|K_reg\| dt` | `2.03` | against the `δ`'s weight of `1` |
| `T` at lowest bin `ω = 0.00488` | `4.10e-07` | `T(0) = 0` exactly; this is the lowest **sampled** bin |

A rigid exchange kernel is `δ(t)`, possibly with a delay: whatever enters leaves
undistorted, and a static signal passes perfectly. The real geometry **blocks a
static signal completely**, and does so *entirely* through the memory term.

**The durable statement is the pair of limits, not the mass.** With the chosen
asymptotic normalisation `T(∞) = 1` while `T(0) = 0`; therefore the transfer
function cannot be a rigid phase or delay, and the regular memory must carry
signed weight `−1`. That is analytic and depends on no measured number.

`∫K dt = T(ω = 0)`, and `T(0) = 0` because the centrifugal barrier reflects
zero-frequency waves completely. With `K = δ + K_reg` that forces

```
∫ K_reg dt = −1     exactly
```

so the memory does not merely modify the instantaneous term — it **cancels it
exactly** at zero frequency. In absolute mass the memory is about twice the
delta — a quantitative extra, not the argument.

**The mass is quoted to three significant figures, which is what it earns.**
Across the subtraction parameter, `(ω_max, count)` and the time quadrature the
worst spread is `~2e-3`; quoting `2.0309` would claim precision the knobs do not
support. See the convergence study below.

**Scope of that claim.** This is a statement about the transfer kernel of a test
scalar on a fixed `D = 5` Tangherlini background, per angular channel. It says
what the causal geometry does. Whether any particular BAM exchange kernel is
*intended* to approximate this object is a modelling question this round does
not settle.

---

## Why the frequency domain works here, having failed twice

PR #270 and PR #274 both failed to converge a quasinormal frequency by shooting
in real `r`. #274 recorded the reason: for `Im ω < 0` the outgoing solution
grows like `e^{|Im ω|R}` and swamps the coefficient being zeroed.

**That objection does not apply to this calculation.** Here `ω` is *real*, so
the **asymptotic** in and out waves are unit-flux and the coefficient being
extracted is not swamped by an exponentially larger companion. This is not a
repaired version of the shoot that failed — it is a different, well-posed
problem that happens to use the same background.

**Scoped honestly:** unit modulus is a property of the asymptotic normalisation,
*not* of the propagation everywhere. Under the barrier `k = √(ω² − V)` is
imaginary and the local propagator is hyperbolic, so one solution does grow
relative to the other across the forbidden region. The correct claim is **well
conditioned on the tested real-frequency range** — evidenced by unitarity and
step refinement — not structural immunity as `ω → 0` or `ℓ` grows.

The evidence is unitarity, which is imposed nowhere and therefore measures the
computation rather than decorating it:

| `ω` | `\|T\|` | `\|R\|` | `\|R\|²+\|T\|²−1` |
|--|--|--|--|
| `0.1` | `0.00079553` | `0.99999968` | `+1.1e-15` |
| `0.7` | `0.19749250` | `0.98030440` | `+8.9e-15` |
| `1.0` | `0.64441842` | `0.76467307` | `+3.0e-14` |
| `1.5` | `0.98880233` | `0.14923117` | `−3.9e-13` |
| `3.0` | `0.99999977` | `0.00068283` | `+1.1e-13` |
| `30` | `1.00000000` | `0.00000000` | `+6.3e-13` |

The solver is a piecewise-constant transfer matrix, vectorised over `ω` so the
spatial loop is shared across all frequencies. It is second order in the step
(`|T|` at `ω = 1`: `0.64442644 → 0.64442002 → 0.64441842`, differences `6.4e-6`,
`1.6e-6`), and it agrees with an independent `solve_ivp` integration to `~1e-6`.

---

## Three exact anchors, derived not recalled

### 1. The barrier peak at `ℓ = 1` is rational

In `x = 1/r²` the potential is a cubic, `V = Lx + (9/4 − L)x² − (9/4)x³`, so the
peak solves `27/4 x² − 2(9/4 − L)x − L = 0`. At `ℓ = 1` the discriminant is a
perfect square (`110.25 = 10.5²`), giving `x = 5/9` and

```
r² = 9/5 ,      V_max = 100/81       (exactly)
```

| `ℓ` | `r²` at peak | `r*` | `V_max` |
|--|--|--|--|
| 0 | `1.605551` | `0.1978` | `0.505383744` |
| 1 | **`1.800000`** | `0.3792` | **`1.234567901`** |
| 2 | `1.893190` | `0.4541` | `2.476709048` |
| 3 | `1.935691` | `0.4862` | `4.223427135` |

`ℓ = 1` is the **only** non-negative integer `ℓ` that is rational, and this is a
**theorem** rather than an inference from three spot checks. With `m = ℓ+1` the
discriminant is `(16m⁴ + 28m² + 73)/4`, so rationality asks when
`16m⁴ + 28m² + 73` is a perfect square. For `m ≥ 4`,

```
(4m² + 3)²  <  16m⁴ + 28m² + 73  <  (4m² + 4)²
```

— the left inequality is `4m² + 64 > 0`, the right is `4m² > 57`, i.e. `m ≥ 4`.
Strictly between consecutive squares, hence never square. That leaves
`m = 1, 2, 3`, giving `117`, `441 = 21²`, `1621` — only `m = 2`. (The `1621` is
the same surd that appears in the `ℓ = 2` row above.) It is also the mode #270's
two codes disagreed on. And the peak
radius increases toward `r² = 2`, the photon sphere PR #274 pinned exactly, which
is the consistency bridge between the two rounds.

### 2. The integrated potential is elementary

`dr* = dr/f` cancels the overall `f` in `V`, so

```
∫ V_ℓ dr*  =  ∫₁^∞ (L/r² + (9/4)/r⁴) dr  =  L + 3/4  =  ℓ(ℓ+2) + 3/2
```

| `ℓ` | exact | quadrature (truncated at `r* = 150`) | predicted missing tail |
|--|--|--|--|
| 0 | `1.500` | `1.4950` | `0.0050` |
| 1 | `4.500` | `4.4750` | `0.0250` |
| 2 | `9.500` | `9.4417` | `0.0583` |
| 3 | `16.500` | `16.3950` | `0.1050` |

The deficit matches the predicted tail `L/r*_out` in **every** row, so the domain
truncation is accounted for rather than assumed small — the concern PR #274's
review raised about the failed shoot, handled quantitatively here.

### 3. Hence the kernel's high-frequency phase is closed-form

The leading eikonal phase through the barrier is `−(1/2ω)∫V dr*`, so

```
T_ℓ(ω) → exp(−i c_ℓ/ω) ,    c_ℓ = (ℓ(ℓ+2) + 3/2)/2 ,    c₁ = 9/4  exactly
```

Verified below against the fitted asymptote of the computed spectrum, which
agrees to `1.2e-4` relative and — with the Jost outer condition — is
**independent of the outer edge**. (Under the earlier plane-wave matching it was not: it drifted
`2.2365 → 2.2427 → 2.2458 → 2.2474` at `r*_out = 150, 300, 600, 1200`, the
deficit halving as the edge doubled, exactly as the truncated `L/(2r*_out)`
predicts. That drift is what the Jost condition removes.)

**This is what makes the kernel computable at all.** `T − 1 ~ −i c/ω` decays too
slowly to transform numerically. Knowing `c` exactly lets it be removed
analytically rather than windowed away:

```
A(ω) = −i c/(ω + ia)      →      −c e^{−at} θ(t)
```

whose only pole is at `ω = −ia`, in the lower half plane, so the subtraction
cannot itself introduce an acausal piece. The remainder decays like `1/ω²` and
transforms cleanly.

**The exact `c` is what gets subtracted, not a fit.** If
`T − 1 = −i c_exact/ω + O(ω⁻²)` and a fitted `c_fit ≠ c_exact` is subtracted, the
remainder retains `−i(c_exact − c_fit)/ω` and still falls only as `1/ω` — the
entire purpose of the subtraction, silently forfeited.

### The estimator, and a correction

An earlier draft measured `c` with `Re[iω(T−1)]` and reported a uniform `0.047%`
shortfall, attributing it to the `1/r⁴` and `1/r⁶` parts of `V` that a purely
centrifugal boundary condition cannot capture. **That attribution was wrong.**
The estimator is biased: even for the *exact* asymptote `T = e^{−ic/ω}`,

```
Re[iω(e^{−ic/ω} − 1)] = ω sin(c/ω) = c − c³/(6ω²) + O(ω⁻⁴)
```

At `c = 2.25` over the sampled band `ω ≈ 36–40` that deterministic bias is
`−1.3e-3` — the entire size of the reported deficit. The unbiased estimator
`c(ω) = −ω arg T(ω)` is exact for the toy asymptote:

| outer edge | `−ω arg T` (unbiased) | dev | `Re[iω(T−1)]` (biased) | dev |
|--|--|--|--|--|
| `150` | `2.25026029` | `+2.60e-04` | `2.24894513` | `−1.06e-03` |
| `300` | `2.25025990` | `+2.60e-04` | `2.24894474` | `−1.06e-03` |
| `600` | `2.25025985` | `+2.60e-04` | `2.24894469` | `−1.06e-03` |

The true residual is `+2.6e-4` — four times smaller and of the **opposite
sign**. No physical attribution is offered for it here.

**The tell was visible in the old data.** The biased estimator's deficit is
*itself* independent of the outer edge, to `4e-9`. A truncation effect moves
with the edge; a fixed-form bias does not. That should have prompted the check.

### And the remainder really does fall like `1/ω²`

Checking `c` is necessary but not sufficient — the property the transform needs
is the remainder's decay. `ω²|S − A|` stays bounded:

| `a` | `ω=5` | `ω=10` | `ω=20` | `ω=40` |
|--|--|--|--|--|
| `0.5` | `1.422` | `1.410` | `1.407` | `1.406` |
| `1.0` | `0.380` | `0.309` | `0.289` | `0.283` |
| `2.0` | `1.842` | `1.934` | `1.960` | `1.967` |

The plateau differs with `a` because a different `A` redistributes the `O(1/ω²)`
tail between the analytic and numerical pieces. The *kernel* must not depend on
`a` at all — checked separately below.

---

## The outer boundary condition is not plane waves

Matching to free `e^{±iωr*}` at a finite outer edge assumes `ωr* ≫ ν`. At the
low-frequency end of a kernel grid that is badly false, and it is the **worst**
place to be sloppy: the lowest bin sets the DC end of the inverse transform, and
therefore the numerical realisation of the `−1` sum rule.

At `ω = 40/(2·4096) ≈ 0.00488`:

| quantity | value |
|--|--|
| outer turning scale `√L/ω` | `396.6` |
| `V` at `r* = 150` | `1.67e-04` |
| `ω²` | `2.38e-05` |

`V` at the edge still **exceeds** `ω²` — that bin sits inside the centrifugal
tail, where plane waves are simply the wrong basis.

The fix is the exact solutions of the tail. Asymptotically the `r*` equation is
`ψ'' + (ω² − L/r*²)ψ = 0` with `L = ν² − ¼`, `ν = ℓ+1`, whose exact solutions are
`√x H^{(1,2)}_ν(x)`, `x = ωr*`. Normalised by `e^{±i(νπ/2 + π/4)}` they reduce to
`e^{±ix}`, so the high-frequency convention is untouched and only the
low-frequency matching moves. Their Wronskian is exactly `−2i`, which the tests
check rather than assume.

| outer edge | `\|T\|` at `ω=0.00488` | `ω=0.02` | `ω=0.1` |
|--|--|--|--|
| `150` | `4.099913e-07` | `1.393649e-05` | `7.955278e-04` |
| `300` | `4.100095e-07` | `1.393644e-05` | `7.955282e-04` |
| `600` | `4.100145e-07` | `1.393644e-05` | `7.955282e-04` |

Relative spread across outer edges `5.6e-05`, `3.6e-06`, `4.3e-07` — **converged**,
where plane-wave matching drifted. And `T(ω→0)` fell from `1.73e-06` to
`4.10e-07`, four times closer to its exact value of zero.

---

## The three gates

### Gate 1 — causal support

`K(t) = 0` for `t < 0`, by construction of the retarded kernel.

| `t` | `K_reg(t)` |
|--|--|
| `−300` | `+1.67e-07` |
| `−100` | `+4.67e-07` |
| `−20` | `−2.13e-06` |
| `−5` | `+8.26e-06` |
| `−1` | `−3.40e-05` |
| `−0.5` | `−1.03e-04` |

Worst `1.03e-04` next to the front, `~9e-07` away from it.

### Gate 2 — flux conservation

`|R|² + |T|² = 1` to `6.3e-13`. See the table above.

### Gate 3 — the published ringdown

Fitting the kernel's *own* ringdown against the external published value
`1.01601691149 − 0.36232802385i`:

| `dt` | window | fitted `ω` | real err | damping err |
|--|--|--|--|--|
| `0.2` | `(3, 14)` | `1.027952 − 0.376104i` | `1.175%` | `3.802%` |
| `0.2` | `(5, 18)` | `1.023143 − 0.363089i` | `0.701%` | `0.210%` |
| `0.3` | `(3, 14)` | `1.017028 − 0.363837i` | `0.100%` | `0.416%` |
| `0.3` | `(5, 18)` | `1.015390 − 0.361124i` | `0.062%` | `0.332%` |

Band across nine settings: real part `0.062%`–`1.17%`, damping `0.11%`–`3.80%`.
Following PR #274, the **band is the honest statement** rather than the best row
— extraction choices move the fit, and reporting only the best would repeat the
error that round measured. The comparison is with the published value, never
with this repository's own solver: scoring a kernel against a frequency
extracted from the same machinery would not be a check.

---

## An independent method reproduces the kernel

Deep inside, the transmitted wave as a function of `v = t + r*` is exactly
`K ⋆ g`. PR #274's evolution is a characteristic null-grid march rather than a
frequency-domain transfer matrix, so this is real cross-validation of the
*propagation* — and it exposed a subtlety about where the pulse may be launched.

It is **not** independent in every sense: both rest on the same `potential` and
tortoise definitions, and this module imports the characteristic solver. The
accurate claim is *independent numerical propagation methods on the same
operator*, which tests the propagation and not the operator.

The incident amplitude is only defined where the wave is free, and the phase it
has *not yet* accumulated is set by the potential remaining beyond the launch
point — which is known exactly, `∫_{r*}^∞ V dr* ≈ L/r*`:

| launch `r*` | `∫V` beyond launch | max diff | rms diff |
|--|--|--|--|
| `100` | `0.0375` | `2.73%` | `0.50%` |
| `200` | `0.0187` | `1.40%` | `0.25%` |
| `400` | `0.0094` | **`0.73%`** | **`0.13%`** |

Successive ratios `1.94`, `1.93` — the residual **halves as the launch radius
doubles**, exactly as `L/r*_launch` predicts. That is what identifies it as
placement rather than a disagreement between methods.

PR #274 launched at `r* = 6`, where `V ≈ 0.1`. That is harmless for extracting a
quasinormal frequency, because a ringdown does not care how it was excited, and
fatal for defining a transmission ratio.

> **An earlier number here was flattered by cancellation.** With plane-wave outer
> matching this check read `0.92%` at `r* = 100`. That was two errors partly
> cancelling — the plane-wave outer condition carried its own error in the
> opposite direction. Under the correct Jost condition the same launch reads
> `2.73%`, and the series converges cleanly as `1/r*_launch`. **The larger number
> is the honest one**, and the convergence is worth more than the smaller
> number was.

---

## What the causality gate caught

Two artefacts, both sitting at exactly the amplitude of the physics being looked
for, and neither visible from the `t > 0` side alone.

**A missing DC cell.** Sampling `ω` at right endpoints left `[0, dω]`
uncovered. Since `T(0) = 0`, `S(0) = −1`, and the omission put a constant
`≈ −dω/π ≈ −1.9e-3` under the *entire* kernel — including `t < 0`, where the
answer is exactly zero. Midpoint sampling fixed it. **Would have been read as:**
a constant offset, i.e. a spurious permanent memory.

**Gibbs ringing from the `1/ω` tail.** Before the analytic subtraction the
kernel sat on a `1.6e-3` plateau at large `|t|` — on *both* sides of the origin.
**Would have been read as:** a late-time tail, which is precisely the quantity
this round would most like to measure.

> **The general point, correctly scoped.** A quantity that is exactly zero on
> part of its domain is a free monitor for the artefacts that violate it — on
> the same run, with no reference value.

**It is not a total error bar.** A numerical error can be perfectly causal and
live only at `t > 0`: outer-boundary error, finite `ω_max`, the arbitrary
subtraction parameter, quadrature bias. This round demonstrated exactly that —
replacing plane-wave outer matching with the Jost condition moved positive-time
results with no matching negative-time signature. So the negative-time residual
bounds one *family* of error and no other.

Accordingly the late-time tail is judged out of reach on **positive-time**
evidence, not from the causality floor:

| setting | `t=40` | `t=60` | `t=100` | `t=150` | `t=200` |
|--|--|--|--|--|--|
| base | `+7.9e-07` | `+7.7e-08` | `+3.5e-07` | `+1.9e-07` | `−3.3e-07` |
| `ω_max=20, N=2048` | `−5.7e-06` | `+6.8e-07` | `−2.3e-06` | `−5.8e-07` | `+7.9e-07` |
| `N=8192` | `+7.9e-07` | `+7.4e-08` | `+3.4e-07` | `+1.8e-07` | `−2.9e-07` |
| `a=0.5` | `+5.4e-06` | `+6.1e-07` | `+2.0e-06` | `+8.4e-07` | `−1.7e-06` |
| `a=2.0` | `−8.7e-06` | `−7.0e-07` | `−3.1e-06` | `−9.8e-07` | `+2.3e-06` |

The spread **exceeds** the values and the sign is not even stable. No exponent
is quoted.

This is the complement to PR #274's lesson. That round needed an external
published spectrum to discover its error floor. Here the structure of the
problem supplied a monitor for one family of error — but only one.

---

## Ledger

| claim | verdict | evidence |
|--|--|--|
| the retarded transfer kernel exists as a computable object | **YES, DELIVERED** | `K = δ(t) + K_reg` from `T(ω)` on a well-conditioned real-frequency problem |
| the frequency domain can be used despite #270's and #274's failures | **YES — DIFFERENT PROBLEM** | those failed for `Im ω < 0`; for real `ω` nothing dominates, unitarity `~1e-13` |
| the kernel is causal | **YES, TO `~1e-6`** | measured directly; bounds **acausal** artefacts, not total error |
| the kernel carries the published ringdown | **YES** | against the external value: real `0.062%`–`1.17%`, damping `0.11%`–`3.80%` |
| an independent propagation method reproduces the kernel | **YES, `0.73%` PEAK / `0.13%` RMS** at `r* = 400` | convolution vs #274's characteristic null-grid march — a different algorithm on the same operator |
| the transfer is rigid / instantaneous | **NO — AND ANALYTICALLY SO** | `T(∞) = 1` while `T(0) = 0` forces `∫K_reg dt = −1`; measured `−0.999996`, mass `2.03` |
| the late-time power-law tail is measured | **NO** | argued at positive times: spread across settings exceeds the values, sign unstable |
| #274's pulse placement was adequate for this object | **NO — AND IT WAS FOR ITS OWN** | `43%` mismatch from launching at `r* = 6`; harmless for the ringdown |

---

## What this round adds to the standing list

**An exactly-zero region is a free monitor for the artefacts that violate it.**
PR #274 needed a published spectrum to find its floor. Here causality supplied a
monitor on the same run, with no external reference, and it caught two separate
*acausal* artefacts at exactly the amplitude of the physics being sought. When a
problem has a region where the answer is known in advance, compute there first.
But scope what it covers: it bounds the errors that violate the known value, not
the ones that respect it.

**Check what an estimator measures before attributing what it returns.**
`Re[iω(T−1)]` carries a deterministic `−c³/(6ω²)` bias, and an earlier draft read
the resulting shortfall as a physical property of the potential's sub-leading
tail. The tell was already in the data: the deficit was independent of the outer
edge, which a truncation effect would not be and a fixed-form bias would.

**A sum rule beats a fit.** `∫K_reg dt = −1` is forced by `T(0) = 0`, not
obtained by integration. The measured `−0.9978` is then a check on the numerics
rather than a result, which is the right way round.

**The same data can be adequate for one question and not another.** #274's pulse
placement gave a correct quasinormal frequency and a 43%-wrong transmission
ratio, from the identical run. What made it adequate there — that a ringdown is
insensitive to its excitation — is exactly what makes it inadequate here.

---

## Still open

**The late-time tail**, which needs dynamic range where this method has none. The
frequency-domain route is excellent near the barrier and poor in the tail — the
opposite of a long time-domain evolution. That is a separate calculation with a
separate method, not a refinement of this one, and no exponent should be quoted
until it is done.

**Whether any BAM exchange kernel is meant to approximate this object.** That is
a modelling question, not a numerical one, and it is the natural next step now
that the geometric side is in hand.

**Still blocked, unchanged and unrelated:** the nonlinear backreaction oracle, by
the `C(v)` / inner-flux issue.
