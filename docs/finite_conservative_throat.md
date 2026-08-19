# A finite conservative throat: a DtN map, a point limit, and a traversal time

**PR #260.** Roadmap: the flux-conserving throat operator that PRs #253–#259 all
listed as still owed.

> **The closure result, stated first.** The tube is conservative and has a real
> traversal delay — but **point-mouth matching makes it dynamically unstable**,
> for every choice of parameters, with a growing mode that belongs to the mouth
> and not to the interior (§8). Nothing here cures it; moving the retarded
> contour above the mode evaluates the correct retarded solution *of an unstable
> system*. Whether a finite-radius mouth or neck geometry removes the mode is
> the open question, and it should be settled before PR #261's stationary action
> or PR #262's backreaction.

> **Scope.** A linear conformally coupled scalar on a **fixed** Einstein static
> universe. The throat now has an interior — a tube of proper length `L`,
> cross-section `𝒜` and interior mass `m` — but that interior is
> **one-dimensional** (one transverse channel), and the mouths are still
> **points** in the ambient. §8 measures exactly what the second approximation
> costs. **No backreaction**, no topology change, and `L`, `𝒜`, `m` are chosen
> rather than derived.

---

## 1. The debt

Every round from #253 to #259 carried the same disclaimer, in the same words:
*the throat is point-supported — no interior, no proper length, no delay*. And
the capstone's ledger was more specific about what was wrong with the object
that stood in for one: a rank-one **mouth-transfer** model — field values only,
no normal-derivative matching, no reflected channel, `1×1` where a conserving
junction needs `2×2`, and **lossy** for `κ < 1`.

This round replaces it.

## 2. What is put in — two lines

A tube of length `L` and cross-section `𝒜` joins the mouths. Its interior
problem is `u'' + k²u = 0` with `k² = λ − m²`, and its **Dirichlet-to-Neumann
map** — the outward flux at each end given the end values — is exact:

```
N(λ)  =  𝒜k · [[ cot kL , −csc kL ] ,
               [ −csc kL ,  cot kL ]]
```

The matching to the ambient is value continuity and flux continuity at the
mouths:

```
q  =  −N(λ) Φ ,          Φ  =  φ^reg + q/(4πa)
```

with `a` the mouth radius (optional; `a → ∞` is the pure point-mouth matching
used throughout unless stated). Everything below follows from those two lines.

Both entries of `N` are **even in `k`**, so `k = √(λ − m²)` needs no branch
choice: `N` is meromorphic in `λ` with no cut, and the square root is cosmetic.
`det N = −(𝒜k)²` exactly, which is why the chart has a closed form:

```
A(λ)  =  −N(λ)⁻¹ − I/(4πa)  =  (1/𝒜k)[[ cot kL , csc kL ] ,
                                       [ csc kL , cot kL ]] − I/(4πa)
```

So the throat's **transmission amplitude** is `β(λ) = csc(kL)/(𝒜k)` and its
**self-energy** is `α(λ) = cot(kL)/(𝒜k) − 1/(4πa)`: one function of the
interior, read at two arguments.

**The boundary condition is now frequency-dependent, and that dependence *is*
the interior.** A point throat is a fixed Hermitian `A`; a finite throat is
`A(λ)`. Everything in this round is a consequence of that one sentence.

## 3. Where the self-adjointness lives

The conservative object is the **enlarged system** — ambient `⊕` tube — with the
*fixed*, `λ`-independent matching of §2. That is one self-adjoint operator on
`L²(S³) ⊕ L²([0,L])`, and it is what conserves flux.

Eliminating the tube leaves an ambient-only problem with a **`λ`-dependent**
boundary condition: `A(λ) = −N(λ)⁻¹` is the **Weyl (`M`-) function** of that
elimination. `A(λ)` is *not* itself a self-adjoint operator on the ambient space
— an energy-dependent boundary condition never is — and nothing here claims it
is. What it is, is a matrix **Nevanlinna** function, monotone in `λ` between its
poles, and that monotonicity is the enlarged system's self-adjointness showing
through after the elimination. It is measured (§8), not assumed.

**A sign convention, since the word cuts both ways.** `A(λ)` here *decreases*
with `λ`, so in the standard Herglotz/Nevanlinna convention — a map of the upper
half-plane to itself — it is `−A` that is Nevanlinna and `A` that is
anti-Nevanlinna. (Read the small-`L` symmetric channel, `A_sym(z) ∼ 2/(𝒜Lz)`,
which sends the upper half-plane to the lower one.) The boundary-triple
convention used here puts the sign on `A`; nothing depends on the choice, but the
monotonicity direction does, so it is stated rather than left to the reader.

What *can* be checked pointwise is that the elimination is faithful: at each
`λ`, the eliminated problem is a maximal self-adjoint boundary condition for the
ambient problem **at that `λ`** — `rank[B|C] = 2` with `BC†` Hermitian. Here
`B = N` and `C = −(I + N/(4πa))`, so `BC† = −N − N²/(4πa)`, Hermitian because
`N` is real symmetric.

| | |
| --- | ---: |
| `‖BC† − CB†‖`, seven `λ` on both sides of zero | **`0.0`** |
| `rank[B\|C]` at every one of them | `2` |
| imaginary part of `[B\|C]` | `0.0` |
| DtN vs its own interior, by Green's identity under quadrature | `1.5e-07` |
| **control** — #255's rank-one transfer model | **`0.30`** |

The control's defect is the size of the coupling itself. That is what "lossy for
`κ < 1`" was, stated as an operator property rather than as bookkeeping.

The Green's-identity check matters because it tests the matrix against the
*interior it claims to summarize* rather than against itself: a wrong sign or a
`cot`/`csc` swap fails it immediately. It is the **sesquilinear** form — `Φ†NΦ`
and `|u'|²` — because that is the one that expresses energy; the bilinear
version agrees for real `k` and real end values, which is why a first draft
passed while checking something weaker. The interior derivative is analytic —
#248's lesson about `np.gradient`'s one-sided end differences applies with
particular force to a boundary-value identity.

## 4. The result: the throat transmits at the traversal time

The measured object is the **two-mouth block's** impulse response — the inverse
transform of `R(ω) = (A(ω) − Γ(ω))⁻¹` along the retarded contour. The source and
observer legs are gone. But `Γ` is the *ambient's* own mouth-to-mouth
propagator and stays in, so this is the **coupled** ambient+tube response and
not the throat in isolation — which is exactly why the second prediction below
reads `min(L, d)` rather than `L`.

Two predictions, and they are different:

| | prediction | measured |
| --- | ---: | ---: |
| `r₁₁` — same mouth in and out | `0` | **`0.0000`** |
| `r₁₂` — opposite mouths | `min(L, d)` | slope `1.007` |

* **A wave that reaches a mouth is partly reflected instantaneously.** `r₁₁`
  starts at `t = 0` and its tube echoes arrive later, at `2L, 4L, …`. The
  point throat has the reflection and none of the echoes.
* **`r₁₂` starts at `min(L, d)`.** The tube's own path takes the traversal time
  `L`. But the *ambient* also connects the two mouths, along a geodesic of
  length `d`, and that path is there **whether or not the mouths are joined**.

| `L` | `σ*` | contour `ε` | onset of `r₁₂` | `min(L, d)` |
| ---: | ---: | ---: | ---: | ---: |
| `0.4` | `1.769` | `2.569` | `0.2586` | `0.4` |
| `0.6` | `1.575` | `2.375` | `0.4601` | `0.6` |
| `0.9` | `1.417` | `2.217` | `0.7622` | `0.9` |
| `1.2` | `1.325` | `2.125` | `1.0643` | `1.2` |
| `2.0` | `1.206` | `2.006` | `1.1375` | `1.3` |
| `3.0` | `1.152` | `1.952` | `1.1375` | `1.3` |

The probe pulse's tail puts every onset early by the same fixed amount, so what
is quoted is the **slope**: `d(onset)/dL = 1.0071` below the ambient path, and
`0` above it — the last two rows agree to `0.0`, exactly.

**That ambient path is #258's cross-mouth channel and #259's `β = 0` control,
now separated in time rather than by rank counting.** Those two rounds had to
argue that the off-diagonal block is not "through the throat"; here the two
contributions arrive at *different times*, and which one is first is decided by
`min(L, d)`.

The point-throat control, `A` frozen at `A(λ₀)`: `r₁₂` starts at **`0.0000`**.
A point throat transmits instantaneously, which is what a point throat is.

### The massless tube's ledger is a derivation

On the retarded contour `Im x > 0`, so both geometric series converge:

```
cot x  =  −i − 2i Σ_{k≥1} e^{2ikx}          csc x  =  −2i Σ_{k≥0} e^{i(2k+1)x}
```

with `x = kL`, checked against the closed forms **on the actual contour** to
`4.5e-16` and `1.7e-15`. Reading the delays off:

| entry | delays |
| --- | --- |
| same mouth (`cot`) | `0, 2L, 4L, …` |
| opposite mouths (`csc`) | `L, 3L, 5L, …` |

**The parities are the physics** — an even number of traversals returns to the
mouth it came from — and the reflected channel is the one the rank-one transfer
model does not have at all.

Two scopes belong on this. It is the **`m = 0`** ledger: there `k = ω`, so
`e^{ikL}` is a pure translation and the interior really does return shifted
copies. With `m ≠ 0`, `k = √(ω²−m²)` and `e^{ikL}` is **dispersive** — the causal
*front* is still at `L`, but the echoes are no longer translates of one shape,
which is §7's evanescent physics seen from the other side. And it is the ledger
of the **tube kernel `A(ω)`**, not of the coupled `R = (A − Γ)⁻¹`, which also
carries the ambient's `d`-paths — which is precisely why §4 reads `min(L, d)`.

## 5. There *is* a point limit — and it is not a finite `A`

Freezing `A` at `A(λ₀)` reproduces the finite throat exactly at `λ₀` and nowhere
else. The band it is right on has width `Δω ∼ 1/L` in **frequency** — in `λ`
that is `Δλ ∼ 2√λ/L`, and the two are not the same:

| `λ/λ₀` | `β` exact | `β` frozen | relative error |
| ---: | ---: | ---: | ---: |
| `1.00` | `0.101589` | `0.101589` | **`0.0`** |
| `1.05` | `0.097446` | `0.101589` | `4.3%` |
| `1.20` | `0.087127` | `0.101589` | `16.6%` |
| `2.00` | `0.058864` | `0.101589` | `72.6%` |
| `3.00` | `0.045947` | `0.101589` | `121%` |

And the two channels fail **differently**, which is the part that matters:

* the **antisymmetric** channel does have a limit. `A_anti = −tan(kL/2)/(𝒜k) →
  −L/(2𝒜)`, with the error falling like `L²` — measured `1.3e-3 · L²` down to
  `1.7e-4 · L²` over a decade in `L`;
* the **symmetric** channel does not. `A_sym = cot(kL/2)/(𝒜k)` diverges like
  `2/(𝒜λL)`, matched to `1.4%` at `L = 0.4` and better below, because **a
  massless tube holds a zero mode and a point cannot**.

A first draft concluded from those two facts that **the limit does not exist**.
That was wrong, and the correction is the more interesting statement.

A boundary pair is defined only up to `(B, C) → (MB, MC)` for invertible `M`, so
a chart matrix running off to infinity means the limit has **left the chart**,
not that it is absent. Row-scaled — the symmetric row by `1`, the antisymmetric
row by `1/N_anti` — the pair converges cleanly:

| `L` | `‖B − P_anti‖`, `‖C + P_sym‖` | rate `/L` | `rank[B\|C]` |
| ---: | ---: | ---: | ---: |
| `0.40` | `2.5473` | `6.3683` | `2` |
| `0.20` | `1.2608` | `6.3042` | `2` |
| `0.10` | `0.6288` | `6.2884` | `2` |
| `0.05` | `0.3142` | `6.2845` | `2` |
| `0.02` | `0.1257` | `6.2834` | `2` |

```
(B, C)  ⟶  (P_anti, −P_sym)  ,   i.e.   Φ_anti = 0   and   q_sym = 0
```

linearly in `L`, at the rate `𝒜λ/2 = 2π` for the working parameters. That is a
**mixed Dirichlet–Neumann** point extension: maximal (`rank[B|C] = 2` the whole
way), self-adjoint (`BC† = 0`), and reachable by **no finite Hermitian `A`**,
since both blocks are singular there. It is exactly the kind of stratum #257's
review insisted the finite-`A` chart does not reach.

Which is also physically what a very short pipe should do: short the two mouths
together (`Φ_anti = 0`) and store nothing (`q_sym = 0`).

So the honest statement is **"no finite-`A` point limit"**, and the constant-`A`
family of #257–#259 is this throat *read at one frequency* — the traversal time
being exactly what a frequency-independent boundary condition cannot carry.

## 6. The static limit is rank one, and #258's tomography breaks on it

The same zero mode empties the symmetric channel at `λ = 0`. The static response
`S = Re R` collapses onto the antisymmetric direction:

| `λ` | `det S` | `det S/λ` | `\|S₁₁ + S₁₂\|/\|S₁₁\|` | `𝒲` |
| ---: | ---: | ---: | ---: | ---: |
| `1e-08` | `1.49e-06` | `149.076` | `0.0` | `−8.84e+06` |
| `1e-06` | `1.49e-04` | `149.076` | `4e-07` | `−8.84e+04` |
| `1e-04` | `1.49e-02` | `149.092` | `4.3e-05` | `−884.2` |
| `1e-02` | `1.507` | `150.715` | `4.3e-03` | `−8.853` |

`det S → 0` **linearly in `λ`**, with the coefficient constant to `1e-3` over
four decades (the last row is kept because it shows the correction turning on,
`1.1%`). So #258's disconnection defect `𝒲 = S₁₂/det S − G₀` **diverges like
`1/λ`**.

**What that falsifies is the generic finite-`A` family** — every member of which
has `rank S = 2` — and **not point-ness.** The tube's own short-tube limit, the
mixed stratum of §5, gives `R → diag(0, −1/Γ_anti)` in the channel basis: rank
one as well, and the tube converges to it. A first draft of this round claimed
the stronger version, that rank one distinguishes a finite throat from a point
one. It does not; it distinguishes both from the chart.

Give the tube an interior mass and the rank comes back:

| `m` | `det S` | `det S/m²` | `𝒲` | `−β(0)` | error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `0.05` | `−0.3724` | `−148.98` | `35.3558` | `35.3558` | **`0.0`** |
| `0.10` | `−1.4869` | `−148.69` | `8.83002` | `8.83002` | **`0.0`** |
| `0.20` | `−5.9012` | `−147.53` | `2.19859` | `2.19859` | **`0.0`** |
| `0.40` | `−22.890` | `−143.06` | `0.540863` | `0.540863` | **`0.0`** |

Same coefficient — `det S ∝ (λ − m²)`, one statement read two ways — and then
the closure: **off the collapse, `𝒲 = −β(λ)` exactly**, to `3.1e-13`, with `β`
the tube's own transmission amplitude. #258's theorem survives the
generalization to a frequency-dependent conservative throat; what it returns is
no longer a constant but the interior's amplitude at that frequency.

## 7. An interior mass gives the channel a mass gap

Below `λ = m²` the wavenumber turns imaginary and

```
β(λ)  =  csc(kL)/(𝒜k)  →  −csch(κL)/(𝒜κ)  ≈  −2e^{−κL}/(𝒜κ)
```

— negative, and exponentially suppressed rather than oscillating. Matched to its
asymptote to `7.6e-05` deep in the evanescent region; suppression `3.1e-03`
across the cutoff.

The discriminator is **monotone decay against oscillation**, not the sign: above
the cutoff `β` changes sign every time `kL` passes a multiple of `π`, so a
single high-`λ` sample proves nothing. Measured across a sweep: **`0` sign
changes below the cutoff, `7` above**, and `|β|` monotone below.

So a throat with a massive interior has a **mass gap**, and below it the channel
is *evanescent*: the mouths look like two ordinary scatterers with a tunnelling
correction, and the same throat above the gap is fully connected. (An earlier
draft called this "low-pass", which is backwards — it is the *low* frequencies
that are suppressed. It is an evanescent cutoff, not a filter passband.) And the
cutoff is where the rank collapses — `k → 0` makes `β → 1/(𝒜k²L)` diverge — so §6 and §7 are one
statement read at `λ = m²`.

Both regimes stay real and self-adjoint: worst imaginary part `0.0`.

## 8. The model fails the stability gate — at the mouth/tube interface

`A(λ)` is decreasing in `λ` and `Γ(λ)` increasing (#257's Gram identity), so
`A − Γ` is strictly monotone between poles and each channel has **at most one
root per pole-interval** — a count, not a scan. In the symmetric channel there
is always exactly one, and it is at `λ < 0`: the tube's zero mode, pushed below
zero by coupling to a point.

Three facts identify what it is. All three are **limits**, and the convergence
is what is measured:

| `𝒜` | `L` | `σ*` | `2√(π/𝒜)` | rel. error | channel split | `σ*L` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `0.2` | `1.5` | `7.926726` | `7.926655` | `9.0e-06` | `1.4e-04` | `11.9` |
| `0.2` | `3.0` | `7.926672` | `7.926655` | `2.2e-06` | `3.5e-05` | `23.8` |
| `0.2` | `6.0` | `7.926672` | `7.926655` | `2.2e-06` | `3.5e-05` | `47.6` |
| `0.5` | `3.0` | `5.014024` | `5.013257` | `1.5e-04` | `1.5e-03` | `15.0` |
| `1.0` | `3.0` | `3.550133` | `3.544908` | `1.5e-03` | `1.1e-02` | `10.7` |

* its rate tends to **`σ* = 2√(π/𝒜)`** — closed form, and **`L` is not in it**;
* its dependence on the mouth separation dies exponentially: at `σ*d ≈ 19` and
  `σ*d ≈ 24`, two throats agree to **`3.9e-09`**;
* the two channels split by **`1.04 · e^{−σ*d}`** — ratio `1.0378` to `1.0448`
  across the five rows with `σ*L > 14` — which is the Euclidean mouth-to-mouth
  propagator itself, not a bound on it.

**So the mode is generated at the point-mouth/tube interface, and in the
`σL, σd ≫ 1` limit it localizes to a single mouth** and stops knowing the tube's
length or the other mouth. Its length scale there is `1/σ* = √(𝒜/4π)` — the
mouth's own radius, which is exactly where "point mouth" stops being an
approximation.

**But the working throat is not in that limit**, and that qualification is
load-bearing. At `𝒜 = 4π` the asymptotic form gives `σ* = 1` while `L = 0.9`
gives `1.417`, and across the lengths this round uses:

| `L` | `0.4` | `0.9` | `1.8` | `3.0` |
| --- | ---: | ---: | ---: | ---: |
| `σ*` | `1.769` | `1.417` | `1.226` | `1.152` |

a spread of `0.617`, **54% of the smallest value**. A first draft wrote "the
mode belongs to the mouth and not the interior"; that is stronger than the data
support. The interface statement plus an asymptotic localization is both correct
and more informative — and it is what makes a finite-radius mouth exactly the
right discriminator. With a finite mouth radius the closed form
generalizes to `σ* = ½[1/a + √(1/a² + 16π/𝒜)]`, which is the same statement with
the mouth's other scale in it.

So the growing mode is the price of gluing a tube of finite area to a *point*,
and it lives entirely in the regime where the gluing is invalid.

**This is the round's closure result.** The candidate throat does not pass the
stability gate #257 set up: it is conservative and it grows. The safe conclusion
is that a 1-D tube is conservative and has a real traversal delay, but that
point-mouth matching makes it dynamically unstable, and the open question is
whether finite-radius mouth or neck geometry removes the mode. That question
should be settled before stationary-action or backreaction work, not after.

Both approximations fail in the same direction and for the same reason: at
`σ*L ≲ 2` the tube is shorter than the mode, `coth(σL/2)` is not yet `1`, and
the closed form is `15%` out. That row is kept in the table.

## 9. So the contour must clear it — and be resolved

The retarded contour `Im ω = ε` must lie above **every** singularity of the
response, and a finite throat puts one at `ω = iσ*`.

| clearance | contour | `÷ 2π/span` | resolved? | onset | pedestal |
| ---: | ---: | ---: | :---: | ---: | ---: |
| `−0.03` | `1.5452` | `−1.43` | no | `0.0000` | **`0.992`** |
| `+0.02` | `1.5952` | `+0.95` | **no** | `0.0000` | **`2.6e-03`** |
| `+0.30` | `1.8752` | `+14.32` | yes | `0.4646` | `1.0e-16` |
| `+0.80` | `2.3752` | `+38.20` | yes | `0.4601` | `1.0e-16` |
| `+1.50` | `3.0752` | `+71.62` | yes | `0.4555` | `1.0e-16` |

Placed just below `σ*`, the inversion returns a field with support **before its
own light cone**: a pedestal at 99% of the peak, arriving at `t = 0` for an
event that cannot begin until `t = 0.6`.

**And clearing the mode is necessary, not sufficient.** A first draft of this
section stated the rule as `ε > σ*`, and the table above contradicts it: the
`+0.02` row *is* above the mode and still has a pedestal of `2.6e-03` with the
onset at `0`. That clearance is `0.95` of the frequency spacing
`2π/span = 0.0209` — the pole is above the contour but **unresolved by the
grid**, which is #259's lesson arriving a second time. The rule has two parts:

```
ε > σ*                 analytic:  the Bromwich contour clears the pole
ε − σ* ≫ 2π/span       numerical: the grid resolves the clearance
```

Both sides have closed forms, so both are checkable before the solve. At
clearances of `14.3`, `38.2` and `71.6` spacings the pedestal is `1.0e-16` and
the recovered onset **converges** — `0.4646`, `0.4601`, `0.4555`, a spread of
`0.0092`, four time steps.

**Clearing the contour stabilizes nothing.** Above `σ*` the inversion returns the
correct retarded solution *of an unstable system*; the growing mode is still
there and dominates every time series at late times. The delay is read from the
causal **onset**, which is a statement about analytic structure and is immune to
what the solution does afterwards — but the instability is a property of the
model, and no contour touches it.

It is the same species of error as #259's under-resolved contour — a
plausible-looking number produced by a contour in the wrong place — and it is
reported the same way: both values and the rule. The difference is that `σ*` has
a closed form, so this time the contour can be placed **before** the solve
rather than diagnosed after it.

## 10. Which frequencies each result uses

The band is not uniform across the round, so it is worth saying which result
lives where.

| result | frequencies | inside `|λ| ≪ σ*²`? |
| --- | --- | :---: |
| the delay (§4), the bounce ledger (§4) | **all** — a causal onset is a UV object | no |
| the short-tube stratum (§5) | fixed `λ₀`, and `L → 0` | yes |
| rank collapse, `𝒲 = −β` (§6) | `λ → 0` | yes |
| the mass gap (§7) | around `λ = m²` | yes |
| monotonicity and the mode count (§8) | all, by construction | — |

The delay is an exact statement about **the analytic structure of this model**,
not a prediction about a resolved physical mouth: the probe pulse that resolves
the onset has width `0.03` and carries content out to `ω ∼ 30`, far above
`σ* ∼ 1.4`. A causal onset cannot be measured any other way — sharpness is a
high-frequency property — so the honest framing is that §4 is a theorem about
the model and §8 is the statement of where the model stops describing a mouth.

Relatedly `𝒜` is a **one-dimensional coupling strength**, not an area with a
radius attached: the interior has one transverse channel and no higher modes.
At the working point `√(𝒜/4π) = 1`, a "mouth radius" of order the whole unit
`S³`, which is the same observation in different words.

## 11. What this closes, and what it does not

**Closed.** The throat operator is now genuinely flux-conserving, `2×2`, with a
reflected channel, a normal-derivative matching and a proper length. The delay
is measured and equals the traversal time. The point family of #257–#259 is
located inside the new one — as a single-frequency reading, with the band
quantified — and the short-tube limit is identified as a *specific* mixed
Dirichlet–Neumann stratum rather than as a failure. #258's `𝒲 = −β` survives,
generalized.

**Falsified.** The point-mouth matching. It is unstable, always, and the
instability is the mouth's rather than the tube's. That is the result to carry
forward.

**Not closed.**

* the mouths are still **points**, and §8 says that costs a growing mode. A
  genuinely finite mouth means solving the ambient outside two small balls, not
  a point interaction with a radius parameter — and whether that removes the
  mode is **the** open question of this round;
* the interior is **one-dimensional** — one transverse channel, so `𝒜` enters as
  a coupling strength and not as a geometry. A real tube has higher transverse
  modes above `ω ∼ 1/√𝒜`, and those are inside the band this round works in;
* `L`, `𝒜`, `m` are **chosen**, not derived. #249 is still the round that would
  fix them from matter;
* **no backreaction.** The throat is a fixed background.

**The next step is not #261.** The roadmap has a gate in it now, and it is not
optional: an action or a backreaction computed on a background with a growing
mode inherits the mode, so #261 and #262 would both be measuring it rather than
the physics they are after. The next construction is **a finite-radius mouth or
neck** — the ambient solved outside two small balls rather than a point
interaction with a radius parameter — and it has one question to answer: *does
the negative mode survive?* If it does, the point-interaction throat family of
#255–#260 is the wrong model of a wormhole mouth and the arc should say so. If it
does not, the mode was an artifact of the matching, the delay and the
conservation law carry over unchanged, and the two steps below resume.

Once that is resolved — #261, the common action and stationary history, now
has an object with a conservation law, a proper length and a delay, which is
what an on-shell action is made of, and it inherits §4's warning that any
quantity integrated over the field has to state which arrivals were present,
because the ledger now has `min(L, d)` and an echo at every `2L`. #262 —
`A`/`B`/`A+B` metric backreaction — inherits #259's warning about *which*
diagnostic to feed (`ΔT`, not `T_A:T_B`), and now also a source with a finite
size and a stated resolution limit.
