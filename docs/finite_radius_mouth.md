# A finite-radius mouth: does the negative mode survive?

**PR #261.** The gate PR #260 set, answered.

> **The answer is no** — and the statement is *structural*, not parametric. At
> `λ < 0` the tube's channel functions are strictly negative and the resolved
> mouth's are strictly positive, so `det(A − 𝒢) = 0` has no solution there for
> **any** choice of `a, d, L, 𝒜, m`. What produced #260's mode was the
> **linearization** of the mouth's self-energy: freezing a *screened* quantity at
> its unscreened leading term. The mode it manufactured sits at `κa ≈ 1`, exactly
> the edge of that approximation's validity.

> **Scope.** A linear conformally coupled scalar on a **fixed** Einstein static
> universe. The mouths are spheres of radius `a` coupled to the tube through
> **one channel each** — monopole matching, so the `ℓ ≥ 1` content on each sphere
> is dropped (§8 quantifies it). They are spheres in a fixed ambient, **not a
> solved neck geometry**. **No backreaction.**

---

## 1. The question

#260 built the conservative finite throat: a tube with an exact
Dirichlet-to-Neumann map, a reflected channel and a real traversal delay. It also
found that with **point** mouths the composite carries an exponentially growing
mode for every choice of parameters, generated at the point-mouth/tube
**interface** and localizing — in the `σL, σd ≫ 1` limit — to a single mouth, at
a rate `σ* = 2√(π/𝒜)` containing neither the tube's length nor the mouths'
separation. One question followed:

> does the negative mode survive a finite-radius mouth?

Both branches were spelled out in advance. If it survives, the point-interaction
throat family of #255–#260 is the wrong model of a wormhole mouth. If it does
not, the mode was an artifact of the matching and the delay carries over.

## 2. What a resolved mouth changes

A point interaction in three dimensions has to subtract the `1/(4πχ)` divergence,
and what it keeps is the **renormalized** self-energy

```
g(λ) = −(ω/4π) cot(πω)  ,     g(0) = −1/4π²  <  0
```

negative at threshold and falling like `−κ/4π` down the imaginary axis. That
negativity is a renormalization artifact — it is the finite part left after an
infinite subtraction — and it is what #260's mode fed on.

A sphere needs no subtraction. Smear the coupling over `∂B_a(c_i)` with a uniform
weight — **the same operator on both sides**, so the composite stays manifestly
self-adjoint — and the ambient block becomes

```
𝒢_self(λ)  = f(a,λ)·G(a,λ)          𝒢_cross(λ) = f(a,λ)²·G(d,λ)
```

with `G` the **unsubtracted** Green function and `f(χ,λ) = sin(ωχ)/(ω sin χ)` the
regular radial solution, `f(0) = 1`. `f` appears once per smearing: it is the
form factor of a uniform shell on emission, and the mean-value factor on
reception.

### These are identities, and they are checked

| | quadrature | prediction | relative |
| --- | ---: | ---: | ---: |
| `⟨G(·,c₂)⟩` over `∂B_a(c₁)` | — | `f(a)G(d)` | **`1.0e-10`** |
| `⟨⟨G⟩⟩` over one sphere, twice | — | `f(a)G(a)` | `4.1e-04` |

The cross identity is the mean-value theorem and lands at quadrature precision.
The self one carries the Green function's integrable singularity at coincidence,
so its residual is **grid-limited** and is reported as quadrature error rather
than as a model error.

## 3. The answer: the signs decide it

At `λ = −σ²` the two sides of the secular equation have opposite signs, and both
are manifest rather than scanned:

* **the tube** gives `−coth(κL/2)/(𝒜κ)` and `−tanh(κL/2)/(𝒜κ)` — strictly
  **negative**. A passive interior supplies no restoring force;
* **the ambient** gives `f·G(a) ± f²·G(d)` — strictly **positive**. Written so
  the positivity is visible by inspection,

```
f·G(a)  = (1 − e^{−2κa})(1 − e^{−2κ(π−a)}) / (8πκ sin²a (1 − e^{−2πκ}))
f²·G(d) = e^{2κa−κd}(1 − e^{−2κa})²(1 − e^{−2κ(π−d)}) / (16πκ² sin²a sin d (1 − e^{−2πκ}))
```

every bracket lies in `(0,1]`, and `2κa − κd < 0` because disjoint mouths need
`a < d/2`. So `𝒢_sym > 0` and `𝒢_anti > 0` for every `κ`.

**A difference of a negative and a positive number has no zero.** No parameter
choice can produce a growing mode, and a sweep confirms what the signs already
settle:

| | |
| --- | ---: |
| samples over `(a, d, L, 𝒜, m, σ)` | **`3078`** |
| sign changes found | **`0`** |
| worst approach to zero | `−5.1e-04` |
| tube side negative at every sample | **yes** |
| ambient side positive at every sample | **yes** |

(The exponentials underflow to exactly zero once `κχ ≳ 745`. That is harmless:
where the ambient term underflows, the difference is the tube's term alone, which
is strictly negative. A test pins that case rather than leaving it to the sweep.)

## 4. And #260's mode was the linearization

That round modelled the mouth by adding a **constant** `1/(4πa)` to the
self-energy — the leading term of `G(a,λ) = 1/(4πa) + g(λ) + O(a)`. The exact
`G(a,−κ²)` is **screened**, `≈ e^{−κa}/(4πa)`, and dies down the imaginary axis.
The constant does not. So the linearized ambient eventually loses to the tube's
`−1/(𝒜κ)` and the difference changes sign, while the exact one never can.

| `a` | linearized root `κ*` | `κ*·a` | exact root |
| ---: | ---: | ---: | :---: |
| `0.02` | `50.02` | **`1.0004`** | none |
| `0.05` | `20.05` | **`1.0025`** | none |
| `0.15` | `6.814` | **`1.0221`** | none |
| `0.35` | `3.220` | **`1.1269`** | none |

**The root sits at `κa ≈ 1` — precisely where `κa ≈ 1` and the linearization is
invalid.** The two models agree to `0.8%` for `κa ≤ 0.1` and differ by `1000%` at
`κa = 3`; they disagree not in magnitude but in **sign**.

A mode that lives at the scale where its own derivation fails is an artifact.
#260 measured that `σ*` loses `L` and `d` asymptotically and located the mode at
the mouth/tube interface, which pointed here — but it could only record the
suspicion. This is the demonstration, and it names the culprit precisely: not the
interface as such, but the *constant* standing in for a screened function
across it.

## 5. Where the mode went: soft, and positive

It did not vanish. The composite has exactly one state below the free ESU gap
`λ = 1`, in the symmetric channel — the tube's zero mode, which a point cannot
hold and two spheres can — and as `a → 0`:

| `a` | `λ₀` | `8πa/(𝒜L)` | ratio |
| ---: | ---: | ---: | ---: |
| `0.200` | `0.35158841` | `0.44444444` | `0.791` |
| `0.050` | `0.10759125` | `0.11111111` | `0.968` |
| `0.020` | `0.04398357` | `0.04444444` | `0.990` |
| `0.005` | `0.01108557` | `0.01111111` | **`0.998`** |

```
λ₀  ⟶  8πa/(𝒜L)
```

two mouth capacitances `4πa` restoring a tube of volume `𝒜L`. The antisymmetric
channel has no bound state at all.

So **the point limit drives this mode to zero from above.** #260 did not get a
rate slightly wrong: it took a mode approaching `0⁺` like `a` and put it on the
other side of zero, at `λ ≈ −1/a²`. That is the whole failure — a sign error
produced by freezing a screened quantity, and its size grows as the mouth
shrinks.

## 6. The delay survives

| `L` | onset of `r₁₂` | `min(L, d)` |
| ---: | ---: | ---: |
| `0.4` | `0.2083` | `0.4` |
| `0.6` | `0.4097` | `0.6` |
| `0.9` | `0.7095` | `0.9` |
| `1.2` | `1.0094` | `1.2` |
| `2.0` | `1.0666` | `1.3` |
| `3.0` | `1.0666` | `1.3` |

`d(onset)/dL = **1.0010**` against a predicted `1`, and the onset stops moving
above the ambient path to **`0.0`** exactly. The mouth adds only a **sub-leading**
`O(a)` shift, measured slope `−0.39`.

A first draft of this section predicted `−2a` — one radius per leg, the wave
reaching the mouth's surface rather than its centre — from a version of the
ambient block with the shell form factor `f(a)` left out. Restoring `f` changes
the answer, because `f` carries an advance of its own that partly cancels the
surface offset. The measured slope is quoted and **the prediction is recorded as
wrong** rather than quietly dropped.

The contour is also easier now. #260 needed `ε > σ* ≈ 2` to clear the growing
mode; with no mode to clear, `ε` is back to #259's single requirement of sitting
above the frequency spacing, and `0.4` is comfortable.

## 7. And so do the static results

| `λ` | `det S` | `det S/λ` | `𝒲` | `−β(λ)` |
| ---: | ---: | ---: | ---: | ---: |
| `1e-08` | `−3.6e-08` | `−3.6379` | `−8841941.294812` | `−8841941.294820` |
| `1e-06` | `−3.6e-06` | `−3.6379` | `−88419.424766` | `−88419.424765` |
| `1e-04` | `−3.6e-04` | `−3.6402` | `−884.206065` | `−884.206065` |

The static response still collapses onto the antisymmetric direction, `det S`
still vanishes linearly in `λ`, and #258's defect is still `𝒲 = −β(λ)` — to
`3.6e-12`. All three came from the **tube's** zero mode, which the mouth does not
touch. Only the coefficient moves (`149.08 → −3.638`, and the sign flips),
because the ambient diagonal is now `+f(a)G(a)` instead of the renormalized
`g₀ < 0`.

## 8. What is still put in

One channel per mouth means only the `ℓ = 0` projection on each sphere is
coupled. The dropped content obeys **#250's screening law**:

| `a` | dipole/monopole | `÷ (a/d)` | dropped power |
| ---: | ---: | ---: | ---: |
| `0.02` | `0.01437` | `0.9343` | `6.9e-05` |
| `0.05` | `0.03594` | `0.9344` | `4.3e-04` |
| `0.15` | `0.10799` | `0.9359` | `3.9e-03` |
| `0.35` | `0.25394` | `0.9432` | `2.2e-02` |

The dipole/monopole ratio is `0.934·(a/d)` across a decade in `a`, so the leading
omission is the dipole at `O(a/d)` and the dropped power fraction is `6.9e-05` at
the working radius.

**And the mouths are spheres in a fixed ambient, not a solved neck.** Nothing
here derives the geometry that would produce them, and there is **no
backreaction**.

## 9. What this ungates

#260 gated stationary action and backreaction on this question, for a concrete
reason: an integral over a solved field on a background with a growing mode
measures the mode. That reason is now removed. The finite-mouth throat is
conservative, has a real traversal delay, and has **no** growing mode — so an
on-shell action, or an `A`/`B`/`A+B` collapse comparison, computed on it measures
the physics.

```
ray closure → field solution → two-wave invariant → finite throat
            → mouth resolved ✓ → neck resolved ✓ → stationary action
            → backreaction → topological branch
```

What is *not* settled here: the neck geometry, the multipole content above
`ℓ = 0`, and backreaction. Those are limitations of this round, not gates on the
next one — and each has a number attached rather than a caveat.

**Two of them are closed by #262** (`docs/finite_radius_neck.md`), which removes
the balls from the ambient and re-asks the question. The answer does not move,
and it becomes a *theorem* rather than a sign on a `2×2`: with nothing
subtracted the energy is a sum of non-negative terms, so `λ > 0` for every
configuration, all multipoles, no sweep. That round also shows the `ℓ ≥ 1`
sectors decouple from a one-channel tube — so the monopole truncation flagged in
§8 was never a *stability* limitation — and it prices the fixed ambient at
`3.8e-03` at the working radius. It corrects one claim above: **§5's "exactly
one state below the gap" holds for `L < π`, not structurally.**
