# The neck: the balls actually removed

**PR #262.** The one limitation #261 named about itself, closed.

> **The answer is still no — and it is now a theorem.** With the balls removed
> from the ambient there is no subtraction anywhere, so the composite energy is a
> sum of manifestly non-negative terms. `E = 0` forces `φ ≡ 0`, matching then
> forces `u(0) = u(L) = 0`, and Poincaré finishes it. Hence `λ > 0` for **every**
> configuration — all multipoles, no truncation, no sweep.

> **Scope.** A linear conformally coupled scalar on a **fixed** round `S³` with
> two balls cut out. The tube has **one transverse channel**, so `𝒜` is a
> coupling and the neck is a quantum-graph edge rather than a solved
> cross-section. The metric is not solved for. **No backreaction.**

---

## 1. What #261 left open about itself

That round answered #260's gate — the growing mode does not survive a
finite-radius mouth — and named its own weakest point in doing so:

> the mouths are spheres in a **fixed ambient**, not a solved neck geometry.

Which is exactly right, and it is load-bearing. #261 smeared the coupling over
`∂B_a` while still using the **whole sphere's** Green function: the balls were
never removed. It also coupled through **one channel** per mouth, so only `ℓ = 0`
talked to the tube. Two things could hide in that gap — a self-energy that is
wrong because the ball is still there, and a multipole that goes soft where the
monopole does not.

This round removes the balls. The ambient is `Ω = S³ ∖ (B_a(c₁) ∪ B_a(c₂))`, the
tube is glued along the boundary spheres, and the binary question is re-asked.

## 2. The theorem

Write the composite energy on the removed-ball geometry:

```
E[φ,u]  =  ∫_Ω (|∇φ|² + φ²) dV  +  𝒜 ∫₀^L (|u'|² + m²|u|²) ds
```

Every term is `≥ 0`. There is **nothing to subtract** — that is the whole
content of removing the ball, and it is why this is a change of footing rather
than a refinement. A point interaction has to renormalize, and what it keeps is
`g(λ) = −(ω/4π)cot(πω)`, negative at threshold. #260's mode fed on that
negativity. Here it has nowhere to enter.

`E = 0` forces `φ ≡ 0` on `Ω`. The matching condition then gives `u(0) = u(L) = 0`,
and for a function vanishing at both ends Poincaré gives

```
𝒜 ∫₀^L |u'|²  ≥  (π/L)² · 𝒜 ∫₀^L |u|²
```

so `u ≡ 0` too. Hence `λ = E/‖·‖² > 0` for every configuration — **all
multipoles, no truncation, no parameter sweep**.

Contrast the two rounds precisely. #261 had to establish a **sign** on a reduced
`2×2`: negative tube against positive ambient, checked over 3078 samples. This
round has **positivity of the form itself**. The first is an argument about a
model; the second is an argument about the problem.

## 3. The object the theorem is about, checked

A two-line argument deserves a check that the object it is about was built the
way the argument assumes.

### The exterior map

Substituting `ψ = v/sin χ` turns the `S³` radial equation for `−∇² + 1` into a
Schrödinger problem with no first derivative,

```
v''  +  [λ − ℓ(ℓ+1)/sin²χ] v  =  0
```

which is also the clearest statement of the physics: at `λ < 0` the bracket is
negative everywhere, `v` cannot oscillate, and the exterior map keeps a sign.
Regularity at the antipode picks the `e^{ℓ+1}` branch (`e = π − χ`), and the
integration starts from it. The map is

```
N_ℓ(λ)  =  −4π sin²a · ψ'(a)/ψ(a)
```

| | |
| --- | ---: |
| shooting vs. the independent `ℓ=0` closed form | **`1.7e-14`** |
| smallest `N_ℓ` over all `ℓ`, `a`, `λ < 0` sampled | `0.253` |
| positive in every channel | **yes** |
| increasing in `ℓ` | **yes** |
| `N₀(0)/4πa` at `a = 0.005` | `1.0016` |

Two of those matter beyond bookkeeping. **`N_ℓ` increases with `ℓ`**, so the
monopole is the softest channel — whatever goes soft first is the one the tube
couples to, which is the one #261 kept. And **`N₀ → 4πa`**, the capacitance of a
sphere, which is what fixes the normalization as physical rather than
conventional.

The `ℓ = 0` closed form is kept as an independent construction on purpose: the
`ℓ ≥ 1` channels are known only through the integrator, so the integrator has to
be checked somewhere it can be.

### The quadratic form

Random trial configurations, evaluated with the honest measures — `4π sin²χ dχ`
outside, `𝒜 ds` inside:

| | |
| --- | ---: |
| trials | `40` |
| smallest Rayleigh quotient | `0.3591` |
| lowest computed mode | `0.10752` |
| every quotient positive | **yes** |
| every quotient above the lowest mode | **yes** |
| purely interior configuration vs. `π²/L²` | **`2e-07`** |

The last row is the theorem's degenerate case checked exactly: `φ ≡ 0` outside
gives `u(0) = u(L) = 0`, whose lowest quotient is `π²/L²`, and it lands.

### And the sweep agrees

| | |
| --- | ---: |
| samples over `(a, d, L, 𝒜, σ)` | **`1197`** |
| sign changes found | **`0`** |
| worst approach to zero | `−1.6e-03` |
| tube side negative at every sample | **yes** |
| exterior side positive at every sample | **yes** |

This reaches #261's conclusion from a **different construction** — the diagonal
is now `1/N₀`, the inverse exterior map, where that round used a smeared average
over an ambient that still contained the balls. Two independent models, one
answer.

## 4. The monopole truncation was never a stability limitation

#261 listed monopole matching as its main open approximation. It is worth being
precise about what that costs, because the answer is: not this.

A tube with **one transverse channel** presents a single number at each mouth, so
it drives the `ℓ = 0` projection of the boundary data and nothing else. The
`ℓ ≥ 1` sectors do not see the tube at all — they are the exterior's own modes,
with the exterior's own boundary condition, and their DtN is positive:

| `ℓ` | `N_ℓ/N₀` at `σ = 0.1` | at `σ = 10` |
| ---: | ---: | ---: |
| `1` | `1.970` | `1.445` |
| `2` | `2.954` | `2.054` |
| `3` | `3.939` | `2.702` |
| `5` | `5.908` | `4.022` |

Every higher channel is **stiffer**, by at least `1.45×`. They can neither be
driven nor go unstable. #250's `(a/d)^ℓ` screening law bounds how much
**response** is missed; for the question this arc gated the roadmap on, the
higher channels were never in play.

## 5. What the fixed ambient cost, in numbers

#261's `f(a)G(a)` against this round's `1/N₀(a,λ)`:

| `a` | `λ` | fixed ambient | balls removed | relative |
| ---: | ---: | ---: | ---: | ---: |
| `0.02` | `0` | `3.954070` | `3.954594` | **`1.3e-04`** |
| `0.02` | `−4` | `3.824388` | `3.826840` | `6.4e-04` |
| `0.05` | `0` | `1.567525` | `1.568812` | `8.2e-04` |
| `0.05` | `−4` | `1.443696` | `1.449166` | `3.8e-03` |
| `0.15` | `−4` | `0.401942` | `0.413552` | `2.8e-02` |
| `0.35` | `−4` | `0.127475` | `0.142798` | **`1.1e-01`** |

The error tracks `a` and `|λ|` together, which is the expected shape: it is the
fraction of the sphere wrongly left in, weighted by how much of the field's
support sits there. So #261 was right where it looked — and a limitation with a
measured size is a different object from one with a caveat.

### And the one approximation this round has *not* removed

The reduced `2×2` keeps a **single-scattering** cross term `f(a)²G(d)`. The exact
two-ball exterior would add the series in which a monopole on one boundary sphere
drives the other and is driven back. Its expansion parameter is `cross/self`,
which the model computes anyway:

| `a` | `d` | `λ` | `cross/self` | `÷ (a/d)` | neglected order |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `0.02` | `1.3` | `0` | `0.01224` | `0.796` | `1.5e-04` |
| `0.05` | `1.3` | `0` | `0.03088` | `0.803` | **`9.5e-04`** |
| `0.15` | `1.3` | `0` | `0.09515` | `0.825` | `9.1e-03` |
| `0.35` | `0.8` | `0` | `0.39429` | `0.901` | `1.6e-01` |

`0.8·(a/d)` across a decade, the same shape as #250's screening law. So the
neglected terms are `~1e-03` at the working point and never exceed `0.16`
anywhere sampled. **A correction that small cannot flip a sign**, which is the
only thing the reduced model is being asked to decide — and the theorem of §2
does not go through the reduced model at all. It is stated here rather than left
as "exact to leading order", because that phrase is where the previous two rounds
each hid something.

## 6. The soft mode is forced, not found

#261 found one state below the free ESU gap and identified it as the tube's zero
mode held up by two mouth capacitances. It found it by scanning. It does not need
scanning: the **two ends of the gap** settle it, by the same style of argument
that kills the growing mode.

* At `λ → 0⁺` the tube's symmetric channel is `cot(kL/2)/(𝒜k) ≈ 2/(𝒜λL) → +∞`
  while the exterior block stays finite. So `F_sym > 0`.
* At `λ → 1⁻` the exterior stiffness **vanishes** — the free ESU threshold is
  exactly where the regular exterior solution goes flat, `ψ ∝ sin(ω(π−χ))/sin χ
  → 1` — so the ambient diagonal `1/N₀ → +∞` and `F_sym → −∞`.

A continuous function running from `+∞` to `−∞` has a zero. The rate at the gap
edge is closed-form: writing `ω = 1 − δ`, first order gives

```
N₀(λ)  ⟶  2π (π − a + sin a cos a) · (1 − λ)
```

| `1 − λ` | `N₀/(1−λ)` | vs. `2π(π − a + sin a cos a)` |
| ---: | ---: | ---: |
| `1e-03` | `19.15186` | `2.97e-02` |
| `1e-04` | `19.67839` | `3.05e-03` |
| `1e-05` | `19.73264` | `3.06e-04` |
| `1e-06` | `19.73808` | **`3.06e-05`** |

First order exactly — the error falls by `0.1000` per decade, so what is reported
is a **rate**, not a tolerance. At `a → 0` the slope is `2π²`. The antisymmetric
channel starts at `−L/(2𝒜) < 0` instead, which is why it has no state — *while
the tube is short enough*, which brings up the next section.

## 7. A correction to #261: "exactly one state" is a statement about `L < π`

#261 reported that the composite has **exactly one** state below the gap, in the
symmetric channel. That is true for every length that round used, and it is
**not structural**.

The channel functions have **poles**. `A_sym = cot(kL/2)/(𝒜k)` blows up at
`kL/2 = jπ` and `A_anti = −tan(kL/2)/(𝒜k)` at `kL/2 = (j−½)π`, i.e. at

```
λ = (2πj/L)²        and        λ = (π(2j−1)/L)²
```

The first of these enters the gap `λ < 1` at exactly `L = π`. Above that, each
tube harmonic that falls in brings a genuine extra state just above it:

| `L` | tube harmonics in the gap | symmetric roots | antisymmetric roots | states |
| ---: | ---: | :--- | :--- | ---: |
| `0.9` | `0` | `0.107516` | — | `1` |
| `3.0` | `0` | `0.031929` | — | `1` |
| `4.0` | `1` | `0.023771` | `0.667445` | `2` |
| `6.0` | `1` | `0.015606` | `0.307777` | `2` |
| `8.0` | `2` | `0.011525`, `0.638527` | `0.179227` | `3` |

**The trap is worth naming, because it is easy to fall into and this round
nearly did.** A pole is a *sign change with no zero*. Counting sign changes
alone reports five states at `L = 8`; two of them are poles sitting on
`(2πj/L)²` and `(π(2j−1)/L)²` to machine precision. The separation is not
subtle once looked for — at a genuine root the residual is `~1e-15`, at a pole
`~1e+15` — so classifying by residual size is a fifteen-order-of-magnitude
discrimination and not a tuned threshold.

None of this touches the stability answer. Every one of these states is
**positive**, as the quadratic form requires. What it corrects is a count that
was reported as structural when it was parametric.

That is the same species of error the ledger already tracks — *a limit written
down as a fact* — in a new dress: **a parameter-range fact written down as a
structural one**. It is recorded rather than quietly patched.

## 8. What this ungates, and what is still put in

#260 gated stationary action and backreaction on whether its growing mode
survives a resolved mouth. #261 released that gate on a reduced model with the
balls left in. This round closes the gap that release stood on, and upgrades the
answer from a sign on a `2×2` to positivity of the energy itself.

```
ray closure → field solution → two-wave invariant → finite throat
            → mouth resolved ✓ → neck resolved ✓ → stationary action
            → backreaction → topological branch
```

What is **not** settled, and each with its shape stated rather than waved at:

* the tube still has **one transverse channel**, so `𝒜` is a coupling and the
  neck is a quantum-graph edge, not a solved cross-section;
* the reduced `2×2`'s cross term is **single-scattering**, with the neglected
  series measured at `9.5e-04` at the working point (§5);
* the ambient metric is **fixed** — a round `S³` with two balls cut out. Nothing
  here derives the geometry that would produce them;
* there is **no backreaction**.

Those are the next round's subject. They are not gates on it: an action computed
on this background measures the physics, because the background's positivity is
now proved rather than sampled.
