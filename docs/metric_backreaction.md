# Metric backreaction: does `A + B` do what rescaling `A` or `B` cannot?

**PR #263.** The gate #260 set is closed, so this is the roadmap's first GR
question — and deliberately *not* "does spacetime pinch off?".

> **Yes.** At the working point `0.921` of the interference response lies
> outside everything rescaling `A` or `B` can produce, and the interference term
> is **comparable in size** to the single-wave responses rather than a
> correction to them. The fraction is **not a universal constant** — it runs
> `0.56–0.99` across carriers and `0.88–1.00` across time windows — so it is
> reported as a range, and the claim is that it is *large everywhere*, not that
> it is any particular number.

> **Scope.** A linear conformally coupled scalar on a **fixed** ESU with #257's
> **point** throat, and the metric response to **first order** in the lowest
> (`n = 3`, homogeneous) transverse-traceless harmonic. Backreaction as a linear
> response, not a solved coupled system. **The mouths are not the resolved ones
> of #261/#262.**

---

## 1. Why this is a linear-algebra question

Three facts compose:

* the field equation is **linear**, so `φ_{A+B} = φ_A + φ_B`;
* `T` is **quadratic**, so `T[A+B] = T[A] + T[B] + ΔT` with `ΔT` bilinear;
* linearized Einstein is **linear in `T`**, so `β[A+B] = β_A + β_B + β[ΔT]`.

Rescaling `A → cA` sends `β_A → c²β_A`. So **everything reachable by rescaling
is the two-parameter family `{c²β_A + d²β_B}`**, and the question is exactly
whether `β[ΔT]` lies in it. That is a projection residual with a number attached.

The split is measured, not assumed: the self term scales as `c²` and the cross
term as `c` to machine precision, and the cross term is **exactly zero** with one
source switched off.

Note also that asking whether the *total* `β[A+B]` is reachable and whether the
*cross term* is are the same question — `β[A+B] − (αβ_A + γβ_B)` equals
`β[ΔT] − ((α−1)β_A + (γ−1)β_B)` — so nothing hinges on which is measured.

## 2. The channel, and why it is the only honest one

The **transverse-traceless** sector, lowest harmonic: the shear of the universe.
The first reason is the one that matters.

**The ESU is held static by a fluid this arc never specifies.** A perfect fluid
carries no anisotropic stress — in an orthonormal frame its
`T_ab = diag(ρ,p,p,p)` *whatever* the anisotropy — so it contributes nothing to
the traceless spatial part. Neither does `Λ`, which contributes `−Λη_ab`. **This
is the only channel whose answer does not depend on what was never put in.** The
scalar sector depends on the fluid's sound speed, and it is also where the
Eddington mode lives, so an `A`/`B`/`A+B` comparison there would be measuring the
background rather than the waves.

Three lesser reasons: TT perturbations are gauge-invariant, so no gauge choice
enters; it is the softest tensor channel, so if anything responds this does; and
it reduces the response to a five-component driven oscillator, which makes the
measurement a statement about the *source* rather than about a PDE solver.

## 3. The response, derived rather than recalled

Linearizing Einstein about the ESU in the homogeneous anisotropy
`a_i = a e^{β_i}`, `Σβ_i = 0`:

```
δG^TT_{ij}  =  β̈_{ij}  +  (8/a²) β_{ij}          ⟹   ω₃ = 2√2 / a
```

The connection is obtained by **solving** the torsion-free condition
`W^a_{cb} − W^a_{bc} = C^a_{bc}` rather than quoting a formula — because a first
draft *did* quote one and produced a `G₀₀` containing `ä`, which is impossible.
Hence three validations against known answers:

| | computed | expected |
| --- | --- | --- |
| round `S³` | `Ric = 2/a²`, `R = 6/a²` | ✓ |
| ESU | `G = diag(3,−1,−1,−1)/a²` | ✓ — and this independently reproduces `two_wave`'s `_EINSTEIN` |
| closed FRW | `G₀₀ = 3(ȧ²+1)/a²` | ✓ |

Two more facts fall out of the linearization itself, and each would have caught a
different error:

* the first-order piece is **automatically traceless** — so the anisotropy
  sources no trace and cannot mix with the fluid's `δρ`, `δp`;
* `δG_{0i} = 0` **identically** — the momentum constraint, so the mode is
  genuinely transverse and does not mix with vector or scalar perturbations.

And `ω₃² = +8/a² > 0`: **the tensor sector of the ESU is stable.** That is the
gate #260 taught this arc to check before computing anything on a background,
and it is why the comparison measures the waves.

## 4. What was actually hard: the quadrature

The source is `S_{ij}(t) = ⟨T^TT_{ij}⟩` over all of `S³`, and the integrand has
`1/χ⁴` singularities at **eight** points — the two sources, the two mouths, and
all four antipodes. Three things were needed.

**A batched solver.** The kernels do not depend on the observation point, only
`χ` and `n̂` do, so one batched FFT replaces `P` separate ones. Verified against
`two_wave.solve_field` to **`6e-16`**. Without it the 30 000-point rule the
convergence control needs is unaffordable.

**Exact angular rules.** Near a source `T ~ n_i n_j/χ⁴`, and the traceless part
of that has **zero** angular average. A rule that integrates `n_i n_j` exactly
kills the `1/χ⁴` exactly; a random rule leaves `O(1/√N)` of it. With exact rules
the traceless angular average is *flat* into the singularity — `2.44e-03` at
`r = 0.005` against `2.64e-03` at `r = 0.32` — so these integrals converge.

**A smooth partition of unity.** Excising the balls leaves the bulk integrand
discontinuous and the rule stalls near `1e-03`; the partition reaches `4e-04` and
keeps falling.

### And the mouths

The failure that took longest to find: the singular set was missing the **two
mouths and their antipodes**. Worse, `two_source.mouth_positions` puts the first
at `(1,0,0,0)` — exactly the axis a product grid naturally uses — so the
quadrature's own pole sat on a singularity, and **refining the grid made the
answer worse**. That is what the magnitude ratios were saying: `1.44`, `1.50`
between levels, growing rather than settling.

## 5. The control, which is the point

> **A first attempt at this round reported `0.982` unreachable. It was pure
> quadrature noise.**

Nothing about that run looked wrong. The tell appeared only on recomputing the
*same* quantity with an independent rule:

| | correlation between independent quadratures |
| --- | ---: |
| `β_A` | **`−0.04`** |
| `β_B` | `−0.13` |
| `β_ΔT` | `0.52` |

Not slightly off — meaningless. A residual of `0.982` is exactly what two
uncorrelated noise vectors give.

With all eight singular points handled and the partition of unity in place, the
same check now reads:

| refinement | `β_A` | `β_B` | `β_ΔT` |
| --- | ---: | ---: | ---: |
| correlation `7802 → 15988` | `0.9961` | `0.9895` | `0.9992` |
| magnitude ratio | `1.0009` | `0.9728` | `0.9946` |

with the residual moving `0.9186 → 0.9215` — a drift of `0.0029`.

These are the numbers under a **deterministic** quadrature basis, and that
mattered. A first version took the tangent basis from `np.linalg.svd`, which
returns *a* valid null-space basis but not a canonical one, so the whole rule
was platform-dependent — CI on another Python read the refined level as *less*
accurate than the coarse one. Replaced by a Householder reflection, the rule is
bit-identical everywhere, and the headline moved by `2.0e-05`: `0.921509` to
`0.921529`. The control came out **cleaner** rather than worse (worst
correlation `0.9704 → 0.9895`, drift `0.0076 → 0.0029`).

The volume check, separately, does **not** descend: the rule reaches `~1e-04`
and wanders there (`1.02e-02, 6.24e-05, 3.32e-04, 9.83e-05` down the ladder). A
`C^∞` partition bump was tried and does not help, so smoothness is not the
cause; the bulk grid spacing `π/26 ≈ 0.12` barely resolves the partition's
`0.267`-wide transition. The operative control is the one above — agreement
between *responses*, which does descend — and the plateau is stated rather than
papered over.

**The rule this round adopts: two refinement levels must agree in correlation
*and* in magnitude before any residual is quoted.**

## 6. The answer

| | |
| --- | ---: |
| unreachable fraction, window `(4,30)` | **`0.9215`** |
| `‖β_ΔT‖ / ‖β_A‖` | `1.026` |
| `|cos(β_A, β_B)|` | `0.171` |
| quadrature points | `15 990` |

The interference response is **comparable in size** to the single-wave ones and
nearly orthogonal to both. The two single-wave responses are themselves
independent (`|cos| = 0.171`), so the span is genuinely two-dimensional and the
residual is not large for a trivial reason.

The residual is reported off the full linear **span**, which strictly contains
the physical cone `{c²β_A + d²β_B}` — so the figure is **conservative**: the true
unreachable fraction is at least this large.

### It is not a universal constant

| window | unreachable | `‖β_ΔT‖/‖β_A‖` |
| --- | ---: | ---: |
| `(4, 12)` | `0.996` | `2.20` |
| `(4, 20)` | `0.964` | `1.28` |
| `(4, 30)` | `0.922` | `1.02` |
| `(8, 45)` | `0.879` | `0.96` |

| carrier | unreachable |
| ---: | ---: |
| `0.700` | `0.936` |
| `1.414` | `0.977` |
| `2.828` | `0.929` |
| `4.000` | `0.565` |

Large everywhere, constant nowhere. The honest headline is the **range**, and
the claim that survives it is that the interference response is *generically*
outside the reachable family — not that `0.92` is a property of nature.

## 7. The resonance test, done on the coupled system — which reverses it

A first version of this section argued the channel was off resonance **by
construction**. The argument: the conformally coupled scalar on the ESU has
spectrum `ω_n = n + 1`, **integers**; the space is compact and static so nothing
decays and the field rings on those modes forever; `T` is quadratic and integers
are closed under sums and differences, so a source built from the free field
rings on integers too — measured peaks at `5.969` and `7.959`, within `0.031`
and `0.041` of `6` and `8`. And `ω₃ = 2√2` is **irrational**, `0.172` from the
nearest integer.

**Every step of that is true, and the conclusion is false**, because all of it
is a statement about the *uncoupled ambient*.

The throat rings where `det(A − Γ(ω))` **vanishes**, and those zeros are not
integers. At the working point:

```
0.875   1.854   1.872   2.706   2.878   3.698   3.825   4.713   4.739 …
```

`ω₃ = 2.8284` sits **`0.050`** from the resonance at `2.878`. Near `ω₃` the
coupled resonances are spaced only about `0.17` apart, so the mode *cannot* be
far from one:

| `d` | boundary condition | nearest coupled resonance | distance to `ω₃` |
| ---: | :--- | ---: | ---: |
| `0.9` | `(0.30, 0.35, 0.06)` | `2.829473` | **`0.00105`** |
| `0.9` | `(0.10, 0.15, 0.06)` | `2.742287` | `0.08614` |
| `1.3` | `(0.30, 0.35, 0.06)` | `2.878498` | `0.05007` |
| `1.3` | `(0.60, 0.70, 0.06)` | `2.845978` | `0.01755` |
| `2.4` | `(0.30, 0.35, 0.20)` | `2.819321` | `0.00911` |

Across sixteen configurations the nearest is **always within `0.086`**, and at
`d = 0.9` with the working boundary condition it is `0.001` — resonant to the
width of the scan. `ω₄ = 3.873` is likewise `0.048` from one.

**So the corrected statement is a working-point one, and it points the other
way:** off resonance with the free ambient, **generically near-resonant with the
coupled system**. That is also the mechanism the first version lacked for §6's
carrier sensitivity — a response driven near a pole is sensitive, and it changes
*sign* across the pole.

### The species of error

An argument was established for one system and carried over to a system that
differs **precisely in the thing being argued about**. Coupling is what this
whole arc is about, and the spectrum is the first thing it moves. The earlier
draft's mistake was not a wrong number; it was not noticing that its premise had
a scope.

*(A still earlier draft said `T` being quadratic puts the source's power at
`2ω₀` and chose the carrier to match. Also wrong, and for a related reason — the
measured peak sits at `5.969` whatever the carrier is, because the ringing is
the background's rather than the pulse's.)*

## 7b. What this channel cannot say

A traceless shear preserves volume **exactly** — `det e^{β} = 1` when `Σβ_i = 0`,
at every amplitude, not just to first order — and the area of a geodesic sphere
**to first order**:

| `ε` | volume ratio | `δA/A` | `÷ ε²` |
| ---: | ---: | ---: | ---: |
| `0.200` | `1.0` | `+1.551e-02` | `0.3876` |
| `0.100` | `1.0` | `+3.963e-03` | `0.3963` |
| `0.050` | `1.0` | `+1.002e-03` | `0.4007` |
| `0.025` | `1.0` | `+2.518e-04` | **`0.4029`** |

The areal change is `O(ε²)` with the coefficient converging to `0.403`, and it
is **positive** — what second-order effect exists *opens* rather than pinches.

**So this channel cannot answer whether the interaction metric moves toward a
neck or away from one.** It distorts the mouth into an ellipse of the same area.
Quoting a shear amplitude as evidence about pinching would be a category error,
not a measurement error — and the areal sector on a resolved neck is a different
module.

## 7c. And neither can the rest of the tower

The natural next move from §7b is to build the inhomogeneous TT harmonics, on
the grounds that the *homogeneous* mode's areal blindness came from its
constancy. It does not. **The whole transverse-traceless sector changes the area
of no surface, at first order, exactly.**

For a surface with unit normal `n̂`, `δA/A = ½⟨tr h − h_nn⟩`, which is
`−½⟨h_nn⟩` for traceless `h`. For a TT plane wave `ε_ij cos(k·x)` the average
`⟨n_i n_j f(k·n̂)⟩` can only be `A δ_ij + B k̂_i k̂_j` by symmetry — and
contracting with `ε` kills the first term by **tracelessness** and the second by
**transversality**. Nothing is left.

| | worst relative `⟨h_nn⟩` |
| --- | ---: |
| TT plane waves, `a ∈ [0.05, 5]`, `\|k\| ∈ {1,4}` | **`4.5e-15`** |
| generic superpositions of 12 TT waves | `5.5e-17` |
| **control: tracelessness kept, transversality dropped** | **`1.5e-01`** |

The control is `3.3e+13` times larger, so the test has power — this is a
cancellation, not a quantity that happens to be small.

**And curvature does not spoil it,** because `S³` is maximally symmetric. The
moments of the outward normal over a geodesic sphere, in the global
left-invariant frame, are isotropic to `6e-16`, and `⟨n_i n_j n_k n_l⟩` matches
`(δδ + δδ + δδ)/15` to `2e-15` — at every radius from `0.05` to `1.2`, about
both a mouth and a source. Isotropic moments against a traceless transverse
tensor vanish on `S³` exactly as in flat space.

So building more harmonics adds more exact zeros. Asking the `n = 3` mode about
pinching was a category error; asking the whole tower is **the same category
error, one harmonic at a time**.

**A signed `δA/A` therefore has to come from the scalar sector** — which is
precisely where the fluid holding the ESU static enters, and where the Eddington
mode lives. §2 chose the TT sector because it is the only one whose answer does
not depend on what was never put in; the price of that choice, now measured, is
that it cannot answer a question about area. **The dependence cannot be worked
around by computing harder in this sector. It has to be named.**

## 8. Which branches were there

#257 measured the same configuration giving an invariant of anything from `0` to
`4` depending on which arrival branches were present, and concluded that any
quantity built by integrating over the field has to state them. Backreaction is
such a quantity, so the throat is switched off as a control:

| | unreachable | `‖β_ΔT‖/‖β_A‖` | `‖β_ΔT‖` |
| --- | ---: | ---: | ---: |
| with throat | `0.929` | `1.02` | `0.1213` |
| no throat | `0.986` | `9.46` | `0.1026` |

The conclusion **survives switching the throat off**, so this is a statement
about two waves rather than about the throat. What the throat changes is the
*single-wave* response — without it `β_A` is much smaller, which is why the ratio
jumps to `9.5` while `‖β_ΔT‖` barely moves.

## 9. What is still put in

* the **`n = 3` harmonic only** — the homogeneous shear, not the full TT tower;
* a **fixed** ESU background: the metric is not solved for;
* **point** sources and #257's **point** throat — *not* the resolved mouths of
  #261/#262. The traceless angular average happens to kill the `1/χ⁴` exactly, so
  the integrals converge, but resolving the sources is the obvious next
  refinement and it is the same finite-support move those rounds made;
* a **linear** response: first-order backreaction, not a solved coupled system.

```
ray closure → field solution → two-wave invariant → finite throat
            → mouth resolved ✓ → neck resolved ✓ → backreaction ✓
            → stationary action → topological branch
```

What this ungates is the **stationary action** round, which now has both a
positivity theorem for its background (#262) and a measured statement that the
two-wave configuration is not a rescaled one-wave configuration. What it does not
settle is whether the response *closes* — a solved coupled system, rather than a
first-order response on a fixed background, is the round after that.
