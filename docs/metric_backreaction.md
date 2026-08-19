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
| correlation `7802 → 15990` | `0.9942` | `0.9704` | `0.9991` |
| magnitude ratio | `0.9927` | `1.0094` | `0.9937` |

with the residual moving `0.9291 → 0.9215` — a drift of `0.008`. A third level
at `32 496` points gives `0.9173`, and correlations that *rise* toward `1`
(`0.9988`, `0.9979`, `0.9995`).

**The rule this round adopts: two refinement levels must agree in correlation
*and* in magnitude before any residual is quoted.**

## 6. The answer

| | |
| --- | ---: |
| unreachable fraction, window `(4,30)` | **`0.9215`** |
| `‖β_ΔT‖ / ‖β_A‖` | `1.021` |
| `|cos(β_A, β_B)|` | `0.173` |
| quadrature points | `15 990` |

The interference response is **comparable in size** to the single-wave ones and
nearly orthogonal to both. The two single-wave responses are themselves
independent (`|cos| = 0.17`), so the span is genuinely two-dimensional and the
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

## 7. The channel is never on resonance

A second-order fact that came out of trying to justify the carrier, and is more
interesting than the justification was.

The conformally coupled scalar on the ESU has spectrum `ω_n = n + 1`:
**integers**. The space is compact and static, so nothing decays and the field
rings on those modes forever. `T` is quadratic, and integers are closed under
sums and differences, so the shear source rings on integers too — measured, its
peaks land on integers to within a grid bin:

| measured peak | nearest integer | offset |
| ---: | ---: | ---: |
| `2.094` | `2` | `0.094` |
| `3.142` | `3` | `0.142` |
| **`5.969`** | `6` | `0.031` |
| `6.912` | `7` | `0.088` |
| **`7.959`** | `8` | `0.041` |

And `ω₃ = 2√2 = 2.8284` is **irrational** — `0.172` from the nearest integer, and
it cannot ever coincide with one. **So on this background the gravitational shear
channel is driven off resonance by construction, whatever the source does.**

### A wrong first draft, recorded

That draft said: `T` is quadratic, so a carrier `ω₀` puts the source's power at
`2ω₀`; therefore pick the carrier to put `2ω₀` at `ω₃`. **It is wrong.** The
measured spectral peak is `5.969` for carriers `0.7`, `1.414`, `2.0` and `2.828`
alike — it does not follow the carrier at all, because the ringing is the
*background's*, not the pulse's. (At `carrier = 4` the envelope's DC lobe takes
over and the largest bin moves to the bottom of the band; the ringing peaks are
still on integers. So "the peak never moves" would also have been too strong, and
what is asserted is the claim that survives: no peak of the source ever lands on
`ω₃`.)

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
