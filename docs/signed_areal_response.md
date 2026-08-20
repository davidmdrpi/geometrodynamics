# The signed `ΔA/A`: the interference metric deforms toward a neck

**Module:** `geometrodynamics/waves/areal.py`
**Probe:** `python -m experiments.closure_ledger.areal_probe` (8/8)
**Tests:** `tests/test_waves_areal.py` (32)
**Renderer:** `python scripts/geometrodynamics_v64_areal.py --still v64.png`

---

## The question, and why it had to move sectors

PR #263 asked whether `A + B` produces a metric response that rescaling `A` or
`B` alone cannot, measured it in the transverse-traceless sector, and got a
large answer — `0.92` of the interference response unreachable. Then, asked for
a *geometric* verdict, it proved that sector cannot give one:

> `δA/A = −½⟨h_nn⟩`, and that average vanishes **identically** for any
> transverse-traceless field. With `⟨n_i n_j f(k·n̂)⟩ = A δ_ij + B k̂_i k̂_j`,
> tracelessness kills the first term and transversality kills the second.
> Building more harmonics adds more exact zeros.

So the question actually asked —

> Does the interaction metric deform **toward** a neck, **away** from one, or
> merely **oscillate**?

— had to be answered in the scalar sector. PR #263 had avoided that sector for
two specific reasons: it depends on a sound speed the arc never specified, and
it is where the Eddington mode lives. Posing it as an **initial-data** problem
rather than an evolution removes both. On a maximal slice (`K̄_ij = 0`) the `K`
terms in the Hamiltonian constraint are *quadratic*, so at the order the field's
stress tensor enters,

    δR⁽³⁾ = 16πG δρ

with no time derivatives in it. A constraint does not have a sound speed.

## The answer

    ΔA/A = ( −2.06 × 10⁻³ ,  −1.88 × 10⁻³ )     in units of 2πG

**Toward a neck. Both mouths close.** Negative in all eight control
combinations — two quadrature levels, two mouth radii, two gluings.

And the mechanism is not the obvious one. The interference energy *alone* would
**open** the mouths: `U(c_j) > 0` at both. What inverts it is the throat. The
contributions to `u` at the first mouth are

| piece | value |
|-------|-------|
| the source `U(c₁)` | `+3.44 × 10⁻³` |
| the two monopole layers | `−3.62 × 10⁻³` |
| the six dipole layers | `−4.2 × 10⁻⁹` |
| the free kernel element | `−3.43 × 10⁻⁴` |
| **total** | **`−5.14 × 10⁻⁴`** |

The throat is a low-impedance load for the constraint. It cannot support the
conformal factor the interference energy piles up around its mouths, and it
responds with negative monopole layers that more than cancel it.

## The setup

Writing `g_ij = ψ⁴ĝ_ij` with `ψ = 1 + u` on the unit `S³` (`R̂ = 6`):

    ∇²u + 3u = −2πG δρ ,        ΔA/A = 4u .

The conformal ansatz is not a restriction: a transverse-traceless piece
contributes nothing to `δR⁽³⁾`, and a longitudinal piece is a diffeomorphism of
a constant-curvature background. The conformal factor is the whole of what this
equation sees.

Two things make it hard, and both were established in
[`initial_data`](../geometrodynamics/waves/initial_data.py):

**The source diverges at the mouths.** `ΔT₀₀` goes as `1/χ⁴` there, with the
radial integrand *growing* into the singularity. PR #263's traceless projection
survived that because the angular average killed it; an energy density has no
such cancellation. Removing the balls removes the singularity from the domain.

**The operator is exactly degenerate.** `∇²u = −3u` is `−∇²u + u = 4u`, and
`4 = (n+1)²` at `n = 1`: the constraint operator sits **on** the ESU's dipole
level. Its kernel on `S³` is the four `k = 1` harmonics `x^A`, and the
interference source is not orthogonal to them.

## The representation, and the identity underneath it

    u = U + Σ_j [ A_j G_⊥(χ_j) + D_j·𝒟_j ] + c·x

with `U = −G_⊥[σ 1_Ω]`, `σ = −2πG δρ`, `𝒟_j` the dipole layer at mouth `j`,
`D_j ⊥ c_j`, and `c ∈ R⁴` free. Twelve unknowns: two monopole strengths, six
dipole components, four kernel coefficients.

The pseudo-inverse Green function satisfies

    L G_⊥ = −δ + (2/π²) cos χ ,     G_⊥(χ) = (π−χ) cos 2χ / (4π² sin χ)

and the trailing term is **exactly** the integral kernel of the orthogonal
projector onto `span{x^A}`: `Σ_A x^A y^A = cos χ(x,y)` and
`‖cos χ‖²_{S³} = π²/2`, so `1/‖cos χ‖² = 2/π²` and nothing else. That is
verified rather than asserted (`residual − tail = 2e-07`,
`1/‖·‖² − 2/π² = 0` to `1e-12`), because it is what turns the solvability
condition into an identity:

    Σ_j A_j c_j + Σ_j D_j = S_σ ,      S_σ = ∫_Ω y σ(y) dV .

Four equations. The remaining eight are the throat.

## The throat is a cavity

A tube of cross-section `𝒜` is a round `S²` of radius `r = √(𝒜/4π)` crossed with
a line, so `R̂ = 2/r² = 8π/𝒜`, and the same reduction gives `∇² + R̂/2`:

| tube channel | operator | behaviour |
|--------------|----------|-----------|
| `ℓ = 0` | `∂_s² + 4π/𝒜` | oscillatory, `k = 1/r` |
| `ℓ = 1` | `∂_s² − 4π/𝒜` | evanescent, **same rate** |

Continuity of `u` and of flux at each end closes the system. The first row is
the important one: **the tube is a resonant cavity for the constraint.** At
`kL = nπ` the response has a pole and the sign of `ΔA/A` flips — the scan finds
flips at `3.133` and `6.260` against closed-form poles at `π` and `2π`, both
within one grid cell. Past the first pole the two mouths can even move in
*opposite* directions.

The working throat (`𝒜 = 4π`, `L = 0.9`, so `k = 1` and `kL = 0.9`) is well
inside the first cell. **A signed answer here is a statement about a throat,
not about the sphere.**

## Four things that were checked, one that was corrected

### The assembly, against exact solves, in both sectors

The 12×12 matching is checked against one-dimensional boundary-value solves on
the punctured sphere that share no code with it. Both sectors, on purpose:
every closed-form check in PR #262 was `ℓ = 0`, which is exactly how an
`ℓ(ℓ+2)` error survived that round's entire suite. A suite that cannot see a
sector cannot defend it.

| sector | `a = 0.2` | `a = 0.1` | `a = 0.05` |
|--------|-----------|-----------|------------|
| `ℓ = 0` | `3.7e-11` | `7.7e-11` | `8.4e-12` |
| `ℓ = 1` | `5.1e-12` | `3.8e-11` | `3.9e-10` |

The error is **flat** in `a`. That is the assembly's numerical floor, not a
truncation error, so no order of convergence is claimed from it.

### The correction worth recording

The first assembly agreed with the reference at `1e-06` and no better — at
every radius, in both sectors. That number did **not move** when the
reference's quadrature tolerance, its finite-difference step, and its
boundary-value tolerance were each tightened by four orders of magnitude.

So it was not the reference's floor. It was a systematic error of the assembly:
the two-point stencil for the mouth-sphere radial derivative, whose *relative*
truncation error on a `1/χ` field is exactly `step²` — `1e-06` at `step = 1e-3`,
which is what was measured. A five-point rule removed four orders.

> **A discrepancy that refuses to move when you refine the other side is the
> other side telling you it is not the problem.**

### The dipole layers: required, and then invisible

`A₁c₁ + A₂c₂` sweeps only the plane of the two mouths, so two monopoles can
meet at most two of the four solvability equations. Measured, the monopole-only
condition fails by **62.5%** of the obstruction: without the `ℓ = 1` layers
there is **no solution at all**. This was expected to be the heart of the round.

It is not. The response to the off-plane obstruction is `6 × 10⁻¹⁷` — zero. The
layers deposit it in the kernel elements `x²` and `x³`, which vanish at both
mouths, because both mouths lie in the `(x⁰, x¹)` plane.

A first draft of the module docstring said the `ℓ = 1` sector *is* the
calculation. It is required for the calculation to **exist** and contributes
nothing to its **value**.

### What the answer is actually made of, and why that is lucky

| variant | `ΔA/A` mouth 1 | vs. full |
|---------|----------------|----------|
| full | `−2.0575e-03` | — |
| obstruction only, no local source data | `−2.0556e-03` | `0.09%` |
| local source data only, no obstruction | `−1.9077e-06` | `0.09%` of full |
| `ℓ = 1` moments zeroed | `−2.0579e-03` | `2.4e-04` |
| `ℓ = 1` moments randomised (±50%) | `−2.0575e-03` | `4e-05` |
| obstruction doubled | `−4.1112e-03` | exactly `×2` |

`ΔA/A` is, to `0.09%`, a **linear functional of the obstruction alone** — of
its in-plane part alone, exactly.

That matters because the `ℓ = 1` source moments are the worst-converged input
the source quadrature produces: **41%** drift between quadrature levels, against
**1.5%** for the obstruction. The signed answer rests on the single
best-converged number available. That is not how these rounds usually go, and it
is worth saying out loud rather than leaving as luck.

## The controls

| radius | points | gluing | `ΔA/A` mouth 1 | `ΔA/A` mouth 2 |
|--------|--------|--------|----------------|----------------|
| 0.05 | 5158 | transported | `−2.1014e-03` | `−1.9148e-03` |
| 0.05 | 5158 | reflected | `−2.1014e-03` | `−1.9148e-03` |
| 0.05 | 12630 | transported | `−2.0575e-03` | `−1.8818e-03` |
| 0.05 | 12630 | reflected | `−2.0575e-03` | `−1.8818e-03` |
| 0.10 | 5158 | transported | `−1.2247e-03` | `−1.0899e-03` |
| 0.10 | 5158 | reflected | `−1.2247e-03` | `−1.0899e-03` |
| 0.10 | 12630 | transported | `−1.1783e-03` | `−1.0547e-03` |
| 0.10 | 12630 | reflected | `−1.1783e-03` | `−1.0547e-03` |

Quadrature spread at fixed radius: `2.1%`. Worst condition number `2.1e+05`,
worst residual `2.8e-15`.

**The mouth radius is not a regulator.** The source goes as `1/χ⁴` at the
mouths; `a` is a parameter of the throat; there is no `a → 0` limit to converge
to — that is precisely the singular point PR #261/#262 removed. Doubling `a`
moves the answer by a factor `1.75`, almost exactly the factor by which the
source's own mouth-weighted moment moves, because a bigger ball excludes more of
the `1/χ⁴` pile-up.

**The gluing cannot reach the answer.** Identifying the two mouths' transverse
frames through the tube is a genuine modelling choice — a handle may be glued
with a twist. A full `2π` sweep of that twist moves the individual dipole
strengths by `1%` and `ΔA/A` by less than `1e-12`, because the areal response
sees the dipoles only through `Σ_j D_j`, which the solvability rows pin.

## What is still put in

* The **fluid** holding the ESU static is held rigid (`δρ_fluid = 0`). That is
  consistent because the scalar's stress tensor is separately conserved, but a
  responsive fluid is the obvious next refinement.
* The **source** is PR #263's: a linear conformally coupled scalar on a *fixed*
  ESU with *point* sources, with the throat entering the field solve as PR
  #257's point throat.
* The source's **channel data** come from the local expansion of `U` about each
  mouth centre — exact through `O(a²)` in `ℓ = 0` (the correction is free:
  `∇²U = −3U − (2/π²) S_σ·x` inside the ball, where the source does not reach)
  and leading order in `ℓ = 1`. Bounded from above by the fact that the answer
  barely depends on those data at all.
* The response is **linear** in `G` and quadratic in the wave amplitude. Note
  that the kernel coefficient comes out at `|c| ≈ 1.7` against a mouth response
  of `2e-03`: the solution is dominated by a nearly-zero mode that is large far
  from the throat and small at it. The perturbative window is set by that
  region, not by the mouths.

## What this closes, and what it opens

The roadmap's line was `mouth resolved ✓ → neck resolved ✓ → backreaction ✓ →
stationary action`. This round adds the geometric verdict that PR #263 could
not give: **the interference deforms the neck toward closing, at the working
throat, off resonance.**

What it opens is the resonance structure. The sign is not universal — it is a
property of the throat's aspect ratio, with poles at `kL = nπ`. Whether a
physical throat can be driven onto one of those poles, and what happens there,
is a question this round can pose precisely and does not answer.
