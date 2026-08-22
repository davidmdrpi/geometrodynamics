# Which throat is physical — and the sign reverses on the one that is forced

**Module:** `geometrodynamics/waves/physical_throat.py`
**Probe:** `python -m experiments.closure_ledger.physical_throat_probe` (7/7)
**Tests:** `tests/test_waves_physical_throat.py` (45)
**Renderer:** `python scripts/geometrodynamics_v65_throat.py --still v65.png`

---

## The question, and why it could not be deferred

PR #264 measured a signed `ΔA/A` on the throat this arc has carried since
PR #261 — a product tube of cross-section `𝒜 = 4π` and length `L = 0.9` — and
found that the interference **closes** both mouths. Then, asked to match the
tube's area to its own mouths', it found the sign **reverses**.

That turned the throat's geometry from decoration into the load-bearing part of
the answer:

> `𝒜` and `L` were free parameters. Which values are physical?

## They were never free

On a maximal slice the **background** Hamiltonian constraint is `R̂ = 16πG ρ̄`.
A profile therefore does not get to choose its matter — **the matter is
whatever the profile implies.** For a rotationally symmetric slice
`ds² + f(s)²dΩ²`,

    R̂ = 2/f² − 4f''/f − 2f'²/f²

which is derived here from the metric (Christoffels contracted, not a
remembered formula) and checked against two cases with known answers: the round
`S³` (`f = sin s`) must give exactly `6`, and the vacuum profile exactly `0`.
Both hold to machine precision.

A product tube of area `𝒜` is `S²(r) × R` with `r = √(𝒜/4π)`, so `R̂ = 8π/𝒜`:

| throat | `𝒜` | `R̂` | `ρ_tube/ρ̄` |
|--------|-----|-----|-------------|
| PR #261–#264, wide | `4π` | `2` | **`1/3`** |
| matched to the mouths (`a = 0.05`) | `4π sin²a` | `800.7` | **`133`** |
| the one area that carries the ambient fluid | `4π/3` | `6` | `1` |

Neither area in use was chosen for a reason, and neither holds the ambient's own
fluid. And the third row is not the answer either: two regions of equal constant
`R̂` joined without a surface layer are *one* region of that `R̂`, and the only
rotationally symmetric one is the round sphere — which has no neck.

## The throat that is forced

Ask instead for a throat that needs **no matter at all** (`R̂ = 0`) and glues on
with **no surface layer**. `R̂ = 0` integrates once:

    f'² = 1 + C/f ,   and a neck f'(f₀) = 0 fixes C = −f₀ ,
    so  f'² = 1 − f₀/f .

Smooth gluing to the round `S³` at mouth radius `a` requires `f' = cos a` where
`f = sin a`, and that forces

    cos²a = 1 − f₀/sin a   ⟹   **f₀ = sin³a** .

Two conditions, two unknowns, **nothing left over.** Length and resistance then
follow in closed form, each verified against quadrature to `1e-12`:

    L = 2[ sin³a · arccosh(1/sin a) + sin a cos a ]   ≈ 2a
    I = ∫ ds/f² = 4 cos a / sin³a

| `a` | `f₀ = sin³a` | `L` | `L/a` | `I` | `(4π/I)/N₀` |
|-----|--------------|-----|-------|-----|-------------|
| 0.02 | `7.998e-06` | `0.040063` | `2.003` | `5.000e+05` | `0.250000` |
| 0.05 | `1.248e-04` | `0.100754` | `2.015` | `3.200e+04` | `0.250000` |
| 0.10 | `9.950e-04` | `0.204629` | `2.046` | `4.000e+03` | `0.250000` |
| 0.20 | `7.841e-03` | `0.425477` | `2.127` | `4.999e+02` | `0.250000` |

Two things worth noticing. **The throat is short** — `L ≈ 2a`, so `0.10` at the
working radius, not the `0.9` carried since PR #261. And the conductance is
*exactly* a quarter of the exterior's own monopole stiffness,
`4π/I = N₀(a,4)/4`, at every mouth radius, with no free parameter in the
relation.

## It is an Einstein–Rosen bridge, and that derives its mass

`R̂ = 0`, `K = 0` and a spherical neck do not merely *permit* an Einstein–Rosen
bridge. They are one.

Write the time-symmetric Schwarzschild slice in proper radial distance:

    ds² = dr²/(1 − 2M/r) + r²dΩ²   ⟹   df/ds = √(1 − 2M/f) ,  f = r .

That is `f'² = 1 − f₀/f` with **`f₀ = 2M`**. The forced neck radius *is* twice a
mass, so the gluing condition `f₀ = sin³a` reads

    **M = sin³a / 2** ,

the throat's mass derived from the size of the excised mouth, with nothing left
to choose. Small mouths give `M ≈ a³/2`; inverting, `a = arcsin((2M)^{1/3})`.

**Three independent quasi-local masses agree on it exactly** (to `1.3e−13`, the
cancellation floor):

| quantity | value |
|----------|-------|
| Schwarzschild parameter `f₀/2` | `M` |
| irreducible mass `√(A_neck/16π)` | `M`, since `A_neck = 16πM²` exactly |
| Hawking mass `(f/2)(1 − f'²)` | `M`, **constant** along the whole vacuum piece |

### And the gluing condition *is* a mass statement

The Hawking mass of a round sphere in `ds² + f²dΩ²` is `(f/2)(1 − f'²)`. On the
throat that is `f₀/2` at every radius. On the round `S³` at geodesic radius `χ`
it is `sin³χ / 2`. Setting them equal at `χ = a` is, term for term, the equation
`f₀ = sin³a`.

> The seam is smooth exactly when the Hawking mass does not jump across it.

That is a much better way to see the result than "solve for `f₀`": the mouth
radius `a` and the throat mass `M` were never two parameters, because the
ambient already assigns a mass to every sphere it contains, and the throat has
to carry the one at the cut.

### Four things this does not say

The claim is strong enough to be worth not overstating, so each of these is
asserted in the tests rather than left as prose.

1. **It is a truncated bridge.** The geometry is Schwarzschild only between `f₀`
   and `sin a`; past that it is the ESU. There is no asymptotic region, so the
   **ADM mass is not defined here.** What is derived is quasi-local — and it is
   unambiguous only because the Hawking mass happens to be *constant* on the
   vacuum piece, so there is no question of where to evaluate it.
2. **The mass is dimensionless.** It is `M/R` against the ESU's curvature
   radius. No absolute unit is claimed, which is what PR #52's scale-modulus
   theorem requires of anything derived in this framework — a ratio is
   available; a unit is not.
3. **It is a handle, not a bridge between universes.** Both ends sew into the
   *same* `S³`, which is Misner's construction rather than the two-sheeted
   Einstein–Rosen picture.
4. **The neck is a marginal sphere of this slice — and that is all.** The
   result is an identity:

   > `H = 0` at the neck, `K_ij = 0` on the slice  ⟹  `θ_+ = θ_− = 0` .

   Both null expansions vanish, so the neck is a minimal surface and a
   marginally outer trapped sphere (a MOTS) **of the initial data**. That is a
   statement about one surface in one slice, and it is the whole of what is
   established here.

   Two stronger things are *not* established, and earlier drafts of this
   document asserted them directly.

   - **It is not shown to be an apparent horizon.** An apparent horizon is the
     **outermost** MOTS of a slice. Deciding that means examining every
     surface enclosing this one — a global condition on the slice, which
     nothing in this round evaluates.
   - **It is not shown to be non-traversable.** Traversability is a property
     of the **Lorentzian development**, and everything here is spatial initial
     data with no lapse chosen. Non-traversability *does* follow if one
     additionally takes the development to be the standard vacuum
     Schwarzschild / Einstein–Rosen one — the natural reading, since the slice
     *is* the time-symmetric Schwarzschild slice — and it is under that added
     assumption, not from the data alone, that this connects to §7's result
     that a traversable connection requires exotic matter.

   The identity is exact and needs no qualification. The conclusions drawn
   from it need the conditions named above, and the difference between the two
   is the correction.

## There is no cavity

The constraint operator is `∇² + R̂/2`. With `R̂ = 0` the tube carries the **plain
Laplacian**, and the `ℓ = 0` channel is

    (f²u')' = 0   ⟹   u = A + B∫ds/f²   — monotone.

No standing wave. No poles. Nothing for a sign to flip across.

PR #264's cavity, its resonances at `kL = nπ`, and the sign flips across them
were **properties of matter in the tube**, not of a throat. That is worth
stating plainly because that round's headline was qualified as *off resonance*,
and on the physical throat there is no resonance to be off.

## The mechanism: one number decides the sign

Any symmetric two-port splits as

    Y = G·[[−1, 1], [1, −1]]  +  shunt·[[1, 0], [0, 1]]

— a **conductance** through the tube and a **shunt** into it. The shunt is
`Y·(1,1)`: the flux a *uniform* potential drives into the throat. For a vacuum
tube `(f²u')' = 0` makes flux conservation an identity, so the shunt vanishes
**exactly** — there is nowhere to put that flux. A tube with matter in it has
somewhere, and shunts.

Scanned separately, the two behave completely differently:

| conductance | shunt | `ΔA/A` mouth 1 | verdict |
|-------------|-------|----------------|---------|
| `3.927e-04` (vacuum) | `0` | `+6.63750e+00` | opens |
| `1.604e+01` (wide) | `0` | `+7.60933e+00` | opens |
| `3.927e-04` (vacuum) | `6.070` | `−1.59264e-03` | closes |
| `1.604e+01` (wide) | `6.070` | `−2.05839e-03` | closes |

The **conductance**, scanned over eight orders of magnitude, never changes the
sign. The **shunt** passes through a pole between `1e-03` and `3e-03` and flips
it. PR #264's tube sat at `6.07` — three orders past.

> The sign of `ΔA/A` is not about how well the throat conducts. It is about
> whether there is anywhere inside it to put monopole flux — which is to say,
> whether there is matter in it.

## The answer

    ΔA/A = ( +6.64 , +8.58 )     in units of 2πG

**Away from a neck. Both mouths open** — the reverse of PR #264's sign on its
own tube.

| control | result |
|---------|--------|
| two quadrature levels | opens, spread `2.4%` |
| two mouth radii (`0.05`, `0.10`) | opens |
| two gluings (transported, reflected) | opens, identical |
| the whole vacuum family, four orders in `f₀` | opens |

That last one matters most. The smooth-gluing condition is what removes the
last freedom, so the answer had better not depend on hitting it exactly. Scanned
from `f₀ = 0.02 sin³a` to `f₀ = 300 sin³a` — with the glued point in the middle
— the sign never changes.

## What it costs, stated plainly

The response is **3000×** PR #264's and grows as **`a⁻³`**, and the matching
system's condition number grows with it (`3.1e+07` at `a = 0.05`, against
`4.5e+05` for the wide tube).

That is not noise, and it is not a numerical complaint. It is the physics of a
throat with zero shunt *by identity*: it does almost nothing to lift the
constraint operator's exact degeneracy, so the linear response sits close to a
mode the operator nearly annihilates. The near-zero mode is the `k = 1` kernel
that `initial_data` identified — `4 = (n+1)²` at `n = 1` — and a vacuum throat
barely separates it from zero.

**The sign is robust. The window in which linearising was legitimate is now the
binding question**, and this round does not answer it. A response of order
`10⁰ × 2πG` against a conformal factor `ψ = 1 + u` is a statement about a
coefficient, not about a geometry that has moved by that much — but it does mean
the amplitude at which the linearisation can be trusted is set far from the
throat, in the region where the near-zero mode is large.

## Release hardening: the two-port in closed form, and where the solve stopped

*Added after the round closed. This is a correctness repair to an implementation
this round shipped, not a new result about the geometry — the geometry is
unchanged and so is the answer above.*

**The static problem on this profile has a closed form.** `f'² = 1 − f₀/f` is
`f = f₀cosh²x` with `ds = 2f dx`, which turns

    (f²u')' = ℓ(ℓ+1) u    into    R_xx + 2 tanh x · R_x − 4ℓ(ℓ+1) R = 0 ,

and the substitution `R = y/cosh x` reduces *that* to `y'' = (2ℓ+1)² y` — a
constant-coefficient equation. The half-length is `X = arcosh(1/sin a)`, whose
exponential is exactly `e^{−X} = tan(a/2)`, so with

    k = 2ℓ+1 ,   q = e^{−2kX} = tan^{2k}(a/2) ,   f_m = sin a

every hyperbolic function collapses to a rational function of `q`
(`coth 2kX = (1+q²)/(1−q²)`, `csch 2kX = 2q/(1−q²)`, `tanh X = cos a`), and the
two-port is

> `D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]` ,
> `C_ℓ = +4π sin a · kq/(1−q²)` ,   `Y_ℓ = [[D, C], [C, D]]` .

That is now `VacuumThroat.admittance`. `ℓ = 0` is special-cased to
`(π sin³a / cos a)·[[−1,1],[1,−1]]` — the same number to `1e-17` relative, but
written that way the zero-shunt identity `Y·(1,1)ᵀ = 0` holds **to the last
bit** instead of to a tolerance, and the general diagonal's own near-cancellation
(`k coth 2kX − tanh X`, two quantities both tending to `1` as `a → 0`) is
avoided.

**Where the Riccati solve stopped resolving.** The solve is retained, as
`admittance_riccati`, and demoted to a validator rather than deleted — the
closed form needs something independent to be checked against. What is wrong
with it is its last step: it forms

    Y₁₂ = ½(s − t)

from the two eigenchannel values. Both are `O(1)` once `ℓ ≥ 1`, and they agree
to more and more digits as the neck closes, so the cross term ends up made
entirely of the cancelling tail. At `a = 0.05`:

| `ℓ` | `C_ℓ` closed form | `C_ℓ` Riccati | rel. error | signs | diagonal |
|--|--|--|--|--|--|
| 0 | `+3.92699e-04` | `+3.92699e-04` | `7.5e-14` | agree | `7.5e-14` |
| 1 | `+4.60578e-10` | `+4.60517e-10` | `1.3e-04` | agree | `9.7e-14` |
| 2 | `+3.00105e-16` | `−1.16573e-14` | `3.98e+01` | **opposite** | `9.4e-15` |
| 3 | `+1.64257e-22` | `−7.54952e-15` | `4.60e+07` | **opposite** | `4.0e-15` |

The **diagonal was never affected** — it is dominated by the mouth term and the
two routes agree to `1e-14` in every channel. The floor is around `1e-12` in
`|Y₁₂|`.

The `39×` at `ℓ = 2` is deliberately **not** pinned by a test. That factor is
one solver's step sequence in one build of SciPy and would move under any of
them. What the tests pin is the boundary: the sign is wrong, the magnitude is
more than an order of magnitude out, and the honest answer lies below the
solver's floor. A difference of two numbers cannot carry information the
numbers have already lost.

**The answer above is unchanged.** `solve_matching` consumes only `ℓ = 0` and
`ℓ = 1`, both above the floor; swapping the closed form in moves `ΔA/A` in the
**thirteenth digit** (`6.637503334434484 → 6.637503334434201`). That is checked
directly rather than asserted — if it had moved, the old number would have been
carried by the broken term.

**What the closed form then shows, for free.** As `a → 0`,

    C_ℓ  →  4π(2ℓ+1) sin a · tan^{4ℓ+2}(a/2)   ~   a^{4ℓ+3} ,

verified to `2e-16` relative already at `a = 0.05`, with fitted exponents
`3.000000`, `7.000004`, `11.000009`, `15.000013`. **Each unit of angular
momentum costs four powers of the mouth radius.**

| `a` | `C₀` | `C₁` | `C₀/C₁` | `C₁/D₁` |
|--|--|--|--|--|
| `0.05` | `3.926992e-04` | `4.605780e-10` | `8.5262e+05` | `7.33e-10` |
| `0.10` | `3.141614e-03` | `5.910168e-08` | `5.3156e+04` | `4.70e-08` |
| `0.20` | `2.513546e-02` | `7.641316e-06` | `3.2894e+03` | `3.03e-06` |
| `0.30` | `8.487017e-02` | `1.327738e-04` | `6.3921e+02` | `3.50e-05` |

That connects to the ESU kernel. The four `n = 1` harmonics `x^A` are degenerate
on the round `S³`, but a throat cut at one point splits them by *local* angular
momentum about that point: `X⁰ = cos χ` is `ℓ = 0` there and the three
`Xⁱ = sin χ n̂ⁱ` are `ℓ = 1`. The kernel therefore splits `1 ⊕ 3`, and the two
pieces cross the throat four powers of `a` apart.

**The statement this supports, and the only one:**

> The static scalar Laplacian on this scalar-flat spatial throat suppresses the
> local `ℓ = 1` mouth-to-mouth channel by `~10⁻⁹` at `a = 0.05`, while
> preserving a much stronger monopole channel.

It is **not** a statement that orientation information cannot cross the throat.
This is one operator (the static scalar Laplacian) on one slice — and, per the
scope note below, `physical_throat` supplies spatial initial data only, with no
lapse chosen, so there is no dynamical problem here to make such a claim about.
The `ℓ = 1` channel is **small, not zero**.

**Two scope corrections to this round's own prose, while here.** `f_min > 0` is
forced *within the class this round works in* — spherically symmetric,
scalar-flat, `K_ij = 0`, `C¹`-matched — not by Einstein's equations in general;
and everything in this file is **spatial initial data**. The dynamic problem is
not well-posed until a lapse is chosen, and none is.

## What is still put in

* The **source** is PR #263's: a linear conformally coupled scalar on a *fixed*
  ESU with *point* sources.
* The ESU's **fluid** is held rigid.
* The **exterior** is the round `S³` with two balls removed. The `C¹` gluing is
  exact for the profile, but the ambient's own response to carrying a handle is
  not modelled.
* A **vacuum** tube is the simplest acceptable matter, not the only one. Any
  `ρ(s) ≥ 0` that glues smoothly is admissible — a function's worth of freedom.
  The vacuum throat is distinguished by needing no matter at all and by having
  no free parameter, and the sign result is stated for it.

## What this closes, and what it opens

The roadmap's line reads `mouth resolved ✓ → neck resolved ✓ → backreaction ✓ →
signed ΔA/A ✓ → which throat ✓ → stationary action`. The throat's cross-section
has been a free parameter since PR #261, carried because nothing measured had
ever depended on its value. PR #264 made something depend on it; this round
removes it.

What it opens is sharper than what it closed, in two directions.

**The perturbative one.** The physical throat barely lifts the constraint's
degeneracy, so the next obstruction is not geometric: at what amplitude is the
linearised Lichnerowicz problem still the right problem? That question has a
definite answer and this arc has not asked it.

**And the one the mass law opens.** `M = sin³a / 2` ties the throat's mass to
the mouth's angular size with no free parameter — a *dimensionless* relation
between the throat and the cosmological radius, which is the only kind PR #52's
scale-modulus theorem permits. Whether this arc's larger framework can accept a
throat that is an Einstein–Rosen bridge — one whose neck is a marginal sphere
already on the initial slice — is a separate question, and it is the one that
decides whether the mass law is a result about the model or only about this
slice.

Discharging it takes two steps, and only the second is about the framework.
First the geometric question this round left open: **is this MOTS the outermost
one, and what is the Lorentzian development?** Neither is settled by the data —
the first is a global condition on the slice, the second needs a lapse. If the
development is the standard vacuum Schwarzschild/Einstein–Rosen one then the
throat does not carry signal and the framework has to accept that; if it is
not, the question stays open. Only then does the second step — checking the
answer against whatever the framework requires of a throat that carries signal
— have something definite to check. Neither step has been done here.
