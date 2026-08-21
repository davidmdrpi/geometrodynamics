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
4. **The neck is an apparent horizon.** It is a minimal surface, and on a
   `K = 0` slice `θ_± = ±H`, so `H = 0` makes both null expansions vanish. The
   forced throat carries a horizon and is therefore **not traversable** — which
   is the vacuum face of §7's result that a traversable connection requires
   exotic matter. The throat with no exotic matter in it is precisely the one
   you cannot get through.

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
throat that is an Einstein–Rosen bridge — carrying an apparent horizon, and so
not traversable — is a separate question, and it is the one that decides whether
the mass law is a result about the model or only about this slice. Discharging
it means checking the horizon against whatever the framework requires of a
throat that carries signal, and that check has not been done here.
