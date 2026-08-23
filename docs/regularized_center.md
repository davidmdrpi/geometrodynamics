# The centre as a finite bearing, not a point

**Module:** `geometrodynamics/viz/regularized_center.py`
**Probe:** `python -m experiments.closure_ledger.regularized_center_probe` (13/13)
**Tests:** `tests/test_viz_regularized_center.py` (110)
**Renderer:** `scripts/geometrodynamics_v69_regularized_center.py`

---

## The proposal, and what it is trying to keep

Every picture in this arc has put a **point** in the middle. That is not
decoration — the point is doing work. Two radial arms

    P_A → O → P_B

can change direction at `O` for free, because at a point there is no angular
direction left to change. So the link's cost does not depend on where the
mouths sit, which is the property the connection is wanted for. And it is
bought with `f = 0`, which is where the geometry stops existing.

Regularise it. Keep the metric rotationally symmetric,

    dℓ² = ds² + f(s)² dΩ² ,

and replace `f(0) = 0` with

> **`f_min = f₀ > 0`.**

In the 2-D cross-section the middle is now a small **circle**. In a
`d`-dimensional bulk it is the space of radial directions, `S^{d−1}` — or
`RP^{d−1}` if the clock hand is an unoriented *axis*, `n ∼ −n`. Nothing is
singular, and three things the point picture forced become free.

The worked example throughout is the scalar-flat profile `f'² = 1 − f₀/f`.
That is not an invention: it is the profile PR #265 **forced**, which is what
makes every closed form below checkable against `physical_throat` instead of
asserted.

## The arms are the repo's own geometry, with the symmetry dropped

The proper distance from the neck out to an end of transverse scale `F` is

    L(F) = √(F(F − f₀)) + f₀ arcosh√(F/f₀) ,

and one arm's share of the resistance is

    I(F) = (2/f₀)√(1 − f₀/F) .

Both check against quadrature to `1e-9`. The check that matters more is the
other one: setting `f_o = f_i = sin a` with `f₀ = sin³a` reproduces
`VacuumThroat.length()` and `.resistance()` **bit for bit** — difference
`0.0e+00`. The regularised centre is not a new geometry. It is the forced
throat's own profile, with the two arms allowed to differ.

**That is an *intrinsic* statement, and it needs a caveat the first draft did
not carry.** Unequal arms are honest truncations of one scalar-flat profile —
nothing about the interior objects. But PR #265's mouths are `C¹`-matched to a
**unit round `S³`**, giving `F_j = sin a_j` and `|f′_j| = cos a_j`; with
`f′² = 1 − f₀/F_j` that forces `f₀ = sin³a_j = F_j³` at each end. One common
`f₀` then fixes `F = f₀^{1/3}` uniquely, so **both matched mouths have the same
scale**. Asymmetric arms therefore require asymmetric ambient attachments —
different curvature scales, a shell or junction, nonzero `K` — and are *not* a
fully matched #265 throat with the symmetry merely removed.

**And they can differ.** There is no requirement that `f_o = f_i` or that the
two proper distances match:

| `f_i` | `L_o` | `L_i` | `L_o/L_i` | `T(1 rad)` | `α²/(2I)` |
|--|--|--|--|--|--|
| `0.002` | `1.003647` | `0.002296` | `437.21` | `1.4461e-04` | `1.4649e-04` |
| `0.01` | `1.003647` | `0.011305` | `88.78` | `1.2718e-04` | `1.2832e-04` |
| `0.1` | `1.003647` | `0.102492` | `9.79` | `1.2430e-04` | `1.2535e-04` |
| `1.0` | `1.003647` | `1.003647` | `1.00` | `1.2403e-04` | `1.2506e-04` |

Nothing pathological happens at a ratio of 437. The old vacuole picture had to
give the inner and outer boundaries **one shared arbitrary radial gap**; a
bearing with two arms has no such constraint.

One correction to how that was first put. `R_inner/R_outer` does *not* stop
carrying physical significance — the next section derives `w_i/w_o = f_i/f_o`,
so the **endpoint scale ratio is load-bearing** for scale transport. What loses
its significance is the vacuole's **arbitrary drawing ratio**, which was never
identified with the endpoint scales at all.

## Scale transport becomes explicit

A feature of angular width `Δθ` has physical width `w(s) = f(s)Δθ`, so along
the route

    w_o = f_o Δθ  →  w_min = f₀ Δθ  →  w_i = f_i Δθ ,

with `w_i/w_o = f_i/f_o` and `f₀` nowhere in it. **A packet does not teleport
from one radius to another at fixed size** — it is squeezed into the bearing
and let out again. What is transported is the *angular* width, which is
constant along the route because the metric's angular part is `f²dΩ²` and the
coordinate does not change; the physical width is not.

One bookkeeping note, because this arc has made the mistake before. `w_i/w_o =
f_i/f_o` is an identity, but it is not asserted as a *bitwise* one: the drawn
ratio is `(f_i Δθ)/(f_o Δθ)`, whose two products round separately, and it sits
one ulp away at the first row of the scan. The identity is exact in the field;
the picture carries the rounding of the widths it actually drew.

## The result: the turn cost is quadratic, not linear

The natural guess is that turning through `α` around a bearing of radius `f₀`
costs its arc, `f₀α`, giving `L_throat(α) ≃ L_o + L_i + f₀α`. **That is exactly
right for one route** — hold direction all the way down, turn on the bearing,
hold it all the way up — and `corner_route_length` computes it. But that route
is not the geodesic.

**Why the corner loses is Pythagoras, not leverage.** The corner pays its
angular displacement as *pure transverse motion* at the neck — `f₀α`, first
order in `α`. The geodesic instead **tilts** motion that is already happening
radially, and a tilt costs `√(ds² + f²dθ²) − ds ≈ ½f²θ′²ds` — *second* order.
Converting a first-order cost into a second-order one is the whole of the
saving.

*Where* the turn goes then follows from minimising `∫½f²θ′²ds`, which gives
Clairaut's `θ′ = h/f²` — and that puts the angular rate **highest where `f` is
smallest**. The geodesic hugs the neck: `46%` of the turn is done by
`f = 1.27 f₀`, `76%` by `f = 2.4 f₀`.

> An earlier draft of this page said the geodesic turns "where the lever arm
> `f` is longer and a given angle costs less arc". That is **backwards on both
> counts** — an angular increment at larger `f` costs *more* transverse arc,
> `f dθ` — and it is withdrawn.

Integrating for the `h` that sweeps a given `α`:

> **`T(α) = α²/(2I₂) + O(α⁴)`**,  with  `I₂ = ∫ ds/f²` .

**Leading order, not an exact law** — the exact object is the
Clairaut-integrated `turn_cost`, which comes to `0.9248` of the quadratic value
at `α = π`. Quadratic in the angle, and `I₂` is **not a new quantity**: it is
the resistance `physical_throat.VacuumThroat.resistance` already computes, so
at angular dimension `q = 2`

    T(α) = α² · (4π/I₂) / (8π) + O(α⁴) ,

and the geometric cost of swinging the clock hands is the electrical cost of
pushing monopole flux through the same tube. (That coincidence is specific to
`q = 2` — see below.) With two long arms `I₂ → 4/f₀` and the law reads
`T → f₀α²/8`.

Checked against the *integrated geodesic*, not against the expansion that
predicts it:

| `f₀` | `α` | `T(α)` | `α²/(2I)` | ratio | `f₀α` | `T/(f₀α)` |
|--|--|--|--|--|--|--|
| `1e-03` | `0.02` | `5.00481e-08` | `5.00483e-08` | `0.999997` | `2.00e-05` | `0.00250` |
| `1e-03` | `0.10` | `1.25110e-06` | `1.25121e-06` | `0.999916` | `1.00e-04` | `0.01251` |
| `1e-03` | `0.60` | `4.49084e-05` | `4.50435e-05` | `0.997002` | `6.00e-04` | `0.07485` |
| `1e-03` | `1.00` | `1.24085e-04` | `1.25121e-04` | `0.991722` | `1.00e-03` | `0.12408` |
| `1e-03` | `π` | `1.14201e-03` | `1.23489e-03` | `0.924785` | `3.14e-03` | `0.36351` |

**So the linear guess is an upper bound, and a loose one.** The geodesic spends
`1.25%` of `f₀α` at `α = 0.1` and `36%` at `α = π`. That fraction is a function
of **α alone** — the same to three figures across two decades of `f₀` — which
is itself the signature of the quadratic law.

## Which means the property the point was wanted for survives

The reason for a point in the middle was that the link's cost then did not care
where the mouths were. A finite bearing charges *something*, so the question is
whether it charges enough to matter. It does not:

| `f₀` | `L_o + L_i` | `T(π)` | `T(π)/(L_o+L_i)` |
|--|--|--|--|
| `1e-02` | `1.394519` | `1.150e-02` | `8.25e-03` |
| `1e-03` | `1.356768` | `1.142e-03` | `8.42e-04` |
| `1e-04` | `1.350907` | `1.141e-04` | `8.45e-05` |

And `π` is **the largest separation there is** — `bearing_distance` reduces any
pair of directions to `[0, π]` — so that row is the worst case over the entire
configuration space, not a sample of one.

> **There is no reachable orientation at which the hinge costs as much as the
> journey.**

An earlier draft instead quoted a "break-even angle" of `104` radians, from
`√(2I₂·L)`. That is **withdrawn**: it extrapolates the small-angle quadratic
law far outside its own domain *and* outside the bearing's configuration space,
so it describes nothing. Only the fact that it exceeds `π` carries content.

**The conclusion the point-centre picture was wanted for survives the
regularisation, and survives it more strongly than proposed.** This is a
correction that runs the helpful way.

And the point model is the **limit**, not a rival: `T` is linear in `f₀` at
fixed `α`, and the arms' own excess over their naive value obeys
`L(F) − F ≃ (f₀/2)[ln(4F/f₀) − 1]`. Both vanish with the bearing, so `f₀ → 0`
returns `L_o + L_i`, independent of where the mouths are.

## What "intersection" becomes

Two fronts no longer have to collide at `r = 0`. They land on the bearing at
angular positions, and the question splits into two halves that behave
completely differently:

- **Whether** they meet is a statement about angles — their angular extents
  overlap, or they do not — and `f₀` does not enter it at all.
- **How big** the meeting is depends on `f₀` — *and on the bearing's
  dimension*. In the drawn 2-D cross-section the overlap is a **length**,
  `f₀ × (overlap angle)`. On an `S²` bearing it is an **area**,
  `~f₀² × (angular area)`; on `S³` a **volume**. The first draft quoted only
  the 2-D length. The yes/no criterion is dimension-free; the size law is not.

That is how the point limit is recovered, and it is not the way it is usually
described:

> **As `f₀ → 0` the angular incidence survives and the physical interaction
> region collapses.** Which directions the fronts arrive from, and whether
> their angular extents overlap, are untouched by `f₀`. The region in which
> they actually share space is `f₀ × (overlap angle)` and goes to zero — as
> does the physical separation of two fronts that miss.

So `f₀ → 0` does **not** make every route meet. It shrinks the overlap **and**
the gap together, and the distinction between meeting and missing survives as a
yes/no while disappearing as a length.

## The drawn circle is honest

Two directions span a totally geodesic 2-plane, so the walk between them is a
great-circle arc whatever the bearing's dimension is: `S¹`, `S²`, `S³` give the
same `T`. Checked by computing the cost from an `S²` dot product and again from
the planar angle in the orthonormal frame the pair spans — the same to `0.0e+00`.
So the 2-D cross-section may be drawn without the answer being an artifact of
having drawn a circle.

The one thing that **does** change the answer is the projective identification
`n ∼ −n`, which replaces `α` by `min(α, π−α)`. Below `π/2` it does nothing at
all; near `α = π` it is worth a factor of **418**. An unoriented clock hand has
almost nothing to pay to reverse.

## What is general, and what is only about this neck

The derivation of `T = α²/(2I)` never used the profile. Minimising
`∫½f²θ'²ds` at fixed `∫θ'ds = α` gives `θ' = h/f²`, hence `α = hI` and an
excess length `½h²I = α²/(2I)` — for **any** rotationally symmetric `f(s)` with
a finite minimum, each with its own `I`. The profile enters only through that
one integral.

Asserting that from the algebra would be cheap, so it is checked on a second
profile with nothing to do with the first: `f(s) = √(f₀² + s²)`, not
scalar-flat and not glued to anything. The law holds there to `1.3e-04`, across
six decades of `f₀`, with its own resistance `I = (1/f₀)∫dx/cosh x → π/f₀`
against the scalar-flat neck's `4/f₀`.

**What does not carry over is the correction beyond leading order:**

| | shape `T/(α²/2I)` at `α = π` |
|--|--|
| scalar-flat, `f'² = 1 − f₀/f` | `0.9250` |
| hyperbolic, `f = √(f₀²+s²)` | `0.8886` |

Same law, different `O(α⁴)` term. So the quadratic law is a statement about
**necks**; the `8–11%` deficit at a half turn is a statement about a
**particular** neck.

## The identity underneath: one Dirichlet form

`I₂` turning up in both the monopole conductance and the hinge cost is not a
coincidence of two calculations. They are **one variational problem** on the
tube:

> minimise `E[φ] = ∫ w(s) φ′² ds` at fixed total increment `Δφ`.

The Euler–Lagrange equation is `(wφ′)′ = 0`, so the current `wφ′` is conserved,
`Δφ = c·∫ds/w`, and the minimum is `Δφ²/∫ds/w`. Read it twice:

| | `φ` | the conserved current | what `I₂` gives |
|--|--|--|--|
| static monopole | the potential `u` | the flux `Φ = 4πf²u'` | conductance `4π/I₂` |
| infinitesimal rotation | the azimuth `θ` | Clairaut's `h = f²θ'` | cost `α²/(2I₂)` |

Same equation, different field — but **not automatically the same weight, and
an earlier draft of this page said it was.** The azimuth's weight is the
*metric* coefficient on the angular direction, `w = f²`, for every bearing
dimension. The monopole's weight is the *volume* element, `w = f^q` for an
`S^q` cross-section, so its resistance is `∫ds/f^q`:

| angular dimension | monopole resistance | same as the hinge's `∫ds/f²`? |
|--|--|--|
| `q = 1` (2-D slice, `S¹`) | `∫ds/f` | no |
| **`q = 2`** (3-D slice, `S²`) | **`∫ds/f²`** | **yes** |
| `q = 3` (4-D slice, `S³`) | `∫ds/f³` | no |

They coincide **exactly at `q = 2`** — the case in play, a three-dimensional
spatial throat with `S²` cross-sections, and what `physical_throat` carries. So
`T(α) = α²·(4π/I₂)/(8π)` is an identity where it is used, but it is a fact
about *that dimension* rather than about the construction. Unlike the
great-circle reduction of the hinge, it is **not** dimension-free.

**The sharpest form of it is not about the numbers but about the profiles.**
Normalised each to its own total, the monopole potential and the geodesic's
azimuth are *the same function of position along the tube*:

| `α` | worst deviation of the two normalised profiles | `/α²` |
|--|--|--|
| `1.0` | `4.855e-03` | `4.855e-03` |
| `0.1` | `4.856e-05` | `4.856e-03` |
| `0.01` | `4.856e-07` | `4.856e-03` |

A clean factor of 100 per decade — the deviation is `O(α²)`, and the ratio
column is flat to `1e-3`. **That is why "infinitesimal" is the right word**: at
finite `α` the geodesic starts to feel `I₄` and the static potential never
does.

## `I₂` is universal; `I₄` is where the neck shape enters

Expanding the geodesic exactly rather than to leading order —
`dℓ/ds − 1 = (1 − h²/f²)^{−1/2} − 1` — gives

```
T = ½h²I₂ + ⅜h⁴I₄ + …          α = hI₂ + ½h³I₄ + …          I_n = ∫ds/fⁿ
```

and eliminating `h` leaves

> **`T(α) = α²/(2I₂) − α⁴I₄/(8I₂⁴) + O(α⁶)`** ,
> equivalently `shape = T/(α²/2I₂) = 1 − α²I₄/(4I₂³)` .

The correction is **negative**, which is why every measured shape sits below
one. And the division of labour is exact:

- **What is universal is the leading functional form** `α²/(2I₂)`: every
  rotationally symmetric neck obeys it, each with its own `I₂`. `I₂` itself is
  *not* universal as a number — `4/f₀` here, `π/f₀` on the hyperbolic neck.
  What `I₂` does have is a second job: at `q = 2` it is also the monopole
  resistance, so the leading term is fixed by a quantity the throat already
  carried.
- **`I₄` is the first *additional independent* moment to enter**, and the first
  place the neck's shape shows up. It is **not** the entire profile dependence
  at large angle: `I₆` and beyond enter at `O(α⁶)`, and by `α = π` they matter
  — the quartic prediction misses the measured shape by `7.2e-03` there,
  against `7.9e-09` at `α = 0.1`.

| profile | `I₂` | `I₄` | shape |
|--|--|--|--|
| scalar-flat, `f'² = 1 − f₀/f` | `4/f₀` | `32/(15f₀³)` | `1 − α²/120` |
| hyperbolic, `f = √(f₀²+s²)` | `π/f₀` | `π/(2f₀³)` | `1 − α²/(8π²)` |

`1/120 = 0.00833` against `1/(8π²) = 0.01267` sets **how fast** the two
separate: the shape difference is `α²(1/79 − 1/120) = α² × 4.33e-03`.

| `α` | scalar-flat | hyperbolic | difference |
|--|--|--|--|
| `0.1` | `0.99991667` | `0.99987336` | `4.33e-05` |
| `1.0` | `0.99174547` | `0.98749421` | `4.25e-03` |
| `π` | `0.92496954` | `0.88959734` | `3.54e-02` |

> An earlier draft said the two "agree to eight digits at `α = 0.1`". **That is
> wrong and is withdrawn.** Eight digits is how well *each profile matches its
> own quartic law* (`7.9e-09`); the two profiles match *each other* only to
> `4.3e-05`. Two different quantities, and the smaller was reported for the
> larger.

The moments have closed forms of their own. Substituting `t = f' = √(1 − f₀/f)`
— the profile's own slope — turns all of them into polynomials:

```
I_n(F) = (2/f₀^{n−1}) ∫₀^T (1−t²)^{n−2} dt ,   T = √(1 − f₀/F)
```

so `I₂ = (2/f₀)T` and `I₄ = (2/f₀³)[T − ⅔T³ + ⅕T⁵]`, both checked against
quadrature to `1e-9`. Only even orders appear, because `(1 − h²/f²)^{−1/2}` is
a series in `h²/f²`.

## A formulation trap, recorded

The first draft computed `T(α)` as `route_length − (L_o + L_i)` with the
Clairaut integrals taken in `f`. Both halves of that are wrong, and both were
caught by the same test — that the shape `T/(α²/2I)` must be one number across
decades of `f₀`.

- **The subtraction.** `T` is as small as `1e-15` while the two lengths are
  `O(1)`. That is the failure mode PR #267 had just finished repairing in
  `physical_throat.admittance`, reproduced one round later in a different
  module. Answers were `2.7×` wrong at `f₀ = 1e-07`.
- **The variable.** Integrating in `f` puts a spike of width `√f₀` in front of
  an adaptive quadrature, which walks past it. The direct integral was *worse*
  than the difference at `f₀ = 1e-07` for this reason alone.

The fix is the substitution that made PR #267's admittance a closed form:
`f = f₀cosh²x`, giving `ds = 2f dx` and

    dθ/dx = 2κ/√(cosh⁴x − κ²) ,
    dT/dx = 2f₀cosh²x · κ²/[R(cosh²x + R)] ,   R = √(cosh⁴x − κ²) ,   h = f₀κ ,

with `cosh²x − R` rewritten as `κ²/(cosh²x + R)` so nothing cancels. Both
integrands are smooth and decay like `e^{−2x}`. The result is then stable over
**seven decades** of `f₀`, and the `tanh X = √(1 − f₀/F)` identity — which at
the symmetric point is PR #267's `tanh X = cos a` — reproduces `T = α²/(2I)`
analytically as well as numerically.

## Scope

This is **geometry**, deliberately: a metric, its geodesics, and the transport
of an angular width along them.

- **No field equation is solved on it.** Nothing here is a wave, and nothing
  here evolves — the arms are static.
- **It does not choose between the three candidate bulks.** The finite bearing
  is worked out far enough to be compared against the finite caustic ring and
  the finite neck with moving attachment maps; it is not shown to be the right
  one. What it now has that the others do not is a measured hinge cost.
- **The scalar-flat profile is a worked example, not a claim.** It was chosen
  because every closed form in it is checkable against `physical_throat`. See
  the next section for exactly how much of the result outlives it.
- **The turn cost is not a dynamical statement.** It is a length. What a wave
  does with that length is the next question, and this round does not ask it.
