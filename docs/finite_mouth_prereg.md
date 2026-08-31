# Pre-registration: the finite-mouth scalar-flat bulk handle

**Status: frozen before any solver exists.** This document is committed on its
own, ahead of the module it constrains. That ordering is the point — the repo
audit (`docs/repo_audit.md`) found 45 probe runs with 45 passes and 0 failures,
because probe checks were authored after the answer was known. These numbers
are written down first so the next computation has a real chance to fail.

---

## The one new assumption

The observed closed `S³` universe is the totally geodesic equator of a round
four-sphere spatial bulk, `Σ_bulk = S⁴_R`:

```
ds₊² = −dt² + R²(dχ² + sin²χ dΩ₃²)
```

It introduces **no new scale and no dimensionless shape parameter**, and it
makes the existing `S³` exactly a reflection-fixed totally geodesic
hypersurface. Everything below is claimed to follow.

Two geodesic four-balls of angular radius `a < π/2` are excised at `P_A, P_B`;
their `S³` boundaries are joined by `[−S,S] × S³` with
`ds₋² = −N(s)²dt² + ds² + f(s)²dΩ₃²`.

---

## Independently verified before freezing

Confirmed symbolically and to `≤1e-38` numerically (`sympy` + `mpmath`, at
`a = 0.05, 0.3, 0.8, 1.2` and `R = 1, 2.5`):

| step | result |
|--|--|
| `⁴R` of `ds²+f²dΩ₃²` | `−6f″/f + 6(1−f′²)/f²` — matches, exactly |
| `⁴R = 0` | `f = √(s²+b²)` gives `⁴R ≡ 0` |
| Darmois matching | `S = R sin a cos a`, `b = R sin²a` |
| Misner–Sharp | `μ_throat = b²`, `μ_ext(χ) = R²sin⁴χ`, equal at `χ = a` |
| general-lapse stress | both `p_s` and `p_Ω` formulas reproduce exactly |
| exterior ultrastatic `S⁴_R` | `8πG₅ρ = 6/R²`, `8πG₅p = −3/R²` |
| static reduction | `u = y/cosh x`, `s = b sinh x` ⟹ `y″ = (ℓ+1)²y`, residual `0` |
| `tanh X = cos a` at `X = arcosh(1/sin a)` | exact |
| `I₃`, `G`, ANEC closed forms | match quadrature to `1e-38` |

**One correction to the proposal.** The neck NEC theorem does *not* require
`N′(0) = 0`. In

```
8πG₅ p_s = 3 f′N′/(fN) − 3b²/f⁴
```

the `N′` term carries a factor `f′(0) = 0`, which is what *makes* `s=0` a neck.
Substituting a deliberately asymmetric lapse `N = c₀ + c₁s + c₂s²` still returns
exactly `−3/b²`. So the result holds for **every** smooth traversable lapse,
symmetric or not — stronger than stated, and it should be recorded that way.

---

## The six frozen predictions

A finite-mouth construction must reproduce all six **before** it is permitted to
discuss transaction closure. Any failure stops the round.

**P1 — geometry is parameter-free.** Given only `R` and `a`:

```
b = R sin²a ,   S = R sin a cos a ,   L = 2S = R sin 2a ,   f_m = R sin a
```

No independent neck radius, throat length, or tube area exists.

**P2 — Misner–Sharp continuity.** `μ = f²(1−f′²)` equals `b² = R²sin⁴a` inside
the throat and `R²sin⁴χ` outside, so it is continuous across the seam. The seam
is smooth exactly when the quasi-local mass parameter does not jump.

**P3 — no Israel shell.** `[h_ab] = 0` and `[K_ab] = 0`, hence `S_ab = 0`, with
`K^t_t = 0` and `K^A_B = (f′/f)δ^A_B → cos a/(R sin a)` from both sides. The
geometry is `C¹` and not `C²`: `f″` jumps, giving a finite step in bulk stress
and **no** delta-function surface layer. The normal pressure must agree
exactly, `p_s^throat = p_n^ext = −3/(8πG₅R²)`.

**P4 — the neck NEC price.** For **every** smooth lapse with `N(0) > 0`:

```
8πG₅(ρ + p_s)|₀ = −3/b² = −3/(R² sin⁴a)
```

with `ρ ≡ 0` for any lapse (Gauss–Codazzi on a time-symmetric slice with
`⁴R = 0`). Smooth **and** traversable ⟹ radial NEC violation at the neck. The
only escape is `N(0) = 0`, the Tangherlini horizon.

**P5 — the static finite-mouth admittance.** With `X = arcosh(1/sin a)`,
`k = ℓ+1`, `F = R sin a`, and flux `q = 2π²f_m³ n^s∂_sφ`:

```
Y_ℓ(0) = 2π²F² [ k coth(2kX) − cos a      −k csch(2kX)     ]
              [ −k csch(2kX)              k coth(2kX) − cos a ]
```

and for the monopole `Y₀(0) = G·[[1,−1],[−1,1]]` with

```
G = π²R²sin⁴a / cos a = 2π²/I₃ ,   I₃ = 2cos a/(R²sin⁴a)
```

Row sums vanish **exactly**: no static monopole shunt.

**P6 — the vacuum control from the same spatial metric.** The Tangherlini
branch must be obtained by choosing the lapse alone,

```
N_vac(s) = |s|/√(s²+b²)
```

on the *identical* spatial profile, with its horizon at the neck `N(0) = 0`.
The two branches of the repo differ **only** in `N`, not in spatial geometry.

---

## Falsifiers

Stated so they can be checked mechanically:

1. Any independent neck/length/area parameter surviving in the construction
   falsifies **P1**.
2. `|μ_in − μ_out|` above solver tolerance at the seam falsifies **P2**.
3. Any nonzero `[K_ab]`, or a normal pressure mismatch, falsifies **P3**.
4. A smooth `N > 0` giving `(ρ+p_s)|₀ ≥ 0` falsifies **P4** — and would be the
   most interesting outcome in the round, since it would open a traversable
   branch with no local NEC price.
5. A numerically built DtN disagreeing with the closed-form `Y_ℓ(0)`, or
   nonvanishing monopole row sums, falsifies **P5**.
6. A vacuum control needing any spatial change falsifies **P6**.

---

## What must stay separate

Per the audit's dependency-ledger recommendation, the metric says nothing about
the discrete BAM identification. These five must remain distinct objects and
must not be folded into `f(s)`, a sign convention, or an antipodal harmonic
phase:

```
Φ_spatial ,  (−1)^ℓ ,  η_orientation ,  η_wrap ,  U_spin
```

They may be multiplied into a boundary condition only after the boundary lift
is mathematically fixed.

---

## Why this replaces the asymptotic S-matrix

There is no infinity in this geometry, so `T_ℓ(ω)` defined against Jost phases
at `s = ±∞` is not the physical observable. The finite-mouth DtN operator
`Y_ℓ(ω)` at the two actual `S³` mouth surfaces gives the field an absolute
physical reference surface, which removes exactly the constant-phase ambiguity
that dissolved PR #276's closure verdict. Closure then becomes

```
det 𝓜_ℓ(ω) = 0 ,   𝓜_ℓ = Y_ℓ^ext + 𝒰_ℓ† Y_ℓ^th 𝒰_ℓ
```

a **finite-boundary determinant, not an asymptotic phase test**. A change of
mode basis conjugates the operator; its determinant zeros cannot appear or
disappear because a constant was added. That is the process failure #276
exposed, closed structurally rather than by convention.

---

## Scope

This document fixes what the next construction must reproduce. It does **not**
claim the `S⁴_R` completion is correct — that is the one new assumption, and it
is falsifiable only by its consequences. Nor does it address the question the
proposal correctly identifies as the remaining hard one: *what classical BAM
degree of freedom, if any, supplies the stress that keeps `N(0) > 0`?* If there
is none, the geometry collapses onto the Tangherlini horizon branch and the MTY
transaction mechanism is unavailable.
