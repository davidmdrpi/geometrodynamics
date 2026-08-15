# Does a two-source invariant distinguish a throat from two scatterers?

`geometrodynamics/waves/two_source.py` · renderer
`scripts/geometrodynamics_v58_two_source.py` · probe
`experiments/closure_ledger/two_source_probe.py`

## 0. The question PR #253 left

Rank counting ended with an explicit statement of what it could not supply:

> a quantity that **vanishes** when a source is removed rather than merely
> becoming underdetermined.

Deleting a scalar equation costs a dimension whatever was deleted — a theorem
about square systems, not about photons. The replacement has to be a **field**
quantity, and PRs #254–#257 built the field it has to be written in: an exactly
solvable point throat with a known positive sector, a known chart, and a known
marginal endpoint.

## 1. The object

Superposition is exact for a linear field, so no *linear* functional carries
two-source information — it is additive. A **quadratic** one does, in its cross
term:

```
𝒞  =  Q[φ_A + φ_B] − Q[φ_A] − Q[φ_B]
```

identically zero if either source is switched off. For static sources `Q` is the
interaction energy and `𝒞` is the throat's Green function between the two source
points:

```
𝒞(y_A, y_B) = G(y_A,y_B) + Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B) ,

R = (C − BΓ)⁻¹B          ( = (A − Γ(λ))⁻¹ in PR #257's finite-A chart )
```

Written that way it is exactly the index PR #255 asked for: a **matrix in a pair
of branches**, the branch being which mouth the field entered and which it left.
The diagonal is enter-and-leave the same mouth; the off-diagonal is *through*;
and there is one extra channel that used neither.

| | |
| --- | ---: |
| largest value with a source removed | **`0.0`** |
| smallest value with both present | `0.0366` |
| bilinear in the source strengths | yes |

Zero, not underdetermined. That is the property the ray rounds could not have.

**Every number below is evaluated at** `A = (α₁, α₂, β) = (0.30, 0.35, 0.06)`
**with** `d = 1.3` — strictly inside PR #257's cone, **Löwner margin `0.323`** —
and the exact antipodal endpoint is tested separately in §8 rather than
approached.

## 2. The throat channel is rank two, at any number of sources

The `N × N` table of throat-mediated cross terms is `Vᵀ S V` with `V` of shape
`2 × N`, so its rank is at most two however many sources are placed.

| | |
| --- | ---: |
| sources | 12 |
| rank of the direct table | **12** |
| rank of the throat table | **2** |

The deficiency is the throat's, not the geometry's.

Off the chart the statement needs care, and the careless version is wrong. The
*complex* response obeys `rank R = rank B`, so PR #257's Dirichlet strata are
rank one. But a pair of **real** static sources sees only `Re R`, and the real
part of a rank-one *complex* Hermitian matrix `c·uu†` is
`c(Re u Re uᵀ + Im u Im uᵀ)` — generically **rank two**. So:

| Dirichlet direction | `rank B` | `rank R` | rank of the static table |
| --- | :---: | :---: | :---: |
| complex (time-reversal breaking) | 1 | 1 | **2** |
| real | 1 | 1 | **1** |

A rank-one boundary condition still fills both channels of the static observable
unless it is time-reversal invariant.

## 3. Two things that look like the signature and are not

**The cross term being nonzero.** So is every interference pattern. It is a
statement about there being two sources, not about there being a throat.

**The interaction being anisotropic.** Hold the geodesic separation fixed and
move one source over the sphere of that radius. A free field on this background
*cannot* move — `𝒞` is a function of `χ_AB` alone — and it does not, to
`8e-17`. The throat's varies by **66%** of its mean.

That is a real effect and a real measurement. It is also **not** a wormhole
signature, because two *disconnected* scatterers — `A = diag(α₁, α₂)`, no
mouth-to-mouth term at all — produce a **69%** variation on the same
configuration. Anisotropy detects structure at the mouths, not a connection
between them.

**And the response matrix having off-diagonal entries** is the same trap one
level down: `Γ` is off-diagonal whatever the boundary data says, because two
unconnected scatterers still talk through the ambient field.

## 4. What does discriminate is a parameter count

The static invariant determines exactly **three** numbers — the entries of the
real symmetric `S = Re R`. Two independent scatterers have **two** knobs. So
their image cannot fill the space, and the surface it does fill has an exact
equation. With `p = α₁ − g₀`, `q = α₂ − g₀`,

```
S = [[q, G₀], [G₀, p]] / (pq − G₀²)   ⟹   det S = 1/(pq − G₀²)
```

and therefore, identically in `α`,

```
S₁₂  =  G₀ · det S .
```

The **disconnection defect**

```
𝒲  =  S₁₂ / det S  −  G₀
```

is zero on that surface and nonzero off it.

| | |
| --- | ---: |
| worst `\|𝒲\|` on 200 disconnected draws | **`1.4e-16`** |
| smallest `\|𝒲\|` on connected draws | `5.0e-3` |
| connected draws detected | 199 / 199 |

## 5. And on real `β` the defect **is** the coupling

Not merely nonzero. Working the algebra out for general complex `β`, with
`P = (α₁−g)(α₂−g) − (Re β − G_d)²`,

```
𝒲  =  −Re β  −  (G_d − Re β)(Im β)² / P
```

so on the time-reversal-invariant slice (`β` real — PR #256's `Im β → −Im β`
symmetry)

```
𝒲  =  −β        exactly,  to 5.0e-16 over 120 random draws.
```

Independent of the self-energies, independent of the mouth separation, and —
the part that matters — **independent of the Löwner margin.**

| Löwner margin | `\|𝒞\|` | `𝒲` |
| ---: | ---: | ---: |
| `0.400` | `0.0510` | `−0.060000` |
| `0.100` | `0.0791` | `−0.060000` |
| `0.020` | `0.1503` | `−0.060000` |
| `0.004` | `0.1960` | `−0.060000` |

That is the answer to PR #255's caution that anything built from a resummed
field measures the pole rather than the source. The raw invariant grows `3.8×`
as the cone's boundary is approached; the discriminator drifts `2.1e-17`. Every
row is strictly inside the cone.

## 6. It is a protocol, not a formula in the boundary data

The discriminator has to be buildable by someone who measures interaction
energies, knows the background, and knows where the mouths are — but is not told
`A`. Each observation gives `vᵀ_A S v_B` after the known free term is
subtracted, linear in the three entries of `S`, so three placements determine it
and more overdetermine it.

| | |
| --- | ---: |
| observations | 24 |
| condition number | `9.9` |
| worst entry error in `S` | `1.3e-15` |
| `𝒲` from the observations vs from `A` | **`1.1e-16`** |

## 7. Against the round: a one-frequency test has a blind spot

`𝒲 = 0` has solutions away from `β = 0`. Setting the general expression to zero,

```
(Im β)²  =  −Re β · P / (G_d − Re β)
```

which has a real root on **two** branches — `Re β < 0`, and `Re β > G_d`.

PR #257's gate separates them cleanly. On the invisibility surface
`det(A − Γ) = P · G_d/(G_d − Re β)`, so the upper branch has a **negative
determinant** and is unstable — excluded. The lower branch is not:

| `α` | `Re β` | `Im β` | `\|β\|` | `𝒲` | Löwner margin |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(0.30, 0.35)` | `−0.05` | `0.2390` | `0.244` | `−2.1e-17` | **`+0.0907`** |
| `(0.50, 0.40)` | `−0.02` | `0.2529` | `0.254` | `0.0` | **`+0.2086`** |
| `(0.25, 0.25)` | `−0.10` | `0.1904` | `0.215` | `6.9e-18` | **`+0.0340`** |

These are **connected throats with couplings larger than their own
self-energies**, sitting strictly inside the stable cone, that the static
two-source invariant cannot distinguish from two disconnected scatterers. Not
fine-tuned, not unstable, and not removed by anything the previous round
established.

A single-frequency two-source test therefore **cannot falsify a throat**. It can
only confirm one.

## 8. The antipodal endpoint, on its own

PR #257 ended by showing `d = π` is a different statement rather than a limit:
`Γ(0)` is negative *semi*definite there, with a zero eigenvalue in the symmetric
channel, so `A = 0` sits **on** the cone's boundary. The consequence here is
immediate — `(A − Γ(0))⁻¹` is singular as `A → 0`, so the invariant **diverges**.

| `ε` (`A = εI`) | `𝒞` | `ε·𝒞` | `𝒲` |
| ---: | ---: | ---: | ---: |
| `1e-2` | `0.391` | `0.003913` | `3.5e-18` |
| `1e-3` | `3.428` | `0.003428` | `−3.5e-17` |
| `1e-4` | `33.80` | `0.003380` | `2.5e-16` |
| `1e-5` | `337.5` | `0.003375` | `3.6e-15` |

A clean `1/ε`, and the discriminator does not move: **`𝒲` stays exactly zero
through four decades of divergence**, because it is algebraic in `G₀`, which is
finite at the antipode. The identity survives the endpoint too — a connected
antipodal throat with `β = 0.06` still gives `𝒲 = −0.060000`.

So the loudest available two-source signal carries no information about whether
the mouths are connected. **Size is not evidence** — and this is the cleanest
instance of it in the arc, because the size is unbounded.

## 9. The repair, measured

`Γ` depends on `λ`, so the blind surface moves with it. One frequency gives
three numbers for four parameters; two give six.

| | |
| --- | ---: |
| frequencies | `λ = 0` and `λ = −1` |
| random throats reconstructed | 6 |
| worst parameter error | **`1.2e-15`** |
| worst residual | `1.8e-15` |
| the blind-family member reconstructed | **`3.9e-16`** |

The boundary matrix comes back outright — the blind family included. The only
thing still not observable is the **sign** of `Im β`, and PR #256 already
established that is a time reversal rather than a gap in the measurement.

## 10. What this closes, and what it does not

**Closes:** the gap PR #253 named. There is now a field quantity that is
identically zero without a second source, and its disconnection defect `𝒲` is a
genuine discriminator with an exact value — `−β` on the time-reversal-invariant
slice — that is independent of the self-energies and of the distance to the
stability boundary. It is a protocol, recoverable from measurements without
knowing the boundary data, and it survives the antipodal endpoint.

**Does not close, and this is the sharper half:**

* a **one-frequency** test cannot falsify a throat. The blind family is
  one-parameter, has `|β| ≈ 0.25`, and lies inside the stable cone;
* the fix is two frequencies, which is a stronger experimental demand than "a
  two-source coincidence" and should be quoted as such;
* everything here is **static-source**, so it is an interaction-energy statement,
  not a scattering or radiation statement;
* the boundary data is still four real numbers **chosen**, not derived. `𝒲`
  measures `β`; nothing here says what `β` should be. `shells/junction.py`
  (PR #249) is still what would fix it from matter;
* the throat remains **point-supported**: no interior, no proper length, no
  delay. No backreaction, no stress tensor, no topology change, no rate.

## 11. What this changes for the next step

The stationary-action round (roadmap step 4) now has a target rather than a
description. Three things it inherits:

* **a quantity to vary.** `𝒞` is a cross term of a quadratic functional, which
  is what an on-shell action is made of. "Is the whole history jointly
  stationary" can be asked about *this* object rather than about a picture;
* **a measured non-locality.** `𝒲 = −β` means the arc now has a number
  attached to the one part of the operator with no local realization on `S³`.
  Whatever the action round produces has to reproduce it;
* **a warning with a size.** The blind family is the concrete demonstration that
  a two-source coincidence is not automatically a wormhole detection, and the
  antipodal endpoint is the concrete demonstration that a large effect is not
  automatically a measurement. Both are cheap to re-check and should be.
