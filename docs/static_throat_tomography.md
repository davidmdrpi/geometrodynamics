# Static two-source throat tomography — not yet the two-wave invariant

`geometrodynamics/waves/two_source.py` · renderer
`scripts/geometrodynamics_v58_two_source.py` · probe
`experiments/closure_ledger/two_source_probe.py`

## 0. What this is not

Roadmap step 3 asked for a **two-wave collision invariant**, the sharp
two-source falsifier the ray rounds could not supply. This round does not
deliver it, and the first draft's framing that it did was wrong.

The object built here is a **static source-interaction kernel** at a fixed
spectral parameter. It carries **no local null momenta**, so it cannot
distinguish equal-energy collinear from counterpropagating waves — which was the
entire load-bearing control behind

```
𝒞 = I_A I_B (k_A · k_B)²
```

The dynamical object, built from `T_A^{μν} T^B_{μν}` and resolved on the
geodesic/winding branches of PRs #253–#255, is still owed. **The roadmap entry
stays open.**

Relatedly: the index pair `(i, j)` used below labels **mouth channels** — which
mouth the field entered and which it left. It is *not* the branch index of
#253–#255, which are geodesic histories with winding numbers and Maslov signs.
The functions are named for mouths.

What *is* delivered is a static inverse result: an exactly solvable tomography
of the two-mouth boundary condition, with one clean identity and a
carefully-bounded negative half.

## 1. The object

PR #253 closed rank counting with an explicit statement of what it could not
supply:

> a quantity that **vanishes** when a source is removed rather than merely
> becoming underdetermined.

Superposition is exact for a linear field, so no *linear* functional carries
two-source information — it is additive. A **quadratic** one does, in its cross
term:

```
𝒞  =  Q[a φ_A + b φ_B] − Q[a φ_A] − Q[b φ_B]  =  2ab · 𝒞(y_A, y_B)
```

For static sources `Q` is the interaction energy and the cross term is the
throat's Green function between the two source points:

```
𝒞(y_A, y_B) = G(y_A,y_B) + Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B) ,

R = (C − BΓ)⁻¹B          ( = (A − Γ(λ))⁻¹ in PR #257's finite-A chart )
```

**Computed the honest way.** A first draft "removed a source" by multiplying the
answer by zero, which proves nothing about the field construction.
`energy_functional` builds `Q` with explicit source strengths and **its own
self-energy terms** — `G^reg_A(y,y) = g₀ + vᵀ_y S v_y`, using the same `g₀` PR
#256 built the boundary condition on — and the cross term comes from three
separate evaluations.

| | |
| --- | ---: |
| `Q[a,b] − Q[a,0] − Q[0,b]` against `ab·𝒞` | **`2.8e-17`** |
| a self-energy term, for scale | `−0.0232` |
| value with source B absent | **`0.0`** |
| smallest cross term with both present | `0.0560` |

The self-energies are present and are not small; they cancel because the
functional is quadratic, not because anything was multiplied by zero.

**Every number below** is evaluated at `A = (α₁, α₂, β) = (0.30, 0.35, 0.06)`
with `d = 1.3` — strictly inside PR #257's cone, **Löwner margin `0.323`** — and
the exact antipodal endpoint is tested separately in §8.

## 2. The throat channel is rank two, at any number of sources

The `N × N` table of throat-mediated cross terms is `Vᵀ S V` with `V` of shape
`2 × N`, so its rank is at most two however many sources are placed.

| | |
| --- | ---: |
| sources | 12 |
| rank of the direct table | **12** |
| rank of the throat table | **2** |

Off the chart the statement needs care. The *complex* response obeys
`rank R = rank B`, so PR #257's Dirichlet strata are rank one. But a pair of
**real** static sources sees only `Re R`, and the real part of a rank-one
*complex* Hermitian matrix `c·uu†` is `c(Re u Re uᵀ + Im u Im uᵀ)` — generically
**rank two**.

| Dirichlet direction | `rank B` | `rank R` | rank of the static table |
| --- | :---: | :---: | :---: |
| complex (time-reversal breaking) | 1 | 1 | **2** |
| real | 1 | 1 | **1** |

## 3. Three things that look like the signature and are not

**The cross term being nonzero.** So is every interference pattern.

**The interaction being anisotropic.** Hold the geodesic separation fixed and
move one source over the sphere of that radius. A free field on this background
*cannot* move — `𝒞` is a function of `χ_AB` alone — and it does not, to
`8e-17`. The throat's varies by **66%** of its mean. Real effect, real
measurement — and two *disconnected* scatterers produce **69%** on the same
configuration. It detects structure at the mouths, not a connection.

**The off-diagonal response block.** `R₁₂ ≠ 0` even for `β = 0`, because `Γ` is
off-diagonal whatever the boundary data says: two unconnected scatterers talk
through the ambient field. So it is a **cross-mouth** channel. Calling it
"through the throat", as the first draft did, contradicts the discriminator
section two paragraphs later.

## 4. What does discriminate is a parameter count

The static invariant determines exactly **three** numbers — the entries of the
real symmetric `S = Re R`. Two independent scatterers have **two** knobs. With
`p = α₁ − g₀`, `q = α₂ − g₀`,

```
S = [[q, G₀], [G₀, p]] / (pq − G₀²)   ⟹   det S = 1/(pq − G₀²)
```

and therefore, identically in `α`,

```
S₁₂  =  G₀ · det S .
```

The **disconnection defect** `𝒲 = S₁₂/det S − G₀` is zero on that surface and
nonzero off it.

| | |
| --- | ---: |
| worst `\|𝒲\|` on 200 disconnected draws | **`1.4e-16`** |
| smallest `\|𝒲\|` on connected draws | `5.0e-3` |

**Scope, stated exactly.** What `𝒲` detects is

> **off-diagonal mouth-boundary mixing, relative to the diagonal two-scatterer
> null model**, inside this point-interaction model.

Not topology, not a traversable interior, not anything the model does not
contain. The first draft's "𝒲 confirms a throat" was broader than what is
proved.

## 5. And on real `β` the defect **is** the coupling

With `P = (α₁−g)(α₂−g) − (Re β − G_d)²`,

```
𝒲  =  −Re β  −  (G_d − Re β)(Im β)² / P
```

so for real `β`

```
𝒲  =  −β        exactly,  to 5.0e-16 over 120 random draws.
```

Independent of the self-energies, the mouth separation, and the Löwner margin.

| Löwner margin | `\|𝒞\|` | `𝒲` |
| ---: | ---: | ---: |
| `0.400` | `0.0510` | `−0.060000` |
| `0.100` | `0.0791` | `−0.060000` |
| `0.020` | `0.1503` | `−0.060000` |
| `0.004` | `0.1960` | `−0.060000` |

That answers PR #255's caution that anything built from a resummed field
measures the pole rather than the source: the invariant grows `3.8×` as the
cone's boundary is approached; the discriminator drifts `2.1e-17`. Every row is
strictly inside the cone.

Real `β` is not a convenient slice — see §7, it is what the field requires.

## 6. It is a protocol, not a formula in the boundary data

Each observation gives `vᵀ_A S v_B` after the known free term is subtracted,
linear in the three entries of `S`, so three placements determine it.

| | |
| --- | ---: |
| observations | 24 |
| condition number | `9.9` |
| worst entry error in `S` | `1.3e-15` |
| `𝒲` from observations vs from `A` | **`1.1e-16`** |

## 7. Which field is being solved — and what that costs the blind family

`𝒲 = 0` has solutions away from `β = 0` **only for complex `β`**. That is not a
free choice.

PR #254 solves a **real** time-domain scalar. A real solution requires the
self-adjoint domain to be invariant under complex conjugation: conjugating
`φ^reg = A q` gives `φ^reg* = A* q*`, which lies in the domain only when
`A* = A`. With Hermiticity that is `A` real symmetric — **`β` real**.

Measured rather than argued:

| `β` | real-field compatible | `max \|Im R\|` | `Im` of the field from a **real** unit source |
| --- | :---: | ---: | ---: |
| `0.06` | yes | `0.0` | **`0.0`** |
| `0.06 + 0.20i` | no | `2.44` | **`−2.4e-3`** |

A complex `β` makes a real static source produce a **complex field**. It is a
complex-scalar, time-reversal-breaking extension, and saying so is the
difference between a limitation of the test and a change of the model.

**So for the arc's field content there is no blind family at all**, and
`𝒲 = −β` settles the question at a single spectral parameter.

## 8. The blind family, scoped

Inside a deliberately complex extension the family is real, and worth stating
precisely because the first draft overstated it. `𝒲 = 0` off `β = 0` has roots
on two branches, and both of the things that remove them are reported.

* **PR #257's gate** removes one. On the invisibility surface
  `det(A − Γ) = P·G_d/(G_d − Re β)`, so `Re β > G_d` has a negative determinant
  and is unstable.
* **Reality of the field** removes the rest, by §7.

| `α` | `Re β` | `Im β` | `\|β\|` | `< min α` | real-field | Löwner margin |
| --- | ---: | ---: | ---: | :---: | :---: | ---: |
| `(0.30, 0.35)` | `−0.05` | `0.2390` | `0.244` | yes | **no** | `+0.0907` |
| `(0.50, 0.40)` | `−0.02` | `0.2529` | `0.254` | yes | **no** | `+0.2086` |
| `(0.25, 0.25)` | `−0.10` | `0.1904` | `0.215` | yes | **no** | `+0.0340` |

The couplings are **comparable to, and smaller than, the self-energies** — the
first draft said larger, which is false for every row.

What survives is a real but narrow statement: *the real-static-source protocol
at a single `λ` is blind on this family, in the complex-scalar model.*

## 9. And even there, the limitation is the protocol

Real static sources see only `Re R` — three numbers for four parameters. That is
where the whole multi-parameter story comes from, and it is a property of the
probe.

A pair of sources with complex strengths `a, b` contributes
`2 Re[a* b · G_A(y_A,y_B)]`, so scanning the relative phase gives both
quadratures, hence the full complex `R`, and then

```
A  =  Γ(λ)  +  R⁻¹
```

at a **single** spectral parameter.

| | |
| --- | ---: |
| placements | 8 |
| condition number | `43` |
| worst error in `R` | `1.8e-14` |
| worst error in `A` | **`3.9e-15`** |

The Hermitian structure matters here: `R₀₁ ≠ R₁₀` for complex `β`, so the
unknowns are `(R₀₀, R₁₁, Re R₀₁, Im R₀₁)` and each complex measurement supplies
two real equations. A symmetric ansatz silently drops `Im R₀₁` — which is
exactly the `Im β` the real protocol could not see.

## 10. The real-static-source reconstruction

Two **spectral parameters**, not "two frequencies": `λ = ω²`, so a negative `λ`
is an *imaginary* frequency. Both defaults are positive and below the free
ground state `λ = 1`, so both are genuinely drivable.

| | |
| --- | ---: |
| spectral parameters | `λ = 0.3` and `λ = 0.8` |
| random throats reconstructed | 6 |
| worst parameter error | **`1.1e-16`** |
| worst residual | `3.6e-15` |
| blind-family member | **`3.9e-16`** |

Solved from several starts: a single start does land in a local minimum on one
draw, with residual `1.5` — caught by the reported residual, which is why it is
reported.

At `λ = 0.8` the invisibility equation has **no real root at all** in this
region, so a second spectral parameter does not merely move the blind set — it
can empty it.

## 11. The antipodal endpoint, on its own

PR #257 showed `d = π` is a different statement rather than a limit: `Γ(0)` is
negative *semi*definite there, so `(A − Γ(0))⁻¹` is singular as `A → 0`.

| `ε` (`A = εI`) | `𝒞` | `ε·𝒞` | `𝒲` |
| ---: | ---: | ---: | ---: |
| `1e-2` | `0.391` | `0.003913` | `3.5e-18` |
| `1e-3` | `3.428` | `0.003428` | `−3.5e-17` |
| `1e-4` | `33.80` | `0.003380` | `2.5e-16` |
| `1e-5` | `337.5` | `0.003375` | `3.6e-15` |

A clean `1/ε`, and the discriminator does not move: `𝒲` stays exactly zero
through four decades of divergence, because it is algebraic in `G₀`, which is
finite at the antipode. `𝒲 = −β` still holds for a connected antipodal throat.

**Size is not evidence** — and this is the cleanest instance in the arc, because
the size is unbounded.

## 12. What this closes, and what it does not

**Closes:** the gap PR #253 named, for *static* probes. There is a field
quantity that is identically zero without a second source, built from a
quadratic functional that carries its own self-energies; its disconnection
defect `𝒲` is exactly `−β` on the real-field sector, independent of the
self-energies and of the distance to the stability boundary; and it is a
protocol, recoverable from measurements without knowing the boundary data.

**Does not close:**

* **the two-wave invariant.** No null momenta, no collinear/head-on control, no
  `T^{μν}`. This is the next PR, not a footnote;
* everything here is **static-source**, so it is an interaction-energy
  statement, not a scattering or radiation statement;
* the boundary data is still four real numbers **chosen**, not derived. `𝒲`
  measures `β`; nothing here says what `β` should be. `shells/junction.py`
  (PR #249) is still what would fix it from matter;
* the throat remains **point-supported**: no interior, no proper length, no
  delay. No backreaction, no stress tensor, no topology change, no rate.

## 13. What this changes for the next step

The next PR is the **dynamical two-wave invariant**, and this round gives it
three things:

* **a control it must pass.** Equal-energy collinear versus counterpropagating.
  A quantity that cannot separate those is not the invariant, whatever else it
  measures — that is the lesson of this round applied to the next one;
* **a number it must reproduce.** `𝒲 = −β` attaches a measured value to the one
  part of the operator with no local realization on `S³`. Whatever the
  dynamical object gives, its static limit has to land here;
* **a habit.** Two false signatures were excluded here by finding a null model
  that reproduces them — anisotropy by two disconnected scatterers, the
  off-diagonal block by `Γ` itself. The two-wave round should name its null
  model before quoting its effect.
