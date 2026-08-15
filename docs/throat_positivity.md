# For which Hermitian `A` is the point-throat operator non-negative?

`geometrodynamics/waves/throat_positivity.py` · renderer
`scripts/geometrodynamics_v57_throat_positivity.py` · probe
`experiments/closure_ledger/throat_positivity_probe.py`

## 0. The question PR #256 left

That round established three things: a point-supported throat on `S³` is a
self-adjoint extension parametrized by a Hermitian `2×2` boundary matrix `A`;
Hermiticity is **exactly** flux conservation; and it does **not** imply
stability — `λ = ω²` is real but need not be positive, and a negative `λ` means
`ω = ±i√|λ|` with a growing mode.

It mapped the stable region on the two-parameter exchange-symmetric slice, by
scanning. Three of the four parameters were left open.

> **For which Hermitian `A` is the operator non-negative?**

The full answer is not a scan. It is one inequality.

## 1. The answer

```
non-negative   ⟺   A ⪰ Γ(0)          (Löwner order)

Γ(0) = [[ g₀, G₀ ], [ G₀, g₀ ]] ,
g₀ = −1/(4π²) = −0.02533030 ,   G₀ = (π−d)/(4π² sin d)
```

for **distinct mouths in the finite-`A` chart** — the scope is not decoration,
and §7 and §10 are what it buys and what it costs.

**Why**, in one line. `Γ(λ)` is real symmetric for real `λ` of either sign, and
`dΓ/dλ ≻ 0` below threshold, so every eigenvalue of `M(λ) = A − Γ(λ)` is
**strictly decreasing** in `λ`. As `λ → −∞`, `Γ → −(σ/4π)I` and both eigenvalues
run to `+∞`. So an eigenvalue crosses zero somewhere below threshold **iff it is
already negative at threshold** — which is the inequality.

That monotonicity is the load-bearing step, and it is a *theorem*, not a sample:
see §9.

Checked, not assumed: for 200 random Hermitian `A` — every one with complex `β`
and unequal mouths, so all four parameters are exercised — the criterion is
compared against an actual negative-`λ` root scan.

| | |
| --- | ---: |
| mismatches | **0 / 200** |
| stable | 19 |
| unstable | 181 |

Both verdicts occur, so the test is not vacuous.

## 2. And the same argument counts

Nothing in the argument is special to `λ* = 0`. For any threshold below the free
ground state `λ = 1`:

```
#{mouth-active eigenvalues < λ*}  =  #{negative eigenvalues of A − Γ(λ*)}
```

a Krein-type **inertia theorem**.

| `λ*` | draws | mismatches |
| ---: | ---: | ---: |
| `−2.0` | 40 | **0** |
| `0.0` | 40 | **0** |
| `0.5` | 40 | **0** |
| `0.9` | 40 | **0** |

Stability is the `λ* = 0` case, and the **count** of growing modes — not just
the yes/no — comes out of it: two negative eigenvalues of `A − Γ(0)` means two
growing modes, verified.

## 3. The geometry: a forward light cone

Hermitian `2×2` matrices are `ℝ⁴` under `A − Γ(0) = x₀I + x·σ`, and positive
semidefiniteness is `x₀ ≥ |x|`. So the stable set is a **forward light cone**:

| | |
| --- | --- |
| **apex** | `A = Γ(0)` — a *doubly* degenerate zero mode |
| **null boundary** `x₀ = \|x\| > 0` | exactly one zero mode, `λ = 0` in the spectrum |
| **interior** `x₀ > \|x\|` | strictly positive, no growing mode |

In coordinates, with `A = [[α₁, β], [β*, α₂]]`:

```
x₀ = (α₁+α₂)/2 − g₀ ,   x₁ = Re β − G₀ ,   x₂ = −Im β ,   x₃ = (α₁−α₂)/2
```

It is a cone in the literal sense, and that is tested: convex, and closed under
positive scaling **from the apex** — not from the origin, which is a different
statement and a false one.

## 4. The boundary is a zero mode, not a convention

On the null surface `A − Γ(0)` is rank one, so `λ = 0` is in the spectrum: a
**static solution supported by the throat**, sitting below the free ground state
`λ = 1`. Its charge pattern is the null vector.

| | |
| --- | ---: |
| secular function at `λ = 0` on the boundary | **`1.8e-17`** |
| marginal mode located by independent root-finding | **`1.4e-14`** |
| strictly growing modes at any boundary point | **0** |
| zero modes at the apex `A = Γ(0)` | **2** |

The second row matters: the zero mode is found by *root-finding on the secular
function*, not read off the eigenvalue that defined the boundary. The apex is
therefore a distinguished point of the parameter space rather than an artefact
of coordinates.

## 5. The instability turns on like a square root

Step past the boundary by a **Löwner margin** `ε = −λ_min(A − Γ(0))`, and the
eigenvalue appears at `λ ≈ −ε/μ′`:

| `ε` | `λ` | `λ/ε` | `σ` | `σ/√ε` |
| ---: | ---: | ---: | ---: | ---: |
| `1e-2` | `−7.820e-2` | `−7.8196` | `0.27964` | `2.7964` |
| `1e-3` | `−7.417e-3` | `−7.4172` | `0.08612` | `2.7235` |
| `1e-4` | `−7.379e-4` | `−7.3787` | `0.02716` | `2.7164` |
| `1e-5` | `−7.375e-5` | `−7.3749` | `0.008588` | `2.7157` |
| `1e-6` | `−7.374e-6` | `−7.3745` | `0.002716` | `2.7156` |

`λ` is **linear** in `ε`, so the growth rate `σ = √|λ|` rises with exponent
**`0.50001`**. And the coefficient is not fitted: `μ′ = d(g+G_d)/dλ` at
threshold is computed independently and predicts `λ/ε = −7.37443`, against the
measured `−7.37448`.

So the boundary is a genuine threshold — the instability is continuous there,
which is what makes "just inside the cone" a meaningful place to work rather
than a knife edge.

## 6. PR #256's wedge is a two-dimensional slice

Setting `α₁ = α₂` and `β` real is exactly `x₂ = x₃ = 0` (verified to `0.0`), and
on that slice the wedge `α ± β ≥ g₀ ± G₀` reproduces the cone at all **143**
sampled points. So the general criterion contains the special one.

Applied to *general* boundary data — by averaging the mouths and dropping
`Im β`, which is the natural way someone would reuse it — the same rule gets
**65 of 400** draws wrong. Two concrete failures:

| `A` | wedge says | truth |
| --- | :---: | :---: |
| `α = 0.05, β = 0.03` | stable | stable |
| `α = 0.05, β = 0.03 + 0.20i` | stable | **unstable** |
| `α₁ = 0.30, α₂ = −0.20, β = 0.03` | stable | **unstable** |

`Im β` and the mouth asymmetry are precisely the two dimensions the wedge cannot
see, and both push `A` out of the cone without changing `α ± Re β`. That is why
the general form was needed rather than a wider scan.

## 7. Where the apex sits

| `d` | `g₀ − G₀` | `g₀ + G₀` | `tr Γ(0)` | indefinite |
| ---: | ---: | ---: | ---: | :---: |
| `0.2` | `−0.40038` | `+0.34972` | `−0.0506606` | yes |
| `0.8` | `−0.10801` | `+0.05735` | `−0.0506606` | yes |
| `1.3` | `−0.07374` | `+0.02308` | `−0.0506606` | yes |
| `2.0` | `−0.05713` | `+0.00647` | `−0.0506606` | yes |
| `3.0` | `−0.05075` | `+0.00008` | `−0.0506606` | yes |
| **`π`** | `−0.05066` | **`0`** | `−0.0506606` | **no** — neg. semidefinite |

Three facts, each with a consequence:

* **`tr Γ(0) = 2g₀ = −1/(2π²)` at every separation** — the mouth distance does
  not enter it at all;
* its eigenvalues are exactly PR #256's two **channel thresholds**, so that
  round's wedge edges are this round's apex spectrum;
* `det Γ(0) = g₀² − G₀² < 0` for `0 < d < π`, so `Γ(0)` is **indefinite** there
  — and therefore **`A = 0` is unstable wherever the mouths are actually
  apart.** No placement short of the antipode fixes it; the throat needs
  boundary data that is not zero.

### The exact antipode is a different statement

`d = π` is not the limit of the rule — it is a case the rule does not cover, and
for `S³` with a through-throat geodesic it is the natural configuration, so it
is worth stating separately. Writing `e = π − d`, the free Green function

```
G(χ, ω) = sin(ω(π−χ)) / (4π sin χ · sin πω)
```

has a **removable** singularity at `χ = π`: numerator and denominator both
vanish linearly in `e`, and the limit is finite,

```
G(π, ω) = ω / (4π sin πω)   ⟹   G₀(d=π) = +1/(4π²) = −g₀ .
```

So at the antipode

```
Γ(0) = g₀ [[1, −1], [−1, 1]] ,   eigenvalues (2g₀, 0) ,
```

**negative semidefinite, not indefinite.** `A = 0` is therefore **marginally
non-negative** at the exact antipode — it sits *on* the cone's boundary, with a
zero mode in the **symmetric** channel — rather than unstable. The earlier code
raised on `χ = π` as though it were a pole; it is not, and `free_green` and
`gamma_at` now take the limit.

Measured: `the_apex_is_indefinite_away_from_the_antipode: True`,
`the_apex_is_negative_semidefinite_at_the_antipode: True`,
`at_the_antipode_A_zero_sits_on_the_boundary: True`,
`the_marginal_channel_is_symmetric: symmetric`.

## 8. How big the region is

A cone is unbounded, so this only means something with its box stated:

| | |
| --- | ---: |
| box | `\|α_j\|, \|Re β\|, \|Im β\| ≤ 0.2` |
| uniform draws | 4000 |
| stable fraction | **`0.083`** |

Under 9% — the positive sector is a genuine restriction, not a formality.

## 9. The monotonicity is a Gram matrix

Everything above rests on `dΓ/dλ ≻ 0` below threshold. That was originally
*sampled* — eigenvalues checked positive at a handful of `λ`. It does not need
to be, and a criterion advertised as exact should not lean on a sample.

Up to a `λ`-independent subtraction (the coincidence-limit renormalization,
which drops out of the derivative),

```
Γ_ij(λ) = ⟨δ_i, (H₀ − λ)⁻¹ δ_j⟩   ⟹   dΓ_ij/dλ = ⟨δ_i, (H₀ − λ)⁻² δ_j⟩ ,
```

which is the **Gram matrix** of the vectors `(H₀ − λ)⁻¹δ_j`. A Gram matrix is
positive semidefinite for free, and positive *definite* exactly when those
vectors are linearly independent — which they are whenever the two mouths are
distinct points. So Löwner monotonicity of `Γ` is an identity, not a numerical
observation, and it holds at every `λ` below threshold and at every separation
including the antipode.

Verified independently: `gram_derivative` builds the sum mode by mode from the
`S³` addition theorem — `Σ_m m²/(2π²(m²−λ)²)` on the diagonal and
`Σ_m m(−1)^{m+1} sin(me)/(2π² sin e · (m²−λ)²)` off it, with the analytic tail
`1/(2π²M)` — and compares it to the closed-form `dΓ/dλ`.

| | |
| --- | ---: |
| worst `\|gram − closed form\|` over 15 `(d, λ)` pairs | **`8.1e-12`** |
| positive definite at every sampled point | **yes** |
| including at `d = π` | **yes** |

The root scans elsewhere in this document are now regression checks on a proved
statement rather than the evidence for it.

## 10. `A ⪰ Γ(0)` is the criterion **in a chart**

`φ^reg = A q` is not the general self-adjoint boundary condition. The general
one is a **pair**,

```
B φ^reg = C q ,   rank[B | C] = 2 ,   B C†  Hermitian ,
```

a `U(2)`'s worth of extensions under `B = i(U−I)`, `C = U+I`. Solving for `A`
needs `B` invertible, so the finite-`A` family is a **chart**, and it misses
exactly the strata where `det B = 0`: **Dirichlet directions**, in which some
combination of charges is forced to zero. Those are reached only as
`‖A‖ → ∞` and are not represented by any finite Hermitian `A` — so "the positive
sector is `A ⪰ Γ(0)`" is a statement about the chart, and saying it about
`U(2)` would be false.

The general criterion drops out of the same argument. Let `P` project onto the
allowed-charge subspace `ker B^⊥`-complement — concretely, `q` is constrained to
a subspace of dimension `k = rank B`, and on it the boundary condition reduces
to an effective Hermitian form `A_eff`. Then

```
non-negative   ⟺   A_eff ⪰ P† Γ(0) P     on the allowed-charge subspace .
```

The chart is the `k = 2` case, `P = I`, `A_eff = A`. The Dirichlet strata are
`k = 1`: one charge direction is frozen, one survives, and the criterion is a
scalar inequality. The free stratum `k = 0` (`q ≡ 0`) has no mouth-active
spectrum at all and is non-negative trivially.

| | |
| --- | ---: |
| chart draws, general form vs the cone | **0 mismatches / 60** |
| `k = 1` stratum draws, general form vs a root scan | **0 mismatches / 60** |
| every stratum carries exactly one Dirichlet direction | yes |
| worst Hermiticity defect of `A_eff` | `1.8e-12` |
| worst row-space defect of the reduction | `9.6e-16` |
| free stratum non-negative | yes |

The last two rows are the reduction checking its own assumptions rather than
asserting them: `A_eff` is verified Hermitian, and `C` is verified to lie in the
row space that the reduction assumes.

## 11. What this closes, and what it does not

**Closes:** the question PR #256 left. For **distinct non-antipodal mouths in
the finite-`A` self-adjoint chart**, the positive sector of the two-mouth family
is exactly `A ⪰ Γ(0)`, a forward light cone with apex at the threshold Krein
matrix; its boundary is exactly where a zero mode enters, the number of growing
modes outside is the inertia of `A − Γ(0)`, the instability turns on
continuously as `√ε`, and the previous round's wedge is the `x₂ = x₃ = 0` slice
of it. The monotonicity behind all of that is a Gram identity (§9), the exact
antipode is the separate marginal case (§7), and off the chart the criterion
reads `A_eff ⪰ P†Γ(0)P` on the allowed-charge subspace (§10).

**Does not close:** the boundary data itself. `A` is still four real numbers
**chosen**, not derived — `shells/junction.py` (PR #249) is what would fix them
from a matter model, and nothing here computes the exotic-matter bill. Which
point *inside* the cone a physical throat corresponds to is exactly as open as
it was. The throat remains **point-supported**: no interior, no proper length,
no delay. No backreaction, no stress tensor, no topology change, no rate.

## 12. What this changes for the next step

The two-source invariant now has a well-posed place to live. Previously "state
the test at stable boundary data" was an instruction with no way to check it on
general `A`; it is now a single matrix inequality, and the failure mode outside
is quantified rather than qualitative — `σ ∝ √ε`, with a count.

Three concrete constraints it hands the next round:

* **quote the point.** Any invariant evaluated on this background must name the
  `A` it used and show `A − Γ(0) ⪰ 0`; the criterion is one line, so there is no
  excuse for leaving it implicit;
* **stay off the boundary.** On the null surface there is a `λ = 0` mode, and a
  static mode contaminates any quantity built by integrating over the field;
  "strictly inside" is the usable region, and the **Löwner (spectral) margin**
  `x₀ − |x| = λ_min(A − Γ(0))` is how far in — quote it with the point;
* **`A = 0` is not a neutral choice** — except at the antipode, where it is
  exactly the marginal one. For `0 < d < π` the apex is indefinite, so the
  tempting "no self-energy, pure transmission" boundary data is *unstable*; at
  `d = π` the same `A = 0` sits on the cone's boundary with a static symmetric
  mode. Both are reasons not to use it as a baseline without saying which case
  you are in, and the antipodal endpoint deserves its own test rather than a
  nearby-`d` stand-in.
