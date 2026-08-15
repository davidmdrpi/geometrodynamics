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

**Why**, in one line. `Γ(λ)` is real symmetric for real `λ` of either sign, and
`dΓ/dλ ≻ 0` below threshold (measured, eigenvalues positive at every sampled
`λ`), so every eigenvalue of `M(λ) = A − Γ(λ)` is **strictly decreasing** in
`λ`. As `λ → −∞`, `Γ → −(σ/4π)I` and both eigenvalues run to `+∞`. So an
eigenvalue crosses zero somewhere below threshold **iff it is already negative
at threshold** — which is the inequality.

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

Step a distance `ε` past the boundary — measured in the smallest eigenvalue of
`A − Γ(0)` — and the eigenvalue appears at `λ ≈ −ε/μ′`:

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

Three facts, each with a consequence:

* **`tr Γ(0) = 2g₀ = −1/(2π²)` at every separation** — the mouth distance does
  not enter it at all;
* its eigenvalues are exactly PR #256's two **channel thresholds**, so that
  round's wedge edges are this round's apex spectrum;
* `det Γ(0) = g₀² − G₀² < 0` everywhere, so `Γ(0)` is **indefinite** — and
  therefore **`A = 0` is unstable at every separation.** No placement of the
  mouths fixes it; the throat needs boundary data that is not zero.

As `d → π` the positive threshold `g₀ + G₀` closes toward zero: antipodal mouths
make the symmetric channel marginally stable at `A = 0`.

## 8. How big the region is

A cone is unbounded, so this only means something with its box stated:

| | |
| --- | ---: |
| box | `\|α_j\|, \|Re β\|, \|Im β\| ≤ 0.2` |
| uniform draws | 4000 |
| stable fraction | **`0.083`** |

Under 9% — the positive sector is a genuine restriction, not a formality.

## 9. What this closes, and what it does not

**Closes:** the question PR #256 left. The positive sector of the two-mouth
self-adjoint family is `A ⪰ Γ(0)`, a forward light cone with apex at the
threshold Krein matrix; its boundary is exactly where a zero mode enters, the
number of growing modes outside is the inertia of `A − Γ(0)`, the instability
turns on continuously as `√ε`, and the previous round's wedge is the
`x₂ = x₃ = 0` slice of it.

**Does not close:** the boundary data itself. `A` is still four real numbers
**chosen**, not derived — `shells/junction.py` (PR #249) is what would fix them
from a matter model, and nothing here computes the exotic-matter bill. Which
point *inside* the cone a physical throat corresponds to is exactly as open as
it was. The throat remains **point-supported**: no interior, no proper length,
no delay. No backreaction, no stress tensor, no topology change, no rate.

## 10. What this changes for the next step

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
  "strictly inside" is the usable region, and `x₀ − |x|` is how far in;
* **`A = 0` is not a neutral choice.** The apex is indefinite at every
  separation, so the tempting "no self-energy, pure transmission" boundary data
  is *unstable* — which is worth knowing before it is used as a baseline.
