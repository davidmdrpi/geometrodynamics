# What does making the throat self-adjoint buy — and what does it not?

`geometrodynamics/waves/throat_operator.py` · renderer
`scripts/geometrodynamics_v56_throat_operator.py` · probe
`experiments/closure_ledger/throat_operator_probe.py`

## 0. What PR #255 owed, and what this round first got wrong

PR #255 solved a mouth relation self-consistently and said plainly what it was
not: a relation between field **values** carried by the free Green function, with
no normal-derivative matching, no reflected channel, a `1×1` mouth object where a
conserving junction needs `2×2`, and `κ²` power throughput.

The upgrade is the right one — but the first version of this round drew a
conclusion from it that is **false**, and the correction is the round's main
content:

> ~~"the spectrum is real for every coupling, so a conserving throat cannot ring
> up"~~ — **withdrawn.** Self-adjointness makes `λ = ω²` real. It does not make
> `λ ≥ 0`, and a negative `λ` is a growing mode.

> **What does self-adjointness buy, then, and where is the family stable?**

Conservation and a real `λ`. Stability is a **separate** condition, and for the
exchange-symmetric family it has a closed form.

## 1. Scope

* a **linear** scalar field on a **fixed** background (`S³ × R`), `c = 1`;
* the throat is **point-supported**: no interior, no proper length, and
  therefore **no delay**. The `Δ` of PRs #251–#255 is not a parameter of a
  self-adjoint point extension and does not survive into one. That is a real
  loss of structure, not a simplification;
* the boundary data is **four real numbers chosen**, not derived.
  `shells/junction.py` (PR #249) is what would fix them from a matter model, and
  nothing here computes the exotic-matter bill;
* **no backreaction, no stress tensor, no topology change, no rate, no
  two-source invariant.**

## 2. The operator, as a *pair* of matrices

Write the field as free plus two point sources,
`φ(x) = φ_in(x) + Σ_j G(χ_j, ω) q_j`, so that near mouth `j` it is
`q_j/(4πχ_j) + φ_j^reg + O(χ_j)` with `φ^reg = φ_in|_mouths + Γ q`. A linear
boundary condition is then a pair:

```
B φ^reg = C q       ⟹       (C − BΓ) q = B φ_in|_mouths
```

so the mouth-active spectrum is `det(C − BΓ) = 0`. The extension is
**self-adjoint** iff `rank[B|C] = 2` and `BC†` is Hermitian. The familiar case
`B = I`, `C = A` recovers `φ^reg = A q` with `A` Hermitian — but the general form
is needed, because **PR #255's relation is not of that shape** (§8).

## 3. It is definable at all

```
G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))
```

real on the real axis, poles exactly at `ω = n+1`, and **finite at the antipode**
— the numerator's zero cancels `sin χ` there, so the antipodal focus is not a
singularity.

| | |
| --- | ---: |
| closed form against PR #255's branch series (`γ → 0`) | **`6.3e-12`** |
| `G(χ) − 1/(4πχ) − g(ω)`, ratio per decade in `χ` | **`10.0`** |
| antipodal limit `ω/(4π sin πω)`, relative | **`< 1e-9`** |

with `g(ω) = −(ω/4π)cot(πω)`. The divergence is the universal Coulomb one, so
the subtraction a point interaction needs is **forced**, not chosen.

## 4. Hermiticity **is** flux conservation

The radial current through a small sphere at mouth `j` is `Im(q_j* φ_j^reg)`,
independent of the sphere's radius, so the total absorbed is `Im(q† A q)` — zero
for **every** `q` iff `A = A†`.

| | |
| --- | ---: |
| worst relative net flux, 200 random Hermitian draws | **`1.8e-16`** |
| pure-transmission throat: one mouth's loss is the other's gain | **`1.7e-16`** |
| a non-Hermitian control, median net flux | **`0.54`** |

**And a correction.** An earlier version read the Cayley entries
`S = (A − ic)(A + ic)⁻¹` as reflection and transmission amplitudes. Unitarity of
`S` is a real fact about the *parametrization* — it is von Neumann's `U(2)` made
concrete, to `2.2e-16` — but the entries' magnitudes depend on the arbitrary
reference scale `c`:

| `c` | `\|S₁₁\|` | `\|S₁₂\|` | column norm |
| ---: | ---: | ---: | ---: |
| `0.05` | `0.955` | `0.296` | `1.000` |
| `0.10` | `0.855` | `0.519` | `1.000` |
| `0.20` | `0.713` | `0.701` | `1.000` |

Same `A`. They are **boundary-mixing coefficients**, not a scattering matrix: a
closed universe has no asymptotic region to normalize flux against. The physical
conservation statement is the identity above.

## 5. Self-adjointness makes `λ` real — and that is all

`Γ(λ)` is real symmetric for real `λ` of **either sign**. For `λ < 0`, with
`ω = iσ`,

```
g = −(σ/4π) coth(πσ) ,      G_d = sinh(σ(π−d)) / (4π sin d · sinh(πσ))
```

both real, both continuous through `λ = 0`. So `det M(λ)` is a real function
(relative imaginary part `2.9e-15`) and its roots are real `λ`. **Nothing forces
them positive**, and a negative `λ` means `ω = ±i√|λ|` with one member of the
pair going like `e^{√|λ| t}`.

Two of the three boundary matrices this module originally advertised do exactly
that:

| `(α₁, α₂, β)` | `σ` | `λ` | growing? |
| --- | ---: | ---: | :---: |
| `(0.05, 0.05, 0.03)` | — | — | no |
| `(0.2, −0.13, 0.15+0.07i)` | **`2.470532`** | `−6.1035` | **yes** |
| `(−0.4, 0.07, −0.09+0.31i)` | **`7.090982`** | `−50.282` | **yes** |

They were missed because the earlier complex-root search seeded only
`Re ω ∈ [1.1, 6.9]` and discarded roots that left that window — a search which
**by construction** could not find a root on the imaginary axis. That is a
scope error in the instrument, not a numerical one, and it is the third time in
this arc that a correct calculation carried a wrong claim.

## 6. What replaces it: the stability region, in closed form

Along the imaginary axis both channel functions are **monotone decreasing** from
their `λ = 0` values to `−∞` (verified on a dense grid), so `g ± G_d = α ± β` has
a negative-`λ` root **iff** the right-hand side falls below the threshold:

```
stable  ⟺  α + β ≥ g₀ + G₀   and   α − β ≥ g₀ − G₀

g₀ = −1/(4π²) = −0.02533030      G₀ = (π−d)/(4π² sin d) = +0.04841232
                                  (d = 1.3)
```

| | |
| --- | ---: |
| symmetric threshold `g₀ + G₀` | **`+0.02308202`** |
| antisymmetric threshold `g₀ − G₀` | **`−0.07374262`** |

Every point of a `13 × 17` grid is also scanned for a negative-`λ` root:

| | |
| --- | ---: |
| grid points | 221 |
| stable | **56** |
| closed-form vs scan mismatches | **0** |

So positivity is a genuine restriction — most of the sampled family is *not*
stable — and it is a condition on the boundary data, not a consequence of
self-adjointness.

## 7. `det(C − BΓ) = 0` is the mouth-active sector, not the spectrum

A two-point interaction is a **rank-two** perturbation. Level `n` on `S³` has
degeneracy `(n+1)²` and only the two combinations with support at the mouths can
move:

| level `n` | degeneracy | mouth-active | untouched |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 1 | 0 |
| 1 | 4 | 2 | 2 |
| 2 | 9 | 2 | 7 |
| 4 | 25 | 2 | **23** |

Twenty-three of twenty-five modes at level 4 never leave the free eigenvalue, and
the secular equation never sees them. Every statement below is scoped to the
sector.

Within it there are three regions, and the first two are invisible to an
`ω`-scan that starts above `1`:

* **`λ < 0`** — the growing modes of §5–6;
* **`0 ≤ λ < 1`** — one mode **below the free ground state**, at `λ = 0.311` for
  the default data;
* **`λ ∈ (m², (m+1)²)`** — two roots per gap, strictly inside.

**And the convenient theorem is false as stated.** "Both channel functions run
`−∞ → +∞` across every gap" fails at the first one: the `n = 0` constant mode has
the same value at both mouths, so it couples only to the symmetric combination
and the antisymmetric channel's pole at `ω = 1` cancels. Its first-gap endpoint
is finite (`−0.0383`), and whether a root exists there depends on `α − β`:

| `α − β` | antisymmetric root in gap 1 |
| ---: | :---: |
| `+0.02` | yes |
| `−0.05` | no |
| `−0.09` | no |

## 8. Where PR #255 sits — and a control that was wrong

That round set `q₁ = 0` (mouth `M⁺` absorbs but does not radiate) and
`q₂ = gain · φ₁^reg`. In the `B φ^reg = C q` form:

```
B = [[0, 0], [gain, 0]] ,   C = I    ⟹    det(C − BΓ) = 1 − gain·G_d
```

which is **exactly** PR #255's own pole condition `1 − L = 0`.

| | |
| --- | ---: |
| embedding against PR #255's own `1 − L` | **`3.5e-18`** |
| `rank[B\|C]` | `2` — maximal |
| `‖BC† − (BC†)†‖` | `1.667` — **not** self-adjoint |

No finite Hermitian `A` reproduces it: it needs the singular `B`.

**The correction.** An earlier version used `A = [[0,0],[1/gain,0]]` with `B = I`
as the control. That gives `det(A − Γ) = g² − G_d² + G_d/gain`, a **different
function** — off by `1.44` at a sample point. The old control was not PR #255's
model, so nothing was concluded from comparing against it, and the conclusion
that was drawn ("#255's instability was its non-conservation") is **withdrawn**:
a self-adjoint throat can have growing modes too, so this is a *classification*
of PR #255's boundary condition, not a diagnosis of its poles.

## 9. What the phase of `β` does

The secular function is `(α₁−g)(α₂−g) − |β − G_d|²` with `G_d` **real**, so it
depends on `Re β` as well as `|β|`: the mouths are joined through the bulk as
well as through the throat, and that fixes the relative phase. So `arg β` is
physical — spread `0.060` across the phase at fixed `|β|`.

But it is invariant under `β → β*`, which is **time reversal** — conjugation
defect exactly `0`. An earlier version called a complex `β` "non-reciprocal",
reading the (scale-dependent) Cayley entries as physical amplitudes. Withdrawn.

## 10. What this closes, and what it does not

**Closes:** that two point removals on `S³` admit a `U(2)` self-adjoint-extension
family; that the regularized Green function and the Krein matrix are well defined
and the point interaction is definable at all; that Hermitian boundary data
conserve the point-boundary flux **exactly**, as an identity; that
`det(C − BΓ) = 0` describes the rank-two mouth-active sector, with the untouched
modes counted; that positivity is a separate condition with a closed form and a
mapped region; and that PR #255's relation embeds exactly as a maximal but
non-self-adjoint boundary condition.

**Does not close:** the boundary data itself; the exotic-matter bill; the
throat's finite size, and with it the delay. And — the one the review forced —
**nothing here explains PR #255's off-axis poles.** A self-adjoint throat can be
unstable too; the two facts are independent.

## 11. What this changes for the next step

The two-source invariant (PR #257) needs a field to be built on, and it now has
one whose stability is a *checkable property of stated boundary data* rather than
an assumption. The practical consequence is a constraint on how that round may be
set up: any invariant evaluated on this background must be evaluated at boundary
data inside the wedge of §6, and the wedge must be quoted, because outside it the
field has a growing mode and every quantity built from it inherits it.

What carries forward from PR #255 is the **branch-pair index**. What does not is
the **delay** — a point extension has no interior — so any two-source statement
that leaned on `Δ` has to be restated in terms of the mouth separation `d` and
the boundary data `(α₁, α₂, β)`.
