# Does a flux-conserving throat behave differently from an effective mouth relation?

`geometrodynamics/waves/throat_operator.py` · renderer
`scripts/geometrodynamics_v56_throat_operator.py` · probe
`experiments/closure_ledger/throat_operator_probe.py`

## 0. What PR #255 owed

That round solved a mouth relation self-consistently and said plainly what it
was not: a relation between field **values** carried by the free Green function,
with no normal-derivative matching, no reflected channel, a `1×1` mouth object
where a flux-conserving junction needs `2×2` unitary, and `κ²` power throughput.
It also found its poles **off** the real axis, and had to separate three
thresholds — existence, convergence, stability — to say what that meant.

> **Does replacing the relation with a genuine flux-conserving boundary operator
> change the answer — and what is the spectrum of the coupled system?**

**Yes**, and the change is the one that mattered: the instability goes away.

## 1. The operator, and why there is no freedom about its form

A point-supported throat is a **self-adjoint extension** of the Laplacian on
`S³ ∖ {M⁺, M⁻}`. Von Neumann's theorem says the extensions are parametrized by
a unitary map between the deficiency subspaces — here `U(2)` — and Krein's
formula makes that concrete as a **Hermitian** `2×2` matrix `A`.

Write the field as free plus two point sources,
`φ(x) = φ_in(x) + Σ_j G(χ_j, ω) q_j`. Near mouth `j` this is
`q_j/(4πχ_j) + φ_j^reg + O(χ_j)`, and the extension is the statement that the
regular parts are a Hermitian image of the charges, `φ^reg = A q`:

```
M(ω) q = φ_in|_mouths ,     M(ω) = A − Γ(ω) ,
Γ(ω) = [[ g(ω), G_d(ω) ], [ G_d(ω), g(ω) ]]
```

`A`'s **diagonal** is each mouth's self-energy — the **reflection** channel PR
#255 had none of. Its **off-diagonal** is the mouth-to-mouth amplitude, the
**transmission** channel, and the only part with no local realization on `S³`:
that non-locality is the wormhole.

## 2. It is definable at all: the Green function has a finite part

Fourier-transforming PR #254's image sum — equivalently, summing PR #255's
branch series in closed form — gives

```
G(χ,ω) = sin(ω(π−χ)) / (4π sin χ · sin(πω))
```

**real** on the real axis, with poles exactly at the free spectrum `ω = n+1`.

| | |
| --- | ---: |
| closed form against PR #255's branch series (`γ → 0`) | **`6.3e-12`** |

Its short-distance expansion is `1/(4πχ) + g(ω) + O(χ)` with

```
g(ω) = −(ω/4π) · cot(πω)
```

and the remainder is **first order in `χ`** — measured convergence ratio `10.0`
per decade. That is what makes a point interaction definable: the divergence is
the universal Coulomb one, so the subtraction is forced rather than chosen.

## 3. The boundary operator is a unitary `2×2`

The Cayley transform `S = (A − ic)(A + ic)⁻¹` is unitary for every Hermitian `A`
and every real reference scale `c`, and inverts back.

| | |
| --- | ---: |
| `‖S†S − I‖` | **`4.4e-16`** |
| `\|r\|² + \|t\|² − 1` at each mouth | **`4.4e-16`** |
| Cayley round-trip back to `A` | **`2.0e-16`** |

Diagonal is reflection, off-diagonal transmission. A **real** `β` is reciprocal
(`t₁₂ = t₂₁`); a complex one is not. And PR #255's model in the same language:

| `κ` | `\|r\|` | `\|t\|` | column norm | in `U(2)` |
| ---: | ---: | ---: | ---: | :---: |
| `0.3` | `0` | `0.3` | `0.09` | no |
| `0.6` | `0` | `0.6` | `0.36` | no |
| `1.0` | `0` | `1.0` | `1.00` | yes |

Outside `U(2)` unless `κ = 1` — and even there it is only half a junction,
because a value relation with no derivative matching is not a boundary
condition.

## 4. Flux conservation **is** Hermiticity

Around mouth `j` the field is `q_j/(4πχ) + φ_j^reg`, and the radial current
integrated over a small sphere is `Im(q_j* φ_j^reg)` — independent of the
sphere's radius. Summing over mouths with `φ^reg = A q`, the total absorbed is

```
Σ_j Im(q_j* φ_j^reg) = Im(q† A q)
```

which is zero for **every** `q` iff `A = A†`. Not on average, not to leading
order — identically.

| | |
| --- | ---: |
| worst relative net flux, 200 random Hermitian draws | **`1.8e-16`** |
| pure-transmission throat: what one mouth absorbs, the other emits | **`1.7e-16`** |
| PR #255's directional relation, median net flux | **`0.54`** |

## 5. And therefore the spectrum is real — for every coupling

`Γ(ω)` is real symmetric on the real axis, so `M = A − Γ` is Hermitian there and
`det M(ω)` is a **real** function of real `ω` (relative imaginary part
`1.5e-15`). Its zeros are the coupled eigenfrequencies.

Newton from a grid of *complex* seeds — deliberately the same method PR #255
used to find its poles, so the comparison is like for like:

| boundary data | roots found | off the axis | worst `\|Im ω\|` | growing |
| --- | ---: | ---: | ---: | ---: |
| Hermitian, `α = (0.05, 0.05)`, `β = 0.03` | 12 | **0** | `6.4e-24` | 0 |
| Hermitian, `α = (0.2, −0.13)`, `β = 0.15+0.07i` | 11 | **0** | `3.9e-18` | 0 |
| Hermitian, `α = (−0.4, 0.07)`, `β = −0.09+0.31i` | 8 | **0** | `4.5e-18` | 0 |
| PR #255 directional, `κ = 0.3` | 9 | **9** | `0.684` | 2 |
| PR #255 directional, `κ = 1.0` | 11 | **11** | `0.357` | 3 |

**PR #255's instability was its own non-conservation.** And it is unstable at
`κ = 1` too, where nothing is lost — so the culprit is the **directionality**,
not the loss. The three thresholds that round had to separate collapse here into
one statement: a conserving throat cannot ring up, at any coupling.

## 6. The coupled spectrum interlaces the free one

For an exchange-symmetric pair (`α₁ = α₂`, `β` real) the secular equation
factorizes into the symmetric and antisymmetric channels:

```
g(ω) + G_d(ω) = α + β        g(ω) − G_d(ω) = α − β
```

Both left-hand sides run **monotonically** from `−∞` to `+∞` across every unit
gap — verified on a dense grid — so each contributes exactly one root there.

| | |
| --- | ---: |
| roots per gap, 8 gaps | **exactly 2** |
| every root strictly between consecutive free frequencies | yes |
| the determinant roots against the channel roots | agree to `1e-8` |

Two coupled frequencies strictly between every pair of consecutive free ones.
Interlacing, not merely shifting — which is what a rank-two perturbation of a
self-adjoint operator does.

## 7. Switching the throat off returns the ESU spectrum

**Off is `‖A‖ → ∞`, not `A → 0`.** The diagonal of `A` is an *inverse*
scattering length, so `α → ∞` is what decouples (`q → 0`) and `α = 0` is a
strongly coupled, resonant throat. Getting this backwards is easy, and the
measurement is written so that it would show.

| `‖A‖` scale | worst distance to the nearest free `ω` | exponent |
| ---: | ---: | ---: |
| `1` | `0.4953` | — |
| `10` | `0.3488` | `0.152` |
| `10²` | `0.0602` | `0.763` |
| `10³` | `0.00621` | `0.987` |
| `10⁴` | `0.000623` | **`0.999`** |

The shift falls like `1/‖A‖`, and the free spectrum `ω = n+1` comes back.

## 8. Where PR #255 sits

Its relation is recovered exactly as the strictly lower-triangular boundary
matrix

```
A(ω) = [[ 0, 0 ], [ 1/(ηκ e^{−iωΔ}), 0 ]]
```

with three defects, each a number:

* **no self-energy**, so no reflection — the diagonal is identically zero;
* **one direction only**, so not Hermitian: the anti-Hermitian part equals the
  whole matrix (`1.667` at `κ = 0.3`);
* **frequency-dependent**, through its phase. A boundary condition is not
  allowed to be: `A` is a constraint on the domain of the operator, not a
  dynamical response.

None of this touches PR #255's resolvent, which was exact for the model it
posed. They are statements about *which model*.

## 9. What this closes, and what it does not

**Closes:** that a two-mouth throat on the ESU can be given a genuine
flux-conserving boundary operator; that flux conservation is exactly Hermiticity
of the boundary matrix, as an identity rather than a limit; that the operator is
a unitary `2×2` with both reflection and transmission; that the coupled spectrum
is real for every coupling, interlaces the free spectrum two per gap, and
returns it on decoupling; and that PR #255's off-axis poles were an artefact of
directionality rather than a property of throats.

**Does not close:** the boundary matrix itself. `A` is four real numbers
**chosen**, not derived — `shells/junction.py` (PR #249) is what would fix them
from a matter model, and nothing here computes the exotic-matter bill. The
throat is **point-supported**, so it has no interior, no proper length, and
therefore **no delay**: the `Δ` that carried PRs #251–#255 is not a parameter of
a self-adjoint point extension and does not survive into one. That is a real
loss of structure relative to those rounds and is stated rather than hidden —
recovering it needs a throat of finite size, which is a different construction.
No backreaction, no stress tensor, no topology change, no rate.

## 10. What this changes for the next step

The two-source invariant (roadmap step 3, PR #257) now has a **self-adjoint**
field to be built on rather than an effective one, which matters for the reason
PR #255 flagged at the end: a quantity built from a resummed field inherits its
poles, and a test stated at "fixed sub-critical gain" is only meaningful if
"sub-critical" is well defined. Here it always is — the spectrum is real, so
there is no critical coupling to stay below, and the invariant can be evaluated
at any `A` without the answer being contaminated by a resonance about to go
unstable.

What carries forward from PR #255 is the **index**: the pair `(a,b)` of branch
labels, on which the theory's conditions live even though its amplitudes
factorize over them. What does not carry forward is the delay, and any
two-source statement that leaned on `Δ` will have to be restated in terms of the
mouth separation and the boundary matrix instead.
