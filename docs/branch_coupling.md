# Is the throat part of the field problem, or a shift applied to its free branches?

`geometrodynamics/waves/branch_coupling.py` · renderer
`scripts/geometrodynamics_v55_branch_coupling.py` · probe
`experiments/closure_ledger/branch_coupling_probe.py`

## 0. The gap

PR #254 got the strong result it was after. On the Einstein static universe the
conformally coupled retarded Green function has **exact** image support, so PR
#253's ray branches are the field's branches, carrying the `1/(4π sin χ)` shell
law and a Maslov sign the ray ledger structurally could not.

It got that result with the throat on the **outside**. The identification

```
φ(M⁺, t) = η · φ(M⁻, t + Δ)
```

was applied *to the free branches after they were computed*: enumerate the
arrivals at `M⁺`, shift them by `Δ`, re-emit from `M⁻`, enumerate again. One
traversal, by construction — a post-processing step has no way to notice that
what it re-emits will come back.

> **Does writing the identification into the field problem change the answer,
> and what is the primitive object once it is?**

**Yes**, and the primitive is indexed by a **pair of branches**.

## 1. Scope, stated before the result

* a **linear scalar field** on a **fixed** background (`S³ × R`), `c = 1`;
* the throat is still an **identification map** with a transmission `κ` put in
  by hand. `shells/junction.py` (PR #249) priced it and the bill is inherited,
  unpaid. What is new is that it enters an equation that is *solved*;
* **no backreaction, no stress tensor, no topology change, no rate**;
* the two-throat cross term in §6 is a **throat–throat** interference. It is
  **not** the two-source invariant of roadmap step 3, and is measured here only
  because it has the same bilinear shape and shows which object carries it;
* the winding sums need a regulator. A damping `γ` per unit path length is used
  throughout — the Abel factor that makes the image series converge — and every
  result is either `γ`-independent or reported *as* a `γ`-scaling.

## 2. The throat as a boundary condition

In frequency space, the amplitude re-emitted at `M⁻` is driven by everything
that reaches `M⁺`, including its own earlier emission after a trip around the
sphere:

```
a(ω) = ηκ e^{−iωΔ} [ S(ω) + T_d(ω) a(ω) ]

a(ω) = ηκ e^{−iωΔ} S(ω) / (1 − L(ω)) ,     L = ηκ e^{−iωΔ} T_d(ω)
```

`L` is the **round-trip gain**. The resolvent `1/(1−L)` is the solve; PR #254's
answer is its `n = 0` term.

| | |
| --- | ---: |
| resolvent against an explicit walk over 400 traversals | **`3.5e-18`** |
| relative error of the one-traversal answer | exactly `\|L\|` |

That last row is an identity, not a fit: `|1/(1−L) − 1| / |1/(1−L)| = |L|`. The
error of post-processing **is** the round-trip gain.

## 3. The branch series has a closed form, and its poles are the spectrum

The short-way images all carry the Maslov factor `+1` and the long-way images
all carry `−1` — that is what PR #254's sign rule says, and it means the winding
sum is two geometric series:

```
Σ_b s_b e^{−u ℓ_b} = ( e^{−uχ} − e^{−u(2π−χ)} ) / ( 1 − e^{−2πu} ) ,   u = γ + iω
```

| | |
| --- | ---: |
| term-by-term sum against the closed form | **`2.7e-15`** |

As `γ → 0` the denominator vanishes at `ω = 1, 2, 3, …` — the **conformal ESU
eigenfrequencies** `ω_n = n+1` — while the numerator's own zero removes the pole
at `ω = 0` that would otherwise be there. The residues are the mode functions
over `2ω`:

| `ω` | residue | `−(mode weight)/2ω` |
| ---: | ---: | ---: |
| 1 | `−0.0253302` | `−0.0253303` |
| 2 | `−0.0135516` | `−0.0135516` |
| 3 | `+0.0180801` | `+0.0180802` |

So the image representation and the mode representation are **one function**,
recovered from either side. That is the strongest statement so far that the
branch labels are a *representation* rather than an approximation.

## 4. The solve adds events, not amplitudes

The sharpest difference. A history word is `(a, c₁ … c_n, b)` — one branch into
`M⁺`, `n` round trips, one branch out of `M⁻` — arriving at

```
t = ℓ_a + Δ + Σ_i (ℓ_{c_i} + Δ) + ℓ_b
```

with amplitude `κ^{n+1} A₁ A₂ A_dⁿ e^{−γΣℓ}` and sign `η^{n+1}` times **every**
Maslov factor in the word.

| | |
| --- | ---: |
| solved waveform against the sum over history words | **`5.4e-06`** relative |
| solved field at an isolated echo time, over the control | **`3.3e+12`** |
| amplitude ladder, `n = 1, 2, 3, 4` | `1`, `0.0483`, `0.00233`, `0.000113` |
| echo signs against their Maslov words | all agree |

The first row matters: the word enumeration is *checked against the waveform*,
not told as a story about it. The second says the echoes are not small
corrections to arrivals PR #254 already had — at `t = 5.4` (`ℓ_a + ℓ_c + ℓ_b +
2Δ`) the one-traversal field is at the level of numerical zero, and the solved
field is not. **These are events at times that ledger does not contain.**

## 5. The primitive is indexed by a pair of branches

```
K_ab(ω) = ηκ · s_a A₁ e^{−u ℓ_a} · e^{−iωΔ} · s_b A₂ e^{−u ℓ_b}
```

Row `a` is a branch of the leg `source → M⁺`; column `b` a branch of
`M⁻ → observer`. `Σ_ab K_ab / (1 − L)` is the solved propagator.

**Closure is broadband coherence.** `K_ab` carries the phase
`e^{−iω(ℓ_a + Δ + ℓ_b)}`, so PR #253's closure condition `ℓ_a + Δ + ℓ_b = 0` is
*exactly* the statement that `K_ab` does not depend on `ω`: the closed pair
contributes with the same phase at every frequency while every other pair winds
and averages away.

| | |
| --- | ---: |
| band coherence of a closed pair | **`1.000`** |
| band coherence of every other pair | **`< 0.091`** |

**And that is why the pair is the primitive.** The amplitude *factorizes* over
the pair index — `K` is an outer product, rank one — but the condition does not.
At `Δ = −(χ₁ + χ₂ + 4π)` the closed set is `{(k, j) : k + j = 2}` on the
short-way branches: **three** pairs, sitting inside the **nine** that any rule
phrased on `a` alone and `b` alone would have to admit. An anti-diagonal, not a
rectangle. No single-index bookkeeping reproduces it.

## 6. Rank counts histories

| configuration | singular values of `K` (normalised) |
| --- | --- |
| one throat | `1`, `1.3e-16`, `4.6e-17`, … |
| two throats | `1`, `0.542`, `1.0e-16`, … |

One throat is rank one; a second throat adds a second outer product and the rank
goes to two. The interference between them is a **full fringe** — visibility
`−1.000` to `+1.000` across the band — and it is bilinear, one factor from each
throat's pair sum, so it is identically zero without either.

That is the same shape as roadmap step 3's invariant and is **not** that
invariant: these are throats, not sources. What it does show is that the object
carrying a "vanishes when one is removed" quantity is `K`, not either leg.

## 7. Where post-processing stops being an approximation

`|T_d|` peaks where `1 − e^{−2πu} ≈ 2πγ`, so the critical coupling
`κ_c = 1/max|T_d|` falls **linearly** in the regulator:

| `γ` | `κ_c` | ratio |
| ---: | ---: | ---: |
| `0.08` | `3.0465` | — |
| `0.04` | `1.5237` | `1.99944` |
| `0.02` | `0.76189` | `1.99987` |
| `0.01` | `0.38095` | `1.99997` |

Fitted exponent **`0.9998`**, and the peak sits exactly on an ESU
eigenfrequency (`ω = 6` for `d = 1.3`, where `sin(6d)` is near its maximum). As
the regulator is removed, **every** coupling is critical at some frequency, and
there the one-traversal answer is the first term of a divergent series rather
than the leading term of a convergent one.

At fixed `κ = κ_c/2`:

| `ω` | round-trip gain | relative miss of one traversal |
| ---: | ---: | ---: |
| `6.000` (on resonance) | `0.500` | `0.500` |
| `6.500` (off) | `0.0177` | `0.0177` |

## 8. What this closes, and what it does not

**Closes:** that the throat identification can be imposed *inside* the field
problem rather than applied to its output, and that doing so is not a
rearrangement — it adds arrivals at times the free-branch ledger does not
contain; that the branch (image) representation and the mode representation are
the same function, poles and residues included; that the primitive object of a
through-throat history is indexed by a pair of branches, because the amplitude
factorizes over that index and PR #253's closure condition does not; and that
the one-traversal treatment has an expansion parameter with no bound as the
regulator is removed.

**Does not close:** the coupling `κ` is still put in by hand and the background
is still fixed and unbackreacted. When `Δ + ℓ_c < 0` the loop is **closed in
time**, and `1/(1−L)` is then a self-consistency condition rather than a
convergent history sum — it has a unique solution exactly when the
branch-resolved loop gain is subcritical. That is a *bound* on `κ`, and PR #251's
fixed point is where it came from; it is not a derivation of `κ`. No stress
tensor, no topology change, no rate, and **no two-source invariant** — roadmap
step 3 is still ahead.

## 9. What this changes for the next step

PR #254 flagged that a two-source invariant `𝒞 = A_A² A_B² (k_A · k_B)²` would
need a branch-resolved definition, because `k = ∇S` is multivalued once the
field arrives on several branches at once.

This round supplies the index it should be written in. The pair `(a, b)` is not
bookkeeping convenience: it is the index on which the *conditions* of the theory
live even though the *amplitudes* factorize over it. A branch-resolved `𝒞` should
therefore be a matrix in the same sense `K` is — and §6 says what to expect from
it, that the rank counts independent histories and the off-diagonal is the part
that vanishes when a source is removed.

The one caution the `κ_c ∝ γ` scaling raises for that step: any quantity built
from a resummed field inherits the resonances, so a "vanishes when a source is
removed" test must be stated at fixed sub-critical gain, or it will be measuring
the pole rather than the source.
