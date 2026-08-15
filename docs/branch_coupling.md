# Is the mouth transfer part of the field problem, or applied to its free branches?

`geometrodynamics/waves/branch_coupling.py` · renderer
`scripts/geometrodynamics_v55_branch_coupling.py` · probe
`experiments/closure_ledger/branch_coupling_probe.py`

## 0. The gap

PR #254 got the strong result it was after. On the Einstein static universe the
conformally coupled retarded Green function has **exact** image support, so PR
#253's ray branches are the field's branches, carrying the `1/(4π sin χ)` shell
law and a Maslov sign the ray ledger structurally could not.

It got that result with the mouth relation on the **outside**. The relation

```
φ(M⁺, t) = η · φ(M⁻, t + Δ)
```

was applied *to the free branches after they were computed*: enumerate the
arrivals at `M⁺`, shift them by `Δ`, re-emit from `M⁻`, enumerate again. One
traversal, by construction — a post-processing step has no way to notice that
what it re-emits will come back.

> **Does writing the relation into the field problem change the answer, and what
> is the primitive object once it is?**

**Yes**, and the primitive is indexed by a **pair of branches**.

## 1. Scope — which model this is

Stated first and plainly, because the resolvent being exact for a model says
nothing about which model it is. What is solved here is a **self-consistent
rank-one mouth-transfer model**:

* a **linear scalar field** on a **fixed** background (`S³ × R`), `c = 1`;
* the mouth relation carries field **values** through the free Green function.
  It is **not a boundary operator and not a quotient of the manifold.** It has
  **no normal-derivative (flux) matching**; **no reflected channel**, so the
  mouth scattering object is `1×1` where a flux-conserving two-mouth junction
  needs at least a `2×2` unitary; and with `κ < 1` it is **lossy** — power out
  over power in is `κ²` exactly, measured in §8. So it is not literally an
  identification. `shells/junction.py` (PR #249) is what would fix `κ` and
  supply the missing channels, and the flux-conserving boundary operator is the
  next step rather than this one;
* **no backreaction, no stress tensor, no topology change, no rate**;
* the two-throat cross term in §6 is a **throat–throat** interference. It is
  **not** the two-source invariant of roadmap step 3, and is measured here only
  because it has the same bilinear shape and shows which object carries it;
* the winding sums need a regulator. A damping `γ` per unit path length is used
  throughout — the Abel factor that makes the image series converge — and every
  result is either `γ`-independent or reported *as* a `γ`-scaling.

What *is* claimed is narrow and holds: the relation is solved self-consistently
rather than applied to already-computed branches, and that changes the answer.

## 2. The relation, solved rather than applied

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

## 6. Rank counts transfer channels, not histories

| configuration | singular values of `K` (normalised) |
| --- | --- |
| one throat | `1`, `1.3e-16`, `4.6e-17`, … |
| two throats | `1`, `0.542`, `1.0e-16`, … |

Rank is **not** a history count. A single throat with `n = 12` branches per leg
already carries `144` distinct `(a,b)` histories, each with its own arrival time,
and `K` is still rank one. What an outer product counts is independent
**separable transfer channels**, and one value-feedback throat supplies exactly
one.

A second throat adds a second channel, so the rank goes to two. For that sum to
mean anything the two matrices must be in the **same basis**, and they are: both
indices are the topological branch label `(long_way, winding)` in `branch_labels`
order, not the leg length — a label that is shared by legs of different `χ`,
which the two throats certainly have. That is checked in the measurement rather
than left to the coincidence that `χ < 2π−χ < χ+2π < …` holds for every
`χ ∈ (0,π)`.

The interference between the two channels is a **full fringe** — visibility
`−1.000` to `+1.000` across the band — and it is bilinear, one factor from each
throat's pair sum, so it is identically zero without either.

That is the same shape as roadmap step 3's invariant and is **not** that
invariant: these are throats, not sources. What it does show is that the object
carrying a "vanishes when one is removed" quantity is `K`, not either leg.

## 7. Existence, convergence and stability are three different conditions

This is the correction that matters most, and the one the round originally got
wrong. Three statements were being run together:

| | condition | |
| --- | --- | --- |
| **existence** | `L ≠ 1` | `1/(1−L)` is fine for `\|L\| > 1` |
| **convergence** | `\|L\| < 1` | the radius of `Σ Lⁿ`, and it is delay-blind |
| **stability** | `Im ω > 0` for every root of `D(ω) = 1 − L(ω)` | phase-sensitive |

`|T_d|` peaks where `1 − e^{−2πu} ≈ 2πγ`, so the **series radius**
`κ_series = 1/max|T_d|` falls linearly in the regulator — fitted exponent
**`0.9998`**, peak exactly on a bare resonance. That says the traversal
expansion's radius of convergence shrinks as the bare poles approach the real
axis. It does *not* say the coupled system becomes unstable, because `γ` is an
Abel regulator, `T_d` carries the *bare* poles, and a finite coupling **moves**
those poles rather than inflating a gain.

Where they move is computable. With `e^{+iωt}` a root at `ω = ω_r + iω_i` rings
as `e^{iω_r t}e^{−ω_i t}`, so `Im ω > 0` is decay. The uncoupled poles sit at
`ω = m + iγ` and the coupling displaces the `m`-th by

```
δ_m = −ηκ e^{−imΔ} sin(md) / (4π² sin d) ,   Im δ_m ∝ sin(md) sin(mΔ)
```

which **changes sign with the mode** — measured against the actual roots to
`2.2e-04`. Stability is therefore phase-sensitive, and no bound on `|L|` can
decide it. Scanning the delay makes the separation concrete, since `T_d` does not
contain `Δ` at all:

| `Δ` | `κ_series` | `κ_stability` | ratio | first mode to go |
| ---: | ---: | ---: | ---: | ---: |
| `1.0` | `0.7619` | `0.7710` | `1.012` | 11 |
| `0.5` | `0.7619` | `0.7798` | `1.024` | 35 |
| `π` | `0.7619` | **`3.0336`** | **`3.98`** | 11 |

At `Δ = π` every first-order displacement is real (`sin(mπ) = 0`), the poles
slide *along* their line instead of off it, and the threshold rises four-fold. In
the gap, at `κ = 1.520`:

| | |
| --- | ---: |
| round-trip gain at the peak | `1.995` |
| the traversal series | `1.3e+119` — diverges |
| the solve | `0.00450 − 0.13033j` — finite |
| least-damped pole | `Im ω = +0.0145` — still stable |

Solving and summing are not the same operation, and this is the coupling that
shows it.

## 8. What the model leaves out, as numbers

| `κ` | power out / power in | `1×1` scattering magnitude | unitary |
| ---: | ---: | ---: | :---: |
| `0.3` | `0.09` | `0.3` | no |
| `0.6` | `0.36` | `0.6` | no |
| `1.0` | `1.00` | `1.0` | yes |

The ratio is `κ²` exactly. Below `κ = 1` the relation is lossy and cannot be an
identification of anything; at `κ = 1` the *value* relation would be one and the
normal derivative would still be unmatched. The missing pieces, each with what it
would cost:

* **no normal-derivative (flux) matching.** A Darmois–Israel junction matches
  both value and derivative; the mismatch is what `shells/junction.py` (PR #249)
  prices.
* **no reflected channel.** Each mouth has one outgoing amplitude, so the
  scattering object is `1×1`. A flux-conserving two-mouth junction needs at least
  a `2×2` unitary — transmit *and* reflect.
* **no self-adjointness**, and therefore no claim that this is a quotient of the
  manifold.

None of this is a defect of the resolvent, which is exact for the model as posed.
It is a statement of *which model*, and the flux-conserving boundary operator is
the next step.

## 9. What this closes, and what it does not

**Closes:** that the mouth relation can be imposed *inside* the field problem
rather than applied to its output, and that doing so is not a rearrangement — it
adds arrivals at times the free-branch ledger does not contain; that the branch
(image) representation and the mode representation are the same function, poles
and residues included; that the primitive object of a through-throat history is
indexed by a pair of branches, because the amplitude factorizes over that index
and PR #253's closure condition does not; and that existence, convergence and
stability are three separate conditions with three different thresholds.

**Does not close:** everything in §1 and §8. This is a **rank-one mouth-transfer
model**, not a throat boundary operator — no flux matching, no reflected channel,
no unitary two-mouth `S`-matrix, and lossy for `κ < 1`. `κ` is put in by hand and
the background is fixed and unbackreacted. When `Δ + ℓ_c < 0` the loop is
**closed in time** and `1/(1−L)` is a self-consistency condition rather than a
convergent history sum; PR #251's fixed point is where that came from. No stress
tensor, no topology change, no rate, and **no two-source invariant** — roadmap
step 3 is still ahead.

## 10. What this changes for the next step

PR #254 flagged that a two-source invariant `𝒞 = A_A² A_B² (k_A · k_B)²` would
need a branch-resolved definition, because `k = ∇S` is multivalued once the
field arrives on several branches at once.

This round supplies the index it should be written in. The pair `(a, b)` is not
bookkeeping convenience: it is the index on which the *conditions* of the theory
live even though the *amplitudes* factorize over it. A branch-resolved `𝒞` should
therefore be a matrix in the same sense `K` is — and §6 says what to expect from
it, that the rank counts independent histories and the off-diagonal is the part
that vanishes when a source is removed.

Two cautions for that step. Any quantity built from a resummed field inherits
the poles, so a "vanishes when a source is removed" test must be stated at a
fixed, explicitly *stable* coupling — and stable means `Im ω > 0` for every root
of `D`, not `|L| < 1`, which §7 shows can differ by a factor of four. And rank,
if it is used as evidence again, counts separable channels rather than histories;
that is a real distinction and it cost this round a claim.
