# Round 9 — can positive classical history counting reach the quantum correlation?

Pre-registered in `docs/positive_counting_prereg.md` (`a987342`) before the open
questions were computed. Module `geometrodynamics/bulk/positive_counting.py`,
tests (15), probe (9/9).

**Verdicts.**
`Q1: NO_BOUND_ALGEBRAIC_MAXIMUM_REACHED` — `sup CHSH = 4`.
`Q2: QUANTUM_LAW_ATTAINABLE_BY_POSITIVE_COUNTING` — witness `Φ(D) = D²(1 − D/5)`.

The round was opened to test a hypothesis: that every mechanism landing *near
but not on* quantum mechanics (`2.1423`, `3.3941`, `3.7712`) reflected a
structural obstruction, and that a no-go theorem was waiting. **The hypothesis
is false, in both directions.**

## The reduction that makes it computable

On the closure circle, with `t = 1 + u·w`,

```
D(x) = t + √(2t) cos ψ,     ψ uniform in arclength
```

as an equality of laws (sorted residual `2.7e-5`, grid-limited, over 16
`(γ, sector)` cases). Since `t_like = 1 − cos γ` and `t_unlike = 1 + cos γ` sum
to `2` identically, every model in the class collapses to **one function of one
variable**:

```
W(t) = ∫₀^{2π} Φ(t + √(2t) cos ψ) dψ,
E(γ) = [W(τ) − W(2−τ)] / [W(τ) + W(2−τ)],    τ = 1 − cos γ.
```

`Φ = |·|` reproduces round 5 and `Φ = id` reproduces round 6, both to `1e-9`.

Two structural consequences fall out immediately. First, `Γ` is invariant under
`x → −x`, which swaps `(s_A,s_B)` with `(−s_A,−s_B)` while preserving `D`, so
**marginals are exactly `1/2` for every member of the class** — no-signalling is
a property of the class and can never be evidence for a member of it. Second,
`E = −cos γ` for all `γ` is equivalent to the functional equation

```
G(t) = W(t)/t   is even about t = 1.
```

## Q1 — positivity bounds nothing

`W(t) > 0` iff `v ∈ [t − √(2t), t + √(2t)]`, and the upper endpoint increases in
`t`, equalling `v = 1 + √2` at `t = 1`. A narrow nonnegative bump there gives
`W = 0` below `t = 1` and `W > 0` above, hence `E(γ) = −sgn(cos γ)` and

```
CHSH = 4.0000000000
```

at the standard angles, for bump widths `0.2` down to `0.02`. A PR box, with
marginals exactly `1/2`. **Positivity imposes no bound whatever.**

## Q2 — and the quantum law is attained

Expanding `Φ = Σ aₙ Dⁿ`, the monomial transform is exact:

```
Dⁿ  ↦  W(t)/2π = Σⱼ C(n,2j) C(2j,j) 2^{−j} t^{n−j}
```

so `Gₙ(u) = Wₙ/(2πt)` in `u = t − 1` has degree `n−1` with leading coefficient
`1`. The quantum condition is that the odd-power coefficients of `Σ aₙ Gₙ`
vanish. `G₁` is constant (the signed `Φ = D`, which is not nonnegative);
`G₂ = u + 2` and `G₃ = u² + 5u + 4`, so the single condition at degree 3 is

```
a₂ + 5a₃ = 0    ⟹    Φ(D) = D² − D³/5 = D²(1 − D/5),
G(u) = 1.2 − 0.2u²,   even about t = 1.
```

Nonnegativity is immediate rather than numerical: `D² ≥ 0` and `1 − D/5 ≥ 1/5`
on `D ≤ 4`. And the transform has a closed form whose second factor is
*manifestly* symmetric,

```
W(t) = (2π/5) · t · (5 + 2t − t²),     5 + 2t − t²  invariant under  t ↦ 2−t,
```

which **proves** `E = −cos γ` with equal sector priors rather than checking it
angle by angle. Measured, as confirmation:

| `γ` | `E` | `−cos γ` | error |
|---|---:|---:|---:|
| `0.3` | `−0.955336489126` | `−0.955336489126` | `1.1e-16` |
| `1.0` | `−0.540302305868` | `−0.540302305868` | `0.0` |
| `2.0` | `+0.416146836547` | `+0.416146836547` | `1.1e-16` |
| `3.0` | `+0.989992496600` | `+0.989992496600` | `0.0` |

`CHSH = 2.828427124746`, equal to `2√2` to `0.0`.

Degree 3 is minimal: degree 2 forces `a₂ = 0` and collapses to the signed
solution. And **no polynomial solution is globally nonnegative** — for a
solution of top degree `N`, the `u^{N−1}` coefficient is exactly `a_N`, which is
an odd-power condition whenever `N` is even, forcing `a_N = 0`; so the top
degree is always odd. Nonnegativity here is necessarily range-limited, which is
legitimate because the model only evaluates `Φ` on `[−1/2, 4]`.

## What this class contains

| weight `Φ(D)` | nonnegative on `[−½,4]` | CHSH |
|---|---|---:|
| `\|D\|` (round 5, and #283's equilibrium limit) | yes | `2.1422831632` |
| `D` (round 6, signed) | **no** | `2.8284271247` |
| `D²(1 − D/5)` (this round) | **yes** | `2.8284271247` |
| threshold bump at `1+√2` | yes | `4.0000000000` |

## What is withdrawn

Round 5's write-up concluded: *"the distance to quantum mechanics is exactly the
distance from `|D|` to `D`: from counting the two closure branches to summing
them with their holonomy."* **That is wrong and is withdrawn.** Counting
suffices. The quantum correlation does not require signed cancellation, the
closure holonomy, or an oriented current; it requires a particular counting
function, and `|D|` is not it.

This also reframes round 7. That round's central result — that the branch
holonomy and the coarea magnitude come from two provably different functionals,
leaving `κ` free — remains correct as stated, but it is no longer the barrier
between BAM and the quantum law. One can bypass the oriented branch entirely.

## Attaining the quantum law does not close the causality gate

The cubic's sector-summed weight has the closed form

```
Σ_s Φ(D_s(x)) = (8/5)[3 − (a·b)(a·x)(b·x)]      (residual 2.7e-15),
```

which is **constant on the circle for orthogonal settings**. That does not
rescue the source-readout hazard of #284: the closure circles for `b = e_x` and
`b = e_y` lie in *different planes*, so the odd readout `F = (x_x + 2x_y)/√5`
still has variances `0.1` and `0.4`. **Reproducing the singlet correlation and
closing the causality gate are independent requirements**, and this round
settles only the first. This is a sharp-closure statement; #285's finite-spread
result does not transfer to the cubic weighting without specifying its coupled
history extension.

## What remains

Everything now rests on a single question, sharper than before: **what selects
`Φ`?** #283's canonical equilibrium selects `|·|`, hence `2.1423`. Nothing in
the closure geometry, in positivity, in no-signalling, or in the marginals
distinguishes `|·|` from `D²(1 − D/5)`. The near-miss is a property of the
mechanism, not of the model class — so a different mechanism is not merely
permitted, it is unconstrained by anything established so far.

Scope: the class is `W_s = ∫_Γ Φ(D_s) dσ` with `Φ ≥ 0` and a paired sector
prior. It excludes weights that are nonlinear functionals of the integral
(round 7's `(∫D)²`), weights depending on `x` other than through `D`, and
sector-dependent `Φ`. Any of those escapes both conclusions.
