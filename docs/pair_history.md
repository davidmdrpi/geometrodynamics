# Two closed histories, sewn at one interaction: what the closure conditions constrain

`geometrodynamics/waves/pair_history.py` · renderer
`scripts/geometrodynamics_v53_pair_history.py` · probe
`experiments/closure_ledger/pair_history_probe.py`

## 0. The question, and why it is the right first one

PR #251 built a single closed history — an expanding leg, a throat that jumps
backwards in coordinate time, a collapsing leg — and showed the *pair* of local
events closes a ledger neither closes alone. PR #252 established that pair
creation is a threshold on an **invariant**, needs **two** independently
propagated waves, and therefore forces a second interaction at the antipode.

The natural next move is to sew two closed histories at one interaction:

```
γ_A + γ_B  ⟶  H_+ + H_-
```

with each `H` an entire closed history rather than a particle trajectory. The
*tempting* next move — attempting topology change — is not one that can be
checked. This round does the prior question, which can:

> **If two closed histories must share their interaction event, is that event
> constrained by the closure conditions, or still put in by hand?**

That is a counting question about a kinematic skeleton, and it has a definite
answer — with a scope that has to be stated first.

## 1. Scope, stated before the result

* a **counting result on a kinematic skeleton**. Fixed round `S³`, `c = 1`;
* the throats are **identification maps** with **given** mouths and **given**
  delays — not solutions of anything. `shells/junction.py` (PR #249) priced a
  connected throat as necessarily exotic; this round puts **two** on the books
  and pays for neither;
* **no action principle, no field equations, no topology change, no dynamics,
  no rate**, and no demonstration that such a configuration is realisable;
* **nothing here derives a worldline.** Whether a particle trajectory is the
  locus where expanding and collapsing components stay mutually consistent is
  untouched — and calling the apparent particles "dragged" by an incoming
  collapsing front would assume a force law this round does not have;
* the conjugacy `Q_A + Q_B = 0` is a **label**, carried and checked, never
  derived. No charge comes out of geometry here.

## 2. The system

An event `C = (c, t_C)` with `c ∈ S³`. Two waves launched from `S_A, S_B` at
times `τ_A, τ_B`; each is null, so its front reaches `c` when geodesic distance
equals elapsed time. Two throats, each a pair of mouths and a delay `Δ < 0`.

A history runs `C → M⁺`, through the throat, then `M⁻ → C`. Every leg is null,
so it **closes in coordinate time** iff

```
d(c, M⁺) + d(M⁻, c) + Δ = 0
```

which **on the principal branch** is a geodesic ellipsoid — the locus whose
summed distance to two foci is the constant `|Δ|`, feasible exactly for
`d(M⁺,M⁻) ≤ |Δ| ≤ 2π − d(M⁺,M⁻)`, verified against 40,000 uniform samples of
`S³` per configuration rather than assumed.

So the global system is

```
|c|² = 1                                normalisation
d(S_A, c) = t_C − τ_A                   C lies on front A
d(S_B, c) = t_C − τ_B                   C lies on front B
d(c, M_A⁺) + d(M_A⁻, c) + Δ_A = 0       history A closes
d(c, M_B⁺) + d(M_B⁻, c) + Δ_B = 0       history B closes
```

**Five equations, five unknowns.** The interaction event is not free.

## 3. The branch scope, which every other claim depends on

`d` is the **principal** geodesic distance, so the closure equation above
describes only short-way, first-pass legs. On a closed `S³` a null leg may also
take the long way (`2π − d`) or wind (`+2πk`), and those are different
constraints. Three things follow, and they are measured rather than argued:

| | |
| --- | ---: |
| feasible branches with `\|Δ\|` inside `[D, 2π − D]` | **1** (principal only) |
| feasible branches sampling the whole delay axis | up to **4** |
| kinds of closure locus found off the principal branch | **sum *and* difference** |

**The prior hides the issue.** `_random_system` draws `|Δ|` inside the principal
band, where the principal branch is the *only* feasible one — so every other
measurement in this round is principal-branch **by construction of its prior**,
not by argument.

**Off it, the locus changes kind.** A *mixed* branch fixes the **difference** of
the two distances rather than their sum: a hyperboloid, not an ellipsoid. So
"closure is a geodesic ellipsoid" is itself a principal-branch statement.

**Discreteness survives per branch.** On any *fixed* branch pair the system is
still five equations in five unknowns, and the roots found are still locally
isolated at full rank — checked on difference-type branches specifically.

That was the state after review. The round is now **branch-complete**; see the
next section.

## 3b. Branch-completeness is exact, not truncated

The winding cannot run away, and that is what makes completeness achievable
rather than a chosen cutoff. A leg length is `ℓ = (d or 2π−d) + 2πk ≥ 2πk`, and
closure demands `ℓ₁ + ℓ₂ = |Δ|`, so

```
k₁ + k₂  ≤  ⌊ |Δ| / 2π ⌋
```

The feasible branch set is therefore **finite and explicitly bounded by the
delay**. Verified by brute enumeration to winding 11 over 400 random
configurations: **no feasible branch ever exceeds the bound** (worst excess `0`).

| `\|Δ\| / 2π` | winding bound | complete branch count |
| ---: | ---: | ---: |
| 0.5 | 0 | 1 |
| 1.0 | 1 | 2 |
| 2.0 | 2 | 4 |
| 3.0 | 3 | 6 |
| 4.0 | 4 | 8 |

**And the three results survive completion.**

| | branch-complete result |
| --- | --- |
| discrete events | **survive** — 18 of 18 roots at full rank 5 |
| shared-throat obstruction | **survives** — 51 distinct branch pairs, **0** restore full rank |
| delay dependence | **survives** — every point closable by choosing `Δ` |

A union of finitely many discrete sets is discrete, so discreteness is now a
property of the **whole problem** rather than of each branch of it. The candidate
**count grows** — that is the honest cost — and the existence rate shifts; the
local structure does not.

The shared-throat check is now **exhaustive rather than scanned**: with the
winding bounded, every distinct branch pair through the one throat can be tried,
so the earlier "no counterexample at winding ≤ 1" caveat is withdrawn for the
delays reachable here.

## 4. The allowed events are discrete, and existence is restrictive

Solved **blind** from random starts — the question is how many *distinct* events
satisfy the system, so the search cannot be given a hint.

| | |
| --- | ---: |
| events found | 12 |
| at full Jacobian rank 5 | **12** |
| configurations tried | 12 |
| configurations admitting any pair-history | **6** |

Every root found is at **full rank 5**. Said precisely: multi-start root-finding
plus a full-rank Jacobian shows **each root found is locally isolated** — *not*
that all roots were found, and *not* that the event is unique. An earlier draft
of this round called that "the event is selected"; **that was an overstatement
and is withdrawn.** The claim is that the allowed events are discrete on the
sampled branch.

Existence is also restrictive: only about half of random configurations *drawn
from this module's prior* admit a closed pair-history at all.

## 5. Removing a wave is a dimensionality control, not physics

This was presented as a falsification. It is not one, and the correction
matters.

| | equations | Jacobian rank | solutions |
| --- | ---: | ---: | --- |
| both waves | 5 | **5** | isolated events |
| wave B removed | 4 | **4** | a **one-parameter family**, ~159 sampled points |

Deleting wave B removes one scalar equation from a square nondegenerate system,
so rank 5 → 4 and a one-parameter family is exactly what the implicit function
theorem predicts — **for any deleted equation**. The measurement therefore also
deletes a **closure** equation instead, and gets the identical drop:

| deleted equation | rank | solutions |
| --- | ---: | --- |
| front of wave B | 4 | one-parameter family |
| **closure of history A** | 4 | one-parameter family |

So this establishes that the system behaves as a nondegenerate square system. It
is **not** evidence that pair creation needs two photons — that content lives in
the invariant `s`, which needs two independent momenta and is measured
separately.

What survives as interesting is only the *direction*: the solutions do not
vanish, they stop being isolated. Dropping the constraint can even **create**
solutions where two waves admitted none.

## 6. In this model the conjugate pair needs two distinct throats

Not something the two-history picture assumed — it falls out of the rank. A
single shared throat fails in **both** senses of traversal:

* **traversed oppositely.** History B sees `Δ_B = −Δ_A > 0`, so its closure
  demands `ℓ₁ + ℓ₂ < 0`. Leg lengths are non-negative **on every branch**, so
  this is infeasible identically — the one conclusion here that does *not*
  depend on the branch scope.
* **traversed the same way, on the same branch.** The two closure equations
  coincide, the rank drops to 4, and the events stop being discrete.

**The second half is scoped.** It is a statement about the *minimal single-pass
model*. Letting the two histories take **different branches** through the one
throat could in principle give two independent equations, so that is **scanned**
rather than argued: at winding `≤ 1` every distinct branch pair either reduces
to the identical constraint or is jointly inconsistent, and **no counterexample
was found**. Not found is not impossible, and a different gluing is not
excluded.

## 7. The non-circularity check, which is the one that matters

If the throat delays were unknowns rather than given data, the system would be
**five equations in seven unknowns**. The nullity is **measured on the actual
5×7 Jacobian** — rank 5, nullity 2 — rather than counted from the shapes. And
*every* event lying on both fronts can then be closed by choosing `Δ`
afterwards: **100% of 345 sampled events**, with feasibility checked for **both**
throats, not just the first.

So the entire result rests on the throat being **given**. A version of this
calculation that solved for the delays would select nothing, and would look
identical from the outside. That is exactly where a circular version of this
would hide, so it is measured rather than asserted.

## 8. Closure constrains *where*; the invariant decides *whether*

Kept strictly apart, because conflating amplitude with invariant is the error PR
#252 unwound. **But the numbers come with two warnings.**

| `E/m` | fraction of selected events clearing `s ≥ 4m²` |
| ---: | ---: |
| 1.0 | **0%** |
| 1.5 | 78% |
| 2.0 | 94% |
| 3.0 | 100% |

**First: the `0%` at `E = m` is forced, not measured.** `s = 2E²(1 − cos θ) ≤
4E²` always, with equality **only** at exactly head-on — a measure-zero set. So
zero is an inequality, not a finding, and the row is a consistency check.

**Second: every fraction here is prior-dependent.** They are conditioned on
`_random_system`'s arbitrary distribution over mouths, delays and launch times.
They are **regression diagnostics, not predictions**, and nothing should be
inferred from their values.

What remains is the structural point: these are two independent conditions, and
the module never lets one stand in for the other.

## 9. What this closes, and what it does not

**The surviving claim, stated in full:** *with fixed mouth data and delays, and
on a fixed propagation branch, intersecting two null fronts with two independent
closure hypersurfaces generically produces locally isolated candidate events;
removing one front constraint restores a continuous degree of freedom.*

That is useful without promoting the control into pair-production dynamics.

**Also closes:** that a conjugate pair cannot ride one shared throat *in this
model* (with the opposite-traversal half holding on every branch); and that all
of it depends on the throat delays being data rather than unknowns.

**Does not close:** everything the proposal is actually about. No topology
change is attempted, so nothing here shows a two-wave encounter *creating* a
two-mouth sector. No action principle is written, so "the whole history is
jointly stationary" remains a description rather than a computation. No
worldline emerges. And the throats remain identification maps whose
exotic-matter bill is inherited twice over.

## 10. Where this leaves rank counting

At its end. The surviving statement — with given throat data, intersecting two
null fronts with two independent closure hypersurfaces produces locally isolated
candidate events, branch-completely — is as much as constraint counting can
supply.

What it **cannot** supply is a quantity that *vanishes* when a source is removed
rather than merely becoming underdetermined. Deleting any scalar equation costs
a dimension; that is a theorem about square systems, not about photons. A
two-wave discriminator has to be something like

```
𝒞(x) = A_A² A_B² (k_A · k_B)²
```

which is zero without a second source rather than under-determined by its
absence. That is a **field** quantity, and it is the next thing to build.

## 11. The lesson

Every correction this round needed was a **scope** correction, not an arithmetic
one: which branch, which prior, which equation was deleted, and what a rank
actually proves. **The numbers never changed.**

That is a different failure mode from the rest of the arc, which kept catching
wrong objects and wrong quantities. Here the object and the quantity were right
and the *claims attached to them* were too wide — which is the failure mode that
survives having good numbers, and the one a probe full of passing checks is
least able to catch by itself.
