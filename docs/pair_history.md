# Two closed histories, sewn at one interaction: is the event selected?

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
> selected by the closure conditions, or still put in by hand?**

That is a counting question about a kinematic skeleton, and it has a definite
answer.

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

which is a **geodesic ellipsoid** — the locus whose summed distance to two foci
is the constant `|Δ|`. Feasible exactly for `d(M⁺,M⁻) ≤ |Δ| ≤ 2π − d(M⁺,M⁻)`,
verified against 40,000 uniform samples of `S³` per configuration rather than
assumed.

So the global system is

```
|c|² = 1                                normalisation
d(S_A, c) = t_C − τ_A                   C lies on front A
d(S_B, c) = t_C − τ_B                   C lies on front B
d(c, M_A⁺) + d(M_A⁻, c) + Δ_A = 0       history A closes
d(c, M_B⁺) + d(M_B⁻, c) + Δ_B = 0       history B closes
```

**Five equations, five unknowns.** The interaction event is not free.

## 3. The event is selected, and existence is restrictive

Solved **blind** from random starts — the question is how many *distinct* events
satisfy the system, so the search cannot be given a hint.

| | |
| --- | ---: |
| events found | 12 |
| at full Jacobian rank 5 | **12** |
| configurations tried | 12 |
| configurations admitting any pair-history | **6** |

Every event found is at **full rank 5**, so isolation is a property of the
system rather than of the solver stopping early. And only about **half** of
random feasible configurations admit a closed pair-history at all — the closure
conditions can forbid a configuration outright, not merely locate an event
inside one.

## 4. Removing a wave removes the *selection*, not the solution

This is the falsification, and it does not land the way the proposal expected.

| | equations | Jacobian rank | solutions |
| --- | ---: | ---: | --- |
| both waves | 5 | **5** | isolated events |
| wave B removed | 4 | **4** | a **one-parameter family**, ~159 sampled points |

The pair-history solution does **not** disappear. There is still a locus of
events closing both histories; there is no longer a *selected* one. Dropping the
constraint can even **create** solutions where two waves admitted none.

So the Breit–Wheeler two-wave requirement appears in this global geometry as
**loss of isolation** — a sharper and weaker statement than nonexistence, and it
is the weaker one that is true. (Separately, with a single front there is no
second independent history to form an opening angle with, so the invariant
`s = 2E₁E₂(1 − cos θ)` is not even defined.)

## 5. And the conjugate pair needs two distinct throats

Not something the two-history picture assumed — it falls out of the rank. A
single shared throat fails in **both** senses of traversal:

* **traversed oppositely.** History B sees `Δ_B = −Δ_A > 0`, so its closure
  demands `b₁ + b₂ = −Δ_B < 0`: a sum of geodesic distances that is negative.
  **Infeasible identically**, with no configuration to search.
* **traversed the same way.** The two closure equations become the *same*
  equation. The rank drops to 4 and the event stops being selected — a family
  again.

## 6. The non-circularity check, which is the one that matters

If the throat delays were unknowns rather than given data, the system would be
**five equations in seven unknowns** with nullity 2, and *every* event lying on
both fronts could be closed by choosing `Δ` afterwards — measured at **100% of
345 sampled events**.

So the entire result rests on the throat being **given**. A version of this
calculation that solved for the delays would select nothing, and would look
identical from the outside. That is exactly where a circular version of this
would hide, so it is measured rather than asserted.

## 7. Closure selects *where*; the invariant decides *whether*

Kept strictly apart, because conflating amplitude with invariant is the error PR
#252 unwound.

| `E/m` | fraction of selected events clearing `s ≥ 4m²` |
| ---: | ---: |
| 1.0 | **0%** |
| 1.5 | 78% |
| 2.0 | 94% |
| 3.0 | 100% |

**Not one** selected event clears threshold at `E = m` — median `2.48`, maximum
`3.97`, just under. The geometry picks the event; the event still has to be paid
for in energy. Two independent conditions, and the module never lets one stand
in for the other.

## 8. What this closes, and what it does not

**Closes:** that requiring two closed histories to share an interaction event is
a *determinate* condition — five equations in five unknowns, isolated
nondegenerate solutions, and roughly half of configurations forbidden outright;
that removing a wave costs isolation rather than existence; that a conjugate
pair cannot ride one shared throat; and that all of this depends on the throat
delays being data rather than unknowns.

**Does not close:** everything the proposal is actually about. No topology
change is attempted, so nothing here shows a two-wave encounter *creating* a
two-mouth sector. No action principle is written, so "the whole history is
jointly stationary" remains a description rather than a computation. No
worldline emerges. And the throats remain identification maps whose
exotic-matter bill is inherited twice over.

## 9. The lesson

The result is a rank, and the interesting part is what happens when a constraint
is removed. *"The solution disappears"* would have been the satisfying answer and
is the wrong one: it survives, and stops being isolated.

Reporting the weaker, truer statement is the discipline this arc keeps
rediscovering — and this time it arrived at the moment the picture was starting
to look like a confirmation.
