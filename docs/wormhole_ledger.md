# One conserved wave, seen in pieces: an S³ wormhole loop and its ledger

`geometrodynamics/viz/wormhole_ledger.py` · renderer
`scripts/geometrodynamics_v51_wormhole_ledger.py` · probe
`experiments/closure_ledger/wormhole_ledger_probe.py`

## 0. Scope, because everything below depends on it

Kinematics and bookkeeping on a **fixed round `S³`**. No Einstein equation, no
backreaction, no throat construction, and **no claim that such an identification
is dynamically realisable**. The wormhole is an *identification map* — a pair of
mouths and a time offset — and `κ` and the delay `Δ` are **parameters**.

`shells/junction.py` (PR #249) priced what such a throat costs: a minimal
surface has `σ = −(D−2)(β₊+β₋)/8πGR < 0` identically, so a connected throat
needs exotic matter. This round **inherits that bill without paying it**.

Amplitudes are scalar, and the "momentum exchange" at the receiver is a
bookkeeping label on the flux, not a computed cross-section. `c = 1`, `R = 1`.

## 1. The scene

An emitter fires a shell. While it expands, a *collapsing* shell sweeps past it
and lands on a receiver, which recoils. Locally that is two unrelated events —
an emission here, an arriving shell from somewhere else there. Globally it is
one object:

```
E ──expand──▶ M_future ──throat, Δt < 0──▶ M_past ──collapse──▶ R
```

The previous movie in this idiom (`antipodal_crossing`, v10) drew worldlines and
their antipodal traces and had **no notion of this topology**: nothing in it
said the incoming and outgoing shells were the same wave. Here they are drawn in
one colour, and the second panel closes the circuit on the one axis where the
two halves are continuous — **path length against coordinate time**. Path length
only ever increases; the throat is a horizontal jump *backwards in time*.

## 2. The staging is geometry, and it is used twice

A geodesic sphere at distance `χ` on `S³` has area `4π sin²χ`, so an
energy-conserving shell has amplitude `∝ 1/sin χ` and **refocuses exactly at
`χ = π`**. Checked against uniform sampling of `S³` through the enclosed volume
— the area's integral, since a thin-band estimator is dominated by noise exactly
where the area is small — and scored against the estimator's own **binomial
standard error**: worst `z = 2.02`.

That single fact places **both** ends of the loop:

| | placed at | because |
| --- | --- | --- |
| future mouth | the emitter's antipode | where `S³` refocuses the emitted shell |
| receiver | the past mouth's antipode | the only place the arriving shell is *collapsing* |

The second row was a correction. A first draft put the receiver at a generic
point and kept the word "collapse" anyway — but "collapsing" is the sign of
`dA/dχ = 4π sin 2χ`, and at a generic `χ < π/2` the very same wave is still
**expanding** when it lands:

| receiver | arrival `χ` | `dA/dχ` on approach | area at arrival |
| --- | ---: | --- | ---: |
| past mouth's antipode | `π` | negative throughout — collapsing | `1.9e-31` |
| displaced | `1.2` | positive throughout — expanding | `10.92` |

A caption is not a measurement; the sign is.

## 3. The content: one amplitude, and linearity is why

A closed loop through a time-displaced throat invites the obvious objection.
For a **linear** wave there is no paradox and no freedom either:

```
A = A_source + κA    ⟹    A = A_source / (1 − κ)
```

a single fixed point, unique for **every** `κ ≠ 1`. Verified against brute
iteration — 4000 round trips, which solves nothing and just keeps sending the
wave around — agreeing to `7.1e-13` even at `κ = 0.99`, where the amplification
is `×100`. The only obstruction in the whole complex plane is the resonance
`κ = 1`, and the fixed point exists even where `|κ| > 1` and the *iteration*
diverges: divergence of a summation method is not absence of a solution.

**And that is a statement about linear equations, not about time travel.** With
a quadratic return `A = S + κA²` the same closed loop has **two** solutions below
a source threshold of `1/4κ` and **none** above it. Demonstrated rather than
asserted, because the uniqueness would otherwise read as a result about
wormholes.

## 4. The picture measures rather than illustrates

The drawn shell is the geodesic level set — every point sits at distance `χ` from
its centre to `6.7e-15` — and, more usefully, its **screen extent is proportional
to `sin χ` with one constant, to `3.6e-16`**. That constant is `√(A/4π)`, so the
apparent size in the figure *is* the area law rather than a depiction of it.

Getting that required changing the projection, which is this round's real
correction.

**A stereographic chart is unbounded at its own pole, and a shell launched from a
point sweeps all of `S³` — so whatever pole is chosen, the shell crosses it
once.** The radius grows as `2/ε` on approach, verified constant to `4e-6` across
four decades of `ε`, and never converges under refinement:

| `ε` | stereographic radius | `radius × ε` | shadow radius |
| ---: | ---: | ---: | ---: |
| `1e-2` | 200.0 | 1.99998 | 0.563 |
| `1e-3` | 2000.0 | 2.00000 | 0.561 |
| `1e-4` | 20000.0 | 2.00000 | 0.561 |
| `1e-5` | 200000.4 | 2.00000 | 0.561 |

The first renderer projected from `q₃ = +1`, which **is** the emitter's own
position: the emitter was a division by zero and never got drawn at all, and the
shells sprawled off the panel and had to be clipped.

The replacement is the **orthographic shadow** — project `R⁴ → R³` along a fixed
direction and keep the discarded coordinate as brightness. It is bounded by the
unit ball everywhere, and a shell of radius `χ` has diameter `≤ 2 sin χ`, so the
refocus is *visible*: the shell shrinks to a point at the antipode.

The price is that the shadow is **2-to-1** — the two halves of `S³` land on the
same ball, separated only by depth. **A crossing on screen is not a crossing on
`S³`**, and no claim in this round rests on one. The projection axis is also
chosen non-orthogonal to every marked point: if `n ⊥ C` the `χ = π/2` shell about
`C` collapses onto a flat disc, a degenerate artefact rather than a shell.

## 5. What closing the ledger is, and is not, evidence for

Flux conservation through the throat is **put in** when the mouths are
identified. So the ledger residual — `1.1e-16` — checks the arithmetic and
nothing else. It is reported that way because a number that small invites the
opposite reading.

What the ledger *does* establish is the identification of the two halves:
neither local event conserves anything alone — the emitter loses flux, the
receiver gains it — and the pair closes exactly. That balance is the only thing
making them one object rather than two, which is the measurable version of "it is
the same conserved wave".

The delay is a dial. It decides the **entire story** — whether the receiver is
struck before the emitter fires flips between `Δ = −2` and `Δ = −12` — and
changes **no conserved quantity at all**: the amplitude spread across delays from
`−2` to `−10⁴` is `0.0`. The renderer uses `Δ = −5`, chosen so both shells are on
screen together, which is the informative movie and not a preference of the
bookkeeping.

## 6. What this closes, and what it does not

**Closes:** that the four apparent objects can be exhibited as one wave with a
single conserved balance; that the antipodal staging is geometry and is what
makes the arriving shell collapsing; that a linear wave on such a loop has
exactly one amplitude with nothing tuned; and that the figure's drawn sizes are
the area law.

**Does not close:** whether the identification is realisable. No Einstein
equation is solved, `κ` and `Δ` are inputs, the throat's exotic-matter bill is
inherited from PR #249 rather than paid, and the "Compton" label on the
receiver's recoil is bookkeeping, not a cross-section.

## 7. The lesson, which repeats

Both of this round's corrections are **representational**, and neither was
numerical:

* a receiver placed where the caption said "collapse" but the geometry said
  "arrive";
* a projection that sent the very point the scene is about to infinity.

Refining a grid would have caught neither. That is the same shape as the arc's
recurring finding — *a converged number is not a correct number* — with the
sharper version this round: **an object drawn correctly can still be the wrong
object, and only an independent construction says which.**
