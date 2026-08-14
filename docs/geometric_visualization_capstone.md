# The geometric-visualization capstone: what twelve rounds settled (PRs #242–#252, plus the current round)

> **Framing.** Every round below is a *representation* question first and a
> physics question second: given a geometry and a field on it, which object
> should be drawn, and what does the choice smuggle in? The arc began by asking
> whether a wave could fold a surface through itself, passed through locating
> exactly which angular modes let two concentric surfaces couple at all, drew
> four apparently separate objects as one conserved wave, and ended by finding
> that the event it had been drawing all along was not the event at all. Each
> step is machine-checked by a probe in
> `experiments/closure_ledger/`.

## 0. The answer, stated first

**Five negative results and five positive ones, and the negatives are the load-bearing half.**

| | question | answer |
| --- | --- | --- |
| 1 | can a scalar height field on a slice wind through a glued bulk? | **no** — a graph has degree 0 identically |
| 2 | does the gluing scale rescue it? | **no** — it sets shape and magnification, not winding |
| 3 | can the converging ring reach across and intersect? | **no** — it reaches, and never crosses |
| 4 | does drawing the *vectors* rather than their tips change that? | **yes** — hundreds of crossings, threshold `ρ = 1/κ` |
| 5 | does focusing reach a singularity? | **no** — it reaches a *caustic*, and passes through |
| 6 | can a detached oppositely-glued shell replace exotic matter? | **no** — connected implies exotic, by identity |
| 7 | where does the two-shell coupling start? | **`ℓ = 2`** — and it is screened as `(b/a)^ℓ` |
| 8 | are the emitted shell, the passing collapsing shell, the receiver's recoil and a past mouth one wave? | **yes** — one conserved balance, and *linearity* is why it costs nothing |
| 9 | is the antipodal caustic a particle-creation event? | **no** — a caustic is a *venue*; the threshold is on an invariant, needs two waves, and then **forces** a second antipodal interaction |
| 10 | do two closed histories *constrain* the event they share? | **yes, discretely** — five equations, five unknowns, rank 5 on a fixed branch; removing *any* one equation costs a dimension |

The through-line: **almost every apparent obstruction turned out to be a
property of the object being drawn, not of the physics** — until the closing
rounds, where it turned out to be an identity no drawing choice can move, and
finally something sharper still: the arc had been drawing the right object and
measuring the *wrong quantity*.

## 1. The graph cannot wind (`docs/circle_slice_bulk.md`, PR #246)

A great-circle slice through both poles reproduces the sphere to `1.4e-14` and
carries both wave lobes on one circle, meeting at `σ = ±π`. Gluing an inner and
outer shell makes the curve's home a torus, so winding is *available* — and the
wave never uses it. The reason is not dynamical: a graph `r = f(σ)` has winding
number 0 **identically**, whatever the wave does. The seam is crossed in pairs,
up to 346 times, and the degree stays zero.

## 2. The gluing scale is a choice, and it does not rescue anything (`docs/seam_scale.md`, PR #246)

Translate (`r → r − gap`) and conformal (`r → r·R_in/R_out`) gluings are both
legitimate and give different geometry: distortion `1.7027` against `1.0000`,
and only the translate gluing lets the radius reach `−0.30`. The conformal
magnification is start-independent to `2.2e-16`. So the scale decides what
winding *would look like* — and does not decide whether it happens.

## 3. The ring reaches across but never crosses (`docs/ring_and_fold.md`, PR #246)

The converging ring's lead reaches `1.413`, comfortably spanning the gap, and
the curve still never self-intersects. A fold requires *tangential* freedom the
graph does not have, and its threshold `λε ≈ 0.385 w²` is **independent of the
gap** — so shrinking the vacuole could never buy an intersection.

## 4. What a drawn wave has to obey (`docs/wave_constraints.md`, PR #247)

Two constraints the earlier pictures violated. The front is causal, but a
Gaussian source carries a permanent monopole `w²/4` that a closed surface can
never shed — `ℓ = 0` has `ω = 0`. Compactly supported monopole-free data
`(1−x²)⁴` fixes both. Spin-2 is immune: DC content `2e-06` against `8.3e-03` for
the scalar.

## 5. Draw the vectors, not their tips (`docs/normal_field.md`, PR #247)

Four rounds of negative results were all about the same object — the graph of
the displacement's tips, embedded by construction. Draw the displacement
*itself* and the obstruction vanishes for a classical reason: neighbouring
normals meet at the centre of curvature. Same wave, same instant: **0
self-intersections as a graph, 520 as a normal field**, with threshold
`ρ = 1/κ` falling `0.1408 → 0.0540` as the ring converges.

## 6. Focusing reaches a caustic, and passes through (`docs/congruence.md`, PR #248)

A singularity is a failure of evolution, not a bright spot, so the object is a
congruence with a deforming cross-section and the quantity is `J = det F`.
Three things that had been sharing one name come apart: ordinary focus, caustic,
curvature singularity — the last **not reachable in this program**, and that is
a fact about the program.

Of passage, singular termination, and finite-radius reconnection, the equations
give **passage**: `J` crosses zero with slope `−17.877`, converged to `−17.836`
at half the timestep, and evolution continues. And the neck is a **ring** of
radius `0.44 w`, never the antipode — an axisymmetric spin-2 refocus cannot
support a rotationally invariant tensor amplitude at the symmetry pole.

## 7. Connected implies exotic (`docs/shell_junction.md`, PR #249)

Deriving the orientation from the gluing rather than setting a sign by hand
gives **four** gluings, not two, and corrects the framing: `η = ε₊ε₋ = −1`
covers two of them with *opposite* forced signs. What is forced is that a
**minimal surface** has `σ = −(D−2)(β₊+β₋)/8πGR < 0` and a **maximal surface**
the same with the other sign — identities, unviolated in 40,000 random bulks
across `D = 4, 5, 6`.

The dichotomy: a detached surface that **connects** to the throat's region does
so through a minimal surface and is necessarily exotic; one that is non-exotic
by its gluing is a maximal surface, which caps off on both sides and shares no
bulk — non-exotic *precisely because* it is disconnected. Within
Einstein–Israel spherical thin shells, exotic matter is relocated, never
removed.

## 8. Where the coupling starts (`docs/multipole_coupling.md`, PR #250)

In the **static Newtonian two-shell model** the monopole mutual stiffness
vanishes while higher angular multipoles couple, with the coupling suppressed
geometrically by separation:

```
∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ + D − 3) · b^ℓ / a^{ℓ+D−3} · κ_ℓ(D)
```

The prefactor is the eigenvalue of the Laplacian on `S^{D−2}`, so the `ℓ = 0`
decoupling *is* that zero eigenvalue — in every dimension. `D = 4` (`9e-06`) and
`D = 5`, the BAM case (`3.3e-04`, `ℓ = 0` to `1.7e-12`), are each checked by
brute force in their own dimension; at `D = 5` the closed form is
`ℓ(ℓ+2)/(ℓ+1) · b^ℓ/a^{ℓ+2}`. **This is the Newtonian analogue of what
`shell_junction` established in GR — Birkhoff remains a GR theorem, imported by
that round and not replaced here.**

The coupling **starts at `ℓ = 2`**, not `ℓ = 1`. A first draft claimed
"everything `ℓ ≥ 1` couples", from checking translation invariance of the
*area*; run on the mutual *energy*, a rigidly displaced sphere leaves it at
exactly `−G m_b m_a / a` (Newton's shell theorem, `1e-15`), so the translation
mode does not couple. The pure `P₁` *shape* mode is a different object and does,
at `1.78e-02`. And the same formula screens the coupling as `(b/a)^ℓ` — `544×`
from `ℓ = 1` to `ℓ = 8` at `b/a = 0.4`.

## 9. One conserved wave, seen in pieces (`docs/wormhole_ledger.md`, PR #251)

Back to drawing, on `S³`, and the question is whether four apparent objects —
an emitted expanding shell, a passing collapsing one, a receiver's recoil, and a
time-displaced past mouth — can be exhibited as **one wave**. They can, and the
measurable version is that neither local event conserves anything alone while
the pair closes to `1.1e-16`.

The staging is geometry and the **same fact is used twice**: `4π sin²χ` puts the
future mouth at the emitter's antipode *and* the receiver at the past mouth's,
which is the only place `dA/dχ = 4π sin 2χ` is negative all the way in. A first
draft placed the receiver generically and kept the word "collapse"; against a
receiver at `χ = 1.2` the same wave is still expanding when it lands.

The actual content is that a **linear** wave on a closed timelike loop has
exactly one amplitude, `A = A_src/(1 − κ)`, unique for every `κ ≠ 1` — so nothing
is tuned and no paradox is available. That is a fact about linear equations, and
it is fenced as one: a quadratic return gives two solutions or none.

Two things are **put in** and labelled: the throat is an identification map, and
flux conservation through it is an assumption, so the closing ledger checks the
arithmetic. The exotic-matter bill from §7 is inherited, not paid.

## 10. Pair creation is a collision, not a focus (`docs/pair_creation.md`, PR #252)

The arc had been drawing one wavefront refocusing at its antipode and calling
the caustic a creation event. The correction is not that the drawing was
inaccurate — it is that **the quantity was wrong**. A caustic is where the
*amplitude* gets large; Breit–Wheeler is a threshold on an *invariant*,
`s = 2E₁E₂(1 − cos θ) ≥ (2m)²`, and `θ` is not something focusing supplies.

So focusing is **neither sufficient** (collinear momenta have `s = 0` after
amplifying by `10¹²`) **nor necessary** (crossed beams need no focus at all).
The honest complication is kept: a converging front does contain opposed rays,
so the distinction is *independence of the sources*, not brightness.

Two sources `δ` apart give `1 − cos θ = (1 − cos δ)/sin²t` — an identity of
geodesic triangles, so the same on `S²` and `S³`, verified to `2.0e-14` against
embedded tangent vectors. The consequence is the round's result: `s(t)` is
**U-shaped**, head-on at *both* ends of the crossing window and minimal at the
equator, so a threshold cuts **two disjoint windows and never one** — and only
the far one is a collision of waves that have actually propagated
independently, by a factor of `9.6` in path length. **The second interaction
has to be antipodal, and that is derived rather than staged.**

## 11. Two closed histories constrain their shared event (`docs/pair_history.md`, current round)

Sewing two `#251` histories at one interaction is a *determinate* condition, not
a picture. Every leg is null, so a history closes — on the **principal branch** —
on a geodesic ellipsoid with its mouths as foci, and the global system is five
equations in five unknowns. Solved blind, every root found is locally isolated at
**full rank 5**.

**The scope is the result here.** `d` is the principal geodesic distance, and the
prior draws `|Δ|` inside the band where that is the only feasible branch — so the
rest is principal-branch by construction. Off it a mixed branch fixes the
*difference* of distances, a hyperboloid rather than an ellipsoid; discreteness
survives per branch, and what branching changes is the candidate count.

What looked like a falsification is a **dimensionality control**: deleting any
one scalar equation from a square nondegenerate system drops the rank by one, and
deleting a *closure* gives the identical result. It is not evidence about
photons. What survives is the direction — solutions do not vanish, they stop
being isolated.

Two things fall out that were not assumed: in this model a conjugate pair
**cannot ride one shared throat** (infeasible on every branch one way,
rank-deficient the other, scanned rather than argued), and the entire result
rests on the throat delays being **given** — with them free the measured nullity
is 2 and 100% of sampled events close.

## 12. What the arc cost in errors, and what caught them

Worth recording, because the failure modes repeat:

* **A circular measurement.** Normalising each pulse by its own threshold made
  all waves identical *by construction*. Caught by asking what a null result
  would look like.
* **A non-periodic derivative.** `np.gradient` uses one-sided end differences,
  so the closing chord faked a self-intersection where the Jacobian was
  positive.
* **A mismatched grid, twice.** `κ` and `∂²_σu` are second derivatives; refining
  the angular sampling alone gave a `ρ_min` **70% wrong**. Documented after the
  first occurrence and walked into again one round later.
* **An error refinement cannot see.** Seeding `ḧ` from `h(−dt) = h(0)` injects a
  spurious `½ḣ(0)` impulse whose `1/dt²` cancels against the update's `dt`, so
  the ladder *reproduces* the error: `0.628` (no caustic) against `−12.84`.
* **An analytic tail mistaken for a numerical one.** A Gaussian deforms the far
  side at `t = 2.00` against a causal bound of `2.77`, **grid-converged**.
* **A plotting artefact reported as physics.** The causal bound `t = |σ| − w` has
  two branches; drawn as one sorted polyline it fills the panel with a
  triangular hatch that looks exactly like a light cone.
* **A free sign standing in for a derivation.** Carrying `ε = ±1` as a flag made
  a theorem conditional; deriving it from the retained branch corrected the
  headline.
* **A probe scale far outside the physical one.** Extrapolating `V'(p)` from
  `p = 1` when `σ₀ ~ 3e-4` let an `O(h·p²)` term dominate — `p₀` wrong by
  `3.4×`.
* **A zero mode that was not one, and then the wrong test for it.** A pure `P₁`
  deformation is not a translation past linear order; its area cost is `4/3`.
  Worse, the first repair checked translation invariance of the **area** when
  the claim was about the **energy** — and the energy test reverses the
  conclusion, moving the onset of coupling from `ℓ = 1` to `ℓ = 2`.
* **A caption standing in for a sign.** A receiver placed at a generic point was
  described as being *collapsed onto*; `dA/dχ` says the shell arriving there is
  still expanding. Nothing numerical was wrong.
* **A chart that removes the subject.** Stereographic projection is unbounded at
  its own pole, and a shell launched from a point sweeps all of `S³`, so no pole
  is safe. The pole chosen was the emitter's own position — the emitter was a
  division by zero and never appeared in the figure at all.
* **The wrong quantity entirely.** Nine rounds drew a caustic and called it a
  particle-creation event. The drawings were accurate and the object was right;
  amplitude simply is not an invariant, and only asking what the threshold is a
  threshold *on* surfaces that.
* **An angle destroyed by its own projection.** The momenta at a crossing are
  exact to `1e-15`, but measured off the *picture* the opening angle is wrong by
  up to `67°` — and by up to `56°` differently at two crossing points whose true
  angle is identical.
* **An exact identity tested through `arccos`.** `1 − cos θ = 2` at the head-on
  ends is exact; converting it to an angle costs eight digits to the square-root
  singularity at `−1`, turning the statement into a tolerance argument.

The recurring lesson is narrow and practical: **a converged number is not a
correct number.** Three of the errors above survived grid refinement, and were caught only by
an independent construction — a closed form against brute
force, an operator form against a difference, an exact surface against a
truncated family.

The last round sharpens it in a different direction. Neither of its two errors
was numerical, and refining anything would have caught neither: **an object
drawn correctly can still be the wrong object, and only an independent
construction says which.** The repair in both cases was to make the picture
carry a measurement — the sign of `dA/dχ`, and a screen extent proportional to
`sin χ` with one constant to `3.6e-16`, which is `√(A/4π)`.

## 13. What is imported rather than derived

* Birkhoff's theorem (`shell_junction`) — a GR result, still relied on there;
  `multipole_coupling` supplies its static Newtonian analogue, not a
  replacement.
* The Darmois–Israel formalism.
* `β²` as a free parameter at the equilibrium.
* The fixed round `S²` background of every wave round — and the fixed round `S³`
  of the last one: curvature `1` everywhere, at every time, with no Einstein
  equation and no backreaction.
* The wormhole identification itself: a pair of mouths, a time offset `Δ`, a
  loop transfer `κ`, and flux conservation through the throat. All inputs.
* The Breit–Wheeler threshold and cross-section — QED, checked against the
  textbook peak but not derived. Rays-as-photons is a correspondence, and no
  rate is computed anywhere in the arc.
* In the closing round, the throat *data* — mouth positions and delays — which
  is where its entire result lives: with the delays free instead of given,
  every event closes and nothing is selected.

## 14. What would come next

The honest next object is not another drawing. Three of the closing results name
their own missing ingredient:

* a curvature singularity needs **metric evolution with backreaction** — the
  congruence round could not have produced one;
* a finite-radius reconnection needs the congruence to **act back** on
  something, which a test field never does;
* an `ℓ ≥ 2` resonance needs a **shear modulus**, which no equation of state
  supplies and which this round carries as an explicit input.

Each is a construction, not a visualization.

The closing round adds two more, and names them precisely because it stopped
short of them: **topology change** — nothing here shows a two-wave encounter
*creating* a two-mouth sector, only that a shared event is determinate once the
sector is assumed — and an **action principle**, without which "the whole
history is jointly stationary" stays a description rather than a computation.
A worldline would follow from the second, not the first.
