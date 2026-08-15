# The geometric-visualization capstone: what thirteen rounds settled (PRs #242–#253, plus the current round)

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

**Five negative results and six positive ones, and the negatives are the load-bearing half.**

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
| 10 | do two closed histories *constrain* the event they share? | **yes, discretely** — five equations, five unknowns, rank 5, branch-completely |
| 11 | does a solved *field* reproduce that ray ledger? | **yes, exactly** — the branches are its exact support, and it carries a Maslov phase the rays could not |

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

## 11. Two closed histories constrain their shared event (`docs/pair_history.md`, PR #253)

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

## 12. A solved field reproduces the ledger and signs it (`docs/field_solve.md`)

The move off rays, and it pays immediately. On the Einstein static universe the
**conformally** coupled massless scalar has `ω = n+1` — integer — so its retarded
Green function is exactly periodic and is a sum of images. The geometric-optics
branches are the field's **exact support**, not a stationary-phase
approximation, and a mode sum matches the image sum to `8.3e-13`.

Two earlier rounds fall out of it rather than being modelling choices: the
branch arrival times of `#253` (peaks land on them, grid-limited at `3.0e-04`)
and the shell law `A ∝ 1/sin χ` of `#251` (`peak × sin χ` constant to `7.0e-16`).

**And one thing appears that rays could not have carried.** Every arrival has a
sign, `(−1)^m` with `m` the number of focal crossings — the Maslov index, 12 of
12. Path lengths give times; only a field gives a phase.

A scope fact also surfaces that the ray picture was structurally unable to
notice: the sharp ledger belongs to the **conformally** coupled field. The
minimally coupled one has irrational frequencies and no images — `63%` of the
peak sits between the arrivals against `4.0e-08`. Rays cannot tell the two
apart, because both have the same geodesics.

## 13. The mouth transfer solved for, not applied (`docs/branch_coupling.md`)

The previous round solved the field but kept the mouth relation on the
**outside**: `φ(M⁺,t) = η φ(M⁻,t+Δ)` was applied *to the free branches after
they were computed*, which gives one traversal by construction. Written into the
field problem instead, the amplitude re-emitted at `M⁻` is driven by everything
reaching `M⁺` including its own return, and the solution carries a resolvent
`1/(1−L)` with `L = ηκ e^{−iωΔ} T_d(ω)` the round-trip gain.

Scope first, because the resolvent being exact for a model says nothing about
which model: this is a **rank-one mouth-transfer model**, not a throat boundary
operator and not a quotient. No flux matching, no reflected channel, a `1×1`
mouth scattering object where a conserving junction needs `2×2` unitary, and
power out over power in is `κ²` — lossy, so not an identification.

**Within that, it is not a rearrangement.** The resolvent equals an explicit walk
over 400 traversals to `3.5e-18`, `#254`'s answer is its `n = 0` term, and the
relative error of that term is *exactly* `|L|` — an identity, not a fit.

**The branch series sums in closed form.** Short-way images all carry Maslov
factor `+1` and long-way images `−1`, so the winding sum is two geometric series
equal to `(e^{−uχ} − e^{−u(2π−χ)})/(1 − e^{−2πu})`, checked term-by-term to
`2.7e-15`. Its poles, as the regulator goes, are the conformal ESU
eigenfrequencies `ω = n+1`, with residues equal to the mode functions over `2ω`.
**The image and mode representations are one function** — the strongest
statement in the arc that the branch labels are a representation rather than an
approximation.

**And the solve adds events, not amplitudes.** The solved waveform *is* the sum
over history words `(a, c₁…c_n, b)` to `5.4e-06`; at echo times
`ℓ_a + Δ + n(ℓ_c+Δ) + ℓ_b` the control is at numerical zero and the solved field
is not (`3.3e+12`), on a `κⁿ` ladder, each echo signed by every Maslov factor in
its word. Those are arrivals at times `#254`'s ledger does not contain.

**The primitive is indexed by a pair of branches**, one per leg, and for a
reason that is not bookkeeping taste. `K_ab` carries the phase
`e^{−iω(ℓ_a + Δ + ℓ_b)}`, so `#253`'s closure condition is exactly the statement
that it is `ω`-independent — closed pairs have band coherence `1.000`, all
others below `0.091`. The **amplitude factorizes** over that index and the
**condition does not**: three pairs close inside the nine any single-index rule
would admit.

**Rank counts separable transfer channels, not histories** — a distinction this
round had to be corrected on. One throat carries `144` distinct `(a,b)` histories
at rank one; a second throat adds a second channel, in a shared topological
branch-label basis rather than by leg length, and the interference between the
two is a full fringe that is bilinear and therefore identically zero without
either.

**Finally, existence, convergence and stability came apart.** `1/(1−L)` exists
for any `L ≠ 1`; `|L| < 1` is only the radius of `Σ Lⁿ` and does not depend on
the delay; stability is `Im ω > 0` for every root of `D = 1 − L` in complex `ω`.
The coupling *displaces* the bare poles `ω = m + iγ` by
`δ_m = −ηκ e^{−imΔ} sin(md)/(4π² sin d)` — matched to `2.2e-04` — whose imaginary
part goes like `sin(md) sin(mΔ)` and changes sign with the mode. So `κ_series` is
`0.762` for every delay while `κ_stability` is `0.771` at `Δ = 1` and `3.034` at
`Δ = π`; at `κ = 1.520` in that gap the series diverges to `1.3e+119` while the
solve is finite and stable. **Solving and summing are not the same operation.**

## 14. A flux-conserving throat cannot ring up (`docs/throat_operator.md`, current round)

The previous round owed a boundary operator and named it as the next
construction. A point-supported throat is a **self-adjoint extension** of the
Laplacian on `S³ ∖ {M⁺, M⁻}`; von Neumann parametrizes those by a unitary
between the deficiency spaces — `U(2)` — equivalently, by Krein's formula, a
Hermitian `2×2` matrix `A`, with `M(ω) = A − Γ(ω)`. Everything follows from that
one substitution.

**It is definable at all.** The free Green function is
`G(χ,ω) = sin(ω(π−χ))/(4π sin χ sin(πω))` — real on the axis, poles exactly at
`ω = n+1`, matching `#255`'s branch series to `6.3e-12` — and its short-distance
split is `1/(4πχ) + g(ω) + O(χ)` with `g = −(ω/4π)cot(πω)`, remainder first
order in `χ`. The divergence is the universal Coulomb one, so the subtraction is
forced.

**The operator is a unitary `2×2` with both channels**, by Cayley transform, to
`4.4e-16`, with `|r|²+|t|² = 1` at each mouth. `#255`'s model has `r = 0` and
`|t| = κ`: outside `U(2)` unless `κ = 1`.

**Flux conservation is exactly Hermiticity.** The current through a small sphere
at a mouth is `Im(q_j* φ_j^reg)`, so the total absorbed is `Im(q† A q)` — zero
for *every* `q` when `A = A†`, at `1.8e-16` over 200 draws, against a median
`0.54` for the directional control.

**And therefore the spectrum is real, for every coupling.** Newton from complex
seeds — the same method `#255` used — converges only onto the axis: `0` off-axis
roots, worst `|Im ω| = 4.5e-18`. The directional control gives every root off
the axis, several growing, and **is unstable at `κ = 1` too**, where nothing is
lost. So the culprit was the *directionality*. `#255`'s three thresholds
collapse into one statement, and its instability is retired as an artefact of
its own non-conservation.

**The coupled spectrum interlaces the free one**, exactly two per gap over eight
gaps, because each channel function `g ± G_d` is monotone from `−∞` to `+∞`
across a gap. Switching the throat off returns `ω = n+1` — off being
`‖A‖ → ∞`, since the diagonal of `A` is an *inverse* scattering length — with
the shift falling like `1/‖A‖`, exponent `0.999`.

**What it costs.** The throat is point-supported: no interior, no proper length,
and therefore **no delay**. The `Δ` that carried `#251`–`#255` is not a
parameter of a self-adjoint point extension and does not survive into one. That
is a real loss of structure, not a simplification, and any statement from those
rounds that leaned on `Δ` has to be restated in terms of the mouth separation
and the boundary matrix.

## 15. What the arc cost in errors, and what caught them

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
* **Three conditions run together as one.** `|L| < 1` was called the condition
  for the coupled system to have a unique solution. It is the radius of
  convergence of `Σ Lⁿ` and nothing else: the solution exists for any `L ≠ 1`,
  and stability is `Im ω > 0` for every root of `1 − L(ω)` in complex `ω`. The
  three thresholds differ — at `Δ = π` by a factor of `3.98` — and the number
  originally quoted was blind to the delay while the real one is not.
* **A regulator artefact promoted to physics.** `κ_c ∝ γ` was read as "every
  coupling is critical as the damping is removed". `γ` is an Abel regulator and
  `T_d` carries the *bare* poles, so `1/max|T_d| → 0` is what a shrinking radius
  of convergence looks like, not an instability. A finite coupling *moves* the
  poles; the honest claim is about where the expansion fails.
* **A rank read as a count of the wrong thing.** `K` being rank one was called
  the number of independent histories the geometry supports. One throat carries
  `144` distinct `(a,b)` histories at rank one — an outer product counts
  separable *channels*. Worse, the two-throat sum was taken over matrices built
  from different leg lengths; it needed an explicit common branch-label basis
  before the rank meant anything at all.
* **A model named for what it was meant to be.** "The throat as a boundary
  condition" described a relation between field *values* with no flux matching,
  no reflected channel and `κ² < 1` power throughput. The solve was exact; the
  name promised a boundary operator the model does not contain.

The recurring lesson is narrow and practical: **a converged number is not a
correct number.** Three of the errors above survived grid refinement, and were caught only by
an independent construction — a closed form against brute
force, an operator form against a difference, an exact surface against a
truncated family.

The last rounds sharpen it in a different direction. Neither of the wave rounds'
errors was numerical, and refining anything would have caught none of them: **an
object drawn correctly can still be the wrong object, and only an independent
construction says which.** The repair was to make the picture carry a
measurement — the sign of `dA/dχ`, and a screen extent proportional to `sin χ`
with one constant to `3.6e-16`, which is `√(A/4π)`.

The coupling round pushes it one step further, and it is worth naming: **every
one of its four errors was a statement about a correct calculation.** The
resolvent, the pole displacement, the rank and the power ratio were all right to
machine precision *before* the corrections; what was wrong was which condition
the number was a condition for, which limit the scaling described, what the rank
counted, and what the model was called. No amount of numerical care reaches any
of that. What reached it was being asked to name the object precisely.

## 16. What is imported rather than derived

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
* In the ray-closure round, the throat *data* — mouth positions and delays —
  which is where its entire result lives: with the delays free instead of given,
  every event closes and nothing is selected.
* **Conformal coupling**, in the field round. `ξ = 1/6` is what makes `ω = n+1`
  and the branch structure sharp; the minimally coupled field has no images at
  all. That the ledger belongs to one and not the other is a *result*, but
  choosing the conformal field is an input.
* **The whole mouth-transfer model**, in the coupling round — not just its
  coupling `κ`. It relates field *values* through the free Green function, with
  no normal-derivative matching, no reflected channel, a `1×1` mouth scattering
  object where a conserving junction needs `2×2` unitary, and `κ²` power
  throughput. Solving it self-consistently fixes what a given `κ` *does*; it
  does not make the relation a boundary operator or a quotient. **Replaced** in
  the operator round by a self-adjoint extension — which promptly showed that
  the model's instability was its own.
* **The boundary matrix `A`**, in the operator round. Self-adjointness fixes the
  *form* — Hermitian `2×2`, four real parameters, a unitary by Cayley — and
  nothing else. Which four numbers is a choice, and `shells/junction.py` is
  still what would derive them from matter. The exotic-matter bill remains
  unpaid; what changed is that the thing being priced is now the right object.
* **The regulator `γ`**, likewise. A damping per unit path length is what makes
  the winding series converge. Every result is either `γ`-independent or reported
  as a `γ`-scaling, but that the physical answer is the `γ → 0` limit of this
  family is an assumption — and the bare poles sit at `Im ω = γ`, so the limit is
  where stability is decided.

## 17. What would come next

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

### And rank counting has now reached its end

The closing round was pushed to **branch-completeness** — the winding is bounded
by the delay (`k₁ + k₂ ≤ ⌊|Δ|/2π⌋`), so the feasible branch set is finite and the
union over all of it is still discrete. That is as much as constraint counting
can give.

What it structurally *cannot* give is a quantity that **vanishes** when a source
is removed rather than merely becoming underdetermined: deleting any scalar
equation costs a dimension, which is a theorem about square systems, not about
photons. The next discriminator has to be a field quantity, e.g.

```
𝒞(x) = A_A² A_B² (k_A · k_B)²
```

zero without a second source rather than under-determined by its absence.

The staged order that follows from this, and the reason it is that order:

```
ray closure → field solution → two-wave invariant
            → stationary action → backreaction → topological branch
```

* **field solution** — ~~done~~, in two stages, with one piece still owed. The
  first (`docs/field_solve.md`) solved the field with the mouth relation applied
  afterwards, and found the branches to be *exact support* rather than
  stationary-phase contributions, with a Maslov phase rays could not carry. The
  second (`docs/branch_coupling.md`) put the relation inside the equation, which
  turned out not to be a rearrangement: it adds arrivals the free-branch ledger
  does not contain, and it names the index the next step has to be written in —
  a **pair of branches**, the index on which the theory's *conditions* live even
  though its *amplitudes* factorize over it. **Still owed:** a genuine
  flux-conserving throat operator. What is solved so far is a rank-one
  *mouth-transfer* model — field values only, no normal-derivative matching, no
  reflected channel, `1×1` where a conserving junction needs `2×2` unitary, and
  lossy for `κ < 1`. That operator, not another visualization, is the immediate
  next construction;
* **two-wave invariant** — `𝒞` above, which is the sharp two-source falsifier
  the ray round could not supply. It now has both an index and a caution. The
  index: it should be a matrix in `(a,b)` in the same sense `K_ab` is, with the
  rank counting independent histories and the off-diagonal carrying the part
  that vanishes when a source is removed. The caution: anything built from a
  resummed field inherits the throat's resonances, and `κ_c ∝ γ` means those are
  unbounded as the regulator goes — so the test must be stated at fixed
  sub-critical gain, or it will be measuring the pole rather than the source;
* **stationary action** — evaluate the on-shell action and ask whether the
  candidate events are stationary. *Not* with Lagrange multipliers imposing this
  round's five equations, which would only rename them. This is where the
  retrocausal language earns its keep or fails: the backward-in-time throat
  contribution should fall out of one stationary solution rather than be
  narrated afterwards;
* **backreaction** — and the first GR question is not "does spacetime pinch
  off?" but whether `A + B` produces a collapse response not reproducible by
  rescaling `A` or `B` alone;
* **topological branch** — the detached resonator, last, and only if
  backreaction produces a finite-radius neck.
