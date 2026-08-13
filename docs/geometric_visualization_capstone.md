# The geometric-visualization capstone: what nine rounds of drawing settled (PRs #242–#250)

> **Framing.** Every round below is a *representation* question first and a
> physics question second: given a geometry and a field on it, which object
> should be drawn, and what does the choice smuggle in? The arc began by asking
> whether a wave could fold a surface through itself, and ended by deriving why
> a spherical model could not couple two surfaces at all. Each step is
> machine-checked by a probe in `experiments/closure_ledger/`.

## 0. The answer, stated first

**Four negative results and three positive ones, and the negatives are the load-bearing half.**

| | question | answer |
| --- | --- | --- |
| 1 | can a scalar height field on a slice wind through a glued bulk? | **no** — a graph has degree 0 identically |
| 2 | does the gluing scale rescue it? | **no** — it sets shape and magnification, not winding |
| 3 | can the converging ring reach across and intersect? | **no** — it reaches, and never crosses |
| 4 | does drawing the *vectors* rather than their tips change that? | **yes** — hundreds of crossings, threshold `ρ = 1/κ` |
| 5 | does focusing reach a singularity? | **no** — it reaches a *caustic*, and passes through |
| 6 | can a detached oppositely-glued shell replace exotic matter? | **no** — connected implies exotic, by identity |
| 7 | does `ℓ ≥ 2` supply the coupling `ℓ = 0` forbade? | **yes**, and it is screened as `(b/a)^ℓ` |

The through-line: **almost every apparent obstruction turned out to be a
property of the object being drawn, not of the physics** — until the last two,
where the obstruction turned out to be an identity that no drawing choice can
move.

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

## 8. Birkhoff is a vanishing eigenvalue (`docs/multipole_coupling.md`, PR #250)

The mutual stiffness of two concentric deformed shells is
`G m_b m_a ℓ(ℓ+1)(b/a)^ℓ/(a(2ℓ+1)²)`, verified to `9e-06` against brute-force
integration. The prefactor is the **Laplacian eigenvalue**, so the previous
round's decoupling is the `ℓ = 0` case of a multipole fact, not an imported
theorem. `ℓ ≥ 2` supplies the coupling — and the same formula screens it as
`(b/a)^ℓ`, `544×` from `ℓ = 1` to `ℓ = 8` at `b/a = 0.4`.

## 9. What the arc cost in errors, and what caught them

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
* **A zero mode that was not one.** A pure `P₁` deformation is not a translation
  past linear order; its area cost is `4/3`.

The recurring lesson is narrow and practical: **a converged number is not a
correct number.** Three of the nine errors above survived grid refinement, and
were caught only by an independent construction — a closed form against brute
force, an operator form against a difference, an exact surface against a
truncated family.

## 10. What is imported rather than derived

* Birkhoff's theorem (`shell_junction`) — and then *not* needed in
  `multipole_coupling`.
* The Darmois–Israel formalism.
* `β²` as a free parameter at the equilibrium.
* The fixed round `S²` background of every wave round: curvature `1` everywhere,
  at every time, with no Einstein equation and no backreaction.

## 11. What would come next

The honest next object is not another drawing. Three of the closing results name
their own missing ingredient:

* a curvature singularity needs **metric evolution with backreaction** — the
  congruence round could not have produced one;
* a finite-radius reconnection needs the congruence to **act back** on
  something, which a test field never does;
* an `ℓ ≥ 2` resonance needs a **shear modulus**, which no equation of state
  supplies and which this round carries as an explicit input.

Each is a construction, not a visualization.
