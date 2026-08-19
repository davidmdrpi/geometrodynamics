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

## 14. Conservation is not stability (`docs/throat_operator.md`)

The previous round owed a boundary operator. A point-supported throat is a
**self-adjoint extension** of the Laplacian on `S³ ∖ {M⁺, M⁻}`, parametrized by
`U(2)`; written as the *pair* `B φ^reg = C q` — general enough to hold `#255`'s
relation, which is not of the form `φ^reg = A q` — the mouth-active spectrum is
`det(C − BΓ) = 0`.

**It is definable at all.** `G(χ,ω) = sin(ω(π−χ))/(4π sin χ sin(πω))`, real on
the axis, poles at `ω = n+1`, finite at the antipode, matching `#255`'s branch
series to `6.3e-12`; short-distance split `1/(4πχ) + g(ω) + O(χ)` with remainder
first order in `χ`, so the subtraction is forced.

**Hermiticity is exactly flux conservation** — `Im(q†Aq) = 0` for every `q`, at
`1.8e-16` over 200 draws, against `0.54` for a non-Hermitian control.

**And that is all it buys.** This is where the round's first version went wrong,
and the correction is its main content. `Γ` is real symmetric for real `λ = ω²`
of *either* sign, so `λ` is real — but nothing forces it positive, and `λ < 0`
means `ω = ±i√|λ|` with a growing mode. Two of the round's own three advertised
boundary matrices have exactly that, at `σ = 2.470532` and `σ = 7.090982`. They
were missed because the root search seeded only `Re ω ∈ [1.1, 6.9]` and
discarded roots leaving that window — a search that by construction could not
find a root on the imaginary axis. **The claim that a conserving throat cannot
ring up is withdrawn.**

**What replaces it is a stability region in closed form.** Both channel
functions fall monotonically along the imaginary axis, so stability is
`α + β ≥ g₀ + G₀` and `α − β ≥ g₀ − G₀` with `g₀ = −1/(4π²)` and
`G₀ = (π−d)/(4π² sin d)` — verified against a negative-`λ` scan at all 221 grid
points with 0 mismatches, and only 56 of them stable. Positivity is a genuine
restriction on the boundary data, separate from self-adjointness.

**Scope.** `det(C − BΓ) = 0` is the **rank-two mouth-active sector**, not the
spectrum: level `n` has degeneracy `(n+1)²` and only two combinations move, so
23 of 25 modes at level 4 never leave the free eigenvalue. Inside the sector
there is a mode *below* the free ground state that an `ω`-scan starting above 1
cannot see; and the antisymmetric channel does not span the first gap, because
the `n = 0` constant mode is equal at both mouths and its pole cancels.

**`#255` embeds exactly** as `B = [[0,0],[gain,0]]`, `C = I`, giving
`det(C − BΓ) = 1 − gain·G_d` — its own `1 − L` to `3.5e-18`. Maximal, not
self-adjoint, and unreachable by any finite Hermitian `A`. That is a
classification of its boundary condition and **not** a diagnosis of its poles: a
self-adjoint throat can be unstable too.

**What it costs.** The throat is point-supported: no interior, no proper length,
and therefore **no delay**. The `Δ` that carried `#251`–`#255` is not a parameter
of a point extension.

## 15. The positive sector is a light cone (`docs/throat_positivity.md`)

The previous round left three of the boundary matrix's four parameters open: it
showed that flux conservation does not give stability, and mapped the stable
region on the exchange-symmetric slice by scanning. The general answer is one
inequality.

```
non-negative  ⟺  A ⪰ Γ(0),   Γ(0) = [[g₀,G₀],[G₀,g₀]],
g₀ = −1/(4π²),  G₀ = (π−d)/(4π² sin d)
```

for **distinct non-antipodal mouths in the finite-`A` chart** — the scope
carries weight, and the two paragraphs before the last one are why.

**Why**: `dΓ/dλ ≻ 0` below threshold, so every eigenvalue of `A − Γ(λ)` falls
with `λ` while both run to `+∞` as `λ → −∞` — one crosses zero below threshold
iff it is already negative at it. Checked against a negative-`λ` root scan on
200 random Hermitian `A`, all with complex `β` and unequal mouths: **0
mismatches**, with both verdicts occurring.

**The monotonicity is an identity, not a sample.** `dΓ_ij/dλ =
⟨δ_i,(H₀−λ)⁻²δ_j⟩` is a **Gram matrix** — positive semidefinite for free,
positive definite exactly when the two mouths are distinct points. Rebuilt mode
by mode from the `S³` addition theorem and matching the closed form to
`8.1e-12` at every sampled `(d, λ)`, antipode included. This matters for what
the rest of the section is allowed to claim: a criterion advertised as exact
cannot rest on eigenvalues spot-checked at a handful of `λ`.

**And the same argument counts.** `#{eigenvalues < λ*} = #{negative eigenvalues
of A − Γ(λ*)}` for any `λ*` below the free ground state — a Krein-type inertia
theorem, 0 mismatches in 160 tests. Stability is its `λ* = 0` case, and the
*number* of growing modes comes out of it.

**The geometry is a forward light cone.** Hermitian `2×2` is `ℝ⁴` under
`A − Γ(0) = x₀I + x·σ`, and PSD is `x₀ ≥ |x|`: apex at `Γ(0)` with a doubly
degenerate zero mode, null boundary carrying exactly one, interior strictly
positive. Tested as a cone — convex, and closed under positive scaling *from the
apex*.

**The boundary is detectable rather than conventional**: the secular function
vanishes there to `1.8e-17` and the marginal mode is located by independent
root-finding at `1.4e-14`. Outside it the instability is continuous — `λ` linear
in the distance `ε`, so `σ` rises with exponent `0.50001`, with the coefficient
predicted from the eigenvalue slope rather than fitted.

**And the previous round's wedge is the `x₂ = x₃ = 0` slice** — exact there, and
wrong on 65 of 400 general draws when reused by averaging the mouths and
dropping `Im β`. Those two dimensions are exactly what a two-parameter scan
cannot see, which is the round's methodological point: the earlier result was
not wrong, it was *lower-dimensional*, and only the closed form showed which.

**Where the apex sits.** `tr Γ(0) = −1/(2π²)` at every mouth separation, its
eigenvalues are the previous round's two channel thresholds, and `det Γ(0) < 0`
for `0 < d < π` — so `Γ(0)` is indefinite there and **`A = 0` is unstable
wherever the mouths are actually apart**, which no placement short of the
antipode repairs.

**The exact antipode is a separate case, and for this geometry the interesting
one.** `G_d` has a *removable* singularity at `d = π`, not a pole:
`G_π(0) = +1/(4π²) = −g₀`, so `Γ(0) = g₀[[1,−1],[−1,1]]` has eigenvalues
`(2g₀, 0)` — negative *semi*definite — and `A = 0` there is **marginally
non-negative**, sitting on the cone's boundary with a static mode in the
symmetric channel. The earlier code raised on `χ = π` as a singularity; it is
one only in the formula, not in the function. This is the arc's recurring
lesson in a new dress: a limit rejected as a pole, where taking it changes the
verdict at the one configuration a through-throat geodesic on `S³` actually
picks out.

**And the inequality is stated on a chart.** `φ^reg = A q` requires `B`
invertible in the general pair `Bφ^reg = Cq`; the `U(2)` strata it misses are
Dirichlet directions, reached only as `‖A‖ → ∞`, so "the positive sector is
`A ⪰ Γ(0)`" is true of the chart and false of the whole family. The general
criterion is `A_eff ⪰ P†Γ(0)P` on the allowed-charge subspace — 0 mismatches
against the cone on 60 chart draws and against a root scan on 60 stratum draws,
with the reduction's own assumptions (Hermiticity of `A_eff`, the row-space
condition) checked rather than asserted.

## 16. Static two-source throat tomography (`docs/static_throat_tomography.md`)

**Not the two-wave invariant**, and the round is filed under what it actually
is. §10 ended rank counting by naming what it could not supply: a quantity that
**vanishes** when a source is removed rather than merely becoming
underdetermined. This round supplies one — for *static* probes — and then
establishes carefully what it does and does not measure.

Superposition decides the shape: every linear functional of a linear field is
additive, so the object has to be **quadratic**, and what carries the two-source
information is its cross term,

```
𝒞(y_A, y_B) = G(y_A,y_B) + Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B),  R = (A − Γ(λ))⁻¹
```

built from a functional that carries its own self-energy terms, so
`Q[a,b] − Q[a,0] − Q[0,b]` is a regression on the field construction rather than
a multiplication by zero. Its `N × N` table has rank two at every source count.

**What it is not** matters as much. It is a *static* kernel at a fixed spectral
parameter, with **no local null momenta**, so it cannot distinguish equal-energy
collinear from counterpropagating waves — the control the whole two-wave idea
rests on. And its index `(i,j)` labels **mouth channels**, not §11's geodesic
branches with their winding numbers and Maslov signs. The dynamical object built
from `T_A^{μν} T^B_{μν}` is still owed.

**Three things that look like the signature and are not**, and separating them
is most of the round. The cross term being nonzero is interference. The
interaction being **anisotropic** — depending on more than the geodesic
separation, which no free field on this background can do at all — is a real
`66%` effect that **two disconnected scatterers reproduce at `69%`**. And the
off-diagonal response block is nonzero for diagonal boundary data too, because
`Γ` couples the mouths through the ambient field; it is a *cross-mouth* channel.

**What discriminates is a parameter count.** The static invariant determines
three numbers — the entries of `S = Re R` — and two independent scatterers have
two knobs, so their image is a surface with the exact equation
`S₁₂ = G₀ det S`. The defect `𝒲 = S₁₂/det S − G₀` is its defining function, and
on the real sector

```
𝒲 = −β        exactly,  to 5e-16
```

independent of the self-energies, the separation, and the **Löwner margin** —
which is §13's resonance caution answered rather than managed. Scoped exactly,
`𝒲` detects **off-diagonal mouth-boundary mixing relative to the diagonal
two-scatterer null model**, inside this point-interaction model. It is also a
**protocol**: recovered by least squares from measured interaction energies by
an observer never told the boundary data.

**Which field is being solved turned out to decide the round's negative half.**
The first draft advertised a one-parameter family of connected throats invisible
to the test. Every such point needs `Im β ≠ 0` — and a **real** scalar, which is
what §12 solves, requires the self-adjoint domain to be conjugation-invariant,
`A = A*`, hence `β` real. Measured, not argued: with complex `β` a real static
source produces a **complex field**. So the blind family belongs to a
deliberately time-reversal-breaking complex extension, not to the arc's field,
and for the arc's field `𝒲 = −β` settles the question at one spectral parameter.
Inside the complex extension §15's gate removes one branch, and even the
remainder is a limit of the **probe** rather than the operator: phase-sensitive
complex sources give the full complex `R`, hence `A = Γ + R⁻¹` at a single `λ`.

**And the antipodal endpoint, tested as itself.** At `d = π` the static response
is singular as `A → 0`, so the invariant **diverges like `1/ε`** — while `𝒲`
stays exactly zero through four decades of it. The loudest available two-source
signal carries no information about whether the mouths are connected: **size is
not evidence**, in its most extreme form, because the size is unbounded.

## 17. The two-wave invariant is branch-resolved (`docs/two_wave_invariant.md`)

§16 built a *static* kernel and said plainly it was not the two-wave invariant:
no local null momenta, so it could not tell equal-energy collinear from
counterpropagating waves. This round solves the **time-dependent** field on the
throated background, builds the improved conformal stress tensors, and applies
that control.

The known WKB result is the **control, not the result**. What the round is about
is the difference between the exact multipath throat-coupled field and that
limit.

**The solver earns it.** The retarded field comes from Krein's resolvent formula
inverted along the contour `Im ω = ε`, which is exact rather than approximate;
every derivative is analytic, closing in form from `∇χ = −(y−(x·y)x)/sin χ` and
`∇∇χ = cot χ(δ_ab − n̂_an̂_b)`. Its free part reproduces §12's closed-form
winding-image sum to `3.3e-16` *including the Maslov signs* — two constructions
sharing no code. The conformal wave equation holds to `4e-16` relative with and
without the throat, and the improved stress tensor is traceless to `1.9e-15`,
which is a real test because `□φ` is taken from the solve rather than substituted
on shell: computed honestly the trace equals `φ(□φ − φ)`.

**The limit is recovered, with rates.** Sources and observer on one great circle
make the arriving directions exactly parallel or antiparallel to `1e-12`, so
WKB's `𝒩 = (1 − n̂_A·n̂_B)²` predicts 0 and 4 by geometry rather than by fitting.
The exact field gives `3.99995` head-on and `1.8e-10` collinear. And the
collinear null turns out to be *stronger* than a leading-order statement: on this
geometry the two wavefronts share their normal exactly, so amplitude gradients
cannot tilt either `k`, and the residue falls faster than any fixed power.

**Which is what makes the result large.** Holding the sources and the
observation point fixed and changing only *which branch has arrived*:

| branch pair | `𝒩` exact | geometry |
| --- | ---: | ---: |
| direct + direct | `1.905e-07` | `0` |
| long-way winding image + direct | `3.99806` | `4` |
| direct + via a mouth | `0.56501` | `0.56367` |

The winding image propagates the other way round the sphere — its phase is
`t + χ`, so the front moves *toward* the source and the arrival direction
reverses. A collinear pair reads head-on. The free control at the same instant
has **no second arrival at all** (`4e-29` against `1.2e-02`): the mouths *create*
the branch rather than bending one.

So the collinear null is not spoiled by curvature corrections, which are `1e-7`
here. It is spoiled by **multipath**, at `O(1)`. The invariant carries the branch
index §13 said it had to, and a single-branch WKB formula is not merely
approximate on this background — it answers a different question.

**The third row, audited.** Writing `j` for the mouth the source drives and `i`
for the mouth the signal leaves from, all four two-leg paths are enumerated
instead of minimised over. The prediction depends on `i` alone, so the four paths
carry **two** values — `0.563669` and `0.651935` — and the field must *pick*
rather than match a number. It picks correctly at each delay, to `8.1e-04`.

**And the control that scopes the whole row.** Re-run with `β = 0`, two
*disconnected* mouths: the invariant does not move. Swept over `β ∈ [0, 0.26]`,
every point inside §15's cone, `𝒩` shifts by `6.2e-07` — `7e-06` of the `0.0883`
separating the two exit mouths — while the channel's weight shifts `0.6%`. It has
to: `𝒩` is amplitude-normalized and a single channel is a single arrival
direction, so `β` rescales the weight without touching the geometry. **This
observable sees structure at the mouths, not the connection between them** — the
dynamical restatement of §16's lesson about the off-diagonal block. What sees the
connection is `𝒲 = −β`, recovered below from this same solve.

**`ΔT_{μν}`, and what it says the invariant does not.** `T` is quadratic, so the
two-wave content of the total stress tensor is the bilinear cross term
`ΔT = T[φ_A+φ_B] − T[φ_A] − T[φ_B]`, built from three evaluations of the same
functional rather than a hand-derived form. It is traceless to `1.8e-15` and
vanishes **exactly** when either source is switched off — §13's missing property,
now carried by a tensor. And it disagrees with the invariant completely:
`ΔT⁰⁰/√(T_A⁰⁰T_B⁰⁰)` is `2.000` collinear, the coherent maximum, against `1.044`
head-on. The interference energy is *largest* precisely where `T_A:T_B` is null.
**A backreaction estimate driven by `𝒞 = T_A:T_B` would look at the collinear
case, see nothing, and be wrong about its own source by the size of the whole
effect.** `ΔT` is what backreaction integrates; `T_A:T_B` is what the collision
invariant measures.

**And the smaller corrections have numbers.** *Tail:* `S³ × R` is conformally
flat, so the conformal scalar obeys Huygens **exactly** — between geometric
arrivals the free field is `1.4e-08` of its peak against `8.1e-02` with the
throat, a factor of `5.7e+06`. Every tail in this model is the throat's ringing.
*Caustic:* geometric optics diverges at the antipode where the exact kernel is
finite and *linear in* `ω`; the exact/WKB ratio is `|sin(ωe)|`, a function of
`ωe` alone, collapsing across three carriers to `6.6e-15`, so the caustic is cut
off at `e* ∼ 1/ω`. *Closure:* the DC content of the solved time series is exactly
§16's static kernel, and running that round's tomography on numbers from the
dynamic solver returns `𝒲 = −0.060010` against `−β = −0.06`.

**No backreaction**: the stress tensor is computed *from* the field and never fed
back. That is the next step, and it now has a concrete object to feed.

## 18. The throat has an interior, and the interior is the delay — and the point mouth is unstable (`docs/finite_conservative_throat.md`, current round)

Every round from §13 to §17 carried the same disclaimer — *point-supported, no
interior, no proper length, no delay* — and §19's ledger was more specific about
what stood in for one: a rank-one **mouth-transfer** model, field values only,
no normal-derivative matching, no reflected channel, `1×1` where a conserving
junction needs `2×2`, and lossy for `κ < 1`. This round replaces it.

**Two lines are put in.** A tube of length `L`, cross-section `𝒜` and interior
mass `m` joins the mouths; its Dirichlet-to-Neumann map is
`N(λ) = 𝒜k[[cot kL, −csc kL],[−csc kL, cot kL]]` with `k² = λ − m²`, and the
matching is value and flux continuity, `q = −NΦ`. Everything else follows.
Both entries are **even in `k`**, so the square root needs no branch choice, and
`det N = −(𝒜k)²` makes the chart closed-form: `A(λ) = −N(λ)⁻¹`, transmission
`β(λ) = csc(kL)/(𝒜k)`, self-energy `α(λ) = cot(kL)/(𝒜k)`.

**The one structural change is that the boundary condition is now
frequency-dependent, and that dependence is the interior.** A point throat is a
fixed Hermitian `A`; a finite throat is `A(λ)`. Every result below is a
consequence of that sentence.

**Where the self-adjointness lives.** The conservative object is the *enlarged*
system, ambient `⊕` tube, with the `λ`-independent matching — one self-adjoint
operator on `L²(S³) ⊕ L²([0,L])`. Eliminating the tube leaves a `λ`-dependent
boundary condition, and `A(λ)` is the **Weyl function** of that elimination: not
itself a self-adjoint operator on the ambient space, which an energy-dependent
boundary condition never is, but a matrix **Nevanlinna** function whose
monotonicity in `λ` *is* the enlarged system's self-adjointness showing through.
What is checked pointwise is that the elimination is faithful: `rank[B|C] = 2`
and `BC† = CB†` to `0.0` at seven values of `λ` on both sides of zero, with the
DtN map checked against the interior it summarizes by the **sesquilinear**
Green's identity to `1.5e-07` — against itself would prove nothing. §15's
`DirectionalThroat`, the control, has defect `0.30`: the size of the coupling
itself, which is what "lossy for `κ < 1`" was.

**The result: the throat transmits at the traversal time.** The measured object
is the **two-mouth block's** impulse response — `R(ω) = (A(ω) − Γ(ω))⁻¹`
inverted along the retarded contour. The source and observer legs are gone, but
`Γ` is the *ambient's* own mouth-to-mouth propagator and stays in, so this is
the coupled response and not the throat alone. Two different predictions, both
met: `r₁₁`, same mouth in and out, starts at **`t = 0`**, because a wave that
reaches a mouth is partly reflected *instantly*; and `r₁₂`, opposite mouths,
starts at **`min(L, d)`**, with `d(onset)/dL = 1.0071` against a predicted `1`
and a spread of `0.0` once `L` exceeds `d`.

That second number is the round's sharpest connection backwards. **The ambient
also connects the two mouths**, along a geodesic of length `d`, whether or not
they are joined — §16's cross-mouth channel and §17's `β = 0` control, which
those rounds could separate only by rank counting and by a sweep. Here they are
**separated in time**, and which arrives first is decided by `min(L, d)`. A
point throat transmits at `0.0`, which is what a point throat is.

The delay ledger is a derivation: on the contour `cot x = −i − 2iΣe^{2ikx}` and
`csc x = −2iΣe^{i(2k+1)x}`, verified to `4.5e-16` and `1.7e-15`, so the
same-mouth entry carries `0, 2L, 4L…` and the cross-mouth entry `L, 3L, 5L…`.
**The parities are the physics** — an even number of traversals returns to the
mouth it entered — and the reflected channel is the one the rank-one model does
not have at all.

**There *is* a point limit; it is not a finite `A`.** Freezing `A` at `A(λ₀)` is
exact at `λ₀` and `121%` wrong at `3λ₀` — a band of width `Δω ∼ 1/L` in
frequency. As `L → 0` the antisymmetric channel converges to `−L/(2𝒜)` while the
symmetric one diverges like `2/(𝒜λL)`, and a first draft concluded from that
that the limit does not exist. It does. A boundary pair is defined up to
`(B, C) → (MB, MC)`, so a diverging *chart matrix* means the limit has **left
the chart**; row-scaled, the pair converges linearly in `L` to
`(P_anti, −P_sym)` — `Φ_anti = 0`, `q_sym = 0` — a **mixed Dirichlet–Neumann**
stratum, maximal throughout and reached by no finite Hermitian `A`. That is
exactly the stratum §15's review said the chart does not cover, and it is also
what a very short pipe should do: short the mouths together and store nothing.
So §15–§17's constant-`A` family is this object read at one frequency.

**And the same zero mode breaks §16's tomography.** At `λ = 0` the static
response collapses onto `[[1,−1],[−1,1]]` to `4e-05` and `det S → 0` *linearly*
in `λ`, coefficient `149.08` constant to `1e-3` over four decades, so `𝒲`
diverges like `1/λ`. **What that falsifies is the generic finite-`A` family**,
every member of which is rank two — and *not* point-ness: the tube's own
short-tube stratum is rank one as well, `R → diag(0, −1/Γ_anti)`, and the tube
converges to it. The first draft claimed the stronger and wrong version. An
interior mass restores the rank, and off the collapse **`𝒲 = −β(λ)` exactly**,
to `3.1e-13`: §16's theorem survives the generalization and returns the
interior's own amplitude rather than a constant.

**And the candidate fails the stability gate.** `A(λ)` decreases and `Γ(λ)`
increases (§15's Gram identity), so `A − Γ` is strictly monotone between poles
and each channel has at most one root — a count, not a scan. The symmetric
channel always has exactly one, and it is at `λ < 0`: **an exponentially growing
mode, for every choice of parameters.** Three facts say whose it is, all limits
with their convergence measured: its rate matches `σ* = 2√(π/𝒜)` to `1.5e-03`
with **no `L` in it**; two mouth separations agree to `3.9e-09`; and the channel
splitting is `1.04·e^{−σ*d}`, the Euclidean propagator between the mouths, which
is the mechanism and not a bound on it. A mode that ignores the tube's length and
the mouths' separation, and does not distinguish the channels, is a
**single-mouth object**: the instability is the **point-mouth matching's**, not
the interior's.

**That is the round's closure result.** The retarded contour must clear `σ*` —
placed `0.03` below it the inversion returns a field with support *before its own
light cone*, a pedestal at 99% of the peak for an event that cannot begin until
`t = 0.6`, against `1.0e-16` placed above, and `σ*`'s closed form means the
contour is placed before the solve rather than diagnosed after it. But clearing
the contour evaluates the correct retarded solution **of an unstable system** and
cures nothing. Whether a finite-radius mouth or neck geometry removes the mode is
open, and it should be settled before §21's stationary action or backreaction —
either of those computed on this background would inherit the mode and measure
it rather than the physics they are after.

**Which frequencies.** The delay and the ledger are statements about the model's
analytic structure at **all** frequencies: a causal onset is a UV object, and
the pulse that resolves it carries content to `ω ∼ 30`, far above `σ* ∼ 1.4`.
They are exact results *about this model*, not predictions about a resolved
physical mouth. The static and low-frequency results sit inside the band.

**Still open:** the mouths are points — and §18 now says that costs a growing
mode — the interior is one-dimensional (so `𝒜` is a coupling and not a radius,
and a real tube's transverse modes above `ω ∼ 1/√𝒜` are inside the working
band), and `L, 𝒜, m` are chosen rather than derived. **No backreaction.**

## 19. What the arc cost in errors, and what caught them

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
* **A theorem's conclusion widened past what it says.** Self-adjointness makes
  the *spatial* eigenvalue `λ = ω²` real. It was reported as "the spectrum is
  real, so a conserving throat cannot ring up" — but `λ < 0` is real and means
  `ω = ±i√|λ|`, one member growing. Two of the round's own three examples were
  unstable, at `σ = 2.4705` and `7.0910`.
* **An instrument that could not see the answer.** The search for those roots
  seeded `Re ω ∈ [1.1, 6.9]` and discarded anything leaving the window, so a
  root on the imaginary axis was outside its reach *by construction*. The
  numbers it returned were all correct.
* **A control that was a different model.** PR #255's relation was represented
  as `A = [[0,0],[1/gain,0]]` with `B = I`, giving `g² − G_d² + G_d/gain`,
  where its actual pole condition is `1 − gain·G_d`. Everything concluded by
  comparing against it was unsupported. The true embedding needs a singular `B`:
  `q₁ = 0` is not of the form `φ^reg = A q` at all.
* **A parametrization read as a measurement.** The Cayley entries
  `(A − ic)(A + ic)⁻¹` were called reflection and transmission. Their magnitudes
  depend on the arbitrary reference scale `c` — the same `A` gives `0.955` and
  `0.713` — so they are boundary-mixing coefficients, and the physical
  conservation statement is the flux identity instead.
* **A removable singularity rejected as a pole.** `G(χ,ω)` was guarded against
  `χ = π` because `sin χ` vanishes there — but so does the numerator, and the
  limit is finite. The guard changed a verdict rather than a formatting detail:
  at the antipode `Γ(0)` is negative *semi*definite, not indefinite, so `A = 0`
  is marginally stable rather than unstable — and the antipode is the one
  configuration a through-throat geodesic on `S³` actually picks out.
* **A vanishing proved by multiplying by zero.** The two-source cross term was
  advertised as "identically zero when a source is removed" and demonstrated by
  computing `0.0 * answer`. That is a tautology, not a regression on the field
  construction. The honest version builds the quadratic functional with its own
  self-energy terms and takes `Q[a,b] − Q[a,0] − Q[0,b]` from three separate
  evaluations, where the large self-energies have to cancel.
* **A result named for the goal instead of for itself.** A static
  source-interaction kernel was presented as the roadmap's two-wave collision
  invariant. It has no local null momenta, so it cannot distinguish collinear
  from counterpropagating waves — the one control the two-wave idea exists to
  apply — and its index labels mouth channels rather than geodesic branches.
  Nothing computed was wrong; the label was, and a label decides which roadmap
  entry gets crossed off.
* **A negative result that depended on an undeclared change of model.** The
  "blind family" of undetectable throats needs `Im β ≠ 0`, and a *real* scalar
  field requires the self-adjoint domain to be conjugation-invariant — `A = A*`,
  hence `β` real. With complex `β` a real static source produces a *complex*
  field. So the round's headline caveat silently assumed a complex,
  time-reversal-breaking scalar it had never adopted.

* **A theorem on a chart, stated globally.** `A ⪰ Γ(0)` was presented as the
  positive sector of the `U(2)` throat family. It is the positive sector of the
  `B`-invertible *chart*; the Dirichlet strata are not reached by any finite
  Hermitian `A`, and on them the criterion is `A_eff ⪰ P†Γ(0)P` instead. The
  inequality was right, its quantifier was not.
* **An effect attributed to the wrong feature, for want of its null control.**
  §17's third branch was first described as arriving "through the throat", which
  credits the *connection* between the mouths for it. Running the obvious control
  — `β = 0`, the mouths disconnected — moves the invariant by `6.2e-07`, five
  orders below the signal it was supposed to explain. The branch is real and the
  number is right; what creates it is the presence of the mouths, and the free
  control alone could never have shown the difference, because it removes both.
  **A control that removes everything cannot tell you which part mattered.**
* **A limit written down as a fact.** §18's three identifying properties of the
  growing mode — the closed form `2√(π/𝒜)`, the blindness to the mouth
  separation, the degeneracy of the two channels — were first asserted with
  thresholds, and all three failed at the working point: `15%`, `0.18` and
  `0.87` against tolerances of `2%`, `1e-6` and `1e-2`. They are true *as
  limits*, in `σ*L` and `σ*d`, and the failures were the measurement pointing
  at the parameter it was being taken in. Restated with the convergence
  measured, the third one improved from a bound into a mechanism: the channel
  splitting is not merely small, it **is** `1.04·e^{−σ*d}`, the Euclidean
  propagator between the mouths. **A threshold that fails is often a limit
  asking to be named.**
* **A coordinate singularity reported as an absence.** §18's first draft said a
  finite tube has **no point limit**, on the evidence that its chart matrix
  `A(λ)` diverges as `L → 0`. But a boundary pair is defined only up to
  `(B, C) → (MB, MC)`, so that divergence says the limit has left the *chart*,
  not that it is missing: row-scaled, the pair converges perfectly well, to the
  mixed Dirichlet–Neumann stratum `Φ_anti = 0`, `q_sym = 0`. The corrected
  statement is strictly stronger — the short-tube limit *selects* a specific
  self-adjoint extension — and it also punctured the claim built on top of it,
  that a rank-one static response distinguishes a finite throat from a point
  one. It does not: the limiting stratum is rank one too, so what rank one
  falsifies is the finite-`A` chart. **§15's own lesson, arriving from the other
  direction: a quantity blowing up in a chart is a fact about the chart.**

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

## 20. What is imported rather than derived

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
* **The boundary data**, in the operator round. Self-adjointness fixes the
  *form* — a maximal pair `(B, C)` with `BC†` Hermitian, four real parameters —
  and nothing else. Which four numbers is a choice, `shells/junction.py` is
  still what would derive them from matter, and the exotic-matter bill remains
  unpaid. **Stability is also not fixed by the form**, though the positivity
  round narrowed the input sharply: the admissible set is now the closed-form
  cone `A ⪰ Γ(0)` rather than an unmapped region, and only `0.083` of a stated
  box lies inside it. *Which* point in the cone is still an input.
* **The regulator `γ`**, likewise. A damping per unit path length is what makes
  the winding series converge. Every result is either `γ`-independent or reported
  as a `γ`-scaling, but that the physical answer is the `γ → 0` limit of this
  family is an assumption — and the bare poles sit at `Im ω = γ`, so the limit is
  where stability is decided.

## 21. What would come next

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
ray closure → field solution → two-wave invariant → finite throat
            → ** RESOLVE THE MOUTH ** → stationary action → backreaction
            → topological branch
```

**The order now has a gate in it, and it is not optional.** §18 built the
conservative finite throat the arc had owed since §11 — and found that with
point mouths it carries an exponentially growing mode for *every* choice of
parameters, at a rate that knows neither the tube's length nor the mouths'
separation. Stationary action and backreaction are both integrals over a solved
field on this background. Run on a background with a growing mode, each would be
measuring the mode: an on-shell action evaluated on a solution whose amplitude
diverges is not stationary in any useful sense, and an `A`/`B`/`A+B` collapse
comparison would report the instability's response rather than the waves'. So
the next construction is **a finite-radius mouth or neck** — the ambient solved
outside two small balls rather than a point interaction with a radius parameter
— and the question it has to answer is single and sharp: *does the negative mode
survive?* If it does, the point-interaction throat family of §15–§18 is the wrong
model of a wormhole mouth and the arc has to say so. If it does not, the mode was
an artifact of the matching, the delay and the conservation law carry over, and
the two steps below resume with an object that has a proper length.

* **field solution** — ~~done~~, in two stages, with one piece still owed. The
  first (`docs/field_solve.md`) solved the field with the mouth relation applied
  afterwards, and found the branches to be *exact support* rather than
  stationary-phase contributions, with a Maslov phase rays could not carry. The
  second (`docs/branch_coupling.md`) put the relation inside the equation, which
  turned out not to be a rearrangement: it adds arrivals the free-branch ledger
  does not contain, and it names the index the next step has to be written in —
  a **pair of branches**, the index on which the theory's *conditions* live even
  though its *amplitudes* factorize over it. **The flux-conserving operator that
  round listed as owed is now built** (§18): a tube with an exact
  Dirichlet-to-Neumann map, whose elimination is faithful at every frequency to
  `0.0` against the transfer model's `0.30`, with the reflected channel, the
  normal-derivative matching and a proper length. What it costs is a
  frequency-dependent boundary condition; what it buys is a traversal delay,
  measured at slope `1.0071`. **And it fails the stability gate** — with
  point mouths it always carries a growing mode, whose rate knows neither the
  tube's length nor the mouths' separation, so the instability is the
  matching's. That is the thing to settle next;
* **two-wave invariant** — ~~done~~, in two rounds and with the second correcting
  the first's scope. §16 built the *static* kernel: zero without a second source,
  rank two at any source count, discriminator `𝒲 = −β`. §17 built the
  **dynamical** object it was explicitly not — time-dependent fields on the
  stable throat background, improved conformal stress tensors, exact
  `T_A^{μν}T^B_{μν}` — and applied the collinear/head-on control the whole idea
  rests on. The known WKB limit is recovered as a limit, with rates; the
  departure that matters is **multipath**, `O(1)` where every other correction is
  `1e-7`; and the invariant is **branch-resolved** in exactly §13's index, giving
  0, 4 or the mouth's angle for the same sources at the same event — the last of
  those audited over all four `(i,j)` channels, whose prediction depends on the
  **exit** mouth alone. Its `β = 0` control scopes the claim: the invariant sees
  the mouths, not the connection between them, which only `𝒲 = −β` sees. Tail,
  caustic and low-frequency closure all have numbers. What is *not* done is
  backreaction — and §17 hands that step a warning as well as an object, since
  the interference tensor `ΔT` peaks exactly where `T_A:T_B` is null;
* **stationary action** — *gated on the mouth above.* Evaluate the on-shell
  action and ask whether the candidate events are stationary. *Not* with Lagrange multipliers imposing the
  ray round's five equations, which would only rename them. This is where the
  retrocausal language earns its keep or fails: the backward-in-time throat
  contribution should fall out of one stationary solution rather than be
  narrated afterwards. It now has an object rather than a description — §17's
  `T_{μν}` is built from a solved field and traceless to `1.9e-15`, which is what
  an on-shell action is made of, and `𝒲 = −β` is a measured number the action
  round has to reproduce. It also inherits a warning: any quantity built by
  integrating over the field has to state which branches were present, because
  §17 measured the same configuration giving anything from 0 to 4 — and §18
  widens that warning, since the arrival ledger now contains `min(L, d)` and the
  throat contributes an echo series at every even multiple of the traversal
  time. The action round also inherits §18's conservation law, which is what
  makes a common action possible at all;
* **backreaction** — *gated on the mouth above.* The first GR question is not
  "does spacetime pinch off?" but whether `A + B` produces a collapse response
  not reproducible by rescaling `A` or `B` alone. It inherits §17's warning about *which* diagnostic
  to integrate — `ΔT`, not `T_A:T_B` — and §18's source with a finite size and a
  delay. It also inherits §18's **blocker**: computed on a background with a
  growing mode, a backreaction estimate measures the mode. The mouth has to be
  resolved first;
* **topological branch** — the detached resonator, last, and only if
  backreaction produces a finite-radius neck.
