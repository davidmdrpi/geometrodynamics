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

## 18. The throat has an interior, and the interior is the delay — and the point mouth is unstable (`docs/finite_conservative_throat.md`)

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
mode, for every choice of parameters.** In the `σL, σd ≫ 1` limit its rate
matches `σ* = 2√(π/𝒜)` to `1.5e-03` with **no `L` in it**, two separations agree
to `3.9e-09`, and the channel splitting is `1.04·e^{−σ*d}` — the Euclidean
propagator between the mouths, the mechanism rather than a bound on it. So the
mode is generated at the **point-mouth/tube interface** and *localizes to a
single mouth in that limit*.

The working throat is **not** in it: at `𝒜 = 4π` the asymptotic form gives
`σ* = 1` while `L = 0.9` gives `1.417`, and `σ*` runs `1.769 → 1.152` across
`L = 0.4 → 3`, a spread of `54%`. A first draft claimed the mode "belongs to the
mouth and not the interior"; the interface statement plus an asymptotic
localization is what the data support, and it is also what makes §19's
finite-radius mouth exactly the right discriminator.

**That is the round's closure result.** The retarded contour must clear `σ*` —
placed `0.03` below it the inversion returns a field with support *before its own
light cone*, a pedestal at 99% of the peak for an event that cannot begin until
`t = 0.6`. **And clearing it is necessary, not sufficient**: at a clearance of
`+0.02` the contour is above the mode and the pedestal is still `2.6e-03`,
because that clearance is `0.95` of the frequency spacing `2π/span` and the grid
does not resolve it — §17's lesson arriving a second time. The rule is
`ε > σ*` *and* `ε − σ* ≫ 2π/span`; both have closed forms, so both are checkable
before the solve, and at `14–72` spacings the pedestal is `1.0e-16` with the
onset converged to four time steps. Neither condition cures anything: above `σ*`
the inversion returns the correct retarded solution **of an unstable system**.
Whether a finite-radius mouth or neck geometry removes the mode was left open
here, and it had to be settled before §24's stationary action or backreaction —
either of those computed on this background would inherit the mode and measure
it rather than the physics they are after. §19 settles the mouth and §20 the
neck; the answer is no, twice, from two constructions.

**Which frequencies.** The delay and the ledger are statements about the model's
analytic structure at **all** frequencies: a causal onset is a UV object, and
the pulse that resolves it carries content to `ω ∼ 30`, far above `σ* ∼ 1.4`.
They are exact results *about this model*, not predictions about a resolved
physical mouth. The static and low-frequency results sit inside the band.

**Still open:** the mouths are points — and §18 now says that costs a growing
mode — the interior is one-dimensional (so `𝒜` is a coupling and not a radius,
and a real tube's transverse modes above `ω ∼ 1/√𝒜` are inside the working
band), and `L, 𝒜, m` are chosen rather than derived. **No backreaction.**

## 19. The negative mode does not survive a finite mouth (`docs/finite_radius_mouth.md`)

§18 ended by gating the roadmap on one question: does its growing mode survive a
finite-radius mouth? **No — and the statement is structural rather than
parametric.**

**What a resolved mouth changes.** A point interaction has to subtract the
`1/(4πχ)` divergence and keeps the *renormalized* self-energy `g(λ)`, which is
**negative** — the finite part left after an infinite subtraction, and what the
mode fed on. A sphere needs no subtraction. Smearing the coupling over `∂B_a`
with the same operator on both sides — so the composite stays manifestly
self-adjoint — replaces `g` by the **unsubtracted** `𝒢_self = f(a)G(a)`,
`𝒢_cross = f(a)²G(d)`, with `f(χ,λ) = sin(ωχ)/(ω sin χ)` the regular radial
solution. Both are mean-value identities, and both are checked against direct
quadrature on `S³` — the cross one to `1.0e-10`, the self one to `4.1e-04`,
grid-limited by the singularity at coincidence and reported as such rather than
as a model error.

**The signs settle it.** At `λ = −σ²` the tube's channel functions are strictly
**negative** — a passive interior supplies no restoring force — and the resolved
mouth's are strictly **positive**, once written so the positivity is visible:
every bracket in `(0,1]`, and `2κa − κd < 0` because disjoint mouths need
`a < d/2`. A difference of a negative and a positive number has no zero, so **no
parameter choice can produce a growing mode**; a sweep of `3078` points finds `0`,
worst approach `−5.1e-04`.

**And §18's mode was the linearization.** That round froze the mouth at the
*constant* `1/(4πa)`, the leading term of `G(a,λ) = 1/(4πa) + g(λ) + O(a)`. The
exact `G(a,−κ²)` is **screened**, `≈ e^{−κa}/(4πa)`; the constant is not, so it
eventually beats the tube's `−1/(𝒜κ)` and crosses — at `κa = 1.0004, 1.0025,
1.0221, 1.1269` for `a = 0.02, 0.05, 0.15, 0.35`. **The root sits at the edge of
its own approximation.** The two models agree to `0.8%` for `κa ≤ 0.1` and differ
by `1000%` at `κa = 3`, disagreeing not in magnitude but in *sign*. §18 measured
that its mode knew only the mouth's scale and suspected exactly this; this is the
demonstration.

**Where the mode went: soft, and positive.** One state below the free gap
`λ = 1`, in the symmetric channel, at `λ₀ → 8πa/(𝒜L)` — two mouth capacitances
`4πa` restoring a tube of volume `𝒜L`, ratio `0.998` at `a = 0.005`. So the point
limit drives it to zero **from above**, and §18 took a mode approaching `0⁺` like
`a` and put it on the other side, at `λ ≈ −1/a²`.

**The good results survive.** The traversal delay keeps slope `1.0010` in `L` and
saturates at the ambient path to `0.0`; the static rank-one collapse and §16's
`𝒲 = −β(λ)` hold to `3.6e-12` — all of them properties of the *tube's* zero mode,
which the mouth does not touch. And the contour is easier: `ε = 0.4` where §18
needed `ε > σ* ≈ 2`.

**Still put in:** one channel per mouth, so only `ℓ = 0` couples — the dropped
multipoles obey §8's screening law, dipole/monopole `= 0.934·(a/d)` across a
decade — and the mouths are spheres in a *fixed* ambient, not a solved neck.
**No backreaction.** §20 closes the first two of those.

## 20. The balls removed, and the answer made a theorem (`docs/finite_radius_neck.md`)

§19 answered §18's gate and named its own weakest point in doing so: its mouths
were **spheres in a fixed ambient**. The balls were never removed — it smeared
the coupling over `∂B_a` while still using the *whole sphere's* Green function —
and only `ℓ = 0` coupled to the tube. This round removes them:
`Ω = S³ ∖ (B_a(c₁) ∪ B_a(c₂))`, tube glued along the boundary spheres.

**The answer is still no, and it is now a theorem.** With the balls removed there
is **no subtraction anywhere**, so

```
E[φ,u] = ∫_Ω (|∇φ|² + φ²) dV  +  𝒜 ∫₀^L (|u'|² + m²|u|²) ds
```

is a sum of non-negative terms. `E = 0` forces `φ ≡ 0` on `Ω`; matching then
gives `u(0) = u(L) = 0`, and Poincaré gives `𝒜∫|u'|² ≥ (π/L)²𝒜∫|u|²`. Hence
`λ > 0` for **every** configuration — all multipoles, no truncation, no sweep.

This is a change of footing rather than a refinement, and §22 records it as its
own error species: **§19 had to establish a sign on a reduced `2×2` and support
it with a `3078`-sample sweep; this round has positivity of the form itself.**
The first is an argument about a model, the second about the problem.
The renormalized `g(λ) < 0` that §18's mode fed on has nowhere to enter, because
nothing is renormalized.

**The object the theorem is about, checked.** The exterior DtN
`N_ℓ(λ) = −4π sin²a·ψ'(a)/ψ(a)` comes from shooting
`v'' + [λ − ℓ(ℓ+2)/sin²χ]v = 0` from the far pole (regularity picks the `e^{ℓ+1}`
branch), and matches an independent `ℓ = 0` closed form to `1.7e-14` — which is
the check the `ℓ ≥ 1` channels, which have no closed form, inherit. It is
positive in every channel and **increasing in `ℓ`**, so the monopole is the
softest; and `N₀ → 4πa`, the capacitance of a sphere, which fixes the
normalization as physical. Explicit trial configurations give Rayleigh quotients
from `0.359` up, all above the computed lowest mode `0.1075`, and the degenerate
purely interior case lands on `π²/L²` to `2e-07`. A `1197`-sample sweep finds `0`
roots, worst approach `−1.6e-03`.

**The monopole truncation was never a stability limitation.** A one-channel tube
presents a single number at each mouth, so it drives `ℓ = 0` and nothing else;
the `ℓ ≥ 1` sectors are the exterior's own modes, `1.45×` stiffer or more, and
can neither be driven nor go unstable. §8's `(a/d)^ℓ` screening law bounds missed
**amplitude**, not the answer — §19 listed this as its main open approximation,
and it turns out to have been open in a direction that did not matter.

**What the fixed ambient cost, priced.** `f(a)G(a)` against `1/N₀`: `1.3e-04` at
`a = 0.02, λ = 0`; `3.8e-03` at the working radius; `11%` at `a = 0.35, λ = −4`.
The error tracks the fraction of the sphere wrongly left in. And the one
approximation this round has *not* removed is priced in the same place: the
reduced `2×2` keeps a **single-scattering** cross term, whose neglected series
expands in `cross/self = 0.8·(a/d)` — `9.5e-04` at the working point, at worst
`0.16` anywhere sampled. Too small to flip a sign, and outside the theorem,
which does not go through the reduced model.

**The soft mode is forced, not found.** The same style of argument that kills the
growing mode produces this one, from the two ends of the gap: `F_sym → +∞` as
`λ → 0⁺`, and `→ −∞` as `λ → 1⁻` because the exterior stiffness **vanishes** at
the free ESU threshold, `N₀ → 2π(π − a + sin a cos a)(1 − λ)`, reproduced to
`3.1e-05` and first order exactly (error `×0.1` per decade). A continuous
function running `+∞ → −∞` has a zero.

**And one correction to §19.** Its "exactly one state below the gap" is a
statement about `L < π`, not a structural one. The channel functions have
**poles** at `λ = (2πj/L)²` and `(π(2j−1)/L)²`; above `L = π` these enter the gap
and each brings a genuine extra state just above it — three at `L = 8`. A pole is
a *sign change with no zero*, so crossing-counting alone reports states that are
not there; the residual separates roots (`~1e-15`) from poles (`~1e+15`) by
fifteen orders of magnitude, so the discrimination is not a tuned threshold. None
of it touches the stability answer: every one of those states is positive, as the
form requires.

**Still put in:** the tube has **one transverse channel**, so `𝒜` is a coupling
and the neck is a quantum-graph edge, not a solved cross-section; the ambient
metric is **fixed**; **no backreaction**.

## 21. `A + B` does what rescaling `A` or `B` cannot (`docs/metric_backreaction.md`)

§18 gated backreaction on a growing mode; §19 and §20 closed that gate. So this
round asks the roadmap's first GR question — and deliberately *not* "does
spacetime pinch off?":

> does `A + B` produce a metric response that rescaling `A` or `B` alone cannot?

**Yes.** At the working point `0.921` of the interference response lies outside
everything rescaling can reach, and the interference term is **comparable in
size** to the single-wave responses rather than a correction to them.

**Why it is a linear-algebra question.** The field equation is linear so the
fields add; `T` is quadratic so `T[A+B] = T[A] + T[B] + ΔT` with `ΔT` bilinear
(measured: self term `∝ c²`, cross term `∝ c`, exactly zero with one source
off); linearized Einstein is linear so the responses add. Rescaling `A → cA`
sends `β_A → c²β_A`, so everything reachable is the two-parameter family
`{c²β_A + d²β_B}` and the question is a projection residual. Reported off the
full linear **span**, which strictly contains that cone — so the figure is
conservative.

**The channel is forced, not chosen.** The ESU is held static by a fluid this
arc never specifies. A perfect fluid carries no anisotropic stress — in an
orthonormal frame `T_ab = diag(ρ,p,p,p)` whatever the anisotropy — so neither it
nor `Λ` touches the traceless spatial part. **The transverse-traceless sector is
the only channel whose answer does not depend on what was never put in.** The
scalar sector depends on the sound speed and carries the Eddington mode.

**The response, derived.** Cartan about the ESU in the homogeneous anisotropy
gives `δG^TT_{ij} = β̈_{ij} + (8/a²)β_{ij}`, so `ω₃ = 2√2` and `ω₃² > 0` — the
tensor sector is **stable**, which is §18's gate applied to this round's own
background. The connection comes from *solving* the torsion-free condition, after
a first draft's remembered formula produced a `G₀₀` containing `ä`. Three
validations pass — round `S³`, ESU (independently reproducing `two_wave`'s
`_EINSTEIN`), closed FRW — the first-order piece is automatically traceless, and
`δG_{0i} = 0` identically.

**And it is not a universal constant.** `0.88–1.00` across time windows,
`0.56–0.99` across carriers, `0.929` with the throat and `0.986` without. Large
everywhere, constant nowhere — so the range is the headline and the claim is
that the interference response is *generically* unreachable.

**The resonance test reverses itself when done on the coupled system.** A first
version argued the channel was off resonance *by construction*: the conformal
scalar on the ESU has integer spectrum `ω_n = n+1`, so a free-field source rings
on integers, and `ω₃ = 2√2` is irrational. True — **of the uncoupled ambient**.
The throat rings where `det(A − Γ(ω))` vanishes, and those zeros are not
integers; the nearest sits `0.050` from `ω₃`, with a local spacing of `0.17`.
Across sixteen configurations it is always within `0.086` and at `d = 0.9` it is
`0.001`. **Off resonance with the free ambient, generically near-resonant with
the coupled one** — which is also the mechanism the first version lacked for the
carrier sensitivity.

**And the whole TT sector is areal-blind, not just the homogeneous mode.** A
traceless shear preserves volume *exactly* and mouth area *to first order*
(second-order coefficient `0.403`, positive). The natural next move — build the
inhomogeneous harmonics — does not help: `δA/A = −½⟨h_nn⟩`, and that average
vanishes **identically** for any transverse-traceless field, because
`⟨n_in_jf(k·n̂)⟩` can only be `Aδ_ij + Bk̂_ik̂_j` and tracelessness kills the
first while transversality kills the second. Measured to `4.5e-15` against a
transversality-dropped control of `1.5e-01`, and it survives curvature because
`S³` is maximally symmetric — the normal's moments are isotropic to `6e-16`. So
§21 cannot say whether the interaction metric moves toward a neck or away, and
**neither can any refinement of it**: a signed areal response has to come from
the scalar sector, where the fluid enters. The price of choosing the only
fluid-free channel, now measured, is that it cannot answer a question about
area.

**The superseded phrasing.** A traceless shear preserves volume *exactly*
and mouth area *to first order*; the second-order coefficient is `0.403` and
**positive**. So §21 cannot say whether the interaction metric moves toward a
neck or away — it distorts the mouth into an equal-area ellipse. That question
needs the areal sector on a resolved neck, which §21 is not.

**Still put in:** the `n = 3` harmonic only; a **fixed** ESU; **point** sources
and §18's **point** throat rather than §19/§20's resolved mouths; and a
**linear** response, not a solved coupled system.

## 22. The interference metric deforms toward a neck (`docs/signed_areal_response.md`)

§21 answered *does `A + B` do something rescaling cannot* with a large yes, and
then proved that its own channel could not answer the question actually asked:

> `δA/A = −½⟨h_nn⟩`, and that average vanishes **identically** for any
> transverse-traceless field. Tracelessness kills the isotropic part of
> `⟨n_i n_j f(k·n̂)⟩` and transversality kills the rest.

So the geometric question moved to the **scalar** sector, posed as an
**initial-data constraint solve** rather than an evolution. That is what removes
the two reasons §21 avoided that sector: on a maximal slice the `K` terms in the
Hamiltonian constraint are quadratic, so `δR⁽³⁾ = 16πG δρ` has no time
derivatives in it, and a constraint has neither a sound speed nor an Eddington
mode.

**Toward a neck — at the wide working throat, off resonance. Both mouths
close.**

    ΔA/A = ( −2.06e−03 , −1.88e−03 )     in units of 2πG

Negative in all eight controls. **The mechanism is not the obvious one:** the
interference energy alone would *open* the mouths — `U(c_j) > 0` at both — and
the throat's own monopole layers overshoot that and invert it. The throat is a
low-impedance load for the constraint; it cannot support the conformal factor
the energy piles up around its mouths.

**The structure.** With `g = ψ⁴ĝ` and `ψ = 1+u`, the constraint is
`∇²u + 3u = −2πG δρ` and `ΔA/A = 4u`. That operator is exactly degenerate on
`S³` (`4 = (n+1)²` at `n = 1`), so it has no solution on the closed sphere at
all; removing the two mouth balls is what makes it solvable. The field is a
source term plus monopole and dipole layers on the two mouth spheres plus a free
kernel element — twelve unknowns, twelve equations — and the solvability
condition `Σ_j A_j c_j + Σ_j D_j = S_σ` is an **identity**, because
`(2/π²) cos χ` is exactly the projector onto `span{x^A}`.

**Two things came out the other way from expectation, and both are recorded as
such.**

*The dipole layers are required, and then invisible.* Two monopoles sweep only
the plane of the two mouths, and the monopole-only solvability condition fails
by `62.5%` of the obstruction — without the `ℓ = 1` layers there is no solution.
That was expected to be the heart of the round. It is not: the response to the
off-plane obstruction is `6e−17`. The layers deposit it in the kernel elements
`x²` and `x³`, which vanish at both mouths because both mouths lie in the
`(x⁰, x¹)` plane. A first draft of the module docstring said the `ℓ = 1` sector
*is* the calculation. It is required for the calculation to **exist** and
contributes nothing to its **value**.

*The answer rests on the best-converged input, not the worst.* `ΔA/A` is, to
`0.09%`, a linear functional of the obstruction alone. The `ℓ = 1` source moments
drift `41%` between quadrature levels, against `1.5%` for the obstruction —
and scaling them by three, by zero, or replacing them with noise moves the answer
by `5e−04`.

**The sign is a statement about a throat, not about the sphere** — and that was
put to the test rather than left as a caveat. The tube's `ℓ = 0` constraint
channel is `∂_s² + 4π/𝒜`: a **cavity**, with poles at `kL = nπ` (the scan finds
flips at `3.133` and `6.260` against closed forms `π` and `2π`). The working
tube carries `𝒜 = 4π` against a mouth sphere of area `4π sin²a` — wider than its
own mouths by `400×` at `a = 0.05`, which is what puts it at `kL = 0.9`, inside
the first cell. **Set the two equal and the sign does not survive:**

| `a` | `kL/π` | `𝒜 = 4π` | matched `𝒜 = 4π sin²a` |
|-----|--------|----------|------------------------|
| 0.05 | 5.732 | closes / closes | **opens / opens** |
| 0.10 | 2.870 | closes / closes | **closes / opens** |

At `a = 0.05` both mouths open; at `a = 0.10` they disagree. The matched throat
is past five poles and past two, under `5%` of its own length from the next one.
The neck's cross-sectional area has been a free parameter since §19, carried
along because nothing measured had ever depended on its value. Something does
now — and §23 answers it, by showing the area and length were never free at
all, and that the sign reverses on the throat the constraint forces.

**What is still put in.** The ESU's fluid is held rigid — consistent because the
scalar's stress tensor is separately conserved, but a responsive fluid is the
next refinement. The source is §21's, on a fixed background with point sources.
And the solution is dominated by a nearly-zero mode with `|c| ≈ 1.7` against a
mouth response of `2e−03`, so the perturbative window is set far from the throat.

## 23. Which throat is physical, and the sign reverses (`docs/which_throat_is_physical.md`)

§22 ended by finding that matching the tube's area to its own mouths flips the
sign of `ΔA/A`. The throat's geometry had been a free parameter since §19,
carried because nothing measured had ever depended on its value. Something did
now, so the question could not be deferred.

**They were never free.** On a maximal slice the background constraint is
`R̂ = 16πGρ̄`, so a profile does not choose its matter — the matter is whatever
the profile implies. A product tube of area `𝒜` is `S²(r) × R` with `R̂ = 8π/𝒜`,
so §19–§22's `𝒜 = 4π` implies `ρ_tube = ρ̄/3`, and the matched tube implies
`133 ρ̄`. Neither is the ambient's own fluid, and neither was chosen for a
reason.

**The throat that is forced.** Ask for one needing no matter (`R̂ = 0`) and
gluing on with no surface layer. `R̂ = 0` integrates once to `f'² = 1 − f₀/f`,
and smooth gluing to the round `S³` at mouth radius `a` forces `f₀ = sin³a`.
Two conditions, two unknowns, **nothing left over**. Length and resistance
follow in closed form — `L = 2[sin³a·arccosh(1/sin a) + sin a cos a] ≈ 2a` and
`I = 4 cos a/sin³a` — each checked against quadrature to `1e-12`, and the
conductance comes out exactly `N₀(a,4)/4` at every radius.

**And the throat has a name, which derives its mass.** `R̂ = 0`, `K = 0` and a
spherical neck do not merely permit an Einstein–Rosen bridge — they are one.
`f'² = 1 − f₀/f` is the time-symmetric Schwarzschild slice written in proper
radial distance, with `f₀ = 2M`. So the forced neck radius is twice a mass:

    M = sin³a / 2

the throat's mass derived from the size of the excised mouth, with nothing left
to choose. The Schwarzschild parameter, the irreducible mass `√(A_neck/16π)` and
the Hawking mass agree to `1.3e−13`, and the neck area is `16πM²` exactly. Better
still, **the gluing condition *is* a mass statement**: the Hawking mass
`(f/2)(1 − f'²)` is `f₀/2` on the throat and `sin³χ/2` on the round `S³`, so the
seam is smooth exactly when it does not jump. The mouth radius and the throat
mass were never two parameters — the ambient already assigns a mass to every
sphere it contains, and the throat has to carry the one at the cut.

Four things it does not say, each asserted in the tests because the claim is
strong enough to be worth not overstating: there is **no asymptotic region**, so
no ADM mass — what is derived is quasi-local, and unambiguous only because the
Hawking mass is constant on the vacuum piece; it is **dimensionless**, `M/R`
against the ESU curvature radius, which is all §52's scale-modulus theorem
allows; both ends sew into the **same** `S³`, a handle of Misner's kind rather
than a bridge between universes; and the neck is a **marginal sphere of this
slice**.

That last one is narrower than it first read here. What is proved is an
identity — `H = 0` at the neck with `K_ij = 0` on the slice gives
`θ_+ = θ_− = 0`, so the neck is a minimal surface and a MOTS of the initial
data. It is **not** shown to be an *apparent horizon*, which is the
**outermost** MOTS and needs a global condition on the slice that nothing here
evaluates; and it is **not** shown to be *non-traversable*, which is a property
of the Lorentzian development and needs a lapse that nothing here chooses.
Non-traversability does follow if the development is additionally taken to be
the standard vacuum Schwarzschild/Einstein–Rosen one, and it is under **that**
added assumption — not from the data — that this becomes the vacuum face of
§7's "connected implies exotic".

**There is no cavity.** The constraint operator is `∇² + R̂/2`, so `R̂ = 0`
leaves the plain Laplacian: `(f²u')' = 0`, solutions `A + B∫ds/f²`, monotone.
§22's resonances at `kL = nπ` and the sign flips across them were **properties
of matter in the tube**. On the physical throat there is no resonance to be off.

**And the sign reverses.** `ΔA/A = (+6.64, +8.58)` in units of `2πG`: both
mouths **open**. The mechanism is a single number. Split a symmetric two-port
as `Y = G·[[−1,1],[1,−1]] + shunt·I`; the shunt is the flux a *uniform*
potential drives into the throat, and `(f²u')' = 0` makes it vanish
**identically** for a vacuum tube. Scanned over eight orders the conductance
never changes the sign; the shunt passes through a pole near `2e−03` and flips
it. §22's tube sat at `6.07`, three orders past.

**What it costs, said out loud.** The response is `3000×` §22's and grows as
`a⁻³`, with the condition number rising alongside. That is the physics of a
throat with zero shunt by identity: it barely lifts the constraint's exact
degeneracy — the `k = 1` kernel at `4 = (n+1)²` — so the linear response sits
close to a mode the operator nearly annihilates. **The sign is robust; the
amplitude at which linearising was legitimate is now the binding question**, and
this round does not answer it. Nor does it discharge the other one the mass law
opens: whether a throat that is an Einstein–Rosen bridge — one whose neck is a
marginal sphere already on the initial slice — is one this arc's larger
framework can accept. That question now has a geometric step in front of it,
which this round also leaves open: whether that MOTS is the outermost one, and
what the Lorentzian development actually is. Together they decide whether
`M = sin³a/2` is a result about the model or only about this slice.

## 24. Two waves connect where one could not (`docs/two_wave_slice.md`)

> **Superseded reading.** §25 shows this section's *construction* was
> ambiguous — it drew two curves and read their overlap as a connection.
> Every number below survives, because the quantity it plots is the
> one-surface deformation; but the two configurations swap names, which
> inverts the headline. Read §25 alongside it.

A revisit rather than a new result, taken before the nonlinear arc starts,
because §4's picture was doing work the rest of the arc had outgrown.

v46 put **one** scalar wave on the great circle through a source and its
antipode, drew it as a radial height in the vacuole, and glued `R_outer` to
`R_inner`. Its finding was negative and sharp: the curve is a **graph**
`r = f(σ)`, so its radial winding number is identically zero — every outward
crossing of the seam is paid for by an inward one. **A height field cannot
wind**, and one wave running to its antipode never meets itself.

Everything since has needed **two**. So: one wave driven outward, one driven
inward, both refocusing at the antipode — do they connect, inner to outer?

**Yes, at the antipode, on the seam, at the refocus.** At threshold one curve
sits at exactly `R_inner` and the other at exactly `R_outer`, which after gluing
is one point.

**And the threshold is not a new number.** A single wave crosses the seam when
`εu = gap/2`; the pair spans `2ε|u|` of the radial circle and touches through it
when that reaches `gap`. The same inequality, so the same gain — `0.220059`
against `0.220059`, differing by zero. v46's *"the wave comes back inside the
circle"* and this round's *"the two pulses connect inner to outer"* are **one
event described twice**, and the check asserts it at zero rather than at a
tolerance, because that is what an identity means.

**What two waves can do that one cannot.** Below threshold nothing connects; at
threshold there is a single tangency; above it that point opens into one arc,
bounded by two genuine crossings, on which the band between the two curves
covers the *entire* radial circle — no radius at those `σ` is outside the pair.
A single wave past its own wrap threshold does nothing of the kind, and the v46
winding result is re-checked here at four gains and still holds. **Two graphs
bound a band, and a band can be radially surjective.** That is the whole
difference, and it is why the question was worth asking twice.

Two details went the other way from the guess. The antipodal refocus is a
**rarefaction**, so it is the *inward*-driven wave that bulges out to `R_outer`
— naming a wave "outward" says which way it is driven, not which way it goes.
And meeting **mid-flight**, the two travelling pulses crossing at the quarter
points, looks like the natural place for a connection and is the **worst** one:
they partially cancel, and it costs `7–9×` more than a refocus.

**What is still put in**, restated because this round is more positive than v46
and therefore easier to overread: the crossing rule is a *representation* choice
and not a derived boundary condition; the field is linear on a fixed background,
so **the two waves do not interact at all** — they are drawn on the same torus
and the question is only whether their images meet; and the gain is a *display*
amplitude. This is not a claim that two physical waves reconnect a throat. It is
that v46's obstruction does not apply to two of them, at an amplitude v46 had
already reported.

### Off the degenerate axis (the offset and the signs)

The co-located pair is the most degenerate configuration available, and saying
so is not a caveat — it is the next experiment. Both wave histories hang off one
antipodal axis, so bringing them to the same pole encourages exact overlap or
exact cancellation and tests neither. What it does not test is whether an
inner-going branch from one axis can meet an outer-going branch that has crossed
the identification and re-entered on another.

Two knobs answer it: the source separation `α`, and the radial sense each wave
is driven in, giving `δ = ε(s_A u_A − s_B u_B)`. **Opposed** signs sum the two
fields; **like** signs difference them.

The first thing measurement said was that the framing needed correcting.
Inner–inner and outer–outer are **one case, not two** — `|δ|` agrees to the bit
as a difference of fields — because flipping both signs is a reflection about
`R_mid`, an isometry of the glued radial circle. The picture cannot tell them
apart, and the reason belongs in the open: *the radial direction carries the
field's amplitude, not its direction of propagation.* A construction in which
those differ has to encode propagation in the curve, and this one does not. A
limitation, stated.

What the picture *can* distinguish carries the answer. `σ = α/2` is equidistant
from both sources, so `u_A = u_B` there identically, at every time and every
amplitude. For a like-signed pair that makes `δ ≡ 0`: the two curves are the
same curve on that axis and no gain however large carries them through the seam.
For an opposed pair the same equality makes `δ = 2εu(α/2)`, the largest it can
be. **The bisector is where one pair is maximally connected and the other is
identically not** — and there are two of them, `α/2` and `α/2 − π`, the far one
cheaper because it sits nearer the antipodal caustic.

So the answer is yes. Driven above threshold, the opposed pair's contact opens
into an arc centred on the bisector **to machine zero**, off both sources and
both antipodes, on which the like-signed pair's contact set is empty at every
offset tested. At `α = 0` the bisector collapses onto the source axis and there
is nothing off-axis to find — the degeneracy, recovered by measurement as a
coordinate fact rather than assumed.

The slider slides it: `α` moves the exclusive connection continuously from the
source axis to the quarter point and raises its price from `0.220` to `1.66`,
reached at the pulse-crossing time `t = α/2` to `0.0031π`. And **exclusive is
not cheap** — from `α = 0.125π` up, the globally cheapest connection sits exactly
on one of the four axes, is available to *both* pairs, and costs `1.7–3.7×` less.
Both numbers are reported, because the interesting one is not the cheapest.

## 25. One field on one surface, and the parity of the antipode (`docs/one_surface.md`)

§24 has a construction problem, and naming it is the round.

v66 drew **two** curves over the circle and asked whether their images meet
through the glued seam. Its own scope note said the two waves do not interact —
but that labels the problem rather than repairing it. Two curves in one frame
are two surfaces, and reading their overlap as a connection is a statement about
a picture, not about a field. The right object is one scalar deformation of one
surface, `u = s_A u_A + s_B u_B`, drawn as the single curve
`r = R_mid + ε u`, and the right question is whether *that* curve reaches
`R_outer` at one `θ` and `R_inner` at another — so the surface passes through the
identification.

**The repair costs nothing numerically and the headline everything.**
`δ = r_A − r_B = ε(s_A u_A − s_B u_B)`, the quantity v66 plotted as a separation,
*is* the one-surface deformation with the second sign flipped — the same array,
to two ulps of `R_mid`, which is the mid-radius's rounding and not the fields'.
So every v66 number survives, with the two configurations swapping names. That
inverts what §24 concluded: its cheapest-when-co-located result belongs to the
**like** pair, and its identically-zero bisector is the **node of the opposed
field**, which is both where it belongs and why it is exact.

With that fixed, the monochromatic reduction is clean and closed-form:
`u = −2A sin(mα/2) sin(mθ − ωt)`, so `B = 2A|sin(mα/2)|`. Three things follow,
all checked rather than asserted. Coincident foci with opposite orientation
cancel **exactly** — `u ≡ 0`, so no amplitude connects them and the required gain
is infinite rather than large. The optimum is `α* = π/m`, **half a wavelength**,
and the antipode is simply the `m = 1` member of that family. And `sin(mπ/2)` is
`±1` for odd `m` and `0` for even `m`, so **the antipode is parity-dependent**:
maximal for odd modes, exactly cancelling for even ones. That last is not an
artefact of the plane-wave reduction — on `S³` it appears as `Z_n(π) = (−1)ⁿ`,
and the antipodal difference field measures `2.0000` for odd `n` and `0.0000`
for even `n`.

Then the measurement disagrees with the derivation about one thing, and the
disagreement is the most useful part. A zonal harmonic is **centred**:
`Z_n(0) = 1` is a global maximum and `|Z_n| ≤ 1`, so `|Z_A − Z_B| ≤ 2` with
equality only where one focus sees `+1` and the other `−1` at the same point —
which is exactly the antipode with odd `n`. So for the real spectrum `α* = π`
for **every** odd `n`, saturating the bound, while half a wavelength reaches only
`1.41` at `n = 1` and `1.10` by `n = 5`; for even `n` the antipode cancels and
nothing reaches the bound at all. A plane wave has no distinguished centre, so
only its wavelength sets a scale; a zonal mode has one. **The parity carries
across the two models exactly. The location of the optimum does not.** The kernel
this programme cares about is `n = 1`, which is odd, so for it the antipode is
both optimal and saturating.

A second departure, and the same species: v46's field is a launched *pulse*, not
a mode — power-weighted mean `n ≈ 10`, fifteen modes for 90% of the power. Two
localized pulses cancel only while they overlap, so the `1/|sin(mα/2)|`
divergence is confined to about one pulse width and past it the threshold
**saturates** at `0.2163` rather than falling to `0.13`. The coincident
cancellation is real for a pulse too; the law governing its approach is not.

Finally the geometry, which is what the round was for. At the optimum the outward
and inward extrema sit `π/m` apart, so the chord between them is
`L = √(D² + 4 R_out R_in sin²(π/2m))` — the law of cosines regrouped so the
purely radial gap is the `Δθ → 0` limit — falling from `2.000` at `m = 1` to
`0.553` at `m = 16`, with limit `D = 0.520`. At fixed display amplitude the span
is flat at `2.0000` across the whole half-wavelength family: **the same
deformation on a progressively shorter connection.** But `E ∝ ω²A²`, so at fixed
*energy* `A ∝ 1/ω` and the span falls exactly as fast as the chord; the highest
mode that still spans the gap is `m = 6`. No favourable frequency is claimed —
that needs an energy normalisation and a packet focusing law this model does not
contain. The narrower claim is the one the visualisation has to respect: a
frequency slider cannot hold displacement fixed and then be read as
constant-energy physics.

**What is still put in** is what was always put in: the crossing rule is a
representation choice, the field is linear on a fixed background, and the gain is
a display amplitude. What changed is not the physics but the *object* — v66 asked
whether two drawn curves meet, which was never a well-posed question about a
field; this asks whether one surface reaches both boundaries, which is.

### Where each front sits on that surface (v68)

The offset question, asked inside the one-surface object rather than alongside
it. Only the **surface** is ever a closed curve in the annulus; the two
contributions are components of its deformation and appear on graphs of field
against `σ`. The annulus panels colour the single curve by which front owns each
arc — and since an inward dent is a negative contribution and an outward one
positive, each front's sign is read straight off the surface. That is the whole
answer to "where do A and B individually lie, and what are their signs".

The sweep says two things. At `α = 0` the surface is a **perfect circle** — not
a small deformation, a circle — because the contributions cancel identically.
And past one pulse width they stop overlapping *at all*: the overlap arc is
`0.000` by `α = 0.25π` and the amplification `max|u| / max|c_A|` sits at
`1.012–1.017` across the rest of the sweep, with the total peaking exactly where
a single contribution peaks.

So the honest reading is less flattering than §24's: **the offset does not turn
interference on, it turns the cancellation off.** What is left at wide offset is
two nearly independent dents in one surface. Interference here is confined to
the arc where both fronts are actually present, and for a localized pulse that
arc closes as soon as the foci clear each other — which is the same pulse-versus-
mode split that §25 found in the threshold law, seen now in the geometry rather
than in a number.

## 26. The centre as a finite bearing, not a point (`docs/regularized_center.md`)

Every section above put a **point** in the middle, and the point was doing work:
two radial arms `P_A → O → P_B` change direction at `O` for free, because at a
point there is no angular direction left to change. That is the property the
connection is wanted for — the link's cost does not care where the mouths are —
and it is bought with `f = 0`, where the geometry stops existing.

**Regularise it.** Keep `dℓ² = ds² + f(s)²dΩ²` and set `f_min = f₀ > 0`. The
middle becomes a circle in the 2-D cross-section, or the space of radial
directions `S^{d−1}` (`RP^{d−1}` for an unoriented axis). Three things §19–§23
had to force become free.

**The arms are this repo's own geometry with the symmetry dropped.** The proper
distance to an end of scale `F` is `L(F) = √(F(F−f₀)) + f₀ arcosh√(F/f₀)` and
one arm's resistance is `I₂(F) = (2/f₀)√(1−f₀/F)`; at `f_o = f_i = sin a`,
`f₀ = sin³a` they reproduce §23's `length()` and `resistance()` **bit for bit**.
So `L_o ≠ L_i` is ordinary — scanned to `437` — and the vacuole's one shared
arbitrary radial gap is gone. Two caveats the first draft missed: that is an
*intrinsic* statement, since `C¹` matching to a unit round `S³` forces
`f₀ = sin³a_j = F_j³` and so `F_o = F_i` for two genuinely matched mouths; and
`R_inner/R_outer` keeps its significance as an *endpoint scale ratio* (the very
next line uses it), what it loses is the vacuole's arbitrary drawing ratio. A
feature of angular width `Δθ` has physical width `f(s)Δθ`, so it is squeezed
into the bearing and let out, `w_i/w_o = f_i/f_o`, with `f₀` nowhere in it.

**The result is a correction to the proposal that prompted it.** Turning through
`α` looks like it should cost the bearing's arc `f₀α` — exact for the
down-turn-up route, and an honest upper bound. But the geodesic spreads the
turn, and **the reason it wins is Pythagoras, not leverage**: the corner pays
its angle as pure transverse motion (first order), the geodesic tilts motion
already radial (second order). So

> `T(α) = α²/(2I₂) + O(α⁴)` , `I₂ = ∫ds/f²` — **quadratic**, leading order, with
> the exact object the integrated `turn_cost` (`0.9248` of the quadratic at `π`).

At angular dimension `q = 2` that `I₂` is also §23's resistance, so
`T(α) = α²(4π/I₂)/(8π)`: the geometric cost of swinging the clock hands and the
electrical cost of pushing monopole flux through the tube are **one integral**
— in that dimension. The geodesic spends `1.25%` of the arc at `α = 0.1` and
`36%` at `π`; a half turn costs `8.4e-04` of the arms, and `π` is the largest
separation there is, so **no reachable orientation makes the hinge cost as much
as the journey**. **The property the point was wanted for survives, and survives
more strongly than proposed** — and the point model is the `f₀ → 0` limit, not
a rival.

**Intersection becomes overlap.** *Whether* two fronts meet on the bearing is a
question about angles with `f₀` nowhere in it; *how big* the meeting is, is
`f₀ ×` that angle. So `f₀ → 0` does **not** make everything meet — it shrinks
the overlap and the gap together, and the distinction survives as a yes/no while
vanishing as a length — or, sharper: **as `f₀ → 0` the angular incidence
survives and the physical interaction region collapses.**

**Underneath it is one Dirichlet form.** Minimising `E[φ] = ∫wφ′²ds` at fixed
increment gives `(wφ′)′ = 0`. At `φ = u` the conserved current is the monopole
flux and the answer is the conductance; at `φ = θ` it is Clairaut's `h` and the
answer is `α²/(2I₂)`. **Static monopole flux and infinitesimal throat rotation
follow the same spatial weighting** — and normalised to their own totals the two
*profiles* are the same function of position, with the deviation falling as
`α²`. But the weights are not automatically equal, and the first draft said they
were: the azimuth's is the metric coefficient `f²` in any dimension, the
monopole's is the volume element `f^q`, so they coincide **only at `q = 2`** —
the physical case here, but a fact about that dimension.

**The moment hierarchy then says what is universal, carefully.**
`T(α) = α²/(2I₂) − α⁴I₄/(8I₂⁴)`, so the **leading functional form** is
universal — `I₂` itself is not (`4/f₀` against `π/f₀`) — while **`I₄` is the
first additional independent moment**, and where the neck's shape first shows:
`1 − α²/120` scalar-flat against `1 − α²/(8π²)` hyperbolic. Not the *whole*
profile dependence: `I₆` and beyond enter at `O(α⁶)` and matter by `π`.

**And the generality is measured rather than asserted.** The law never used the
profile and holds to `1.3e-04` on an unrelated one; the `O(α⁴)` correction does
not carry over (`0.9250` scalar-flat against `0.8886` hyperbolic at `α = π`,
separating as `α² × 4.33e-03` — `4.3e-05` at `α = 0.1`, *not* the "eight
digits" the first draft claimed). The quadratic law is about necks; the deficit
at a half turn is about a particular neck.

**Scope.** Geometry only — no field equation, nothing evolving — and it does not
choose between the three candidate bulks. It works the finite bearing out far
enough to be compared with the other two, which is what it now has that they do
not: a measured hinge cost.

## 27. What higher dimension does to the picture (`docs/hyperspherical.md`)

§26 left the finite bearing with a measured hinge cost and a `f^q` scoping that
review had forced on the overlap *size* law. Following that scoping to its
conclusion changes what the picture is a picture of.

**The collapse is `fⁿ`, not `f`.** The transverse measure of an angular patch in
`ds² + f²dΩ_n²` is `fⁿdΩ_n`, so a squeeze of `1e-03` costs `1e-03` on a circle,
`1e-06` on `S²`, `1e-09` on `S³`. So: *the angular overlap can stay finite while
the physical overlap collapses as `f₀ⁿ`* — not two thick ribbons squeezing until
they touch, but **large angular structure → tiny proper measure → large angular
structure**, with the angular labels intact throughout. The yes/no criterion is
untouched; it was always angular. Which has a practical consequence for the
drawing: two fronts with finite angular overlap `ΔΩ` keep it constant all the way
down while `V_overlap ∝ f(s)ⁿ`, so the picture need not force a dramatic
macroscopic crossing — the intersection is real in the topology of the coordinate
mapping and extraordinarily small in proper measure.

**And which `n` is physical depends on which object is drawn** — `n` is the
dimension of that object's own transverse sphere, which is a fact about the
object, not a modelling choice. PR #265's spatial throat has an `S²`
cross-section: `n = 2`, its neck area the measured `4πf₀²`, its understatement
against the 2-D drawing a **thousand**. The `S³` that gives the *millionfold*
figure is a bearing in a four-spatial-dimensional embedding, a different object.
Kept explicit because otherwise the same `f^n` law migrates between objects that
do not share an `n`, and a figure derived for one gets quoted for the other.

**The finite centre is a routing manifold, not a hub.** Directional capacity is
**dimensionless** — `f₀` never enters it. At `20°` an `S³` bearing distinguishes
`113.5` directions and an `S²⁰` one `2.2e+10`, and a neck scan over six decades
returns the *identical float* while the proper measure runs `1.974e-02 →
1.974e-20`. So the singular centre is not obtained because every direction
becomes equivalent there; it is obtained because an entire finite direction space
is compressed to zero proper measure with its angular structure intact. The
`f₀ → 0` limit separates three things that had been run together: **angular
incidence survives, the overlap verdict survives, only the proper interaction
measure collapses**. It merges the routes' sizes, not their labels.

**Almost all of a sphere is at the equator of any chosen point** — the shell
measure `sin^{n−1}χ` has band width `1/√n`, measured as `std(χ)·√n → 1.000000`.
**So the antipodal relation is vanishingly non-generic:** random directions pile
up at `π/2`, and the near-antipodal fraction runs `3.2e-04` on `S²` to zero by
`n = 10`. The identification picks a measure-zero relation out of an
overwhelmingly generic alternative — which removes a bland reading of it without
making it correct.

**Orientability flips with parity, and the two quotients this arc uses are
always opposite.** `ℝP^n` is orientable iff `n` is odd. The spatial quotient
`ℝP^d` and the two-body exchange space `ℝP^{d−1}` are one apart, so at `d = 3`
the spatial `ℝP³` is orientable while the exchange `ℝP²` — where §7's Pin⁻
structure lives — is not. Raising the spatial dimension swaps them.

**And `S³` is exceptional** (`SU(2)`, parallelizable, Hopf), which is a standing
argument against reading any of this as a trend.

**Scope.** Measure, orientation and frames only; no field equation, nothing
evolving. The reading of the bearing as the blown-up embedding-centre direction
space is marked as an interpretation. And the parity section says what would
have to be re-derived *if* the dimension moved, not that it should.

## 28. What the arc cost in errors, and what caught them

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
* **An instability manufactured by a linearization.** §18 reported that its
  throat carries a growing mode for every choice of parameters, and correctly
  identified the mode's scale as the mouth's own. What it did not see is that
  the mode was *produced* by how the mouth had been written: a **constant**
  `1/(4πa)` standing in for the screened `G(a,λ)`. The exact self-energy decays
  like `e^{−κa}` and the constant does not, so the approximation fails exactly
  where the mode lives — `κa ≈ 1`, measured to `1.0004` — and the two versions
  disagree in **sign** rather than in magnitude. §19 shows the resolved mouth has
  no such mode at all, and that the true one is soft and *positive*,
  `λ₀ → 8πa/(𝒜L)`. **A quantity frozen at its leading term is a model, not a
  limit, and its failures can have the wrong sign — which is why §18's own rule,
  that a result living at the edge of its derivation is suspect, was the right
  instinct and needed one more round to cash.**
* **A parameter-range fact written down as a structural one.** §19 reported that
  the composite has **exactly one** state below the free gap, in the symmetric
  channel. True at every length that round used, and not structural: the channel
  functions have poles at `λ = (2πj/L)²` and `(π(2j−1)/L)²`, the first entering
  the gap at exactly `L = π`, and each tube harmonic that falls in brings
  another state — three at `L = 8`. The mechanism that hid it is worth naming on
  its own, because it is a counting trap rather than a physics one: **a pole is
  a sign change with no zero**, so bracketing sign changes reports states that do
  not exist. Classifying by the residual separates the two by fifteen orders of
  magnitude (`1e-15` against `1e+15`), which is a discrimination rather than a
  threshold. The same shape as the entry three above — a claim true in a regime,
  stated without the regime — but reached by counting rather than by asserting.
* **A sign on a model, mistaken for the strength of the claim.** Not an error in
  §19's *answer* — that answer survives — but in what kind of statement it was.
  That round established a sign on a reduced `2×2` and swept `3078` samples to
  support it. §20 removes the balls and the same conclusion becomes a two-line
  theorem: with nothing subtracted the energy is a sum of non-negative terms, so
  `λ > 0` for every configuration, all multipoles, no sweep. **When a sweep is
  doing the work of an argument, the model usually still contains the thing that
  made the argument unavailable** — here, the balls that were never removed.

* **A number produced by noise that looked exactly like signal.** §21's first
  answer was that `0.982` of the interference response was unreachable. It was
  **pure quadrature noise** — recomputing the same quantity with an independent
  rule gave `corr(β_A) = −0.04`. Nothing about the run looked wrong: the
  measurement ran, the figures were plausible, the residual was in the range one
  would expect a real answer to occupy. That is the point. Two uncorrelated
  noise vectors give a residual of `≈ 1`, so **the failure mode and the desired
  result are the same number**, and only an independent recomputation
  distinguishes them. Every other error in this ledger announced itself as a
  wrong value; this one announced itself as a right one. **When a measurement's
  failure mode is indistinguishable from its success, the control is not
  optional and it is not a formality.**
* **A quadrature whose own grid sat on a singularity.** The reason the above
  failed: the singular set was missing the two **mouths** and their antipodes,
  and `two_source.mouth_positions` puts the first at `(1,0,0,0)` — exactly the
  axis a product rule naturally uses. So the rule's own pole coincided with a
  `1/χ⁴` divergence, and **refining the grid made the answer worse** (magnitude
  ratios `1.44`, `1.50` between levels, growing instead of settling). The tell
  was there in the convergence table before it was understood, which is an
  argument for tabulating the ratio rather than eyeballing the result.
* **A mechanism asserted from a plausible argument rather than measured.** §21's
  first draft explained its carrier choice by "`T` is quadratic, so a carrier
  `ω₀` puts the source's power at `2ω₀`". The argument is fine in flat space and
  **wrong here**: the measured spectral peak is `5.969` for carriers `0.7`,
  `1.414`, `2.0` and `2.828` alike, because on a compact static space the field
  rings on the *background's* integer spectrum rather than on the pulse's
  carrier. Correcting it produced a better result than the one it replaced —
  `ω₃ = 2√2` is irrational, so the gravitational channel is off resonance by
  construction. **A reasonable-sounding mechanism is a hypothesis, and this arc's
  own rule is that hypotheses get measured.**
* **An argument carried across the very boundary it was about.** §21 proved the
  tensor mode `ω₃ = 2√2` incommensurate with the ESU's integer scalar spectrum
  and concluded the channel was off resonance **by construction**. Every step
  was true — *of the uncoupled ambient*. The coupled system rings where
  `det(A − Γ(ω))` vanishes, and those zeros are not integers: the nearest sits
  `0.050` from `ω₃`, `0.001` at one configuration, and never more than `0.086`
  across sixteen. **The conclusion is reversed** — generically *near*-resonant —
  and it supplies the mechanism §21 was missing for its own carrier
  sensitivity. The species: an argument established for one system and carried
  over to a system that differs **precisely in the thing being argued about**.
  Coupling is what this whole arc is about, and the spectrum is the first thing
  it moves. **Not a wrong number — an unnoticed scope on a premise**, which no
  amount of numerical care reaches.
* **A channel that cannot answer the question asked of it.** §21 measures the
  homogeneous transverse-traceless mode, and a traceless shear preserves volume
  *exactly* and area *to first order* — second-order coefficient `0.403`, and
  positive. So it distorts the mouth into an equal-area ellipse and can say
  nothing about pinching. Reading a shear amplitude as evidence about necks
  would have been a **category error rather than a measurement error**, and the
  only defence against that species is asking what the instrument is *able* to
  register before asking what it registered.
* **A check stronger than the claim it defends — three times in one round.**
  Distinct from the threshold-for-a-limit species above, and it kept recurring
  because the stronger proposition was always the more convenient one to code.
  §21 asserted that the nearest source peak sits more than `3` grid bins from
  `ω₃` (measured `2.99`); that the spectral peak *never moves* with the carrier
  (it does, at `carrier = 4`, where the envelope's DC lobe takes over); and that
  *no* coupled resonance lands near an integer (one sits `0.018` from one, while
  the mean offset is `0.264`). **In all three cases the prose claim was true and
  the assertion was not.** The load-bearing statements were narrower: that the
  mode is resolvably off every peak; that the dominant peak is the background's
  rather than the pulse's; and that the coupled resonance is *closer to `ω₃`
  than the nearest integer is* — `0.086` against `0.172`, which is the entire
  point. **Assert the proposition being defended, not a stronger one that
  happens to be easier to write.**
* **And the threshold-for-a-limit species, a fourth time.** The replacement
  check demanded that the nearest source peak sit more than `3` grid bins from
  `ω₃`; the measurement is `2.99` bins, so it failed by `0.4%`. The exact
  statement was never the bin count — it is that `2√2` is **irrational** and can
  never equal an integer, which needs no tolerance at all. §18 learned this with
  its three identifying properties, §20 with its "exactly one state", and it
  recurred here: **when a check fails narrowly, look for the exact statement it
  is a numerical shadow of** rather than retuning the number.

* **A discrepancy blamed on the wrong side.** The mouth-matching assembly agreed
  with its independent reference solve at `1e-06` and no better, at every radius
  and in both angular sectors. Flat error against a *reference* reads as the
  reference's own floor, and it was written up that way in a first draft. It was
  not: tightening the reference's quadrature tolerance, its finite-difference
  step and its boundary-value tolerance by four orders each moved the number by
  **zero**. The fault was the assembly's two-point radial stencil, whose
  *relative* truncation error on a `1/χ` field is exactly `step²` — `1e-06` at
  the step in use, which is precisely what was measured. A five-point rule
  removed four orders. **A discrepancy that refuses to move when you refine the
  other side is the other side telling you it is not the problem.**
* **A prose claim that ran ahead of its own measurement, again.** The module
  docstring asserted that the `ℓ = 1` sector "is not a refinement of this
  calculation; without it there is no calculation" — written before the
  measurement that would have tested it. Half of it is right and quantified
  (`62.5%` shortfall in the monopole-only solvability condition). The other half
  is false: the off-plane obstruction those layers carry produces an areal
  response of `6e−17`. This is the same species as §21's carrier rationale and
  §21's incommensurability argument — **a mechanism asserted from a plausible
  argument rather than measured** — and it is now three rounds running, which
  suggests the fix is procedural: *write the docstring after the control, not
  before it.*
* **An order of convergence claimed from a floor.** Once the stencil was fixed
  the reference agreement sat at `4e−10` and stopped falling with `a`. The
  measurement's summary key still asserted `both_sectors_converge` from a fitted
  order — a claim the data could no longer support, since a flat residual has no
  order. Replaced by a tolerance assertion plus an explicit statement that the
  flatness *is* the finding.

* **A singular system read out as though it meant something.** The `ℓ = 1`
  rows of the matching were a `cosh`/`sinh` transfer matrix across the tube,
  which costs a condition number of `e^{2κL}`. At the working area that is
  `e^{1.8} = 6`, so nothing about it was visible for the whole round. Asked
  for the matched-area case — the same model with **one parameter moved** — it
  becomes `e^{36} = 4.4e+15`, and the first run duly reported a condition
  number of `2.9e+15` **and an answer**, in a format indistinguishable from
  the good ones. Every digit was noise. The repair was to carry the tube's two
  end amplitudes as unknowns instead of eliminating them, which never forms
  `e^{+κL}`; the reference solves then reproduce to the *same* `4e-10`,
  digit for digit, so it was a change of form and not of content. Two things
  are worth keeping. **A model parameter moved by a factor of four hundred is
  not a perturbation of a formulation, it is a test of one** — and a condition
  number is not a diagnostic unless something actually reads it, which is why
  it is now asserted in the tests rather than merely reported.

* **A refactor validated against itself.** Rewriting the throat as an
  *admittance* was a change of form, so the test is that it reproduces the
  previous round's numbers — and the first version agreed to `1e-04` at
  `a = 0.05` and `2e-03` at `a = 0.10`. That looks like agreement, and it was
  nearly written up as one. It is not: two algebraically identical
  formulations of a well-conditioned linear system must agree to *machine
  precision*, and anything else is a bug. It was one — a sign in the
  evanescent channel's admittance. Fixed, the agreement is `1e-15`. **For a
  refactor, "close" is a failing result; only exact is passing.**
* **A tolerance that was absolute where the claim was relative, again.** The
  vacuum throat's monopole admittance is singular by an identity, and the test
  asserted `|det Y| < 1e-20`. That passes at `a = 0.05`, where `‖Y‖ ~ 4e-04`,
  and fails at `a = 0.20`, where `‖Y‖ ~ 2.5e-02` and machine precision on the
  determinant is `9e-20`. The claim was never about an absolute size; it was
  that `det` vanishes *relative to the matrix*. Third occurrence of this
  species in three rounds.
* **A label that named the wrong object.** A row of the mechanism table was
  called "the wide tube (PR #264)" while actually holding the `ℓ = 1`
  admittance fixed at the vacuum throat's — the controlled experiment, but not
  that tube. It differed from #264 by `4e-04`, which a test caught only
  because the test demanded machine precision. The controlled comparison and
  the reproduction of the earlier round are two different rows and are now
  reported as two.

* **Refining the wrong axis, and calling the survival evidence.** The
  bisector-threshold sweep showed a small turn-over at the symmetric endpoint
  `α = π`. To test whether it was sampling, the *time* grid was refined
  fourfold; the dip did not shrink, and it went into a draft as measured
  structure. It did not shrink because time was never the axis it lived on. The
  bisector was being evaluated at the nearest point of the **`σ`** grid, and at
  `α = 0.958π` the bisector falls exactly halfway between two samples. Evaluated
  off-grid at the angle it actually has, the dip is one fifth the size — and
  still real, which is the part that makes this an error rather than an
  artefact. **Refining an axis a discrepancy does not live on will always leave
  it standing, and leaving it standing is not evidence.** This is the mirror
  image of the earlier "a discrepancy that refuses to move when you refine the
  other side is the other side telling you it is not the problem": the same
  observation, and the opposite conclusion, decided entirely by which axis was
  refined.
* **An exactness asserted about the wrong quantity.** `(out,out)` and `(in,in)`
  give the same `|δ|` by an isometry, so the test asserted it at zero — and
  failed, by `2.2e-16`. The identity was never in doubt. What was in doubt was
  which quantity carried it: as a difference of *fields* the agreement is
  bit-exact, but through the *drawn radii* `(R_mid + εu_A) − (R_mid + εu_B)` and
  `(R_mid − εu_A) − (R_mid − εu_B)` round differently, and one ulp of `R_mid` is
  exactly what comes out. Both are now reported, because they are two different
  claims and only one of them is exact.
* **A claim that named too few objects.** "The cheapest connection stays on a
  source axis" was written against `A`'s axis alone and the sweep falsified it
  immediately: above `α = 0.66π` the winner alternates between `A`'s axis and
  `B`'s, the two being degenerate by symmetry and separated only by which grid
  index `argmax` reaches first. There are **four** axes in this configuration,
  not one, and the claim is true of the set and false of any member.

* **A drawing that dropped a rule the construction supplies.** The one-surface
  rig plotted the radius raw, without v46's crossing rule, so an excursion past
  a boundary stuck out past the dashed ring instead of re-entering at the other
  one. That is not a cosmetic loss: with the rule gone, two features on opposite
  sides of the gap can never meet, and the picture silently answers "no" to a
  question it was never asked. Caught by a reader driving the slider, not by any
  test — the module's numbers were all correct, because the fault was entirely
  in the renderer. **A measurement can be right while the picture of it is
  answering a different question.**
* **A symmetric field drawn asymmetrically.** `max u = −min u` to the last bit,
  but `r = R_mid + εu` in polar gives `out/mid ≠ mid/in`, so the inward
  excursion is squeezed onto a shorter arc and comes to a sharp tip while the
  outward one stays round — a `13%` effect at the drawing gain and growing.
  The repository already contained the fix and had tested it a dozen rounds
  earlier, as the `translate` versus `conformal` seam; the one-surface rig
  simply did not reach for it. **The asymmetry looked like a result and was a
  coordinate.**
* **A readout computed in the wrong place.** The seam-crossing count was
  computed inside the draw call while the readout re-derived the state from a
  fresh call that never received it, so it printed `0` at every gain — a plain
  `undefined || 0`. It looked exactly like the correct answer for a surface that
  never crosses, which is what made it worth catching rather than shrugging at.

* **A caveat mistaken for a repair.** v66's scope note said plainly that the two
  waves do not interact and that the question was only whether their *images*
  meet. That is accurate, and it is not a defence: naming an ambiguity does not
  remove it, and a reader is entitled to take a picture's central object
  seriously. The fix was not a better caveat but a different object — one field
  on one surface, where "do they meet" becomes "does it reach both boundaries",
  which is well posed. **A scope note can only bound a claim, never repair a
  construction.**
* **A factor of two in a closed form, caught by hand-checked values.** The span
  is `4A|sin(mα/2)| = 2A·(amplitude factor)`, so the required amplitude is
  `gap/(2·factor)`; the first draft wrote `gap/factor` and was exactly twice too
  large everywhere. Nothing internal to the module could catch it — every
  consumer used the same wrong function. What caught it was a test pinning four
  values (`0.130`, `0.184`, `0.260`, `0.502`) computed independently by
  bisection on the measured amplitude. **A closed form needs at least one value
  checked against something that does not use it.**
* **A dropped factor of two in the *test*, this time.** The `S³` zonal Laplacian
  is `f'' + 2 cot(χ) f'`, and the eigenvalue check wrote `cot(χ) f'` once. It
  failed cleanly with ratio exactly `1.000000` against `−n(n+2)` — the residue
  `−(n²+2n−1)` versus `−n(n+2)` — which is what an *exactly* wrong constant
  looks like as opposed to a noisy one. Worth pairing with the entry above: the
  same slip landed in module and test on the same day, and only the one in the
  module would have reached a reader.
* **A guard that fired at the wrong endpoint.** `Z_n(χ) = sin[(n+1)χ]/[(n+1)sinχ]`
  needs limits taken at `χ = 0` *and* `χ = π`, because `sin χ` vanishes at both.
  A first pass guarded only the small-χ series, so at the antipode it returned
  the Taylor polynomial: `|Z_8| ≈ 131` against a true bound of `1`, and the
  parity result it was written to test came out as "the even orders do not
  cancel". The failure was loud enough to catch immediately; the lesson is that
  it destroyed *exactly* the finding it was instrumenting, which is the way an
  endpoint bug usually presents.
* **A stall diagnosed from an instrument, not a system.** A redirected `pytest`
  run appeared frozen at `53%` for minutes at a time, so it was killed and
  restarted — three times. It was never frozen: the output was block-buffered,
  so the log only showed what had last flushed while the process ran on at
  `130%` CPU. The one genuinely slow file (`test_waves_backreaction.py`)
  dominates the suite and always had. Two of the kills also used
  `pkill -f "pytest tests/"`, whose pattern matches the killing shell's own
  command line, so a run survived and competed with its replacement for four
  cores — which then made the *next* run look slow for real. Same species as
  refining the wrong axis: **the reading came from the instrument's artefact
  rather than the system, and every action taken on it made the measurement
  worse.** The suite was finally verified in three chunks covering all 49 files.
* **A "first optimum" finder that returned the left edge.** Searching for the
  first grid point exceeding `2 − 1e-9` on a grid that never reaches it makes
  `argmax` return index `0`, so the measured optimum was the smallest sampled
  `α` and the error read `2.09`. Fixed by searching one period of
  `|sin(mα/2)|` rather than thresholding. **`argmax` on an all-false mask is a
  silent zero, and it looks like an answer.**

* **A production quantity built out of a cancellation.** `VacuumThroat.
  admittance` formed the mouth-to-mouth cross term as `½(s − t)` from the two
  eigenchannel values of a Riccati solve. Both are `O(1)` once `ℓ ≥ 1` and they
  agree to more and more digits as the neck closes, so the answer was made
  entirely of the cancelling tail. It survived every check it had, because the
  only channels anything consumed — `ℓ = 0` and `ℓ = 1` — happen to sit above
  the solver's floor. At `ℓ = 2, a = 0.05` it returns `−1.17e-14` for a true
  `+3.00e-16`: **wrong sign, and larger than the answer.** Caught by deriving
  the closed form, which reaches the same number as a *product* of small
  factors. **A difference of two numbers cannot carry information the numbers
  have already lost** — and a formulation with no cancellation in it is not a
  refinement of one that has, it is a different instrument.
* **The same failure mode, in the verification of the repair.** Checking the
  closed form, the first pass computed `C₂` as `½(ν_s − ν_a)` from two `O(1.26)`
  values and got `3.4879e−16` against the `sinh`-form's `3.0011e−16` — `14%`
  out. That is the *original* bug, reproduced by hand while auditing it, and it
  is worth recording because the wrong route looked like the natural one both
  times.
* **A retraction that never reached the repo.** An interpretation offered in
  discussion — that the low modes carry a `400×` WKB amplitude enhancement —
  was wrong and is withdrawn. It was never committed, so there is nothing in
  the code or docs to correct; it is recorded here because a claim made out
  loud and then quietly dropped is the one that comes back. The `400×` that
  *does* appear in §22 is unrelated and stands: it is the working tube's area
  against its own mouths', `1/sin²a` at `a = 0.05`.
* **An identity and its two corollaries, run together as one claim.** `H = 0`
  with `K_ij = 0` gives `θ_± = 0`, so the neck is a marginal sphere of the
  slice. That much is exact. §23 then wrote it as *"the neck is an apparent
  horizon, and therefore the throat is not traversable"* — two further steps,
  neither of them checked. An apparent horizon is the **outermost** MOTS, a
  global condition on the slice; non-traversability is a property of the
  **Lorentzian development**, which needs a lapse the round never chose. The
  second does follow if the development is taken to be the standard vacuum
  Schwarzschild/ER one, which is a natural reading and still an *assumption*.
  Nothing numerical was wrong — the identity holds to the last bit — and that
  is what let it pass: **an exact result carried two unexamined inferences
  downstream because they were written in the same sentence as the thing that
  was proved.**
* **Two claims stated wider than what was shown.** `f_min > 0` was written as
  though forced by Einstein's equations; it is forced within the class actually
  worked in — spherically symmetric, scalar-flat, `K_ij = 0`, `C¹`-matched.
  And `physical_throat` was discussed as if it posed a dynamical problem; it
  supplies **spatial initial data only**, and the dynamic problem is not
  well-posed until a lapse is chosen. Both were corrected in review, and the
  pattern is the same one this section keeps finding: **the calculation was
  right and the sentence around it was not.**

* **The cancellation lesson, ignored one round after writing it down.** §26's
  first draft of the turn cost computed `T(α)` as `route_length − (L_o + L_i)`
  — a `1e-15` quantity from two `O(1)` lengths. That is precisely the failure
  mode the entry three above describes, reproduced in a different module the
  week after it was diagnosed and repaired in `physical_throat`. Answers were
  `2.7×` wrong at `f₀ = 1e-07`. **Writing a failure mode into a ledger does not
  inoculate against it**; what caught this one was a test demanding that the
  shape `T/(α²/2I)` be one number across decades of `f₀`, which no incorrect
  formulation can satisfy.
* **And the obvious repair was worse than the bug.** Rewriting `T` as a direct
  integral — no subtraction anywhere — was *less* accurate at `f₀ = 1e-07`,
  because integrating in `f` puts a spike of width `√f₀` in front of an
  adaptive quadrature, which walks straight past it. Two independent defects
  with opposite signatures, and fixing only the one that had a name would have
  left a number that looked better and was not. The fix was the *variable*:
  `f = f₀cosh²x`, the same substitution that made PR #267's admittance a closed
  form, after which the result is stable over seven decades. **"Remove the
  cancellation" is not the lesson; "check the answer is independent of
  something it must be independent of" is.**

* **Six overstatements in one round, none of them numerical.** §26 shipped with
  a correct set of geodesics and a prose layer that outran them, and review
  caught all of it: an `O(α⁴)` law written as an equality; a "break-even angle"
  of `104` rad extrapolated outside both its domain and the configuration space
  (`π` is the largest separation there is); asymmetric arms called a matched
  #265 throat when matching to a unit round `S³` forces `F_o = F_i`; a scale
  ratio declared meaningless one paragraph before being used; a Dirichlet
  identity stated dimension-free when the monopole weight is `f^q` and only
  `q = 2` matches the hinge's `f²`; and `I₂` called universal when it is the
  *functional form* that is. **The calculations were right in every case.**
  What the round lacked was a pass asking, of each sentence, whether the
  measurement underneath it actually says that.
* **And one that was backwards rather than merely wide.** The explanation of
  *why* the geodesic beats the corner — "it turns where the lever arm `f` is
  longer, so a given angle costs less arc" — is wrong twice over: an angular
  increment at larger `f` costs *more* arc (`f dθ`), and `θ' = h/f²` puts the
  turn where `f` is *smallest* (76% of it inside `f < 2.4 f₀`). The real reason
  is Pythagoras: the corner pays its angle as pure transverse motion, first
  order; the geodesic tilts motion already radial, second order. **A correct
  number with a fabricated mechanism attached is the failure mode that survives
  every numerical check there is** — the tests all passed, because none of them
  tested the sentence.
* **A number reported for the wrong quantity.** "The two profiles agree to eight
  digits at `α = 0.1`" — `8` digits is how well *each profile matches its own
  quartic law* (`7.9e-09`); the two match *each other* only to `4.3e-05`. Two
  quantities computed in adjacent lines, and the smaller was quoted for the
  larger. It propagated to six files before review caught it.

* **A maximiser reported at the edge of its own search range.** §27's first
  scan of where the ball volume peaks looked over `d < 80` and reported `79` for
  `R = 4`. That is not a peak, it is a range: widened, the answer is `100`. The
  same species as the `argmax`-on-an-all-false-mask entry above — **a maximiser
  that returns its own boundary looks exactly like an answer** — and the scan
  now flags any peak that touches the ceiling.
* **And a tie resolved by floating point.** In the same scan the sphere measure
  at `R = ½` is *exactly* tied between `d = 2` and `d = 3`, since `2πR = 4πR²`
  there. A bare `max` picked a side, and the side it picked *changed* when the
  computation moved into log space — which is how it was noticed. Third
  appearance of the maximiser-degeneracy failure mode in this arc, after the
  even-`n` harmonic peak and the `argmax` mask. Ties are now reported rather
  than resolved.
* **An exponent quoted without its object.** §27's headline "the 2-D drawing
  understates it by a million" reached three files — the doc, the README and the
  renderer — with no statement of *which object* has `n = 3`. The number is
  right for an `S³` bearing and wrong by three decades for the `S²` throat of
  PR #265, which is the object the rest of the repository means by "throat".
  Nothing numerical was wrong; a law was stated without the scope that makes it
  a law about anything. **A dimensionless exponent will migrate between objects
  unless each statement of it names the object it belongs to** — the catalogue
  is now measured (`measure_which_n_is_physical_for_which_object`, checked
  against `physical_throat`'s own `4πf₀²` neck area) and a test fails if any
  file quotes the figure without naming the object within six lines.
* **An overdetermined boundary, caught by the constraint and nothing else.** The
  first evolution froze `φ(v, R)` at the outer edge, which looks like the obvious
  "no incoming radiation" condition. It is not available: the characteristic
  hierarchy spends its three integration constants on the gauge and on central
  regularity, so `ψ(v, R)` is already determined and imposing anything on it
  overdetermines the slice. Every diagnostic that looked at the *solution* was
  clean — the field stayed smooth, the amplitude was stable to four digits across
  resolutions — and the `vv` residual at the outer edge sat at `O(1)` and got
  **worse** under refinement, diverging like `h⁻²`. **A residual that diverges
  under refinement is not noise, it is a boundary condition that should not be
  there**, and it was the only instrument in the round that could see it.
* **Two converged numbers that disagree.** Two independent horizon-penetrating
  ringdown codes, each stable and each flat to four digits under a four-fold
  refinement and across three extraction windows, give damping rates `37%` apart.
  There is no resolution study that resolves this, because both *are* converged.
  The response was to report no frequency and name the first thing to chase —
  the Kerr–Schild operator fails a flat-space exactness test at its inner cut,
  error flat at `1.07e-02` across a four-fold refinement, which points at the
  excision boundary rather than at the operator. **Convergence establishes that a
  code computes something well; it says nothing about what.**
* **A wrong operator that no test protected.** The radial scalar potential was
  short of the minimally coupled scalar master form by `3A²/(4r²)` for the whole
  life of the module. Flipping it to the correct operator broke **2 tests out of
  1582**, both of them the previous round's own bookkeeping. The `γ` sums, the
  `R_OUTER` fixed point and the `1.054` factor — all load-bearing — are not
  regression-locked anywhere in the suite; they live in prose and probe
  narratives. **A green suite is evidence that the code does what the tests say,
  not that the tests say anything about the physics.** The bug was caught by a
  derivation done for an unrelated round, and nothing in the tree would have
  caught it otherwise.
* **A sub-percent residual that was not small.** The `γ` lock was reported as a
  `−2.2%` or `−0.21%` near-equality, which reads like a rounding-level detail.
  Passing three geometries through the locked lepton chain gives a secant
  `d ln m_μ / d ln γ = −16.6`, and a local derivative at the lock of `−17.5`:
  the corrected `−0.75%` residual is a **measured `15.2%`** muon error. **The
  size of a residual means nothing until its elasticity is measured** — and the
  elasticity had never been measured, so a number quoted as reassuring was
  actually the most sensitive input in the chain.
* **The same residual is worth different amounts in different sectors.** The
  quark ladder reads the *same* barrier sum, and its elasticity is `+4.8` —
  3.6× softer — so a comparable geometric agreement is a `15%` miss in one
  sector and a `1.8%` miss in the other. The follow-on lesson to the one above:
  elasticity is a property of the *consumer*, not of the geometry, and must be
  measured once per chain rather than inherited.
* **Three residuals individually right can be jointly wrong.** Each quark
  residual sits within `0.7%` of its lock, and substituting all three at once
  misses the ladder by `3.8%` — worse than the legacy set, whose larger
  individual errors partially cancel. **Per-knob agreement is not system
  agreement**, and only the per-knob version had ever been checked.

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

And the dynamics round supplies the sharpest instance of the standing lesson yet,
because this time **both** numbers converged. Two horizon-penetrating time-domain
constructions of the same ringdown — a Kerr–Schild slicing and a tortoise
evolution — are each stable under refinement and each flat to four digits, and
they disagree on the damping rate by `37%`. Neither converged more; both
converged, to different answers. That is the cleanest possible statement of *a
converged number is not a correct number*, and the only honest response was to
report no frequency at all rather than pick the one that looked better.

## 29. The first evolved Einstein equations (`docs/tangherlini_dynamics.md`)

Everything §1–§27 did with gravity was stationary, weak-field or linearized. This
round evolves the Einstein equations in time for the first time in the tree, at
the highest-symmetry 4+1 problem: `D = 5`, spherical symmetry, one minimally
coupled massless scalar, in horizon-penetrating ingoing Eddington–Finkelstein
coordinates. Vacuum is not an option — Birkhoff in `D` dimensions makes
Tangherlini unique, so the scalar is the dynamical content.

**The system is derived for general `n` and self-checks at `n = 2`**, the known
`D = 4` case. The `vr` equation turns out to be an exact **quadrature**,
`(r^{n−1}e^δA)' = (n−1)r^{n−2}e^δ`, not an ODE — which is what makes the geometry
machine-precise on Tangherlini (`1.6e-15`).

**The Einstein equation the code never solves converges at second order**,
`1.989 → 1.997 → 1.999`. The hierarchy solves `rr` and `vr` on every slice, so
their residuals are identically zero and testing them would be circular; `vv` is
the one component left over and it carries `∂_v A`, which the code never
otherwise forms. Named as the characteristic-scheme *analogue* of a
Hamiltonian/momentum constraint test, not as one.

**A regular centre forbids a trapped surface.** The quadrature's right-hand side
is a positive integrand over a positive interval, so `A > 0` strictly for `r > 0`
identically, and no trapped surface can sit on a regular-centred ingoing null
slice. Four profile families pushed to `min A = 5.63e-03` never cross. So horizon
*formation* is not observable in this gauge, and the criterion is the loss of
central regularity rather than `A` changing sign — a statement about the chart,
not the physics, and the reason production characteristic codes use outgoing
cones or excise the centre.

**Two targets missed and reported as missed:** the perturbation spectrum and the
retarded outer→inner transfer function. See the entry added to §28.

**And one discrepancy found in passing.** The master potential for a minimally
coupled massless scalar with `ψ = r^{n/2}φ` differs from
`tangherlini.radial.V_tangherlini` by exactly `3A²/(4r²)`, with the flat-limit
Bessel form settling which is which to `4.3e-16`. Nothing was changed —
`V_tangherlini` is consumed by roughly fifty probes and several derived
constants, so acting on it is a decision about the repository's published
numbers, not a side effect of a dynamics round.

## 30. The radial scalar operator, corrected (`docs/scalar_operator_audit.md`)

§29 discovered that `tangherlini.radial.V_tangherlini` was short of the minimally
coupled massless scalar master potential by an `ℓ`-independent `3A²/(4r²)`, and
deliberately changed nothing. This round corrects it and audits what moved.

**The correction is exact and triple-confirmed** — the gap matches its closed
form to `2.4e-15`, carries no `ℓ` to `2.4e-15`, and the flat limit reproduces the
Bessel form to `2.2e-16` while the legacy operator fails that same test.

**Three verdicts, not two.** The cross-`ℓ` operator `V_{ℓ+2} − V_ℓ` is
*exactly invariant* (`3.6e-15`); eigenvalues, actions, ratios and matrix elements
are *numerically shifted*; and the claim that the `ℓ = 0` channel closes the `γ`
gap is *interpretation changed* — withdrawn, not replaced.

**And the narrow re-derivation inverted its own question.** Passing three
geometries through the locked lepton Hamiltonian showed **B and C bit-identical**,
because the locked block discards `r_outer` and sees only the scalar `γ`. So
`γ = 22.5` is the selector and `R_OUTER` is downstream of it; fixing `R_OUTER`
breaks the ladder at the 15–21% level where the legacy geometry missed by 3.8%.
The correction therefore *weakens* the geometry-supplies-`γ` story even while
improving the `1..5` residual, because `d ln m_μ / d ln γ = −17.5` at the lock
(secant `−16.6` across the two corrected geometries — two quantities, not one)
makes a sub-percent geometric residual anything but small.

## 31. The quark residual sector, re-derived (`docs/quark_residual_reaudit.md`)

§30 classified the downstream chains and deliberately re-derived none of them.
This round takes the quark residual sector — three knobs of the frozen v3 lock,
all derived from the same eigensolver — and **reverses §30's expectation**.

**All three residuals move toward their locked values**: `pinhole −1.09% →
+0.36%` (the residual changes sign), `transport +0.88% → +0.70%`, `resistance
+0.49% → −0.02%`. The lepton sector's closure broke; the quark sector's
improves.

**And the reason is not geometric.** One barrier feeds both sectors and the
correction moves it once — `Σ V_max[1..5]` at `R = 1.26` goes `22.008 →
22.331`, and that single number *is* the lepton `γ` and *is* the quark
`pinhole`. What differs is how hard each Hamiltonian leans on it:
`d ln m_μ/d ln γ = −17.5` against `d ln m_s/d ln pinhole = +4.8`, a factor of
3.6. The lepton's `−0.75%` residual is a measured `15.2%` muon error; the
quark's `+0.36%` is a measured `1.79%` strange error.

> **A percentage agreement between a geometric quantity and a fitted knob means
> nothing until multiplied by the elasticity of what it feeds.**

**One claim gets weaker on inspection, under both operators.** Substituting all
three derived residuals with nothing retuned gives `3.44%` (legacy) and `3.78%`
(corrected) against the fitted lock's `1.61%` — and the *ordering reverses*,
because the legacy triple's errors partially cancel. "Each knob is derived to
within 1%" and "the derived knobs reproduce the ladder" are different claims and
only the first was ever established. The correction did not create that gap; it
exposed it.

**One thing genuinely improves.** Read as demands on `R_OUTER`, the legacy
operator has both sectors asking for more than `1.26`, putting the canonical
value outside the bracket they define; the corrected operator has them
**straddle** it (`1.25645` quark, `1.26788` lepton). Different evidence from the
single-sector fixed point §30 reopened — and weak, since a `0.91%` window admits
anything inside it. Recorded as suggestive, not as a restored derivation.

## 32. What is imported rather than derived

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

## 32. What would come next

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
            → mouth resolved ✓ → neck resolved ✓ → backreaction ✓
            → signed ΔA/A ✓ → which throat ✓
            → stationary action → topological branch
```

**The gate is answered, and the reason it was there is worth keeping.** §18 built the
conservative finite throat the arc had owed since §11 — and found that with
point mouths it carries an exponentially growing mode for *every* choice of
parameters, at a rate that knows neither the tube's length nor the mouths'
separation. Stationary action and backreaction are both integrals over a solved
field on this background. Run on a background with a growing mode, each would be
measuring the mode: an on-shell action evaluated on a solution whose amplitude
diverges is not stationary in any useful sense, and an `A`/`B`/`A+B` collapse
comparison would report the instability's response rather than the waves'. §19 answered it: **the mode does not survive**, it was an artifact of
linearizing a screened self-energy, and the delay and conservation law carry
over. §20 then closed the gap that answer stood on — §19's mouths were spheres
in a *fixed* ambient with monopole coupling, so the balls were never removed —
and with them removed the same conclusion becomes a **theorem**: nothing is
subtracted, so the energy is a sum of non-negative terms and `λ > 0` for every
configuration, all multipoles, no sweep. So the two steps below resume, with an
object that has a proper length, a delay, a conservation law and positivity
that is proved rather than sampled. What is still owed to them is the neck's
**cross-section** — the tube carries one transverse channel, so `𝒜` is a
coupling rather than a solved geometry — and the ambient metric, which is fixed.
Both are limitations with numbers attached rather than gates.

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
* **stationary action** — *ungated by §19, secured by §20, and now the next step.* Evaluate the on-shell
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
* **backreaction** — ~~done~~ (§21), to first order and in one channel. The
  first GR question is not "does spacetime pinch off?" but whether `A + B`
  produces a response not reproducible by rescaling `A` or `B` alone, and the
  answer is **yes**: `0.921` of the interference response is unreachable at the
  working point, with the interference term *comparable in size* to the
  single-wave ones. §17's warning about *which* diagnostic to integrate — `ΔT`,
  not `T_A:T_B` — turned out to be the whole structure of the measurement, since
  `β[ΔT]` is exactly the object the reachable family fails to contain. The
  §18 **blocker** was cleared first: §19 resolved the mouth and §20 removed the
  balls, and the second proves positivity rather than sampling it. What §21 does
  *not* do is close the loop — it is a linear response on a fixed background in
  the lowest TT harmonic, with point sources and §18's point throat, so a
  **solved coupled system** is still owed;
* **topological branch** — the detached resonator, last, and only if
  backreaction produces a finite-radius neck.
