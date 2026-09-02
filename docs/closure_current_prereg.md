# Pre-registration: positive coarea count or holonomy-weighted coarea current?

**Status: frozen before any code.** Sixth round of the finite-mouth chain,
following `docs/closure_measurement_dependence_prereg.md` (`1b0144e`) and
`5ecfed4`, whose interpretation this round may reverse. Success criteria may
change afterwards only by an explicit correction note.

## 0. Four corrections to `5ecfed4`, recorded first

**C1 — the signed density is not an inserted quasi-probability.** On the
closure set `N = 0`, `Ω = 2 atan2(0, D)` is `0` for `D > 0` and `±2π` for
`D < 0`, so the spin-½ closure holonomy is

```
e^{iΩ/2} = sgn D          (away from the measure-zero points D = 0)
```

and the "signed variant" `D/(2|u×v|)` is exactly `e^{iΩ/2} × (coarea density)`:
the coarea measure weighted by the actual holonomy of the `0`/`π` closure
branch. The sentence "adopting `D` imports the quantum rule" is withdrawn.
What is unproved is whether physical histories combine with their closure
holonomy or as positive counts. That is this round's question.

**C2 — scope.** `history/closure.py` states four conditions — topological
closure, conservation, phase closure, stationarity — and its `total_phase`
sums event and worldline phases plus `π/2` per throat transport; the branch
`0`/`1` is recorded but the Gaussian weight discards its sign; stationarity
is named in the docstring and implemented nowhere. `5ecfed4` used the
phase-closure component alone, on a geodesic-realignment loop. It is to be
described as *the phase-closure component of BAM under the
geodesic-realignment loop ansatz*, not as "the repository's axiom with
nothing imported".

**C3 — the sector prior.** `closure_weights` normalises the four coarea
volumes with equal prior weight `π_{s_A s_B}` on the four disconnected
outcome sectors. Coarea fixes the measure within a sector, not across
sectors. `π₊₊ = π₋₋`, `π₊₋ = π₋₊` (from the `(x,u,v) ↦ (−x,−u,−v)`
symmetry) keeps the marginals at `1/2` and still lets `E` move. The equal
prior is the counting measure on sectors: **chosen**, entered in the ledger.

**C4 — no-signalling, qualified.** Proved: detector no-signalling,
`P(A|a,b) = P(A|a) = 1/2`. Not proved: operational no-signalling to the
past. `ρ(x|a,b)` depends on the future settings, maximally for
non-coplanar pairs, so a source observer able to read `x` would learn the
later analyzer plane. The model has no source read-out coupling; whether
the dynamics makes `x` inaccessible is open and is an explicit falsifier
below.

## 1. Results established before freezing (R1, R2)

* **R1 — the loop is derived, not inserted.** In the quaternion model of
  the previous rounds, geodesic transports are right multiplications by the
  minimal-rotation lifts `g(p→q) = cos(θ/2) + sin(θ/2) (p×q)/|p×q|`, and
  the mouth transport is `J = L_{−j}`. For the history
  `source → A → J → B → source`,

  ```
  q₅ = J · [ J · q₀ g₁⁻¹ ] g₂⁻¹ g₃⁻¹ = −q₀ (g₃g₂g₁)⁻¹
  ```

  because left and right multiplications commute and `J² = −1`. So the two
  throat transports contribute the spinor sign and nothing else; the
  holonomy is fibre-independent; and with `R(−u→v) = R(u→−v)`,
  `R(v→−x) = R(−v→x)`, the frame holonomy is that of the geodesic triangle
  `x → u → −v → x`: a rotation about `x` by `Ω(x; u, −v)` (checked to
  `5e-15`), with SU(2) lift `G = cos(Ω/2) + sin(Ω/2) x` exactly, sign
  included (`2e-15`). On the closure circle `G = sgn D(x, u, −v)`
  (`2e-14`). **The Pin-derived loop is the triangle of `5ecfed4` with the
  partner sign `v → −v` (singlet) and the closure holonomy `−sgn D`.** The
  triangle was the right object.
* **R2 — the holonomy-weighted coarea current is the quantum joint law,
  without projectors.** With `W̃(s) ∝ ∫_Γ D_s dσ / (2|u×v|)` and
  `∫_Γ x dσ = 0`, `∫_Γ D_s dσ = L_Γ (1 + u·v)`, so after normalisation

  ```
  P₊₊ = P₋₋ = (1 + cos γ)/4 = ½ cos²(γ/2),    P₊₋ = P₋₊ = ½ sin²(γ/2),
  E = cos γ  (triplet sign)  →  −cos γ with the partner sign of R1,   S_max = 2√2.
  ```

  The Bargmann identity `D/4 = Re Tr(P_x P_u P_v)` is then a representation
  of what is already in the closure geometry, not its source.

## 2. The fork, with the rule fixed in advance

Two candidate measures on the closure set, both parameter-free:

* `POSITIVE_COAREA` — histories are ordinary classical alternatives and
  closure conditions them with `1/|∇Φ|`: `E(γ) = [c²(π−γ) − s²γ]/[2 sin γ + c²(π−γ) + s²γ]`,
  `S_max = 2.1423`.
* `HOLONOMY_WEIGHTED_COAREA` — the closure current retains its derived
  spin holonomy `e^{iΦ}`, so the `0` and `π` branches enter with opposite
  orientation: `E = ∓cos γ`, `S_max = 2√2`.
* `NEITHER` — a fuller dynamics (field, action, stationarity) selects a
  different measure.
* `FORK_UNRESOLVED_BY_CURRENT_STRUCTURES` — permitted as the narrower
  verdict: neither the orientation/Pin structure nor any implemented
  classical dynamics in the repository selects between the two.

**The rule:** `HOLONOMY_WEIGHTED_COAREA` may be adopted only if the
orientation/Pin structure or a stationary classical field equation is
shown to *require* the oriented sum. It may not be adopted because it gives
`2√2`.

## 3. Targets, with falsifiers

* **R3 — the sector prior.** Compute `E` under `π_like/π_unlike = r` for
  `r ∈ {0.5, 1, 2}`: marginals stay `1/2`, `E` moves (pre-computed at
  `γ = 1` from the closed-form weights `W_like = 2 + c(π−γ)/s`,
  `W_unlike = 2 + sγ/c`: `E = (rW_l − W_u)/(rW_l + W_u)`, i.e.
  `0.0947, 0.3985, 0.5695`). No symmetry of the model relates like and
  unlike sectors at fixed settings (`v → −v` exchanges them but changes
  `b`); the equal prior is chosen. *Falsifier:* a symmetry that fixes `r`.
* **R4 — stationarity.** The repository does not implement it (grep). The
  stationary set of the closure phase `Ω(x)` on `S²` (`∇Ω = 0`) is
  computed and shown to be distinct from the closure set `N = 0` (isolated
  points versus a circle), so imposing both would select finitely many
  source directions and no continuous measure. Stationarity as stated in
  the repository does **not** decide the fork. *Falsifier:* the stationary
  set coincides with the closure set.
* **R5 — what the Pin structure does and does not force.** The Pin
  transport supplies the branch holonomy `−sgn D` geometrically (R1); a
  holonomy is a *label* on each closed history. Whether the label enters
  the measure as a weight is not a property of the bundle. Prediction:
  no computation available in the repository's classical structures
  selects a weighting; the choice is the quantization step, now localised
  to "positive count versus oriented sum over the `0`/`π` closure
  branches". *Falsifier:* a derivation, from the Pin structure or the
  classical equations, that the closure current is oriented.
* **R6 — operational no-signalling to the past.** The model as defined
  has no source coupling and cannot answer it. Recorded as the open
  falsifier for a dynamical round: a source read-out whose statistics
  depend on `(a, b)` would falsify the model as a physical theory.

## Correction notes (post-implementation, from review)

**First correction — R4 substituted the wrong notion of stationarity, and the
result is stronger than predicted.** The repository's fourth closure
condition is *stationarity: the history has extremal action*. No action
functional exists in `history/closure.py`, so replacing that condition with
stationarity of the geometric phase, `∇Ω = 0`, introduces a **new
assumption**; it is now labelled a **proxy** and does not test the
repository's condition. Analytically, with `N = x·q`, `q = u×v`,
`D = A + x·p`, `A = 1 + u·v`, `p = u+v`:

```
∇Ω = 2(D ∇N − N ∇D)/(D² + N²),     and on N = 0,     ∇_{S²}Ω = 2q/D = 2(u×v)/D
```

(`q` is already tangent there), of norm `2|u×v|/|D| > 0`. So for
non-collinear analyzers, **sharp closure and phase-stationarity are disjoint**,
not intersecting in isolated points as R4 predicted; the only points with
`D = 0` are `x = −u` and `x = −v`, where `D = N = 0` and the `arg` chart is
singular — a chart singularity, not stationarity. Verified by a
branch-wrapped finite difference to `3.4e-7` (`Ω` itself jumps by `4π`
across the closure circle where `D < 0`, while the holonomy
`e^{iΩ/2} = −1` is continuous). This strengthens the fork diagnosis: **no
variational principle available in the repository can choose between
positive and holonomy-weighted coarea.** The earlier grid-search version of
R4 and its Lexell level-set check are superseded by this identity.

**Second correction — R3's pre-computed values.** `0.0947, 0.3985, 0.5695`
were mis-evaluated when the document was written; the closed form
`E = (rW_l − W_u)/(rW_l + W_u)` with `W_l = 2 + c(π−γ)/s = 5.9201590573`,
`W_u = 2 + sγ/c = 2.5463024898` at `γ = 1` gives
`0.075145, 0.398497, 0.646018`. The criterion (marginals fixed at `1/2`,
`E` moves) is unchanged and holds.

**Addendum, not a criterion — the oriented current has non-negative sector
integrals.** Although `D` changes sign pointwise, `∫_Γ x dσ = 0` gives
`∫_Γ D_s dσ = 2π(1 + u·v) ≥ 0` for every outcome pair. The
holonomy-weighted construction therefore uses destructive cancellation
*internally* and yields non-negative normalised sector weights: the
analogy is classical wave interference, not a negative-probability
distribution. Recorded, and it does not decide the fork.

**Addendum, not a criterion — the sharper form of the fork.** Whether the
Pin/Hopf data make the closure locus naturally an *oriented current with
local coefficients*, whose observables are integrals of sections (in which
case the sign is geometrically mandatory), or whether the physical object
is a *measure on histories* (in which case positivity forces `|D|`). Added
to the implementation audit as an open direction; establishing it would
resolve the fork under the rule of §2.

**Third correction (review of PR #280) — the headline was one item too
narrow, and R2 tested the wrong object.** Four findings, none changing a
criterion of §2, all narrowing what may be claimed:

1. *The equal outcome-sector prior is load-bearing on the oriented branch
   too.* R3 recorded it only for the positive branch. The oriented sector
   integrals are proportional to `1 ± s_A s_B cos γ`, so under
   `r = π_like/π_unlike` the correlation is
   `E_r = [r(1+cos γ) − (1−cos γ)]/[r(1+cos γ) + (1−cos γ)]` (triplet) and
   its analogue for the singlet, equal to `±cos γ` **only at `r = 1`**:
   `0.25243, 0.54030, 0.74031` and `−0.74031, −0.54030, −0.25243` at
   `γ = 1` for `r = 0.5, 1, 2`, with marginals `1/2` throughout. So the
   fork is **not** a single binary `|D|` versus `D` choice; the headline is
   narrowed to *branch aggregation + sector coefficients + readout remain
   underived*.
2. *R2 did not test the derived loop.* It computed `closure_weights(γ,
   "signed")`, the triplet `D(x, u, v)`, and reached the singlet by a verbal
   sign substitution. The derived object `x → u → −v → x` with current
   `−D(x, u, −v)/(2|u×v|)` is now integrated directly:
   `P(s_A,s_B) = (1 − s_A s_B cos γ)/4`, `E = −cos γ` to `≤ 1e-17`.
3. *"The loop is derived" was too strong.* Derived: the quaternionic
   reduction of the **chosen** itinerary, under the **chosen** geodesic
   realignment, to the triangle, plus its holonomy. The itinerary and the
   realignment remain model choices, as the dependency ledger already said.
4. *Non-negative sector integrals are not yet an event-frequency law.* The
   wave analogy is qualified: a classical detector normally responds to a
   quadratic functional of an amplitude, so why frequencies would be
   **linear** in the integrated current is a third open item, and forcing
   the sign by local coefficients would not by itself complete the
   derivation.

## 4. Expected verdict, stated in advance

`FORK_UNRESOLVED_BY_CURRENT_STRUCTURES`, with the two candidates' predictions
recorded side by side and the quantization gap located at the weighting of
the closure branches. If R4 or R5 produces a derivation, the verdict changes
accordingly; the rule of §2 governs.

## 5. Dependency ledger to be printed

```
loop                =  loop( source pair [x, −x via J], geodesic realignment at A, B [chosen],
                             mouth transport J [derived: L_{−j}] )  →  triangle x → u → −v → x [derived]
closure holonomy    =  −sgn D(x, u, −v)  [derived from J² = −1 and the lift G]
positive coarea     =  |D|/(2|u×v|)      [derived conditioning; counting measure on sectors: chosen]
oriented current    =  e^{iΦ} × coarea  [derived label; adopting it as weight: the open step]
sector prior        =  counting          [chosen]
```
