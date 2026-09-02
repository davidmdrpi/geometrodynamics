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
