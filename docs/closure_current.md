# Positive coarea count or holonomy-weighted coarea current?

*Sixth round of the finite-mouth chain. Pre-registered in
`docs/closure_current_prereg.md` (`f954e3d`), which also records four
corrections to `5ecfed4`; module `geometrodynamics/bulk/closure_current.py`;
tests `tests/test_closure_current.py`; probe `closure_current_probe.py`.*

## Verdict

`FORK_UNRESOLVED_BY_CURRENT_STRUCTURES`

Two candidate measures on the closure set, both parameter-free, differ only
in whether the two closure branches (`Ω/2 = 0` and `Ω/2 = π`) are counted or
summed with their holonomy:

| candidate | density on `Γ` | `E(γ)` | `S_max` |
|--|--|--|--|
| `POSITIVE_COAREA` | `|D| / (2|u×v|)` | `[c²(π−γ) − s²γ]/[2 sin γ + c²(π−γ) + s²γ]` | `2.1423` |
| `HOLONOMY_WEIGHTED_COAREA` | `e^{iΩ/2}|D| / (2|u×v|) = D/(2|u×v|)` | `∓cos γ` | `2√2` |

Nothing classical in the repository selects between them, and the
pre-registered rule forbids selecting the second because it gives `2√2`.

## 1. The corrections to `5ecfed4`, now incorporated

* **C1.** On the closure set `e^{iΩ/2} = sgn D` (checked to `1e-16`), so the
  "signed" density is the coarea density times the branch holonomy, not an
  inserted quasi-probability. The sentence "adopting `D` imports the quantum
  rule" is withdrawn.
* **C2.** The previous round used the phase-closure component of BAM under a
  geodesic-realignment loop ansatz. `history/closure.py` names stationarity
  and implements it nowhere (identifier audit), records the branch and
  discards its sign in the Gaussian weight.
* **C3.** The equal prior on the four outcome sectors is the counting
  measure on sectors: chosen. With `π_like/π_unlike = r ∈ {0.5, 1, 2}` the
  marginals stay `1/2` and `E(1) = 0.075, 0.398, 0.646` (the pre-registration's `0.095, 0.570` were mis-evaluated; corrected there). No symmetry of the
  model at fixed settings relates like and unlike sectors.
* **C4.** Detector no-signalling is proved; operational no-signalling to the
  past is open (the model has no source read-out).

## 2. What is now derived (R1, R2)

**R1 — the reduction of the loop** *(scope corrected after review: what is
derived is the reduction and its holonomy, not the itinerary; the
`source → A → J → B → source` history and the geodesic realignment at the
detectors remain model choices).* In the quaternion model, frame rotation `R = Ad_g` is
`q ↦ q g⁻¹` and the mouth transport is `J = L_{−j}`; the two commute. For
the history `source → A → J → B → source`,

```
q₅ = J·[J·q₀ g₁⁻¹]·g₂⁻¹·g₃⁻¹ = −q₀ (g₃g₂g₁)⁻¹          (4e-16)
```

so the throat transports contribute the spinor sign `J² = −1` and nothing
else; the holonomy is fibre-independent (`2e-15`); and since
`R(−u→v) = R(u→−v)`, `R(v→−x) = R(−v→x)`, the frame holonomy is that of the
geodesic triangle `x → u → −v → x`: a rotation about `x` by `Ω(x; u, −v)`
(`2e-15`), with SU(2) lift `G = cos(Ω/2) + sin(Ω/2) x` exactly, sign included
(`9e-16`), equal to `sgn D(x, u, −v)` on the closure circle (`4e-16`). **The
chosen itinerary reduces to the triangle of the previous round with the
partner sign** `v → −v` (singlet), **and its closure holonomy is `−sgn D`.**
What this removes is the *ansatz status of the triangle given the itinerary*;
the itinerary itself and the geodesic realignment remain chosen.

**R2 — the oriented current is the quantum joint law, without projectors.**
*(Tested on the derived loop itself after review finding 2; see R2b.)*
`W̃(s) ∝ ∫_Γ D_s dσ / (2|u×v|)` and `∫_Γ x dσ = 0` (`3e-17`) give
`∫_Γ D_s dσ = L_Γ(1 + u·v)`, hence

```
P₊₊ = P₋₋ = (1 + cos γ)/4,     P₊₋ = P₋₊ = (1 − cos γ)/4,     E = cos γ
```

(deviation from the quadrature `6e-17`), `S = 2√2`. The Bargmann identity
`D/4 = Re Tr(P_x P_u P_v)` represents this; it does not produce it.

**R2b — the derived singlet loop, computed directly.** A first version of
this document obtained the singlet by appending "`−cos γ` with the partner
sign" to the triplet result, which left the strongest physical claim one
verbal substitution beyond the tested function. The derived object
`x → u → −v → x`, with current `−D(x, u, −v)/(2|u×v|)` on the same closure
circle, is now integrated directly: the common `−1` from `J²` and the common
`1/(2|u×v|)` cancel in the normalisation, and

```
P(s_A, s_B) = (1 − s_A s_B cos γ)/4 ,        E = −cos γ
```

to `≤ 1e-17` at `γ = 0.3, 1, 2`, with all four normalised weights positive.

## 3. What does not decide the fork (R3–R5)

* **R3, the sector prior.** Chosen; moves `E`; fixes nothing about the
  branch weighting. **And it is load-bearing on the *oriented* branch too**
  (review finding 1): the oriented sector integrals are proportional to
  `1 ± s_A s_B cos γ`, so under `r = π_like/π_unlike`

  ```
  E_r^triplet = [r(1+cos γ) − (1−cos γ)] / [r(1+cos γ) + (1−cos γ)]
  E_r^singlet = [r(1−cos γ) − (1+cos γ)] / [r(1−cos γ) + (1+cos γ)]
  ```

  giving `0.25243, 0.54030, 0.74031` and `−0.74031, −0.54030, −0.25243` at
  `γ = 1` for `r = 0.5, 1, 2`: the quantum law only at `r = 1`, with
  marginals `1/2` throughout, so no-signalling does not constrain it. **The
  holonomy weighting alone therefore does not deliver the quantum joint
  law.**
* **R4, stationarity — corrected after review.** The repository's condition
  is *extremal action*, and no action functional exists in
  `history/closure.py`. Substituting stationarity of the geometric phase is
  a **proxy**, labelled as such, and it does not test the repository's
  condition. The proxy is moreover analytically incompatible with sharp
  closure: on `N = 0` the tangential gradient is

  ```
  ∇_{S²}Ω |_{N=0} = 2q/D = 2(u×v)/D,     |∇Ω| = 2|u×v|/|D| > 0
  ```

  so closure and phase-stationarity are **disjoint** for non-collinear
  analyzers (minimum gradient norm `0.5107` on the closure set;
  branch-wrapped finite-difference residual `3.4e-7`). The only points with
  `D = 0` are `x = −u` and `x = −v`, where `D = N = 0` and the `arg` chart
  is singular — not stationary points. *This is stronger than the
  pre-registered "isolated points", and it strengthens the diagnosis: no
  variational principle available in the repository can choose between the
  two measures.* (`Ω` jumps by `4π` across the closure circle where `D < 0`;
  the holonomy `e^{iΩ/2} = −1` is continuous there, which is why the
  finite difference must be branch-wrapped.)
* **R5, the Pin structure.** It supplies the closure holonomy `−sgn D` as
  a **label** on each closed history (R1). Whether a label enters the
  measure as a weight is not a property of a bundle. The repository's
  classical field equations are linear and carry no measure on histories.
  No computation available here selects a weighting.

## 3b. The oriented current does not give negative probabilities

Although `D` changes sign pointwise on the closure circle (on an arc of
fractional length `γ/2π` for like outcomes), the integrated sector weights
are non-negative: with `∫_Γ x dσ = 0`,

```
∫_Γ D_s dσ = 2π(1 + u·v) ≥ 0        (quadrature residual 1.8e-15)
```

for every outcome pair, vanishing only at `u = −v`. The holonomy weighting
therefore uses destructive cancellation **internally** and yields
non-negative normalised sector weights. This removes the naive
"negative probabilities" objection and is what makes the oriented branch a
candidate at all rather than an obvious non-starter.

**The wave analogy stays qualified** (review finding 4). Non-negativity of
the integrated current is not yet an event-frequency law. If `D` is an
amplitude or a current, a classical detector normally responds to a
quadratic (energy) functional of it; if it is an event measure, the signed
cancellation needs a reason. So "the Pin/Hopf data make an oriented current
mandatory" would not by itself complete the probability derivation — the
linear-versus-quadratic readout is a separate open item.

## 4. Where the quantization gap is — three inputs, not one binary choice

Not the spin structure (derived in the spin-frame rounds). Not the
reduction of the chosen itinerary to the triangle (derived, R1). Not Bell:
both candidates violate CHSH by measurement dependence with exact detector
no-signalling. *Corrected after review:* it is **not** a single binary
question either. Three inputs remain underived:

1. **Branch aggregation.** Are observed event frequencies the **positive
   count** of closed histories, `Σ_branches |D|` (`S_max = 2.1423`), or their
   **oriented sum with the closure holonomy**,
   `Σ_branches e^{iΩ/2}|D| = Σ D` (`2√2`)?
2. **The relative sector coefficients.** The equal outcome-sector prior is a
   chosen counting measure, it moves the correlation on **both** branches,
   and the quantum law appears only at `r = 1` (R3, R3b). No symmetry of the
   model at fixed settings fixes it, and no-signalling does not either.
3. **The readout.** Even granting an oriented current, why would observed
   frequencies be **linear** in the integrated current rather than quadratic
   in it — the usual classical detector response to an amplitude — or given
   by some other functional? (Review finding 4.)

So the honest headline is: **the derived Pin/spin-frame holonomy supplies
exactly the sign that turns positive coarea into the Born-shaped oriented
sector integrals, but current BAM selects none of (1), (2), (3).** That is a
sharper localisation of the quantization gap than the repository had before
this round, and it is narrower than "one binary choice".

## 5. Open falsifiers carried forward

* **R6, operational no-signalling to the past.** `ρ(x|a,b)` depends on the
  future settings. A source coupling whose read-out statistics depend on
  `(a, b)` would falsify the model as a physical theory; the model has no
  such coupling yet.
* **A derivation of the sector coefficients**, from a symmetry or a chain
  argument that fixes `π_like = π_unlike` independently of the answer.
* **A derivation of the weighting.** From the Pin structure, from a
  stationary classical field equation, or from a counting argument that
  distinguishes the two branches physically. Any of these would resolve
  the fork; by the rule, only such a derivation may.
* **The sharper form of the fork, recorded as an open audit item and not a
  criterion.** Do the Pin/Hopf data make the closure locus naturally an
  *oriented current with local coefficients*, whose observables are
  integrals of sections? Then the sign `e^{iΩ/2}` is geometrically
  mandatory and the oriented branch is forced. Or is the physical object a
  *measure on histories*? Then positivity forces `|D|`. Neither is
  established here, and this is the form in which the question should be
  attacked next. Establishing it would settle (1) but **not** (2) or (3).

## 6. Dependency ledger

```
loop                =  loop( source pair [x, −x via J], geodesic realignment at A, B [chosen],
                             mouth transport J [derived: L_{−j}] )  →  triangle x → u → −v → x [derived]
closure holonomy    =  −sgn D(x, u, −v)  [derived from J² = −1 and the lift G]
positive coarea     =  |D|/(2|u×v|)      [derived conditioning; counting measure on sectors: chosen]
oriented current    =  e^{iΦ} × coarea  [derived label; adopting it as weight: the open step]
sector prior        =  counting          [chosen]
```
