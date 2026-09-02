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

**R1 — the loop.** In the quaternion model, frame rotation `R = Ad_g` is
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
Pin-derived loop is the triangle of the previous round with the partner
sign** `v → −v` (singlet), **and its closure holonomy is `−sgn D`.** The
loop is no longer an ansatz; the geodesic realignment at the detectors is.

**R2 — the oriented current is the quantum joint law, without projectors.**
`W̃(s) ∝ ∫_Γ D_s dσ / (2|u×v|)` and `∫_Γ x dσ = 0` (`3e-17`) give
`∫_Γ D_s dσ = L_Γ(1 + u·v)`, hence

```
P₊₊ = P₋₋ = (1 + cos γ)/4,     P₊₋ = P₋₊ = (1 − cos γ)/4,     E = cos γ
```

(deviation from the quadrature `6e-17`), `−cos γ` with the partner sign,
`S = 2√2`. The Bargmann identity `D/4 = Re Tr(P_x P_u P_v)` represents this;
it does not produce it.

## 3. What does not decide the fork (R3–R5)

* **R3, the sector prior.** Chosen; moves `E`; fixes nothing about the
  branch weighting.
* **R4, stationarity.** Named in the repository, implemented nowhere. And it
  could not help even if implemented as stationarity of the closure phase
  in the source direction: `Ω(x; u, v)` has **no critical points** on `S²`
  minus the two singular points `−u, −v` where `D = N = 0` (`min |∇Ω| = 0.51`
  on a `721 × 1441` grid). Its level sets are Lexell circles through `−u`
  and `−v` (a sampled level set is coplanar to `2e-3` and its plane contains
  both to `< 1e-2`). *Correction to the pre-registration, which predicted isolated
  points; the conclusion is stronger, not weaker.*
* **R5, the Pin structure.** It supplies the closure holonomy `−sgn D` as
  a **label** on each closed history (R1). Whether a label enters the
  measure as a weight is not a property of a bundle. The repository's
  classical field equations are linear and carry no measure on histories.
  No computation available here selects a weighting.

## 4. Where the quantization gap is

Not the spin structure (derived in the spin-frame rounds). Not the loop
(derived, R1). Not Bell: both candidates violate CHSH by measurement
dependence with exact detector no-signalling. Not the relative outcome law:
the oriented current gives it analytically. It is one binary question,
stated in the repository's own objects:

> Are observed event frequencies the **positive count** of closed histories,
> `Σ_branches |D|`, or their **oriented sum with the closure holonomy**,
> `Σ_branches e^{iΩ/2}|D| = Σ D`?

The first is classical event counting and predicts `S_max = 2.1423`. The
second is an amplitude-like signed history measure and predicts quantum
mechanics. The distance between BAM and quantum mechanics, for this
degree of freedom, is exactly that choice, and it is made nowhere in the
repository.

## 5. Open falsifiers carried forward

* **R6, operational no-signalling to the past.** `ρ(x|a,b)` depends on the
  future settings. A source coupling whose read-out statistics depend on
  `(a, b)` would falsify the model as a physical theory; the model has no
  such coupling yet.
* **A derivation of the weighting.** From the Pin structure, from a
  stationary classical field equation, or from a counting argument that
  distinguishes the two branches physically. Any of these would resolve
  the fork; by the rule, only such a derivation may.

## 6. Dependency ledger

```
loop                =  loop( source pair [x, −x via J], geodesic realignment at A, B [chosen],
                             mouth transport J [derived: L_{−j}] )  →  triangle x → u → −v → x [derived]
closure holonomy    =  −sgn D(x, u, −v)  [derived from J² = −1 and the lift G]
positive coarea     =  |D|/(2|u×v|)      [derived conditioning; counting measure on sectors: chosen]
oriented current    =  e^{iΦ} × coarea  [derived label; adopting it as weight: the open step]
sector prior        =  counting          [chosen]
```
