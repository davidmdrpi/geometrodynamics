# Sharp phase closure of a closed history: measurement dependence at the source, without an imported measure or signalling

*Fifth round of the finite-mouth chain. Pre-registered in
`docs/closure_measurement_dependence_prereg.md` (`1b0144e`); module
`geometrodynamics/bulk/closure_measurement.py`; tests
`tests/test_closure_measurement.py`; probe `closure_measurement_probe.py`.*

## Verdict

`CLOSURE_INDUCES_SETTING_DEPENDENT_SOURCE_MEASURE_NO_SIGNALLING_NOT_BORN`

The repository's closed-history axiom, applied to the classical state the
spin-frame rounds supply and with no amplitude, no Gaussian and no width,
does three things and fails to do a fourth:

1. it makes the admissible source ensemble depend on **both** future
   analyzer settings (the source direction is confined to the great circle
   through `a` and `b`, with a density that depends on the pair);
2. it does so with **exactly** no-signalling marginals, `1/2` to `1e-12`;
3. it **violates the CHSH bound**: `S = 2.1423` at the standard angles,
   which is also the maximum over settings, below `2√2`;
4. it does **not** reproduce the quantum correlation: `E(γ)` is a
   different, closed-form function of the angle. The quantum correlation is
   what one gets from the **holonomy-weighted** closure density
   `e^{iΩ/2} |D| = D` (see the correction below).

*Correction after review (`docs/closure_current_prereg.md`, C1–C4).* A first
version of this document called the signed density an inserted
quasi-probability and said adopting it "imports the quantum rule". On the
closure set `e^{iΩ/2} = sgn D`, so the signed density is the coarea density
times the actual holonomy of the `0`/`π` closure branch; the open question
is whether histories combine as positive counts or with their holonomy,
and that is the subject of the next round. Three further corrections: the
model uses the phase-closure component of BAM under a geodesic-realignment
loop ansatz, not the full closure axiom (stationarity is unimplemented and
the branch sign is discarded by the repository's Gaussian weight); the
equal prior over the four outcome sectors is a chosen counting measure,
entered in the ledger; and what is proved is detector no-signalling, not
operational no-signalling to the past.

This is outcome **D** of the Born-rule round's typology, obtained by
calculation rather than asserted. It evades Bell's theorem by violating
measurement independence, at the standard cost, and it stops short of
quantum mechanics by a computable margin.

## 1. What was audited first

`history/closure.py` is the repository's two-boundary principle. Its Bell
computation weights each outcome branch by `|amplitude|² × closure_weight`,
with the amplitude imported from `bell.analyzers` and the closure weight a
Gaussian in the phase mismatch of chosen width `σ = 0.6`. That derives
nothing about measurement dependence; the width control below shows how
much the answer moves with `σ` (`E(1) = 0.170` at `σ = 0.6` against the
sharp `0.398`).

## 2. The model, with its choices named

* Source variable `x ∈ S²`, the created pair's frame direction, Haar prior
  (invariant; *chosen: physical rather than gauge*). The partner carries
  `−x`.
* Closed history `source → A → B → source` as a loop of frame directions
  `x → u → v → x`, `u = s_A a`, `v = s_B b`, geodesic legs. *Chosen:*
  detection realigns the frame to the outcome direction along the geodesic;
  the outcome signs are boundary data of the history.
* Closure phase: the spin-½ geometric phase, half the solid angle
  `Ω = 2 atan2(N, D)`, `N = x·(u×v)`, `D = 1 + x·u + u·v + v·x`. *Derived:*
  it is the spin-frame holonomy of the loop.
* Closure: `Ω/2 ≡ 0 or π (mod 2π)` ⟺ `N = 0` ⟺ `x ∈ Γ(a,b)`, the great
  circle through both settings. *Repository axiom.*
* Conditioning variable: the **phase**. `window_monte_carlo` samples the
  window `|Ω mod 2π| < ε`, and its `ε → 0` limit is the coarea measure *with
  respect to `θ`*, density `|D|/(2|u×v|)`. *The density is derived once the
  variable is fixed; the variable is **chosen**.* Conditioning on `|N| < ε`
  instead has the same support `Γ` and a different limit — uniform in
  arclength, `E = 0` — because `|∇N| = |u×v|` is constant on `Γ` while
  `|∇θ| = |u×v|/|D|` is not, and `|N|` is identical in all four outcome
  sectors. Conditioning on a measure-zero set is not determined by the set
  (Borel–Kolmogorov). Justified by the axiom being stated on phase, not
  derived from the closure set. See `docs/conditioning_variable_prereg.md`.
* `P(s_A, s_B | a, b) ∝ ∫_Γ |D| dσ / (2|u×v|)`. Not a deterministic
  detector: several histories close for a given `x`, and the outcome is
  which closed history is realised.

## 3. Results, all pre-computed in closed form and then reproduced

With `γ = ∠(a,b)`, `c = cos(γ/2)`, `s = sin(γ/2)`:

```
W(+,+) = W(−,−) = 2 + c(π−γ)/s        W(+,−) = W(−,+) = 2 + sγ/c
E(γ) = [c²(π−γ) − s²γ] / [2 sin γ + c²(π−γ) + s²γ]
```

| `γ` | `E` (quadrature) | `E` (closed form) | `cos γ` |
|--|--|--|--|
| `0.3` | `0.82095` | `0.82095` | `0.95534` |
| `π/4` | `0.53557` | `0.53557` | `0.70711` |
| `1` | `0.39850` | `0.39850` | `0.54030` |
| `π/2` | `0` | `0` | `0` |
| `2` | `−0.30350` | `−0.30350` | `−0.41615` |

`E(π−γ) = −E(γ)`; marginals `1/2` exactly for every variant and angle.
`S(0, π/2, π/4, −π/4) = 4E(π/4) = 2.1423`, and a scan over settings finds
no larger value: `S_max = 2.1423 < 2.8284`.

**Independent construction (P7).** Uniform `x` on `S²` with the window
`|Ω mod 2π| < ε`: `E(1) = 0.3729, 0.3923, 0.3973, 0.3972` at
`ε = 0.4, 0.2, 0.1, 0.05` (`2×10⁶` samples), converging to the coarea value
`0.3985`. This confirms the positive density `|D|`, not the signed `D`.

**The holonomy-weighted variant is quantum mechanics (P5, corrected).**
With `D` in place of `|D|`, `E = (c²−s²)/(c²+s²) = cos γ` exactly and
`S = 2√2`. Since `e^{iΩ/2} = sgn D` on the closure set, `D/(2|u×v|)` is the
coarea density weighted by the closure holonomy of each branch, the `π`
branch entering with opposite orientation on the arc where `D < 0`
(fractional length `γ/2π` for like outcomes; measured `0.0796, 0.1592,
0.3183` at `γ = 0.5, 1, 2`). The identity `D/4 = Re Tr(P_x P_u P_v)`
(checked to `2e-16`) is then a representation of what the closure geometry
already contains, not its source: with `∫_Γ x dσ = 0` the oriented current
gives `P_like = (1+cos γ)/4`, `P_unlike = (1−cos γ)/4` with no projectors.
What this identifies precisely is *where* the quantum rule sits relative to
the classical closure model: the same closure set, with the branches
summed with their holonomy instead of counted.

**The strict-zero variant (P6).** Dropping the `π` branch of the axiom
keeps only the `D > 0` arc: `E(1) = 0.46495`, `S = 2.4649`. Still not
`cos γ`.

## 4. Measurement dependence, stated correctly

*Correction to the pre-registration.* P1 said the supports of `ρ(x|a,b)`
and `ρ(x|a′,b′)` are distinct great circles, so `TV = 1`. That holds for
**non-coplanar** settings (computed: `TV = 1`). In the standard Bell
configuration every analyzer lies in one plane, every pair shares the same
great circle, and the dependence is in the **density** on it: the
outcome-summed conditioned density `∝ Σ_s |D_s(x)|/(2|u×v|)` changes with
the pair (total variation `0.018` between the pairs `(0, 1)` and
`(0.5, 1)`, `0.063` between `(0, π/4)` and `(π/2, −π/4)`). Either way
measurement independence fails; the correction is to the mechanism of
failure in the coplanar case, not to its existence.

## 5. Controls

| control | result |
|--|--|
| loop topology: `x → u → x` (does not link both detectors) | solid angle `3e-15`; closure automatic; no conditioning |
| local detectors (previous round's C5 on the same source, independent `y_A, y_B`) | `S = 0.94 ≤ 2` |
| Gaussian width (the repository's `σ = 0.6`, then `0.3, 0.1, 0.03`) | `E(1) = 0.170, 0.293, 0.380, 0.395 → 0.398`: depends on `σ`, sharp limit is the coarea value |
| sign (`v → −v`) | `E → −E`, singlet versus triplet; `S` unchanged |
| symmetry | marginals `1/2` at every angle; `E(π−γ) = −E(γ)` |

## 6. Dependency ledger

```
ρ(x|a,b)  =  ρ( Haar on S² [invariant prior; pair direction physical: chosen],
                closure axiom Ω ≡ 0 mod 2π [repository axiom],
                geodesic-realignment detection model [chosen],
                conditioning variable = phase, not N [chosen; justified by
                    the axiom being stated on phase],
                coarea density given that variable [derived: 1/|grad theta|] )
E(γ)      =  E( ρ(x|a,b) [above], outcome signs as history boundary data [chosen: D-type] )
oriented current  =  e^{iΩ/2} × coarea = D/(2|u×v|)  [derived label; adopting it as weight: the open step]
sector prior      =  counting on the four outcome sectors  [chosen]
```

## 7. What this means for the program

* The question the previous round posed is answered: the BAM two-boundary
  principle **does** make the source ensemble depend on both settings, with
  no signalling and no imported measure beyond an invariant prior and a
  parameter-free conditioning. The mechanism is the closure constraint
  linking the source direction to both analyzers through one loop; a loop
  that does not link both does nothing.
* It is therefore not in Bell's local class, and it violates CHSH. This is
  the first place in the repository where a Bell violation is *computed*
  from a classical global constraint rather than from an inserted singlet.
* Under the positive count it is not quantum mechanics. The correlation
  function is `[c²(π−γ) − s²γ]/[2 sin γ + c²(π−γ) + s²γ]`, and the maximal
  violation is `2.14`, not `2.83`. The distance to quantum mechanics is
  exactly the distance from `|D|` to `D`: from counting the two closure
  branches to summing them with their holonomy.
* Two modelling choices carry the result and are open: that detection is
  geodesic realignment to the outcome direction, and that outcomes are
  history boundary data rather than functions of the source state. A
  different detection model changes the loop and hence `E`; the next
  question is whether any classical detection dynamics selects the loop
  and the measure, and whether the closure principle can be stated for the
  field rather than for the ray of frame directions.
