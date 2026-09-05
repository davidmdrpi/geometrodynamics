# Pre-registration: is the `|D|` density forced by the closure set?

**Status: frozen before any correction.** Eighth round of the finite-mouth
chain, following PR #281 (`92a915b`). A **correction round**, prompted by an
external audit of that snapshot dated 5 September 2026, which found that the
repository conditions on a *phase* window while describing the result as Haar
conditioned on `N = 0`. Success criteria may change afterwards only by an
explicit correction note appended to this file.

## 0. What prompted it

`geometrodynamics/bulk/closure_measurement.py:15` states

> Haar conditioned on `N = 0` is the coarea measure, density `|D|/(2|u×v|)`

and `docs/closure_measurement_dependence.md:150` ledgers `coarea conditioning`
as **`[derived; window limit]`**. The same claim appears at `README.md:4831`
and `docs/closure_measurement_dependence.md:70`. The `|D|` density is the input
every later result rests on: round 5's positive branch, round 6's oriented
current, round 7's Morse–Bott component masses.

The audit's counter-observation is that `window_monte_carlo` windows on
`|Ω mod 2π| < ε` — the **phase** — not on `|N| < ε`.

## 1. Results established before freezing (analytic oracle)

Verified numerically before this file was committed; frozen as oracles, so a
disagreement in the round is a bug in the round, not a result. With
`u = s_A a`, `w = −s_B b`, `q = a×b`, `N = x·(u×w)`, `D = 1 + x·u + u·w + w·x`,
`θ = atan2(N, D)`:

* **O1.** `u×w = −s_A s_B (a×b)`, so `|N|` is **identical in all four outcome
  sectors**. A `|N| < ε` window therefore selects the same set of `x` in every
  sector, and the four counts are exactly equal (measured: 5956 each at
  `γ = 1`, `ε = 0.01`, `n = 5×10⁵`).
* **O2.** On `Γ`, `x ⊥ q`, so `∇_{S²}N = P_x(q) = q` and `|∇_{S²}N| = |a×b| =
  sin γ` is **constant**. The coarea density of Haar with respect to `N` is
  therefore **uniform in arclength**, not `|D|`.
* **O3.** On `Γ`, `∇_{S²}θ = q/D`, so `|∇_{S²}θ| = |q|/|D|` and the coarea
  density with respect to `θ` is proportional to **`|D|`**. This is where the
  repository's density comes from.
* **O4.** Consequently the `N`-window with equal sector priors gives
  `P(s_A,s_B) = 1/4` and `E = 0` exactly, at every angle.
* **O5.** The phase window converges to the round-5 law: `E(1)` measured
  `0.3971, 0.3955, 0.3983, 0.3953, 0.4009` at `ε = 0.05 … 0.002` against the
  closed form `0.398497`.

## 2. The question, and the rule fixed in advance

Conditioning a measure on a measure-zero set is not determined by the set
(Borel–Kolmogorov). O2 and O3 exhibit two limiting families with the same
support `Γ` and different limits. So:

**Is the `|D|` density forced by the closure set, or is the conditioning
variable an additional input?**

Permitted verdicts:

* `CONDITIONING_VARIABLE_IS_FORCED_BY_THE_CLOSURE_SET`
* `CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_JUSTIFIED_BY_THE_PHASE_AXIOM`
* `CONDITIONING_VARIABLE_IS_A_CHOSEN_INPUT_WITH_NO_JUSTIFICATION_IN_THE_REPOSITORY`
* or a narrower correction.

**Rules.**

1. The phase window may **not** be preferred because it gives the interesting
   answer. A justification must be a structure already in the repository,
   named by file and line, that privileges the phase as the conditioning
   variable — not an appeal to the resulting correlation.
2. Equally, the `N` window may not be preferred because it is "simpler". If
   both are unjustified, the third verdict is the honest one.
3. **No round 5–7 number may change.** This is a correction to the *status* of
   an input, not to any computation. A test asserts the round-5 correlations,
   the round-6 oriented law and the round-7 masses are bit-unchanged. If any
   number moves, the round has found something larger than a mislabelled
   ledger entry and must say so.

## 3. Two further narrowings from the same audit

**B — `κ`'s status.** Round 7's prose says the round "adds a fourth (`κ`)",
which reads as a fourth *universal* underived input alongside branch
aggregation, sector coefficients and readout. `κ` is the normalisation of
`e^{iκS_H}` and exists only inside the holonomy-trace route that round 7
closed. To be narrowed to a parameter of that attempted route.

**C — velocity-field non-uniqueness in `docs/born_rule_equivariance.md`.**
Continuity does not determine a unique velocity field: if `∇·K = 0` then
`J + K` satisfies the same continuity equation and `v' = (J+K)/ρ` is another
equivariant flow wherever `ρ > 0`. To be demonstrated with an explicit
divergence-free `K` rather than asserted, and the cited Goldstein–Struyve
uniqueness result recorded as assuming Bohmian dynamics plus locality
conditions rather than deriving them.

## 4. Expected verdict, stated in advance

The second: `history/closure.py:12` states the repository's third closure
condition as *"Phase closure: total phase around every loop ≡ 0 or π
(mod 2π)"* — a condition on **phase**. If that is judged to privilege the
phase as the conditioning variable, the choice is justified but remains a
choice, and the ledger entry must move from `derived` to `chosen` with the
justification named. Obtaining this is not a success; the round's job is to
record the status correctly.

## 5. Dependency ledger to be printed

The corrected round-5 ledger, with the conditioning variable as its own line.
