# The self-consistent network loop (PR #217)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #216 built the two-port throat and
> showed the network's projection is advanced — but treated the loop at
> first order, and postulated the I± assembly with the network weight.
> This PR constructs the **full two-mouth transfer system**, solves the
> loop **self-consistently**, and derives the effective Green function
> **before** any comparison with I±. The companion probe machine-checks
> every claim (~60 s).

## 0. The self-consistent constraint

The returned confirmation re-enters the network — the closed universe
refocuses it onto mouth A again — so the field at the crossing obeys
the CTC wave constraint

```
F = F₀ + Λ·F        ⟹        G_eff(ω) = g / (1 − Λ(ω))
```

with everything in Λ derived. Solving this honestly forced three
corrections/completions of #216, each machine-checked.

## 1. Clock-rate-correct traversal

The throat is static in its own proper time, so the ports and the
interior loop phase run at the throat frequency ω_τ = ω/rate_A, and the
emergent global frequency is **ω·rate_B/rate_A** — a slow-clock
(deep-well) exit mouth *redshifts* the wave, as climbing out of a well
must. This **corrects an inverted rate ratio shipped in #216**
(exercised only at equal rates there; the unit test asserted blueshift
4.8 where physics demands redshift 1.2). The exit time follows the
exact clock composition t_exit = [ρ_A(t_entry − o_A) + τ_local]/ρ_B +
o_B.

**The elastic selection rule.** Unequal rates return the confirmation
at ω′ ≠ ω — state closure fails at the crossing. Elastic confirmation
*requires* rate-matched mouths at traversal epoch: the differential
aging that builds Δ_BA must have happened earlier (MTY history), with
only the offset surviving. #216's equal-rate analysis is thereby
*derived*, not assumed.

## 2. The transfer system and the ring-validated eigenvalue

The full signal-flow system — both barriers explicit, interior loops
and windings as separate cycles — resolves to the closed form:

- **matrix resolvent** (I − M)⁻¹ = 1/(1 − Λ) to 9×10⁻¹⁵;
- **nested double resummation** (interior loops × windings) to 10⁻¹⁵.

**The value-transport eigenvalue.** The self-consistency of the field
is governed by the returned *value* per emitted value. Summing the echo
train with its k-dependent arrival displacements:

```
Λ(ω) = t_net(ω_τ) · deco · e^{+iωD_loop} ,    D_loop = d_A + d_B + τ_glob + Δ_BA
```

Validated **from first principles against a ring spectrum**: a ring
(free arc + the two-barrier throat section) has transfer-matrix modes
at ω = 2.7429/3.6000, and the Λ = 1 comb predicts 2.7390/3.6071 — 0.2%,
the residual being the O(|r|) backscatter splitting. The opposite phase
convention misses by 0.16 and is **excluded**.

**The #216 convention correction.** #216's closure comb tracked the
engine's time-phase loop amplitude (`retro_phase_match`'s bookkeeping),
which coincides with value transport on the tower's exterior legs
(e^{−2πin} = 1) but not for the transit phase — the offset between the
two is exactly e^{−iω(d_A+d_B+τ)} (machine-checked to 10⁻¹⁶). At time
closure the carrier actually closes on **the throat's own scattering
phase**, arg t_net = 0 (mod 2π). The deform knob (φ = ωδ) and all
magnitude/advanced-projection results of #216 are unchanged.

## 3. Wigner-correct closure and the transaction points

Group closure demands D_loop = −τ_W (the composite Wigner delay);
carrier closure demands arg Λ = 0 (mod 2π). Solving **both
simultaneously**:

- **Above the barrier** (ω₀ = 3, a tower mode): the throat phase is
  nearly flat, so the mouths' geometric transfer phase θ* supplies the
  carrier tuning, and Δ* = −(d_A+d_B+τ) − τ_W. Result:
  Λ = 0.999998 + 0i exactly at the point, and the **Wigner-corrected
  packet lands ON the crossing**: return peak 0.006 (versus 0.175 = τ_W
  uncorrected), phase-aligned to 0.005 — the fully closed transaction,
  envelope and carrier together.
- **Below the barrier** (ω₀ = 0.5): both closure equations are solved
  by the throat alone — τ* = 5.1105 with residual 2×10⁻¹² — and the
  solution sits **on the interior Fabry–Pérot resonance**
  (|t_net|² = 0.969, storage delay 31). **Completed transactions live
  on the throat's resonance comb** — #216's "resonant confirmation"
  observation upgraded to a solved fixed point.

## 4. Source and mouth state evolution

- **Input–output rates**: κ_tot = (T_A+T_B)/2τ* = 0.0652 matches the
  derived linewidth of |t_net|² (0.052, within the CMT tolerance at
  finesse ~8) and the resonant storage time 2/κ_tot = 30.7 matches the
  Wigner delay 31.1 (1.4%) — the mouth-state dynamics is the coupled-
  mode reading of the derived scattering.
- **The source fixed point**: iterating x ← 1 + Λx converges to G_eff
  at rate |Λ| exactly (measured 0.23239907 vs |Λ| = 0.23239900) — the
  winding/generation picture of the self-consistent constraint.
  At a completed transaction |Λ| → 1: convergence is marginal — the
  transaction is a *persistent* self-consistent history.
- **The cavity ledger**: the interior amplitude builds monotonically to
  its steady state through the transfer-matrix iteration (a genuinely
  long-lived interior transient — spectral radius 0.9991).

## 5. The effective Green function — then I±

With G_eff = g/(1 − Λ) in hand, its structure is measured first:

- **the quartic line**: at a completed transaction, group closure makes
  the loop phase *stationary* — the resonance is anomalously flat:
  |1−Λ|² ≈ (1−|Λ|)² + (|Λ|θ″δω²/2)², FWHM = 2√(2(1−|Λ|)/|θ″|) —
  measured 0.0116 vs predicted 0.0118. Completed transactions are
  **robust to detuning**;
- **the Lorentzian line** at carrier-only closure: FWHM = 2(1−|Λ|)/|θ′|
  — measured 0.0439 vs 0.0456;
- **passivity**: max|Λ(ω)| = 0.9999993 ≤ 1 over the sweep — two-port
  unitarity makes every self-consistent solution stable or marginal.
  **The Novikov fixed point cannot run away**: no paradoxical
  self-amplification, existence and uniqueness of the CTC field for
  Λ ≠ 1.

Only now, the comparison with the orderings:

- the **O(Λ) truncation** of the winding series is exactly #216's
  assembly K₁ = I₊ + ΛI₋;
- the **full resummation** renormalizes the confirmation weight to
  **Λ/(1−Λ)** (series identity to 10⁻¹¹ at a partial-confirmation
  point); at completion the weight reaches ~5×10⁵ — **the divergence is
  the transaction pole**;
- the **ε unification**: the confirmation deficit 1−|Λ| → 0 at
  completion (2×10⁻⁶ vs 0.11 partial) — the completed transaction is
  the #213 coherent ε → 0 limit, and the deficit is #214's absorber
  damping in its transactional form. Three PRs' ε's, one number.

## 6. Honest scope

- The #216 convention correction affects the closure *phase*
  bookkeeping only; magnitudes, the advanced projection, and the deform
  knob are unchanged. The engine weld statement of #216 remains true of
  the engine's own convention; migration remains staged.
- Continuous-ω analysis interpolates between physical tower modes; the
  above-barrier carrier tuning uses the mouths' geometric transfer
  phase (a free parameter of the frozen network), while the
  below-barrier point is solved by the throat alone.
- CMT rates are a moderate-finesse consistency check, not the
  derivation.
- The winding resummation assumes a frozen, linear network per winding;
  back-reaction of stored cavity energy on the throat is out of scope.
- Classical, zonal scalar, ℓ = 0 greybody; network history (MTY aging)
  posited.

## 7. What would falsify this

- Ring modes not matching the Λ = 1 comb — the eigenvalue convention
  would be wrong. (Checked: 0.2%; opposite convention excluded at 40σ
  of the match.)
- A packet missing the crossing after the Wigner correction. (Checked:
  0.006 vs grid 0.001.)
- A transaction point off the interior resonance below the barrier.
  (Checked: |t_net|² = 0.969 at the root.)
- |Λ| > 1 anywhere — a runaway time machine. (Checked: ≤ 0.9999993.)
- The fixed-point rate ≠ |Λ|, or the resummed weight ≠ Λ/(1−Λ).
  (Checked: 8×10⁻⁸ and 10⁻¹¹.)

## 8. Companion probe

`experiments/closure_ledger/self_consistent_network_loop_probe.py`
(T1–T8, ~60 s), plus module upgrades in
`geometrodynamics/transaction/network.py` (clock-rate-correct
`t_AB`/`traverse_throat`/`emergent_frequency`, the value-transport
`loop_eigenvalue`, `effective_green`) and 16 unit tests.

**Verdict:**
`THE_LOOP_SOLVED_SELF_CONSISTENTLY_G_EFF_DERIVED_COMPLETED_TRANSACTIONS_SIT_ON_THE_THROAT_RESONANCES_AND_THE_NOVIKOV_FIXED_POINT_IS_PASSIVE`
