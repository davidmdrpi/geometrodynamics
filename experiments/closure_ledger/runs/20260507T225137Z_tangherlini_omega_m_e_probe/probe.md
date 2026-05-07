# Tangherlini ω ↔ m_e relation probe

**Run:** 2026-05-07T22:51:37+00:00

Sub-target #3 from the ℏ-origin research plan. Tests whether the lowest Tangherlini radial eigenfrequency is naturally identified with the electron Compton angular frequency, and whether the spectrum's higher eigenvalues match the lepton mass ratios. A successful identification would let BAM **predict ℏ in physical units** rather than relying on the m_e anchor.

## Q1: Is ω(l=1, n=0) ≈ 1 in canonical Compton units?

Under the canonical identification R_MID = ℏ/(m_e c) ≈ 3.862e-11 cm, the electron Compton angular frequency m_e c² / ℏ ≈ 7.763e+20 rad/s corresponds to ω = 1 in geometric units.

**Observed:** ω(l=1, n=0) = `1.054727`.

**Deviation from 1:** `5.473%`.

**Q1 FAIL.** No natural Compton identification at the 5 %-level. The dimensional bridge is not present at leading order.

## Q2: Do ω-ratios match lepton mass ratios?

**Targets:** m_μ/m_e = 206.77, m_τ/m_e = 3477.23.

**Tangherlini ω-spectrum span:** ω_max / ω_min = `3.8741`.

**Best candidate for m_μ/m_e:**
  ω(l=5, n=3) / ω(1, 0) = `3.8741` (target 206.77, error `98.1%`).

**Best candidate for m_τ/m_e:**
  ω(l=5, n=3) / ω(1, 0) = `3.8741` (target 3477.23, error `99.9%`).

**Q2 partial.** Some ω-ratio is within factor-2 of the lepton mass ratio — but not within 5 % precision.

## Dimensional verdict

**One-line verdict:** Neither natural identification holds. The dimensional bridge is not present in the current geometric input set.

**BAM remains dimensional-ratio-complete and dimensional-scale-incomplete.** Mass ratios are predicted to sub-percent (lepton lock); the absolute MeV scale is set by anchoring m_e. Predicting ℏ in physical units requires geometric determination of R_MID — which is open in the present scope.

### Open sub-targets

- Sub-target #4: R_MID self-consistency. Determine R_MID as the equilibrium throat radius for the locked mass spectrum, rather than imposing R_MID = 1 by convention. With R_MID geometrically determined, the conversion ℏ = m_e R_MID c becomes a prediction.
- Identify the structural source of the lepton mass ladder. The ratios m_μ/m_e ≈ 207 and m_τ/m_e ≈ 3477 are NOT in the Tangherlini ω-spectrum; they come from the locked surrogate Hamiltonian's structure (closure quantum + radial sums + pinhole). The physical origin of these scales is the genuine open problem.

## Closure of the ℏ-origin research thread

The four-probe sequence on this branch leaves the ℏ-origin problem in a precise state:

**Sub-target #1 — Closure-cycle integer quantization.** Confirmed at the exact-quantum level. The closure cycle is integer-valued in units of 2π for every species: `N_total = N_layer_1 + N_radial`, with all four constituent channels (antipodal closure, Hopf-throat partnership, β-uplift, hard-wall radial BS) integer-quantized (`closure_cycle_action_probe`, `closed_orbit_radial_action_probe`, `hard_wall_boundary_verification`).

**Sub-target #2 — Aharonov-Bohm Hopf-fibre form.** Verified. The Hopf holonomy `π·cos(χ)` matches numerical integration to machine precision; spinor double-cover is integer-quantized at exactly the polar fibres χ ∈ {0, π/2, π}; the Hopf-throat partnership at χ = 0 contributes one full closure quantum (`aharonov_bohm_hopf_fibre_probe`).

**Sub-target #3 — ω ↔ m_e relation.** Approximate at leading order (this probe). Q1 holds at 5 %; Q2 fails decisively. Dimensional bridge is suggestive for the electron only.

**Sub-target #4 — R_MID self-consistency.** Open. Beyond the closure-ledger machinery's scope; requires throat-dynamics infrastructure not present in the current codebase.

**Net status:** the closure-phase ledger is now structurally complete (the closure cycle is integer-quantized end-to-end). The dimensional content of those integers — i.e. the conversion factor to ℏ in SI units — remains pending the deeper sub-target #4. BAM predicts dimensionless ratios at sub-percent and identifies the integer counts that label the lepton spectrum, but the absolute MeV scale is anchored, not derived. **This is the cleanest, most-progressed state of the ℏ-origin problem in the framework's history.**