# Compton-scattering antipodal-transaction kinematics probe

**Run:** 2026-05-16T01:37:16+00:00

Opens the QFT-event-reinterpretation thread (`docs/qft_event_reinterpretation_research_plan.md`). Tests whether BAM's antipodal-closure picture is kinematically consistent with standard Compton scattering `γ + e → γ + e`, and identifies which proper-time skew (inter-mouth, throat-pinch, closure-cycle) carries the BAM-specific dynamical content.

**Reinterpretation:** Local QED 4-vertex (γ + e → γ + e) ↦ antipodal transaction with front-vertex at X and back-vertex at the `S³` antipode X̃, linked by throat worldlines. Antipodal map on 4-momentum: (E, p) → (E, −p).

**Units:** c = ℏ = m_e = 1; R_S3 = 1

## Predictions and results

| # | Prediction | Key metric | Value | PASS? |
|---|---|---|---:|---|
| P1 | `P1_closure_no_overconstraint` | max back-vertex residual | 1.67e-16 | **PASS** |
| P2 | `P2_inter_mouth_skew_vanishes` | max inter-mouth γ skew | 0.00e+00 | **PASS** |
| P3 | `P3_throat_pinch_skew` | soft-ω power-law exponent at θ=π | 1.885 | **PASS** |
| P4 | `P4_crossing_via_antipodal` | max u-invariance residual | 0.00e+00 | **PASS** |

## P1_closure_no_overconstraint

For every (ω, θ), the back-vertex 4-momentum balance (under antipodal flip of all 4-momenta) holds whenever the front vertex does. Tests that BAM closure does not add constraints beyond standard QED Compton kinematics.

Sample row (first of 63):

```
  omega: 0.01
  theta: 0
  front_residual_inf: 0
  back_residual_inf: 8.67362e-18
  stu_residual: 0
```

## P2_inter_mouth_skew_vanishes

The proper-time skew between the two mouths of the same particle vanishes identically: γ_front = γ_back because the antipodal map flips 3-momentum while preserving energy, so the boost factor γ = E/m is identical at both mouths.

**Interpretation.** Inter-mouth skew is NOT the BAM-specific dynamical variable — it is identically zero by the antipodal momentum map. Throat-pinch skew (P3) is the next candidate.

Sample row (first of 5):

```
  omega: 0.01
  theta: 0
  particle: p_e_in
  gamma_front: 1
  gamma_back: 1
  inter_mouth_skew: 0
```

## P3_throat_pinch_skew

Throat-pinch proper-time skew Δτ_throat = (π·R_S3/c)·(1/γ_out − 1/γ_in) — the integrated information-propagation delay between front and back mouth during the scattering, and the BAM-specific dynamical quantity. Tests its scaling with photon energy and identifies whether it is structurally a recoil-induced quantity or a topological invariant.

**Interpretation.** For ω → 0 (Thomson limit), Δτ_throat → 0 quadratically in ω (electron recoil is `O(ω²/m_e²)` to leading non-trivial order). This identifies Δτ_throat as a *recoil-induced* dynamical quantity, not a topological invariant. The closure-quantum reading is ruled out at this level; the BAM-specific content here is the specific propagation-delay formula tied to the `S³` antipodal traverse time π·R_S3/c.

Sample row (first of 63):

```
  omega: 0.01
  theta: 0
  omega_over_me: 0.01
  gamma_e_out: 1
  delta_tau_throat: 0
  delta_tau_over_omega_over_me: 0
```

## P4_crossing_via_antipodal

Tests that the BAM antipodal map on the outgoing-photon 4-momentum leaves the Mandelstam-u invariant in the electron rest frame. This is the relevant statement for s↔u crossing symmetry in the BAM antipodal-transaction reinterpretation.

**Interpretation.** The antipodal map on a single 4-momentum (E, p) → (E, −p) preserves p · q for any 4-momentum q with q_spatial = 0 (i.e. for a rest-frame particle). In the electron rest frame, all Mandelstam invariants involving p_e_in are invariant under the antipodal map on any other 4-momentum. This makes the BAM antipodal action structurally compatible with crossing, but it does NOT directly implement crossing (which would require a full p → −p map across all particles, not just the outgoing photon).

Sample row (first of 5):

```
  omega: 0.01
  theta: 0
  s: 1.02
  u_original: 0.98
  u_after_antipodal: 0.98
  u_invariance_residual: 0
```

## Verdict

**PASS — BAM antipodal-closure picture is kinematically consistent with standard Compton scattering. The predictions are: (P1) closure adds no constraints beyond the standard ones; (P2) inter-mouth proper-time skew vanishes identically; (P3) throat-pinch skew Δτ_throat = (π·R_S3/c)·(1/γ_out − 1/γ_in) is the BAM-specific dynamical quantity, scaling as ω² at small recoil (Thomson limit) — this is a recoil-induced effect tied to the `S³` antipodal traverse time, not a topological closure quantum; (P4) the antipodal map on the outgoing photon preserves the Mandelstam-u invariant in the electron rest frame, structurally compatible with crossing symmetry though not implementing it directly. The thread is cleared to proceed to the amplitude / cross-section probe.**

## What this leaves open

- **Amplitude algebra.** This probe is kinematic only. The cross-section probe (Klein-Nishina at tree level, Thomson in the soft-photon limit) requires extending the existing `transaction/handshake.py` amplitude machinery from GW-triggered transactions to QED scattering. That is the designated follow-on probe.
- **Off-rest-frame antipodal action.** The P4 invariance is established only in the electron rest frame. The general covariance question — does the antipodal action commute with Lorentz boosts in any natural sense? — is open.
- **Photon throat-pair structure.** The probe treats the photon as a null 4-momentum without explicit throat-pair representation (the boost factor diverges for photons). BAM's photon structure (Hopf-fibre excitation? massless throat pair?) needs explicit modelling for the amplitude probe.
- **Throat-pinch skew physical interpretation.** P3 shows Δτ_throat scales as `ω²/m_e²` at low energies, identifying it as a recoil-induced quantity tied to the `S³` antipodal traverse time. Whether this quantity enters a measurable observable (a phase shift in the amplitude, a delay in the transaction completion) is the next physical question.
