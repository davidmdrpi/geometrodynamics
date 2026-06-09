# Finite-momentum charge form factor on the antipodal cavity (PR #146)

**Run:** 2026-06-09T23:53:27+00:00

Extends the #145 zero-momentum charge (F(0) = c₁, Z₁ = Z₂) to finite momentum transfer. The Bethe sum rule — current conservation at every q — protects the form factor's normalisation; the dressed charge density integrates to c₁ exactly; and the charge radius is GEOMETRIC: the bare cavity mode profile, with the one-loop cloud a ~1e-4 correction. The throat charge is spread over the cavity, finite, with no UV divergence. *(QFT on the classical throat, not quantum gravity.)*

- **Form factor**: F(q) = ∫ρ(x)e^{iq(x−x̄)}dx; F(0) = c₁ (#145 anchor)
- **Finite-q Ward**: Bethe sum rule Σ(E_m−E_n)|⟨m|e^{iqx}|n⟩|² = q² (~1e-4, q ∈ [0.5,10])
- **Z₂ cross-check**: dressed-state Z₂ = Dyson Z₂ = 0.985543 (machine precision)
- **Charge radius**: r_c = 0.2649 (tortoise units) — geometric; cloud shift ~1e-4
- **Contrast**: charge-violating cloud ⟹ total ≠ c₁; conservation is the protection (#58/#142)
- **Open**: recoil; F₁/F₂ (g−2, #62); higher loops; normalisation (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | extend the #145 q = 0 charge to finite momentum transfer | **PASS** |
| T2 | `T2_bare_charge_density_and_form_factor` | bare form factor: F(0) = c₁, monotone fall-off (spatial structure) | **PASS** |
| T3 | `T3_finite_q_ward_bethe_sum_rule` | Bethe sum rule = q² across q ∈ [0.5, 10] (finite-q Ward identity) | **PASS** |
| T4 | `T4_dressed_density_exact_charge` | dressed-state Z₂ = Dyson Z₂ (#145); total dressed charge = c₁ exact | **PASS** |
| T5 | `T5_geometric_charge_radius` | charge radius geometric: variance = small-q F; cloud shift ~1e-4 | **PASS** |
| T6 | `T6_charge_violation_counterfactual` | charge-violating cloud ⟹ total ≠ c₁ — conservation is the protection | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: derived sum rule/anchor/radius vs modelled coupling vs input α | **PASS** |
| T8 | `T8_assessment` | CHARGE_FORM_FACTOR_FINITE_Q_WARD_ANCHORED_GEOMETRIC_CHARGE_RADIUS | **PASS** |

## The finite-q Ward identity (Bethe sum rule = q²)

| q | Bethe sum | q² | ratio |
|---:|---:|---:|---:|
| 0.5 | 0.249992 | 0.25 | 0.999969 |
| 1.0 | 0.999968 | 1.0 | 0.999968 |
| 2.0 | 3.999858 | 4.0 | 0.999964 |
| 5.0 | 24.998521 | 25.0 | 0.999941 |
| 10.0 | 99.985644 | 100.0 | 0.999856 |

## The form factor falls off (the charge has spatial structure)

| q | F(q) |
|---:|---:|
| 0.0 | 1.0 |
| 0.5 | 0.991255 |
| 1.0 | 0.965387 |
| 2.0 | 0.867254 |
| 4.0 | 0.55029 |
| 8.0 | 0.028812 |

## The charge radius is geometric

| quantity | value |
|---|---:|
| r_c bare (mode profile) | 0.265 |
| r_c dressed (one loop) | 0.2649 |
| r_c from small-q F(q) | 0.2649 |
| relative cloud shift | 9e-05 |
| cavity width (tortoise) | 1.4663 |
| healing length √(2·rs·ε) | 0.2 |

The radius from the density variance matches the small-q fall-off of F(q); the one-loop cloud moves it only at the 1e-4 level — the charge radius is the bare cavity mode profile (classical geometry), not a cloud effect.

## Verdict

**CHARGE_FORM_FACTOR_FINITE_Q_WARD_ANCHORED_GEOMETRIC_CHARGE_RADIUS.** THE CHARGE FORM FACTOR ON THE ANTIPODAL CAVITY IS ANCHORED AT F(0) = c₁ (#145), PROTECTED AT EVERY MOMENTUM TRANSFER BY THE BETHE SUM RULE, AND FALLS OFF WITH A FINITE, GEOMETRIC CHARGE RADIUS — THE THROAT CHARGE IS SPREAD OVER THE CAVITY, NOT POINTLIKE, WITH NO UV DIVERGENCE. PR #145 fixed the dressed charge at q = 0; this probe supplies its q ≠ 0 structure.

THE FORM FACTOR. The dressed charge density ρ(x) integrates to c₁; F(q) = ∫ρ e^{iq(x−x̄)} dx is its transform about the centroid, with F(0) = c₁ and the small-q expansion defining the charge radius r_c² = Var_ρ(x).

THE FINITE-q WARD IDENTITY. The Bethe sum rule Σ(E_m−E_n)|⟨m|e^{iqx}|n⟩|² = q² holds at every momentum transfer (verified to ~1e-4 across q ∈ [0.5, 10]): the double commutator [e^{−iqx},[H,e^{iqx}]] = 2q² is V-independent, so current conservation constrains the cavity at EVERY q — the finite-q generalization of the #144 TRK sum rule (its q² → 0 limit) and the finite-q face of the #142 Ward identity.

THE DRESSED DENSITY. The one-loop dressed state reproduces the #145 Dyson Z₂ exactly (1/(1 + Σa²) = 1/(1 − Σ'), machine precision — the two one-loop pictures agree), and its real-space charge density integrates to c₁ EXACTLY at every coupling: the #145 anchor, now in real space.

THE CHARGE RADIUS IS GEOMETRIC. r_c from the density variance equals r_c from the small-q fall-off; the one-loop cloud moves it only at the 1e-4 level (× coupling²). The radius is the bare cavity mode profile — an O(cavity-scale) geometric length in throat units. The throat charge is not pointlike: it is spread over the cavity, finite, with no UV divergence — the form-factor face of the #55 finite self-energy — and its size is set by the CLASSICAL GEOMETRY, with the QFT dressing a small correction on top (geometry → fields, the program's arrow).

THE COUNTERFACTUAL. A cloud carrying c' ≠ c₁ shifts the total dressed charge away from c₁: the anchor — and the whole normalisation of F — rests on exact charge conservation at the unitary throat (Σc₁ = 0, #58/#141/#142/#145); the absorbing throat leaks charge and loses it.

SCOPE. One loop, radial reduction: no relativistic recoil, no F₁/F₂ decomposition (the magnetic form factor / g−2 is #62's territory). Modelled cubic coupling (the #136 posture); higher loops, the absolute normalisation (#133), and the flavor residuals (#134) stand; α(μ₀) stays the one EM input (#143).
