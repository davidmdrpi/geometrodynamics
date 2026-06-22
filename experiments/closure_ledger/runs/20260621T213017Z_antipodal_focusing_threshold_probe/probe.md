# Antipodal wave-packet focusing threshold (PR #166)

**Run:** 2026-06-21T21:30:17+00:00

Computes the antipodal focusing the THESIS asserts but never simulated: a conformal wave packet on the closed S³ reconverges EXACTLY at the antipode at t = πR, the geometric trigger for throat nucleation. *(QFT on the classical throat, not quantum gravity.)*

- **Refocus**: exact at t=πR (machine precision); full revival at 2πR
- **Conformal**: sharp focus requires conformal coupling (minimal dephases)
- **Caustic**: 1/sin²χ density gain; diffuse S³ wave → throat scale (R/R_MID)
- **Threshold**: focused energy ≥ E(R*); pair (Σc₁=0) → 2 m_e c² (PR #58)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_and_scope` | goal, framing, honest scope (computed vs inherited) | **PASS** |
| T2 | `T2_zonal_reduction` | exact S³ zonal → 1D-string reduction; 1/sin χ focusing | **PASS** |
| T3 | `T3_exact_antipodal_refocus` | exact antipodal refocus at t=πR; full revival at 2πR | **PASS** |
| T4 | `T4_conformal_coupling_required` | sharp focus requires conformal coupling (minimal dephases) | **PASS** |
| T5 | `T5_caustic_focusing_gain` | 1/sin²χ caustic; gain sharpens with ℓ_max ~ R/R_MID | **PASS** |
| T6 | `T6_nucleation_threshold` | threshold: focused energy → 2 m_e c² (pair, Σc₁=0) | **PASS** |
| T7 | `T7_honesty_and_scope` | honesty/scope: trigger computed, nonlinear throat not simulated | **PASS** |
| T8 | `T8_assessment` | ANTIPODAL_FOCUS_EXACT_AT_PI_R_TRIGGERS_NUCLEATION | **PASS** |

## The focus, quantified

| quantity | value |
|---|---:|
| refocus time | πR |
| refocus identity error | 3e-15 |
| amplitude recovery (conformal) | ×1.0000 |
| revival error at 2πR | 0e+00 |
| recovery, conformal vs minimal | ×1.000 vs ×0.877 |

## Verdict

**ANTIPODAL_FOCUS_EXACT_AT_PI_R_CAUSTIC_TRIGGERS_NUCLEATION_AT_2MEC2.** THE ANTIPODAL FOCUS IS REAL, EXACT, AND CONFORMAL — THE GEOMETRIC TRIGGER FOR PAIR NUCLEATION.

THE FOCUS. A conformal wave packet on the closed S³ does not dissipate: it refocuses EXACTLY at the antipode at t = πR (half the great-circle period), to machine precision (3e-15), amplitude recovered to ×1.0000; at t = 2πR it fully revives (the sub-threshold focus passes through and re-disperses). The zonal sector reduces exactly to a 1D string, and the physical field ψ = f/sin χ carries the geometric focusing factor.

CONFORMAL REQUIRED. The sharp focus needs the conformal tower ω_ℓ = (ℓ+1)/R (recovery ×1.0000); the minimally-coupled tower dephases (×0.8771). The same conformal coupling that makes the S³ vacuum tower equally spaced (PR #165) makes the antipodal caustic sharp.

THE CAUSTIC. The energy density ∝ 1/sin²χ diverges as the wavefront converges on the antipode — a true caustic, regularized by ℓ_max ~ R/R_MID. It lets a delocalized, S³-wide wave reconcentrate onto the throat scale: the dynamical bridge from a diffuse wave to a local nucleation density.

THE THRESHOLD. The focus is the trigger; nucleation requires the focused energy reach the lowest stable throat E(R*) = m_e c², and the C-conjugate pair (Σc₁=0) makes the threshold 2 m_e c² = 1.022 MeV (inherited, PR #58). The disperse-below / persist-above dichotomy is the bubble barrier; the nonlinear throat formation is named, not simulated.
