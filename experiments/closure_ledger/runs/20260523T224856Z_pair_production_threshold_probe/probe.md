# Pair-production / throat-nucleation threshold probe

**Run:** 2026-05-23T22:48:56+00:00

Derives the pair-production threshold as 2× the lowest stable throat configuration, forced into a C-conjugate pair by Hopf-charge / antipodal-Z₂ conservation.

- **Threshold**: 2 m_e c² = 1.022 MeV
- **Mechanism**: C-conjugate throat pair (Σc₁=0); lowest stable R*
- **Dynamics**: nucleation barrier (disperse/persist); Schwinger field
- **B4 caveat**: factor 2 derived; absolute scale = single anchor m_e c² = ℏc/R_MID

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_lowest_stable_configuration` | E(R*) global min = ground-state throat (electron) | **PASS** |
| T2 | `T2_pair_from_conservation` | |c₁|=1; C-conjugate pair Σc₁=0; single forbidden | **PASS** |
| T3 | `T3_threshold_2mc2` | E_thr = 2 m_e c² = 1.0220 MeV | **PASS** |
| T4 | `T4_nucleation_barrier` | R_c=2σ/ρ; disperse below / persist above | **PASS** |
| T5 | `T5_schwinger_critical_field` | E_S=1.323e+18 V/m; e E_S R_MID=m_e c² | **PASS** |
| T6 | `T6_subthreshold_dispersal` | below 2 m_e c² → disperses; above → persists | **PASS** |
| T7 | `T7_b4_accounting` | factor 2 derived; scale = single anchor | **PASS** |
| T8 | `T8_assessment` | threshold = 2× lowest stable throat | **PASS** |

## T1: Lowest stable configuration

- R* = 3.8616e-13 m (numeric min 3.8616e-13, rel err 8.0e-06)
- E(R*) = 4.4808e-16 J; stable minimum: True

## T2: Charge/topology conservation → pairs

- throat Hopf charge |c₁| = 1
- orientations: -1.0000, +1.0000 (throat / antithroat)
- pair charge sum Σc₁ = 0; vacuum = 0
- single-throat creation allowed: False; pair creation allowed: True

## T3: Threshold = 2 m_e c²

- single throat rest energy = 0.510999 MeV
- pair threshold = 2 m_e c² = 1.021998 MeV (expected 1.021998)

## T4: Nucleation barrier (disperse vs persist)

- σ = 1.0, ρ = 1.5; critical radius R_c = 2σ/ρ = 1.3333
- below R_c: dE/dR > 0 (disperse): True
- above R_c: dE/dR < 0 (persist): True
- barrier height E_c = (16π/3)σ³/ρ² = 7.4467

## T5: Schwinger critical field

- E_S = m_e²c³/(eℏ) = 1.3233e+18 V/m (expected 1.3233e+18)
- e E_S · R_MID / (m_e c²) = 1.000000 (R_MID = λ_C = 3.8616e-13 m)

## T6: Sub-threshold dispersal

| E_in (MeV) | above threshold? | outcome |
|---:|:---:|---|
| 0.300 | False | sub-critical: disperses to vacuum |
| 0.600 | False | sub-critical: disperses to vacuum |
| 1.000 | False | sub-critical: disperses to vacuum |
| 1.022 | True | real pair persists |
| 1.500 | True | real pair persists |
| 2.000 | True | real pair persists |

## T7: B4 accounting

- m_e c² from bridge ℏc/R_MID = 8.1871e-14 J (actual 8.1871e-14); consistent: True
- derived/dimensionless: {'pair_factor': 2, 'charge_sum': 0, 'schwinger_relation': 'e E_S R_MID / (m_e c²) = 1'}
- absolute scale: single anchor m_e c² = ℏc/R_MID

## T8: Assessment

- **Threshold**: 2 m_e c² = 1.022 MeV
- **Lowest stable config**: ground-state throat R* (electron)
- **Pair mechanism**: C-conjugate (Σc₁=0); single-throat creation forbidden
- **Dynamics**: nucleation barrier R_c=2σ/ρ; disperse below / persist above
- **Field connection**: Schwinger E_S: e E_S R_MID = m_e c²
- **Remaining**: full instanton rate; heavier-lepton thresholds 2m_μc², 2m_τc²

## Verdict

**PAIR_THRESHOLD_DERIVED.** PAIR THRESHOLD DERIVED. The pair-production threshold falls out as 2× the lowest stable throat configuration.

LOWEST STABLE CONFIGURATION. The self-energy E(R)=A/R+B·R² (EM repulsion + brane tension, PRs #55–#57) has a unique global minimum at R*=(A/2B)^(1/3); the ground-state throat has rest energy E(R*) = m_e c² (the single anchor) — this is the electron.

PAIRS FROM CONSERVATION. A throat carries one unit of Hopf charge (|c₁|=1, compute_c1); the vacuum has c₁=0. The two Hopf orientations (±1) are the throat and its C-conjugate antithroat (the inner/outer swap / antipodal Z₂). Conservation Σc₁=0 forbids single-throat creation and forces C-conjugate PAIR creation.

THRESHOLD. By C-symmetry both partners have rest energy m_e c², so the pair-production threshold is E_thr = 2 m_e c² = 1.022 MeV — the minimal on-shell pair.

DYNAMICS. Creating a throat from the vacuum is bubble nucleation, E_nuc(R)=4πσR²−(4/3)πρR³, with critical radius R_c=2σ/ρ: below R_c the configuration shrinks (the antipodal focus disperses → vacuum), above R_c it grows to the stable R* (a real throat persists) — exactly the THESIS dichotomy. The Schwinger critical field E_S=m_e²c³/(eℏ)≈1.32×10¹⁸ V/m is where the field work over a reduced Compton wavelength (= R_MID) equals m_e c²: e E_S·R_MID = m_e c², tying the throat scale to the threshold.

B4. The factor 2, the pair structure (Σc₁=0), the disperse/persist dichotomy, and e E_S R_MID = m_e c² are all derived / dimensionless; the absolute threshold 2 m_e c² rides on the single anchor m_e c² = ℏc/R_MID, consistent with the B4 scale-modulus theorem. Remaining: the full instanton/tunneling rate and the heavier-lepton thresholds 2 m_μ c², 2 m_τ c² (the excited radial throats).

## What this leaves open

- **The anchor scale.** m_e c² = ℏc/R_MID still rides on the one dimensionful input (the bulk gravitational scale, PR #57).
- **The full dynamical nucleation rate.** The static critical-bubble barrier here; the tunneling/instanton rate (Schwinger exponential e^{−πm²c³/(eEℏ)}) from the BAM throat action is the follow-on.
- **Heavier-lepton thresholds.** 2 m_μ c², 2 m_τ c² from the excited radial throats (the closure ladder).
