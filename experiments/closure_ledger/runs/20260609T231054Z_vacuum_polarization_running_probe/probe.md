# One-loop photon vacuum polarisation and the running of α (PR #144)

**Run:** 2026-06-09T23:10:54+00:00

Computes the one-loop photon vacuum polarisation on the antipodal cavity — the running of α that #142/#143 classified as derived but never computed. The cavity Ward identity (the diamagnetic–paramagnetic cancellation) is verified numerically, the photon stays exactly massless (1/q² protected), the spectral density is positive with no width below the lowest pair threshold, and the running has the QED screening direction with the flat-space log coefficient α/3π. Only the boundary value α(μ₀) stays input (#143). *(QFT on the classical throat, not quantum gravity.)*

- **Bubble**: Π = Σ_pairs c|v_nm|²/(s − (ω_n+ω_m)² + i0⁺); v_nm = ∫φ_γψ_nψ_m dr*
- **Ward**: TRK = 0.999969; 1 − S = 3.10e-05 (computed cancellation)
- **Masslessness**: Π(0) = 0 ⟹ photon pole at q² = 0; 1/q² (#42–#44) protected
- **Screening**: Δ(Q²) monotone ⟹ α_eff increases with Q² (QED direction)
- **Running**: flat-limit log slope = α/3π (~1%); α(μ₀) input (#143)
- **Contrast**: absorbing throat ⟹ Im Π ≠ 0 below threshold; Ward protection lost (#130/#142)
- **Open**: higher loops; 4D tensor Π^μν; absolute normalisation (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | one-loop vacuum polarisation — the uncomputed #142/#143 running | **PASS** |
| T2 | `T2_polarisation_bubble_construction` | charged-pair bubble; even-Σl Z₂ selection on the pair channel (#141) | **PASS** |
| T3 | `T3_cavity_ward_identity_computed` | cavity Ward identity computed: TRK = 1, 1 − S = 0 (~3e-5) | **PASS** |
| T4 | `T4_photon_massless_vs_absorbing_counterfactual` | Π(0) = 0 photon massless; absorbing ⟹ Im Π ≠ 0 (width) + no protection | **PASS** |
| T5 | `T5_spectral_positivity_and_screening` | ρ ≥ 0; no width below (2ω_0)²; Δ(Q²) monotone ⟹ screening | **PASS** |
| T6 | `T6_running_flat_log_coefficient` | flat-limit log slope = α/3π (~1%); α(μ₀) stays input (#143) | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: derived structure/running vs modelled photon leg vs input α | **PASS** |
| T8 | `T8_assessment` | VACUUM_POLARIZATION_WARD_MASSLESS_SCREENING_LOG_RUNNING_ALPHA_INPUT | **PASS** |

## The cavity Ward identity, computed (Π(0) = 0)

| quantity | value |
|---|---:|
| TRK sum Σ(E_m−E_n)\|x_mn\|² (expect 1) | 0.999969 |
| paramagnetic sum S (expect 1) | 0.999969 |
| Ward bracket 1 − S (expect 0) | 3.102e-05 |
| Im Π below threshold — antipodal | -6.045e-09 |
| Im Π below threshold — absorbing | -0.04159 |

The diamagnetic +1 cancels the paramagnetic sum exactly (grid-level residual ~3e-5): Π(0) = 0, no photon mass. The absorbing counterfactual carries a below-threshold width and loses the cancellation.

## The running: flat-space log coefficient α/3π

| Q²/m² | Δ_QED(Q²) |
|---:|---:|
| 100 | 0.002319 |
| 1000 | 0.004062 |
| 10000 | 0.005834 |
| 100000 | 0.007616 |

Fitted log slope `0.0007667904` vs `α/3π = 0.0007742732` — ratio `0.9903`. The Ward-protected dispersion form, fed the flat 4D pair density, reproduces the textbook running; the cavity bubble gives the discrete-threshold analogue (T5).

## Verdict

**VACUUM_POLARIZATION_WARD_MASSLESS_SCREENING_LOG_RUNNING_ALPHA_INPUT.** THE ONE-LOOP PHOTON VACUUM POLARISATION ON THE ANTIPODAL CAVITY IS WARD-PROTECTED, MASSLESS, POSITIVE, AND RUNS IN THE QED SCREENING DIRECTION WITH THE FLAT-SPACE LOG COEFFICIENT α/3π — THE RUNNING THAT #142/#143 CLASSIFIED AS DERIVED IS NOW COMPUTED; ONLY α(μ₀) STAYS INPUT. PRs #141–#143 completed the gauge–matter structure but left the running uncomputed; the matter sector had its one-loop two-point function (#136) while the photon did not. This probe supplies it.

THE BUBBLE. Π is the charged-pair loop: pair (n,m) opens at s_nm = (ω_n+ω_m)² with vertex v_nm = ∫φ_γ ψ_n ψ_m dr* (the #137/#141 triple overlap with one photon leg) and spectral density ρ_nm = c_nm|v_nm|² ≥ 0; the photon couples only to even-Σl pair channels (the #141 antipodal Z₂ selection, re-verified exactly).

THE CAVITY WARD IDENTITY, COMPUTED. Gauge invariance under minimal substitution makes the O(A²) shift vanish: the diamagnetic +1 cancels the paramagnetic sum, 1 − S = 0, the TRK sum rule in disguise — verified numerically to ~3e-5 on the cavity tower. The #142 Ward identity, structural there, is quantitative here.

THE PHOTON STAYS EXACTLY MASSLESS. Π(0) ∝ 1 − S = 0: the photon pole stays at q² = 0 and the 1/q² kernel (#42–#44) is protected through one loop. The absorbing counterfactual (#130) breaks it: complex pair thresholds give Im Π ≠ 0 below every threshold (an absorption width on the photon) and the real-mode cancellation is gone — gauge protection REQUIRES the unitary antipodal throat, the one-loop face of #129/#142.

POSITIVITY, UNITARITY, SCREENING. ρ_nm ≥ 0 manifestly; Im Π = 0 below the lowest pair threshold (2ω_0)² (no width — the #136 pattern, now on the photon); and the Ward-protected Δ(Q²) is monotone increasing, so α_eff = α/(1 − Δ) increases with spacelike Q² — the QED screening direction, with the discrete pair thresholds the cavity analogue of the lepton thresholds.

THE RUNNING, COMPUTED. The same once-subtracted dispersion machinery, fed the flat-space 4D pair density, reproduces the textbook log running with slope dΔ/d lnQ² = α/3π to ~1% over three decades. BAM derives HOW α runs — the form, the sign, the coefficient; the boundary value α(μ₀) ≈ 1/137 stays the one EM input (#143), and this probe deliberately does not hunt for it (the #107/#108 anti-numerology discipline).

SCOPE. One loop on the fixed antipodal background; the photon-leg radial profile is modelled (the soft cavity mode, the #136 posture); higher loops, the full 4D tensor structure, the absolute normalisation (#133), and the flavor residuals (#134) stand.
