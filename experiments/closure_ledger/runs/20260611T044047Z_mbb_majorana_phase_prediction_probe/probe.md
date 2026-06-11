# Majorana phase and m_ββ prediction from the PMNS flavor ensemble (PR #154)

**Run:** 2026-06-11T04:40:47+00:00

Completes the neutrino-sector prediction card: the effective Majorana mass extracted exactly (m_ββ = |M_fl,ee|) from the derived #151–#153 flavor ensemble. The prediction is a few meV, EXACTLY independent of the charged-side rotation (the e-row argument, machine zero), with generic Majorana phases and a negligible lightest-state contribution — the program's earlier ≲ 8 meV claim sharpened into a falsifiable distribution. *(QFT on the classical throat, not quantum gravity.)*

- **Shortcut**: m_ββ = |(W M W^T)_ee| exact; Takagi Σt_i cross-check ~1e-12
- **Invariance**: m_ββ exactly φ_ℓ-independent (e-row argument; machine zero)
- **Prediction**: median ≈ 3 meV; 68% [1.5, 5.9]; 95% [0.5, 8.7] (self-consistent)
- **Falsification**: P(>10 meV) ~ few %; P(>20 meV) = 0 — detection above ~10 meV falsifies
- **Phases**: Φ₂₃ broad (P(>π/2) ≈ 69%) — generic Majorana CP; m₁ negligible
- **Open**: sharpen O_geom e-row; CKM analogue; joint neutrino-sector test

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | supply #153's Majorana/m_ββ open item from the same ensemble | **PASS** |
| T2 | `T2_construction_takagi_consistency` | m_ββ = |M_fl,ee| exact; Takagi Σt_i consistency ~1e-12 | **PASS** |
| T3 | `T3_exact_phi_ell_invariance` | m_ββ exactly φ_ℓ-invariant (machine zero) — e-row argument | **PASS** |
| T4 | `T4_mbb_prediction` | median ≈ 3 meV; schemes agree; conditioning robust | **PASS** |
| T5 | `T5_falsification_structure` | P(>20 meV) = 0; detection above ~10 meV falsifies | **PASS** |
| T6 | `T6_majorana_phase_structure` | Φ₂₃ broad (generic Majorana CP); m₁ contribution negligible | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: exact structure + statistical distribution; budget unchanged | **PASS** |
| T8 | `T8_assessment` | MBB_FEW_MEV_PHI_ELL_EXACT_INVARIANT_GENERIC_MAJORANA_PHASES | **PASS** |

## The m_ββ prediction

| scheme | median (meV) | 68% | 95% |
|---|---:|---|---|
| self-consistent (draw masses → m₃) | 3.2 | [1.45, 5.85] | [0.52, 8.72] |
| data-anchored (m₂, m₃ from splittings) | 2.93 | [1.26, 5.99] | [0.49, 11.3] |
| conditioned (r₃₂ ∈ [4, 8], n = 690) | 3.07 | — | — |

## The falsification card

| probability | value |
|---|---:|
| P(m_ββ < 1 meV) — cancellation | 7.9% |
| P(m_ββ > 6 meV) | 14.9% |
| P(m_ββ > 10 meV) | 0.5% |
| P(m_ββ > 20 meV) | 0.0% |

A detection above ~10 meV would falsify the ensemble; ton-scale experiments are predicted to see nothing or a floor-level signal.

## The Majorana-phase structure

Two-term interference |t₂ + t₃| (m₁ median 0.074 meV — negligible); relative phase broad: P(|Φ₂₃| > π/2) = 69% — generic Majorana CP.

## Verdict

**MBB_FEW_MEV_PHI_ELL_EXACT_INVARIANT_GENERIC_MAJORANA_PHASES.** m_ββ IS PREDICTED AT A FEW meV — MEDIAN 3.2 meV, 68% [1.45, 5.85], 95% [0.52, 8.72] — EXACTLY INDEPENDENT OF THE CHARGED-SIDE ROTATION, WITH GENERIC MAJORANA PHASES AND A NEGLIGIBLE LIGHTEST-STATE CONTRIBUTION. #153 assembled the PMNS and left the Majorana sector open; this probe completes the neutrino-sector prediction card.

THE EXACT SHORTCUT. m_ββ = |(W M W^T)_ee| — the (e,e) element of the flavor-basis Majorana matrix: no mixing-matrix approximation, all phases included. The Takagi decomposition reproduces it term by term (Σt_i = M_fl,ee to ~1e-12).

THE EXACT INVARIANCE. The charged-side μ–τ rotation never touches the e-row, and M_fl,ee depends only on the e-row of W: m_ββ is EXACTLY φ_ℓ-independent (machine zero across 3000 draws). The one modelled O(1) angle in the #153 assembly drops out entirely — m_ββ is more robust than the mixing angles, inheriting only the geometric e-row overlap and the derived ensemble.

THE PREDICTION. Self-consistent and data-anchored normalizations agree (medians 3.2 vs 2.93 meV); conditioning on data-compatible spreads gives 3.07 meV — robust. The few-meV scale is structural: the light m₁ (#151/#152) makes m_ββ a two-term interference |t₂ + t₃| of comparable terms.

THE FALSIFICATION CARD. P(m_ββ < 1 meV) = 7.9% (anarchic-phase cancellation uncommon); P(> 10 meV) = 0.5%; P(> 20 meV) = 0.0% — a detection above ~10 meV would falsify the ensemble, and ton-scale experiments are predicted to see nothing or a floor-level signal. The earlier ≲ 8 meV program claim is sharpened into a distribution.

GENERIC MAJORANA PHASES. P(|Φ₂₃| > π/2) = 69% — the relative Majorana phase is broad, the Majorana-sector face of the #153 generic Dirac CP; no alignment predicted or needed; m₁ (median 0.074 meV) contributes negligibly.

THE NEUTRINO-SECTOR CARD, COMPLETE. Normal ordering (#113), m₁ ≈ 0.05 meV and Σm_ν ≈ 58.8 meV (#151/#152), angles anarchy-natural (#153), CP generic — Dirac (#153) and Majorana (here) — and m_ββ ≈ 1.5–6 meV (68%). SCOPE: the O_geom e-row is the one systematic; NMEs are experimental overlay; no new input (#150 budget unchanged).
