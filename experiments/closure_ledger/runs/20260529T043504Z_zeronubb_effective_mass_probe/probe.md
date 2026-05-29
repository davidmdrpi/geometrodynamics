# Neutrinoless double-beta (0νββ) effective-mass probe (PR #95)

**Run:** 2026-05-29T04:35:04+00:00

Turns the neutrino arc (#85–#94) into a concrete, falsifiable prediction for the effective Majorana mass `m_ββ = |Σ U_ei² m_i|` measured in 0νββ. **0νββ occurs** (the neutrino is Majorana, `c₁=0`, PR #86); BAM selects the **normal-ordering** band (PR #91), populated by the **anarchic Majorana phases** (PR #94); and at the **light scale** (PR #90) `m_ββ ≲ 8 meV` — below current reach and a falsifiable target for next-gen experiments.

- **Identification**: 0νββ occurs (Majorana, c₁=0); BAM normal-ordering band populated by anarchic Majorana phases; at the light scale m_ββ ≲ 8 meV — below current reach, falsifiable
- **Occurs because**: neutrino Majorana ⟸ c₁=0 (PR #86)
- **Ordering**: normal (PR #91) — NO band, below IO floor ~19 meV
- **Phases**: anarchic (PR #94) — full band incl. cancellation to ~0
- **Scale**: light, ~few meV (PR #90) — m_ββ ≲ 8 meV
- **Falsifier**: discovery with m_ββ ≳ 19 meV ⟹ IO/degenerate

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_setup` | m_ββ = |Σ U_ei² m_i|; needs Majorana+ordering+phases+scale | **PASS** |
| T2 | `T2_zeronubb_occurs_majorana` | 0νββ occurs ⟸ neutrino Majorana ⟸ c₁=0 (PR #86) | **PASS** |
| T3 | `T3_normal_ordering_band` | normal ordering ⟹ NO band (below IO floor ~19 meV) | **PASS** |
| T4 | `T4_anarchic_phases_populate_band` | anarchic phases populate full band incl. cancellation→~0 | **PASS** |
| T5 | `T5_bam_light_scale_mbb` | BAM light scale ⟹ m_ββ ≲ 8 meV | **PASS** |
| T6 | `T6_experimental_comparison_falsifiable` | below current (28–122 meV) & next-gen (~9–20); falsifiable | **PASS** |
| T7 | `T7_honest_scope` | qualitative firm; exact m_ββ a band (residuals) | **PASS** |
| T8 | `T8_assessment` | ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV | **PASS** |

## m_ββ band (normal ordering)

| m_lightest (meV) | m_ββ min (meV) | m_ββ max (meV) |
|---:|---:|---:|
| 0 | 1.45 | 3.68 |
| 2 | 0.15 | 5.11 |
| 3 | 0.02 | 5.88 |
| 5 | 0.02 | 7.49 |

The anarchic phases populate the whole band, with a cancellation trough down to `0.02 meV` (m_ββ → ~0) around m_lightest ~ 3–5 meV.

## BAM (normal) vs inverted ordering, and experiment

| | m_ββ (meV) |
|---|---|
| **BAM (normal, light scale)** | ≲ 7 (with cancellation to ~0) |
| inverted-ordering band | 19–48 (floor ~19, no cancellation) |
| current bound (KamLAND-Zen) | 28–122 |
| next-gen reach (LEGEND-1000 / nEXO) | ~9–20 |

BAM predicts `m_ββ ≲ 8 meV` — below current bounds (so the null result is expected) and largely below next-gen reach. **Falsifier:** a 0νββ discovery with `m_ββ ≳ 19 meV` (the IO floor) would imply inverted ordering or a quasi-degenerate scale, contradicting BAM.

## Verdict

**ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV.** 0νββ OCCURS, IN NORMAL ORDERING, WITH m_ββ ≲ 8 meV — BELOW CURRENT REACH AND FALSIFIABLE. The neutrino arc (#85–#94) fixed the structure of the BAM neutrino sector; this probe turns it into a prediction for the effective Majorana mass m_ββ = |Σ U_ei² m_i| measured in 0νββ.

0νββ OCCURS. 0νββ requires a Majorana neutrino (ΔL=2). In BAM the neutrino is Majorana because it is chargeless (c₁=0, C-invariant, PR #86), so 0νββ occurs — a qualitative prediction. A Dirac neutrino would forbid it (m_ββ ≡ 0).

NORMAL ORDERING ⟹ THE NO BAND. BAM predicts normal ordering (PR #91: generations = cavity overtones, m_ν ∝ m_D), selecting the NO m_ββ band — m_ββ ≈ 1.5–3.7 meV at zero lightest mass. The inverted-ordering band (the contrast) sits at ~19–48 meV with a hard floor ~19 meV and no cancellation to zero — entirely above the BAM band.

ANARCHIC PHASES POPULATE THE BAND. The anarchic Majorana phases (PR #94, uniform) populate the WHOLE NO band, including the cancellation trough where the three terms partially cancel and m_ββ → ~0 (around m_lightest ~ 3–5 meV).

LIGHT SCALE ⟹ m_ββ ≲ 8 meV. At the BAM light scale (PR #90: lightest ~ few meV) the NO band is m_ββ ≈ 0–8 meV. This is below the current bound (KamLAND-Zen, m_ββ ≲ 28–122 meV) — consistent with the null result — and largely below even next-generation reach (LEGEND-1000 / nEXO, ~9–20 meV). It is a sharp falsifier: a 0νββ discovery with m_ββ ≳ 19 meV would imply inverted ordering or a quasi-degenerate scale, contradicting the BAM normal-ordering + light-scale prediction.

HONEST SCOPE. ESTABLISHED (BAM-native): 0νββ occurs (the neutrino is Majorana, PR #86); BAM selects the normal-ordering band (PR #91), below the IO floor; the anarchic Majorana phases (PR #94) populate the full band (cancellation to ~0 allowed); at the BAM light scale (PR #90) m_ββ ≲ 8 meV — below current bounds and a target/falsifier for next-gen experiments. NOT established: a single m_ββ value — it is a band, because the lightest neutrino mass is unmeasured and the Majorana phases are anarchic (uniform); the exact spectrum (the PR #91 χ_n-corrected ratios, the absolute scale) and the specific phases are the residuals.

## What this leaves open

- **A single m_ββ value** — it is a band, because the lightest neutrino mass is unmeasured and the Majorana phases are anarchic (uniform).
- **The exact spectrum** (PR #91 `χ_n`-corrected ratios; the absolute scale) and the **specific phases**.
