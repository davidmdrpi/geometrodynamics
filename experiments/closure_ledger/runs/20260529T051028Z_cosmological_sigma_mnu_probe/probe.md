# Cosmological Σm_ν prediction probe (PR #96)

**Run:** 2026-05-29T05:10:28+00:00

The same BAM neutrino spectrum behind the 0νββ prediction (PR #95) also fixes the sum of neutrino masses `Σm_ν = m1 + m2 + m3` — the quantity probed by cosmology (CMB lensing + LSS/BAO). BAM fixes the **normal ordering** (PR #91) and a **light scale** (PR #90), giving `Σm_ν ≈ 59–65 meV` — pinned at the normal-ordering floor, right at the DESI DR2 + CMB frontier.

- **Identification**: BAM light, normal-ordered spectrum ⟹ Σm_ν ≈ 59–65 meV (at the NO floor, below IO ~99 meV); consistent with Planck, at the DESI-DR2 + CMB frontier; falsifiable
- **NO floor**: 58.7 meV (IO floor 99.2 meV)
- **BAM band**: 59–65 meV
- **Cosmo frontier**: DESI DR2 + CMB ~60–64 meV
- **Consistency**: same spectrum as the PR #95 0νββ prediction (m_ββ ≲ 8 meV)
- **Falsifiers**: Σm_ν < 58.7 meV ⟹ NO excluded; Σm_ν ≳ 100 meV ⟹ not light

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_setup` | Σm_ν = m1+m2+m3; needs ordering + scale + Δm² | **PASS** |
| T2 | `T2_normal_ordering_floor` | normal ordering ⟹ NO floor ≈ 58.7 meV (IO floor ≈ 99) | **PASS** |
| T3 | `T3_light_scale_sigma` | light scale ⟹ Σm_ν ≈ 59–65 meV (near floor) | **PASS** |
| T4 | `T4_cosmological_comparison` | below Planck (<120), at DESI DR2+CMB frontier (~60–64) | **PASS** |
| T5 | `T5_falsifiability` | falsifiable: <58.7 ⟹ NO excluded; ≳100 ⟹ not light | **PASS** |
| T6 | `T6_consistency_with_0nubb` | one spectrum: Σm_ν + m_ββ (PR #95) cross-checkable | **PASS** |
| T7 | `T7_honest_scope` | pinned near NO floor; exact value the absolute-scale residual | **PASS** |
| T8 | `T8_assessment` | SIGMA_MNU_AT_NORMAL_ORDERING_FLOOR_60MEV | **PASS** |

## Σm_ν vs the lightest mass (normal ordering)

| m_lightest (meV) | Σm_ν (meV) |
|---:|---:|
| 0 | 58.7 |
| 2 | 60.9 |
| 5 | 65.2 |
| 10 | 74.2 |

## BAM vs cosmology

| | Σm_ν (meV) |
|---|---|
| **BAM (normal, light scale)** | ≈ 59–65 (NO floor 59) |
| inverted-ordering floor (contrast) | ≈ 99 |
| Planck 2018 + BAO (95%) | < 120 |
| DESI DR1 + CMB (95%) | < 72 |
| DESI DR2 + CMB (95%, tightest) | ≲ 62 |

BAM sits **at the normal-ordering floor**, below the IO floor, consistent with Planck and right at the DESI DR2 + CMB frontier. **Falsifier:** a robust Σm_ν below the NO floor (~58.7 meV) excludes normal ordering — contradicting BAM (and the oscillation Δm² themselves); a quasi-degenerate Σm_ν ≳ 100 meV contradicts the light scale.

## Verdict

**SIGMA_MNU_AT_NORMAL_ORDERING_FLOOR_60MEV.** Σm_ν ≈ 59–65 meV — AT THE NORMAL-ORDERING FLOOR, AT THE COSMOLOGICAL FRONTIER. The same BAM neutrino spectrum behind the 0νββ prediction (PR #95) also fixes the sum of neutrino masses Σm_ν = m1 + m2 + m3 — the quantity probed by cosmology (CMB lensing + large-scale structure / BAO).

NORMAL ORDERING ⟹ THE NO FLOOR. BAM predicts normal ordering (PR #91: generations = cavity overtones, m_ν ∝ m_D), so the floor (lightest mass → 0) is Σm_ν = √Δm²_21 + √Δm²_31 ≈ 8.7 + 50 ≈ 58.7 meV. The inverted-ordering floor (the contrast) is ≈ 99 meV — well above.

LIGHT SCALE ⟹ Σm_ν ≈ 59–65 meV. At the BAM light scale (PR #90, lightest ~ few meV) the sum rises only slightly above the floor: Σm_ν ≈ 59 meV (m_lightest = 0), 61 meV (2 meV), 65 meV (5 meV). So BAM predicts Σm_ν ≈ 59–65 meV — pinned near the NO floor, not in the quasi-degenerate regime.

COSMOLOGICAL COMPARISON. This sits exactly where cosmology is now probing: comfortably below Planck 2018 + BAO (< 120 meV), just inside DESI DR1 + CMB (< 72 meV), and right at the DESI DR2 + CMB frontier (~60–64 meV). It is a sharp, near-term falsifier: if cosmology robustly pushes Σm_ν below ~58.7 meV (the NO floor), normal ordering is excluded and the BAM prediction fails — such a result would also be in tension with the oscillation Δm² themselves, making it a deep consistency test. Conversely a quasi-degenerate detection (Σm_ν ≳ 100 meV) would contradict the BAM light scale.

ONE SPECTRUM, TWO OBSERVABLES. Σm_ν (~59–65 meV) and the 0νββ effective mass m_ββ (≲ 8 meV, PR #95) are both consequences of the SAME light, normal-ordered, Majorana spectrum — a joint, cross-checkable prediction.

HONEST SCOPE. ESTABLISHED (BAM-native): the light, normal-ordered spectrum gives Σm_ν ≈ 59–65 meV — at the NO floor, below the IO floor (~99 meV), consistent with Planck and at the DESI-DR2 + CMB frontier; a near-term-falsifiable target from the same spectrum as the PR #95 0νββ prediction. NOT established: the exact Σm_ν — a narrow band (59–65 meV) because the lightest mass is unmeasured (the BAM light scale keeps the band narrow); the absolute scale is the PR #90 residual.

## What this leaves open

- **The exact Σm_ν** — a narrow band (59–65 meV) because the lightest neutrino mass is unmeasured; the absolute scale is the PR #90 residual (the BAM light scale keeps the band narrow).
