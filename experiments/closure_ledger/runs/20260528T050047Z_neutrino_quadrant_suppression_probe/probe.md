# Neutrino-quadrant suppression: Majorana seesaw from `c₁ = 0` (PR #86)

**Run:** 2026-05-28T05:00:47+00:00

PR #85 found the `(k=0, n<3)` neutrino quadrant gives the lightest states but ~10⁵–10⁶ too heavy for neutrinos. This probe identifies the BAM-native suppression: the neutrino is **Majorana** (because `k=0 ⟹ c₁=0 ⟹ C-invariant`), and the **Majorana seesaw** `m_ν = m_D²/M_R` is the mechanism.

- **Identification**: neutrino suppression = Majorana seesaw, forced by k=0 ⟹ c₁=0 ⟹ C-invariance; available only to the chargeless sector (charged leptons are Dirac); seesaw scale M_R ~ TeV is a new open heavy input
- **Mechanism**: Majorana seesaw m_ν = m_D²/M_R
- **Why only neutrinos**: only c₁=0 (k=0) sector is C-invariant ⟹ Majorana
- **Open**: M_R (lepton-number-violating / throat↔antithroat scale)
- **B4 caveat**: m_D from cavity floor + electron anchor; M_R open; absolute m_ν needs M_R + L_eff unification + B4 anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_suppression_gap` | bare ν/charged ≈ 0.07; need ~10⁵–10⁶ suppression | **PASS** |
| T2 | `T2_neutrino_is_majorana_from_c1_zero` | k=0 ⟹ c₁=0 ⟹ C-invariant ⟹ Majorana (BAM-native) | **PASS** |
| T3 | `T3_charged_lepton_is_dirac_no_seesaw` | charged lepton k≠0 ⟹ Dirac ⟹ no seesaw (why only ν) | **PASS** |
| T4 | `T4_bare_dirac_masses_cavity_floors` | bare Dirac masses = cavity floors ~43,80,118 keV | **PASS** |
| T5 | `T5_required_seesaw_scale` | required seesaw M_R = m_D²/m_ν ≈ 0.3–1.8 TeV | **PASS** |
| T6 | `T6_bam_scale_candidates_for_M_R` | no current BAM scale at ~TeV → M_R is open input | **PASS** |
| T7 | `T7_honest_scope` | mechanism BAM-native; scale M_R open | **PASS** |
| T8 | `T8_assessment` | NEUTRINO_SUPPRESSION_IS_MAJORANA_SEESAW_SCALE_OPEN | **PASS** |

## T2–T3: Why only neutrinos are suppressed

| sector | k | c₁ | under C (c₁→−c₁) | nature | seesaw? |
|---|---:|---:|---|---|:---:|
| neutrino | 0 | 0 | 0 → 0 (invariant) | **Majorana** | ✓ |
| charged lepton | 2g−1 | ±1 | +1 → −1 (e⁻→e⁺) | Dirac | ✗ |

Only the chargeless (`c₁=0`) sector is C-invariant, hence Majorana, hence gets the ΔL=2 seesaw. Charged leptons are Dirac and keep their full winding mass `β·k²`.

## T4–T5: Bare Dirac masses and required seesaw scale

| gen | n | m_D (keV) | m_ν obs (eV) | M_R (GeV) |
|---:|---:|---:|---:|---:|
| 1 | 0 | 42.9 | 0.001 | 1836 |
| 2 | 1 | 80.2 | 0.009 | 715 |
| 3 | 2 | 117.6 | 0.05 | 277 |

Required seesaw scale `M_R ≈ 1836–277 GeV` (order TeV), a single new heavy scale.

## T6: BAM heavy-scale candidates for `M_R`

| candidate scale | value (eV) | matches ~TeV? |
|---|---:|:---:|
| pair-production 2 m_e c² | 1.02e+06 | ✗ |
| leptoquark gen-1 (raw operator) | 5.32e+05 | ✗ |
| bulk/cosmological ℏc/R_cosmo (Hubble) | 1.00e-33 | ✗ |
| required seesaw M_R | 1.00e+12 | ✓ |

No existing BAM scale lands at ~TeV. `M_R` — the lepton-number-violating (throat↔antithroat) scale — is a new heavy input, not yet BAM-derivable.

## Verdict

**NEUTRINO_SUPPRESSION_IS_MAJORANA_SEESAW_SCALE_OPEN.** NEUTRINO SUPPRESSION IS THE MAJORANA SEESAW; SCALE OPEN. PR #85 found the (k=0, n<3) neutrino quadrant gives the lightest states but with a mass ratio m_ν/m_charged ≈ 0.07 — ~10⁵–10⁶ too heavy for neutrinos. This probe identifies the BAM-native suppression mechanism.

THE STRUCTURAL KEY: k = 0 ⟹ c₁ = 0 ⟹ MAJORANA. The throat winding number k IS the Hopf charge: a winding mode carries c₁ = ±1, a non-winding mode carries c₁ = 0. Under charge conjugation C (the inner/outer swap c₁ → −c₁, PR #63), the neutrino quadrant (k=0, c₁=0) maps 0 → −0 = 0 — invariant. The neutrino is its own C-conjugate, hence NECESSARILY MAJORANA. This is not assumed; it follows from the chargeless (k=0) sector being the C-invariant sector.

THE SUPPRESSION: MAJORANA SEESAW. A Majorana mass violates lepton number by ΔL = 2 and admits the seesaw m_ν = m_D²/M_R, where m_D is the bare neutrino-quadrant Dirac mass (the cavity floor) and M_R is the heavy lepton-number-violating scale — in BAM, the throat ↔ antithroat coupling scale (a Majorana mass flips throat to antithroat, ΔL = 2). Because M_R ≫ m_D, the seesaw produces m_ν ≪ m_D automatically: the anomalous lightness of neutrinos is GENERIC to the mechanism.

WHY ONLY NEUTRINOS ARE SUPPRESSED. The seesaw is available ONLY to the Majorana (C-invariant, c₁=0) sector. Charged leptons (k≠0, c₁=±1) are Dirac — under C, e⁻ → e⁺ (distinct), so they cannot have a ΔL=2 Majorana mass, get no seesaw, and keep their full winding mass β·k². This is the BAM-native explanation of why the neutrino is anomalously light and the charged lepton is not — the two differ precisely by their Hopf charge (winding) c₁.

NUMBERS. With the electron anchored at (k=1, n=0), the bare neutrino Dirac masses are the cavity floors m_D ≈ 43, 80, 118 keV (generations 1, 2, 3 at n = 0, 1, 2). For observed neutrino masses ~0.001–0.05 eV, the required seesaw scale is M_R = m_D²/m_ν ≈ 0.3–1.8 TeV — a single new heavy scale, roughly generation-uniform.

THE SCALE M_R IS OPEN. No BAM scale in the current catalog lands at ~TeV: pair-production ~1 MeV, leptoquark sub-MeV (raw operator units), bulk/cosmological ~10⁻³³ eV. So the seesaw scale M_R is a NEW heavy input — the lepton-number-violating (throat↔antithroat / B−L breaking) scale — not yet BAM-derivable. The TeV value is read off from the observed neutrino mass, not predicted.

HONEST SCOPE. ESTABLISHED (BAM-native): the neutrino is Majorana because the chargeless k=0 quadrant is the C-invariant sector; the Majorana seesaw is the natural suppression mechanism; and the seesaw is available only to the c₁=0 sector, explaining why neutrinos (not charged leptons) are anomalously light. NOT established: the value of M_R (a new heavy input, ~TeV from observed m_ν, not predicted; no current BAM scale matches), and absolute neutrino masses (need M_R + the L_eff unification still open from PR #83 + the B4 anchor).

## What this leaves open

- **The seesaw scale `M_R`** — the lepton-number-violating (throat↔antithroat / B−L) scale; ~TeV from observed `m_ν`, not predicted; no current BAM scale matches.
- **Absolute neutrino masses** — need `M_R` + the L_eff unification (PR #83 open) + the B4 anchor.
- **The neutrino mass ordering / mixing (PMNS)** — beyond this probe's scope.
