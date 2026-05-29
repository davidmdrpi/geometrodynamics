# CP / Majorana phase probe (PR #94)

**Run:** 2026-05-29T04:25:02+00:00

PRs #92–#93 fixed the PMNS mixing angles. This probe gives the BAM-native statements about the three CP-violating phases. **CP violation is generic** (the winding amplitudes carry the complex Hopf holonomy, so the PMNS is generically complex); the **Jarlskog dichotomy** mirrors the angle dichotomy (PMNS anarchic / large CP violation, CKM aligned / suppressed); and **two Majorana phases exist** because the neutrino is Majorana (`c₁=0`, PR #86).

- **Identification**: CP violation generic (Hopf-complex winding amplitudes; Jarlskog dichotomy PMNS-anarchic / CKM-aligned); two Majorana phases exist because the neutrino is Majorana (c₁=0); specific values anarchic, not pinned
- **CP source**: Hopf holonomy e^{ikχ} (PR #60) ⟹ complex PMNS
- **Jarlskog dichotomy**: PMNS anarchic (large J), CKM aligned (suppressed J)
- **Majorana phases**: two exist ⟸ neutrino Majorana ⟸ c₁=0 (PR #86); 0νββ
- **Open**: specific phase values (anarchic, not pinned); δ_CP poorly measured

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_phases_remain` | angles done (#92–#93); 3 phases (δ_CP + 2 Majorana) remain | **PASS** |
| T2 | `T2_cp_violation_generic_hopf_phase` | CP generic: winding amplitudes Hopf-complex; CP-cons. measure-zero | **PASS** |
| T3 | `T3_delta_cp_anarchic_uniform` | δ_CP anarchic (uniform); observed consistent | **PASS** |
| T4 | `T4_jarlskog_pmns_typical_of_anarchy` | |J_PMNS| ≈ 0.026 typical of anarchy (52nd/82nd pct) | **PASS** |
| T5 | `T5_jarlskog_dichotomy_ckm_aligned` | |J_CKM| ≈ 3e-5 extremely atypical (~0.1th pct) = aligned | **PASS** |
| T6 | `T6_two_majorana_phases_from_c1_zero` | two Majorana phases exist ⟸ c₁=0 (Majorana); 0νββ | **PASS** |
| T7 | `T7_honest_scope` | CP generic + Majorana existence firm; values anarchic | **PASS** |
| T8 | `T8_assessment` | CP_GENERIC_MAJORANA_PHASES_EXIST_VALUES_ANARCHIC | **PASS** |

## T4–T5: The Jarlskog dichotomy

`|J| ≤ 1/(6√3) ≈ 0.0962`; anarchic median `≈ 0.0255`.

| | observed \|J\| | anarchy percentile (μ=0) | residual (μ≈3) | reading |
|---|---:|---:|---:|---|
| **PMNS** | 0.0260 | 51th | 81th | typical → large CP violation |
| **CKM** | 3.08e-05 | 0.08th | — | extremely atypical → aligned/suppressed |

The Jarlskog mirrors the angle dichotomy (PR #92): PMNS anarchic (cross-coordinate, generic CP violation), CKM aligned (intra-coordinate, suppressed). CP conservation (`J=0`) is measure-zero — CP is generically violated.

## T6: Two Majorana phases exist because c₁=0

- Dirac neutrino: 1 physical phase (δ_CP), no Majorana phases.
- Majorana neutrino: 3 physical phases ⟹ **2 Majorana phases**.
- In BAM the neutrino is Majorana ⟸ c₁=0 (C-invariant, PR #86), so the two Majorana phases EXIST — observable in 0νββ (neutrinoless double beta decay).

## Verdict

**CP_GENERIC_MAJORANA_PHASES_EXIST_VALUES_ANARCHIC.** CP VIOLATION IS GENERIC, AND TWO MAJORANA PHASES EXIST BECAUSE THE NEUTRINO IS MAJORANA. PRs #92–#93 fixed the PMNS mixing angles (anarchic, with θ13 suppressed by a residual nearest-neighbour alignment). This probe addresses the remaining sector — the three CP-violating phases.

CP VIOLATION IS GENERIC (THE HOPF PHASE). CP conservation requires a real PMNS matrix (up to rephasing) — a measure-zero condition. In BAM the winding (charged-lepton) amplitudes carry the Hopf holonomy e^{ikχ} (PR #60: the throat Berry phase ∮A = π cos χ), so they are intrinsically COMPLEX; the cross-channel overlaps that build the PMNS matrix (PR #92) are therefore generically complex, and δ_CP ≠ 0, π with probability 1. No BAM symmetry forces real amplitudes, so CP is generically violated.

THE JARLSKOG DICHOTOMY. The rephasing-invariant measure of Dirac CP violation is J = Im(U_e1 U_μ2 U*_e2 U*_μ1), with |J| ≤ 1/(6√3) ≈ 0.096 and anarchic median ≈ 0.025. The observed |J_PMNS| ≈ 0.026 (δ ~ 230°) is the ~52nd percentile of anarchy (μ=0), ~82nd with the PR #93 residual alignment — typical of anarchy: large, generic CP violation. The observed |J_CKM| ≈ 3.1×10⁻⁵ is the ~0.1th percentile — extremely atypical = aligned ⟹ CP violation suppressed. So the Jarlskog mirrors the angles (PR #92): PMNS anarchic (cross-coordinate, large CP violation), CKM aligned (intra-coordinate, suppressed).

TWO MAJORANA PHASES EXIST BECAUSE c₁=0. A Dirac neutrino has one physical phase (δ_CP) and no Majorana phases; a Majorana neutrino has two extra. In BAM the neutrino is Majorana precisely because it is chargeless (c₁=0, C-invariant, PR #86), so the EXISTENCE of two physical Majorana phases is a firm BAM prediction — CP-violating phases of the ΔL=2 throat↔antithroat sector (the bounce of PRs #87–#90), observable in 0νββ. Were the neutrino Dirac, there would be none.

HONEST SCOPE. ESTABLISHED (BAM-native): CP violation is generic (the winding amplitudes are Hopf-complex; CP conservation is measure-zero); the Jarlskog dichotomy mirrors the angle dichotomy (PMNS typical of anarchy, CKM extremely atypical = aligned/suppressed); and two physical Majorana phases EXIST because the neutrino is Majorana (c₁=0, PR #86) — a firm prediction for 0νββ. NOT established: the specific values of δ_CP and the two Majorana phases — they are anarchic (uniform), set by the Hopf phases of the cross-channel overlaps and the throat↔antithroat tunnelling, and not pinned (δ_CP is itself poorly measured, consistent with the uniform anarchic expectation).

## What this leaves open

- **The specific phase values** (δ_CP and the two Majorana phases) — anarchic (uniform), set by the Hopf phases of the cross-channel overlaps and the throat↔antithroat tunnelling; not pinned.
- **δ_CP measurement** — itself poorly measured; the observed value is consistent with the uniform anarchic expectation.
