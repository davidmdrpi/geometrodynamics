# The analytic Berger-Dirac ladder - companion probe (PR #197)

**Run:** 2026-07-02T15:43:59+00:00

The deliverable is `docs/berger_dirac_analytic_ladder.md` - the closed-form Dirac spectrum on the Berger family applied to the odd-k ladder. This probe verifies every identity to machine precision. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | item 2: the analytic upgrade; the document is the argument | **PASS** |
| T2 | `T2_closed_form_vs_operator` | closed forms = assembled operator (1e-15, j <= 3, four lambda) | **PASS** |
| T3 | `T3_checkpoints` | round +-(3/2+n) mult (n+1)(n+2) exact; collapse -> S2(1/2) | **PASS** |
| T4 | `T4_lichnerowicz_and_asymmetry` | Lichnerowicz: all zeros at lam > 2; eta asymmetry off round | **PASS** |
| T5 | `T5_ladder_protection` | ladder ordered+gapped; O(1) sensitivities (no fine-tuning) | **PASS** |
| T6 | `T6_boundaries_closed_form` | boundaries exact: sqrt(2k+4); 5.668/8.035/9.851 (e first) | **PASS** |
| T7 | `T7_k5_question_and_scope` | k5: NO spectral cutoff at any lambda (gap independent of k) | **PASS** |
| T8 | `T8_assessment` | spectral fact with located boundaries; cutoff dynamical | **PASS** |

## The boundaries (closed form vs numeric roots)

| k | lam_x closed | lam_x numeric | lam* closed | lam* numeric |
|---:|---:|---:|---:|---:|
| 1 | 2.44949 | 2.44949 | 5.667849 | 5.667849 |
| 3 | 3.162278 | 3.162278 | 8.034777 | 8.034777 |
| 5 | 3.741657 | 3.741657 | 9.850411 | 9.850411 |

## The round-point sensitivities (metric-soft)

| k | numeric | closed form (1/2-k)/(k+1/2) |
|---:|---:|---:|
| k=1 | -0.333333 | -0.333333 |
| k=3 | -0.714286 | -0.714286 |
| k=5 | -0.818182 | -0.818182 |

## Verdict

**ANALYTIC_BERGER_DIRAC_ODD_K_LADDER_PROTECTED_WITH_ALL_CROSSINGS_LOCATED_IN_CLOSED_FORM_NO_SPECTRAL_K5_CUTOFF_AT_ANY_LAMBDA.** CLOSED-FORM SPECTRAL GEOMETRY (the argument is in docs/berger_dirac_analytic_ladder.md; this probe checks its identities).

THE SPECTRUM. The Berger Dirac operator reduces by Peter-Weyl to exact 2x2 blocks: family A (the winding tower) a_j = (2j+1)/lam + lam/2 and family B b+- = lam/2 +- 2 sqrt((j+1/2)^2 + m'^2(lam^-2 - 1)) - validated against the assembled operator to 1e-15, against the round spectrum +-(3/2+n) with EXACT multiplicities (n+1)(n+2), against the lam->0 collapse limit (the S2(1/2) Dirac spectrum), and against Lichnerowicz (every zero lies in the scal < 0 regime lam > 2).

THE LADDER. The odd-k sector grounds are m_k(lam) = k/lam + lam/2 with UNIFORM gaps 2/lam: ordered and gapped at every lambda, absolutely protected on the squash side, smooth with O(1) sensitivities (-1/3, -5/7, -9/11) at the round point - the infinitesimal-squash refutation is analytically excluded, and the #192 surrogate's fine-tuning has no counterpart in the true spinor spectrum. Every crossing is located in closed form: character changes at lam_x(k) = sqrt(2k+4) (sqrt6 first, a 145% stretch margin around round), harmonic-spinor masslessness at lam*(k) = 5.668, 8.035, 9.851 - the ELECTRON SECTOR COLLAPSES FIRST under extreme stretch.

THE k5 QUESTION. No spectral counterpart at any lambda: the gap 2/lam is independent of k, the tower continues 7, 9, ... with identical spacing, and the collapse boundaries grow with k (stretch removes sectors from the bottom, never truncates the top). The three-generation cutoff is DYNAMICAL - the phase budget - now as a closed-form statement on the whole Berger family.
