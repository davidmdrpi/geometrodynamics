# Is 832 = N_q + ΔN an independent scale-ratio selection, or recycled n_part? (PR #107)

**Run:** 2026-05-31T02:16:49+00:00

Tests, skeptically, the tempting candidate derivation of PR #106's underived lepton/QCD scale ratio: `N_q + ΔN = 2N_q − N_lepton = 832 ≈ √σ/m_e ≈ 830` (0.2%). **Answer: we are recycling `n_part`.** 832 is built from the phenomenological compensator; it §8-drifts 764–920 (so the match is a baseline coincidence); and no independent bulk shell-stress integral selects ~466/832. The channel-normalisation derivation fails; `√σ/m_e` stays underived.

- **Answer**: recycling n_part — NOT an independent selection
- **Decisive test**: 4·n_part−100 drifts 764–920 (±9%) while 830 is fixed
- **Independent integrals**: O(10–70) (Σω²≈70, Σ(n+1)π≈47); never ~466/832
- **Ledger**: PR #106 unchanged — √σ/m_e underived; one scale G + one open ratio

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_the_832_observation` | N_q+ΔN = 2N_q−N_lepton = 832 ≈ 830 (0.2%) — tempting | **PASS** |
| T2 | `T2_832_built_from_npart` | 832 = 4·n_part − 4·k_5² — built from n_part | **PASS** |
| T3 | `T3_npart_is_phenomenological_compensator` | n_part is the compensator (§8-drifts 216–255, PR #76/#97) | **PASS** |
| T4 | `T4_s8_drift_decisive_test` | DECISIVE: 4·n_part−100 drifts 764–920 (±9%) vs fixed 830 | **PASS** |
| T5 | `T5_no_independent_integral_selects_832` | independent shell integrals O(10–70), never ~466/832 | **PASS** |
| T6 | `T6_circularity` | circular: n_part fit to the spectrum it would "predict" | **PASS** |
| T7 | `T7_ledger_unchanged` | ledger unchanged — √σ/m_e underived; one scale G + one ratio | **PASS** |
| T8 | `T8_assessment` | RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT | **PASS** |

## T4: the decisive §8-drift test

- observed `√σ/m_e = 830.3` (FIXED)
- `4·n_part − 100` across §8 ablations: `764`–`920` (span 156, ≈ ±9%)
- the baseline `n_part = 233` lands it at 832 (0.2% from 830) — but a quantity that drifts ±9% under ablations is **not** a stable selection of a fixed ratio. **Baseline coincidence.**

## T5: no independent integral selects 832

| bulk shell-stress integral | value |
|---|---:|
| Σ ω²(l=1, n=3..5) | 69.8 |
| Σ (n+1)π (Bohr–Sommerfeld closure) | 47.1 |

All `O(10–70)`. The number 466 enters **only** through the v3 closure count `4β_quark/(2π) = 2·n_part` — the fit. No independent integral yields ~466/832.

## Verdict

**RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT.** WE ARE RECYCLING n_part — 832 IS NOT AN INDEPENDENT SELECTION OF THE LEPTON/QCD SCALE RATIO. PR #106 left √σ/m_e ≈ 830 as the one underived ratio whose derivation would reduce BAM to a single anchor. The candidate N_q + ΔN = 832 ≈ 830 (0.2%) is seductive — and false.

832 IS BUILT FROM n_part. N_q + ΔN = 2·N_q − N_lepton = 2·(2·n_part) − 4·k_5² = 4·n_part − 4·k_5² = 4·233 − 100 = 832. So 832 is a linear function of n_part — the quark closure integer that PR #76/#97 established is a PHENOMENOLOGICAL COMPENSATOR: it absorbs the quark flavor puzzle, drifts 216–255 across the quark_axioms §8 ablations, and only its parity is invariant. 832 inherits that non-derived status.

THE DECISIVE §8-DRIFT TEST. If 832 independently selected the FIXED observed ratio 830, it would be §8-stable. It is not: propagating n_part ∈ {216..255} through 4·n_part − 100 gives [764, 920] — a span of 156, about ±9%. So the quantity drifts nearly ±10% while √σ/m_e = 830.3 is fixed; the baseline n_part = 233 merely happens to land it at 832 (0.2% from 830). That is a BASELINE COINCIDENCE — the same kind as 50π·k_5 = 785 (PR #106) and F_13 = 233 (PR #76) — not a stable selection.

NO INDEPENDENT BULK SHELL-STRESS INTEGRAL SELECTS 832. The genuine bulk shell-stress integrals over the Tangherlini geometry are O(10–70), nowhere near 466/832: the sum of shell eigenvalues Σω²(l=1,n=3..5) ≈ 69.8, the Bohr–Sommerfeld closure sum Σ(n+1)π ≈ 47.1. The number 466 enters ONLY through the v3 Hamiltonian closure count 4β_quark/(2π) = 2·n_part — i.e. through the fit. There is no independent integral that yields ~466 or ~832.

THE CIRCULARITY. n_part was FIT to reproduce the quark spectrum, which already encodes the physical scales. Recovering a scale ratio from n_part is therefore circular — you get back (an unstable version of) what was put in. So 832 ≈ 830 is the compensator echoing the spectrum it was fit to, not a derivation of the lepton/QCD hierarchy.

LEDGER UNCHANGED. The channel-normalisation derivation via N_q + ΔN fails. PR #106 stands: √σ/m_e ≈ 830 remains an UNDERIVED open dimensionless residual, and BAM remains at one foundational scale (G) + one open ratio. A genuine derivation would need an INDEPENDENT bulk shell-stress integral (not the v3-fit closure count) that selects ~830 AND is §8-stable — which is not available.

## What this leaves open

- **A genuine derivation of `√σ/m_e ≈ 830`** — would require an INDEPENDENT bulk shell-stress integral (not the v3-fit closure count) that selects ~830 AND is §8-stable. Not available; the PR #106 status (underived; one scale G + one open ratio) stands.
