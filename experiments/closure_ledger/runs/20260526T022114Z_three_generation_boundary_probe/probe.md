# Three-generation boundary: β-uplift cutoff + throat-shell availability

**Run:** 2026-05-26T02:21:14+00:00

Pins the sharp k ≤ 5 charged-lepton boundary as the join of the β-uplift quadratic closure-quantum growth (within-lepton scaling) and the throat-shell mode availability cutoff (#68), where the B2 mapping sends the would-be 4th generation to the shell/QCD sector (#69).

- **Boundary**: k ≤ 5 (charged-lepton sector); k=7 (and beyond) → shell/QCD sector
- **Within-lepton scaling**: β-uplift = 50π·(k−3)² → 0,0,100,400,900 (quadratic)
- **Sector boundary**: throat ↔ shell at n=3 (#68); B2: k=7 → n=3 (shell, QCD #69)
- **Closure arithmetic**: closes for any k (recap #67) — cutoff is NOT arithmetic
- **B4 caveat**: k dimensionless integer; boundary topological/structural

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_beta_uplift_quadratic_growth` | β-uplift = 50π·(k−3)²: 0, 0, 100, 400, 900 (quadratic) | **PASS** |
| T2 | `T2_closure_arithmetic_holds_for_any_k` | Φ_avail ≡ 0 mod 2π for any k (cutoff not arithmetic) | **PASS** |
| T3 | `T3_k_to_n_mapping_B2` | B2: n=(k−1)/2; (k,n)=(1,0),(3,1),(5,2),(7,3),… | **PASS** |
| T4 | `T4_throat_shell_availability` | n=0,1,2 throat (lepton); n≥3 shell (QCD, #69) | **PASS** |
| T5 | `T5_would_be_4th_generation_in_shell` | k=7 → n=3 (shell/QCD); no 4th charged lepton | **PASS** |
| T6 | `T6_combined_mechanism_pins_k_le_5` | β-uplift + throat-shell → k ≤ 5 (3 generations) | **PASS** |
| T7 | `T7_falsification_b4` | 3 generations observed; 4th would falsify | **PASS** |
| T8 | `T8_assessment` | three throat-localized odd-k fermionic modes | **PASS** |

## T1: β-uplift quadratic growth

| k | β-uplift / 2π | activates? |
|---:|---:|:---:|
| 1 | 0 | False |
| 3 | 0 | False |
| 5 | 100 | True |
| 7 | 400 | True |
| 9 | 900 | True |
| 11 | 1600 | True |

## T2: Closure arithmetic closes for any k (recap #67)

| k | N_total | closes mod 2π |
|---:|---:|:---:|
| 1 | 2 | True |
| 3 | 4 | True |
| 5 | 106 | True |
| 7 | 408 | True |
| 9 | 910 | True |
| 11 | 1612 | True |

## T3: k ↔ n mapping (B2)

| k | n = (k−1)/2 | species |
|---:|---:|---|
| 1 | 0 | e |
| 3 | 1 | μ |
| 5 | 2 | τ |
| 7 | 3 | (would-be n=3) |
| 9 | 4 | (would-be n=4) |

## T4: Throat-shell availability cutoff at n=3 (recap #68 + #69)

| n | sector |
|---:|---|
| 0 | throat-localized (lepton, #68) |
| 1 | throat-localized (lepton, #68) |
| 2 | throat-localized (lepton, #68) |
| 3 | shell-saturated (QCD, #69) |
| 4 | shell-saturated (QCD, #69) |
| 5 | shell-saturated (QCD, #69) |

## T6: Combined mechanism (the sharp k ≤ 5 boundary)

| k | n | β-uplift / 2π | throat-localized? | charged lepton? |
|---:|---:|---:|:---:|:---:|
| 1 | 0 | 0 | True | True |
| 3 | 1 | 0 | True | True |
| 5 | 2 | 100 | True | True |
| 7 | 3 | 400 | False | False |

Boundary at k ≤ 5: True

## Verdict

**THREE_GENERATIONS_PINNED.** THREE GENERATIONS PINNED. The sharp k ≤ 5 charged-lepton boundary is the join of the β-uplift quadratic closure-quantum growth and the throat-shell mode availability cutoff (#68), where the B2 mapping places the would-be 4th generation in the shell/QCD sector (#69).

β-UPLIFT (the within-lepton scaling). The locked lepton surrogate uses Φ_avail(k) = 2π(k+1) + 50π·max(0,k−3)² (the odd-k lemma); the second term gives 25·(k−3)² closure quanta — 0, 0, 100, 400, 900, … — activating at k=5 (τ) and exploding quadratically beyond. This is the within-lepton mass-scaling structure that fits the muon and tau masses correctly; it quantifies the steep closure-quantum cost of any would-be higher generation but is NOT itself the sector cutoff (closure arithmetic alone closes for every integer k, #67).

THROAT-SHELL AVAILABILITY (the hard cutoff). The throat-localized fermionic ladder has exactly three modes (n=0,1,2 = e,μ,τ, #68); from n=3 the modes are shell-saturated (participation → 2/3, the QCD/quark sector with the structural Z₂ partition, #69). The B2 reading n=(k−1)/2 maps the would-be 4th generation k=7 to n=3 — which is in the shell/QCD sector, NOT a throat-localized charged lepton. So there is no 4th charged lepton.

JOIN. The three observed generations are exactly the three throat-localized odd-k fermionic modes (e,μ,τ at k=1,3,5 = n=0,1,2); higher k modes are in the shell/QCD sector. The β-uplift quantifies the within-lepton scaling; the throat-shell availability provides the hard cutoff. Together they pin k ≤ 5 sharply. Closure arithmetic alone does not cut (#67); the cutoff is the mode-availability boundary.

HONEST SCOPE. The β_lepton = 50π itself is the locked surrogate's structure (closure-quantum origins in docs/hbar_origin_status.md: 8π=4·(2π) transport, 7π/100 resistance); a first-principles derivation of the specific 50π factor and why exactly three throat-localized modes fit in the cavity is future work. B4: k is a dimensionless integer; the boundary is topological/structural — scale-independent.

## What this leaves open

- **First-principles β_lepton = 50π.** The closure-quantum origins (`8π = 4·(2π)` transport, `7π/100` resistance) are established in `docs/hbar_origin_status.md`; a clean derivation of the specific 50π factor sharpens further.
- **Why exactly 3 throat-localized modes.** The cavity geometry (R_OUTER − R_MID, the centrifugal barrier) sets it; a deeper structural argument is future work.
