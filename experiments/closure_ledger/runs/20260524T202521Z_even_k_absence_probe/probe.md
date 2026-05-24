# Even-k absence in the charged-lepton throat sector

**Run:** 2026-05-24T20:25:21+00:00

Classifies why even-k modes are not charged leptons: a spin-statistics selection rule. k mod 2 is the orientability grading of the throat closure (T^k off-diagonal/odd = orientation-reversing spin-½ fermion; diagonal/even = orientation-preserving boson); charged leptons (spin-½ Dirac fermions, #59–#66) are the odd class.

- **Grading**: k mod 2 = orientability (T^k off-diagonal odd / diagonal even)
- **Rule**: charged lepton = spin-½ fermion ⟺ orientation-reversing ⟺ odd k
- **Even k**: orientable double-cover / bosonic sector (not a charged lepton)
- **Not arithmetic**: Φ_avail(k) ≡ 0 mod 2π for every integer k
- **B4 caveat**: k dimensionless integer; selection topological, scale-independent

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_monodromy_Tk` | T^k off-diag (odd) / diag (even) = k mod 2 | **PASS** |
| T2 | `T2_orientation_classes` | odd = orientation-reversing fermion; even = boson | **PASS** |
| T3 | `T3_selection_rule` | spin-½ ⟹ odd k; (1,3,5) all odd | **PASS** |
| T4 | `T4_even_k_classified` | even k = orientable bosonic sector (not a lepton) | **PASS** |
| T5 | `T5_selection_not_arithmetic` | Φ_avail ≡ 0 mod 2π for any k (not arithmetic) | **PASS** |
| T6 | `T6_unification_with_spin_arc` | same T²=−I structure as #60/#61/#65/#66 | **PASS** |
| T7 | `T7_falsification_b4` | even-k bosons would falsify; BAM odd-k fermions | **PASS** |
| T8 | `T8_assessment` | spin-statistics selection rule | **PASS** |

## T1: Monodromy T^k

T = iσ_y, T² = −I: True

| k | parity | T^k off-diagonal | Z₂ class |
|---:|---|:---:|---|
| 1 | odd | True | opposite |
| 2 | even | False | same |
| 3 | odd | True | opposite |
| 4 | even | False | same |
| 5 | odd | True | opposite |
| 6 | even | False | same |

## T2: Orientation classes

- odd k = orientation-reversing (non-orientable) fermion: True
- even k = orientation-preserving (orientable) boson: True

## T3: The selection rule

- charged leptons are spin-½ fermions (#59–#66): True
- lepton depths [1, 3, 5] all odd: True; fermionic closure: True
- even k excluded from charged-lepton sector: True

## T4: Even-k classified

| k | T^k diagonal | sector | charged lepton? |
|---:|:---:|---|:---:|
| 2 | True | orientable double-cover (bosonic) | False |
| 4 | True | orientable double-cover (bosonic) | False |
| 6 | True | orientable double-cover (bosonic) | False |

## T5: Not arithmetic (the ledger closes for any k)

| k | Φ_avail/2π | closes mod 2π |
|---:|---:|:---:|
| 1 | 2.0 | True |
| 2 | 3.0 | True |
| 3 | 4.0 | True |
| 4 | 30.0 | True |
| 5 | 106.0 | True |
| 6 | 232.0 | True |

Even k also closes: True — the odd-only selection is spin-statistics, not arithmetic.

## T6: Unification with the spin / discrete-symmetry arc

Shared structure: T² = −I (the non-trivial RP³ spin structure, B2)

  - **B2_topological_sector**: T = iσ_y, T² = −I (the spin structure)
  - **PR60_wigner_rotation**: spinor double cover (−1 under 2π)
  - **PR61_g_factor**: g = 2 (spin-½ Pauli term)
  - **PR65_cpt**: T² = −I in the CPT operator
  - **PR66_throat_dirac_spinor**: the fermionic 4-spinor
  - **this_probe**: odd-k = the fermionic (non-orientable) closure

## T7: Falsification / B4

- leptons odd and fermionic: True
- even k would be bosonic (would falsify if leptons): True
- k dimensionless/topological: True

## T8: Assessment

- grading: k mod 2 = orientability (T^k off/on diagonal)
- rule: charged lepton = spin-½ fermion ⟺ orientation-reversing ⟺ odd k
- even k: orientable double-cover / bosonic sector (not a charged lepton)
- not arithmetic: Φ_avail(k) ≡ 0 mod 2π for every k
- unification: same T²=−I fermionic structure as #60/#61/#65/#66
- remaining: the even-k (bosonic) spectrum; why exactly 3 generations (k≤5)

## Verdict

**EVEN_K_EXCLUDED_BY_SPIN_STATISTICS.** EVEN-K EXCLUDED BY SPIN-STATISTICS. The even-k absence from the charged-lepton throat sector is classified as a spin-statistics selection rule, upgrading the odd-k closure lemma from a "choice of sector" to a genuine rule.

THE GRADING. Each throat pass applies T = iσ_y (T²=−I, B2); the spinor monodromy T^k after k passes is off-diagonal for odd k (mapping the spinor to the opposite Z₂ class — the antipodal/deck partner p~−p on RP³, the orientation-reversing closure across the non-orientable throat = the non-trivial spin structure = a spin-½ fermion) and diagonal for even k (same Z₂ class — orientation-preserving on the orientable double cover S³ = bosonic/integer class). So k mod 2 IS the orientability (= spin-statistics) grading of the closure.

THE RULE. The charged lepton is a Dirac spin-½ fermion (established across the recent arc: throat Dirac spinor #66, g=2 #61, Wigner rotation #60, CPT T²=−I #65), which requires the orientation-reversing (non-orientable) closure ⟹ odd k. The lepton depths (1,3,5) are all odd and fermionic; even k (orientation-preserving, bosonic) is excluded from the charged-lepton sector. The exclusion is NOT arithmetic — the Layer-1 ledger Φ_avail(k) = 2π(k+1) + 50π·max(0,k−3)² ≡ 0 mod 2π for every integer k (even k closes just as well) — it is spin-statistics / topology.

WHAT EVEN-K IS. Even-k modes are not forbidden; they are the orientation-preserving closures on the orientable double cover S³, the bosonic/integer-class sector of the same geometry — simply not charged leptons. The even-k absence is the spin-statistics face of the fermionic throat: the same T²=−I structure as #60/#61/#65/#66, seen from the generation-counting side. B4: k is a dimensionless integer; the selection is topological, scale-independent. Remaining: the even-k (bosonic) spectrum, and why the charged leptons stop at k=5 (three generations).

## What this leaves open

- **The even-k (bosonic) spectrum.** Whether the orientable even-k closures host a physical (integer-spin) spectrum is not worked out — only that they are not charged leptons.
- **Generation count (why exactly 3).** The odd-k rule allows k = 1,3,5,7,…; why the charged leptons stop at k=5 (three generations) is a separate question (the β-uplift / closure cutoff).
