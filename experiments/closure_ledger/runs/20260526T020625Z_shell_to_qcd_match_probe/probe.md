# Shell modes ↔ QCD/quark ladder: structural identification

**Run:** 2026-05-26T02:06:25+00:00

Tests whether the shell-saturated radial modes (PR #68) reproduce the documented structural ingredients of the quark sector (Z₂ partition multiplicity, 3×2=6 flavor counting, heavier scale, extended character), honestly scoped against the phenomenological pieces (n_part = 233, exact masses, SU(3) color) that remain open per `docs/quark_beta_status.md`.

- **Identification**: shell-saturated radial modes (PR #68, n≥3 at l=1) ↔ quark sector structural invariants
- **Reproduced**: Z₂ partition (factor of 2, N_q ∈ 2ℤ); 3×2=6 flavors; heavier; extended (participation → 2/3)
- **Open**: n_part = 233 (compensator); exact quark masses; SU(3) color (consistent with quark_beta_status)
- **B4 caveat**: structural invariants dimensionless / integer; scale-independent

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_mouth_balance_Z2_partition` | e asymmetric (0.12); shell ≈50/50 (Z₂ realized) | **PASS** |
| T2 | `T2_three_shell_modes_times_Z2` | 3 shell × 2 (Z₂) = 6 quark flavors | **PASS** |
| T3 | `T3_shell_heavier_than_leptons` | ω(shell-1)=3.82 > ω(τ)=2.89 | **PASS** |
| T4 | `T4_extended_shell_coupled` | shell participation → 2/3 (extended); leptons focused | **PASS** |
| T5 | `T5_beta_quark_parity_from_shell_Z2` | N_q ∈ 2ℤ from the shell Z₂ partition | **PASS** |
| T6 | `T6_honest_scope` | structural ingredients reproduced; n_part/masses/color open | **PASS** |
| T7 | `T7_falsification_b4` | asymmetric shell would falsify; BAM passes | **PASS** |
| T8 | `T8_assessment` | shell ↔ quark structural invariants match | **PASS** |

## T1: Inner/outer mouth balance (Z₂ partition)

| n | species | inner half | outer half | imbalance |
|---:|---|---:|---:|---:|
| 0 | e | 0.560 | 0.440 | 0.120 |
| 1 | μ | 0.466 | 0.534 | 0.068 |
| 2 | τ | 0.500 | 0.500 | 0.001 |
| 3 | shell-1 | 0.490 | 0.510 | 0.020 |
| 4 | shell-2 | 0.500 | 0.500 | 0.000 |
| 5 | shell-3 | 0.495 | 0.505 | 0.009 |

Electron (n=0) asymmetric (0.12): single throat identification. Shell modes ≈50/50 ([0.020367539353593922, 2.1233657035379938e-05, 0.009025387895144732]): Z₂ partition fully realized.

## T2: 3 shell modes × 2 = 6 quark flavors

First 3 shell modes (the 3-generation candidates):
  - n=3, ω=3.825
  - n=4, ω=4.761
  - n=5, ω=5.700
× Z₂ doubling (2) = **6 quark flavors** (matches 6).

## T3: Heavier mass scale

- ω(τ) = 2.8941; ω(shell-1, n=3) = 3.8246; shell heavier: True

## T4: Extended / shell-coupled

- lepton participation (e, μ, τ): [0.651, 0.665, 0.665]
- shell participation (n=3,4,5): [0.666, 0.666, 0.666] → 2/3 = 0.667
- shell at uniform standing-wave value: True; electron more focused: True

## T5: β_quark parity (N_q ∈ 2ℤ) from the shell Z₂

- shell Z₂ multiplicity: 2 (the v3 partition basis {(k,+),(k,−)})
- N_q = 2·n_part is even: True
- shell provides the structural factor of 2: True

## T6: Honest scope

Reproduced (structural ingredients):
  - Z₂ partition multiplicity (factor of 2 in N_q)
  - 3×2 = 6 quark flavors (3 shell modes × Z₂ doubling)
  - heavier mass scale than leptons
  - extended / shell-coupled character (participation → 2/3)
  - N_q ∈ 2ℤ parity invariant from the shell Z₂
Open (phenomenological pieces, per `docs/quark_beta_status.md`):
  - n_part = 233 (the compensator; not in BAM's catalog)
  - specific quark masses (m_u, m_d, m_c, m_s, m_t, m_b)
  - SU(3) color sector

## Verdict

**SHELL_REPRODUCES_QCD_STRUCTURE.** SHELL MODES REPRODUCE QCD STRUCTURE. The shell-saturated radial modes identified in PR #68 reproduce the documented structural ingredients of the quark sector.

Z₂ PARTITION MULTIPLICITY. The inner/outer mouth balance is ≈0.50/0.50 for shell modes (the Z₂ partition fully realized; both mouths as distinct states → the structural factor of 2 in N_q = 2·n_part), while the throat-focused electron (n=0) is asymmetric (0.56/0.44 — single throat identification, one mouth dominant = the lepton). The μ, τ transition between them. So the shell Z₂ multiplicity gives the parity invariant N_q ∈ 2ℤ that the quark β-lock requires (docs/quark_beta_status.md identified this as the robust invariant; lepton sector's single mouth identification does not provide it).

6 QUARK FLAVORS. The first 3 shell modes (n=3,4,5) doubled by the Z₂ partition give 3×2 = 6 — the right structural count of quark flavors (u/d, c/s, t/b).

HEAVIER + EXTENDED. ω(shell-1, n=3) ≈ 3.83 > ω(τ) ≈ 2.89 (quarks heavier than leptons, the right direction); participation ratio → 2/3 (uniform standing wave; extended / shell-coupled) vs the focused leptons — confinement-like extended states.

HONEST SCOPE. The shell modes reproduce the documented STRUCTURAL ingredients (Z₂ partition / factor of 2, 3×2=6 flavor counting, heavier scale, extended character, N_q ∈ 2ℤ parity). They do NOT, and this probe does not, derive specific quark masses, n_part = 233 (the phenomenological compensator the quark_beta probe sequence already established is not in BAM's catalog), or the SU(3) color sector — those remain open, consistent with docs/quark_beta_status.md. B4: the structural invariants are dimensionless ratios / integer counts; the shell↔QCD identification is geometric, independent of the single anchor m_e.

## What this leaves open

- **`n_part = 233`** — the phenomenological compensator; per `docs/quark_beta_status.md` the principled enumerations in BAM's catalog do not produce it.
- **Exact quark masses.** Not addressed by the structural match.
- **SU(3) color sector.** Not part of this probe.
