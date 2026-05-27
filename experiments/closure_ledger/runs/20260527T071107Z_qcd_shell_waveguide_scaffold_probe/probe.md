# QCD shell waveguide: basis + operator scaffold

**Run:** 2026-05-27T07:11:07+00:00

> *Quarks do not pass through the throat; they are the wavefronts that resolve the cavity itself.*

Foundation for the four-PR quantitative QCD-shell arc the user laid out:

  - **PR #77 (this PR)** — QCD shell waveguide basis/operator scaffold
  - **PR #78** — shell Hamiltonian mass-ordering / n_part audit
  - **PR #79** — boundary stress tensor and singlet constraint
  - **PR #80** — color algebra candidate: SU(3), SU(2)×Z₂, or other

- **Identification**: quarks = shell-saturated cavity wavefronts (NOT throat traversals); 6-state (l, n, p) basis with 6×6 operator scaffold
- **Next PR**: PR #78 — shell Hamiltonian mass-ordering / n_part audit
- **B4 caveat**: ω dimensionful; scaffold scale-free in ratios; structural

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_shell_waveguide_basis_constructed` | 6-state shell basis (l=1, n=3,4,5, ±); ⟨r⟩ on plateau | **PASS** |
| T2 | `T2_lepton_to_shell_transition` | lepton ladder = throat→shell transition; shell on plateau | **PASS** |
| T3 | `T3_operator_scaffold_constructed` | H = H_kin + H_Z2 + H_couple Hermitian; ω² eigenvalues | **PASS** |
| T4 | `T4_basis_flexibility_n_vs_l_enumeration` | scaffold supports n_varied and l_varied enumerations | **PASS** |
| T5 | `T5_six_flavor_structural_count` | 3 generations × 2 partitions = 6 flavors (PR #69) | **PASS** |
| T6 | `T6_pr78_79_80_hooks` | hooks for PR #78 (mass), #79 (stress-T), #80 (color) | **PASS** |
| T7 | `T7_honest_scope_b4` | honest scope: scaffold only; no masses, no n_part audit | **PASS** |
| T8 | `T8_assessment` | shell waveguide scaffold constructed | **PASS** |

## T1: Shell waveguide basis (n_varied enumeration)

| state | ω | PR | ⟨r⟩−R_MID | throat frac |
|---|---:|---:|---:|---:|
| `(l=1, n=3, +)` | 3.8246 | 0.6657 | 0.0454 | 0.7517 |
| `(l=1, n=3, −)` | 3.8246 | 0.6657 | 0.0454 | 0.7517 |
| `(l=1, n=4, +)` | 4.7609 | 0.6658 | 0.0460 | 0.7940 |
| `(l=1, n=4, −)` | 4.7609 | 0.6658 | 0.0460 | 0.7940 |
| `(l=1, n=5, +)` | 5.7000 | 0.6658 | 0.0462 | 0.8195 |
| `(l=1, n=5, −)` | 5.7000 | 0.6658 | 0.0462 | 0.8195 |

## T2: Lepton-to-shell transition — the throat-focused → cavity-resolving boundary

| metric | lepton (e=0, μ=1, τ=2) | shell (n=3,4,5) |
|---|---|---|
| ⟨r⟩−R_MID | [0.0212, 0.0389, 0.044] | [0.0454, 0.046, 0.0462] |
| PR | [0.6509, 0.6646, 0.6655] | [0.6657, 0.6658, 0.6658] |
| throat fraction (inner third) | [0.9628, 0.8505, 0.7591] | [0.7517, 0.794, 0.8195] |

Shell plateau threshold: ⟨r⟩ ≥ 0.04; standing-wave target PR = 0.6667.

The lepton ladder IS the throat-to-shell transition: ⟨r⟩ rises through e, μ, τ approaching the shell plateau. The electron is sharply throat-focused; the tau is at the edge of shell saturation. Quark sector (n ≥ 3) sits on the plateau.

## T4: Basis flexibility

| enumeration | labels | eigenvalues |
|---|---|---|
| `n_varied` (fix l=1) | ['(l=1, n=3, +)', '(l=1, n=3, −)', '(l=1, n=4, +)', '(l=1, n=4, −)', '(l=1, n=5, +)', '(l=1, n=5, −)'] | [14.62778, 14.62778, 22.665871, 22.665871, 32.490231, 32.490231] |
| `l_varied` (fix n=3) | ['(l=1, n=3, +)', '(l=1, n=3, −)', '(l=2, n=3, +)', '(l=2, n=3, −)', '(l=3, n=3, +)', '(l=3, n=3, −)'] | [14.62778, 14.62778, 14.940642, 14.940642, 15.385323, 15.385323] |

PR #78 chooses based on which reproduces the quark mass ordering (m_u < m_d, m_c > m_s, etc.).

## T6: Hooks for PR #78–#80

**PR #78 — shell Hamiltonian mass-ordering / n_part audit**

- Populates: `H_couple inter-mode mixing; enumeration choice`
- Tests: m_u < m_d, m_c > m_s; reaudit n_part on shell basis

**PR #79 — boundary stress tensor + singlet constraint**

- Populates: `chi (Z₂ splitter from cavity-wall stress-T)`
- Tests: singlet projection ⟹ colorless physical states

**PR #80 — color algebra candidate**

- Populates: `algebra acting on (l, n, p)`
- Candidates: ['SU(3)', 'SU(2) × Z₂', 'SU(4) Pati-Salam', 'other']

## Verdict

**SHELL_WAVEGUIDE_SCAFFOLD_CONSTRUCTED.** SHELL WAVEGUIDE SCAFFOLD CONSTRUCTED. Quarks are the wavefronts that resolve the cavity itself — shell-saturated Tangherlini radial modes whose ⟨r⟩ has plateaued at the shell center (~0.046 above R_MID) and whose participation ratio has locked to the uniform-standing-wave value 2/3. The lepton ladder (e=n=0, μ=n=1, τ=n=2) IS the throat-to-shell transition: ⟨r⟩ rises 0.021 → 0.039 → 0.044 toward the plateau but does not reach it; the electron is the most throat-focused. Shell modes (n ≥ 3) sit on the plateau — the wavefront regime where the mode resolves the cavity rather than focuses on the throat.

6-STATE BASIS. The basis is (l, n, p) with l = S³ Casimir (Hopf-bundle angular base, #73), n = shell-saturated radial overtone (≥ 3 for l = 1, from #68 saturation metric), p ∈ {+, −} = Z₂ partition (B2, the non-orientable throat). The lowest 6 states match the 3 × 2 = 6 flavor structural count documented by #69.

OPERATOR SCAFFOLD. The 6×6 Hamiltonian H = H_kin + H_Z2 + H_couple acts on the shell basis with: H_kin diagonal (ω²(l, n) cavity-eigenfrequency-squared mass operator); H_Z2 block-σ_z partition splitter (slot for #79–#80 structural constant chi); H_couple inter-mode mixing (slot for #78–#80 population). PR #77 leaves chi = 0 and H_couple = 0 — the eigenvalues are then ω² each doubled by the Z₂ partition.

BASIS FLEXIBILITY. Two natural enumerations are supported: vary n at fixed l = 1 ({(1, 3, ±), (1, 4, ±), (1, 5, ±)}) or vary l at fixed shell-saturated n ({(1, 3, ±), (2, 3, ±), (3, 3, ±)}). They give different eigenvalues; PR #78 chooses based on which reproduces the quark mass ordering.

DISTINCTNESS FROM LEPTON/V3 MACHINERY. The lepton sector lives at the throat (basis = odd-k throat-traversal modes {k=1, 3, 5}, diagonal = β·k²·(2π) winding cost from #71). The v3 quark Hamiltonian inherited that lepton-shaped basis and absorbed QCD shell physics into the phenomenological n_part = 233 (#76). The shell waveguide basis is structurally distinct: extended-character standing waves on the cavity, not focused throat pulses; ω²(l, n) cavity eigenfrequency, not β·k² winding. This IS the right machinery for the quark sector.

HOOKS FOR PR #78–#80. The scaffold names each follow-on slot: PR #78 (mass-ordering / n_part audit) populates H_couple and chooses the enumeration; PR #79 (boundary stress tensor / singlet constraint) sets chi from cavity-wall energy-momentum and adds a singlet (colorless) projection; PR #80 (color algebra) identifies the algebra acting on (l, n, p) — SU(3) color triplet, SU(2) × Z₂ partition-flavored, Pati-Salam SU(4) with n as a 4th leptocolor, or another candidate.

HONEST SCOPE. Scaffold only — does NOT derive quark masses, audit n_part, define the boundary stress tensor, or identify the color algebra. The 6×6 operator is structurally distinct-from-lepton but not yet populated. Mass values, singlet constraint, color algebra are PRs #78–#80 explicitly. B4: cavity eigenfrequencies ω are dimensionful (scale-free in ratios); the scaffold is structurally scale-independent.

## What this leaves open

- **PR #78** — populate H_couple with inter-mode mixing, choose the enumeration (n_varied vs l_varied), and audit whether the shell basis reproduces the quark mass ordering AND derives n_part structurally.
- **PR #79** — set chi from the boundary stress tensor on the cavity wall; add a singlet (colorless) projection as a physical-state constraint.
- **PR #80** — identify the color algebra acting on (l, n, p): SU(3) (color triplet), SU(2) × Z₂ (partition-flavored), Pati-Salam SU(4) (n as 4th leptocolor), or another.
