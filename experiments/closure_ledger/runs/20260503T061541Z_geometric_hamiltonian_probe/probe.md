# Geometric-Hamiltonian probe — Tangherlini matrix elements as the lepton H

**Run:** 2026-05-03T06:15:41+00:00
**Grid:** N = 80

Tests whether a 3×3 lepton Hamiltonian whose entries come directly from Tangherlini radial matrix elements (no fitted parameters) can reproduce the observed lepton mass ratios AND, via its eigenvectors, close the closure-phase ledger.

## Anchors

**Observed lepton masses (PDG)**

| species | mass (MeV) | ratio to electron |
|---|---:|---:|
| electron | 0.510999 | 1.0000 |
| muon | 105.658374 | 206.7683 |
| tau | 1776.860000 | 3477.2283 |

**Tangherlini ground-mode eigenfrequencies (l=1, 3, 5; n=0)**

| l | ω(l, 0) | ω² | ratio ω/ω(1) |
|---|---:|---:|---:|
| 1 | 1.054727 | 1.112449 | 1.000000 |
| 3 | 1.219083 | 1.486163 | 1.155828 |
| 5 | 1.395973 | 1.948742 | 1.323540 |

Tangherlini ω(5)/ω(1) = 1.3235; observed m_τ/m_e = 3477.23. Dynamic range mismatch is ~2627×.

**Locked-surrogate reference (`_build_generation_block`)**

| eigenvalue index | eigenvalue | √eigenvalue | ratio to lowest |
|---|---:|---:|---:|
| 0 | 0.1996 | 0.4468 | 1.0000 |
| 1 | 41.2600 | 6.4234 | 14.3763 |
| 2 | 694.9832 | 26.3625 | 59.0025 |

The locked surrogate's mass-ratio dynamic range comes predominantly from the diagonal `β · max(0, k − 3)²` term (β = 50π), which contributes 200π ≈ 628 to the τ row only. No Tangherlini matrix element is anywhere near that scale.

## Per-variant results

| variant | eigenvalues | √eigenvalue ratios (μ/e, τ/e) | log₁₀ mass-ratio error (μ, τ) | mass match? | B1 spread (rad) | D1 spread (rad) |
|---|---|---|---|---|---:|---:|
| `GH_A_diagonal_omega_squared` | [1.1124, 1.4862, 1.9487] | (1.156, 1.324) | (-2.25, -3.42) | no | 0.3814 | 0.0000 |
| `GH_B_potential_average` | [0.0000, 0.0284, 1.6659] | (24.007, 183.942) | (-0.94, -1.28) | no | 0.0898 | 1.8194 |
| `GH_C_induced_hamiltonian` | [0.6137, 1.2176, 2.7161] | (1.409, 2.104) | (-2.17, -3.22) | no | 0.2559 | 1.7474 |
| `GH_D_omega_squared_overlap` | [-1.8292, -1.3197, 7.6963] | (0.000, 0.000) | (+inf, +inf) | no | 0.1882 | 1.7561 |
| `GH_E_symmetric_momentum` | [0.0011, 0.0881, 2.8901] | (8.801, 50.396) | (-1.37, -1.84) | no | 0.1331 | 1.8897 |

### Per-variant detail

#### `GH_A_diagonal_omega_squared`

H_ij = δ_ij · ω_i². Eigenvalues = ω_i², eigenvectors = standard basis. Sanity baseline.

- Eigenvalues: ['1.1124', '1.4862', '1.9487']
- √eigenvalue ratios (predicted m_e:m_μ:m_τ = 1.0000 : 1.1558 : 1.3235); observed: 1 : 206.77 : 3477.23
- B1-with-probe-eigenvectors closure: residues ['0.8819', '0.7707', '0.7605'] π, circular spread 0.381385 rad
- D1-with-probe-eigenvectors closure: residues ['0.0000', '0.0000', '0.0000'] π, circular spread 0.000000 rad  ⚠ **trivial closure** (identity eigenvectors × D1's zero diagonal → Φ_radial = 0 for every species)

#### `GH_B_potential_average`

H_ij = ⟨u_i | V̄ | u_j⟩ with V̄ = (V_1+V_3+V_5)/3. Common potential operator projected onto the basis.

- Eigenvalues: ['0.0000', '0.0284', '1.6659']
- √eigenvalue ratios (predicted m_e:m_μ:m_τ = 1.0000 : 24.0070 : 183.9422); observed: 1 : 206.77 : 3477.23
- B1-with-probe-eigenvectors closure: residues ['0.7875', '0.8094', '0.8161'] π, circular spread 0.089811 rad
- D1-with-probe-eigenvectors closure: residues ['1.9105', '1.7552', '0.3343'] π, circular spread 1.819366 rad

#### `GH_C_induced_hamiltonian`

H_ii = ω_i²; H_ij = (ω_j²−ω_i²)·⟨u_i|u_j⟩ mirrored. The exact projection of H_l_j − H_l_i onto the {u_i} basis (provably equivalent to the D1 transport matrix element).

- Eigenvalues: ['0.6137', '1.2176', '2.7161']
- √eigenvalue ratios (predicted m_e:m_μ:m_τ = 1.0000 : 1.4086 : 2.1038); observed: 1 : 206.77 : 3477.23
- B1-with-probe-eigenvectors closure: residues ['0.8519', '0.7705', '0.7907'] π, circular spread 0.255934 rad
- D1-with-probe-eigenvectors closure: residues ['1.7772', '1.8893', '0.3335'] π, circular spread 1.747417 rad

#### `GH_D_omega_squared_overlap`

H_ii = ω_i²; H_ij = π · ⟨u_i|u_j⟩. Eigenfrequency diagonal with closure-quantum-scaled overlap coupling.

- Eigenvalues: ['-1.8292', '-1.3197', '7.6963']
- √eigenvalue ratios (predicted m_e:m_μ:m_τ = 0.0000 : 0.0000 : 0.0000); observed: 1 : 206.77 : 3477.23
- B1-with-probe-eigenvectors closure: residues ['0.8361', '0.7763', '0.8007'] π, circular spread 0.188151 rad
- D1-with-probe-eigenvectors closure: residues ['1.8656', '1.7877', '0.3467'] π, circular spread 1.756145 rad

#### `GH_E_symmetric_momentum`

H_ij = ⟨u_i|√max(ω̄²−V̄, 0)|u_j⟩, ω̄/V̄ symmetric averages. Diagonal recovers WKB action; off-diagonal couples through symmetric momentum kernel (the D2 form interpreted as H).

- Eigenvalues: ['0.0011', '0.0881', '2.8901']
- √eigenvalue ratios (predicted m_e:m_μ:m_τ = 1.0000 : 8.8008 : 50.3963); observed: 1 : 206.77 : 3477.23
- B1-with-probe-eigenvectors closure: residues ['0.7846', '0.8270', '0.8015'] π, circular spread 0.133116 rad
- D1-with-probe-eigenvectors closure: residues ['1.9081', '1.7452', '0.3467'] π, circular spread 1.889659 rad

## Verdict

**One variant closes D1 trivially (identity eigenvectors × D1's zero diagonal → Φ_radial = 0 for every species), but no variant closes the ledger non-trivially or matches the observed mass ladder.** Tangherlini ground-mode matrix elements alone cannot supply the lepton sector — the gap is structural, not parametric.

### Why the dynamic range is missing

All Tangherlini matrix elements between the (l, 0) ground modes live within an O(1) range: ω² spans 1.11 → 1.95 (factor 1.8), potential matrix elements ⟨u_i|V_j|u_j⟩ are O(0.1), overlap elements ⟨u_i|u_j⟩ are between 0.96 and 1.0 (the eigenfunctions live in nearly the same region of the tortoise grid). Diagonalizing any 3×3 matrix with entries in this range produces eigenvalues within an O(1) window — far short of the 5-orders-of-magnitude span of m_e² : m_τ².

The locked surrogate gets its dynamic range from the **integer closure quantum** uplift: 4·β = 100·(2π) for the lepton sector, applied as β·max(0, k−3)² to the τ row only. This term is geometric in origin (the integer 100 is the antipodal closure quantum count, not a fit parameter) but it is **not** a Tangherlini matrix element — it is a constraint from the antipodal-cavity closure condition operating in parallel to the radial sector.

**Operator structure of D1 is provably the projected radial Hamiltonian.** The identity

```
⟨u_i | V_j − V_i | u_j⟩  =  (ω_j² − ω_i²) · ⟨u_i | u_j⟩
```

(verified numerically to 1e-12 in this run) shows that D1's off-diagonal element is exactly the projection of the differential operator H_{l_j} − H_{l_i} onto the basis state (u_i, u_j). Variant GH_C uses this identity to build a closed-form geometric Hamiltonian whose entries are determined purely by the radial spectrum — no closure-quantum input. Its eigenvectors give the same B1-style closure result as C1 (because the GH_C eigenvalues are dominated by the ω² diagonal), confirming that the radial channel alone reaches the same ~0.3 rad floor without ever seeing the mass ladder.

### Conclusion

The lepton mass ladder and the closure ledger do **not** emerge together from Tangherlini ground-mode matrix elements. The radial spectrum is too gentle (factor 1.8 in ω² over l = 1, 3, 5) to support either the ~3500× mass span of m_τ/m_e or a residue redistribution that closes mod 2π. The locked surrogate's mass ladder relies on the antipodal-closure-quantum uplift β·max(0, k−3)², which is **not** a single matrix element of the radial operator — it is a separate geometric constraint. Closing the ledger and explaining the mass ratios at the same time requires both channels (radial AND closure-quantum), not either alone.