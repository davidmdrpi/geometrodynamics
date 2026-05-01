# Closure-phase ledger — candidate comparison

**Run:** 2026-05-01T06:09:22+00:00
**Transport convention:** `T2_sign_flip`
**χ (Hopf fibre):** 0.0

## Status table

| candidate | layer | result | per-lepton mod 2π (units of π) | linear spread | circular spread | universal value |
|---|---|---|---|---|---|---|
| `none` | layer1 | **PASS** | [0.000000, 0.000000, 0.000000] | 1.421e-14 | 1.421e-14 | 0.000000π |
| `A_lowest_radial_per_l` | layer2 | **FAIL** | [0.881876, 1.652620, 0.413097] | 3.894e+00 | 3.862e+00 | — |
| `B1_single_angular_mode` | layer2 | **FAIL** | [0.881876, 0.770744, 0.760477] | 3.814e-01 | 3.814e-01 | — |
| `B2_single_radial_excitation` | layer2 | **FAIL** | [0.881876, 1.994352, 0.999437] | 3.495e+00 | 3.158e+00 | — |
| `C1_eigenvector_weighted_B1` | layer2 | **FAIL** | [0.864195, 0.788396, 0.760506] | 3.257e-01 | 3.257e-01 | — |
| `C2_eigenvector_weighted_B2` | layer2 | **FAIL** | [1.059227, 1.817711, 0.998727] | 2.573e+00 | 2.573e+00 | — |
| `C1_maslov_standard` | layer2 | **FAIL** | [0.364195, 0.288396, 0.260506] | 3.257e-01 | 3.257e-01 | — |
| `B2_maslov_standard` | layer2 | **FAIL** | [0.381876, 1.994352, 0.999437] | 5.066e+00 | 3.158e+00 | — |
| `C2_maslov_standard` | layer2 | **FAIL** | [0.638760, 1.738289, 0.998617] | 3.454e+00 | 3.454e+00 | — |
| `D0_overlap_phase` | layer2 | **FAIL** | [1.775308, 0.283574, 0.941118] | 4.686e+00 | 3.663e+00 | — |
| `D1_potential_difference_phase` | layer2 | **FAIL** | [0.097827, 1.914264, 1.987910] | 5.938e+00 | 5.767e-01 | — |
| `D2_symmetrized_momentum_phase` | layer2 | **FAIL** | [0.533647, 0.093822, 0.320892] | 1.382e+00 | 1.382e+00 | — |

## Mode-map comparison

| candidate | k=1 | k=3 | k=5 |
|---|---|---|---|
| `A_lowest_radial_per_l` | {(1,0)} | {(1,0), (3,0)} | {(1,0), (3,0), (5,0)} |
| `B1_single_angular_mode` | {(1,0)} | {(3,0)} | {(5,0)} |
| `B2_single_radial_excitation` | {(1,0)} | {(1,1)} | {(1,2)} |
| `C1_eigenvector_weighted_B1` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |
| `C2_eigenvector_weighted_B2` | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} |
| `C1_maslov_standard` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |
| `B2_maslov_standard` | {(1,0)} | {(1,1)} | {(1,2)} |
| `C2_maslov_standard` | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} |
| `D0_overlap_phase` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |
| `D1_potential_difference_phase` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |
| `D2_symmetrized_momentum_phase` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |

## Per-candidate one-liners

- `none`: Layer 1 PASS: lepton ledger universal mod 2π at 0.000π. Layer 2 DISABLED (sk_candidate='none').
- `A_lowest_radial_per_l`: Layer 2 FALSIFIES candidate 'A_lowest_radial_per_l': lepton ledger does NOT close universally mod 2π (spread 3.894e+00).
- `B1_single_angular_mode`: Layer 2 FALSIFIES candidate 'B1_single_angular_mode': lepton ledger does NOT close universally mod 2π (spread 3.814e-01).
- `B2_single_radial_excitation`: Layer 2 FALSIFIES candidate 'B2_single_radial_excitation': lepton ledger does NOT close universally mod 2π (spread 3.495e+00).
- `C1_eigenvector_weighted_B1`: Layer 2 FALSIFIES candidate 'C1_eigenvector_weighted_B1': lepton ledger does NOT close universally mod 2π (spread 3.257e-01).
- `C2_eigenvector_weighted_B2`: Layer 2 FALSIFIES candidate 'C2_eigenvector_weighted_B2': lepton ledger does NOT close universally mod 2π (spread 2.573e+00).
- `C1_maslov_standard`: Layer 2 FALSIFIES candidate 'C1_maslov_standard': lepton ledger does NOT close universally mod 2π (spread 3.257e-01).
- `B2_maslov_standard`: Layer 2 FALSIFIES candidate 'B2_maslov_standard': lepton ledger does NOT close universally mod 2π (spread 5.066e+00).
- `C2_maslov_standard`: Layer 2 FALSIFIES candidate 'C2_maslov_standard': lepton ledger does NOT close universally mod 2π (spread 3.454e+00).
- `D0_overlap_phase`: Layer 2 FALSIFIES candidate 'D0_overlap_phase': lepton ledger does NOT close universally mod 2π (spread 4.686e+00).
- `D1_potential_difference_phase`: Layer 2 FALSIFIES candidate 'D1_potential_difference_phase': lepton ledger does NOT close universally mod 2π (spread 5.938e+00).
- `D2_symmetrized_momentum_phase`: Layer 2 FALSIFIES candidate 'D2_symmetrized_momentum_phase': lepton ledger does NOT close universally mod 2π (spread 1.382e+00).

## Interpretation

- Layer 1 (`sk_candidate='none'`): closed-form topological ledger is universal mod 2π across leptons. PASS.
- Layer 2A (`A_lowest_radial_per_l`, WKB convention): cumulative odd-l ground modes do not preserve closure mod 2π. FAIL — rejects this radial-mode interpretation of lepton depth.
- Layer 2B1 (`B1_single_angular_mode`, WKB convention): single l=k ground mode per generation. The per-row Φ values are exactly the candidate-A summands; if A's cumulative residues are non-universal, the same per-mode numbers cannot be universal alone unless they happen to coincide mod 2π. FAIL.
- Layer 2B2 (`B2_single_radial_excitation`, WKB convention): single l=1 ladder. Asymptotic WKB gives Φ(1, n) → (n + 1) π, which produces a parity pattern across generations rather than a single universal value. FAIL.
- Layer 2C1 (`C1_eigenvector_weighted_B1`, WKB convention): B1 modes weighted by |v_species,i|² from the locked lepton generation block. The eigenvector mixing tightens the spread relative to B1 (the {1,3} mixing pulls e and μ residues toward each other) but does NOT close them to a universal value. FAIL — under this WKB convention, the surrogate Hamiltonian's own eigenvectors do not supply the missing bridge.
- Layer 2C2 (`C2_eigenvector_weighted_B2`, WKB convention): B2 modes weighted by the same eigenvectors. The (n+1)π spacing between B2 modes dominates the eigenvector mixing, leaving residues distributed across [0, 2π). FAIL.
- Layer 2 C1+Maslov (`C1_maslov_standard`, Bohr-Sommerfeld convention): same C1 mode set and weights, but each integrated phase is shifted by −π/2 per detected classical turning point (sign change of ω² − V_eff inside the tortoise grid). When the turning-point count is uniform across the B1 ground modes the Maslov shift is a uniform offset and the spread inherits from C1; when the count varies across modes the correction can in principle redistribute residues. The radial-detail table in each sub-run reports the per-mode `n_turning_points` so the regime can be read directly.
- Layer 2 B2+Maslov (`B2_maslov_standard`, Bohr-Sommerfeld convention): single-mode B2 ladder (l = 1, n = (k−1)/2) with the standard −π/2 turning-point shift. The B2 ladder is the natural test for differential Maslov: the (l=1, n=0) ground mode has a centrifugal-barrier soft turning point while the n ≥ 1 excitations sit above the barrier and have no interior turning points. Whether the differential shift collapses the (n+1)π parity pattern is read directly off the per-row Φ + Δ_maslov in the breakdown table.
- Layer 2 C2+Maslov (`C2_maslov_standard`, Bohr-Sommerfeld convention): C2's eigenvector-weighted B2 ladder with the same standard Maslov shift. Combines C-family weighting (lifts the B2 single-mode degeneracy) with the per-mode turning-point correction. Differential N_turning across the B2 ladder lets the eigenvector weights redistribute the Maslov shift across species — the most aggressive composition of the two known structural levers.
- Layer 2 D-family (`D0_overlap_phase`, `D1_potential_difference_phase`, `D2_symmetrized_momentum_phase`): operator-valued radial phase. A 3×3 Hermitian Φ matrix is built on the depth-basis B1 modes from the Tangherlini eigenfunctions on the tortoise grid, then contracted as v_species^T Φ v_species using the same locked lepton eigenvectors that drive C1/C2. D0 uses Φ_ij = π·⟨u_i|u_j⟩ (L²-normalized overlap with a symmetric phase scale); D1 uses the canonical quark transport matrix element Φ_ij = ⟨u_i|V_j−V_i|u_j⟩ with the upper triangle mirrored to make Φ Hermitian (so Φ_ij = Φ_ji = the transport coupling, diagonal = 0); D2 uses the symmetrized-momentum kernel Φ_ij = ⟨u_i|√max(ω̄²−V̄, 0)|u_j⟩ with ω̄² = (ω_i² + ω_j²)/2 and V̄ = (V_i + V_j)/2. The off-diagonal entries are the new degree of freedom relative to scalar-per-mode candidates: cross-terms 2·Σ_{i<j} v_i v_j Φ_ij can redistribute phase across species when the eigenvector signs disagree.
- Reading the D-family results: the residues mod 2π may straddle the wrap (e.g. one near 0, two near 2π), in which case the linear max−min `spread` over-states the actual closure error. The `circular spread` column is the honest closure metric — 2π minus the largest gap between consecutive residues on the circle — and is what `universal` is keyed off.
