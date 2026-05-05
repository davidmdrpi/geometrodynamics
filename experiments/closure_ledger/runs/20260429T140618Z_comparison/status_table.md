# Closure-phase ledger — candidate comparison

**Run:** 2026-04-29T14:06:18+00:00
**Transport convention:** `T2_sign_flip`
**χ (Hopf fibre):** 0.0

## Status table

| candidate | layer | result | per-lepton mod 2π (units of π) | spread | universal value |
|---|---|---|---|---|---|
| `none` | layer1 | **PASS** | [0.000000, 0.000000, 0.000000] | 1.421e-14 | 0.000000π |
| `A_lowest_radial_per_l` | layer2 | **FAIL** | [0.881876, 1.652620, 0.413097] | 3.894e+00 | — |
| `B1_single_angular_mode` | layer2 | **FAIL** | [0.881876, 0.770744, 0.760477] | 3.814e-01 | — |
| `B2_single_radial_excitation` | layer2 | **FAIL** | [0.881876, 1.994352, 0.999437] | 3.495e+00 | — |
| `C1_eigenvector_weighted_B1` | layer2 | **FAIL** | [0.864195, 0.788396, 0.760506] | 3.257e-01 | — |
| `C2_eigenvector_weighted_B2` | layer2 | **FAIL** | [1.059227, 1.817711, 0.998727] | 2.573e+00 | — |

## Mode-map comparison

| candidate | k=1 | k=3 | k=5 |
|---|---|---|---|
| `A_lowest_radial_per_l` | {(1,0)} | {(1,0), (3,0)} | {(1,0), (3,0), (5,0)} |
| `B1_single_angular_mode` | {(1,0)} | {(3,0)} | {(5,0)} |
| `B2_single_radial_excitation` | {(1,0)} | {(1,1)} | {(1,2)} |
| `C1_eigenvector_weighted_B1` | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} | {(1,0), (3,0), (5,0)} |
| `C2_eigenvector_weighted_B2` | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} | {(1,0), (1,1), (1,2)} |

## Per-candidate one-liners

- `none`: Layer 1 PASS: lepton ledger universal mod 2π at 0.000π. Layer 2 DISABLED (sk_candidate='none').
- `A_lowest_radial_per_l`: Layer 2 FALSIFIES candidate 'A_lowest_radial_per_l': lepton ledger does NOT close universally mod 2π (spread 3.894e+00).
- `B1_single_angular_mode`: Layer 2 FALSIFIES candidate 'B1_single_angular_mode': lepton ledger does NOT close universally mod 2π (spread 3.814e-01).
- `B2_single_radial_excitation`: Layer 2 FALSIFIES candidate 'B2_single_radial_excitation': lepton ledger does NOT close universally mod 2π (spread 3.495e+00).
- `C1_eigenvector_weighted_B1`: Layer 2 FALSIFIES candidate 'C1_eigenvector_weighted_B1': lepton ledger does NOT close universally mod 2π (spread 3.257e-01).
- `C2_eigenvector_weighted_B2`: Layer 2 FALSIFIES candidate 'C2_eigenvector_weighted_B2': lepton ledger does NOT close universally mod 2π (spread 2.573e+00).

## Interpretation

- Layer 1 (`sk_candidate='none'`): closed-form topological ledger is universal mod 2π across leptons. PASS.
- Layer 2A (`A_lowest_radial_per_l`, WKB convention): cumulative odd-l ground modes do not preserve closure mod 2π. FAIL — rejects this radial-mode interpretation of lepton depth.
- Layer 2B1 (`B1_single_angular_mode`, WKB convention): single l=k ground mode per generation. The per-row Φ values are exactly the candidate-A summands; if A's cumulative residues are non-universal, the same per-mode numbers cannot be universal alone unless they happen to coincide mod 2π. FAIL.
- Layer 2B2 (`B2_single_radial_excitation`, WKB convention): single l=1 ladder. Asymptotic WKB gives Φ(1, n) → (n + 1) π, which produces a parity pattern across generations rather than a single universal value. FAIL.
- Layer 2C1 (`C1_eigenvector_weighted_B1`, WKB convention): B1 modes weighted by |v_species,i|² from the locked lepton generation block. The eigenvector mixing tightens the spread relative to B1 (the {1,3} mixing pulls e and μ residues toward each other) but does NOT close them to a universal value. FAIL — under this WKB convention, the surrogate Hamiltonian's own eigenvectors do not supply the missing bridge.
- Layer 2C2 (`C2_eigenvector_weighted_B2`, WKB convention): B2 modes weighted by the same eigenvectors. The (n+1)π spacing between B2 modes dominates the eigenvector mixing, leaving residues distributed across [0, 2π). FAIL.
