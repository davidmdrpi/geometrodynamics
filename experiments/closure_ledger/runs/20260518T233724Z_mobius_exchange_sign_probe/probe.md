# Möbius throat exchange-sign probe

**Run:** 2026-05-18T23:37:24+00:00

Derives the QED Bhabha / Møller interference signs from BAM's non-orientable antipodal throat transport T = iσ_y, closing the last Fermi-statistics gap identified in PR #42 and partially resolved in PR #43.

## Key geometric identification

```
T = iσ_y = [[0, 1], [−1, 0]] = SU(2) ε tensor
```

## Derivation chain

`T = ε  →  antisymmetric singlet  →  swap eigenvalue = −1  →  one transposition Möbius sign = −1  →  Bhabha Wick + Møller Pauli derived`

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_throat_transport_from_repo` | 7/7 repo properties verified | **PASS** |
| T2 | `T2_T_equals_SU2_epsilon` | ||T − ε|| = 0.00e+00 | **PASS** |
| T3 | `T3_two_spinor_swap_operator` | eigenvalues: 3 × (+1), 1 × (−1) | **PASS** |
| T4 | `T4_singlet_via_epsilon_tensor` | singlet from T residual = 0.00e+00; antisym check = 0.00e+00 | **PASS** |
| T5 | `T5_mobius_exchange_sign` | 1 transposition → -1.0000 (target -1.0); 2 transpositions → 1.0000 | **PASS** |
| T6 | `T6_bhabha_s_t_diagrammatic_exchange` | BAM Möbius = -1, QED Wick = -1 | **PASS** |
| T7 | `T7_moller_t_u_diagrammatic_exchange` | BAM Möbius = -1, QED Pauli = -1 | **PASS** |
| T8 | `T8_bhabha_end_to_end` | max rel diff = 0.00e+00 | **PASS** |
| T9 | `T9_moller_end_to_end` | max rel diff = 0.00e+00 | **PASS** |

## T1: T = iσ_y from `geometrodynamics.embedding.transport`

```
  +0 +1
  -1 +0
```

Repo properties:
  - `T²=−I` ✓ (residual: 0.00e+00)
  - `T†T=I` ✓ (residual: 0.00e+00)
  - `det=1` ✓ (residual: 0.00e+00)
  - `TσzT†=−σz` ✓ (residual: 0.00e+00)
  - `TσxT†=−σx` ✓ (residual: 0.00e+00)
  - `T|↑⟩=−|↓⟩` ✓ (residual: 0.00e+00)
  - `T|↓⟩=+|↑⟩` ✓ (residual: 0.00e+00)

## T2: T = SU(2) antisymmetric ε tensor

||T − ε_SU(2)|| = `0.00e+00` (machine precision). `T + T^T = `0.00e+00` (exact antisymmetry).

## T3: Two-spinor swap operator P_12

Eigenvalues (sorted): `[-1.0, 1.0, 1.0, 1.0]`. Triplet (symmetric, +1) dimension = 3; singlet (antisymmetric, −1) dimension = 1.

## T4: Singlet via T = ε

Reshape T → 4-vector, normalise: matches canonical singlet (|↑↓⟩ − |↓↑⟩)/√2 to residual 0.00e+00.

## T5: Möbius exchange sign from T

One transposition → **-1.0000** (target -1.0). Two transpositions → **+1.0000** (target 1.0).

## T6: Bhabha s↔t diagrammatic exchange

s-channel pairing: `((1, 2), (3, 4))`; t-channel pairing: `((1, 3), (2, 4))`; transposition: `(p_2 ↔ p_3)` (1 swap). BAM Möbius sign = **-1**; QED Wick sign = **-1**.

## T7: Møller t↔u diagrammatic exchange

t-channel pairing: `((1, 3), (2, 4))`; u-channel pairing: `((1, 4), (2, 3))`; transposition: `(p_3 ↔ p_4)` (1 swap of identical fermions). BAM Möbius sign = **-1**; QED Pauli sign = **-1**.

## T8: Bhabha end-to-end |M̄|²/(8e⁴)

Möbius sign = `-1` (one transposition).

| θ | s | t | u | BAM | QED | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 391.7307 | 391.7307 | 0.00e+00 |
| 60 | 4.00 | -1.000 | -3.000 | 21.1250 | 21.1250 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 4.5000 | 4.5000 | 0.00e+00 |
| 120 | 4.00 | -3.000 | -1.000 | 2.3472 | 2.3472 | 0.00e+00 |
| 150 | 4.00 | -3.732 | -0.268 | 2.0193 | 2.0193 | 0.00e+00 |

## T9: Møller end-to-end |M̄|²/(8e⁴)

Möbius sign = `-1` (one transposition).

| θ | s | t | u | BAM | QED | rel diff |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 4.00 | -0.268 | -3.732 | 450.0000 | 450.0000 | 0.00e+00 |
| 60 | 4.00 | -1.000 | -3.000 | 37.5556 | 37.5556 | 0.00e+00 |
| 90 | 4.00 | -2.000 | -2.000 | 18.0000 | 18.0000 | 0.00e+00 |
| 120 | 4.00 | -3.000 | -1.000 | 37.5556 | 37.5556 | 0.00e+00 |
| 150 | 4.00 | -3.732 | -0.268 | 450.0000 | 450.0000 | 0.00e+00 |

## Verdict

**SIGNS_DERIVED.** SIGNS DERIVED FROM MÖBIUS THROAT TRANSPORT. The QED Wick / Pauli interference signs for Bhabha (s↔t) and Møller (t↔u) emerge automatically from BAM's T = iσ_y non-orientable antipodal throat transport:
  (1) T = iσ_y = [[0, 1], [−1, 0]] is the SU(2) antisymmetric ε tensor (T2);
  (2) the antisymmetric two-spinor singlet state — the Fermi state — is generated directly by ε (T4);
  (3) the Möbius exchange sign for one transposition is −1, derived as the eigenvalue of the swap operator on the antisymmetric subspace (T5);
  (4) Bhabha s↔t and Møller t↔u differ by exactly one fermion-leg transposition, giving σ_M = −1 in both cases (T6, T7), matching the QED Wick and Pauli rules respectively.
Combined with PR #43 (BAM SU(2) Pauli/Weyl traces giving the diagonal-channel numerators (s²+u²), (u²+t²), (s²+t²)), the end-to-end Bhabha and Møller scalar intensities |M̄|²/(8e⁴) are reproduced from BAM geometric ingredients alone, to machine precision (T8, T9). The "Fermi-statistics overlay" identified in PR #42 is not foreign to BAM — it is the natural action of the non-orientable throat transport.

## What this closes

- **The Bhabha/Møller derivation thread (PRs #42, #43, this)**: all of Bhabha and Møller tree-level scalar intensity now derives from BAM-geometric ingredients (Hopf SU(2) Pauli traces for diagonals; non-orientable T = iσ_y throat transport for interference signs), with no remaining Fermi-statistics overlay.

## What this leaves open

- **Virtual-photon propagator beyond 1/q²**: still an ansatz; a BAM throat-fibre propagator derivation is a separate probe target.
- **Loop corrections**: tree-level only.
