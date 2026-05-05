# Dynamic-phase probe — closure ledger

**Run:** 2026-05-03T05:59:48+00:00

Question: can the missing radial-channel residual for the best current S(k) candidates (C1, D1) be supplied by one of the natural BAM loop phases — moving throat transport, Hopf fibre loop, antipodal out-and-back, or worldline crossing? The probe scans each mechanism's natural parameter grid and reports whether any (mechanism, parameters) tightens the circular spread mod 2π.

## Mechanisms tested

| name | module | Δ(k) | parameter combinations |
|---|---|---|---:|
| `M1_moving_throat` | `geometrodynamics.embedding.transport` | `Δ(k) = m · k · (π/2)` | 7 |
| `M2_hopf_fibre_loop` | `geometrodynamics.hopf.connection` | `Δ(k) = m · k · π · cos(χ)` | 25 |
| `M3_antipodal_out_and_back` | `geometrodynamics.transaction.handshake` | `Δ(k) = α · k · 2π` | 5 |
| `M4_worldline_crossings` | `geometrodynamics.embedding.transport` | `Δ(k) = m · k · (π/2)   [same family as M1]` | 7 |
| `M5_antipodal_pair_with_hopf` | `geometrodynamics.transaction.handshake + geometrodynamics.hopf.connection` | `Δ(k) = m · k · (π/2) + n · π · cos(χ)` | 42 |

## Candidate: `C1_eigenvector_weighted_B1`

Baseline circular spread (no dynamic phase): **0.325749 rad** = 0.103689 π

### Best per mechanism

| mechanism | best parameters | Δ(k) (units of π) | new residues (units of π) | spread (rad) | helps? | closes? |
|---|---|---|---|---:|---|---|
| `M1_moving_throat` | m=2 | [+1.000, +3.000, +5.000] | [1.864, 1.788, 1.761] | 0.325749 | no | no |
| `M2_hopf_fibre_loop` | m=2, chi=1.0471975511965976 | [+1.000, +3.000, +5.000] | [1.864, 1.788, 1.761] | 0.325749 | no | no |
| `M3_antipodal_out_and_back` | alpha=1.0 | [+2.000, +6.000, +10.000] | [0.864, 0.788, 0.761] | 0.325749 | no | no |
| `M4_worldline_crossings` | m=2 | [+1.000, +3.000, +5.000] | [1.864, 1.788, 1.761] | 0.325749 | no | no |
| `M5_antipodal_pair_with_hopf` | m=2, n=1, chi=0.0 | [+2.000, +4.000, +6.000] | [0.864, 0.788, 0.761] | 0.325749 | no | no |

**Overall best for C1_eigenvector_weighted_B1**: `M2_hopf_fibre_loop` with params {'m': 2, 'chi': 1.0471975511965976} → spread = 0.325749 rad. (Helps: False; closes: False.)

## Candidate: `D1_potential_difference_phase`

Baseline circular spread (no dynamic phase): **0.576680 rad** = 0.183563 π

### Best per mechanism

| mechanism | best parameters | Δ(k) (units of π) | new residues (units of π) | spread (rad) | helps? | closes? |
|---|---|---|---|---:|---|---|
| `M1_moving_throat` | m=-2 | [-1.000, -3.000, -5.000] | [1.098, 0.914, 0.988] | 0.576680 | no | no |
| `M2_hopf_fibre_loop` | m=2, chi=1.0471975511965976 | [+1.000, +3.000, +5.000] | [1.098, 0.914, 0.988] | 0.576680 | no | no |
| `M3_antipodal_out_and_back` | alpha=0.0 | [+0.000, +0.000, +0.000] | [0.098, 1.914, 1.988] | 0.576680 | no | no |
| `M4_worldline_crossings` | m=-2 | [-1.000, -3.000, -5.000] | [1.098, 0.914, 0.988] | 0.576680 | no | no |
| `M5_antipodal_pair_with_hopf` | m=-2, n=-1, chi=0.0 | [-2.000, -4.000, -6.000] | [0.098, 1.914, 1.988] | 0.576680 | no | no |

**Overall best for D1_potential_difference_phase**: `M2_hopf_fibre_loop` with params {'m': 2, 'chi': 1.0471975511965976} → spread = 0.576680 rad. (Helps: False; closes: False.)

## Verdict

**No natural loop phase closes or tightens** the residual for either C1 or D1 within the scanned parameter grid. All four mechanisms (moving throat, Hopf fibre, antipodal out-and-back, worldline crossing) reduce to either k-independent shifts (universal across species, do not affect spread) or to multiples of k · π/2 (which give residues with spread π for k ∈ {1, 3, 5}, wider than C1's 0.326 rad). The radial-channel residual cannot be absorbed into any existing geometric channel — it must come from a new piece of physics.

### Why the natural mechanisms cannot fit

- **M1/M4 (spin-½ passes; T eigenvalue arg per crossing)**: both give Δ(k) = m · k · π/2 for integer m. For k ∈ {1, 3, 5}, the k·π/2 pattern mod 2π takes only the values [π/2, 3π/2, π/2] (m = 1) or [π, π, π] (m = 2) etc. The first has spread π; the second is universal (spread 0) and only shifts the universal value. Neither matches C1's residue pattern [0.864, 0.788, 0.761] π.
- **M2 (Hopf fibre loop)**: Δ(k) = m·k·π·cos(χ). At χ = 0 this is m·k·π, which mod 2π is universal for odd k. At χ = π/2 it vanishes entirely. Intermediate χ scales the universal piece without breaking the k-symmetry.
- **M3 (antipodal out-and-back)**: Δ(k) = α·k·2π. mod 2π this is α·k·2π mod 2π. For α = 1/2 (spin-½): k·π mod 2π = π for odd k (universal). For α = 1: 2π·k mod 2π = 0. Always universal.
- **General constraint**: any loop phase that scales as k·c with constant c gives k-dependent Δ only through the modular arithmetic, and the resulting patterns over k = (1, 3, 5) live on a finite lattice of multiples of π. C1's residues do not lie on this lattice.

Conclusion: the BAM loop-phase channels already in the ledger (antipodal closure k·2π, Hopf holonomy π·cos(χ), throat T² = π) exhaust what natural integer-quantized loop phases can supply. The remaining ~0.3 rad C1 residue is too small and the wrong shape to come from another instance of the same mechanisms.