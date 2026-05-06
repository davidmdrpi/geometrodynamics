# Quark β: sub-block stability sub-probe

**Run:** 2026-05-06T04:35:13+00:00
**k_5 = 5**, **N_c = 3**, **baseline N_q = 466** = 100 + 350 + 15 + 1.

Asks the focused question: **does any sub-block of the baseline decomposition remain perturbation-stable across the §8 ablation table?** Two test families are run: (a) constant-shift subtraction (sanity-check that subtracting a baseline-locked constant cannot reduce variance — confirmed below), and (b) modular-arithmetic invariants (which IS a non-trivial test, since residue-class membership can survive ablations even when the absolute value drifts).

## Constant-shift subtraction (null-hypothesis check)

Each row subtracts a baseline sub-block S from N_q across every ablation and reports the residual envelope. Constant-shift subtraction cannot reduce variance (Var(N − c) = Var(N)), so this is a sanity check that no sub-block is differentially stable in the trivial sense.

| sub-block S | formula | S value | residual range | width | std dev |
|---|---|---:|---|---:|---:|
| `S_zero` | `0 (no subtraction)` | 0 | [432, 510] | 78 | 21.87 |
| `S_lepton` | `(k_5−1)·k_5²` | 100 | [332, 410] | 78 | 21.87 |
| `S_lepton_plus_boundary` | `(k_5−1)·k_5² + (N_c−2)` | 101 | [331, 409] | 78 | 21.87 |
| `S_shell` | `2·k_5²(k_5+2)` | 350 | [82, 160] | 78 | 21.87 |
| `S_color_x_k5` | `N_c·k_5` | 15 | [417, 495] | 78 | 21.87 |
| `S_color_plus_boundary` | `N_c·k_5 + (N_c−2)` | 16 | [416, 494] | 78 | 21.87 |
| `S_lepton_plus_shell` | `(k_5−1)·k_5² + 2·k_5²(k_5+2)` | 450 | [-18, 60] | 78 | 21.87 |
| `S_lepton_plus_color` | `(k_5−1)·k_5² + N_c·k_5` | 115 | [317, 395] | 78 | 21.87 |
| `S_shell_plus_color` | `2·k_5²(k_5+2) + N_c·k_5` | 365 | [67, 145] | 78 | 21.87 |
| `S_all_but_boundary` | `(k_5−1)·k_5² + 2·k_5²(k_5+2) + N_c·k_5` | 465 | [-33, 45] | 78 | 21.87 |
| `S_full_baseline` | `(k_5−1)·k_5² + 2·k_5²(k_5+2) + N_c·k_5 + (N_c−2)` | 466 | [-34, 44] | 78 | 21.87 |

As expected: **all 11 sub-block subtractions give the same drift width (78 units)**. No sub-block is differentially stable in the constant-shift sense — they all carry the same drift envelope. This rules out the simple reading 'subtract the lepton sub-block to expose a tighter compensator.'

## Modular-arithmetic invariants

Tests whether N_q ≡ baseline_residue (mod m) across all 12 logged ablations. Each row is one modulus m, with the structural interpretation of m. A modular invariance that holds across every ablation is a genuine sub-block stability — even if the absolute value of N_q wanders.

| modulus m | residues across 12 ablations | invariant? | interpretation |
|---:|---|:---:|---|
| 2 | [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] | **✓** | Z₂ partition class — two partition sectors (p = ±) contribute equally to β |
| 3 | [1, 1, 1, 2, 0, 0, 2, 0, 2, 2, 2, 0] | — | SU(3) color triplet (N_c = 3) |
| 4 | [2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 0, 2] | — | spinor double cover or 4-fold winding |
| 5 | [1, 1, 1, 1, 4, 4, 2, 2, 4, 4, 0, 0] | — | k_5 modulus (heaviest shell label) |
| 6 | [4, 4, 4, 2, 0, 0, 2, 0, 2, 2, 2, 0] | — | (k_5 + 1) modulus — closure-quantum on shell+1 |
| 10 | [6, 6, 6, 6, 4, 4, 2, 2, 4, 4, 0, 0] | — | 2 · k_5 (Z₂ × k_5 combined) |
| 15 | [1, 1, 1, 11, 9, 9, 2, 12, 14, 14, 5, 0] | — | N_c · k_5 |
| 25 | [16, 16, 16, 1, 24, 24, 7, 7, 19, 19, 15, 10] | — | k_5² |
| 26 | [24, 24, 24, 8, 6, 6, 14, 16, 0, 0, 24, 16] | — | k_5² + 1 |
| 100 | [66, 66, 66, 76, 74, 74, 82, 32, 94, 94, 40, 10] | — | lepton closure quantum |

### Preserved modular invariants

- **N_q ≡ 0 (mod 2)** holds across all 12 logged ablations. Interpretation: Z₂ partition class — two partition sectors (p = ±) contribute equally to β.

## Half-partition diagnostic

**N_q is even across all 12 logged §8 ablations.** Per-partition closure-quantum count n_part = N_q / 2:

- baseline n_part = 233
- §8 envelope n_part ∈ [216, 255], drift ±22 units

The factor-of-2 reflects the v3 ansatz's two-partition Hamiltonian basis {(k, +), (k, −)}: each Z₂ partition sector contributes equally to the total β, forcing N_q to be twice an integer count. The two-ness IS the topological invariant; the per-partition count itself drifts under perturbations as the fit compensator.

## Verdict

**The only perturbation-stable structural property of the baseline decomposition is parity** (`N_q ≡ 0 (mod 2)`, i.e. `mod 2` invariance). All other modular structures and all constant-shift sub-block subtractions vary across §8 ablations.

Structurally: the two-Z₂-partition multiplicity is the topological invariant; the specific N_q value within the even integers is the compensator. This refines the decomposition probe's reading from

  `N_q = ((k_5−1)·k_5 + 2·k_5(k_5+2) + N_c)·k_5 + (N_c−2)`

to the more honest

  `N_q = 2 · n_part`,   where `n_part` is a phenomenological per-partition closure-quantum count whose specific value (233 at baseline) is NOT topologically locked.

The 2× partition factor is the sole structural piece; the rest is fit. This is the cleanest summary of the §8 compensator behavior under a sub-block-stability lens.