# Multi-throat mechanics & the exchange kernel from the GR ψ–Φ–q soliton (PR #185)

**Run:** 2026-06-29T06:40:47+00:00

Derives the two-throat exchange kernel from GR — a geometric overlap `K(R)` from the actual #180 self-gravitating throat-soliton times the topological sign `−1` (Pin⁻ geon statistics) — and the multi-throat consequences (Pauli exclusion, the Fermi `n^{5/3}` EoS). *(QFT on the classical throat, not quantum gravity.)*

- **Kernel**: K_exchange(R) = (−1) × K(R) : Pin⁻ sign × GR soliton overlap
- **Spatial**: K(R) from the #180 soliton, decays over the soliton size (range)
- **Sign**: −1 (fermion), from the Pin⁻ geon statistics (exchange ≃ 2π rotation = T²=−I)
- **Fermi pressure**: the exchange hole ∝ K(R)²; N throats give E∝N^{5/3}, Γ=5/3 (the #172 EoS)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | derive the two-throat exchange kernel from the GR soliton | **PASS** |
| T2 | `T2_exchange_operator` | the exchange operator P (P²=1, eigenvalues ±1) | **PASS** |
| T3 | `T3_spatial_exchange_kernel` | the spatial kernel K(R) from the #180 soliton (overlap, range) | **PASS** |
| T4 | `T4_exchange_sign_pin_minus` | the exchange sign −1 from the Pin⁻ geon statistics | **PASS** |
| T5 | `T5_pauli_exclusion_at_coincidence` | Pauli: the antisymmetric two-throat state vanishes at coincidence | **PASS** |
| T6 | `T6_exchange_hole_and_fermi_pressure` | the exchange hole + the Fermi pressure (E∝N^{5/3}, Γ=5/3, #172) | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | MULTI_THROAT_EXCHANGE_KERNEL_FROM_GR | **PASS** |

## The spatial exchange kernel K(R) (two #180 throat-solitons)

| R | K̂(R) |
|---:|---:|
| 0.0 | 1.0 |
| 1.0 | 0.7911 |
| 2.0 | 0.4093 |
| 3.0 | 0.1529 |
| 4.0 | 0.0454 |
| 6.0 | 0.0026 |

Decays over the throat-soliton size (RMS ≈ 1.274) — the GR-geometric exchange range — then multiplied by the Pin⁻ sign `−1` for the full antisymmetric kernel.

## Verdict

**MULTI_THROAT_EXCHANGE_KERNEL_FROM_GR_OVERLAP_TIMES_PIN_MINUS_SIGN_ANTISYMMETRIC_FERMI_PRESSURE.** DERIVED — THE TWO-THROAT EXCHANGE KERNEL FROM GR. K_exchange(R) = (−1) × K(R): a geometric overlap times the Pin⁻ sign.

SPATIAL KERNEL. K(R), the overlap of two actual #180 self-gravitating throat-solitons, decays smoothly from K̂(0) = 1 over the soliton size (RMS ≈ 1.274) — a GR-geometric exchange range, not a postulated form factor (K̂: {'0.0': 1.0, '2.0': 0.4093, '4.0': 0.0454}).

SIGN. The swap large-diffeomorphism is homotopic to a 2π rotation of one throat (geon statistics); on the Pin⁻ throat that is T² = −I (½ tr T² = -1), so the exchange sign is -1 — the geometry selects the antisymmetric Fermi sector.

PAULI. The antisymmetric two-throat state vanishes at coincidence (max|Ψ₋(z,z)| = 0e+00) — two identical throats cannot occupy the same state — while the boson state does not.

FERMI PRESSURE. The exchange term ∝ K(R)² carves an exchange hole of GR range; the exclusion fills a degenerate Fermi tower, E ∝ N^1.667 ⟹ P ∝ n^{5/3} (Γ = 1.6667) — the Fermi EoS measured in #172. The GR-derived exchange kernel is the microscopic origin of the Fermi pressure of throat matter.
