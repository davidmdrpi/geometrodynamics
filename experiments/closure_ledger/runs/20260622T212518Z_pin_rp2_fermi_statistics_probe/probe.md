# Pin⁻ on the throat's RP² mouth: the exchange sign and the Fermi equation of state (PR #170)

**Run:** 2026-06-22T21:25:18+00:00

Takes the Pin⁻ structure on the non-orientable throat mouth (#169) and shows it **delivers** the −1 exchange sign and the Fermi equation of state — the calculation that makes the topology matter. *(QFT on the classical throat, not quantum gravity.)*

- **Pin structure**: RP² admits Pin⁻ only (w₁=w₂=a; Spin/Pin⁺ excluded)
- **Exchange sign**: −1 (spinor 2π = −1; exchange ≃ 2π rotation, Finkelstein–Rubinstein)
- **Fermi EoS**: P=⅔u Γ=5/3 (NR); P=⅓u Γ=4/3 (UR); T=0 degeneracy pressure > 0
- **Scope**: computed: Pin⁻ class, spinor sign, EoS; cited: the FR exchange↔rotation homotopy

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | goal: deliver the exchange sign + the Fermi EoS | **PASS** |
| T2 | `T2_rp2_is_pin_minus` | RP² carries Pin⁻ only (SW: w₁=w₂=a) | **PASS** |
| T3 | `T3_spinor_2pi_is_minus_one` | the Pin⁻ spinor is spin-½ (2π = −1, 4π = +1) | **PASS** |
| T4 | `T4_exchange_sign_minus_one` | exchange sign = −1 (Finkelstein–Rubinstein); antisymmetry | **PASS** |
| T5 | `T5_pauli_exclusion_fermi_dirac` | Pauli exclusion: occupation n_p ∈ {0,1} | **PASS** |
| T6 | `T6_fermi_equation_of_state` | Fermi EoS: P=⅔u Γ=5/3 (NR); P=⅓u Γ=4/3 (UR) | **PASS** |
| T7 | `T7_honesty_and_scope` | honest scope (computed vs the cited FR homotopy) | **PASS** |
| T8 | `T8_assessment` | PIN_MINUS_ON_RP2_DELIVERS_FERMI_EXCHANGE_SIGN_AND_EOS | **PASS** |

## The Fermi equation of state (T = 0)

| regime | P/u | polytropic Γ = d ln P / d ln n |
|---|---:|---:|
| non-relativistic (ε = p²/2m) | 0.6667 (= 2/3) | 1.6667 (= 5/3) |
| ultra-relativistic (ε = pc) | 0.3333 (= 1/3) | 1.3333 (= 4/3) |
| Bose gas at T=0 (contrast) | — | pressure 0 (no degeneracy) |

## Verdict

**PIN_MINUS_ON_RP2_DELIVERS_FERMI_EXCHANGE_SIGN_AND_EOS.** DELIVERED. The Pin⁻ structure on the throat's RP² mouth yields the −1 exchange sign and the Fermi equation of state — the deferred calculation, done.

THE SPINOR STRUCTURE. RP² has Stiefel–Whitney classes w₁ = w₂ = a, so it admits no Spin and no Pin⁺ structure — only Pin⁻ (w₂ + w₁² = 0). The throat mouth therefore has a unique, definite spinor structure, the non-orientable analogue of Spin.

THE EXCHANGE SIGN. The Pin⁻ spinor is spin-½: a 2π rotation acts as R(2π) = −I (and 4π = +I). By Finkelstein–Rubinstein the two-throat exchange is homotopic to a 2π rotation of one throat, so the exchange sign is −1 and the two-throat wavefunction is ANTISYMMETRIC — the spin-statistics connection realised by the SAME holonomy that gives 2π = −1.

THE EQUATION OF STATE. Antisymmetry ⟹ Pauli exclusion (occupation 0 or 1) ⟹ filling the Fermi sphere gives the degenerate EoS: P = ⅔u, Γ = 1.667 = 5/3 (non-relativistic) and P = ⅓u, Γ = 1.333 = 4/3 (ultra-relativistic), with a strictly positive T=0 degeneracy pressure — the support of white dwarfs and neutron stars — that a Bose gas (which collapses to p=0, P=0) lacks.

SCOPE. Computed: the Pin⁻ classification, the spinor 2π sign, and the Fermi EoS. Cited: the Finkelstein–Rubinstein exchange↔rotation homotopy — the one configuration-space theorem linking the throat's internal Pin holonomy to the physical exchange.
