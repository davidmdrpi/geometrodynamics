# Two-throat Hartree–Fock sandbox: direct plus exchange terms (PR #187)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The two-throat Hartree–Fock energy

PR #185 gave the two-throat exchange kernel `K_exchange(R) = (−1)·K(R)`; PR
#186 hardened its overlap kernels — the direct density-overlap `D(R)` and the
exchange amplitude-overlap `K(R)` — and separated the channels. This probe
convolves them with an interaction `V` to build the actual two-throat
**Hartree–Fock energy**, with both the direct (Hartree) and exchange terms:

```
E±(R) = J(R) ± K_ex(R) ,
  J(R)    = ∫∫ ρ_a(r₁) ρ_b(r₂) V(r₁−r₂)            (DIRECT / Hartree)
  K_ex(R) = ∫∫ τ(r₁) τ(r₂) V(r₁−r₂),  τ = φ_a φ_b   (EXCHANGE)
```

The orbitals are two rigid #180 throat-solitons at separation `R`; the
interaction `V` is a screened-photon (Yukawa) stand-in for the BAM
throat-fibre exchange (the unscreened Coulomb/photon `1/q²` is the #42–#44
result; the screening is a numerical regulator). The GR-selected Pin⁻ sign
`−1` (#185) puts the throats in the antisymmetric sector, so their physical
energy is `E₋ = J − K_ex`.

## The energies (3D-FFT Coulomb on the #180 orbitals)

| R | direct `J` | exchange `K_ex` | `K_ex/J` |
|---:|---:|---:|---:|
| 0.0 | 0.0393 | 0.0393 | 1.00 |
| 1.0 | 0.0310 | 0.0238 | 0.77 |
| 2.0 | 0.0173 | 0.0058 | 0.34 |
| 3.0 | 0.0086 | 0.0007 | 0.08 |
| 4.0 | 0.0044 | 0.0001 | 0.01 |
| 6.0 | 0.0013 | 0.0000 | 0.00 |

Both positive (repulsive `V`) and decaying; the **direct dominates** — the
exchange-to-direct ratio falls from 1 at contact to ~0 far apart (the exchange
has the shorter, overlap-set range).

## The exchange splitting — the fermion branch is lower

| R | `E₊ = J+K_ex` (boson) | `E₋ = J−K_ex` (fermion) | split `2 K_ex` |
|---:|---:|---:|---:|
| 0.5 | 0.0715 | 0.0024 | 0.0692 |
| 1.0 | 0.0548 | 0.0073 | 0.0475 |
| 2.0 | 0.0231 | 0.0115 | 0.0116 |
| 3.0 | 0.0093 | 0.0079 | 0.0014 |

The boson branch `E₊` lies **above** the fermion branch `E₋` by `2 K_ex`
everywhere — the exchange hole reduces the interaction energy of the
antisymmetric state. So the GR-selected Pin⁻ `−1` places the two throats in
the **lower-energy** branch; the splitting has a GR range set by the soliton
overlap (it dies as the throats separate).

## Pauli exclusion, exact

At coincidence (`R = 0`) the two orbitals are identical, so
`ρ_a = ρ_b = τ` and the direct and exchange energies are **equal**
(`J = K_ex = 0.0393`); the fermion energy

```
E₋ = J − K_ex = 0   (exactly)
```

— two identical throats at the same point have **zero interaction energy**,
the Pauli hole removing it entirely, while the boson has `E₊ = 2J` (bunching).
For a **contact** interaction `V = g δ` the cancellation is exact at all `R`:
`J = K_ex = g·D(R)` (the hardened #186 direct overlap), so `E₋ = 0` everywhere
and `E₊ = 2g·D(R)`. The `−1` the geometry selects is exactly the Pauli
exclusion of two throats — now at the level of the interaction energy.

## Controls + convergence

- **Far separation** (`R = 6`): both energies vanish (`J = 0.0013`,
  `K_ex ≈ 10⁻⁷ → 0`) — widely separated throats are distinguishable, no
  interaction, no exchange.
- **Grid convergence**: under 3D-grid refinement (`N = 64 → 80 → 96`) the
  direct `J(2)` and exchange `K_ex(2)` vary by `< 0.1%` — the FFT energies are
  grid-stable. (The precise values still carry the #186 soliton-profile ~3%
  uncertainty.)

## Honest scope

A **sandbox**: rigid #180 orbitals (not relaxed in each other's presence — the
self-consistent two-throat HF solve is a follow-up), a screened-photon
(Yukawa) regulated stand-in for the BAM Coulomb/photon exchange, and the
spatial (orbital) exchange only — the spin/statistics factor is the separate
Pin⁻ `−1` (#185). The energies are in code units (the interaction strength is
a scale, not calibrated to α). The qualitative HF structure — direct +
exchange, the `2 K_ex` splitting, the fermion-lower ordering, and `E₋ = 0` at
coincidence — is robust; the precise numbers carry the #186 soliton-profile
~3% uncertainty. Weak-field / semi-dynamical soliton.

## Reproduce

```bash
python -m experiments.closure_ledger.two_throat_hartree_fock_probe
# Verdict: TWO_THROAT_HARTREE_FOCK_DIRECT_PLUS_EXCHANGE_FERMION_BRANCH_LOWER_PAULI_ZERO_AT_COINCIDENCE
```
