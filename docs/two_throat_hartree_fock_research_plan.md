# Two-throat Hartree–Fock sandbox: direct plus exchange terms (PR #187)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The two-throat Hartree–Fock energy (overlap-normalized)

PR #185 gave the two-throat exchange kernel `K_exchange(R) = (−1)·K(R)`; PR
#186 hardened its overlap kernels — the direct density-overlap `D(R)` and the
exchange amplitude-overlap `K(R)` — and separated the channels. This probe
convolves them with an interaction `V` to build the actual two-throat
**Hartree–Fock energy**, with both the direct (Hartree) and exchange terms.

Two displaced throats are **non-orthogonal** — their orbital overlap
`S(R) = ⟨φ_a|φ_b⟩ ≠ 0` — so the properly normalized two-body energy carries
the overlap normalization `(1 ± S²)`:

```
E±(R) = (J(R) ± K_ex(R)) / (1 ± S²) ,
  S(R)    = ⟨φ_a|φ_b⟩                            (orbital overlap)
  J(R)    = ∫∫ ρ_a(r₁) ρ_b(r₂) V(r₁−r₂)          (DIRECT numerator)
  K_ex(R) = ∫∫ τ(r₁) τ(r₂) V(r₁−r₂),  τ = φ_a φ_b (EXCHANGE numerator)
```

`J` and `K_ex` are the **unnormalized HF numerators**; the physical energies
are `E±` above. The orbitals are two rigid #180 throat-solitons at separation
`R`; `V` is a screened-photon (Yukawa) stand-in for the BAM throat-fibre
exchange (the unscreened Coulomb/photon `1/q²` is #42–#44; the screening is a
regulator). The GR-selected Pin⁻ sign `−1` (#185) puts the throats in the
antisymmetric sector, so their physical energy is `E₋ = (J − K_ex)/(1 − S²)`.

## The overlap, the numerators, and the normalized energies

| R | `S(R)` | `J` (direct num.) | `K_ex` (exch. num.) | `E₊=(J+K_ex)/(1+S²)` | `E₋=(J−K_ex)/(1−S²)` |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.942 | 0.0369 | 0.0346 | 0.0379 | 0.0211 |
| 1.0 | 0.791 | 0.0310 | 0.0238 | 0.0337 | 0.0194 |
| 1.5 | 0.597 | 0.0238 | 0.0130 | 0.0271 | 0.0168 |
| 2.0 | 0.409 | 0.0173 | 0.0058 | 0.0198 | 0.0138 |
| 3.0 | 0.153 | 0.0086 | 0.0007 | 0.0091 | 0.0081 |
| 4.0 | 0.045 | 0.0044 | 0.0001 | 0.0044 | 0.0044 |

The overlap `S` decays from 1 (coincidence) to ~0 (far apart, orthogonal);
`J, K_ex` are positive (repulsive `V`) and decaying, the **direct dominating**
(`K_ex/J` from 1 at contact to ~0 far apart).

## Fermion lower — for a repulsive interaction

For the **repulsive** screened `V`, the antisymmetric (fermion) branch `E₋`
sits **below** the symmetric (boson) branch `E₊` at every finite separation
(see the table) — the exchange hole reduces the interaction energy of the
antisymmetric state, so the GR-selected Pin⁻ `−1` places the two throats in
the **lower** branch. The gap closes as the overlap dies (`S → 0 ⟹ E₊ ≈ E₋`,
distinguishable, nearly degenerate).

This ordering is **scoped to a repulsive `V`**: with an attractive interaction
the exchange term flips the comparison. (Note the normalization itself works
*against* it — dividing `E₋` by `1 − S² < 1` raises it and `E₊` by `1 + S² > 1`
lowers it — yet the fermion branch remains lower across the tested range.)

## Pauli at coincidence — the zero vector (forbidden)

As `R → 0` the two orbitals become identical (`S → 1`), and the antisymmetric
combination

```
Ψ₋ = (φ_a φ_b − φ_b φ_a) / √(2(1 − S²))
```

has **both** a vanishing numerator (`J − K_ex → 0`) **and** a vanishing
normalization (`1 − S² → 0`): it is the **zero vector**. Two identical
fermions cannot occupy the same orbital, so the antisymmetric state is
**Pauli-forbidden** — *not* a state with zero interaction energy. The
symmetric (boson) state survives (`E₊ = (J + K_ex)/(1 + S²)`, bunching).

For a **contact** interaction `V = g δ` the numerator `J − K_ex = 0` at all
`R` (`J = K_ex = g·D(R)`, the hardened #186 direct overlap), so the
antisymmetric energy `E₋ = 0` at every finite separation (the exchange
exactly cancels the direct — the Pauli hole removes the contact interaction),
the state being forbidden only at exact coincidence.

## Controls + convergence

- **Far separation** (`R = 6`): `S = 0.003 → 0` (orthogonal, distinguishable
  throats) and both numerators vanish (`J = 0.0013`, `K_ex ≈ 10⁻⁷`), so
  `E₊ ≈ E₋` (no interaction, no exchange splitting).
- **Grid convergence**: under 3D-grid refinement (`N = 64 → 80 → 96`) the
  direct `J(2)` and exchange `K_ex(2)` vary by `< 0.1%`. (The precise values
  still carry the #186 soliton-profile ~3% uncertainty.)

## Honest scope

A **sandbox**: rigid #180 orbitals (not relaxed in each other's presence — the
self-consistent two-throat HF solve is a follow-up), a screened-photon
(Yukawa) regulated stand-in for the BAM Coulomb/photon exchange, and the
spatial (orbital) exchange only — the spin/statistics factor is the separate
Pin⁻ `−1` (#185). The energies are in code units (the interaction strength is
a scale, not calibrated to α). The overlap-normalized structure — `S(R)`, the
direct/exchange numerators, the `(1 ± S²)` normalization, the fermion-lower
ordering **for a repulsive `V`**, and the **forbidden (zero-vector)**
antisymmetric state at coincidence — is robust; the precise numbers carry the
#186 soliton-profile ~3% uncertainty. Weak-field / semi-dynamical soliton.

## Reproduce

```bash
python -m experiments.closure_ledger.two_throat_hartree_fock_probe
# Verdict: TWO_THROAT_HARTREE_FOCK_OVERLAP_NORMALIZED_DIRECT_PLUS_EXCHANGE_FERMION_LOWER_FOR_REPULSIVE_V_ANTISYM_FORBIDDEN_AT_COINCIDENCE
```
