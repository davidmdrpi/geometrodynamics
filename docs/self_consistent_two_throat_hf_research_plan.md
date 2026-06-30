# Self-consistent two-throat Hartree–Fock relaxation (PR #189)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Relaxing the rigid orbitals of #187

PR #187 assembled the two-throat Hartree–Fock energy with **rigid** orbitals
(two fixed #180 throat-solitons) and flagged the self-consistent solve as a
follow-up. This probe does that follow-up: it **relaxes** the two orbitals in
each other's **direct (Hartree) + exchange (Fock)** field by a genuine
self-consistent-field (SCF) iteration, so the throats deform in the mean field
and the energy drops to its variational minimum.

## The Fock operator

Each orbital is relaxed in

```
F = h + V_H − K ,
  h   = −½∇² + V_ext       (kinetic + the confining well of each throat),
  V_H = ∫ ρ(x') V(x−x')    (the DIRECT / Hartree field of the other throat),
  K   = the non-local Fock EXCHANGE operator,
        (Kφ)(x) = Σ_j φ_j(x) ∫ φ_j(x') V(x−x') φ(x') .
```

Two same-spin fermions (the two throats, whose spatial state is antisymmetric
— the Pin⁻ sector of #185/#188) occupy the two lowest orbitals; the SCF
relaxes them (orthonormal, by imaginary-time gradient descent) to the
self-consistent ground state.

## Convergence

The imaginary-time HF relaxation lowers the energy **monotonically**:

```
E : −3.808 → −3.905 ,   final ΔE ~ 10⁻⁸ → 0 .
```

It settles to a self-consistent fixed point — the orbitals come to rest in the
mean field they themselves produce (the defining HF condition). The monotone
descent (no oscillation) confirms a genuine variational relaxation; a naive
diagonalization-SCF *oscillates* for these near-degenerate bonding/antibonding
orbitals (eigenstate-swapping), which is why the imaginary-time gradient
descent is used.

## The relaxation (energy lowering)

| | energy |
|---|---:|
| rigid (unrelaxed orbitals, full HF energy) | −3.808 |
| relaxed (self-consistent) | **−3.905** |
| lowering | **2.54 %** |

The rigid orbitals — evaluated with the full HF energy including the
interaction — give `E_rigid`; relaxing them self-consistently lowers the
energy by 2.54%. In #187 the orbitals were held rigid; here they deform in the
direct + exchange field to lower the energy (the variational principle the SCF
realizes).

## The orbitals deform

The two-throat density polarizes/spreads in the mean field:

| quantity | rigid | relaxed |
|---|---:|---:|
| density RMS width | 2.646 | 3.055 |
| density fidelity (rigid vs relaxed) | — | 0.978 |

The throats respond to each other's field — a real deformation of the
orbitals, the content the rigid #187 sandbox omitted.

## The exchange field matters

| SCF | energy |
|---|---:|
| full Hartree–Fock (direct − exchange) | −3.905 |
| Hartree-only (direct, no exchange) | −3.288 |
| **exchange lowering** | **0.616** |

For the two **same-spin** throats the non-local Fock exchange **substantially
lowers** the energy: the exchange hole (#186/#187) keeps the like throats
apart, reducing the repulsive direct energy. The `−1` of #185/#188 is doing
real work in the self-consistent mean field, not only in the rigid kernel.

## Honest scope

A **sandbox SCF**, in 1D: the confinement is an external double well (a
stand-in for the throats' self-binding — the #180 self-gravity); the
interaction is a screened-photon (Yukawa) stand-in for the BAM throat-fibre
exchange; the Hartree–Fock is spatial-orbital (same-spin / unrestricted), in
one dimension for tractability. The **SCF itself is genuine** — an
imaginary-time relaxation to self-consistency, monotone and machine-converged
— and the qualitative physics (convergence, the variational energy lowering,
the orbital deformation, the exchange lowering) is robust and standard. The
full 3D self-gravitating two-throat SCF — relaxing actual #180 ψ–Φ–q
throat-solitons in each other's direct + exchange field — is the follow-up.
Weak-field, code units.

## Reproduce

```bash
python -m experiments.closure_ledger.self_consistent_two_throat_hf_probe
# Verdict: SELF_CONSISTENT_TWO_THROAT_HF_RELAXATION_CONVERGES_LOWERS_ENERGY_ORBITALS_DEFORM_IN_DIRECT_PLUS_EXCHANGE_FIELD
```
