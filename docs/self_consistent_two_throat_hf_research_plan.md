# Self-consistent two-throat Hartree–Fock relaxation (PR #189)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Relaxing the rigid orbitals of #187

PR #187 assembled the two-throat Hartree–Fock energy with **rigid** orbitals
(two fixed #180 throat-solitons) and flagged the self-consistent solve as a
follow-up. This probe does that follow-up: it **relaxes** the two orbitals in
each other's **direct (Hartree) + exchange (Fock)** field by a genuine
self-consistent-field (SCF) iteration, so the throats deform in the mean field
and the energy drops to its variational fixed point.

## The Fock operator (orbital-specific, self-interaction-free)

Each orbital is relaxed in

```
F_i = h + J_{≠i} − K_{≠i} ,
  h      = −½∇² + V_ext            (kinetic + the confining well),
  J_{≠i} = ∫ |φ_{≠i}|² V(x−x')     (the DIRECT / Hartree field of the OTHER
                                    throat — self-interaction excluded),
  K_{≠i} = the non-local Fock EXCHANGE with the other orbital.
```

Excluding self-interaction makes `F_i` the **exact variational derivative** of
the reported energy `E = Σ_i ⟨i|h|i⟩ + (J₀₁ − K₀₁)` — the same functional for
the operator and the energy, in **both** the full-HF and the Hartree-only
control. Two same-spin fermions (the two throats, whose spatial state is
antisymmetric — the Pin⁻ sector of #185/#188) occupy two orthonormal orbitals,
relaxed by imaginary-time gradient descent.

## Convergence and robustness

The imaginary-time HF relaxation lowers the energy **monotonically**:

```
E : −3.808 → −3.905 ,   final ΔE ~ 10⁻⁸ → 0 .
```

It settles to a self-consistent fixed point — the orbitals come to rest in the
mean field they produce. The fixed point is **robust across seeded restarts**:
five random localized initial orbitals converge to

```
[−3.897, −3.905, −3.905, −3.905, −3.905] ,   spread ~ 8×10⁻³ ,
```

so the relaxed state is a self-consistent **variational fixed point** (robustly
reached, **not certified the global ground state**). The monotone descent (no
oscillation) confirms a genuine variational relaxation; a naive
diagonalization-SCF *oscillates* for these near-degenerate bonding/antibonding
orbitals (eigenstate-swapping), which is why the imaginary-time gradient
descent is used.

## The relaxation (energy lowering)

| | energy |
|---|---:|
| rigid **1D #187-style reference** (unrelaxed orbitals, full HF energy) | −3.808 |
| relaxed (self-consistent) | **−3.905** |
| lowering | **2.54 %** |

The rigid 1D reference — the unrelaxed orbitals evaluated with the full HF
energy, the **1D analogue** of #187's rigid-orbital evaluation (not the 3D
number) — is lowered by 2.54% when the orbitals relax self-consistently. In
#187 the orbitals were held rigid; here they deform in the direct + exchange
field to lower the energy.

## The orbitals deform

| quantity | rigid | relaxed |
|---|---:|---:|
| density RMS width | 2.646 | 3.055 |
| density fidelity (rigid vs relaxed) | — | 0.978 |

The two-throat density polarizes/spreads in the mean field — a real
deformation of the orbitals, the content the rigid #187 sandbox omitted.

## The exchange field matters (consistent control)

The Hartree-only control uses the **same self-interaction-free variational
functional**, just with the exchange dropped from **both** the Fock operator
(`F_i = h + J_{≠i}`) and the energy (`E = Σ⟨i|h|i⟩ + J₀₁`) — so its operator
and energy are the same functional.

| SCF | energy |
|---|---:|
| full Hartree–Fock (direct − exchange) | −3.905 |
| Hartree-only (direct, consistent control) | −3.338 |
| **exchange lowering** | **0.567** |

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
imaginary-time relaxation to self-consistency, monotone and machine-converged,
with an orbital-specific **self-interaction-free** Fock operator consistent
with the reported energy — and the qualitative physics (convergence, the
variational energy lowering, the orbital deformation, the exchange lowering)
is robust and standard. The relaxed state is a **self-consistent variational
fixed point** (robust across seeded restarts, not certified the global ground
state), and the "rigid" reference is the **1D #187-style** evaluation (not the
3D number). The full 3D self-gravitating two-throat SCF — relaxing actual #180
ψ–Φ–q throat-solitons in each other's direct + exchange field — is the
follow-up. Weak-field, code units.

## Reproduce

```bash
python -m experiments.closure_ledger.self_consistent_two_throat_hf_probe
# Verdict: SELF_CONSISTENT_TWO_THROAT_HF_RELAXATION_CONVERGES_LOWERS_ENERGY_ORBITALS_DEFORM_IN_DIRECT_PLUS_EXCHANGE_FIELD
```
