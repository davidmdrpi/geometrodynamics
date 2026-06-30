# Rigid soliton exchange-kernel hardening (PR #186)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## From a promising kernel to a trustworthy one

PR #185 derived the two-throat exchange kernel `K_exchange(R) = (−1)·K(R)`
from the rigid overlap of two #180 self-gravitating throat-solitons. This
probe hardens it (as #177 hardened #176) with normalization, convergence, and
direct-term controls — the benchmark the #187 Hartree–Fock sandbox needs.

## Normalization

| check | value |
|---|---|
| orbital norm `∫|φ|² d³r` | `1.000000` |
| self-overlap `K(0)` | `1.001` (reproduces the norm to 0.1%) |
| parity `K(2)` vs `K(−2)` | `0.40963` = `0.40963` (exact) |
| Cauchy–Schwarz `K(R) ≤ K(0)` | `K(2) = 0.41 < 1` ✓ |

The exchange amplitude is a proper normalized overlap; the `K(0) = 1.001`
residual is the overlap-quadrature discretization (T3).

## Convergence

**Overlap-integral grid** — refining the angular/radial quadrature `(Nr, Nu)`:

| grid | 160×80 | 320×160 | 640×320 | 1000×500 |
|---|---:|---:|---:|---:|
| `K(2)` | 0.40977 | 0.40963 | 0.40960 | 0.40959 |

The production grid (320×160) differs from the finest by **0.01%** — the
overlap integral is converged to < 0.1%.

**Soliton grid** — rebuilding the #180 throat-soliton at `N = 240 → 320`:

| N | 240 | 320 |
|---|---:|---:|
| `K(2)` | 0.4096 | 0.4215 |

A **~2.9%** shift — *not* the overlap integral but the soliton profile's core
grid-sensitivity (the documented #180 caveat: the sharp core carries ~10%
per-point grid uncertainty, integral overlaps a few %). The kernel's shape and
scale are trustworthy to a few %, with the soliton profile the limiting
factor — an inherited, known limitation, not a flaw in the kernel computation.

## Direct vs exchange — two distinct GR-geometric kernels

The DIRECT density-overlap `D(R) = ∫ ρ_a ρ_b d³x` (the Hartree channel,
sign-independent) and the EXCHANGE amplitude-overlap `K(R)` (the ±-carrying
channel) are genuinely different objects:

| R | direct `D(R)` | exchange `X = K²` |
|---:|---:|---:|
| 0.0 | 0.065 | 1.00 |
| 1.0 | 0.038 | 0.63 |
| 2.0 | 0.008 | 0.17 |
| 3.0 | 0.001 | 0.023 |
| 4.0 | 0.000 | 0.002 |
| 6.0 | 0.000 | 0.000 |

Both decay monotonically to zero at large R, but the **direct decays faster**
(it is the product of densities; the exchange is the amplitude overlap, which
carries the statistics).

## Direct-term controls

1. **Far separation** (`R = 6`): both channels vanish (`D = 0`, `X = 0`) — widely
   separated throats are distinguishable (no overlap, no exchange, no
   interaction).
2. **Sign-independence**: `D(R)` is a positive integral with no sign — identical
   for the boson (+) and fermion (−) sectors — while the exchange carries the
   Pin⁻ `−1`. So the `−1` the geometry selects lives **purely in the exchange
   channel**; the direct is the control that isolates it.
3. **The energy split**: consequently the two-body energy is
   `E = E_direct ∓ E_exchange` — the direct the common part, the exchange the
   sign-dependent splitting. The #187 Hartree–Fock sandbox evaluates this
   against an interaction `V`.

## Honest scope

This hardens the **rigid** soliton-overlap kernel (two rigid copies of the
#180 radial throat-soliton at separation `R`). The overlap integral is
benchmark-converged (< 0.1%); the residual ~3% is inherited from the soliton
profile's core grid-sensitivity (#180), the honest dominant uncertainty. The
direct and exchange OVERLAP kernels are separated and controlled here;
convolving them with an interaction `V` to get the Hartree (direct) and
exchange ENERGIES — the two-throat Hartree–Fock sandbox — is **PR #187**.
Relaxing the orbitals self-consistently in each other's presence (beyond the
rigid approximation) is a further follow-up. Weak-field / semi-dynamical
soliton.

## Reproduce

```bash
python -m experiments.closure_ledger.rigid_soliton_exchange_kernel_hardening_probe
# Verdict: RIGID_SOLITON_EXCHANGE_KERNEL_HARDENED_NORMALIZED_OVERLAP_CONVERGED_DIRECT_TERM_CONTROLLED
```
