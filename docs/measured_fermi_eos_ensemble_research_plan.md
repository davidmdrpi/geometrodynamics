# The measured Fermi equation of state: a many-throat ensemble (PR #172)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Why this exists (the second closing option, same branch as #171)

PR #170 **assumed** antisymmetry and read the index 5/3 off the analytic
Fermi integral; PR #171 derived the −1 exchange sign **topologically** (the
Pin⁻ geon). This probe is the third leg: **simulate** a many-throat
ensemble and **measure** the equation of state, rather than assuming the
occupation and reading off 5/3. It lives on the **same branch as #171** so
the three routes — assumed-analytic (#170), topological-sign (#171),
measured (here) — can be compared directly.

## The simulation

`N` identical throats are free fermions in a cubic box of side `L`. The −1
exchange sign (the Pin⁻ result of #170/#171) is realised concretely as
**Pauli single-occupancy**: each single-particle box mode `(n_x, n_y, n_z)`,
`n_i = 1, 2, 3, …`, holds at most `g = 2` fermions. The many-body ground
state is the **filled Fermi sea** (a Slater determinant) — built by filling
the `N` lowest of ~1.3 million enumerated modes. No occupation distribution
is assumed; the ground state is *constructed* by level-filling.

The equation of state is then **measured** from the simulated ground-state
energies:

- **P/u** from the box volume derivative `P = −dE/dV` (the virial
  relation): `E = K/L²` (NR) or `K/L` (UR) with `K` independent of `L`, so
  `P = (q/3)·E/V`.
- **Γ = d ln P / d ln n** from the local log-slope of the filled-mode energy
  sum `K(N)`, finite-size-extrapolated to the thermodynamic limit via the
  Weyl surface correction `Γ(N) = Γ∞ − a·N^{−1/3}`.

## The result (measured, not assumed)

| regime | P/u (measured) | Γ (measured) |
|---|---:|---:|
| non-relativistic (`ε ∝ p²`) | 0.6667 = 2/3 | **1.6665 ≈ 5/3** (0.01%) |
| ultra-relativistic (`ε ∝ p`) | 0.3333 = 1/3 | **1.3332 ≈ 4/3** (0.01%) |
| Bose control (all `N` in ground mode) | — | **1** (no degeneracy) |

The local slopes converge monotonically toward the targets and the
`N^{−1/3}` extrapolation lands on them. The **Bose control** (Γ = 1, with
the `T=0` degeneracy pressure vanishing as the mode energy `∝ 1/L² → 0`)
shows the stiffening is a measured **consequence of the −1 exchange sign**,
not a universal property of the box.

## The three routes agree

| regime | assumed (#170) | topological (#171) | measured (#172) |
|---|---|---|---:|
| non-relativistic | 5/3 | −1 ⇒ Fermi ⇒ 5/3 | 1.6665 |
| ultra-relativistic | 4/3 | −1 ⇒ Fermi ⇒ 4/3 | 1.3332 |

The measured indices reproduce #170's assumed-analytic values and confirm
the equation of state implied by #171's topological −1 sign.

## Honest scope

- **Measured**: the box spectrum (enumerated), the Pauli-filled
  ground-state energy `K(N)`, the `P/u` virial ratio, the polytropic index
  Γ (finite-size-extrapolated), and the Bose control.
- **Input, not re-derived**: the −1 exchange sign itself — realised as Pauli
  single-occupancy — which is the Pin⁻/geon result of #170/#171. This probe
  measures the *equation of state that the −1 produces*; it does not
  re-establish the sign.
- **Idealizations** (the standard degenerate-gas ones): free
  (non-interacting) throats, `T = 0`, and the cubic box as the confining
  volume. Interactions and finite temperature add the usual corrections
  without changing the leading degeneracy index.

## Reproduce

```bash
python -m experiments.closure_ledger.measured_fermi_eos_ensemble_probe
# Verdict: MEASURED_FERMI_EOS_FROM_PAULI_FILLING_GAMMA_5_3_AND_4_3
```
