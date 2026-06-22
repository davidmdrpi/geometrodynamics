# Pin⁻ on the throat's RP² mouth: the exchange sign and the Fermi equation of state (PR #170)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The deferred calculation, done

PR #169 established the throat mouth is the non-orientable `RP²` and noted —
as a *remark* — that it admits a Pin structure carrying the spin-½
character. This probe stops deferring: it takes the Pin⁻ structure and
shows it **delivers** the two things that make it matter — the −1 exchange
sign and the Fermi equation of state.

## The chain

```
RP² mouth  →  Pin⁻ structure (the only one RP² admits)
           →  spin-½ spinor: 2π rotation acts as −1
           →  exchange of two throats ≃ 2π rotation  (Finkelstein–Rubinstein)
           →  antisymmetric two-throat wavefunction (exchange sign −1)
           →  Pauli exclusion: occupation n_p ∈ {0, 1}
           →  the Fermi equation of state (degeneracy pressure):
                P = ⅔u, P ∝ n^{5/3}  (non-relativistic, Γ = 5/3)
                P = ⅓u, P ∝ n^{4/3}  (ultra-relativistic, Γ = 4/3)
```

## RP² carries Pin⁻ — and only Pin⁻

The total Stiefel–Whitney class is `w(RP²) = (1+a)³ = 1 + a + a²`, so
`w₁ = a` (non-orientable → no Spin) and `w₂ = a²`. The admissibility
conditions then give:

| structure | condition | RP²? |
|---|---|---|
| Spin | `w₁ = w₂ = 0` | **no** |
| Pin⁺ | `w₂ = 0` | **no** (`w₂ = a² ≠ 0`) |
| Pin⁻ | `w₂ + w₁² = 0` | **yes** (`a² + a² = 0`) |

So the throat mouth has a **unique, definite** spinor structure — Pin⁻, the
non-orientable analogue of Spin.

## The exchange sign = −1

The Pin⁻ spinor is spin-½: `R(2π) = exp(−iπσ_z) = −I` (and only `R(4π) = +I`).
By the **Finkelstein–Rubinstein** construction, exchanging two identical
throats is homotopic to a 2π rotation of one (the orientation-entanglement /
belt trick; `π₁` of the two-particle configuration space in ≥3D is `ℤ₂`,
with the exchange generator identified with the 2π-rotation generator). The
Pin⁻ spinor assigns that loop the value −1, so the two-throat wavefunction
is **antisymmetric**. The spin-statistics connection is *realised*, not
assumed: the same spin-½ holonomy gives both `2π = −1` and `exchange = −1`.

## The Fermi equation of state

Antisymmetry ⟹ Pauli exclusion (`n_p ∈ {0,1}`) ⟹ filling the Fermi sphere:

| regime | P/u | polytropic Γ |
|---|---:|---:|
| non-relativistic (`ε = p²/2m`) | 2/3 | **5/3** |
| ultra-relativistic (`ε = pc`) | 1/3 | **4/3** |
| Bose gas at `T=0` (contrast) | — | pressure **0** |

The T=0 **degeneracy pressure is strictly positive** — the pressure that
holds up white dwarfs and neutron stars — and it exists *only* because of
the exclusion: a Bose gas at T=0 collapses all quanta to `p=0` with zero
pressure. The Fermi equation of state is delivered by the −1 exchange sign
of the Pin⁻ mouth.

## Honest scope

- **Computed** here: the Pin⁻ classification (Stiefel–Whitney classes), the
  spinor `2π = −1` sign, and the Fermi-gas EoS (the `P/u` ratios and
  polytropic indices from the T=0 momentum integrals).
- **Cited**, not re-derived: the Finkelstein–Rubinstein homotopy that
  identifies the two-particle exchange loop with a 2π rotation — the one
  configuration-space theorem linking the throat's internal Pin holonomy to
  the physical exchange. The spinor sign and the EoS are not citations.

This is the same orientability grading already in BAM as the C-swap
(`C = iσ_y`, `T² = −1`; PR #63) and the even-`k` absence (PR #67), now
carried through to the statistics and the equation of state.

## Reproduce

```bash
python -m experiments.closure_ledger.pin_rp2_fermi_statistics_probe
# Verdict: PIN_MINUS_ON_RP2_DELIVERS_FERMI_EXCHANGE_SIGN_AND_EOS
```
