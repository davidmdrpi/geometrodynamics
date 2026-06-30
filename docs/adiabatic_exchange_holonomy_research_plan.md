# Adiabatic two-throat exchange holonomy: measure the Pin⁻ sign along a swap path (PR #188)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Measuring the exchange sign, not asserting it

PR #185 derived the two-throat exchange sign `−1` from the abstract identity
`T² = −I` (the Pin⁻ monodromy; the swap `≃` a 2π rotation by geon statistics).
This probe makes that **operational**: it transports the throat's spin-½ state
**adiabatically** along an explicit two-throat **swap path** and **measures**
the accumulated holonomy — which comes out exactly `−1` (a π Berry phase).

## The swap path and the spin-statistics connection

The relative-coordinate space of two identical throats is
`(ℝ³∖0)/ℤ₂ ≃ RP² × ℝ₊`, with the angular factor `S²/antipodal = RP²` — the
**same antipodal closure** that makes the throat (#169/#170). The exchange
moves the relative coordinate `r → −r`, a path from a point to its antipode on
`S²` that is **closed in RP²**: the generator of `π₁(RP²) = ℤ₂`. By the
**Finkelstein–Rubinstein / Friedman–Sorkin** spin-statistics theorem for
solitons/geons, this swap loop is **homotopic to a 2π rotation** of one throat
(the belt-trick / tether twist). So the exchange holonomy = the spin-½
holonomy of a 2π rotation.

## The holonomy measured

Path-ordering the spin connection along the swap (2π) loop,

```
dU/ds = −i (ω·σ / 2) U ,   [ω]× = Ṙ Rᵀ   (the loop's angular velocity),
```

gives the adiabatic holonomy

```
Hol = −I   (‖Hol + I‖ ~ 10⁻⁶) ,
```

so the throat's spin-½ state returns to **minus itself**:

| quantity | value |
|---|---|
| measured exchange sign `⟨ψ\|Hol\|ψ⟩` | **−1** |
| Berry phase (arg of the holonomy eigenvalue) | **π** |

The exchange of two throats is a π Berry phase — a *measured* `−1`, not an
assumed one.

## Topological (path-independent)

The holonomy is the `ℤ₂` homotopy class of the loop, not the path:

| loop | holonomy | measured sign |
|---|---|---:|
| swap (2π, the exchange) | `−I` | **−1** |
| double-swap (4π, two exchanges) | `+I` | +1 |
| contractible (no exchange) | `+I` | +1 |

A swap path with a **wandering** rotation axis (a non-commuting, genuinely
path-ordered transport) gives the same `−I`, converging as the transport is
refined (`N = 250 → 4000` steps). Any way of carrying out the exchange gives
the same `−1`. The controls confirm the `−1` is the **single-swap (odd)**
class: a double-swap (two exchanges, 4π) gives `+1` (two fermion exchanges =
a boson) and a contractible loop gives `+1` (no exchange, no sign).

## The Pin⁻ identification

The measured `−1` **is** the Pin⁻ monodromy `T = iσ_y`, `T² = −I`
(`½ tr T² = −1`; #170/#174/#183) — the Pin⁻ structure on the non-orientable
`RP²` closure (the unique spin structure `RP²` admits), which makes the throat
a **spinor** (spin-½). The 2π/swap holonomy of a spin-`j` object is
`(−1)^{2j}`:

| object | `(−1)^{2j}` |
|---|---:|
| scalar `j=0` | +1 |
| **spinor `j=1/2` (the throat)** | **−1** |
| vector `j=1` | +1 |

So the spin-½ throat gives `−1` (fermion); a scalar/bosonic throat would give
`+1` along the **same** swap path. The adiabatic holonomy `−I` is exactly
`T² = −I`, now transported along the swap rather than read off the algebra.

## Honest scope

This **operationalizes** the Finkelstein–Rubinstein / geon-statistics result:
the holonomy is computed exactly (a path-ordered SU(2) transport) and is
**topological** (the `ℤ₂` class), so the measured exchange sign `−1` is exact.
The **swap path** is the reduced relative-coordinate / frame model — the loop
in the two-throat configuration space (`RP²`) lifted to the spin frame
bundle — and the spin-statistics **connection** (exchange `≃` 2π rotation) is
the FR theorem, cited rather than re-derived; the throat's spinor (Pin⁻)
nature is the #170 result (`RP²` admits only Pin⁻). The **adiabatic** limit
(slow transport) is assumed. No #180 soliton dynamics is invoked — this is the
statistics / holonomy layer, complementary to the #185–#187 spatial exchange
kernel and Hartree–Fock energies.

## Reproduce

```bash
python -m experiments.closure_ledger.adiabatic_exchange_holonomy_probe
# Verdict: ADIABATIC_TWO_THROAT_EXCHANGE_HOLONOMY_MEASURES_THE_PIN_MINUS_SIGN_MINUS_ONE_ALONG_THE_SWAP_PATH
```
