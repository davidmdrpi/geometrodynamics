# Dynamic two-throat exchange path with back-reaction (PR #191)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question

PR #188 measured the exchange sign as the **adiabatic** holonomy of the swap
loop — the Pin⁻ monodromy `T² = −I`, the geometric `−1`, transported
infinitely slowly. PR #189 relaxed two throats in each other's static field.
This probe asks the dynamical question those two leave open:

**What happens to the exchange `−1` when the swap is performed at finite
speed, with a field that back-reacts on the moving throat?**

Is the `−1` a robust dynamical output, or only an adiabatic idealization? And
what is the *cost* of a fast swap?

## The model (an effective dynamical model)

The throat's internal Pin spinor `ψ` is driven around the swap loop
`n̂(s) = (cos 2πs, sin 2πs, 0)`, `s = t/T`, over a finite duration `T` (so the
swap speed is `1/T`), under

```
H(s, x) = (Δ/2) n̂(s)·σ  +  g·x·σ_z ,
```

the second term coupling the spinor to a **back-reacting field** coordinate
`x(t)` obeying

```
ẍ + γẋ + ω²x = −κ ⟨σ_z⟩ ,
```

so the moving throat **sources** the field (`⟨σ_z⟩` drives `x`) and the field
**acts back** on the spinor (`g·x·σ_z` in `H`) — a genuine two-way coupling,
not a one-way drive. The spinor is integrated with RK4 and the field with a
velocity-Verlet/Euler step. The geometric phase is extracted by continuously
factoring out the dynamical phase
`φ_dyn(t) = ∫₀ᵗ ⟨ψ|H|ψ⟩ dt'`:

```
φ_geo = arg( ⟨ψ(0)|ψ(T)⟩ · e^{+i φ_dyn} ) ,
```

robust because `e^{iφ_dyn}` is periodic (factoring it out cannot lose
windings), and `P_exc = 1 − |⟨ψ_lower(H_f)|ψ(T)⟩|²` measures the
non-adiabatic leakage out of the instantaneous lower state.

Params: `Δ = 1.0`, `ω = 2.0`, `κ = 0.6`, `γ = 0.3`, `dt = 0.02`.

## What the probe establishes

- **The invariant.** The swap loop's exact adiabatic Berry phase is
  `−½ · 2π = −π` — the exchange sign `−1` (the #188 holonomy `e^{iπ} = −1`),
  a geometric invariant independent of speed and back-reaction.
- **Adiabatic recovery.** As the swap slows (`T → ∞`), the *dynamical*
  geometric phase `φ_geo → −π` (deviation `≈ 0.02` at `T = 1000`) and the
  non-adiabatic excitation `P_exc → 0`. The full dynamics reproduce the
  exchange `−1` in the adiabatic limit.
- **The non-adiabatic cost.** At finite speed `φ_geo` deviates from `−π` by
  `O(1/T)` — `deviation × T ≈ 29.7` (constant across `T = 10…300`) — and
  `P_exc` grows with speed. This is the quantified price of a fast swap: the
  throat cannot follow it, and the geometric phase is not yet the clean `−π`.
- **Back-reaction.** The moving throat sources the field — peak field energy
  `≈ 0.028` at `T = 20`, falling to `≈ 1.5e-4` when the swap is slow. Energy
  is exchanged with the field and it acts back on the spinor, but the
  adiabatic limit (the `−1`) is **unchanged**. The dynamics add a quantified,
  adiabatically-vanishing cost; they do **not** alter the exchange sign.

## Honest scope

This is an **effective** dynamical model: the throat's internal Pin spinor on
the swap loop, plus one back-reacting field mode. It is **not** a
field-resolved real-time solve of two physical throats moving through the BAM
photon field of #190 — that (the orbitals translating, the Coulomb field
re-solved every step, the geon swap geometry resolved in 3D) is the follow-up.
What this probe shows is the *mechanism*: the exchange `−1` is the adiabatic
limit of a real dynamical process, recovered by slowing the swap, with
non-adiabatic and back-reaction corrections that are explicit, `O(1/T)`, and
vanish adiabatically. Code units; weak-field.

## Reproduce

```bash
python -m experiments.closure_ledger.dynamic_two_throat_exchange_probe
# Verdict: DYNAMIC_TWO_THROAT_EXCHANGE_RECOVERS_THE_ADIABATIC_MINUS_ONE
#          _WITH_NON_ADIABATIC_AND_BACK_REACTION_CORRECTIONS_VANISHING_ADIABATICALLY
```
