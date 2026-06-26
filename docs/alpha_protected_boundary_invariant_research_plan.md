# α as a protected boundary invariant, not a continuous tuning parameter (PR #184)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Reframing the 137 problem

PR #105/#143 classified the EM coupling α: BAM **derives** the structure
around it — the charge **quantum** (`|c₁| = 1`, the integer Hopf number;
charge quantization is topological), the `1/2π` closure **loop measure** (the
`2π` in the Schwinger anomaly `a = α/2π`), and the **running** — but **not**
the **value** `α⁻¹ ≈ 137` (the residual "137 problem"; the fit-independent
scans against the closure numbers failed).

The roadmap has a sharper lens. #181/#182 showed the order-field **winding**
is a topological charge that survives smooth evolution and changes only at
`|q| = 0`; #183 showed the odd-k **generation** sector survives smooth bulk
deformation and changes only at a topology change. This probe applies the
**same test** to α: is the derived structure a **protected boundary
invariant** — topologically robust to smooth deformation — or just another
continuous knob?

## The charge quantum is a boundary Chern number

The BAM charge quantum is the first Chern number of the Hopf / spin-½ monopole
over the boundary S² — the Gauss-law charge `c₁ = (1/2π)∮_{S²} F`. Computed by
the exactly-quantized **Fukui–Hatsugai–Suzuki** lattice method on the lower
band of `H = d·σ` (`d = (sinθ cosφ, sinθ sinφ, cosθ + m)`, at `m = 0`):

```
c₁ = +1     (|c₁| = 1, the integer Hopf number; sign is orientation)
```

An exact integer — charge quantization is a boundary integral, not an input.

## Protected vs continuous tuning (the discriminator)

Under **30 smooth diffeomorphisms** of the boundary geometry
(`θ → θ + a sinθ + b sin2θ`, monotone):

| quantity | response |
|---|---|
| charge quantum `c₁` (Chern number) | **invariant** — same integer, max move `5×10⁻⁷` |
| a continuous coupling functional `⟨A_φ⟩` | **drifts** — `6.8%` mean (`15.8%` max) |

A tuning parameter responds continuously to the geometry; a protected boundary
invariant does not. The charge quantum is on the **protected** side — so α's
derived structure should be **tested as** a protected boundary invariant, not
fit as a continuous knob (the failure mode of the 137-problem scans).

## The loop measure, and changes only at a topology change

The total boundary flux is quantized in units of the closure quantum `2π`:

```
∮ F = 2π · c₁ = 2π·(integer) ,
```

so the `2π` in the Schwinger anomaly `a = α/2π` is the protected loop measure,
leaving only the α prefactor as the residual.

The charge integer changes **only** when the Berry **gap closes**. Sweeping
the gap parameter `m` (which moves the degeneracy `d = 0` relative to the
boundary):

| m | 0.0 | 0.5 | 0.9 | 1.0 | 1.1 | 1.5 | 2.0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `C(m)` | 1 | 1 | 1 | **0** | 0 | 0 | 0 |
| gap `min|d|` | 1.0 | 0.5 | 0.1 | **0.001** | 0.1 | 0.5 | 1.0 |

`C` stays `1` while the gap is open and jumps to `0` exactly at `m = 1`, where
`min|d| → 0` — the degeneracy crosses the boundary. This is the EM-boundary
analog of the order-field `|q| = 0` (#182) and the spin-closure `½ tr T² = 0`
(#183): the protected invariant moves only through the singular /
topology-change event, never by smooth deformation.

## Unity with #181/#182/#183

The EM charge quantum is to the **boundary** what the order-field winding is to
the **soliton** (#181/#182) and the generation sector is to the **bulk**
(#183): a protected topological charge robust to smooth deformation, changing
only at a topology-change event (`|q| = 0` / `½ tr T² = 0` / the Berry gap
closing). The program's protected-invariant picture now covers the EM sector.

## Honest scope

This does **not** derive the value `α⁻¹ ≈ 137` — that residual input (the 137
problem) stands, exactly as #105/#143 found. What is established is that α's
**quantized structure** (the charge quantum + the `1/2π` loop measure) is a
**protected boundary invariant** (topological), so α decomposes as
protected-boundary-structure × one residual scale, and should be tested that
way rather than fit as a continuous tuning family. The Chern number here is the
spin-½/Hopf monopole boundary invariant (the BAM charge quantum); the
deformations are diffeomorphisms of the boundary 2-sphere. This refines the
#105/#143 ledger — the derived EM structure is specifically protected.

## Reproduce

```bash
python -m experiments.closure_ledger.alpha_protected_boundary_invariant_probe
# Verdict: ALPHA_CHARGE_QUANTUM_AND_LOOP_MEASURE_ARE_PROTECTED_BOUNDARY_INVARIANTS_NOT_TUNING_THE_VALUE_REMAINS_RESIDUAL
```
