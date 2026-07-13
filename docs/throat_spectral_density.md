# The throat's spectral density (PR #215)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #214 promoted the absorber to a
> degree of freedom but modelled its spectrum as a flat oscillator bank
> — a shape chosen by hand, and the one thing about the absorber still
> not read off the geometry. This PR retires the bank: the spectral
> density is the throat's own greybody transmission, with both wings
> pinned by theorems. The companion probe machine-checks every claim
> (~15 s).

## 0. The named successor

After #213 (the propagator from complete histories) and #214 (the
absorber as dynamics), one modeling choice remained: the absorber's
spectral shape J(ν). The physical absorber is the throat itself — a 5D
Tangherlini mouth (f(r) = 1 − r_h²/r², the #168 Killing-horizon
identification) whose horizon is a classical one-way membrane with a
*computable* frequency response. Deriving that response and replacing
the flat bank removes the last hand-chosen element of the absorber
sector.

## 1. The scattering problem

A massless scalar mode φ(r)Y_ℓ(S³)e^{−iωt} on the Tangherlini exterior
reduces, with ψ = r^{3/2}φ and the tortoise coordinate dx = dr/f, to

```
ψ_xx + [ω² − V]ψ = 0 ,     V(r) = f·[ (ℓ(ℓ+2) + 3/4)/r² + (9/4)r_h²/r⁴ ] ,
x(r)  = r + (r_h/2)·ln((r − r_h)/(r + r_h))          (analytic).
```

Ingoing data at the horizon, Hankel least-squares matching at large r
(the far radius scales with ω to keep the metric's long-range phase
correction below the matching tolerance) gives the greybody
transmission T_ℓ(ω). Validation: flux conservation |R + T − 1| ≤ 2×10⁻⁵
across the working band (≤ 1.2×10⁻³ at the highest frequencies), and
horizon-regulator independence to 4×10⁻⁶.

## 2. The IR wing is the area theorem

The universal low-frequency absorption theorem (valid in any dimension)
states σ_s-wave(ω→0) = A_horizon. For 5D Tangherlini A_h = 2π²r_h³.
Measured: σ_s/A_h = **1.012** at ωr_h = 0.04, decreasing monotonically
toward 1, extrapolating to **0.983**; and T₀ ∝ ω³ with log-log slope
**3.02**. So the spectral density's IR wing is

```
J(ω → 0)  ∝  A_h ω³ / 4π
```

— *the horizon area*, not a model parameter. The universe's low tower
modes are nearly free: **IR transparency** (everyday fields propagate;
the throat is invisible to wavelengths that don't resolve it).

## 3. The UV wing is the photon sphere

Null geodesics on Tangherlini: the impact-parameter function
b(r)² = r²/f maximizes at the photon sphere **r_c = √2 r_h** with
critical impact parameter **b_c = 2 r_h** (both machine-checked to
10⁻⁶). The transmission's half-crossings sit at the eikonal frequencies
ω(T = ½) ≈ (ℓ+1)/b_c — measured ratios 1.04 / 1.02 / 1.01 for
ℓ = 1, 2, 3, with the WKB ratio increasing toward 1 as WKB demands —
and above the barrier T → 1 (0.9999 at ωr_h = 3): **UV-black**. The
throat is a *perfect* absorber for every mode that resolves it.

## 4. The horizon is the continuum

#214's finite bank revived at the Poincaré time, and ε → 0 required the
continuum limit N → ∞ *taken first*. On the throat this is not a limit
anyone takes — the geometry has taken it: the tortoise depth of a
regulated interior grows as

```
L(δ) = const − (r_h/2)·ln δ        (fitted coefficient −0.5005 r_h)
```

— divergent at the horizon, so the level spacing π/L → 0 and the
recurrence time 2L → ∞. Classically the throat never revives. The
infinite tortoise throat **is** the N → ∞ bank.

## 5. The weld: a parameter-free transit law

Terminate a cavity (the universe stand-in) with the throat and find its
complex quasimodes ω_q − iγ/2 (ingoing at the horizon, node at the
wall, secant iteration in complex ω). The prediction, with **no free
parameter**:

```
γ = T(ω_q) / (2 L_cav) ,      L_cav = x_wall − x_barrier-peak
```

(the sub-barrier region is the drain, not the cavity — the wave
reflects at the barrier's outer face and T is the tunneling probability
into the horizon per round trip). Measured down the ladder:

| ω_q | γ measured | T/(2L_cav) predicted | ratio |
|---:|---:|---:|---:|
| 0.23188 | 4.117×10⁻⁴ | 4.143×10⁻⁴ | **0.994** |
| 0.33472 | 1.476×10⁻³ | 1.457×10⁻³ | **1.013** |
| 0.43635 | 3.913×10⁻³ | 3.623×10⁻³ | **1.080** |

with ladder spacing π/L_cav to 4%. (The slow drift of the ratio up the
ladder is the WKB transit picture degrading as ω_q approaches the
barrier top — expected, and it brackets the law's validity domain.) The
mode's ε is now fully geometric: **ε = ω·T(ω)/(2L_cav)**.

## 6. The tower density

Applied to the S³ tower through the #213 refocusing transit (the
antipodal caustic delivers the zonal front onto the throat once per
half-period πR — the #214 impedance matching): γₙ = T(ωₙr_h)/(πR).
At r_h/R = 0.1:

| n | 1 | 2 | 4 | 8 | 12 | 16 | 24 |
|---|---:|---:|---:|---:|---:|---:|---:|
| per-pass absorption T | 0.0017 | 0.015 | 0.16 | 0.85 | 0.99 | 0.999 | 1.000 |

IR-transparent as n³ (measured ratio 9.1 vs 8 exact-cubic — the excess
is the known finite-ω correction to the area law), crossing at
ωr_h ≈ 1, UV-black above the photon-sphere edge. With the #210 anchor
r_s ~ αλ̄_C, laboratory modes sit far into the transparent wing; the
throat blackens only at its own Compton-edge scale.

**After this PR, nothing about the absorber is chosen: its spectrum
(this PR), its address (#214: the antipode), and its ε (#214, now with
the geometric spectrum) are all read off the frozen bulk.**

## 7. Honest scope

- Classical scalar greybody: no Hawking flux; no spin/tensor channels
  (the photon's greybody differs at O(1) — follow-up).
- The tower rate uses the zonal 1D reduction (one antipodal delivery
  per half-period) that the cavity weld validates exactly; the full
  S³-with-throat matched asymptotics (3D flux spreading, ℓ-mixing at
  the caustic) is the named successor.
- r_h/R is a parameter here (0.1 in the table); the #210 anchor
  supplies the physical value.
- What the horizon does with the energy (the interior, re-emission, the
  antipodal mouth identification) is #168/#200 territory, untouched.
- Matching accuracy ~10⁻³ at the highest frequencies (the long-range
  metric phase); UV numbers carry that error bar.
- Frozen geometry, no backreaction.

## 8. What would falsify this

- σ_s(ω→0)/A_h converging to anything but 1 — the solver or the
  identification would be wrong. (Checked: 1.012 → 0.983 extrapolated.)
- Half-crossings away from the eikonal (ℓ+1)/b_c — the barrier would
  not be the photon sphere. (Checked: 1–4%.)
- Quasimode widths violating γ = T/(2L_cav) — the greybody would not be
  the operative spectral density. (Checked: 0.6–8% down the ladder.)
- A tortoise depth finite at the horizon — recurrence would survive and
  the continuum claim fail. (Checked: log-divergent, coefficient
  −0.5005 vs −1/2 exact.)

## 9. Companion probe

`experiments/closure_ledger/throat_spectral_density_probe.py` (T1–T8;
~15 s) machine-checks: flux conservation and regulator independence;
the area-theorem ratio, its monotone approach, and the ω³ slope; the
photon-sphere algebra, the eikonal crossings, the WKB ordering, and
T → 1; the tortoise log law; the quasimode ladder against the transit
law and the assembled tower table.

**Verdict:**
`THE_BANK_RETIRED_THE_SPECTRAL_DENSITY_IS_THE_THROATS_GREYBODY_IR_TRANSPARENT_BY_THE_AREA_THEOREM_UV_BLACK_ABOVE_THE_PHOTON_SPHERE_AND_THE_HORIZON_IS_THE_CONTINUUM`
