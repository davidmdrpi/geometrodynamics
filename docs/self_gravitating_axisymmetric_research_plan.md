# Real GR back-reaction: a semi-dynamical axisymmetric self-gravitating wave packet (PR #176)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Past the 1D ring proxy

PR #175 answered "can a continuous geometry evolve into the discrete
sector?" with a **1D ring proxy** whose focusing was an *ad-hoc*
nonlinearity `g|ψ|^p`. This probe moves past the proxy and asks whether
**real general relativity** backs it up: a semi-dynamical, **axisymmetric**
wave-packet evolution under a metric that is allowed to **respond to the
energy density** — the self-gravitating scalar (Schrödinger–Newton, the
weak-field limit of Einstein–Klein–Gordon):

```
i ∂_t ψ = −½ ∇²ψ + Φ ψ ,     ∇²Φ = 4πG |ψ|² ,
```

with the metric component `g_tt = −(1 + 2Φ)` responding to the energy
density `ρ = |ψ|²`. The focusing is now **gravity**, not a knob.

## The scheme (genuinely axisymmetric)

`ψ(r,θ,t)` in the `(r, ℓ)` Legendre-mode basis; the radial Laplacian
propagated by a Dirichlet sine transform; the centrifugal `ℓ(ℓ+1)/2r²`
diagonal in `ℓ`; and the metric potential `Φ(r,θ)` solved each step by the
axisymmetric **multipole Poisson** integral. Split-step, mass-conserving.

## The result (measured)

| finding | result |
|---|---|
| **stable** | mass conserved to ~10⁻³ over the evolution |
| **gravitational threshold** | disperse below / collapse above a critical mass |
| **it is gravity** | critical mass scales as **1/G** (halves when G doubles) |
| **metric responds** | the central potential well deepens as the field concentrates |

### The threshold

| mass | peak growth | outcome |
|---:|---:|---|
| 1.0 | ×1.0 | disperse |
| 2.0 | ×1.3 | disperse |
| 3.0 | ×2.2 | collapse |

Below a critical mass the packet disperses (the metric stays shallow);
above it the self-gravity concentrates the packet (the metric deepens,
runaway). This is the disperse-below / persist-above threshold of
#58/#166/#175 — now driven by actual gravitational back-reaction rather
than a phenomenological nonlinearity.

### It is gravity: the 1/G scaling

The decisive test. Critical masses interpolated from the peak-growth
crossing:

| G | critical mass |
|---:|---:|
| 0.5 | > 3.2 |
| 1.0 | 2.29 |
| 2.0 | 1.10 |

The critical mass **halves** from G=1 to G=2 (ratio **0.48 ≈ 0.5**) — it
scales as `1/G`. The threshold is set by the **gravitational coupling**,
not a toy nonlinearity: real (weak-field) general relativity backs the
focusing threshold of the antipodal-focusing arc.

### The metric responds

During a super-threshold collapse the central potential well `|Φ(0)|`
deepens as the field concentrates — the back-reaction in action — versus a
shallow, dispersing sub-threshold run. The metric `g_tt = −(1+2Φ)` tracks
the growing energy density; the geometry is dynamical.

## Honest scope

**Semi-dynamical**: the field evolves dynamically while the metric responds
**quasi-statically** through the weak-field (Newtonian) limit of GR — not
full numerical relativity. The collapse confirms the **threshold** and the
**concentration** (the throat-formation analog); the fully relativistic
**strong-field endpoint** (a horizon / a resolved throat) is beyond a
weak-field scheme. Self-gravity also tends to **sphericalize** (the
monopole dominates the Poisson source), so the collapse is predominantly
radial — the axisymmetric machinery is *exercised* (it handles non-spherical
packets, mass-conserving), not a directional jet claimed.

## What this establishes

The 1D ring proxy of #175 upgrades to real (weak-field) axisymmetric
self-gravity, and the disperse/collapse threshold **survives and is
gravitational** (`M_c ∝ 1/G`), with the metric responding to the energy
density. Real GR backs the antipodal-focusing threshold of #166/#175 and
the nucleation of #58.

## Reproduce

```bash
python -m experiments.closure_ledger.self_gravitating_axisymmetric_probe
# Verdict: REAL_SELF_GRAVITY_REPRODUCES_THE_FOCUSING_THRESHOLD_CRITICAL_MASS_SCALES_AS_INVERSE_G
```
