# The nonlinear antipodal focusing PDE sandbox (PR #175)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question

**Can a continuous time-dependent geometry actually evolve into the
discrete sector?** PR #166 simulated the *linear* antipodal focusing (a
conformal packet refocuses exactly at the antipode) and explicitly deferred
the nonlinear throat formation; PR #174 showed the discrete sector (the
odd-k winding) sits *outside* the continuous deformation manifold. This
probe is the dynamical capstone: a nonlinear antipodal-focusing PDE
sandbox.

## The sandbox

A complex field `ψ(χ,t)` on the antipodal ring (the great circle of S³),
evolved by a focusing nonlinear Schrödinger equation

```
i ∂_t ψ = −∂_χχ ψ − g |ψ|^p ψ      (p = 2 cubic, p = 4 critical),
```

split-step Fourier (mass-conserving). The **discrete sector** is the winding
number `Q = (1/2π) ∮ d(arg ψ)` — an integer topological charge, the proxy
for the throat's quantized `k`.

## The answer (measured)

| finding | result |
|---|---|
| **smooth evolution conserves Q** | Q = 1.00000 → 1.00000 while \|ψ\|>0 (mass conserved) — the discrete sector is **locked out** of continuous evolution (#174, dynamically) |
| **the gate is at the antipode** | changing Q forces an amplitude-zero node (min\|ψ\|→0 along the interpolation), located **exactly at χ=π** — the focus |
| **the focusing threshold** | mass <~1.6 disperses (continuous); mass >~1.9 concentrates toward the core (nucleation) — the #58/#166 threshold, now nonlinearly simulated |
| **the jump is quantized** | crossing the antipodal node changes Q by exactly **±1** — a discrete response to the smooth focusing drive |

### Winding conserved (the discrete sector is locked out)

A smooth Q=1 field (`|ψ|>0` everywhere) keeps its winding *exactly* under
the nonlinear evolution, with mass conserved to ~1e-6 and `min|ψ|` staying
well above zero. A continuous geometry cannot smoothly deform into a
different discrete sector — the dynamical confirmation of #174's topological
separation.

### The topological gate, at the antipodal focus

The winding is a homotopy invariant of maps to `ℂ∖{0}`, so Q can change only
across an amplitude-zero **node**. Interpolating `ψ_s = (1−s)·1 + s·e^{iχ}`
between Q=0 and Q=1, the minimum amplitude is forced to `0` at `s=½`,
located **exactly at the antipode χ=π** — the focus. The discrete sector is
gated by a singular core at the antipodal caustic, and the antipodal
focusing (#166) is precisely the dynamics that drives the field there.

### The focusing threshold (critical-NLS mass scan)

| mass | peak growth | outcome |
|---:|---:|---|
| 0.99 | ×1.0 | disperse |
| 1.59 | ×1.0 | disperse |
| 1.94 | ×1.5 | concentrate |
| 2.53 | ×5.1 | concentrate |

Below the critical mass the field disperses (stays smooth, Q frozen,
continuous); above it concentrates toward the core (the nucleation onset).
This is the disperse-below / persist-above threshold of #58/#166, now
actually simulated nonlinearly.

## The synthesis

A continuous geometry reaches the discrete sector **only by developing a
focusing singularity at the antipode — never by smooth deformation**. Smooth
evolution conserves the winding; the gate is an amplitude-zero core forced
at the antipodal focus; the nonlinear focusing reaches that core only above
a critical mass; and the winding jump there is quantized ±1. The sandbox
confirms #174's topological separation dynamically and realizes #166's
threshold and #58's nucleation.

**Answer: yes, but only through the caustic.**

## Honest scope

A reduced 1D zonal/ring model: `Q` is the winding proxy for the discrete
`k`, the collapse core is the proxy for throat nucleation (not the full
4D/5D GR throat), and the critical-NLS collapse is marginal (the threshold
is resolved, the singular core is not fully resolved). The conceptual
answer — smooth conserves, the gate is at the antipodal node, the threshold
sets reachability, the jump is quantized — is robust; the specific numbers
are model-dependent.

## Reproduce

```bash
python -m experiments.closure_ledger.nonlinear_antipodal_focusing_pde_probe
# Verdict: CONTINUOUS_REACHES_DISCRETE_SECTOR_ONLY_THROUGH_ANTIPODAL_FOCUSING_SINGULARITY
```
