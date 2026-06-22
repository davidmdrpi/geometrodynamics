# Antipodal wave-packet focusing threshold (PR #166)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The gap this closes

The THESIS section "Why antipodal focusing matters" asserts that a
wavefront on a closed S³ does not dissipate but **reconverges at the
antipode**, and that a strong enough antipodal caustic nucleates a throat
— *"the focus is the trigger; the particle is the persistent topological
response."* The threshold **energy** (`2 m_e c²`, twice the lowest stable
throat) was derived *statically* from the self-energy functional
`E(R)=A/R+B·R²` in `pair_production_threshold_probe` (PR #58). But the
antipodal focusing **itself** — whether a packet actually refocuses, when,
how sharply, with what gain — was asserted, never computed. This probe
computes it.

## The exact reduction

The zonal sector of S³ (fields of the polar angle `χ∈[0,π]`, antipode at
`χ=π`) reduces **exactly** to a 1D wave on a string: with the reduced
field `f(χ,t)=sin(χ)·ψ(χ,t)`, the modes are `sin((ℓ+1)χ)` with Dirichlet
ends, and for the **conformal** scalar the frequencies are the equally
spaced tower `ω_ℓ=(ℓ+1)/R` (cf. PR #165). The physical field `ψ=f/sinχ`
carries the **geometric focusing factor** `1/sinχ`: as a wavefront
converges on the antipode (`sinχ→0`) its amplitude is amplified — the
caustic.

## What is computed (new) vs inherited

| | result |
|---|---|
| **Exact refocus** at `t=πR` | identity `ψ(χ,πR)=−ψ(π−χ,0)` to **3e-15**; amplitude recovery ×1.0000 |
| **Full revival** at `t=2πR` | `0e+00` — the sub-threshold focus passes through and re-disperses |
| **Conformal required** | conformal recovery ×1.000 vs minimal-coupling ×0.877 (dephased) |
| **Caustic** | density `∝1/sin²χ`; antipodal-point gain rises steeply with `ℓ_max ~ R/R_MID` |
| **Threshold (inherited)** | focused energy ≥ `E(R*)=m_e c²`; pair (`Σc₁=0`) → `2 m_e c² = 1.022 MeV` (PR #58) |

The **conformal** coupling that makes the S³ vacuum tower equally spaced
(PR #165) is the very same coupling that makes the antipodal caustic
sharp. The caustic `1/sin²χ` gain — regularized by the spectral cutoff
`ℓ_max ~ R/R_MID` (the finest wavefront feature is the throat scale) —
lets a **delocalized, S³-wide (diffuse) wave reconcentrate its energy onto
the throat scale**, the concentration factor scaling with `R/R_MID`. This
is the dynamical bridge from a spread wave to a local nucleation density.

## Honest scope (no rigging)

- **Computed** from first principles (linear conformal wave propagation):
  the exact refocusing, its timing `t=πR`, the `2πR` revival, the
  conformal requirement, and the caustic gain.
- **Not simulated**: the *nonlinear* throat formation — linear waves
  cannot nucleate a topology change. The probe maps the **trigger** and
  applies the threshold criterion; it does not watch a throat form.
- **Inherited, not re-derived**: the `2 m_e c²` value and the
  disperse-below / persist-above barrier (PR #58). No constant is fit
  here; the only free choices are the packet's launch angle and width,
  which set the focus location and sharpness, not the threshold.

## Reproduce

```bash
python -m experiments.closure_ledger.antipodal_focusing_threshold_probe
# Verdict: ANTIPODAL_FOCUS_EXACT_AT_PI_R_CAUSTIC_TRIGGERS_NUCLEATION_AT_2MEC2
```
