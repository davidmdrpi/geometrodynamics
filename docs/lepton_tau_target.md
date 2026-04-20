# Next target analysis: lifting the k=5 (tau) sector

## Context

The current max-winding strong-mixing branch can recover exact 
\(\mu/e = 206.768282987\) roots, but tau remains far below 1776.86 MeV.

## Scan used for this analysis

Using `calibrate_grid(...)` with:

- `depth_cost_mode="tunnel_only"`
- `winding_mode="max"`
- `phase_steps=9`, `transport_steps=9`, `pinhole_steps=9`
- `transport ∈ [20, 30]`, `pinhole ∈ [0, 30]`, `resistance ∈ [0.2, 0.5]`

## Results

- Total candidates: **729**
- Exact mu/e roots: **23**
- Best exact-root tau (minimum tau relative error): **195.31 MeV**
- Tau range across exact roots: **1.13 → 195.31 MeV**
- Mean tau over exact roots: **71.21 MeV**

Implied uplift required for the best exact-root candidate:

\[
\text{uplift} = \frac{1776.86}{195.31} \approx 9.10
\]

So the next target is not finding mu/e roots (already solved), but adding a mechanism that increases
k=5 by roughly **one order of magnitude** while keeping k=3 anchored.

## Sensitivity hints from exact-root set

- Correlation(`hard_pinhole_gamma`, `tau_pred`) ≈ **+0.49** (moderately positive)
- Correlation(`transport_strength`, `tau_pred`) ≈ **−0.19** (weakly negative)

Interpretation: stronger pinhole/crossing contributions appear to help tau uplift more than transport
in this regime.

## Recommended immediate implementation target

Add a **k-selective diagonal uplift term** that is minimal at k=1 and k=3 but stronger at k=5, e.g.

\[
\Delta E_k = \beta \,(k-3)^2\,\Theta(k-3)
\]

or a topological variant tied to crossing combinatorics that activates predominantly for k≥5.

This is the cleanest next step because it directly addresses the observed failure mode (k=5 under-lift)
without sacrificing exact mu/e roots already available in the current architecture.

## k-selective uplift experiment (implemented)

Applied uplift term in generation diagonal:

\\[
\\Delta E_k = \\beta \\max(0, k-3)^2
\\]

with `k_uplift_beta = 80` in max-winding tunnel-only scan (`7×7×7`):

- exact mu/e roots still present (**9** exact roots),
- tau improved from ~170–229 MeV to a range including **~1332 MeV** best exact-root candidate,
- best exact-root tau relative error reduced to **~25.0%**.

Conclusion: the k-selective uplift is moving the model in the correct direction while preserving the
mu/e manifold; next step is refining `k_uplift_beta` (and possibly mild k≥5 nonlinear terms) to close
the remaining ~25% tau gap.

## High-density beta sweep (exact mu/e enforced)

Using `scripts/sweep_k_uplift_beta.py` with:

- `beta ∈ [60, 140]` (21 points),
- `phase_steps=5`, `transport_steps=5`, `pinhole_steps=5`,
- `depth_cost_mode=\"tunnel_only\"`, `winding_mode=\"max\"`.

Result:

- exact mu/e roots for all tested betas in this window,
- best beta found: **β ≈ 72.0**,
- `tau_pred ≈ 1783.26 MeV`, tau relative error **~0.36%**,
- corresponding root remains on stable manifold:
  - `phase_per_pass ≈ 0.2947743`,
  - `transport_strength = 25`,
  - `hard_pinhole_gamma = 0`,
  - `resistance_scale ≈ 0.2055224`.

This closes the previously identified tau gap while preserving exact mu/e calibration.
