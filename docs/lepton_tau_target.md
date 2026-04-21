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

## Basin-of-attraction map (local gradient probe)

Using `scripts/map_basin_k_uplift.py` initialized at the best exact-root point:

- Optimizer converged back to essentially the same solution:
  - `phase≈0.2947743`, `transport≈25`, `pinhole≈0`,
  - `resistance≈0.2055224`, `beta≈72`,
  - `action_base≈6.6643244`,
  - `mu/e` exact within numerical tolerance,
  - `tau≈1783.26 MeV` (0.36% error).

1% perturbation sensitivity (holding other variables fixed) is extremely sharp:

- `action_base` ×0.99 → mu error ~237%, tau ~6034 MeV
- `action_base` ×1.01 → mu error ~41%, tau ~1046 MeV
- `resistance` ×0.99 → mu error ~35%, tau ~2411 MeV
- `resistance` ×1.01 → mu error ~20.5%, tau ~1416 MeV

Interpretation: this exact solution currently sits on a **narrow attractor needle**, not a broad basin.
The next physical step is to widen this basin (e.g., constrain action/base-resistance relation from a
geometric identity) so the root is structurally robust, not just numerically reachable.

## Hard S^3 lock experiment (`action_base = 2π`)

Implemented a geometric lock that removes `action_base` as a free optimizer slider:

- `action_base` is fixed to the S^3 circumference invariant `2π ≈ 6.283185307180`.
- calibration and sweep CLIs now pass this value explicitly unless overridden.

With this lock enabled:

- exact mu/e roots still exist in the tested `beta ∈ [40, 140]` window,
- best coarse-grid exact-root point in that window is near:
  - `beta ≈ 140`,
  - `phase ≈ 0.001`,
  - `transport ≈ 25`,
  - `hard_pinhole_gamma ≈ 22.5`,
  - `resistance ≈ 0.2154574`,
  - `tau ≈ 1607.9 MeV` (tau error ~9.51%).

Local 1% resistance perturbation around that locked exact-root solution:

- `resistance × 0.99` → mu error ~7.84%, tau ~1735.17 MeV (tau error ~2.35%)
- `resistance × 1.01` → mu error ~6.75%, tau ~1498.31 MeV (tau error ~15.68%)

Compared to the pre-lock needle behavior (20–35% mu error for ±1% resistance shifts),
this indicates a meaningfully wider local resistance basin under the hard topological lock.
The remaining task is to recover the sub-percent tau fit while preserving this broader basin.

## Dense locked local scan near β≈140 (beta window 130→150)

Executed dense local locked scan (`scripts/refine_locked_tau.py`) over:

- `beta ∈ [130, 150]` with fine stepping (tested down to `Δbeta=0.5` in code path, and
  smoke/quick reporting at `Δbeta=1.0`),
- constrained small cavity phase window `phase ∈ [0.001, 0.01]`,
- transport and pinhole windows centered near the prior locked optimum.

Best local exact-root candidate found in this window:

- `beta ≈ 150.0`,
- `phase ≈ 0.001`,
- `transport ≈ 24.9`,
- `hard_pinhole_gamma ≈ 21.75`,
- `resistance ≈ 0.22338`,
- `tau ≈ 1732.31 MeV` (tau error **~2.51%**),
- local resistance basin remained broader than the pre-lock needle
  (`mu/e` error ~8.34% and ~7.12% at ±1% resistance).

Interpretation: the lock preserves the wider basin trend, but this particular beta window did not
yet recover sub-percent tau. A natural next extension is to scan a slightly higher beta shoulder
(e.g. `150→165`) and test integer-anchored geometric beta families.

Geometric beta anchor helper was added:

- `derive_geometric_beta(...) = 4π (R_MID / l_core) * winding_integer * scale`.
- For `winding_integer=5`, this gives `beta ≈ 133.286`, which lies in the explored local window.

## Higher beta shoulder scan (`150→165`) + integer families

Extended local locked scan to `beta ∈ [150, 165]` (small-phase window `0.001→0.01`):

- best exact-root candidate in this shoulder:
  - `beta ≈ 157`,
  - `phase ≈ 0.001`,
  - `transport ≈ 25.1`,
  - `hard_pinhole_gamma ≈ 22.5`,
  - `resistance ≈ 0.217875`,
  - `tau ≈ 1778.89 MeV`, tau error **~0.114%** (sub-percent recovered),
  - resistance basin check remained broad relative to the pre-lock needle
    (`mu/e` error ~7.90% / ~6.80% at ±1% resistance).

Integer-anchored geometric beta family checks (nearest scanned beta in this shoulder):

- `n=5`: `beta_geom≈133.286` → nearest scanned `beta=150`, tau error ~2.51%
- `n=6`: `beta_geom≈159.944` → nearest scanned `beta=160`, tau error ~1.86%
- `n=7`: `beta_geom≈186.601` → nearest scanned `beta=165`, tau error ~4.78%

Interpretation: the shoulder scan succeeded in recovering sub-percent tau while preserving the broader
basin trend. The integer-family anchor near `n=6` sits directly inside this high-performing shoulder.

## Hard lock test: `beta = 50π` (exact `200π` tau uplift)

Added `scripts/lock_beta_50pi_probe.py` to freeze:

- `action_base = 2π` (existing lock),
- `beta = 50π` exactly, so for `k=5`: `4*beta = 200π` and uplift quanta = `100 × (2π)`.

With beta hard-locked to `50π`, optimizing only `(phase, transport, pinhole, resistance)` yielded:

- `phase ≈ 0.001000001`,
- `transport ≈ 25.10006`,
- `hard_pinhole_gamma ≈ 22.49962`,
- `resistance ≈ 0.217869`,
- `mu/e` error ~`0.000001%`,
- `tau ≈ 1779.72 MeV`, tau error **~0.161%**.

So the strict topological `200π` lock remains in the sub-percent regime and keeps the broad basin behavior,
but did not hit literal `0.000%` tau in this run.
