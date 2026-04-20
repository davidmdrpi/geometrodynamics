# Lepton instanton ablation: results and interpretation

## Scan configuration

Used the current calibration machinery with:

- `phase_steps=6` on `[1e-3, π/8]`
- `transport_strength ∈ [20, 30]` with 6 points
- `hard_pinhole_gamma ∈ [0, 30]` with 6 points
- `resistance_scale ∈ [0.2, 0.5]` (solved via root/minimization)
- depths `(1, 3, 5)` and electron calibration

Evaluated three ablation modes:

- `both` (depth cost on diagonal and tunneling)
- `diag_only`
- `tunnel_only`

## Best candidate per ablation mode

| mode | best mu/e | mu/e rel err | tau_pred (MeV) | tau rel err |
|---|---:|---:|---:|---:|
| both | 11.6455 | 94.37% | 12.9808 | 99.27% |
| diag_only | 11.6473 | 94.37% | 12.9808 | 99.27% |
| tunnel_only | 24.7666 | 88.02% | 14.0017 | 99.21% |

## Boundary behavior

For all three modes, top-ranked points concentrate at search boundaries:

- `resistance_scale` at lower bound (0.2) for all top-20 entries
- `hard_pinhole_gamma` at upper bound (30) for all top-20 entries
- `transport_strength` frequently at upper edge (30)

This indicates the current objective is **still capacity-limited by model structure**, not by local optimizer noise.

## Interpretation

1. **Tunneling-side depth cost carries more leverage** than diagonal depth cost.
   The `tunnel_only` mode almost doubles best mu/e relative to `both/diag_only`.
2. **Diagonal-only cost mostly shifts absolute scale**, not generation separation.
3. **Crossing/pinhole term is being used maximally** by the optimizer, implying topological interaction terms are beneficial but still insufficient in current form.
4. **No exact mu/e roots** were found in this regime, and tau remains strongly underpredicted.

## Recommended immediate follow-up

- Keep `depth_cost_mode="tunnel_only"` as the primary branch for next experiments.
- Expand winding mode from `delta_k` to `max(i+1, j+1)` and compare envelopes.
- Add reporting of effective off-diagonal amplitudes `exp(-S_ij)` to confirm whether suppression is still too flat between k=1 and k=3 sectors.

## Winding-mode envelope comparison (new)

Using `depth_cost_mode="tunnel_only"` and the same coarse grid:

- `winding_mode="delta"` best mu/e ≈ 24.77 with k=3 envelope
  `min≈2.5e-24, mean≈1.5e-07, max≈2.3e-07`.
- `winding_mode="max"` best mu/e ≈ 23.56 with k=3 envelope
  `min≈7.6e-52, mean≈7.1e-25, max≈2.1e-24`.

Interpretation: `max` mode strongly suppresses off-diagonals (as intended), but the current
observable/selection pipeline still maps to almost identical mu/e outcomes. Next step is to make
the spectrum observable explicitly sensitive to these envelope differences.

## Strong-mixing update (linear diagonal + reduced suppression)

After flattening the diagonal and using stronger off-diagonal mixing:

- `winding_mode="max"` no exact mu/e roots on the coarse grid, best mu/e ≈ 10.54.
- `winding_mode="delta"` found an **exact mu/e root** on the same grid:
  - `mu/e = 206.768282987...`
  - `tau_pred ≈ 270.48 MeV` (still under target by ~84.8%).
  - `k=3` envelope moved to order `1e-1` (`mean≈3.89e-1`), i.e. no longer in weak-mixing limit.

Interpretation: the weak-mixing bottleneck is broken (off-diagonals now active). The remaining
issue is shaping the **k=5 sector** without collapsing the successful mu/e fit.

## Quadratic-diagonal + rescaled-max update

With generation-block diagonal updated to `action_base + resistance_scale * k^2` and a softer
suppression factor for `winding_mode="max"`:

- coarse scan (`7×7×7`) found **11 exact mu/e roots** under `winding_mode="max"`.
- best exact-root candidates predict `tau ≈ 170–229 MeV` (still far below 1776 MeV).
- `k=3` envelope remains active (`mean≈1.43e-1`), so mixing is healthy.

Interpretation: the `max` mode can now recover the mu/e root, but the model still under-lifts the
`k=5` sector; additional k-dependent self-energy shaping is required.
