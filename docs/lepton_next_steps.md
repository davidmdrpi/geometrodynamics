# Lepton instanton surrogate: recommended next steps

## 1) Diagnose why mu/e saturates near ~10-12

The calibration scans repeatedly hit a boundary optimum (high transport, low resistance) while
mu/e remains far below 206.8. This usually means the objective is constrained by model structure,
not search range.

**Immediate check:**
- Add a diagnostic sweep that records `eig_k` and off-diagonal amplitudes `exp(-S_ij)` per depth.
- Verify whether couplings are effectively frozen (near-zero or near-constant) across the scan.

## 2) Separate diagonal growth from tunneling suppression explicitly

Right now both the diagonal term and tunneling action carry depth costs. This can lock ratios into
narrow bands.

**Recommended experiment:**
- Introduce a controlled ablation mode:
  - (A) depth cost only on diagonal,
  - (B) depth cost only in tunneling action,
  - (C) both.
- Compare mu/e sensitivity surfaces for A/B/C.

## 3) Promote winding cost option from `delta_k` to `max_k`

Current action uses `delta_k = |i-j|`. A higher-generation transition may need full higher-k cost.

**Recommended experiment:**
- Add a switch `winding_mode in {"delta", "max"}` where
  - `delta`: `delta_k = |i-j|`
  - `max`: `delta_k = max(i+1, j+1)`
- Re-run the same grid and compare reachable mu/e envelope.

## 4) Tie crossing cost to identified-node count, not just gamma magnitude

Crossing effects are currently parameterized by `hard_pinhole_gamma`; include combinatorics from
actual identified crossings so cost scales with topology, not a single knob.

**Recommended experiment:**
- Add `n_identified`-dependent term in action, e.g. proportional to observed crossing multiplicity.

## 5) Upgrade search strategy

Boundary-hugging solutions indicate grid limits or non-convexity.

**Recommended experiment:**
- Keep coarse grid for broad map, then run local optimizer on best seeds.
- Report Pareto front: `(mu/e error, tau error, action regularity)`.

## 6) Add one script-level integration test

The CLI is central to iteration speed.

**Recommended test:**
- Run a tiny deterministic grid and assert stable ordering + finite outputs.

## 7) Decision gate

After the above diagnostics, choose one path:
- **Path A:** identified-node (pinhole/crossing) dominant mechanism.
- **Path B:** instanton action dominant mechanism.

Gate criterion:
- If mu/e > 100 becomes reachable under physically constrained settings, continue that path.
- Otherwise, revise action decomposition before further parameter tuning.
