# Explicit Hopf-transport derivation of φ_h = π/k₅ (PR #159)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The CP phase
> scale is computed by explicit parallel transport on the Hopf bundle.

PR #158 relocated the quark CP phase to the Hopf-fiber transport of the
same-partition shell couplings and flagged φ_h = π/k₅ as a candidate (five
CP observables, zero parameters) pending an independent transport
derivation. This PR supplies it — the #152 modelled→derived path — and
**returns the #156-consumed input** to the budget.

## The derivation chain

| ingredient | status | content |
|---|---|---|
| the rate ½ | derived | A_φ(χ=0) = ½ — the spin-½ factor (`hopf/connection.py`) |
| the sign ± | derived | the two Z₂ partition classes traverse the fiber oppositely (#63 C-swap; explicit opposite transport conjugates — verified) |
| the winding content dk = max(k,k′) | locked | the SAME `winding_mode='max'` rule the **mass** calibration locked — phase and magnitude rules share their dk (independent corroboration) |
| the arc 2π/k₅ | identification | one winding-sector of the fiber (capacity k₅ = 5, #73/#126) — **data-selected over all principled alternatives** (below) |

    phase = ± dk · (½) · (2π/k₅) = ± dk · π/k₅

**Explicit path-ordered transport** (module connection, 20k steps): full
circuit at k = 1 → exactly π (the spinor sign flip — the module's own
consistency anchor); sector arc → k·π/k₅ for k = 1, 3, 5, **exact to
1e-15**.

## The exclusion scan (observed: J/t = 1.0, β = 22.2°, γ = 65.9°, sin δ = 0.887)

| candidate | J/target | β | γ | sin δ | passes? |
|---|---|---|---|---|---|
| **π/k₅ = π/5** | **0.969** | **22.8** | **63.5** | **0.888** | **✓ (all five)** |
| π/3 (generations) | −0.020 | −0.9 | −1.3 | −0.022 | ✗ (CP killed) |
| π/4 | 0.661 | 20.8 | 42.2 | 0.668 | ✗ |
| π/6 | 1.140 | 23.2 | 78.2 | 0.971 | ✗ |
| 2π/5 | −0.499 | −17.7 | −32.8 | −0.538 | ✗ (sign flip) |
| π/10 | 1.120 | 19.0 | 114.7 | 0.901 | ✗ |

π/k₅ is the **unique survivor** — the identification is data-selected among
principled candidates: the anti-numerology discipline in its positive mode.

## The derived prediction (no calibration anywhere)

J at 0.969 of target · (β, γ, α) = (22.8°, 63.5°, 93.8°) vs (22.2°, 65.9°,
91.9°) · sin δ = 0.888 vs 0.887 · masses shifted 0.09% · V_cb untouched
(0.0377). **Quark CP is derived.**

## Budget impact

The #156-consumed input (the quark CP phase content) is **returned**. The
flavor card's last open row closes (quark CP: derived, sector identification
documented). Net flavor-arc bookkeeping #149–#159: **inputs +0** (one
consumed in #156, returned here), **modelling knobs −1** (the β
interpolation, #152). A postscript is added to the THESIS flavor section.

## Scope

The final geometric mile — deriving the hop arc (one capacity sector) from
the explicit shell wavefunctions rather than from the capacity structure +
data selection — is flagged, the same status the #152 saddle had before its
controlled model. The soft V_us direction stands.

## Reproduce

```bash
python -m experiments.closure_ledger.phi_h_transport_derivation_probe
# Verdict: PHI_H_PI_OVER_K5_DERIVED_TRANSPORT_EXACT_ALTERNATIVES_EXCLUDED
```
