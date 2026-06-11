# Flavor-sector synthesis capstone (PR #157)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

The synthesis capstone for the flavor arc #149–#156, in the #150/#131
convention: one keystone re-verified from every arc member, the complete
flavor card assembled, the bookkeeping audited, and the card added to
`docs/THESIS.md` ("The flavor sector, assembled").

## The arc

bracket the residual (#149) → test the mixing/anarchy hypothesis (#151) →
derive the channel-dominant saddle (#152, the β knob retired) → extract both
mixing matrices (#153 PMNS, #155 CKM) → complete CP in both sectors (#154
Majorana, #156 quark).

## Keystones re-run together

| arc member | keystone | re-verified value |
|---|---|---|
| #149 | required ε₃/ε₂ (pure-bounce inversion) | 1.435 |
| #151 | observed r₃₂ percentile in the ensemble | ~77th |
| #152 | exact 2-channel element follows the max rule | t/(Δ_max/2) ≈ 0.2; geo fails ×~9 |
| #153 | PMNS angle percentiles (s12/s23/s13) | natural |
| #154 | m_ββ median; exact φ_ℓ invariance | ≈ 3.2 meV; 0.0 |
| #155 | V_cb from the mass-locked blocks | 0.0377 |
| #156 | ceiling identity; calibration point | exact; J = target |

## The flavor card

Masses (normal ordering derived; m₁ predicted light; Σm_ν ≈ 58.8 meV
falsifiable), the ε_n spread (derived — β retired), the three PMNS angles
(anarchy-natural), CP generic in the lepton sector (Dirac and Majorana),
m_ββ ≈ 3.2 meV (falsifiable), the CKM (out-of-sample, zero inputs, V_cb/V_ts
stiff at 10%), and quark CP (calibrated, with β = 22° the Hopf-phase
acceptance test).

## The bookkeeping

Eight probes: **one input consumed** (the quark CP phase content, #156 — the
flavor puzzle's CP entry made explicit in the #150 budget) and **one
modelling knob retired** (the #151 β interpolation, derived by #152). The
net modelled-assumption count *decreased* while the sector was assembled.

## The falsifiable targets

1. Σm_ν 58.8 vs 61.1 meV — cosmology at ~1–2 meV precision discriminates.
2. An m_ββ detection above ~10 meV falsifies the ensemble.
3. β = 22° — the acceptance test for the Hopf-connection φ_q(k).
4. The J ceiling must rise to 3.5×10⁻⁵ when the soft V_us/V_ub land.
5. V_cb = 0.038, stiff at 10%.

## Scope

The capstone consolidates; the remaining residuals (the anarchic draw, the
CP phase content, the soft V_us/V_ub direction, the O_geom e-row) are
listed with their falsifiable exits, not removed.

## Reproduce

```bash
python -m experiments.closure_ledger.flavor_sector_synthesis_probe
# Verdict: FLAVOR_SECTOR_ASSEMBLED_ONE_INPUT_ONE_KNOB_RETIRED_FIVE_TARGETS
```
