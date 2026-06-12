# Flavor phase addendum: Hopf CP derivation and full CKM realization (PR #162)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. A
> consolidating addendum in the #131/#150 convention.

The CP-phase arc #156→#161 ended in a very different place from where it
started. This PR consolidates it: one keystone from every step re-verified
in a single run (`flavor_phase_addendum_probe`), and the full account
written into `docs/THESIS.md` — the interim #158–#159 postscript replaced by
the complete addendum subsection, and the #157 card's quark-CP row updated
from "calibrated; shape open" to "**derived** (φ_h = π/k₅); full dataset
realized".

## The keystones, re-run together

| arc step | keystone | re-verified |
|---|---|---|
| #158 relocation | exact unitarity + quartet-consistent J at the Hopf phases; the excluded route's unitarized J | 7e-16 / 2e-18 / J ≈ 0 |
| #159 transport scale | sector-arc phases k·π/k₅ (unit-circle at the k = 5 branch point) | exact to 1e-12 |
| #160 algebra | the Weyl commutator e^{2πi/k₅} on the capacity-k₅ space | exact |
| #161 realization | the down-dominant solution re-solved: nine observables | all ≤ 1%, masses 1e-14 |

## The consolidated chain (written to the thesis)

| ingredient | value | status |
|---|---|---|
| rate ½ | A_φ(χ=0) — the spin-½ factor | derived (connection) |
| sign ± | Z₂ partition orientation | derived (#63 C-swap) |
| winding dk | max(k, k′) | locked (the v3 mass calibration) |
| arc 2π/k₅ | the Weyl commutator quantum | derived (#160) |
| **φ_h = π/k₅** | 0.6283 | **derived** (#159; alternatives excluded) |

## The bookkeeping (written to the thesis)

Thirteen probes #149–#161: **net zero new inputs** (the #156 input consumed,
then returned by the #159 derivation) and **one modelling knob retired**
(the #151 β interpolation, #152). The quark flavor-CP sector stands as a
consistency statement: locked masses + derived CP phase + the realized
nine-observable dataset + complete re-lock targets.

## Remaining (unchanged from #161)

The knob-level v3+CP re-lock against the tabulated targets, and the lepton
sector's anarchic draw.

## Reproduce

```bash
python -m experiments.closure_ledger.flavor_phase_addendum_probe
# Verdict: FLAVOR_PHASE_ADDENDUM_ALL_KEYSTONES_REVERIFIED_THESIS_CONSOLIDATED
```
