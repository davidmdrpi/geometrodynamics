# The spectral deformation test: the {1,3,5} ladder on the Berger-squashed S³ (PR #192)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Upgrading #183 from algebra to spectrum

PR #183 proved the odd-k {1,3,5} generation sector is protected by
metric-independent **algebra** — the deck determinant, the Pin⁻
spin-closure class `½ tr T² = −1`, the odd-parity grading. But algebra is
not spectrum: the claim the thesis actually needs is that the **ladder** —
three positive, ordered levels with the observed mass ratios — survives a
metric deformation. #183 cannot distinguish "the topology protects the
spectrum" from "the round metric was doing spectral work the topology
can't protect."

This probe distinguishes them, with machinery the repo already has: the
#165 R-unification audit's **genuine SU(2) Berger-sphere spectrum**
(`Δ(j,m) = 4j(j+1) + 4m²(λ⁻²−1)`, imported, validated to collapse to the
round tower at λ=1). The locked lepton Hamiltonian's geometric ingredients
are rebuilt on the Berger-squashed S³_λ (fiber × λ, base S² round) and the
ladder is tracked as λ moves off 1.

- **Pass:** the {1,3,5} structure and mass ratios deform smoothly over a
  finite λ window.
- **Fail:** the ladder breaks at infinitesimal squash — the round metric
  was doing spectral work the topology can't protect.

Either outcome is a real result. The measured outcome is **both**, cleanly
separated.

## The ingredient map (declared before the sweep)

The Berger squash multiplies the Hopf-fiber length by λ and leaves the
base S² and the Hopf **connection** untouched. Each locked ingredient of
the generation block is classified accordingly:

| ingredient | locked value | rides on | under S³_λ |
|---|---|---|---|
| `action_base` | 2π | the fiber circumference | 2πλ |
| `action_slope` | 0.5 | per-winding tunnel path | 0.5λ |
| `resistance_scale` | 0.2179 | transport along the winding | ×λ |
| `k_uplift_beta` | 50π | fiber 2π-quanta | ×λ |
| `phase_per_pass` | 0.001 | the Hopf **holonomy** (connection) | fixed |
| `hard_pinhole_gamma` | 22.5 | the round base S² | fixed |
| `transport_strength` | 25.1 | attempt-frequency prefactor | fixed |

The two ambiguous assignments (resistance, uplift) are **flipped** in a
robustness control, plus a minimal map where only the unambiguous 2π fiber
circumference rides λ. The conclusion survives all four maps.

## The pass: no infinitesimal breakdown

- At λ=1 the deformed Hamiltonian **is** the locked baseline
  (`solved_lepton_masses_mev` reproduced to machine precision).
- The {1,3,5} structure — three positive, ordered levels — persists over
  the finite window **λ ∈ (0.986, ≥3]** (dense check, 81 points on the
  stretch side).
- The response at the round point is **linear**: the central-difference
  slope d(μ/e)/dλ converges to a finite limit (−1.47×10⁴, stable to <0.1%
  between h=10⁻⁴ and 10⁻⁵). No jump, no divergence — an ordinary
  differentiable spectral flow. The #183 fail-mode is ruled out.

## The discovery: the e–μ hierarchy is metric-fine-tuned

The squash-side boundary is the **electron level's zero crossing** at

```
λ_break = 0.98598   (a 1.4% fiber squash — finite, but close)
```

The electron eigenvalue at the round point is a **near-zero**: 0.1996 in
action units against μ = 41.26 and τ = 694.98 — the celebrated μ/e ≈ 207
is the ratio of a fine-tuned near-cancellation to a normal level. The
Berger deformation exposes it quantitatively:

```
d ln(μ/e)/dλ = −70.9   ≈   −1/(1−λ_break) = −71.3
```

— the sensitivity **is** the inverse distance to the spectral boundary,
and the identity holds in **every** ingredient map (71/71, 64/64, 70/70,
31/31). By contrast τ/μ is **gentle**: d ln(τ/μ)/dλ = +0.82 (O(1) or less
in all maps).

| map | λ_break | d ln(μ/e)/dλ | 1/(1−λ_break) | d ln(τ/μ)/dλ |
|---|---:|---:|---:|---:|
| default (all fiber) | 0.98598 | −70.92 | 71.34 | +0.822 |
| flip resistance → blind | 0.98435 | −63.61 | 63.89 | +0.895 |
| flip uplift → blind | 0.98573 | −69.68 | 70.06 | −0.081 |
| minimal (only 2π base) | 0.96823 | −31.32 | 31.47 | −0.143 |

## What this means for the thesis claim

The topology guarantees **three generations** — the {1,3,5} structure
survives a finite metric window with smooth, linear response, so the round
metric is not smuggling in the sector structure. But the topology does
**not** guarantee the mass hierarchy: the electron's near-zero is
round-metric spectral work, 1.4% of squash from a spectral catastrophe,
and τ/μ is metric-robust while μ/e is not. The protection claim and the
fine-tuning are now **separated by measurement** — which sharpens the
thesis rather than weakening it.

## Honest scope

- The lepton Hamiltonian is the repo's locked instanton-transition
  **surrogate** (the LEPTON_BASELINE constants), not a first-principles
  wave operator on S³_λ; what is deformed is its geometric ingredients
  under the declared fiber/base map.
- The fiber/base classification of two ingredients is ambiguous; the
  conclusion is invariant under flipping them (T7), but the map itself is
  a modeling choice.
- The Berger spectrum used to justify the fiber scaling is the genuine
  SU(2) result (#165); the surrogate couples to it only through the
  ingredient map, not through a mode-by-mode expansion — the
  field-theoretic ladder on S³_λ (the actual deformed wave operator) is
  the natural follow-up.
- Action units; the overall mass scale is electron-calibrated at each λ,
  so only ratios are meaningful.

## Reproduce

```bash
python -m experiments.closure_ledger.odd_k_ladder_spectral_deformation_probe
# Verdict: ODD_K_LADDER_SPECTRUM_DEFORMS_SMOOTHLY_OVER_A_FINITE_BERGER_WINDOW
#          _BUT_THE_MU_E_HIERARCHY_IS_METRIC_FINE_TUNED
```
