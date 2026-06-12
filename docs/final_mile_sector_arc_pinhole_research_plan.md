# The final geometric mile: the Hopf sector arc + the pinhole refinement (PR #160)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. Part A is
> fiber algebra; Part B is a mass-preserving refinement audit.

The two items flagged at the close of #159, taken together.

## Part A: the sector arc is the Weyl quantum

A capacity-k₅ winding space (k ∈ {1,3,5} capped at k₅ = 5, #73/#126) carries
the canonical clock–shift (Weyl) pair, and on **any** k₅-dimensional space

    U V U† V† = e^{2πi/k₅} · 1     (machine-exact, verified)

— the fiber position dual to the winding label is a k₅-site lattice
θ_n = 2πn/k₅, and winding-changing transport (the shell hop) is the shift V,
whose minimal step is **one site: arc 2π/k₅**. The #159 "identification" is
the Weyl commutator quantum of the capacity-bounded fiber — algebra, not
choice, with no radial-profile model anywhere. Composed with the
connection's ½: **φ_h = π/k₅, algebraic end-to-end** (re-verified by
transport integration; the k = 5 case compared on the unit circle to avoid
the ±π branch point).

## Part B: the pinhole refinement endgame

- **Exclusion 1 — pinhole single-knob**: pinhole* = 20.77 lands V_us = 0.225
  but shifts m_s by **−22.5%** (≫ the 1.6% accuracy). V_us and m_s ride the
  same d–s direction.
- **Exclusion 2 — transport rescale**: the exact (t, α) map that doubles the
  dk = 3 element while fixing dk = 5 raises V_us only to 0.133 with m_s
  +50% — level repulsion self-defeats (the 2×2 invariant
  `sin 2θ = 2|H_ds|/Δλ`, verified).
- **The mass-preserving family**: eigenvector rotations at exactly fixed
  eigenvalues (masses to 1e-15). The joint (V_us = 0.225, β = 22.2°)
  solution at (δθ_u, δθ_d) = (−5.2°, +9.9°) lands **seven of eight**
  observables:

  | observable | refined/observed |
  |---|---|
  | V_us | ×1.00 (exact) |
  | V_cb / V_ts | ×0.90 / ×0.89 |
  | V_ub | ×1.10 |
  | V_td | ×1.19 |
  | J | **×1.05** |
  | β | 22.2° (exact-fit) |
  | γ | **104° vs 65.9° — the single misfit** |

- **The J-ceiling lock, verified**: #156/#158 predicted the Jarlskog ceiling
  rises to the observed 3.5×10⁻⁵ when the soft elements land. At the refined
  point: ceiling = 0.99 of observed, J = ×1.05 — a prediction made about a
  state that did not yet exist, now checked in that state and **passed**.
- **The re-lock targets**: the rotated minus-block elements (H_ds ×1.84-class
  with %-level diagonal compensation) are tabulated as the precise targets
  for the next v3+CP joint re-lock.

## The residual after the final mile

The flavor sector's remaining residual is **one number — the γ angle (~38°
misfit)** — plus the knob-level realization of the re-lock targets. No new
inputs consumed; the #150 budget unchanged.

## Reproduce

```bash
python -m experiments.closure_ledger.final_mile_sector_arc_pinhole_probe
# Verdict: SECTOR_ARC_WEYL_DERIVED_SOFT_DIRECTION_REDUCED_TO_GAMMA_CEILING_VERIFIED
```
