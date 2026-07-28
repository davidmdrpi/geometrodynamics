# What selects the twist sector? (PR #230)

_Run 2026-07-28T01:18:14.934093+00:00 · 5/5 PASS_

## The twist costs energy — and the sign is the statistics

| C | statistics | `E0` untwisted | `E0` twisted | `ΔE` | favoured |
|---:|---|---:|---:|---:|---|
| 1.0 | bosonic | -0.523599 | +0.261799 | +0.785398 | untwisted |
| 1.0 | fermionic | +0.523599 | -0.261799 | -0.785398 | twisted |
| 2.0 | bosonic | -0.261799 | +0.130900 | +0.392699 | untwisted |
| 2.0 | fermionic | +0.261799 | -0.130900 | -0.392699 | twisted |
| 3.5 | bosonic | -0.149600 | +0.074800 | +0.224399 | untwisted |
| 3.5 | fermionic | +0.149600 | -0.074800 | -0.224399 | twisted |

## But it cannot relax: the gate is an amplitude zero

| s | closing amplitude | `b₁` | holonomy |
|---:|---:|---:|---|
| 0.0 | +1.00 | 1 | +1 |
| 0.25 | +0.50 | 1 | +1 |
| 0.5 | +0.00 | 0 | UNDEFINED |
| 0.75 | -0.50 | 1 | -1 |
| 1.0 | -1.00 | 1 | -1 |

## The holonomy is deformation-invariant

  - untwisted: lowest eigenvalue ≤ 1.6e-15 (machine zero) over 200 random re-weightings
  - twisted: lowest eigenvalue ≥ 4.82e-02 — the gap never closes

## Both Möbius sectors, from one Casimir sign

  - **QCD Moebius (#100/#103)** — flux-tube transverse displacement (phonon) (bosonic) ⟹ untwisted favoured. presents the Moebius states as an extra tower above the orientable one: +pi*sigma for the glueball tower (#100), +2 sqrt(sigma) = 849 MeV for the heavy baryons (#103).
  - **BAM matter (#170/#195/#202)** — throat modes are Pin- spinors (fermionic) ⟹ twisted favoured. #202's Pin parity selects the electron k = 1 mode ODD at the neck, and #227's T6 finds constructive antipodal self-return in the twisted odd-l sector.

## Verdict

**THE_TWIST_SECTOR_IS_SELECTED_ENERGETICALLY_WITH_A_SIGN_SET_BY_THE_LINK_FIELD_STATISTICS_AND_THEN_TOPOLOGICALLY_FROZEN_AT_NUCLEATION**
