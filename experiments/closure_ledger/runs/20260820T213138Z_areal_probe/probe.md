# Areal probe — the signed ΔA/A on the resolved neck

**Question.** does the interaction metric deform toward a neck, away from one, or merely oscillate?

**Answer.** toward a neck at the wide working throat: dA/A < 0 at both mouths -- but the sign is a property of the throat, and matching the tube's area to the mouths' flips it

ΔA/A = (-2.0575e-03, -1.8818e-03) in units of 2 pi G, with the source being PR #263's interference dT00.

Throat: area 12.5664, length 0.9, wavenumber 1.0000 — phase 0.90, inside the first cavity cell.

**9/9 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | the constraint and the areal response | PASS |
| T2 | the degeneracy, and both sectors of it | PASS |
| T3 | the projector is the Green function's own tail | PASS |
| T4 | the assembly against exact solves, both sectors | PASS |
| T5 | the dipole layers: required, then invisible | PASS |
| T6 | what the answer is made of | PASS |
| T7 | *** the signed areal response *** | PASS |
| T8 | the throat is a resonant cavity | PASS |
| T9 | *** the sign does not survive a matched tube *** | PASS |

## The controls

| radius | points | gluing | ΔA/A mouth 1 | ΔA/A mouth 2 |
|--------|--------|--------|--------------|--------------|
| 0.05 | 5158 | transported | -2.1014e-03 | -1.9148e-03 |
| 0.05 | 5158 | reflected | -2.1014e-03 | -1.9148e-03 |
| 0.05 | 12630 | transported | -2.0575e-03 | -1.8818e-03 |
| 0.05 | 12630 | reflected | -2.0575e-03 | -1.8818e-03 |
| 0.10 | 5158 | transported | -1.2247e-03 | -1.0899e-03 |
| 0.10 | 5158 | reflected | -1.2247e-03 | -1.0899e-03 |
| 0.10 | 12630 | transported | -1.1783e-03 | -1.0547e-03 |
| 0.10 | 12630 | reflected | -1.1783e-03 | -1.0547e-03 |

Sign agrees in every variant: True. Quadrature spread at fixed radius: 2.13%. Worst condition number: 1.45e+04, worst residual: 4.3e-16.

**What it means.** toward a neck AT THIS THROAT — the conformal factor falls at both mouths, so both mouth areas contract.  The interference energy alone would open them (U(c_j) > 0 at both); the throat's monopole layers overshoot that and invert it.  Matching the tube's area to the mouths' flips the sign, so this is a property of the throat and not of the source.
