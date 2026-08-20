# Areal probe — the signed ΔA/A on the resolved neck

**Question.** does the interaction metric deform toward a neck, away from one, or merely oscillate?

**Answer.** toward a neck: dA/A < 0 at both mouths — ΔA/A = (-2.0575e-03, -1.8818e-03) in units of 2 pi G, with the source being PR #263's interference dT00.

Throat: area 12.5664, length 0.9, wavenumber 1.0000 — phase 0.90, inside the first cavity cell.

**8/8 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | the constraint and the areal response | PASS |
| T2 | the degeneracy, and both sectors of it | PASS |
| T3 | the projector is the Green function's own tail | PASS |
| T4 | the assembly against exact solves, both sectors | PASS |
| T5 | the dipole layers: required, then invisible | PASS |
| T6 | what the answer is made of | PASS |
| T7 | *** the signed areal response *** | PASS |
| T8 | the sign is a statement about a throat | PASS |

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

Sign agrees in every variant: True. Quadrature spread at fixed radius: 2.13%. Worst condition number: 2.12e+05, worst residual: 2.8e-15.

**What it means.** toward a neck — the conformal factor falls at both mouths, so both mouth areas contract.  The interference energy alone would open them (U(c_j) > 0 at both); the throat's monopole layers overshoot that and invert it.
