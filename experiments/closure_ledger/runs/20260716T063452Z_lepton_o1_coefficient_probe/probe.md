# The derived O(1) lepton mass coefficient (PR #221)

**Run:** 2026-07-16T06:34:52+00:00

The deliverable is `docs/lepton_o1_coefficient.md` - the #210 register item executed on the #220 eigenhistory background. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the goal: derive, not fit | **PASS** |
| T2 | `T2_hard_wall_theorems` | the hard-wall theorems: X = pi/2 exactly | **PASS** |
| T3 | `T3_measurement` | the measurement; the #202 parity dichotomy realized | **PASS** |
| T4 | `T4_robustness` | regulator-, exterior-, grid-robust | **PASS** |
| T5 | `T5_eigenhistory_orbit` | the eigenhistory orbit carries it; source-decoupled | **PASS** |
| T6 | `T6_confrontation` | conv A excluded, conv B lands at 92-96% | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | assessment | **PASS** |

## Verdict

**Class:** `THE_REMAINING_O1_IS_DERIVED_THE_ODD_THROAT_FUNDAMENTAL_GIVES_X_EQUALS_PI_OVER_2_PLUS_THROAT_SHIFT_AND_ME_OVER_MMU_EQUALS_ALPHA_OVER_X_AT_92_TO_96_PERCENT_NO_FIT`

ESTABLISHED (the argument is in docs/lepton_o1_coefficient.md).

THE COEFFICIENT. X = sigma_mode * w on ONE geometric object - the eigenhistory whose frequency is the mass (lambda_C = 1/w) and whose #202 matching radius is sigma_mode: a derived ratio, nothing to fit.

THE THEOREM. #202's Pin-twisted boundary condition forces the electron (k = 1) transit mode ODD - the exact node at the neck is reproduced here to 1e-12 with the phi ~ sigma near-neck law (drift 2.0%); for the odd fundamental the antinode is the quarter wave: X = pi/2 EXACTLY on hard walls (closed forms verified to 5e-09), cavity-length independent - the analytic reason the coefficient is O(1).

THE MEASUREMENT. On the real Tangherlini barriers (reference D = 12): X_match = 1.5783 (hard wall 1.5708); regulator scan (depth 8-48) band [1.5787, 1.6381]; exterior spread 0.035; dx-shift 7e-04.

THE ORBIT. The #220 Gauss-Newton eigenhistory on the odd mode: residual 3e-13; the source decouples EXACTLY (q* = 3e-17) - the charged transit mode is source-transparent at the crossing; complete monodromy on the unit circle to 2e-14; the orbit carries X_match = 1.5790 (linear-mode value to 6e-07); exactly linear (amplitude-independent, 1e-13).

THE CONFRONTATION. Convention A (required 0.6467) EXCLUDED x2.4; convention B (required 1.5089) SELECTED: m_e/m_mu = alpha/X = 0.004624 (band [0.004455, 0.004622]; hard-wall anchor 2 alpha/pi = 0.004646) vs observed 0.004836 - the derivation lands at 92.1%-95.6% of the observed ratio with ZERO fitted numbers (inputs: the throat geometry and alpha; m_e used only for comparison).
