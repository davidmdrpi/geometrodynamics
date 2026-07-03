# The 5D throat-core solve - companion probe (PR #202)

**Run:** 2026-07-03T03:08:27+00:00

The deliverable is `docs/throat_core_5d_solve.md` - the exact suppression law on the J-quotiented Tangherlini bridge, deriving #201's fitted exponent. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the register item: the 5D core solve | **PASS** |
| T2 | `T2_bridge_and_twisted_bc` | the bridge; odd-k Dirichlet at the cross-cap (quotient equality) | **PASS** |
| T3 | `T3_exact_radial_law` | the exact law: phi = sigma (k=1); exponents = k; tail theorem | **PASS** |
| T4 | `T4_suppression_law` | eps_k = (rs/sigma)^k e^{c0}; c0 = {0, -0.405, -0.783} | **PASS** |
| T5 | `T5_naturalness_refined` | physical sensitivity = -k = -1 exactly (74.7 -> 4.48 -> 1) | **PASS** |
| T6 | `T6_inversion_and_impossibility` | inversion exact: sigma/rs = 88.6; eps3/eps1 ~ 8e-5 | **PASS** |
| T7 | `T7_zero_winding_contrast_and_scope` | only the zero-winding channel touches the cross-cap; scope | **PASS** |
| T8 | `T8_assessment` | one number from an m_e/m_mu prediction, with an exact law | **PASS** |

## The law

- slopes d S/d ln(sigma): {1: 1.0, 3: 2.99962, 5: 4.99925} (= k)
- matching constants c0(k): {1: 0.0, 3: -0.4054, 5: -0.7827} (c0(1) = 0 exact)
- physical sensitivities: {1: 1.0, 3: 2.9999, 5: 4.9997} (electron = 1)
- inversion: sigma_mode/rs = {'A': 88.61, 'B': 206.77} ; eps3/eps1 = 1.91e-04 (needed for pairing hierarchy: 88.6)

## Verdict

**FIVE_D_THROAT_CORE_SOLVED_ODD_K_NODE_AT_THE_CROSS_CAP_SUPPRESSION_IS_A_POWER_LAW_IN_THE_SCALE_RATIO_ELECTRON_SENSITIVITY_ONE.** THE FITTED EXPONENT IS NOW A LAW (the argument is in docs/throat_core_5d_solve.md; this probe verifies every step).

THE GEOMETRY. On the t=0 slice of the J-quotiented Tangherlini throat the k-winding problem is 1D on the bridge rho = sqrt(rs^2 + sigma^2), and the deck parity (-1)^k forces the Pin-twisted boundary condition: ODD-k MODES HAVE A NODE AT THE CROSS-CAP (quotient-spectrum equality to 1e-5) - the geometric realization of #195's forbidden one-mouth lift; even-k modes are Neumann (the cross-cap repels fermions, admits bosons - contrast measured).

THE LAW. At E~0 (the physical regime; tail theorem verified at 2%) the regular solution is exactly phi = sigma for k=1 and sigma^k in general: the suppression is eps_k = (rs/sigma_mode)^k e^(c0) with slopes {1: 1.0, 3: 2.99962, 5: 4.99925} and computed constants c0 = {1: 0.0, 3: -0.4054, 5: -0.7827} - #201's e^(-kc) DERIVED, with c = ln(sigma_mode/rs).

THE CONSEQUENCES. The electron mass is a LINEAR readout of the throat-to-mode hierarchy: physical-parameter sensitivity {1: 1.0, 3: 2.9999, 5: 4.9997} - the naturalness chain closes 74.7 -> 4.48 -> 1. The inversion is exact (c0(1)=0): sigma_mode/rs = 88.61 (conv A) - the throat core is ~1% of the wave extent; ONE dimensionless number, governed by an exact law, now separates the program from an outright m_e/m_mu prediction (the coupled 5D+soliton solve). And the impossibility bound hardens: eps3/eps1 = 2e-04 vs the 88.6 a pairing hierarchy would need.
