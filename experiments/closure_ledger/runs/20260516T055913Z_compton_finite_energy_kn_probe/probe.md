# Finite-energy Klein-Nishina recoil probe

**Run:** 2026-05-16T05:59:13+00:00

Follow-on to PR #28 (Thomson-limit KN restoration). Tests whether the natural BAM construction — antipodal `S³` Green-function propagator + photon transverse polarization + scalar electron — also reproduces Klein-Nishina at finite `ω/m_e`, or whether the Thomson match conceals a leading-order failure.

**Construction:**

```
BAM: |M|² ∝ (1+cos²θ) · (1/ψ_s + 1/ψ_u)²  with ψ_x = (x−m²)/(2m²). KN: |M|² ∝ x²·(x + 1/x − sin²θ) with x = ω'/ω.
```

## Test summary

| # | Test | Key metric | Value | PASS? |
|---|---|---|---:|---|
| T1 | `T1_thomson_sanity` | max |f_BAM − f_KN| at ω/m=1e-4 | 6.00e-04 | **PASS** |
| T2 | `T2_leading_order_discrepancy` | fitted leading order in ε (predicted 1.0) | 0.999, 0.994, 0.991 | **PASS** |
| T3 | `T3_forward_backward_asymmetry` | max |A_BAM(ε) − A_KN(ε)| | 1.6059 | **PASS** |
| T4 | `T4_compton_edge` | ε at which BAM/KN diverge by 10 % | 0.0178 | **PASS** |
| T5 | `T5_angular_fit_at_finite_omega` | Δc0, Δc1, Δc2, Δc3 (BAM − KN) at ω/m=0.5 | +0.697, -0.608, +0.683, -0.761 | **PASS** |

## T1_thomson_sanity

At ω/m = 1e-4 the BAM and KN normalised angular distributions agree across θ ∈ [0, π].

## T2_leading_order_discrepancy

Lowest non-vanishing order in ε = ω/m of f_BAM(ε, θ) − f_KN(ε, θ) at generic θ. Predicted fitted slope = 1.0 (O(ε)).

Per-θ fits:

| θ/π | fitted order in ε |
|---:|---:|
| 0.2500 | 0.9992 |
| 0.5000 | 0.9941 |
| 0.7500 | 0.9909 |

## T3_forward_backward_asymmetry

Forward-backward asymmetry A(ε) = (f(0)−f(π))/(f(0)+f(π)) compared between BAM and KN across ε ∈ [0.01, 1].

| ε | f_BAM(π) | f_KN(π) | A_BAM | A_KN | ΔA |
|---:|---:|---:|---:|---:|---:|
| 0.010 | 1.0202 | 0.9614 | -0.0100 | +0.0197 | -0.0297 |
| 0.050 | 1.1052 | 0.8302 | -0.0500 | +0.0928 | -0.1428 |
| 0.100 | 1.2225 | 0.7060 | -0.1001 | +0.1723 | -0.2724 |
| 0.200 | 1.5034 | 0.5394 | -0.2011 | +0.2992 | -0.5003 |
| 0.500 | 3.0568 | 0.3125 | -0.5070 | +0.5238 | -1.0308 |
| 1.000 | 23.5061 | 0.1852 | -0.9184 | +0.6875 | -1.6059 |

## T4_compton_edge

Compton backscatter (θ = π) comparison vs ε on log scale. Identifies the ε at which BAM and KN diverge by > 10 %.

Smallest ε with > 10 % BAM/KN divergence: **0.01778279410038923**

| ε | f_BAM(π) | f_KN(π) | rel. diff |
|---:|---:|---:|---:|
| 1.0000e-04 | 1.0002e+00 | 9.9960e-01 | 0.0006 |
| 3.1623e-04 | 1.0006e+00 | 9.9874e-01 | 0.0019 |
| 1.0000e-03 | 1.0020e+00 | 9.9601e-01 | 0.0060 |
| 3.1623e-03 | 1.0063e+00 | 9.8749e-01 | 0.0191 |
| 1.0000e-02 | 1.0202e+00 | 9.6136e-01 | 0.0612 |
| 3.1623e-02 | 1.0653e+00 | 8.8624e-01 | 0.2020 |
| 1.0000e-01 | 1.2225e+00 | 7.0602e-01 | 0.7315 |
| 3.1623e-01 | 1.9402e+00 | 4.2122e-01 | 3.6060 |
| 1.0000e+00 | 2.3506e+01 | 1.8519e-01 | 125.9328 |
| 3.1623e+00 | 1.6746e+00 | 6.9536e-02 | 23.0822 |
| 1.0000e+01 | 1.1152e+00 | 2.3864e-02 | 45.7330 |

## T5_angular_fit_at_finite_omega

Fit |M|²(θ) to {1, cos θ, cos²θ, cos³θ} at ω/m = 0.5 for both BAM and KN. Compare coefficients.

BAM coefficients (c0, c1, c2, c3): +0.9408, -0.5141, +1.0799, -0.5143
KN coefficients  (c0, c1, c2, c3): +0.2442, +0.0941, +0.3964, +0.2462

Max pointwise |M_BAM − M_KN| at ω/m=0.5: **2.7443**

## Verdict

**PARTIAL_MATCH_O1_THOMSON_ONLY.** PARTIAL MATCH — KN reproduced at O(1) in ε = ω/m_e (Thomson limit, PR #28), but the natural BAM construction fails to reproduce KN at leading order O(ε). The discrepancy f_BAM − f_KN scales linearly in ω/m_e (T2 fitted slopes near 1.0). The forward-backward asymmetry (T3) and the Compton edge at θ = π (T4) show the structural signature: BAM's per-channel propagator sum (G_s + G_u)² gives a (1 + 1/x)² energy dependence, while QED's vertex-factor algebra gives x² · (x + 1/x − sin²θ). The missing BAM ingredient is the per-channel kinematic weighting — the QED vertex factors that contract photon polarization with momentum (ε·k structures), which the present BAM construction does not have natively. This locates the next structural piece: a vertex-coupling derivation from the Hopf connection or throat-transport algebra.

## What this leaves open

- **Vertex-factor BAM derivation.** The probe locates the finite-ω gap in the missing per-channel kinematic weighting (ε·k vertex structure in QED). Whether a BAM-derived vertex factor — from explicit Hopf-connection coupling, throat-transport algebra, or another natural construction — closes the gap is the natural follow-on probe.
- **Photon as Hopf-fibre excitation.** The current photon model treats polarization as flat-space tangent-plane vectors. An `S³`-native photon (Hopf-fibre excitation, `hopf/connection.py`) would automatically carry the connection coupling that supplies ε·k vertex factors. Sketching this is the natural design step.
- **Electron spin at finite energy.** The Thomson sanity used a scalar electron; finite-ω Compton in QED has spin-½ corrections that BAM's previous probe identified as spurious at Thomson but should reappear at higher orders.
