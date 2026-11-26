# Alice-handle exchange: V-curve, aperture scaling and envelope control (freeze `291ab5f`)

**Verdict: COUNTEREXAMPLE** (G4)

Largest |u + v − 1|, sin², N = 2000, K = 150. Each cell is measured (predicted).

| cell | η | 3.21 | 7.05 | 14.75 | 31.42 | 63.49 | 127.61 | 255.87 |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.5 | 8.40e-03 (8.0e-03) | 7.25e-04 (6.8e-04) | 2.05e-04 (1.2e-04) | 1.66e-04 (1.1e-04) | 2.18e-04 (2.1e-04) | 3.82e-04 (4.6e-04) | 8.20e-04 (1.0e-03) |
| 0.5 | 0.9 | η unreachable (max 0.873) | 1.85e-03 (1.8e-03) | 5.94e-04 (4.8e-04) | 6.95e-04 (6.6e-04) | 1.25e-03 (1.3e-03) | 2.52e-03 (2.7e-03) | 5.34e-03 (5.8e-03) |
| 0.35 | 0.5 | 8.26e-03 (8.0e-03) | 5.93e-04 (6.5e-04) | 6.96e-05 (7.2e-05) | 1.45e-05 (1.2e-05) | 1.14e-05 (1.1e-05) | 1.76e-05 (2.1e-05) | 3.32e-05 (4.2e-05) |
| 0.35 | 0.9 | η unreachable (max 0.874) | 1.51e-03 (1.7e-03) | 1.82e-04 (2.0e-04) | 5.09e-05 (5.2e-05) | 6.49e-05 (6.9e-05) | 1.23e-04 (1.3e-04) | 2.45e-04 (2.7e-04) |

| gate | result |
|---|---|
| G1 predictor within ×3 | pass (measured/predicted 0.79–1.76) |
| G2 V shape | pass |
| G3 aperture scaling | pass (η=0.5: 1.66e-04 → 1.14e-05, ×14.6; η=0.9: 5.94e-04 → 5.09e-05, ×11.7) |
| G4 spectral nulls | FAIL (6.413: min at 6.4, depth 0.29; 7.695: min at 7.7, depth 0.31; 8.978: min at 9.0, depth 0.40; 10.260: min at 10.3, depth 0.50) |
| G5 envelope | pass (τ=8 sin4/sin² 0.136 vs 0.117; τ=8 bump/sin² 1.760 vs 1.225; τ=16 sin4/sin² 0.151 vs 0.061; τ=16 bump/sin² 3.399 vs 3.631) |
| G6 numerics | pass (K True, N True, balance 5.1e-12, undetected 1.7e-11) |

Run time 1259 s.
