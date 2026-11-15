# Normal reversal and circulating real waves

Public freeze: `533f9257fcb28f015740b9f08557959bc3c1a83b`.

Numerical validity: **PASS**.

| Check | Pass |
|---|---|
| geometry | True |
| spectral_selection | True |
| continuum_convergence | True |
| angular_symmetry | True |
| susceptibility | True |
| seam | True |
| momentum | True |
| mouth_pair | True |
| energy | True |
| response_scaling | True |
| independent_integrator | True |

| N | n | m | amplitude | normal parity | max opening | relative L error | energy excursion |
|---|---|---|---|---|---|---|---|
| 16 | 1 | 0 | 0.05 | -1 | 0.000749715 | 6.05e-15 | 1.83e-08 |
| 16 | 1 | 0 | 0.1 | -1 | 0.00299584 | 2.87e-15 | 5.67e-08 |
| 16 | 1 | 0 | 0.2 | -1 | 0.0119323 | 2.07e-15 | 2.23e-07 |
| 16 | 1 | 0 | 0.4 | -1 | 0.0468293 | 3.5e-15 | 8.54e-07 |
| 16 | 1 | 1 | 0.05 | 1 | 0.000999395 | 1.85e-15 | 2.61e-08 |
| 16 | 1 | 1 | 0.1 | 1 | 0.00398932 | 2.53e-15 | 1.03e-07 |
| 16 | 1 | 1 | 0.2 | 1 | 0.015838 | 1.52e-15 | 4.05e-07 |
| 16 | 1 | 1 | 0.4 | 1 | 0.0614643 | 3.88e-15 | 1.52e-06 |
| 16 | 2 | 0 | 0.05 | 1 | 0.000749706 | 2.08e-15 | 2.17e-08 |
| 16 | 2 | 0 | 0.1 | 1 | 0.00299558 | 2.91e-15 | 5.56e-08 |
| 16 | 2 | 0 | 0.2 | 1 | 0.0119293 | 5.69e-15 | 2.09e-07 |
| 16 | 2 | 0 | 0.4 | 1 | 0.0468605 | 1.8e-15 | 8.11e-07 |
| 16 | 2 | 1 | 0.05 | -1 | 0.000999328 | 4.11e-15 | 4.38e-08 |
| 16 | 2 | 1 | 0.1 | -1 | 0.00398836 | 1.53e-15 | 1.03e-07 |
| 16 | 2 | 1 | 0.2 | -1 | 0.0158136 | 2.7e-15 | 3.52e-07 |
| 16 | 2 | 1 | 0.4 | -1 | 0.0611418 | 4.11e-15 | 1.31e-06 |
| 24 | 1 | 0 | 0.05 | -1 | 0.000749727 | 2.07e-15 | 1.83e-08 |
| 24 | 1 | 0 | 0.1 | -1 | 0.00299589 | 1.91e-15 | 5.68e-08 |
| 24 | 1 | 0 | 0.2 | -1 | 0.0119314 | 1.43e-15 | 2.24e-07 |
| 24 | 1 | 0 | 0.4 | -1 | 0.0467982 | 7.96e-16 | 8.54e-07 |
| 24 | 1 | 1 | 0.05 | 1 | 0.000999347 | 2.36e-15 | 2.61e-08 |
| 24 | 1 | 1 | 0.1 | 1 | 0.00398994 | 1.35e-15 | 1.03e-07 |
| 24 | 1 | 1 | 0.2 | 1 | 0.0158404 | 3.87e-15 | 4.05e-07 |
| 24 | 1 | 1 | 0.4 | 1 | 0.0614693 | 3.03e-15 | 1.52e-06 |
| 24 | 2 | 0 | 0.05 | 1 | 0.000749692 | 2.64e-15 | 2.17e-08 |
| 24 | 2 | 0 | 0.1 | 1 | 0.00299554 | 1.39e-15 | 5.56e-08 |
| 24 | 2 | 0 | 0.2 | 1 | 0.0119287 | 5.69e-15 | 2.09e-07 |
| 24 | 2 | 0 | 0.4 | 1 | 0.046846 | 4.16e-15 | 8.11e-07 |
| 24 | 2 | 1 | 0.05 | -1 | 0.000999337 | 2.11e-15 | 4.38e-08 |
| 24 | 2 | 1 | 0.1 | -1 | 0.00398835 | 2.35e-15 | 1.03e-07 |
| 24 | 2 | 1 | 0.2 | -1 | 0.0158175 | 2.23e-15 | 3.52e-07 |
| 24 | 2 | 1 | 0.4 | -1 | 0.0612038 | 2.82e-15 | 1.31e-06 |

Angular momentum is that of one real harmonic. Signed loads are waist
loads, not a demonstrated bulk momentum balance or electric charge.

Independent integrator:
```json
{
  "duration": 40.0,
  "coarse_final_error": 0.0016859669517852122,
  "fine_final_error": 0.00042147072473773696,
  "improvement": 4.000199427455647,
  "dop853_momentum_error": 5.330992904385433e-10
}
```

Post-freeze transfer exploration:
```json
{
  "scope": "POST_FREEZE_EXPLORATION_NOT_A_FROZEN_GATE",
  "max_local_normal_flux": 0.0007263151854960083,
  "max_net_normal_flux": 9.215718466126788e-19,
  "max_face_momentum_change": 1.3877787807814457e-17,
  "control_initial_net_flux": 0.12000000000000002,
  "control_face_momentum_change": 0.04736196447487004,
  "interpretation": "CIRCULATION_SURVIVES_BUT_EXACT_DECK_SYMMETRY_FORBIDS_NET_NORMAL_ANGULAR_TRANSFER"
}
```
