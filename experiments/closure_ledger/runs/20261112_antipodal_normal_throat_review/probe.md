# Four-cell review: circulation and transfer are different constraints

Public amendment: `31cfc67e5ef92286c3893f2a8667845cf59b544d`. Pattern already supplied by review.

Verdict: **PASS**.

Scale is the initial squared field-phase-space norm (0.16 for generic runs).
Each row covers eight draws; face change is measured relative to its initial value.

| N | n | cell | max total L / scale | min face change / scale | max face change / scale | max net flux / scale |
|---|---|---|---:|---:|---:|---:|
| 16 | 1 | shift_kept | 0.19212 | 0.076958 | 0.2635 | 0.092601 |
| 16 | 1 | antipodal | 0.15351 | 5.5972e-16 | 1.7781e-15 | 3.3153e-16 |
| 16 | 1 | reflective | 1.8093e-15 | 3.6185e-16 | 1.4034e-15 | 4.3554e-16 |
| 16 | 1 | reflect_reversed | 2.6563e-15 | 0.058403 | 0.21283 | 0.10329 |
| 24 | 2 | shift_kept | 0.26451 | 0.094527 | 0.34981 | 0.12382 |
| 24 | 2 | antipodal | 0.13928 | 1.6155e-15 | 4.1335e-15 | 9.297e-16 |
| 24 | 2 | reflective | 4.9873e-15 | 1.7781e-15 | 4.391e-15 | 1.1357e-15 |
| 24 | 2 | reflect_reversed | 6.3751e-15 | 0.24298 | 0.7501 | 0.15252 |

| Check | Pass |
|---|---|
| spectral_selection | True |
| angular_symmetry | True |
| seam | True |
| momentum | True |
| mouth_pair | True |
| response_scaling | True |
| independent_integrator | True |
| coordinate_geometry | True |
| one_face_susceptibility | True |
| phi_symbol_convergence | True |
| original_runs_tightened_energy | True |
| generic_coverage | True |
| generic_multimode_and_moving_opening | True |
| generic_seam_and_momentum | True |
| generic_energy | True |
| four_cell_circulation | True |
| four_cell_transfer | True |
| strong_response | True |
| strong_independent_integrator | True |
| validation_negative_controls | True |

Strong response (one independent opening coordinate):
```json
{
  "initial_total": 0.37903469074266727,
  "max_total": 0.37903469074266843,
  "momentum_drift": 1.1862974370951239e-15,
  "face_change": 5.931487185475619e-16,
  "face_max": 0.18951734537133422,
  "net_flux": 0.0,
  "seam_error": 3.129846246676283e-15,
  "energy_excursion": 1.0423669469960979e-05,
  "outside_strongest_phi_bin": 0.0,
  "initial_mouth_momentum": 0.0,
  "initial_mouth_position": 0.0,
  "max_opening": 1.064924436958468,
  "amplitude": 2.5,
  "duration": 40.0,
  "relative_momentum_drift": 3.1297859168793608e-15,
  "coarse_final_error": 0.0014707198365964175,
  "fine_final_error": 0.000367828123318142,
  "improvement": 3.9983887673655722,
  "independent_opening_coordinates": 1
}
```

The orientation/transfer correspondence holds for this four-cell Hamiltonian family.
It is not a general non-orientability obstruction.
