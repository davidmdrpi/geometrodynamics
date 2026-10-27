# Conservative coherent recoil controls

Source: `experiments/closure_ledger/coherent_recoil_probe.py`.
Interpretation and equations: `docs/coherent_mouth_recoil.md`.

- `spectrum.json`: 9,100 normalized all-walk neck-normal pressures, now with a one-sided source. Axis order and the 3D impedance law are explicit.
- `exposure.json`: corresponding externally fed incident exposures; not stored energy gain.
- `spectral_pressure_per_source.json`: pressure per source energy alongside the normalized ratio.
- `source_controls.json`, `uniform_throat_phase.json`: asymmetric-source checks and the old symmetric-source degeneracy.
- `birth_flash_spectra.json`, `birth_weighted_pressure.json`: actual flash energy weights and 192 physical-length controls; incident port projections remain assumptions.
- `stress_projection_controls.json`: zero centred-mode tangent resultant and local charged cross-stress integral q E.
- `review_summary.json`: review-control counts and residuals.
- `closed_paths.json`: 560 finite-seed fixed-geometry cases, each coherent and dephased.
- `pulse_controls.json`, `duration_controls.json`: 30 width/delay and 10 duration checks.
- `clock_controls.json`: 650 stationary boundary checks, including advanced directed delays.
- `finite_mass.json`: 315 individual energy/momentum-conserving recoil completions.
- `geometric_sizes.json`: repository vacuum-bridge ratios under both 2D rim and 3D area impedance laws.
- `summary.json`: counts, residuals and declared-grid sign fractions.
- `sha256.json`: repository-relative source/data/figure hashes.

The centred m=0 mode has zero tangential force by symmetry. Its negative stress
is neck-normal; fixed supports absorb it. Time-domain bins are positive-frequency
carrier envelopes, not real broadband flashes. The finite-mass encounter completion is tested separately;
no moving history, binding, collapse, or shrinking throat is simulated here.
No old trajectory experiment was rerun. All original negative evidence remains.
