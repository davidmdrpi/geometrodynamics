# Source-readout and grounding probe

Public preregistration: `d258bb14e73674fd6ecd7ff6f5d2a46edf20eeab`.

These are conditional laws and added classical couplings. The operational BAM gate remains open.

| beta | choice | mean F | variance F | P(F>0) | P(abs(Y)>0.6), noise=0.15 |
|---|---|---:|---:|---:|---:|
| None | 0 | 3.3923309e-18 | 0.1000000000 | 0.5000000000 | 0.0528578345 |
| None | 1 | 6.7846618e-18 | 0.4000000000 | 0.5000000000 | 0.5136058265 |
| 64.0 | 0 | 1.419474e-18 | 0.1335456636 | 0.5000000000 | 0.1231586349 |
| 64.0 | 1 | -7.3619759e-18 | 0.3904155247 | 0.5000000000 | 0.4897143360 |

The finite-pulse pointer reads q_osc^2, not the spin-frame variable.

```json
{
  "period": 2.299760742237193,
  "source_q0": -0.04614619141493452,
  "pointer_initial_momentum": 0.0,
  "baseline_closure_error": 2.758904216193514e-14,
  "source_endpoint_change": 2.2776919239575477e-13,
  "source_trajectory_change": 6.037725874818989e-12,
  "pointer_shift": 0.0004258902142597511,
  "pointer_momentum_change": 0.0,
  "energy_range": 9.997003225237222e-12,
  "observable": "q_osc^2, not a function of the spin-frame x"
}
```

| criterion | pass |
|---|---|
| A1 mean beta=None choice=0 | True |
| A1 sign beta=None choice=0 | True |
| A1 exact variance choice=0 | True |
| A2 noiseless antiderivative choice=0 | True |
| A1 mean beta=None choice=1 | True |
| A1 sign beta=None choice=1 | True |
| A1 exact variance choice=1 | True |
| A2 noiseless antiderivative choice=1 | True |
| A1 variance gap beta=None | True |
| A2 noisy tail gap beta=None | True |
| A1 mean beta=64.0 choice=0 | True |
| A1 sign beta=64.0 choice=0 | True |
| A1 refine normal choice=0 | True |
| A1 refine azimuth choice=0 | True |
| A1 mean beta=64.0 choice=1 | True |
| A1 sign beta=64.0 choice=1 | True |
| A1 refine normal choice=1 | True |
| A1 refine azimuth choice=1 | True |
| A1 variance gap beta=64.0 | True |
| A2 noisy tail gap beta=64.0 | True |
| A3 periodic source orbit | True |
| A3 zero-P source preserved | True |
| A3 nonzero pointer record | True |
| A3 zero-P energy and momentum | True |
| A3 finite-P recoil | True |
| B1 finite-difference parent | True |
| B1 existing DtN static limit | True |
| B3 Gaussian momentum integral | True |
| B3 quartic momentum integral | True |
| B3 variable inertia volume | True |
| B3 anisotropic unit volume | True |
| B3 nonuniform prior mean | True |
| C quaternion_phase | True |
| C holonomy_trace | True |
| C positive_mass | True |
| C oriented_mass | True |
| C positive_joint | True |
| C oriented_joint | True |
| C phase_gradient | True |
| C frame_hessian | True |
| C source_density | True |
| C source_parity | True |
| C source_normalization | True |
| C sector_prior | True |
| C normal_model_E | True |
| C normal_window_E | True |

Passed 46/46 criteria.

The JSON includes all cross-round residuals, masses, recoil, and field-parent controls.
