# Classical closure equilibrium probe

Local preregistration: `76ed50e`. Checks: **19/19**.

GitHub publication after calculation: `d83d46aab8cfe87b0a0adc7f9674401be455d74f`.

Conditional model: round rotor with identical sector inertia; canonical equilibrium preparation; specified classical frame-restoring energy; equal sector priors except P6; geodesic triangle with singlet partner sign.

| criterion | passed | error | bound |
|---|---|---|---|
| P1 quaternion frame identity and invariances | True | 5.273559366969494e-16 | 1e-11 |
| P2 limiting joint law | True | 0.0001906017568821805 | 0.002 |
| P2 scaled partition masses | True | 0.0020911258282508838 | 0.01 |
| P2 independent coordinate refinements | True | 1.6687925680214377e-07 | 0.0001 |
| P2 independent full-sphere Monte Carlo (standard errors) | True | 0.7360391882498649 | 6.0 |
| P2 Monte Carlo correlation (standard errors) | True | 0.5561263171578296 | 6.0 |
| P3 normal residual exact Gaussian | True | 4.107825191113079e-15 | 1e-12 |
| P3 normal residual zero correlation | True | 0.0 | 1e-12 |
| P4 covariance of finite-temperature partition | True | 6.938893903907228e-18 | 1e-12 |
| P4 covariance of pointwise energy | True | 1.1102230246251565e-16 | 1e-12 |
| P4 fixed-stiffness analytic mass | True | 1.7763568394002505e-15 | 1e-10 |
| P4 fixed-stiffness limiting joint law | True | 4.033218242416314e-05 | 0.002 |
| P4 fixed stiffness changes the ensemble | True |  |  |
| P5 quartic limiting joint law | True | 0.00019005175201153302 | 0.002 |
| P5 quartic scaled partition masses | True | 0.002812895436047458 | 0.01 |
| P5 quartic changes finite-temperature response | True |  |  |
| P2 puncture excision asymptotic | True | 0.0008330556051556748 | 0.001 |
| P6 paired half marginals | True | 5.551115123125783e-17 | 1e-12 |
| P6 prior dependence survives | True |  |  |

Full-sphere partition integrals; beta = K/(k_B T).

| gamma | beta | E | max joint error to limit | max scaled mass relative error |
|---|---|---|---|---|
| 0.78539816 | 16 | -0.4105985761 | 3.124e-02 | 2.649e-01 |
| 0.78539816 | 64 | -0.4944317382 | 1.028e-02 | 1.034e-01 |
| 0.78539816 | 256 | -0.5239631784 | 2.902e-03 | 3.109e-02 |
| 0.78539816 | 1024 | -0.5325540359 | 7.542e-04 | 8.236e-03 |
| 0.78539816 | 4096 | -0.5348083838 | 1.906e-04 | 2.091e-03 |
| 1.00000000 | 16 | -0.3086809320 | 2.245e-02 | 1.847e-01 |
| 1.00000000 | 64 | -0.3704913366 | 7.001e-03 | 6.420e-02 |
| 1.00000000 | 256 | -0.3908711388 | 1.906e-03 | 1.810e-02 |
| 1.00000000 | 1024 | -0.3965421444 | 4.886e-04 | 4.687e-03 |
| 1.00000000 | 4096 | -0.3980052277 | 1.229e-04 | 1.183e-03 |
| 2.35619449 | 16 | 0.4105985761 | 3.124e-02 | 2.649e-01 |
| 2.35619449 | 64 | 0.4944317382 | 1.028e-02 | 1.034e-01 |
| 2.35619449 | 256 | 0.5239631784 | 2.902e-03 | 3.109e-02 |
| 2.35619449 | 1024 | 0.5325540359 | 7.542e-04 | 8.236e-03 |
| 2.35619449 | 4096 | 0.5348083838 | 1.906e-04 | 2.091e-03 |

Independent Monte Carlo: E = -0.309181587 ± 0.000900255 (one standard error; n = 120000, seed = 20260905). This integrates an equilibrium law; it does not simulate equilibration.

Standard-angle CHSH: limit 2.1422831632; beta = 4096: 2.1392335351. No global maximum asserted.

Open: BAM origin of the coupling; equilibration and local detector implementation; sector prior; composition and event readout; source-readout causality.

The positive measure follows from the stated energy and preparation. The closure locus alone does not fix it. Detailed controls are in probe.json.
