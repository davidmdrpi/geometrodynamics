# TT field-to-triangle-rotor test

Unconstrained linear homogeneous ESU TT mode and its proposed uniaxial cone; not a no-go for every BAM field or driven solution.

Public freeze: `0eb684b86a57be9cfbd2d41d89377cc11f6b76cf`.

Physical verdict: **FREE_ROTATING_UNIAXIAL_TT_FAMILY_NOT_INVARIANT**.

| time | distance to full uniaxial cone | distance / time squared |
|---:|---:|---:|
| 0.04 | 1.81361262086e-06 | 0.00113350788804 |
| 0.02 | 4.52761061856e-07 | 0.00113190265464 |
| 0.01 | 1.13150364829e-07 | 0.00113150364829 |
| 0.005 | 2.82851009989e-08 | 0.00113140403996 |

Predicted coefficient: 0.0011313708499.

| criterion | pass |
|---|---|
| Q1 independent ADM round and linear anchors | True |
| Q1 ADM kinetic and Cartan frequency agree | True |
| Q1 matrix velocity and restricted action | True |
| Q1 finite-difference pullback metric | True |
| Q1 relational rotation covariance and antipodal identity | True |
| Q2 omitted equation matches full tensor residual | True |
| Q2 normal projector matches explicit biaxial basis | True |
| Q2 free normal norm identity | True |
| Q2 restricted radial and angular equations vanish | True |
| Q2 primary full-field trajectory leaves the cone | True |
| Q2 second-order departure in time | True |
| Q2 predicted small-time coefficient | True |
| Q2 minimization includes signed amplitude and all axes | True |
| Q2 stationary-director invariant control | True |
| Q2 departure persists at linear field order | True |
| Q2 every nonzero test speed has nonzero normal residual | True |
| Q2 independent unconstrained DOP853 flow | True |
| Q2 conserved free energy and rotation charge | True |
| Q2 manufactured forced rotation solves full equation | True |
| Q2 manufactured rotation requires radial and normal stress | True |
| Q2 frozen stationary amplitude also needs support | True |
| Q3 repository source normalization and response | True |
| R1 exact full-period spectrum and distance | True |
| R2 periodic returns do not imply interval invariance | True |
| R3 a continuously advancing eigenline exists | True |
| R4 nearest-uniaxial axis switches between distinct optima | True |

Passed 26/26 checks (22 frozen, 4 post-review).

A passing obstruction check is not a successful rotor reduction. Phi selection and operational source readout remain unproved.
