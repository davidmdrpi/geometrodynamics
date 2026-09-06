# Finite pointer-spread experiment

Specified coupled rotor history laws; no operational BAM field or intervention derivation.

Public freeze: `baf856df06b26de6c9750697ba5df5ad521000c5`.

| sigma P | fixed-marginal contrast | scramble SE | joint-conditioned contrast | frozen-posterior contrast | momentum recoil RMS |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.34297949 | 0.002 | 0.34297949 | 0.34297949 | 6.02517e-16 |
| 0.01 | 0.34242324 | 0.002 | 0.34242330 | 0.34279825 | 0.00816216 |
| 0.025 | 0.33978305 | 0.0018 | 0.33978510 | 0.34185074 | 0.0204052 |
| 0.05 | 0.33096737 | 0.0015 | 0.33098460 | 0.33852179 | 0.040809 |
| 0.1 | 0.29958141 | 0.0012 | 0.29967452 | 0.32603166 | 0.0816067 |
| 0.25 | 0.16903365 | 0.00043 | 0.17317721 | 0.26396515 | 0.203819 |
| 0.5 | 0.02372069 | 0.0008 | 0.05027816 | 0.15976592 | 0.40624 |

Scramble SE is a numerical diagnostic, not a rigorous confidence interval.

| criterion | pass |
|---|---|
| P4 rotation-covariant blind instrument | True |
| P4 constant response erases information | True |
| P1 quaternion product | True |
| P1 stationary triangle bridge | True |
| P1 stationary pointer bridge | True |
| P5 ODE at P=-0.3 | True |
| P5 work at P=-0.3 | True |
| P5 second-order at P=-0.3 | True |
| P5 ODE at P=0.0 | True |
| P5 work at P=0.0 | True |
| P5 ODE at P=0.2 | True |
| P5 work at P=0.2 | True |
| P5 second-order at P=0.2 | True |
| P1 independent sphere probability bridge | True |
| P2 primary contrast exceeds 0.1 in every scramble | True |
| P5 flow and transport invariants | True |
| P4 weight and record parity | True |
| P4 fixed Gaussian pointer moments | True |
| P4 odd means and smoothed signs | True |
| P5 nonzero primary source recoil | True |
| P5 zero-P source motion is free | True |
| resolution primary steps | True |
| resolution primary Hermite | True |
| resolution primary source | True |
| resolution wide Hermite | True |

Passed 25/25 criteria.

The JSON retains all preparation-specific records, posterior P moments, path diagnostics and independent refinements.
