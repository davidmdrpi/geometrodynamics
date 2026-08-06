# Probe L — the minimal missing interaction, tested directly (PR #239)

_Run 2026-08-06T06:12:03.500976+00:00 · 8/8 PASS_

**Q. #238 said the Bell chain needs a winding-sector-mixing operation. Can it be had, and what does it cost?**

**A. Yes — uniquely, and it costs charge conservation.**

## It is unique at first order

| harmonic `m` | `\|⟨+1\|V\|−1⟩\|` |
|---:|---:|
| 1 | 0.0000 |
| 2 | 0.4082 |
| 3 | 0.0000 |
| 4 | 0.0000 |

Restricted to the qubit it is **pure `σ_x`**, coefficient **0.408248** — `I`, `σ_y`, `σ_z` coefficients all zero to 1e-16. Exactly the generator #238 found missing.

## It works — with the *derived* detector

| CHSH with the minimal mixing | **2.828427** |
|---|---:|
| CHSH with no mixing (#238) | 2.000000 |

## Cost 1 — charge

`‖[H, K]‖ = 6.6e-16` for the committed dynamics, but `‖[V₂, K]‖ = 1.414`. Winding **is** charge (#42–#44), and the fiber U(1) that protects it is exactly the symmetry the mixing term must break.

## Cost 2 — leakage, and the loophole-free window

| `\|t\|` span | CHSH `S` | `η` | `4/η − 2` | margin |
|---:|---:|---:|---:|---:|
| 0.05 | 2.002219 | 0.9995 | 2.0019 | +0.000350 |
| 0.1 | 2.008855 | 0.9981 | 2.0075 | +0.001370 |
| 0.2 | 2.035079 | 0.9925 | 2.0300 | +0.005031 |
| 0.3 | 2.077659 | 0.9833 | 2.0680 | +0.009639 |
| 0.4 | 2.134938 | 0.9704 | 2.1220 | +0.012973 |
| 0.5 | 2.204661 | 0.9540 | 2.1927 | +0.011964 |
| 0.7 | 2.369812 | 0.9113 | 2.3892 | -0.019371 |
| 0.9 | 2.545816 | 0.8565 | 2.6702 | -0.124407 |
| 1.2 | 2.760876 | 0.7553 | 3.2958 | -0.534970 |
| 1.5 | 2.827939 | 0.6378 | 4.2717 | -1.443798 |
| 1.9 | 2.828028 | 0.4697 | 6.5162 | -3.688143 |

A loophole-free window **exists** — best margin **+0.0130** at `S = 2.1349`, `η = 0.9704` — but **Tsirelson saturation is not in it**.

## Verdict

**THE_MINIMAL_MISSING_INTERACTION_IS_THE_SECOND_CHI_HARMONIC_AT_A_MOUTH_WHICH_IS_EXACTLY_PURE_SIGMA_X_ON_THE_QUBIT_AND_RESTORES_CHSH_TO_2_83_BUT_IT_BREAKS_WINDING_CHARGE_CONSERVATION_AND_ONLY_S_2_13_SURVIVES_THE_LEAKAGE_LOOPHOLE_FREE**
