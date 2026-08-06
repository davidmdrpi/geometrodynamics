# Probe L — the minimal missing interaction, tested directly (PR #239)

_Run 2026-08-06T06:33:25.399371+00:00 · 8/8 PASS_

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

## Cost 1 — charge: **retracted**

The prescribed term breaks the winding of the throat subsystem (`1.414` against `6.6e-16`) — but that is not the loss of total charge. Include a **winding-2 carrier** and `‖[H_ext, K_total]‖ = 0.0e+00` exactly, with the mean-field limit reproducing the prescribed coupling. *The cost is a carrier requirement, not a broken conservation law.*

## Cost 2 — the operational Bell test (**replaces the proxy**)

| `\|t\|` span | `S` operational | `S` post-selected | `P(both click)` |
|---:|---:|---:|---:|
| 0.2 | **2.032289** | 2.035101 | 0.986741 |
| 0.4 | **2.117169** | 2.135025 | 0.947838 |
| 0.6 | **2.222481** | 2.284225 | 0.885839 |
| 0.8 | **2.306045** | 2.458690 | 0.804739 |
| 1.0 | **2.330905** | 2.628578 | 0.709620 |
| 1.2 | **2.330891** | 2.761508 | 0.606212 |
| 1.5 | **2.328751** | 2.825902 | 0.448355 |
| 1.9 | **2.326141** | 2.828028 | 0.259820 |

Leakage is a genuine **no-click outcome**, not an efficiency proxy; probabilities normalize to 1.000000000000 and both fixed assignments agree. Operational maximum **2.330905** — above the local bound at *every* span tested, against Tsirelson 2.828427. *The first version's narrow `S = 2.13` window was an artifact of combining two different experiments.*

## Verdict

**THE_MINIMAL_MISSING_INTERACTION_IS_THE_SECOND_CHI_HARMONIC_AND_IT_WORKS_BUT_BOTH_COSTS_WERE_OVERSTATED_TOTAL_CHARGE_IS_CONSERVED_ONCE_THE_WINDING_TWO_CARRIER_IS_INCLUDED_AND_THE_OPERATIONAL_BELL_VALUE_IS_2_33_BROADLY_NOT_A_NARROW_2_13_WINDOW**
