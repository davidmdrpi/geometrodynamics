# Pinhole-origin probe

**Run:** 2026-05-05T02:08:07+00:00
**Targets:** γ_lepton = 22.5, γ_quark = 22.25

Asks whether the lepton pinhole γ ≈ 22.5 (active at k = 3, 5 in `_build_generation_block`) has a natural origin in the Tangherlini machinery. Four categories of candidates: barrier-maximum sums, non-ground eigenvalues ω(l, n)², bound-state matrix elements ⟨u_l|O|u_l⟩, and smaller closure quanta. A candidate ‘explains’ the pinhole if its value is within ≤ 5 % of either γ_lepton or γ_quark.

## Best per category

| category | best candidate | formula | value | %Δ vs γ_lep | %Δ vs γ_qrk |
|---|---|---|---:|---:|---:|
| `barrier_sum` | `Sum_l=1..5_V_max[chebyshev_N80]` | `Σ_{l=1..5} V_max(l)` | 22.0082 | -2.19% | -1.09% |
| `non_ground_mode` | `omega_sq_l=1_n=4` | `ω(l=1, n=4)²` | 22.6666 | +0.74% | +1.87% |
| `matrix_element` | `<u_5|V_5|u_5>` | `⟨u_l|V_l|u_l⟩, l=5` | 0.7523 | -96.66% | -96.62% |
| `closure_quantum` | `beta_lepton_over_7` | `50π / 7` | 22.4399 | -0.27% | +0.85% |

## Candidates within 5% of γ_lepton = 22.5

| candidate | category | formula | value | %Δ vs γ_lep | %Δ vs γ_qrk |
|---|---|---|---:|---:|---:|
| `Sum_l=1..5_V_max[chebyshev_N80]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 22.0082 | -2.19% | -1.09% |
| `Sum_l=1..5_V_max[dense_r]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 23.1872 | +3.05% | +4.21% |
| `Sum_l=1..5_V_max[uniform_rstar]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 22.0082 | -2.19% | -1.09% |
| `omega_sq_l=1_n=4` | non_ground_mode | `ω(l=1, n=4)²` | 22.6666 | +0.74% | +1.87% |
| `omega_sq_l=2_n=4` | non_ground_mode | `ω(l=2, n=4)²` | 22.9784 | +2.13% | +3.27% |
| `omega_sq_l=3_n=4` | non_ground_mode | `ω(l=3, n=4)²` | 23.4195 | +4.09% | +5.26% |
| `Sum_l=1..5_omega_sq_n=1` | non_ground_mode | `Σ_{l=1..5} ω(l, 1)²` | 23.3541 | +3.80% | +4.96% |
| `7_pi` | closure_quantum | `7 · π` | 21.9911 | -2.26% | -1.16% |
| `beta_lepton_over_7` | closure_quantum | `50π / 7` | 22.4399 | -0.27% | +0.85% |
| `15_half_pi` | closure_quantum | `15/2 · π` | 23.5619 | +4.72% | +5.90% |

## Candidates within 5% of γ_quark = 22.25

| candidate | category | formula | value | %Δ vs γ_lep | %Δ vs γ_qrk |
|---|---|---|---:|---:|---:|
| `Sum_l=1..5_V_max[chebyshev_N80]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 22.0082 | -2.19% | -1.09% |
| `Sum_l=1..5_V_max[dense_r]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 23.1872 | +3.05% | +4.21% |
| `Sum_l=1..5_V_max[uniform_rstar]` | barrier_sum | `Σ_{l=1..5} V_max(l)` | 22.0082 | -2.19% | -1.09% |
| `omega_sq_l=1_n=4` | non_ground_mode | `ω(l=1, n=4)²` | 22.6666 | +0.74% | +1.87% |
| `omega_sq_l=2_n=4` | non_ground_mode | `ω(l=2, n=4)²` | 22.9784 | +2.13% | +3.27% |
| `Sum_l=1..5_omega_sq_n=1` | non_ground_mode | `Σ_{l=1..5} ω(l, 1)²` | 23.3541 | +3.80% | +4.96% |
| `7_pi` | closure_quantum | `7 · π` | 21.9911 | -2.26% | -1.16% |
| `beta_lepton_over_7` | closure_quantum | `50π / 7` | 22.4399 | -0.27% | +0.85% |

## Mass sensitivity to γ

Plug each candidate γ into the locked lepton block, diagonalize, anchor the lightest eigenvalue to m_e = 0.511 MeV, and read off the predicted m_μ and m_τ. The locked surrogate's claimed agreement (errors ≤ 0.2%) requires a γ that hits within roughly a percent of 22.5 — the radial geometric origin pins the SCALE but does not by itself pin the value to mass-fit precision.

| candidate γ | value | m_μ (MeV) | m_τ (MeV) | err μ | err τ |
|---|---:|---:|---:|---:|---:|
| `locked_baseline_22.5` | 22.5000 | 105.613 | 1778.938 | 0.04% | 0.12% |
| `gamma_quark_22.25` | 22.2500 | 131.380 | 2223.483 | 24.34% | 25.14% |
| `Sum_l=1..5_V_max[chebyshev_N80]` | 22.0082 | 173.049 | 2942.211 | 63.78% | 65.58% |
| `Sum_l=1..5_V_max[dense_r]` | 23.1872 | 69.547 | 1156.359 | 34.18% | 34.92% |
| `omega_sq_l=1_n=4` | 22.6666 | 93.591 | 1571.480 | 11.42% | 11.56% |
| `omega_sq_l=2_n=4` | 22.9784 | 77.422 | 1292.363 | 26.72% | 27.27% |
| `omega_sq_l=3_n=4` | 23.4195 | 62.592 | 1036.194 | 40.76% | 41.68% |

## Verdict

**The SCALE of γ has a natural Tangherlini origin** but the PRECISE value 22.5 is calibrated. Three independent natural candidates land within ~3% of γ_lepton:

- `Σ_{l=1..5} V_max(l)` = 22.0082 on the canonical Chebyshev grid (-2.19% vs γ_lepton). This is the same formula the QCD residual-sector pinhole already uses — locked there at γ_quark = 22.25 with 1.1% offset from the same Chebyshev value. The lepton and quark pinholes are the SAME radial geometric quantity, sampled with small (~1–3%) calibration offsets.
- `ω(l=1, n=4)²` ≈ 22.67 (+0.74%): a non-ground Tangherlini eigenvalue. Cleaner numerical match than the barrier sum, but the (l=1, n=4) mode has no obvious physical role in the locked surrogate — coincidence or hidden structure.
- `7π` ≈ 21.99 (−2.26%): a small closure-quantum integer. Suggestive but not directly tied to a known BAM channel; the match is comparable to the barrier-sum offset and could be incidental to both falling near the same scale.

**Sensitivity caveat: the muon-mass match is sharper than the γ-origin candidates.** Plugging the geometric γ = Σ V_max = 22.008 into the locked block gives m_μ off by ~64% (see Mass sensitivity table); the locked γ = 22.5 hits 0.04%. The lepton ladder requires γ accurate to better than 1%, while the radial barrier-sum origin only fixes γ to ~3%. Either (a) there is a small calibration term on top of the geometric base that the present probe has not isolated, or (b) the precise lepton-ladder agreement is itself tuned to ~1% of a geometric anchor that this probe cannot distinguish from the family of nearby candidates.

### Bottom line

γ ≈ 22.5 has a clean radial-geometric origin at the SCALE of Σ_{l=1..5} V_max(l) = 22.0; this is the same operator the QCD residual sector locks at 22.25. The remaining ~2.2% gap between the geometric value and the locked γ_lepton = 22.5 is small in absolute terms but determines the lepton ladder agreement at the fraction-of-percent level — meaning the pinhole is approximately, but not exactly, a single Tangherlini matrix-element quantity. Closing the residual 2.2% is the next concrete probe target — candidate refinements include a different sampling (dense_r gives 23.19, +3% on the other side; the locked value is intermediate), an additive correction from the sub-leading ω(l, 1)² spectrum (Σ ω(l, 1)² ≈ 23.35 within 3.8%), or a small turning-point evaluation rather than V_max.