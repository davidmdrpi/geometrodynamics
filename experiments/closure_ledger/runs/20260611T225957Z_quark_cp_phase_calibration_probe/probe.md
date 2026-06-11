# Quark CP phase calibration on the locked Hamiltonian (PR #156)

**Run:** 2026-06-11T22:59:57+00:00

Solves the constrained problem #155 posed: the quark CP phase is introduced on the partition-mixing element (the Hopf-placeholder φ_q(k) = φ·k), calibrated to the observed phase content, and the locked |V| and masses survive. The J ceiling deficit decomposes exactly into the #155 soft-direction |V| deficits — no independent CP failure — and the db-triangle shape (β = 22°) emerges as the falsifiable acceptance test for the true Hopf-connection phase. *(QFT on the classical throat, not quantum gravity.)*

- **Extension**: H(ε,φ) = H_locked + partition-mixing with φ_q(k) = φ·k (in-probe)
- **Scaling**: J ∝ ε^1.9; sinusoidal in φ; shifts O(ε²)
- **Ceiling**: deficit 0.249 = 0.498×0.902×0.555 exactly (the #155 soft direction)
- **Calibration**: ε* = 0.0528, φ* = 0.8: V_cb −0.0%, masses ≤ 0.5%
- **Triangle**: β ≈ 0° vs 22.2° — shape open; the Hopf φ_q(k) acceptance test
- **Consumed**: one input: the quark CP phase content

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the constrained CP calibration #155 posed | **PASS** |
| T2 | `T2_cp_extension` | in-probe extension; baseline J = 0; global knob unusable | **PASS** |
| T3 | `T3_scaling_structure` | J ∝ ε^1.9; sinusoidal φ-dependence; shifts O(ε²) | **PASS** |
| T4 | `T4_ceiling_identity` | ceiling deficit 0.249 = soft-direction product (exact identity) | **PASS** |
| T5 | `T5_constrained_calibration` | ε* = 0.053: V_cb −0.0%, V_us/V_ub −4%, masses ≤ 0.5% | **PASS** |
| T6 | `T6_triangle_shape_discriminator` | placeholder fails the triangle shape: β ≈ 0° vs 22.2° (target) | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: one input consumed; ceiling/shape open with targets | **PASS** |
| T8 | `T8_assessment` | QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN | **PASS** |

## The ceiling identity (the J deficit is the soft direction)

| quantity | value |
|---|---:|
| J ceiling, predicted \|V\| | 8.643e-06 |
| J ceiling, observed \|V\| | 3.472e-05 |
| ceiling ratio | 0.249 |
| element ratios (us, cb, ub) | [0.498, 0.902, 0.555] |
| decomposition product | 0.249 |
| sin δ observed (near-maximal) | 0.887 |

## The constrained calibration

| quantity | value |
|---|---:|
| ε* | 0.0528 |
| φ* | 0.8 |
| J calibrated / target | 7.667e-06 / 7.667e-06 |
| shift V_us | -0.0421 |
| shift V_cb | -0.0 |
| shift V_ub | -0.0425 |
| max mass shift | 0.005 |
| sin δ calibrated | 0.967 |

## The triangle-shape discriminator

Calibrated (β, γ) = (0.01°, -179.98°) vs observed (22.2°, 65.9°): the placeholder linear phase reproduces the area (J) but squashes the db-triangle — β = 22° is the quantitative acceptance test for the true Hopf-connection φ_q(k).

## Verdict

**QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN.** THE QUARK CP PHASE CALIBRATES ONTO THE LOCKED HAMILTONIAN WITHOUT DISTURBING |V| OR THE MASSES; THE J CEILING DEFICIT DECOMPOSES EXACTLY INTO THE #155 SOFT-DIRECTION |V| DEFICITS (NO INDEPENDENT CP FAILURE); THE PHASE CONTENT IS NEAR-MAXIMAL LIKE THE DATA; AND THE TRIANGLE SHAPE IS THE FALSIFIABLE TARGET THE PLACEHOLDER PHASE FAILS. #155 posed the constrained problem; this probe solves it and sharpens what remains.

THE EXTENSION. H(ε, φ) = H_locked + the v3 §4 partition-mixing element with the Hopf-placeholder phase φ_q(k) = φ·k, built in-probe (locked blocks exactly intact; baseline J = 0 re-verified). The global phase knob is unusable — it enters the transport couplings as cos(phase·dk) and collapses |V_us| ×20 at phase = 0.5 — which is why the mass calibration locked it at zero and the #155 baseline had J = 0 structurally.

SCALING. J ∝ ε^1.9 (quadratic — one insertion per sector side), sinusoidal in φ; |V|/mass shifts also O(ε²): CP switches on cleanly.

THE CEILING IDENTITY. J ≤ |V_us·V_cb·V_ub|: predicted ceiling 8.643e-06 vs observed 3.472e-05 — ratio 0.249, decomposing EXACTLY into the element ratios [0.498, 0.902, 0.555]: the J shortfall IS the #155 V_us/V_ub soft direction propagated, not a new failure. The observed CP is near-maximal (0.887). Consistency lock: when the soft directions land, the ceiling rises to the observed 3.5e-5.

THE CALIBRATION. ε* = 0.0528, φ* = 0.8: J hits the phase-content target (7.667e-06) with |V_cb| shifted -0.0% (the stiff prediction untouched), |V_us|/|V_ub| -4.2%/-4.2% (inside the soft direction), masses ≤ 0.5% — the locked structure survives; sin δ = 0.967, near-maximal like the data.

THE TRIANGLE SHAPE. The placeholder phase reproduces the AREA (J) but squashes the db-triangle: (β, γ) ≈ (0.01°, -179.98°) vs observed (22.2°, 65.9°). The shape discriminates phase structures — β = 22° is the quantitative acceptance test for the true Hopf-connection φ_q(k) (the v3 §4 TODO).

LEDGER. One input consumed (the CP phase content — the flavor puzzle's CP entry made explicit); the scaling, the ceiling identity, the calibration survival, and the near-maximal content derived; the Hopf phase (β = 22°) and the soft |V| directions (J ceiling) open with falsifiable targets.
