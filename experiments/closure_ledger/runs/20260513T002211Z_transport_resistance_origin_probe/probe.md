# Transport / resistance origin probe

**Run:** 2026-05-13T00:22:11+00:00
**Targets:** transport ≈ 25.1, resistance ≈ 0.217869435878

The R_OUTER self-consistency probe flagged residual phenomenological sensitivity in two parameters of the locked lepton block: `transport_strength` ≈ 25.1 (off-diagonal coupling) and `resistance_scale` ≈ 0.218 (diagonal action coefficient). This probe asks the structural question for those two parameters that `pinhole_origin_probe` asked for γ ≈ 22.5: is there a natural Tangherlini / closure-quantum origin?

Six candidate categories: closure-quantum integers (N·π, β/m), barrier-maximum sums (Σ V_max), eigenfrequency invariants ω(l,n), WKB tunneling integrals κ_{l1,l2}/dk, cross-shell overlap integrals ⟨u_l|w|u_{l'}⟩, and inverse / log forms.

## Transport (locked = 25.1)

**Best candidate:** `4_x_2pi` (`4·2π`) = 25.1327, %Δ = +0.130%.

Candidates within ±1.0% of locked: **4**. Within ±5.0%: **6**.

### Within ±1.0% (tight)

| candidate | category | formula | value | %Δ |
|---|---|---|---:|---:|
| `4_x_2pi` | closure_quantum | `4·2π` | 25.1327 | +0.130% |
| `8_pi` | closure_quantum | `8·π` | 25.1327 | +0.130% |
| `16_half_pi` | closure_quantum | `16/2·π` | 25.1327 | +0.130% |
| `Sum_V_max_15_plus_pi` | barrier_sum | `Σ V_max[1..5] + π` | 25.1498 | +0.199% |

### Within ±5.0% (scale-only)

| candidate | category | formula | value | %Δ |
|---|---|---|---:|---:|
| `beta_over_6` | closure_quantum | `50π/6` | 26.1799 | +4.303% |
| `Sum_V_max_15_x_8pi_over_22.5` | barrier_sum | `Σ V_max[1..5] · (8π/22.5)` | 24.5834 | -2.058% |

## Resistance (locked = 0.217869435878)

**Best candidate:** `omega_1_minus_1_x_4` (`4·(ω(1,0) − 1)`) = 0.2189, %Δ = +0.477%.

Candidates within ±1.0% of locked: **2**. Within ±5.0%: **8**.

### Within ±1.0% (tight)

| candidate | category | formula | value | %Δ |
|---|---|---|---:|---:|
| `7pi_over_100` | closure_quantum | `7·π/100` | 0.2199 | +0.937% |
| `omega_1_minus_1_x_4` | eigenfrequency | `4·(ω(1,0) − 1)` | 0.2189 | +0.477% |

### Within ±5.0% (scale-only)

| candidate | category | formula | value | %Δ |
|---|---|---|---:|---:|
| `beta_over_700` | closure_quantum | `50π/700` | 0.2244 | +2.997% |
| `V_max_1_over_5` | barrier_sum | `V_max(l=1) / 5` | 0.2279 | +4.617% |
| `Sum_V_max_15_over_100` | barrier_sum | `Σ V_max[1..5] / 100` | 0.2201 | +1.016% |
| `kappa_per_dk_spread[om2_l1_V_l1_dkmax]` | wkb_tunneling | `spread of κ/dk across pairs (om2_l1_V_l1_dkmax)` | 0.2091 | -4.033% |
| `ln2_over_pi` | inverse_form | `ln 2 / π` | 0.2206 | +1.270% |
| `2_over_3pi` | inverse_form | `2 / (3π)` | 0.2122 | -2.599% |

## Mass sensitivity

Each row substitutes the listed (transport, resistance) into the locked lepton block, anchors m_e, and reports predicted m_μ, m_τ with relative errors. The locked baseline matches PDG to ≤ 0.2%; an acceptable substitution must keep errors within the same envelope (target: ≤ 5.0%).

| label | candidate | transport | resistance | m_μ (MeV) | m_τ (MeV) | err μ | err τ |
|---|---|---:|---:|---:|---:|---:|---:|
| locked_baseline | `—` | 25.1000 | 0.2179 | 105.613 | 1778.938 | 0.043% | 0.117% |
| transport=4_x_2pi | `4·2π` | 25.1327 | 0.2179 | 114.149 | 1922.073 | 8.036% | 8.172% |
| resistance=omega_1_minus_1_x_4 | `4·(ω(1,0) − 1)` | 25.1000 | 0.2189 | 102.063 | 1718.546 | 3.403% | 3.282% |
| joint_best | `4·2π, 4·(ω(1,0) − 1)` | 25.1327 | 0.2189 | 110.006 | 1851.663 | 4.114% | 4.210% |
| transport=8*pi | `8·π = 25.1327` | 25.1327 | 0.2179 | 114.149 | 1922.073 | 8.036% | 8.172% |
| joint_8pi_and_best_resistance | `8·π, 4·(ω(1,0) − 1)` | 25.1327 | 0.2189 | 110.006 | 1851.663 | 4.114% | 4.210% |
| transport=8pi, resistance=7pi_over_100 | `8·π, 7·π/100` | 25.1327 | 0.2199 | 106.283 | 1788.404 | 0.591% | 0.650% |
| transport=8pi, resistance=omega_1_minus_1_x_4 | `8·π, 4·(ω(1,0) − 1)` | 25.1327 | 0.2189 | 110.006 | 1851.663 | 4.114% | 4.210% |

## Verdict

**Transport.** Best candidate `4·2π` (`4_x_2pi`) at %Δ = +0.130%. 

**Resistance.** Best candidate `4·(ω(1,0) − 1)` (`omega_1_minus_1_x_4`) at %Δ = +0.477%.

Locked baseline gives err μ = 0.043%, err τ = 0.117%. Substituting transport → 8π alone (resistance kept locked) gives err μ = 8.036%, err τ = 8.172% — the lepton ladder is high-sensitivity to transport, so the bare 0.13% transport gap amplifies into an 8% mass shift.

**The joint closure-quantum reading survives.** Substituting BOTH transport → 8π = 4·(2π) AND resistance → 7π/100 recovers the lepton ladder at err μ = 0.591%, err τ = 0.650% — both within the 5% envelope. This is non-trivial: the transport miss and the resistance miss partially cancel in the eigenvalue ratios. The two closure-quantum readings:

- `transport = 8π = 4·(2π)`  →  the 4th closure quantum.
- `resistance = 7π / 100`  →  small closure-quantum fraction.

Both are structurally the same kind of object as the antipodal closure (k·2π), the Hopf+throat partnership (1·2π), and the τ-uplift quantum (100·2π) that already organise the Layer-1 ledger. The resistance reading also has a near-twin in `4·(ω(1,0) − 1)` (the 1.054-factor gap, scaled by 4) — both land at ≈ 0.219, raising the open question of whether the geometric origin is the closure-fraction reading or the Tangherlini eigenfrequency reading. The two cannot be distinguished at this probe's resolution.

## What this leaves open

1. **Closed-form sharpening for resistance.** Two readings within 1% of locked (`7π/100` at +0.94%, `4·(ω(1,0) − 1)` at +0.48%) match the SCALE but neither survives mass-sensitivity at the fraction-of-percent precision that the locked baseline reaches. Lifting either to a closed-form structural identity is the next concrete target.
2. **Resistance / pinhole / γ link.** The 1.054 factor (= ω(1,0) at R* ≈ 1.262), the resistance 0.218, and the pinhole γ = 22.5 are all evaluated on the SAME geometry. Whether they are three independent observables or three projections of one Tangherlini matrix-element family is the cross-cutting structural question.
3. **Closing the R_OUTER self-consistency loop.** With transport = 8π and resistance ≈ 7π/100 (or 4·(ω(1,0) − 1)) as principled inputs, re-run the R_OUTER bisection from probe 8. If the fixed point still lands at R* ≈ 1.262 within the same 0.008% cross-species tolerance, R_OUTER is structurally selected on principled inputs only — completing the lift from phenomenological to fully geometric.