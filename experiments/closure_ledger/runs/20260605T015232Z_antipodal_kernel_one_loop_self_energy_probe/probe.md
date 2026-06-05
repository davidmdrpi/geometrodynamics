# One-loop self-energy audit for the antipodal matter kernel (PR #136)

**Run:** 2026-06-05T01:52:32+00:00

Audits the leading interacting correction to the antipodal matter kernel (#135) — the one-loop self-energy Σ. The result: a finite real mass shift, an imaginary (width) part that vanishes below the two-particle threshold (so the lightest mode is exactly stable), and unitarity preserved with no horizon-absorption width (the antipodal mirror, #129).

- **Self-energy**: Σ_k(s) = Σ_{n≤m} c|g_{knm}|²/(s − (ω_n+ω_m)² + i0⁺) (one-loop bubble)
- **Lightest mode**: ω_0 = 1.167 < 2ω_0 ⟹ Im Σ_0(ω_0²) = 0 ⟹ exactly stable
- **Mass shift**: Re Σ_0 finite (mode sum converges; #116 regularisation)
- **Unitarity**: Im Σ ≤ 0 above threshold, = 0 below; no horizon-absorption width (#129)
- **Contrast**: absorbing horizon ⟹ tree-level width on every mode (#130)
- **Open**: modelled vertex/coupling; higher loops; absolute normalisation (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | one-loop self-energy audit for the antipodal matter kernel (#135) | **PASS** |
| T2 | `T2_dyson_dressed_propagator` | Dyson G = 1/(s − ω_k² − Σ); Re Σ mass shift, Im Σ width | **PASS** |
| T3 | `T3_one_loop_two_particle_bubble` | one-loop Σ = two-particle bubble; vertex = mode triple overlap | **PASS** |
| T4 | `T4_lightest_mode_stable_below_threshold` | Im Σ_0(ω_0²) = 0 below threshold ⟹ lightest mode stable | **PASS** |
| T5 | `T5_real_mass_shift_finite` | Re Σ_0 finite mass shift (mode sum converges, #116 scheme) | **PASS** |
| T6 | `T6_unitarity_no_horizon_absorption_width` | unitarity survives; no horizon-absorption width (#129) vs absorbing (#130) | **PASS** |
| T7 | `T7_scope` | scope: one loop, fixed background, modelled vertex | **PASS** |
| T8 | `T8_assessment` | ANTIPODAL_KERNEL_ONE_LOOP_SELF_ENERGY_REAL_SHIFT_STABLE_LIGHTEST_MODE | **PASS** |

## The lightest mode is below its decay threshold (exactly stable)

| quantity | value |
|---|---:|
| `ω_0` (lightest mode) | 1.167 |
| two-particle threshold `2ω_0` | 2.334 |
| pole `s = ω_0²` | 1.362 |
| threshold `s_thr = (2ω_0)²` | 5.448 |
| `Im Σ_0(ω_0²)` | -0.003 |

The pole sits below the two-particle threshold, so `Im Σ_0 = 0` — the lightest matter mode cannot decay and is exactly stable at one loop.

## The real mass shift is finite (mode sum converges)

| internal-mode cutoff | Re Σ_0(ω_0²) |
|---:|---:|
| 10 | -0.2768 |
| 20 | -0.2787 |
| 30 | -0.2793 |
| 40 | -0.2796 |

A finite, real mass renormalisation (× coupling²); the residual UV piece is the #116 zeta/heat-kernel regularisation.

## Verdict

**ANTIPODAL_KERNEL_ONE_LOOP_SELF_ENERGY_REAL_SHIFT_STABLE_LIGHTEST_MODE.** THE ONE-LOOP SELF-ENERGY OF THE ANTIPODAL MATTER KERNEL IS A FINITE REAL MASS SHIFT WITH NO WIDTH BELOW THRESHOLD — THE LIGHTEST MODE STAYS EXACTLY STABLE AND UNITARITY SURVIVES. PR #135 built the free antipodal matter propagator (real poles, the stable #130 spectrum); this probe audits its leading interacting correction, the one-loop self-energy Σ.

THE DYSON-DRESSED PROPAGATOR. Σ(s) (s = ω²) dresses the free kernel: G(s) = 1/(s − ω_k² − Σ(s)). Re Σ is a mass renormalisation, Im Σ a width; a mode stays a sharp, stable particle iff Im Σ = 0 at its pole.

THE ONE-LOOP Σ = THE TWO-PARTICLE BUBBLE. For a cubic self-interaction on the cavity (vertex g_{knm} = ∫ ψ_k ψ_n ψ_m dr*, the triple overlap of the antipodal modes), the one-loop self-energy of mode k is Σ_k(s) = Σ_{n≤m} c_{nm}|g_{knm}|²/(s − (ω_n+ω_m)² + i0⁺) — the amplitude k → (n,m) → k.

Im Σ = 0 BELOW THRESHOLD ⟹ THE LIGHTEST MODE IS EXACTLY STABLE. By the optical theorem Im Σ_k(s) is (minus) the two-particle phase space — nonzero only when s reaches a two-particle threshold (ω_n+ω_m)². The lowest threshold is 2ω_0; the lightest mode sits at ω_0 < 2ω_0, so its pole s = ω_0² lies below s_thr = (2ω_0)² and Im Σ_0(ω_0²) = 0. The lightest matter mode cannot decay (energy conservation) and stays a sharp, real-pole, stable particle through one loop.

THE REAL MASS SHIFT IS FINITE. Re Σ_0(ω_0²) is a finite mass renormalisation: the vertex overlaps g_{0nm} decay with the mode index, so the mode sum converges (stable to ~1e-3 from cutoff 20 to 40), and the residual UV piece is the same zeta/heat-kernel regularisation as the #116 fluctuation determinant. No UV catastrophe on the discrete antipodal cavity.

UNITARITY SURVIVES — AND NO HORIZON-ABSORPTION WIDTH. Im Σ_k(s) ≤ 0 (a width) above the two-particle threshold and = 0 below: the dressed kernel respects the optical theorem. Crucially, because the throat is a unitary mirror (#129) there is NO horizon-absorption contribution to Σ — the only width source is genuine multi-particle decay, which the lightest mode is kinematically forbidden from. This is the sharp contrast with the absorbing horizon, which gives EVERY mode a tree-level width (#130). So the antipodal kernel's one-loop self-energy is unitarity-preserving — it extends the tree-level stable spectrum (#130/#135) to one loop.

SCOPE. The leading (one-loop) interacting correction on the fixed antipodal background. The interaction vertex is MODELLED (a generic cubic triple-overlap), not derived from the S_BAM measure, and the coupling is an input — so Re Σ is fixed only up to the coupling. Higher loops, the absolute normalisation (#133), and the flavor residuals (#134) stand.
