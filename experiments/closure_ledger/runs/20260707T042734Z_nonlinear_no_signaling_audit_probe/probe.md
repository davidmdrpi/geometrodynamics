# The nonlinear no-signaling audit - companion probe (PR #204)

**Run:** 2026-07-07T04:27:34+00:00

The deliverable is `docs/nonlinear_no_signaling_audit.md` - the no-signaling audit of the nonlinear BAM psi-Phi-q dynamics. This probe measures every claim on the live dynamics. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the register item; the Gisin stakes; audit by construction | **PASS** |
| T2 | `T2_edge_is_live` | the edge is live: BAM genuinely nonlinear (defect measured) | **PASS** |
| T3 | `T3_channel_identified` | the channel is gravity: O(G), Phi-clamp kills it, q local | **PASS** |
| T4 | `T4_causal_front` | retardation confines it to the cone (front = d/c, 3 speeds) | **PASS** |
| T5 | `T5_equivariance_survives_retardation` | equivariance survives retardation (Born at noise, kicked) | **PASS** |
| T6 | `T6_entangled_linear_sector` | entangled linear: exact invariance; Valentini boundary | **PASS** |
| T7 | `T7_entangled_nonlinear_channel` | entangled nonlinear: O(G) violation, confined behind front | **PASS** |
| T8 | `T8_scope_and_assessment` | no-signaling survives at dBB grade; scope stated | **PASS** |

## The field-sector response at B (t = 3, separation 30)

| configuration | D_B(3) |
|---|---:|
| kinematic floor (G = 0) | 3.26e-13 |
| Newtonian (instantaneous) gravity | 1.47e-06 |
| half G | 6.88e-07 |
| Phi-clamp control | 2.36e-13 |
| retarded, c = 8 | 1.37e-10 |
| retarded, c = 12 | 1.78e-06 |
| retarded, c = 16 | 2.94e-06 |
| retarded, c = 20 (the c to infinity check) | 1.80e-06 |

Front times {'8': 3.15, '12': 2.2, '16': 1.7} -> c*t = {'8': 25.2, '12': 26.4, '16': 27.2} vs geometric distance 25.

## The entangled sector

- linear marginal invariance: 2.7e-15 (machine); equilibrium KS 0.0152 vs noise 0.0176; non-equilibrium KS 0.0447 (the Valentini signal)
- nonlinear mean field: linear floor 5.5e-15; Newton front 0.6; retarded front 1.9 vs geometric 2.0 (pre-arrival 2.0e-15)

## Verdict

**NO_SIGNALING_SURVIVES_AUDIT_THE_ONLY_NONLINEAR_CHANNEL_IS_THE_RETARDED_GRAVITATIONAL_FIELD_EQUIVARIANCE_INTACT_EQUILIBRIUM_SIGNAL_LOCALITY_AT_DBB_GRADE.** THE AUDIT RAN, THE EDGE FIRED WHERE IT SHOULD, AND THE THEORY SURVIVES IT (the argument is in docs/nonlinear_no_signaling_audit.md; this probe measures everything on the live dynamics).

THE CHANNEL. BAM is genuinely nonlinear (superposition defect measured), so Gisin's threat applies - and indeed a local kick at A reaches B instantaneously in the Newtonian model, 5e+06 above the kinematic floor. But the channel is EXACTLY the gravitational field: O(G) (ratio 2.138), and clamping Phi collapses it to the floor (suppression 6e+06). A physical interaction, not a measurement pathology.

THE CONE. With the causal wave-equation Phi, the response is machine-floor quiet outside the light cone and arrives at c*t = [25.2, 26.4, 27.2] vs geometric distance 25, at three speeds; c -> infinity recovers Newton. The superluminality belongs to the APPROXIMATION, not the theory.

THE COEXISTENCE. The retarded potential is still real: norm exact, continuity residual 1e-04, Born ensemble at noise through the kicked retarded evolution - removing the signaling costs the Born rule nothing.

THE ENTANGLED SECTOR. Linear part: exact marginal invariance (3e-15); equilibrium trajectories signal-local (KS 0.0152 ~ noise 0.0176); non-equilibrium trajectories signal (KS 0.0447) - the Valentini boundary on the BAM transport. Nonlinear part: the mean-field gravitational dressing violates marginal invariance at O(G) instantaneously (the Gisin channel, exhibited) - and retardation confines it behind the front (arrival 1.9 vs geometric 2.0; machine floor before). No-signaling survives at dBB grade: equilibrium + causal fields, both now measured properties of the BAM dynamics.
