# The SN-phenomenology audit - companion probe (PR #205)

**Run:** 2026-07-07T08:24:29+00:00

The deliverable is `docs/sn_phenomenology_audit.md` - the lab phenomenology of the #204-committed Phi[rho], with the sourcing ambiguity adjudicated by the committed dynamics itself. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the commitment (#204); the pick #198 left open | **PASS** |
| T2 | `T2_two_readings` | the two readings; why the dynamics can adjudicate | **PASS** |
| T3 | `T3_beamsplitter_pick` | the beamsplitter: whole-body transport -> conditional | **PASS** |
| T4 | `T4_branch_attraction_rate` | branch attraction real, exactly Newtonian; weight f^2 | **PASS** |
| T5 | `T5_classical_channel_cannot_entangle` | classical channel cannot entangle (BMV null) | **PASS** |
| T6 | `T6_lab_confrontation` | lab confrontation: nothing excluded; two nulls predicted | **PASS** |
| T7 | `T7_consistency_backward` | consistent with #198/#199/#204 | **PASS** |
| T8 | `T8_scope_and_assessment` | register update: the nearest-term falsifiers | **PASS** |

## The beamsplitter sweep (R, T by launch velocity)

| v | BAM (R, T) | linear control (R, T) |
|---:|---|---|
| 0.3 | 0.9999, 0.0001 | -, - |
| 0.4 | 1.0, 0.0 | 0.8211, 0.1299 |
| 0.5 | 0.9088, 0.0912 | -, - |
| 0.6 | 0.0002, 0.9998 | 0.6915, 0.1951 |
| 0.7 | 0.0, 1.0 | -, - |
| 0.8 | 0.0, 1.0 | 0.4677, 0.3633 |
| 1 | 0.0, 0.9998 | 0.2928, 0.5197 |
| 1.2 | 0.0, 0.9992 | 0.175, 0.6375 |

(whole-body at 7/8; bracket (0.5, 0.6); leakage 6e-05; lab kinetic/binding 5e-13)

## The branch-attraction rate (measured / Newtonian)

| t | d measured | d predicted | ratio |
|---:|---:|---:|---:|
| 1.0 | 11.799 | 11.7988 | 1.0 |
| 2.0 | 11.194 | 11.1935 | 1.0 |
| 3.0 | 10.180 | 10.179 | 1.0001 |
| 4.0 | 8.753 | 8.7469 | 1.0007 |
| 5.0 | 6.922 | 6.8851 | 1.0054 |

## The discriminators

- classical channel: S_max = 0.0e+00 vs quantized comparator S = 0.1484 - BMV null predicted (witness phase would be 0.791 rad)
- SN scale: m* = {'100nm': 5670000000.0, '500nm': 3310000000.0, '1um': 2630000000.0} amu; omega_SN(Si) = 0.0505 1/s - BAM predicts null; existing record sits 1e+05 below m*

## Verdict

**PHI_SOURCES_THE_ACTUAL_CONFIGURATION_SN_SIGNATURES_PREDICTED_NULL_BMV_ENTANGLEMENT_PREDICTED_NULL_EXISTING_BOUNDS_UNTOUCHED_TWO_NEAR_TERM_DISCRIMINATORS.** THE PROGRAM PICKED, IN WRITING - AND THE DYNAMICS DID THE PICKING (the argument is in docs/sn_phenomenology_audit.md; this probe runs the adjudication on the committed structure).

THE PICK. The beamsplitter experiment: the bound throat-soliton transmits or reflects WHOLE (7/8 velocities at max(R,T) >= 0.95; sub-threshold leakage 6e-05; co-occupation confined to the bracket (0.5, 0.6)), while the linear control co-occupies branches across the whole sweep. The lab regime (kinetic/binding ~ 5e-13) sits eleven orders below even the sandbox threshold: the mass-carrying field never occupies both arms. Effective-level sourcing is CONDITIONAL - #198's aside, measured.

THE CHANNEL AND ITS WEIGHT. The branch attraction is real and exactly the field-equation Newtonian rate (ratio 1.000-1.005 to merger) - the SN signatures scale as f^2, and BAM's f is exponentially zero in the lab regime: SN SIGNATURES PREDICTED NULL.

THE ENTANGLEMENT DISCRIMINATOR. The classical Phi cannot entangle (S = 0e+00 vs 0.1484 for the quantized comparator at the same coupling): BMV WITNESS PREDICTED NULL where quantized gravity gives 0.791 rad.

THE CONFRONTATION. Existing data exclude nothing (SN phase at the interferometry record 5e-17 rad; mass margin 1e+05 below m*): the audit ends in outcome (3), sharpened - BAM's nearest-term falsifiers are two NULLS (SN-scale signatures; the BMV witness) that active experimental programs are trying to violate. The register gains its first living-experimentalist channel.
