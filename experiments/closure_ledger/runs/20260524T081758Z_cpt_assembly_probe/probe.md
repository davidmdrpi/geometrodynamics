# Geometric CPT assembly for BAM throat histories

**Run:** 2026-05-24T08:17:58+00:00

Assembles C (inner/outer swap, #63), P (S³ reflection), and T (iσ_y, B2) into the antiunitary CPT symmetry on throat histories, guaranteed by the throat's local Lorentz invariance (#59–#60), mapping a throat to the antithroat run backwards (Stückelberg = pair production, #58).

- **CPT**: `C·P·T → q→−, p→+, x→−, s→−, t→−, E→+`
- **Operations**: C = inner/outer swap (#63); P = S³ reflection; T = iσ_y (B2)
- **Signatures**: C²=+1, P²=+1, T²=−I
- **Throat histories**: CPT(throat fwd) = antithroat bwd (Stückelberg = #58)
- **Theorem**: local Lorentz invariance (#59–#60) ⟹ CPT; global S³ breaking suppressed
- **B4 caveat**: discrete geometric operations; scale-independent

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_three_operations` | C=swap (#63), P=S³ reflection, T=iσ_y (B2) | **PASS** |
| T2 | `T2_involution_signatures` | C²=+1, P²=+1, T²=−I: True | **PASS** |
| T3 | `T3_cpt_transformation_table` | C·P·T → q−,p+,x−,s−,t−,E+: True | **PASS** |
| T4 | `T4_bam_realizations` | C↔inner/outer swap (#63); T↔iσ_y (B2) | **PASS** |
| T5 | `T5_stuckelberg_pair_production` | CPT(throat fwd)=antithroat bwd (#58); ΣQ=0 | **PASS** |
| T6 | `T6_cpt_theorem_from_local_lorentz` | local Lorentz ⟹ CPT; violation ~8e-78 | **PASS** |
| T7 | `T7_falsification_b4` | O(1) violation would falsify; BAM suppressed | **PASS** |
| T8 | `T8_assessment` | discrete-symmetry sector unified | **PASS** |

## T1: The three operations

- **C_charge_conjugation**: inner/outer swap S: r ↦ 2R_MID − r (#63) — c₁ → −c₁ (throat → antithroat)
- **P_parity**: spatial S³ reflection x → −x — p → −p; spin (axial) P-even
- **T_time_reversal**: iσ_y K (B2, antiunitary) — t → −t, s → −s, E → +E

## T2: (Anti)involution signatures

- C² = +1: True; P² = +1: True
- T² = (iσ_y K)² = [[-1.0, 0.0], [0.0, -1.0]] = −I: True (fermionic; the RP³ spin structure, B2)

## T3: CPT transformation table

| observable | C | P | T | CPT |
|---|---:|---:|---:|---:|
| q | -1 | +1 | +1 | -1 |
| p | +1 | -1 | -1 | +1 |
| x | +1 | -1 | +1 | -1 |
| s | +1 | +1 | -1 | -1 |
| t | +1 | +1 | -1 | -1 |
| E | +1 | +1 | +1 | +1 |

Matches the standard CPT (q→−, p→+, x→−, s→−, t→−, E→+): True.

## T4: BAM realizations

- C = inner/outer swap involution: True; c₁ → −c₁ (PR #63)
- T = iσ_y, T²=−I: True (B2)

## T5: Stückelberg / pair production

- throat forward: {'q': 1, 't_direction': 1}
- CPT image: {'q': -1, 't_direction': -1} (antithroat running backward: True)
- pair total charge: 0 (the #58 throat–antithroat pair)

## T6: CPT theorem from local Lorentz invariance

- local Lorentz invariance (#59–#60): True
- global Lorentz broken by S³: True
- CPT-violation suppression (R_MID/R_cosmo)² = 7.977e-78 (unobservably small)

## T7: Falsification / B4

- O(1) CPT violation would falsify: True
- BAM violation suppressed (8.0e-78): True
- operations dimensionless/geometric: True

## T8: Assessment

- CPT: q→−, p→+, x→−, s→−, t→−, E→+
- signatures: C²=+1, P²=+1, T²=−I
- throat histories: CPT(throat fwd) = antithroat bwd (Stückelberg = #58)
- theorem: local Lorentz invariance (#59–#60) ⟹ CPT
- remaining: full CPT operator on the throat spinor from S_BAM; P vs antipodal Z₂; observable bounds

## Verdict

**CPT_ASSEMBLED.** CPT ASSEMBLED. The three BAM discrete symmetries compose to the antiunitary CPT symmetry on throat histories, unifying the discrete-symmetry sector.

THE OPERATIONS. C = charge conjugation = the inner/outer swap (S: r ↦ 2R_MID − r, c₁ → −c₁, PR #63); P = parity = spatial S³ reflection (x→−x, p→−p); T = time reversal = iσ_y K (B2, antiunitary, t→−t, s→−s, E→+E). Their signatures are C²=+1, P²=+1, and T²=−I — the fermionic spin structure (the non-trivial RP³ spin structure, B2), with (iσ_y K)²=−I verified.

THE COMPOSITION. The sign tables compose to CPT: q→−q, p→+p, x→−x, s→−s, t→−t, E→+E — a particle (q,p,s,E>0) mapped to an antiparticle (−q,p,−s,E>0) with x,t reversed (E>0 preserved by T's antiunitarity).

THROAT HISTORIES. A throat going forward with c₁=+1 maps under CPT to an antithroat (c₁=−1) running backwards — the Feynman–Stückelberg antiparticle. This IS the pair-production structure (PR #58): a throat–antithroat pair is one worldline turning around in time at the nucleation point, the two arms related by CPT (total charge 0).

THE THEOREM. CPT is guaranteed for any local, Lorentz-invariant theory (Lüders–Pauli). The throat carries LOCAL Lorentz invariance (PRs #59–#60), so CPT is exact locally; the closed S³ breaks GLOBAL Lorentz invariance (a preferred frame, #59), so any CPT violation is suppressed by (R_MID/R_cosmo)² ~ 10⁻⁷⁸ — calculable, unobservably small. An O(1) violation would falsify; BAM passes. B4: C, P, T, CPT are discrete geometric operations (c₁ a topological integer, T²=−1 a group fact) — scale-independent. Remaining: the full CPT operator on the throat Dirac spinor from S_BAM, disentangling P from the antipodal Z₂ (B2), and observable CPT bounds.

## What this leaves open

- **The full CPT operator on the throat Dirac spinor** from S_BAM (the explicit Θ=CPT matrix and Θ²), beyond the sign table.
- **P vs the antipodal Z₂.** Disentangling spatial parity from the antipodal deck transformation of RP³ = S³/Z₂ (B2).
- **Observable CPT bounds.** Mapping (R_MID/R_cosmo)² to specific CPT-violation observables.
