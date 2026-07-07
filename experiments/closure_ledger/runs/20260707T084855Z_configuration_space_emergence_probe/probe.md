# Configuration-space emergence: entanglement is bridge topology - companion probe (PR #206)

**Run:** 2026-07-07T08:48:55+00:00

The deliverable is `docs/configuration_space_emergence.md` - the derivation of the effective two-throat state from the universal 3-space field plus the bulk identification. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the sharp target: Bell's theorem vs a 3-space-local field | **PASS** |
| T2 | `T2_no_go_and_gap` | LHV max = 2 (enumerated); Tsirelson 2*sqrt2 (the gap) | **PASS** |
| T3 | `T3_identification_on_lattice` | the identification on the lattice (conjugation/locking/holonomy/cut) | **PASS** |
| T4 | `T4_emergence_lemmas` | the emergence lemmas: (I x T)|Phi+> = singlet; two paths agree | **PASS** |
| T5 | `T5_er_epr_curve` | the ER=EPR curve: CHSH = 2*sqrt(1+C^2) from bridge content | **PASS** |
| T6 | `T6_no_signaling_consistency` | no-signaling: marginals invariant; the budget closes | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | classical ER=EPR quantitative; condition 2 split | **PASS** |

## The ER=EPR curve (bridge content -> Bell violation)

| beta | sin 2beta | C extracted | S_ent | CHSH | grid check | 2*sqrt(1+C^2) |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.0 | 0.0 | -0.0 | 2.0 | 2.0 | 2.0 |
| pi/12 | 0.5 | 0.5 | 0.2458 | 2.2361 | 2.2356 | 2.2361 |
| pi/8 | 0.7071 | 0.7071 | 0.4165 | 2.4495 | 2.4489 | 2.4495 |
| pi/6 | 0.866 | 0.866 | 0.5623 | 2.6458 | 2.645 | 2.6458 |
| pi/4 | 1.0 | 1.0 | 0.6931 | 2.8284 | 2.8284 | 2.8284 |

(bridge cut: CHSH = 2.0, entropy -0e+00; singlet fidelity at beta = pi/4: 1.0)

## Verdict

**CONFIGURATION_SPACE_EMERGES_FROM_THE_BRIDGE_SCHMIDT_WEIGHTS_ARE_BRIDGE_MODE_AMPLITUDES_CHSH_2_WITHOUT_2SQRT2_WITH_CLASSICAL_ER_EPR_QUANTITATIVE.** ENTANGLEMENT IS BRIDGE TOPOLOGY - MEASURED (the argument is in docs/configuration_space_emergence.md; this probe checks every step).

THE TARGET. Bell's theorem, enumerated exactly: local strategies cap at CHSH = 2; the quantum sector needs 2.8284. A 3-space-local classical field cannot cross that gap - the entangled sector must come from the bridge, the one nonlocal element BAM owns.

THE IDENTIFICATION. On the lattice, the bulk gluing delivers exact charge conjugation (purity 1.0), rigid phase locking (slope 1.0), the pi-holonomy Bell-state selector, and one-object readout - all destroyed by cutting the bridge.

THE EMERGENCE. One shared fiber read from two mouths embeds isometrically as a maximally correlated pair state; the derived transport T = i sigma_y (T^2 = -1, the Pin- sign) makes it THE SINGLET - bell.bulk_identity's postulated state, now derived (fidelity 1.000000), its E(a,b) = -cos(a-b) matching the module to 2e-16.

THE LAW. Schmidt weights = bridge-mode amplitudes (concurrence tracks the preparation to 2e-3); CHSH = 2*sqrt(1 + C^2) exactly across the sweep; bridge cut -> CHSH = 2.0000; symmetric bridge -> 2*sqrt(2). Marginals stay setting-independent (the budget: Bell correlation without a telegraph, #204). Classical ER=EPR is quantitative, and #198's condition 2 splits: the Bell-sector half discharged, the dynamical half sharply bounded.
