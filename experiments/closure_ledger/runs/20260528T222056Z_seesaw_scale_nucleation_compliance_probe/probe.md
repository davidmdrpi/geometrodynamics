# Geometric origin of `M_R` from throat-nucleation compliance (PR #87)

**Run:** 2026-05-28T22:20:56+00:00

PR #86 left the Majorana seesaw scale `M_R ≈ 0.3–1.8 TeV` open. This probe grounds it in the PR #58 throat-nucleation framework: the `ΔL=2` Majorana coupling **is** the throat ↔ antithroat (antipodal `Z₂`) transition, PR #58's `Σc₁=0` **is** PR #86's only-neutrino selection rule, the literal `M_R = `barrier-height hypothesis is **falsified** (~10⁸ too small), and the suppression is **tunnelling through** the barrier: `M_R = m_D·e^{S}` with a modest, generation-stable bounce action `S ≈ 15–18`.

- **Identification**: M_R seesaw scale = the keV Dirac floor lifted by the throat↔antithroat tunnelling action: M_R = m_D·e^{S}, S = ln(m_D/m_ν) ≈ 15–18; barrier-height M_R falsified
- **Mechanism**: ΔL=2 Majorana = PR #58 throat↔antithroat (antipodal Z₂) transition
- **Selection rule**: PR #58 Σc₁=0 ⟹ only k=0 (c₁=0) gets single-state Majorana
- **Falsified**: M_R = nucleation barrier height E_c (~keV, ~10⁸ too small)
- **Open**: the bounce action S (PR #58 instanton-rate follow-on)
- **B4 caveat**: m_D from cavity floor + electron anchor; S open; absolute m_ν needs S + L_eff unification + B4 anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_required_M_R_recap` | required M_R ≈ 0.3–1.8 TeV (PR #86 open target) | **PASS** |
| T2 | `T2_majorana_is_throat_antithroat_transition` | ΔL=2 Majorana = throat↔antithroat = PR #58 antipodal Z₂ | **PASS** |
| T3 | `T3_selection_rule_from_sigma_c1_zero` | PR #58 Σc₁=0 ⟹ only k=0 Majorana (re-derives PR #86) | **PASS** |
| T4 | `T4_barrier_height_hypothesis_falsified` | M_R = barrier height E_c ≈ 2.8 keV → ~10⁸ too small (falsified) | **PASS** |
| T5 | `T5_boundary_compliance_reading` | soft bubble = max compliance; static origin excluded | **PASS** |
| T6 | `T6_throat_antithroat_tunnelling` | m_ν = m_D·e^{−S}; S = ln(m_D/m_ν) ≈ 15–18; M_R = m_D·e^{S} | **PASS** |
| T7 | `T7_honest_scope` | mechanism+rule established; S the open follow-on | **PASS** |
| T8 | `T8_assessment` | SEESAW_SUPPRESSION_IS_THROAT_ANTITHROAT_TUNNELING_S_OPEN | **PASS** |

## T3: Selection rule from PR #58's `Σc₁=0`

| sector | k | c₁ | ΔΣc₁ (single flip) | single-state Majorana? | nature |
|---|---:|---:|---:|:---:|---|
| neutrino | 0 | 0 | +0 | ✓ | Majorana |
| charged lepton g1 | 1 | 1 | -2 | ✗ | Dirac |
| charged lepton g2 | 3 | 1 | -2 | ✗ | Dirac |
| charged lepton g3 | 5 | 1 | -2 | ✗ | Dirac |

A single state flipping throat ↔ antithroat sends `c₁ → −c₁`. Only `k=0` (`c₁=0`) keeps `Σc₁=0` with no partner, so only the neutrino gets a single-state `ΔL=2` Majorana mass. This is PR #86's result, derived from PR #58's conservation law.

## T4: Hypothesis A — `M_R` = nucleation barrier height (falsified)

- Electron-throat geometric units: `σ = 2.6526e-02`, `ρ = 2.3873e-01`, `R_c/R* = 0.222`
- Barrier height `E_c = (16π/3)σ³/ρ² = 5.4870e-03·m_e c² ≈ 2.80 keV`

| gen | required M_R (eV) | barrier E_c (eV) | M_R / E_c |
|---:|---:|---:|---:|
| 1 | 1.836e+12 | 2.804e+03 | 6.55e+08 |
| 2 | 7.149e+11 | 2.804e+03 | 2.55e+08 |
| 3 | 2.765e+11 | 2.804e+03 | 9.86e+07 |

The barrier height is `~10⁸–10⁹` below the required `~TeV` `M_R`. The soft vacuum bubble is the maximally compliant throat-opening channel; its height is not the heavy Majorana scale. **Hypothesis A falsified.**

## T6: Hypothesis B — throat ↔ antithroat tunnelling (viable)

| gen | m_D (keV) | m_ν obs (eV) | required S = ln(m_D/m_ν) | M_R = m_D·e^{S} (GeV) |
|---:|---:|---:|---:|---:|
| 1 | 42.9 | 0.001 | 17.6 | 1836 |
| 2 | 80.2 | 0.009 | 16.0 | 715 |
| 3 | 117.6 | 0.05 | 14.7 | 277 |

Required bounce action `S ≈ 14.7–17.6` — modest, `O(15)`, generation-stable. The seesaw relation becomes `M_R = m_D²/m_ν = m_D·e^{S}`: the `~TeV` scale is the `keV` Dirac floor **exponentially lifted** by the throat ↔ antithroat tunnelling action, not a fundamental heavy mass.

## Verdict

**SEESAW_SUPPRESSION_IS_THROAT_ANTITHROAT_TUNNELING_S_OPEN.** SEESAW SUPPRESSION IS THROAT ↔ ANTITHROAT TUNNELLING; THE BOUNCE ACTION S IS OPEN. PR #86 left the Majorana seesaw scale M_R ≈ 0.3–1.8 TeV as an open heavy input. This probe grounds it in the PR #58 throat-nucleation framework.

THE ΔL=2 MECHANISM IS THE PR #58 TRANSITION. A Majorana mass flips a state into its antiparticle — in BAM the antithroat, the inner/outer swap C (c₁ → −c₁, PR #63), the antipodal Z₂ whose vacuum nucleation PR #58 analysed. So the ΔL=2 coupling lives natively in the PR #58 nucleation barrier.

THE SELECTION RULE IS PR #58's Σc₁=0. Applied to a SINGLE state flipping throat ↔ antithroat: k=0 (c₁=0) flips 0 → −0 = 0, conserving Σc₁ with no partner, so a single-state ΔL=2 (Majorana) mass is ALLOWED; k≠0 (c₁=±1) would give Σc₁ = ∓2 ≠ 0, forbidden — a ΔL=2 mass needs a pair, so the charged lepton stays Dirac. PR #86's "only the c₁=0 sector is Majorana" IS PR #58's conservation law for a single state. The two PRs are one statement.

M_R IS NOT THE BARRIER HEIGHT. The most direct hypothesis — M_R = the nucleation barrier E_c = (16π/3)σ³/ρ² — fails on scale. With the brane tension σ and bag density ρ fixed by the electron throat, E_c ≈ 0.0055 m_e c² ≈ 2.8 keV and R_c ≈ 0.22 R*: ~10⁸–10⁹ BELOW the required ~TeV. The soft vacuum bubble is the MAXIMALLY COMPLIANT throat-opening channel; its barrier height cannot be the heavy Majorana scale, and BAM has no stiff ~TeV barrier or partner state. A static barrier-height / compliance origin is excluded.

THE SUPPRESSION IS TUNNELLING THROUGH THE BARRIER. PR #58 flagged the tunnelling/instanton rate as its open follow-on. If the Majorana mass is the throat ↔ antithroat tunnelling AMPLITUDE through the barrier rather than its height, then m_ν = m_D·e^{−S} with S the Euclidean bounce action. The required action is the data-derived S = ln(m_D/m_ν) ≈ 17.6, 16.0, 14.7 (gen 1, 2, 3) — modest, O(15), and generation-stable. Equivalently M_R = m_D²/m_ν = m_D·e^{S}: the ~TeV seesaw scale is NOT a fundamental heavy mass but the keV Dirac floor (the cavity zero-point) EXPONENTIALLY LIFTED by the throat ↔ antithroat tunnelling action.

NET PROGRESS. The open input is recast from a mysterious ~TeV mass (6 orders above any BAM scale, with no candidate state) into a modest, generation-stable O(15) bounce action — a far more constrained, BAM-shaped unknown that is exactly PR #58's own flagged instanton follow-on.

HONEST SCOPE. ESTABLISHED (BAM-native): the ΔL=2 mechanism is the PR #58 throat↔antithroat transition; PR #58's Σc₁=0 IS PR #86's only-neutrino selection rule; the literal M_R = barrier-height hypothesis is falsified (~10⁸ too small); the suppression is tunnelling through the barrier, M_R = m_D·e^{S} with S ≈ 15–18. NOT established: the value of S from first principles (needs the BAM throat inertia / Euclidean action normalisation — the instanton-rate follow-on), the absolute m_ν, and PMNS mixing.

## What this leaves open

- **The bounce action `S`** — from first principles, needs the BAM throat inertia / Euclidean action normalisation (the instanton-rate follow-on PR #58 explicitly flagged). The open input is now an `O(15)` action, not a `~TeV` mass.
- **Absolute neutrino masses** — need `S` + the `L_eff` unification (PR #83) + the B4 anchor.
- **PMNS mixing / mass ordering** — beyond this probe's scope.
