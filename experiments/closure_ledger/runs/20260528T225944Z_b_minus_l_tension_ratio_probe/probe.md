# ΔL=2 / B−L throat-tension ratio: derive or constrain (PR #89)

**Run:** 2026-05-28T22:59:44+00:00

PR #88 localised the neutrino mass's open input to the ΔL=2 (B−L) throat-tension ratio `t ≈ 6–12`. This probe constrains `t` from BAM structure: it is bracketed, parameter-free, by the **closure quantum `2π`** (minimal orientation reversal) and the **winding action `k_5√(2π) = √β_lepton`** (full winding) — `t ∈ [6.28, 12.53]`, exactly PR #88's required band.

- **Identification**: ΔL=2/B−L tension ratio t bracketed by the closure quantum 2π (lower) and the winding action k_5√(2π)=√β_lepton (upper): t ∈ [6.28, 12.53], matching PR #88's required 6–12
- **Lower bound**: 2π (closure quantum, minimal orientation reversal)
- **Upper bound**: k_5√(2π) = √β_lepton (full throat winding)
- **Residual**: where in the window = the boundary compliance ε
- **B4 caveat**: t constrained parameter-free to [2π, k_5√(2π)]; unique value needs ε; bounce-normalisation model dependence

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_t_epsilon_degeneracy` | recap t≈6–12; (t,ε) degeneracy t²·L*(ε)=S_req/P0 | **PASS** |
| T2 | `T2_b_minus_l_global_closure_enhancement` | t = B−L global-closure enhancement (orientation reversal is global) | **PASS** |
| T3 | `T3_lower_bound_closure_quantum_2pi` | lower bound = closure quantum 2π ≈ 6.28 | **PASS** |
| T4 | `T4_upper_bound_winding_action_sqrt_beta` | upper bound = winding action k_5√(2π) = √β ≈ 12.53 | **PASS** |
| T5 | `T5_window_brackets_required_t` | window [2π, k_5√(2π)] brackets required 6–12 | **PASS** |
| T6 | `T6_residual_is_compliance` | residual = compliance ε; m_charged/m_D ≈ √β cross-check | **PASS** |
| T7 | `T7_honest_scope` | t constrained parameter-free; residual = ε | **PASS** |
| T8 | `T8_assessment` | B_MINUS_L_TENSION_BRACKETED_BY_CLOSURE_AND_WINDING_ACTIONS | **PASS** |

## The closure-to-winding window

| bound | BAM scale | value | physical meaning |
|---|---|---:|---|
| lower | closure quantum `2π` | 6.283 | single great-circle orientation reversal (minimal global flip) |
| upper | winding action `k_5√(2π) = √β_lepton` | 12.533 | full throat winding to the antipode |

PR #88's required `t` (for `S=16` at sane compliance `ε ∈ [1e-6, 1e-2]`) is `[6.41, 12.05]` — **inside** the BAM window `[6.28, 12.53]`. The `6–12` band was not a fit but the BAM closure-to-winding window.

## Where in the window? → the compliance ε

| t (window edge) | BAM scale | compliance ε |
|---:|---|---:|
| 6.283 | closure quantum 2π | 5.85e-07 |
| 12.533 | winding action k_5√(2π) | 1.33e-02 |

Fixing `t` fixes `ε`: the closure-quantum edge is near-rigid (`ε≈6e-7`), the winding edge more compliant (`ε≈1.3e-2`). Cross-check: the winding/cavity-floor mass ratio `m_charged/m_D ≈ 11.9 ≈ √β_lepton` lands at the winding edge, consistent with the flip borrowing the winding channel.

## Verdict

**B_MINUS_L_TENSION_BRACKETED_BY_CLOSURE_AND_WINDING_ACTIONS.** THE ΔL=2 / B−L TENSION RATIO IS BRACKETED BY THE CLOSURE QUANTUM AND THE WINDING ACTION. PR #88 localised the neutrino mass's open input to a single dimensionless number — the ΔL=2 (B−L) throat tension must be a factor t ≈ 6–12 stiffer than the EM-throat tension. This probe constrains t from BAM structure.

THE (t, ε) DEGENERACY. The bounce S = t²·P0·L*(ε) means the real constraint is t²·L*(ε) = S_req/P0: the tension ratio t and the boundary compliance ε trade off. Pinning t needs an independent argument.

t IS A B−L GLOBAL-CLOSURE ENHANCEMENT. The EM-throat tension is a LOCAL surface tension (σ·Area, PR #56). The ΔL=2 Majorana flip reverses the throat's orientation (c₁→−c₁, PR #63) — a GLOBAL operation on S³, impossible locally — so its tension is the local tension times a global-closure factor. t is that factor: a B−L-breaking enhancement, not a free coupling. BAM has exactly two fundamental action scales for such a closure, and they bracket t.

LOWER BOUND = CLOSURE QUANTUM 2π. The cheapest global orientation reversal is a single great-circle traversal of the antipodal identification: one closure quantum 2π = action_base (the same 2π in Hopf holonomy, throat dwell, the PR #74 loop measure, the PR #83 throat closure quantum). You cannot pay less than one closure quantum for a global flip, so t ≥ 2π ≈ 6.28.

UPPER BOUND = WINDING ACTION k_5√(2π). The most expensive route is a full throat winding through the k_5 structure to reach the antipode: the winding action √β_lepton = √(k_5²·2π) = k_5·√(2π) ≈ 12.53 (PR #71). Beyond a full winding there is no costlier lepton-sector deformation, so t ≤ k_5√(2π).

THE WINDOW MATCHES THE REQUIREMENT. Hence t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53] — bracketed, parameter-free, by the closure quantum and the winding action, the two most basic BAM action scales. This is exactly PR #88's required t ≈ 6–12 (the computed band [6.41, 12.06] at sane compliance ε∈[1e-6,1e-2] sits INSIDE the window): the band was not a fit but the BAM closure-to-winding window. The winding/cavity-floor mass ratio m_charged/m_D ≈ 11.9 ≈ √β_lepton independently lands at the winding edge, consistent with the flip borrowing the winding channel.

WHERE IN THE WINDOW = THE COMPLIANCE. Fixing t in the window fixes the compliance ε: t=2π ⟹ ε≈6e-7 (a near-rigid throat), t=k_5√(2π) ⟹ ε≈1.3e-2 (more compliant). The neutrino's actual t sits between these, fixed by how much winding the ΔL=2 flip requires — i.e. by the compliance. The open input is sharpened from "an O(10) tension ratio" to "where in the [2π, k_5√(2π)] window," a single residual number.

PROGRESSIVE LOCALISATION. ~TeV mass (PR #86) → O(15) action S (PR #87) → O(10) tension ratio t (PR #88) → the BAM closure-to-winding window [2π, k_5√(2π)] + a compliance residual (PR #89).

HONEST SCOPE. ESTABLISHED (BAM-native): t is a B−L-breaking global-closure enhancement; it is bracketed, parameter-free, by the closure quantum 2π (minimal orientation reversal) and the winding action k_5√(2π) (full winding) — exactly PR #88's required 6–12. NOT established: a UNIQUE t (the (t,ε) degeneracy remains — the window edges map to ε ≈ 6e-7 … 1.3e-2), the compliance ε from first principles, and the bounce-normalisation model dependence (the t² scaling, the P0 prefactor). So this is a CONSTRAINT + identification, not a closed derivation; the residual open number is the compliance ε within the window.

## What this leaves open

- **A unique `t`** — the `(t, ε)` degeneracy remains; the window edges map to `ε ≈ 6×10⁻⁷ … 1.3×10⁻²`. The residual open number is now "where in the `[2π, k_5√(2π)]` window," i.e. the compliance `ε`.
- **The compliance `ε`** — a sub-throat length, not yet derived from the bulk.
- **Bounce-normalisation model dependence** — the `t²` scaling and the `P0` prefactor of the reduced bounce.
