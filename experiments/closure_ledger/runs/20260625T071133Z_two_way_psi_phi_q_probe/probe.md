# Two-way ψ–Φ–q evolution: the self-consistent matter–metric–order system (PR #179)

**Run:** 2026-06-25T07:11:33+00:00

Closes PR #178's one-way coupling into the full TWO-WAY self-consistent system of three co-evolving fields — the matter wave `ψ`, the gravitational potential `Φ`, and the throat-order field `q` — all from one energy functional. *(QFT on the classical throat, not quantum gravity.)*

- **System**: `∂_τψ=½∇²ψ−Φψ+½gq²ψ; ∂_τq=κ∇²q−(a₀−g|ψ|²)q−λq³; ∇²Φ=4πG(|ψ|²+μq²)`
- **Self-consistent**: the coupled flow converges to a throat-soliton fixed point
- **Back-reaction**: throat deepens the well (~5%) and densifies the core (~12%); two-way
- **Saturation vs runaway**: stable bound soliton (sub-critical μ) vs runaway (super-critical μ)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | close #178's one-way coupling into two-way ψ–Φ–q | **PASS** |
| T2 | `T2_coupled_system_one_functional` | the coupled system from one energy functional; four channels | **PASS** |
| T3 | `T3_self_consistency_converges` | self-consistency: the coupled flow converges to a fixed point | **PASS** |
| T4 | `T4_two_way_back_reaction` | back-reaction: throat deepens the well / densifies the core | **PASS** |
| T5 | `T5_saturation_vs_collapse` | saturation (stable soliton) vs runaway (q self-gravity) | **PASS** |
| T6 | `T6_threshold_continuity` | threshold continuity: sub = pure SN; super = throat-soliton | **PASS** |
| T7 | `T7_honest_scope` | honest scope (weak-field, spherical, effective constants) | **PASS** |
| T8 | `T8_assessment` | TWO_WAY_PSI_PHI_Q_SELF_CONSISTENT_THROAT_SOLITON | **PASS** |

## The two-way back-reaction (M = 3, super-threshold)

| quantity | pure Schrödinger–Newton (q=0) | two-way ψ–Φ–q | change |
|---|---:|---:|---:|
| well depth Φ(0) | -3.0274 | -3.1779 | 4.97% deeper |
| core density ρ_peak | — | — | 13.44% denser |

## Saturation vs collapse

- **stable** (sub-critical self-gravity): max\|q\| [0.046, 0.218, 0.414, 0.471, 0.487, 0.493, 0.494] — plateaus (a self-consistent throat-soliton)
- **collapse** (super-critical self-gravity): max\|q\| → 31.17, Φ(0) → -252.0, residual 7e+03 — no weak-field fixed point

## Verdict

**TWO_WAY_PSI_PHI_Q_CONVERGES_TO_A_SELF_CONSISTENT_THROAT_SOLITON_BACK_REACTION_REAL.** LOOP CLOSED — A SELF-CONSISTENT TWO-WAY THROAT-SOLITON. Closing #178's one-way coupling into the full ψ–Φ–q system, the throat-order field back-reacts both ways.

ONE FUNCTIONAL. The whole coupled system — ∂_τψ = ½∇²ψ − Φψ + ½g q²ψ, ∂_τq = κ∇²q − (a₀−g|ψ|²)q − λq³, ∇²Φ = 4πG(|ψ|²+μq²) — descends from one energy functional, so the four channels (ψ↔Φ, ψ→q, q→ψ, q→Φ) are consistently coupled (the ordering and binding terms share the same g).

SELF-CONSISTENT. The coupled flow converges: the energy plateaus (ΔE = -6.00e-04) and the q residual drops to 1.1e-04 — a self-consistent throat-soliton exists.

TWO-WAY BACK-REACTION. Versus the pure Schrödinger–Newton soliton, the ordered throat core deepens the well by 5.0% and densifies the core by 13.4% — the throat traps the wave, which concentrates it, which strengthens the order.

SATURATION vs RUNAWAY. The binding feedback saturates (λq⁴) into a stable bound soliton (|q| plateaus); q's self-gravity, pushed past a coupling threshold, drives runaway collapse (|q| climbs without bound) — the strong-field endpoint.

CONTINUOUS. Sub-threshold the order field vanishes and the system reduces exactly to the Schrödinger–Newton soliton of #176/#177; the #176 → #178 → #179 arc is one continuous system.

SCOPE. Weak-field, semi-dynamical, spherically reduced; the constants are effective (the structure is the result); the stable soliton is sub-critical and the strong-field runaway is for full NR.
