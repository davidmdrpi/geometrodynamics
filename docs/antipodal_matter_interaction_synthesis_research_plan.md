# Antipodal matter interaction ledger synthesis (PR #139)

This is the **capstone** of the antipodal matter-interaction arc — PRs #129–#138,
on the cavity operator of #116 and the measure of #115–#122. That arc built the
BAM matter interaction **layer by layer** on the antipodal horizon boundary
data: the boundary condition (#129), the spectrum (#130), the free propagator
(#135), the one-loop self-energy (#136), and the cubic and quartic vertices
(#137/#138). This PR introduces **no new physics**; it (a) re-verifies, in one
place, a keystone from each arc member (a cross-arc consistency check), (b) shows
the whole arc is governed by **two threads from a single postulate**, and (c)
lays out the honest epistemic ledger.

## The arc, layer by layer

| PR | layer |
|---|---|
| **#129** | the antipodal l-parity BC (even-l Neumann, odd-l Dirichlet), a unitary mirror, from `Y_l(−x) = (−1)^l Y_l` |
| **#130** | the spectrum: antipodal ⟹ real, undamped (stable matter); absorbing ⟹ complex ringdown |
| **#135** | the free propagator: the cavity resolvent with the antipodal data — reciprocal, unitary, parity-graded; a mode-sum over the stable modes |
| **#136** | the one-loop self-energy: finite real mass shift, lightest mode exactly stable (`Im Σ = 0`), unitarity preserved (no horizon-absorption width) |
| **#137** | the cubic vertex: angular selection rule `Σl` even (antipodal Z₂) + triangle; geometric shape derived, coupling input |
| **#138** | the quartic vertex + bounded interaction: same Z₂; `∫ψ⁴ > 0` ⟹ bounded below ⟹ stable vacuum = the #122 measure-convergence condition |

## The keystones, re-verified together

| PR | keystone | value |
|---|---|---|
| #129 | `Y_l` antipodal parity | `[1, −1, 1, −1] = (−1)^l` |
| #130 | antipodal fundamental (real) | `≈ 1.17` |
| #130 | absorbing fundamental (complex) | `≈ 1.89 − 1.16i` |
| #135 | kernel reciprocity `\|G−Gᵀ\|/\|G\|` | `~1e-14` |
| #136 | lightest-mode `Im Σ` (stable) | `≈ 0` |
| #138 | `∫ψ⁴` overlap (bounded vacuum) | `1.03 > 0` |

All mutually consistent in one run.

## Two threads, one postulate

The whole arc is two threads, both flowing from the single antipodal
identification (#128/#129):

  **Thread A — the antipodal Z₂ `(−1)^l`.** The C-swap inversion `x → −x` (#63)
  carries `Y_l → (−1)^l Y_l`. This one Z₂ fixes the boundary condition (#129),
  grades the exchange kernel (#135), and selects the cubic (#137) and quartic
  (#138) vertices (`Σl` even). One parity threading the entire interaction
  structure.

  **Thread B — unitarity / stability.** The antipodal BC is a unitary mirror
  (#129), giving a real, stable spectrum (#130), a unitary reciprocal propagator
  (#135), a unitarity-preserving self-energy with an exactly-stable lightest mode
  (#136), and a bounded-below interacting vacuum (#138, = the #122 measure
  convergence). BAM matter is stable at every order because the throat reflects
  antipodally.

|  | Thread A (Z₂ selection) | Thread B (unitarity/stability) |
|---|---|---|
| #129 BC | even-l N / odd-l D (parity) | unitary mirror (zero flux) |
| #130 | — | real stable spectrum |
| #135 | kernel parity grading | reciprocal unitary propagator |
| #136 | — | stable lightest mode |
| #137/#138 | Σl-even vertex selection | bounded vacuum (#122) |

The two threads are **not independent**: the real l-parity boundary condition IS
both the Z₂ grading (Thread A) and the unitary mirror (Thread B). **One
postulate, two faces.**

## The epistemic ledger

  - **Derived (given the antipodal BC):** the Z₂ selection structure of the
    propagator and the vertices; the unitary, reciprocal, real-pole propagator;
    the unitarity-preserving self-energy and stable lightest mode; the
    bounded-below stable vacuum.
  - **Postulated (BAM's axiom):** the antipodal identification itself (#128) —
    the throat glued antipodally, not absorbing. The arc shows it is
    **self-consistent** (a unitary, stable, bounded interacting theory), **not
    forced**.
  - **Input:** the coupling magnitudes `λ_3, λ_4` (the sign `λ_4 > 0` is required
    by #122 convergence; the magnitudes are not derived from S_BAM).
  - **Open:** the S_BAM generation of the vertices; higher loops / higher
    vertices; the bulk-scale (#133) and flavor (#134) residuals.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | synthesise the antipodal matter-interaction arc (#129–#138) |
| T2 | arc layer by layer | #129 BC → #130 spectrum → #135 propagator → #136 Σ → #137/#138 vertices |
| T3 | keystones together | cross-arc consistency (all reproduce in one run) |
| T4 | Thread A | antipodal Z₂ `(−1)^l`: BC, kernel grading, vertex selection |
| T5 | Thread B | unitarity/stability: spectrum → propagator → Σ → vacuum |
| T6 | one postulate | the real l-parity BC is both the Z₂ grading and the mirror |
| T7 | epistemic ledger | derived / postulated / input / open |
| T8 | assessment | `ANTIPODAL_MATTER_INTERACTION_SYNTHESIS_Z2_SELECTION_AND_UNITARY_STABILITY` |

## Established and open

  - **Established (BAM-native):** the antipodal matter-interaction arc is two
    threads from one postulate — the antipodal Z₂ `(−1)^l` selecting the
    propagator and vertices, and the unitary mirror giving a stable spectrum,
    propagator, self-energy, and bounded vacuum; the keystones are mutually
    consistent.

  - **Does not / open:** a synthesis/consistency capstone — it re-verifies and
    organises the arc; it does **not** add new derivations, remove any open item,
    or strengthen the antipodal postulate from "self-consistent" to "forced".
    The couplings are input; the S_BAM vertex generation, higher loops/vertices,
    and the bulk-scale (#133) and flavor (#134) residuals stand.

## Cross-references

  - `docs/null_throat_boundary_conditions_research_plan.md` — #129 (BC)
  - `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #130 (spectrum)
  - `docs/antipodal_horizon_exchange_kernel_research_plan.md` — #135 (propagator)
  - `docs/antipodal_kernel_one_loop_self_energy_research_plan.md` — #136 (self-energy)
  - `docs/cubic_vertex_ledger_research_plan.md` — #137 (cubic vertex)
  - `docs/quartic_vertex_bounded_interaction_research_plan.md` — #138 (quartic + bounded)
  - `docs/bam_factorized_sector_sum_research_plan.md` — #122 (measure convergence)
  - `docs/charge_conjugation_swap_research_plan.md` — #63 (C-swap inversion)

## Run

```
python -m experiments.closure_ledger.antipodal_matter_interaction_synthesis_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_antipodal_matter_interaction_synthesis_probe/`.
Expected verdict:
`ANTIPODAL_MATTER_INTERACTION_SYNTHESIS_Z2_SELECTION_AND_UNITARY_STABILITY`, 8/8 PASS.
