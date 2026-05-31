# Heavy Möbius baryon: a sharper LHCb / Belle II search table (PR #114)

PRs #103/#109/#110 predicted the heavy Möbius/hybrid baryon (ground heavy
baryon + `2√σ`), worked out its twist-unwinding decays and the hybrid
selection rule, and compiled the non-orientable sector into one note. This
PR **converts that physics into a sharper, actionable LHCb / Belle II search
table**: for each target, the concrete reconstruction chain, the facility
and production mode, the primary discovery channel with its kinematic
handle, the discriminator, the constraint status, and a tiered priority.
Every mass is a pushforward of the single scale `√σ`. The compiled table is
committed standalone as `docs/heavy_mobius_baryon_search_table.md`.

## The new quantitative handle: the dipion endpoint

The preferred twist-unwinding channel is `Λ_Q(ππ)_S`. Its dipion
invariant-mass endpoint — the maximum `m(ππ)`, reached when the heavy baryon
recoils at rest — is

```
m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV,
```

**flavor-independent**: the same 849 MeV endpoint above the charm and the
bottom ground baryon, with the spectrum peaking **high** (coherent isoscalar
S-wave, like `ψ(2S) → J/ψ ππ`). This is sharper than a Q-value: it is a
fixed edge in a directly-plotted observable, identical for c and b — a
single overlay that tests the whole framework.

## The search table (tiered)

| Tier | Target (mass MeV) | Reconstruction chain | Facility / production | Constraint status |
|---|---|---|---|---|
| 1 | Λ_c Möbius (3135) | `Λ_c⁺π⁺π⁻, Λ_c⁺ → pK⁻π⁺` | LHCb (prompt) + Belle II | above Λ_c(2940) |
| 1 | Λ_b Möbius (6469) | `Λ_b⁰π⁺π⁻, Λ_b⁰ → Λ_c⁺π⁻` | LHCb · from b-decays | above Λ_b(6152) |
| 2 | Ω_b Möbius (6894) | `Ω_b⁻π⁺π⁻, Ω_b⁻ → Ω_c⁰π⁻` | LHCb · from b-decays (rare) | unexplored |
| 2 | Ξ_cc Möbius (4471) | `Ξ_cc⁺⁺π⁺π⁻, Ξ_cc⁺⁺ → Λ_c⁺K⁻π⁺π⁺` | LHCb · prompt (rare) | unexplored |
| 3 | Ω_c Möbius (3544) | `Ω_c⁰π⁺π⁻` | LHCb · prompt / b-decays | above 2017 Ω_c (≤3120) |

  - **Tier 1 — the discovery pair:** highest yield × the cross-flavor
    clincher. Λ_c has the statistics-richest channel (golden `pK⁻π⁺`); Λ_b
    is its cross-flavor partner. Together they carry the same dipion endpoint
    (849 MeV) and Q-values (`ππ` 569, `η` 301 MeV) above both — a
    simultaneous two-channel fit.
  - **Tier 2 — the clean frontier:** Ξ_cc and Ω_b have no measured
    excitation spectrum at all (PR #103's "freest of the free"), so a clean
    bump = discovery — but doubly-heavy / Ω_b production is rare
    (rate-limited).
  - **Tier 3 — calibratable:** Ω_c sits above the well-mapped 2017 LHCb Ω_c
    excitations (≤3120), which calibrate the search but mean it is not virgin
    territory.

Priority is a transparent sum of four factors (each 0–2): production yield,
reconstruction cleanliness, freedom from existing-data confusion
(unexplored), and cross-flavor leverage.

## The discriminators (per channel)

  - **Suppressed single-π-to-ground** — the hybrid selection rule forbids
    decay to the ground baryon + one S-wave π; an ordinary radial excitation
    does the opposite. The branching PATTERN distinguishes them.
  - **Dipion endpoint + shape** — a fixed 849 MeV edge, peaking high; the
    same overlay for c and b.
  - **Cross-flavor Q-match** — `Λ_Q ππ` at 569 and `Λ_Q η` at 301 MeV
    identical for charm and bottom; a joint constraint, not a single bump.

## Tests

| # | test | finding |
|---|---|---|
| T1 | scope | convert #109/#110 → tiered actionable table |
| T2 | targets/chains | five targets, golden exclusive reconstruction modes |
| T3 | facility | LHCb prompt/b-decay; Belle II charm threshold |
| T4 | dipion endpoint | `m(ππ)_max = 2√σ ≈ 849 MeV`, flavor-independent (new handle) |
| T5 | discriminators | selection rule; endpoint; cross-flavor Q-match |
| T6 | priority | Tier 1 Λ_c+Λ_b; Tier 2 Ξ_cc+Ω_b; Tier 3 Ω_c |
| T7 | honest scope | masses ±band, broad, BFs/J^P open |
| T8 | assessment | `HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED` |

## Established and open

  - **Established (BAM-native):** PRs #109/#110 converted into a tiered,
    actionable LHCb / Belle II search table — reconstruction chains,
    facility/production, the flavor-independent dipion endpoint (849 MeV) as
    the new sharp handle, per-channel discriminators, and explicit priority.

  - **Open (unchanged):** exact masses (±lattice hybrid gap ~0.8–1.3 GeV),
    broad widths (~tens–150 MeV ⟹ amplitude analyses), absolute branching
    fractions, and the `J^P`. A prioritization deliverable, not new physics.

## Cross-references

  - `docs/heavy_mobius_baryon_search_table.md` — the compiled table (the
    deliverable).
  - `docs/heavy_mobius_baryon_decay_research_plan.md` — PR #109, the
    twist-unwinding decays and the hybrid selection rule.
  - `docs/bam_nonorientable_experimental_note.md` /
    `docs/nonorientable_experimental_note_research_plan.md` — PR #110, the
    sector note this sharpens for the heavy-baryon channels.
  - `docs/heavy_mobius_baryon_research_plan.md` — PR #103, the masses.

## Run

```
python -m experiments.closure_ledger.heavy_mobius_baryon_search_table_probe
```

Writes `probe.json`, `probe.md`, and the standalone `search_table.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_heavy_mobius_baryon_search_table_probe/`.
Expected verdict: `HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED`, 8/8 PASS.
