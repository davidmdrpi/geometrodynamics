# Sharpening the meV-scale neutrino predictions (PR #111)

The BAM neutrino sector (PRs #87–#96) predicts a light, normal-ordered,
Majorana spectrum (generations = cavity overtones; `m_ν ∝ m_D`; `c₁ = 0 ⟹`
Majorana; anarchic phases). PR #96 left the headline meV-scale observable as
a **band** — `Σm_ν ≈ 59–65 meV` — because the lightest mass `m₁` was an open
absolute-scale residual (~few meV, PR #90). This PR **sharpens** the
meV-scale predictions by: (i) updating the oscillation inputs to the latest
global fit (**NuFIT 6.0**, 2024), (ii) using the 2025 **DESI DR2 + CMB**
cosmology to corner `m₁` against the normal-ordering floor, and (iii)
reading off the full pinned spectrum — `Σm_ν`, the three masses, the β-decay
effective mass `m_β`, and the 0νββ effective mass `m_ββ` — with an honest
statement of which are reachable.

## The cornering

NuFIT 6.0 (NO) fixes two of the three masses outright:

```
√Δm²_21 = 8.65 meV,   √Δm²_31 = 50.34 meV,
Σm_ν|_floor (m₁→0) = 59.0 meV.
```

DESI DR2 + CMB (2025) bounds `Σm_ν ≲ 60–64 meV (95%)`, which pins `m₁` to
the low end of "few meV":

```
m₁ ≲ 3 meV   ⟹   Σm_ν ∈ [59.0, 62.6] meV,
```

a sharpening of the PR #96 band (59–65) toward the floor. The sharpest
single-number BAM statement is the hierarchical limit `m₁ → 0`, the natural
endpoint of the light cavity-overtone spectrum and the point cosmology is
cornering us toward.

## The pinned spectrum (hierarchical limit)

```
m₁ ≈ 0  (≲ 3 meV),   m₂ = 8.65 meV,   m₃ = 50.34 meV.
```

From this the two laboratory effective masses follow:

  - **β-decay (KATRIN):** `m_β = √(Σ|U_ei|² m_i²) ≈ 8.8 meV` (rising to
    ~9.3 meV at `m₁ = 3 meV`).
  - **0νββ:** `m_ββ = |Σ U_ei² m_i|`. Crucially, in normal ordering the
    contributions **cannot fully cancel** at `m₁ → 0`, because
    `s12²c13² m₂ = 2.60 meV > s13² m₃ = 1.10 meV`, leaving a **nonzero
    floor**

    ```
    m_ββ ∈ [1.5, 3.7] meV   (m₁ → 0, over the Majorana phases),
    ```

    widening to ~`[0, 5.9]` meV as `m₁ → 3 meV`. So BAM predicts a small but
    nonvanishing 0νββ rate in the hierarchical limit.

## Reachability (the honest other half)

| observable | BAM (NO, light) | testable by | reach vs BAM |
|---|---|---|---|
| **Σm_ν** | **59.0–62.6 meV** | **DESI DR2 + CMB** | **at the floor NOW** |
| m_β | 8.8–9.3 meV | KATRIN / Project 8 | ~4–5× below |
| m_ββ | 1.5–3.7 meV (floor) | LEGEND-1000 / nEXO | ~3–10× below |

The spectrum is pinned, but only `Σm_ν` is near-term testable — and it is
being tested at the floor right now. The Majorana and absolute-scale
predictions (`m_ββ`, `m_β`) are sharp yet below foreseeable sensitivity — an
honest readout, not a promise of imminent discovery.

## Falsifiers (sharpened)

  - A robust cosmological `Σm_ν < 59.0 meV` (below the NO floor) excludes
    normal ordering ⟹ BAM fails (and would clash with the oscillation `Δm²`
    themselves — a deep consistency break). **Flag:** some 2025 DESI + CMB
    analyses already prefer central `Σm_ν` at or below the floor (the data
    wants less lensing than the floor supplies); if that hardens it is
    tension for ALL normal-ordered models, BAM included.
  - A quasi-degenerate detection (`Σm_ν ≳ 100 meV`) contradicts the BAM
    light scale.
  - A 0νββ signal implying `m_ββ ≫ 4 meV` with confirmed normal ordering
    would sit above the BAM hierarchical band.

## Tests

| # | test | finding |
|---|---|---|
| T1 | setup | sharpen the Σm_ν band → full pinned spectrum |
| T2 | NuFIT 6.0 | `√Δm²_21` 8.65, `√Δm²_31` 50.34 meV ⟹ NO floor 59.0 |
| T3 | cornering | DESI DR2 ⟹ `m₁ ≲ 3 meV` ⟹ `Σm_ν ∈ [59.0, 62.6]` |
| T4 | pinned spectrum | `m = (≲3, 8.65, 50.34) meV` |
| T5 | lab masses | `m_β` 8.8–9.3; `m_ββ` nonzero floor [1.5, 3.7] meV |
| T6 | reachability | only `Σm_ν` near-term testable (DESI); `m_β`, `m_ββ` below |
| T7 | falsifiers/scope | sharpened falsifiers + 2025 sub-floor tension flag |
| T8 | assessment | `NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE` |

## Established and open

  - **Established (BAM-native + data):** updating to NuFIT 6.0 and cornering
    `m₁` with DESI DR2 pins the meV-scale spectrum — `m = (≲3, 8.65, 50.34)
    meV`, `Σm_ν ∈ [59.0, 62.6]` (sharpened toward the floor), `m_β ≈
    8.8–9.3 meV`, `m_ββ` a nonzero floor `[1.5, 3.7] meV`. Only `Σm_ν` is
    near-term testable (DESI, at the floor now).

  - **Open:** `m₁` within its cornered band (0–3 meV) — the absolute-scale
    residual; and the anarchic Majorana phases, which set the exact `m_ββ`
    within the floor band (the universal flavor puzzle).

## Cross-references

  - `docs/cosmological_sigma_mnu_research_plan.md` — PR #96, the `Σm_ν`
    band this PR sharpens.
  - `docs/cp_majorana_phase_research_plan.md` — PR #93/#94, the anarchic
    Majorana phases that set `m_ββ` within the floor band.
  - `docs/seesaw_scale_nucleation_compliance_research_plan.md` /
    `docs/majorana_bounce_action_research_plan.md` — PRs #87/#88, the
    seesaw + Majorana origin.
  - `docs/neutrino_quadrant_suppression_research_plan.md` — PR #90, the
    light absolute-scale residual now cornered.

## Run

```
python -m experiments.closure_ledger.neutrino_mev_scale_sharpening_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_neutrino_mev_scale_sharpening_probe/`.
Expected verdict:
`NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE`, 8/8 PASS.
