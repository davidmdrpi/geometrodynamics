# BAM-specific baryonic exotics: classification + experimental constraints (PR #102)

PR #101 matched BAM's non-orientable (Möbius) topology to the OBSERVED
mesonic exotics (1-+ hybrids, multiquark zoo) and flagged the Möbius
baryon as a BAM-specific prediction. This probe classifies the
BAM-specific baryonic exotics and — as requested — identifies which
channels are experimentally constrained.

## The crucial subtlety: no exotic J^P for baryons

For mesons, the Möbius flux tube produced a SMOKING-GUN exotic quantum
number, `J^PC = 1-+`, forbidden to ordinary qq̄ by the `C = (−1)^{L+S}`
constraint. **Baryons have no such constraint.** A qqq baryon has
`P = (−1)^L`, `S ∈ {½, 3/2}`, and no good `C`, so every half-integer
`J^P` is reachable by an ordinary baryon — there is NO forbidden (exotic)
`J^P`. Consequently a BAM Möbius / hybrid baryon carries ORDINARY quantum
numbers: it is a **supernumerary** state, identifiable only by COUNTING
(an extra resonance beyond the qqq quark-model spectrum), not by a
smoking-gun `J^P`. This makes baryonic exotics the hardest to identify,
and — because they must hide inside an already densely-measured spectrum
— the MOST experimentally constrained corner of BAM's non-orientable
predictions (the opposite end from the unobserved glueballs of PR #100).

## The BAM-specific baryonic exotics

  - **Hybrid baryon** (excited Y-junction glue): the baryon analogue of
    the hybrid meson — the Y-junction flux in its first excited/twisted
    mode. Mass ≈ nucleon + `2√σ` ≈ 0.94 + 0.85 ≈ **1.79 GeV**; the Δ-based
    partner ≈ **2.08 GeV** (the same `2√σ` flux-tube quantum as PR #101).
  - **Möbius baryon** (non-orientable Y-network, `make_mobius_baryon`):
    the Z₂-twisted partner of a baryon; ordinary `J^P`, supernumerary.

Both land in the LIGHT N*/Δ* region — the densest, best-measured part of
the baryon spectrum.

## Experimental-constraint ranking of channels

| channel | data density | constraint on a BAM exotic |
|---|---|---|
| light N*/Δ* (< 2.5 GeV) | ~40 PDG states; "missing-resonance" problem | **strongest** |
| strange hyperons (Λ*, Σ*) | well measured | strong |
| charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*) | sparse (LHCb handful) | **weakest** (BAM freest) |

The lightest BAM baryonic exotics (~1.8–2.1 GeV) sit squarely in the
MOST-constrained light sector (`N(1710)`, `N(1875)`, `N(1900)`, …;
`Δ(1900)`, `Δ(1950)`, `Δ(2000)`, …).

## The test / tension

BAM's Möbius topology doubles the spectrum (a Z₂-twisted partner per
state). The dense light sector cannot absorb arbitrary extras, so the
Möbius/hybrid baryons must either (i) COINCIDE with observed-but-
unexplained resonances (filling "missing resonances" the quark model
under-predicts), or (ii) DECOUPLE from the dominant `πN` production/decay
(the standard missing-resonance mechanism). Over-prediction without
decoupling would be excluded. So this is BAM's most exposed topological
prediction — a genuine, near-term-testable constraint.

## Tests

| # | test | finding |
|---|---|---|
| T1 | setup | baryonic exotics from non-orientable topology (Möbius/hybrid) |
| T2 | no exotic J^P | every half-integer `J^P` ordinary for qqq ⟹ supernumerary |
| T3 | masses | hybrid N ≈ 1.79 GeV, hybrid Δ ≈ 2.08 GeV (base + `2√σ`) |
| T4 | constraint ranking | light N*/Δ* > strange hyperons > heavy baryons |
| T5 | most constrained | lightest exotic in the dense N*/Δ* region |
| T6 | test/tension | Möbius doubling: coincide or decouple, else excluded |
| T7 | vs glueballs | opposite extreme (most vs least constrained) |
| T8 | assessment | `BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED` |

## Established and open

  - **Established (BAM-native):** baryonic exotics have no exotic-`J^P`
    smoking gun (any `J^P` is ordinary for qqq), so BAM's Möbius/hybrid
    baryons are supernumerary ordinary-`J^P` states identifiable only by
    counting; they land in the light N*/Δ* region (~1.8–2.1 GeV, the
    `2√σ` gap), the MOST experimentally constrained channel; the
    constraint ranking is light N*/Δ* > strange hyperons > heavy-quark
    baryons (the freest). This is the opposite extreme from the unobserved
    glueballs (PR #100).

  - **Open:** a specific Möbius-baryon state confirmed in data — without a
    smoking-gun `J^P` the prediction is a counting one, testable only by
    matching the dense spectrum or demonstrating decoupling; and the
    heavy-quark baryon exotics (the least-constrained channel, where a
    clean new state is most likely findable).

## Cross-references

  - `geometrodynamics/qcd/topology.py` — `make_mobius_baryon_y_network`,
    `make_mobius_baryon_v12`.
  - `docs/mobius_exotic_sector_research_plan.md` — PR #101, the mesonic
    exotics (with the smoking-gun 1-+) this contrasts against.
  - `docs/glueball_closed_flux_loop_research_plan.md` — PR #100, the
    least-constrained pole (unobserved glueballs).

## Run

```
python -m experiments.closure_ledger.baryonic_exotics_classification_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_baryonic_exotics_classification_probe/`.
Expected verdict: `BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED`, 8/8 PASS.
