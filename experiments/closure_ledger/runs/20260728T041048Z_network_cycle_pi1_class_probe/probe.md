# Probe E — does the #224 cycle generate π₁(RP³)? (PR #230)

_Run 2026-07-28T04:10:48.177905+00:00 · 5/5 PASS_

## The lift test

| network | arcs | total exterior | lift | deck class | generator? |
|---|---:|---:|---|---:|---|
| #224 two-throat ring (build_two_throat) | 2 | 1.999π | closes | 0 | no |
| #223 single two-mouth bridge ring | 1 | 1.000π | antipodal | 1 | **YES** |

## Consequence

| probe | labelled by | bears on the #224 twist? |
|---|---|---|
| B (rp3_dirac_eta_probe) | the deck Z2 | **no** |
| C (bulk_aps_spin_structure_probe) | the deck Z2 (and, for ABK, the mouth Pin- structure) | **no** |
| D (twist_sector_einstein_dirac_probe) | the deck Z2 | **no** |

What would bear: a cycle with nontrivial deck class -- e.g. the #223 single-bridge ring, whose single pi-arc IS the generator (T3)

## Verdict

**THE_224_CYCLE_DOES_NOT_GENERATE_PI1_RP3_ITS_DECK_CLASS_IS_TRIVIAL_AND_ITS_HOLONOMY_LIVES_ON_A_HANDLE_GENERATOR_SO_PROBES_B_TO_D_DO_NOT_BEAR_ON_IT**
