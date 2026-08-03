# Probe G — the freeze/deck parity obstruction (PR #233)

_Run 2026-08-03T06:08:59.973803+00:00 · 6/6 PASS_

## The two parity laws

| `N` throats | arcs | deck class | labels coincide? | doublet? | full freeze? |
|---:|---:|---:|---|---|---|
| 1 | 1 | 1 | yes | no | no |
| 2 | 2 | 0 | no | yes | **yes** |
| 3 | 3 | 1 | yes | yes | no |
| 4 | 4 | 0 | no | yes | **yes** |
| 5 | 5 | 1 | yes | yes | no |
| 6 | 6 | 0 | no | yes | **yes** |

**No `N` has both** — freeze ⟺ `N` even, deck class 1 ⟺ `N` odd.

## Measured on real N-throat rings

| `N` | gap untwisted | gap twisted | frozen? |
|---:|---:|---:|---|
| 2 | 2.082e-03 | 1.995e-11 | **YES** |
| 3 | 3.958e-01 | 3.953e-01 | no |
| 4 | 4.592e-05 | 5.664e-12 | **YES** |

## Consequences

  - *Probe F: the phenomenon and the coupling never co-occur* — **UPGRADED**
  - *#231's recommended next construction: an interior channel (D > 0) on the #223 ring, or some other geometry, carrying BOTH the doublet and deck class 1* — **WITHDRAWN**
  - *the exchange-freeze sector is cut off from RP^3 spin-structure data* — **STRENGTHENED**

## Verdict

**THE_FREEZE_REQUIRES_EVEN_N_AND_DECK_CLASS_ONE_REQUIRES_ODD_N_SO_THEY_ARE_MUTUALLY_EXCLUSIVE_BY_PARITY_AND_THE_231_RECOMMENDED_CONSTRUCTION_CANNOT_EXIST**
