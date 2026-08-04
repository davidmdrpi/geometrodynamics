# Probe H — one throat, two mouth states, no twist control (PR #235)

_Run 2026-08-04T03:02:05.645591+00:00 · 7/7 PASS_

**Q. Can a single throat with a finite interior channel produce two mouth-localized states while retaining deck class 1?**

**A. Yes — and the twist still cannot touch them.**

## The states exist

| bump `h` | deck | splitting | single-well localization |
|---:|---:|---:|---:|
| 2.0 | 1 | 1.281e-01 | 0.9781 |
| 4.0 | 1 | 5.955e-02 | 0.9928 |
| 8.0 | 1 | 1.959e-02 | 0.9983 |
| 14.0 | 1 | 5.849e-03 | 0.9991 |

## But the twist does not act on them

| bump `h` | split untwisted | split twisted | ratio | `|t_bump|` | `|t_arc|` |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 2.7441e-01 | 2.7442e-01 | 1.0000 | 1.37e-01 | 2.26e-06 |
| 1.0 | 2.0752e-01 | 2.0752e-01 | 1.0000 | 1.04e-01 | 2.03e-06 |
| 2.0 | 1.2811e-01 | 1.2812e-01 | 1.0001 | 6.41e-02 | 2.05e-06 |
| 3.0 | 8.5069e-02 | 8.5078e-02 | 1.0001 | 4.25e-02 | 2.16e-06 |
| 4.0 | 5.9552e-02 | 5.9561e-02 | 1.0002 | 2.98e-02 | 2.27e-06 |
| 6.0 | 3.2543e-02 | 3.2552e-02 | 1.0003 | 1.63e-02 | 2.45e-06 |
| 8.0 | 1.9585e-02 | 1.9596e-02 | 1.0005 | 9.80e-03 | 2.58e-06 |
| 10.0 | 1.2560e-02 | 1.2571e-02 | 1.0009 | 6.28e-03 | 2.70e-06 |
| 14.0 | 5.8489e-03 | 5.8604e-03 | 1.0020 | 2.93e-03 | 2.87e-06 |

`|t_arc|` is independent of the bump height — the two-path picture confirmed — but the exterior path crosses **two** mouth barriers, so the interference never balances.

## The mechanism

| walls | `‖[H,S]‖` | `S²=W`? | `R²=+1`? | twisted gap | frozen? |
|---|---:|---|---|---:|---|
| mirror | 4.9e+01 | yes | yes | 5.299e-01 | no |
| translate | 0.0e+00 | yes | yes | 2.954e-12 | **YES** |

Same lengths, same deck class 1, same two mouths — only the wall orientation differs, and it decides the freeze completely. A single throat's two walls both face the same interior, so they **mirror**; and `R² = +1` in both twist sectors, always.

## Consequences

  - *#233 (corrected): 'the one-throat ring with D > 0 gives a nondegenerate interior ladder, so no doublet'* — **TOO NARROW**
  - *the integer governing the freeze is the throat count N* — **SHARPENED**
  - *#231's construction (a geometry carrying the doublet and deck class 1)* — **STILL EXISTS**
  - *#234's gate set could be moved onto a one-throat network* — **NO**

## Verdict

**ONE_THROAT_GIVES_TWO_MOUTH_LOCALIZED_STATES_AT_DECK_CLASS_ONE_BUT_NEVER_A_TWIST_CONTROLLED_DOUBLET_BECAUSE_ITS_TWO_MOUTH_WALLS_MIRROR_RATHER_THAN_TRANSLATE_SO_THE_RING_CARRIES_A_REFLECTION_WITH_R_SQUARED_PLUS_ONE**
