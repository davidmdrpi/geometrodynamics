# Probe H — one throat, two mouth states, no twist control (PR #235)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question

#233 (corrected) showed the exchange freeze and deck class 1 are
**independent**, and exhibited a two-throat ring carrying both. For the
one-throat ring it checked only the interior channel, found a
nondegenerate box ladder, and concluded "no doublet". But #224's doublet
is a **mouth** doublet, and a single throat already **has** two mouths. So
the sharper question:

> **Can a single throat with a finite interior channel produce two
> mouth-localized states while retaining deck class 1?**

**Yes to the states. No to the doublet that matters.** The gap between
those two answers is the result.

## 1. The states exist — #233's "no doublet" was too narrow

Put a barrier in the middle of the interior channel and each sub-well
abuts one **mouth**. With one exterior π-arc the ring keeps deck class 1.
The symmetric/antisymmetric pair rotates (the #234 basin construction) to
two states each living in a single well:

| bump `h` | deck class | splitting | single-well localization |
|---:|---:|---:|---:|
| 2 | 1 | 1.281e-01 | 0.9781 |
| 4 | 1 | 5.955e-02 | 0.9928 |
| 8 | 1 | 1.959e-02 | 0.9983 |
| 14 | 1 | 5.849e-03 | 0.9991 |

So #233 was right about the channel **ladder** and wrong to generalize
from it: a *structured* channel gives a mouth **pair**, on a
deck-class-1 ring, with the splitting tunable over nearly two decades.

## 2. But the twist has no purchase on them

The two wells are joined by **two paths** — through the central bump, or
the long way round through the exterior arc — and only the second carries
the holonomy. So `split = 2|t_bump + W·t_arc|`, and the two twist sectors
separate the paths:

| bump `h` | split untwisted | split twisted | ratio | `\|t_bump\|` | `\|t_arc\|` |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 2.7441e-01 | 2.7442e-01 | 1.0000 | 1.37e-01 | 2.26e-06 |
| 2 | 1.2811e-01 | 1.2812e-01 | 1.0001 | 6.41e-02 | 2.05e-06 |
| 4 | 5.9552e-02 | 5.9561e-02 | 1.0002 | 2.98e-02 | 2.27e-06 |
| 8 | 1.9585e-02 | 1.9596e-02 | 1.0005 | 9.80e-03 | 2.58e-06 |
| 14 | 5.8489e-03 | 5.8604e-03 | 1.0020 | 2.93e-03 | 2.87e-06 |

`|t_arc|` is **independent of the bump height** (2.0–2.9e-6 while
`|t_bump|` falls 47×) — the two-path decomposition confirmed, not assumed.
But the exterior path crosses **two mouth barriers** and is 3–5 orders of
magnitude weaker, so the interference never balances. The doublet is
real; it is **not twist-controlled**.

## 3. The resonant two-basin doublet is no better

The other way to build a pair: tune the channel length so a **channel**
level crosses an **arc** level. Those two basins *are* joined by two mouth
barriers of equal strength — the most favourable case for the twist. An
avoided crossing does appear near `D ≈ 8.25` (gap 3.2e-3, with the lower
state's channel fraction swapping across it), and the twist moves it by
only ~10%. Section 4 says why.

## 4. The mechanism: a reflection is not a translation

The freeze is the `S^k = W` mechanism — a cyclic **translation** whose
`k`-th power is the holonomy, with `k` even so that no eigenvalue is real.
**A one-throat ring has no such translation.** Its two mouth walls both
face the *same* interior channel, so they are **mirror images**, and the
ring carries a **reflection** `R` instead. And

```
R² = +1   in BOTH twist sectors, always
```

A reflection has no holonomy-dependent square, so it can never do what
`S² = W` does. Demonstrated by building the same ring both ways — same
lengths, same basin profiles, same deck class 1, same two mouths, the
**only** difference being wall orientation:

| walls | `‖[H,S]‖` | `‖[H,R]‖` | `S²=W` | `R²=+1` | twisted gap | frozen? |
|---|---:|---:|---|---|---:|---|
| **mirror** (physical) | 4.90e+01 | **0.00e+00** | yes | yes | 5.299e-01 | no |
| **translate** (counterfactual) | **0.00e+00** | 4.90e+01 | yes | yes | 2.954e-12 | **YES** |

The symmetry swaps exactly, and the freeze follows the *translation*, not
the ring, the lengths or the deck class.

**The counterfactual is not realizable by one throat.** Translate-oriented
walls would need the second wall to face *away* from the interior — which
is precisely what a **second throat** provides.

## 5. What this does to the arc

| claim | status |
|---|---|
| #233 (corrected): "one throat with `D > 0` gives a nondegenerate ladder, so no doublet" | **too narrow** — true of the ladder; a structured channel does give a mouth pair (§1). The conclusion survives for the better reason in §4 |
| the integer governing the freeze is the throat count `N` | **sharpened** — it is the **order of the cyclic translation**. That equals the throat count *because* each throat contributes one interior basin bounded by a mirrored pair of walls, so the shortest translation period is one throat |
| #231's construction (doublet **and** deck class 1) | **still exists** — #233's quarter-circle two-throat ring. Two throats are now the minimum for a structural reason, not a numerological one |
| #234's gate set could move onto a one-throat network | **no** — and not because the doublet is missing. It is there (§1); the **twist** has nothing to act on (§2), so the memory/gate switch has nothing to switch |

The through-line of #233 → #235: the freeze was first attributed to
**arc counting** (wrong), then to **channel structure** (too narrow), and
is now attributed to the property that actually controls it — whether a
cyclic translation commutes with `H`. That last one is checkable on *any*
network, which the previous two were not.

## 6. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | two mouth-localized states at deck class 1 | **yes**, localization > 0.97 |
| T3 | the twist does not act on them | ratio > 0.999; both paths extracted |
| T4 | the resonant two-basin doublet | avoided crossing, still not frozen |
| T5 | mirror vs translate walls | `‖[H,S]‖` 4.9e+01 vs 0; freeze follows S |
| T6 | consequences for #231, #233, #234 | 1 narrowed, 1 sharpened, 2 standing |
| T7 | assessment | 7/7 |

## 7. Open

  - Whether any BAM geometry realizes translate-oriented walls on a single
    throat. It needs the second wall to face away from the interior, which
    is what a second throat provides.
  - The selection question is still untouched — nothing selects the twist
    on any network examined so far.
  - Whether #233's quarter-circle two-throat ring supports #234's dressing
    calibration is uncomputed.

## 8. Reproduce

```bash
python -m experiments.closure_ledger.single_throat_mouth_doublet_probe
```

Expected verdict:
`ONE_THROAT_GIVES_TWO_MOUTH_LOCALIZED_STATES_AT_DECK_CLASS_ONE_BUT_NEVER_A_TWIST_CONTROLLED_DOUBLET_BECAUSE_ITS_TWO_MOUTH_WALLS_MIRROR_RATHER_THAN_TRANSLATE_SO_THE_RING_CARRIES_A_REFLECTION_WITH_R_SQUARED_PLUS_ONE`, 7/7 PASS.

A discretization note, continuing #233's: sampling is **cell-centred**, so
that the reflection and the translation are exact symmetries of the
*discrete* operator. On vertex-centred grids both pick up a one-point
residual that swamps the effect being measured.

## 9. Cross-references

  - `docs/freeze_deck_parity_obstruction_research_plan.md` — #233, whose
    one-throat conclusion this narrows and whose freeze law this sharpens
  - `docs/rp3_spin_structure_on_223_ring_research_plan.md` — #231, whose
    construction this leaves standing
  - `docs/physical_gate_set.md` — #234, whose twist switch this shows
    cannot be carried onto a one-throat network
  - `docs/mouth_exchange_dynamics.md` — #224, the mouth doublet being
    generalized
