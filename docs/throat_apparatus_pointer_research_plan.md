# Probe M — the charge-conserving apparatus and the complete pointer statistics (PR #240)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The open item

#239 closed with an explicit instruction: **construct a charge-conserving
throat–apparatus interaction and evaluate the complete multi-outcome
pointer statistics before deciding whether the winding carrier must be
replaced by the spin doublet.**

Both halves mattered. #239's two reported costs had already been retracted
once under review, and the remaining one — a *"2.33 ceiling"* from leakage
— rested on a treatment of the leaked population that had never been
checked against what the apparatus actually does.

**It doesn't survive the check either. The winding carrier stands.**

## 1. The apparatus — and it conserves the *right* charge

Two pieces:

  - **Setting interaction.** The field's winding moves by `Δk = +2` while
    one quantum is absorbed from a **carrier** holding 2 units of winding
    each.
  - **Measurement interaction.** The committed winding Stern–Gerlach,
    `g·K_field ⊗ p̂`, which is winding-**diagonal** and therefore conserves
    charge identically.

| `n_max` | `‖[H, e^{2πiK_total/8}]‖` | `‖[H, K_integer]‖` |
|---:|---:|---:|
| 8 | **1.1e-15** | 9.24 |
| 16 | **1.1e-15** | 13.06 |
| 32 | **6.9e-15** | 18.48 |

> **A correction to #239 falls out.** Winding on a discrete `N_χ = 8`
> fiber is a **Z₈** charge, not an integer one, so the conserved operator
> is `exp(2πiK_total/8)` — **not** the integer-valued `K`. #239's T5
> obtained `0.0e+00` only because it *excluded the wrap-around
> transitions*; with them included the integer test fails outright. **The
> conclusion — total charge is conserved — survives; the operator that
> expresses it does not.**

## 2. Every channel is detected — there is no "no-click"

Under `Δk = ±2` the orbit of `k = +1` closes on the **4-cycle
`{+1, +3, −3, −1}`**. So the population #239 called *leakage* sits in
`k = ±3` and nowhere else.

A winding Stern–Gerlach deflects **proportionally to `k`**, so those
channels land at their own pointer positions rather than failing to
register:

| `g_p` | pointer centres | min separation | branch overlap |
|---:|---|---:|---:|
| 0.5 | 0.5, 1.5, −1.5, −0.5 | 1.0σ | 6.2e-01 |
| 1.0 | 1, 3, −3, −1 | 2.0σ | 3.2e-01 |
| **2.0** | 2, 6, −6, −2 | **4.0σ** | **4.6e-02** |

> **The apparatus has four outcomes per mouth and no no-click outcome at
> all.** #239's three-outcome model was *wrong about the detector*, not
> conservative about it.

## 3. The ceiling disappears

Binning the complete statistics with the natural setting-independent rule
`sign(k)` — positive winding → +1, negative → −1:

| `\|t\|` | `sign(k)` binning | #239's no-click binning | best of all 16 binnings |
|---:|---:|---:|---:|
| 0.4 | **2.134677** | 2.116652 | 2.134677 |
| 0.8 | **2.457397** | 2.305325 | 2.457397 |
| 1.2 | **2.758947** | 2.330532 | 2.758947 |
| 1.5 | **2.824508** | 2.327297 | 2.824508 |
| **1.9** | **2.828028** | 2.331306 | **2.828028** |

Tsirelson is **2.828427** — so `sign(k)` saturates it to grid accuracy,
and is itself the best of the 16 setting-independent binnings.

> #239's rule — send `k = ±3` to a fixed outcome as though undetected — is
> **one of those 16**, and on the same data it yields only 2.33.
> **The "2.33 leakage ceiling" was not a property of the winding carrier;
> it was the cost of throwing away channels the apparatus resolves.**

## 4. Finite-carrier back-action, computed exactly

The mean-field setting is *not* assumed: the carrier is kept and traced
out exactly, so its back-action on the field is included.

| `n̄` | max CHSH | coherent-state truncation |
|---:|---:|---:|
| 1 | 2.290317 | 0.0e+00 |
| 4 | 2.601231 | 1.0e-15 |
| 16 | 2.760053 | 2.5e-13 |
| ∞ (mean field) | 2.8245 | — |

Monotone convergence. A moderately populated carrier (`n̄ ≈ 16`) already
recovers most of the violation.

> *A bug worth recording:* an earlier run of this calculation showed a
> spurious **non-monotonic drop** at large `n̄` — the cause was a coherent
> state truncated *below its own mean*. The truncation figures above are
> the fix and the check.

## 5. No-signaling on the complete statistics

With four outcomes per mouth there is more room for a marginal to move, so
it is worth checking. Alice's complete outcome distribution is invariant
under every one of Bob's settings to **2.2e-16**.

## 6. The decision

| claim | status |
|---|---|
| #239 v1: the interaction is charge-non-conserving | **already retracted** (in #239, under review) — a prescribed external potential |
| #239 v2: leakage caps the violation at 2.33 | **retracted here** — the leaked channels are detected; `sign(k)` gives 2.828 |
| #239's charge-test operator | **corrected** — winding is Z₈; the integer-`K` test fails once wrap-around is included |
| **the winding carrier must be replaced by the spin doublet** | **no** |
| what actually remains | **one structural requirement** — a winding-2 carrier at the mouth, populated to `n̄ ≈ 16` |

The relocation question is settled in the negative, and the pattern is
worth stating plainly: across #239 and this probe, **three separate
"costs" of the winding carrier were reported and all three dissolved on
closer calculation** — a prescribed potential mistaken for broken charge
conservation, an efficiency proxy mistaken for a Bell test, and detected
channels mistaken for lost ones.

What survives is a single requirement *on the geometry* rather than an
argument for moving the chain. The spin doublet may still be worth
exploring on its own merits, but **not as a rescue**.

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | the open item | posed |
| T2 | apparatus + the Z₈ correction | `‖[H,Z₈]‖ = 1.1e-15`; integer test fails |
| T3 | complete outcome set | 4-cycle, all resolved at 4.0σ |
| T4 | **the ceiling disappears** | `sign(k)` → **2.828028** vs #239's 2.33 |
| T5 | finite-carrier back-action | monotone 2.29 / 2.60 / 2.76 |
| T6 | no-signaling | 2.2e-16 |
| T7 | the decision | relocation **not** required |
| T8 | assessment | 8/8 |

## 8. Open

  - **Whether any BAM geometry supplies a winding-2 carrier at a mouth,
    populated to `n̄ ≈ 16`** — now the only surviving requirement, and
    untested.
  - The spin doublet remains unexplored on its own merits, but is no
    longer needed as a rescue for the winding carrier.
  - The multipartite chain (#207/#208) under the same apparatus.
  - A spatially resolved pointer rather than the orthogonal-branch
    idealization used for the outcome statistics here.

## 9. Reproduce

```bash
python -m experiments.closure_ledger.throat_apparatus_pointer_probe
```

Expected verdict:
`THE_APPARATUS_IS_CHARGE_CONSERVING_UNDER_THE_Z8_OPERATOR_AND_THE_COMPLETE_POINTER_STATISTICS_REACH_TSIRELSON_BECAUSE_THE_LEAKED_CHANNELS_ARE_DETECTED_SO_THE_2_33_CEILING_WAS_AN_ARTIFACT_AND_THE_WINDING_CARRIER_NEED_NOT_BE_REPLACED`, 8/8 PASS.

## 10. Cross-references

  - `docs/minimal_mixing_interaction_research_plan.md` — #239, whose
    ceiling and charge-test operator this corrects
  - `docs/detector_algebra_research_plan.md` — #238, the gap being filled
  - `docs/measurement_sector.md` — #209, the winding Stern–Gerlach whose
    `k`-proportional deflection is the whole point of §2
