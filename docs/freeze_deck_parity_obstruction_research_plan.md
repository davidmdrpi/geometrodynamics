# Probe G — the freeze and deck class 1 are **independent** (PR #233)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

> **Correction notice.** The first version of this document claimed a
> **parity obstruction** making the exchange freeze and deck class 1
> mutually exclusive, and on that basis withdrew #231's recommended next
> construction as impossible. **That claim was wrong.** It rested on an
> unstated assumption — that every exterior arc has length exactly π —
> and the geometry it ruled out exists. Everything below is the corrected
> account, with the retraction spelled out in §5 rather than quietly
> rewritten.

## 0. The question

Probe F (#231) **observed** that the two properties needed to make the
exchange-freeze result a statement about the spin structure never occur
together:

| network | throats | deck class | labels coincide? | doublet? |
|---|---:|---:|---|---|
| #223 | 1 | 1 | **yes** | no |
| #224 | 2 | 0 | no | **yes** |

and named "an interior channel (`D > 0`) on the #223 ring, or some other
geometry carrying both" as the obvious next construction. An observation
across the two networks the repo happens to contain is not an obstruction.
**Is it structural, or an accident of which two rings were built?**

**An accident.**

## 1. What stands: the freeze law

The cyclic translation `S` advancing the ring by one throat satisfies
`S^N = W`. The freeze is the statement that **every** level is forced into
a degenerate multiplet, which happens iff `S` has **no real eigenvalue** —
a real operator pairs complex eigenvalues, and only those. For `W = −1`
the eigenvalues are the `N`-th roots of `−1`, `e^{iπ(2k+1)/N}`, and one is
real iff `(2k+1)/N` is an integer for some `k` — possible **iff `N` is
odd**. Hence

```
full freeze  ⟺  N EVEN
```

Checked three ways this time, not one:

| check | result |
|---|---|
| root-of-unity argument | freeze predicted ⟺ `N` even, for `N = 2…6` |
| measured on real rings (max intra-pair gap, 16 lowest levels) | `N=2` **2.2e-12** frozen · `N=3` **4.0e-01** not frozen · `N=4` **1.7e-12** frozen |
| **operator level** | `‖[H,S]‖ < 1e-8` and `S^N = W` **exactly**, every ring |

A numerical trap worth recording: the throat cell must be an **integer**
number of grid points, or `S` is not a symmetry of the *discrete* operator
and the freeze spuriously appears to fail. On an incommensurate grid the
counterexample ring below reports a twisted gap of `6e-4` instead of
`3e-12`.

## 2. What was wrong: the deck-class law

Probe E's criterion is that the `S³` lift closes iff the **total exterior
arc length** is an even multiple of π:

```
deck class = (Σ arc lengths) / π   mod 2
```

**That formula does not mention `N`.** The earlier claim — "an `N`-throat
ring carries `N` exterior π-arcs, so deck class = `N mod 2`" — silently
assumed every arc has length π. But an arc of length π joins a point to
its **antipode**, so assuming it for *every* arc of the ring forces all
`N` throats to share the **same antipodal pair of mouth locations**. That
is #224's configuration, not the general one. Where the mouths sit is free
geometric data.

| ring | arcs/π | Σ/π | deck class | `N mod 2` | old law |
|---|---:|---:|---:|---:|---|
| 1 throat, arc π (#223) | 1.0 | 1.0 | 1 | 1 | ✔ |
| 2 throats, arcs π (#224) | 1.0 | 2.0 | 0 | 0 | ✔ |
| 3 throats, arcs π | 1.0 | 3.0 | 1 | 1 | ✔ |
| **2 throats, arcs π/2** | **0.5** | **1.0** | **1** | 0 | ✘ |
| **2 throats, arcs 3π/2** | **1.5** | **3.0** | **1** | 0 | ✘ |
| 2 throats, arcs 2π | 2.0 | 4.0 | 0 | 0 | ✘ |
| **4 throats, arcs 3π/4** | **0.75** | **3.0** | **1** | 0 | ✘ |

The old law is exactly the all-arcs-π restriction — correct on the three
rings that satisfy it, wrong everywhere else.

## 3. The counterexample

Take `N` **even** so the freeze law applies, then place the mouths so the
arcs sum to an **odd** multiple of π. The simplest such ring: **two
throats whose mouth pairs are rotated a quarter circle apart**, so each
exterior arc is π/2 and the total is π.

| ring | Σ/π | deck | doublet | gap untwisted | gap twisted | frozen? |
|---|---:|---:|---|---:|---:|---|
| `N=2`, arcs π/2 | 1.0 | **1** | **yes** | 1.908e-02 | **8.482e-13** | **YES** |
| `N=2`, arcs 3π/2 | 3.0 | **1** | **yes** | 8.266e-03 | **3.297e-12** | **YES** |
| `N=4`, arcs 3π/4 | 3.0 | **1** | **yes** | 4.705e-04 | **2.093e-12** | **YES** |

And the freeze there is the **same mechanism** as #224's, not a numerical
coincidence: on each of these rings `‖[H,S]‖ < 1e-8` and `S^N = W` holds
exactly.

**So the geometry #231 asked for exists.** #224's own configuration is the
degenerate case in which the two throats happen to share one antipodal
mouth pair; nothing required that.

## 4. The other evasion classes — computed, not asserted

| class | what was done | result |
|---|---|---|
| **tuned / non-π arcs** | §3 | the two properties **co-occur** |
| **one throat, finite channel** (`D > 0`) | full ring operator at `D = 0.5, 4, 8, 16` | interior modes are a **nondegenerate ladder** (`D=8`: 0.3785, 0.7568, 1.1350, 1.5131, 1.8910); twisted and untwisted spectra agree to 5 digits; no freeze (gaps 1.9e-1 … 8.5e-1) |
| **inhomogeneous / internally structured** | unequal arcs (π, 2π); unequal channels (4, 5); internal bump in one throat | the freeze is **destroyed** even at even `N` (gaps 2.5e-1, 3.8e-1, 4.3e-1) |
| **non-cyclic** | theta network: 3 throats between two tri-mouth junctions, `b₁ = 2`, four twist sectors | **no full freeze in any sector**; the degeneracies present are `S₃` edge-permutation degeneracies, which appear in the **untwisted** sector too |

Two of these deserve comment.

  - The **one-throat `D > 0`** case is #231's own recommendation, and the
    earlier version of this probe dismissed it **by assertion** ("adding an
    interior channel does not add a throat, so it cannot create a mouth
    doublet"). Computed, the assertion holds — the interior channel gives a
    box ladder, not a pair, and the localized states never reach the
    holonomy. But asserting it was the error even though the answer was
    right, and it is moot regardless: §3 supplies the co-occurring geometry
    a different way.
  - The **fragility** result cuts against the freeze, and is the honest
    caveat that replaces the earlier version's (wrong) one. The freeze
    needs the cyclic symmetry *exactly*; any detuning of arc lengths,
    channel lengths or internal structure destroys it. A physical reading
    of the freeze has to explain why that exact symmetry would hold.

The theta network is **one representative** non-cyclic network, not a
classification. `S^N = W` has no analogue on a graph with no translation,
so both halves of any obstruction argument would have to be re-derived
there.

## 5. What is retracted, restored and still standing

| claim | status |
|---|---|
| full freeze ⟺ `N` even (`S^N = W`) | **stands** — now verified at operator level too |
| deck class = `N mod 2` | **corrected** — true only when every arc is π |
| freeze and deck class 1 are mutually exclusive by parity | **retracted** — they are *independent*; §3 exhibits rings carrying both |
| #231's recommended next construction | **restored, in corrected form** — it exists, by *moving the mouths*, not by adding an interior channel to the one-throat ring |
| the freeze sector is cut off from RP³ spin-structure data "permanently, not contingently" | **retracted** — on the quarter-circle ring the cycle **is** the π₁(RP³) generator |
| Probe F (#231): the phenomenon and the coupling never co-occur | **demoted** to an observation about the two rings the repo contains — correct as reported there, not a theorem |

The net effect is to **reopen** rather than close. Probes B–D (η = ±1/4,
equal `|det D|`, `h = 0`, `Δ⟨T_AB⟩ = 0`) can now be applied to a network
that **actually freezes**, which was the whole point of looking for one.
They still find nothing selecting the twist — but on the quarter-circle
ring that null result is at last a statement *about the freeze sector*,
which it never was before.

## 6. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal, and the claim under test | posed |
| T2 | freeze law: `S^N = W`, freeze ⟺ `N` even | argued + measured + operator-level |
| T3 | deck class is set by mouth placement, not `N` | old law fails on 4 of 7 rings |
| T4 | **the counterexample** | 3 rings carry freeze **and** deck class 1 |
| T5 | the remaining evasion classes | one-throat, inhomogeneous, non-cyclic — all computed |
| T6 | retracted / restored / standing | 1 stands, 1 corrected, 2 retracted, 1 restored, 1 demoted |
| T7 | assessment | 7/7 |

## 7. Open

  - **The selection question is untouched.** On the quarter-circle ring
    Probes B–D now apply to a network that freezes, and still nothing
    selects the twist.
  - **Why the exact cyclic symmetry would hold.** §4 shows small
    departures destroy the freeze. This is now the sharpest constraint on
    reading the freeze physically.
  - **Non-cyclic networks in general.** One theta network is a data point,
    not a classification.
  - **Is the quarter-circle ring preferred or merely permitted?** Nothing
    here says the mouths *would* sit a quarter circle apart rather than
    coincident as in #224.

## 8. Reproduce

```bash
python -m experiments.closure_ledger.freeze_deck_parity_obstruction_probe
```

Expected verdict:
`THE_FREEZE_NEEDS_EVEN_N_BUT_THE_DECK_CLASS_IS_SET_BY_MOUTH_PLACEMENT_NOT_BY_N_SO_THEY_ARE_INDEPENDENT_AND_A_TWO_THROAT_RING_WITH_QUARTER_CIRCLE_ARCS_CARRIES_BOTH`, 7/7 PASS.

## 9. Cross-references

  - `docs/rp3_spin_structure_on_223_ring_research_plan.md` — #231, Probe F,
    whose recommended construction this restores
  - `docs/rp3_spin_structure_selection_research_plan.md` — Probes A–E;
    Probe E supplies the deck-class criterion used in §2
  - `docs/mouth_exchange_dynamics.md` — #224, whose `build_two_throat` is
    generalised here to `N` throats **and** to arbitrary arc lengths
