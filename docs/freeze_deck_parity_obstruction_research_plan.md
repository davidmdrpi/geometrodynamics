# Probe G — the freeze/deck parity obstruction (PR #233)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question

Probe F (#231) **observed** that the two properties needed to make the
exchange-freeze result a statement about the spin structure never occur
together:

| network | throats | deck class | labels coincide? | doublet? |
|---|---:|---:|---|---|
| #223 | 1 | 1 | **yes** | no |
| #224 | 2 | 0 | no | **yes** |

and named "an interior channel (`D > 0`) on the #223 ring, or some other
geometry carrying both" as the obvious next construction. But an
observation across the two networks the repo happens to contain is not an
obstruction. **Is it structural, or an accident of which two rings were
built?**

## 1. The answer: structural — a parity obstruction

Both properties are governed by the **same integer `N`** with **opposite
parities**.

### 1.1 Deck class = `N mod 2`

A cyclic ring of `N` throats carries exactly `N` exterior π-arcs (one
between each consecutive pair of mouths), so its total exterior length is
`Nπ` and the `S³` lift closes iff `N` is even:

```
deck class = N mod 2      →   class 1  ⟺  N ODD
```

This reproduces Probe E's two data points as the `N = 1` and `N = 2` cases
of one law.

### 1.2 Full freeze ⟺ `N` even

The cyclic translation `S` advancing the ring by one throat satisfies
`S^N = W`. The freeze is the statement that **every** level is forced into
a degenerate multiplet, which happens iff `S` has **no real eigenvalue** —
a real operator pairs complex eigenvalues, and only those. For `W = −1`
the eigenvalues are the `N`-th roots of `−1`, `e^{iπ(2k+1)/N}`, and one is
real iff `(2k+1)/N` is an integer for some `k` — possible **iff `N` is
odd**. Hence:

```
full freeze  ⟺  N EVEN
```

### 1.3 Therefore mutually exclusive

```
freeze       ⟺  N even
deck class 1 ⟺  N odd
```

**No `N` satisfies both**, for any ring size.

## 2. Verified, not just argued

The root-of-unity argument is checked against **actual** `N`-throat ring
operators, built by extending #224's own construction (`N` channels, `2N`
barriers, `N` arcs), measuring the max intra-pair gap over the 16 lowest
levels:

| `N` | gap untwisted | gap twisted | frozen? |
|---:|---:|---:|---|
| 2 | 2.082e-03 | **1.995e-11** | **YES** |
| 3 | 3.958e-01 | **3.953e-01** | **no** |
| 4 | 4.592e-05 | **5.664e-12** | **YES** |

`N = 3` is the decisive case: the smallest ring with **both** several
throats **and** deck class 1 — exactly the combination #231 was looking
for — and it simply **does not freeze**.

## 3. What this does to PR #231

| claim | status |
|---|---|
| Probe F: the phenomenon and the coupling never co-occur | **UPGRADED** — observed across two networks, now proven for all cyclic throat rings |
| #231's recommended next construction (a geometry carrying both) | **WITHDRAWN** — no such geometry exists |
| the freeze sector is cut off from RP³ spin-structure data | **STRENGTHENED** — permanently, not contingently |

The withdrawal matters most. #231 ended by naming a construction as the
obvious next step; this rules it out. Adding an interior channel to the
`N = 1` ring does not add a throat, so it cannot create a mouth doublet;
and any ring with enough throats to have a doublet **and** deck class 1
has odd `N`, which §2 measures as not freezing. **The construction cannot
succeed, and the reason is parity, not effort.**

Closing a line of attack is worth more than another inconclusive probe
along it.

## 4. Consequence for the selection question

Every freezing network has even `N`, hence deck class 0, hence a holonomy
that is a **free-`Z` handle character** rather than a torsion class. So
`η` invariants, APS inflow and bordism are the wrong tools **by
construction**, on every such network — not merely on #224.

What remains open is unchanged in content but sharper in direction: a
selection argument for the freeze sector must address a free-`Z` handle
character, and no torsion-flavoured invariant will reach it.

## 5. Tests

| # | test | outcome |
|---|---|---|
| T1 | is the non-co-occurrence structural? | posed |
| T2 | deck class = `N mod 2` | law holds for `N = 1…5` |
| T3 | full freeze ⟺ `N` even | argued (roots of unity) + measured (`N = 2,3,4`) |
| T4 | mutually exclusive by parity | no `N` has both |
| T5 | consequences for #231 | 1 upgraded, 1 withdrawn, 1 strengthened |
| T6 | assessment | 6/6 |

## 6. Open

  - A selection argument for the freeze sector must address a **free-`Z`
    handle character**; no torsion invariant reaches it.
  - **The parity law is derived for *cyclic* rings of identical throats.**
    Non-cyclic or inhomogeneous networks — trees with several independent
    cycles, tri-mouth junctions — are **not covered** and could in
    principle evade it. That, not the #231 construction, is where to look
    next.
  - Whether the freeze itself survives on a non-cyclic network is
    **untested**: the `S^N = W` argument needs the cyclic symmetry, so
    both halves of the obstruction would have to be re-derived there.

## 7. Reproduce

```bash
python -m experiments.closure_ledger.freeze_deck_parity_obstruction_probe
```

Expected verdict:
`THE_FREEZE_REQUIRES_EVEN_N_AND_DECK_CLASS_ONE_REQUIRES_ODD_N_SO_THEY_ARE_MUTUALLY_EXCLUSIVE_BY_PARITY_AND_THE_231_RECOMMENDED_CONSTRUCTION_CANNOT_EXIST`, 6/6 PASS.

## 8. Cross-references

  - `docs/rp3_spin_structure_on_223_ring_research_plan.md` — #231, Probe F,
    whose observation this upgrades and whose next step this withdraws
  - `docs/rp3_spin_structure_selection_research_plan.md` — Probes A–E
  - `docs/mouth_exchange_dynamics.md` — #224, whose `build_two_throat` is
    generalised here to `N` throats
