# Probe F — the spin-structure machinery on the #223 ring (PR #231)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. Why this ring

Probe E (#230) showed the two Z₂ labels — the B2 spin structure in
`H¹(RP³;Z₂)` and the network holonomy `W` in `H¹(Γ;Z₂)` — coincide only on
a cycle of nontrivial deck class, and located exactly one such network:

| network | exterior arcs | deck class | labels coincide? |
|---|---:|---:|---|
| **#224** two-throat ring | 2 × π = 2π | **0** | no |
| **#223** single-bridge ring | 1 × π | **1** | **yes** |

So #223 is the one place the composition rule `θ = 0 iff W·η = +1` is
justified — and therefore the one place Probes B–D bear on the **network
twist** rather than merely on the spin structure. This pushes the
machinery through it.

## 1. What `W` is on this ring

`bridge_dressing_network_probe.build_ring` solves a **half** ring by
shooting from the neck out to the source point at `χ = π`, with a parity
choice at each end:

| parity | about | meaning |
|---|---|---|
| `odd_neck` | the neck | #202's Pin parity; selects the electron `k = 1` mode |
| `end` ∈ {`deriv`, `val`} | the **source point** (`χ = π`) | the two ways to close the deck-class-1 loop — **this is `W`** |

The source point *is* the antipode, so closing the ring there with Neumann
(`deriv`) versus Dirichlet (`val`) is precisely the two completions of the
loop whose deck class Probe E measured as 1.

Half-ring length `L = 5.1400`, closure comb `π/L = 0.6112`.

## 2. The measurement: the moding shift is **real** here

The composition rule predicts that flipping `W` shifts the mode comb by
**half a spacing**. Measured on the real construction at the electron
(odd-neck) parity:

| `r_s` | `W` = `val` (Dirichlet) | `W` = `deriv` (Neumann) | interleave fractions |
|---:|---|---|---|
| 0.3 | 1.3892, 2.1127, 2.7359, 3.3655 | 0.9374, 1.7722, 2.4268, 3.0498 | 0.529, 0.504, 0.499, 0.506, 0.498 |
| 0.4 | 1.3846, 2.1004, 2.7152, 3.3369 | 0.9347, 1.7644, 2.4099, 3.0257 | 0.531, 0.503, 0.499, 0.506, 0.498 |

The two `W` sectors **interleave**, with each `deriv` mode sitting at
fractional position → **0.5** between consecutive `val` modes, converging
as the modes climb out of the neck region where the bridge potential
distorts the low comb.

**The moding shift that the earlier arc could only assert *conditionally*
is realized here** — on the network where Probe E showed the composition
rule is justified. This is the first time in the whole investigation that
`W` and the spin structure have been shown to be the same label on an
actual construction.

## 3. And yet nothing selects — now legitimately

With the labels genuinely identified, the Probe B–D results apply to the
**network twist** and not merely to the spin structure:

| `W` | `η(0)` | `h` | `\|det D\|` | `arg det D` |
|---:|---:|---:|---:|---:|
| `+1` | +1/4 | 0 | 0.8963003229 | −π/8 |
| `−1` | −1/4 | 0 | 0.8963003229 | +π/8 |

| mechanism | selects? | why |
|---|---|---|
| energetics / vacuum energy | **no** | identical `\|spec\|` mode by mode ⟹ `ΔE_vac = 0` exactly |
| anomaly inflow / APS | **no** | `Ω₃^Spin = 0`, both extend |
| semiclassical back-reaction | **no** | `Δ⟨T_AB⟩ = 0` pointwise |

**The selection question is therefore answered for the one network where
it is well posed: nothing selects the twist.** This is no longer a scoping
caveat — the three mechanisms are now applied to the *right* label, and
they still fail.

## 4. The structural finding: the phenomenon and the coupling never co-occur

| network | deck class | labels coincide? | mouth doublet? | exchange freeze? |
|---|---:|---|---|---|
| **#223** | 1 | **yes** | **no** | **no** |
| **#224** | 0 | no | yes | **yes** |

  - **#223 couples the labels but has no doublet.** The repo says so in its
    own words: the exterior is one connected cavity (the two brane arcs
    join at the antipode), and the pure ultrastatic bridge "has no interior
    state to exchange" (#224's honest-scope section). There is no exchange
    freeze here to explain.
  - **#224 has the doublet and the freeze but decouples the labels.** Its
    deck class is 0, so its `W` is a handle character the spin-structure
    results do not reach (Probe E).

**This is why the selection question kept slipping: every mechanism was
being computed on one ring and applied to the other.** Naming the
obstruction is the useful part — a selection argument for the
exchange-freeze sector cannot come from RP³ spin-structure data *at all*,
because the network exhibiting the phenomenon does not carry that label.

## 5. Tests

| # | test | outcome |
|---|---|---|
| T1 | the one ring where the labels coincide | stated |
| T2 | `W` here is the source-point parity; deck class 1 | `L = 5.1400`, comb `0.6112` |
| T3 | **the moding shift, measured** | interleave fractions → 0.5 |
| T4 | B–D applied legitimately | identical `\|det\|`, `h = 0`, nothing selects |
| T5 | phenomenon and coupling never co-occur | verified across both rings |
| T6 | assessment | 6/6 |

## 6. Open

  - **Why matter sits in the twisted sector** — still open, and now known
    to be *unanswerable from RP³ spin-structure data* for the #224
    exchange-freeze sector.
  - **What a selection argument for a handle character would even look
    like.** The #224 `W` is a free-`Z` label, not a torsion class, so
    bordism invariants and `η` are the wrong tools. This is a genuine
    change in what to look for.
  - ~~**The one geometry that could join them.**~~ **WITHDRAWN (PR #233,
    Probe G): no such geometry exists.** A cyclic ring of `N` throats has
    `N` exterior π-arcs, so deck class `= N mod 2`; and the freeze requires
    the cyclic translation `S^N = W` to have no real eigenvalue, which for
    `W = −1` happens iff `N` is **even**. So freeze ⟺ `N` even and deck
    class 1 ⟺ `N` odd — **mutually exclusive by parity**, verified on real
    `N = 2, 3, 4` rings (`N = 3` twisted gap 3.953e-01: does not freeze).
    Adding an interior channel to the `N = 1` ring does not add a throat,
    so it cannot create a doublet either. The construction recommended
    above cannot succeed, and the reason is parity, not effort. §5's
    "never co-occur" is thereby **upgraded from observation to theorem**.
    See `docs/freeze_deck_parity_obstruction_research_plan.md`.

## 7. Reproduce

```bash
python -m experiments.closure_ledger.rp3_spin_structure_on_223_ring_probe
```

Expected verdict:
`ON_THE_223_RING_THE_LABELS_COINCIDE_AND_THE_MODING_SHIFT_IS_REAL_BUT_NOTHING_SELECTS_THE_TWIST_AND_THE_EXCHANGE_FREEZE_LIVES_ON_THE_OTHER_RING`, 6/6 PASS.

## 8. Cross-references

  - `docs/rp3_spin_structure_selection_research_plan.md` — Probes A–E; E
    located this ring
  - `docs/bridge_dressing_network.md` — #223, the ring solved here
  - `docs/mouth_exchange_dynamics.md` — #224, whose own honest-scope
    section supplies the "no interior state to exchange" finding
  - `docs/twist_sector_selection_research_plan.md` — the topological
    freeze, still untouched by any of this
