# Generation-dependent ε_n and the neutrino hierarchy spread (PR #113)

PR #112 derived the **smallness** of the neutrino mass from a sub-throat,
bulk-geometric healing length `ε ~ R_c³` (the bounce `S = t²·P0·L*(ε)` with
the meV scale an output), but flagged that a **uniform** `ε` gives a uniform
action, hence `m_ν ∝ m_D` — a mere ×2.7 generation spread, far short of the
observed hierarchy. PR #91 suggested the fix: the generations are cavity
radial overtones `n`, and the overtone boundary stress `χ_n` (PR #79)
**decreases** with `n` (0.304, 0.097, 0.039), so higher-overtone necks are
more compliant — a generation-dependent `ε_n` that widens the spread "in the
right direction." This PR makes that quantitative and tests it honestly.

## The mechanism and its direction (right)

Compliance is the inverse of stiffness, so the natural law is

```
ε_n ∝ 1/χ_n     (more compliant neck ⟹ larger healing length).
```

With `χ_n` decreasing, `ε_n` increases with `n` ⟹ `L*(ε_n)` decreases ⟹
`S(n)` decreases ⟹ less suppression for higher overtones ⟹ `m_ν` increases
with `n`. With `m_D` also increasing with `n`, this gives **normal
ordering**, untuned — the **direction is correct**.

## The magnitude (overshoots badly)

But the observed hierarchy needs only a **gentle** `ε_n` variation.
Anchoring gen 1 at `ε_1 = R_c³` (the PR #112 value, `m_ν1 ≈ 2.08 meV`) and
demanding the observed `m_2 = 8.65`, `m_3 = 50.34 meV`:

| gen | χ_n | required ε_n | required ratio | χ-driven ε_n ratio | χ-driven m_ν (meV) |
|---|---:|---:|---:|---:|---:|
| 1 | 0.304 | 0.0110 | 1.00 | 1.00 | 2.1 |
| 2 | 0.097 | 0.0130 | 1.18 | 3.13 | 1038 |
| 3 | 0.039 | 0.0172 | 1.57 | 7.79 | 167650 |

The required `ε_n` rises only ×1.57 over three generations; the principled
`ε_n ∝ 1/χ_n` rises ×7.79, giving

```
m_ν = (2.1, 1038, 167650) meV  ⟹  m_ν3/m_ν2 ≈ 162  vs observed 5.85
                               ⟹  ×28 overshoot (orders of magnitude in absolute mass).
```

## Why: the steep bounce

The culprit is the bounce **steepness** established in PR #112,
`∂ln m_ν/∂ln ε ≈ 4.8`: the factor-~8 variation in `χ_n` is amplified into ~4
orders of magnitude in mass. The power in `ε_n ∝ χ_n^{−p}` that would
reproduce the data is not the principled `p = 1` but an inconsistent
fractional `p ≈ 0.15` (gen 1→2) to `0.31` (gen 2→3) — no single clean law
fits both ratios. So a generation-dependent `ε_n` can **accommodate** the
spread (by fitting a gentle profile) but cannot **predict** it from `χ_n`.

## The honest verdict

The hierarchy **direction** is derived (overtone compliance ⟹ normal
ordering, untuned — a sharpening of PR #91), but its **magnitude** is not:
the natural `χ_n` driver overshoots by ~28×, and the same bounce steepness
that made `ε`'s absolute value a residual (PR #112) now blocks the natural
overtone variation from setting the spread. The neutrino hierarchy spread
stays a **residual** — `ε_n` accommodates it but does not derive it — and it
plausibly belongs to the **mixing / anarchy sector** (PR #92) rather than a
generation-dependent healing length.

## Tests

| # | test | finding |
|---|---|---|
| T1 | setup | uniform ε ⟹ ×2.7 spread; PR #91 proposed χ_n-driven ε_n |
| T2 | direction | `ε_n ∝ 1/χ_n` ⟹ normal ordering (direction right, untuned) |
| T3 | required gentle | required ε_n ratios (1, 1.18, 1.57) to hit observed m_2, m_3 |
| T4 | χ overshoots | χ-driven m_3/m_2 = 162 vs 5.85 (×28) |
| T5 | steepness | bounce ×4.8 amplifies ×8 χ_n; power p 0.15→0.31 (≠1) |
| T6 | no clean law | accommodates (fit), does not predict |
| T7 | honest scope | direction derived; spread residual (mixing/anarchy, PR #92) |
| T8 | assessment | `HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** the hierarchy DIRECTION — a
    generation-dependent `ε_n ∝ 1/χ_n` from the overtone boundary stress
    gives normal ordering, untuned (sharpens PR #91's "right direction").

  - **Open:** the hierarchy MAGNITUDE — the natural `χ_n` driver overshoots
    by ×28; the required `ε_n` is gentle (×1.57) and not a principled
    function of `χ_n`. The same bounce steepness (`m_ν ∝ ε^{4.8}`, PR #112)
    amplifies the overtone variation past the data. The spread stays a
    residual, plausibly the mixing/anarchy sector (PR #92).

## Cross-references

  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, the bounce
    steepness (`m_ν ∝ ε^{4.8}`) that both made `ε`'s value residual and now
    blocks the spread.
  - `docs/boundary_compliance_bulk_geometry_research_plan.md` — PR #90, the
    `ε` healing length and the uniform-`S` ×2.7 spread.
  - `docs/generation_spread_pmns_mixing_research_plan.md` — PR #91, the
    overtone `χ_n` "right direction" claim this sharpens.
  - `docs/cross_channel_pmns_overlap_research_plan.md` — PR #92, the
    mixing/anarchy sector the spread plausibly belongs to.

## Run

```
python -m experiments.closure_ledger.generation_dependent_eps_n_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_generation_dependent_eps_n_probe/`.
Expected verdict:
`HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL`,
8/8 PASS.
