# m_e and √σ: independent anchors, or one bulk-gravity scale G? (PR #106)

PR #105 found `G` the dimensionful anchor and the root the #104 sector
anchors (`m_e`, `√σ`) descend from, and left the sharp question: are
`m_e` and `√σ` two INDEPENDENT dimensionful inputs, or readouts of a
single bulk-gravity scale `G` times dimensionless ratios? This probe
answers it.

## Common gravitational origin

Both #104 anchors are brane/geometric scales of the SAME bulk geometry,
descending from the bulk gravity (PR #57): the throat size `R_MID`
(⟹ `m_e = ℏc/R_MID`) is set by the Randall–Sundrum-tuned brane tension
`λ_crit = √(6|Λ₅|)/κ₅²`, and the confinement tension `σ ∝ √|Λ₅|/κ₅²`. So
`m_e` and `√σ` are NOT two independent KINDS of input — both are the one
bulk-gravity scale `G` read out in two channels (the throat-winding
lepton channel and the cavity-confinement QCD channel, PR #83).

## But the ratio is not derived

The dimensionless ratio is `√σ/m_e ≈ 0.424 GeV / 0.511 MeV ≈ 830`, the
lepton-throat (electron Compton) to QCD-confinement length hierarchy. If
the geometry fixed it, it would be a clean closure number. It is not:

| candidate | value | off |
|---|---:|---:|
| `50π·k_5` | 785 | −5.4% |
| `k_5⁴` | 625 | −24.7% |
| `k_5⁴ + k_5²` | 650 | −21.7% |
| `(2π)³` | 248 | −70.1% |
| `2·k_5⁴` | 1250 | +50.6% |

The nearest is `50π·k_5 = 785` (5.4% off) — a near-coincidence in the
spirit of `F_13 = 233` (PR #76), not a derivation. So the ~830 ratio is
**underived**: the geometry does not yet fix the relative normalisation
of the throat-winding and cavity-confinement channels.

## The verdict and the honest bookkeeping

`m_e` and `√σ` are NOT independent dimensionful anchors — they share the
one bulk-gravity scale `G`. Equivalently, the program has **one
foundational dimensionful anchor (`G`) plus one open dimensionless ratio**
(`m_e/√σ ≈ 1/830`). This reduces the dimensionful-anchor count from two to
one.

But it is a **repackaging, not a free-lunch reduction**: a dimensionful
anchor has been converted into a dimensionless residual, so the TOTAL
count of irreducible inputs is unchanged (was: 2 dimensionful anchors;
now: 1 dimensionful anchor + 1 dimensionless ratio). What it buys is
conceptual cleanliness — the sole fundamental SCALE is the gravitational
foundation `G`, with everything else dimensionless, exactly the
GR-foundational posture. The new ratio joins `ε`, `n_part`, and `α` as the
program's open dimensionless residuals; its smallness (the electron being
anomalously light relative to `Λ_QCD`) overlaps the flavor puzzle, though
it is not identical to it.

## Tests

| # | test | finding |
|---|---|---|
| T1 | the question | independent anchors, or one `G` + ratios? |
| T2 | common origin | both brane scales of one bulk (PR #57) ⟹ not independent kinds |
| T3 | ratio | `√σ/m_e ≈ 830` (lepton-throat / QCD-confinement hierarchy) |
| T4 | underived | no clean closure match (nearest `50π·k_5=785`, 5.4% off) |
| T5 | verdict | NOT independent — one `G` + one open ratio; dimensionful 2→1 |
| T6 | bookkeeping | a repackaging; total inputs unchanged; cleaner "one scale `G`" |
| T7 | ledger | one anchor `G`; `m_e/√σ` a new dimensionless residual |
| T8 | assessment | `M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO` |

## Established and open

  - **Established (BAM-native):** `m_e` and `√σ` both descend from the
    single bulk-gravity scale `G` (PR #57), so they are not independent;
    the dimensionful-anchor count reduces 2→1. But the dimensionless ratio
    `m_e/√σ ≈ 1/830` is underived (no clean closure match), so it becomes
    a new open dimensionless residual — the repackaging leaves the total
    irreducible-input count unchanged.

  - **Open:** a derivation of the ~830 lepton/QCD scale hierarchy (from
    the throat-winding vs cavity-confinement channel normalisation).
    Fixing it would genuinely reduce BAM to a single irreducible input
    (the gravitational scale `G`).

## Cross-references

  - `docs/alpha_G_ledger_classification_research_plan.md` — PR #105, where
    `G` was placed as the root anchor and the scale-count question raised.
  - `docs/brane_tension_tuning_research_plan.md` — PR #57, the bulk-gravity
    relation `λ_crit = √(6|Λ₅|)/κ₅²` both scales descend from.
  - `docs/throat_shell_mass_operator_unification_research_plan.md` — PR #83,
    the throat-winding vs cavity-resolving channels.
  - `docs/quark_npart_origin_research_plan.md` — PR #76, the `F_13 = 233`
    near-coincidence precedent.

## Run

```
python -m experiments.closure_ledger.scale_count_anchors_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_scale_count_anchors_probe/`.
Expected verdict: `M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO`, 8/8 PASS.
