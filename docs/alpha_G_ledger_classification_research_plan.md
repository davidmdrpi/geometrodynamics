# Where do α and G sit in the #104 epistemic ledger? (PR #105)

PR #104 classified the mass-sector results into five tiers and reduced
BAM's dimensionful content to two B4 anchors (`m_e`, `√σ`). This probe
answers a sharper question: under that ledger, are the fine-structure
constant `α` and Newton's constant `G` **derived**, **residual**, or
**anchors**? The answer also sharpens #104.

## The fundamental-constant baseline

| constant | tier | basis |
|---|---|---|
| **c** | unit convention | sets length = time |
| **ℏ** | DERIVED (geometric) | the closure quantum `2π`; `ℏ = m_e·R_MID·c` |
| **G** | DIMENSIONFUL ANCHOR (B4) | the GR-foundational scale; root of `m_e`, `√σ` |
| **α** | UNIVERSAL RESIDUAL | value (1/137) a free input; only the running derived |

## G is an ANCHOR (the gravitational, foundational one)

BAM is GR-foundational: the throat is a gravitational wormhole, and its
size — the invariant bulk length `ΔR = R_OUTER − R_INNER` (equivalently
`R_MID`), the ONE dimensionful input the B4 theorem (PR #52/#53) says is
mandatory — is set by the bulk gravity sector. The brane tension that
fixes the throat equilibrium is the Randall–Sundrum-tuned
`λ_crit = √(6|Λ₅|)/κ₅²` (PR #57), so the throat scale descends from the
5D gravitational coupling `κ₅` (∝ G) and the bulk `Λ₅`. So **G is the
dimensionful anchor** — the GR-foundational scale, not derivable within a
theory whose foundation is gravity (it sets the units of the geometry).

It is moreover the **root** anchor: both #104 sector anchors descend from
it — `m_e = ℏc/R_MID` (R_MID set by the gravity-tuned brane tension) and
`√σ` (with `σ ∝ √|Λ₅|/κ₅²`). Whether `m_e`, `√σ`, `G` are three
independent scales or one gravitational scale plus dimensionless ratios
is the open "how many scales" question — but in every case G is Tier 2
(anchor), not derived.

## α is a RESIDUAL (universal)

Throughout BAM, `α` is a NUMERICAL INPUT: the EM self-energy
`A_EM = α·ℏc/2`, the capped self-energy `U_EM/(mc²) = α/2` (PR #55), and
the one-loop anomaly `a = α/2π` (PR #62). BAM derives the STRUCTURE
around `α` — the charge UNIT (`|c₁| = 1`, charge quantization), the loop
measure `1/2π` (PR #74), the `α/2` self-energy form — but NEVER the VALUE
1/137. As in the Standard Model and every current framework, only `α`'s
RUNNING is derived (the QFT β-function); its value is a free input — the
"137 problem". So `α` is a dimensionless **residual**, and a UNIVERSAL
one: not a BAM-specific shortfall but the same open input every theory
takes. It sits with the flavor puzzle, not with the BAM-specific
residuals (`ε`, `n_part`). BAM's geometric aspiration to fix `α` from the
S³/Hopf structure is open and unfulfilled.

## The refinement of #104

PR #104's "derived geometry" tier used `α` as a silent input (`a = α/2π`
derives the `1/2π` GIVEN `α`), so `α` is properly a **residual input** to
that tier, not itself derived. And the two #104 anchors (`m_e`, `√σ`)
both descend from the gravitational scale `G`, the **root** anchor.
Sharpened ledger: **G — the foundational dimensionful anchor; α — a
universal dimensionless residual (the value, not the running); ℏ —
geometric (the closure quantum); c — units.**

## Tests

| # | test | finding |
|---|---|---|
| T1 | the question | classify `α` and `G` under the #104 ledger |
| T2 | baseline | `c` = units; `ℏ` = closure quantum (`ℏ=m_e·R_MID·c`), geometric |
| T3 | `G` = anchor | throat size (B4 length) set by bulk gravity (`λ_crit=√(6\|Λ₅\|)/κ₅²`, PR #57) |
| T4 | `α` = input | `A_EM=α·ℏc/2`, `a=α/2π`; structure derived, value not |
| T5 | `α` = universal residual | running derived; value (137) a free input, as in the SM |
| T6 | refines #104 | `α` a residual input; `G` the root anchor |
| T7 | classification table | ℏ geometric, c units, G anchor, α universal residual |
| T8 | assessment | `G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** `G` is the dimensionful (gravitational,
    foundational) anchor — the GR scale the B4 length and the #104 sector
    anchors descend from; `α` is a universal dimensionless residual — BAM
    derives the charge unit, the EM structure, and `α`'s running, but the
    value 1/137 is a free input as in every framework; `ℏ` is geometric
    (the closure quantum); `c` is units.

  - **Open:** the geometric derivation of `α` (the "137 problem"); a
    first-principles `G` (not derivable within a gravity-foundational
    theory); and how many truly-independent dimensionful scales there are
    (whether `m_e`, `√σ`, `G` collapse toward one gravitational anchor).

## Cross-references

  - `docs/program_synthesis_research_plan.md` — PR #104, the five-tier
    ledger this places `α` and `G` into.
  - `docs/hbar_origin_status.md`, `docs/hbar_origin_note.md` — the ℏ /
    B4-anchor origin (`ℏ = m_e·R_MID·c`).
  - `docs/brane_tension_tuning_research_plan.md` — PR #57, the bulk-gravity
    relation `λ_crit = √(6|Λ₅|)/κ₅²`.

## Run

```
python -m experiments.closure_ledger.alpha_G_ledger_classification_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_alpha_G_ledger_classification_probe/`.
Expected verdict: `G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL`, 8/8 PASS.
