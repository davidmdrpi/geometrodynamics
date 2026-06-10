# Residual-bracket synthesis and input-budget ledger (PR #150)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The ledger
> accounts for the inputs of the matter QFT the classical geometry
> reconstructs.

The synthesis capstone for the program's input accounting. The five-tier
budget (#104), the constants placement (#105/#106), the anti-numerology
negative results (#107/#108), the APS matter-partition ledger (#123–#125),
the α ledger (#143), the bulk-scale ledger and audit (#133/#148), and the
flavor audits (#113/#134/#149) each categorized one piece. This PR assembles
them into **one categorized input budget**, re-verifies a keystone from every
category (the #131 capstone convention), adds the consolidated table to
`docs/THESIS.md`, and makes the program's central bookkeeping claim
checkable: **BAM is not accumulating loose knobs.**

## The categorized budget

| category | item | status | source PRs |
|---|---|---|---|
| Anchor (dimensionful) | `G` (→ `ΔR = 0.52·R_MID`) | mandatory (B4), relocatable | #52/#53/#57/#106/#133 |
| Fixed tuning | `√6` (RS flatness) | derived constant, not a knob | #57 |
| Universal residual | `α ≈ 1/137` | structure/running derived; value scan-excluded | #74/#141–#147; #143 |
| Universal residual | `√σ/m_e ≈ 830` | one-G repackaging derived; value scan-excluded | #106; #107/#108 |
| Program residual | `n_part = 233` | doubling topological (APS); value compensator | #97/#123/#125 |
| Program residual | `ε` (ν compliance) | order-of-mag derived; window `[2π, k₅√(2π)]` | #89/#112 |
| Bracketed sub-residual | `k·r_s` | `(0, 0.0064–0.070]` two-sided | #133/#148 |
| Bracketed sub-residual | `ε_n` spread | `[1.32, 1.44]`/step, ~0.3%; power laws excluded | #113/#149 |
| Universal open problem | flavor puzzle | RG-invariant ⟹ not running | #97/#107/#108/#134 |
| No residual (contrast) | lepton `N = 4k₅² = 100` | structure AND value derived | #124 |

## Keystones re-verified (one per category)

- **Anchor:** `ΔR/R_MID = 0.52` exactly; `√6 = 2.449490` fixed.
- **Universal residuals:** the #143 α scan and #108 √σ/m_e scan re-run — best
  *principled* candidates `k₅³ + 2π` (−4.2%) and `2π·k₅³` (−5.4%); every
  sub-% match needs an ad-hoc term.
- **Program residuals:** the #107 circularity re-derived — `n_part` drift
  216–255 ⟹ `4n_part − 100 ∈ [764, 920]` (±9%) against the fixed 830; the
  lepton contrast `N = 4k₅² = 100`, 3 generations, fully derived.
- **Bracketed sub-residuals:** the #148 spectrum sensitivity re-checked with
  a light eigensolve (`c = 9.86`, matching #148); the #149 required-profile
  inversion re-derived (`ε₃/ε₂ = 1.435` pure-bounce).

## The no-loose-knobs claim, checkable

| PR | inputs added |
|---|---|
| #144 vacuum polarisation / running | 0 |
| #145 Z₁ = Z₂ | 0 |
| #146 charge form factor | 0 |
| #147 F₁/F₂ EM-arc capstone | 0 |
| #148 k·r_s bracket | 0 |
| #149 ε_n bracket | 0 |

Six probes, zero inputs. The budget today is the **same** budget as
#104/#125 — while the derived ledger grew by the full one-loop EM sector and
two bracket audits.

## Scope

A consolidating synthesis: it organizes and re-verifies, it does **not**
remove residuals (#125's honesty). What would change the ledger: deriving α
or √σ/m_e (the scans make closure-number routes unlikely); deriving
`n_part`'s dynamics (the flavor puzzle); pinning ε or the `ε_n` profile
(mixing/anarchy sector, #92); deriving `k·r_s`'s value (the absolute
normalisation, #112).

## Deliverables

- `experiments/closure_ledger/residual_bracket_synthesis_probe.py` (+ run)
- The consolidated table appended to `docs/THESIS.md` ("The categorized
  input budget (PR #150)")
- README claims-table row

## Reproduce

```bash
python -m experiments.closure_ledger.residual_bracket_synthesis_probe
# Verdict: INPUT_BUDGET_ONE_ANCHOR_FOUR_RESIDUALS_TWO_BRACKETED_ZERO_NEW_KNOBS
```
