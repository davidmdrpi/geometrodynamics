# Program-wide synthesis: the BAM input budget and epistemic ledger (PR #104)

A capstone meta-probe. After the lepton/QED foundations (B1–B5,
throat-as-particle), the neutrino arc (#85–#96), the quark arc
(#76–#98), and the QCD confinement / exotic arc (#99–#103), this probe
classifies every major result by EPISTEMIC STATUS and counts BAM's
actual free/anchored inputs.

## The five epistemic tiers

| tier | what | count |
|---|---|---:|
| 1. derived geometry | parameter-free / reconstructed structural results | ~22 |
| 2. **dimensionful anchors (B4)** | the mandatory scale inputs | **2** |
| 3. open dimensionless residuals | constrained but not derived, localized | 2 |
| 4. the flavor puzzle | the Yukawa hierarchy — universal, not BAM-specific | 1 |
| 5. topological predictions | non-orientable, BAM-specific, falsifiable | 6 |

## The input budget (headline)

  - **2 dimensionful anchors** — `m_e = ℏc/R_MID` (QED/lepton scale) and
    `√σ ≈ Λ_QCD` (confinement scale). The B4 scale-modulus theorem
    (PR #52) requires one dimensionful input per sector, so **two is the
    irreducible minimum — BAM sits at it.**
  - **~2 open dimensionless residuals** — the neutrino boundary compliance
    `ε` (the seesaw/bounce residual, itself bracketed to `[2π, k_5√(2π)]`,
    PR #88–#90) and the quark `n_part = 233` (a compensator for the flavor
    puzzle, PR #76/#97).
  - **1 universal open problem** — the flavor puzzle (the quark Yukawa
    hierarchy: RG-invariant ratios, irregular, derivable by no current
    theory, PR #98). Shared with all of physics; not BAM-specific.
  - Everything else is **derived geometry** (~22 results) or a
    **topological prediction** (~6).

## Tier 1 — derived geometry (the bulk)

`|c₁| = 1` charge quantization; spin-½ from Hopf holonomy; the Coulomb
law from throat flux; `g = 2`; the one-loop Schwinger `a = α/2π` (0.15%);
`C` = inner/outer swap, CPT, the throat Dirac spinor; the pair-threshold
factor 2 (`Σc₁=0`); `k_5 = dim(S³)+2 = 5`; `β_lepton = k_5²·2π = 50π`;
`#generations = (k_5+1)/2 = 3`; the unified Bohr-Sommerfeld mass operator;
the neutrino Majorana selection rule (`c₁=0`); PMNS anarchy / CKM
alignment; generic CP + two Majorana phases; six quarks = 3×2 + Z₂; the
Cornell form (flux-tube bridge + throat/gluon exchange); string-breaking
= Schwinger (`eE→σ`); the Regge slope `α′ = 1/(2πσ)`; the exotic-`J^PC`
bookkeeping; the no-exotic-`J^P`-for-baryons result.

## Tier 5 — topological predictions (the distinctive output)

Spanning the full testability gamut:

  - **MATCHED** — the mesonic `1-+` hybrids (`π₁`, `η₁`) (PR #101);
  - **FALSIFIABLE** — neutrino normal ordering, `m_ββ ≲ 8 meV`,
    `Σm_ν ≈ 59–65 meV` (PR #95/#96);
  - **ACCOMMODATED** — the multiquark exotic zoo (`X, Z_c, T_cc, P_c`)
    (PR #101);
  - **CONSTRAINED** — the light baryonic exotics, a counting test
    (PR #102);
  - **FINDABLE** — the heavy Möbius baryon at the flavor-independent `2√σ`
    gap (PR #103);
  - **FREE** — the Möbius glueball tower (glueballs unobserved, PR #100).

## In one line

Two mandatory B4 anchors + a couple of localized open dimensionless
residuals + the universal flavor puzzle, with the rest derived geometry
and a set of falsifiable non-orientable topological predictions.

## Honest scope

This is a CLASSIFICATION, not a new derivation — it inventories and ranks
the program's results by epistemic status. The "derived" tier groups
parameter-free results with reconstructed/fitted structural matches (the
charged-lepton ladder, the Schwinger 0.15% match), without relitigating
each. The two anchors and the flavor puzzle are genuinely open inputs
(B4-mandatory and universal respectively); the localized residuals
(`ε`, `n_part`) are the program's own remaining homework.

## Run

```
python -m experiments.closure_ledger.program_synthesis_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_program_synthesis_probe/`.
Expected verdict: `BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS`, 8/8 PASS.
