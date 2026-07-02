# The geon-statistics adjudication — companion probe (PR #196)

**Run:** 2026-07-02T04:37:32+00:00

The deliverable is `docs/geon_statistics_adjudication.md` — mathematics against the Sorkin-school literature on the actual BAM topology. This probe machine-checks the document's finite/algebraic steps. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the adjudication; the document is the argument | **PASS** |
| T2 | `T2_lemma1_orientation_arithmetic` | Lemma 1: the throat prime is the RP³ geon (orientation arithmetic) | **PASS** |
| T3 | `T3_lemma2_hypotheses` | Lemma 2: prime ✓ non-chiral ✓ abelian ✓ (SSC hypotheses hold) | **PASS** |
| T4 | `T4_pin_lift_of_rotation` | Lemma 3/5: bare rotation trivial (cited); pin lift = −I (computed) | **PASS** |
| T5 | `T5_two_throat_mcg_sectors` | Lemma 4: G = (Z₂∗Z₂)⋊S₂; 4 scalar + indefinite sectors | **PASS** |
| T6 | `T6_soh_selection` | SOH: SSC ⟹ bare GR → Bose; Pin⁻ framing → Fermi (= #188) | **PASS** |
| T7 | `T7_cobordism_ledger` | #58 cobordism ledger: 3 closed, 1 honest open construction | **PASS** |
| T8 | `T8_assessment` | label change: 'Pauli from GR + the forced Pin⁻ framing' | **PASS** |

## The sector structure (two throats)

| slides | exchange | sector |
|---:|---:|---|
| +1 | +1 | Bose |
| +1 | -1 | Fermi |
| -1 | +1 | Bose |
| -1 | -1 | Fermi |
| (2-dim family) | mixed | indefinite (Sorkin–Surya anomalous) |

## The selection

- SSC hypotheses: {'prime': True, 'non_chiral': True, 'abelian': True} → the correlation is a **theorem**
- bare metric GR: rotation +1 → **Bose**
- Pin⁻-framed GR: rotation -1 → **Fermi** (consistent with #188: True)

## The #58 cobordism ledger

- **mirror_pair_identity** [✓]: amphichirality (T3): the created mirror partner is diffeomorphic to the throat — identical geons
- **bordism_existence** [✓]: Omega_3^SO = 0 (every closed orientable 3-manifold bounds): 4-manifolds S3 -> S3 # RP3 # RP3 exist (cited)
- **spin_structure_extension** [✓]: RP3 is spin (orientable Lie group SO(3), parallelizable); Omega_3^Spin = 0: the pin/spin framing extends over a bounding 4-manifold (cited)
- **explicit_BAM_4manifold** [OPEN]: the Dowker-Sorkin construction is over the R3 background; the S3 transplant (connected-sum locality) with their causal/Morse conditions has NOT been written down - the honest open construction, with no topological obstruction found

## Verdict

**SSC_THEOREM_HYPOTHESES_HOLD_SIGN_IS_FRAMING_DEPENDENT_BARE_GR_SELECTS_BOSE_PIN_MINUS_FRAMING_SELECTS_FERMI.** THE ADJUDICATION (the argument is in docs/geon_statistics_adjudication.md; this probe checks its arithmetic).

WHAT IS A THEOREM. The BAM throat prime is RP³ — the RP³ geon, from the repo's own #169 quotient (orientation arithmetic verified) — and it passes ALL THREE hypotheses of the Dowker–Sorkin spin-statistics theorem: PRIME (elliptic, irreducible), NON-CHIRAL (q² ≡ −1 mod 2; the explicit reflection descends), ABELIAN (π₁ = Z₂). The spin-statistics CORRELATION is therefore a theorem for pair-created BAM throats — the strongest topological proof point available.

WHAT THE SIGN DEPENDS ON. RP³ is NON-SPINORIAL in bare Diff (Hendriks: cyclic-π₁ primes; the #170/#171 'spinorial' sentence is hereby corrected): in bare metric GR the SSC selects BOSE, and in frozen-topology canonical GR the statistics is a sector choice among Bose, Fermi, and a continuum of indefinite sectors (G = (Z₂∗Z₂)⋊S₂ — four scalar sectors + the anomalous 2-dim family, all verified in a faithful model). The −1 is NOT a bare-GR theorem. In PIN⁻-FRAMED GR — forced by the non-orientable RP² slice (#169), the uniqueness of Pin⁻ on it (#170), and the pin-spinor matter content (#195) — the trivialized rotation lifts to −I (the SU(2) lift of the π₁(SO(3)) generator, computed), the throat is spinorial WITH FRAMING, and the same SSC theorem selects FERMI — the #188 holonomy, correctly reinterpreted as the pin lift.

THE #58 CHANNEL. Mirror-pair identity ✓ (amphichirality), bordism existence ✓ (Ω₃^SO = 0), framing extension ✓ (RP³ spin, Ω₃^Spin = 0); the explicit BAM 4-manifold (the Dowker–Sorkin construction transplanted from ℝ³ to the S³ background) is the honest open construction — no obstruction found. LABEL CHANGE: 'Pauli from GR' → 'Pauli from GR + the forced Pin⁻ framing.'
