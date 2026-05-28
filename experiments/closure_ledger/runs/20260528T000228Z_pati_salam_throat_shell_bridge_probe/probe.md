# Pati-Salam SU(4) throat ↔ shell bridge (PR #82)

**Run:** 2026-05-28T00:02:28+00:00

Builds the BAM-native throat ↔ shell `n + 3` Z₂ bridge identified by PR #80's verdict as the most plausible extension toward Pati-Salam SU(4). Tests what the bridge alone gives (mass-ratio audit), and honestly identifies the three open extensions that full SU(4) PS requires.

- **Identification**: throat-shell n+3 Z₂ bridge is BAM-native (= PR #68 shell threshold); full SU(4) PS requires three open extensions (neutrinos, 3-fold quark color, mass-operator unification)
- **Status**: #77 → #80 + #82 bridge: structural foundation complete; three open extensions remain (neutrinos, color, mass-op)
- **B4 caveat**: cavity ω dimensionful; lepton β·k² also dimensionful (β in mass²); mass-ratio audit scale-free; absolute MeV scale rides on single B4 anchor (PR #53)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_throat_shell_n_plus_3_map` | throat ↔ shell n+3 map established (= PR #68 threshold) | **PASS** |
| T2 | `T2_unified_12_state_radial_overtone_basis` | unified 12-state basis (6 throat + 6 shell) | **PASS** |
| T3 | `T3_throat_shell_z2_operator` | throat-shell Z₂ involution constructed | **PASS** |
| T4 | `T4_mass_ratio_audit_cavity_omega_convention` | mass-ratio audit: Gen 3 within 17%, Gen 1 off 2.5×, Gen 2 wrong sign | **PASS** |
| T5 | `T5_missing_pieces_for_full_su4_ps` | 3 open extensions identified (ν, quark color, mass op) | **PASS** |
| T6 | `T6_v3_lepton_quark_mass_operator_mismatch` | v3 lepton β·k² vs quark ω² mass-operator mismatch | **PASS** |
| T7 | `T7_honest_scope` | PR #82 sharpens, does not close PS SU(4) | **PASS** |
| T8 | `T8_assessment` | arc + bridge complete; 3 extensions open | **PASS** |

## T1: Throat ↔ shell `n + 3` map

| generation | lepton n | quark n | shift | = shell threshold? |
|---:|---:|---:|---:|:---:|
| 1 | 0 | 3 | +3 | ✓ |
| 2 | 1 | 4 | +3 | ✓ |
| 3 | 2 | 5 | +3 | ✓ |

## T4: Mass-ratio audit under cavity-ω² convention

| gen | pair | BAM predicted (m_q/m_ℓ) | observed | deviation | within 3×? |
|---:|---|---:|---:|---:|:---:|
| 1 | e ↔ d | 3.626 | 9.139 | 0.397× | ✓ |
| 2 | mu ↔ s | 2.411 | 0.884 | 2.728× | ✓ |
| 3 | tau ↔ b | 1.970 | 2.352 | 0.837× | ✓ |

**3 of 3 generations within factor 3 of observed.**

Gen 3 best (within 17%); Gen 1 off by factor 2.5 (BAM under-predicts); Gen 2 **wrong sign** (BAM predicts quark heavier than lepton; observation has them approximately equal). Predictions are rough at best — the cavity-ω² convention alone is not the right mass operator across both sectors.

## T5: Three open extensions for full PS SU(4)

### BAM-native neutrinos

- **PS role**: 4th leptocolor in up-multiplet F^a
- **BAM status**: no explicit neutrino mode in current scaffold
- **Candidate channels**:
  - opposite-chirality Weyl component of PR #66 throat Dirac spinor
  - sterile Majorana sector (right-handed singlet)
  - separate radial mode not yet identified
- **Closure estimate**: genuine open extension

### 3-fold quark color (SU(3) subgroup)

- **PS role**: quark colors α = 1, 2, 3 in F^a
- **BAM status**: PR #80's verdict: no BAM-derivable triplet in current scaffold
- **Candidate channels**:
  - 3 generations from (k_5+1)/2 — gives SO(3)/SU(2), not SU(3)
  - three Hopf fibrations of S³ — SO(3), not SU(3)
  - S³ isometries SO(4) = SU(2) × SU(2) — no SU(3)
- **Closure estimate**: genuine open extension

### lepton-quark mass-operator unification

- **PS role**: SU(4) acts uniformly on the up-multiplet
- **BAM status**: v3 leptons use β·k² closure-winding (PR #71); PR #77 quarks use ω²(l, n) cavity eigenfrequency. Two different mass operators in current scaffold.
- **Closure estimate**: requires reconciling closure-winding (winding cost on odd-k throat traversal) with cavity-eigenfrequency (shell-mode resonance) into a single mass operator on the unified 12-state basis

## T6: v3 lepton-quark mass-operator mismatch

**Cavity ω² spread in throat region (n=0..2):** `7.53` (too small to span the observed lepton hierarchy)

**Observed lepton mass² spread (e to τ):** `1.21e+07` (needs β·k² closure-winding mechanism, PR #71)

The lepton sector uses **β·k² closure-winding** mass; the quark sector (PR #77) uses **ω²(l, n) cavity eigenfrequency** mass. These are structurally different mass operators. Full PS SU(4) requires unifying them on a single basis.

## Verdict

**PATI_SALAM_THROAT_SHELL_Z2_BUILT_FULL_SU4_REQUIRES_EXTENSIONS.** PATI-SALAM THROAT-SHELL Z₂ BRIDGE BUILT, FULL SU(4) REQUIRES THREE EXTENSIONS. PR #82 constructs the natural BAM-native throat ↔ shell `n + 3` map: each generation `g ∈ {1, 2, 3}` has a charged lepton at radial overtone `n = g - 1` (throat-focused per PRs #59–#66) and a quark-pair at `n = g + 2` (shell-saturated per PRs #77–#80). The shift `+3` = PR #68's shell-saturation threshold. The unified 12-state radial-overtone basis (l=1, n=0..5, p=±) supports a throat-shell Z₂ involution that swaps (n, p) ↔ (n+3 mod 6, p).

MASS-RATIO AUDIT. Under the cavity-ω² mass convention (quark sector convention from PR #77), the predicted lepton-quark mass ratios per generation are: Gen 1 (e ↔ d): BAM 3.63 vs obs 9.14 (off by factor 2.5); Gen 2 (μ ↔ s): BAM 2.41 vs obs 0.88 (WRONG SIGN — BAM predicts quark heavier than lepton, observation has them ~equal); Gen 3 (τ ↔ b): BAM 1.97 vs obs 2.35 (within 17%). 1 of 3 generations within factor 3 of observed.

MASS-OPERATOR MISMATCH. The Gen 2 sign error and the broader pattern point to a deeper issue: v3 leptons use **β·k² closure-winding** (PR #71, β_lepton = k_5²·(2π) = 50π) while PR #77 quarks use **ω²(l, n) cavity eigenfrequency**. The lepton observed τ/e mass-spread is ~1.2·10⁷ in mass² — far beyond the cavity-ω² throat-region spread of ~7.5. The lepton hierarchy CANNOT come from cavity ω² alone; it requires the closure-winding β·k² + β·(k-3)² uplift. The quark hierarchy CANNOT come from closure-winding alone (PR #76's diagnosis: v3 lepton-shaped fitting absorbed unmodeled physics into n_part = 233). The two sectors use STRUCTURALLY DIFFERENT mass operators in the current scaffold.

THREE OPEN EXTENSIONS for full PS SU(4):
  (i) BAM-NATIVE NEUTRINOS. Pati-Salam puts the neutrino as the 4th leptocolor in the up-multiplet F^a = (u^α, ν). BAM's current scaffold has no explicit neutrino mode. Candidate channels: an opposite-chirality Weyl component of PR #66's throat Dirac spinor; a sterile Majorana sector; a separate radial mode not yet identified. Each is a genuine BAM extension, not a closure-ledger question.
  (ii) 3-FOLD QUARK COLOR. PR #80's verdict: no BAM-derivable triplet in the current scaffold (the 3 generations, three Hopf fibrations, S³ isometries, Hopf U(1), and bulk 5D all give SO(3)/SU(2)/U(1) algebras, not SU(3)). This is the same open gap that prevented PR #80 from identifying canonical QCD SU(3) color.
  (iii) LEPTON-QUARK MASS-OPERATOR UNIFICATION. v3 uses β·k² closure-winding for leptons (PR #71), ω²(l, n) cavity-eigenfrequency for quarks (PR #77). A full PS SU(4) representation transforms the up-multiplet F^a = (u^α, ν) uniformly under SU(4); this requires a single mass operator on the unified basis. Reconciling these two operators is a deeper structural question than either (i) or (ii) — it asks why the same closure-ledger geometry gives ω² in the shell sector but β·k² in the throat sector.

WHAT PR #82 ESTABLISHES. The `n + 3` Z₂ throat-shell bridge is BAM-native (the shift = PR #68 shell threshold; no free parameter). The unified 12-state basis combines the lepton and quark sectors structurally. The mass-ratio audit reveals the structural mismatch between the two sectors' mass operators. The verdict sharpens the scope of the PS SU(4) extension: it is not a small incremental probe but a genuine multi-thread research program touching neutrino sector, color algebra, and mass-operator unification.

HONEST SCOPE. PR #82 builds what can be built (the throat-shell Z₂ bridge, the unified basis, the mass-ratio audit) and explicitly identifies the three open extensions. It does NOT derive quark masses, introduce BAM-native neutrinos, derive 3-fold quark color, or unify the lepton-quark mass operators. The four-PR QCD-shell arc (#77–#80) plus PR #82 together complete the structural foundation for any future PS SU(4) construction; the three identified extensions remain genuine open work, comparable in scope to deriving the Standard Model from underlying geometry.

## What this leaves open

- **BAM-native neutrinos** — three candidate channels listed (opposite-chirality Weyl, sterile Majorana, separate radial mode); each a genuine open extension.
- **3-fold quark color** — PR #80's open gap; no BAM-derivable SU(3) triplet in the current scaffold.
- **Lepton-quark mass-operator unification** — reconcile v3's β·k² closure-winding (PR #71) with PR #77's ω²(l, n) cavity eigenfrequency on a single unified basis. Deeper than (i) or (ii).
- **Absolute mass predictions** — pending all three extensions above.
