# Pati-Salam SU(4) throat ↔ shell bridge (PR #82)

The natural extension flagged by PR #80's verdict: unify the lepton
sector (throat traversal modes, PRs #59–#66) and the quark sector
(shell waveguide modes, PRs #77–#80) under a single SU(4) algebra in
which lepton number is a 4th leptocolor.

## What Pati-Salam SU(4) needs

Standard Pati-Salam has 8 fermions per generation per chirality:

  - Up-multiplet `F^a = (u^α, ν)` — 3 quark colors α + 1 lepton =
    4-plet of SU(4).
  - Down-multiplet `F^a' = (d^α, e)` — 3 quark colors α + 1 charged
    lepton = 4-plet of SU(4).

Per chirality across 3 generations: 24 fermions; both chiralities:
48 fermions.

## What BAM currently has

  - 3 charged leptons (e, μ, τ) at throat-traversal modes `n = 0, 1,
    2` (PRs #59–#66).
  - 6 quarks (u, d, s, c, b, t) at shell-saturated cavity modes
    `n = 3, 4, 5` with Z₂ partition (PRs #77–#80).

Total: 9 states. **No explicit neutrinos.** **No 3-fold quark
color** (PR #80's verdict).

## What this PR DOES build

  - **The throat ↔ shell `n + 3` map.** Each generation
    `g ∈ {1, 2, 3}` has a lepton at `n = g - 1` and a quark-pair at
    `n = g + 2`. The shift `+3` = PR #68's shell-saturation
    threshold. **No free parameter** once the radial-overtone ladder
    is computed.

  - **The unified 12-state radial-overtone basis.** `(l = 1, n =
    0..5, p = ±)`. 6 throat states + 6 shell states.

  - **The throat-shell Z₂ involution** on the unified basis:
    `(n, p) ↔ (n + 3 mod 6, p)`. Permutation, involution, swaps
    throat ↔ shell sectors.

  - **Mass-ratio audit** under the cavity-ω² convention.

## Mass-ratio audit results

Under the assumption that lepton and quark masses both come from
cavity `ω²(l, n)` (the quark convention from PR #77):

| generation | pair | BAM (m_q/m_ℓ) | observed | deviation | within 3×? |
|---:|---|---:|---:|---:|:---:|
| 1 | e ↔ d   | 3.63 | 9.14 | 2.52× under | ✓ |
| 2 | μ ↔ s   | 2.41 | 0.88 | wrong sign  | ✗ |
| 3 | τ ↔ b   | 1.97 | 2.35 | within 17%  | ✓ |

**1 of 3 generations within factor 3 of observed.** Gen 3 best; Gen
2 has the wrong sign (BAM predicts quark heavier than lepton;
observation has them ~equal).

The Gen 2 sign error is not a small numerical issue. It points to a
deeper structural mismatch.

## The mass-operator mismatch

v3 leptons use **`β·k²` closure-winding** mass (PR #71, with
`β_lepton = k_5²·(2π) = 50π` derived structurally):

```
m²_lepton ~ baseline + β·k² + β·(k-3)²·(uplift for k=5)
```

PR #77 quarks use **`ω²(l, n)` cavity-eigenfrequency** mass:

```
m²_quark ~ ω²(l, n) + χ_n(p) + H_couple corrections
```

The lepton observed `(τ/e)²` mass-spread is ~`1.2·10⁷`. The
cavity-ω² spread in the throat region is only ~`7.5`. **Cavity ω²
alone cannot give the lepton hierarchy** — the lepton sector
fundamentally requires closure-winding `β·k²`.

So BAM's lepton and quark sectors use **structurally different**
mass operators in the current scaffold. A full Pati-Salam SU(4)
representation transforms the up-multiplet `F^a = (u^α, ν)`
uniformly under SU(4); this requires a single mass operator on the
unified basis. That requires reconciling closure-winding (winding
cost on odd-`k` throat traversal) with cavity-eigenfrequency
(shell-mode resonance) into a single mass operator.

## Three open extensions for full PS SU(4)

  1. **BAM-native neutrinos.** Pati-Salam puts the neutrino as the
     4th leptocolor in the up-multiplet. BAM has no explicit
     neutrino mode in the current scaffold. Candidate channels:

       - The opposite-chirality Weyl component of PR #66's throat
         Dirac spinor (the SUSY factorization gives two
         SUSY-partner sectors per mouth — one might be the
         neutrino sector).
       - A sterile Majorana sector (right-handed singlet).
       - A separate radial mode not yet identified.

     Each is a genuine BAM extension, not a closure-ledger
     machinery problem.

  2. **3-fold quark color.** PR #80's verdict: no BAM-derivable
     triplet in the current scaffold. All natural triplet
     candidates (3 generations, three Hopf fibrations, S³
     isometries, Hopf U(1), bulk 5D) give SO(3)/SU(2)/U(1)
     algebras, not SU(3).

  3. **Lepton-quark mass-operator unification.** Reconcile the v3
     lepton mass operator `β·k²` (PR #71) with the quark mass
     operator `ω²(l, n)` (PR #77) into a single operator on the
     unified 12-state basis. This is structurally deeper than (1)
     or (2): it asks why the same closure-ledger geometry gives
     `ω²` in the shell sector but `β·k²` in the throat sector.

## What PR #82 establishes (and does not)

**Establishes:**

  - The throat-shell `n + 3` Z₂ bridge as a BAM-native structural
    ingredient (no free parameter).
  - The unified 12-state radial-overtone basis.
  - Explicit identification of the three open extensions any full
    PS SU(4) construction must close.
  - Mass-ratio audit revealing the structural mismatch.

**Does NOT establish:**

  - Quark masses or absolute predictions.
  - BAM-native neutrinos.
  - 3-fold quark color.
  - Lepton-quark mass-operator unification.

The verdict is therefore **`PATI_SALAM_THROAT_SHELL_Z2_BUILT_FULL_
SU4_REQUIRES_EXTENSIONS`** — PR #82 sharpens the scope of the
Pati-Salam extension; it does not close it. The QCD-shell arc
(#77–#80) plus PR #82 together complete the structural foundation;
the three identified extensions remain genuinely open and are each
substantial research programs in their own right.

## Honest scope

  - **Is:** the throat-shell `n + 3` Z₂ bridge; the unified 12-state
    basis; the cavity-ω² mass-ratio audit; explicit identification
    of three open extensions.

  - **Is not:** a derivation of quark masses; an introduction of
    BAM-native neutrinos; a derivation of 3-fold quark color; a
    unification of the two mass operators; absolute mass
    predictions.

## B4 accounting

Cavity ω is dimensionful (`1/length`); the lepton `β` is also
dimensionful (mass² units). Mass ratios are scale-free. The mass
operators (closure-winding vs cavity ω²) both rely on the single
B4 anchor (PR #53). Different functional forms; same dimensional
content.

## Tests

  T1. Throat-shell `n + 3` map established; shift equals PR #68's
      shell threshold for every generation.
  T2. Unified 12-state radial-overtone basis (l=1, n=0..5, p=±)
      constructed; 6 throat + 6 shell.
  T3. Throat-shell Z₂ involution operator on the unified basis;
      verified as permutation, involution, swaps throat ↔ shell.
  T4. Mass-ratio audit under cavity-ω² convention: Gen 3 within
      17%, Gen 1 off by factor 2.5, Gen 2 wrong sign.
  T5. Three open extensions identified explicitly: (i) BAM-native
      neutrinos, (ii) 3-fold quark color, (iii) mass-operator
      unification.
  T6. v3 lepton-quark mass-operator mismatch: explicit
      demonstration that cavity-ω² spread in throat region (~7.5)
      cannot give observed `(τ/e)²` mass-spread (~1.2·10⁷); leptons
      structurally require `β·k²` closure-winding.
  T7. Honest scope: PR #82 sharpens, does not close.
  T8. Assessment.

## Verdict structure

  - **`PATI_SALAM_THROAT_SHELL_Z2_BUILT_FULL_SU4_REQUIRES_EXTENSIONS`**
    (expected): T1–T7 pass; the `n + 3` Z₂ bridge is built; three
    open extensions identified.

  - **`PATI_SALAM_INCONCLUSIVE`**: a structural test fails;
    investigate before relying on the bridge in downstream work.

## What this leaves open

  - **BAM-native neutrinos** — three candidate channels listed;
    each a genuine open extension.
  - **3-fold quark color** — PR #80's open gap.
  - **Lepton-quark mass-operator unification** — `β·k²` vs `ω²`.
  - **Absolute mass predictions** — pending all three above.

## Cross-references

  - `docs/color_algebra_shell_research_plan.md` — PR #80, which
    flagged Pati-Salam SU(4) as the natural extension.
  - `docs/throat_to_shell_transition_research_plan.md` — PR #68,
    the shell-saturation threshold = `n + 3` shift.
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    lepton β·k² closure-winding mass operator.
  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77,
    the quark ω²(l, n) cavity mass operator.
  - `docs/throat_dirac_spinor_research_plan.md` — PR #66, the
    SUSY factorization with two SUSY-partner sectors (candidate
    for neutrino channel).
  - `experiments/closure_ledger/pati_salam_throat_shell_bridge_probe.py`
    — this probe.
