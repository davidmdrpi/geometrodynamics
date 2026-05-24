# Even-k absence in the charged-lepton throat sector

The THESIS-flagged "highest-leverage near-term result." The charged
leptons sit at odd depths `k ∈ {1, 3, 5}` (e, μ, τ); even `k` is absent.
The odd-k closure lemma (`docs/odd_k_closure_lemma.md`) already showed
(i) the Layer-1 ledger closes mod 2π for *any* integer `k`, so the
selection is **not** arithmetic, and (ii) odd `k` is the
orientation-reversing closure across the non-orientable throat, framed as
a "choice of sector." This probe **upgrades that to a classification /
selection rule**: even-`k` modes are not charged leptons because the
charged leptons are **spin-½ fermions** (established across PRs #59–#66),
and only odd `k` realizes the fermionic (non-orientable, orientation-
reversing) closure. Even `k` is the orientation-preserving (orientable /
bosonic) sector.

## The monodromy classification

Each throat pass applies the transport `T = iσ_y` (`T² = −I`, B2). After
`k` passes the spinor monodromy is `T^k`, with period-4 pattern

```
k:    1        2     3        4     5        6
T^k:  iσ_y    −I    −iσ_y    +I    iσ_y    −I
```

The decisive feature is **off-diagonal vs diagonal**, i.e. the Z₂
partition class:

  - **Odd `k` → `T^k` off-diagonal** (`±iσ_y`): the spinor is mapped to
    the **opposite** Z₂ class — the antipodal/deck partner (`p ~ −p` on
    `RP³ = S³/Z₂`). This is the **orientation-reversing** closure across
    the non-orientable throat — the non-trivial spin structure
    (`T² = −I`, B2), i.e. a genuine **spin-½ fermion**.

  - **Even `k` → `T^k` diagonal** (`±I`): the spinor returns to the
    **same** Z₂ class. This is the **orientation-preserving** closure on
    the orientable double cover `S³` — the trivial spin structure, an
    integer-class / **bosonic** closure.

So `k mod 2` is exactly the orientability (= spin-statistics) grading of
the closure.

## The selection rule

The charged lepton is a **Dirac spin-½ fermion** — established
geometrically across the recent arc: the throat Dirac 4-spinor
(`throat_dirac_spinor_probe`, #66), `g = 2` (#61), the Wigner rotation
(#60), `T² = −I` / CPT (#65). A spin-½ fermion **requires** the
non-orientable, orientation-reversing closure (the `T² = −I` spin
structure, B2). Therefore:

```
charged-lepton throat  ⟺  spin-½ fermion  ⟺  orientation-reversing
                       ⟺  T^k off-diagonal  ⟺  k odd .
```

Even `k` (`T^k` diagonal, orientation-preserving, bosonic) is **excluded
from the charged-lepton sector** — not arithmetically (the ledger closes
for any `k`), but by **spin-statistics**: an even-`k` closure is not a
spin-½ fermion. The lepton depths `(1, 3, 5)` are all odd; even `k` is
absent.

## What even-k modes *are*

Even-`k` modes are not geometrically forbidden — they are the
orientation-preserving closures on the orientable double cover `S³`, the
integer-class / bosonic sector of the same geometry. They simply are not
charged leptons (which are the non-orientable spin-½ throats). The
absence is a **classification**: charged leptons populate the odd
(fermionic) Z₂ class only.

## Unification with the spin / discrete-symmetry arc

Odd-`k` (non-orientable, `T² = −I`) is the *same* fermionic structure
that runs through the whole recent arc: the spinor double cover
(#60 Wigner rotation), `g = 2` (#61), the CPT `T² = −I` (#65), and the
throat Dirac spinor (#66). The even-`k` absence is the **spin-statistics
face** of the fermionic throat — one structure (`T² = −I`, B2) seen from
the generation-counting side.

## B4 accounting

`k` is a **dimensionless integer** (a winding / depth); the classification
is **topological** (the Z₂ orientability grading), independent of the
single anchor `m_e`. The mass *values* at each odd `k` carry the scale;
the odd-only *selection* does not.

## Tests

  T1. **Monodromy `T^k`.** `T = iσ_y`, `T² = −I`; `T^k` off-diagonal
      (odd) vs diagonal (even) — the Z₂ orientability grading `k mod 2`.
  T2. **Odd k = non-orientable fermionic closure.** off-diagonal `T^k`
      ⟹ opposite Z₂ class ⟹ orientation-reversing ⟹ spin-½ (`T²=−I`,
      B2). Even k = orientation-preserving / bosonic.
  T3. **Selection rule.** charged lepton = spin-½ fermion (#59–#66) ⟹
      orientation-reversing ⟹ odd `k`; the depths `(1,3,5)` all odd.
  T4. **Even-k classified.** even `k` = orientable double-cover /
      bosonic sector — exists, but not a charged lepton.
  T5. **Not arithmetic (ledger closes for any k).** `Φ_avail(k) =
      2π(k+1) + 50π·max(0,k−3)² ≡ 0 mod 2π` for every integer `k` — the
      selection is spin-statistics, not closure arithmetic.
  T6. **Unification.** odd-k `T²=−I` = the same fermionic structure as
      #60/#61/#65/#66 — the spin-statistics face of the fermionic throat.
  T7. **Falsification / B4.** even-k leptons (bosons) would falsify;
      BAM passes (leptons spin-½, odd-k). `k` dimensionless/topological.
  T8. **Assessment.**

## Verdict structure

  - **EVEN_K_EXCLUDED_BY_SPIN_STATISTICS** (expected): even-`k` modes are
    absent from the charged-lepton sector by a spin-statistics selection
    rule. `k mod 2` is the orientability grading of the throat closure
    (`T^k` off-diagonal/odd = orientation-reversing/spin-½ fermion;
    diagonal/even = orientation-preserving/bosonic); charged leptons,
    being spin-½ Dirac fermions (#59–#66), are the odd class. The ledger
    closes for any `k`, so the selection is topological/spin-statistical,
    not arithmetic. Even-`k` is the orientable bosonic sector, not a
    charged lepton. The result unifies with the fermionic-throat arc.

  - **CLASSIFICATION_INCOMPLETE**: the monodromy grading or the
    spin-statistics link fails.

## What this leaves open

  - **The even-k (bosonic) spectrum.** Whether the orientable even-`k`
    closures host a physical (e.g. integer-spin) spectrum is not worked
    out — only that they are not charged leptons.
  - **Generation count (why exactly 3).** The odd-k rule allows
    `k = 1, 3, 5, 7, …`; why the charged leptons stop at `k = 5` (three
    generations) is a separate question (the β-uplift / closure cutoff).

## Cross-references

  - `docs/odd_k_closure_lemma.md` — the closure arithmetic + the
    sector-choice framing this probe upgrades.
  - `docs/throat_dirac_spinor_research_plan.md` — the throat Dirac spinor
    (#66).
  - `docs/cpt_dirac_operator_research_plan.md` — `T² = −I` (#65).
  - `docs/topological_discrete_sector_research_plan.md` — the antipodal
    `Z₂` / `T = iσ_y` (B2).
  - `geometrodynamics/embedding/transport.py` — `derive_throat_transport`
    (`T = iσ_y`).
  - `experiments/closure_ledger/even_k_absence_probe.py` — this probe.
