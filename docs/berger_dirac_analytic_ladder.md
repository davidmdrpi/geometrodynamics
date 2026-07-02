# The analytic Berger–Dirac ladder: closed-form spectral geometry for the {1,3,5} sectors (PR #197)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This document is **closed-form
> analysis**, not simulation: the Dirac spectrum on the Berger-squashed
> 3-sphere is derived exactly (Peter–Weyl representation theory; the
> spectrum itself is classical — Hitchin 1974, Bär's harmonic-spinor
> work — re-derived here self-contained and checkpointed), and the fate
> of the odd-k ladder is read off analytically: which levels persist,
> where every crossing occurs, and whether the k₅ = 5 cutoff has a
> spectral counterpart at λ ≠ 1. The companion probe verifies each
> identity to machine precision.

## 0. Results, up front

For the Berger sphere S³_λ (unit base, Hopf-fiber length ∝ λ; the
#165/#193 conventions), with the throat's actual field content — the
**spinor** (#170/#195/#196) — the full Dirac spectrum consists of two
closed-form families (§2). Sectored by total fiber momentum k = 2m′:

- **The odd-k ladder is the family-A winding tower**

  ```
  m_k(λ) = k/λ + λ/2 ,        k = 1, 3, 5, …
  ```

  with **uniform gaps m_{k+2} − m_k = 2/λ**, positive for every
  λ ∈ (0, ∞): the ladder never reorders and never collapses at any
  squash. Round values: 3/2, 7/2, 11/2.
- **No infinitesimal-squash breakdown — analytically.** The sector
  grounds are smooth in λ with **O(1) sensitivities** at the round
  point: d ln m_k/dλ|₁ = (½−k)/(k+½) = −1/3, −5/7, −9/11. The
  metric-softness of the true spinor ladder (contrast the surrogate's
  −71, PR #192/#194) is now a closed-form fact.
- **Every crossing is located exactly.** On the stretch side two
  boundaries occur, both in closed form:
  - at **λ×(k) = √(2k+4)** (= √6, √10, √14 ≈ 2.449, 3.162, 3.742) the
    sector's lowest level changes character, from the family-A winding
    state to a family-B level descending toward zero;
  - at **λ*(k) = [8(k+1) + 2√(16(k+1)² + k²)]^{1/2}** (≈ 5.668, 8.035,
    9.851) the sector becomes **massless** — a harmonic spinor, the
    Hitchin phenomenon — with the **electron sector collapsing first**.
  The ladder's ordering m₁ < m₃ < m₅ holds on the whole window
  (0, λ*₁ ≈ 5.668); around the round point there is a huge protected
  window (0, √6): every sector ground is the pure winding state.
- **The k₅ = 5 cutoff has NO spectral counterpart at any λ.** The gap
  2/λ is independent of k: the tower continues 7, 9, … with identical
  spacing on the entire Berger family. The three-generation cutoff is
  dynamical (the phase budget, #122/#136), not spectral — now known in
  closed form off the round point, not just at it.

**Adjudication of the refutation edge:** the clean-failure scenario —
"the ladder structure collapses at infinitesimal squash" — is
analytically excluded (smooth closed forms, O(1) derivatives). The clean
pass is delivered with more than asked: the flagship protection claim is
now **spectral fact with exactly located boundaries**, not algebraic
protection (#183) or discretized-eigensolver evidence (#192/#193).

## 1. Setup and conventions

S³ = SU(2), left-invariant frame E₁, E₂, E₃ with [Eᵢ, Eⱼ] = 2ε_{ijk}E_k
(so the bi-invariant metric with Eᵢ orthonormal is the round S³ of
radius 1, scalar curvature 6). The Berger metric g_λ makes
X₁ = E₁, X₂ = E₂, X₃ = E₃/λ orthonormal — fiber length ∝ λ, volume
λ·2π², and the scalar Laplacian is the #165/#193 spectrum
Δ = 4j(j+1) + 4m²(λ⁻²−1) (checked: the conventions match exactly).
Scalar curvature: **scal(λ) = 8 − 2λ²** (round: 6 ✓; positive iff λ < 2).

Spinors: S³ is simply connected — the spin structure is unique. Left
trivialization: ψ: SU(2) → ℂ²; Peter–Weyl decomposes L² into blocks
V_j ⊗ V_j* (dim (2j+1)²), j = 0, ½, 1, …; left-invariant fields act on
the right factor as π_j(Eᵢ) with Casimir −4j(j+1). The Dirac operator of
g_λ in this trivialization (Levi-Civita spin connection of a
left-invariant metric; Koszul/Milnor connection coefficients
Γ₁₂₃ = Γ₂₃₁ = λ, Γ₃₁₂ = 2/λ − λ):

```
𝒟_λ  =  (σ₊ J₋ + σ₋ J₊)  +  (2/λ) σ₃ J₃  +  (λ/2 + 1/λ) ,
```

acting on V_j ⊗ ℂ² for each j (with the left factor contributing
multiplicity 2j+1); σ± = σ₁ ± iσ₂, J the spin-j angular momentum on the
right factor. The overall sign/constant conventions are pinned by the
round checkpoint below.

## 2. The spectrum in closed form

𝒟_λ commutes with the **total fiber momentum** m′ = m + s₃ (the good
sector label; k = 2m′ is the winding of #193 carried by the full spinor
mode). Each m′-sector of V_j ⊗ ℂ² is at most 2-dimensional:

**Family A** — the extreme states m = ±j, spin aligned (m′ = ±(j+½)),
on which the ladder terms annihilate:

```
a_j(λ) = (2j+1)/λ + λ/2 ,      multiplicity 2(2j+1) ,   j = 0, ½, 1, …
```

**Family B** — the interior 2×2 blocks, basis {|m′−½⟩⊗↑, |m′+½⟩⊗↓},

```
block = (λ/2 + 1/λ) − (1/λ) + (2m′/λ)·τ₃ + 2√((j+½)² − m′²)·τ₁
```

(τ = Pauli matrices on the block), giving

```
b±(j, m′; λ) = λ/2 ± 2√((j+½)² + m′²(λ⁻² − 1)) ,
multiplicity (2j+1),   j ≥ ½,   m′ ∈ {−(j−½), …, j−½} .
```

**Verification (machine-checked in the companion probe):**

1. *Algebra:* assembling 𝒟_λ as an explicit matrix on V_j ⊗ ℂ² for all
   j ≤ 3 at λ ∈ {0.6, 1, 2.7, 5}: agreement with the closed forms to
   ~10⁻¹⁵.
2. *Round checkpoint:* at λ = 1 the two families assemble to exactly
   ±(3/2 + n) with multiplicities (n+1)(n+2) — family A supplies
   2(n+1) of the +(3/2+n) level, family B's plus branch n(n+1), the
   minus branch all the negatives. Verified through n = 7, exact
   counts.
3. *Collapse limit* (λ → 0, Gromov–Hausdorff collapse to the base
   S²(½)): all m′ ≠ 0 modes diverge as 2|m′|/λ (the KK fiber momenta),
   and the m′ = 0 survivors converge to ±2(l+1) — **precisely the Dirac
   spectrum of the round S² of radius ½**, the Ammann–Bär collapsing
   circle-bundle picture. Verified.
4. *Lichnerowicz consistency:* 𝒟² ≥ scal/4 forbids zero modes while
   scal(λ) = 8 − 2λ² > 0, i.e. for all λ < 2. Every zero of the closed
   forms indeed lies at λ > 2 (the first at λ = 4). The formulas had no
   right to pass this by accident.
5. *Spectral asymmetry:* for λ ≠ 1 the spectrum is asymmetric (family A
   is entirely positive; B± is centered at λ/2) — the nonzero eta
   invariant of Berger spheres, as it must be.

## 3. The odd-k ladder off the round metric

The winding-k sector (k = 2m′ odd; on the RP³ quotient these are the
antipodally-odd, Pin-twisted modes — the deck map is a fiber
translation, an isometry of every Berger metric, so the grading is
λ-independent exactly as in #193). The sector's lowest |eigenvalue|:

```
m_k(λ) = k/λ + λ/2                                   for  λ ≤ λ×(k) ,
       = 2√((k+2)²/4 + (k²/4)(λ⁻²−1)) − λ/2          for  λ×(k) ≤ λ ≤ λ*(k) ,
```

- **λ×(k) = √(2k+4)**: setting k/λ + λ/2 = |b₋((k+1)/2, k/2)| gives
  λ² = 2k+4 exactly — the character change. Below it (including the
  whole squash side and a wide stretch margin) the sector ground is the
  pure **winding state** with the Kaluza–Klein linear momentum k/λ (the
  first-order counterpart of the scalar (k/λ)² of #193 — the #83 throat
  winding term, again derived).
- **λ*(k)² = 8(k+1) + 2√(16(k+1)² + k²)**: the harmonic-spinor point
  where m_k → 0. Numbers: λ*(1) ≈ 5.668, λ*(3) ≈ 8.035,
  λ*(5) ≈ 9.851 — under extreme stretch **the electron sector is the
  first to become massless**, with μ and τ following in order. (The
  very first harmonic spinor on the family sits in an even sector at
  λ = 4(j+½) = 4 — Hitchin's phenomenon, located.)

**Persistence and crossings, answered exactly:**

| question | answer |
|---|---|
| which levels persist? | the family-A winding tower persists on all of (0,∞): a_j > 0, no crossings within A |
| where do crossings occur? | sector-minimum character changes at λ×(k) = √6, √10, √14; zero crossings (harmonic spinors) at λ = 4 (even sector) and λ*(k) ≈ 5.67, 8.03, 9.85 (odd sectors); none for λ < 2 (Lichnerowicz) |
| infinitesimal squash? | analytically smooth; sensitivities (½−k)/(k+½) = −0.33, −0.71, −0.82 — the true spinor ladder is metric-soft (the #192 surrogate's Δ ≈ 71 has no counterpart here, confirming #194/#195: the fine-tuning was dynamics, not spectral geometry) |
| ordering window? | m₁ < m₃ < m₅ on the whole (0, λ*₁ ≈ 5.668); all three grounds are pure winding states on (0, √6 ≈ 2.449) |
| squash side? | absolutely protected: m_k → k/λ, gaps 2/λ grow |

## 4. The k₅ = 5 question

Does the three-generation cutoff have a spectral counterpart at λ ≠ 1?
**No.** The family-A gap is

```
m_{k+2}(λ) − m_k(λ) = 2/λ    — independent of k, for every λ,
```

so nothing in the spectrum distinguishes k = 5 from k = 7 at any squash:
the tower continues 7, 9, … with identical spacing on the whole Berger
family. The k ≤ k₅ = 5 cutoff is **dynamical** — the phase budget
Φ_avail(k) of #122/#136 — and this is now a closed-form statement for
every λ, not a round-point observation. (Even the collapse boundaries
λ*(k) grow monotonically with k — extreme stretch removes sectors from
the bottom, never truncates from the top.)

## 5. Relation to the arc, and honest scope

- **Upgrades:** #183 (algebraic protection) → #192 (surrogate window,
  numeric) → #193 (scalar operator, closed form) → **this: the actual
  spinor field content, closed form, whole family, all crossings
  located.** The flagship generation claim is now spectral fact.
- **What stays dynamical:** the round-point ratios are 7/3 and 11/7 —
  O(1), as for the scalar. The spectrum supplies grading, ordering, and
  protection; the mass hierarchy remains the instanton dynamics
  (#193–#195). Nothing here changes that division of labor.
- **Scope:** the analysis is on the S³ cover with its unique spin
  structure; on the RP³ quotient the odd-k sectors are the Pin-twisted
  modes and the statements descend unchanged (the deck map is a fiber
  translation — an isometry of every Berger metric). The Berger family
  is the one-parameter fiber/base-separating deformation axis; λ×, λ*
  are boundaries within that family. The spectrum itself is classical
  (Hitchin; Bär); the self-contained derivation, the checkpoint chain,
  and the application to the BAM ladder questions are what this PR adds.

## References

- N. Hitchin, *Harmonic spinors*, Adv. Math. 14 (1974) 1–55. [Dirac
  spectrum on Berger spheres; the first harmonic spinors from squashing.]
- C. Bär, *Metrics with harmonic spinors*, GAFA 6 (1996) 899–942. [Zero
  modes on Berger spheres; surgery constructions.]
- B. Ammann, C. Bär, *The Dirac operator on nilmanifolds and collapsing
  circle bundles*, Ann. Glob. Anal. Geom. 16 (1998) 221–253. [The λ→0
  collapse limit checked in §2.]
- T. Friedrich, *Der erste Eigenwert des Dirac-Operators…*, Math. Nachr.
  97 (1980) 117. [The round-point first-eigenvalue bound 3/2.]
- A. Lichnerowicz, *Spineurs harmoniques*, C. R. Acad. Sci. Paris 257
  (1963) 7–9. [𝒟² ≥ scal/4: the λ < 2 exclusion of zeros.]
- J. S. Dowker, *Fermions on the squashed 3-sphere* (and related).
  [Physics treatments of the same spectrum.]

## Reproduce (the machine-checkable identities)

```bash
python -m experiments.closure_ledger.berger_dirac_analytic_ladder_probe
# Verdict: ANALYTIC_BERGER_DIRAC_ODD_K_LADDER_PROTECTED_WITH_ALL_CROSSINGS
#          _LOCATED_IN_CLOSED_FORM_NO_SPECTRAL_K5_CUTOFF_AT_ANY_LAMBDA
```
