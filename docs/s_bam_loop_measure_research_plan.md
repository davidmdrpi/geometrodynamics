# First-principles `S_BAM` loop measure: `1/(2π) = ` BAM closure quantum

Attempts the hardest standing piece of the spin-sector arc — the
**origin of `1/(2π)` in the Schwinger anomaly `a = α/(2π)` from a
`S_BAM` loop measure** (the open prize flagged by PR #62, which
*reconstructed* Schwinger using tree-normalized BAM primitives without
deriving the `1/(2π)` measure factor).

## Honest target

A full rigorous derivation of the `(2π)^d` covariant Fourier measure
from the `S_BAM` path integral on the throat configuration space is
genuine open work (the user explicitly flagged this as "likely an honest
structural sketch + identified gap").

This probe delivers the **structural identification** that the `1/(2π)`
in `a = α/(2π)` is the **same `2π`** that underlies BAM's closure
quantum (`action_base = 2π`, the S³ great-circle quantum). The QFT loop
measure `d^d k/(2π)^d` and the BAM closure ledger's `Φ_avail(k) = 2π(k+1)
+ …` share the same primitive — the Fourier-conjugate quantum of a
closed cycle of length `2π`. Each loop in a BAM diagram contributes one
factor of `1/(2π)` per momentum dimension from this closure-cycle
measure.

## The structural correspondence

**Closed cycle ↔ Fourier measure (the same `2π`):**

Take a closed cycle of length `L`. Momenta are quantized in units of
`2π/L`, with density of states `L/(2π)` per momentum interval. Sums
over modes become integrals with the measure `dk/(2π)·L`:

```
Σ_n  f(k_n)   →   (L / 2π)  ∫ dk  f(k)        (one closed cycle, length L)
```

For `L = 2π` (the BAM S³ great-circle quantum, `action_base = 2π` in
natural units), the density of states is `2π/(2π) = 1`, and the loop
integration becomes `∫dk/(2π)`. **Each momentum dimension contributes
one factor of `1/(2π)` from the closure-cycle Fourier measure** — the
same `2π` as BAM's closure quantum.

## Where the `2π` appears (across all BAM sectors)

| Sector | `2π` appearance | source |
|---|---|---|
| Closure ledger | `Φ_avail(k) = 2π(k+1) + 50π·max(0,k−3)²` | `docs/odd_k_closure_lemma.md` |
| Action base | `action_base = 2π` (S³ great circle) | `docs/hbar_origin_status.md` |
| Hopf holonomy | `∮A = π cos χ` (half-cycle); 2π at full | `geometrodynamics/hopf/connection.py` |
| β_lepton | `β = k_5²·(2π) = 50π` | PR #71 |
| Throat dwell | `τ(ω) = π/ω` (half-cycle); 2π at full | PR #51, K-channel |
| ε / closure-quantum integers | `4β/(2π) = 100 = 4·k_5²` | PRs #71, hbar_origin |
| **Schwinger loop measure** | **`a = α/(2π)`** | this PR (identification with the above) |
| QFT Fourier measure | `d^d k/(2π)^d` per loop | standard QFT |

The Schwinger `1/(2π)` is structurally **the same `2π`** as all the
other BAM appearances — the Fourier-conjugate quantum of the closed
S³ great circle, the foundational BAM closure quantum.

## What this identifies (and what it doesn't)

  - **Identifies:** the `1/(2π)` in `a = α/(2π)` as the BAM
    closure-quantum loop measure factor; the same `2π` primitive that
    underlies the closure ledger, the action base, β_lepton, the Hopf
    holonomy, and the throat dwell. One closure quantum per loop momentum
    dimension. Each `n`-loop correction picks up `(α/(2π))^n` ×
    dimensionless coefficient, consistent with the structural reading.

  - **Does not (and this probe does not claim to):** rigorously derive
    the full `(2π)^d` covariant Fourier measure from a written-out
    `S_BAM` path integral on the throat configuration space. That is
    genuine open work — the explicit covariant loop measure requires
    constructing the BAM path integral, gauge fixing, the appropriate
    Jacobians, etc., which would be a substantial undertaking and is
    beyond this probe's scope.

  - **What this DOES advance over PR #62:** #62 reconstructed `a = α/(2π)`
    using tree-normalized BAM primitives (S³ Green function as virtual
    photon, throat-pinch vertex). The `1/(2π)` there was *inherited from
    the tree normalization* — silent. This PR identifies it as the BAM
    closure quantum, the same primitive that appears throughout the BAM
    structural arc, giving the loop measure an explicit BAM-native
    origin.

## B4 accounting

`2π` is **dimensionless** (radians / phase / closure quantum). The
identification is structural/topological — the same Fourier-conjugate
quantum across all BAM sectors. Independent of the dimensionful anchor.

## Tests

  T1. **`a = α/(2π) ≈ 0.0011614`.** the Schwinger value (recap #62);
      vs measured `a_e ≈ 0.00115965`.
  T2. **The closed-cycle ↔ Fourier-measure correspondence.** length-L
      closed cycle ⟹ density of states `L/(2π)`; for `L = 2π` (the
      S³ great circle) the loop integration measure is `dk/(2π)`.
  T3. **BAM's `action_base = 2π` is the closure quantum.** the same
      `2π` as the QFT Fourier measure denominator — Fourier-conjugate
      quantum of the S³ great circle.
  T4. **The same `2π` across all BAM sectors.** closure ledger
      (`Φ_avail = 2π(k+1)+…`), action_base, β_lepton = k_5²·(2π),
      Hopf holonomy, throat dwell, ε's `4β/(2π) = 100`, Schwinger
      `1/(2π)`. One primitive.
  T5. **One loop ↔ one closure quantum (the `n`-loop expansion).**
      `(α/(2π))^n` per n-loop contribution. Check the leading `α/(2π)`
      coefficient + the standard `(α/π)²` two-loop coefficient
      structure.
  T6. **Honest scope.** structural identification of `1/(2π)` = BAM
      closure quantum (same primitive across all sectors); does NOT
      rigorously derive the full `(2π)^d` covariant Fourier measure
      from `S_BAM` (genuine open work).
  T7. **Falsification / B4.** if BAM's closure quantum were not `2π`,
      or if Schwinger's `1/(2π)` had a different origin, the
      identification would fail. BAM passes (same primitive). `2π`
      dimensionless; structural.
  T8. **Assessment.**

## Verdict structure

  - **LOOP_MEASURE_IDENTIFIED** (expected): the `1/(2π)` in
    `a = α/(2π)` is identified as the **BAM closure-quantum loop
    measure factor** — the same `2π` that underlies BAM's `action_base`,
    closure ledger, β_lepton, Hopf holonomy, throat dwell, ε integer,
    and every other BAM "2π" appearance. One closed cycle of length
    `2π` (the S³ great circle) gives a Fourier-conjugate measure of
    `1/(2π)` per loop momentum dimension — the standard QFT loop
    measure and the BAM closure quantum are the same primitive. This
    closes the structural piece of the open follow-on (origin of
    `1/(2π)`); a fully rigorous covariant `S_BAM` path-integral
    derivation of the `(2π)^d` measure remains future work.

  - **NO_IDENTIFICATION**: the `2π` factors do not match across the
    sectors (structurally).

## What this leaves open

  - **A full covariant `S_BAM` path-integral derivation** of the
    `(2π)^d` Fourier measure on the throat configuration space.
    Requires writing out the BAM path integral, gauge fixing, the
    appropriate Jacobians — substantial work outside this probe's scope.
  - **Higher-loop coefficients of `a_e`** (the α², α³ corrections beyond
    the leading α/(2π)): the structural `(α/(2π))^n` scaling matches
    one closure quantum per loop, but the dimensionless multiplicative
    coefficients (1, 0.328…, 1.181…, ...) are separate calculations
    not addressed here.

## Cross-references

  - `docs/throat_vertex_loop_research_plan.md` — PR #62 (reconstructed
    `a = α/(2π)`; flagged the open `1/(2π)` derivation).
  - `docs/odd_k_closure_lemma.md` — `Φ_avail(k) = 2π(k+1) + …` (the
    closure ledger's `2π` appearance).
  - `docs/hbar_origin_status.md` — `action_base = 2π` (the BAM closure
    quantum).
  - `docs/beta_lepton_derivation_research_plan.md` — `β_lepton =
    k_5²·(2π) = 50π` (PR #71, the same `2π`).
  - `geometrodynamics/hopf/connection.py` — `hopf_holonomy(χ) = π cos χ`
    (half-cycle).
  - `experiments/closure_ledger/s_bam_loop_measure_probe.py` — this
    probe.
