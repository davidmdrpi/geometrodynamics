# BAM throat-vertex loop probe — the Schwinger anomalous moment

The capstone of the spin sector. PR #61 derived `g = 2` at tree level
from the throat's Pauli/SU(2) + Hopf-monopole structure and recorded the
Schwinger anomaly `a = (g−2)/2 = α/2π` as the known one-loop correction.
This probe constructs the **one-loop throat-vertex correction** in
BAM-native terms — the throat emits and reabsorbs one virtual photon (the
S³ Green-function propagator, PRs #45–#46) at the throat-pinch vertex
(PRs #38–#41) — and shows it reproduces the Schwinger value.

## Honest scope (stated upfront)

This is a one-loop calculation, harder than the tree probes. The honest
result is a **reconstruction**, not an independent first-principles
derivation of `1/2π`:

  - The BAM tree thread (PRs #35–#46) normalized its primitives — the
    S³ Green-function photon propagator (`1/q²`), the throat-pinch vertex,
    the Hopf coupling `α` — to QED at tree level.
  - The one-loop **integrand** is therefore the QED vertex-correction
    integrand expressed in those primitives; the probe verifies the loop
    **integral** gives the Schwinger coefficient.
  - What is BAM-native: the loop's *pieces* (virtual photon = S³
    exchange; vertex = throat pinch) and the *structure* (the throat
    dressing its own moment by one self-exchange).
  - What is **not** independently derived: the `1/2π` from a pure-geometry
    loop measure derived from `S_BAM` — that needs the full covariant
    throat loop, the standing follow-on.

So the verdict is "the BAM throat-vertex loop *reproduces* Schwinger,"
with the coefficient inherited from the tree-normalized primitives.

## The loop

The QED vertex-correction triangle, in BAM terms:

  - **External throat lines** — the in/out electron (throat) at the
    pinch vertex.
  - **Virtual photon** — an S³ Green-function exchange. The S³ scalar
    Green function `G(ψ) = ((π−ψ)cot ψ − ½)/(4π²R)` has the flat limit
    `G → 1/(4π d)` (Coulomb), i.e. the `1/q²` photon propagator
    (PRs #45–#46).
  - **Vertex** — the throat-pinch `F²` vertex (PRs #38–#41), the Hopf
    coupling `α`.

The loop expansion parameter is `α/π` (coupling `α`, loop phase space
`1/π`); the anomalous moment is the `F₂(0)` form factor.

## The Schwinger coefficient

After the loop integral, the anomalous form factor is

```
F₂(0) = (α/2π) · I ,
   I = ∫dx dy dz δ(x+y+z−1) · 2m²z(1−z) / [m²(1−z)²]
     = ∫₀¹ 2z dz = 1 ,
```

so `a = F₂(0) = α/2π ≈ 0.0011614`, matching the measured
`a_e = 0.00115965` (the leading term; higher orders `α², α³, …` make up
the small remainder). The Feynman-parameter integral `I = 1` is verified
numerically.

## B4 accounting

`a = α/2π` is **dimensionless**; `α` is the EM coupling (an input,
related to the Hopf-charge structure), and the absolute mass scale is the
single anchor `m` (`m_e c² = ℏc/R_MID`). `g = 2(1 + a)` at one loop. The
anomaly is a pure number; the scale is the one anchor — consistent with
the whole arc.

## Tests

  T1. **BAM one-loop vertex structure.** Throat emits/reabsorbs one
      virtual photon (S³ propagator) at the throat-pinch vertex — the
      BAM image of the QED vertex-correction triangle; loop parameter
      `α/π`.
  T2. **Virtual photon = S³ exchange.** `G(ψ) → 1/(4π d)` (Coulomb /
      `1/q²`) in the flat limit (PRs #45–#46).
  T3. **Schwinger integral.** `I = ∫dxdydz δ(x+y+z−1) 2z(1−z)/(1−z)² =
      ∫₀¹2z dz = 1` ⟹ `F₂(0) = α/2π`.
  T4. **a = α/2π vs experiment.** `a = 0.0011614` vs `a_e = 0.00115965`.
  T5. **One-loop g.** `g = 2(1 + a) = 2.00232…`.
  T6. **Honest scope.** Reconstruction (tree-normalized primitives), not
      an independent `1/2π` from pure geometry; the full `S_BAM` loop
      measure is the open piece.
  T7. **B4 accounting.** `a`, `g` dimensionless; `α` the coupling; scale
      is the single anchor.
  T8. **Assessment.**

## Verdict structure

  - **SCHWINGER_RECONSTRUCTED** (expected): the BAM throat-vertex loop —
    the throat dressing its moment by one S³-photon self-exchange at the
    throat-pinch vertex — reproduces the Schwinger anomaly
    `a = α/2π` (the Feynman-parameter integral `I = 1`), matching
    `a_e` at the leading order. The loop's pieces and structure are
    BAM-native; the coefficient is inherited from the tree-normalized
    primitives. This is a reconstruction, not an independent
    first-principles derivation of `1/2π` (the full `S_BAM` loop measure
    is the open follow-on).

  - **LOOP_FAILS**: the integral does not give `α/2π`, or the loop
    pieces do not map to the BAM primitives.

## What this leaves open

  - **`1/2π` from a first-principles `S_BAM` loop measure.** The covariant
    throat loop derived from the action, rather than the QED-normalized
    integrand — the genuine derivation.
  - **Higher-order `a_e`.** The `α²`, `α³`, … terms (two- and three-loop).
  - **The explicit throat spinor / vertex from `S_BAM`** (shared with
    #59–#61).

## Cross-references

  - `docs/gyromagnetic_ratio_research_plan.md` — tree `g = 2` (#61).
  - `docs/tree_qed_status.md` — the tree-QED primitives (#35–#46).
  - `docs/bam_exchange_kernel_research_plan.md`,
    `docs/hopf_vector_exchange_kernel_research_plan.md` — the photon
    propagator from the S³ Green function (#45–#46).
  - `geometrodynamics/transaction/s3_geometry.py` — `s3_green_potential`.
  - `experiments/closure_ledger/throat_vertex_loop_probe.py` — this probe.
