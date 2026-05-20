# Covariant BAM effective-action scaffold — research plan

The deepest open problem of the BAM programme, flagged since PR #35
and restated in `docs/tree_qed_status.md` ("a first-principles BAM
*Lagrangian* for QED ... is not claimed") and `docs/THESIS.md`: is
there a **single covariant variational principle** from which the
three distinct BAM structures derive?

  1. the Compton vertex factor `F²(x, c)` (QFT amplitude level);
  2. the Hopf-bundle U(1) connection `A_φ(χ) = ½cos χ` (gauge geometry);
  3. the 5D Tangherlini bulk boundaries `ΔR = R_OUTER − R_INNER`
     (bulk geometric scale).

This probe builds the **scaffold** — writes the candidate action,
performs the sector-by-sector variation, records which targets derive
cleanly, and **rigorously identifies the mismatch terms that prevent
full closure**. It is explicitly *not* a claim that the action closes;
it is a map of where it does and where it does not.

## The candidate covariant action

On the 5D Tangherlini bulk `M₅` with a non-orientable wormhole throat
and an `S³` spatial slice carrying the Hopf fibration:

```
S_BAM[g, A, ψ, Φ]
  = ∫_{M₅} d⁵x √(−g₅) [
        (1/2κ₅)(R₅ − 2Λ₅)              # (G) 5D gravity
        − ¼ F_{MN} F^{MN}              # (A) U(1) gauge
        + ψ̄ (i Γ^M D_M − m) ψ          # (D) Dirac fermions
        + L_throat[Φ, g, A]           # (T) throat / closure scalar
    ]
  + S_∂[hard walls]                   # (B) topological boundary conditions
  + S_closure[∮A = 2πn]               # (C) closure-quantum constraint
```

The three targets attach to three sectors:

| target | sector | variation |
|---|---|---|
| `A_φ = ½cos χ` | (A) U(1) gauge | `δS/δA = 0` on `S³` |
| `ΔR = R_OUTER − R_INNER` | (G) gravity + (B) BCs | `δS/δg = 0` + boundary data |
| `F²(x, c)` | (T) throat + (C) closure | `δS/δΦ = 0` + `∮A = 2πn` |

## Sector status (anticipated)

  - **Gauge (A) — closes.** The Hopf connection `A_φ = ½cos χ` is the
    homogeneous (symmetric-critical) solution of source-free Maxwell
    on the Hopf bundle over `S³`; its first Chern number is exactly
    `c₁ = 1` (verified to `< 1e-9` in `geometrodynamics.hopf.chern`).
    This sector derives from `δS/δA = 0` with no external input.

  - **Throat (T) — closes only with the closure constraint (C).** The
    vertex factor `F² = K(x)²·Q(x, c)` derives from the throat action
    (PR #41) given the closure quantum `S = 2π` and the antipodal Z₂.
    Both are *constraints*, not local Lagrangian densities.

  - **Gravity (G) — partial.** The Tangherlini metric
    `f(r) = 1 − (r_s/r)²` solves the 5D vacuum Einstein equation
    (the Hayward interior emerges from the throat density,
    `geometrodynamics.blackhole.derivation`), but the radial domain
    boundaries `R_INNER, R_OUTER` are fixed by boundary data, not by
    the bulk equation alone.

## The mismatch terms (barriers to full closure)

The scaffold exposes a structured set of barriers. The probe
enumerates and, where possible, quantifies each.

  - **B1 — closure quantum is topological.** The `S = 2π`
    (`action_base`) that selects `F²` is a Wilson-loop / holonomy
    constraint `∮A = 2πn`, not a local density. Local variation of
    the throat action leaves a one-parameter family (the closure
    constant); the `2π` selects one point. *Mismatch: one global
    topological input per closed orbit.*

  - **B2 — antipodal Z₂ is a discrete quotient.** The equal-action
    splitting (giving `K` and `Q`) requires the orientation-reversing
    isometry `σ: p → −p` (`T = iσ_y`, `T² = −I`). A smooth covariant
    action carries continuous symmetries; this discrete Z₂ must be
    imposed as an identification `S³` ↔ antipode (or a boundary
    condition). *Mismatch: one discrete symmetry not encoded in the
    smooth Lagrangian.*

  - **B3 — boundary conditions are topologically imposed.** The
    hard-wall (Dirichlet) condition at the throat comes from the
    spinor double cover `T² = −I` forcing `ψ = 0` at the throat fixed
    point — a topological consequence of the non-orientable throat,
    not derivable from `δS/δg = 0`. *Mismatch: boundary data imposed
    by topology.*

  - **B4 — dimensional bridge / scale incompleteness.** The action
    fixes ratios and dimensionless structure; the absolute mass scale
    requires the external anchor `ℏ = m_e · R_MID · c`. The action is
    *scale-complete modulo `m_e`* (consistent with
    `docs/hbar_origin_status.md`). *Mismatch: one external scale.*

  - **B5 — 5D → 4D reduction not constructed.** `F²` is a 4D-effective
    vertex factor; the bulk action is 5D. The Kaluza–Klein-like
    integration over the radial channel that produces `F²` from the
    5D action is not written. *Mismatch: the dimensional-reduction map
    is a structural gap, not a number.*

## Tests

  T1. **Action scaffold structure.** State the action and the
      target → sector → variation map (informative).

  T2. **Gauge sector closes.** `A_φ = ½cos χ` from
      `geometrodynamics.hopf.connection`; verify it is the
      homogeneous Maxwell solution and `c₁ = 1` to `< 1e-9`
      (`geometrodynamics.hopf.chern.compute_c1`). No external input.

  T3. **Throat sector closes with constraint (C).** Recompute
      `K(x) = 2x/(1+x)` and `Q(x, c) = x² + x(1−x)²/(1+c²)` from the
      equal-action splitting and verify `F² = K²·Q`. Flag the closure
      quantum (B1) and antipodal Z₂ (B2) as imposed.

  T4. **Gravity sector partial.** Verify `f(r) = 1 − (r_s/r)²` is the
      Tangherlini form on the radial grid; note that `R_OUTER` is
      fixed by the cross-species γ-lock fixed point and the throat
      `R_MID = 1` is a gauge choice — boundary data, not bulk-equation
      output (B3).

  T5. **Dimensional bridge barrier (B4).** Show the action is
      scale-free in the relevant sectors; the absolute scale enters
      only through `m_e` (`ℏ = m_e R_MID c`).

  T6. **5D → 4D reduction barrier (B5).** Document that the map from
      the 5D bulk to the 4D `F²` vertex is not constructed; identify
      what it would require (radial-mode integration + vertex
      projection).

  T7. **Barrier ledger / closure assessment.** Tabulate the five
      mismatch terms, their type (topological / discrete / boundary /
      scale / reduction), and which sectors close. Verdict.

## Verdict structure

  - **SCAFFOLD_WITH_BARRIERS** (expected): the action unifies the
    three targets structurally — the gauge sector closes outright, the
    throat sector closes given the closure constraint, the gravity
    sector is partial — but full closure is blocked by the five
    identified mismatch terms. The barriers are structured (two
    topological, one boundary, one scale, one reduction), which scopes
    the remaining work rather than leaving it open-ended.

  - **UNEXPECTED_FULL_CLOSURE**: if every target derived from `δS = 0`
    with no imposed constraint. Not anticipated; would be a major
    result.

  - **SCAFFOLD_FAILS**: if a sector that should close (gauge) does not
    reproduce its target — would indicate an error in the action
    ansatz.

## What this is and is not

**Is:** an honest map of a candidate unifying action — which pieces
derive, which are imposed, and the precise nature of each imposed
input.

**Is not:** a derivation of QED + gravity + gauge from one Lagrangian.
The barriers are real. The value is in their structure: they are not
five unrelated patches but a recognizable set (topological
quantisation, discrete antipodal identification, the dimensional
anchor) that recurs throughout the BAM programme.

## Cross-references

  - `geometrodynamics/hopf/connection.py`, `geometrodynamics/hopf/chern.py`
    — gauge sector (A).
  - `geometrodynamics/tangherlini/radial.py`,
    `geometrodynamics/blackhole/derivation.py` — gravity sector (G).
  - PR #41 `throat_action_derivation_probe` — throat sector (T).
  - PR #44 `mobius_exchange_sign_probe` — antipodal Z₂ (B2).
  - `docs/hbar_origin_status.md` — dimensional bridge (B4).
  - `docs/tree_qed_status.md` — the tree-QED result the F² target
    summarises.
  - `experiments/closure_ledger/bam_effective_action_scaffold_probe.py`
    — this probe.
