# B4 dimensional-bridge audit + Maslov closure-ledger probe

Targets the **last surviving scaffold barrier, B4** — the dimensional
bridge `ℏ = m_e·R_MID·c`, the single external `m_e` anchor — and
formalizes the **Maslov index** of the closure-ledger Bohr–Sommerfeld
quantization. The two are connected: the Maslov machinery is what makes
the ledger *scale-free*, and scale-freeness is exactly *why* B4 is
irreducible.

## Part A — the Maslov closure-ledger

The closure ledger quantizes radial modes by Bohr–Sommerfeld with a
Maslov correction at each boundary of the classically allowed region:

```
∮ p dr*  =  2π·(n + μ/4)        μ = Maslov index = Σ_boundaries (phase/(π/2))
```

with the textbook per-boundary phase deficits

| boundary | reflection phase | Maslov β (= phase / (π/2)·¼·4) | μ contribution |
|---|---|---|---|
| soft turning point (Airy) | π/2 | 1/4 | 1 |
| hard wall (Dirichlet) | π | 1/2 | 2 |

The Tangherlini radial cavity `r ∈ [R_MID, R_OUTER]` has a **hard wall
at each end**: the throat (Dirichlet, B3, forced by `T² = −I`) and the
outer wall. So `μ = 4`, and

```
∮ p dr* / 2π  =  n + 4/4  =  n + 1
```

— exactly the closure-ledger Layer-2 integer `S_full(l,n)/2π → (n+1)`.
**The "+1" per radial mode is the Maslov index of the doubly-Dirichlet
throat cavity.** Its throat half (`μ = 2`, reflection phase `π`) is the
B3 hard wall; and that `π` is simultaneously

  - the **closure half-quantum** (the `2π` closure split by the
    antipodal `Z₂`), and
  - the **throat dwell phase** `ω·τ(ω) = π` with `τ(ω) = π/ω` (the
    K-channel impedance, PR #51).

Three readings of one `π`: Dirichlet reflection = closure half-quantum
= dwell phase. The per-species ledger integers `(3, 6, 109)` are then
sums of dimensionless quanta (antipodal `k`, Hopf+throat `1`, β-uplift,
radial Maslov `n+1`).

## Part B — the B4 dimensional-bridge audit

**Claim audited:** can B4 be derived away, i.e. can the absolute MeV
scale (or `ℏ` in SI) be produced from BAM geometry without an external
anchor?

**Finding: no — and provably so.** Every quantity the closure-ledger +
Maslov machinery produces is **dimensionless**:

  - winding integers `k, n`, Maslov index `μ`;
  - action ratios `S/2π`; the closure quantum `2π`;
  - eigenfrequencies `ω(l,n)` in geometric units (`R_MID = 1`);
  - mass **ratios** `m_μ/m_e`, `m_τ/m_e`;
  - the geometric residuals `R_OUTER`, `ε`, transport, resistance — all
    closure-quantum invariants (ℏ-origin thread);
  - the `1.054` factor.

A theory whose every output is dimensionless **cannot** yield a
dimensionful scale (dimensional analysis). Exactly **one** dimensionful
anchor is mathematically required, and `m_e` (equivalently `R_MID` via
`ℏ = m_e·R_MID·c`) is the minimal choice. This is demonstrated
concretely by **scale invariance**: rescaling `R_MID → λ·R_MID` leaves
every dimensionless ledger output unchanged (`ω·R_MID`, `S/2π`, `μ`,
mass ratios) and shifts only the absolute scale. B4 is therefore
**irreducible by dimensional necessity** — not an unsolved gap but a
structural feature, the single mandatory unit any scale-free theory
must import (cf. SI fixing `c, ℏ, e` by convention).

## Tests

  T1. **Maslov dictionary** — `{soft TP: ¼, hard wall: ½}` reproduces
      the eigenvalues of the analytic cavities: the infinite square
      well (hard+hard → `∫k dx = π(n+1)`) and the harmonic oscillator
      (soft+soft → `∫k dx = π(n+½)`).
  T2. **Tangherlini cavity is doubly-Dirichlet (μ=4)** — the one-sided
      WKB action `∫√(ω²−V) dr* → π(n+1)` (so `∮/2π → n+1`); identify
      the `+1` as the `μ=4` Maslov count. (Reuses `sk_bridge`.)
  T3. **Throat π is three things at once** — Dirichlet reflection phase
      `π` (B3) = closure half-quantum `π` = dwell phase `ω·τ = π`
      (`τ = π/ω`).
  T4. **Ledger integers are dimensionless quanta** — `(3, 6, 109)` from
      `k + 1 + uplift + (n+1)`; each an integer.
  T5. **Scale invariance (core audit)** — `R_MID → λ·R_MID` leaves
      `ω·R_MID`, `S/2π`, `μ`, and mass ratios invariant to machine
      precision; only the absolute scale changes. The machinery is
      scale-free.
  T6. **Dimensionless residuals closed** — `R_OUTER ≈ 1.262`,
      `ε = 7π/(100·5⁴)`, transport `= 8π`, resistance `= 7π/100`; all
      dimensionless closure-quantum invariants (ℏ-origin thread).
  T7. **1.054 does not close B4** — recap the clean negative (no closed
      form within 0.01 %); reframe: `1.054` is dimensionless, so even a
      closed form would not supply the dimensionful scale.
  T8. **B4 irreducibility** — exactly one dimensionful input (`m_e`);
      `ℏ = m_e·R_MID·c` is dimensionally consistent; all else derived
      and dimensionless.
  T9. **Scaffold final assessment** — B1, B2, B3, B5 closed; B4
      irreducible-by-necessity. The programme is complete in the precise
      sense that the only remaining input is the single mandatory
      dimensionful unit.

## Verdict structure

  - **B4_IRREDUCIBLE** (expected): the Maslov closure-ledger is
    formalized (the radial `+1` is the `μ=4` Maslov index; the throat
    `π` is the B3 reflection = closure half-quantum = dwell phase), and
    the B4 audit shows the machinery is scale-free, so exactly one
    dimensionful anchor (`m_e`) is mathematically required. B4 is
    irreducible by dimensional necessity, not an open gap. The scaffold
    is closed: four barriers derived, the fifth shown mandatory.

  - **B4_REDUCIBLE** (not expected): some dimensionful output emerges
    from the dimensionless machinery — would contradict scale
    invariance.

  - **AUDIT_INCONSISTENT**: a Maslov count or the scale-invariance
    identity fails numerically.

## What this leaves open

  - **`m_e` itself.** The single anchor is not derived — by the audit,
    it *cannot* be from scale-free geometry. Deriving `m_e` would
    require new dimensionful physics outside the closure-ledger
    scaffold (e.g. a coupling to a fixed external length).
  - **First-principles internal action.** The Maslov indices are read
    from the WKB/boundary structure; an explicit covariant action whose
    stationary phase reproduces them is the standing follow-on (shared
    with the master-integral residual).

## Cross-references

  - `docs/bam_scaffold_status.md` — barrier ledger (B1–B5); B4 the
    surviving barrier.
  - `docs/hbar_origin_status.md` — the ℏ-origin thread: dimensionless
    residuals (`R_OUTER`, `ε`, transport, resistance) closed; `1.054`
    clean negative.
  - `docs/master_integral_research_plan.md` — B5 closed (the master
    integral); B4 noted orthogonal/dimensional.
  - `experiments/closure_ledger/sk_bridge.py` — WKB radial action +
    Maslov policy (`−π/2` per soft turning point).
  - `experiments/closure_ledger/closed_orbit_radial_action_probe.py` —
    `S_full(l,n)/2π → (n+1)` (hard-wall BS).
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`,
    tortoise coordinate.
  - `experiments/closure_ledger/maslov_dimensional_bridge_probe.py` —
    this probe.
