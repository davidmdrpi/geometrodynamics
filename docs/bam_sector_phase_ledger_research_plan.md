# The BAM sector-phase ledger: continuous η-phases vs discrete Z₂ topology (PR #121)

PRs #117–#120 built and lattice-validated the `det'(∂_τ)` η-invariant
machinery. This PR converts it into a BAM **sector-phase ledger** that
separates the two — and only two — sources of phase/sign in the loop
measure, and proves they are **not double-counted**:

1. **continuous η-phases** from the U(1) holonomy (the Hopf/Wilson line);
2. **discrete Z₂ signs** from the Möbius / odd-k orientation (the
   non-orientable closure).

## The two structures (independent geometric data)

A BAM closure loop carries two independent twists:

  - a **U(1) holonomy** `a ∈ [0,1)` — the connection / Wilson line `∮A =
    e^{2πia}` (the Hopf holonomy, `a = kχ/2π`). The bundle's **connection**
    (`π₁` of the fibre), continuous.
  - an **orientation class** — orientable vs non-orientable (Möbius),
    captured by the winding parity `k mod 2` (the odd-k closure lemma) `=
    w₁ ∈ H¹(loop; Z₂)`, the bundle's **orientability**. Discrete.

## (1) Continuous η-phase

The twisted determinant (PR #119/#120) is `det P_a = 2 sin(πa)·e^{iθ(a)}`
with

```
θ(a) = (π/2)·η_A(0) = (π/2)(1 − 2a),      a ∈ (0,1).
```

Because the twisted operator has no zero mode, `ζ(0) = 0`, so the phase is
**purely** the η-invariant piece. As `a` sweeps `(0,1)`, `θ(a)` sweeps
`(−π/2, +π/2)`: the η-phase is confined to the **open right half-circle**
(`Re > 0`), reaching `+1` only at `a = 1/2` and **never `−1`**.

| a | θ = (π/2)(1−2a) | e^{iθ} | Re |
|---:|---:|---|---:|
| →0 | +π/2 | +0.000+1.000i | 0 |
| 1/4 | +π/4 | +0.707+0.707i | +0.707 |
| 1/3 | +π/6 | +0.866+0.500i | +0.866 |
| 1/2 | 0 | +1.000 | +1.000 |
| 2/3 | −π/6 | +0.866−0.500i | +0.866 |
| 3/4 | −π/4 | +0.707−0.707i | +0.707 |

## (2) Discrete Z₂ sign

The Möbius / non-orientable closure contributes the orientation sign of the
odd-k closure lemma (PR #115/#118):

```
e^{ikπ} = (−1)^k :  +1 (even k, orientable / torus cover),
                    −1 (odd k, non-orientable / Möbius half-twist).
```

A discrete element of `{±1}` — the `w₁` holonomy — independent of `a`.

## No double-counting (the proof)

Three independent reasons:

  - **(a) Different groups.** The η-phase is `U(1)`-valued and continuous in
    `a`; the Z₂ sign is a discrete element of `{±1}`.
  - **(b) Different geometry.** The η-phase comes from the **connection**
    (holonomy `a`, spectral asymmetry); the Z₂ sign from the **orientation**
    (`w₁`, orientability) — distinct invariants.
  - **(c) No collision on the nontrivial element.** `θ(a) ∈ (−π/2, +π/2)`
    for `a ∈ (0,1)`, so the η-phase lies in the open right half-circle and
    can **never** equal `−1` (at `θ = ±π`). The Möbius `−1` is therefore
    inaccessible to the continuous η-phase — it is purely the discrete Z₂.
    At `a = 1/2` the η-phase is exactly `+1` (`det P_{1/2} = 2`, real), so
    all of the Möbius character there is carried by the separate `(−1)^k`.

Numerically (dense `a`-grid): the η-phase's minimum real part is `> 0`, and
its closest approach to `−1` is `≈ √2` (far) — it never collides with the
Z₂ sign.

## The sector-phase ledger

| contribution | source | object | group | value | type |
|---|---|---|---|---|---|
| magnitude | twisted spectrum | `\|det P_a\| = 2 sin(πa)` | ℝ₊ | — | continuous |
| **η-phase** | U(1) holonomy `a` (Hopf/Wilson) | spectral asymmetry `η(0)=1−2a` | U(1) | `e^{i(π/2)(1−2a)}` | **CONTINUOUS** |
| **Z₂ sign** | Möbius / odd-k orientation | `w₁ / e^{ikπ}` | Z₂ | `(−1)^k` | **DISCRETE** |
| (scaling ζ(0)) | local heat-kernel | `ζ(0)` | — | absorbed (`ζ(0)=0` twisted) | removed |

## The factorized measure phase

```
det_full = |det P_a| · e^{i(π/2)(1−2a)} · (−1)^k
         = [magnitude]  · [continuous η-phase] · [discrete Z₂ sign],
```

each factor counted exactly once. At `a = 1/2`, `k = 1`: `det_full = 2·(+1)·
(−1) = −2` (real) — the antiperiodic det's Möbius sign is purely the
`(−1)^k`.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | convert η-machinery → sector-phase ledger |
| T2 | two structures | U(1) holonomy (connection) vs Z₂ (`w₁`) |
| T3 | continuous η-phase | `θ(a) ∈ (−π/2,+π/2)`, right half-circle, never −1 |
| T4 | discrete Z₂ | `(−1)^k`: +1 torus, −1 Möbius |
| T5 | ledger | the sector-phase ledger table |
| T6 | no double-count | groups, geometry, no collision on −1 |
| T7 | factorization | `det_full = |det|·η-phase·Z₂`, each once |
| T8 | assessment | `BAM_SECTOR_PHASE_LEDGER_CONTINUOUS_ETA_TIMES_DISCRETE_Z2_NO_DOUBLE_COUNT` |

## Established and open

  - **Established (BAM-native):** the BAM loop-measure phase factorizes as a
    continuous η-phase `e^{i(π/2)(1−2a)}` (U(1) holonomy, open right
    half-circle, never −1) times a discrete Z₂ sign `(−1)^k` (Möbius/odd-k
    orientation, `w₁`); independent (different groups, geometry, no collision
    on −1) ⟹ no double-counting; the measure phase is the product, each
    factor once.

  - **Open (unchanged):** the analytic open pieces of the measure arc
    (absolute `Z` normalization `κ₅²/Λ₅`, multi-loop). This ledger organizes
    the phase structure, it does not close those.

## Cross-references

  - `docs/detprime_dtau_eta_invariant_phase_research_plan.md` — PR #119, the
    η-phase framework.
  - `docs/lattice_validation_research_plan.md` — PR #120, the lattice
    validation (incl. generic holonomy).
  - `docs/diff_s1_first_order_ghost_audit_research_plan.md` /
    `docs/odd_k_closure_lemma.md` — PR #118/the odd-k Z₂ orientation sign.

## Run

```
python -m experiments.closure_ledger.bam_sector_phase_ledger_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_bam_sector_phase_ledger_probe/`.
Expected verdict:
`BAM_SECTOR_PHASE_LEDGER_CONTINUOUS_ETA_TIMES_DISCRETE_Z2_NO_DOUBLE_COUNT`, 8/8 PASS.
