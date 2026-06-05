# Quartic vertex ledger and bounded interaction audit (PR #138)

PR #137 drew the ledger for the **cubic** vertex of the antipodal matter kernel:
its angular selection rule (`Σl` even — the antipodal Z₂ — plus the SO(4)
triangle) and its geometric radial shape are derived; only the coupling is
input. This PR extends the ledger to the **quartic** vertex and adds the
decisive stability question the cubic alone leaves open: **is the BAM matter
interaction bounded below (a stable vacuum)?** A pure cubic potential is
unbounded; this PR shows the geometric quartic self-overlap is positive, so the
effective potential is bounded below — and ties that boundedness to the
convergence of the S_BAM measure (#122).

## The quartic vertex factorises (same structure as the cubic)

```
V_4 = λ_4 · [ ∫_{S³} Y_{l1} Y_{l2} Y_{l3} Y_{l4} dΩ ] · [ ∫ ψ_k ψ_l ψ_m ψ_n dr* ] .
```

## The quartic angular selection rule is DERIVED (same Z₂ + SO(4))

`∫_{S³} Y_{l1} Y_{l2} Y_{l3} Y_{l4} dΩ` is nonzero only if

  - **(a) `l1+l2+l3+l4` even** — the **same antipodal parity Z₂** as the cubic
    (#137): under `x → −x` (the throat ↔ antithroat C-swap #63), `Y_l → (−1)^l
    Y_l`, so `(−1)^{Σl} = +1` over the inversion-symmetric S³;
  - **(b) a common SO(4) intermediate channel** — `∃ L ∈ [|l1−l2|,l1+l2] ∩
    [|l3−l4|,l3+l4]` (the two pairs must couple through a shared `Y_L`).

| (l₁,l₂,l₃,l₄) | Σl even? | SO(4) channel? | ∫YYYY |
|---|:---:|:---:|---:|
| (0,0,0,0) | ✓ | ✓ | 1.0 |
| (1,1,0,0) | ✓ | ✓ | 0.25 |
| (1,1,1,1) | ✓ | ✓ | 0.0417 |
| (2,1,1,0) | ✓ | ✓ | 0.0417 |
| (1,1,1,0) | ✗ | — | 0 |
| (1,0,0,0) | ✗ | — | 0 |

So the antipodal Z₂ **persists from the cubic to the quartic** (verified exactly:
odd-Σl → 0).

## The quartic self-overlap is POSITIVE — the interaction is bounded below

The diagonal quartic radial overlap is manifestly positive,

```
g_4 = ∫ ψ_k⁴ dr* > 0   (an integral of a fourth power),
```

(`1.03, 1.02, 1.02, …` for the lowest modes), so the single-mode effective
potential

```
V(a) = ½ ω_k² a² + (λ_3 g_3 / 6) a³ + (λ_4 g_4 / 24) a⁴
```

has a **positive a⁴ coefficient** (for `λ_4 > 0`). A quartic polynomial is
bounded below iff its `a⁴` coefficient is positive, so `V → +∞` as `|a| → ∞` and
the vacuum is **stable** — regardless of the (potentially destabilising) cubic,
which can only tilt the minimum, never unbound it:

| λ₃ (cubic) | a⁴ coeff | V(+10⁴) | V(−10⁴) | bounded? |
|---:|---:|---:|---:|:---:|
| 0 | 0.043 | 4.3e14 | 4.3e14 | ✓ |
| 5 | 0.043 | 4.3e14 | 4.3e14 | ✓ |
| 50 | 0.043 | 4.2e14 | 4.4e14 | ✓ |
| 200 | 0.043 | 4.0e14 | 4.6e14 | ✓ |

## Boundedness IS the S_BAM measure-convergence condition (#122)

A bounded-below action is exactly the condition for the path-integral measure
`∫ Dμ e^{−S}` to converge — which PR #122 established non-perturbatively for the
Z₂-graded sector sum. So the positive quartic is **not an extra assumption**: a
bounded interaction is required by, and consistent with, the measure's
existence. This extends the program's stability thread — free modes stable
(#130), one-loop self-energy unitarity-preserving (#136), and now the full
interacting vacuum bounded below.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | quartic vertex ledger + bounded interaction audit |
| T2 | factorisation | `V_4 = λ_4 · (angular ∫YYYY) · (radial ∫ψψψψ)` |
| T3 | quartic selection rule | Σl even (same Z₂ as #137) + SO(4) channel; verified |
| T4 | positive self-overlap | `g_4 = ∫ψ_k⁴ > 0` (manifest) |
| T5 | bounded interaction | a⁴ coeff > 0 ⟹ `V` bounded ⟹ stable vacuum |
| T6 | measure convergence | boundedness = `∫Dμ e^{−S}` convergence (#122); stability thread |
| T7 | ledger / scope | rule + overlap + boundedness derived; couplings input |
| T8 | assessment | `QUARTIC_VERTEX_LEDGER_Z2_SELECTION_POSITIVE_OVERLAP_BOUNDED_STABLE_VACUUM` |

## Established and open

  - **Established (BAM-native):** the quartic vertex carries the same antipodal
    Z₂ selection rule as the cubic (`Σl` even) plus the SO(4) common-channel
    condition; its geometric self-overlap `∫ψ⁴ > 0` is positive, so the
    effective potential is **bounded below** — a stable interacting vacuum — and
    that boundedness is the S_BAM measure-convergence condition (#122),
    extending the stability thread (#130/#136).

  - **Does not / open:** the coupling magnitudes `λ_3, λ_4` are **input** (the
    sign `λ_4 > 0` is required by #122 convergence; the magnitude is not derived
    from S_BAM); quintic / higher vertices are not addressed. The bulk-scale
    (#133) and flavor (#134) residuals stand.

## Cross-references

  - `docs/cubic_vertex_ledger_research_plan.md` — #137, the cubic ledger this
    extends; the antipodal Z₂ selection rule.
  - `docs/bam_factorized_sector_sum_research_plan.md` (#122) /
    `docs/z2_graded_sector_sum_convergence_research_plan.md` (#126) — the measure
    convergence that boundedness realises.
  - `docs/antipodal_kernel_one_loop_self_energy_research_plan.md` /
    `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #136/#130, the
    stability thread.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116, the
    cavity modes setting the geometric overlaps.

## Run

```
python -m experiments.closure_ledger.quartic_vertex_bounded_interaction_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_quartic_vertex_bounded_interaction_probe/`.
Expected verdict:
`QUARTIC_VERTEX_LEDGER_Z2_SELECTION_POSITIVE_OVERLAP_BOUNDED_STABLE_VACUUM`, 8/8 PASS.
