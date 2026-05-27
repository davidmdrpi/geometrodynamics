# Quark `n_part = 233` origin: phenomenological compensator

Attempts the hardest open piece of the quark sector arc — the origin
of `n_part = 233` in the locked-quark closure integer
`N_q = 2 · n_part = 466` (so `β_quark = N_q · π/2 = 233π`). Lands on
the structural reading that has been implicit since
`docs/quark_beta_status.md`'s headline finding: **`n_part = 233` is a
phenomenological compensator at the v3 baseline, not a topological
invariant**.

## Honest target

The user explicitly flagged this as the hardest open piece. The
deliverable here is **structural classification**, not a
first-principles derivation:

  - Extend the prior candidate catalog (S³/S² harmonics, SU(3)
    representations, torus-knot crossings, Tangherlini barrier sums)
    with Fibonacci, Lucas, Padovan, Perrin, tribonacci, color × flavor
    × generation, QCD β₀ combinations, Tangherlini QCD-shell mode
    counts. **No exact-match enumeration survives `docs/quark_axioms.md`
    §8 drift.**
  - Identify `n_part = 233 = F_13` and `9·k_5² + k_5 + 3 = 233` as
    **baseline coincidences** (the §8-drift values `{216, 220, 237,
    238, 241, 247, 255}` are neither Fibonacci nor of the same `k_5`
    polynomial form).
  - Reframe the question **structurally**: the v3 quark Hamiltonian
    uses the lepton-shaped closure-quantum basis `{(k=1,±), (k=3,±),
    (k=5,±)}`, but the quark sector lives in the QCD shell channel
    (#68–#69). The v3 Hamiltonian is the **wrong machinery** for the
    quark sector; `n_part = 233` absorbs the unmodeled QCD physics
    (confinement, αs running, color) as a phenomenological
    compensator.
  - Identify the **right derivation route** — quantitative development
    of #68–#69 (throat-to-shell + shell↔QCD match) into a full
    QCD-shell model — as the genuine open work, and note this is
    structurally outside the closure-ledger machinery's scope.

## What survives the §8 audit

The prior `quark_beta_subblock_stability_probe` established that across
all 12 logged `docs/quark_axioms.md` §8 ablations, the closure integer
takes values

```
N_q ∈ {466, 466, 466, 476, 474, 474, 482, 432, 494, 494, 440, 510}
n_part = N_q / 2 ∈ {233, 233, 233, 238, 237, 237, 241, 216, 247, 247,
                    220, 255}
```

with min 216, max 255, span 39, mean 236.4. The only invariant is
**parity** (`N_q ≡ 0 mod 2`), interpreted as the Z₂ partition
multiplicity from the v3 basis `{(k, +), (k, −)}` (B2's
non-orientable-throat structure). `n_part = 233` itself is just one
realization of a fit compensator that absorbs the per-species
perturbations.

## The extended candidate catalog (this probe)

Beyond the prior probe's catalog, this probe scans **129 candidates
across 9 families**:

| family | entries | exact n_part match? |
|---|---:|:---:|
| Fibonacci `F_n` | 18 | `F_13 = 233` (coincidence) |
| Lucas `L_n` | 18 | none |
| Padovan `P_n` | 24 | none |
| Perrin `R_n` | 23 | none |
| tribonacci `T_n` | 12 | none |
| color × flavor × generation | 11 | none |
| QCD β₀ = 7 combinations | 5 | none |
| Tangherlini QCD-shell modes | 12 | none |
| `k_5` polynomials | 6 | `9·k_5² + k_5 + 3 = 233` (coincidence) |

Two baseline coincidences (`F_13 = 233`, `9·k_5² + k_5 + 3 = 233`); no
enumeration survives §8 drift because all are parameter-free in the
ablation variables.

## Why `F_13 = 233` is a coincidence

The 13th Fibonacci number `F_13 = 233` happens to equal the v3 baseline
`n_part`. But the §8 drift values are

```
216, 220, 237, 238, 241, 247, 255
```

None of these are Fibonacci (the relevant `F_n` are `…, 144, 233, 377,
…`, with gaps of 89 and 144 around 233). Only the baseline value lands
on `F_n`; the drift values do not. `F_13 = 233` is therefore a
**baseline coincidence**, not a structural invariant.

## The structural reading

| sector | lives in | basis | closure integer |
|---|---|---|---:|
| lepton | throat (odd-`k` modes) | `{(k=1,±), (k=3,±), (k=5,±)}` | `4·k_5² = 100` (clean closure quantum, PR #71) |
| quark | QCD shell channel (#68–#69) | (same v3 basis = wrong machinery) | `466` (phenomenological) |

The v3 quark Hamiltonian fits the quark spectrum on a **lepton-shaped**
basis (the 6 odd-`k` throat modes — the same modes that give
`β_lepton = k_5²·(2π) = 50π` via `4·k_5² = 100`). But per PRs #68–#69,
the quark mass sector is **delocalized into the QCD shell channel**:
extended-character wavefronts reproducing the documented quark-sector
structural invariants (Z₂ partition, `3 × 2 = 6` flavors, heavier mass
scale).

The v3 ansatz is therefore fitting QCD-confined quarks on
closure-quantum throat basis vectors — the wrong machinery — and
`n_part = 233` is the empirical price of that sector mismatch: it
absorbs unmodeled QCD physics (confinement, αs running, color) that the
closure-ledger primitives are designed for the lepton-throat sector,
not for the shell channel.

## The right derivation route

Quantitative development of #68–#69 (throat-to-shell transition +
shell↔QCD structural match) from "structural match" into a full
QCD-shell model — with confinement, αs running, and color sector
explicit — is the structurally honest path to `n_part`. This is:

  - **Outside the closure-ledger machinery's scope.** The
    closure-ledger primitives (`action_base = 2π`, `transport = 8π`,
    `resistance = 7π/100`, `pinhole γ = Σ V_max[1..5]`, `ε = 7π/(100·k_5⁴)`)
    are lepton-throat-sector primitives.
  - **A substantial research program in its own right** — comparable
    in scope to deriving lattice QCD's hadron spectrum from underlying
    geometric principles.
  - **Not the next-most-tractable work** in the current BAM framework.
    It should not be pursued by further enumeration on the v3
    Hamiltonian.

## Honest scope

  - **Is:** structural classification of `n_part = 233` as a
    phenomenological compensator at the v3 baseline; extension of the
    candidate catalog with eight new families (zero surviving exact
    matches); identification of `F_13 = 233` and
    `9·k_5² + k_5 + 3 = 233` as baseline coincidences; structural
    framing of *why* (lepton/shell sector mismatch); identification of
    the right derivation route (quantitative #68–#69 development) and
    its scope (substantial, outside closure-ledger machinery).
  - **Is not:** a first-principles derivation of `n_part = 233`. No
    such derivation exists in the current catalog. The prior verdict
    from `docs/quark_beta_status.md` ("phenomenological compensator")
    is upheld and sharpened, not closed.

## B4 accounting

`n_part` is a dimensionless integer; scale-independent. The
identification is structural (sector-mismatch reading) — independent
of the dimensionful anchor `m_e`.

## Tests

  T1. **Parity invariance recap.** `N_q ≡ 0 (mod 2)` across all 12 §8
      ablations (Z₂ partition multiplicity; from
      `quark_beta_subblock_stability`).
  T2. **Extended candidate catalog.** Fibonacci, Lucas, Padovan,
      Perrin, tribonacci, color × flavor × generation, QCD β₀,
      Tangherlini QCD-shell modes — 129 candidates, 9 families.
  T3. **F_13 = 233 baseline coincidence.** §8 drift values are not
      Fibonacci → coincidence, not invariant.
  T4. **v3 Hamiltonian is lepton-shaped.** Basis `{(k=1,±), (k=3,±),
      (k=5,±)}` = the same odd-`k` throat modes that give the lepton
      ladder.
  T5. **PRs #68–#69 = right route.** Quark sector lives in the QCD
      shell channel; quantitative #68–#69 development is the genuine
      derivation route — outside closure-ledger scope.
  T6. **§8 drift.** `n_part` ranges 216–255 (span 39, mean 236.4) —
      not a topological invariant.
  T7. **Honest scope / B4.** Structurally classifies; does NOT derive
      `n_part = 233` first-principles. `n_part` dimensionless integer;
      structural.
  T8. **Assessment.**

## Verdict structure

  - **N_PART_IS_PHENOMENOLOGICAL_COMPENSATOR** (expected): `n_part =
    233` at the v3 baseline is one realization of a phenomenological
    fit compensator that drifts across §8 ablations (216–255, span 39);
    the only topological invariant is parity (Z₂ partition
    multiplicity). No principled enumeration in the extended catalog
    reproduces `n_part` across §8. The v3 Hamiltonian is lepton-shaped
    machinery (basis = 6 odd-`k` throat modes giving the lepton
    ladder), but the quark sector lives in the QCD shell channel
    (#68–#69); `n_part` absorbs unmodeled QCD physics (confinement, αs
    running, color) as a phenomenological compensator. The right
    derivation route is quantitative #68–#69 development — outside
    closure-ledger scope, a substantial research program in its own
    right.

  - **N_PART_DERIVED**: a previously-untested candidate enumeration
    reproduces `n_part = 233` exactly AND survives §8 drift.

## What this leaves open

  - **A first-principles derivation of `n_part = 233`.** Genuine open
    work; no principled enumeration in the current catalog reproduces
    it across §8 drift.
  - **A quantitative QCD-shell model** (extending #68–#69 from
    "structural match" to "quantitative computation"). This is the
    structurally honest route; substantial research program; outside
    closure-ledger machinery's scope.
  - **Heavy-lepton thresholds** (`2 m_μ c²`, `2 m_τ c²`) and the
    analogous quark-sector pair-production thresholds — related open
    work flagged by PR #58.

## Cross-references

  - `docs/quark_beta_status.md` — the prior five-probe summary
    establishing the parity invariance and the phenomenological reading.
  - `docs/quark_axioms.md` §8 — the original N-ablation evidence with
    the 12 logged perturbations.
  - `docs/throat_to_shell_transition_research_plan.md` — PR #68,
    higher excitations delocalize into the QCD shell channel.
  - `docs/shell_to_qcd_match_research_plan.md` — PR #69, shell modes
    reproduce the documented quark-sector structural invariants.
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    lepton-sector analogue `β_lepton = k_5²·(2π) = 50π` with the
    clean closure-quantum integer `4·k_5² = 100`.
  - `docs/odd_k_closure_lemma.md` — the lepton-sector closure lemma.
  - `experiments/closure_ledger/quark_beta_origin_probe.py` — the
    prior catalog (S³/S²/SU(3)/torus-knot/Tangherlini).
  - `experiments/closure_ledger/quark_beta_subblock_stability_probe.py`
    — establishes the parity invariance.
  - `experiments/closure_ledger/quark_npart_origin_probe.py` — this
    probe.
