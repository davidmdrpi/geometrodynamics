# Shell modes ↔ QCD/quark ladder: structural identification

Tests whether the shell-saturated radial modes identified in PR #68
reproduce the **structural ingredients** of the QCD/quark ladder. The
shell↔QCD match was flagged as future work in #68. Honest scope: test
the documented *structural* invariants of the quark sector — not the
phenomenological pieces (specific quark masses, `n_part = 233`) which
the quark β-status probes already showed are not derivable from BAM's
current catalog.

## What the quark sector requires (per `docs/quark_beta_status.md`)

The probe sequence that closed the quark β-status established **one
robust structural invariant**:

```
N_q ≡ 0 (mod 2),   i.e.   N_q = 2 · n_part .
```

This **factor of 2** is the **Z₂ partition multiplicity** of the
non-orientable throat (the v3 ansatz basis `{(k, +), (k, −)}`) and
survives every §8 ablation. Together with the standard QCD facts, the
structural ingredients to reproduce are:

  1. **Z₂ partition multiplicity** — the factor of 2 (two states per
     mode index): the symmetric realization of both mouths/partition
     classes.
  2. **Six quark flavors = 3 generations × 2** — the 3×2 counting (u/d,
     c/s, t/b).
  3. **Higher mass scale than leptons.**
  4. **Extended / shell-coupled character** (vs pointlike leptons) —
     the confinement-like nature of quark states.
  5. **β_quark = N·π/2 closure-quantum integer, with even N** — the
     parity invariant.

What is **not** part of the structural target (open per
`docs/quark_beta_status.md`): the specific value `n_part = 233`, exact
quark masses, the SU(3) color sector. Those remain phenomenological.

## What this probe finds

The shell modes (`n ≥ 3` at `l=1`, the saturated channel from PR #68)
reproduce the documented structural ingredients:

  - **Z₂ partition multiplicity (50/50 mouth balance).** For shell modes,
    the inner-half vs outer-half probability is `≈ 0.50 / 0.50` — both
    mouths realized symmetrically, the Z₂ partition fully active. The
    throat-focused electron (n=0) is `0.56 / 0.44` (asymmetric, single
    throat identification — the lepton has ONE mouth dominant). The
    transition through the μ, τ leads to the symmetric shell. The
    structural factor of 2 (`N_q = 2·n_part`) is **provided** by the
    shell mouth-symmetry.

  - **3 generations × 2 = 6 quark flavors.** The first 3 shell modes
    (`n = 3, 4, 5`) × the Z₂ doubling = **6** — the right structural
    count of quark flavors (u/d, c/s, t/b).

  - **Heavier than leptons.** `ω(shell)` starts at `≈ 3.83`, above the
    heaviest lepton `ω(τ) ≈ 2.89` — quarks heavier than leptons (the
    correct direction).

  - **Extended / shell-coupled** (#68): participation ratio `→ 2/3`
    (uniform shell standing wave) vs lepton `≈ 0.65` (focused) —
    confinement-like extended states.

  - **Even-N (parity) closure invariant.** `N_q ∈ 2ℤ` follows from the
    shell Z₂ partition — the structural invariant the quark β-sequence
    isolated, now provided geometrically by the shell modes.

## What this is and is not

  - **Is:** a structural identification of the shell-saturated modes
    with the documented quark-sector invariants: the Z₂ partition
    multiplicity (factor of 2 / `N_q ∈ 2ℤ`), the 3×2=6 flavor counting,
    the heavier-than-lepton mass scale, the extended (shell-coupled)
    character.
  - **Is not:** a derivation of specific quark masses, of `n_part = 233`,
    or of the SU(3) color sector. Those remain phenomenological /
    open, consistent with `docs/quark_beta_status.md`.

## B4 accounting

The structural invariants tested (mouth-balance, mode counting,
participation ratio) are **dimensionless ratios / integer counts**;
the identification is geometric, independent of the single anchor.
Mass values carry the scale; mass ratios and the structural counting
do not.

## Tests

  T1. **Inner/outer mouth balance: lepton asymmetry vs shell symmetry.**
      Lepton end (n=0=e) `0.56/0.44` (asymmetric, single throat
      identification); shell `≈ 0.50/0.50` (Z₂ partition fully realized,
      both mouths as distinct states → factor of 2).
  T2. **3 shell modes × 2 = 6 quark flavors.** the first 3 shell modes
      (`n=3,4,5`) × Z₂ doubling = 6 — the right structural quark count.
  T3. **Heavier mass scale.** `ω(shell) ≥ ω(n=3) ≈ 3.83 > ω(τ) ≈ 2.89`.
  T4. **Extended / shell-coupled** (recap #68): participation `→ 2/3`
      vs lepton focused — confinement-like extended states.
  T5. **β_quark parity (`N_q ∈ 2ℤ`) from the shell Z₂.** the
      shell mouth-symmetry provides the documented structural invariant
      of the quark β-lock.
  T6. **Honest scope.** structural ingredients reproduced;
      `n_part = 233`, exact masses, SU(3) color open.
  T7. **Falsification / B4.** asymmetric shell modes (no Z₂
      multiplicity) would falsify; BAM shows 50/50 shell. Metrics
      dimensionless/scale-independent.
  T8. **Assessment.**

## Verdict structure

  - **SHELL_REPRODUCES_QCD_STRUCTURE** (expected): the shell-saturated
    radial modes (PR #68, `n ≥ 3` at `l=1`) reproduce the documented
    structural ingredients of the quark sector: the **Z₂ partition
    multiplicity** (`≈ 0.50/0.50` inner/outer mouth balance, providing
    the structural factor of 2 in `N_q = 2·n_part`), the **3×2 = 6**
    flavor counting from the first 3 shell modes doubled, the
    **heavier-than-lepton mass scale**, and the **extended /
    shell-coupled** character. The phenomenological pieces (`n_part`,
    specific quark masses, SU(3) color) remain open, consistent with
    `docs/quark_beta_status.md`.

  - **NO_STRUCTURAL_MATCH**: shell modes do not show the Z₂ multiplicity
    or the right counting — the shell↔QCD identification fails.

## What this leaves open

  - **`n_part = 233`** (the phenomenological compensator) — the quark
    β-status doc already established that none of the principled
    enumerations in the catalog produces it; remains open.
  - **Exact quark masses.** The structural ladder match here does not
    address the specific `m_u, m_d, …` values.
  - **SU(3) color.** Not part of this probe.

## Cross-references

  - `docs/throat_to_shell_transition_research_plan.md` — the shell
    channel (#68).
  - `docs/quark_beta_status.md` — the documented quark-sector
    invariants (`N_q ∈ 2ℤ`, the structural factor of 2).
  - `docs/quark_axioms.md` — the §8 ablation evidence.
  - `experiments/closure_ledger/shell_to_qcd_match_probe.py` — this
    probe.
