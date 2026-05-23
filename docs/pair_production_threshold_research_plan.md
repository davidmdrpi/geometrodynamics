# Pair-production / throat-nucleation threshold probe

The last open follow-on of the B4-anchor arc (PRs #53–#57): the THESIS
target *"the pair-production threshold falling out as the lowest stable
configuration."* In BAM a particle is a throat — a stable equilibrium
`R*` of the self-energy functional (PR #55), with rest energy
`E(R*) = m_e c²`. This probe derives the **pair-production threshold** as
`2 m_e c²`: twice the lowest stable throat, forced into a C-conjugate
pair by charge/topology conservation, with the nucleation barrier giving
the *disperse-below / persist-above* dichotomy and the Schwinger critical
field connecting the throat scale to the threshold.

## The picture

  - **Lowest stable configuration.** The self-energy `E(R) = A/R + B·R²`
    (EM repulsion + brane tension, PRs #55–#57) has its global minimum at
    `R*`, with `E(R*) = m_e c²` — the ground-state throat is the
    electron.

  - **Pairs from charge/topology conservation.** A throat carries one
    unit of Hopf charge (`|c₁| = 1`, `geometrodynamics/hopf/chern.py`).
    The vacuum has `c₁ = 0`. Conservation forces creation as a
    **C-conjugate pair**: throat (`c₁ = +1`) + antithroat (`c₁ = −1`,
    the inner/outer swap / antipodal `Z₂`, B2). Single-throat creation is
    forbidden; `Σ c₁ = 0`.

  - **Threshold.** By C-symmetry both partners have rest energy `E(R*) =
    m_e c²`, so the pair-production threshold is

    ```
    E_thr = 2 E(R*) = 2 m_e c² = 1.022 MeV .
    ```

  - **Nucleation barrier (disperse vs persist).** Creating a throat from
    the vacuum is a bubble-nucleation problem: a surface (brane-tension)
    cost competes with a volume energy gain,
    `E_nuc(R) = 4πσ R² − (4/3)π ρ R³`, with a critical radius
    `R_c = 2σ/ρ`. Below `R_c` the configuration shrinks back to vacuum
    (the antipodal focus disperses); above `R_c` it grows to the stable
    `R*` (a real throat persists) — exactly the THESIS dichotomy.

  - **Schwinger critical field.** Pair production in an external field
    requires the field to do work `2 m_e c²` to separate the pair. The
    Schwinger critical field `E_S = m_e²c³/(eℏ)` is where the work over a
    reduced Compton wavelength (`= R_MID`) equals `m_e c²`:
    `e E_S · R_MID = m_e c²`. This ties the throat scale `R_MID = λ_C` to
    the threshold.

## B4 accounting

The threshold `2 m_e c²` rides on the single anchor
`m_e c² = ℏc/R_MID` (the rest energy at `R*`). The **factor 2**, the
pair structure (`Σ c₁ = 0`), the disperse/persist dichotomy, and the
Schwinger-field relation (`e E_S R_MID = m_e c²`) are all derived /
dimensionless; the absolute scale is the one anchor — consistent with the
B4 scale-modulus theorem (PR #52).

## Tests

  T1. **Lowest stable configuration.** `E(R) = A/R + B·R²` has its global
      minimum at `R*` with `E(R*) = m_e c²` (the ground-state throat).
  T2. **Charge/topology conservation → pairs.** `|c₁| = 1` (the throat),
      `c₁ = 0` (vacuum); creation as a C-conjugate pair, `Σ c₁ = 0`;
      single-throat creation forbidden.
  T3. **Threshold = 2 m_e c².** `E_thr = 2 E(R*) = 2 m_e c² = 1.022 MeV`
      — the minimal on-shell pair energy.
  T4. **Nucleation barrier.** `E_nuc(R) = 4πσR² − (4/3)πρR³`: critical
      radius `R_c = 2σ/ρ`; shrinks below `R_c` (disperses), grows above
      (persists).
  T5. **Schwinger critical field.** `E_S = m_e²c³/(eℏ) ≈ 1.32×10¹⁸ V/m`;
      `e E_S R_MID = m_e c²` ties the throat scale to the threshold.
  T6. **Sub-threshold dispersal.** Below `2 m_e c²` no real pair forms
      (sub-critical fluctuation; the focus relaxes to vacuum).
  T7. **B4 accounting.** Factor 2 + pair structure + critical-field
      relation derived; absolute `2 m_e c²` rides on the single anchor.
  T8. **Assessment.**

## Verdict structure

  - **PAIR_THRESHOLD_DERIVED** (expected): the pair-production threshold
    is `2 m_e c²` — twice the lowest stable throat rest energy, forced
    into a C-conjugate pair by Hopf-charge / antipodal-`Z₂` conservation
    (`Σ c₁ = 0`). The nucleation barrier gives the disperse-below /
    persist-above dichotomy, and the Schwinger critical field
    (`e E_S R_MID = m_e c²`) ties the throat scale to the threshold. The
    factor 2 and the structure are derived; the absolute scale rides on
    the single anchor (B4-consistent).

  - **THRESHOLD_FAILS**: the threshold is not `2 m_e c²`, the charge
    bookkeeping does not close, or the nucleation barrier has no critical
    radius.

## What this leaves open

  - **The bulk-gravity / anchor scale.** `m_e c² = ℏc/R_MID` still rides
    on the one dimensionful input (the bulk gravitational scale, PR #57).
  - **The full dynamical nucleation rate.** The barrier here is the
    static critical-bubble picture; the tunneling/instanton rate
    (Schwinger exponential `e^{−πm²c³/(eEℏ)}`) from the BAM throat action
    is the follow-on.
  - **Heavier-lepton thresholds.** `2 m_μ c²`, `2 m_τ c²` from the
    excited radial throats (the closure ladder).

## Cross-references

  - `docs/self_consistent_throat_radius_research_plan.md` — `E(R*)`, the
    stable throat (#55).
  - `docs/cohesive_tension_derivation_research_plan.md` — the brane
    tension `σ` (#56).
  - `docs/brane_tension_tuning_research_plan.md` — the bulk-gravity
    tuning (#57).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `geometrodynamics/hopf/chern.py` — `compute_c1` (the Hopf charge).
  - `experiments/closure_ledger/throat_nucleation_caustic_derivation_probe.py`
    — the throat-nucleation caustic (earlier thread).
  - `experiments/closure_ledger/pair_production_threshold_probe.py` —
    this probe.
