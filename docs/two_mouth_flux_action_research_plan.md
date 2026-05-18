# Two-mouth flux/action derivation of the Padé factor — research plan

Follow-on to PR #38 (geometric decomposition of `F²(x, c)`). PR #38
identified the closed-form Compton vertex factor as

    F²(x, c) = K(x)² · Q(x, c)

with `K(x) = 2x/(1+x)` (the "caustic Padé") and verified that the
harmonic-mean interpretation is consistent with the algebra. This
probe goes further: derives `K(x)` **from a concrete BAM throat
model** rather than postulating it.

## The model

A BAM transaction between two antipodal mouths on `S³`:

```
  mouth 1                      mouth 2
  ω₁ = ω  ◀── throat conduit ──▶  ω₂ = ω'
   |                                  |
   └──── closed orbit (2π action) ────┘
```

  - Two mouths at antipodal points on `S³`, connected by a throat.
  - Each mouth couples to a local cavity-mode oscillator (the
    `CavityMode` of `geometrodynamics.transaction.cavity`, with
    eigenfrequency `ω`).
  - A photon worldline traverses the closed orbit: mouth 1 → throat →
    mouth 2 → throat → mouth 1.
  - **Action at each mouth**: `S_i = ℏ·ω_i·τ_i` where `τ_i` is the
    proper time spent at mouth `i`.

## The principle: closure quantum splits equally

The BAM closure constant `action_base = 2π` (verified in
`experiments.closure_ledger.ledger`) is the action of one S³
great-circle. For one full closed orbit through both mouths, the
accumulated action must equal one closure quantum:

    S_total = S_1 + S_2 = 2π·ℏ

The **two-mouth flux-continuity postulate** is the natural
strengthening: action accumulates *equally* at each mouth (otherwise
one mouth would dominate the transaction). Equivalently, the flux
of action through the throat conduit is continuous:

    S_1 = S_2 = π·ℏ
    ω_1·τ_1 = ω_2·τ_2 = π          (each mouth carries half the closure quantum)

This is the conservation of action-flux at a bottleneck — the same
principle that governs any flux through a varying-cross-section
conduit (acoustic horns, hydraulic flow, transmission lines).

## The derivation

With equal-action splitting:

    τ_i = π / ω_i
    T_period = τ_1 + τ_2 = π·(1/ω_1 + 1/ω_2) = π·(ω_1 + ω_2)/(ω_1·ω_2)

The effective angular frequency of the closed orbit is

    ω_eff = 2π / T_period = 2·ω_1·ω_2 / (ω_1 + ω_2) = H(ω_1, ω_2)

— **exactly the harmonic mean** of the mouth frequencies. Normalising
to the incoming photon frequency `ω_1`:

    K(x) = ω_eff / ω_1 = 2·ω_2/(ω_1 + ω_2) = 2x/(1+x)   with x = ω_2/ω_1

The Padé factor `K(x) = 2x/(1+x)` is **derived**: the harmonic-mean
of the two mouth frequencies, normalised to the incoming frequency.
The squared factor `K(x)²` in `F²` arises because two amplitudes
(one per mouth) combine multiplicatively.

## Equivalent classical derivations (sanity)

The harmonic mean is the unique flux-continuous combination for two
series elements. Three convergent classical pictures:

  - **Series impedance**: two impedances `Z_i = 1/ω_i` in series →
    `Z_total = 1/ω_1 + 1/ω_2`. Effective admittance =
    `1/Z_total = ω_1·ω_2/(ω_1+ω_2)`. The flux through the series
    chain is the harmonic-mean rate.

  - **Acoustic / hydraulic horn**: a fluid passing through a
    varying-cross-section conduit has time-averaged velocity equal
    to the harmonic mean of the segment velocities (under
    incompressible-flux continuity).

  - **Wheeler–Feynman handshake**: offer wave from mouth 1
    (frequency ω_1) and confirm wave from mouth 2 (frequency ω_2)
    overlap at the rate of their harmonic mean — the unique rate
    consistent with both endpoint phase conditions.

## Tests

  T1. **BAM closure quantum value**: verify `action_base = 2π`
      from `experiments.closure_ledger.ledger.compute_lepton_ledger`
      (or `S3_ACTION_BASE` fallback).

  T2. **Equal-action splitting algebraic derivation**: starting
      from `ω_1·τ_1 = ω_2·τ_2 = π`, derive `K(x) = 2x/(1+x)`.

  T3. **Numerical two-segment orbit simulation**: integrate the
      piecewise oscillator phase from t = 0 to one full period,
      using equal-action splitting at the two mouths. Measure the
      effective period and verify `K(x) = ω_eff/ω_1 = 2x/(1+x)` at
      a grid of `(ω_1, ω_2)` ratios.

  T4. **Alternative splitting postulates rejected**: try
      `τ_1 = τ_2` (equal-time), `τ_1·ω_1² = τ_2·ω_2²`
      (equal-energy-squared), and `τ_1/τ_2 = ω_2/ω_1` (linear-ratio
      swap). Show none reproduce `2x/(1+x)`.

  T5. **Series-impedance equivalence**: solve the classical
      series-impedance problem with `Z_i = 1/ω_i` and verify the
      effective rate equals the harmonic mean.

  T6. **K² in F²**: cross-check that the derived `K(x)` reproduces
      the Padé factor in `F²(x, c) = K(x)² · Q(x, c)` to machine
      precision.

  T7. **Cross-process analytic continuation**: under crossing
      `x → x_⊗ < 0`, the derived `K(x_⊗) = 2x_⊗/(1+x_⊗)` is the
      analytic continuation; verify it has a pole at `x_⊗ = -1`
      (the throat-closure breakdown locus) and otherwise agrees
      with the closed-form Padé factor.

## Verdict structure

  - **PADÉ_DERIVED**: T1, T2, T3, T5, T6 all pass at machine
    precision; T4 rejects all alternative splittings; T7 confirms
    the analytic continuation. The Padé factor is derived from
    BAM closure-quantum splitting + equal-action throat flux.

  - **DERIVATION_AMBIGUOUS**: T4 fails to reject some alternative
    splitting; the equal-action postulate is not unique.

  - **DERIVATION_BROKEN**: T2 or T6 fails; the model does not
    produce the Padé factor.

## What this leaves open

  - **Why exactly equal action?** The equal-action postulate is
    *physically* the most natural (flux continuity), but a deeper
    BAM derivation would tie it to a specific S³ throat action or
    Lagrangian. The probe verifies the postulate's consequences,
    not its first-principles necessity.

  - **Loop corrections**: still tree-level.

  - **Q(x, c) derivation**: PR #38 identified the Hopf-fibre
    helicity-transport meaning of the `Q` factor; this probe only
    derives the `K` factor. A complementary "Hopf-fibre helicity
    transport" probe (mirror of this one) would close the
    derivation thread.

## Cross-references

  - PR #38: `throat_nucleation_caustic_derivation_probe.py` —
    F² = K(x)²·Q(x, c) decomposition.
  - `geometrodynamics/transaction/cavity.py` — `CavityMode`,
    `closure_check`.
  - `experiments/closure_ledger/ledger.py` — `action_base = 2π`.
  - `experiments/closure_ledger/two_mouth_flux_action_probe.py`
    — this probe.
