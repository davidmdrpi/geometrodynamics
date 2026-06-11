# Quark CP phase calibration on the locked Hamiltonian (PR #156)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The CP phase
> rides the partition-mixing element of the mass-locked Hamiltonian.

PR #155 extracted the CKM |V| (zero new inputs) and posed the constrained CP
problem: the baseline has J = 0 exactly, the observed J = 3.08×10⁻⁵, and any
phase calibration must reproduce J **without disturbing the fixed |V| or the
mass calibration**. This PR performs the calibration — and finds sharper
structure than expected.

## The CP extension (locked blocks preserved)

The global `phase` knob **cannot** be used: it enters the same-partition
transport couplings as `cos(phase·dk)` — setting it to 0.5 collapses |V_us|
×20 — which is why the mass calibration locked it at zero and why #155's
baseline had J = 0 structurally. The v3 §4 partition-mixing element carries
its own (Hopf-placeholder) phase `φ_q(k) = φ·k`, and the locked params have
`partition_mixing = 0`, so CP switches on cleanly in-probe:

    H(ε, φ) = H_locked + Σ_k [ −ε e^{iφk} |k,+⟩⟨k,−| + h.c. ]

with the CKM read through the shell-partner charged current on the full 6×6.

## Results

- **Scaling derived**: J ∝ ε^1.9 (≈ quadratic — one insertion per sector
  side), sinusoidal sign-changing φ-dependence; |V|/mass shifts also O(ε²).
- **The ceiling identity**: J ≤ |V_us·V_cb·V_ub|. Predicted ceiling
  8.64×10⁻⁶ vs observed 3.47×10⁻⁵ — ratio **0.249 = 0.498 × 0.902 × 0.555
  exactly** (the per-element #155 soft-direction ratios). The J shortfall
  *is* the V_us/V_ub soft direction propagated — **no independent CP
  failure**. The observed CP is near-maximal (J_obs/ceiling_obs = 0.887).
  Consistency lock: when the soft directions land on data, the ceiling rises
  to the observed value.
- **The constrained calibration**: targeting the observed *phase content*
  (sin δ = 0.887) on the predicted |V| (J_target = 7.67×10⁻⁶):
  **ε\* = 0.0528, φ\* = 0.80** — |V_cb| shifts −0.0% (the stiff prediction
  untouched), |V_us|/|V_ub| −4.2%/−4.3% (inside the soft direction), masses
  ≤ 0.5% (inside the 1.6% calibration accuracy). **The locked structure
  survives.** Calibrated sin δ = 0.967 — near-maximal, like the data.
- **The triangle-shape discriminator (the honest sharp edge)**: the
  calibrated placeholder phase reproduces the *area* (J) but squashes the
  db-unitarity triangle — (β, γ) ≈ (0°, 180°) vs observed (22.2°, 65.9°).
  The linear `φ_q(k) = φ·k` puts the CP in the wrong quartet orientation.
  **β = 22° is the quantitative acceptance test for the true
  Hopf-connection phase** (the v3 §4 TODO tracked since the handoff).

## Ledger

- **Consumed**: one new input — the quark CP phase content (like α:
  structure derived, value input; the flavor puzzle's CP entry made
  explicit in the #150 budget).
- **Derived**: the quadratic switching, the ceiling identity with its exact
  soft-direction decomposition, the survival of the locked structure, the
  near-maximal phase content.
- **Open, with falsifiable targets**: the Hopf-connection `φ_q(k)` (target:
  β = 22°); the soft V_us/V_ub directions (target: J ceiling → 3.5×10⁻⁵).

## Reproduce

```bash
python -m experiments.closure_ledger.quark_cp_phase_calibration_probe
# Verdict: QUARK_CP_CALIBRATED_J_CEILING_DECOMPOSES_TRIANGLE_SHAPE_OPEN
```
