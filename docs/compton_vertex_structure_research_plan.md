# Compton vertex-structure probe — research plan

Follow-on to the finite-energy Klein-Nishina probe (PR #29, merged).
That probe established that the natural BAM amplitude (`(G_s + G_u)²
· (1 + cos²θ)`) reproduces KN at Thomson but fails at leading order
in ε = ω/m_e. The structural gap was localised to the missing
**per-channel kinematic weighting** — QED's vertex factors that
contract photon polarization with momentum (ε·k structures).

The previous probe also gave the precise analytic discrepancy:

    f_BAM − f_KN  ≈  +ε·(1 − cos θ)·(1 + 3cos²θ)/2  +  O(ε²)

at fixed θ in the small-ε expansion. The vertex modification must
subtract this specific θ-dependence at O(ε).

## The question

Does a natural BAM-derivable vertex factor — built from
combinations of `(ε·ε'*)`, `(ε·k̂')`, `(ε'*·k̂)`, channel-kinematic
weights — close the finite-energy gap and reproduce KN at all
orders in ε?

If yes (and the optimal coupling is a clean BAM-derived value): the
QED Compton amplitude is reproducible from BAM ingredients.

If no (or the optimal coupling is fine-tuned to a non-natural
value): the structural gap is deeper than vertex factors and the
thread needs another piece (Dirac spinor structure, second-order
connection coupling, etc.).

## Parametric vertex ansatz family

The probe tests four families of vertex modification to the baseline
`V₀ = ε(k)·ε'*(k')`:

**Family A — additive ε·k coupling**

    V_A(α)  =  ε(k)·ε'*(k')  +  α · (ε(k)·k̂') · (ε'*(k')·k̂)

This is the natural QED-style vertex contraction. For physical
photon polarisations summed over,

    Σ_λ |ε^(λ)(k)·k̂'|²  =  sin²θ
    Σ_λ Re[ε^(λ)·ε^(λ')* · (ε^(λ)·k̂')(ε^(λ')*·k̂)*]  =  −cos θ · sin²θ

so this modification adds `O(sin²θ)` corrections to the angular
factor.

**Family B — angular-modulation of V₀**

    V_B(β, γ)  =  ε(k)·ε'*(k') · (1 + β·sin²θ + γ·(1−cos θ))

Tests whether a multiplicative angular correction (not necessarily
from explicit ε·k contractions) closes the gap.

**Family C — per-channel kinematic weighting**

    M_s  =  G_S3(ψ_s)^p · exp(iφ_s) · V₀
    M_u  =  G_S3(ψ_u)^p · exp(iφ_u) · V₀

with p ∈ {0.5, 1, 1.5, 2}. The natural BAM construction uses p = 1;
the QED scaling `|M_s|² ∝ 1/(s−m²)` suggests the BAM analog might
need p = 0.5 (since `|M|²` gets squared).

**Family D — combined**

The most promising candidate from A/B/C, combined: per-channel
kinematic weight AND vertex factor. Tests whether the gap closes
when both pieces are present.

## Tests

### T1 — Thomson preservation

Whatever vertex modification is tested, the Thomson limit
(1 + cos²θ)/2 must be reproduced at ε → 0. Any ansatz that breaks
Thomson is rejected outright.

### T2 — Family A coupling sweep

Sweep α ∈ [−2, 2], compute max KN residual at finite ε (e.g.
ε = 0.1). Report:
  - α_opt: minimizing α
  - residual_at_opt: KN match quality
  - whether α_opt has a clean value (0, ±1, ±½)

### T3 — Family B coupling sweep

2D sweep over (β, γ) ∈ [−2, 2]². Same outputs as T2.

### T4 — Family C power sweep

Test p ∈ {0.5, 1.0, 1.5, 2.0}. Find optimal p.

### T5 — Family D best-combined

Take the best from A/B/C, combine if applicable, evaluate KN
residual on a full (ε, θ) grid. Report whether the combined ansatz
achieves < 1 % residual at all sampled (ε, θ).

### T6 — Optimal coupling interpretation

For each family, identify whether the optimal coupling has a
natural BAM derivation (Hopf connection-mediated coupling,
throat-transport-mediated, etc.) or is empirically fitted.

## Predicted outcomes

**Best case (FULL_CLOSURE)**: Family A with α = ±1 gives < 1 % KN
residual. Identifies the ε·k vertex coupling as the missing
ingredient with a clean unit-strength value, possibly derivable
from Hopf-connection coupling.

**Likely case (PARTIAL_CLOSURE)**: Some family gives ≤ 10 %
residual at finite ε, but exact KN requires a combination or higher-
order corrections. The probe identifies which family is closest and
quantifies remaining gap.

**Pessimistic case (NO_CLOSURE)**: No natural vertex ansatz brings
the residual below the baseline (PR #29 value, ~50–100 %). The
structural gap is deeper — likely requires explicit Dirac spinor
structure or a non-trivial photon throat-pair representation.

## What this leaves open even if the probe finds clean closure

  - **Derivation of the optimal coupling from BAM principles.**
    Empirical fit ≠ derivation. Even with α = 1 giving exact KN,
    the question "why α = 1?" needs answering from Hopf-connection
    or throat-transport algebra.
  - **Loop corrections.** Vertex, self-energy, vacuum polarization
    — still require BAM's bulk radial channel.
  - **Electron spin at finite ω.** The probe still uses scalar
    electron (Thomson choice). Spin-½ Dirac structure at finite ω
    is a separate target.
  - **Lorentz covariance.** Probe is set up in the electron rest
    frame.

## Cross-references

- `experiments/closure_ledger/compton_finite_energy_kn_probe.py`
  — predecessor establishing the O(ε) gap (PR #29).
- `experiments/closure_ledger/compton_photon_structure_probe.py`
  — predecessor establishing Thomson match (PR #28).
- `geometrodynamics/hopf/connection.py` — natural source of
  Hopf-mediated vertex coupling.
- `geometrodynamics/embedding/transport.py` — throat-transport
  source of `T = iσ_y`; vertex-coupling derivation might involve
  T-conjugation algebra.
- `experiments/closure_ledger/compton_vertex_structure_probe.py`
  — this probe.
