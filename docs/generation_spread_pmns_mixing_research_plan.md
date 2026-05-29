# Generation spread + PMNS mixing sector (PR #91)

Follows PR #90, which closed the neutrino-mass *scale* (`m_ν ~ few meV`)
from bulk throat geometry but left two residuals: a generation-uniform
bounce gives `m_ν ∝ m_D` (the cavity-floor ratios, a ×2.7 spread), and
the large PMNS mixing was untouched. This probe addresses both with
BAM-native structure.

## Generations = radial overtones; the bare prediction

In the unified operator (PR #83) the neutrino generations are the cavity
radial overtones `n = 0, 1, 2` (`k = 0`). A generation-uniform bounce
gives `m_ν,g ∝ m_D,g = √(ω²(0, n))`, the cavity floors:

```
m_ν,1 : m_ν,2 : m_ν,3  =  1 : 1.87 : 2.74   (normal ordering).
```

A falsifiable BAM prediction — normal ordering, masses in the
cavity-floor ratios. (Only `Δm²` are measured; the absolute scale is
unknown, so this is a prediction, not a fit.)

## The spread: overtone-dependent neck coupling

The bounce suppression `S_n` grows with how strongly overtone `n` couples
to the throat neck (PR #88). That coupling is exactly PR #79's
Z₂-antisymmetric boundary stress

```
χ_n = T_odd(n) = (T_inner − T_outer)/2,
```

which **decreases** monotonically with `n`:

| n (gen) | cavity floor m_D | χ_n |
|---:|---:|---:|
| 0 (1) | 1.055 | 0.304 |
| 1 (2) | 1.974 | 0.097 |
| 2 (3) | 2.894 | 0.039 |

Higher overtones fill the cavity more uniformly and feel the mouth
asymmetry less, so higher-`n` neutrinos are **less throat-coupled** ⟹ a
more compliant effective neck ⟹ a **smaller bounce** `S_n` ⟹ **less
suppression** ⟹ relatively **heavier**. This widens the spread in the
right direction: it lifts `m₃` relative to `m₂`, pushing the bare
`m₂:m₃ = 1:1.47` toward the `Δm²`-implied `1:5.8`. The mechanism and
direction are BAM-native; the exact magnitude needs the `ε_n(χ_n)`
relation (an `O(1)` coefficient, not derived).

## PMNS large vs CKM small: cross-channel vs intra-channel (headline)

The PMNS matrix is the overlap of the charged-lepton mass basis with the
neutrino mass basis. In BAM these live in **different channels** of the
PR #83 unified operator:

| sector | basis 1 | basis 2 | channels | mixing |
|---|---|---|---|---|
| **leptons (PMNS)** | charged: throat-winding (`k≠0`) | neutrino: cavity-resolving (`k=0`) | **different** | **large** |
| **quarks (CKM)** | up: cavity-shell (`k=0`) | down: cavity-shell (`k=0`) | **same** | **small** |

  - **Leptons** mix *across* the throat-winding / cavity-resolving divide
    (the two channels of PR #83), related only by the throat↔shell Z₂ /
    the +3 shift (PR #82). Two bases built from different quantum numbers
    (`k` vs `n`) are generically strongly misaligned ⟹ **large PMNS**.
  - **Quarks** (up-type and down-type) are *both* cavity-shell modes
    (`k=0`, `n≥3`; PR #85 quadrant, PR #82 quark pair) ⟹ same channel ⟹
    nearly aligned bases ⟹ **small CKM**.

So the long-standing puzzle — why is lepton mixing large and quark mixing
small — is the BAM **cross-channel** (leptons) vs **intra-channel**
(quarks) distinction.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | uniform `S` ⟹ `m_ν ∝ m_D` (×2.7); PMNS large |
| T2 | overtone prediction | `m_ν ∝ m_D` cavity floors `1:1.87:2.74` (normal) |
| T3 | spread mechanism | `χ_n` ↓ with `n` ⟹ higher-`n` less suppressed ⟹ heavier |
| T4 | PMNS large | cross-channel (charged `k≠0` × neutrino `k=0`) |
| T5 | CKM small | intra-channel (up & down both `k=0` shell) |
| T6 | data context | only `Δm²` measured; spread is a prediction |
| T7 | honest scope | direction + dichotomy structural; angles open |
| T8 | assessment | `PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING` |

## Established and open

  - **Established (BAM-native):** generations are radial overtones, so
    the bare prediction is normal ordering with `m_ν ∝ m_D`
    (`1:1.87:2.74`); the spread is widened in the right direction by the
    overtone-dependent neck coupling (PR #79 `χ_n` decreasing with `n` ⟹
    higher-`n` less suppressed ⟹ heavier); and large PMNS vs small CKM is
    the cross-channel (leptons: throat-winding × cavity-resolving) vs
    intra-channel (quarks: shell × shell) distinction — the BAM-native
    reason `PMNS ≫ CKM`.

  - **Open:** the precise mass spectrum (the `ε_n(χ_n)` coefficient is
    `O(1)`, not derived; the absolute scale is unmeasured — only `Δm²`),
    the explicit PMNS angles (need the cross-channel overlap integrals;
    large, but not computed to specific angles here), and the CP /
    Majorana phases.

## Run

```
python -m experiments.closure_ledger.generation_spread_pmns_mixing_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_generation_spread_pmns_mixing_probe/`.
Expected verdict: `PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING`, 8/8 PASS.
