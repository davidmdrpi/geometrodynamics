# Lepton Sector Topological Axioms (Locked Baseline)

This document records the locked axioms now used as the default lepton baseline.

## 1) Locked geometric invariants

- `action_base = 2π` (`S3_ACTION_BASE`).
- `k_uplift_beta = 50π` (`TAU_BETA_50PI`).
- For `k=5`, uplift contribution is:
  - `ΔE_5 = beta * (5-3)^2 = 4*beta = 200π`,
  - i.e. exactly `100 × (2π)` topological quanta.

## 2) Matrix structure

### Diagonal sector

In the generation block, diagonal entries include a dominant `O(k^2)` term:

- `resistance_scale * (k**2)`,
- plus resistance profile, optional pinhole barrier, and uplift term.

This captures increasing geometric writhe/tension cost with generation complexity.

### Tunneling sector

- Default tunneling mode is `winding_mode="max"`.
- This scales off-diagonal suppression by crossing-number complexity of the deeper branch.

### Uplift sector

- Uplift is hard-coded by default as:
  - `ΔE_k = 50π * max(0, k-3)^2`.

## 3) Baseline solver anchors

The default baseline starting point for local optimization/probing is:

- `phase_per_pass ≈ 0.001`,
- `transport_strength ≈ 25.1`,
- `hard_pinhole_gamma ≈ 22.5`,
- `resistance_scale ≈ 0.217869`.

These values are exposed as baseline constants in the lepton module and used by
the basin/probe CLIs as defaults.

## 4) Downstream immutable mass vector

The API exports `solved_lepton_masses_mev()` which returns a read-only NumPy array:

- `[m_e, m_mu, m_tau]` in MeV from the locked baseline.

This gives downstream sectors a stable anchor for comparative fits and coupling
constraints without re-running scans.

## 5) Geometric origin of the locked anchors (post-hoc identifications)

After the closure-ledger probes (`experiments/closure_ledger/`), three
of the four solver anchors have geometric counterparts — see
`docs/odd_k_closure_lemma.md` for the depth-count and Layer-1
universality, and the `gamma_offset_probe` and `pinhole_origin_probe`
archives for the pinhole identification:

- `4 · k_uplift_beta = 100 · (2π)` is an integer count of antipodal
  closure quanta (the closure quantum). Geometric.
- `hard_pinhole_gamma ≈ 22.5` matches `Σ_{l=0..5} V_max(l)` ≈ 22.453
  on the canonical Chebyshev tortoise grid (within −0.21%, recovering
  the muon mass within 3.8% under the locked block). The `l = 0`
  channel is the 5D-specific centrifugal-free `3·rs²/r⁴` barrier
  unique to the Tangherlini metric. Geometric.
- `action_base = 2π` is the S³ great-circle action. Geometric.
- `transport_strength`, `resistance_scale`, `phase_per_pass` retain
  their phenomenological status — no geometric identification yet.

Together with the odd-k closure lemma, this leaves the lepton
diagonal fully geometric:

| row | source | value |
|---|---|---|
| `e (k=1)` | radial baseline + Hopf + throat | small (~6.9) |
| `μ (k=3)` | barrier-spectrum sum `Σ_{l=0..5} V_max(l)` | ≈ 22.45 |
| `τ (k=5)` | closure quantum `4β = 100·(2π)` | ≈ 628 |

The remaining sub-percent gap between the geometric μ-row anchor
(22.45) and the locked value (22.5) is recorded as a calibration
caveat in the probe archives.
