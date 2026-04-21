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
