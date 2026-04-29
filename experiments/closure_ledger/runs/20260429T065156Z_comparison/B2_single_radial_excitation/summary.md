# Closure-phase ledger — run summary

**Run:** 2026-04-29T06:51:56+00:00
**Experiment:** closure_ledger.layer2
**Transport convention:** `T2_sign_flip`
**χ (Hopf fibre):** 0.0
**S(k) candidate:** `B2_single_radial_excitation`

**Overall status:** Layer 2 FALSIFIES candidate 'B2_single_radial_excitation': lepton ledger does NOT close universally mod 2π (spread 3.495e+00).

## Constants used

| name | value | source |
|---|---|---|
| `action_base` | 6.283185 (2.000π) | repo |
| `beta_lepton` | 157.079633 (50.000π) | repo |
| `beta_quark` | 731.991088 (233.000π) | repo |
| `lepton_quanta` (4β/2π) | 100 | derived |
| `quark_quanta` (4β/2π) | 466 | derived |

## Per-lepton ledger

Phase contributions in units of π. `mod 2π` column shows the available-term sum reduced to [0, 2π).

| lepton | k | antipodal | Hopf | throat | uplift | radial | moving | total | mod 2π | status |
|---|---|---|---|---|---|---|---|---|---|---|
| electron | 1 | 2.000π | 1.000π | 1.000π | 0.000π | 0.882π | — | 4.882π | 0.881876π | partial_fails_to_close |
| muon | 3 | 6.000π | 1.000π | 1.000π | 0.000π | 1.994π | — | 9.994π | 1.994352π | partial_fails_to_close |
| tau | 5 | 10.000π | 1.000π | 1.000π | 200.000π | 2.999π | — | 214.999π | 0.999437π | partial_fails_to_close |

## Radial bulk channel — per-mode breakdown

| lepton | k | (l, n) | ω(l, n) | Φ(l, n) | status |
|---|---|---|---|---|---|
| electron | 1 | (l=1, n=0) | 1.054727 | 0.881876π | computed |
| muon | 3 | (l=1, n=1) | 1.974381 | 1.994352π | computed |
| tau | 5 | (l=1, n=2) | 2.894070 | 2.999437π | computed |

## Universality check (Layer 2)

- Per-lepton totals mod 2π (in units of π): ['0.881876', '1.994352', '0.999437']
- Spread across leptons: 3.49e+00
- Universal mod 2π (within 1e-9): **False**
- Universal value: N/A (not universal)

## Quark sector (structural)

- Lepton lock quanta: 100
- Quark lock quanta: 466
- Lock quanta gap: 366

_Both β-locks are integer-compatibility conditions for the closure ledger: each multiplier ensures the heaviest-shell uplift contributes 0 mod 2π. The lepton lock fits at 100 quanta of 2π (channels 2 and 3 nearly inactive); the quark lock fits at 466 (channels 2 and 3 active). The gap of 366 is a hypothesis about what the bulk-coupling channel contributes across the quark spectrum, contingent on S(k) being defined._

## Layer 2 blocker — S(k) bridge

**Verdict:** Φ_radial(k) is wired in the lepton sector via candidate `B2_single_radial_excitation` using `geometrodynamics.tangherlini.radial.solve_radial_modes` plus a Bohr-Sommerfeld integration of √max(ω² − V_eff, 0) over the tortoise grid. The remaining candidates are kept as open thesis-level alternatives. The quark-sector S(k) bridge is still undefined, so the 366-quanta gap remains S(k)-conditional.

### Evidence
- compute_knotted_lepton_spectrum accepts `l`, `n_points`, `rs`, `r_outer` for API compatibility but explicitly does not use them in the surrogate Hamiltonian.
- The locked lepton diagonal H_kk = action_base + resistance_scale·k² + res_diag(k) + pinhole(k∈{3,5}) + β·max(0, k−3)² has no eigenmode index and no Tangherlini potential evaluation.
- The quark residual sector (transport, pinhole, resistance) DOES read scalars off the tortoise grid — but those are single integrated quantities, not a per-generation S(k) → (l, n) mapping.

### Candidate S(k) maps

**A_lowest_radial_per_l**

- Formula: `S(k) = { (l, n=0) : l = 1, 3, ..., k }`
- Physical picture: Generation k is the sum of odd-l radial ground states up to angular harmonic l = k. Ties odd-k closure to odd-l angular content.
- Advantages: Consistent with the non-orientable transport rule: closure must flip the Z₂ partition class, and odd-l modes are the natural candidates for partition-flipping angular harmonics.
- Open questions: Result with WKB radial-action convention: FALSIFIES universal closure mod 2π. The cumulative-odd-ground-mode interpretation of lepton depth is rejected; whether a different phase convention (Maslov, Bohr-Sommerfeld with two soft turning points) revives it is open.

**B1_single_angular_mode**

- Formula: `S(k) = { (l = k, n = 0) }`
- Physical picture: Generation k is a single angular harmonic (l = k) in its radial ground state. One mode per generation, indexed by the angular quantum number alone.
- Advantages: Cleanest single-mode-per-generation interpretation: one angular eigenstate per lepton family. Easy to falsify because the per-row Φ is a single integrated quantity.
- Open questions: Does Φ(l = k, n = 0) close to a universal value mod 2π under the chosen WKB convention? The B1 phases are the same numbers that already appear in the candidate-A decomposition, so candidate-A's failure mode constrains B1 directly.

**B2_single_radial_excitation**

- Formula: `S(k) = { (l = 1, n = (k − 1) / 2) }`
- Physical picture: All generations share the lowest angular harmonic l = 1; depth labels successive radial excitations n = 0, 1, 2. Lepton depth is reinterpreted as radial-mode quantum number rather than angular content.
- Advantages: Uses one fixed angular sector and a single radial ladder, matching the surrogate's `β · k²` uplift if and only if ω²(1, n) ≈ ω²(1, 0) + β · n² (testable directly). Maximally falsifiable: a single mode integral per row.
- Open questions: Does Φ(l = 1, n) close to a universal value mod 2π? Asymptotically Φ(1, n) → (n + 1) π by WKB, so for high n the residues converge to a parity pattern {π, 0, π, 0, …} rather than a single universal value — universality requires a Maslov-shifted convention or a convention that absorbs the (n + 1) π structure.

**C_eigenvector_weighted**

- Formula: `S(k) = eigenvector / generation-block decomposition of the surrogate Hamiltonian onto Tangherlini (l, n) modes`
- Physical picture: The map S(k) is not imposed by hand. Instead it is derived: project the depth-k eigenvector of the lepton instanton surrogate onto Tangherlini radial modes, and use the resulting weights to build a coherent superposition phase Φ_radial(k).
- Advantages: Most physically principled: the bridge from instanton surrogate to Tangherlini operator becomes a derived quantity. Eliminates the by-hand selection of (l, n) shells that A, B1, and B2 all rely on.
- Open questions: Requires explicitly mapping `compute_knotted_lepton_spectrum`'s eigenvectors onto an (l, n) basis — currently the surrogate Hamiltonian is depth-only and uses no radial basis. The mapping has to be defined before the candidate is computable.

### Next steps
- Read the universality_check on the lepton ledger: candidate `B2_single_radial_excitation` either preserves closure mod 2π (supports the candidate) or breaks it (falsifies the candidate; rerun with another wired candidate).
- Define a quark-sector S(k) bridge. The lepton implementation operates on (l, n) Tangherlini modes; the quark sector needs its own per-generation mode assignment before P3 (the 366-quanta gap) becomes testable end-to-end.
- Investigate whether the candidate's per-mode Φ(l, n) values match the surrogate's β·k² uplift coefficients — this is the open-question line in candidate A's record.
- Update THESIS.md to reflect whichever way the universality check went: either close the ℏ open-problem bullet to a derived dimensionless invariant, or report which channel was overcredited.

### Downgraded predictions
- P3 (366-quanta gap as integrated bulk-mode contribution): downgraded from near-term falsification test to S(k)-conditional hypothesis. The number 366 = 466 − 100 is a structural prediction about what the bulk-coupling channel must contribute across the quark ladder relative to the lepton ladder, contingent on S(k) being defined for both.
