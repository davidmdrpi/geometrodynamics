# Closure-phase ledger — run summary

**Run:** 2026-04-29T14:06:17+00:00
**Experiment:** closure_ledger.layer2
**Transport convention:** `T2_sign_flip`
**χ (Hopf fibre):** 0.0
**S(k) candidate:** `C1_eigenvector_weighted_B1`

**Overall status:** Layer 2 FALSIFIES candidate 'C1_eigenvector_weighted_B1': lepton ledger does NOT close universally mod 2π (spread 3.257e-01).

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
| electron | 1 | 2.000π | 1.000π | 1.000π | 0.000π | 0.864π | — | 4.864π | 0.864195π | partial_fails_to_close |
| muon | 3 | 6.000π | 1.000π | 1.000π | 0.000π | 0.788π | — | 8.788π | 0.788396π | partial_fails_to_close |
| tau | 5 | 10.000π | 1.000π | 1.000π | 200.000π | 0.761π | — | 212.761π | 0.760506π | partial_fails_to_close |

## Radial bulk channel — per-mode breakdown

| lepton | k | (l, n) | ω(l, n) | Φ(l, n) | weight | weight·Φ | status |
|---|---|---|---|---|---|---|---|
| electron | 1 | (l=1, n=0) | 1.054727 | 0.881876π | 0.840935 | 0.741600π | computed |
| electron | 1 | (l=3, n=0) | 1.219083 | 0.770744π | 0.158673 | 0.122296π | computed |
| electron | 1 | (l=5, n=0) | 1.395973 | 0.760477π | 0.000393 | 0.000299π | computed |
| muon | 3 | (l=1, n=0) | 1.054727 | 0.881876π | 0.158844 | 0.140081π | computed |
| muon | 3 | (l=3, n=0) | 1.219083 | 0.770744π | 0.841087 | 0.648263π | computed |
| muon | 3 | (l=5, n=0) | 1.395973 | 0.760477π | 0.000069 | 0.000052π | computed |
| tau | 5 | (l=1, n=0) | 1.054727 | 0.881876π | 0.000221 | 0.000195π | computed |
| tau | 5 | (l=3, n=0) | 1.219083 | 0.770744π | 0.000241 | 0.000185π | computed |
| tau | 5 | (l=5, n=0) | 1.395973 | 0.760477π | 0.999538 | 0.760126π | computed |

## Universality check (Layer 2)

- Per-lepton totals mod 2π (in units of π): ['0.864195', '0.788396', '0.760506']
- Spread across leptons: 3.26e-01
- Universal mod 2π (within 1e-9): **False**
- Universal value: N/A (not universal)

## Quark sector (structural)

- Lepton lock quanta: 100
- Quark lock quanta: 466
- Lock quanta gap: 366

_Both β-locks are integer-compatibility conditions for the closure ledger: each multiplier ensures the heaviest-shell uplift contributes 0 mod 2π. The lepton lock fits at 100 quanta of 2π (channels 2 and 3 nearly inactive); the quark lock fits at 466 (channels 2 and 3 active). The gap of 366 is a hypothesis about what the bulk-coupling channel contributes across the quark spectrum, contingent on S(k) being defined._

## Layer 2 blocker — S(k) bridge

**Verdict:** Φ_radial(k) is wired in the lepton sector via candidate `C1_eigenvector_weighted_B1` using `geometrodynamics.tangherlini.radial.solve_radial_modes` plus a Bohr-Sommerfeld integration of √max(ω² − V_eff, 0) over the tortoise grid. The remaining candidates are kept as open thesis-level alternatives. The quark-sector S(k) bridge is still undefined, so the 366-quanta gap remains S(k)-conditional.

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

**C1_eigenvector_weighted_B1**

- Formula: `Φ_radial(k) = Σ_i |v_species(k),i|² · Φ(l = k_i, n = 0), weights from the locked lepton generation block`
- Physical picture: The depth basis {1, 3, 5} is shared with the instanton surrogate; each species' radial phase is the squared-amplitude weighted sum of the B1 ground modes (l = k_i, n = 0) over the depth basis. The B1 hand-imposed single-mode bridge is the |v|² = δ_ij limit of this map.
- Advantages: Derives weights from the existing lepton Hamiltonian without introducing fitted parameters. The eigenvector mixing lifts B1's degeneracy and produces tighter residues than any prior hand-imposed candidate.
- Open questions: Does the lifting close the residues to a single value mod 2π? If not, the eigenvector mixing alone does not supply the missing bridge under the WKB convention; test whether a different phase convention (Maslov, Bohr-Sommerfeld) revives universality given the tightened spread.

**C2_eigenvector_weighted_B2**

- Formula: `Φ_radial(k) = Σ_i |v_species(k),i|² · Φ(l = 1, n = (k_i−1)/2), weights from the locked lepton generation block`
- Physical picture: Same eigenvector weights as C1 but selecting B2's single-l radial-excitation ladder. Each species mixes across radial excitations n ∈ {0, 1, 2} of l = 1.
- Advantages: Tests whether the eigenvector mixing of B2's ladder (which under WKB asymptotes to (n+1)π) collapses the parity pattern across generations to a single value.
- Open questions: Same convention-dependence as C1; weight rows are the same in both candidates, only the per-mode Φ values differ. The two are independent tests of which mode ladder (B1 angular or B2 radial) the lepton Hamiltonian's eigenvectors actually align with.

### Next steps
- Read the universality_check on the lepton ledger: candidate `C1_eigenvector_weighted_B1` either preserves closure mod 2π (supports the candidate) or breaks it (falsifies the candidate; rerun with another wired candidate).
- Define a quark-sector S(k) bridge. The lepton implementation operates on (l, n) Tangherlini modes; the quark sector needs its own per-generation mode assignment before P3 (the 366-quanta gap) becomes testable end-to-end.
- Investigate whether the candidate's per-mode Φ(l, n) values match the surrogate's β·k² uplift coefficients — this is the open-question line in candidate A's record.
- Update THESIS.md to reflect whichever way the universality check went: either close the ℏ open-problem bullet to a derived dimensionless invariant, or report which channel was overcredited.

### Downgraded predictions
- P3 (366-quanta gap as integrated bulk-mode contribution): downgraded from near-term falsification test to S(k)-conditional hypothesis. The number 366 = 466 − 100 is a structural prediction about what the bulk-coupling channel must contribute across the quark ladder relative to the lepton ladder, contingent on S(k) being defined for both.
