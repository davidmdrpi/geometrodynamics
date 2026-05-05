# Closure-phase ledger — run summary

**Run:** 2026-04-29T06:12:22+00:00
**Experiment:** closure_ledger.layer1
**Transport convention:** `T2_sign_flip`
**χ (Hopf fibre):** 0.0

**Overall status:** Layer 1 PASS: lepton ledger universal mod 2π at 0.000π. Layer 2 BLOCKED: S(k) bridge undefined.

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
| electron | 1 | 2.000π | 1.000π | 1.000π | 0.000π | — | — | 4.000π | 0.000000π | partial_closes_mod_2pi |
| muon | 3 | 6.000π | 1.000π | 1.000π | 0.000π | — | — | 8.000π | 0.000000π | partial_closes_mod_2pi |
| tau | 5 | 10.000π | 1.000π | 1.000π | 200.000π | — | — | 212.000π | 0.000000π | partial_closes_mod_2pi |

## Universality check (Layer 1)

- All available-term sums mod 2π: ['0.000000', '0.000000', '0.000000']
- Spread across leptons: 1.42e-14
- Universal mod 2π (within 1e-9): **True**
- Universal value: 0.000000π

## Quark sector (structural)

- Lepton lock quanta: 100
- Quark lock quanta: 466
- Lock quanta gap: 366

_Both β-locks are integer-compatibility conditions for the closure ledger: each multiplier ensures the heaviest-shell uplift contributes 0 mod 2π. The lepton lock fits at 100 quanta of 2π (channels 2 and 3 nearly inactive); the quark lock fits at 466 (channels 2 and 3 active). The gap of 366 is a hypothesis about what the bulk-coupling channel contributes across the quark spectrum, contingent on S(k) being defined._

## Layer 2 blocker — S(k) bridge

**Verdict:** Φ_radial(k) is structurally blocked for the lepton sector: lepton_spectrum.py operates on depth labels k ∈ {1, 3, 5} via an instanton-transition surrogate, with no map from generation depth to Tangherlini eigenmodes (l, n).

### Evidence
- compute_knotted_lepton_spectrum accepts `l`, `n_points`, `rs`, `r_outer` for API compatibility but explicitly does not use them in the surrogate Hamiltonian.
- The locked lepton diagonal H_kk = action_base + resistance_scale·k² + res_diag(k) + pinhole(k∈{3,5}) + β·max(0, k−3)² has no eigenmode index and no Tangherlini potential evaluation.
- The quark residual sector (transport, pinhole, resistance) DOES read scalars off the tortoise grid — but those are single integrated quantities, not a per-generation S(k) → (l, n) mapping.

### Candidate S(k) maps

**A_lowest_radial_per_l**

- Formula: `S(k) = { (l, n=0) : l = 1, 3, ..., k }`
- Physical picture: Generation k is the sum of odd-l radial ground states up to angular harmonic l = k. Ties odd-k closure to odd-l angular content.
- Advantages: Consistent with the non-orientable transport rule: closure must flip the Z₂ partition class, and odd-l modes are the natural candidates for partition-flipping angular harmonics.
- Open questions: Does this map reproduce the locked lepton masses when Φ_radial(k) is added back into the diagonal Hamiltonian? If yes, the lepton instanton surrogate is a depth-block effective theory of this microscopic mode set.

**B_fixed_total_quantum_number**

- Formula: `S(k) = { (l, n) : 2n + l = k }`
- Physical picture: Generation k is a single shell in the joint (l, n) spectrum, with each generation mapping to one mode rather than a sum.
- Advantages: Clean single-mode-per-generation interpretation. Easier to falsify: each generation's Φ_radial value is a single eigenmode integral, with no sum-over-modes ambiguity.
- Open questions: Does the surrogate's β·k² term then read as the eigenfrequency ω(l, n) for the selected single mode? If so, this candidate gives a sharper microscopic interpretation than candidate A.

**C_closure_coherent_superposition**

- Formula: `S(k) = { (l, n) : ω(l, n) satisfies S³ closure at depth k }`
- Physical picture: The mode set is determined by which (l, n) combinations actually close on the antipodal cavity at depth k — the S³ closure condition selects the membership of S(k).
- Advantages: Most physically principled: makes S(k) emergent from the antipodal closure condition rather than imposed. Connects channel 1 (closure) and channel 3 (eigenmodes) directly.
- Open questions: Requires the S³ closure phase condition to be spelled out as an explicit equation on ω(l, n) before the membership of S(k) becomes computable. This is the highest-leverage candidate but also the highest-effort to define.

### Next steps
- Define S(k) for the lepton sector. Three candidates listed in `candidates`; selection is thesis-level, not implementation-level.
- Once S(k) is fixed, the radial phase extractor becomes well-posed: Φ_radial(k) = Σ_(l,n)∈S(k) ∫ k_local(r*) dr*.
- Test the three falsification predictions (P1 lepton universality, P2 quark universality, P3 the 366-quanta gap) against the implemented S(k).
- Update THESIS.md to reflect the result: either close the ℏ open-problem bullet to a derived dimensionless invariant, or report which channel was overcredited.

### Downgraded predictions
- P3 (366-quanta gap as integrated bulk-mode contribution): downgraded from near-term falsification test to S(k)-conditional hypothesis. The number 366 = 466 − 100 is a structural prediction about what the bulk-coupling channel must contribute across the quark ladder relative to the lepton ladder, contingent on S(k) being defined for both.
