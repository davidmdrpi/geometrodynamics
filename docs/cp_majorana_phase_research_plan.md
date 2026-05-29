# CP / Majorana phase probe (PR #94)

Follows PRs #92–#93, which fixed the PMNS mixing-angle structure
(anarchic, with θ13 suppressed by a residual nearest-neighbour
alignment). The remaining sector is the three CP-violating phases: the
Dirac phase `δ_CP` (seen in oscillations) and the two Majorana phases
`α21, α31` (seen only in lepton-number-violating processes, e.g. 0νββ).
This probe gives the BAM-native statements about them.

## CP violation is generic (the Hopf phase)

CP conservation requires the PMNS matrix to be real (up to rephasing) —
a **measure-zero** condition. In BAM the winding (charged-lepton)
amplitudes carry the Hopf holonomy `e^{ikχ}` (PR #60: the throat Berry
phase `∮A = π cos χ`), so they are intrinsically **complex**; the
cross-channel overlaps that build the PMNS matrix (PR #92) are therefore
generically complex, and `δ_CP ≠ 0, π` with probability 1. No BAM
symmetry forces real amplitudes, so CP is generically violated.

## The Jarlskog dichotomy (CP analogue of the angle dichotomy)

The rephasing-invariant measure of Dirac CP violation is the Jarlskog
invariant `J = Im(U_e1 U_μ2 U*_e2 U*_μ1)`, with `|J| ≤ 1/(6√3) ≈ 0.096`.
For an anarchic (Haar-random) `U(3)`, `|J|` has median ≈ 0.025.

| | observed \|J\| | anarchy percentile (μ=0) | residual (μ≈3) | reading |
|---|---:|---:|---:|---|
| **PMNS** | 0.026 | 51st | 81st | typical → large CP violation |
| **CKM** | 3.1×10⁻⁵ | 0.08th | — | extremely atypical → aligned/suppressed |

So the Jarlskog mirrors the angles (PR #92): PMNS is anarchic
(cross-coordinate, large CP violation), CKM is aligned (intra-coordinate,
suppressed CP violation). CP conservation (`J = 0`) is measure-zero — CP
is generically violated.

## The two Majorana phases exist because c₁ = 0

A Dirac neutrino has one physical phase (`δ_CP`) and no Majorana phases;
a Majorana neutrino has two extra physical phases. In BAM the neutrino is
Majorana precisely because it is chargeless (`c₁ = 0`, C-invariant,
PR #86), so the **existence** of two physical Majorana phases is a firm
BAM prediction — CP-violating phases of the ΔL=2 throat↔antithroat sector
(the bounce of PRs #87–#90), observable in 0νββ. Were the neutrino Dirac,
there would be none.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | angles done (#92–#93); 3 phases (δ_CP + 2 Majorana) remain |
| T2 | CP generic | winding amplitudes Hopf-complex (`e^{ikχ}`); CP-cons. measure-zero |
| T3 | δ_CP | anarchic (uniform for Haar U(3)); observed consistent |
| T4 | Jarlskog PMNS | \|J_PMNS\| ≈ 0.026 typical of anarchy (51st/81st pct) |
| T5 | Jarlskog CKM | \|J_CKM\| ≈ 3e-5 extremely atypical (~0.1th pct) = aligned |
| T6 | Majorana phases | two exist ⟸ c₁=0 (Majorana, PR #86); observable in 0νββ |
| T7 | honest scope | CP generic + Majorana existence firm; values anarchic |
| T8 | assessment | `CP_GENERIC_MAJORANA_PHASES_EXIST_VALUES_ANARCHIC` |

## Established and open

  - **Established (BAM-native):** CP violation is generic (the winding
    amplitudes are Hopf-complex; CP conservation is measure-zero); the
    Jarlskog dichotomy mirrors the angle dichotomy (PMNS CP violation
    typical of anarchy, CKM extremely atypical = aligned/suppressed); and
    two physical Majorana phases EXIST because the neutrino is Majorana
    (`c₁ = 0`, PR #86) — a firm prediction for 0νββ.

  - **Open:** the specific phase VALUES. `δ_CP` and the two Majorana
    phases are anarchic (uniform), set by the Hopf phases of the
    cross-channel overlaps and the throat↔antithroat tunnelling — not
    pinned. (`δ_CP` is itself poorly measured; the observed value is
    consistent with the uniform anarchic expectation.)

## Run

```
python -m experiments.closure_ledger.cp_majorana_phase_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_cp_majorana_phase_probe/`.
Expected verdict: `CP_GENERIC_MAJORANA_PHASES_EXIST_VALUES_ANARCHIC`, 8/8 PASS.
