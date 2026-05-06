# Quark β status — summary of the closure-ledger probe sequence

The locked quark sector uses

  `β_quark = N · π/2,  N = 466`

with the closure-quantum integer count `4β / (2π) = N`. Per
`docs/quark_axioms.md` §8 ("decisive N-ablation"), `N = 466` is **not
a topological invariant**: it drifts by tens of units under per-species
mass perturbations and anchor-species changes, while staying invariant
under uniform-mass scaling. The repo's pre-probe verdict was that `N`
is a fit compensator with no clean topological reading.

The closure-ledger experiment (initially built for the lepton-sector
S(k) bridge) led to a five-probe sequence that **narrowed the
structural status of `β_quark`**. This document records that
sequence's headline findings and the irreducible structural piece
that survives all of `quark_axioms` §8 ablations.

## Probe sequence

Each probe is a Python module under `experiments/closure_ledger/`,
with JSON + markdown archives under
`experiments/closure_ledger/runs/<timestamp>_<probe_name>/`. The
chronological order narrows the claimed structural content at each
step:

| # | probe                                | claim                                          |
|---|--------------------------------------|------------------------------------------------|
| 1 | `quark_beta_origin_probe`            | No principled enumeration produces N = 466 exactly. Suggestive: `(k_5²+1)·m ± 2`. |
| 2 | `quark_beta_boundary_probe`          | Cleaner reading: `N_q = m · k_5 + 1` with `m = 93, +1 = N_c − 2`. |
| 3 | `quark_beta_decomposition_probe`     | `m_q = (k_5−1)·k_5 + 2·k_5(k_5+2) + N_c = 93`; lepton block ⊂ quark formula. |
| 4 | `quark_beta_robustness_audit`        | The +1 invariance holds for only 33 % of §8 ablations; the structural reading is descriptive of the baseline, not predictive. |
| 5 | `quark_beta_subblock_stability`      | **Only `N_q ≡ 0 (mod 2)` is preserved across all 12 §8 ablations.** All sub-block subtractions and all other modular structures fail. |

## Headline finding

`N_q ≡ 0 (mod 2)` holds across every logged §8 ablation. Concretely,
the 12 §8 N values are

```
[466, 466, 466, 476, 474, 474, 482, 432, 494, 494, 440, 510]
```

— all even. The per-partition count `n_part = N_q / 2` takes the
values `[233, 233, 233, 238, 237, 237, 241, 216, 247, 247, 220, 255]`,
drifting `±22` units around the baseline `233`.

The structural reading is therefore

  `N_q  =  2 · n_part`

with the **factor of 2 topological** (the Z₂ partition multiplicity
from the v3 ansatz Hamiltonian basis `{(k, +), (k, −)}`) and `n_part`
**phenomenological** (it carries the entire §8 perturbation drift as
the fit compensator).

## What's principled-bounded vs derived

- **Principled-bounded**: yes. `N ∈ 2ℤ` — the parity constraint comes
  from the Z₂ partition class structure of the non-orientable throat.
  This survives every §8 ablation.

- **Derived**: no. `n_part = 233` at the baseline does not have a
  first-principles derivation in the present catalog. None of the
  scanned principled-enumeration families (S³ harmonics, S² harmonics,
  S³ angular eigenvalues, SU(3) representation dimensions, torus-knot
  invariants, Tangherlini barrier-derived counts) produces `233` or
  `466` in a way that survives the §8 robustness audit.

## How this refines `docs/quark_axioms.md` §8

§8's verdict was: "N is a compensator, not a topological invariant."
The probe sequence agrees, with one refinement: **the parity of N
is structural, the rest is fit.** Specifically:

- "Uniform mass scaling leaves N exactly at 466" → all §8 N values
  are even, including under uniform scaling. Consistent.
- "Per-species mass perturbations shift N by tens of units" → §8
  perturbed values [432, 494, 510, 440] remain even. The shift respects
  the parity constraint.
- "Anchor-species choice shifts N by 8–16 units" → [474, 476, 474,
  482] all even.

The Z₂ partition multiplicity is preserved across every documented
perturbation. No probe found a stronger invariant.

## Implications for the closure ledger

- The 366-quanta lepton/quark gap (`P3` in the Layer-2 blocker) is
  even (`366 = 2 · 183`), consistent with the parity constraint
  applied to both sectors. But `183` itself is not a documented
  topological invariant — the gap is parity-locked, not value-locked.

- The closure-quantum lock `4β_lepton = 100·(2π)` has `100` even, so
  the parity constraint is automatic for the lepton sector too. The
  quark sector merely extends the same structural lock at the
  coarsest level (parity), without inheriting the full integer-100
  topological reading.

- Closing the closure-phase ledger at Layer 2 (the unsolved
  S(k)-bridge problem from the original closure_ledger experiment)
  does not depend on the precise value of `N_q`. It depends on the
  per-mode radial phases on the Tangherlini grid, which the C1/D1
  candidates approached but did not close.

## What is now closed in `docs/THESIS.md`

The research target

  "Quark β = N · π/2 with N = 466 derived or principled-bounded.
   Either via torus-knot closure numbers on the non-orientable S³
   throat, or via representation dimensions in the throat condensate's
   color sector, or by some other finite enumeration the framework
   already supports. A clean negative result — 466 is not in the
   principled list — is also progress, and would refocus the search."

is closed at the **principled-bounded** level (parity), with the
**clean negative result** for the listed enumerations recorded by
probe #1. The full derivation of `n_part = 233` remains open.

## What's next

The natural follow-up is **not** to push further on `β_quark` — the
probe sequence has localized the irreducible structural piece, and
deriving `n_part` first-principles would require either an
enumeration framework not present in BAM (e.g. throat-condensate
color-sector representation theory developed independently) or a
non-perturbative QCD calculation that is well outside the BAM
machinery's current scope.

The larger open problem flagged in `docs/THESIS.md` "Open problems"
section is **the absolute action scale** — where does ℏ enter? The
recovered ratios are dimensionless; the absolute MeV scale is
imposed by anchoring `m_e`. A first-principles derivation of the
action quantum from the full closure cycle (antipodal closure × Hopf
holonomy × throat T² × radial bulk) would simultaneously close the
closure-phase ledger AND give a dimensionful prediction. That is the
natural next major work after this PR.

## Cross-references

- `experiments/closure_ledger/quark_beta_origin_probe.py` and archives
- `experiments/closure_ledger/quark_beta_boundary_probe.py` and archives
- `experiments/closure_ledger/quark_beta_decomposition_probe.py` and archives
- `experiments/closure_ledger/quark_beta_robustness_audit.py` and archives
- `experiments/closure_ledger/quark_beta_subblock_stability.py` and archives
- `docs/quark_axioms.md` §8 — original N-ablation evidence
- `docs/THESIS.md` — research roadmap and open problems
- `docs/odd_k_closure_lemma.md` — companion lemma for the lepton sector
