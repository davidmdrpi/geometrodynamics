# The v4 library migration: the flavor-CP lock lands in the library (PR #164)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. This PR
> performs the migration #163 staged; no new physics, no new inputs.

PR #163 verified the v4 candidate lock as probe-local code and flagged the
library migration as the successor, in four steps. This PR performs all
four — additively, so nothing is re-baselined.

## What moved into `geometrodynamics/qcd/quark_spectrum.py`

1. **Six new `QuarkParams` fields**, all default-off:
   - `phi_h` — the Hopf transport holonomy, applied to every
     same-partition off-diagonal as `e^{±i·σ(p)·φ_h·dk}`;
   - `eta_k1k3_plus`, `eta_k1k3_minus`, `eta_k1k5_minus` — the three
     targeted couplings completing the #161 target state (subtract
     convention, like the existing `eta_k3k5_minus`);
   - `diag_shift_plus`, `diag_shift_minus` — the per-shell diagonal
     retunes (indexed by `PASS_COUNTS = (1, 3, 5)`).
2. **The complexified same-partition transport**: `_offdiag_same_partition`
   now returns the holonomy-multiplied element. At `φ_h = 0` it is exactly
   real — the v3 spectrum is reproduced bit-for-bit.
3. **`LOCKED_QUARK_PARAMS_V4`** — `replace(LOCKED_QUARK_PARAMS, …)` at the
   derived `φ_h = π/k₅`, beside the **frozen** v3 lock.
4. **`extract_ckm_matrix(params)`** — diagonalizes the two partition blocks
   and returns `V = U₊† U₋`.

## The two views: mass vs mixing (the #158 relocation, in code)

The holonomy is a **pure phase**: `|transport|` is the closure cost that
sets the mass, `arg(transport)` is the Hopf holonomy that sets the mixing.
So:

- `extract_physical_spectrum` **strips** `φ_h` and reads the real
  closure-cost spectrum — the v4 lock **inherits the v3 masses exactly**
  (max relative drift ~3×10⁻⁹ across all six species).
- `extract_ckm_matrix` **keeps** `φ_h` and reads the misalignment of the
  two partition eigenbases — the CKM matrix, carrying the physical
  Jarlskog invariant.

## The unmixed reference (the species-labeling regime)

The three targeted couplings are **mixing** couplings — they generate the
CKM rotation — so they are treated exactly like `transport`/`γ_q`/
`partition_mixing`: `_unmixed_params` zeroes them, and the adiabatic
continuation ramps them on together from the unmixed reference. This keeps
the species-labeling reference cleanly block-diagonal (the three v4
off-diagonals — u–c, d–s, d–b — vanish exactly) and the ramp path
continuous from it. The structural level-repulsion `eta_k3k5_minus` is
**not** a mixing knob (it fixes the s/t masses across an enormous diagonal
gap); it stays on in both the reference and every ramp step, so the v3
mass-extraction path is bit-for-bit untouched.

## What the library now realizes (verified in the probe and the tests)

| property | result |
|---|---|
| v3 Hamiltonian at the default | exactly real (max \|Im\| = 0) |
| v3 CKM at the default | real rotation, J = 0 (no CP) |
| v4 masses vs v3 | inherited to ~3×10⁻⁹ |
| \|V_us\|, \|V_cb\|, \|V_ub\|, \|V_td\|, \|V_ts\|, J | all ≤ 1% |
| (β, γ, α) | (22.3, 65.9, 91.8)° vs (22.2, 65.9, 91.9)° |
| sin δ | 0.889 vs 0.887 |
| unitarity | dev ~7×10⁻¹⁶ |

## The additive regression decision (step 4)

The v3 `LOCKED_QUARK_PARAMS` is kept **frozen**; `LOCKED_QUARK_PARAMS_V4`
sits beside it. Every PR #155–#162 probe pins to the v3 lock, so the
migration needs **no re-baseline** — those probes stay bit-reproducible,
and the v4 lock is pinned by `tests/test_quark_v4_lock.py` (19 tests,
including the unmixed-reference regression). The full suite is **245
passed, 1 xfailed**.

## The counting (unchanged from #163)

+3 parameters buy +5 independent observables (V_us, V_cb, V_ub, β, γ — the
other four follow from unitarity and the derived phase): net predictive
surplus **+2**. The entire CP sector costs **zero** parameters
(φ_h = π/k₅ derived, #158–#160). The #150 input budget is unchanged.

## Reproduce

```bash
python -m experiments.closure_ledger.v4_library_migration_probe
# Verdict: V4_LOCK_MIGRATED_LIBRARY_DEFAULT_OFF_MASSES_INHERITED_NINE_OBSERVABLES
pytest tests/test_quark_v4_lock.py -q
```
