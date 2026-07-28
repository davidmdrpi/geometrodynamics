# The BAM Pin⁻ spinor spectrum and effective action on the cycle (PR #230)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. Why this exists: it corrects its own PR

`docs/twist_sector_selection_research_plan.md` (same PR) argued that the
sign of the twist selection is **the statistics of the link field**, using
the generic rule that the zero point is `+½Σω` for a boson and `−½Σω` for a
fermion, and concluded

```
fermionic link field  →  ΔE = −π/(4C)  →  TWISTED favoured
```

then read that as confirming that BAM matter sits in the twisted sector.

**That is wrong.** It *asserted* the statistics of the link field instead
of computing BAM's actual spinor — and in doing so silently gave the
spinor the same moding as the scalar. It does not have it.

## 1. The moding swaps

For a spinor the network holonomy `W` composes with the **intrinsic spin
structure** `η`, and BAM's B2 fixes the non-trivial one — `T² = −I`,
"antiperiodic 4π-spinors" — so `η = −1`. The boundary phase is `θ = 0`
when `W·η = +1` and `θ = π` otherwise:

| field | `η` | `W` | `W·η` | `θ` | moding |
|---|---:|---:|---:|---:|---|
| scalar | `+1` | `+1` | `+1` | 0 | integer |
| scalar | `+1` | `−1` | `−1` | π | half-integer |
| **BAM Pin⁻ spinor** | `−1` | `+1` | `−1` | π | **half-integer** |
| **BAM Pin⁻ spinor** | `−1` | `−1` | `+1` | 0 | **integer** |

The assignment is **exactly reversed**. That single omission is what made
the earlier sign wrong.

## 2. The spectrum

Computed through the repo's own SUSY framing: the throat Dirac operator is
the square root of the second-order operator
(`docs/throat_dirac_spinor_research_plan.md`), so on the cycle `D²` is the
ring Laplacian at the **spinor's** boundary phase and
`spec(D) = ±√(spec(Laplacian_θ))`.

Convergence to the analytic `|k_m| = |θ + 2πm|/C`, at `n = 200…1600`:

| sector | max error `n=200` | `n=1600` | successive ratios |
|---|---:|---:|---|
| `θ=0` (integer) | 6.98e-03 | 1.09e-04 | 4.00, 4.00, 4.00 |
| `θ=π` (half-integer) | 1.11e-02 | 1.73e-04 | 4.00, 4.00, 4.00 |

Second order, as the 3-point Laplacian requires.

> **Method note — a real trap.** A direct first-order central-difference
> Dirac operator was tried first and is **unusable here**: fermion doubling
> puts a spurious second zero mode at the Brillouin-zone edge, giving
> `dim ker D = 2` at `θ = 0` and destroying precisely the index statement
> §3 rests on. The SUSY route is doubling-free.

## 3. The index: the zero mode is in the **twisted** sector

The swap has a sharper consequence than any vacuum-energy sign. Integer
moding carries an exact Dirac zero mode, and for the BAM spinor that is
the **twisted** sector:

| `W` | `θ` | `dim ker D` | `det D = 2sin(θ/2)` | lowest `\|k\|` |
|---:|---:|---:|---:|---:|
| `+1` | π | **0** | 2.000000 | 3.141592 (`= π/C`) |
| `−1` | 0 | **1** | **0.000000** | 0.000036 |

A zero mode is **topological, not energetic** — it cannot be traded
against a `π/(4C)` vacuum-energy difference. This is #195's index
mechanism (Pin⁻ Dirac zero modes pinned by Atiyah–Singer, "not by
tuning") seen from the network side.

## 4. The effective action — and the correction

| field | `C` | `E₀` untwisted | `E₀` twisted | `ΔE` | favoured |
|---|---:|---:|---:|---:|---|
| scalar (boson) | 1.0 | −0.523599 | +0.261799 | **+0.785398** | untwisted |
| scalar (boson) | 2.0 | −0.261799 | +0.130900 | +0.392699 | untwisted |
| **BAM Pin⁻ spinor** | 1.0 | −0.261799 | +0.523599 | **+0.785398** | **untwisted** |
| **BAM Pin⁻ spinor** | 2.0 | −0.130900 | +0.261799 | +0.392699 | untwisted |

With the moding computed rather than assumed, the statistics flip and the
spin-structure flip **cancel**. Both fields give `ΔE = +π/(4C)`:
**energetics favour the untwisted sector universally.** The magnitude
`π/(4C)` per degree of freedom is unchanged (1e-12).

So there is **no "sign = statistics" rule**, and the earlier cross-check —
which read the fermionic sign as explaining why matter is twisted — does
not stand.

## 5. The corrected picture: two mechanisms, not one sign

| sector | link field | selected by | outcome |
|---|---|---|---|
| **QCD Möbius** (#100/#103) | flux-tube phonon | **ENERGY** | untwisted cheaper by `π/(4C)` ⟹ Möbius states are excitations *above* the orientable ground — how #100/#103 present them (`+πσ`, `+2√σ`) |
| **BAM matter** (#170/#195/#202) | Pin⁻ throat spinor | **INDEX** | energetics favour untwisted here **too**, so energy cannot explain matter; the twisted sector is selected because it carries the Dirac zero mode — matching #202/#227 |

Both sectors still come out right, but for **different and now correctly
identified** reasons. This is a better answer than the one it replaces: it
explains *why* the two sectors behave differently, instead of attributing
it to a sign that does not exist.

**Untouched:** the topological freeze (`docs/twist_sector_selection_research_plan.md`
§1.2) — changing `W` still requires severing the cycle, and `W` is still
exactly deformation-invariant. Those never depended on statistics.

## 6. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal: compute the spinor, do not assert its statistics | stated |
| T2 | spectrum converges to `±\|k_m\|` | `O(1/n²)`, ratios 4.00 |
| T3 | the moding swaps | assignment reversed |
| T4 | the index | `dim ker D` = 1 twisted / 0 untwisted; `det D` = 0 / 2 |
| T5 | **the correction** | `ΔE = +π/(4C)` for both — untwisted favoured |
| T6 | two mechanisms, not one sign | energy vs index |
| T7 | assessment | 6/6 |

## 7. Open

  - This is the **free massless cycle mode**. The full throat spinor
    carries the radial/Dirac and SU(2) factors — a degeneracy, which
    scales the magnitude but changes neither the sign nor the index.
  - The index statement here is the **1D cycle** zero mode; #195 derives
    the `k`-dependent count on the monopole base. Relating the two is not
    attempted.

## 8. Cross-references

  - `docs/twist_sector_selection_research_plan.md` — the same-PR claim this
    corrects
  - `docs/throat_dirac_spinor_research_plan.md` — the SUSY square-root
    framing used to compute the spectrum
  - `docs/k1_zero_mode_index_mechanism_research_plan.md` — #195, the index
    mechanism this meets from the network side
  - `docs/pin_rp2_fermi_statistics_research_plan.md` — #170, Pin⁻ on the
    mouth; `docs/bam_scaffold_status.md` — B2, the `T² = −I` spin structure

## 9. Reproduce

```bash
python -m experiments.closure_ledger.bam_spinor_spectrum_effective_action_probe
```

Expected verdict:
`THE_BAM_SPINOR_SWAPS_THE_MODING_SO_ENERGETICS_FAVOUR_UNTWISTED_IN_BOTH_SECTORS_AND_MATTER_IS_SELECTED_BY_THE_ZERO_MODE_INDEX_NOT_BY_ENERGY`, 6/6 PASS.
