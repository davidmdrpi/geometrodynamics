# Berger deformation audit of the R-unification assumption (PR #165)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. This is an
> **audit**, not a demonstration; **negative results are reported first**.

## What this is — and is not

Berger deformation is **not** a quantum test, **not** a throat-formation
test, **not** a wave-propagation test. It is an **R-unification test**.

BAM's unified Bohr–Sommerfeld mass operator (THESIS PR #83)

```
m²(k, n) = (k·2π / L_throat)² + ((n+1)·π / L_cavity)²,   L_throat = √(2π)/k₅
```

rides the *throat* scale (the Hopf-fiber winding) and the *cavity* scale
(the radial/base resolution) on a **single** S³ radius R — "everything
rides on one R". The Berger sphere S³_λ squashes the Hopf fiber by λ while
leaving the base S² round: it is the one geometric move that pulls those
two scales apart, so it directly tests whether the global cosmic-cavity
vacuum energy and the local throat self-energy are one dynamical object on
one R, or whether "one R" is merely scale-free bookkeeping.

## Enforcement guardrails (anti-rigging)

1. **No derived inversions** — `L_throat = √(2π)/k₅` is used as the
   THESIS-locked value; the Casimir uses the genuine SU(2) Berger spectrum
   `Δ(j,m) = 4j(j+1) + 4m²(λ⁻²−1)`; no constant is reverse-fit to the
   target, and no fitted number is relabelled as an arbitrary π-multiple
   (enforced by a source scan in T7).
2. **No hidden imports** — the Born rule and singlet state
   (`geometrodynamics.bell`) are never imported; `_forbid_quantum_inputs`
   raises an explicit exception if a step ever requests them.
3. **No false victories** — the `A/R + B·R²` well's stability is computed
   and then **explicitly discounted**: it is an artifact of a hardcoded
   potential and is not counted as evidence.

## The computation

- **Global Casimir `E_cav(λ)`** — the zeta-regularized conformal-scalar
  vacuum energy on the unit Berger sphere. Method: subtract the *analytic*
  Weyl growth `a₃(a)·(n+1)³` with `a₃ = ½[√(1+a)+asinh(√a)/√a]`
  (`a=λ⁻²−1`), fit the residual tail, and add back the zeta values
  `a₃ζ(−3)+b₁ζ(−1)+b₀ζ(0)`. **Validated**: at λ=1 it recovers the exact
  closed form `E_cav(1) = 1/240` to 1e-5, anomaly-free.
- **Local self-energy `λ_min(λ)`** — the mass-operator ground state
  `√((2π/(λ·L_throat))² + ω₀²)`; the winding term tracks the squashed
  fiber (∝ 1/λ), the radial Tangherlini floor `ω₀` is fiber-blind. It
  **moves** (×1.99 across λ ∈ [0.7, 1.4]).
- **`ρ(λ) = E_cav(λ)/λ_min(λ)`** — parameter-free, and **not flat**.

## The result (negative first): a clean failure

R-unification asserts both energies equal `(pure number)/R` on one R, so
`ρ = E_cav/E_self` must be a pure number — computed here as

```
ρ(1) = 3.31e-04          (the geometric "one-R" prediction)
```

But physically the cavity is the global **cosmic** vacuum (`~1/R_Hubble`)
and the throat is the local **particle** self-energy (`~1/λ_Compton`),
whose measured ratio is

```
ρ_measured = λ_C / R_Hubble = 2.97e-39
```

The one-R prediction overshoots the measured global/local ratio by **~35
orders of magnitude**. The cosmic cavity and the local throat **cannot
both ride on one R** — the global-Casimir vs local-self-energy mismatch is
the cosmological-constant problem, rediscovered as the failure of the
geometric shorthand.

## What survives (not a victory)

Within geometric units the structure is internally consistent: `ρ(λ)` is
parameter-free, `E_cav(1)` is validated against `1/240R`, and `λ_min(λ)`
moves dynamically. But `ρ(λ)` is **not flat** — the cavity and the throat
respond differently to the *same* deformation — so they are not one
dynamical object even in shape. R-unification is therefore a valid
**scale-free bookkeeping** device (consistent with the B4 single-anchor
audit) but a **failed physical single-R identification**. The boundary is
mapped: the shorthand breaks at the global-Casimir vs local-self-energy
line, by tens of orders of magnitude, under un-dialed conditions.

## Reproduce

```bash
python -m experiments.closure_ledger.berger_r_unification_audit_probe
# Verdict: R_UNIFICATION_BREAKS_GLOBAL_CASIMIR_LOCAL_SELFENERGY_MISMATCH
```
