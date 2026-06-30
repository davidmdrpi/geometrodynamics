# Multi-throat mechanics & the exchange kernel from the GR ψ–Φ–q soliton (PR #185)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The exchange kernel of two GR throat-solitons

The arc built a single self-gravitating ψ–Φ–q throat-soliton (#176–#180) and
showed its charges are protected topological invariants (#181–#184). This
probe takes **two** of those solitons — multi-throat mechanics — and derives
the **exchange kernel** from GR: the amplitude/operator for swapping two
identical throats. It factorizes into a GR-geometric **spatial** part and a
**topological sign**, both derived (no postulated statistics):

```
K_exchange(R) = (−1)  ×  K(R) .
                 sign     overlap
```

## The two-throat configuration and the exchange operator

Two identical throats live on the two-mouth configuration space (its `π₁` was
computed in #171). The exchange operator `P` swaps them; since swapping twice
is the identity, `P² = 1`, so `P` has eigenvalues `±1` — the symmetric (boson)
and antisymmetric (fermion) sectors. **Which sector a pair of throats occupies
is not a free choice**: it is fixed by the GR geometry of the swap.

## The spatial exchange kernel K(R) (from the actual #180 soliton)

`K(R)` is the overlap of two #180 self-gravitating throat-solitons separated
by `R`:

```
K(R) = ∫ φ(|x|) φ(|x − Rẑ|) d³x = 2π ∫ r² dr ∫ du φ(r) φ(√(r²+R²−2rRu)) ,
```

with `φ(r)` the normalized single-throat orbital (the #180 soliton). Measured:

| R | 0.0 | 1.0 | 2.0 | 3.0 | 4.0 | 6.0 |
|---|---:|---:|---:|---:|---:|---:|
| `K̂(R)` | 1.00 | 0.79 | 0.41 | 0.15 | 0.045 | 0.003 |

It decays smoothly and monotonically over the throat-soliton size
(RMS ≈ 1.27) — a GR-geometric exchange **range**, not a postulated form
factor. Two throats exchange strongly only when they overlap; far apart, the
exchange kernel is exponentially small (they are effectively distinguishable).

## The exchange sign −1 (Pin⁻ geon statistics)

The sign is `−1` — fermionic — derived from GR. The large diffeomorphism of
the spatial 3-geometry that **swaps** two throats is homotopic to a **2π
rotation** of one throat (the Friedman–Sorkin / Dowker–Sorkin spin-statistics
theorem for geons: *exchange ≃ rotation*). On the non-orientable throat (the
antipodal RP² closure, Pin⁻) a 2π rotation is the monodromy `T² = −I`, so
`½ tr T² = −1` ⟹ the exchange phase is `−1`. The GR geometry therefore
**selects** the antisymmetric (Fermi) eigenvalue of `P` — a boson would
require the orientable (`T² = +I`) closure the throat does not have
(#170/#174/#183). The full kernel is `K_exchange(R) = (−1)·K(R)`.

## Pauli exclusion at coincidence

With two throats in orbitals `φ_a, φ_b`, the two-body spatial state is
`Ψ∓(z₁,z₂) = φ_a(z₁)φ_b(z₂) ∓ φ_a(z₂)φ_b(z₁)`. At coincidence (`z₁ = z₂`) the
antisymmetric (fermion) state vanishes **identically** —
`max|Ψ₋(z,z)| = 0` (machine zero, the determinant of two equal rows) — so two
identical throats **cannot occupy the same state**. The symmetric (boson)
state does not vanish (in fact enhanced — bunching). The `−1` the geometry
selects is exactly the Pauli exclusion of two throats.

## The exchange hole and the Fermi pressure

The exchange term `∝ K(R)²` carves an **exchange hole** (suppressed
coincidence): `0.89, 0.63, 0.17, 0.002` at `R = 0.5, 1, 2, 4` — the soliton
overlap squared, so the hole has a GR range = the soliton size.
Macroscopically, the exclusion (one throat per state) fills a degenerate
**Fermi tower**: with the 3D density of states `g(E) ∝ √E`,

```
N(E_F) ∝ E_F^{3/2} ,   E(E_F) ∝ E_F^{5/2}   ⟹   E ∝ N^{5/3}
⟹  P = (2/3)(E/V) ∝ n^{5/3} ,   Γ = 5/3 ,
```

exactly the Fermi EoS **measured** in #172. The GR-derived exchange kernel is
the microscopic origin of the Fermi pressure of throat matter.

## Honest scope

- The exchange **sign** is exact / topological — the Pin⁻ geon statistics, a
  representation of the GR large-diffeomorphism (mapping-class) group, with
  exchange ≃ 2π rotation = `T² = −I`.
- The **spatial** kernel is the soliton-**overlap** model: two rigid copies of
  the #180 radial throat-soliton at separation `R` as the single-particle
  orbitals — a real GR object, but the full two-body GR problem (the actual
  two-throat metric, the gravitational direct/Hartree interaction alongside the
  exchange, and the dynamical swap with back-reaction) is a follow-up.
- The Fermi EoS index `5/3` is the standard degenerate-gas result, here
  attributed to the GR-derived exchange kernel (and matching the independently
  measured #172 EoS).
- Still weak-field / semi-dynamical (the #180 soliton).

## Reproduce

```bash
python -m experiments.closure_ledger.multi_throat_exchange_kernel_probe
# Verdict: MULTI_THROAT_EXCHANGE_KERNEL_FROM_GR_OVERLAP_TIMES_PIN_MINUS_SIGN_ANTISYMMETRIC_FERMI_PRESSURE
```
