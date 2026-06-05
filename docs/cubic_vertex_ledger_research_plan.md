# Cubic vertex ledger for the antipodal matter kernel (PR #137)

PR #136 computed the one-loop self-energy of the antipodal matter kernel using a
cubic vertex `g_{knm} = ∫ ψ_k ψ_n ψ_m`, and flagged it as **modelled** — a
generic triple overlap, not derived from the S_BAM measure. This PR is the
**ledger** for that vertex: it separates what the antipodal structure **derives**
about the cubic vertex (its selection rules, its geometric shape, its symmetry)
from what stays an **input** (the overall coupling strength, and whether S_BAM
generates the cubic term at all).

## The vertex factorises

A cubic vertex of three matter modes (each a radial profile times an S³ harmonic
`Y_l`) factorises into an angular integral, a radial overlap, and a coupling:

```
V = λ · [ ∫_{S³} Y_{l1} Y_{l2} Y_{l3} dΩ ] · [ ∫ ψ_k ψ_n ψ_m dr* ] .
```

## The angular selection rule is DERIVED (antipodal parity + SO(4))

The S³ harmonic triple integral is nonzero only if

  - **(a) `l1 + l2 + l3` is EVEN** — the **antipodal parity** rule: under the
    inversion `x → −x` (the throat ↔ antithroat C-swap, #63), `Y_l → (−1)^l Y_l`,
    so the integrand over the inversion-symmetric S³ must be even,
    `(−1)^{l1+l2+l3} = +1`. This is the **same Z₂** that fixed the antipodal
    boundary condition (#129), graded the kernel (#135), and sorted the flavor
    sectors (#134);
  - **(b) the triangle inequality** `|l1−l2| ≤ l3 ≤ l1+l2` — SO(4)
    angular-momentum addition on S³.

Verified exactly via the S³ monomial integral:

| (l₁,l₂,l₃) | Σl even? | triangle? | allowed? | ∫YYY | nonzero? |
|---|:---:|:---:|:---:|---:|:---:|
| (0,0,0) | ✓ | ✓ | ✓ | 1.0 | ✓ |
| (1,1,0) | ✓ | ✓ | ✓ | 0.25 | ✓ |
| (1,1,2) | ✓ | ✓ | ✓ | 0.0417 | ✓ |
| (2,2,2) | ✓ | ✓ | ✓ | 0.0052 | ✓ |
| (1,1,1) | ✗ | ✓ | ✗ | 0 | ✗ |
| (3,1,1) | ✗ | ✗ | ✗ | 0 | ✗ |
| (1,1,4) | ✓ | ✗ | ✗ | 0 | ✗ |

So the antipodal structure **derives which cubic vertices exist**: odd-Σl
vertices are forbidden by the Z₂ parity, triangle-violating ones by angular
momentum.

## The radial overlap is geometric (DERIVED shape, symmetric, real)

The radial factor `∫ ψ_k ψ_n ψ_m dr*` is fixed by the antipodal cavity modes
(#116/#135/#136): a definite geometric number for each `(k,n,m)`, **totally
symmetric** in its indices (Bose symmetry, verified to ~1e-14) and **real**
(Hermitian theory). The vertex **shape** is derived; only its overall scale is
free.

## What stays INPUT

The overall coupling strength `λ` (dimensionless) is **not** derived — it is the
input the #136 self-energy set to 1. And whether the S_BAM measure (#115–#122)
actually generates a cubic term (the existence/normalisation of the vertex) is
modelled, not derived. The ledger isolates these as the residual.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | cubic vertex ledger for the antipodal matter kernel (#136 vertex) |
| T2 | factorisation | `V = λ · (angular ∫YYY) · (radial ∫ψψψ)` |
| T3 | angular rule DERIVED | nonzero ⟺ Σl even (Z₂) ∧ triangle (SO(4)); verified exactly |
| T4 | parity = arc Z₂ | the Σl-even rule is the `(−1)^l` of #129/#134/#135 (#63 C-swap) |
| T5 | radial geometric | `∫ψ_kψ_nψ_m` definite, totally symmetric, real (#116/#136) |
| T6 | ledger | derived (rule, shape, symmetry, reality) vs input (coupling, S_BAM gen) |
| T7 | scope | vertex structure audited; coupling / higher vertices open |
| T8 | assessment | `CUBIC_VERTEX_LEDGER_ANTIPODAL_PARITY_SELECTION_GEOMETRIC_SHAPE_COUPLING_INPUT` |

## Established and open

  - **Established (BAM-native):** the antipodal structure derives the cubic
    vertex's angular **selection rule** (`Σl` even — the antipodal Z₂ — plus the
    SO(4) triangle) and its geometric, symmetric, real radial **shape**; the
    parity rule is the same `(−1)^l` that runs through #129/#134/#135. The
    vertex structure is BAM-native.

  - **Does not / open:** the overall coupling `λ` (dimensionless) is **input**,
    not derived from S_BAM; whether the S_BAM measure generates the cubic term
    (existence/normalisation) is modelled; the quartic / higher vertices are not
    addressed. The bulk-scale (#133) and flavor (#134) residuals stand.

## Cross-references

  - `docs/antipodal_kernel_one_loop_self_energy_research_plan.md` — #136, the
    self-energy that modelled this vertex.
  - `docs/null_throat_boundary_conditions_research_plan.md` /
    `docs/flavor_hierarchy_log_bounce_audit_research_plan.md` /
    `docs/antipodal_horizon_exchange_kernel_research_plan.md` — #129/#134/#135,
    the antipodal Z₂ `(−1)^l` the selection rule reproduces.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116, the
    cavity modes setting the geometric radial overlap.
  - `docs/charge_conjugation_swap_research_plan.md` — #63, the C-swap inversion.

## Run

```
python -m experiments.closure_ledger.cubic_vertex_ledger_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_cubic_vertex_ledger_probe/`.
Expected verdict:
`CUBIC_VERTEX_LEDGER_ANTIPODAL_PARITY_SELECTION_GEOMETRIC_SHAPE_COUPLING_INPUT`, 8/8 PASS.
