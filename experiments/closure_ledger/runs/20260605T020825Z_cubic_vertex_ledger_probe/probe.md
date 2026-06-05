# Cubic vertex ledger for the antipodal matter kernel (PR #137)

**Run:** 2026-06-05T02:08:25+00:00

Ledger for the cubic vertex `g_{knm} = ∫ ψ_k ψ_n ψ_m` the #136 self-energy modelled. It separates what the antipodal structure DERIVES (the angular selection rule, the geometric radial shape, the symmetry) from what stays INPUT (the coupling strength, the S_BAM generation).

- **Factorisation**: V = λ · (angular ∫YYY) · (radial ∫ψψψ)
- **Angular rule**: DERIVED: Σl even (antipodal Z₂, #63 C-swap) + triangle (SO(4))
- **Radial overlap**: DERIVED shape: geometric (#116 modes), totally symmetric, real
- **Parity Z₂**: same (−1)^l as #129 BC / #134 flavor / #135 kernel grading
- **Input**: coupling λ (dimensionless); whether S_BAM generates the cubic term
- **Open**: coupling λ not derived; quartic/higher vertices; S_BAM cubic generation

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | cubic vertex ledger for the antipodal matter kernel (#136 vertex) | **PASS** |
| T2 | `T2_vertex_factorisation` | factorisation V = λ · (angular ∫YYY) · (radial ∫ψψψ) | **PASS** |
| T3 | `T3_angular_selection_rule_derived` | angular selection rule DERIVED: Σl even (Z₂) + triangle (SO(4)) | **PASS** |
| T4 | `T4_parity_rule_is_the_arc_z2` | parity rule IS the arc Z₂ (−1)^l (#129/#134/#135, #63 C-swap) | **PASS** |
| T5 | `T5_radial_overlap_geometric_symmetric` | radial overlap geometric, totally symmetric, real (#116/#136) | **PASS** |
| T6 | `T6_cubic_vertex_ledger` | ledger: derived (rule, shape, symmetry) vs input (coupling, S_BAM gen) | **PASS** |
| T7 | `T7_scope` | scope: vertex structure audited; coupling / higher vertices open | **PASS** |
| T8 | `T8_assessment` | CUBIC_VERTEX_LEDGER_ANTIPODAL_PARITY_SELECTION_GEOMETRIC_SHAPE_COUPLING_INPUT | **PASS** |

## The angular selection rule (Σl even ∧ triangle), verified exactly

| (l₁,l₂,l₃) | Σl even? | triangle? | allowed? | ∫YYY | nonzero? |
|---|:---:|:---:|:---:|---:|:---:|
| (0,0,0) | ✓ | ✓ | ✓ | 1.0 | ✓ |
| (1,1,0) | ✓ | ✓ | ✓ | 0.25 | ✓ |
| (1,1,2) | ✓ | ✓ | ✓ | 0.041667 | ✓ |
| (2,2,0) | ✓ | ✓ | ✓ | 0.041667 | ✓ |
| (2,2,2) | ✓ | ✓ | ✓ | 0.005208 | ✓ |
| (1,1,1) | ✗ | ✓ | ✗ | 0.0 | ✗ |
| (1,0,0) | ✗ | ✗ | ✗ | 0.0 | ✗ |
| (3,1,1) | ✗ | ✗ | ✗ | 0.0 | ✗ |
| (1,1,4) | ✓ | ✗ | ✗ | 0.0 | ✗ |

Nonzero exactly when `Σl` is even (the antipodal Z₂ parity, the #63 C-swap inversion) **and** the SO(4) triangle holds — the derived cubic-vertex selection rule.

## The radial overlap is geometric and totally symmetric

| (k,n,m) | ∫ψ_kψ_nψ_m dr* | symmetry err |
|---|---:|---:|
| (0,0,0) | -0.9958 | 0.0 |
| (0,0,1) | -0.2031 | 0.0 |
| (0,1,2) | 0.2403 | 0.0 |
| (1,1,2) | 0.4813 | 0.0 |

## Verdict

**CUBIC_VERTEX_LEDGER_ANTIPODAL_PARITY_SELECTION_GEOMETRIC_SHAPE_COUPLING_INPUT.** THE ANTIPODAL STRUCTURE DERIVES THE CUBIC VERTEX'S SELECTION RULE AND GEOMETRIC SHAPE; ONLY THE COUPLING REMAINS INPUT. PR #136 modelled the cubic vertex g_{knm} = ∫ ψ_k ψ_n ψ_m for the one-loop self-energy; this ledger separates its derived structure from its input magnitude.

THE VERTEX FACTORISES. A cubic vertex of three matter modes (each a radial profile times an S³ harmonic Y_l) is V = λ · [∫_{S³} Y_{l1}Y_{l2}Y_{l3} dΩ] · [∫ ψ_k ψ_n ψ_m dr*] — an angular integral, a radial overlap, and a coupling.

THE ANGULAR SELECTION RULE IS DERIVED. ∫_{S³} Y_{l1}Y_{l2}Y_{l3} is nonzero only if (a) l1+l2+l3 is EVEN — the antipodal parity rule: under x → −x (the throat ↔ antithroat C-swap, #63), Y_l → (−1)^l Y_l, so the integrand over the inversion-symmetric S³ must be even, (−1)^{Σl} = +1 — and (b) the triangle inequality |l1−l2| ≤ l3 ≤ l1+l2 (SO(4) angular-momentum addition). Odd-Σl vertices are forbidden by the Z₂ parity, triangle-violating ones by angular momentum (verified exactly via the S³ monomial integral).

THE PARITY RULE IS THE ARC'S Z₂. The Σl-even rule is the same antipodal parity (−1)^l that fixed the boundary condition (#129), graded the kernel (#135), and sorted the flavor sectors (#134). The cubic vertex respects the same Z₂; the #136 self-energy bubble connects only even-Σl mode triples.

THE RADIAL OVERLAP IS GEOMETRIC. ∫ ψ_k ψ_n ψ_m dr* is a definite geometric number set by the antipodal cavity modes (#116/#135/#136), totally symmetric in (k,n,m) (Bose symmetry) and real (Hermitian theory). The vertex SHAPE is derived; only its overall scale is free.

WHAT STAYS INPUT. The overall coupling strength λ (dimensionless) is NOT derived — it is the input the #136 self-energy set to 1 — and whether the S_BAM measure (#115–#122) actually generates a cubic term (the existence/normalisation of the vertex) is modelled, not derived. The ledger isolates these as the residual.

SCOPE. Audits the cubic-vertex STRUCTURE: the derived angular selection rule, the geometric radial shape, the symmetry and reality. It does NOT derive the coupling λ from S_BAM, nor the quartic / higher vertices; the bulk-scale (#133) and flavor (#134) residuals stand.
