# Compton vertex analytic derivation — probe plan

Follow-on to the vertex-structure scan (PR #30, merged). That probe
empirically identified clean values `(β, γ) = (−1/2, −1)` for the
Family B vertex modification that reduces the O(ω/m) KN residual by
12.6× at ε ≤ 0.2 (down to 7.6 %). This probe asks whether those
clean values can be **derived analytically** from the small-ε
expansion of the BAM amplitude, and whether the gap can be closed
exactly at O(ε) by any natural combination of BAM vertex families.

## Setup — small-ε expansion

With `c ≡ cos θ`, `ε ≡ ω/m`, and `x ≡ ω'/ω = 1/(1 + ε(1−c))`:

**Klein-Nishina at O(ε)** (deriving by Taylor expansion of
`x²(x + 1/x − sin²θ)` and normalising at θ=0):

    f_KN(ε, θ)/f_KN(ε, 0)
       ≈ (1+c²)/2  −  ε·(1−c)·c²  +  O(ε²)

**BAM baseline at O(ε)** (with the photon-structure probe's
`(1 + cos²θ)·(G_s + G_u)²`):

    f_BAM(ε, θ)/f_BAM(ε, 0)
       ≈ (1+c²)/2  +  ε·(1−c)·(1+c²)/2  +  O(ε²)

**Required correction**:

    Δ(θ) = f_KN − f_BAM
         = −ε·(1−c)·c²  −  ε·(1−c)·(1+c²)/2
         = −ε·(1−c)·(1+3c²)/2.

## Vertex modification at O(ε)

The probe in PR #30 used three families, all scaled by `ε_kin = ε`
to preserve Thomson:

**A**: `V_A = ε·ε'* + ε·α·(ε·k̂')(ε'*·k̂)`
**B**: `V_B = ε·ε'* · (1 + ε·(β·sin²θ + γ·(1−c)))`

For the polarization-summed amplitude with these modifications, the
O(ε) correction beyond baseline (after extensive but elementary
algebra) is

    Δ_mod(θ; α, β, γ)
       =  ε · {  (1+c²)·[β·(1−c²) + γ·(1−c)]
              −  α·c·(1−c²) }

Matching to the required correction `Δ_required = −ε·(1−c)·(1+3c²)/2`
gives the equation

    (1+c²)·[β·(1−c²) + γ·(1−c)] − α·c·(1−c²) = −(1−c)·(1+3c²)/2.

Factor `(1−c)` (valid for `c ≠ 1`):

    (1+c²)·[β·(1+c) + γ] − α·c·(1+c) = −(1+3c²)/2.

Expanding both sides in `c`:

|        | c⁰      | c¹       | c²       | c³ |
|--------|---------|----------|----------|-----|
| LHS    | β+γ     | β−α      | β+γ      | β−α  |
| LHS… cont'd       | (from γ + β(1+c) terms) | (cross terms in (1+c²)(1+c)β...) | ... | ... |
| RHS = −(1+3c²)/2  | −1/2 | 0 | −3/2 | 0 |

(See probe code for the full expansion.) The resulting linear system
matches each coefficient in `{1, c, c², c³}`.

## Sub-targets

### T1. Analytic expansion verification

Verify numerically that the analytic small-ε expansions of f_KN and
f_BAM are correct (computed Taylor coefficients of (f − (1+c²)/2)/ε
at small ε for each test family, compared to closed-form
predictions).

### T2. Linear system for coefficient matching

Compute the linear system relating (α, β, γ) to coefficients of
{c⁰, c¹, c², c³} in the LHS, and check whether it admits any
solution matching the RHS coefficients (−1/2, 0, −3/2, 0).

### T3. Over-determination diagnostic

If T2 produces no solution, the system is over-determined. Compute
the least-squares optimum and verify it matches the numerical scan
result from PR #30 — checking the analytic predicts the empirical
`(β, γ) ≈ (−1/2, −1)`.

### T4. Identify the missing structural piece

The over-determination implies a specific residual structure
(typically a c³-like term that no natural combination of Family A
or B can produce). Identify what natural BAM ingredient could
provide this term. Candidates:

  - Second-order Hopf-connection coupling (`hopf/connection.py`):
    the connection `A_φ = ½cos(χ)dφ` has cosine structure that
    might generate cos³θ terms through repeated coupling.
  - Per-channel kinematic weights at quadratic order (Family C
    generalization).
  - Explicit electron-spinor structure (Dirac trace algebra would
    naturally produce cos^n θ terms at all orders).

### T5. Numerical confirmation with extended family

Add a "Family E" with explicit cubic angular structure
(`δ·(ε·cos θ · (1−cos θ))` or `δ·(ε·sin²θ · cos θ)`) and verify
that a clean δ value closes the gap when combined with optimal
Family B.

## Expected outcome

**CLEAN_DERIVATION** if the analytic linear system has a solution
matching PR #30's empirical `(β, γ) = (−1/2, −1)` (with α free or
forced to a specific value), and the residual structure is fully
explained by the analytic small-ε expansion.

**OVER_DETERMINED** if the system has no exact solution — meaning
no natural polynomial vertex (Family A or B) can close the O(ε)
gap exactly, and the empirical best fit is a least-squares
compromise. This would identify the residual cubic structure as
the next BAM target.

**INCONSISTENCY** if the analytic expansion disagrees with PR #30's
numerical result — would indicate a bug in either the analytic
work or the probe construction.

The most likely outcome is OVER_DETERMINED — the small-ε system has
four constraints (c⁰, c¹, c², c³) but only three free parameters
(α, β, γ), so generically over-constrained. The probe's value is
in **deriving** this analytically (rather than the empirical
least-squares fit) and identifying the cubic residual term.

## Cross-references

- PR #25–#30 in the QFT-event-reinterpretation thread.
- `experiments/closure_ledger/compton_vertex_structure_probe.py` —
  PR #30 empirical predecessor.
- `geometrodynamics/hopf/connection.py` — candidate source of
  natural cubic angular structure.
- `experiments/closure_ledger/compton_vertex_derivation_probe.py` —
  this probe.
