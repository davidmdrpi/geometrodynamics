# BAM Compton amplitude — structural reproduction probe

**Run:** 2026-05-16T01:58:42+00:00

Tests whether a BAM transaction amplitude built from natural ingredients reproduces structural features of the QED Compton amplitude — propagator denominators and Klein-Nishina angular dependence.

**Construction:**

```
M_BAM(θ, ω) = M_s + M_u with M_x = G_S3(ψ_x) · exp(i·φ_x) · T² where ψ_x = (x − m²)/(2m²), φ_x = π ± θ/2, T² = −I.
```

## Test summary

| # | Test | Key metric | Value | PASS? |
|---|---|---|---:|---|
| T1 | `T1_propagator_pole_reproduction` | fitted pole exponent (expected 1.0) | 1.0002 | **PASS** |
| T2 | `T2_thomson_angular` | fitted cos²θ coefficient C (KN: 0.5) | -0.0000 | **FAIL** |
| T3 | `T3_ansatz_sensitivity` | closest ansatz residual to KN | 0.5000 | **FAIL** |
| T4 | `T4_photon_factor` | baseline residual to KN | 1.0000 | **FAIL** |

## T1_propagator_pole_reproduction

Verify |M_BAM| ∝ 1/(s − m²) as s → m² (ω → 0). The fitted log-log slope should be −1; equivalently, the pole exponent is +1.

Fitted pole exponent: **1.000210** (expected 1.0). Residue: 2.7503e-01.

Sampled (ω, |M_BAM|):

```
  ω = 1.00e-02,  |M_BAM| = 1.375143e+01
  ω = 1.00e-03,  |M_BAM| = 1.378008e+02
  ω = 1.00e-04,  |M_BAM| = 1.378291e+03
  ω = 1.00e-05,  |M_BAM| = 1.378319e+04
  ω = 1.00e-06,  |M_BAM| = 1.378322e+05
```

## T2_thomson_angular

Compute |M_BAM(θ)|² in the Thomson limit (ω → 0) and fit A + B·cos θ + C·cos²θ. Klein-Nishina: (½, 0, ½).

Fitted: A = +0.5001, B = +0.5000, C = -0.0000.
Klein-Nishina: A = +0.5000, B = +0.0000, C = +0.5000.

**Dominant angular form:** `cos θ-dominant (single-cosine spin-½ structure)`

**Max pointwise residual to KN (normalised):** 1.0000

## T3_ansatz_sensitivity

Repeat T2 with three spin-phase ansätze and identify which (if any) approaches Klein-Nishina.

| ansatz | A | B | C | residual to KN |
|---|---:|---:|---:|---:|
| `spin_half` | +0.5001 | +0.5000 | -0.0000 | 1.0000 |
| `spin_one` | +0.0000 | -0.0001 | +1.0001 | 0.5000 |
| `only_u` | +0.5001 | +0.5000 | -0.0000 | 1.0000 |

**Closest ansatz:** `spin_one` with residual 0.5000.

## T4_photon_factor

Multiply M_BAM by exp(iθ) (a putative photon spin-1 contribution) and check whether Klein-Nishina structure is restored.

Baseline fit: A = +0.5001, B = +0.5000, C = -0.0000.
With photon factor: A = +0.5001, B = +0.5000, C = -0.0000.

Residual to KN: baseline = 1.0000, with photon factor = 1.0000 (improvement +0.0000).

**Photon factor restores Klein-Nishina:** NO — does not restore

## Verdict

**PARTIAL_MATCH.** PARTIAL STRUCTURAL MATCH — the BAM amplitude reproduces the propagator pole structure (T1 PASS: `G_S3(ψ) ∼ 1/ψ` with `ψ ∝ s − m²` gives the correct 1/(s − m²) leading divergence) but NOT the Klein-Nishina angular dependence (T2 FAIL). The dominant angular form is `cos θ-dominant (single-cosine spin-½ structure)` rather than KN`s cos²θ-dominant structure. The diagnostic T4 identifies whether an extra photon-side phase factor restores KN — if yes, the missing structure is the photon throat-pair representation; if no, the gap is deeper. This is the expected outcome and identifies the next structural piece.

## What this leaves open

- **Klein-Nishina from polarization sum.** The QED angular factor (1 + cos²θ) arises from the polarization sum of two transverse photon polarisations |ε·ε'|² (averaged over incoming, summed over outgoing). BAM does not yet have explicit photon polarization machinery; reproducing this sum requires either (a) explicit photon throat-pair representation (T4 diagnostic), (b) a Hopf-fibre photon model coupled to the connection (see `geometrodynamics/hopf/connection.py`), or (c) an `S³` vector Green function rather than the scalar one used here.
- **Energy-dependent recoil enhancement.** Klein-Nishina at high energy `ω ~ m_e` differs from Thomson by the recoil factor `(ω'/ω)²·[ω'/ω + ω/ω' − sin²θ]`. The BAM amplitude inherits energy-dependence only through `ψ_s, ψ_u` near the pole. Whether the full recoil enhancement is reproduced is a separate test (not run here).
- **Loop corrections.** Even if tree-level structure is reproduced, one-loop corrections (vertex, self-energy) would require BAM's bulk radial channel — the missing radial-bulk phase from the closure-ledger framework. Cross-connection between threads is open.
- **Mandelstam → ψ mapping uniqueness.** The choice `ψ = (s − m²)/(2m²)` is the simplest natural choice but not the only one. Alternative mappings (e.g. `ψ = arctan((s − m²)/Λ²)`) may give different sub-leading behaviour. The pole reproduction is robust to choice; finer structure is not.
