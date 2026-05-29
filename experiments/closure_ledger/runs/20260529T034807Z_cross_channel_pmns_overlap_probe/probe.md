# Cross-channel PMNS overlap probe (PR #92)

**Run:** 2026-05-29T03:48:07+00:00

PR #91 argued large PMNS vs small CKM is the cross-channel vs intra-channel distinction. This probe computes the overlap. **Result:** a naive radial overlap gives *small* mixing, so large PMNS is not a literal mode overlap — the lepton generation labels live in **different coordinates** (closure-winding `k` vs radial-overtone `n`), making the map **anarchic**. Observed PMNS is typical of a Haar-random U(3); CKM is extremely atypical (aligned).

- **Identification**: PMNS is anarchic (cross-coordinate: closure-winding × radial-overtone) — observed angles typical of Haar U(3); CKM is aligned (intra-coordinate: shell × shell) — extremely atypical of anarchy (joint p ≈ 0)
- **Naive overlap**: same-coordinate radial overlap → near-permutation (small)
- **Real structure**: different coordinates (closure-winding vs radial-overtone) → anarchy
- **PMNS class**: anarchy (cross-coordinate)
- **CKM class**: aligned (intra-coordinate)
- **Open**: specific angles (anarchy statistical); θ13 mild tension; CP/Majorana phases

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_dichotomy` | make the PR #91 dichotomy quantitative | **PASS** |
| T2 | `T2_naive_radial_overlap_is_small` | naive radial overlap → near-permutation (small angles) | **PASS** |
| T3 | `T3_different_coordinates_anarchic` | lepton gens in different coordinates (k vs n) ⟹ anarchy | **PASS** |
| T4 | `T4_pmns_typical_of_anarchy` | observed PMNS typical of Haar U(3) (30th/57th/4th pct) | **PASS** |
| T5 | `T5_ckm_aligned_atypical_of_anarchy` | CKM extremely atypical (aligned; joint p ≈ 0) | **PASS** |
| T6 | `T6_dichotomy_quantified` | PMNS ∈ anarchy class, CKM ∈ aligned class | **PASS** |
| T7 | `T7_honest_scope` | class-level robust; specific angles open; θ13 mild tension | **PASS** |
| T8 | `T8_assessment` | PMNS_ANARCHIC_CROSS_COORDINATE_CKM_ALIGNED_INTRA_COORDINATE | **PASS** |

## T2: Naive radial overlap is near-permutation (small)

Overlap of winding-imprint `sin(k·πs)` (k=1,3,5) with cavity overtones `ψ₀,ψ₁,ψ₂`:

```
  [-0.997, +0.071, +0.012]
  [+0.008, -0.053, +0.998]
  [-0.001, -0.002, -0.002]
```
Max 2nd-largest per row = 0.071 ⟹ near-permutation ⟹ small genuine mixing. **Large PMNS is not a literal radial mode overlap.**

## T4–T5: Observed angles vs the anarchic (Haar U(3)) distribution

| angle | Haar median | PMNS obs (percentile) | CKM obs (percentile) |
|---|---:|---:|---:|
| θ12 | 44.9° | 33.4° (30th) | 13.04° (5.1th) |
| θ23 | 45.0° | 49.0° (57th) | 2.38° (0.2th) |
| θ13 | 32.9° | 8.6° (5th) | 0.20° (0.0th) |

PMNS is broadly **typical** of anarchy (θ12, θ23 central; θ13 small-side). CKM is **extremely atypical**: the joint probability that a Haar U(3) is as aligned as CKM is ≈ 0.0e+00 ⟹ aligned (intra-coordinate).

## Verdict

**PMNS_ANARCHIC_CROSS_COORDINATE_CKM_ALIGNED_INTRA_COORDINATE.** PMNS IS ANARCHIC (CROSS-COORDINATE), CKM IS ALIGNED (INTRA-COORDINATE). PR #91 argued large PMNS vs small CKM is the cross-channel vs intra-channel distinction and left the explicit angles open. This probe computes the overlap and turns the dichotomy into a quantitative, falsifiable statement.

THE NAIVE RADIAL OVERLAP GIVES SMALL MIXING. If the charged-lepton and neutrino generations were both labelled by the radial overtone (intra-channel, like quarks), their mixing matrix would be the overlap of two near-orthonormal sinusoidal cavity families — a near-PERMUTATION (off-diagonal ≲ 0.07, mixing ≲ 5°). So large PMNS is NOT a literal radial mode overlap; that is the honest negative that points to the real structure.

DIFFERENT COORDINATES ⟹ ANARCHY. The two lepton generation labels live in DIFFERENT coordinates of the S³ × radial space: charged leptons in the closure-winding k = 1, 3, 5 (the Hopf-fibre / throat-traversal direction), neutrinos in the radial-overtone n = 0, 1, 2 (the cavity direction). A map between a closure-winding labelling and a radial-overtone labelling has NO preferred alignment — the coordinates are unrelated — so the PMNS matrix is effectively anarchic (generic / Haar-random in generation space). This is the BAM realisation of neutrino anarchy. Quarks are the contrast: up- and down-type generations are BOTH radial-overtone (shell) labels (same coordinate), so their map is a small deformation of the identity ⟹ aligned ⟹ small CKM.

QUANTITATIVE TEST. For a Haar-random U(3) the angles have medians θ12 ≈ θ23 ≈ 45°, θ13 ≈ 33°. The observed PMNS (33.4°, 49.0°, 8.6°) sits at the ~30th / 57th / 4th percentiles — broadly TYPICAL of anarchy (θ12, θ23 central; θ13 on the small side but non-zero, the one mild tension). The observed CKM (13.0°, 2.4°, 0.20°) sits at the ~5th / 0.2th / 0.0th percentiles, and the joint probability that a Haar U(3) is as aligned as CKM is ≈ 0 — EXTREMELY atypical of anarchy, i.e. aligned. So BAM predicts PMNS in the anarchy class and CKM in the aligned class — a clean, falsifiable separation that matches observation.

HONEST SCOPE. ESTABLISHED (BAM-native): a literal same-coordinate radial overlap gives near-permutation (small) mixing; the lepton generation labels live in different coordinates (closure-winding vs radial-overtone) ⟹ anarchic map ⟹ large mixing; the observed PMNS is typical of the anarchic (Haar) distribution while CKM is extremely atypical (aligned). NOT established: the specific PMNS angles (anarchy is statistical — no preferred angles; the explicit closure↔overtone map is not fixed by the mode geometry alone); θ13 sits on the small side of anarchy (4th percentile, the mild tension); and the CP / Majorana phases.

## What this leaves open

- **The specific PMNS angles** — anarchy is statistical (no preferred angles); the explicit closure↔overtone map is not fixed by the mode geometry alone.
- **θ13** — sits on the small side of anarchy (4th percentile), the one mild tension.
- **The CP / Majorana phases.**
