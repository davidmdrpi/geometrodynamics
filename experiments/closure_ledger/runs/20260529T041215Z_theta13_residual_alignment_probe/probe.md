# θ13 suppression / residual alignment probe (PR #93)

**Run:** 2026-05-29T04:12:15+00:00

PR #92 found PMNS broadly anarchic (θ12, θ23 typical of a Haar-random U(3)) but θ13 = 8.6° at the 4th percentile — the one mild tension. This probe explains it: θ13 = |U_e3| is the corner (most coordinate-distant, **two-hop**) element, so a residual **nearest-neighbour** alignment of the channels (the throat↔shell coupling is local in the (k,n) lattice) suppresses it — making θ13 the smallest angle and moving the observed value from the 4th to ~21st percentile.

- **Identification**: θ13 (the two-hop corner U_e3) suppressed by a residual nearest-neighbour alignment of the closure-winding × radial-overtone channels; observed θ13 moves from 4th to ~21st percentile, and θ13 becomes the smallest angle
- **Mechanism**: corner U_e3 = two-hop amplitude (throat↔shell coupling local in (k,n))
- **Fiducial μ**: 3.0
- **Open**: exact θ13 (μ one parameter); θ13 median saturates ~14–16°; CP phases

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_theta13_tension` | θ13 at 4th percentile of pure anarchy (PR #92 tension) | **PASS** |
| T2 | `T2_theta13_is_corner_two_hop` | θ13 = corner U_e3 = two-hop (gap |g−i|=2); θ12,θ23 adjacent | **PASS** |
| T3 | `T3_residual_nearest_neighbour_alignment` | residual = nearest-neighbour coupling (throat↔shell local) | **PASS** |
| T4 | `T4_model_theta13_shifts_down` | μ≈3: θ13 median 33°→~16°, θ12/θ23 stay; θ13 smallest 0.50→0.72 | **PASS** |
| T5 | `T5_theta13_tension_resolved` | observed θ13 4th→~21st pct; θ12,θ23 stay typical | **PASS** |
| T6 | `T6_theta13_smallest_prediction` | θ13 robustly the smallest angle (observed hierarchy) | **PASS** |
| T7 | `T7_honest_scope` | mechanism robust; μ one param; θ13 median saturates ~14–16° | **PASS** |
| T8 | `T8_assessment` | THETA13_SUPPRESSED_BY_RESIDUAL_NEAREST_NEIGHBOUR_ALIGNMENT | **PASS** |

## T4–T5: Pure anarchy (μ=0) vs residual alignment (μ≈3)

| | θ12 median | θ23 median | θ13 median | θ13 smallest (frac) |
|---|---:|---:|---:|---:|
| μ=0 (anarchy) | 44.8° | 44.8° | 32.9° | 0.50 |
| μ=3 (residual) | 37.0° | 36.6° | 15.5° | 0.72 |

**Observed-angle percentiles:**

| angle | obs | pure anarchy | residual (μ≈3) |
|---|---:|---:|---:|
| θ13 | 8.6° | 4th | 21th |
| θ12 | 33.4° | 30th | 44th |
| θ23 | 49.0° | 57th | 70th |

The residual alignment moves θ13 from the 4th to ~21st percentile (tension resolved) while θ12, θ23 stay typical, and makes θ13 robustly the smallest angle.

## Verdict

**THETA13_SUPPRESSED_BY_RESIDUAL_NEAREST_NEIGHBOUR_ALIGNMENT.** θ13 IS SUPPRESSED BY A RESIDUAL NEAREST-NEIGHBOUR ALIGNMENT; THE PR #92 TENSION IS RESOLVED. PR #92 found the PMNS matrix broadly anarchic, with θ12 and θ23 typical of a Haar-random U(3), but θ13 = 8.6° sitting at the 4th percentile (Haar median 33°) — the one mild tension. This probe explains the θ13 suppression.

θ13 IS THE MOST COORDINATE-DISTANT ELEMENT. In the standard parametrisation θ13 = |U_e3| connects the electron flavour (charged lepton generation 1, the LOWEST winding k=1) to the heaviest neutrino mass eigenstate (overtone n=2, the HIGHEST) — the corner of the generation/channel lattice, a gap |g−i|=2. θ12 and θ23 are adjacent (gap 1).

RESIDUAL ALIGNMENT = NEAREST-NEIGHBOUR COUPLING. The two channels are not perfectly unrelated: the throat↔shell coupling (the PR #82 +3 shift, the PR #83 unified Bohr-Sommerfeld operator) is LOCAL in the (k,n) lattice — it links a winding to a nearby overtone. So adjacent generations still mix anarchically (a single channel-hop, unsuppressed), but reaching the g=1↔g=3 extreme requires TWO channel-hops, so the corner amplitude U_e3 is suppressed (a two-hop amplitude, as in a tight-binding model). This makes θ13 generically the SMALLEST angle and pulls its distribution below pure anarchy.

QUANTITATIVE. With a modest residual-alignment strength μ≈3 (μ=0 being pure anarchy), the θ13 distribution shifts down (median 33°→~16°) while θ12, θ23 stay large (~37°); θ13 becomes robustly the smallest angle (fraction θ13<θ12,θ23 rises from 0.50 to ~0.72); and the observed θ13=8.6° moves from the 4th percentile (pure anarchy, the tension) to the ~21st (comfortable), while θ12=33.4° (~44th) and θ23=49° (~70th) stay typical. So a modest nearest-neighbour residual alignment resolves the θ13 tension AND explains why θ13 is the smallest mixing angle — both as consequences of the corner being a two-hop amplitude.

HONEST SCOPE. ESTABLISHED (BAM-native): θ13=|U_e3| is the most coordinate-distant (two-hop) element; a residual nearest-neighbour alignment (the throat↔shell coupling is local in the (k,n) lattice) suppresses it relative to the adjacent θ12, θ23, robustly making θ13 the smallest angle and moving the observed value from the 4th to ~21st percentile — resolving the tension while keeping θ12, θ23 typical. NOT established: the exact θ13 (μ is one parameter, not derived; the θ13 median saturates at ~14–16° under this mechanism, so observed 8.6° is on the low-typical side), and the BAM origin of the nearest-neighbour locality (Bohr-Sommerfeld / +3 shift) is identified but not fully derived.

## What this leaves open

- **The exact θ13** — μ is one residual-alignment parameter (not derived); the θ13 median saturates at ~14–16° under this mechanism, so observed 8.6° is on the low-typical side.
- **The BAM origin of the nearest-neighbour locality** — identified (Bohr-Sommerfeld / +3 shift), not fully derived.
- **The CP / Majorana phases.**
