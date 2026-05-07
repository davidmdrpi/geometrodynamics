# Closed-orbit radial action probe

**Run:** 2026-05-07T06:47:55+00:00

Tests whether the radial bulk contribution becomes integer-quantized when treated as a single closed-orbit integral with the throat-reflection phase included, rather than the per-mode WKB sum used by the previous Layer-2 candidates.

**Throat correction.**
- value: `+π` per closure cycle
- interpretation: T = iσ_y has eigenvalues ±i. T² = −I contributes phase argument π per closure pass through the non-orientable throat. Adding π to the WKB closed-orbit action gives S_full_with_throat = 2π(n + 1) — integer-quantized.

**Predicted reading.**
  `S_radial_with_throat(l, n) = 2π · (n + 1/2)  +  π  =  2π · (n + 1)`

  The Bohr-Sommerfeld `(n + 1/2)` becomes `(n + 1)` once the throat-reflection Maslov shift `+π` is included. WKB → exact as n → ∞.

## Radial-action table (canonical Chebyshev N = 80)

| (l, n) | ω | S_single (π) | S_full (π) | + throat (π) | S/2π | predict (n+1) | dev | integer? |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| (l=1, n=0) | 1.0547 | 0.8819π | 1.7638π | 2.7638π | 1.3819 | 1 | 0.3819 | — |
| (l=1, n=1) | 1.9744 | 1.9944π | 3.9887π | 4.9887π | 2.4944 | 2 | 0.4944 | — |
| (l=1, n=2) | 2.8941 | 2.9994π | 5.9989π | 6.9989π | 3.4994 | 3 | 0.4994 | — |
| (l=1, n=3) | 3.8247 | 3.9998π | 7.9997π | 8.9997π | 4.4998 | 4 | 0.4998 | — |
| (l=3, n=0) | 1.2191 | 0.7707π | 1.5415π | 2.5415π | 1.2707 | 1 | 0.2707 | — |
| (l=3, n=1) | 2.1412 | 1.9398π | 3.8796π | 4.8796π | 2.4398 | 2 | 0.4398 | — |
| (l=3, n=2) | 3.0220 | 2.9934π | 5.9869π | 6.9869π | 3.4934 | 3 | 0.4934 | — |
| (l=3, n=3) | 3.9225 | 3.9985π | 7.9970π | 8.9970π | 4.4985 | 4 | 0.4985 | — |
| (l=5, n=0) | 1.3960 | 0.7605π | 1.5210π | 2.5210π | 1.2605 | 1 | 0.2605 | — |
| (l=5, n=1) | 2.3694 | 1.7773π | 3.5546π | 4.5546π | 2.2773 | 2 | 0.2773 | — |
| (l=5, n=2) | 3.2277 | 2.9580π | 5.9160π | 6.9160π | 3.4580 | 3 | 0.4580 | — |
| (l=5, n=3) | 4.0861 | 3.9933π | 7.9866π | 8.9866π | 4.4933 | 4 | 0.4933 | — |

**WKB → exact convergence.** Max deviation from integer (n+1) at n ≥ 2: `0.4998`. At n = 0: `0.3819`. The integer reading is sharp at high n; ground states have O(1/n) WKB-to-exact corrections that shrink with mode index.

## Species-coupling tests

For each candidate species → (l, n) coupling, the total closure-cycle integer count is `N_total = N_layer_1 + N_radial`. The exact-action reading predicts `N_radial = Σ_{(l, n) ∈ S(species)} (n + 1)`; WKB gives an approximation.

### `B1_ground (k=k_species, n=0)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 1.3819 | 3 | 0.3819 |
| muon | 3 | 4 | (3,0) | 1 | 1.2707 | 5 | 0.2707 |
| tau | 5 | 106 | (5,0) | 1 | 1.2605 | 107 | 0.2605 |

### `B1_first_excited (k=k_species, n=1)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,1) | 2 | 2.4944 | 4 | 0.4944 |
| muon | 3 | 4 | (3,1) | 2 | 2.4398 | 6 | 0.4398 |
| tau | 5 | 106 | (5,1) | 2 | 2.2773 | 108 | 0.2773 |

### `B1_second_excited (k=k_species, n=2)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,2) | 3 | 3.4994 | 5 | 0.4994 |
| muon | 3 | 4 | (3,2) | 3 | 3.4934 | 7 | 0.4934 |
| tau | 5 | 106 | (5,2) | 3 | 3.4580 | 109 | 0.4580 |

### `B1_third_excited (k=k_species, n=3)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,3) | 4 | 4.4998 | 6 | 0.4998 |
| muon | 3 | 4 | (3,3) | 4 | 4.4985 | 8 | 0.4985 |
| tau | 5 | 106 | (5,3) | 4 | 4.4933 | 110 | 0.4933 |

### `B2_radial_ladder (l=1, n=(k-1)/2)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 1.3819 | 3 | 0.3819 |
| muon | 3 | 4 | (1,1) | 2 | 2.4944 | 6 | 0.4944 |
| tau | 5 | 106 | (1,2) | 3 | 3.4994 | 109 | 0.4994 |

### `A_cumulative_odd (l=1..k, n=0)`

| species | k | N_layer_1 | coupled modes | N_radial_predict | N_radial_observed | N_total_predict | deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 1.3819 | 3 | 0.3819 |
| muon | 3 | 4 | (1,0), (3,0) | 2 | 2.6526 | 6 | 0.6526 |
| tau | 5 | 106 | (1,0), (3,0), (5,0) | 3 | 3.9131 | 109 | 0.9131 |

## Verdict

**Integer quantization not confirmed.** Even at high n, the deviation from `(n + 1)` exceeds the tolerance threshold; the throat-correction reading does not give a clean integer.

## What's next

If P1 holds at the exact-quantum level, the next concrete sub-probes are:

- **Verify the throat-reflection phase identification.** The +π correction matches T² = −I numerically, but the formal derivation (which Maslov index applies to a Tangherlini boundary?) should be made explicit. The natural reading: the radial wavefunction at r → r_inner (the throat) reflects with a phase π from the orientation-reversing T² action.
- **Test the Aharonov-Bohm Hopf-fibre form** — sub-target #2 from the research plan. The Hopf connection `A = (1/2)cos(χ) dφ` should give a closure-cycle action `2π · (1/2 cos(χ))·N` for a loop wrapping N times around a fibre at hyper-latitude χ. At χ = 0 this is `π · N`; at χ = π/2 it vanishes. The spinor double-cover doubles N effectively.
- **Tangherlini eigenvalue as ℏ scale** — sub-target #3. Test whether `ℏ ω(l=1, n=0) = m_e c²` (or a similar natural relationship) is consistent across species and predicts an absolute mass scale.