# Closed-orbit radial action probe

**Run:** 2026-05-07T06:51:04+00:00

Tests whether the radial bulk contribution becomes integer-quantized when treated as a single closed-orbit integral on the canonical Tangherlini grid.

## Exact quantum reading

  `S_radial(l, n) = (n + 1) · 2π   [hard-wall BS]`

The eigensolver imposes Dirichlet boundary conditions at both grid endpoints. For hard-wall BS, ∮p·dq = 2π·N with N = 1, 2, 3, … Mapping to the eigensolver's 0-indexed n via N = n + 1 gives S_radial = (n + 1) · 2π.

WKB → exact: the WKB closed-orbit integral approaches `(n + 1) · 2π` cleanly at high n; at low n there are O(1/n) WKB-to-exact corrections, with deviation increasing slightly with l (because the centrifugal barrier narrows the classically-allowed region at fixed n).

## Radial-action table (canonical Chebyshev N = 80)

| (l, n) | ω | S_single (π) | S_full (π) | S_full / 2π (WKB) | exact (n+1) | WKB dev | integer? |
|---|---:|---:|---:|---:|---:|---:|:---:|
| (l=1, n=0) | 1.0547 | 0.8819π | 1.7638π | 0.8819 | 1 | 0.1181 | — |
| (l=1, n=1) | 1.9744 | 1.9944π | 3.9887π | 1.9944 | 2 | 0.0056 | **✓** |
| (l=1, n=2) | 2.8941 | 2.9994π | 5.9989π | 2.9994 | 3 | 0.0006 | **✓** |
| (l=1, n=3) | 3.8247 | 3.9998π | 7.9997π | 3.9998 | 4 | 0.0002 | **✓** |
| (l=3, n=0) | 1.2191 | 0.7707π | 1.5415π | 0.7707 | 1 | 0.2293 | — |
| (l=3, n=1) | 2.1412 | 1.9398π | 3.8796π | 1.9398 | 2 | 0.0602 | — |
| (l=3, n=2) | 3.0220 | 2.9934π | 5.9869π | 2.9934 | 3 | 0.0066 | **✓** |
| (l=3, n=3) | 3.9225 | 3.9985π | 7.9970π | 3.9985 | 4 | 0.0015 | **✓** |
| (l=5, n=0) | 1.3960 | 0.7605π | 1.5210π | 0.7605 | 1 | 0.2395 | — |
| (l=5, n=1) | 2.3694 | 1.7773π | 3.5546π | 1.7773 | 2 | 0.2227 | — |
| (l=5, n=2) | 3.2277 | 2.9580π | 5.9160π | 2.9580 | 3 | 0.0420 | — |
| (l=5, n=3) | 4.0861 | 3.9933π | 7.9866π | 3.9933 | 4 | 0.0067 | **✓** |

**Convergence.** Max WKB-to-exact deviation: `0.0420` at n ≥ 2 (essentially exact); `0.2395` at n = 0 (significant WKB correction).

## Species-coupling tests

Under the **exact reading**, every species → (l, n) coupling gives an integer total cycle count `N_total_exact = N_layer_1 + Σ (n_i + 1)`. The WKB approximation deviates by the per-mode WKB error; for ground-state couplings (n = 0), the deviation matches the C1 / B1 residues observed in the earlier closure-ledger sweep.

### `B1_ground (k=k_species, n=0)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 0.8819 | **3** | 0.1181 |
| muon | 3 | 4 | (3,0) | 1 | 0.7707 | **5** | 0.2293 |
| tau | 5 | 106 | (5,0) | 1 | 0.7605 | **107** | 0.2395 |

### `B1_first_excited (k=k_species, n=1)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,1) | 2 | 1.9944 | **4** | 0.0056 |
| muon | 3 | 4 | (3,1) | 2 | 1.9398 | **6** | 0.0602 |
| tau | 5 | 106 | (5,1) | 2 | 1.7773 | **108** | 0.2227 |

### `B1_second_excited (k=k_species, n=2)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,2) | 3 | 2.9994 | **5** | 0.0006 |
| muon | 3 | 4 | (3,2) | 3 | 2.9934 | **7** | 0.0066 |
| tau | 5 | 106 | (5,2) | 3 | 2.9580 | **109** | 0.0420 |

### `B1_third_excited (k=k_species, n=3)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,3) | 4 | 3.9998 | **6** | 0.0002 |
| muon | 3 | 4 | (3,3) | 4 | 3.9985 | **8** | 0.0015 |
| tau | 5 | 106 | (5,3) | 4 | 3.9933 | **110** | 0.0067 |

### `B2_radial_ladder (l=1, n=(k-1)/2)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 0.8819 | **3** | 0.1181 |
| muon | 3 | 4 | (1,1) | 2 | 1.9944 | **6** | 0.0056 |
| tau | 5 | 106 | (1,2) | 3 | 2.9994 | **109** | 0.0006 |

### `A_cumulative_odd (l=1..k, n=0)`

| species | k | N_layer_1 | coupled modes | N_radial_exact | N_radial_WKB | N_total_exact | WKB deviation |
|---|---:|---:|---|---:|---:|---:|---:|
| electron | 1 | 2 | (1,0) | 1 | 0.8819 | **3** | 0.1181 |
| muon | 3 | 4 | (1,0), (3,0) | 2 | 1.6526 | **6** | 0.3474 |
| tau | 5 | 106 | (1,0), (3,0), (5,0) | 3 | 2.4131 | **109** | 0.5869 |

## Verdict

Integer quantization not confirmed at high n; the (n + 1) reading does not hold within the tolerance threshold.

## What's next

If the exact-reading P1 holds, the next sub-probes are:

- **Identify which species → (l, n) coupling is physical.** Multiple couplings give integer N_total but predict different values (e.g. electron N_total = 3 for (1, 0), 4 for (1, 1), 5 for (1, 2), …). The natural next probe constrains which coupling reproduces an independent observable (e.g. the lepton-quark mass-ratio chain or the QCD-pinhole γ identification).
- **Verify hard-wall BS for the actual eigenproblem.** The eigensolver uses Dirichlet boundary conditions at the grid endpoints, but the physical content (bounded antipodal cavity + throat) might want a different boundary structure (e.g. throat reflection at r → R_MID with a non-trivial phase). The (n + 1) integer pattern would shift to `(n + 1/2)` or `(n + 3/2)` if the boundaries are mixed, reading a different Maslov index off the geometry.
- **Aharonov-Bohm form** (sub-target #2). The Hopf connection should give `N · π · cos(χ)` action per fibre loop, with spinor double-cover doubling N. Compute and verify the consistency with the closure-cycle reading.