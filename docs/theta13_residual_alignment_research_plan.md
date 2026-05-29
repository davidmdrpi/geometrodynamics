# θ13 suppression / residual alignment probe (PR #93)

Follows PR #92, which found the PMNS matrix broadly anarchic
(cross-coordinate: charged-lepton closure-winding × neutrino
radial-overtone), with θ12, θ23 typical of a Haar-random `U(3)` — but
θ13 = 8.6° sitting at the **4th percentile** (Haar median ≈ 33°), the one
mild tension. This probe explains the θ13 suppression as a **residual
alignment** between the two coordinate channels.

## θ13 is the most coordinate-distant element

θ13 = `|U_e3|` connects the electron flavour (charged lepton generation
1, the **lowest** winding `k = 1`) to the heaviest neutrino mass
eigenstate (overtone `n = 2`, the **highest**). In the
generation/channel lattice this is the **corner** — the most
coordinate-distant pair, generation gap `|g − i| = 2`. θ12 and θ23 are
adjacent (gap 1).

## Residual alignment = nearest-neighbour channel coupling

The two channels are not perfectly unrelated: the throat↔shell coupling
(the PR #82 `+3` shift, the PR #83 unified Bohr-Sommerfeld operator) is
**local** in the `(k, n)` lattice — it links a winding to a *nearby*
overtone. So adjacent generations still mix anarchically (a single
channel-hop, unsuppressed), but reaching the `g = 1 ↔ g = 3` extreme
requires **two** channel-hops, so the corner amplitude `U_e3` is
suppressed (a two-hop amplitude, as in a tight-binding model). This makes
θ13 generically the **smallest** angle and pulls its distribution below
pure anarchy.

## Model and result

Structured ("nearest-neighbour") anarchy: a complex-Gaussian matrix with
element variance 1 for `|g − i| ≤ 1` and `exp(−μ)` for the corner
`|g − i| = 2`, projected to the nearest unitary. `μ = 0` reproduces
PR #92 pure anarchy; `μ` is the residual-alignment strength. With a
modest `μ ≈ 3`:

| | θ12 median | θ23 median | θ13 median | θ13 smallest (frac) |
|---|---:|---:|---:|---:|
| μ=0 (pure anarchy) | 44.8° | 44.8° | 32.9° | 0.50 |
| μ=3 (residual) | 37.0° | 36.6° | 15.5° | 0.72 |

**Observed-angle percentiles:**

| angle | obs | pure anarchy | residual (μ≈3) |
|---|---:|---:|---:|
| θ13 | 8.6° | 4th | 21st |
| θ12 | 33.4° | 30th | 44th |
| θ23 | 49.0° | 57th | 70th |

So a modest nearest-neighbour residual alignment (i) shifts the θ13
distribution down (median 33° → ~16°) while θ12, θ23 stay large; (ii)
makes θ13 robustly the **smallest** angle (frac 0.50 → 0.72); and (iii)
moves the observed θ13 = 8.6° from the 4th to the ~21st percentile —
**resolving the tension** — while θ12, θ23 stay typical. The θ13
suppression and the observed hierarchy θ13 < θ12, θ23 are both
consequences of the corner being a two-hop amplitude.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | θ13 at the 4th percentile of pure anarchy (PR #92 tension) |
| T2 | corner element | θ13 = `U_e3` = two-hop (gap 2); θ12, θ23 adjacent (gap 1) |
| T3 | residual alignment | nearest-neighbour coupling (throat↔shell local in `(k,n)`) |
| T4 | model shift | μ≈3: θ13 median 33°→~16°, θ12/θ23 stay; θ13 smallest 0.50→0.72 |
| T5 | tension resolved | observed θ13 4th→~21st percentile; θ12, θ23 typical |
| T6 | prediction | θ13 robustly the smallest angle (observed hierarchy) |
| T7 | honest scope | mechanism robust; μ one param; θ13 median saturates ~14–16° |
| T8 | assessment | `THETA13_SUPPRESSED_BY_RESIDUAL_NEAREST_NEIGHBOUR_ALIGNMENT` |

## Established and open

  - **Established (BAM-native):** θ13 = `|U_e3|` is the most
    coordinate-distant (two-hop) element; a residual nearest-neighbour
    alignment (the throat↔shell coupling is local in the `(k,n)` lattice)
    suppresses it relative to the adjacent θ12, θ23, robustly making θ13
    the smallest angle and moving the observed value from the 4th to
    ~21st percentile — resolving the PR #92 tension while keeping θ12,
    θ23 typical.

  - **Open:** the exact θ13 (μ is one residual-alignment parameter, not
    derived; the θ13 median saturates at ~14–16° under this mechanism —
    the corner cannot be driven fully to zero by Gaussian suppression +
    unitarity — so observed 8.6° is on the low-typical side); the BAM
    origin of the nearest-neighbour locality (Bohr-Sommerfeld / `+3`
    shift) is identified, not fully derived; and the CP / Majorana phases.

## Run

```
python -m experiments.closure_ledger.theta13_residual_alignment_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_theta13_residual_alignment_probe/`.
Expected verdict: `THETA13_SUPPRESSED_BY_RESIDUAL_NEAREST_NEIGHBOUR_ALIGNMENT`, 8/8 PASS.
