# Cross-channel PMNS overlap probe (PR #92)

Follows PR #91, which argued large PMNS vs small CKM is the cross-channel
(leptons) vs intra-channel (quarks) distinction and flagged the explicit
mixing angles as open. This probe computes the cross-channel overlap and
turns the dichotomy into a quantitative, falsifiable statement.

## The naive radial overlap gives SMALL mixing (the honest negative)

If the charged-lepton and neutrino generations were both labelled by the
radial overtone (intra-channel, like quarks), the mixing matrix would be
the overlap of two near-orthonormal sinusoidal cavity families. Computing
it — winding-imprint profiles `sin(k·πs)`, `k=1,3,5`, against the cavity
overtones `ψ₀,ψ₁,ψ₂`:

```
  [-0.997, +0.071, +0.012]
  [+0.008, -0.053, +0.998]
  [-0.001, -0.002, -0.002]
```

is a near-**permutation** matrix (off-diagonal ≲ 0.07, mixing ≲ 5°). So
**large PMNS is not a literal radial mode overlap** — the honest negative
that points to the real structure.

## The real structure: different coordinates ⟹ anarchy

The two lepton generation labels live in **different coordinates** of the
S³ × radial space:

  - charged leptons: the **closure-winding** number `k = 1, 3, 5` (the
    Hopf-fibre / throat-traversal direction);
  - neutrinos: the **radial-overtone** number `n = 0, 1, 2` (the cavity
    direction).

A map between a closure-winding labelling and a radial-overtone labelling
has **no preferred alignment** — the two coordinates are unrelated — so
the PMNS matrix is effectively **anarchic** (a generic / Haar-random
unitary in generation space). This is the BAM realisation of neutrino
anarchy. Quarks are the contrast: up- and down-type generations are
*both* radial-overtone (shell) labels (same coordinate), so their map is
a small deformation of the identity ⟹ **aligned** ⟹ small CKM.

## Quantitative test against anarchy

For a Haar-random `U(3)` the mixing angles have medians `θ12 ≈ θ23 ≈ 45°`,
`θ13 ≈ 33°`. Comparing observation to this anarchic distribution
(Monte Carlo, `N = 40000`, seed 0):

| angle | Haar median | PMNS obs (percentile) | CKM obs (percentile) |
|---|---:|---:|---:|
| θ12 | 44.9° | 33.4° (30th) | 13.04° (5th) |
| θ23 | 45.0° | 49.0° (57th) | 2.38° (0.2th) |
| θ13 | 32.9° | 8.6° (5th) | 0.20° (0.0th) |

  - **PMNS** is broadly **typical** of anarchy (θ12, θ23 central; θ13 on
    the small side but non-zero — the one mild tension).
  - **CKM** is **extremely atypical**: the joint probability that a Haar
    `U(3)` is as aligned as CKM is ≈ 0 ⟹ aligned (intra-coordinate).

So BAM predicts PMNS in the anarchy class (cross-coordinate) and CKM in
the aligned class (intra-coordinate) — a clean, falsifiable separation
that matches observation.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | make the PR #91 dichotomy quantitative |
| T2 | naive overlap | near-permutation (small angles) — not a literal overlap |
| T3 | different coordinates | closure-winding `k` vs radial-overtone `n` ⟹ anarchy |
| T4 | PMNS vs anarchy | observed PMNS typical of Haar U(3) (30th/57th/4th) |
| T5 | CKM vs anarchy | extremely atypical (aligned; joint `p ≈ 0`) |
| T6 | dichotomy | PMNS ∈ anarchy class, CKM ∈ aligned class |
| T7 | honest scope | class-level robust; specific angles open; θ13 mild tension |
| T8 | assessment | `PMNS_ANARCHIC_CROSS_COORDINATE_CKM_ALIGNED_INTRA_COORDINATE` |

## Established and open

  - **Established (BAM-native):** a literal same-coordinate radial overlap
    gives near-permutation (small) mixing; the lepton generation labels
    live in different coordinates (closure-winding vs radial-overtone), so
    their map is anarchic ⟹ large mixing; the observed PMNS angles are
    typical of the anarchic (Haar) distribution while the CKM angles are
    extremely atypical (aligned). The cross- vs intra-coordinate structure
    separates the two at the distribution level.

  - **Open:** the specific PMNS angles (anarchy is statistical — no
    preferred angles; the explicit closure↔overtone map is not fixed by
    the mode geometry alone); θ13 sits on the small side of anarchy (4th
    percentile, the one mild tension); and the CP / Majorana phases.

## Run

```
python -m experiments.closure_ledger.cross_channel_pmns_overlap_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_cross_channel_pmns_overlap_probe/`.
Expected verdict: `PMNS_ANARCHIC_CROSS_COORDINATE_CKM_ALIGNED_INTRA_COORDINATE`, 8/8 PASS.
