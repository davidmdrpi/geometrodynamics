# PMNS angle extraction from mouth-localized cross-channel overlaps (PR #153)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The mixing
> matrix is assembled from the classical cavity's overtone eigenvectors and
> mouth overlaps.

PRs #151/#152 derived the neutrino-side structure (the channel-dominant
anarchic seesaw) and localized the cross-channel conversion vertex at the
cavity mouths. This PR assembles the PMNS matrix itself:

    U_PMNS = R₂₃(φ_ℓ) · O_geom · U_ν

- **U_ν** — eigenvectors of the derived channel-dominant complex anarchic
  ensemble (#151/#152; 3000 seeded draws);
- **O_geom** — the orthogonal polar factor of the winding(k = 1,3,5) ↔
  overtone(n = 0,1,2) mouth overlaps, computed near-diagonal (misalignments
  ~8–13°);
- **φ_ℓ** — one charged-side μ–τ rotation: the single O(1) rotation the
  winding hierarchy permits (m_μ/m_τ = 0.060, ×12 less hierarchical than
  m_e/m_μ = 0.0048).

## The two failure modes bracket the structure

| O structure | sin²θ₁₂ pct | sin²θ₂₃ pct | sin²θ₁₃ pct |
|---|---|---|---|
| near-diagonal (O_geom) | 62% ✓ | **98% ✗ (too small)** | 27% ✓ |
| fully anarchic (Haar) | ~20% | ~16% | **≤7% ✗ (too large)** |

The observed PMNS sits between the limits — the data select a specific
intermediate charged-side structure.

## The exact resolution: one hierarchy-permitted rotation

A left μ–τ rotation leaves sin²θ₁₂ and sin²θ₁₃ **exactly invariant**
(machine zero — it never touches the e-row) and moves only sin²θ₂₃. So the
data demand exactly **one** charged-side rotation — and it is exactly the one
the winding hierarchy permits: an O(m_τ) off-diagonal gives an O(1) left μ–τ
rotation while the steeper e-hierarchy suppresses the 1–2/1–3 rotations the
data happen not to need. **The e-row is hierarchy-protected — the BAM reason
θ₁₃ is small while θ₂₃ is large.** The natural window is broad
(φ_ℓ ∈ ~[25°, 65°], ~45% of a uniform O(1) draw): no fine-tuning.

## The assembled PMNS (φ_ℓ = 45°)

| angle | ensemble median | observed | percentile |
|---|---|---|---|
| sin²θ₁₂ | 0.200 | 0.304 | **62%** |
| sin²θ₂₃ | 0.407 | 0.450 | **56%** |
| sin²θ₁₃ | 0.050 | 0.0224 | **27%** |

All three natural — the full observed point anarchy-typical. **CP generic**:
median |J| = 0.015 (data |J|_max ≈ 0.033), P(|J| > 0.01) = 61% — the
README "generic CP" claim quantified at the PMNS level.

## Predictions

- The e-row protection keeps θ₁₃ not-too-small (P(sin²θ₁₃ < 0.002) small —
  the historic anarchy success preserved);
- no strong θ₂₃ octant preference;
- the mass-side m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV prediction (#151/#152) rides
  the same ensemble, unchanged.

## Ledger and scope

- **Derived/structural:** the ν-side anarchy; the exact R₂₃ invariance (why
  one rotation suffices); the hierarchy permission (why μ–τ and only μ–τ);
  the e-row protection.
- **Statistical:** the specific angle values (the anarchic draw — the
  localized flavor residual).
- **Modelled:** the winding-profile shape in O_geom (results need only its
  near-diagonality); φ_ℓ as an O(1) draw on the permitted block (broad
  window — not a tuned knob).
- **Open:** deriving the charged-side matrix; Majorana phases / m_ββ
  refinement; the CKM intra-channel analogue (#91). No new input; the #150
  budget unchanged.

## Reproduce

```bash
python -m experiments.closure_ledger.pmns_angle_extraction_probe
# Verdict: PMNS_NATURAL_ONE_HIERARCHY_PERMITTED_CHARGED_ROTATION_E_ROW_PROTECTED
```
