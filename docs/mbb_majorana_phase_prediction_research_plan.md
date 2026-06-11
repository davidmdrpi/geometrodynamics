# Majorana phase and m_ββ prediction from the PMNS flavor ensemble (PR #154)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The effective
> Majorana mass is read off the flavor-basis seesaw matrix the classical
> cavity structure assembles.

PR #153 assembled the PMNS matrix and left "Majorana phases / m_ββ
refinement" as its lead open item. This PR supplies it, completing the
neutrino-sector prediction card.

## The exact shortcut

`m_ββ = |(W M W^T)_ee|` — the (e,e) element of the flavor-basis Majorana
matrix, with `W = R₂₃(φ_ℓ)·O_geom` (#153) and M the derived channel-dominant
anarchic ensemble (#151/#152, 3000 seeded draws, each rescaled to the
measured m₃ = 50.14 meV). No mixing-matrix approximation; all phases
included. The Takagi factorization (`X^T M_fl X = diag(m)`) gives the
per-eigenstate terms with `Σt_i = M_fl,ee` verified to **~1e-12**.

## The exact invariance

The charged-side μ–τ rotation never touches the e-row, and `M_fl,ee` depends
only on the e-row of W: **m_ββ is exactly φ_ℓ-independent** (machine zero
across the ensemble). The one modelled O(1) angle in the #153 assembly drops
out entirely — m_ββ is *more* robust than the mixing angles, inheriting only
the geometric e-row overlap and the derived ensemble.

## The prediction

| scheme | median (meV) | 68% | 95% |
|---|---|---|---|
| self-consistent | 3.2 | [1.5, 5.9] | [0.5, 8.7] |
| data-anchored | 2.9 | [1.3, 6.0] | [0.5, 11.3] |
| conditioned (r₃₂ ∈ [4,8]) | 3.1 | — | — |

The few-meV scale is structural: the light m₁ (#151/#152; ensemble median
0.074 meV, contribution negligible) makes m_ββ a **two-term interference**
`|t₂ + t₃|` of comparable terms.

## The falsification card

- P(m_ββ < 1 meV) = 7.9% — anarchic-phase cancellation uncommon;
- P(m_ββ > 6 meV) = 14.9%; P(> 10 meV) = 0.5%; **P(> 20 meV) = 0.0%**;
- a detection above ~10 meV would **falsify** the ensemble; ton-scale
  experiments (~5–20 meV sensitivity) are predicted to see nothing or a
  floor-level signal;
- the program's earlier claim m_ββ ≲ 8 meV (zeronubb probe) is **sharpened
  into a distribution** (95% upper edge ≈ 9 meV — consistent).

## Generic Majorana phases

The relative Majorana phase Φ₂₃ = arg(t₂/t₃) is broad — P(|Φ₂₃| > π/2) =
69% — generic Majorana CP, the Majorana-sector face of the #153 generic
Dirac CP. No alignment predicted or needed.

## The neutrino-sector card, complete

Normal ordering (#113) · m₁ ≈ 0.05 meV, Σm_ν ≈ 58.8 meV (#151/#152) ·
angles anarchy-natural (#153) · CP generic, Dirac and Majorana (#153/#154) ·
m_ββ ≈ 1.5–6 meV (68%, #154).

## Ledger and scope

- **Exact/structural:** the |M_fl,ee| identity; the φ_ℓ invariance; the
  two-term interference with negligible m₁.
- **Statistical:** the distribution (the anarchic draw — the localized
  flavor residual).
- **Modelled:** the O_geom e-row (winding-profile shape) — the prediction's
  one systematic; nuclear matrix elements are an experimental overlay.
- **Open:** sharpening the e-row overlap; the CKM intra-channel analogue;
  the joint Σm_ν + m_ββ + oscillation test. No new input; the #150 budget
  unchanged.

## Reproduce

```bash
python -m experiments.closure_ledger.mbb_majorana_phase_prediction_probe
# Verdict: MBB_FEW_MEV_PHI_ELL_EXACT_INVARIANT_GENERIC_MAJORANA_PHASES
```
