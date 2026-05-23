# S³ cosmological expansion / ΔR scale-modulus probe

**Run:** 2026-05-23T16:41:04+00:00

**Question:** Is ΔR invariant under S³ cosmological expansion?

Follows the B4 audit (PR #52): BAM needs exactly one external dimensionful anchor. This probe tests whether the invariant bulk separation ΔR can supply it — i.e. whether ΔR is a proper (fixed) or comoving (co-expanding) length.

- **ΔR**: `ΔR = R_OUTER − R_INNER = 0.52·R_MID`
- **Answer**: invariant (proper length; static bound vacuole)
- **Relocated anchor**: `m_e = f_closure·ℏ/(ΔR·c), f_closure = 0.52`
- **Theorem**: scale-modulus theorem satisfied (one anchor), not evaded

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_bound_system_discrete_spectrum` | discrete spectrum → intrinsic proper scale | **PASS** |
| T2 | `T2_comoving_vs_proper_redshift` | comoving excluded (mass shift ~z vs bound 1e-05) | **PASS** |
| T3 | `T3_staticity_vacuum` | vacuum throat static (comoving needs source) | **PASS** |
| T4 | `T4_scale_separation_einstein_straus` | ΔR/R_cosmo ~ 1e-39 (decoupled) | **PASS** |
| T5 | `T5_relocated_anchor` | m_e=f·ℏ/(ΔR c) (rel err 2e-16) | **PASS** |
| T6 | `T6_falsifiable_prediction` | local ratios stable (drift 2e-11) | **PASS** |
| T7 | `T7_two_s3_distinction` | internal vs cosmo S³ ~39 orders | **PASS** |
| T8 | `T8_assessment` | ΔR invariant; anchor relocated, value underived | **PASS** |

## T1: Bound system / discrete spectrum

| l | spectrum ω(l,n) | level spacings | discrete |
|---:|---|---|:---:|
| 1 | 1.055, 1.974, 2.894, 3.825, 4.761 | 0.920, 0.920, 0.931, 0.936 | True |
| 3 | 1.219, 2.141, 3.022, 3.922, 4.839 | 0.922, 0.881, 0.900, 0.917 | True |
| 5 | 1.396, 2.369, 3.228, 4.086, 4.971 | 0.973, 0.858, 0.858, 0.885 | True |

## T2: Comoving-vs-proper redshift (observational discriminator)

| z | a=1/(1+z) | ω comoving | mass ratio (1+z) | shift comoving | shift proper | excluded |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.5 | 0.6667 | 1.5821 | 1.5000 | 0.5000 | 0.0000 | True |
| 1.0 | 0.5000 | 2.1094 | 2.0000 | 1.0000 | 0.0000 | True |
| 2.0 | 0.3333 | 3.1642 | 3.0000 | 2.0000 | 0.0000 | True |
| 3.0 | 0.2500 | 4.2189 | 4.0000 | 3.0000 | 0.0000 | True |

Observational bound on mass drift: 1e-05. Comoving excluded: True.

## T3: Staticity / vacuum

- comoving ∂_t rs = H·rs = 8.467e-31 m/s (≠ 0 → needs a stress-energy source)
- proper (vacuum) ∂_t rs = 0.0 m/s (static → fixed proper rs)
- vacuum throat is static/proper: True

## T4: Scale separation / Einstein–Straus decoupling

- ΔR (proper) = 2.008e-13 m
- R_cosmo (Hubble radius) = 1.367e+26 m
- ΔR/R_cosmo = 1.469e-39
- tidal effect ~ (ΔR/R_cosmo)² = 2.157e-78
- decoupled from Hubble flow: True

## T5: Relocated anchor

- f_closure = ΔR/R_MID = 0.52
- ΔR (proper) = 2.0080e-13 m
- m_e predicted = 9.109384e-31 kg (actual 9.109384e-31; rel err 1.9e-16)
- ℏ bridge reconstructed: rel err 0.0e+00

Honest: the same single anchor re-expressed — relocation to a geometric invariant, not a derivation of the value.

## T6: Falsifiable prediction

- lepton mass ratio ω(3,0)/ω(1,0): a=1 → 1.155828, a=2 → 1.155828 (drift 1.5e-11)
- local throat ratios (a-independent):
  - deltaR_over_R_MID = 0.52
  - R_INNER_over_R_MID = 0.74
  - R_OUTER_over_R_MID = 1.26
- cosmological ratio ΔR/R_cosmo: a=1 → 1.469e-39, a=2 → 7.343e-40 (drifts ∝ 1/a)

**Prediction:** local ratios (lepton mass ratios) constant in cosmic time; only ΔR/R_cosmo drifts. A detection of mass-ratio drift correlated with the cosmic ratio would falsify.

## T7: Two-S³ distinction

- internal/particle S³ radius ~ 3.862e-13 m
- cosmological S³ radius ~ 1.367e+26 m
- separation ~ 39 orders of magnitude
- closure ledger uses internal S³ only: True → cosmological expansion does not enter the spectrum (invariance either way)

## T8: Assessment

- ΔR invariant: True
- supplies B4 anchor as a geometric invariant: True
- derives the value: False
- remaining prize: pin ΔR to a second fixed scale (e.g. Planck via closure quantum)

## Verdict

**DELTA_R_INVARIANT.** ΔR INVARIANT under S³ cosmological expansion. The bulk separation ΔR = R_OUTER − R_INNER = 0.52·R_MID is a PROPER, cosmologically invariant length:

  (1) The throat is a BOUND SYSTEM — a discrete spectrum with O(1) level spacings (hard walls from B3/T²=−I) — so it carries an intrinsic proper scale (as an atom does), not a comoving one.
  (2) COMOVING co-expansion is observationally EXCLUDED: rs ∝ a would give ω ∝ (1+z), redshifting particle masses by ~100–300 % to z~1–3, against quasar bounds ≲10⁻⁵ (~5 orders).
  (3) A comoving throat rs(t)=a(t)rs₀ has ∂_t rs = H·rs ≠ 0, needing a stress-energy source; the vacuum Tangherlini throat is static (∂_t rs = 0) → fixed proper rs (Einstein–Straus vacuole: static interior, comoving boundary).
  (4) Scale separation is overwhelming: ΔR/R_cosmo ~ 10⁻³⁹, tidal ~ 10⁻⁷⁸ — the bound throat is decoupled from Hubble flow.

So ΔR can supply the single B4 anchor as a GEOMETRIC INVARIANT: the dimensional bridge becomes m_e = f_closure·ℏ/(ΔR·c) with f_closure = ΔR/R_MID = 0.52, making the electron mass a consequence of a fixed bulk length, and predicting that local throat ratios (lepton mass ratios included) are constant in cosmic time while only ΔR/R_cosmo(t) ∝ 1/a drifts — consistent with the observed constancy of dimensionless constants.

HONEST CAVEAT: this does NOT evade the scale-modulus theorem (PR #52). ΔR is still ONE external dimensionful input — the bridge is the same single anchor re-expressed (ΔR, R_MID, m_e differ only by the dimensionless 0.52). The gain is identifying that anchor as a cosmologically-invariant geometric length rather than a particle property; ΔR's value itself is still not derived. Pinning ΔR to a second fixed scale (e.g. a closure-quantum relation to the Planck length) is the remaining prize.

## What this leaves open

- **ΔR's value.** Still one external dimensionful number. Pinning it to a second fixed scale (e.g. a closure-quantum relation to the Planck length) would be a genuine derivation rather than a relocation.
- **The internal-vs-cosmological S³ link.** Whether BAM's internal S³ is dynamically tied to the cosmological S³ is an assumption stated, not derived; invariance holds under either reading, but a derived link would sharpen the picture.
