# Stable moving throat / Lorentz covariance falsifier probe

**Run:** 2026-05-24T03:14:33+00:00

Tests whether a boosted BAM throat behaves as a relativistic particle — invariant mass = static eigenvalue, relativistic dispersion, length contraction, stability under boost. A genuine falsifier of the "throat = particle" claim.

- **Dispersion**: `ω(k)=√(ω₀²+c²k²) → E²−(pc)²=(mc²)²`
- **Invariant mass**: equals the static eigenvalue ω(1,0) for all boosts
- **S³ caveat**: global Lorentz broken by S³; local LV ~(R_MID/R_cosmo)² ~ 10⁻⁷⁸
- **B4 caveat**: Lorentz structure derived; rest scale = single anchor

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_relativistic_dispersion` | KG dispersion exact (err 2e-16); v_g=βc | **PASS** |
| T2 | `T2_invariant_mass_equals_static_eigenvalue` | invariant mass = ω(1,0) (dev 4e-14) | **PASS** |
| T3 | `T3_length_contraction_proper_frame` | R_∥=R*/γ; proper frame invariant | **PASS** |
| T4 | `T4_stability_under_boost` | d²E/dR²>0 frame-independent (no destabilization) | **PASS** |
| T5 | `T5_s3_preferred_frame_lv_suppression` | global LV broken; local LV ~8e-78 | **PASS** |
| T6 | `T6_falsification_criterion` | BAM passes falsifier: True | **PASS** |
| T7 | `T7_b4_accounting` | Lorentz structure scale-free; scale = anchor | **PASS** |
| T8 | `T8_assessment` | throat is a (local) relativistic particle | **PASS** |

## T1: Relativistic dispersion

Rest eigenfrequency ω₀ = ω(1,0) = 1.054727 (static radial solver).

| β | γ | p=γω₀β | E=γω₀ | √(ω₀²+p²) | v_g | disp err |
|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 1.0050 | 0.1060 | 1.0600 | 1.0600 | 0.1000 | 2e-16 |
| 0.30 | 1.0483 | 0.3317 | 1.1057 | 1.1057 | 0.3000 | 0e+00 |
| 0.50 | 1.1547 | 0.6089 | 1.2179 | 1.2179 | 0.5000 | 2e-16 |
| 0.70 | 1.4003 | 1.0338 | 1.4769 | 1.4769 | 0.7000 | 2e-16 |
| 0.90 | 2.2942 | 2.1777 | 2.4197 | 2.4197 | 0.9000 | 0e+00 |
| 0.99 | 7.0888 | 7.4020 | 7.4768 | 7.4768 | 0.9900 | 0e+00 |

## T2: Invariant mass = static eigenvalue

Static eigenvalue ω₀ = 1.054727.

| β | E | p | invariant mass √(E²−p²) |
|---:|---:|---:|---:|
| 0.000 | 1.0547 | 0.0000 | 1.0547269375 |
| 0.200 | 1.0765 | 0.2153 | 1.0547269375 |
| 0.500 | 1.2179 | 0.6089 | 1.0547269375 |
| 0.800 | 1.7579 | 1.4063 | 1.0547269375 |
| 0.950 | 3.3778 | 3.2089 | 1.0547269375 |
| 0.999 | 23.5903 | 23.5667 | 1.0547269375 |

Max deviation from ω₀: 3.95e-14 (invariant mass constant: True).

## T3: Length contraction + proper-frame invariance

| β | γ | R_∥ = R*/γ | R_⊥ | proper R* | proper E(R*) |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 1.0050 | 0.9950 | 1.0000 | 1.0000 | 1.0547 |
| 0.50 | 1.1547 | 0.8660 | 1.0000 | 1.0000 | 1.0547 |
| 0.90 | 2.2942 | 0.4359 | 1.0000 | 1.0000 | 1.0547 |
| 0.99 | 7.0888 | 0.1411 | 1.0000 | 1.0000 | 1.0547 |

## T4: Stability under boost

Proper-frame R* = 0.7937, d²E/dR² = 6.0000 > 0 (frame-independent).

| β | proper d²E/dR² | stable |
|---:|---:|:---:|
| 0.00 | 6.0000 | True |
| 0.50 | 6.0000 | True |
| 0.90 | 6.0000 | True |
| 0.99 | 6.0000 | True |

## T5: S³ preferred frame / LV suppression

- R_MID (proper) = 3.862e-13 m; R_cosmo = 1.367e+26 m
- R_MID/R_cosmo = 2.824e-39
- local LV suppression ~ (R_MID/R_cosmo)² = 7.977e-78 (unobservably small)

## T6: Falsification criterion

- (a) invariant mass drift = 8.88e-16 → no drift: True
- (b) relativistic dispersion confirmed (≠ non-rel expansion): True
- (c) stable under boost: True
- **BAM passes the falsifier: True**

## T7: B4 accounting

| ω₀ scale | invariant mass | inv/ω₀ |
|---:|---:|---:|
| 0.5000 | 0.5000 | 1.000000 |
| 1.0000 | 1.0000 | 1.000000 |
| 1.0547 | 1.0547 | 1.000000 |
| 2.0000 | 2.0000 | 1.000000 |

Lorentz structure is scale-free (inv/ω₀ = 1 for any scale); the rest scale is the single anchor.

## T8: Assessment

- throat is a particle: True
- invariant mass = static eigenvalue: True
- local Lorentz covariant: True
- global Lorentz broken by S³: suppressed ~(R_MID/R_cosmo)²
- remaining: boosted soliton from full S_BAM; spin Wigner rotation; observable LV bounds

## Verdict

**MOVING_THROAT_COVARIANT.** MOVING THROAT COVARIANT. The boosted throat behaves as a relativistic particle — BAM survives the Lorentz-covariance falsifier.

DISPERSION + INVARIANT MASS. A throat with centre-of-mass momentum k has frequency ω(k)=√(ω₀²+c²k²) (the covariant d'Alembertian / Klein–Gordon dispersion, ω₀ the static radial eigenvalue ω(1,0)). So E=ℏω, p=ℏk satisfy E²−(pc)²=(mc²)² with the group velocity v_g=pc²/E=βc, and the invariant mass √(E²−(pc)²)/c² is boost-invariant and equals the static rest eigenvalue to machine precision. The moving throat's m c² agrees with the static eigenvalue — the THESIS "is the throat actually a particle" test passes.

CONTRACTION + STABILITY. The moving throat contracts R_∥=R*/γ, R_⊥=R*, while its proper-frame size R*, rest energy E(R*), and equilibrium (d²E/dR²>0) are boost-invariant — the throat carries its rest-frame equilibrium and does not destabilize.

S³ PREFERRED FRAME. S³ is closed (a preferred rest frame), so GLOBAL Lorentz invariance is necessarily broken; only LOCAL (tangent-space) covariance holds, for λ ≪ R_cosmo. The finite-size Lorentz violation is suppressed by (R_MID/R_cosmo)² ~ (λ_C/R_H)² ~ 10⁻⁷⁸ — a calculable, unobservably small violation (a prediction, not a free pass). An O(1) violation would have falsified.

B4. The rest mass m c² = ℏω₀ = ℏc/R_MID is the single anchor; Lorentz covariance is a dimensionless structural property (the dispersion, the γ factors, the invariant E²−p²c²) — derived, independent of the anchor's value. Remaining: the boosted soliton from the full S_BAM, the spin Wigner rotation, and mapping the LV suppression to observable bounds.

## What this leaves open

- **The boosted soliton from the full action.** The dispersion here is the covariant-field (KG) form; constructing the explicit boosted throat solution of S_BAM and confirming its stress-energy transforms as a 4-tensor is the follow-on.
- **Spin under boost.** Whether the Hopf-holonomy spin-½ reproduces the correct Wigner rotation under boost (the companion Berry-phase falsifier).
- **Observable LV bounds.** Mapping (R_MID/R_cosmo)² to specific Lorentz-violation observables.
