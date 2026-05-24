# Stable moving throat / Lorentz covariance falsifier probe

The THESIS open problem **"Stable moving throats"**: *a boosted throat
solution must remain self-consistent; the "is the throat actually a
particle" test is whether `m c²` for a moving throat agrees with the
static eigenvalue.* This probe is a genuine **falsifier**: if BAM throats
do not transform as relativistic particles — invariant mass drifting with
velocity, a non-relativistic dispersion, or destabilization under boost —
BAM fails the "throat = particle" claim.

## The test

In BAM a particle is a throat: a static radial eigenmode `ω₀ = ω(1,0)`
(the rest energy) sitting at the self-energy equilibrium `R*` (PRs
#55–#58). For the throat to be a *particle*, a boosted throat must obey:

  - **Relativistic dispersion.** A throat with centre-of-mass momentum
    `k` has frequency `ω(k) = √(ω₀² + c²k²)` (the covariant
    d'Alembertian / Klein–Gordon dispersion), so `E = ℏω`, `p = ℏk`
    satisfy `E² − (pc)² = (m c²)²` with `m c² = ℏω₀`, and the group
    velocity `v_g = dω/dk = pc²/E = βc` is the particle velocity.

  - **Invariant mass = static eigenvalue.** `E² − (pc)² = (ℏω₀)²` for
    *all* boosts: the invariant mass equals the static rest eigenvalue —
    the THESIS test verbatim.

  - **Length contraction + proper-frame invariance.** The moving throat
    contracts, `R_∥ = R*/γ`, `R_⊥ = R*`, while its *proper-frame* size
    `R*`, equilibrium, and rest energy `E(R*)` are boost-invariant — the
    throat carries its equilibrium with it.

  - **Stability under boost.** The equilibrium and stability (`d²E/dR² >
    0`) are proper-frame conditions, hence frame-independent: a boosted
    throat does not destabilize.

## Local vs global Lorentz: the S³ preferred frame

`S³` is a closed space with a preferred rest frame, so **global** Lorentz
invariance is necessarily broken. Only **local** (tangent-space) Lorentz
covariance can hold, valid for wavelengths `λ ≪ R_cosmo`. The
finite-size-on-`S³` Lorentz violation is suppressed by

```
(R_MID / R_cosmo)²  ~  (λ_C / R_H)²  ~  10⁻⁷⁸ ,
```

the same overwhelming scale separation as the ΔR probe (#53) — a
calculable, unobservably small violation. This is a *prediction* (tiny
LV at the cosmological scale), not a free pass: an `O(1)` violation would
falsify.

## B4 accounting

The rest mass `m c² = ℏω₀ = ℏc/R_MID` is the single anchor (PRs
#55–#58). Lorentz covariance is a **dimensionless structural** property
(the dispersion relation, the `γ` factors, the invariant `E²−p²c²`) —
derived, independent of the anchor's value. Consistent with the whole
arc: structure derived, scale is the one anchor.

## Tests

  T1. **Relativistic dispersion.** `ω(k)=√(ω₀²+c²k²)` (ω₀ the radial
      eigenvalue) → `E²−(pc)²=(mc²)²`; group velocity `v_g = βc`.
  T2. **Invariant mass = static eigenvalue.** `E²−(pc)²=(ℏω₀)²` constant
      across all `β` — the "throat is a particle" test.
  T3. **Length contraction + proper-frame invariance.** `R_∥=R*/γ`,
      `R_⊥=R*`; proper `R*`, `E(R*)` boost-invariant.
  T4. **Stability under boost.** `d²E/dR²>0` is a proper-frame
      condition → frame-independent; no boost destabilization.
  T5. **S³ preferred frame / LV suppression.** Global Lorentz broken;
      local LV `~(R_MID/R_cosmo)² ~ 10⁻⁷⁸`, unobservably small.
  T6. **Falsification criterion.** `O(1)` mass drift / non-relativistic
      dispersion / destabilization would falsify; BAM passes.
  T7. **B4 accounting.** Lorentz structure derived/dimensionless; rest
      scale is the single anchor.
  T8. **Assessment.**

## Verdict structure

  - **MOVING_THROAT_COVARIANT** (expected): a boosted throat obeys the
    relativistic dispersion `E²−(pc)²=(mc²)²` with the invariant mass
    equal to the static rest eigenvalue (to machine precision), contracts
    as `R*/γ` with a boost-invariant proper frame, and remains stable.
    The throat is a relativistic particle (locally). Global Lorentz is
    broken by `S³` only at the suppressed level `(R_MID/R_cosmo)² ~
    10⁻⁷⁸`. BAM survives the falsifier.

  - **LORENTZ_FALSIFIED**: the invariant mass drifts with velocity, the
    dispersion is non-relativistic, or the throat destabilizes under
    boost — the throat is not a particle.

## What this leaves open

  - **The boosted soliton from the full action.** The dispersion here is
    the covariant-field (KG) form; constructing the explicit boosted
    throat solution of the full `S_BAM` and confirming its stress-energy
    transforms as a 4-tensor is the follow-on.
  - **Spin under boost.** Whether the Hopf-holonomy spin-½ (the Berry
    phase) reproduces the correct Wigner rotation under boost — the
    companion "Berry phases under motion" falsifier.
  - **Observable LV bounds.** Mapping `(R_MID/R_cosmo)²` to specific
    Lorentz-violation observables.

## Cross-references

  - `docs/self_consistent_throat_radius_research_plan.md` — `E(R*)`, the
    rest energy (#55).
  - `docs/pair_production_threshold_research_plan.md` — the pair
    threshold (#58).
  - `docs/delta_r_scale_modulus_research_plan.md` — `R_MID/R_cosmo` (#53).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `geometrodynamics/tangherlini/radial.py` — `solve_radial_modes`
    (the static eigenvalue ω(1,0)).
  - `experiments/closure_ledger/stable_moving_throat_probe.py` — this
    probe.
