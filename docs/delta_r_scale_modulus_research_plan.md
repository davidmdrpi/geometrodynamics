# S³ cosmological expansion / ΔR scale-modulus probe

Follows the B4 audit (PR #52), which proved BAM cannot derive an
absolute unit from scale-free topology alone: exactly **one** external
dimensionful anchor is mathematically required. This probe asks whether
that anchor can be supplied **geometrically** by the **invariant bulk
separation ΔR** — and tests the prerequisite: *is ΔR invariant under
S³ cosmological expansion?*

## The question

Define the bulk separation between the two Tangherlini walls,

```
ΔR ≡ R_OUTER − R_INNER = 0.52 · R_MID      (R_INNER=0.74, R_OUTER=1.26)
```

For ΔR to serve as nature's one dimensionful input, it must be a
**proper** length — fixed in cosmic time — not a **comoving** length
that co-expands with the spatial S³ (scale factor `a(t)`,
`R_cosmo(t) = a(t)·R_0`). The B4 scale-modulus theorem says the closure
ledger cannot decide proper-vs-comoving on its own (both preserve every
dimensionless output). The decider is **dynamics**: is the throat a
locally-bound static structure, or is it dragged by the Hubble flow?

## The argument: ΔR is a proper (invariant) length

  1. **The throat is a bound system.** The Tangherlini radial cavity has
     a discrete spectrum (hard walls from `T²=−I`, B3). A bound system
     with a discrete spectrum has an intrinsic *proper* length scale —
     exactly as an atom (Bohr radius `a₀ = ℏ/(m_e c α)`, set by local
     EM) does not expand with the universe.

  2. **Comoving co-expansion is observationally excluded.** If the
     throat co-expanded (`rs ∝ a`), the radial spectrum would redshift,
     `ω ∝ 1/a = (1+z)`, so particle masses would scale as `(1+z)` — a
     ~100–300 % shift to `z ~ 1–3`. Quasar atomic spectra bound any such
     drift to `≲ 10⁻⁵`. Comoving is falsified by ~5 orders of magnitude.

  3. **A comoving throat is not a vacuum solution.** `rs(t)=a(t)·rs₀`
     gives `∂_t rs = H·rs ≠ 0` — a time-dependent metric requiring a
     stress-energy source. The Tangherlini throat is vacuum + a
     dimensionless topological BC, with no such source: it is static
     (`∂_t rs = 0`), i.e. an Einstein–Straus vacuole — fixed proper `rs`,
     comoving matching boundary.

  4. **Scale separation is overwhelming.** `ΔR ~ λ_C ~ 2×10⁻¹³ m`,
     `R_cosmo ~ R_H ~ 1.4×10²⁶ m`; `ΔR/R_cosmo ~ 10⁻³⁹`. The tidal
     effect of expansion on the bound throat scales as `(ΔR/R_cosmo)² ~
     10⁻⁷⁸`. The throat is decoupled from Hubble flow to absurd
     precision.

So **ΔR is invariant** under S³ cosmological expansion — a proper
constant of nature.

## What this buys (and what it does not)

**Buys.** ΔR is a legitimate geometric anchor: a cosmologically-fixed
bulk *length*, rather than a particle property. The dimensional bridge
becomes

```
m_e = f_closure · ℏ / (ΔR · c) ,     f_closure = ΔR/R_MID = 0.52 ,
```

so the electron mass is a *consequence* of a fixed bulk separation. It
also yields a falsifiable prediction: local throat ratios
(`ΔR/R_MID`, `R_INNER/R_MID`, lepton mass ratios) are `a`-independent
→ **constant in cosmic time** (consistent with the observed constancy of
dimensionless constants); only a quantity coupling to the drifting ratio
`ΔR/R_cosmo(t) ∝ 1/a` would vary.

**Does not.** This does **not** evade the scale-modulus theorem. ΔR is
still *one* external dimensionful input — the bridge `m_e = f·ℏ/(ΔR c)`
is the same single anchor re-expressed (`ΔR`, `R_MID`, and `m_e` differ
only by the dimensionless `0.52`). The gain is *identifying* that anchor
as a cosmologically-invariant geometric length; the value of ΔR itself
is still not derived. Pinning ΔR to a second fixed scale (e.g. a
closure-quantum relation to the Planck length) is the remaining prize.

## Two S³'s

The internal/particle S³ (Hopf bundle, hosts the throat, radius `~R_MID
~ λ_C`) and the cosmological S³ (spatial universe, radius `~R_H`) differ
by ~39 orders of magnitude. The closure ledger uses only the **internal**
S³, so cosmological expansion of the spatial S³ does not enter the
spectrum. If the two are the *same* manifold, the Einstein–Straus
vacuole argument (above) still gives invariance. Either way ΔR is
invariant; the probe states the conditional explicitly.

## Tests

  T1. **Bound system / discrete spectrum** — the throat cavity has a
      discrete spectrum with O(1) level spacings → an intrinsic proper
      scale (atom analogy).
  T2. **Comoving-vs-proper redshift (observational discriminator)** —
      comoving (`rs ∝ a`) gives `ω ∝ (1+z)` → masses redshift by ~z;
      proper gives constant masses. Quasar bound `≲10⁻⁵` excludes
      comoving → ΔR proper.
  T3. **Staticity / vacuum** — comoving `rs(t)=a(t)rs₀` has
      `∂_t rs = H·rs ≠ 0` (needs a source); the vacuum Tangherlini
      throat is static → fixed proper `rs`.
  T4. **Scale separation / Einstein–Straus decoupling** —
      `ΔR/R_cosmo ~ 10⁻³⁹`; tidal `~10⁻⁷⁸`. Bound throat decoupled.
  T5. **Relocated anchor** — `m_e = f_closure·ℏ/(ΔR·c)`,
      `f_closure = 0.52`; verify consistency with `ℏ = m_e·R_MID·c`.
      Honest: relocation, not derivation.
  T6. **Falsifiable prediction** — local ratios `a`-independent (stable
      constants); only `ΔR/R_cosmo(t)` drifts (`∝1/a`).
  T7. **Two-S³ distinction** — internal vs cosmological S³; conditional
      clarity; invariance holds either way.
  T8. **Assessment** — ΔR invariant (static bound vacuole; comoving
      observationally excluded); supplies the B4 anchor as a geometric
      invariant; does not derive its value (theorem intact).

## Verdict structure

  - **DELTA_R_INVARIANT** (expected): ΔR is a proper, cosmologically
    invariant length — the throat is a static bound vacuole (discrete
    spectrum + vacuum + dimensionless BC), and comoving co-expansion is
    observationally excluded. ΔR can supply the single B4 anchor as a
    geometric invariant, `m_e = f_closure·ℏ/(ΔR·c)`, with the constancy
    of local ratios as a falsifiable prediction. The scale-modulus
    theorem is satisfied (one anchor), not evaded; ΔR's value remains
    the one external number.

  - **DELTA_R_COMOVING** (not expected): the throat co-expands; ΔR is
    not an anchor. Would require the throat to be conformally dragged and
    would predict redshifting masses — observationally excluded.

  - **INDETERMINATE**: a test fails to discriminate.

## Cross-references

  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 audit /
    scale-modulus theorem (PR #52).
  - `docs/bam_scaffold_status.md` — barrier ledger; B4 the anchor.
  - `docs/hbar_origin_status.md` — dimensionless residuals closed;
    `R_MID = λ_C` at the Compton bridge.
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`, tortoise.
  - `geometrodynamics/constants.py` — `R_MID, R_INNER, R_OUTER, DELTA`.
  - `experiments/closure_ledger/delta_r_scale_modulus_probe.py` — this
    probe.
