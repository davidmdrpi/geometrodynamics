# Probe L — the minimal missing interaction, tested directly (PR #239)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. Why test it before relocating

#238 showed BAM's Bell chain needs exactly one undelivered thing: an
operation or observable at a mouth that coherently connects distinct
winding sectors. With the derived setting (a fiber U(1) rotation) and the
derived detector (the `σ_z` winding Stern–Gerlach) the generated
`*`-algebra collapses to the abelian `span_C{I, σ_z}` and CHSH is
exactly 2.

The obvious next move is to relocate the chain onto the transported-frame
/ spin doublet, where the qubit is not the winding charge. **Before doing
that the missing interaction should be tested directly** — assuming it
cannot be had would repeat the error #233 made about the one-throat ring.

**It can be had. It is unique. It works. And it costs charge.**

## 1. It is unique at first order

Among mouth-localized fiber-angle potentials
`λ cos(2πmχ/N_χ + φ)`:

| harmonic `m` | `\|⟨+1\|V\|−1⟩\|` |
|---:|---:|
| 1 | 0.0000 |
| **2** | **0.4082** |
| 3 | 0.0000 |
| 4 | 0.0000 |

The reason is elementary and worth stating: the qubit is `k = ±1`, so the
coupling must carry `Δk = 2`, and the `m`-th harmonic carries `Δk = m`.

> **The minimal missing interaction is the second χ-harmonic at a mouth** —
> physically, a **fiber-angle anisotropy** of the mouth: a mouth that is
> not round in the fiber direction.

## 2. It is exactly the missing generator

Restricted to the qubit that term is **pure `σ_x`** with coefficient
**0.408248** — the coefficients on `I`, `σ_y` and `σ_z` are all zero to
1e-16. It does not commute with the derived detector
(`‖[V, σ_z]‖ = 0.4082`), so conjugating `σ_z` by settings it generates
sweeps the x–z circle instead of standing still.

That is precisely the non-abelian element #238 found to be absent, and it
is supplied by *one specific term* rather than by an open-ended search.

## 3. It works — with the derived detector

| | max CHSH |
|---|---:|
| derived `σ_z` detector, settings from the second harmonic | **2.828427** |
| derived `σ_z` detector, no mixing term (#238) | 2.000000 |

No transverse readout is assumed anywhere. The minimal interaction
converts #238's abelian `span{I, σ_z}` into a non-abelian algebra and
recovers Tsirelson. **So #238's gap is not a no-go.** The question becomes
what filling it costs.

## 4. Cost 1 — charge conservation

Winding **is** charge in this program (#42–#44, the KK gauge coupling),
and the committed dynamics conserves it exactly:

```
committed H (handle cut):   ‖[H, K]‖ = 6.6e-16
the minimal mixing term:    ‖[V₂, K]‖ = 1.414
```

> **The interaction that rescues the Bell chain on the winding carrier is
> a charge-non-conserving mouth term.**

This is not a technicality to be absorbed later. The fiber U(1) that
protects winding is *exactly* the symmetry the mixing term must break —
the same generator cannot both conserve the charge and rotate the qubit.
They cannot both hold at the same mouth.

## 5. Cost 2 — leakage, and the loophole-free window

The term also drives population out of the encoded `k = ±1` sector (into
`k = ±3`), so the settings that rotate the detector also deplete the
qubit. Treating that leakage conservatively as detection inefficiency and
applying the standard CHSH bound `S > 4/η − 2`:

| `\|t\|` span | CHSH `S` | `η` | `4/η − 2` | margin |
|---:|---:|---:|---:|---:|
| 0.05 | 2.002219 | 0.9995 | 2.0019 | **+0.000350** |
| 0.20 | 2.035079 | 0.9925 | 2.0300 | **+0.005031** |
| **0.40** | **2.134938** | **0.9704** | **2.1220** | **+0.012973** |
| 0.50 | 2.204661 | 0.9540 | 2.1927 | +0.011964 |
| 0.70 | 2.369812 | 0.9113 | 2.3892 | −0.019371 |
| 1.20 | 2.760876 | 0.7553 | 3.2958 | −0.534970 |
| 1.90 | 2.828028 | 0.4697 | 6.5162 | −3.688143 |

> **A loophole-free window exists** — so this is not a no-go either — **but
> it is narrow**, positive only for `|t| ≲ 0.6` and peaking at
> `S = 2.1349` with `η = 0.9704`. **Tsirelson saturation is not in it**: at
> the span reaching `S = 2.828` the retention has fallen to `0.470` and
> the threshold is `6.52`.

**Honest caveat.** The leaked population sits in *other winding channels*
and is in principle detectable there, so reading it as pure inefficiency
is the conservative choice, not the only one. A detector resolving
`k = ±3` would recover some of it — untested here.

## 6. What this settles

| claim | status |
|---|---|
| #238's gap is a no-go | **no** — the generator exists, is unique at first order, and restores Tsirelson |
| the Bell chain can stay on the winding carrier | **only at the price of charge** — `‖[V₂, K]‖ = 1.414` against the committed 6.6e-16 |
| granting the term, the violation is as strong as the formalism allows | **no** — loophole-free only to `S ≈ 2.13` |
| relocate the chain to the spin doublet | **recommended, and now quantitatively** |

Testing beat assuming: the interaction is real and it works. But both
costs are structural rather than incidental, and the charge one is
decisive — which is a clean argument for moving the chain to a carrier
whose qubit is *not* the charge. That is what the transported-frame / spin
doublet (#192/#197) offers, and the recommendation now carries numbers
instead of a preference.

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | uniqueness | only `m = 2`; 0.4082 vs exactly 0 |
| T3 | it is the missing generator | pure `σ_x`, 0.408248, others 1e-16 |
| T4 | it restores the violation | 2.000000 → **2.828427** |
| T5 | cost 1: charge | `‖[V₂,K]‖ = 1.414` vs 6.6e-16 |
| T6 | cost 2: leakage and the window | window exists, best margin +0.0130; Tsirelson outside |
| T7 | consequences | 1 refuted, 1 priced, 1 bounded, 1 recommended |
| T8 | assessment | 8/8 |

## 8. Open

  - Whether a detector resolving `k = ±3` recovers the leaked population
    and widens the window. The conservative reading is a choice.
  - **Whether the transported-frame / spin doublet supplies a rotation
    generator without breaking any conservation law** — the relocation
    this probe was run to inform, still untested and now the clear next
    step.
  - Whether any BAM geometry realizes a fiber-anisotropic mouth, and what
    else that would break besides charge.
  - The multipartite chain (#207/#208) under the same term.

## 9. Reproduce

```bash
python -m experiments.closure_ledger.minimal_mixing_interaction_probe
```

Expected verdict:
`THE_MINIMAL_MISSING_INTERACTION_IS_THE_SECOND_CHI_HARMONIC_AT_A_MOUTH_WHICH_IS_EXACTLY_PURE_SIGMA_X_ON_THE_QUBIT_AND_RESTORES_CHSH_TO_2_83_BUT_IT_BREAKS_WINDING_CHARGE_CONSERVATION_AND_ONLY_S_2_13_SURVIVES_THE_LEAKAGE_LOOPHOLE_FREE`, 8/8 PASS.

## 10. Cross-references

  - `docs/detector_algebra_research_plan.md` — #238, whose gap this fills
    and prices
  - `docs/measurement_sector.md` — #209, the winding Stern–Gerlach used
    unchanged as the detector here
  - `docs/configuration_space_emergence.md` — #206, the pair state
