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

**It can be had. It is unique. It works.** And the two costs this document
first reported were both **overstated** — review caught them, recomputing
corrected them, and the corrections are §4 and §5 below rather than a
footnote.

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

## 4. Cost 1 — charge: **retracted**

The first version of this document reported

```
committed H (handle cut):   ‖[H, K]‖  = 6.6e-16
the minimal mixing term:    ‖[V₂, K]‖ = 1.414
```

and concluded that the interaction is *"a charge-non-conserving mouth
term"*. **That overstates what the number shows.** A *prescribed external
potential* breaks the winding of the **throat subsystem treated as
closed**; it does not show that **total** charge must be abandoned — any
more than an external magnetic field means angular momentum is not
conserved. The compensating charge sits in whatever sources the potential.

Tested directly. Extend the model with a **carrier** mode holding 2 units
of winding per quantum, and let the mixing move winding between field and
carrier. Then `K_total = k_field + 2·n_carrier` is conserved **exactly**:

| truncation `n_max` | `‖[H_ext, K_total]‖` |
|---:|---:|
| 6 | **0.0e+00** |
| 12 | **0.0e+00** |
| 20 | **0.0e+00** |

and in the large-amplitude (mean-field) limit the induced field operator
reproduces the prescribed second harmonic with the **same coefficient
0.408248** at every `n̄` tested.

> **So the prescribed term is the mean-field limit of a charge-conserving
> interaction.** The real cost is a *structural requirement* — a winding-2
> carrier must exist at the mouth and be in a large-amplitude state — not
> the loss of a conservation law.

## 5. Cost 2 — the Bell window: **retracted and replaced**

The first version combined an **ideal projected-qubit CHSH** (computed
with no leakage in it at all) with a **separate leakage proxy**, fed the
two into `S > 4/η − 2`, and reported *"a narrow loophole-free window
peaking at S = 2.13"*. **Those are two different experiments, and the
combination is not a Bell test.**

Computed properly: each mouth has **three outcomes** — a click in `k=+1`,
a click in `k=−1`, or **no click** when the excitation has left the
encoded sector. The full winding space is evolved so the leakage is in the
dynamics rather than bolted on; no-click is assigned to a fixed outcome
(both assignments tested, and they agree exactly); probabilities are
verified to normalize to **1.000000000000**.

| `\|t\|` span | `S` **operational** | `S` post-selected | `P(both click)` |
|---:|---:|---:|---:|
| 0.2 | **2.032289** | 2.035101 | 0.986741 |
| 0.4 | **2.117169** | 2.135025 | 0.947838 |
| 0.6 | **2.222481** | 2.284225 | 0.885839 |
| 0.8 | **2.306045** | 2.458690 | 0.804739 |
| **1.0** | **2.330905** | 2.628578 | 0.709620 |
| 1.2 | **2.330891** | 2.761508 | 0.606212 |
| 1.5 | **2.328751** | 2.825902 | 0.448355 |
| 1.9 | **2.326141** | 2.828028 | 0.259820 |

> **The detection-loophole-free violation is real and broad** — above the
> local bound at *every* span tested — **and stronger than first
> reported**: 2.3309, not a marginal 2.13. The first version's "narrow
> window" was an artifact of combining two different calculations.

What survives is only a **ceiling**: leakage caps the operational
violation at **2.3309** against Tsirelson **2.8284**, conceding roughly
half the available margin. That is the one cost of the two that review
left standing.

## 6. What this settles

| claim | status |
|---|---|
| #238's gap is a no-go | **no** — unique generator, restores Tsirelson |
| *this document's first claim:* "a charge-non-conserving mouth term" | **retracted** — total charge conserved exactly with a winding-2 carrier |
| *this document's first claim:* "a narrow window peaking at S = 2.13" | **retracted** — the operational value is **2.3309**, and it violates everywhere tested |
| leakage caps the violation below Tsirelson | **stands** — 2.3309 vs 2.8284, about half the margin |
| relocate the chain to the spin doublet | **weakened but still reasonable** |

Testing beat assuming — the interaction is real and does the whole job.
But **neither cost is what this probe first claimed**: charge conservation
is *relocated into a carrier* rather than lost, and the Bell violation is
*broad* rather than marginal. The case for moving the chain to the spin
doublet is correspondingly **weaker** than the first version argued: it
now rests on a carrier requirement and a 2.33 ceiling, not on a broken
conservation law.

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | uniqueness | only `m = 2`; 0.4082 vs exactly 0 |
| T3 | it is the missing generator | pure `σ_x`, 0.408248, others 1e-16 |
| T4 | it restores the violation | 2.000000 → **2.828427** |
| T5 | **charge: retraction** | subsystem 1.414, but `‖[H_ext,K_total]‖ = 0.0e+00` |
| T6 | **operational Bell test** | `S = 2.330905`, violates at every span |
| T7 | consequences after both corrections | 2 retracted, 1 standing, 1 weakened |
| T8 | assessment | 8/8 |

## 8. Open

  - Whether a **winding-2 carrier** at a mouth is realizable in any BAM
    geometry — that, not charge conservation, is what the construction
    actually requires.
  - Whether a detector resolving `k = ±3` raises the 2.33 ceiling by
    recovering the leaked population.
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
`THE_MINIMAL_MISSING_INTERACTION_IS_THE_SECOND_CHI_HARMONIC_AND_IT_WORKS_BUT_BOTH_COSTS_WERE_OVERSTATED_TOTAL_CHARGE_IS_CONSERVED_ONCE_THE_WINDING_TWO_CARRIER_IS_INCLUDED_AND_THE_OPERATIONAL_BELL_VALUE_IS_2_33_BROADLY_NOT_A_NARROW_2_13_WINDOW`, 8/8 PASS.

## 10. Cross-references

  - `docs/detector_algebra_research_plan.md` — #238, whose gap this fills
    and prices
  - `docs/measurement_sector.md` — #209, the winding Stern–Gerlach used
    unchanged as the detector here
  - `docs/configuration_space_emergence.md` — #206, the pair state
