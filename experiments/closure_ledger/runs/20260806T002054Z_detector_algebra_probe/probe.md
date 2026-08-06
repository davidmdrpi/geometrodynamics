# Probe K — the detector algebra and the marginals (PR #238)

_Run 2026-08-06T00:20:54.338786+00:00 · 8/8 PASS_

**Q. After #237, correlator coverage is not the problem. What is?**

**A. Which effects are dialable — and the fully derived pairing does not violate Bell at all.**

## The trilemma

| setting | detector | both derived? | dialable span | max CHSH |
|---|---|---|---:|---:|
| fiber U(1) [derived] | sigma_z winding SG [derived] | **yes** | 1 | **2.000000** |
| fiber U(1) [derived] | sigma_x transverse [assumed by #236/#237] | no | 2 | 2.828427 |
| y-rotation [in #209 code] | sigma_z winding SG [derived] | no | 2 | 2.828427 |

The generated algebra is all of `M₂(C)` (dimension 4), so *the algebra is not the restriction* — the dialable set is.

## The marginals follow the same fork

| `β` | `\|cos 2β\|` | fiber×σ_z | fiber×σ_x | y-rot×σ_z |
|---:|---:|---:|---:|---:|
| 0.2618 | 0.8660 | 0.8660 | 0.0e+00 | 0.8660 |
| 0.3927 | 0.7071 | 0.7071 | 0.0e+00 | 0.7071 |
| 0.5236 | 0.5000 | 0.5000 | 0.0e+00 | 0.5000 |
| 0.7854 | 0.0000 | 0.0000 | 0.0e+00 | 0.0000 |

So #237's *identically zero marginals* holds only for the `fiber × σ_x` pairing and is **retracted as a general claim**.

## Both missing ingredients are the same thing

`σ_x` is off-diagonal between `k = ±1`, and the y-rotation mixes those sectors: both require **coherence between distinct winding sectors**. Imposing winding superselection gives max CHSH = **2.000000** — exactly the local bound.

## Verdict

**THE_ALGEBRA_IS_FULL_BUT_THE_DIALABLE_SET_IS_NOT_AND_THE_ONLY_PAIRING_OF_A_DERIVED_SETTING_WITH_A_DERIVED_DETECTOR_GIVES_CHSH_EXACTLY_TWO_SO_THE_VIOLATION_RESTS_ON_WINDING_SECTOR_COHERENCE_THAT_IS_NOWHERE_DERIVED**
