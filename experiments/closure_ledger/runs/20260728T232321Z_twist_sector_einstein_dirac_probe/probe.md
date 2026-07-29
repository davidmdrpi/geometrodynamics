# Probe D — semiclassical back-reaction (PR #230)

_Run 2026-07-28T23:23:21.315650+00:00 · 6/6 PASS_

Precondition (Probe C): both sectors admissible — **met**.

## The comparison

| quantity | differs? | why |
|---|---|---|
| radial semiclassical equations | no | identical source <T_AB>_ren (T3), identical background |
| regularity of the solution | no | same equations, same boundary data |
| on-shell action MODULUS | no | |det D| identical (Probe B T5, zeta'(0) = 0.2189594807) |
| on-shell action PHASE | **YES** | arg det D = -/+ pi/8 from eta = +/-1/4 (Probe B T6) |
| negative-mode spectrum / stability | no | fluctuation operator built from the same real system |

`ΔE_vac = 0` exactly (`E_vac = 0.0088541667`); `Δ⟨T_AB⟩ = 0` pointwise by homogeneity.

## Not done

  - the exterior used is the ROUND RP^3; BAM's actual exterior is the Tangherlini throat geometry, which is NOT homogeneous
  - on a non-homogeneous background, equal integrated energy does NOT force equal pointwise <T_AB> -- the T3 argument fails there
  - a genuine solve would need the mode functions on the actual radial background, the renormalized coincidence limit of the Green function DIFFERENCE between the two spin structures, and then the coupled radial system; none is attempted
  - the fermion is taken massless and free; a mass or a coupling would break the exact mode-by-mode |spec| coincidence that T2 rests on

## Verdict

**BACK_REACTION_IS_EXACTLY_DEGENERATE_BETWEEN_THE_TWIST_SECTORS_SO_GRAVITY_DOES_NOT_SELECT_EITHER_AND_ONLY_THE_ETA_PHASE_DISTINGUISHES_THEM**
