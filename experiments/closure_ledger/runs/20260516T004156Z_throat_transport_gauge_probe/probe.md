# Gauge-sensitive throat-transport falsifier

**Run:** 2026-05-16T00:41:56+00:00

Tests whether the BAM combination `(A_BAM = ½cos(χ)dφ, ψ_BAM, T = iσ_y)` predicts a physically distinct phase from the standard Bloch combination `(A_Bloch = ½(cos(χ)−1)dφ, ψ_Bloch, T_Bloch)` after all allowed gauge freedoms are accounted for. The relevant gauge transformation is the SU(2) spinor-frame rotation

```
g(χ, φ) = e^{−iφ/2}                (scalar U(1), multi-valued)
ψ_BAM = g · ψ_Bloch                (spinor sign flip on 2π loop)
A_BAM − A_Bloch = +½ (per dφ)      (connection-level relation)
```

Pulling the common factor `e^{−iφ/2}` out of the BAM spinor `(cos(χ/2)e^{−iφ/2}, sin(χ/2)e^{+iφ/2})` recovers the Bloch spinor `(cos(χ/2), sin(χ/2)e^{+iφ})`. The relating factor is a scalar U(1) phase, not an SU(2) frame rotation. It is multi-valued: under `φ → φ+2π`, `g → −g`. This is the Wu-Yang transition function for the spin-½ representation of the Hopf bundle.

## Scalar tests

| # | test | metric | residual | PASS? |
|---|---|---|---:|---|
| 1 | `connection_difference_is_half` | A_BAM(χ) − A_Bloch(χ) = ½ for all χ | 5.55e-17 | **PASS** |
| 2 | `spinor_frame_transformation` | ψ_BAM(χ,φ) = g(χ,φ) · ψ_Bloch(χ,φ) with g = e^{−iφ/2} (scalar U(1) phase, multi-valued under φ → φ+2π) | 2.48e-16 | **PASS** |
| 3 | `throat_conjugation_preserves_invariants` | T_Bloch(φ_in, φ_out) preserves Tr = 0 universally; for instantaneous throat (φ_in = φ_out) also preserves det = 1 and T² = −I — same SU(2) invariants as iσ_y | 4.44e-16 | **PASS** |
| 4 | `relative_phase_gauge_invariance` | Δγ(θ_a, θ_b) for the detector-holonomy phase agrees between BAM and Bloch gauges (the offset cancels) | 8.88e-16 | **PASS** |
| 5 | `single_valued_gauge_cannot_remove_pi` | Any single-valued (periodic) gauge function Λ(φ) has ∮dΛ = 0 around a 2π loop | 1.96e-15 | **PASS** |

## Closed-loop holonomy tests

Each loop is checked against the BAM↔Bloch operator relation `U_loop^BAM = (−1)^n_winding · U_loop^Bloch` and against the SO(3) invariant `|Tr U|²` (which must agree regardless of spinor sign).

| loop | crosses throat | n_winding | (−1)^n | Tr U_BAM | Tr U_Bloch | |Tr|² residual | op residual | PASS? |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `equatorial_2pi_chi_pi_over_3` | False | 1 | -1 | (-0.0000+2.0000j) | (+0.0000-2.0000j) | 0.00e+00 | 9.96e-17 | **PASS** |
| `polar_cap_with_throat_at_phi_pi` | True | 1 | -1 | (+0.0000+0.0000j) | (+0.0000+0.0000j) | 0.00e+00 | 3.74e-16 | **PASS** |
| `equatorial_throat_sweep` | True | 1 | -1 | (+0.0000+0.0000j) | (+0.0000+0.0000j) | 0.00e+00 | 1.52e-16 | **PASS** |

## Verdict

**PASS — BAM symmetric gauge is gauge-equivalent to Bloch gauge under the multi-valued scalar U(1) transformation g(χ,φ) = e^{−iφ/2}. The two combinations (A_BAM, ψ_BAM, T_BAM=iσ_y) and (A_Bloch, ψ_Bloch, T_Bloch) make identical predictions for every SO(3)-level observable. They differ by (−1)^n per n-winding closed loop at the SU(2) level — the Wu-Yang transition function that IS the spin-½ double cover. This sign is observable only in interference between coherent spinor branches (Bell singlets), where BOTH gauges produce the same observable phenomena because the relative phase cancels the gauge offset. No experiment can distinguish the two gauges; the choice is computational, not physical.**

Specifically, the probe established:

- `A_BAM − A_Bloch = ½ dφ` at every sampled point of the base — the connection-level statement of the gauge transformation is exact.
- `ψ_BAM = g · ψ_Bloch` with the scalar U(1) gauge function `g = e^{−iφ/2}` at every sampled (χ,φ) — the spinor relation is a multi-valued scalar phase, not an SU(2) frame rotation.
- `T_Bloch(φ_in, φ_out)` defined by gauge conjugation of `T_BAM = iσ_y` preserves Tr = 0 universally and, in the instantaneous-throat case (φ_in = φ_out), also preserves `det = 1` and `T² = −I`. The throat transport iσ_y in BAM becomes the scaled `e^{i(φ_out−φ_in)/2}·iσ_y` in Bloch — a clean scalar conjugation.
- Every closed loop's holonomy satisfies `U_loop^BAM = (−1)^n_winding · U_loop^Bloch` to numerical precision. The SO(3)-level invariant `|Tr U|²` (which discards the spinor sign) agrees exactly. So the two gauges describe the same SO(3) rotation but lift it to opposite signs in SU(2) for every odd-winding loop.
- Detector-holonomy phase differences `Δγ(θ_a, θ_b)` are identical in the two gauges — the π/loop offset cancels in every relative measurement. Bell experiments using either gauge produce the same CHSH = 2√2.
- No single-valued U(1) gauge transformation can remove the π/loop BAM−Bloch offset (any periodic Λ contributes zero around a closed loop). The π offset is a true Wu-Yang transition function — the spinor double cover content — present in any gauge that uses spinor wavefunctions.

Interpretation. The BAM symmetric gauge and the standard Bloch gauge are two sections of the SAME Hopf bundle, related by the multi-valued scalar U(1) gauge function `g = e^{−iφ/2}`. The multi-valuedness `g(φ+2π) = −g(φ)` is exactly the spinor double-cover transition function; it is unobservable in any projective (ray) Hilbert-space measurement and cancels in every relative-phase observable. Every physical prediction agrees between the two gauges. BAM's derivational economy — `T = iσ_y` falls out directly from the orientation-reversing isometry of S³ — is real but represents a choice of computational gauge, not an experimentally falsifiable prediction.

What this rules out:

- Claims that BAM's symmetric gauge predicts a *novel* phase or interference pattern not present in standard spin-½ quantum mechanics. It does not.

What this leaves intact:

- BAM's **derivation** of `T = iσ_y` from the unique orientation-reversing Hopf-preserving S³ isometry — a geometrically forced spinor map. Standard quantum mechanics POSTULATES this transformation as part of the spinor representation; BAM DERIVES it from the throat's non-orientability. The derivation has content independent of any gauge choice.
- BAM's identification of the equator `χ = π/2` as the natural throat location and the zero-self-energy stable orbit. This is a structural / geometric statement, not a gauge artifact.
- The closure-ledger predictions — BAM uses the symmetric gauge because it places the closure-quantum integers in their most natural form (`π·cos(χ)` per loop, `2π` per pair traversal). The Bloch gauge would produce the same integer ledger after the appropriate offset.

## What this leaves open

- **Sharp derivation vs. sharp prediction.** BAM's symmetric-gauge derivation of `T = iσ_y` may have content beyond what this gauge-equivalence test sees — in particular, the derivation argument in `embedding/transport.py` invokes the UNIQUE orientation-reversing isometry of S³ that preserves the Hopf bundle. Whether that uniqueness statement has additional consequences beyond the spinor frame choice (e.g. a uniqueness theorem for the throat sector at the classical level) is open.
- **Non-Abelian gauge sector.** This probe analyses only the U(1) ↔ SU(2) doubling. If BAM's throat carries additional non-Abelian internal structure (e.g. for the quark sector's colour content via the QCD pinhole), the corresponding gauge-equivalence test would need an extended formulation.
- **Time-dependent (Aharonov-Anandan) corrections.** The present test is purely adiabatic. Whether the gauge equivalence persists under finite-velocity transport is a separate target.
