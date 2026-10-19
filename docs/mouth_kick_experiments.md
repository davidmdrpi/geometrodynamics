# Single wave exchanges between antipodal mouths: kicks, throat flux and clock offsets

Classical scalar experiments on #338's glued sphere, building on #344. These are
leading-order kicks with the geometry held fixed: the mouths do not yet move.
"Exchange" means a classical wave that leaves one mouth and reaches the other.
No quantization rule enters anywhere.

**Questions.**
1. **Kick scan.** Can one wave exchange start to move two antipodal mouths toward each other? Through the throat, over the sphere, or both?
2. **Throat flux.** With a static flux threading the throat, one mouth is a source and the other a sink. Does the kick then gain a part that depends on the relative sign of the wave and the flux?
3. **Clock offset.** With a mouth clock offset (the classical Morris–Thorne–Yurtsever ingredient), does a forward-propagating wave re-emerge from the far mouth earlier in external time? What happens to the kicks?

**Answers in brief.**
1. **Yes, the exchange alone moves the mouths together.** A wave released next to mouth A gives A and B equal kicks in the same direction, which is pure approach. Most of it goes over the sphere. The throat carries about 7.5% of the sideways kick, and the handle's orientation matters only there.
2. **Yes, the cross term is sign-dependent and large.** It is linear in the throat flux and in the wave, and at equal energies it is as large as the wave's own kick or larger.
   - When the wave and flux signs are aligned, B follows A much more strongly: follow is 2.6–3.0 for the release, against 1.0 without the flux.
   - When they are opposed, B's following is cancelled or reversed.
3. **Yes, the wave re-emerges from B earlier by exactly the offset.** It reaches A's rim at t = 0.47 for every offset, and leaves B at t = 1.09 + Δ.
   - **Before entry:** for Δ < −0.61 it leaves B before it enters A.
   - **Before release:** for Δ < −1.09 it leaves B before it is released.
   - **Toward the horizon:** the external-time period of the A → throat → B → sphere → A loop is 3.15 + Δ. It shrinks to zero at the chronology horizon Δ = −3.11.

## Model
- **Surface:** `viz.klein_handles.build_surface` (#338). A unit sphere carries mouths of angular radius 0.3 at A = +x and B = −x.
  - The mouths are joined by the catenoid neck r(s) = √(s² + b²), with waist b = 0.087 and length 0.565.
  - The neck is glued as a Klein handle (rim angle ψ → π − ψ, orientation-reversing) or a torus handle (ψ → ψ + π).
- **Field:** u_tt = Δu − m²u, on the intrinsic cotangent mesh with lumped mass. Two values of m² are used:
  - m² = ¼, the #339/#344 conformal operator, for the scan;
  - m² = 0 for the flux and clock-offset parts, because a period around the handle is incompatible with a mass term. At m² = 0 the scan's kicks change by 0.7%.
- **Kick on a mouth:** the angular momentum the sphere band hands to its rim.
  - On the band, rotations are Killing fields, so ℓ = −u_t ξ·∇u is conserved locally.
  - A smooth rim weight turns the rim flux into area integrals (`geometrodynamics/waves/rim_kicks.py`).
  - Each mouth is its own body; neck stiffness is not modelled.
- **Waves:** a single-signed bump of width 0.25 and unit energy, starting 1.0 rad from A on side u(ψ₀), where ψ₀ is measured from the handle's mirror axis e₁. Three launches are used:
  - **release:** the bump starts at rest and spreads both into A and over the sphere;
  - **out:** A emits a beam toward B over the sphere, and the launch recoil is booked to A;
  - **in:** a wave arriving from elsewhere is sent into A.

**How "toward each other" is judged at antipodes.** Let a and b be the ambient directions in which A and B are pushed. The separation is π − |a + b| to first order, so any kick that is not a rigid rotation shortens it, and repulsion is impossible. The real test is whether B follows A, measured by **follow = b·a/|a|²**:
- follow > 0 means both slide toward the same side, so they close in;
- follow < 0 means the pair mostly rotates.

## 1. Kick scan
Two handles × three launches × eight entry sides, each to one full loop (2π + ½), on a 24,000-vertex mesh. Kicks are ×10⁻³ per unit wave energy. "First-pass plateau" means t = 3.0–3.7, after the first passes over both mouths and before the second; the kicks are steady there.

| launch | handle | entry side | A kick | B kick along A | follow (plateau) | follow (one loop) |
|---|---|---|---|---|---|---|
| release | Klein or torus | 0° | 44.8 | +45.5 | **+1.02** | +1.47 |
| release | Klein | 90° | 31.7 | +33.0 | **+1.04** | +0.19 |
| release | torus | 90° | 44.6 | +45.2 | **+1.01** | +1.46 |
| out (A emits) | either | any | 642–679 | +104 to +105 | +0.15 to +0.16 | +0.15 to +0.21 |
| in (throat) | torus (any side); Klein at 0° | — | 104 | +14.4 | +0.14 | −0.41 |
| in (throat) | Klein | 45° | 105 | −2.3 | −0.02 | −0.66 |
| in (throat) | Klein | 90° | 104 | −18.3 | **−0.18** | −1.16 |

![Exchange kicks](images/mouth_kicks/exchange_kicks.png)

**Mechanism.**
- **Emitter and receiver kicks align.** A wave that leaves A over the +u side pushes A toward −u. It reaches B travelling −u and pushes B toward −u as well, so the two antipodal mouths close in on the −u side.
- **In flat space the same exchange would separate them.** On the closed sphere, separating along the exchange path is approaching around the back.
- **A wave that passes over both mouths is different.** An outside wave that pushes A and carries on to B pushes both along its path, and at antipodes that is a rigid rotation. That is why the "in" launch turns into rotation over the loop.

**Throat channel.**
- **The narrow waist reflects most of the sideways part.** The waist (radius 0.087) is narrower than the wavelength, so it reflects most of the m = 1 part of the wave. B's exit kick is 7.56% of A's entry kick, and the same at 48,000 vertices. A wider packet started farther out gives 6.2%.
- **The two handles map that part differently.**
  - The torus map sends cos(ψ − ψ₀) to −cos(ψ − ψ₀), so B follows A whatever the entry side.
  - The Klein map sends it to −cos(ψ + ψ₀), so B's kick scales as cos 2ψ₀: B follows A along the mirror axis, gets nothing at 45°, and is pushed the opposite way at 90°.
- **The absolute sign is a wave effect.** It comes from the relative phase of the transmitted m = 1 and m = 0 parts, and ray optics would predict the opposite sign for the torus. The cos 2ψ₀ relation between the two handles is a symmetry statement.

## 2. Flux through the throat
**The static field.** A massless field admits a static harmonic field whose potential gains a period once around the handle loop:
- on the band, φ = c log tan(θ_A/2);
- in the neck, φ = c asinh(s/b).

This is Wheeler's charge without charge, for a scalar potential.
- **Source and sink.** Field lines leave A into the sphere and return through the neck, so A is a source and B a sink.
- **Period.** The period is 4c log cot(a/2), and the neck carries exactly half of it, because asinh(cot a) = log cot(a/2).
- **Same field for both handles.** The field is axisymmetric, so the Klein and torus handles carry the same one.

**The cross term.** With the field u = φ + w, the kick gains a term bilinear in φ and the wave w (`RimKicks.cross_rates`). Flipping either sign flips it exactly.
- **Size.** At equal (unit) wave and flux energy it is comparable to the wave's own kick, or larger.
- **No plateau.** Unlike the quadratic kick, it does not settle within the first-pass plateau, because it follows the passing wave's sign. The table gives plateau means with ranges.

| launch | handle | quadratic A, B | cross A | cross B (range) | follow, signs aligned | follow, signs opposed |
|---|---|---|---|---|---|---|
| release | torus, or Klein at 0° | 44.8, +45.4 | −8.2 | +50 (+35 to +71) | **+2.64** (+2.10 to +3.57) | **−0.08** (−0.51 to +0.21) |
| release | Klein 90° | 31.7, +32.8 | −4.4 | +50 (+35 to +71) | **+3.03** (+2.52 to +3.77) | **−0.48** (−1.14 to −0.05) |
| in | torus, or Klein at 0° | 104.5, +13.5 | +19.3 | +67 (+26 to +115) | +0.66 (+0.33 to +1.05) | −0.62 (−1.26 to −0.12) |
| in | Klein 90° | 105.0, −19.4 | +24.9 | +67 (+26 to +115) | +0.36 (+0.07 to +0.68) | −1.10 (−1.89 to −0.51) |
| out | either | 682, +106 | −12.8 | +15 (0 to +34) | +0.18 | +0.13 |

![Throat flux](images/mouth_kicks/flux_cross.png)

**Reading.**
- **"Aligned" means a compression pulse released next to the source mouth.** The large cross term on B comes from the sphere-route wave passing through B's strong flux field near its rim. This is the classical force on a charge from a passing field, linear in both.
- **Compression and rarefaction now differ.** The wave's own (quadratic) kick cannot tell them apart; the cross term can.
- **This is the sign-dependent ingredient that radiation pressure lacks.** With aligned signs the mouths close in strongly. With opposed signs, B's following is cancelled or reversed.
- **What a single handle cannot test.** Its two mouths are always opposite (source and sink), so this does not yet compare like with unlike charges. That needs two handles, or a reversed flux at one of two pairs.

## 3. Clock offsets (classical MTY)
#338's `HandleWaveSim` can glue the waist with a clock offset Δ: a signal crossing from A to B emerges Δ later in external time, or earlier when Δ < 0.
- **The spacetime.** It is stationary, and chronology-safe while the clock tilt stays below 1.
- **The bookkeeping still applies.** The sphere band is static in external time, so the same kick bookkeeping works on the solver's external-time frames. At Δ = 0, frames 0.0025 apart reproduce the direct leapfrog kicks to 0.5%.
- **The runs:** the release next to A, on the massless field, on a 12,000-vertex mesh. The bump is the solver's own compact, zero-mean bump of width 0.25, not the scan's Gaussian. That is why follow at Δ = 0 here (+2.6) differs from the scan's +1.0.

| offset Δ | clock tilt | field first at rim A | field first at rim B | B − A | follow, T = π + ½ (Klein 0°) |
|---|---|---|---|---|---|
| −2.0 | 0.77 | 0.479 | **−0.927** | −1.406 | −0.34 |
| −1.5 | 0.58 | 0.477 | −0.424 | −0.901 | −0.20 |
| −1.0 | 0.39 | 0.474 | +0.079 | −0.395 | +1.12 |
| −0.5 | 0.19 | 0.471 | +0.581 | +0.110 | +1.87 |
| 0 | 0 | 0.473 | +1.086 | +0.613 | +2.62 |
| +0.5 | 0.19 | 0.475 | +1.576 | +1.101 | +2.33 |
| +1.0 | 0.39 | 0.433 | +1.614 | +1.181 | +2.24 |

Klein and torus handles give the same arrival times to 0.012. The entry sides 0° and 90° agree to 0.01, except rim A at Δ = +1 (0.433 at 0°, 0.480 at 90°).

**Rim times** are the first external time the field at the rim exceeds 0.5% of the launch amplitude.
- **Re-emission at B tracks Δ with slope 1.** The entry at A does not move.
- **Past re-emission.** At Δ = −2 the wave leaves B at t = −0.93, before the release at t = 0 and 1.4 before it enters A.
- **No closed causal curves, and the B rim time stays put at large positive Δ.** The clock tilt stays at or below 0.81, so no closed causal curve exists. At Δ = +1 the sphere-route wave reaches B's rim before the throat wave, so the B rim time stays at about 1.61.

![Clock offsets](images/mouth_kicks/clock_offsets.png)

**Faster and faster exchanges.**
- **The loop period is one way plus the return.** The one-way interval from A's entry to B's exit is 0.613 + Δ. Adding the 2.54 return over the sphere, rim to rim, gives an external-time loop period of 3.15 + Δ.
- **It reaches zero at the chronology horizon.** The horizon is Δ = −3.11, where the clock gain equals the loop length. In this model, repeated exchange at an ever faster rate is what approaching the chronology horizon looks like.
- **Both the period and the horizon are derived.** The period comes from the measured one-way timing and the geodesic return. The runs themselves stop at Δ = −2, where the tilt is 0.77–0.81.

**Kicks with an offset.**
- **The kick history reshuffles.** The throat-transmitted wave now leaves B earlier and meets the sphere-route waves at different times.
- **The kicks at a fixed external time change sign.** At T = π + ½, follow drops from +2.6 at Δ = 0 to about −0.3 at Δ = −2.
- **That is not yet a like-for-like comparison.** Successive passes move relative to any fixed external time, so a fair comparison across Δ needs a horizon tied to the exchange events themselves.

## Checks
| check | result |
|---|---|
| bare sphere: kick on mouths that do not exist | 0.37% of the packet's angular momentum |
| ring symmetric about the A–B axis: sideways kick | 0.12% of the release kick |
| two different rim weights | agree to 3.8% |
| mesh, 24,000 vs 48,000 vertices | 1.3% (throat ratio 0.0756 vs 0.0756) |
| compression vs rarefaction, quadratic kick | identical to the last digit; the cross term flips exactly |
| m² = 0 vs ¼ | 0.7% |
| energy drift | below 8.5e-5 |
| clock-offset frames vs direct leapfrog at Δ = 0 | 0.3–0.5% |
| tests | `tests/test_rim_kicks.py`, 9 fast checks on a coarse mesh |

## Limits and next steps
- **Fixed geometry.** These are first kicks with the geometry held fixed. The mouths do not move, and neck stiffness is not included.
- **Antipodal mouths only.** `build_surface` places mouth pairs antipodally. The separation dependence ("exchanging faster as the mouths get nearer") needs non-antipodal mouth pairs, a change to the mesh builder. The clock offset is so far the only knob that shortens the exchange interval.
- **No settling.** In the closed, lossless sphere the waves keep passing, so kicks are read on the first-pass plateau. The cross term never settles.
- **Single handle.** A single handle always has opposite mouths, so the like-versus-unlike comparison needs two handles.
- **Collapse and release are not connected yet.** Neither connects to the neck collapse and release of #342 and #344. The natural link is the kick and energy ledger here feeding the neck's own energy balance.

## Reproduction
```bash
python experiments/closure_ledger/mouth_kick_probe.py scan      # about 9 min on 4 cores
python experiments/closure_ledger/mouth_kick_probe.py flux
python experiments/closure_ledger/mouth_kick_probe.py offsets   # about 13 min
python experiments/closure_ledger/mouth_kick_probe.py checks
python experiments/closure_ledger/mouth_kick_report.py          # tables, figures, summary.json, sha256.json
python -m pytest -q tests/test_rim_kicks.py
```
Results are in `experiments/closure_ledger/runs/20261019_mouth_kicks/`.
