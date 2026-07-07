# The measurement sector: pointer outcomes for the entangled sector (PR #209)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR closes the last standing open
> of the entangled-sector thread: the spatial/measurement sector. The
> companion probe measures each link of the chain (~6 min).

## 0. The missing chain

#206–#208 derived the entangled sector's *states* from bridge topology
— the singlet, swapping, GHZ — but their operational content (CHSH,
Mermin) rested on "Born statistics for internal states at dBB grade",
while #198's equivariance theorem covers **spatial** transport only.
What was missing is the measurement sector:

> internal (fiber/spin) state → spatial pointer branches → position
> beables → Born statistics.

Three links, each derived or measured here, plus the spatial sector
itself (the pointer *is* spatial structure; positional EPR; and the
#205 guiding-without-gravitating split realized where it matters).

## 1. The coupling exists in the committed structure

Winding couples minimally to the fiber connection — winding = charge
(#42–#44), the KK gauge coupling. On the lattice this is exact: the
χ-hop with connection θ(x) leaves the winding channels decoupled with

```
V_k(x) = −2 t_χ cos(2πk/N_χ − θ(x)) ,
```

whose **k-odd part** is −2t_χ sin(2πk/N_χ) sin θ — verified as an
identity to 10⁻¹⁶. A connection-*gradient* region therefore exerts
opposite forces on opposite windings: **the winding Stern–Gerlach**,
which is nothing exotic — it is charge measurement by deflection in a
gauge-field gradient. Live, from the raw dispersion (no linearization):
a k = +1 packet crosses the ramp (final centroid +28) while k = −1 is
turned back (−23) — opposite windings, opposite sides: a pointer. (For
the transported-frame/spin doublet — the #208 GHZ carrier — the
analogous device is a fiber-*geometry* gradient: the Berger-squash
region of #192/#197.)

## 2. Fiber-integrated equivariance (the #198 extension)

For the multichannel wave — channels are orthogonal internal states,
per-channel potentials real — the fiber-integrated density and current
close the continuity equation exactly, so a throat transported by
v = J/ρ (fiber-integrated) is equivariant **through the measurement
interaction**. Verified on the live Stern–Gerlach evolution: continuity
residual ~10⁻⁴ (integrator error); a 20 000-throat Born ensemble stays
at sampling noise through branch separation. This is the measurement
theorem: spatial equivariance, which #198 proved, delivers
internal-state statistics once the committed coupling correlates
channel with position.

## 3. Born statistics for internal states, measured

The winding superposition cos β|+⟩ + sin β|−⟩ through the device;
outcomes = final beable positions (branches separated by 7σ):

| β | P(+) measured | cos²β | diff |
|---|---:|---:|---:|
| π/8 | 0.8516 | 0.8536 | −0.0020 |
| π/6 | 0.7477 | 0.7500 | −0.0024 |
| π/4 | 0.4974 | 0.5000 | −0.0027 |
| π/3 | 0.2500 | 0.2500 | −0.0000 |

**P(k) = |a_k|² to ≤ 0.003 (sampling noise)** — the Born rule for the
internal sector is not a new postulate; it is #198's spatial
equivariance routed through the committed coupling.

**Pointer permanence:** deleting the empty branch changes the guidance
of throats in the occupied branch by ~10⁻⁴ during separation, decaying
with the Gaussian branch overlap to 7×10⁻⁹ by the end — the empty branch
becomes dynamically irrelevant as the pointer forms: **effective
collapse, from geometry, without collapse**. (Permanence against
deliberate recombination is amplification/irreversibility — scope §6.)

## 4. The operational closing: CHSH = 2√2 from pointer positions

The full loop, with every ingredient derived upstream: the
**#206-derived singlet** (the bridge state), local setting rotations
(the fiber-frame rotation before the device), Stern–Gerlach branches at
both wings (exact accelerated Gaussians), dBB transport of throat pairs
on (x₁, x₂), outcomes = sign(x):

- sanity: E(0,0) = **−1.0000** (perfect anticorrelation);
- the four CHSH correlators land within 0.012 of −cos(a−b);
- **CHSH = 2.824** (Tsirelson 2.828) — **Bell violation as pointer
  statistics of classical position beables**;
- operational no-signaling: Alice's outcome marginal shifts by ≤ noise
  when Bob changes his setting — the #204/#206 marginal theorems, now
  at the pointer level.

State: #206 topology. Coupling: §1. Statistics: §2 (#198). Equilibrium:
#198/#204. No imported quantum rule anywhere in the chain.

## 5. The spatial sector

- **The pointer is the spatial sector**: in a measurement, the spatial
  part of the pair wavefunction is exactly the branch structure of
  §3–§4 — the open's content is the measurement chain itself.
- **Positional EPR from nucleation kinematics**: the #58/#200
  C-conjugate pair is born at a *local* event with conserved total
  momentum — anticorrelated momenta, co-located birth: the pair's
  spatial wave is the EPR state. On the effective Gaussian pair state
  (σ₋ = 0.5, σ₊ = 4): Var(x₁−x₂) + Var(p₁+p₂) = **0.53 < 2** — the
  Duan–Simon bound for separable states is violated: spatially
  entangled, with the correlation *structure* symmetry-derived and the
  widths as inputs.
- **The #205 split, realized in measurement**: the empty pointer branch
  *guides* until separation (§3: guidance is fiber-integrated over both
  branches while they overlap) and *never gravitates* (#205:
  conditional sourcing, measured); gravitational back-action on
  outcomes is the #205 size, ~10⁻¹⁷ — nil.
  Guiding-without-gravitating is not a paradox; it is the division of
  labor between the phase (guidance) and the mass (sourcing), and the
  measurement context is exactly where it operates.

## 6. What is and is not established (honest scope)

**Established:** the pointer coupling from the committed KK structure
(identity + live deflections); fiber-integrated equivariance (residual
10⁻⁴; ensemble at noise); Born for internal states (≤ 0.003); pointer
formation and the decay of empty-branch influence; the operational Bell
test (CHSH 2.824, marginals flat); positional EPR from conservation.

**Not established (the conditions):**
1. The SG force window is the co-moving description of a
   connection-gradient transit (the raw-dispersion device runs live in
   §1; §3–§4 use the leading linear KK coupling).
2. **Registration/irreversibility**: branch separation makes the
   pointer; permanence against deliberate recombination requires
   amplification/radiative decoherence — the machinery of #204's
   dissipative controls, not modeled here.
3. Equilibrium hypothesis (#198/#204), as throughout.
4. The nucleation spatial **widths** are inputs (structure from
   conservation; the 5D profile underived).
5. Program-wide remaining opens: the 5D pants nucleation (#208),
   W-class reachability (#208), and the strong-field NR core
   contraction (#203's target, the mass thread).

## References

- D. Bohm, Phys. Rev. 85 (1952) 166, 180 [the measurement account];
  D. Dürr, S. Goldstein, N. Zanghì, J. Stat. Phys. 67 (1992) 843.
- L.-M. Duan, G. Giedke, J. I. Cirac, P. Zoller, PRL 84 (2000) 2722;
  R. Simon, PRL 84 (2000) 2726. [The separability criterion.]
- A. Einstein, B. Podolsky, N. Rosen, Phys. Rev. 47 (1935) 777.
- The BAM chain: #42–#44 (winding = charge; the KK coupling), #198
  (equivariance), #204 (no-signaling; equilibrium), #205 (conditional
  sourcing), #206–#208 (the entangled states from topology).

## Reproduce

```bash
python -m experiments.closure_ledger.measurement_sector_probe
# Verdict: THE_MEASUREMENT_CHAIN_CLOSES_THE_KK_COUPLING_MAKES_THE_POINTER
#          _FIBER_INTEGRATED_EQUIVARIANCE_MAKES_BORN_OPERATIONAL_CHSH_2SQRT2
#          _FROM_BEABLE_POSITIONS
```
