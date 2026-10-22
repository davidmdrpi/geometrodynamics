# Coherent neck stress and the missing tangential recoil channel

**Corrected after the [review of 11bf1d3](pr348_review_response.md).** The neutral,
centred single-mode calculation supports interference-dependent **neck-normal
stress**, not orbital attraction. Its tangential resultant is exactly zero by
azimuthal symmetry. The original symmetric spectral source also made normalized
pressure blind to the surface phase. The repaired source controls now see that
phase. These corrections supersede the initial candidate-attraction reading.

**This is not yet self-consistent motion of the two mouths.** We derive a local
stress kernel, complete individual collisions with finite-mass energy/momentum
balance, and compare multiple paths on fixed geometry. For a centred m=0 mode,
the tangential radial flux integrates as `∫(cosφ,sinφ)dφ=0`; the remaining
resultant runs along the neck axis. The fixed-geometry cumulative impulses must not be applied wholesale
to moving mouths. No binding, collapse, shrinking throat or QED derivation is
claimed. Earlier negative results remain intact.

![Conservative recoil mechanism comparison](images/history_throat/coherent_recoil.png)

## 1. Neutral junction stress and the separate charged cross term

Let the energy-normalized incident amplitudes be `a=(a_s,a_t)`, with incident
longitudinal directions `(+1,-1)` and outgoing directions `(-1,+1)`. For a
single matched transverse mode, equal wave speeds and positive characteristic
admittance ratio `z=Z_throat/Z_surface`, continuity of value and flux gives

\[
 b=S a,\qquad S=\begin{pmatrix}r&t\\t&-r\end{pmatrix},\qquad
 r=\frac{1-z}{1+z},\quad t=\frac{2\sqrt z}{1+z},\quad S^\dagger S=I.
\]

The real coefficients describe a reciprocal, lossless junction; they are not
derived from a resolved MTY mouth. The wave energy is `E=|a_s|²+|a_t|²`.
With `c=1`, incoming minus outgoing momentum gives

\[
 \Delta p=2r^2(E_s-E_t)+4rt\sqrt{E_s E_t}\cos\phi,
 \qquad \phi=\arg a_t-\arg a_s.
\]

The first term is ordinary reflection pressure from unequal incoming fluxes.
The second is coherent interference. Multiplying **all** amplitudes by −1
changes neither term. Flipping one channel changes only the interference term.
A single incident wave at an **uncharged** passive junction is unaffected by
its compression sign. This does not exclude a sign-linear force on a charged
flux mouth. With a static field E0 and wave field e, the electric stress cross
term (in units ε=1) is

\[
 T_{cross}=E_0\otimes e+e\otimes E_0-I(E_0\cdot e).
\]

It is linear in the wave and reverses with its sign. The implementation checks
this against `T(E0+e)-T(E0)-T(e)`. Integrating it over a small circle in a local
tangent patch, with `E0=q e_r/(2πr)` and a uniform external e, gives **q e**.
The reported quadrature residual is `6.9e-15`. This permits a tangential force
when the incident field has the required directional content; it is excluded
from the neutral centred m=0 junction, not from charged mouths in general.

Thus #345's sign-linear rule can be motivated by static–wave interference.
Its normalization by the front's angular momentum is still a stated rule:
the physical electric impulse is `q∫e dt`, with magnetic and motion corrections
where applicable. This work does not derive that field integral from the
topological mouth or from its global advanced boundary conditions.
Identifying a normal-compression sign with the relevant tangential electric
field also requires the field/mode projection; the local identity alone does
not show that a Klein flip makes the orbital interaction attractive.

There is an independent virtual-work check. For a translated junction,
`S(x)=exp(iωJx) S(0) exp(iωJx)`, with `J=diag(1,-1)`,

\[
 {-i\over\omega}S^\dagger\partial_x S\big|_{x=0}
 =J+S^\dagger J S=K,\qquad \Delta p=a^\dagger K a.
\]

The eigenvalues of `K` are `±2|r|`, so `|Δp| ≤ 2|r| E`. This is a bound per
encounter, not a bound on repeated pressure over an arbitrary time.

For equal incoming energies, `Δp/E=2rt cosφ`. Its uniform phase average is zero.
For unequal inputs the reflection bias remains after phase averaging. A throat-
dominated incident field can give negative pressure without coherence; that
does not explain how such an incident field is generated or maintained.

## 2. Partial reflection and sizing

For **balanced input energies**, maximum phase-sensitive recoil occurs at
`R=r²=T=t²=1/2`. The radius giving that admittance depends on the space:

\[
 z=(b/a)^{D-1},\qquad
 (b/a)_{D=2}=3-2\sqrt2=0.1715728753,\qquad
 (b/a)_{D=3}=\sqrt2-1=0.4142135624.
\]

The 2D surface reduction uses rim circumference; a 3D space uses area. The
original 9,100-case grid remains explicitly a **3D area-law control**, while
the added physical-length, flash-weighted controls evaluate both laws. These
are optima of stated equal-speed, single-mode impedance assumptions, **not a
predicted BAM size or a universal optimum**. For arbitrary input energy ratios
the optimization changes. With `R=T=1/2` and `E_t≤E_s`, a negative single-
encounter impulse is possible only when `E_t/E_s > 3-2√2 ≈ 0.1716`, and the
relative phase is sufficiently close to π. Multiple returns can change that
input ratio, so it is not an initial-seed threshold for the entire network.

The independent Gaussian focus-overlap control gives
`η=1-exp[-a²/(2w²)]`: 0.118, 0.393, 0.865 and 0.989 for `a/w=0.5,1,2,3`.
Unmatched energy stays in the bulk. This assumes a circular Gaussian transverse
profile; it is not an aperture diffraction derivation. Spatial overlap,
temporal overlap and polarization/mode overlap are all needed for the cross term.

The repository's `physical_throat.VacuumThroat` supplies a useful geometry
control. Its scalar-flat, time-symmetric bridge has `b/a_areal=sin²(a_angular)`.
Passing those actual ratios through the assumed area-impedance model gives:

| Angular mouth radius | Neck / mouth | T, 3D area | Max balanced `|Δp|/E`, 3D | T, 2D rim | Max balanced `|Δp|/E`, 2D |
| --- | --- | --- | --- | --- | --- |
| 0.05 | 0.002498 | 0.00002496 | 0.00999 | 0.00994 | 0.19842 |
| 0.16 | 0.025382 | 0.002574 | 0.10133 | 0.09656 | 0.59073 |
| 0.30 | 0.087332 | 0.03005 | 0.34144 | 0.29547 | 0.91251 |
| 0.50 | 0.229849 | 0.19065 | 0.78562 | 0.60785 | 0.97646 |

The suppression inferred from the area law must not be transferred to the 2D
surface model: at angular radius 0.3 its mixing is already near optimal.
Both laws suppress mixing in the limit of a vanishing neck, at different rates.
This bridge
is spatial initial data, not an established traversable dynamical throat;
its ratios must not silently be imposed on an MTY solution.

## 3. Multiple paths, finite inventory and Klein reversal

Two fixed antipodal nodes are connected by a surface path of length π and a
throat path of length L. Each has the junction above. Surface interactions occur
only at the next antipodal node; throat interactions occur only at its ends.
This is an aligned-node reduction of #347, not a continuous local surface force.

The surface factor is either −i (post-caustic capture) or 1 (pre-caustic capture),
retaining #347's explicit open-boundary choice. A throat traversal carries the
**assumed compression holonomy** −1 for Klein and +1 for torus. This changes a
mixed surface/throat loop phase. It does not implement the tube's full momentum
map: the latest #345 review correctly distinguishes compression sign, spin and
the two orbital components. In particular a planar orbital reversal is present
for both of its tube maps.

Two complementary computations avoid confusing a driven resonance with free
energy:

1. **Closed finite pulse:** lossless delay lines and unitary scattering evolve
   one initial energy inventory of 1. There is no later injection or gain.
   The bins are **complex envelopes of a positive-frequency carrier**, not a
   real broadband flash. A real caustic has multiplier `−i sgn(ω)` and reshapes
   the pulse by a Hilbert transform. The envelope experiment does not compute
   that reshaping or its change in overlap. Zero-frequency and integration-
   window conventions must also be specified before using pulse area to infer
   an impulse; the `Re(phase)` bookkeeping alone does not settle them.
   An intensity-only `|S|²` control removes coherence while retaining energy.
   We scan 560 parameter/seed combinations, each with both evolutions, plus 30
   pulse-width/delay controls and 10 longer-duration controls.
2. **All-walk boundary susceptibility:** `a=b+ρ P S a` is solved as
   `(I-ρPS)^{-1}b`, and independently compared with its explicit walk sum.
   Here `ρ=0.8` is amplitude survival; escaped power is
   `(1-ρ²)||Sa||²`. This calculation has a source and is explicitly externally
   fed. Its pressure is divided by total incident exposure. There are 9,100
   combinations of frequency, length, size, handle and capture convention.
   `ρ=1` is rejected: a closed pole needs its homogeneous boundary solution.

### Source degeneracy and repaired spectral controls

The original source `(1,0,−1,0)` restricts the solution to `a_B=−a_A` at zero
clock offset. Writing its reduced loop factors as `D=diag(d_s,d_t)`,
`a_A ∝ (1+r d_t, t d_t)`: its direction depends on the throat loop only.
The surface factor changes exposure but cancels from impulse/exposure.
The old before/after-caustic agreement and Klein half-period shift were thus
structural consequences of the source, not independent evidence about the
surface path or topology selection.

The main grid now uses `(1,0,0,0)` at strictly positive frequencies; it does not
assign the envelope caustic factor to a static/DC mode. Twenty-four added controls compare it with
a normalized biased pair `(1,0,−0.35 exp(0.37i),0)` and the normalized old
antisymmetric pair. The old source's maximum caustic difference is `5.9e-16`;
the one-sided source's is about **0.975** in these controls. Tests independently
check the analytic throat-only reduction and require the asymmetric source to
respond to the surface phase.

Over a full uniform throat phase, the old symmetric control gives identical
negative fractions for both handles: approximately 0, 0.105 and 0.242 at 3D
radius ratios 0.05, 0.2 and √2−1. The former 14.9% versus 15.6% comparison was
a frequency-grid sampling artifact and is withdrawn as a topology comparison.

Each spectral file now reports the source explicitly. Alongside the normalized
ratio, `spectral_pressure_per_source.json` reports impulse per source energy and
`exposure.json` reports the incident exposure. For spectral energy weights w,
the relevant aggregate is `Σ w P / Σ w exposure`; it is generally different
from `Σ w (P/exposure)`. Both are saved to make that distinction reviewable.

### Actual flash energy weights and the model's neck length

The 192 additional controls use `L=sin(2a)` at angular radii 0.05, 0.16, 0.3
and 0.5, with both dimensional admittance laws and balanced/vacuum geometry
ratios. The weights come from the actual `PairBirth.flash` at speed 1/3,
footprint 0.12 and cutoff 32. Source-copy fields are summed **before** squaring,
retaining their interference in the total modal energy. The extracted
frequencies remain the Maxwell values `sqrt(l(l+1))`, not `l+1/2`.

This supplies the birth flash's energy spectrum, **not its missing projection
into neck modes**. All three incident-source choices remain explicit controls.
The exact-focus caustic factor is still an imposed routing assumption; it is
not the exact propagator of those Maxwell modes. Frequency-weighted stationary
stress also does not reconstruct a real flash's time-resolved overlap.

For example, at `a=0.3`, `L=sin(0.6)`, the vacuum radius ratio, 2D rim law and
one-sided source:

| Capture convention | Handle | Weighted impulse / source energy | Weighted exposure | Exposure-weighted normal pressure |
| --- | --- | --- | --- | --- |
| Before caustic | Torus | 0.87954 | 1.08901 | 0.80766 |
| Before caustic | Klein | 0.84235 | 0.81394 | 1.03490 |
| After caustic | Torus | 6.88062 | 5.42781 | 1.26766 |
| After caustic | Klein | 9.46489 | 6.58097 | 1.43822 |

All four are positive normal pressures. Neither this example nor any negative
normal-pressure case establishes a tangential separation force.

The 32π fixed-geometry results are:

| Initial throat share | Initial relative phase | Negative cumulative port impulse / cases |
| --- | --- | --- |
| 0 | No initial throat field | 1 / 140 |
| 0.1 | π | 69 / 140 |
| 0.5 | Both 0 and π tested | 147 / 280 |

These counts are retained for provenance, **not as evidence of attraction**.
The surface-seed exception (`b/a=0.8`, `L=π`, Klein, pre-caustic) has normalized
normal impulse −0.00703 at 32π, +0.000189 at 64π, −0.000189 at 128π and
−0.000188 at 256π. It oscillates near zero at the 10^-4 level in the later
samples; it is not a sustained negative-stress mechanism.

Prepared balanced opposite-phase returns can retain a negative port stress:
at `b/a=√2-1`, `L=π/4`, post-caustic capture, its normalized cumulative value
changes from about −0.834 at 32π to −0.828 at 256π. This is a fixed-support
cavity stress, **not a moving-mouth trajectory**. The supports take the impulse;
they do no work because they do not move. Integrating that stress into a moving
mouth while leaving the same wave inventory and phases unchanged would violate
the required backreaction calculation.

No relaxation or attractor is present in the closed unitary network to choose
the desired coherent preparation automatically.

## 4. Advanced paths are boundary data, not an energy supply

For stationary mouths with `σ_A=t`, `σ_B=t+δ`, the directed throat delays are
`L-δ` and `L+δ`. The added 650 frequency-domain controls include `δ>L` and
`δ=π+L`, where one mixed geometric route has zero total coordinate-time delay.
These delays are checked against #347's actual clock-map inverses.

The propagation matrix remains unitary for real frequency, including negative
delays. An advanced phase changes interference but contributes no extra energy
factor. The regularized, externally fed resolvent is a **boundary diagnostic**,
not a causal initial-value evolution or proof of stability near a closed null
route. Advanced events are never inserted into the forward closed delay-line
simulation. A time-orientation reversing attachment, as opposed to a clock
offset, also requires opposite-frequency sectors and is not included in this
scalar fixed-frequency closure.

## 5. Finite-mass recoil: a local conservation step

For each encounter, boost both incident packet energies into the mouth's
instantaneous rest frame, apply `S`, then rescale the outgoing energies by a
common positive `k`. If `E,p_in` are the total incident rest-frame energy and
momentum, and `p_out` is the unrescaled outgoing momentum, impose

\[
 (M+E-kE)^2-(p_{in}-kp_{out})^2=M^2.
\]

The physical root is evaluated without subtractive cancellation and the result
is boosted back. Both energy and momentum close; the mouth remains on its mass
shell. Perfect reflection reproduces the independent photon/mirror recoil
redshift `E_out=ME/(M+2E)` for a stationary mouth. Matched transmission produces
no impulse. Tests also cover moving mouths and complex inputs.

The common energy rescaling is an explicit **constitutive completion**, not a
derived Doppler-resolved scattering spectrum. It does not yet tell us the
coherence of packets at the next encounter. This is why the finite-mass result
is not silently chained into the fixed-delay interference calculation.

Across 315 collisions (`v=-0.3,0,0.3`, wave/rest-energy ratios 0.001 to 1,
seven sizes and three phases), worst energy and momentum residuals are
`6.7e-16` and `4.5e-16`. An illustrative equal-input, opposite-phase encounter
at the optimal balanced size, `v=0.3`, and incident `E/M=0.0063`, changes mouth
**neck-axis** momentum by `−0.0086364 M` from an initial `0.314485 M`.
Outgoing wave energy grows from `0.0063 M` to `0.0088585 M`, paid by mouth
kinetic energy along that axis. The orbital velocity on the sphere is
perpendicular to this axis: **this number is not orbital braking**. It is
neither a replay of #345's flash nor evidence that the required paired returns
occur on a moving worldtube.

## 6. Mechanism assessment and the remaining closure

| Mechanism | What this work establishes | Assessment |
| --- | --- | --- |
| Compression flip of an isolated front at an uncharged junction | No change in quadratic stress at fixed energy | Does not exclude the static–wave cross term of a charged mouth |
| Unequal counter-incident fluxes | Neck-axis reflection pressure has either sign depending on flux direction | No tangential resultant for the centred m=0 mode |
| Coherent returns with partial reflection | Normal stress, axial translation work, individually conservative axial recoil | Relevant to normal deformation; not a demonstrated orbital attraction mechanism |
| Repeated paths / resonance | Alters occupation and phase; finite seed energy is conserved | Can maintain fixed-geometry stress, not unlimited acceleration |
| Advanced clock paths | Shift boundary phases without increasing flux norm | May affect coherence; no energy-generating mechanism established |
| Charged static–wave cross stress | Local tangent-patch integral gives q E and changes sign with the wave | Candidate tangential interaction; global source, mode projection and advanced-history recoil remain unresolved |

The concrete orbital-force extension needs directional aperture content, such
as m=±1 from an off-centre focus or a momentum-carrying wave. The centred m=0
junction excludes it. An arbitrary dipole mixing coefficient would simply add
a new force assumption, so none is fitted here. A charged static–wave cross
term gives a physically grounded local target for the missing coupling.

Normal stress may instead act on a neck deformation. To claim a force on a
radius a requires the shape-work projection
`Q_a = ∫ traction · (∂X/∂a) dA`. Axial stress alone is not transverse neck
tension and does not establish shrinkage, barrier crossing or collapse.

The next closure is to obtain a frequency- and position-dependent
multi-component mouth scattering matrix from the same bulk/throat field action,
including finite aperture modes and the actual Klein ray/normal transport.
Its translation/rotation derivatives and stress flux must give the same orbital,
spin and radial work. Couple those forces to the mouth histories, update the
packet frequencies and antipodal interception times, and solve the advanced
branches as boundary conditions. Then test whether the **unprepared** birth
field selects a negative separation force, maintains overlap and pays for throat
shrinkage. A prescribed favorable relative phase does not close that argument.

## Reproduction and provenance

```bash
python -m experiments.closure_ledger.coherent_recoil_probe
python -m pytest -q tests/test_coherent_mouth_recoil.py tests/test_antipodal_throat_transport.py tests/test_history_throat.py
```

70 tests pass (26 recoil tests and 44 preceding regressions). Maximum closed-network
energy error is `7.6e-14`, including the 256π duration controls;
explicit-walk/resolvent agreement is `2.0e-15`, including advanced clock controls.
Integer-delay bin rescaling is an exact identity and has been removed from the
validation gates; no continuum convergence claim is made.

All grid data, pulse/duration/clock controls, geometry controls and hashes are in
[`runs/20261022_coherent_recoil`](../experiments/closure_ledger/runs/20261022_coherent_recoil).
The scalar candidate's failure of sustained approach/collapse and #345's signed-
exchange negative binding result are preserved. Its 22 trajectories were not
rerun or relabelled as outputs of this model.

The external motivation is methodological, not evidence for BAM: classical
optical binding through mode scattering/interference is demonstrated by
[Bykov et al., arXiv:1709.04053](https://arxiv.org/abs/1709.04053), and scattering-
matrix position derivatives are used by
[Ambichl et al., arXiv:1703.07250](https://arxiv.org/abs/1703.07250).
Neither supplies this throat's constitutive law or establishes its attraction.
For the charged-field qualification see the
[MIT Maxwell-stress notes](https://ocw.mit.edu/courses/6-642-continuum-electromechanics-fall-2008/resources/lec06_f08-1/).
The real-pulse versus analytic-signal distinction follows the
[Hilbert-transform convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html);
no broadband time-domain transform is implemented in the envelope network.
