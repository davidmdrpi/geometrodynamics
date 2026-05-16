# BAM Compton amplitude — structural reproduction probe plan

Follow-on to the Compton-kinematics probe
(`compton_antipodal_kinematics_probe.py`, PR #25). That probe
established that BAM's antipodal-closure picture is kinematically
consistent with standard Compton scattering, identified the
throat-pinch proper-time skew as a recoil-induced (∝ ω²) kinematic
quantity, and cleared the thread to attempt amplitude-level
reproduction.

## The question

Build a BAM transaction amplitude from the natural ingredients

  - **antipodal propagation** — `S³` Green function `G_S3(ψ)` on
    geodesic distance `ψ`, with `G_S3 ~ 1/(4πψ)` near the source;
  - **throat transport** — `T = iσ_y` at each vertex (one factor per
    mouth traversal), with `T² = −I`;
  - **closure phase** — Hopf holonomy `π·cos(χ)` per fibre plus
    SU(2) spin transport `α_spin · θ_transport = θ/2` for spin-½

and ask: does `|M_BAM(θ)|²` reproduce structural features of the
Compton amplitude?

Two layers of structural feature:

  1. **Propagator denominators.** The Compton amplitude has poles at
     `s = m²` (s-channel) and `u = m²` (u-channel). Does the BAM
     amplitude, with `ψ_s, ψ_u` mapped from the Mandelstam invariants,
     have the same pole structure?

  2. **Klein-Nishina angular dependence.** In the Thomson limit
     `ω → 0`, the differential cross section is
     `dσ/dΩ = (r_e²/2)·(1 + cos²θ)`. The angular factor
     `(1 + cos²θ)/2` is the spin-summed, polarization-averaged
     factor from the QED vertex. Can the BAM amplitude, with its
     natural spin-½ phase structure, reproduce this angular
     dependence?

The expected outcome is *partial* structural match. The propagator
pole reproduction is plausible (the `S³` Green function has the
right `1/ψ` singularity to map onto `1/(s−m²)` under
`ψ ∝ s−m²`). The Klein-Nishina angular factor is harder —
it comes from a vector-current polarization sum which BAM does not
have explicit machinery for (photons in BAM are not yet given a
throat-pair representation, and the existing `bell/hopf_phases.py`
treats them only at the phase level).

The probe will land in one of three regimes:

  - **(A) Full structural match.** Both propagator poles AND
    angular dependence reproduced from BAM ingredients with no
    additional inputs. Surprising and would indicate strong BAM
    progress.
  - **(B) Partial match (propagator only).** BAM reproduces the
    pole structure but the angular dependence is wrong (e.g.,
    `1 + cos θ` instead of `1 + cos²θ`). This is the expected
    outcome and identifies the missing structure: explicit photon
    throat-pair representation for the polarization sum.
  - **(C) No structural match.** BAM amplitude has neither correct
    poles nor correct angular dependence. Would falsify the
    amplitude-level reinterpretation thread.

## Concrete construction

### Mandelstam → ψ mapping

Mapping the s-channel propagator pole to the `S³` Green-function
pole: at the pole, `s = m²` should correspond to `ψ_s = 0` (the
S³ Green function diverges like `1/ψ`). The natural choice is

    ψ_s(s)  =  (s − m²) / (2 m²)

so that near threshold

    G_S3(ψ_s)  ≈  1 / (4π·R·ψ_s)  =  2 m² / (4π·R·(s − m²))

reproducing the propagator pole `1/(s − m²)` up to a constant.
Similarly `ψ_u(u) = (u − m²) / (2m²)`.

The s and u channel mapping is not symmetric in sign — at the
physical point `u < m²` typically, so `ψ_u < 0` formally. We take
the absolute value `|ψ_u|` and absorb the sign into the channel
phase. This is the simplest natural choice; alternatives (e.g.
`ψ = arctan(...)`) are tested as sensitivity probes.

### Throat transport per vertex

Each scattering event involves two throat traversals (one at the
initial vertex, one at the final vertex). Each contributes a factor
of `T = iσ_y`. In the scalar amplitude (treating the spin trace
externally), this is a phase factor `i` per traversal, so `T·T`
contributes `i² = −1` to the channel amplitude.

### Closure phase

The closure phase per channel has two contributions:

  - **Hopf holonomy** at the electron's locked Hopf angle `χ = 0`:
    `π·cos(0) = π`, giving a phase factor `exp(iπ) = −1`.
  - **SU(2) spin transport** along the geodesic of scattering angle
    `θ` for spin-½: `α_spin · θ = θ/2`.

For the s-channel, the spin transport runs forward; for the u-channel
it runs backward (the electron's worldline crosses itself in the
crossed diagram). So
    `φ_s = π + θ/2`,
    `φ_u = π − θ/2`,
giving relative phase factor `exp(+iθ/2)` vs `exp(−iθ/2)` between
the two channels.

### Total amplitude

    M_BAM(θ) = G_S3(ψ_s) · exp(i·φ_s) · T² · η_s
             + G_S3(ψ_u) · exp(i·φ_u) · T² · η_u

with `η_s, η_u` channel signs (the u-channel formally has `ψ_u`
negative, accounting for the crossing sign).

## Tests

### T1. Propagator-pole reproduction

Verify that, in the on-shell limits

    s → m²:  |M_BAM| → ∞  with leading divergence `1/(s − m²)`;
    u → m²:  same with `1/(u − m²)`,

the BAM amplitude exhibits the correct propagator-pole structure
inherited from `G_S3(ψ) ∼ 1/ψ`. Fit the leading divergence and
report the residue.

### T2. Thomson-limit angular dependence

In the limit `ω → 0` at fixed `θ`, evaluate `|M_BAM(θ)|²` over
`θ ∈ [0, π]`. Normalize and compare to the Klein-Nishina /
Thomson prediction

    F_KN(θ)  =  (1 + cos²θ) / 2.

Compute the residual

    R(θ)  =  |M_BAM(θ)|² / |M_BAM(0)|²  −  F_KN(θ) / F_KN(0).

Fit a one-parameter family `(1 + a·cos²θ + b·cos θ)` to
`|M_BAM(θ)|²` and report the fit coefficients. Klein-Nishina has
`a = 1, b = 0`. Pure `(1 + cos θ)` would give `a = 0, b = 1` (a
forward-peaked single-cosine angular pattern).

### T3. Ansatz sensitivity

Vary the spin-phase ansatz: try
  - `φ_s = θ/2, φ_u = −θ/2`  (current natural choice for spin-½),
  - `φ_s = θ, φ_u = −θ`      (spin-1 lift, photon-like phase),
  - `φ_s = 0, φ_u = θ`       (only u-channel carries θ-dependence),

and report which (if any) reproduces (1 + cos²θ).

### T4. Photon-throat-pair effect

If the photon is given an explicit throat-pair representation (its
own `T_γ`), an additional factor enters the amplitude. As a
diagnostic, try multiplying by an extra spin-1 factor
`exp(i·1·θ_transport) = exp(iθ)` and check whether this restores
Klein-Nishina structure. If yes, it identifies the missing piece:
photons need explicit throat-pair representation in BAM.

## Predicted verdict structure

The honest expected outcome is **(B) partial match**:

  - T1 PASS — propagator poles are reproduced. The structural
    statement "the antipodal `S³` Green function plays the role of
    the QED propagator" survives at the kinematic-pole level.

  - T2 likely FAIL with the spin-½-only construction — the BAM
    amplitude gives a single-cosine `(1 + cos θ)` angular pattern,
    not the spin-1 polarization-summed `(1 + cos²θ)`.

  - T4 may identify the missing structure as the photon's
    spin-1 internal phase, currently absent from BAM. This is the
    follow-on target: give the photon explicit BAM machinery
    (Hopf-fibre excitation, throat-pair, or vector connection on
    `S³`) and re-run the angular test.

The verdict line in the probe captures whether the match is full,
partial, or absent, and identifies the next structural piece needed.

## Cross-references

- `experiments/closure_ledger/compton_antipodal_kinematics_probe.py`
  — the kinematic predecessor.
- `geometrodynamics/transaction/handshake.py` — existing
  amplitude algebra (offer × confirm × phase-closure) on `S³`.
  The probe builds a focused Compton-amplitude version using the
  same ingredients but tailored to scattering rather than
  GW-triggered transactions.
- `geometrodynamics/transaction/s3_geometry.py` — `s3_green_potential`
  is the natural antipodal propagator.
- `geometrodynamics/hopf/connection.py` — Hopf holonomy
  `π·cos(χ)` provides the closure phase contribution.
- `geometrodynamics/embedding/transport.py` — `T = iσ_y` throat
  transport.
- `experiments/closure_ledger/compton_amplitude_structure_probe.py`
  — first probe in this sub-thread.
