# Photon-structure probe — research plan

Follow-on to the Compton amplitude structure probe (PR #26, merged).
That probe established a **partial structural match**: BAM's
antipodal `S³` Green function with `ψ ∝ (s − m²)` reproduces the QED
propagator pole exactly (fitted exponent 1.0002), but the natural
spin-½ ansatz gives Thomson angular dependence `(1 + cos θ)/2`, not
the Klein-Nishina `(1 + cos²θ)/2`. The diagnostic T4 confirmed that
an overall photon phase factor `exp(iθ)` cannot restore KN — a
polarization-vector structure is needed.

This probe gives the photon explicit polarization machinery and tests
whether the structural gap closes.

## The missing structure

The QED Thomson angular factor `(1 + cos²θ)` comes from the
polarization sum

```
Σ_{λ, λ'} | ε^(λ)(k) · ε^(λ')*(k') |²  =  1 + cos²θ
```

over two transverse polarisations per photon (the completeness
relation `Σ_λ ε^(λ)_i ε^(λ)_j = δ_ij − k̂_i k̂_j` projects onto the
plane transverse to propagation). BAM's previous amplitude treats
the photon as a phase factor with no vector structure, so this sum
is absent and `(1 + cos²θ)` cannot arise.

## What the photon-structure probe tests

Build the polarization-resolved BAM amplitude

```
M^{(λ, λ')}_BAM(θ, ω)  =  ( ε^(λ)(k) · ε^(λ')*(k') )
                          ·  [ G_S3(ψ_s) · exp(iφ_s)
                             + G_S3(ψ_u) · exp(iφ_u) ]
                          ·  T²
```

with two photon polarisations per particle. The natural choice for
the polarization-vector basis: pick the scattering plane (e.g. xz),
let `ε^(1)` lie in the plane (the "parallel" polarization), `ε^(2)`
perpendicular to it. Outgoing photon at angle θ has its `ε'^(1)`
rotated by θ in the plane and `ε'^(2)` unchanged.

Then

  - `ε^(1)(k) · ε^(1)'*(k') = cos θ`
  - `ε^(2)(k) · ε^(2)'*(k') = 1`
  - mixed `(1, 2)` and `(2, 1)` give 0
  - sum of squared: `cos²θ + 1 = 1 + cos²θ`  ✓

The polarization-summed cross section is

```
|M_total|² = Σ_{λ, λ'} | M^{(λ, λ')}_BAM |²
```

## Predictions

### T1. Polarization sum identity

Sanity check: `Σ_{λ, λ'} | ε^(λ)(k) · ε^(λ')*(k') |² = 1 + cos²θ`
exactly across all sampled angles. Verifies the polarization-vector
implementation.

### T2. Klein-Nishina restored with photon polarization + scalar electron

Build `M^{(λ, λ')}_BAM` with:
  - photon polarization factor `ε · ε'*`,
  - propagator factor `G_s + G_u`,
  - scalar electron (no θ-dependent spin phase: φ_s = φ_u = π),
  - throat factor T² = −I.

Compute `|M_total|²` summed over polarisations, fit
`A + B cos θ + C cos²θ` in the Thomson limit. **Predicted:**
`(A, B, C) ∝ (1, 0, 1)` — exact Klein-Nishina.

### T3. With electron spin-½ phases re-included as a comparison

Repeat T2 with the previous probe's spin-½ ansatz
`φ_s = π + θ/2, φ_u = π − θ/2`. **Predicted:** the
`(1 + cos θ)` factor from the spin-½ phases multiplies the photon
polarization sum, giving `(1 + cos θ)(1 + cos²θ)` — a four-term
angular polynomial that is not Klein-Nishina. This identifies the
spin-½ phases as a spurious contribution at the Thomson limit, where
electron spin does not enter the leading QED amplitude.

### T4. Propagator pole structure robust to new photon structure

Repeat the T1 propagator-pole test from PR #26 with the new
amplitude. At a fixed off-axis angle, the polarization factor is just
a constant, so the `1/(s − m²)` pole should remain reproduced.
**Predicted:** PASS, pole exponent ≈ 1.0.

### T5. Polarization-basis invariance

Replace linear polarisations with circular polarisations
`ε^(R) = (ε^(1) + i ε^(2))/√2`, `ε^(L) = (ε^(1) − i ε^(2))/√2`.
**Predicted:** total `|M_total|²` is invariant under polarization
basis choice (basis-completeness).

## Expected verdict

**FULL_MATCH** if T1, T2, T4, T5 all pass and T3 produces the
predicted spurious convolution. The Compton amplitude's tree-level
structural skeleton would then be reproduced from the BAM
ingredient set

  - antipodal `S³` propagation  (propagator pole)
  - photon transverse polarization in tangent bundle  (KN angular factor)
  - throat transport T = iσ_y  (overall sign / closure)
  - Hopf-holonomy closure phase  (the χ = 0 lock)

with the electron treated as a scalar charge at the Thomson limit
(spin effects appearing only at higher orders in ω/m_e).

If T2 fails to give Klein-Nishina, the structural gap is deeper than
polarization machinery and the thread closes with a deeper
falsification.

## What this leaves open even if T2 passes

  - **High-energy Klein-Nishina recoil** `(ω'/ω)² · [ω'/ω + ω/ω' −
    sin²θ]` for `ω ~ m_e` — a separate probe.
  - **Loop corrections** — vertex, self-energy, vacuum polarization;
    would require BAM's bulk radial channel.
  - **General-frame Lorentz covariance** of the antipodal +
    polarization construction.
  - **Pair production, electron self-energy** — other QED events to
    test the same machinery against.

## Cross-references

- `experiments/closure_ledger/compton_amplitude_structure_probe.py`
  — the partial-match predecessor (PR #26).
- `experiments/closure_ledger/compton_antipodal_kinematics_probe.py`
  — kinematic foundation (PR #25).
- `geometrodynamics/transaction/s3_geometry.py` — `s3_green_potential`.
- `experiments/closure_ledger/compton_photon_structure_probe.py`
  — this probe.
