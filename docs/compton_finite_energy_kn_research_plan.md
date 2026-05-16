# Finite-energy Klein-Nishina recoil — probe plan

Follow-on to the Thomson-limit Klein-Nishina restoration probe
(PR #28, merged). That probe established that the BAM Compton
amplitude with explicit photon transverse polarization and scalar
electron reproduces the Klein-Nishina angular factor
`(1 + cos²θ)/2` exactly at ω → 0. The natural next test: does BAM
also reproduce the **finite-energy** Klein-Nishina recoil structure,
or does the Thomson-limit match conceal a leading-order failure at
ω ~ m_e?

## The finite-energy KN amplitude

The full unpolarized Compton differential cross section in the
electron rest frame is

    dσ/dΩ  ∝  (ω'/ω)² · [ω'/ω + ω/ω' − sin²θ]

with the recoil relation `ω' = ω / (1 + (ω/m)·(1 − cos θ))`. Let
`x ≡ ω'/ω`. The squared-amplitude angular structure is

    |M_KN|²(x, θ)  ∝  x² · (x + 1/x − sin²θ)
                    =  x³ + x − x² sin²θ.

At Thomson (x → 1, ω → 0): `|M_KN|² → 1 + 1 − sin²θ = 1 + cos²θ`.
Two finite-ω features that distinguish KN from a pure (1 + cos²θ):

  - **Forward-backward asymmetry from the (ω'/ω)² factor.** At
    θ = π (backscatter), x = 1/(1 + 2ω/m) < 1, so the cross section
    is suppressed by `x²`. At θ = 0, x = 1 (no recoil) and the cross
    section is unchanged. The result is a forward peak that grows
    with ω/m.
  - **Compton edge.** The maximum photon energy loss occurs at
    θ = π and equals `2ω²/(m + 2ω)`. At ω ~ m, the backscattered
    photon energy collapses to `m/2` regardless of incoming ω.

## What the natural BAM construction predicts

The BAM amplitude from the photon-structure probe (scalar electron,
photon polarization summed) gives

    |M_BAM|²(ω, θ)  ∝  (1 + cos²θ) · (G_s + G_u)²
                    ∝  (1 + cos²θ) · (1/ψ_s + 1/ψ_u)²
                    =  (1 + cos²θ) · (m/ω + m/ω')²
                    =  (1 + cos²θ) · (m/ω)² · (1 + 1/x)²

Normalising at θ = 0 (where x = 1 for all ω) gives

    f_BAM(x, θ)  ≡  |M_BAM|² / |M_BAM(θ=0)|²
                 =  ((x + 1)² / (4 x²)) · ((1 + cos²θ) / 2).

The KN normalised form

    f_KN(x, θ)  ≡  |M_KN|² / |M_KN(θ=0)|²
                 =  (x³ + x − x² sin²θ) / 2.

At Thomson (x → 1) both reduce to `(1 + cos²θ)/2` exactly. At finite
ω they differ:

    f_BAM − f_KN  ≠  0   for x ≠ 1.

## Predictions

### T1. Thomson sanity

At ω/m_e = 1e-4: `|f_BAM(x, θ) − f_KN(x, θ)| < 1e-3` for all sampled
θ. Reproduces PR #28's main result; sanity check that the new
finite-energy machinery still recovers Thomson.

### T2. Leading-order discrepancy

Expand `f_BAM − f_KN` in ε = ω/m_e at fixed θ. Compute the leading
non-vanishing coefficient at each θ. Predicted structural result:
the O(ε) discrepancy is non-zero at θ ≠ 0, π, indicating the natural
BAM construction reproduces KN at zeroth order in ω/m but disagrees
at first order in ω/m.

### T3. Forward-backward asymmetry at finite ω

KN predicts a specific forward peak that grows with ω. Compute the
asymmetry `A(ω) = (σ(θ=0) − σ(θ=π)) / (σ(θ=0) + σ(θ=π))` for both
KN and BAM at several ω values. **Predicted:** KN gives a non-zero,
ω-dependent asymmetry; BAM (with x = 1 at θ=0 and x = 1/(1+2ε) at
θ=π) gives a different ω-dependence.

### T4. Compton edge at backscatter

KN at θ = π: `f_KN(x, π) = (x³ + x)/2` with `x = 1/(1 + 2ε)`. BAM at
θ = π: `f_BAM(x, π) = (x + 1)²/(4 x²)`. Plot both as functions of ε
on a log scale and identify the regime where they diverge by >10%.

### T5. Full angular fit at finite ω

At ω/m_e = 0.5 (relativistic Compton regime), fit `|M_BAM|²(θ)` and
`|M_KN|²(θ)` to the same angular polynomial basis. Compare fit
coefficients. **Predicted:** BAM and KN disagree in their c_2 / c_0
ratios at this energy, with the difference traceable to the missing
"vertex factor" structure (ε·k terms in QED that BAM doesn't have).

## Expected verdict

**PARTIAL_MATCH at finite energy.** The Thomson limit (PR #28) is
recovered cleanly; the natural BAM construction fails to reproduce
finite-ω KN at leading order in ε = ω/m_e. The discrepancy is
structural: BAM's per-channel propagator sum `(1/ψ_s + 1/ψ_u)`
gives `(x + 1)²/(4x²)` energy dependence, while QED's vertex-factor
algebra gives `x²(x + 1/x)`. The missing BAM ingredient is the
**per-channel kinematic weighting** — in QED, the s- and u-channel
amplitudes are weighted by momentum-dependent vertex factors (ε·k
type contractions) that the natural BAM construction lacks.

The probe's value is not in claiming a full KN reproduction but in
**locating the precise structural gap** at finite ω. With this
identified, the natural follow-on is to test whether a BAM-derived
vertex-factor structure (e.g. from explicit Hopf-connection coupling
in `hopf/connection.py`) closes the gap, mirroring how the
photon-polarization probe closed the Thomson-limit gap.

## Cross-references

- `experiments/closure_ledger/compton_photon_structure_probe.py`
  — predecessor establishing Thomson-limit KN match (PR #28).
- `experiments/closure_ledger/compton_amplitude_structure_probe.py`
  — partial-match precursor (PR #26).
- `experiments/closure_ledger/compton_antipodal_kinematics_probe.py`
  — kinematic foundation (PR #25).
- `experiments/closure_ledger/compton_finite_energy_kn_probe.py`
  — this probe.
