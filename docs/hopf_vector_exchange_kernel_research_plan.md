# BAM Hopf/vector exchange-kernel — research plan

Natural follow-on from PR #45 (BAM scalar exchange kernel → `1/q²`).
PR #45 derived the QED photon propagator structure `1/q²` from the
flat-space limit of the **scalar** S³ Green function. That worked for
spin-summed Bhabha/Møller amplitudes because fermion currents are
conserved, so the gauge-fixing structure of the photon propagator
drops out and only the scalar `1/q²` factor survives.

This probe targets the **full Lorentz-vector structure** of the
photon propagator: in Feynman gauge `D_F^{μν}(q) = −g^{μν}/q²`, in
Lorenz gauge `D_L^{μν}(q) = (−g^{μν} + (1−ξ)q^μq^ν/q²)/q²`, in
transverse gauge `D_T^{μν}(q) = P_T^{μν}(q)/q²` with
`P_T^{μν} = −g^{μν} + q^μq^ν/q²`. These differ by gauge-mode
contributions that decouple from conserved-current amplitudes.

BAM's Hopf bundle is the natural geometric structure carrying U(1)
gauge degrees of freedom. The photon is a Hopf-fibre connection (PR #38
T4 verified that the spin-1 photon helicity sum
`(1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2) = Σ_λ|d¹_{1,λ}(θ)|²` is the
Hopf-bundle helicity transport). This probe constructs the
Hopf-bundle vector Green function and shows:

  - It factors as `D^{μν}(q) = (Lorentz tensor structure)·D_scalar(q²)`
    where `D_scalar(q²) = 1/q²` is the flat-limit of PR #45.
  - In Feynman gauge the Lorentz tensor is `−g^{μν}`.
  - The Hopf-bundle naturally selects the 2 transverse polarizations
    (the photon is helicity-±1, not 0).
  - Spin-summed Bhabha/Møller via vector exchange reduce to PR #45's
    scalar exchange (gauge-equivalent via Ward identity).

## What this resolves

  - **Hopf-bundle vector kernel**: explicitly identified as the
    photon propagator structure, with the scalar PR #45 result as
    the gauge-equivalent simplification.
  - **Gauge invariance of physical amplitudes**: verified via the
    Ward identity at the fermion-current level
    `q_μ Tr[γ^μ p̸_1 γ^ν p̸_2] = 0` for `q = p_2 − p_1`.
  - **Photon polarization**: the 2-polarization structure of the
    photon (transverse helicity-±1) is the Hopf-fibre helicity
    structure from PR #38, not a separate gauge-fixing requirement.

## Tests

  T1. **Vector propagator setup**: define Feynman, Lorenz, transverse
      vector propagator forms; verify the standard Lorentz tensor
      structure.

  T2. **Feynman-gauge factorization**: verify that the Feynman vector
      propagator factors as `D_F^{μν}(q) = −g^{μν}·D_scalar(q²)`
      with `D_scalar(q²) = 1/q²` (the PR #45 result).

  T3. **Transverse projector decomposition**:
      `P_T^{μν}(q) = −g^{μν} + q^μq^ν/q²`. Verify properties:
        - Idempotent: `P_T · P_T = P_T` (projector)
        - Transverse: `q_μ P_T^{μν} = 0`
        - Trace: `Tr[P_T] = 3` (3 transverse modes in 4D, including
          longitudinal-time mode; physical photon has 2 of these
          after gauge fixing)

  T4. **Ward identity (current conservation)**: for a conserved
      fermion current `J^μ ∝ ū(p_2)γ^μ u(p_1)` with massless on-shell
      fermions, `q_μ J^μ = 0` where `q = p_2 − p_1`. Verified at the
      spin-summed level via `q_μ q_ν Tr[γ^μ p̸_1 γ^ν p̸_2] = 0`.

  T5. **Gauge equivalence for physical amplitudes**: contract the
      Feynman and transverse propagators with two conserved fermion
      currents; verify the resulting physical amplitudes agree
      (Ward identity kills the gauge difference).

  T6. **Hopf-bundle 2-polarization structure**: the photon has 2
      physical helicity states (Hopf-fibre helicity ±1). Recap the
      PR #38 T4 result: `(1+c²)/2 = Σ_λ|d¹_{1,λ}(θ)|²` gives the
      Hopf-fibre helicity transport sum.

  T7. **End-to-end Bhabha/Møller via vector exchange**: contract
      fermion currents with the Feynman-gauge vector propagator;
      verify the spin-summed `|M̄|²` matches PR #45 and QED.

  T8. **S³ curvature corrections (vector case)**: the gauge-covariant
      Laplacian on S³ has a Ricci mass term `R_μν A^ν ~ (2/R²) A_μ`
      from the constant Ricci curvature `R_μν = (2/R²)g_μν`. In the
      flat limit `R → ∞` this vanishes. Verify the leading
      curvature correction is O(1/R²) and matches the analytic
      prediction.

## Verdict structure

  - **VECTOR_KERNEL_FROM_HOPF_BUNDLE**: all tests pass. The
    Hopf-bundle vector Green function in Feynman gauge factors as
    `−g^{μν}·G_scalar(q²)`; in the flat limit `G_scalar → 1/q²`
    (PR #45). The transverse projector `P_T^{μν}` and the
    polarization sum from PR #38 (Hopf-fibre helicity) together
    describe the 2 physical photon polarizations. Gauge differences
    between Feynman / Lorenz / transverse drop from physical
    amplitudes via the Ward identity. End-to-end Bhabha/Møller via
    vector exchange match QED at machine precision.

  - **STRUCTURAL_GAP**: a non-trivial mismatch between the BAM
    Hopf-bundle vector propagator and the QED photon propagator
    needs investigation.

## What this closes vs leaves open

**Closes**:
  - The full Lorentz-vector photon propagator structure is now
    Hopf-bundle-derivable, not just the scalar gauge-equivalent
    factor.
  - The Hopf-bundle 2-helicity polarization structure (PR #38) is
    explicitly connected to the gauge-fixing of the photon
    propagator.

**Leaves open**:
  - **Loop corrections**: still tree-level only.
  - **Non-Abelian extension**: BAM's QCD program (separate thread)
    would need a Hopf-bundle structure on a non-Abelian group;
    not addressed here.
  - **Off-shell virtual modes**: the Ward identity argument works
    on-shell; off-shell virtual photon modes (in loops, for example)
    would require explicit Faddeev-Popov ghost analysis on S³.

## Cross-references

  - PR #38: `throat_nucleation_caustic_derivation_probe.py` — Hopf
    helicity sum `(1+c²)/2 = Σ|d¹_{1,λ}|²`.
  - PR #43: `dirac_trace_geometry_probe.py` — SU(2) Pauli traces.
  - PR #44: `mobius_exchange_sign_probe.py` — `T = iσ_y = ε`.
  - PR #45: `bam_exchange_kernel_probe.py` — scalar Green function → `1/q²`.
  - `geometrodynamics/hopf/connection.py` — Hopf connection
    `A_φ(χ) = ½cos(χ)`.
  - `geometrodynamics/transaction/s3_geometry.py:s3_green_potential` —
    scalar Green function.
  - `experiments/closure_ledger/hopf_vector_exchange_kernel_probe.py`
    — this probe.
