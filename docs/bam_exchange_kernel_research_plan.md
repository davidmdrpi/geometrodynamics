# BAM exchange-kernel derivation — research plan

Closes the last identified gap from the Bhabha/Møller derivation
thread (PRs #42–#44): the virtual-photon propagator `1/q²`. PRs
#42–#44 derived the full 4-fermion tree QED scalar intensities for
Bhabha and Møller from BAM SU(2) Hopf-bundle spinor traces (diagonal
numerators) and the non-orientable `T = iσ_y` throat transport
(Fermi-statistics signs), but the propagator was still the QED
ansatz `1/q²`. This probe derives `1/q²` directly from BAM
throat-fibre exchange geometry — **without invoking virtual photons**.

## The conceptual shift

In QED, `1/q²` is the propagator of an off-shell intermediate
photon — a "virtual particle" with momentum `q` that is not on the
mass shell `q² = 0`. The notion of virtual particles is a
calculation tool; their physical existence is a contentious
philosophical question.

In BAM, there are **no virtual photons**. Two fermions on `S³` at
points `x_A` and `x_B` interact via a **throat-fibre exchange**:
the throat itself carries the gauge structure (Hopf-bundle U(1) for
electromagnetism), and the amplitude for the exchange is a
**Green function on `S³`** — a purely geometric kernel.

The claim under test: the S³ scalar Green function (already in the
repo, `geometrodynamics.transaction.s3_geometry.s3_green_potential`)

```
G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)
```

reduces in the flat-space limit (ψ → 0, R → ∞ with `d = ψ·R` fixed)
to the Coulomb potential

```
G(d) → 1 / (4π d)
```

whose 3D Fourier transform is **exactly** `1/q²`, the QED photon
propagator in Feynman gauge — **without any virtual particle**.

## Tests

  T1. **S³ Green function from repo**: load `s3_green_potential(ψ)`
      and verify the form `((π-ψ)cot ψ - ½)/(4π²R)`.

  T2. **Flat-space limit** (Coulomb potential): for small ψ,
      `G(ψ) = 1/(4π d) + zero-mean offset + O(d/R²)`. Verify
      numerically that the singular part matches `1/(4π d)`.

  T3. **Fourier transform → 1/q²**: numerically compute the 3D
      Fourier transform of `G(d) = 1/(4π d)` and verify it equals
      `1/q²` (Coulomb-to-momentum-space identity).

  T4. **Spectral representation on S³**: write `G(ψ)` as a sum
      over S³ scalar harmonics `Y_{n,l,m}(ψ, θ, φ)` with Laplacian
      eigenvalues `λ_n = n(n+2)/R²`. Verify the spectral sum
      reproduces the position-space `G(ψ)`.

  T5. **Curvature corrections**: identify the `O(1/R²)` corrections
      to the BAM propagator relative to QED `1/q²`. Verify they
      vanish in the flat-space limit and are suppressed in the
      kinematic regime where flat-space QED is tested.

  T6. **End-to-end Bhabha**: combine PR #43 SU(2) Pauli-trace
      diagonals + PR #44 Möbius interference sign + this probe's
      `1/q²` propagator. Verify `|M̄|²_Bhabha/(8e⁴)` matches QED
      textbook formula to machine precision.

  T7. **End-to-end Møller**: same combination, verify Møller.

  T8. **Verdict**: BAM throat-fibre exchange recovers QED `1/q²`
      from geometric Green function, no virtual photons needed.

## Verdict structure

  - **PROPAGATOR_FROM_GEOMETRY**: T1–T7 all pass. The QED photon
    propagator `1/q²` is the flat-space limit of the BAM S³ Green
    function — a geometric kernel for two-point exchange, with no
    notion of off-shell virtual particles. Combined with PRs #43
    (Pauli traces) and #44 (Möbius signs), all of tree-level QED
    Bhabha and Møller derives from BAM geometric ingredients alone.

  - **PROPAGATOR_GAP**: any test fails. Indicates a non-trivial
    mismatch between the BAM S³ Green function and the QED photon
    propagator that needs further investigation.

## What this closes

If all tests pass, the BAM derivation of tree-level QED two-to-two
scattering processes is **complete**:

  - **Compton/BW/annihilation** (PRs #35–#41): closed-form F²
    from BAM throat action.
  - **Bhabha/Møller diagonals** (PR #43): SU(2) Pauli traces.
  - **Bhabha/Møller interference signs** (PR #44): T = iσ_y =
    antisymmetric ε.
  - **Propagator structure 1/q²** (this PR): S³ Green function in
    flat-space limit.

No remaining "QED overlay" — every ingredient derives from
BAM-geometric primitives.

## What this leaves open

  - **Loop corrections**: tree-level only. Vertex/self-energy/
    vacuum polarisation would require integration over closed
    throat-fibre loops on S³, coupling to the bulk radial channel.
  - **Hopf-bundle photon vs scalar Green function**: this probe
    uses the scalar Green function. Strictly speaking, the QED
    photon has spin-1; the Hopf-bundle generalisation requires
    accounting for the photon's vector structure. In Feynman gauge,
    however, the vector-photon propagator is `−g^{μν}/q²`, and
    the structure reduces to the scalar case under the gauge
    choice. The Hopf-fibre extension is a follow-on.
  - **Curvature corrections at high energy**: as `q² ∼ 1/R²`,
    the BAM propagator deviates from `1/q²` by S³ curvature
    effects. For laboratory energies and `R ∼ Hubble scale`,
    corrections are negligible; at Planck scale they are not.

## Cross-references

  - PR #42: identified the propagator gap.
  - PR #43: SU(2) Pauli traces for diagonal numerators.
  - PR #44: Möbius `T = iσ_y` interference signs.
  - `geometrodynamics/transaction/s3_geometry.py:s3_green_potential`
    — the S³ Green function used here.
  - `experiments/closure_ledger/bam_exchange_kernel_probe.py`
    — this probe.
