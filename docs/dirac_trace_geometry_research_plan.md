# BAM Dirac-trace geometry probe — research plan

Follow-on to PR #42 (Bhabha/Møller two-channel interference): PR #42
identified that the BAM scalar-intensity Compton kernel cannot
reproduce the QED 4-fermion structures `(s²+u²)`, `(u²+t²)`,
`(s²+t²)` because they come from Dirac γ-matrix trace algebra. This
probe tests whether those structures emerge from **geometric traces
over BAM throat spinor transport** — specifically, the SU(2) Pauli
matrices on the Hopf bundle (from `geometrodynamics.hopf.spinor`,
`geometrodynamics.embedding.transport`).

## The claim under test

BAM's natural spinor representation — SU(2) Pauli matrices on the
Hopf-fibre — is the standard chiral (Weyl) decomposition of QED
Dirac spinors. The parity-symmetric Weyl trace product

```
T_BAM(p_a, p_b, p_c, p_d)
  =  Σ_{μν} η_μ η_ν · (2·Re[Tr_2x2[σ^μ σ̄·p_a σ^ν σ̄·p_b]])
                    · (2·Re[Tr_2x2[σ^μ σ̄·p_c σ^ν σ̄·p_d]])
```

equals the QED Dirac trace product
`Tr_4D[γ^μ p̸_a γ^ν p̸_b] · Tr_4D[γ_μ p̸_c γ_ν p̸_d]
   = 32·[(p_a·p_c)(p_b·p_d) + (p_a·p_d)(p_b·p_c)]`

The 4-fermion numerator structures emerge as **specific kinematic
pairings** of the four momenta:

  - Bhabha t-channel `(p_1, p_3)(p_2, p_4)` → `8·(s² + u²)`
  - Bhabha s-channel `(p_1, p_2)(p_3, p_4)` → `8·(u² + t²)`
  - Møller t-channel `(p_1, p_3)(p_2, p_4)` → `8·(s² + u²)`
  - Møller u-channel `(p_1, p_4)(p_2, p_3)` → `8·(s² + t²)`

## What this resolves vs leaves open

**Resolves** (if successful): the diagonal-channel numerators
`(s²+u²)`, `(u²+t²)`, `(s²+t²)` from PR #42 are reproduced by BAM
SU(2) Hopf-bundle spinor traces. The "Dirac-trace structure" is not
a foreign QED imposition — it's the natural BAM spinor algebra.

**Leaves open**: the interference sign / magnitude `2u²/(s·t)` for
Bhabha and `2s²/(t·u)` for Møller. These come from a **single**
8-γ trace combining all four momenta, with the relative sign
determined by Fermi statistics (Wick contractions for Bhabha;
Pauli antisymmetrisation for Møller).

The probe tests:

  - **Diagonal traces** (likely all PASS via the Pauli trace identity).
  - **Interference single trace** (likely matches magnitude but
    sign requires additional Fermi-statistics input).

## Tests

  T1. **Pauli σ matrices and slashed momenta**: setup the BAM
      spinor objects and verify the basic identity
      `Tr_2x2[σ^μ σ̄^ν] = 2·η^{μν}`.

  T2. **Parity-symmetric Weyl trace product formula**: numerically
      verify that the parity-symmetric Pauli trace product equals
      the QED Dirac trace `32·[(A·C)(B·D) + (A·D)(B·C)]` for
      arbitrary 4-momenta.

  T3. **Bhabha t-channel diagonal**: at sample CM angles, verify
      `T_BAM(p_1, p_3, p_2, p_4) = 8·(s² + u²)` to machine precision.

  T4. **Bhabha s-channel diagonal**: verify
      `T_BAM(p_1, p_2, p_3, p_4) = 8·(u² + t²)`.

  T5. **Møller t-channel diagonal**: verify
      `T_BAM(p_1, p_3, p_2, p_4) = 8·(s² + u²)`.

  T6. **Møller u-channel diagonal**: verify
      `T_BAM(p_1, p_4, p_2, p_3) = 8·(s² + t²)`.

  T7. **Bhabha interference single trace**: compute the BAM
      single-trace structure with all four momenta inserted.
      Compare magnitude to QED `8·u²` and identify the sign source.

  T8. **Møller interference single trace**: similar for Møller,
      with QED `8·s²` magnitude and Pauli sign.

  T9. **End-to-end |M̄|² reconstruction**: combine diagonal Pauli
      traces with the interference trace, verify the full QED
      Bhabha and Møller scalar intensities are reproduced (modulo
      the relative interference sign, which requires explicit
      Fermi/Pauli input).

  T10. **Verdict**: identify which QED structures are derivable
       from BAM Pauli/Hopf-bundle spinor traces alone, and which
       still require an explicit Fermi-statistics rule.

## Verdict structure

  - **DIAGONALS_DERIVED_INTERFERENCE_MAGNITUDE_OK**: T1–T8 all
    pass; T9 reproduces QED magnitudes; the interference sign is
    the only remaining input (and it is well-defined from Fermi
    statistics, not arbitrary). The BAM SU(2) Hopf-bundle spinor
    traces are sufficient for all 4-fermion kinematic structures
    in QED Bhabha/Møller.

  - **DIAGONALS_DERIVED_INTERFERENCE_FAILS**: diagonal traces match
    QED but the interference single trace doesn't match — would
    indicate a non-trivial mismatch between BAM Pauli traces and
    QED 8-γ Dirac traces.

  - **DIAGONALS_FAIL**: even the parity-symmetric Pauli trace
    doesn't reproduce QED `(s²+u²)` etc. — would mean BAM's
    SU(2) Hopf-bundle spinor structure is NOT the right spinor
    representation for QED.

## What this leaves open

  - **Fermi-statistics sign rule**: the relative sign between
    interfering diagrams (Bhabha s↔t, Møller t↔u) comes from Wick
    contractions or Pauli antisymmetrisation. BAM's `T = iσ_y`
    non-orientable throat transport (README channel 2) is a
    candidate origin for the spin-½ antisymmetrisation; testing
    whether it predicts the correct interference signs is a
    follow-on Möbius-sign probe.

  - **Virtual-photon propagator**: still `1/q²` ansatz; deriving
    it from BAM throat-fibre dynamics is a separate probe.

  - **Loop corrections**: tree-level only.

## Cross-references

  - PR #35 / PR #38–#41: Compton/F² derivation thread.
  - PR #42: `bhabha_moller_interference_probe.py` — identified the
    Dirac-trace gap.
  - `geometrodynamics/hopf/spinor.py` — SU(2) spinor transport.
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y`
    antipodal throat transport (non-orientable, Pauli-sign origin).
  - `experiments/closure_ledger/dirac_trace_geometry_probe.py`
    — this probe.
