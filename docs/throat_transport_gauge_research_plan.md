# Gauge-sensitive throat-transport falsifier — research plan

Opens the natural follow-up thread to the moving-mouth Berry-phase
test (PR #23, merged f0d4b0d). That test established that the BAM
symmetric Hopf connection `A_BAM = ½cos(χ)dφ` is internally
consistent: every closed loop in the (χ,φ) base accumulates exactly
`∮A_BAM` as its Berry phase, both via the line integral of A and via
the link-variable formula on the explicit BAM spinor
`(cos(χ/2)e^{−iφ/2}, sin(χ/2)e^{+iφ/2})`. The "What this leaves open"
section of that probe flagged a **gauge-experimental discrimination**
question, which is the subject of this thread.

## The question

The BAM and Bloch gauges differ by

    A_BAM − A_Bloch  =  ½ dφ

i.e. they are related by the (multi-valued) U(1) gauge function
`Λ(χ, φ) = φ/2`. Under any 2π φ-loop, `Λ → Λ + π`, so the gauge
transformation `g = exp(iΛ)` returns to `−1`, not `+1`. The two gauges
therefore differ by a Wu-Yang transition function on the equatorial
overlap.

The throat transport `T = iσ_y` is derived in `embedding/transport.py`
in the **BAM-frame** spinor representation. The natural question is
whether the combination

    (BAM connection A_BAM)  +  (BAM-frame throat T_BAM = iσ_y)

predicts a physically distinct observable from the gauge-conjugate
combination

    (Bloch connection A_Bloch)  +  (Bloch-frame throat T_Bloch)

where `T_Bloch = g^{−1}(φ_out) · T_BAM · g(φ_in)` with the *scalar*
multi-valued U(1) gauge function `g(χ,φ) = e^{−iφ/2}`. (Pulling this
common phase out of `ψ_BAM = (cos(χ/2)e^{−iφ/2}, sin(χ/2)e^{+iφ/2})`
recovers `ψ_Bloch = (cos(χ/2), sin(χ/2)e^{+iφ})`, so the relation
between the two spinor parametrisations is a single scalar phase,
not an SU(2) frame rotation.)

If yes → BAM is making a sharper-than-Bloch prediction; an experiment
can in principle distinguish them.

If no → the choice between gauges is a computational convenience; BAM's
"symmetric gauge + T = iσ_y" is a particular section of the same
SU(2) bundle that Bloch describes asymmetrically. The falsifier returns
the verdict that BAM survives this test by being gauge-equivalent to
the textbook treatment, not by predicting something extra.

Either outcome is informative: it sharpens what BAM's geometric
content actually buys, beyond just internal consistency.

## What "all allowed gauge freedoms" means

The candidate freedoms acting on the spin bundle over S²:

  - **Single-valued U(1) gauge transformations.** `g(χ,φ) = exp(iΛ(χ,φ))`
    with `Λ(χ, φ+2π) = Λ(χ, φ) (mod 2π)`. These preserve the bundle's
    transition functions on all overlaps and leave the wavefunction
    single-valued.
  - **Multi-valued U(1) transformations with consistent transition
    functions** (Wu-Yang patches). Allowed iff the multi-valuedness on
    every overlap is compensated by a corresponding change in the
    bundle structure (chern number unchanged).
  - **SU(2) frame rotations** acting on the spinor as `ψ → U·ψ`
    with `U(χ,φ) ∈ SU(2)`. These don't change the Hopf bundle's U(1)
    connection but can rotate the spinor basis. For the BAM↔Bloch
    relation, however, only a SCALAR phase (proportional to identity)
    is needed; no genuine SU(2) frame rotation is involved.

The BAM↔Bloch transformation sits in the SECOND class: `g = e^{−iφ/2}`
is a multi-valued scalar U(1) gauge function with `g(φ+2π) = −g(φ)`.
It is allowed for spinor wavefunctions because the spin-½ '
representation is itself a representation of `SU(2) = Spin(3)`, in
which the `−1` ambiguity on a 2π rotation is part of the structure.
It is NOT allowed for charge-1 wavefunctions (vectors), where it
would generate a genuine bundle change.

## Falsifiable predictions

### (1) Connection-level gauge transformation is exact

The line-integral identity

    A_BAM(χ,φ)  =  A_Bloch(χ,φ)  +  ½

(per `dφ`) must hold pointwise on the base. Numerical check at all
points along a sampled trajectory; failure would indicate a
mis-statement of the gauge relation.

### (2) Spinor relation is a scalar U(1) phase

`ψ_BAM(χ,φ) = g(χ,φ) · ψ_Bloch(χ,φ)` must hold pointwise, with the
SCALAR (proportional to identity) multi-valued phase
`g(χ,φ) = e^{−iφ/2}`. Failure ⇒ either spinor parametrisation is
wrong, or the relation is not in fact a pure U(1) gauge transformation
(would indicate an SU(2) frame-rotation component that needs separate
analysis).

### (3) The throat transport conjugates consistently

If we postulate that the throat lives at the equator `χ = π/2` and is
entered at φ_in, exited at φ_out (both on the equatorial S¹), then by
scalar gauge conjugation

    T_Bloch(φ_in, φ_out)
        =  g(π/2, φ_out)^{−1} · T_BAM · g(π/2, φ_in)
        =  e^{+iφ_out/2} · iσ_y · e^{−iφ_in/2}
        =  e^{i(φ_out − φ_in)/2} · iσ_y

A scalar phase times iσ_y. Tr is preserved universally (= 0); for the
instantaneous-throat case (φ_in = φ_out), det = 1 and T² = −I are
preserved as well. For a separated-endpoints throat, the phase
e^{i(φ_out − φ_in)/2} is the gauge-transition cost of the spinor
relocation — the same physics that BAM bakes into the constant `iσ_y`
by using a different section.

### (4) The full-loop holonomy operator differs only by the spinor sign

For any closed loop γ with net φ-winding n (in units of 2π), the
full evolution operator must satisfy

    U_full^{BAM}(γ)  =  (−1)^n · U_full^{Bloch}(γ)

— a pure scalar sign, the value of the multi-valued gauge function
`g(end)/g(start) = e^{−iπ·n}` around the loop. The two operators are
SO(3)-equivalent (same effective rotation angle θ_loop) but differ in
their SU(2) lift by `(−1)^n`. The SO(3)-invariant

    |Tr U_full(γ)|²  =  4·cos²(θ_loop/2)

must agree exactly between BAM and Bloch.

### (5) Relative observables are gauge-invariant

For a Bell-pair-type observable like `phase_spin(θ_a, θ_b)` in
`bell/hopf_phases.py`, the BAM and Bloch formulas

    γ_BAM(θ) = π·cos(θ),       Δγ_BAM = π(cos θ_a − cos θ_b)
    γ_Bloch(θ) = −π(1 − cos θ), Δγ_Bloch = π(cos θ_a − cos θ_b)

give IDENTICAL relative phases. So any relative measurement (CHSH,
two-path interferometry on the same bundle) cannot distinguish the
two gauges. The π/loop offset cancels in every relative observable.

### (6) Distinguishing experiments — what would falsify gauge-equivalence

The only observables that distinguish the two gauges are:

  (a) Absolute closed-loop phase, mod 2π — but this is itself
      gauge-ambiguous (depends on the global section choice), so
      isn't really observable except indirectly through (b).
  (b) The Wu-Yang transition function on patch overlaps. In a
      hypothetical experiment that compared a wavefunction transported
      across the BAM-string region (equator) to one that avoided it,
      the relative phase would be the integral of the transition
      function over the avoided region. For a spinor wavefunction this
      is `−1`, which is observable only via interference with a
      reference that knows the section choice.

The probe enumerates these and computes, for each, whether BAM and
Bloch agree. The expected outcome under standard gauge-bundle
mathematics is that they agree on all (a)–(b) observables: the two
gauges describe the same SU(2) bundle, related by a spinor-frame
rotation, and all physical observables (which are SU(2)/Z₂-invariant
because we measure projection probabilities) coincide.

## Stopping condition

The thread closes when:

  (a) **Gauge equivalence verified.** All six tests above PASS:
      the BAM and Bloch combinations make identical predictions for
      every observable, with the differences absorbed by the spinor
      frame conjugation. Conclusion: BAM's "symmetric gauge + iσ_y
      throat transport" is not a physically distinct theory from the
      standard Bloch+Wu-Yang treatment — it is a particularly
      symmetric SECTION of the same SU(2) Hopf bundle. The
      derivational economy (T = iσ_y falls out directly from the
      orientation-reversing isometry of S³) is BAM's only advantage
      here, and it is real but not falsifiable.

  (b) **Gauge distinguishability found.** Some observable differs
      between BAM and Bloch under all single-valued U(1) gauge
      transformations AND all SU(2) frame rotations. This would
      indicate BAM is making a sharper prediction than Bloch, and
      that prediction can in principle be tested experimentally.
      Outcome (b) is a much stronger claim and would warrant a
      dedicated follow-up probe on the specific distinguishing
      experiment.

## What this probe is NOT

This is a falsifier for **gauge distinguishability**, not for the
throat transport itself. The Berry-phase probe already verified the
internal consistency of `(A_BAM, T = iσ_y, ψ_BAM)`. This probe asks
whether the same physics, expressed in the Bloch gauge with the
gauge-conjugated `T_Bloch`, predicts something different.

It is also not a CHSH / Bell test — the existing `bell/hopf_phases.py`
already uses the BAM gauge and reproduces CHSH = 2√2. The point here
is that the Bloch gauge with consistent T_Bloch should give
**identical** CHSH = 2√2, because the relative phase is the only
input and the gauge offset cancels.

## Cross-references

- `experiments/closure_ledger/moving_mouth_berry_phase_probe.py` —
  the predecessor probe; established BAM-gauge internal consistency.
- `geometrodynamics/embedding/transport.py` — `T = iσ_y` derivation
  from the orientation-reversing Hopf-preserving S³ isometry. This
  derivation is BAM-frame; the probe constructs the Bloch-frame
  analogue.
- `geometrodynamics/hopf/connection.py` — BAM connection
  `A_φ = ½cos(χ)` and full-fibre holonomy `π·cos(χ)`.
- `geometrodynamics/bell/hopf_phases.py` — uses BAM-derived phase
  formula in the Bell analysis; this probe verifies the Bloch-gauge
  equivalent gives the same CHSH.
- `experiments/closure_ledger/throat_transport_gauge_probe.py` — the
  first probe in this thread.
