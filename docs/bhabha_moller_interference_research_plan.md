# Bhabha / Møller two-channel interference — research plan

Next-target probe following the closed Compton/F² derivation thread
(PRs #35 → #41). The Compton kernel and its crossings (Compton, BW,
annihilation) handle all 2-fermion-2-photon tree-level QED diagrams.
Bhabha (`e⁺e⁻ → e⁺e⁻`) and Møller (`e⁻e⁻ → e⁻e⁻`) are the natural next
step: 4-fermion external content, two crossed channels per process,
and interference between them with a Fermi-statistics-determined sign.

## The question

Can the BAM tree kernel reproduce tree-level QED scalar-intensity
structure when **two crossed kernels must be coherently summed with
the correct relative sign / phase**?

The Compton/BW/ann triangle covers a single class of QED diagrams
(2f + 2γ). Bhabha and Møller are 4f processes that sit OUTSIDE
this triangle — they don't differ from Compton by simple Mandelstam
crossing because the external leg content is different. So the BAM
F² closed form does not give them directly; we need either an
extension or a composite construction.

The honest expectation: the BAM scalar-intensity kernel **cannot**
reproduce the Fermi-statistics interference signs natively, because
those signs come from Wick contractions of fermion operators
(Bhabha s ↔ t exchange) or Pauli antisymmetrisation of identical
fermions (Møller t ↔ u exchange) — both of which require complex
amplitudes with explicit Dirac structure, not real positive scalars.

The probe makes this falsifiable: it tests several natural BAM
constructions, characterises exactly which residuals survive, and
identifies the minimal additional structure that BAM would need.

## QED reference

### Bhabha (massless e)

```
|M̄|²_Bhabha / (8e⁴)  =  (s² + u²)/t²  +  (u² + t²)/s²  +  2·u²/(s·t)
                         t-channel²       s-channel²     interference
```

In CM: `s = 4E²`, `t = -s(1-cosθ)/2`, `u = -s(1+cosθ)/2`. The
interference `2u²/(s·t)` is **negative** (since `t < 0`, `s > 0`).

### Møller (massless e)

```
|M̄|²_Møller / (8e⁴)  =  (s² + u²)/t²  +  (s² + t²)/u²  +  2·s²/(t·u)
                         t-channel²       u-channel²     interference
```

The interference `2s²/(t·u)` is **positive** (since `t < 0`, `u < 0`).

## BAM constructions tested

### Construction I: Naïve √-channel positive sum

Postulate `M_channel = √(|M_channel_diag|²)` real positive for each
diagonal channel of QED, summed coherently with positive sign:

```
M_total = M_s + M_t
|M_total|² = |M_s|² + |M_t|² + 2·M_s·M_t  (positive interference)
```

Falsifies on **sign**: for Bhabha, the actual interference is negative.
Falsifies on **magnitude** too: `2·M_s·M_t ≠ 2u²/(s·t)` in general.

### Construction II: Fermi-statistics signs added by hand

Same as I but with explicit relative sign:
  - Bhabha: `M_total = M_s − M_t` (or `M_t − M_s`)
  - Møller: `M_total = M_t − M_u`

Now the sign is right (by construction). Test whether the magnitudes
match: still no, because `M_s · M_t` from the √-channel ansatz
does not equal `u²/(s·t)`.

### Construction III: BAM Compton-F² inspired vertex

Each channel modelled as `M_channel = K(x_q)·√Q(x_q, c_q) / q²` with
appropriate kinematic identification of `x_q`, `c_q`. Test whether
this reproduces the QED diagonal terms.

If even the diagonals fail, that's the most basic mismatch:
Bhabha/Møller don't sit in the BAM Compton kernel's analytic
continuation class.

## Tests

  T1. **QED Bhabha scalar intensity**: compute `|M̄|²_Bhabha / (8e⁴)`
      at sample CM angles; tabulate diagonal vs interference parts.

  T2. **QED Møller scalar intensity**: same for Møller.

  T3. **Interference sign**: confirm Bhabha interference is negative,
      Møller is positive across the angular range.

  T4. **Construction I (positive-sum)**: compute `|M_s + M_t|²` with
      `M_channel = √(QED diagonal)`. Compare to QED `|M|²`. Identify
      sign + magnitude failures.

  T5. **Construction II (Fermi sign)**: compute `|M_s ± M_t|²` with
      the correct relative sign by hand. Test if the magnitude now
      matches.

  T6. **Construction III (Compton-F²-vertex)**: try
      `M_channel = K(x_q)·√Q(x_q, c_q) / q²` with `x_q = 1` (massless,
      no recoil within a single vertex) and other natural choices.
      Test diagonal-term match.

  T7. **Residual analysis**: identify exactly which QED structure
      cannot be produced by any real-positive scalar BAM construction.

  T8. **Required BAM extension**: name the missing ingredient — a
      complex amplitude with Dirac-trace phase, or an explicit
      antisymmetrisation rule from non-orientable throat topology, or
      a virtual-photon propagator with non-trivial Hopf-fibre phase.

## Verdict structure

  - **BHABHA_MOLLER_DERIVED**: Construction I, II, or III (with
    natural ingredients) reproduces both QED `|M̄|²_Bhabha` and
    `|M̄|²_Møller` (including signs) to machine precision. Would mean
    BAM is interference-correct natively. **Not expected.**

  - **DIAGONAL_OK_INTERFERENCE_FAILS**: BAM constructions reproduce
    the diagonal channel terms but the interference sign / magnitude
    doesn't match — confirming the Fermi-statistics gap. The most
    likely outcome.

  - **DIAGONAL_FAILS_TOO**: even the diagonal channel terms don't
    match — Bhabha/Møller live outside the Compton-F² analytic
    continuation class entirely, requiring a fundamentally new BAM
    kernel.

## What this leaves open

  - **Dirac-trace structure in BAM**: capturing the QED γ-matrix
    algebra that produces the `(s² + u²)`, `(u² + t²)` numerators
    requires either explicit Dirac spinors at vertices or a clever
    geometric proxy. Open.

  - **Non-orientable throat → Pauli signs**: BAM's `T = iσ_y`
    antipodal transport (see README, channel 2) is the geometric
    origin of spin-½ and could in principle drive Pauli signs in
    multi-fermion processes. Whether it gives the *correct* signs
    for Møller `(t ↔ u)` interference is a separate (possibly
    deeper) probe target.

  - **Bhabha s-channel virtual-photon propagator**: in BAM, an
    internal photon propagator is presumably a throat-fibre
    connecting two pinches. Its analog of `1/q²` and any phase
    factors is undefined in the current BAM construction.

## Cross-references

  - PR #35 / PR #38–#41: Compton/F² derivation thread.
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y` antipodal
    transport (non-orientable throat, candidate Pauli-sign origin).
  - `experiments/closure_ledger/bhabha_moller_interference_probe.py`
    — this probe.
